///////////////////////////////////////////////////////////////////////////////
//
// File MMFNeuralEP.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// License for the specific language governing rights and limitations under
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: MMFNeuralEP.
//
///////////////////////////////////////////////////////////////////////////////
#include <iomanip>
#include <iostream>

#include <boost/algorithm/string.hpp>
#include <boost/core/ignore_unused.hpp>
#include <boost/math/special_functions/spherical_harmonic.hpp>

#include <DiffusionSolver/EquationSystems/MMFNeuralEP.h>

#include <CardiacEPSolver/Filters/FilterCellHistoryPoints.h>
#include <CardiacEPSolver/Filters/FilterCheckpointCellModel.h>

#include <SolverUtils/Driver.h>
#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>

#include <LibUtilities/BasicUtils/Timer.h>
#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>
#include <LibUtilities/BasicUtils/SessionReader.h>


using namespace std;
using namespace Nektar::SolverUtils;
using namespace Nektar;

namespace Nektar
{
string MMFNeuralEP::className =
    SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
        "MMFNeuralEP", MMFNeuralEP::create, "MMFNeuralEP equation.");

MMFNeuralEP::MMFNeuralEP(const LibUtilities::SessionReaderSharedPtr &pSession,
                         const SpatialDomains::MeshGraphSharedPtr &pGraph)
    : UnsteadySystem(pSession, pGraph), MMFSystem(pSession, pGraph)
{
}

void MMFNeuralEP::v_InitObject(bool DeclareFields)
{
    UnsteadySystem::v_InitObject(DeclareFields);

    int nq   = GetTotPoints();

    // Conductance parameters
    m_session->LoadParameter("Chi", m_chi, 28.0);
    m_session->LoadParameter("Cm", m_capMembrane, 0.125);

    // Helmsolver parameter
    m_session->LoadParameter("Helmtau", m_Helmtau, 1.0);
    
    // Resting potential
    m_session->LoadParameter("urest", m_urest, 0.0);

    // NeuralEP paramter on temperature
    m_session->LoadParameter("Temperature", m_Temperature, 24.0);
    m_session->LoadParameter("diameter", m_diameter, 0.025);

    m_session->LoadParameter("NumelemperNode", m_numelemperNode, 4);

    m_session->LoadParameter("AnisotropyStrength", m_AnisotropyStrength, 4.0);

    // Define ProblemType
    if (m_session->DefinesSolverInfo("NeuralEPType"))
    {
        std::string NeuralEPTypeStr;
        NeuralEPTypeStr = m_session->GetSolverInfo("NEURALEPTYPE");
        for (int i = 0; i < (int)SIZE_NeuralEPType; ++i)
        {
            if (boost::iequals(NeuralEPTypeMap[i], NeuralEPTypeStr))
            {
                m_NeuralEPType = (NeuralEPType)i;
                break;
            }
        }
    }
    else
    {
        m_NeuralEPType = (NeuralEPType)0;
    }

    // Define SovlerSchemeType
    if (m_session->DefinesSolverInfo("SolverSchemeType"))
    {
        std::string SolverSchemeTypeStr;
        SolverSchemeTypeStr = m_session->GetSolverInfo("SolverSchemeType");
        for (int i = 0; i < (int)SIZE_SolverSchemeType; ++i)
        {
            if (boost::iequals(SolverSchemeTypeMap[i], SolverSchemeTypeStr))
            {
                m_SolverSchemeType = (SolverSchemeType)i;
                break;
            }
        }
    }
    else
    {
        m_SolverSchemeType = (SolverSchemeType)0;
    }

    // TimeMap: Parameters
    m_session->LoadParameter("TimeMapStart", m_TimeMapStart, 0.0);
    m_session->LoadParameter("TimeMapEnd", m_TimeMapEnd, 10000.0);
    if (m_session->DefinesSolverInfo("TimeMapType"))
    {
        std::string TIMEMAPTYPEStr;
        TIMEMAPTYPEStr = m_session->GetSolverInfo("TimeMapType");
        for (int i = 0; i < (int)SIZE_TimeMapType; ++i)
        {
            if (boost::iequals(TimeMapTypeMap[i], TIMEMAPTYPEStr))
            {
                m_TimeMap = (TimeMapType)i;
                break;
            }
        }
    }
    else
    {
        m_TimeMap = (TimeMapType)0;
    }

    switch (m_NeuralEPType)
    {
        case eNeuralEPPT:
        {
            std::string vNeuronModel;
            m_session->LoadSolverInfo("NEURONMODEL", vNeuronModel,
                                      "FrankenHuxley");

            ASSERTL0(vNeuronModel != "", "Neuron Model not specified.");

            m_neuron = GetNeuronModelFactory().CreateInstance(
                vNeuronModel, m_session, m_fields[0]);

            // Ranvier node zone: 0: Myelin, 1: node
            m_nfibers  = 1;
            m_NodeZone = Array<OneD, Array<OneD, int>>(m_nfibers);
            for (int i = 0; i < m_nfibers; i++)
            {
                m_NodeZone[i] = Array<OneD, int>(nq, 1);
            }

            // Stimulus
            m_stimulus = Stimulus::LoadStimuli(m_session, m_fields[0]);
            break;
        }

        case eNeuralEP1D:
        {
            std::string vNeuronModel;
            m_session->LoadSolverInfo("NEURONMODEL", vNeuronModel,
                                      "FrankenHuxley");

            ASSERTL0(vNeuronModel != "", "Neuron Model not specified.");

            m_neuron = GetNeuronModelFactory().CreateInstance(
                vNeuronModel, m_session, m_fields[0]);

            // Node and Myelen elements range
            m_session->LoadParameter("ElemNodeEnd", m_ElemNodeEnd, 4);
            m_session->LoadParameter("ElemMyelenEnd", m_ElemMyelenEnd, 1000);

            // Relative Extracellular resistance: 1 < \beta < 10
            m_session->LoadParameter("ratio_re_ri", m_ratio_re_ri, 1.0);

            // Ranvier node zone: 0: Myelin, 1: node
            m_nfibers  = 1;
            m_NodeZone = Array<OneD, Array<OneD, int>>(m_nfibers);
            for (int i = 0; i < m_nfibers; i++)
            {
                m_NodeZone[i] = Array<OneD, int>(nq);
            }

            m_NodeZone[0] = IndexNodeZone1D(m_fields[0], m_ElemNodeEnd, m_ElemMyelenEnd);

            // Constrct m_NeuralCm: node: 1/Cn, Myelin: 1/Cm
            const NekDouble Rf = m_neuron->GetRecistanceValue();
            const NekDouble Cm = m_neuron->GetCapacitanceValue(0);
            const NekDouble Cn = m_neuron->GetCapacitanceValue(1);

            std::cout << "Cm = " << Cm << ", Cn = " << Cn << ", Rf = " << Rf
                      << ", Cm * Rf = " << Cn * Rf << std::endl;
                      
            int index;
            int cntm = 0, cntn = 0, cnte = 0;

            m_NeuralCm    = Array<OneD, Array<OneD, NekDouble>>(1);
            m_NeuralCm[0] = Array<OneD, NekDouble>(nq);
            for (int i = 0; i < m_fields[0]->GetExpSize(); ++i)
            {
                for (int j = 0; j < m_fields[0]->GetTotPoints(i); ++j)
                {
                    index = m_fields[0]->GetPhys_Offset(i) + j;

                    // Ranvier node zone
                    if (m_NodeZone[0][index] >= 0)

                    {
                        m_NeuralCm[0][index] = 1.0 / Cn;
                        cntm++;
                    }

                    // Myelin zone
                    else if (m_NodeZone[0][index] == -1)
                    {
                        m_NeuralCm[0][index] = 1.0 / Cm;
                        cntn++;
                    }

                    // Extracellular space: \sigma_i = m_ratio_re_ri * \sigma_e
                    else
                    {
                        m_NeuralCm[0][index] = 1.0 / Cm / m_ratio_re_ri;
                        cnte++;
                    }
                }
            }

            std::cout << "v_InitObject: Node = " << cntn
            << ", Myelin = " << cntm << ", extracell = " << cnte
            << std::endl;

            // Stimulus
            m_stimulus = Stimulus::LoadStimuli(m_session, m_fields[0]);
            break;
        }

        case eNeuralHelmTest:
        case eNeuralEP2Dmono:
        case eNeuralEP2Dbi:
        case eNeuralEP2DEmbed:
        {
            std::string vNeuronModel;
            m_session->LoadSolverInfo("NEURONMODEL", vNeuronModel,
                                      "HodgkinHuxley");

            ASSERTL0(vNeuronModel != "", "Neuron Model not specified.");

            m_neuron = GetNeuronModelFactory().CreateInstance(
                vNeuronModel, m_session, m_fields[0]);

            // Node and Myelen elements range
            m_session->LoadParameter("ElemNodeEnd", m_ElemNodeEnd, 0);
            m_session->LoadParameter("ElemMyelenEnd", m_ElemMyelenEnd, 0);

            // Relative Extracellular resistance: 1 < \beta < 10
            m_session->LoadParameter("ratio_re_ri", m_ratio_re_ri, 1.0);

            // Ranvier node zone: 0>: Myelin, -1: node, -2: Extracellular space
            m_NodeZone = Array<OneD, Array<OneD, int>>(1);
            m_NodeZone[0] = IndexNodeZone2D(m_fields[0], m_ElemNodeEnd, m_ElemMyelenEnd);

            // Constrct m_NeuralCm: node: 1/Cn, Myelin: 1/Cm
            const NekDouble Rf = m_neuron->GetRecistanceValue();
            const NekDouble Cm = m_neuron->GetCapacitanceValue(0);
            const NekDouble Cn = m_neuron->GetCapacitanceValue(1);

            std::cout << "Cm = " << Cm << ", Cn = " << Cn << ", Rf = " << Rf
                      << ", Cm * Rf = " << Cn * Rf << std::endl;

            int index;
            int cntm = 0, cntn = 0, cnte = 0;

            m_NeuralCm    = Array<OneD, Array<OneD, NekDouble>>(1);
            m_NeuralCm[0] = Array<OneD, NekDouble>(nq);
            for (int i = 0; i < m_fields[0]->GetExpSize(); ++i)
            {
                for (int j = 0; j < m_fields[0]->GetTotPoints(i); ++j)
                {
                    index = m_fields[0]->GetPhys_Offset(i) + j;

                    // Ranvier node zone
                    if (m_NodeZone[0][index] >= 0)

                    {
                        m_NeuralCm[0][index] = 1.0 / Cn;
                        cntm++;
                    }

                    // Myelin zone
                    else if (m_NodeZone[0][index] == -1)
                    {
                        m_NeuralCm[0][index] = 1.0 / Cm;
                        cntn++;
                    }

                    // Extracellular space: \sigma_i = m_ratio_re_ri * \sigma_e
                    else
                    {
                        m_NeuralCm[0][index] = 1.0 / Cm / m_ratio_re_ri;
                        cnte++;
                    }
                }
            }

            std::cout << "v_InitObject: Node = " << cntn
                      << ", Myelin = " << cntm << ", extracell = " << cnte
                      << std::endl;

            // Stimulus
            m_stimulus = Stimulus::LoadStimuli(m_session, m_fields[0]);
            break;
        }

        default:
            break;
    }

    if (m_session->DefinesSolverInfo("MEDIUMTYPE"))
    {
        std::string MediumTypeStr;
        MediumTypeStr = m_session->GetSolverInfo("MEDIUMTYPE");
        for (int i = 0; i < (int)SIZE_MediumType; ++i)
        {
            if (boost::iequals(MediumTypeMap[i], MediumTypeStr))
            {
                m_MediumType = (MediumType)i;
                break;
            }
        }
    }
    else
    {
        m_MediumType = (MediumType)0;
    }

    // Derive AnisotropyStrength.
    switch (m_NeuralEPType)
    {
        case eNeuralEPPT:
        case eNeuralEP1D:
        {
            Array<OneD, Array<OneD, NekDouble>> AniStrength(m_expdim);
            for (int j = 0; j < m_expdim; ++j)
            {
                AniStrength[j] = Array<OneD, NekDouble>(nq, 1.0);
            }

            // Ratio between myelinated fiber and nodal fiber.
            const NekDouble Cm = m_neuron->GetCapacitanceValue(0);
            const NekDouble Cn = m_neuron->GetCapacitanceValue(1);

            // AnisotropyStrength and AniStrength
            Array<OneD, NekDouble> ones(nq, 1.0);

            // AniStrength: Node = Cn/Cm, Myelined = 1.0
            m_AnisotropyStrength = Cn / Cm;

            if (m_MediumType == eAnisotropy)
            {
                for (int j = 0; j < m_expdim; ++j)
                {
                    Vmath::Smul(nq, Cn, &m_NeuralCm[0][0], 1, &AniStrength[j][0], 1);
                    Vmath::Vsqrt(nq, &AniStrength[j][0], 1, &AniStrength[j][0], 1);
                }
            }

            std::cout << "Max Anistrength = "
                      << Vmath::Vmax(nq, AniStrength[0], 1)
                      << ", Min Anistrength = "
                      << Vmath::Vmin(nq, AniStrength[0], 1) << std::endl;

            MMFSystem::MMFInitObject(AniStrength);
            break;
        }

        case eNeuralHelmTest:
        case eNeuralEP2Dmono:
        {
            // Ratio between myelinated fiber and nodal fiber.
            const NekDouble Cm = m_neuron->GetCapacitanceValue(0);
            const NekDouble Cn = m_neuron->GetCapacitanceValue(1);

            Array<OneD, Array<OneD, NekDouble>> AniStrength(m_expdim);
            for (int j = 0; j < m_expdim; ++j)
            {
                AniStrength[j] = Array<OneD, NekDouble>(nq, 1.0);
            }

            // AniStrength: Node = Cn/Cm, Myelined = 1.0
            m_AnisotropyStrength = Cn / Cm;

            if (m_MediumType == eAnisotropy)
            {
                for (int j = 0; j < m_expdim; ++j)
                {
                    Vmath::Smul(nq, Cn, &m_NeuralCm[0][0], 1, &AniStrength[j][0], 1);
                    Vmath::Vsqrt(nq, &AniStrength[j][0], 1, &AniStrength[j][0],
                                 1);
                }
            }

            std::cout << "Max Anistrength_1  = "
                      << Vmath::Vmax(nq, AniStrength[0], 1)
                      << ", Anistrength_2 = "
                      << Vmath::Vmax(nq, AniStrength[1], 1)
                      << ", Min Anistrength 1 = "
                      << Vmath::Vmin(nq, AniStrength[0], 1)
                      << ", Anistrength 2 = "
                      << Vmath::Vmin(nq, AniStrength[1], 1) << std::endl;

            // Create MMF init object
            MMFSystem::MMFInitObject(AniStrength);
            
            // Set up for phie Poisson solver
            std::string phieMMFdirStr = "LOCAL";
            m_session->LoadSolverInfo("phieMMFDir", phieMMFdirStr, "LOCAL");

            Array<OneD, Array<OneD, NekDouble>> phieAniStrength(m_expdim);
            m_phieNeuralCm = Array<OneD, Array<OneD, NekDouble>>(m_expdim);
            for (int j = 0; j < m_expdim; ++j)
            {
                phieAniStrength[j] = Array<OneD, NekDouble>(nq, 1.0);
            }
            break;
        }

        case eNeuralEP2Dbi:
        case eNeuralEP2DEmbed:
        {
            // Ratio between myelinated fiber and nodal fiber.
            const NekDouble Cm = m_neuron->GetCapacitanceValue(0);
            const NekDouble Cn = m_neuron->GetCapacitanceValue(1);

            Array<OneD, Array<OneD, NekDouble>> AniStrength(m_expdim);
            for (int j = 0; j < m_expdim; ++j)
            {
                AniStrength[j] = Array<OneD, NekDouble>(nq, 1.0);
            }

            // AniStrength: Node = Cn/Cm, Myelined = 1.0
            m_AnisotropyStrength = Cn / Cm;

            if (m_MediumType == eAnisotropy)
            {
                for (int j = 0; j < m_expdim; ++j)
                {
                    Vmath::Smul(nq, Cn, &m_NeuralCm[0][0], 1, &AniStrength[j][0], 1);
                    Vmath::Vsqrt(nq, &AniStrength[j][0], 1, &AniStrength[j][0],
                                 1);
                }
            }

            std::cout << "Max Anistrength_1  = "
                      << Vmath::Vmax(nq, AniStrength[0], 1)
                      << ", Anistrength_2 = "
                      << Vmath::Vmax(nq, AniStrength[1], 1)
                      << ", Min Anistrength 1 = "
                      << Vmath::Vmin(nq, AniStrength[0], 1)
                      << ", Anistrength 2 = "
                      << Vmath::Vmin(nq, AniStrength[1], 1) << std::endl;

            // Create MMF init object
            MMFSystem::MMFInitObject(AniStrength);
            
            // Set up for phie Poisson solver
            std::string phieMMFdirStr = "LOCAL";
            m_session->LoadSolverInfo("phieMMFDir", phieMMFdirStr, "LOCAL");

            Array<OneD, Array<OneD, NekDouble>> phieAniStrength(m_expdim);
            m_phieNeuralCm = Array<OneD, Array<OneD, NekDouble>>(m_expdim);
            for (int j = 0; j < m_expdim; ++j)
            {
                phieAniStrength[j] = Array<OneD, NekDouble>(nq, 1.0);
            }

            // Create UnitMovingFrames
            std::cout << std::endl;
            std::cout << "Constructing phieunitMF "
                        "================================================"
                    << std::endl;
            m_phieMMFdir = FindMMFdir(phieMMFdirStr);
            SetUpMovingFrames(m_phieMMFdir, phieAniStrength,
                            m_unitmovingframes);
            CheckMovingFrames(m_unitmovingframes);

            // Create Phiemovingframes
            std::cout << std::endl;
            std::cout << "Constructing phieMF "
                        "================================================"
                    << std::endl;

            m_phiemovingframes =
                Array<OneD, Array<OneD, NekDouble>>(m_mfdim);
            NekDouble Helmfactor =
                sqrt((1.0 + m_ratio_re_ri) / m_ratio_re_ri);

            for (int j = 0; j < m_expdim; ++j)
            {
                Vmath::Smul(nq, Helmfactor, &phieAniStrength[j][0], 1,
                            &phieAniStrength[j][0], 1);
            }

            SetUpMovingFrames(m_phieMMFdir, phieAniStrength,
                                m_phiemovingframes);
            CheckMovingFrames(m_phiemovingframes);

            break;
        }

        default:
        {
            Array<OneD, Array<OneD, NekDouble>> AniStrength(m_expdim);
            for (int j = 0; j < m_expdim; ++j)
            {
                AniStrength[j] = Array<OneD, NekDouble>(nq, 1.0);
            }
            MMFSystem::MMFInitObject(AniStrength);
            break;
        }
    }

    if (m_explicitDiffusion)
    {
        m_ode.DefineImplicitSolve(&MMFNeuralEP::DoNullSolve, this);
        m_ode.DefineProjection(&MMFNeuralEP::DoOdeProjection, this);
    }

    else
    {
        switch (m_NeuralEPType)
        {
            case eNeuralEPPT:
            case eNeuralEP1D:
            {
                // ComputeVarCoeff1D(m_movingframes, m_varcoeff);
                ComputeVarCoeff2D(m_movingframes, m_varcoeff);
                break;
            }
            
            case eNeuralHelmTest:
            case eNeuralEP2Dmono:
            {
                std::cout << "Compute Varcoeff for phim ====================== " << std::endl;
                ComputeVarCoeff2D(m_movingframes, m_varcoeff);

                break;
            }

            case eNeuralEP2Dbi:
            case eNeuralEP2DEmbed:
            {
                std::cout << "Compute Varcoeff for phim ====================== " << std::endl;
                ComputeVarCoeff2D(m_movingframes, m_varcoeff);

                std::cout << "Compute Varcoeff for phie ====================== " << std::endl;
                ComputeVarCoeff2D(m_phiemovingframes, m_phievarcoeff);
                break;
            }

            default:
                break;
        }

        switch (m_NeuralEPType)
        {
            case eNeuralEPPT:
            case eNeuralEP1D:
            {
                m_ode.DefineImplicitSolve(
                    &MMFNeuralEP::DoImplicitSolveNeuralEP1D, this);
                break;
            }

            case eNeuralHelmTest:
            case eNeuralEP2Dmono:
            {
                m_ode.DefineImplicitSolve(
                    &MMFNeuralEP::DoImplicitSolveNeuralEP2Dmono, this);
                break;
            }

            case eNeuralEP2Dbi:
            {
                m_ode.DefineImplicitSolve(
                    &MMFNeuralEP::DoImplicitSolveNeuralEP2Dbi, this);
                break;
            }

            case eNeuralEP2DEmbed:
            {
                m_ode.DefineImplicitSolve(
                    &MMFNeuralEP::DoImplicitSolveNeuralEP2DEmbed, this);
                break;
            }

            default:
                break;
        }
    }

    switch (m_NeuralEPType)
    {
        case eNeuralEPPT:
        {
            m_ode.DefineOdeRhs(&MMFNeuralEP::DoOdeRhsNeuralEPPT, this);
            break;
        }

        case eNeuralEP1D:
        {
            m_ode.DefineOdeRhs(&MMFNeuralEP::DoOdeRhsNeuralEP1D, this);
            break;
        }

        case eNeuralHelmTest:
        case eNeuralEP2Dmono:
        {
            m_ode.DefineOdeRhs(&MMFNeuralEP::DoOdeRhsNeuralEP2Dmono, this);
            break;
        }

        case eNeuralEP2Dbi:
        {
            m_ode.DefineOdeRhs(&MMFNeuralEP::DoOdeRhsNeuralEP2Dbi, this);
            break;
        }

        case eNeuralEP2DEmbed:
        {
            m_ode.DefineOdeRhs(&MMFNeuralEP::DoOdeRhsNeuralEP2DEmbed, this);
            break;
        }

        default:
            break;
    }

       // Test Helm 2D Solver 
    if(m_NeuralEPType==eNeuralHelmTest)
    {
        int nvar = m_fields.size();

        std::cout << std::endl;

        std::cout << "Initiating eNeuralHelmTest =========================================, nvar = " << nvar << std::endl;

        Array<OneD, int> Oneindex(nq, 1);
        Array<OneD, NekDouble> outputarray(nq);

        Array<OneD, Array<OneD, NekDouble>> inputarray(nvar);
        GetFunction("HelmForcingFunction")->Evaluate(m_session->GetVariables(), inputarray);

        std::cout << "inarray = " << RootMeanSquare(inputarray[0]) << std::endl;

        Array<OneD, Array<OneD, NekDouble>> ExactSoln(nvar);
        GetFunction("HelmExactSolution")->Evaluate(m_session->GetVariables(), ExactSoln);

        std::cout << "ExactSoln = " << RootMeanSquare(ExactSoln[1]) << std::endl;

        SolveHelmholtzDiffusion(Oneindex, inputarray[0], m_unitmovingframes, m_phievarcoeff, outputarray);

        Vmath::Vsub(nq, ExactSoln[1], 1, outputarray, 1, outputarray, 1);
        std::cout << "SolveHelmholtzDiffusion Error = " << RootMeanSquare(outputarray) << std::endl;

        wait_on_enter();
    }

}

/**
 *
 */
MMFNeuralEP::~MMFNeuralEP()
{
}

void MMFNeuralEP::CheckNodeZoneMF(
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    const Array<OneD, const Array<OneD, int>> &NodeZone,
    const Array<OneD, const NekDouble> &inarray)
{
    int nq = GetTotPoints();

    Array<OneD, NekDouble> x0(nq);
    Array<OneD, NekDouble> x1(nq);
    Array<OneD, NekDouble> x2(nq);

    m_fields[0]->GetCoords(x0, x1, x2);

    int i, j, index, npts;;
    NekDouble xp, yp, e1mag, e2mag, inarrayavg;
    NekDouble dx, dy, dist;
    for (i = 0; i < m_fields[0]->GetExpSize(); ++i)
    {
        npts = m_fields[0]->GetTotPoints(i);
        xp = 0.0;
        yp = 0.0;
        e1mag = 0.0;
        e2mag = 0.0;
        inarrayavg = 0.0;
        for (j = 0; j < npts; ++j)
        {
            index = m_fields[0]->GetPhys_Offset(i) + j;
            inarrayavg += inarray[index];
            xp += x0[index];
            yp += x1[index];
            e1mag = e1mag +
                    (movingframes[0][index] * movingframes[0][index] +
                     movingframes[0][nq + index] * movingframes[0][nq + index]);
            e2mag = e2mag +
                    (movingframes[1][index] * movingframes[1][index] +
                     movingframes[1][nq + index] * movingframes[1][nq + index]);
        }

        dx = x0[m_fields[0]->GetPhys_Offset(i)] - x0[m_fields[0]->GetPhys_Offset(i)+npts-1];
        dy = x1[m_fields[0]->GetPhys_Offset(i)] - x1[m_fields[0]->GetPhys_Offset(i)+npts-1];
        dist = sqrt(dx*dx+dy*dy);

        e1mag = sqrt(e1mag / npts);
        e2mag = sqrt(e2mag / npts);
        xp = (xp / npts);
        yp = (yp / npts);
        inarrayavg = inarrayavg/npts;

        std::cout << "Elemid = " << i << ", Nodeid = " << NodeZone[0][index]
                << ", x = " << xp << ", y = " << yp << ", dist = " << dist
                << ", e1mag = " << e1mag << ", e2mag = " << e1mag
                << ", inarray = " << inarrayavg << std::endl;
    }
}

// void MMFNeuralEP::ImportFiberXml(
//     const int nfibers,
//     Array<OneD, MultiRegions::ExpListSharedPtr> &fiberfields,
//     Array<OneD, Array<OneD, Array<OneD, int>>> &fiberindex,
//     Array<OneD, LibUtilities::SessionReaderSharedPtr> &fibersession,
//     Array<OneD, SpatialDomains::MeshGraphSharedPtr> &fibergraph)
// {
//     int argc1D = 2;
//     char * argv1D[2];

//     string cmd = "../MMFNeuralEP";
//     char* tmp;
//     tmp = static_cast<char*>(malloc(cmd.length()+1));

//     strcpy(tmp,cmd.c_str());

//     fibersession = Array<OneD,
//     LibUtilities::SessionReaderSharedPtr>(nfibers); fibergraph = Array<OneD,
//     SpatialDomains::MeshGraphSharedPtr>(nfibers);

//     // Index between 1D fiber and 2D plane
//     int nq = m_fields[0]->GetNpoints();

//     Array<OneD, NekDouble> x0(nq);
//     Array<OneD, NekDouble> x1(nq);
//     Array<OneD, NekDouble> x2(nq);

//     m_fields[0]->GetCoords(x0, x1, x2);

//     NekDouble Tol=1.0e-12;

//     for (int nfib=0; nfib<nfibers; ++nfib)
//     {
//         const char* add = "_fiber";
//         const char* xml = ".xml";
//         char* newadd;
//         newadd = static_cast<char*>(malloc(strlen(add)+10+strlen(xml)));
//         char* nfibstr;
//         nfibstr = static_cast<char*>(malloc(1));
//         strcpy(newadd,"_fiber");
//         sprintf(nfibstr,"%d",nfib);
//         strncat(newadd,nfibstr,5);
//         strncat(newadd,xml,5);

//         char* newname;
//         newname =
//         static_cast<char*>(malloc(m_sessionName.length()+strlen(add)+1));

//         strcpy(newname,m_sessionName.c_str());
//         strcat(newname,newadd);

//         argv1D[0] = tmp;
//         argv1D[1] = newname;

//         for (int j=0; j<argc1D; ++j)
//         {
//             std::cout << "j = " << j << ", argv1D = " << argv1D[j] <<
//             std::endl;
//         }
//         std::cout << std::endl;

//         // Create session reader.
//         fibersession[nfib] =
//         LibUtilities::SessionReader::CreateInstance(argc1D, argv1D);

//         // Create MeshGraph
//         fibergraph[nfib] =
//         SpatialDomains::MeshGraph::Read(fibersession[nfib]);

//         // Create field for ith fiber
//         fiberfields[nfib] =
//         MemoryManager<MultiRegions::DisContField1D>::AllocateSharedPtr(fibersession[nfib],
//         fibergraph[nfib], fibersession[nfib]->GetVariable(0));

//         // Check the locations of the fibers points
//         int fnq = fiberfields[nfib]->GetNpoints();

//         Array<OneD, NekDouble> fx0(fnq);
//         Array<OneD, NekDouble> fx1(fnq);
//         Array<OneD, NekDouble> fx2(fnq);
//         m_fiberfields[0]->GetCoords(fx0, fx1, fx2);

//         fiberindex[nfib] = Array<OneD, Array<OneD, int>>(fnq);

//         int indx;
//         NekDouble fxp, fyp, fzp, px, py, pz, dist;
//         for (int k=0; k<fnq; ++k)
//         {
//             fiberindex[nfib][k] = Array<OneD, int>(10,0);

//             fxp=fx0[k];
//             fyp=fx1[k];
//             fzp=fx2[k];
//             indx=0;

//             for (int i=0;i<nq;++i)
//                 {
//                     px=x0[i];
//                     py=x1[i];
//                     pz=x2[i];

//                     dist = sqrt( (px-fxp)*(px-fxp) + (py-fyp)*(py-fyp) +
//                     (pz-fzp)*(pz-fzp) ); if(dist<Tol)
//                     {
//                         fiberindex[nfib][k][indx] = i;
//                         indx=indx+1;
//                     }
//                 }
//         }
//         // for (int k=0; k<fnq; ++k)
//         // {
//         //     std::cout << "k = " << k << ", pt = ( " << fx0[k] << " , " <<
//         fx1[k] << " , " << fx2[k] << " ) "
//         //     << ", index = ( " << fiberindex[nfib][k][0] << " , " <<
//         fiberindex[nfib][k][1] << " , " << fiberindex[nfib][k][2] << " , "
//         //     << fiberindex[nfib][k][3] << " , " <<  fiberindex[nfib][k][4]
//         << " , " << fiberindex[nfib][k][5] << " , " << fiberindex[nfib][k][6]
//         << " ) " << std::endl;
//         // }
//     }
// }

// Constrcuct Cm vector: 1.0/Cn if node. 1.0/Cm if myeline.
Array<OneD, int> MMFNeuralEP::IndexNodeZone1D(
    const MultiRegions::ExpListSharedPtr &field, const int ElemNodeEnd,
        const int ElemMyelenEnd)
{
    boost::ignore_unused(ElemNodeEnd,ElemMyelenEnd);

    int fnq = field->GetNpoints();

    Array<OneD, NekDouble> x0(fnq);
    Array<OneD, NekDouble> x1(fnq);
    Array<OneD, NekDouble> x2(fnq);

    field->GetCoords(x0, x1, x2);

    int index, npts;
    int Nelem = m_fields[0]->GetExpSize();

    Array<OneD, int> outarray(fnq, 0);
    int cnn=0;
    int cnm=0;
    for (int i = 0; i < Nelem; ++i)
    {
        npts = m_fields[0]->GetTotPoints(i);

        for (int j = 0; j < npts; ++j)
        {
            index = m_fields[0]->GetPhys_Offset(i) + j ;

            if(i <=1)
            {
                outarray[index] = i;
                cnn++;
            }

            else if (i % 2 == 1)
            {
                outarray[index] = i/2+1;
                cnn++;
            }

            else 
            {
                outarray[index] = -1;
                cnm++;
            }
        }
    }

    std::cout << "No.: Node = " << cnn << ", Myelin = " << cnm << std::endl;

    return outarray;
}

// Constrcuct Cm vector: 1.0/Cn if node. 1.0/Cm if myeline.
Array<OneD, int> MMFNeuralEP::IndexNodeZone2D(
    const MultiRegions::ExpListSharedPtr &field, const int ElemNodeEnd,
    const int ElemMyelenEnd)
{
    int fnq = field->GetNpoints();

    int index, npts;
    int Nelem = m_fields[0]->GetExpSize();

    Array<OneD, int> outarray(fnq, 0);
    for (int i = 0; i < Nelem; ++i)
    {
        npts = m_fields[0]->GetTotPoints(i);
        for (int j = 0; j < npts; ++j)
        {
            index = m_fields[0]->GetPhys_Offset(i) + j ;

            // First and last element is all node for easier excitation
            if (i <= ElemNodeEnd)
            {
                outarray[index] = i / m_numelemperNode;
            }

            else if (i <= ElemMyelenEnd)
            {
                outarray[index] = -1;
            }

            else
            {
                outarray[index] = -2;
            }
        }
    }

    return outarray;
}

void MMFNeuralEP::DoOdeProjection(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time)
{
    // Counter variable
    int i;
    int nq    = GetNpoints();
    int nvar = inarray.size();

    // Set the boundary conditions
    SetBoundaryConditions(time);

    switch (m_projectionType)
    {
        case MultiRegions::eDiscontinuous:
        {
            // Just copy over array
            for (i = 0; i < nvar; ++i)
            {
                Vmath::Vcopy(nq, inarray[i], 1, outarray[i], 1);
            }
            break;
        }
        case MultiRegions::eGalerkin:
        case MultiRegions::eMixed_CG_Discontinuous:
        {
            Array<OneD, NekDouble> coeffs(m_fields[0]->GetNcoeffs());

            for (i = 0; i < nvar; ++i)
            {
                m_fields[i]->FwdTrans(inarray[i], coeffs);
                m_fields[i]->BwdTrans(coeffs, outarray[i]);
            }
            break;
        }
        default:
        {
            ASSERTL0(false, "Unknown projection scheme");
            break;
        }
    }
}

void MMFNeuralEP::v_DoSolve()
{
    switch (m_SolverSchemeType)
    {
        case eMMFFirst:
        {
            DoSolveMMFFirst();
            break;
        }

        case ePointWise:
        {
            DoSolvePoint();
            break;
        }

        default:
        {
            DoSolveMMFZero();
            break;
        }
    }
}

void MMFNeuralEP::DoSolveMMFZero()
{
    ASSERTL0(m_intScheme != 0, "No time integration scheme.");

    int i, nchk = 1;
    int nq               = GetTotPoints();
    int ncoeffs          = GetNcoeffs();
    int nvariables       = 0;
    int nfields          = m_fields.size();
    std::string fulltext = ""; // initiate fulltext

    if (m_intVariables.empty())
    {
        for (i = 0; i < nfields; ++i)
        {
            m_intVariables.push_back(i);
        }
        nvariables = nfields;
    }
    else
    {
        nvariables = m_intVariables.size();
    }

    // Set up wrapper to fields data storage.
    Array<OneD, Array<OneD, NekDouble>> fields(nvariables);
    Array<OneD, Array<OneD, NekDouble>> fields_old(nvariables);

    // Order storage to list time-integrated fields first.
    for (i = 0; i < 1; ++i)
    {
        fields[i] = m_fields[m_intVariables[i]]->GetPhys();
        m_fields[m_intVariables[i]]->SetPhysState(false);

        fields_old[i] = Array<OneD, NekDouble>(nq, 0.0);
    }

    // Initialise time integration scheme
    m_intScheme->InitializeScheme(m_timestep, fields, m_time, m_ode);

    // Check uniqueness of checkpoint output
    ASSERTL0((m_checktime == 0.0 && m_checksteps == 0) ||
                 (m_checktime > 0.0 && m_checksteps == 0) ||
                 (m_checktime == 0.0 && m_checksteps > 0),
             "Only one of IO_CheckTime and IO_CheckSteps "
             "should be set!");

    LibUtilities::Timer timer;
    bool doCheckTime  = false;
    int step          = 0;
    NekDouble intTime = 0.0;
    NekDouble cpuTime = 0.0;
    NekDouble elapsed = 0.0;

    Array<OneD, NekDouble> tmpc(ncoeffs);

    Array<OneD, NekDouble> velmag(nq, 0.0);
    Array<OneD, NekDouble> velocity(m_spacedim * nq);

    Array<OneD, NekDouble> x0(nq);
    Array<OneD, NekDouble> x1(nq);
    Array<OneD, NekDouble> x2(nq);

    m_fields[0]->GetCoords(x0, x1, x2);

    // Aligh Moving Frames along the velocit vector
    Array<OneD, Array<OneD, NekDouble>> MF1st(m_spacedim);
    for (int i = 0; i < m_spacedim; ++i)
    {
        MF1st[i] = Array<OneD, NekDouble>(m_spacedim * nq);
        Vmath::Smul(m_spacedim * nq, 1.0, &m_movingframes[i][0], 1,
                    &MF1st[i][0], 1);
    }

    Array<OneD, int> phimhistory(nq, 0.0);
    while (step < m_steps || m_time < m_fintime - NekConstants::kNekZeroTol)
    {
        timer.Start();
        fields = m_intScheme->TimeIntegrate(step, m_timestep, m_ode);
        timer.Stop();

        m_time += m_timestep;
        elapsed = timer.TimePerTest(1);
        intTime += elapsed;
        cpuTime += elapsed;

        if (m_session->GetComm()->GetRank() == 0 && !((step + 1) % m_infosteps))
        {
            // Print out at every info step
            std::cout << "Steps: " << std::setw(8) << std::left << step + 1
                      << " "
                      << "Time: " << std::setw(12) << std::left << m_time
                      << std::endl;

            std::stringstream ss;
            ss << cpuTime / 60.0 << " min.";
            std::cout << " CPU Time: " << std::setw(8) << std::left << ss.str()
                      << std::endl << std::endl;

            // fulltext.append("\n");
            // fulltext.append("Steps: " + std::to_string((step+1)));
            // fulltext.append("\n");
            // fulltext.append("Time: " + std::to_string(m_time));
            // fulltext.append("\n");

            // fulltext.append("CPU Time: " + std::to_string(cpuTime / 60.0) + " min.");
            // fulltext.append("\n");

            cpuTime = 0.0;
        }

        // Write out checkpoint files
        if ((m_checksteps && step && !((step + 1) % m_checksteps)) ||
            doCheckTime)
        {
            int Iummax = Vmath::Iamax(nq, fields[0], 1);
            std::cout << "um_max = " << Vmath::Vamax(nq, fields[0], 1)
                      << " at x = " << x0[Iummax] << ", y = " << x1[Iummax] << ", z = " << x2[Iummax]
                      << std::endl;

            if( (m_NeuralEPType==eNeuralEP2Dbi) || (m_NeuralEPType==eNeuralEP2DEmbed))
            {
                int Iuemax = Vmath::Iamax(nq, m_fields[1]->GetPhys(), 1);
                std::cout << "ue_max = " << Vmath::Vamax(nq, m_fields[1]->GetPhys(), 1)
                        << " at x = " << x0[Iuemax] << ", y = " << x1[Iuemax] << ", z = " << x2[Iuemax]
                        << std::endl;
            }

            // NekDouble umax = Vmath::Vamax(nq, fields[0], 1);
            // fulltext.append("u_max = " + std::to_string(umax));
            // fulltext.append(", x = " + std::to_string(x0[Iumax]));
            // fulltext.append(", y = " + std::to_string(x1[Iumax]));
            // fulltext.append("\n");

            // CheckNodeZoneMF(m_movingframes, m_NodeZone, fields[0]);

            // Print phim and phie at each node
            if( m_expdim>1 )
            {
                DisplayatNode();
            }
            
            Checkpoint_Output(nchk++);
            doCheckTime = false;
        }

        ++step;
    } // namespace Nektar

    // Print out summary statistics
    if (m_session->GetComm()->GetRank() == 0)
    {
        std::cout << "Time-integration  : " << intTime << "s" << std::endl;
    }

    for (i = 0; i < 1; ++i)
    {
        m_fields[m_intVariables[i]]->SetPhys(fields[i]);
        m_fields[m_intVariables[i]]->SetPhysState(true);
    }

    for (i = 0; i < nvariables; ++i)
    {
        m_fields[i]->FwdTrans(m_fields[i]->GetPhys(),
                              m_fields[i]->UpdateCoeffs());
    }
} // namespace Nektar


void MMFNeuralEP::DoSolvePoint()
{
    ASSERTL0(m_intScheme != 0, "No time integration scheme.");

    int i, nchk = 1;
    int nq               = GetTotPoints();
    int nvariables       = 0;
    int nfields          = m_fields.size();
    std::string fulltext = ""; // initiate fulltext

    if (m_intVariables.empty())
    {
        for (i = 0; i < nfields; ++i)
        {
            m_intVariables.push_back(i);
        }
        nvariables = nfields;
    }
    else
    {
        nvariables = m_intVariables.size();
    }

    // Set up wrapper to fields data storage.
    Array<OneD, Array<OneD, NekDouble>> fields(nvariables);

    // Order storage to list time-integrated fields first.
    for (i = 0; i < 1; ++i)
    {
        fields[i] = m_fields[m_intVariables[i]]->GetPhys();
        m_fields[m_intVariables[i]]->SetPhysState(false);
    }

    // Initialise time integration scheme
    m_intScheme->InitializeScheme(m_timestep, fields, m_time, m_ode);

    // Check uniqueness of checkpoint output
    ASSERTL0((m_checktime == 0.0 && m_checksteps == 0) ||
                 (m_checktime > 0.0 && m_checksteps == 0) ||
                 (m_checktime == 0.0 && m_checksteps > 0),
             "Only one of IO_CheckTime and IO_CheckSteps "
             "should be set!");

    LibUtilities::Timer timer;
    bool doCheckTime  = false;
    int step          = 0;
    NekDouble intTime = 0.0;
    NekDouble cpuTime = 0.0;
    NekDouble elapsed = 0.0;

    Array<OneD, NekDouble> x0(nq);
    Array<OneD, NekDouble> x1(nq);
    Array<OneD, NekDouble> x2(nq);

    m_fields[0]->GetCoords(x0, x1, x2);

    int totsteps = (m_steps + 1) / m_checksteps;

    Array<OneD, NekDouble> timevec(totsteps);
    Array<OneD, NekDouble> fieldu(totsteps);
    Array<OneD, NekDouble> fieldm(totsteps);
    Array<OneD, NekDouble> fieldn(totsteps);
    Array<OneD, NekDouble> fieldh(totsteps);
    Array<OneD, NekDouble> fieldp(totsteps);

    while (step < m_steps || m_time < m_fintime - NekConstants::kNekZeroTol)
    {
        timer.Start();
        fields = m_intScheme->TimeIntegrate(step, m_timestep, m_ode);
        timer.Stop();

        m_time += m_timestep;
        elapsed = timer.TimePerTest(1);
        intTime += elapsed;
        cpuTime += elapsed;

        if (m_session->GetComm()->GetRank() == 0 && !((step + 1) % m_infosteps))
        {
            // Print out at every info step
            std::cout << "Steps: " << std::setw(8) << std::left << step + 1
                      << " "
                      << "Time: " << std::setw(12) << std::left << m_time
                      << std::endl;

            std::stringstream ss;
            ss << cpuTime / 60.0 << " min.";
            std::cout << " CPU Time: " << std::setw(8) << std::left << ss.str()
                      << std::endl << std::endl;

            cpuTime = 0.0;
        }

        // Write out checkpoint files
        if ((m_checksteps && step && !((step + 1) % m_checksteps)) ||
            doCheckTime)
        {
            std::cout << "time = " << m_time << ", y = " << x1[0] << ", u = " << fields[0][0] 
            << ", m = " << (m_fields[1]->GetPhys())[0] << ", n = " << (m_fields[2]->GetPhys())[0] 
            << ", h = " << (m_fields[3]->GetPhys())[0] << ", p = " << (m_fields[4]->GetPhys())[0] << std::endl;
            
            timevec[nchk] = m_time;
            fieldu[nchk] = fields[0][0];
            fieldm[nchk] = (m_fields[1]->GetPhys())[0];
            fieldn[nchk] = (m_fields[2]->GetPhys())[0];
            fieldh[nchk] = (m_fields[3]->GetPhys())[0];
            fieldp[nchk] = (m_fields[4]->GetPhys())[0];

            Checkpoint_Output(nchk++);
            doCheckTime = false;
        }

        ++step;
    } // namespace Nektar

    // Print out summary statistics
    if (m_session->GetComm()->GetRank() == 0)
    {
        std::cout << "Time-integration  : " << intTime << "s" << std::endl;
    }

    for (i = 0; i < 1; ++i)
    {
        m_fields[m_intVariables[i]]->SetPhys(fields[i]);
        m_fields[m_intVariables[i]]->SetPhysState(true);
    }

    // Output 

    std::cout << "time: ===============================" << std::endl;
    for (i=0;i<totsteps; ++i)
    {
        std::cout << timevec[i] << " , ";
    }
    std::cout << std::endl << std::endl;

    std::cout << "fieldu: ===============================" << std::endl;
    for (i=0;i<totsteps; ++i)
    {
        std::cout << fieldu[i] << " , ";
    }
    std::cout << std::endl << std::endl;

    std::cout << "fieldm: ===============================" << std::endl;
    for (i=0;i<totsteps; ++i)
    {
        std::cout << fieldm[i] << " , ";
    }
    std::cout << std::endl << std::endl;

    std::cout << "fieldn: ===============================" << std::endl;
    for (i=0;i<totsteps; ++i)
    {
        std::cout << fieldn[i] << " , ";
    }
    std::cout << std::endl << std::endl;

    std::cout << "fieldh: ===============================" << std::endl;
    for (i=0;i<totsteps; ++i)
    {
        std::cout << fieldh[i] << " , ";
    }
    std::cout << std::endl << std::endl;

    std::cout << "fieldp: ===============================" << std::endl;
    for (i=0;i<totsteps; ++i)
    {
        std::cout << fieldp[i] << " , ";
    }
    std::cout << std::endl << std::endl;

    for (i = 0; i < nvariables; ++i)
    {
        m_fields[i]->FwdTrans(m_fields[i]->GetPhys(),
                              m_fields[i]->UpdateCoeffs());
    }
} // namespace Nektar



void MMFNeuralEP::DisplayatNode()
{
    int nvar = m_fields.size();

    if(nvar==1)
    {
        DisplayatNodevar1();
    }

    else if(nvar==2)
    {
        DisplayatNodevar2();
    }
}

void MMFNeuralEP::DisplayatNodevar1()
{
    // Print phim and phie at each node
    int index, Rnodeid = 0;
    NekDouble locphimsum;

    Array<OneD, NekDouble> phimavg(m_ElemNodeEnd);
    for (int i = 0; i < m_ElemNodeEnd; ++i)
    {
        Rnodeid = i / m_numelemperNode;

        if( (i-Rnodeid*m_numelemperNode)==0 )
        {
            locphimsum = 0.0;

            for (int j = 0; j < m_fields[0]->GetTotPoints(i); ++j)
            {
                index =  m_fields[0]->GetPhys_Offset(i) + j;
                locphimsum = locphimsum + m_fields[0]->GetPhys()[index];
            }

            phimavg[Rnodeid] = locphimsum / m_fields[0]->GetTotPoints(i);
        }
    }

    std::cout << " " << std::endl;
    std::cout << "(Nodeid,phim): ";
    for (int i = 0; i < Rnodeid + 1; ++i)
    {
        std::cout << "(" << i << "," << phimavg[i] << "), ";
    }
    std::cout << " " << std::endl << std::endl;
}

void MMFNeuralEP::DisplayatNodevar2()
{
    int nq = GetTotPoints();

    // Print phim and phie at each node
    int index, Rnodeid = 0;
    NekDouble locphimsum, locphiesum;

    Array<OneD, NekDouble> phimavg(m_ElemNodeEnd);
    Array<OneD, NekDouble> phieavg(m_ElemNodeEnd);
    for (int i = 0; i < m_ElemNodeEnd; ++i)
    {
        Rnodeid = i / m_numelemperNode;

        locphimsum = 0.0;
        locphiesum = 0.0;
        for (int j = 0; j < m_fields[0]->GetTotPoints(i); ++j)
        {
            index =  m_fields[0]->GetPhys_Offset(i) + j;

            locphimsum = locphimsum + m_fields[0]->GetPhys()[index];
            locphiesum = locphiesum + m_fields[1]->GetPhys()[index];
        }

        phimavg[Rnodeid] = locphimsum / m_fields[0]->GetTotPoints(i);
        phieavg[Rnodeid] = locphiesum / m_fields[1]->GetTotPoints(i);
    }

    std::cout << " " << std::endl;
    std::cout << "(Nodeid,phim,phie): ";
    for (int i = 0; i < Rnodeid + 1; ++i)
    {
        std::cout << "(" << i << "," << phimavg[i] << "," << phieavg[i] << "), ";
    }
    std::cout << " " << std::endl << std::endl;
}


// void MMFNeuralEP::DisplayConductionVelocity(
//     const Array<OneD, const NekDouble> &field,
//     const Array<OneD, const NekDouble> &dudt)
// {
//     int nq         = GetTotPoints();

//     NekDouble phim_th = 100.0;

//     Array<OneD, NekDouble> x0(nq);
//     Array<OneD, NekDouble> x1(nq);
//     Array<OneD, NekDouble> x2(nq);

//     m_fields[0]->GetCoords(x0, x1, x2);

//     for (int i=0;i<nq; ++i)
//     {
//         if( (dudt[i]>0) && (field[i]>phim_th) )
//         {
//             if (phimhistory[i]==0)
//             {
//                 frontloc = x1[i];
//             }

//             phimhistory[i] = 1;
//         }
//     }

//     // Compute conduction velocity
//     condvel = 1000.0 * (frontloc - frontloc_old)/(m_time - time_old);

//     if(condvel >0.00001)
//     {
//     std::cout << "frontloc = " << frontloc << ", cond. vel. = " << condvel <<
//     " m/s" << std::endl;
//     }
// }

// void MMFNeuralEP::Plotphimphie(const Array<OneD, const NekDouble> &phim,
//                              const Array<OneD, const NekDouble> &phie,
//                             const Array<OneD, const NekDouble> &Exactphie,
//                              const int nstep)
// {
//     int nvar    = 4;
//     int ncoeffs = m_fields[0]->GetNcoeffs();
//     int nq         = GetTotPoints();

//     std::string outname1 = m_sessionName + "_" +
//                            boost::lexical_cast<std::string>(nstep) + ".chk";

//     std::vector<Array<OneD, NekDouble>> fieldcoeffs(nvar);
//     for (int i = 0; i < nvar; ++i)
//     {
//         fieldcoeffs[i] = Array<OneD, NekDouble>(ncoeffs);
//     }

//     std::vector<std::string> variables(nvar);
//     variables[0] = "um";
//     variables[1] = "ue";
//     variables[2] = "Exactue";
//     variables[3] = "ue_err";

//     Array<OneD, NekDouble> err(nq);
//     Vmath::Vsub(nq, &phie[0], 1, &Exactphie[0], 1, &err[0], 1);

//     m_fields[0]->FwdTrans(phim, fieldcoeffs[0]);
//     m_fields[0]->FwdTrans(phie, fieldcoeffs[1]);
//     m_fields[0]->FwdTrans(Exactphie, fieldcoeffs[2]);
//     m_fields[0]->FwdTrans(err, fieldcoeffs[3]);

//     WriteFld(outname1, m_fields[0], fieldcoeffs, variables);
// }

void MMFNeuralEP::DoSolveMMFFirst()
{
    ASSERTL0(m_intScheme != 0, "No time integration scheme.");

    int i, nchk = 1;
    int nq         = GetTotPoints();
    int ncoeffs    = GetNcoeffs();
    int nvariables = 0;
    int nfields    = m_fields.size();

    if (m_intVariables.empty())
    {
        for (i = 0; i < nfields; ++i)
        {
            m_intVariables.push_back(i);
        }
        nvariables = nfields;
    }
    else
    {
        nvariables = m_intVariables.size();
    }

    // Set up wrapper to fields data storage.
    Array<OneD, Array<OneD, NekDouble>> fields(nvariables);
    Array<OneD, Array<OneD, NekDouble>> fieldsold(nvariables);

    // Order storage to list time-integrated fields first.
    for (i = 0; i < nvariables; ++i)
    {
        fields[i] = m_fields[m_intVariables[i]]->GetPhys();
        m_fields[m_intVariables[i]]->SetPhysState(false);

        fieldsold[i] = Array<OneD, NekDouble>(nq);
    }

    // Initialise time integration scheme
    m_intScheme->InitializeScheme(m_timestep, fields, m_time, m_ode);

    // Check uniqueness of checkpoint output
    ASSERTL0((m_checktime == 0.0 && m_checksteps == 0) ||
                 (m_checktime > 0.0 && m_checksteps == 0) ||
                 (m_checktime == 0.0 && m_checksteps > 0),
             "Only one of IO_CheckTime and IO_CheckSteps "
             "should be set!");

    LibUtilities::Timer timer;
    bool doCheckTime  = false;
    int step          = 0;
    NekDouble intTime = 0.0;
    NekDouble cpuTime = 0.0;
    NekDouble elapsed = 0.0;

    Array<OneD, NekDouble> tmpc(ncoeffs);

    Array<OneD, NekDouble> velmag(nq, 0.0);
    Array<OneD, NekDouble> velocity(m_spacedim * nq);

    // Aligh Moving Frames along the velocit vector
    Array<OneD, Array<OneD, NekDouble>> MF1st(m_spacedim);
    Array<OneD, Array<OneD, NekDouble>> MF1sttmp(m_spacedim);
    Array<OneD, Array<OneD, NekDouble>> MF1stAligned(m_spacedim);
    Array<OneD, Array<OneD, NekDouble>> TimeMapMF(m_spacedim);
    for (int i = 0; i < m_spacedim; ++i)
    {
        MF1st[i]    = Array<OneD, NekDouble>(m_spacedim * nq);
        MF1sttmp[i] = Array<OneD, NekDouble>(m_spacedim * nq);

        MF1stAligned[i] = Array<OneD, NekDouble>(m_spacedim * nq, 0.0);

        TimeMapMF[i] = Array<OneD, NekDouble>(m_spacedim * nq, 0.0);

        Vmath::Smul(m_spacedim * nq, 1.0, &m_movingframes[i][0], 1,
                    &MF1st[i][0], 1);
    }

    // Connection 1-form
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> MF1stConnection(m_mfdim);
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> TMMFConnection(m_mfdim);

    Array<OneD, Array<OneD, NekDouble>> MF1stCurvature(m_mfdim);
    Array<OneD, Array<OneD, NekDouble>> TMMFCurvature(m_mfdim);
    for (int i = 0; i < m_mfdim; i++)
    {
        MF1stCurvature[i] = Array<OneD, NekDouble>(nq, 0.0);
        TMMFCurvature[i]  = Array<OneD, NekDouble>(nq, 0.0);

        TMMFConnection[i]  = Array<OneD, Array<OneD, NekDouble>>(m_mfdim);
        MF1stConnection[i] = Array<OneD, Array<OneD, NekDouble>>(m_mfdim);
        for (int j = 0; j < m_mfdim; j++)
        {
            TMMFConnection[i][j]  = Array<OneD, NekDouble>(nq, 0.0);
            MF1stConnection[i][j] = Array<OneD, NekDouble>(nq, 0.0);
        }
    }

    Array<OneD, Array<OneD, NekDouble>> Relacc(m_shapedim);
    Array<OneD, Array<OneD, NekDouble>> TMRelacc(m_shapedim);
    for (int j = 0; j < m_shapedim; ++j)
    {
        Relacc[j]   = Array<OneD, NekDouble>(nq, 0.0);
        TMRelacc[j] = Array<OneD, NekDouble>(nq, 0.0);
    }

    Array<OneD, int> ActivatedPre(nq, 0);
    Array<OneD, int> Activated(nq, 0);
    Array<OneD, int> ActivatedHistory(nq, 0);

    Array<OneD, NekDouble> VelmagHistory(nq, 0.0);
    Array<OneD, NekDouble> fieldHistory(nq, 0.0);

    Array<OneD, NekDouble> DivDiff(nq, 0.0);

    Array<OneD, NekDouble> dudtval(nq);
    Array<OneD, NekDouble> dudtvalHistory(nq, 0.0);
    Array<OneD, NekDouble> NoBoundaryZone(nq, 1.0);

    Array<OneD, int> dudt(nq);
    Array<OneD, int> APindex(nq, 1);

    Array<OneD, int> NewValidTimeMap(nq, 1);

    Array<OneD, NekDouble> Laplacian(nq);
    Array<OneD, NekDouble> LaplacianNew(nq);

    Array<OneD, Array<OneD, NekDouble>> qfield(m_expdim);
    Array<OneD, Array<OneD, NekDouble>> qfieldNew(m_expdim);

    Array<OneD, NekDouble> TimeMap(nq, 0.0);
    Array<OneD, NekDouble> IappMap(nq, 0.0);
    Array<OneD, NekDouble> UnitVelMap(m_spacedim * nq, 0.0);

    int totsteps = (m_steps + 1) / m_checksteps;
    Array<OneD, NekDouble> uval(totsteps, 0.0);
    Array<OneD, NekDouble> mval(totsteps, 0.0);
    Array<OneD, NekDouble> nval(totsteps, 0.0);
    Array<OneD, NekDouble> hval(totsteps, 0.0);
    Array<OneD, NekDouble> pval(totsteps, 0.0);

    Array<OneD, NekDouble> fieldoldchk(nq, 0.0);
    Array<OneD, NekDouble> fieldchkdiff(nq, 0.0);
    while (step < m_steps || m_time < m_fintime - NekConstants::kNekZeroTol)
    {
        // Initialize Activated
        Activated = Array<OneD, int>(nq, 0);

        // Save fields into fieldsold
        for (i = 0; i < nvariables; ++i)
        {
            Vmath::Vcopy(nq, &fields[i][0], 1, &fieldsold[i][0], 1);
        }

        // field time integration
        timer.Start();
        fields = m_intScheme->TimeIntegrate(step, m_timestep, m_ode);
        timer.Stop();

        m_time += m_timestep;
        elapsed = timer.TimePerTest(1);
        intTime += elapsed;
        cpuTime += elapsed;

        // Compute TimeMap
        // dudtsign: wavefront = -1.0, waveback = 1.0
        //  dudt = Computedudt(m_uTol, fields[0], fieldsold[0]);
        Vmath::Vsub(nq, fields[0], 1, fieldsold[0], 1, dudtval, 1);
        Vmath::Smul(nq, 1.0 / m_timestep, dudtval, 1, dudtval, 1);

        // Smoothing dudt map for a smooth time map
        // HelmSolveSmoothing(m_TimeMapSmoothL, dudtval);
        // Vmath::Vmul(nq, m_ValidTimeMap, 1, dudtval, 1, dudtval, 1);

        // Compute Proper Time Map by weight integration of field.
        // TMmode = 0 (Gradient-weighted time)
        // Output = Propertimemap: time when the cell is excited.
        //          fieldHistory: sum of field is updated

        if ((m_TimeMapStart <= m_time) && (m_TimeMapEnd >= m_time))
        {
            ComputeTimeMap(m_time, m_urest, fields[0], dudtval, m_ValidTimeMap,
                           dudtvalHistory, IappMap, TimeMap);
        }

        // Aligning moving frames along the velocity vector

        // For multiple waves, if dudt changes from positive to
        // negative, it is a peak to distinguish WB from WF. if dudt
        // changes from negative to a negligible magnitude or positive,
        // then it is another peak to change index of AP.
        // ComputedudtHistory(dudt, fields[0], dudtHistory, APindex);

        // vector = the gradient of u
        velocity = ComputeDirectionVector(m_movingframes, fields[0], dudtval);

        // Compute the magnitude of velocityls
        velmag = ComputeVelocityMag(velocity);

        // Activated = 1 only where u > m_uTol. is rad >
        // m_NoAlignInitRadius.
        ActivatedPre =
            ComputeZoneActivation(m_uTol, fields[0], m_NoAlignInitRadius);

        // Elementwise activation: Activate when velmag is larger than
        // VATol
        m_fields[0]->ElementWiseActivation(-1, velmag, m_VelActivationTol,
                                           ActivatedPre);

        // Align MF to Velocity vector if Activated is on.
        // Input: Activated, velocity, MF1st_old (movingframes)
        // Output: MF1st
        AlignMFtoVelocity(ActivatedPre, velocity, m_movingframes, MF1sttmp);

        // MF1sttmp has the same magnitude with m_movingframes which may
        // have anisotropy

        // Compute the difference of the divergence of the gradient
        WeakDGMMFirstLaplacian(0, m_movingframes, fields[0], Laplacian);
        WeakDGMMFirstLaplacian(0, MF1sttmp, fields[0], LaplacianNew);

        DivDiff = ComputeLaplacianDiff(Laplacian, LaplacianNew);

        // Validate New frames to modify MF1st and Activated
        // Elementwise activation: Activate when DivDiff is smaller than
        // AdaptNewFramesTol
        Vmath::Vcopy(nq, ActivatedPre, 1, Activated, 1);
        m_fields[0]->ElementWiseActivation(1, DivDiff, m_AdaptNewFramesTol,
                                           Activated);

        // Align MF to Velocity vector if Activated is on.
        // Input: Activated, velocity, MF1st_old (movingframes)
        // Output: MF1st
        AlignMFtoVelocity(Activated, velocity, m_movingframes, MF1st);

        // For Weighted integration
        // ================================================== Input:
        // Activated, MF1st Output: ActivatedHistory: 0 or 1. 1 is
        // activated.
        //         MF1stAligned: MF1st is updated and stored
        //         VelmagHistory: sum of velmag is updated and stored
        UpdateMF1st(Activated, MF1st, velmag, VelmagHistory, MF1stAligned,
                    ActivatedHistory);

        if (m_session->GetComm()->GetRank() == 0 && !((step + 1) % m_infosteps))
        {
            std::cout << "Steps: " << std::setw(8) << std::left << step + 1
                      << " "
                      << "Time: " << std::setw(12) << std::left << m_time
                      << std::endl;

            std::cout << "DivDiff = " << Vmath::Vamax(nq, DivDiff, 1)
                      << ", ActivatedHistory = "
                      << CountActivated(ActivatedHistory) << " ( "
                      << 100.0 * CountActivated(ActivatedHistory) / nq
                      << " % ) " << std::endl;

            std::stringstream ss;
            ss << cpuTime / 60.0 << " min.";
            std::cout << " CPU Time: " << std::setw(8) << std::left << ss.str()
                      << std::endl;

            if (CountActivated(ActivatedHistory) > 0)
            {
                // Check the Curvature 2-form of the aligned moving frames
                Compute2DConnectionCurvature(MF1stAligned, MF1stConnection,
                                             MF1stCurvature);

                // Test moving frames Connection whenever it is possible
                Test2DConnectionCurvature(m_Initx, m_Inity, m_Initz,
                                          ActivatedHistory, MF1stAligned,
                                          MF1stConnection, MF1stCurvature);
            }

            cpuTime = 0.0;
        }

        // Write out checkpoint files
        if ((m_checksteps && step && !((step + 1) % m_checksteps)) ||
            doCheckTime)
        {
            // NekDouble dudtpros, dudtneg;
            // dudtpros = Computedudtpercent(1, dudt);
            // dudtneg  = Computedudtpercent(-1, dudt);

            // NekDouble udiff;
            // udiff = Vmath::Vmax(nq, fields[0], 1) - Vmath::Vmin(nq,
            // fields[0], 1);

            Array<OneD, NekDouble> x0(nq);
            Array<OneD, NekDouble> x1(nq);
            Array<OneD, NekDouble> x2(nq);

            m_fields[0]->GetCoords(x0, x1, x2);

            int Iumax = Vmath::Iamax(nq, fields[0], 1);
            std::cout << "u_max= " << Vmath::Vamax(nq, fields[0], 1)
                      << " at x = " << x0[Iumax] << ", y = " << x1[Iumax]
                      << std::endl;

            int Ivelmax = Vmath::Iamax(nq, velmag, 1);
            std::cout << "vel_max= " << Vmath::Vamax(nq, velmag, 1)
                      << " at x = " << x0[Ivelmax] << ", y = " << x1[Ivelmax]
                      << std::endl;

            if (CountActivated(ActivatedHistory) > 0)
            {
                ComputeRelacc(MF1stAligned, Relacc);

                PlotTrajectoryMF(ActivatedHistory, fields[0], MF1stAligned,
                                 MF1stConnection, Relacc, NoBoundaryZone, nchk);

                // Plot relative acceleration and conduction block zone
                // PlotRelacc2D(ActivatedHistory, MF1stConnection,
                // MF1stCurvature, Relacc, RelaccOmega, nchk);
            }

            if (m_TimeMap == eActivated)
            {
                // Compute velocity field
                // VelocityMap = ComputeVelocityField(m_ValidTimeMap, TimeMap);
                ComputeMFTimeMap(m_ValidTimeMap, TimeMap, NewValidTimeMap,
                                 TimeMapMF);

                // Compute2DConnectionCurvature(TimeMapMF, TMMFConnection,
                // TMMFCurvature);

                Compute2DConnection1form(TimeMapMF, TMMFConnection);

                // Compute Relative Acceleration
                ComputeRelacc(TimeMapMF, TMRelacc);

                // PlotTimeMap(TimeMap, IappMap, TimeMapMF, nchk);

                // PlotTimeMapMF(NoBoundaryZone, TimeMap, TimeMapMF, MF1stAligned,
                //               TMMFConnection, TMRelacc, nchk);

                std::cout << "Time Map: Max = " << Vmath::Vmax(nq, TimeMap, 1)
                          << ", Min = " << Vmath::Vmin(nq, TimeMap, 1)
                          << std::endl;
            }

            Checkpoint_Output(nchk++);
            doCheckTime = false;
        }

        // Step advance
        ++step;
    } // namespace Nektar

    // Print out summary statistics
    if (m_session->GetComm()->GetRank() == 0)
    {
        std::cout << "Time-integration  : " << intTime << "s" << std::endl;
    }

    for (i = 0; i < nvariables; ++i)
    {
        m_fields[m_intVariables[i]]->SetPhys(fields[i]);
        m_fields[m_intVariables[i]]->SetPhysState(true);
    }

    for (i = 0; i < nvariables; ++i)
    {
        m_fields[i]->FwdTrans(m_fields[i]->GetPhys(),
                              m_fields[i]->UpdateCoeffs());
    }
} // namespace Nektar

void MMFNeuralEP::PlotFHIonCurrent(const Array<OneD, const NekDouble> &inarray,
                                   const int nstep)
{
    const int nvar = m_neuron->GetNumNeuronVariables();

    int nq      = m_fields[0]->GetTotPoints();
    int ncoeffs = m_fields[0]->GetNcoeffs();

    std::string outname;
    outname = m_sessionName + "_FHIC_" +
              boost::lexical_cast<std::string>(nstep) + ".chk";

    std::vector<Array<OneD, NekDouble>> fieldcoeffs(nvar);
    for (int i = 0; i < nvar; ++i)
    {
        fieldcoeffs[i] = Array<OneD, NekDouble>(ncoeffs);
    }

    std::vector<std::string> variables(nvar);
    variables[0] = "V";
    variables[1] = "m";
    variables[2] = "n";
    variables[3] = "h";
    variables[4] = "p";

    Array<OneD, NekDouble> tmp(nq);
    Vmath::Vcopy(nq, &inarray[0], 1, &tmp[0], 1);
    m_fields[0]->FwdTrans(tmp, fieldcoeffs[0]);

    for (int i = 1; i < nvar; ++i)
    {
        tmp = m_neuron->GetNeuronSolution(i);
        m_fields[0]->FwdTrans(tmp, fieldcoeffs[i]);
    }

    WriteFld(outname, m_fields[0], fieldcoeffs, variables);
}

Array<OneD, NekDouble> MMFNeuralEP::ComputeLaplacianDiff(
    const Array<OneD, const NekDouble> &Laplacian,
    const Array<OneD, const NekDouble> &LaplacianNew)
{
    int nq = GetTotPoints();

    Array<OneD, NekDouble> DivDiff(nq, 0.0);

    NekDouble Tol = 1.0e-10;
    for (int i = 0; i < nq; ++i)
    {
        if (fabs(Laplacian[i]) > Tol)
        {
            DivDiff[i] = fabs((Laplacian[i] - LaplacianNew[i]) / Laplacian[i]);
        }

        else
        {
            DivDiff[i] = fabs(Laplacian[i] - LaplacianNew[i]);
        }
    }

    return DivDiff;
}

// Implicit solve for NeuralEP solver
void MMFNeuralEP::DoImplicitSolveNeuralEP1D(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time,
    const NekDouble lambda)
{
    boost::ignore_unused(time);

    int nvar = m_fields.size();
    int nq = m_fields[0]->GetNpoints();

    const NekDouble R_f = m_neuron->GetRecistanceValue();
    const NekDouble C_n = m_neuron->GetCapacitanceValue(1);

    StdRegions::ConstFactorMap factors;
    factors[StdRegions::eFactorTau] = m_Helmtau;

    factors[StdRegions::eFactorLambda] = C_n * R_f / lambda;
    if(nvar==1)
    {
        NekDouble betaratio = (m_ratio_re_ri + 1.0) / m_ratio_re_ri;
        factors[StdRegions::eFactorLambda] = C_n * R_f * betaratio / lambda;
    }

    // We solve ( \nabla^2 - HHlambda ) Y[i] = rhs [i]
    // inarray = input: \hat{rhs} -> output: \hat{Y}
    // outarray = output: \hat{Y} where \hat = modal coeffs

    SetBoundaryConditions(time);

    // For the variable of membrane potential: Multiply 1.0/timestep
    Vmath::Smul(nq, -factors[StdRegions::eFactorLambda], inarray[0], 1,
                m_fields[0]->UpdatePhys(), 1);

    // Solve a system of equations with Helmholtz solver and transform
    // back into physical space.

    m_fields[0]->HelmSolve(m_fields[0]->GetPhys(), m_fields[0]->UpdateCoeffs(),
                           factors, m_varcoeff);
    m_fields[0]->BwdTrans(m_fields[0]->GetCoeffs(), outarray[0]);
    m_fields[0]->SetPhysState(true);
}

// Implicit solve for NeuralEP 2D solver
void MMFNeuralEP::DoImplicitSolveNeuralEP2Dmono(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time,
    const NekDouble lambda)
{
    boost::ignore_unused(time);

    int nq   = m_fields[0]->GetNpoints();

    // Set up factors for Helmsolve
    const NekDouble R_f = m_neuron->GetRecistanceValue();
    const NekDouble C_n = m_neuron->GetCapacitanceValue(1);

    StdRegions::ConstFactorMap factors;
    factors[StdRegions::eFactorTau] = m_Helmtau;

    // m_ratio_re_ri determines the conductivity only when the variable is one
    NekDouble betaratio = (m_ratio_re_ri + 1.0) / m_ratio_re_ri;
    factors[StdRegions::eFactorLambda] = C_n * R_f * betaratio / lambda;

    SetBoundaryConditions(time);

    // Multiply 1.0/timestep
    Vmath::Smul(nq, -factors[StdRegions::eFactorLambda], inarray[0], 1,
                m_fields[0]->UpdatePhys(), 1);

    m_fields[0]->HelmSolve(m_fields[0]->GetPhys(), m_fields[0]->UpdateCoeffs(),
                           factors, m_varcoeff);
    m_fields[0]->BwdTrans(m_fields[0]->GetCoeffs(), outarray[0]);
    m_fields[0]->SetPhysState(true);
}

// Implicit solve for NeuralEP 2D solver
void MMFNeuralEP::DoImplicitSolveNeuralEP2Dbi(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time,
    const NekDouble lambda)
{
    boost::ignore_unused(time);
    int nq   = m_fields[0]->GetNpoints();

    // Set up factors for Helmsolve
    const NekDouble R_f = m_neuron->GetRecistanceValue();
    const NekDouble C_n = m_neuron->GetCapacitanceValue(1);

    StdRegions::ConstFactorMap factors;
    factors[StdRegions::eFactorTau] = m_Helmtau;

    // m_ratio_re_ri determines the conductivity only when the variable is one
    factors[StdRegions::eFactorLambda] = C_n * R_f / lambda;

    SetBoundaryConditions(time);
    // SetMembraneBoundaryCondition(time);

    // Multiply 1.0/timestep
    Vmath::Smul(nq, -factors[StdRegions::eFactorLambda], inarray[0], 1,
                m_fields[0]->UpdatePhys(), 1);

    m_fields[0]->HelmSolve(m_fields[0]->GetPhys(), m_fields[0]->UpdateCoeffs(),
                           factors, m_varcoeff);
    m_fields[0]->BwdTrans(m_fields[0]->GetCoeffs(), outarray[0]);
    m_fields[0]->SetPhysState(true);
}


// Implicit solve for NeuralEP 2D solver
void MMFNeuralEP::DoImplicitSolveNeuralEP2DEmbed(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time,
    const NekDouble lambda)
{
    boost::ignore_unused(time);

    int nq = m_fields[0]->GetNpoints();

    // Set up factors for Helmsolve
    const NekDouble R_f = m_neuron->GetRecistanceValue();
    const NekDouble C_n = m_neuron->GetCapacitanceValue(1);

    StdRegions::ConstFactorMap factors;
    factors[StdRegions::eFactorTau] = m_Helmtau;

    NekDouble Cv                       = C_n * R_f;
    factors[StdRegions::eFactorLambda] = Cv / lambda;

    // SetBoundaryConditions(time);
    SetMembraneBoundaryCondition();

    // Multiply 1.0/timestep
    Vmath::Smul(nq, -factors[StdRegions::eFactorLambda], inarray[0], 1,
                m_fields[0]->UpdatePhys(), 1);

    m_fields[0]->HelmSolve(m_fields[0]->GetPhys(), m_fields[0]->UpdateCoeffs(),
                           factors, m_varcoeff);

    m_fields[0]->BwdTrans(m_fields[0]->GetCoeffs(), outarray[0]);
    m_fields[0]->SetPhysState(true);
}


//   void MMFNeuralEP::SetAxonWallBoundaryConditions(Array<OneD, Array<OneD,
//   NekDouble> > &inarray)
//   {
//       std::string varName;
//       int nvariables = m_fields.size();
//       int cnt = 0;
//       int nTracePts  = GetTraceTotPoints();

//       // Extract trace for boundaries. Needs to be done on all processors to
//       avoid
//       // deadlock.
//       Array<OneD, Array<OneD, NekDouble> > Fwd(nvariables);
//       for (int i = 0; i < nvariables; ++i)
//       {
//           Fwd[i] = Array<OneD, NekDouble>(nTracePts);
//           m_fields[i]->ExtractTracePhys(inarray[i], Fwd[i]);
//       }

//       // loop over Boundary Regions
//       for(int n = 0; n < m_fields[0]->GetBndConditions().size(); ++n)
//       {
//           // Wall Boundary Condition
//           if
//           (boost::iequals(m_fields[0]->GetBndConditions()[n]->GetUserDefined(),"AxonWall"))
//           {
//               AxonWallBoundary2D(n, cnt, Fwd, inarray);
//           }

//           cnt += m_fields[0]->GetBndCondExpansions()[n]->GetExpSize();
//       }
//   }

//     void MMFNeuralEP::AxonWallBoundary2D(
//         int bcRegion,
//         int cnt,
//         Array<OneD, Array<OneD, NekDouble> > &Fwd,
//         Array<OneD, Array<OneD, NekDouble> > &physarray)
//     {
//         int nvariables = physarray.size();

//         // Adjust the physical values of the trace to take
//         // user defined boundaries into account
//         int i, e, id1, id2, npts;

//         for(e = 0; e <
//         m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExpSize(); ++e)
//         {
//             npts =
//             m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExp(e)->GetNumPoints(0);
//             id1  =
//             m_fields[0]->GetBndCondExpansions()[bcRegion]->GetPhys_Offset(e)
//             ; id2  =
//             m_fields[0]->GetTrace()->GetPhys_Offset(m_fields[0]->GetTraceMap()->GetBndCondCoeffsToGlobalCoeffsMap(cnt+e));

//             switch(m_expdim)
//             {
//             case 1:
//                 {
//                 // negate the forward flux
//                 Vmath::Neg(npts,&Fwd[1][id2],1);
//                 }
//                 break;
//             case 2:
//                 {
//                 Array<OneD, NekDouble> tmp_n(npts);
//                 Array<OneD, NekDouble> tmp_t(npts);

//                 Vmath::Vmul(npts,&Fwd[1][id2],1,&m_traceNormals[0][id2],1,&tmp_n[0],1);
//                 Vmath::Vvtvp(npts,&Fwd[2][id2],1,&m_traceNormals[1][id2],1,&tmp_n[0],1,&tmp_n[0],1);

//                 Vmath::Vmul(npts,&Fwd[1][id2],1,&m_traceNormals[1][id2],1,&tmp_t[0],1);
//                 Vmath::Vvtvm(npts,&Fwd[2][id2],1,&m_traceNormals[0][id2],1,&tmp_t[0],1,&tmp_t[0],1);

//                 // negate the normal flux
//                 Vmath::Neg(npts,tmp_n,1);

//                 // rotate back to Cartesian
//                 Vmath::Vmul(npts,&tmp_t[0],1,&m_traceNormals[1][id2],1,&Fwd[1][id2],1);
//                 Vmath::Vvtvm(npts,&tmp_n[0],1,&m_traceNormals[0][id2],1,&Fwd[1][id2],1,&Fwd[1][id2],1);

//                 Vmath::Vmul(npts,&tmp_t[0],1,&m_traceNormals[0][id2],1,&Fwd[2][id2],1);
//                 Vmath::Vvtvp(npts,&tmp_n[0],1,&m_traceNormals[1][id2],1,&Fwd[2][id2],1,&Fwd[2][id2],1);
//                 }
//                 break;
//             case 3:
//                 ASSERTL0(false,"3D not implemented for Shallow Water
//                 Equations"); break;
//             default:
//                 ASSERTL0(false,"Illegal expansion dimension");
//             }

//             // copy boundary adjusted values into the boundary expansion
//             for (i = 0; i < nvariables; ++i)
//             {
//                 Vmath::Vcopy(npts,&Fwd[i][id2],
//                 1,&(m_fields[i]->GetBndCondExpansions()[bcRegion]->UpdatePhys())[id1],1);
//             }
//         }
//     }

// Implicit solve for NeuralEP solver
// void MMFNeuralEP::DoImplicitSolveNeuralEP2p1Dfiber(
//     const Array<OneD, const Array<OneD, NekDouble>> &inarray,
//     const Array<OneD, const Array<OneD, NekDouble>> &MF1st,
//     Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time,
//     const NekDouble lambda)
// {
//     boost::ignore_unused(MF1st, time);

//     for (int nfib=0; nfib<m_nfibers; ++nfib)
//     {
//         int fnq   = m_fiberfields[nfib]->GetNpoints();

//         const NekDouble R_f = m_fiberneurons[nfib]->GetRecistanceValue();
//         const NekDouble C_n = m_fiberneurons[nfib]->GetCapacitanceValue(1);

//         // std::cout << "inarray = " << Vmath::Vmax(fnq, inarray[nfib], 1) <<
//         // ", outarray = " << Vmath::Vmax(fnq, outarray[nfib], 1)  <<
//         std::endl;

//         StdRegions::ConstFactorMap factors;
//         factors[StdRegions::eFactorTau] = m_Helmtau;

//         // // factors[StdRegions::eFactorLambda] = 1.0 / lambda;
//         // // Cm dVm/dt = (1/rf) \nabla^2 Vm
//         // // Vm^{n+1} = Vm^{n} + \Delta t/(rf*Cm) \nabla^2 Vn
//         // // m_beta = Relative extracellular resistance = r_ex / r_f =
//         \sigma_f / \sigma_ex

//         // NekDouble Cv = C_n * R_f * ( (m_beta_e + 1.0)/m_beta_e );
//         // NekDouble Cv = C_n * R_f * (m_ratio_re_ri + 1.0);
//         NekDouble Cv = C_n * R_f * (m_ratio_re_ri + 0.0);
//         factors[StdRegions::eFactorLambda] = Cv / lambda ;

//         // We solve ( \nabla^2 - HHlambda ) Y[i] = rhs [i]
//         // inarray = input: \hat{rhs} -> output: \hat{Y}
//         // outarray = output: \hat{Y} where \hat = modal coeffs

//         // For the variable of membrane potential: Multiply 1.0/timestep
//         Vmath::Smul(fnq, -factors[StdRegions::eFactorLambda], inarray[nfib],
//         1,
//                     m_fiberfields[nfib]->UpdatePhys(), 1);

//         // SetBoundaryConditions(time);

//         // Solve a system of equations with Helmholtz solver and transform
//         // back into physical space.
//         m_fiberfields[nfib]->HelmSolve(m_fiberfields[nfib]->GetPhys(),
//         m_fiberfields[nfib]->UpdateCoeffs(),
//                             NullFlagList, factors, m_fibervarcoeff[nfib]);

//         m_fiberfields[nfib]->BwdTrans(m_fiberfields[0]->GetCoeffs(),
//         outarray[nfib]); m_fiberfields[nfib]->SetPhysState(true);
//     }
// }

// // Implicit solve for NeuralEP solver
// void MMFNeuralEP::DoImplicitSolveNeuralEP2p1D(
//     const Array<OneD, const Array<OneD, NekDouble>> &inarray,
//     const Array<OneD, const Array<OneD, NekDouble>> &MF1st,
//     Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time,
//     const NekDouble lambda)
// {
//     boost::ignore_unused(MF1st, time, lambda);

//     int nq   = m_fields[0]->GetNpoints();
//     int nvar = inarray.size();

//     for (int i = 0; i < nvar; ++i)
//     {
//         Vmath::Vcopy(nq, &inarray[i][0], 1, &outarray[i][0], 1);
//     }
// }

void MMFNeuralEP::UpdateFibertoField(
    const int nfib, const Array<OneD, const NekDouble> &fiberinarray,
    Array<OneD, NekDouble> &outarray)
{
    int fnq = m_fiberfields[nfib]->GetNpoints();
    // int nq         = m_fields[0]->GetNpoints();

    int index;
    int maxindx = m_fiberindex[0][0].size();

    // Use the average of the vertex values for the corresponding fiber
    for (int i = 0; i < fnq; ++i)
    {
        for (int k = 0; k < maxindx; ++k)
        {
            index = m_fiberindex[nfib][i][k];
            if (index > 0)
            {
                outarray[index] = fiberinarray[i];
            }
        }
    }
}

Array<OneD, NekDouble> MMFNeuralEP::ExtractFiberValue(
    const int nfib, const Array<OneD, const NekDouble> &inarray)
{
    int fnq = m_fiberfields[nfib]->GetNpoints();

    Array<OneD, NekDouble> fiberoutarray(fnq, 0.0);

    int index, n = 0;
    int maxindx = m_fiberindex[0][0].size();

    NekDouble sum = 0.0;
    // Use the average of the vertex values for the corresponding fiber
    for (int i = 0; i < fnq; ++i)
    {
        sum = 0.0;
        n   = 0;
        for (int k = 0; k < maxindx; ++k)
        {
            index = m_fiberindex[nfib][i][k];
            if (index > 0)
            {
                sum = sum + inarray[index];
                n++;
            }
        }

        if (n > 0)
        {
            fiberoutarray[i] = 1.0 * sum / n;
        }
    }

    return fiberoutarray;
}

// We Return Y[i] = rhs [i] without no Helomsolver
void MMFNeuralEP::DoNullSolve(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time,
    const NekDouble lambda)
{
    boost::ignore_unused(lambda, time);

    int nvariables = inarray.size();
    int nq         = m_fields[0]->GetNpoints();

    for (int i = 0; i < nvariables; ++i)
    {
        Vmath::Vcopy(nq, &inarray[i][0], 1, &outarray[i][0], 1);
    }
}

void MMFNeuralEP::DoOdeRhsNeuralEPPT(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time)
{
    int nvar = m_fields.size();
    int nq   = m_fields[0]->GetNpoints();

    // Compute the reaction function divided by Cm or Cn.
    m_neuron->TimeIntegrate(m_NodeZone[0], inarray[0], outarray[0], time, m_diameter, m_Temperature);

    for (int i=1; i < nvar; ++i)
    {
        Vmath::Vcopy(nq, m_neuron->GetNeuronSolution(i), 1, m_fields[i]->UpdatePhys(), 1);
    }

    Array<OneD, Array<OneD, NekDouble>> RHSstimulus(nvar);
    for (int i = 0; i < nvar; ++i)
    {
        RHSstimulus[i] = Array<OneD, NekDouble>(nq, 0.0);
    }

    for (unsigned int i = 0; i < m_stimulus.size(); ++i)
    {
        m_stimulus[i]->Update(RHSstimulus, time);
    }

    // ONLY simulation at node zone: No excitation at myelinated region
    StimulusAtNode(RHSstimulus[0]);

    // Add it to the RHS
    Vmath::Vadd(nq, RHSstimulus[0], 1, outarray[0], 1, outarray[0], 1);
}


void MMFNeuralEP::DoOdeRhsNeuralEP1D(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time)
{
    int nvar = m_fields.size();
    int nq   = m_fields[0]->GetNpoints();
    NekDouble Rf = m_neuron->GetRecistanceValue();
    // NekDouble Cn = m_neuron->GetCapacitanceValue(1);

    // Compute the reaction function divided by Cm or Cn.
    m_neuron->TimeIntegrate(m_NodeZone[0], inarray[0], outarray[0], time, m_diameter, m_Temperature);

    Array<OneD, Array<OneD, NekDouble>> RHSstimulus(nvar);
    for (int i = 0; i < nvar; ++i)
    {
        RHSstimulus[i] = Array<OneD, NekDouble>(nq, 0.0);
    }

    for (unsigned int i = 0; i < m_stimulus.size(); ++i)
    {
        m_stimulus[i]->Update(RHSstimulus, time);
    }

    // ONLY simulation at node zone: No excitation at myelinated region
    StimulusAtNode(RHSstimulus[0]);

    // Add it to the RHS
    Vmath::Vadd(nq, RHSstimulus[0], 1, outarray[0], 1, outarray[0], 1);

        switch (m_projectionType)
        {
            case MultiRegions::eDiscontinuous:
            {
                std::string diffName;

                // Do not forwards transform initial condition
                m_homoInitialFwd = false;

                m_session->LoadSolverInfo("DiffusionType", diffName, "LDG");
                m_diffusion = SolverUtils::GetDiffusionFactory().CreateInstance(
                    diffName, diffName);
                m_diffusion->SetFluxVector(&MMFNeuralEP::GetFluxVector, this);
                m_diffusion->InitObject(m_session, m_fields);
                break;
            }

            case MultiRegions::eGalerkin:
            case MultiRegions::eMixed_CG_Discontinuous:
            {
                if (m_explicitDiffusion)
                {
                    ASSERTL0(false, "Explicit Galerkin diffusion not set up.");
                }
            }
        }

    // Multiply by 1/Cm for myeline or 1/Cm for node
    if (m_explicitDiffusion)
    {
        Array<OneD, Array<OneD, NekDouble>> outarrayDiff(nvar);
        for (int i = 0; i < nvar; ++i)
        {
            outarrayDiff[i] = Array<OneD, NekDouble>(nq, 0.0);
        }

        m_diffusion->Diffuse(nvar, m_fields, inarray, outarrayDiff);                         

        for (int i = 0; i < nvar; ++i)
        {
            Vmath::Smul(nq, 1.0/Rf, &outarrayDiff[i][0], 1, &outarrayDiff[i][0], 1);
            Vmath::Vmul(nq, &m_NeuralCm[0][0], 1, &outarrayDiff[i][0], 1, &outarrayDiff[i][0], 1);

            Vmath::Vadd(nq, &outarrayDiff[i][0], 1, &outarray[i][0], 1, &outarray[i][0], 1);
        }
    }

    // if (m_explicitDiffusion)
    // {
    //     int nq = m_fields[0]->GetNpoints();

    //     // Laplacian only to the first variable
    //     Array<OneD, NekDouble> Laplacian(nq);
    //     WeakDGMMFDiffusion(0, inarray[0], Laplacian, time);
    //     // WeakDGMMFNeuralEP(0, inarray[0], Laplacian, time);

    //     Vmath::Smul(nq, 1.0 / (Cn * Rf), Laplacian, 1, Laplacian, 1);

    //     Vmath::Vadd(nq, &Laplacian[0], 1, &outarray[0][0], 1, &outarray[0][0],
    //                 1);
    // }
}

void MMFNeuralEP::GetFluxVector(
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    const Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &qfield,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &viscousTensor)
{
    boost::ignore_unused(inarray);

    unsigned int nDim              = qfield.size();
    unsigned int nConvectiveFields = qfield[0].size();
    unsigned int nPts              = qfield[0][0].size();

    for (unsigned int j = 0; j < nDim; ++j)
    {
        for (unsigned int i = 0; i < nConvectiveFields; ++i)
        {
            // Vmath::Smul(nPts, m_epsilon[j], qfield[j][i], 1, viscousTensor[j][i],
            //             1);
            Vmath::Vcopy(nPts, qfield[j][i], 1, viscousTensor[j][i], 1);
        }
    }
}


void MMFNeuralEP::DoOdeRhsNeuralEP2Dmono(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time)
{
    int nvar = m_fields.size();
    int nq   = m_fields[0]->GetNpoints();

    NekDouble Rf = m_neuron->GetRecistanceValue();
    NekDouble Cn = m_neuron->GetCapacitanceValue(1);

    for (int i=0; i<nvar; ++i)
    {
        outarray[i] = Array<OneD, NekDouble>(nq, 0.0);
    }

    // Compute the reaction function divided by Cm or Cn.
    m_neuron->TimeIntegrate(m_NodeZone[0], inarray[0], outarray[0], time, m_diameter, m_Temperature);

    Array<OneD, Array<OneD, NekDouble>> RHSstimulus(nvar);
    for (int i = 0; i < nvar; ++i)
    {
        RHSstimulus[i] = Array<OneD, NekDouble>(nq, 0.0);
    }

    for (unsigned int j = 0; j < m_stimulus.size(); ++j)
    {
        m_stimulus[j]->Update(RHSstimulus, time);
    }

    // ONLY simulation at node zone: No excitation at myelinated region
    StimulusAtNode(RHSstimulus[0]);

    // Add it to the RHS
    Vmath::Vadd(nq, RHSstimulus[0], 1, outarray[0], 1, outarray[0], 1);

    if (m_explicitDiffusion)
    {
        // Laplacian only to the first variable
        Array<OneD, NekDouble> Laplacian(nq);
        WeakDGMMFDiffusion(0, inarray[0], Laplacian, time);

        Vmath::Smul(nq, 1.0 / (Cn * Rf), Laplacian, 1, Laplacian, 1);
        Vmath::Vadd(nq, &Laplacian[0], 1, &outarray[0][0], 1, &outarray[0][0], 1);
    }
}


void MMFNeuralEP::DoOdeRhsNeuralEP2Dbi(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time)
{
    int nvar = m_fields.size();
    int nq   = m_fields[0]->GetNpoints();

    NekDouble Rf = m_neuron->GetRecistanceValue();
    NekDouble Cn = m_neuron->GetCapacitanceValue(1);

    for (int i=0; i<nvar; ++i)
    {
        outarray[i] = Array<OneD, NekDouble>(nq, 0.0);
    }

    // Compute the reaction function divided by Cm or Cn.
    m_neuron->TimeIntegrate(m_NodeZone[0], inarray[0], outarray[0], time, m_diameter, m_Temperature);

    Array<OneD, Array<OneD, NekDouble>> RHSstimulus(nvar);
    for (int i = 0; i < nvar; ++i)
    {
        RHSstimulus[i] = Array<OneD, NekDouble>(nq, 0.0);
    }

    for (unsigned int j = 0; j < m_stimulus.size(); ++j)
    {
        m_stimulus[j]->Update(RHSstimulus, time);
    }

    // ONLY simulation at node zone: No excitation at myelinated region
    StimulusAtNode(RHSstimulus[0]);

    // Add it to the RHS
    Vmath::Vadd(nq, RHSstimulus[0], 1, outarray[0], 1, outarray[0], 1);

    // Compute phi_e to satisfy the following equation
    // \nabla \cdot ( (\signa_e + \sigma_i) \nabla \phi_e) = - \nabla \cdot
    // (\sigma_i \nabla \phi_m)
    Array<OneD, NekDouble> phie(nq,0.0);
    SolveHelmholtzDiffusion(m_NodeZone[0], inarray[0], m_unitmovingframes, m_phievarcoeff, phie);

    // Add the current changes by the external current
    Array<OneD, NekDouble> extcurrent(nq,0.0);
    extcurrent = ComputeMMFDiffusion(m_unitmovingframes, m_fields[1]->GetPhys());
    Vmath::Smul(nq, 1.0 / (Cn * Rf), extcurrent, 1, extcurrent, 1);

    // Let the extcurrent be zero at Myeline nodes (-1).
    OnlyValideinNode(m_NodeZone[0], extcurrent);

    // add divergence of phie to the current
    Vmath::Vadd(nq, &extcurrent[0], 1, &outarray[0][0], 1, &outarray[0][0], 1);

    if (m_explicitDiffusion)
    {
        // Laplacian only to the first variable
        Array<OneD, NekDouble> Laplacian(nq);
        WeakDGMMFDiffusion(0, inarray[0], Laplacian, time);

        Vmath::Smul(nq, 1.0 / (Cn * Rf), Laplacian, 1, Laplacian, 1);
        Vmath::Vadd(nq, &Laplacian[0], 1, &outarray[0][0], 1, &outarray[0][0], 1);
    }
}

void MMFNeuralEP::DoOdeRhsNeuralEP2DEmbed(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time)
{
    int nvar = m_fields.size();
    int nq   = m_fields[0]->GetNpoints();

    // Compute the reaction function divided by Cm or Cn.
    m_neuron->TimeIntegrate(m_NodeZone[0], inarray[0], outarray[0], time,
                            m_Temperature);

    Array<OneD, Array<OneD, NekDouble>> RHSstimulus(nvar);
    for (int i = 0; i < nvar; ++i)
    {
        RHSstimulus[i] = Array<OneD, NekDouble>(nq, 0.0);
    }

    for (unsigned int j = 0; j < m_stimulus.size(); ++j)
    {
        m_stimulus[j]->Update(RHSstimulus, time);
    }

    // ONLY simulation at node zone:
    // No excitation at myelinated region
    // const NekDouble Cm = m_neuron->GetCapacitanceValue(0);

    const NekDouble Rf = m_neuron->GetRecistanceValue();
    const NekDouble Cn = m_neuron->GetCapacitanceValue(1);
    for (int k = 0; k < nq; ++k)
    {
        // Exict only the first and the second node.
        if ((m_NodeZone[0][k] >= 0) && (m_NodeZone[0][k] <= 1))
        {
            outarray[0][k] = outarray[0][k] + RHSstimulus[0][k] / Cn;
        }
    }

    // Compute phi_e to satisfy the following equation
    // \nabla \cdot ( (\signa_e + \sigma_i) \nabla \phi_e) = - \nabla \cdot
    // (\sigma_i \nabla \phi_m)
    Array<OneD, NekDouble> phie(nq);
    SolveHelmholtzDiffusion(m_NodeZone[0], inarray[0], m_unitmovingframes, m_phievarcoeff, phie);

    // Add the current changes by the external current
    Array<OneD, NekDouble> extcurrent;
    extcurrent =
        ComputeCovariantDiffusion(m_unitmovingframes, phie);
    Vmath::Smul(nq, 1.0 / (Cn * Rf), extcurrent, 1, extcurrent, 1);

    // Let the extcurrent be zero at Myeline nodes (-1).
    OnlyValideinNode(m_NodeZone[0], extcurrent);

    // add divergence of phie to the current
    Vmath::Vadd(nq, &extcurrent[0], 1, &outarray[0][0], 1, &outarray[0][0], 1);

    // Multiply by 1/Cm for myeline or 1/Cm for node
    if (m_explicitDiffusion)
    {
        // Laplacian only to the first variable
        Array<OneD, NekDouble> Laplacian(nq);
        // WeakDGMMFNeuralEP(0, inarray[0], Laplacian, time);

        Vmath::Vmul(nq, &m_NeuralCmRf[0], 1, &Laplacian[0], 1, &Laplacian[0],
                    1);
        Vmath::Vadd(nq, &Laplacian[0], 1, &outarray[0][0], 1, &outarray[0][0],
                    1);
    }
}

void MMFNeuralEP::StimulusAtNode(Array<OneD, NekDouble> &outarray)
{
    int nq   = m_fields[0]->GetNpoints();

    const NekDouble Cn = m_neuron->GetCapacitanceValue(1);
    for (int k = 0; k < nq; ++k)
    {
        if ( (m_NodeZone[0][k] == 0) || (m_NodeZone[0][k] == 1) )
        {
            outarray[k] = outarray[k] / Cn;
        }

        else
        {
            outarray[k] = 0.0;
        }
    }
}

void MMFNeuralEP::OnlyValideinNode(const Array<OneD, const int> &NodeZone,
                                   Array<OneD, NekDouble> &outarray)
{
    int nq = m_fields[0]->GetNpoints();

    for (int i = 0; i < nq; ++i)
    {
        if (NodeZone[i] < 0)
        {
            outarray[i] = 0.0;
        }
    }
}

// Compute phi_e from the given distribution of phi_m
// \nabla \cdot ( (1 + \rho) \mathbf{e}_1 + \mathbf{e}_2 ) ( \nabla \phi_e ))
//                         = - \nabla \cdot \mathbf{e}_1 \nabla \phi_m
void MMFNeuralEP::SolveHelmholtzDiffusion(
    const Array<OneD, const int> &NodeZone,
    const Array<OneD, const NekDouble> &phim,
    const Array<OneD, const Array<OneD, NekDouble>> &DiffusionMF,
    StdRegions::VarCoeffMap &Helmvarcoeff,
    Array<OneD, NekDouble> &outarray)
{
    boost::ignore_unused(DiffusionMF);

    int nq = m_fields[0]->GetNpoints();

    // Solve the Poisson equation: \nabla (\sigma_e + \sigma_i ) phi_e = \nabla
    // \sigma_i \nabla phi_m
    StdRegions::ConstFactorMap phiefactors;
    phiefactors[StdRegions::eFactorTau]    = m_Helmtau;
    phiefactors[StdRegions::eFactorLambda] = 0.0;

    // // Compute \nabla \sigma_i \nabla phi_m and use it as point sources for
    // phi_e.
    // // This is equivalently achieved by removing all the point sources in
    // myelinnated fiber region.
    Array<OneD, NekDouble> phimLaplacian(nq,0.0);
    phimLaplacian = ComputeEuclideanDiffusion(phim);
    // phimLaplacian = ComputeMMFDiffusion(DiffusionMF, phim);

    // Only nonzero for node.
    OnlyValideinNode(NodeZone, phimLaplacian);
    NekDouble phimavg = -1.0 * AvgInt(phimLaplacian);
    Vmath::Sadd(nq, phimavg, phimLaplacian, 1, phimLaplacian, 1);
    
    Vmath::Smul(nq, -1.0, phimLaplacian, 1, m_fields[1]->UpdatePhys(), 1);
    // Compute phie distribution
    // SetMembraneBoundaryCondition();
    SetBoundaryConditions(0.0);

    m_fields[1]->HelmSolve(m_fields[1]->GetPhys(), m_fields[1]->UpdateCoeffs(), phiefactors, Helmvarcoeff);
    m_fields[1]->BwdTrans(m_fields[1]->GetCoeffs(), m_fields[1]->UpdatePhys());
    m_fields[1]->SetPhysState(true);

    NekDouble phieavg = -1.0 * Average(m_fields[1]->GetPhys());
    Vmath::Sadd(nq, phieavg, m_fields[1]->GetPhys(), 1, m_fields[1]->UpdatePhys(), 1);

    outarray = m_fields[1]->GetPhys();
}

// Array<OneD, NekDouble> MMFNeuralEP::Computephie(
//            const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
//            const Array<OneD, const NekDouble> &phim)
// {

//     int nq = m_fields[0]->GetNpoints();
//     int ncoeffs         = GetNcoeffs();

//     // Solve the Poisson equation: \nabla (\sigma_e + \sigma_i ) phi_e =
//     \nabla \sigma_i \nabla phi_m StdRegions::ConstFactorMap phiefactors;
//     phiefactors[StdRegions::eFactorTau] = m_Helmtau;
//     phiefactors[StdRegions::eFactorLambda] = 0.0 ;

//     // // Compute \nabla \sigma_i \nabla phi_m and use it as point sources
//     for phi_e.
//     // // This is equivalently achieved by removing all the point sources in
//     myelinnated fiber region. Array<OneD, NekDouble> phimdist(nq,0.0);
//     phimdist = ComputeCovariantDiffusion(movingframes, phim, 1);

//     Vmath::Sadd(nq, -1.0 * AvgInt(phimdist), phimdist, 1, phimdist, 1);
//     Vmath::Neg(nq, phimdist, 1);

//     // Compute phie distribution
//     Array<OneD, NekDouble> outarray(nq,0.0);
//     Array<OneD, NekDouble> tmpcoeff(ncoeffs);

//     m_fields[0]->HelmSolve(phimdist, tmpcoeff, NullFlagList, phiefactors,
//     m_varcoeff);
//     // m_contField->HelmSolve(phimdist, tmpcoeff, NullFlagList, phiefactors,
//     m_varcoeff);

//     m_fields[0]->BwdTrans(tmpcoeff, outarray);

//     return outarray;
// }

// void MMFNeuralEP::DoOdeRhsNeuralEP2p1Dfiber(
//     const Array<OneD, const Array<OneD, NekDouble>> &inarray,
//     const Array<OneD, const Array<OneD, NekDouble>> &MF1st,
//     Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time)
// {
//     boost::ignore_unused(MF1st);
//    // int nvar = m_fields.size();

//     for (int nfib=0; nfib<m_nfibers; nfib++)
//     {
//         int fnq = m_fiberfields[nfib]->GetNpoints();
//         // int fncoeffs = m_fiberfields[nfib]->GetNcoeffs();

//         // Compute the reaction function divided by Cm or Cn.
//         m_fiberneurons[nfib]->TimeIntegrate(m_NodeZone[nfib], inarray[nfib],
//         outarray[nfib], time, m_Temperature);

//         Array<OneD, Array<OneD, NekDouble>> RHSstimulus(1);
//         for (int i=0; i<1; ++i)
//         {
//             RHSstimulus[i] = Array<OneD, NekDouble>(fnq, 0.0);
//         }

//         m_fiberstimulus[nfib]->Update(RHSstimulus, time);

//         // ONLY simulation at node zone:
//         // No excitation at myelinated region
//         const NekDouble Cn = m_fiberneurons[nfib]->GetCapacitanceValue(1);
//         const NekDouble Rf = m_fiberneurons[nfib]->GetRecistanceValue();

//         // Array<OneD, NekDouble> MMFSystem::WeakDGMMFLDGFiber(
//         //     MultiRegions::ExpListSharedPtr &fiberfield,
//         //     const Array<OneD, const NekDouble> &inarray,
//         //     const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
//         //     const Array<OneD, const Array<OneD, NekDouble>> &ncdotMFFwd,
//         //     const Array<OneD, const Array<OneD, NekDouble>> &ncdotMFBwd,
//         //     const Array<OneD, const Array<OneD, NekDouble>> &DivMF,
//         //     const NekDouble time)

//         Array<OneD, NekDouble> phie(fnq,0.0);
//         Array<OneD, NekDouble> d2phie(fnq,0.0);

//         phie = ExtractFiberValue(nfib, m_fields[0]->GetPhys());

//         m_fiberfields[nfib]->PhysDeriv(1, phie, d2phie);
//         m_fiberfields[nfib]->PhysDeriv(1, d2phie, d2phie);
//         Vmath::Smul(fnq, m_AnisotropyStrength, d2phie, 1, d2phie, 1);

//         std::cout << "phie = " << Vmath::Vmax(fnq, phie, 1) << ", d2phie = "
//         << Vmath::Vmax(fnq, d2phie, 1) << std::endl;

//         // Array<OneD, NekDouble> Laplacian(fnq,0.0);
//         // Array<OneD, NekDouble> tmpc(fncoeffs);

//         // tmpc = WeakDGMMFLDGFiber(
//         //     m_fiberfields[nfib], phie, m_fibermovingframes[nfib],
//         //     m_ncdotfiberMFFwd[nfib], m_ncdotfiberMFBwd[nfib],
//         m_DivfiberMF[nfib], time);

//         // m_fiberfields[nfib]->MultiplyByElmtInvMass(tmpc, tmpc);
//         // m_fiberfields[nfib]->BwdTrans(tmpc, Laplacian);

//         // std::cout << "phie = " << Vmath::Vmax(fnq, phie, 1)
//         // << ",  tmpc = " << Vmath::Vmax(fncoeffs, tmpc, 1)
//         // << ",  Laplacian = " << Vmath::Vmax(fnq, Laplacian, 1) <<
//         std::endl;

//         for (int i=0;i<fnq; ++i)
//         {
//             if(m_NodeZone[nfib][i]==1)
//             {
//                 outarray[nfib][i] = outarray[nfib][i] + (
//                 RHSstimulus[nfib][i] + d2phie[i]/Rf  ) / Cn;
//             }
//         }

//         // std::cout << "max RHSstimulus = " << Vmath::Vmax(fnq,
//         RHSstimulus[0], 1) <<
//         //         ", inarray = " << Vmath::Vmax(fnq, inarray[nfib], 1)  <<
//         // ", outarray = " << Vmath::Vmax(fnq, outarray[nfib], 1)  <<
//         std::endl;
//     }

//     // Multiply by 1/Cm for myeline or 1/Cm for node
//     if (m_explicitDiffusion)
//     {
//         int nq = m_fields[0]->GetNpoints();

//         // Laplacian only to the first variable
//         Array<OneD, NekDouble> Laplacian(nq);
//         // WeakDGMMFNeuralEP(0, inarray[0], Laplacian, time);

//         Vmath::Vmul(nq, &m_NeuralCmRf[0], 1, &Laplacian[0], 1, &Laplacian[0],
//         1); Vmath::Vadd(nq, &Laplacian[0], 1, &outarray[0][0], 1,
//         &outarray[0][0], 1);
//     }
// }

// void MMFNeuralEP::DoOdeRhsNeuralEP2p1D(
//     const Array<OneD, const Array<OneD, NekDouble>> &inarray,
//     const Array<OneD, const Array<OneD, NekDouble>> &MF1st,
//     Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time)
// {
//     boost::ignore_unused(inarray, MF1st, time);

//     int nvar = m_fields.size();
//     int nq = m_fields[0]->GetNpoints();

//     for (int i=0; i<nvar; ++i)
//     {
//         outarray[i] = Array<OneD, NekDouble>(nq,0.0);
//     }
// }

void MMFNeuralEP::v_SetInitialConditions(NekDouble initialtime,
                                         bool dumpInitialConditions,
                                         const int domain)
{
    boost::ignore_unused(domain, dumpInitialConditions);

    int nq = GetTotPoints();

    switch (m_NeuralEPType)
    {
        case eNeuralEPPT:
        case eNeuralEP1D:
        case eNeuralEP2Dmono:
        case eNeuralEP2Dbi:
        {
            m_neuron->Initialise();

            // Read initial condition from xml file
            EquationSystem::v_SetInitialConditions(initialtime, false);

            Array<OneD, Array<OneD, NekDouble>> tmp(1);
            tmp[0] = Array<OneD, NekDouble>(nq);

            Array<OneD, NekDouble> initialcondition(nq);
            Vmath::Vcopy(nq, m_fields[0]->GetPhys(), 1, tmp[0], 1);
            Vmath::Vcopy(nq, tmp[0], 1, initialcondition, 1);
            for (unsigned int i = 0; i < m_stimulus.size(); ++i)
            {
                m_stimulus[i]->Update(tmp, initialtime);
                StimulusAtNode(tmp[0]);
                m_fields[0]->SetPhys(tmp[0]);
            }
            
            // Check NodeZone and moving frames
            CheckNodeZoneMF(m_movingframes, m_NodeZone, tmp[0]);

            break;
        }

            // case eNeuralEP2p1D:
            // {
            //     for (int nfib=0; nfib<m_nfibers; nfib++)
            //     {
            //         m_fiberneurons[nfib]->Initialise();
            //     }

            //     // Read initial condition from xml file
            //     EquationSystem::v_SetInitialConditions(initialtime, false);

            //     Array<OneD, Array<OneD, NekDouble>> tmp(1);
            //     tmp[0] = Array<OneD, NekDouble>(nq);

            //     Array<OneD, NekDouble> initialcondition(nq);
            //     Vmath::Vcopy(nq, m_fields[0]->GetPhys(), 1, tmp[0], 1);
            //     Vmath::Vcopy(nq, tmp[0], 1, initialcondition, 1);
            //     for (unsigned int i = 0; i < m_stimulus.size(); ++i)
            //     {
            //         m_stimulus[i]->Update(tmp, initialtime);
            //         m_fields[0]->SetPhys(tmp[0]);
            //     }

            //     // Only the excited regions are considered for m_InitExcitation
            //     Vmath::Vsub(nq, tmp[0], 1, initialcondition, 1, initialcondition, 1); 
            //     m_ValidTimeMap = ComputeTimeMapInitialZone(initialcondition);
            //
            //     Array<OneD, NekDouble> x0(nq);
            //     Array<OneD, NekDouble> x1(nq);
            //     Array<OneD, NekDouble> x2(nq);

            //     m_fields[0]->GetCoords(x0, x1, x2);
            // }
            // break;

        default:
        {
            EquationSystem::v_SetInitialConditions(initialtime, false);
            break;
        }
    }
    
    std::cout << "Initial: max um = "
              << Vmath::Vmax(nq, m_fields[0]->GetPhys(), 1) << std::endl;

    if (m_fields.size() > 1)
    {
        std::cout << "Initial: max ue = "
                  << Vmath::Vmax(nq, m_fields[1]->GetPhys(), 1) << std::endl;
    }

    if (dumpInitialConditions)
    {
        std::string outname;
        outname = m_sessionName + "_initial.chk";

        WriteFld(outname);
    }
}

void MMFNeuralEP::SetMembraneBoundaryCondition(const NekDouble time)
{
    std::string varName;
    int cnt        = 0;
    int nvariables = m_fields.size();
    int nTracePts  = GetTraceTotPoints();
    int nq         = GetTotPoints();

    // Extract trace for boundaries. Needs to be done on all processors to avoid
    // deadlock.
    Array<OneD, Array<OneD, NekDouble>> inarray(nvariables);
    Array<OneD, Array<OneD, NekDouble>> Fwd(nvariables);
    for (int i = 0; i < nvariables; ++i)
    {
        inarray[i] = Array<OneD, NekDouble>(nq);
        Fwd[i]     = Array<OneD, NekDouble>(nTracePts);

        Vmath::Vcopy(nq, &m_fields[i]->GetPhys()[0], 1, &inarray[i][0], 1);
        m_fields[i]->ExtractTracePhys(inarray[i], Fwd[i]);
    }

    // loop over Boundary Regions
    for (int n = 0; n < m_fields[0]->GetBndConditions().size(); ++n)
    {
        // Wall Boundary Condition
        if (boost::iequals(m_fields[0]->GetBndConditions()[n]->GetUserDefined(),
                           "Membrane"))
        {
            MembraneBoundary2D(n, cnt, Fwd, inarray);
        }

        else
        {
            for (int i = 0; i < nvariables; ++i)
            {
                varName = m_session->GetVariable(i);
                m_fields[i]->EvaluateBoundaryConditions(time, varName);
            }
        }

        cnt += m_fields[0]->GetBndCondExpansions()[n]->GetExpSize();
    }
}

// TO DO: IMPLEMENT Nonhomogeneous Neurann boundary conditions
//----------------------------------------------------
/**
 * @brief Wall boundary condition.
 */
void MMFNeuralEP::MembraneBoundary2D(
    int bcRegion, int cnt, Array<OneD, Array<OneD, NekDouble>> &Fwd,
    Array<OneD, Array<OneD, NekDouble>> &physarray)
{
    int nvariables = physarray.size();

    // Adjust the physical values of the trace to take
    // user defined boundaries into account
    int id1, id2, npts;

    for (int e = 0;
         e < m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExpSize(); ++e)
    {
        // npts = m_fields[0]->GetBndCondExpansions()[bcRegion]->
        //     GetExp(e)->GetTotPoints();
        // id1  = m_fields[0]->GetBndCondExpansions()[bcRegion]->
        //     GetPhys_Offset(e);
        // id2  = m_fields[0]->GetTrace()->GetPhys_Offset(
        //             m_fields[0]->GetTraceMap()->
        //                         GetBndCondCoeffsToGlobalCoeffsMap(cnt+e));

        npts = m_fields[0]
                   ->GetBndCondExpansions()[bcRegion]
                   ->GetExp(e)
                   ->GetNumPoints(0);
        id1 = m_fields[0]->GetBndCondExpansions()[bcRegion]->GetPhys_Offset(e);
        id2 = m_fields[0]->GetTrace()->GetPhys_Offset(
            m_fields[0]->GetTraceMap()->GetBndCondIDToGlobalTraceID(cnt + e));

        // copy boundary adjusted values into the boundary expansion
        for (int i = 0; i < nvariables; ++i)
        {
            Vmath::Vcopy(npts, &Fwd[i][id2], 1,
                         &(m_fields[i]
                               ->GetBndCondExpansions()[bcRegion]
                               ->UpdatePhys())[id1],
                         1);
        }
    }
}

// Array<OneD, NekDouble> MMFNeuralEP::PlanePhiWave()
// {
//     int nq = GetTotPoints();
//     Array<OneD, NekDouble> outarray(nq, 0.0);

//     Array<OneD, NekDouble> x(nq);
//     Array<OneD, NekDouble> y(nq);
//     Array<OneD, NekDouble> z(nq);

//     m_fields[0]->GetCoords(x, y, z);

//     NekDouble xmin, ymin, xmax;

//     xmin = Vmath::Vmin(nq, x, 1);
//     xmax = Vmath::Vmax(nq, x, 1);
//     ymin = Vmath::Vmin(nq, y, 1);

//     NekDouble xp, yp, xp2;
//     for (int i = 0; i < nq; i++)
//     {
//         switch (m_InitWaveType)
//         {
//             case eLeft:
//             {
//                 NekDouble radiusofinit;
//                 NekDouble frontstiff;

//                 radiusofinit = 3.0;
//                 frontstiff   = 0.5;

//                 xp = x[i] - xmin;
//                 outarray[i] =
//                     1.0 / (1.0 + exp((xp - radiusofinit) / frontstiff));
//             }
//             break;

//             case eBothEnds:
//             {
//                 NekDouble radiusofinit = 3.0;
//                 NekDouble frontstiff   = 0.1;

//                 xp  = x[i] - xmin;
//                 xp2 = x[i] - xmax;

//                 outarray[i] =
//                     1.0 / (1.0 +
//                            exp((sqrt(xp * xp) - radiusofinit) / frontstiff))
//                            +
//                     1.0 / (1.0 +
//                            exp((sqrt(xp2 * xp2) - radiusofinit) /
//                            frontstiff));
//             }
//             break;

//             case eCenter:
//             {
//                 NekDouble radiusofinit = 6.0;
//                 NekDouble frontstiff   = 0.1;

//                 // NekDouble xc = 0.5*(Vmath::Vmax(nq, x, 1) +
//                 // Vmath::Vmin(nq, x, 1));

//                 xp = x[i] - xmin;
//                 outarray[i] =
//                     1.0 / (1.0 + exp((xp - radiusofinit) / frontstiff));
//             }
//             break;

//             case eLeftBottomCorner:
//             {
//                 NekDouble radiusofinit = 6.0;
//                 NekDouble frontstiff   = 0.1;
//                 NekDouble bs           = 2.0;

//                 xp = x[i] - xmin;
//                 yp = y[i] - ymin;
//                 outarray[i] =
//                     1.0 /
//                     (1.0 + exp((sqrt(xp * xp + yp * yp) / bs - radiusofinit)
//                     /
//                                frontstiff));
//             }
//             break;

//             case ePoint:
//             {
//                 NekDouble xloc, yloc, zloc, rad;
//                 NekDouble radiusofinit = 5.0;

//                 xloc = x[i] - m_Initx;
//                 yloc = y[i] - m_Inity;
//                 zloc = z[i] - m_Initz;

//                 rad = sqrt(xloc * xloc + yloc * yloc + zloc * zloc);

//                 xloc = xloc / radiusofinit;
//                 yloc = yloc / radiusofinit;
//                 zloc = zloc / radiusofinit;

//                 if (rad < radiusofinit)
//                 {
//                     outarray[i] =
//                         exp(-(1.0 / 2.0) *
//                             (xloc * xloc + yloc * yloc + zloc * zloc));
//                 }

//                 else
//                 {
//                     outarray[i] = 0.0;
//                 }
//             }
//             break;

//             case eSpiralDock:
//             {
//                 NekDouble radiusofinit = 3.0;
//                 NekDouble frontstiff   = 0.1;
//                 xp                     = x[i] - 4.0;
//                 yp                     = y[i];
//                 outarray[i] =
//                     (1.0 / (1.0 + exp(2.0 * yp))) *
//                     (1.0 / (1.0 + exp(-2.0 * xp))) *
//                     (1.0 / (1.0 + exp((xp - radiusofinit) / frontstiff)));
//             }
//             break;

//             default:
//                 break;
//         } // namespace Nektar
//     }

//     return outarray;
// }

void MMFNeuralEP::v_EvaluateExactSolution(unsigned int field,
                                          Array<OneD, NekDouble> &outfield,
                                          const NekDouble time)
{
    EquationSystem::v_EvaluateExactSolution(field, outfield, time);
}

// Compute \int \nabla u \cdot e^{dir}
// void MMFNeuralEP::WeakDGDirectionalDeriv(
//     const int direction, const Array<OneD, const Array<OneD, NekDouble>>
//     &MF1st, const Array<OneD, const NekDouble> &InField, Array<OneD,
//     NekDouble> &OutField)
// {
//     int ncoeffs         = GetNcoeffs();
//     int nTracePointsTot = GetTraceNpoints();
//     int nq              = GetNpoints();

//     Array<OneD, NekDouble> physfield(nq);

//     // Get the variables in physical space
//     // already in physical space
//     Vmath::Vcopy(nq, InField, 1, physfield, 1);

//     Array<OneD, NekDouble> WeakDeriv(ncoeffs, 0.0);
//     Array<OneD, NekDouble> tmp(nq);

//     // Directional derivation with respect to the j'th moving frame
//     // tmp[j] = \nabla \physfield[i] \cdot \mathbf{e}^j
//     // Implemented at TriExp::v_IProductWRTDirectionalDerivBase_SumFa
//     m_fields[0]->IProductWRTDirectionalDerivBase(MF1st[direction],
//     physfield,
//                                                  WeakDeriv);

//     // if the NumericalFluxs function already includes the normal in the
//     output Array<OneD, NekDouble> Fwd(nTracePointsTot); Array<OneD,
//     NekDouble> Bwd(nTracePointsTot);

//     Array<OneD, NekDouble> flux(nTracePointsTot, 0.0);
//     Array<OneD, NekDouble> fluxFwd(nTracePointsTot);
//     Array<OneD, NekDouble> fluxBwd(nTracePointsTot);

//     // Evaluate numerical flux in physical space which may in
//     // general couple all component of vectors
//     m_fields[0]->GetFwdBwdTracePhys(physfield, Fwd, Bwd);

//     // evaulate upwinded m_fields[i]
//     Array<OneD, NekDouble> traceVn(nTracePointsTot, 0.0);
//     Array<OneD, NekDouble> tmptrace(nTracePointsTot);
//     for (int i = 0; i < m_spacedim; ++i)
//     {
//         Vmath::Vcopy(nq, &MF1st[direction][i * nq], 1, &tmp[0], 1);
//         m_fields[0]->ExtractTracePhys(tmp, tmptrace);
//         Vmath::Vvtvp(nTracePointsTot, m_traceNormals[i], 1, tmptrace, 1,
//                      traceVn, 1, traceVn, 1);
//     }

//     m_fields[0]->GetTrace()->Upwind(traceVn, Fwd, Bwd, flux);

//     Array<OneD, Array<OneD, NekDouble>> ncdotMFFwd;
//     Array<OneD, Array<OneD, NekDouble>> ncdotMFBwd;

//     ComputencdotMF(MF1st, ncdotMFFwd, ncdotMFBwd);

//     OutField = Array<OneD, NekDouble>(ncoeffs, 0.0);
//     // calculate numflux = (n \cdot MF)*flux
//     Vmath::Vmul(nTracePointsTot, &flux[0], 1, &ncdotMFFwd[direction][0],
//     1,
//                 &fluxFwd[0], 1);
//     Vmath::Vmul(nTracePointsTot, &flux[0], 1, &ncdotMFBwd[direction][0],
//     1,
//                 &fluxBwd[0], 1);

//     // FwdBwdtegral because generallize (N \cdot MF)_{FWD} \neq -(N \cdot
//     // MF)_{BWD}
//     Vmath::Neg(ncoeffs, WeakDeriv, 1);
//     m_fields[0]->AddFwdBwdTraceIntegral(fluxFwd, fluxBwd, WeakDeriv);
//     m_fields[0]->SetPhysState(false);

//     Vmath::Vadd(ncoeffs, &WeakDeriv[0], 1, &OutField[0], 1, &OutField[0],
//     1);
// }

void MMFNeuralEP::v_GenerateSummary(SolverUtils::SummaryList &s)
{
    int nvar = m_fields.size();

    MMFSystem::v_GenerateSummary(s);

    SolverUtils::AddSummaryItem(s, "NeuralEPType",
                                NeuralEPTypeMap[m_NeuralEPType]);
    SolverUtils::AddSummaryItem(s, "SolverSchemeType",
                                SolverSchemeTypeMap[m_SolverSchemeType]);
    SolverUtils::AddSummaryItem(s, "TimeMap", TimeMapTypeMap[m_TimeMap]);
    SolverUtils::AddSummaryItem(s, "TimeMapStart", m_TimeMapStart);
    SolverUtils::AddSummaryItem(s, "TimeMapEnd", m_TimeMapEnd);
    SolverUtils::AddSummaryItem(s, "ElemNodeEnd", m_ElemNodeEnd);
    SolverUtils::AddSummaryItem(s, "ElemMyelenEnd", m_ElemMyelenEnd);

    SolverUtils::AddSummaryItem(s, "Temperature", m_Temperature);
    SolverUtils::AddSummaryItem(s, "diameter", m_diameter);
    SolverUtils::AddSummaryItem(s, "Helmtau", m_Helmtau);
    if(nvar==1)
    {
        SolverUtils::AddSummaryItem(s, "ratio_re_ri", m_ratio_re_ri);
    }
}
} // namespace Nektar

int main(int argc, char *argv[])
{
    LibUtilities::SessionReaderSharedPtr session;
    SpatialDomains::MeshGraphSharedPtr graph;

    LibUtilities::SessionReaderSharedPtr session1D;
    SpatialDomains::MeshGraphSharedPtr graph1D;

    std::string vDriverModule;
    DriverSharedPtr drv;

    try
    {
        // Create session reader.
        session = LibUtilities::SessionReader::CreateInstance(argc, argv);

        // Create MeshGraph
        graph = SpatialDomains::MeshGraph::Read(session);

        // Create driver
        session->LoadSolverInfo("Driver", vDriverModule, "Standard");
        drv = GetDriverFactory().CreateInstance(vDriverModule, session, graph);

        // Execute driver
        drv->Execute();

        // Finalise session
        session->Finalise();
    }

    catch (const std::runtime_error &e)
    {
        return 1;
    }
    catch (const std::string &eStr)
    {
        std::cout << "Error: " << eStr << std::endl;
    }

    return 0;
}
