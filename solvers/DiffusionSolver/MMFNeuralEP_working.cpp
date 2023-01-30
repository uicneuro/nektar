///////////////////////////////////////////////////////////////////////////////
//
// File: MMFNeuralEP.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
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

#include <DiffusionSolver/EquationSystems/MMFNeuralEP.h>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>
#include <SolverUtils/Driver.h>

#include <boost/math/special_functions/spherical_harmonic.hpp>
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
    int nvar = m_fields.size();

    // Conductance parameters
    m_session->LoadParameter("Chi", m_chi, 28.0);
    m_session->LoadParameter("Cm", m_capMembrane, 0.125);

    // Helmsolver parameter
    m_session->LoadParameter("Helmtau", m_Helmtau, 1.0);

    // NeuralEP paramter on temperature
    m_session->LoadParameter("Temperature", m_Temperature, 24.0);
    m_session->LoadParameter("NumelemperNode", m_numelemperNode, 4);

    m_session->LoadParameter("TimeMapStart", m_TimeMapStart, 0.0);
    m_session->LoadParameter("TimeMapEnd", m_TimeMapEnd, 10000.0);

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



    // TimeMap ?
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
        case eNeuralEP1D:
        {
            std::string vNeuronModel;
            m_session->LoadSolverInfo("NEURONMODEL", vNeuronModel,
                                      "FrankenHuxley");

            ASSERTL0(vNeuronModel != "", "Neuron Model not specified.");

            m_neuron = GetNeuronModelFactory().CreateInstance(
                vNeuronModel, m_session, m_fields[0]);

            // Node and Myelen elements range
            m_session->LoadParameter("ElemNodeEnd", m_ElemNodeEnd, 0);
            m_session->LoadParameter("ElemMyelenEnd", m_ElemMyelenEnd, 0);

            // Relative Extracellular resistance: 1 < \beta < 10
            m_session->LoadParameter("ratio_re_ri", m_ratio_re_ri, 1.0);

            // Create capacitance vector
            // Parameter 1: m_Rnodelength: Ravier node elementwise length
            // Parameter 2: m_Rnodelength: Ravier node elementwise length
            m_session->LoadParameter("RanvierNodeLength", m_Rnodelength, 1);
            m_session->LoadParameter("RanvierNodeGap", m_Rnodegap, 11);

            // Ranvier node zone: 0: Myelin, 1: node
            m_nfibers  = 1;
            m_NodeZone = Array<OneD, Array<OneD, int>>(m_nfibers);
            for (int i = 0; i < m_nfibers; i++)
            {
                m_NodeZone[i] = Array<OneD, int>(nq);
            }

            m_NodeZone[0] = IndexNodeZone1D(m_fields[0], m_Rnodegap);

            // Constrct m_NeuralCm: node: 1/Cn, Myelin: 1/Cm
            const NekDouble Rf = m_neuron->GetRecistanceValue();
            const NekDouble Cm = m_neuron->GetCapacitanceValue(0);
            const NekDouble Cn = m_neuron->GetCapacitanceValue(1);

            std::cout << "Cm = " << Cm << ", Cn = " << Cn << ", Rf = " << Rf
                      << ", Cm * Rf = " << Cn * Rf << std::endl;

            m_NeuralCm    = Array<OneD, Array<OneD, NekDouble>>(1);
            m_NeuralCm[0] = Array<OneD, NekDouble>(nq);
            for (int i = 0; i < nq; ++i)
            {
                if (m_NodeZone[0][i] == 0)
                {
                    m_NeuralCm[0][i] = 1.0 / Cm;
                }

                else if (m_NodeZone[0][i] == 1)
                {
                    m_NeuralCm[0][i] = 1.0 / Cn;
                }
            }

            // Stimulus
            m_stimulus = Stimulus::LoadStimuli(m_session, m_fields[0]);
        }
        break;

        case eNeuralEP2D:
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
            m_NodeZone[0] =
                IndexNodeZone2D(m_fields[0], m_ElemNodeEnd, m_ElemMyelenEnd);

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

                    // Myelin zone
                    if (m_NodeZone[0][index] >= 0)

                    {
                        m_NeuralCm[0][index] = 1.0 / Cm;
                        cntm++;
                    }

                    // Ranvier node zone
                    else if (m_NodeZone[0][index] == -1)
                    {
                        m_NeuralCm[0][index] = 1.0 / Cn;
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
        }
        break;

        default:
            break;
    }

    // Derive AnisotropyStrength.
    switch (m_NeuralEPType)
    {
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
                    Vmath::Smul(nq, Cn, &m_NeuralCm[0][0], 1,
                                &AniStrength[j][0], 1);
                    Vmath::Vsqrt(nq, &AniStrength[j][0], 1, &AniStrength[j][0],
                                 1);
                }
            }

            std::cout << "Max Anistrength = "
                      << Vmath::Vmax(nq, AniStrength[0], 1)
                      << ", Min Anistrength = "
                      << Vmath::Vmin(nq, AniStrength[0], 1) << std::endl;

            MMFSystem::MMFInitObject(AniStrength);
        }
        break;

        case eNeuralEP2D:
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
                    Vmath::Smul(nq, Cn, &m_NeuralCm[0][0], 1,
                                &AniStrength[j][0], 1);
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

            // Check NodeZone and moving frames
            CheckNodeZoneMF(m_movingframes, m_NodeZone);

            // Set up for phie Poisson solver
            std::string phieMMFdirStr = "LOCAL";
            m_session->LoadSolverInfo("phieMMFDir", phieMMFdirStr, "LOCAL");

            Array<OneD, Array<OneD, NekDouble>> phieAniStrength(m_expdim);
            m_phieNeuralCm = Array<OneD, Array<OneD, NekDouble>>(m_expdim);
            for (int j = 0; j < m_expdim; ++j)
            {
                phieAniStrength[j] = Array<OneD, NekDouble>(nq, 1.0);
            }

            if (nvar == 2)
            {
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
            }
        }
        break;

        default:
        {
            Array<OneD, Array<OneD, NekDouble>> AniStrength(m_expdim);
            for (int j = 0; j < m_expdim; ++j)
            {
                AniStrength[j] = Array<OneD, NekDouble>(nq, 1.0);
            }
            MMFSystem::MMFInitObject(AniStrength);
        }
        break;
    }

    if (m_explicitDiffusion)
    {
        // m_ode.DefineImplicitSolve(&MMFNeuralEP::DoNullSolve, this);
        // m_ode.DefineProjection(&MMFNeuralEP::DoOdeProjection, this);
    }

    else
    {
        switch (m_NeuralEPType)
        {
            case eNeuralEP1D:
            {
                ComputeVarCoeff1D(m_movingframes, m_varcoeff);
            }
            break;

            case eNeuralEP2D:
            case eNeuralEP2DEmbed:
            {
                ComputeVarCoeff2D(m_movingframes, m_varcoeff);
                if (nvar == 2)
                {
                    ComputeVarCoeff2D(m_phiemovingframes, m_phievarcoeff);
                }
            }
            break;

            default:
                break;
        }

        // Test Helmsolve
        int ncoeffs          = GetNcoeffs();

        Array<OneD, NekDouble> testphys(nq,0.0);
        Array<OneD, NekDouble> testcoeffs(ncoeffs,0.0);

        StdRegions::ConstFactorMap factors;
        factors[StdRegions::eFactorTau] = 1.0;

        factors[StdRegions::eFactorLambda] = 1.0;

        std::cout << "HelmSolve test starts" << std::endl;
        m_fields[0]->HelmSolve(testphys, testcoeffs, factors, m_varcoeff);
        std::cout << "HelmSolve test ends, testcoeffs = " << RootMeanSquare(testcoeffs) << std::endl;

        switch (m_NeuralEPType)
        {
            case eNeuralEP1D:
            {
                m_ode.DefineImplicitSolve(
                    &MMFNeuralEP::DoImplicitSolveNeuralEP1D, this);
            }
            break;

            case eNeuralEP2D:
            {
                m_ode.DefineImplicitSolve(
                    &MMFNeuralEP::DoImplicitSolveNeuralEP2D, this);
            }
            break;

            // case eNeuralEP2DEmbed:
            // {
            //     m_ode.DefineImplicitSolve(
            //         &MMFNeuralEP::DoImplicitSolveNeuralEP2DEmbed, this);
            // }
            // break;

            default:
                break;
        }
    }

    switch (m_NeuralEPType)
    {
        case eNeuralEP1D:
        {
            m_ode.DefineOdeRhs(&MMFNeuralEP::DoOdeRhsNeuralEP1D, this);
        }
        break;

        case eNeuralEP2D:
        {
            m_ode.DefineOdeRhs(&MMFNeuralEP::DoOdeRhsNeuralEP2D, this);
        }
        break;

        // case eNeuralEP2DEmbed:
        // {
        //     m_ode.DefineOdeRhs(&MMFNeuralEP::DoOdeRhsNeuralEP2DEmbed, this);
        // }
        // break;

        default:
        break;
    }
}

void MMFNeuralEP::DoOdeRhsNeuralEP1D(
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

    for (unsigned int i = 0; i < m_stimulus.size(); ++i)
    {
        m_stimulus[i]->Update(RHSstimulus, time);
    }

    // ONLY simulation at node zone:
    // No excitation at myelinated region
    const NekDouble Cn = m_neuron->GetCapacitanceValue(1);
    for (int k = 0; k < nq; ++k)
    {
        // if( (m_NodeZone[0][k]>=0) && (m_NodeZone[0][k]<=1) )
        if (m_NodeZone[0][k] >= 0)
        {
            outarray[0][k] = outarray[0][k] + RHSstimulus[0][k] / Cn;
        }
    }

    // Multiply by 1/Cm for myeline or 1/Cm for node
    if (m_explicitDiffusion)
    {
        int nq = m_fields[0]->GetNpoints();

        // Laplacian only to the first variable
        Array<OneD, NekDouble> Laplacian(nq);
        // WeakDGMMFNeuralEP(0, inarray[0], Laplacian, time);

        Vmath::Vmul(nq, &m_NeuralCmRf[0], 1, &Laplacian[0], 1, &Laplacian[0],
                    1);
        Vmath::Vadd(nq, &Laplacian[0], 1, &outarray[0][0], 1, &outarray[0][0],
                    1);
    }
}

void MMFNeuralEP::DoOdeRhsNeuralEP2D(
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
                std::cout << "DoOdeRhs: HERE 5, nvar = " << nvar << std::endl;

    if (nvar == 2)
    {
        // Compute phi_e to satisfy the following equation
        // \nabla \cdot ( (\signa_e + \sigma_i) \nabla \phi_e) = - \nabla \cdot
        // (\sigma_i \nabla \phi_m)
        Updatephie(m_NodeZone[0], m_phiemovingframes, inarray[0]);
                std::cout << "DoOdeRhs: HERE 6" << std::endl;

        // Add the current changes by the external current
        Array<OneD, NekDouble> extcurrent;
        extcurrent = ComputeCovariantDiffusion(m_unitmovingframes,
                                               m_fields[1]->GetPhys());
        Vmath::Smul(nq, 1.0 / (Cn * Rf), extcurrent, 1, extcurrent, 1);
                std::cout << "DoOdeRhs: HERE 7" << std::endl;

        // Let the extcurrent be zero at Myeline nodes (-1).
        OnlyValideinNode(m_NodeZone[0], extcurrent);
                std::cout << "DoOdeRhs: HERE 8" << std::endl;

        // add divergence of phie to the current
        Vmath::Vadd(nq, &extcurrent[0], 1, &outarray[0][0], 1, &outarray[0][0],
                    1);
    }

    // Multiply by 1/Cm for myeline or 1/Cm for node
    if (m_explicitDiffusion)
    {
        int nq = m_fields[0]->GetNpoints();

        // Laplacian only to the first variable
        Array<OneD, NekDouble> Laplacian(nq);
        // WeakDGMMFNeuralEP(0, inarray[0], Laplacian, time);

        Vmath::Vmul(nq, &m_NeuralCmRf[0], 1, &Laplacian[0], 1, &Laplacian[0],
                    1);
        Vmath::Vadd(nq, &Laplacian[0], 1, &outarray[0][0], 1, &outarray[0][0],
                    1);
    }
}

// Implicit solve for NeuralEP solver
void MMFNeuralEP::DoImplicitSolveNeuralEP1D(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time,
    const NekDouble lambda)
{
    boost::ignore_unused(time);

    int nq = m_fields[0]->GetNpoints();

    const NekDouble R_f = m_neuron->GetRecistanceValue();
    const NekDouble C_n = m_neuron->GetCapacitanceValue(1);

    StdRegions::ConstFactorMap factors;
    factors[StdRegions::eFactorTau] = m_Helmtau;

    // // factors[StdRegions::eFactorLambda] = 1.0 / lambda;
    // // Cm dVm/dt = (1/rf) \nabla^2 Vm
    // // Vm^{n+1} = Vm^{n} + \Delta t/(rf*Cm) \nabla^2 Vn
    // // m_beta = Relative extracellular resistance = r_ex / r_f = \sigma_f /
    // \sigma_ex

    // NekDouble Cv = C_n * R_f * ( (m_beta_e + 1.0)/m_beta_e );
    // NekDouble Cv = C_n * R_f * (m_ratio_re_ri + 1.0);

    NekDouble Cv                       = C_n * R_f * m_ratio_re_ri;
    factors[StdRegions::eFactorLambda] = Cv / lambda;

    // We solve ( \nabla^2 - HHlambda ) Y[i] = rhs [i]
    // inarray = input: \hat{rhs} -> output: \hat{Y}
    // outarray = output: \hat{Y} where \hat = modal coeffs

    // For the variable of membrane potential: Multiply 1.0/timestep
    Vmath::Smul(nq, -factors[StdRegions::eFactorLambda], inarray[0], 1,
                m_fields[0]->UpdatePhys(), 1);

    SetBoundaryConditions(time);

    // Solve a system of equations with Helmholtz solver and transform
    // back into physical space.
    m_fields[0]->HelmSolve(m_fields[0]->GetPhys(), m_fields[0]->UpdateCoeffs(),
                           factors, m_varcoeff);

    m_fields[0]->BwdTrans(m_fields[0]->GetCoeffs(), outarray[0]);
    m_fields[0]->SetPhysState(true);
}

// Implicit solve for NeuralEP 2D solver
void MMFNeuralEP::DoImplicitSolveNeuralEP2D(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time,
    const NekDouble lambda)
{
    boost::ignore_unused(time);

    int nvar = m_fields.size();
    int nq   = m_fields[0]->GetNpoints();

    // Set up factors for Helmsolve
    const NekDouble R_f = m_neuron->GetRecistanceValue();
    const NekDouble C_n = m_neuron->GetCapacitanceValue(1);

    StdRegions::ConstFactorMap factors;
    factors[StdRegions::eFactorTau] = m_Helmtau;

    // NekDouble Cv = C_n * R_f * (m_ratio_re_ri + 0.0);
    // factors[StdRegions::eFactorLambda] = C_n * R_f / lambda ;
    NekDouble Cv = C_n * R_f;

    if (nvar == 1)
    {
        Cv = Cv * m_ratio_re_ri;
    }

    factors[StdRegions::eFactorLambda] = Cv / lambda;

    std::cout << "factors_Lambda = " << factors[StdRegions::eFactorLambda] << std::endl;

    SetBoundaryConditions(time);
    // SetMembraneBoundaryCondition(time);
                std::cout << "DoImplicit: HERE 4" << std::endl;

    // Multiply 1.0/timestep
    Vmath::Smul(nq, -factors[StdRegions::eFactorLambda], inarray[0], 1,
                m_fields[0]->UpdatePhys(), 1);
                std::cout << "DoImplicit: HERE 5" << std::endl;

    m_fields[0]->HelmSolve(m_fields[0]->GetPhys(), m_fields[0]->UpdateCoeffs(),
                           factors, m_varcoeff);

                std::cout << "DoImplicit: HERE 6" << std::endl;

    m_fields[0]->BwdTrans(m_fields[0]->GetCoeffs(), outarray[0]);
    m_fields[0]->SetPhysState(true);
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
void MMFNeuralEP::Updatephie(
    const Array<OneD, const int> &NodeZone,
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    const Array<OneD, const NekDouble> &phim)
{
    boost::ignore_unused(movingframes);

    int nq      = m_fields[0]->GetNpoints();
    int ncoeffs = GetNcoeffs();

    // Solve the Poisson equation: \nabla (\sigma_e + \sigma_i ) phi_e = \nabla
    // \sigma_i \nabla phi_m
    StdRegions::ConstFactorMap phiefactors;
    phiefactors[StdRegions::eFactorTau]    = m_Helmtau;
    phiefactors[StdRegions::eFactorLambda] = 0.0;

    // // Compute \nabla \sigma_i \nabla phi_m and use it as point sources for
    // phi_e.
    // // This is equivalently achieved by removing all the point sources in
    // myelinnated fiber region.
    Array<OneD, NekDouble> phimLaplacian(nq);
    phimLaplacian = ComputeCovariantDiffusion(m_unitmovingframes, phim);

    // Only nonzero for node.
    OnlyValideinNode(NodeZone, phimLaplacian);

    Vmath::Smul(nq, -1.0, phimLaplacian, 1, m_fields[1]->UpdatePhys(), 1);

    // Compute phie distribution
    SetMembraneBoundaryCondition();

    Array<OneD, NekDouble> tmpcoeff(ncoeffs);
    m_fields[1]->HelmSolve(m_fields[1]->GetPhys(), tmpcoeff, phiefactors,
                           m_phievarcoeff);
    m_fields[1]->BwdTrans(tmpcoeff, m_fields[1]->UpdatePhys());
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

// Constrcuct Cm vector: 1.0/Cn if node. 1.0/Cm if myeline.
Array<OneD, int> MMFNeuralEP::IndexNodeZone1D(
    const MultiRegions::ExpListSharedPtr &field, const int Nodegap)
{
    // int nq         = fields[0].size();
    int fnq = field->GetNpoints();

    Array<OneD, NekDouble> x0(fnq);
    Array<OneD, NekDouble> x1(fnq);
    Array<OneD, NekDouble> x2(fnq);

    field->GetCoords(x0, x1, x2);

    int index;
    NekDouble dx, dy, dz;

    Array<OneD, int> outarray(fnq, 0);
    Array<OneD, NekDouble> dist(m_fields[0]->GetExpSize(), 0.0);

    NekDouble npts,ioffset;
    for (int i = 0; i < m_fields[0]->GetExpSize(); ++i)
    {
        npts = m_fields[0]->GetTotPoints(i);
        ioffset = m_fields[0]->GetPhys_Offset(i);
        for (int j = 0; j < npts; ++j)
        {
            index = ioffset + j;

            // First and last element is all node for easier excitation
            if (i == 0 || i == 1 || ((i - 1) % Nodegap == 0) ||
                i == (m_fields[0]->GetExpSize() - 1) || i == (m_fields[0]->GetExpSize() - 2))
            {
                outarray[index] = 1;
            }
        }

        dx = x0[ioffset] - x0[ioffset+npts-1];
        dy = x1[ioffset] - x1[ioffset+npts-1];
        dz = x2[ioffset] - x2[ioffset+npts-1];

        dist[i] = sqrt(dx * dx + dy * dy + dz * dz);
    }

    return outarray;
}

// Constrcuct Cm vector: 1.0/Cn if node. 1.0/Cm if myeline.

Array<OneD, int> MMFNeuralEP::IndexNodeZone2D(
    const MultiRegions::ExpListSharedPtr &field, const int ElemNodeEnd,
    const int ElemMyelenEnd)
{
    // int nq         = fields[0].size();
    int fnq = field->GetNpoints();

    int index;
    int cntn = 0, cntm = 0, cnte = 0;

    // std::cout << "Nelemtj = " << Nelemtj
    // << ", ElemNodeEnd = " << ElemNodeEnd << ", ElemMyelenEnd = " <<
    // ElemMyelenEnd << std::endl;

    Array<OneD, int> outarray(fnq, 0);
    for (int i = 0; i < m_fields[0]->GetExpSize(); ++i)
    {
        for (int j = 0; j < m_fields[0]->GetTotPoints(i); ++j)
        {
            index = m_fields[0]->GetPhys_Offset(i) + j ;

            // First and last element is all node for easier excitation
            if (i <= ElemNodeEnd)
            {
                outarray[index] = i / m_numelemperNode;
                cntn++;
            }

            else if (i <= ElemMyelenEnd)
            {
                outarray[index] = -1;
                cntm++;
            }

            else
            {
                outarray[index] = -2;
                cnte++;
            }
        }
    }
}

void MMFNeuralEP::DoOdeProjection(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time)
{
    // Counter variable
    int i;
    int npoints    = GetNpoints();
    int nVariables = inarray.size();

    // Set the boundary conditions
    SetBoundaryConditions(time);

    // Switch on the projection type (Discontinuous or Continuous)
    for (i = 0; i < nVariables; ++i)
    {
        Vmath::Vcopy(npoints, inarray[i], 1, outarray[i], 1);
    }
}

void MMFNeuralEP::v_DoSolve()
{
    switch (m_SolverSchemeType)
    {
        case eMMFFirst:
        {
            DoSolveMMFFirst();
        }
        break;

        default:
        {
            DoSolveMMFZero();
        }
        break;
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
    // NekDouble frontloc=0.0, frontloc_old=0.0, time_old=0.0;
    while (step < m_steps || m_time < m_fintime - NekConstants::kNekZeroTol)
    {
        // field time integration
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
                      << std::endl;

            // fulltext.append("******************************");
            // fulltext.append("\n");
            // fulltext.append("Steps: " + std::to_string((step+1)));
            // fulltext.append("\n");
            // fulltext.append("Time: " + std::to_string(m_time));
            // fulltext.append("\n");

            // std::stringstream ss;
            // fulltext.append("CPU Time: " + ss.str());
            // fulltext.append("\n");

            cpuTime = 0.0;
        }

        // Write out checkpoint files
        if ((m_checksteps && step && !((step + 1) % m_checksteps)) ||
            doCheckTime)
        {
            // Array<OneD, NekDouble> dudt(nq);
            // Vmath::Vsub(nq, fields[0], 1, fields_old[0], 1, dudt, 1);
            // DisplayConductionVelocity(fields[0],dudt);

            // Print phim and phie at each node


            DisplayphiatNode(fields[0]);

            // Array<OneD, NekDouble> PoissonRHS(nq,0.0);
            // Array<OneD, NekDouble> PoissonLHS(nq,0.0);
            // PoissonLHS = ComputeCovariantDiffusion(m_phiemovingframes,
            // m_fields[1]->GetPhys());
            // // PoissonRHS = ComputeCovariantDiffusion(m_unitmovingframes,
            // m_fields[0]->GetPhys());

            // Array<OneD, NekDouble> Poissonerr(nq);
            // Vmath::Vadd(nq, PoissonRHS, 1, PoissonLHS, 1, Poissonerr, 1);

            // OnlyValideinNode(m_NodeZone[0], Poissonerr);

            // std::cout << "PoissonLHS = " << RootMeanSquare(PoissonLHS) << ",
            // PoissonRHS = " << RootMeanSquare(PoissonRHS)
            // << ", Poisson Error = " << RootMeanSquare(Poissonerr) <<
            // std::endl;

            Checkpoint_Output(nchk++);
            doCheckTime = false;
        }

        for (i = 0; i < 1; ++i)
        {
            Vmath::Vcopy(nq, &fields[i][0], 1, &fields_old[i][0], 1);
        }

        // Step advance
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

void MMFNeuralEP::DisplayphiatNode(const Array<OneD, const NekDouble> &field)
{
    int nq = GetTotPoints();

    // Print phim and phie at each node
    Array<OneD, NekDouble> phie(nq);
    Vmath::Vcopy(nq, m_fields[1]->GetPhys(), 1, phie, 1);

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

            locphimsum = locphimsum + field[index];
            locphiesum = locphiesum + phie[index];
        }

        phimavg[Rnodeid] = locphimsum / m_fields[0]->GetTotPoints(i);
        phieavg[Rnodeid] = locphiesum / m_fields[0]->GetTotPoints(i);
    }

    std::cout << " " << std::endl;
    std::cout << "(Nodeid,phim,phie): ";
    for (int i = 0; i < Rnodeid + 1; ++i)
    {
        std::cout << "(" << i << "," << phimavg[i] << "," << phieavg[i]
                  << "), ";
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
            ComputeTimeMap(m_time, fields[0], dudtval, m_ValidTimeMap,
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

            // Array<OneD, Array<OneD, NekDouble>> stimulusstrength(nvariables);
            // for (unsigned int i = 0; i < m_stimulus.size(); ++i)
            // {
            //     for (int j=0; j<nvariables; ++j)
            //     {
            //         stimulusstrength[j] = Array<OneD, NekDouble>(nq, 0.0);
            //     }

            //     m_stimulus[0]->Update(stimulusstrength, m_time);

            //     if(Vmath::Vmax(nq, stimulusstrength[0], 1)>0.01)
            //     {
            //         std::cout << " =================================== " <<
            //         std::endl; std::cout << "i = " << i << ", Stimulus: = "
            //         << Vmath::Vmax(nq, stimulusstrength[0], 1) << std::endl;
            //     std::cout << " =================================== " <<
            //     std::endl;
            //     }
            // }

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

                PlotTimeMapMF(NoBoundaryZone, TimeMap, TimeMapMF, MF1stAligned,
                              TMMFConnection, TMRelacc, nchk);

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


    // std::cout << "IndexNodeZone2D: cntn = " << cntn << ", cntm =" << cntm <<
    // ", cnte = " << cnte << std::endl;

    // Check Node Zone
    // int nq         = GetTotPoints();

    // Array<OneD, NekDouble> x0(nq);
    // Array<OneD, NekDouble> x1(nq);
    // Array<OneD, NekDouble> x2(nq);

    // m_fields[0]->GetCoords(x0, x1, x2);

    // NekDouble yavg, e1mag, e2mag;
    // for (int i = 0; i < Nelemtj; ++i)
    // {
    //     yavg=0.0;
    //     e1mag=0.0;
    //     e2mag=0.0;
    //     for (int j = 0; j < nptsj; ++j)
    //     {
    //         index = EWIndex[i][j];
    //         yavg = yavg + x1[index];
    //     }

    //     yavg = yavg/nptsj;
    //     if(outarray[index]>=0)
    //     {
    //         std::cout << "Elemid = " << i << ", Nodeid = " << outarray[index]
    //         << ", y = " << yavg << std::endl;
    //     }
    // }

    return outarray;
}

void MMFNeuralEP::CheckNodeZoneMF(
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    const Array<OneD, const Array<OneD, int>> &NodeZone)
{
    int nq = GetTotPoints();

    int i, j, index;
    NekDouble e1mag, e2mag;
    for (i = 0; i < m_fields[0]->GetExpSize(); ++i)
    {
        e1mag = 0.0;
        e2mag = 0.0;
        for (j = 0; j < m_fields[0]->GetTotPoints(i); ++j)
        {
            index = m_fields[0]->GetPhys_Offset(i) + j;
            e1mag = e1mag +
                    (movingframes[0][index] * movingframes[0][index] +
                     movingframes[0][nq + index] * movingframes[0][nq + index]);
            e2mag = e2mag +
                    (movingframes[1][index] * movingframes[1][index] +
                     movingframes[1][nq + index] * movingframes[1][nq + index]);
        }
        e1mag = sqrt(e1mag / m_fields[0]->GetTotPoints(i));
        e2mag = sqrt(e2mag / m_fields[0]->GetTotPoints(i));

        std::cout << "Elemid = " << i << ", Nodeid = " << NodeZone[0][index]
                  << ", e1mag = " << e1mag << ", e2mag = " << e1mag
                  << std::endl;
    }
}
// void MMFNeuralEP::v_InitObject(bool DeclareFields)
// {
//     UnsteadySystem::v_InitObject(DeclareFields);

//     int nq    = m_fields[0]->GetNpoints();
//     int nvar  = m_fields.size();
//     int MFdim = 3;


//     // Diffusivity coefficient for e^j
//     m_epsilon = Array<OneD, NekDouble>(MFdim);
//     m_session->LoadParameter("epsilon0", m_epsilon[0], 1.0);
//     m_session->LoadParameter("epsilon1", m_epsilon[1], 1.0);
//     m_session->LoadParameter("epsilon2", m_epsilon[2], 1.0);

//     // Diffusivity coefficient for u^j
//     m_epsu = Array<OneD, NekDouble>(nvar + 1);
//     m_session->LoadParameter("epsu0", m_epsu[0], 1.0);
//     m_session->LoadParameter("epsu1", m_epsu[1], 1.0);

//     m_session->LoadParameter("InitPtx", m_InitPtx, 0.0);
//     m_session->LoadParameter("InitPty", m_InitPty, 0.0);
//     m_session->LoadParameter("InitPtz", m_InitPtz, 0.0);

//     int shapedim = m_fields[0]->GetShapeDimension();
//     Array<OneD, Array<OneD, NekDouble>> Anisotropy(shapedim);
//     for (int j = 0; j < shapedim; ++j)
//     {
//         Anisotropy[j] = Array<OneD, NekDouble>(nq, 1.0);
//         Vmath::Fill(nq, sqrt(m_epsilon[j]), &Anisotropy[j][0], 1);
//     }

//     MMFSystem::MMFInitObject(Anisotropy);
//     ComputeVarCoeff2D(m_movingframes,m_varcoeff);


//     int ncoeffs = m_fields[0]->GetNcoeffs();

//         Array<OneD, NekDouble> testphys(nq,0.0);
//         Array<OneD, NekDouble> testcoeffs(ncoeffs,0.0);

//         StdRegions::ConstFactorMap factors;
//         factors[StdRegions::eFactorTau] = 1.0;

//         factors[StdRegions::eFactorLambda] = 1.0;

//             std::cout << "HelmSolve test 1:" << std::endl;
//         m_fields[0]->HelmSolve(testphys, testcoeffs, factors, m_varcoeff);
//         std::cout << "testcoeffs = " << RootMeanSquare(testcoeffs) << std::endl;


//     // Define ProblemType
//     if (m_session->DefinesSolverInfo("TESTTYPE"))
//     {
//         std::string TestTypeStr = m_session->GetSolverInfo("TESTTYPE");
//         int i;
//         for (i = 0; i < (int)SIZE_TestType; ++i)
//         {
//             if (boost::iequals(TestTypeMap[i], TestTypeStr))
//             {
//                 m_TestType = (TestType)i;
//                 break;
//             }
//         }
//     }
//     else
//     {
//         m_TestType = (TestType)0;
//     }

//     if (m_session->DefinesSolverInfo("INITWAVETYPE"))
//     {
//         std::string InitWaveTypeStr = m_session->GetSolverInfo("INITWAVETYPE");
//         for (int i = 0; i < (int)SIZE_TestType; ++i)
//         {
//             if (boost::iequals(InitWaveTypeMap[i], InitWaveTypeStr))
//             {
//                 m_InitWaveType = (InitWaveType)i;
//                 break;
//             }
//         }
//     }
//     else
//     {
//         m_InitWaveType = (InitWaveType)0;
//     }


//           std::cout << "HelmSolve test 2:" << std::endl;
//         m_fields[0]->HelmSolve(testphys, testcoeffs, factors, m_varcoeff);
//         std::cout << "testcoeffs = " << RootMeanSquare(testcoeffs) << std::endl;

//     if (!m_explicitDiffusion)
//     {
//         m_ode.DefineImplicitSolve(&MMFNeuralEP::DoImplicitSolve, this);
//     }

//             std::cout << "HelmSolve test 3:" << std::endl;
//         m_fields[0]->HelmSolve(testphys, testcoeffs, factors, m_varcoeff);
//         std::cout << "testcoeffs = " << RootMeanSquare(testcoeffs) << std::endl;

//     m_ode.DefineOdeRhs(&MMFNeuralEP::DoOdeRhs, this);
// }

/**
 *
 */
MMFNeuralEP::~MMFNeuralEP()
{
}

/**OdeRhs
 * @param   inarray         Input array.
 * @param   outarray        Output array.
 * @param   time            Current simulation time.
 * @param   lambda          Timestep.
 */
void MMFNeuralEP::DoImplicitSolve(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time,
    const NekDouble lambda)
{
    int nvariables = inarray.size();
    int nq         = m_fields[0]->GetNpoints();

    StdRegions::ConstFactorMap factors;
    factors[StdRegions::eFactorTau] = 1.0;

    Array<OneD, Array<OneD, NekDouble>> F(nvariables);
    factors[StdRegions::eFactorLambda] = 1.0 / lambda;
    F[0] = Array<OneD, NekDouble>(nq * nvariables);
    for (int n = 1; n < nvariables; ++n)
    {
        F[n] = F[n - 1] + nq;
    }

    // We solve ( \nabla^2 - HHlambda ) Y[i] = rhs [i]
    // inarray = input: \hat{rhs} -> output: \hat{Y}
    // outarray = output: nabla^2 \hat{Y}
    // where \hat = modal coeffs
    SetBoundaryConditions(time);

    for (int i = 0; i < nvariables; ++i)
    {
        factors[StdRegions::eFactorLambda] = 1.0 / lambda / m_epsu[i];

        // Multiply 1.0/timestep
        Vmath::Smul(nq, -factors[StdRegions::eFactorLambda], inarray[i], 1,
                    F[i], 1);
        m_fields[i]->HelmSolve(F[i], m_fields[i]->UpdateCoeffs(), factors,
                               m_varDiffcoeff);
        m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(), outarray[i]);
    }
}

/**
 *
 */
void MMFNeuralEP::DoOdeRhs(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time)
{
    int nq = GetTotPoints();

    switch (m_TestType)
    {
        case eTestPlane:
        {

            Array<OneD, NekDouble> x(nq);
            Array<OneD, NekDouble> y(nq);
            Array<OneD, NekDouble> z(nq);

            m_fields[0]->GetCoords(x, y, z);

            for (int k = 0; k < nq; k++)
            {
                outarray[0][k] = (m_epsilon[0] + m_epsilon[1] - 1.0) * m_pi *
                                 m_pi * exp(-1.0 * m_pi * m_pi * time) *
                                 sin(m_pi * x[k]) * cos(m_pi * y[k]);
            }
        }
        break;

        case eTestCube:
        {

            Array<OneD, NekDouble> x(nq);
            Array<OneD, NekDouble> y(nq);
            Array<OneD, NekDouble> z(nq);

            m_fields[0]->GetCoords(x, y, z);

            for (int k = 0; k < nq; k++)
            {
                outarray[0][k] =
                    (m_epsilon[0] + m_epsilon[1] + m_epsilon[2] - 1.0) * m_pi *
                    m_pi * exp(-1.0 * m_pi * m_pi * time) * sin(m_pi * x[k]) *
                    sin(m_pi * y[k]) * sin(m_pi * z[k]);
            }
        }
        break;

        case eTestLinearSphere:
        {
            Array<OneD, NekDouble> temp(nq);

            NekDouble A = 2.0;
            NekDouble B = 5.0;

            NekDouble m_a, m_b, m_c, m_d;
            m_a = B - 1.0;
            m_b = A * A;
            m_c = -1.0 * B;
            m_d = -1.0 * A * A;

            temp = Array<OneD, NekDouble>(nq, 0.0);
            Vmath::Svtvp(nq, m_a, &inarray[0][0], 1, &temp[0], 1, &temp[0], 1);
            Vmath::Svtvp(nq, m_b, &inarray[1][0], 1, &temp[0], 1,
                         &outarray[0][0], 1);

            temp = Array<OneD, NekDouble>(nq, 0.0);
            Vmath::Svtvp(nq, m_c, &inarray[0][0], 1, &temp[0], 1, &temp[0], 1);
            Vmath::Svtvp(nq, m_d, &inarray[1][0], 1, &temp[0], 1,
                         &outarray[1][0], 1);
        }
        break;

        case eTestNonlinearSphere:
        {
            NekDouble A = 2.0;
            NekDouble B = 5.0;

            Array<OneD, NekDouble> Aonevec(nq, A);

            // cube = phys0*phys0*phy1
            Array<OneD, NekDouble> cube(nq);
            Vmath::Vmul(nq, &inarray[0][0], 1, &inarray[0][0], 1, &cube[0], 1);
            Vmath::Vmul(nq, &inarray[1][0], 1, &cube[0], 1, &cube[0], 1);

            // outarray[0] = A - B*phy0 + phy0*phy0*phy1 - phy0
            NekDouble coeff = -1.0 * B - 1.0;
            Array<OneD, NekDouble> tmp(nq);
            Vmath::Svtvp(nq, coeff, &inarray[0][0], 1, &cube[0], 1, &tmp[0], 1);
            Vmath::Vadd(nq, &Aonevec[0], 1, &tmp[0], 1, &outarray[0][0], 1);

            // outarray[1] = B*phys0 - phy0*phy0*phy1
            Vmath::Svtvm(nq, B, &inarray[0][0], 1, &cube[0], 1, &outarray[1][0],
                         1);
        }
        break;

            // case eFHNStandard:
            // {
            //     // \phi - \phi^3/3 - \psi
            //     NekDouble a  = 0.12;
            //     NekDouble b  = 0.011;
            //     NekDouble c1 = 0.175;
            //     NekDouble c2 = 0.03;
            //     NekDouble d  = 0.55;

            //     Array<OneD, NekDouble> tmp(nq);

            //     // Reaction for \phi = c1 \phi ( \phi - a)*(1 - \phi) - c2 v
            //     Vmath::Smul(nq, -1.0 * c1, inarray[0], 1, outarray[0], 1);
            //     Vmath::Sadd(nq, -1.0 * a, inarray[0], 1, tmp, 1);
            //     Vmath::Vmul(nq, tmp, 1, inarray[0], 1, outarray[0], 1);
            //     Vmath::Sadd(nq, -1.0, inarray[0], 1, tmp, 1);
            //     Vmath::Vmul(nq, tmp, 1, outarray[0], 1, outarray[0], 1);

            //     Vmath::Smul(nq, -1.0 * c2, inarray[1], 1, tmp, 1);
            //     Vmath::Vadd(nq, tmp, 1, outarray[0], 1, outarray[0], 1);

            //     // Reaction for \psi = b (\phi - d \psi )
            //     Vmath::Svtvp(nq, -1.0 * d, inarray[1], 1, inarray[0], 1,
            //                  outarray[1], 1);
            //     Vmath::Smul(nq, b, outarray[1], 1, outarray[1], 1);
            // }
            // break;

            // case eFHNRogers:
            // {
            //     NekDouble a  = 0.13;
            //     NekDouble b  = 0.013;
            //     NekDouble c1 = 0.26;
            //     NekDouble c2 = 0.1;
            //     NekDouble d  = 1.0;

            //     Array<OneD, NekDouble> tmp(nq);

            //     // Reaction for \phi = c1 \phi ( \phi - a)*(1 - \phi) - c2 u
            //     v Vmath::Smul(nq, -1.0 * c1, inarray[0], 1, outarray[0], 1);
            //     Vmath::Sadd(nq, -1.0 * a, inarray[0], 1, tmp, 1);
            //     Vmath::Vmul(nq, tmp, 1, outarray[0], 1, outarray[0], 1);
            //     Vmath::Sadd(nq, -1.0, inarray[0], 1, tmp, 1);
            //     Vmath::Vmul(nq, tmp, 1, outarray[0], 1, outarray[0], 1);

            //     Vmath::Vmul(nq, inarray[0], 1, inarray[1], 1, tmp, 1);
            //     Vmath::Smul(nq, -1.0 * c2, tmp, 1, tmp, 1);
            //     Vmath::Vadd(nq, tmp, 1, outarray[0], 1, outarray[0], 1);

            //     // Reaction for \psi = b (\phi - d \psi )
            //     Vmath::Svtvp(nq, -1.0 * d, inarray[1], 1, inarray[0], 1,
            //                  outarray[1], 1);
            //     Vmath::Smul(nq, b, outarray[1], 1, outarray[1], 1);
            // }
            // break;

            // case eFHNAlievPanf:
            // {

            //     NekDouble a   = 0.15;
            //     NekDouble c1  = 8.0;
            //     NekDouble c2  = 1.0;
            //     NekDouble c0  = 0.002;
            //     NekDouble mu1 = 0.2;
            //     NekDouble mu2 = 0.3;

            //     Array<OneD, NekDouble> tmp(nq);

            //     // Reaction for \phi = c1 \phi ( \phi - a)*(1 - \phi) - c2 u
            //     v Vmath::Smul(nq, -1.0 * c1, inarray[0], 1, outarray[0], 1);
            //     Vmath::Sadd(nq, -1.0 * a, inarray[0], 1, tmp, 1);
            //     Vmath::Vmul(nq, tmp, 1, outarray[0], 1, outarray[0], 1);
            //     Vmath::Sadd(nq, -1.0, inarray[0], 1, tmp, 1);
            //     Vmath::Vmul(nq, tmp, 1, outarray[0], 1, outarray[0], 1);

            //     Vmath::Vmul(nq, inarray[0], 1, inarray[1], 1, tmp, 1);
            //     Vmath::Smul(nq, -1.0 * c2, tmp, 1, tmp, 1);
            //     Vmath::Vadd(nq, tmp, 1, outarray[0], 1, outarray[0], 1);

            //     // Reaction for \psi = (c0 + (\mu1 \psi/(\mu2+\phi) )
            //     )*(-\psi - c1
            //     // * \phi*(\phi - a - 1) )

            //     Vmath::Smul(nq, mu1, inarray[1], 1, outarray[1], 1);
            //     Vmath::Sadd(nq, mu2, inarray[0], 1, tmp, 1);
            //     Vmath::Vdiv(nq, outarray[1], 1, tmp, 1, outarray[1], 1);
            //     Vmath::Sadd(nq, c0, outarray[1], 1, outarray[1], 1);

            //     Vmath::Sadd(nq, (-a - 1.0), inarray[0], 1, tmp, 1);
            //     Vmath::Vmul(nq, inarray[0], 1, tmp, 1, tmp, 1);
            //     Vmath::Smul(nq, c1, tmp, 1, tmp, 1);
            //     Vmath::Vadd(nq, inarray[1], 1, tmp, 1, tmp, 1);
            //     Vmath::Neg(nq, tmp, 1);

            //     Vmath::Vmul(nq, tmp, 1, outarray[1], 1, outarray[1], 1);
            // }
            // break;

        default:
            break;
    }
}

/**
 *
 */
void MMFNeuralEP::v_SetInitialConditions(NekDouble initialtime,
                                          bool dumpInitialConditions,
                                          const int domain)
{
    boost::ignore_unused(domain);

    int nq = GetTotPoints();

    EquationSystem::v_SetInitialConditions(initialtime, false);

    switch (m_TestType)
    {
        case eTestPlane:
        {
            Array<OneD, NekDouble> u(nq);

            TestPlaneProblem(initialtime, u);
            m_fields[0]->SetPhys(u);
        }
        break;

        case eTestCube:
        {
            Array<OneD, NekDouble> u(nq);

            TestCubeProblem(initialtime, u);
            m_fields[0]->SetPhys(u);
            /*for (int k=0; k<nq; ++k)
            {
                //for (int j=0; j<m_spacedim; ++j)
                //{
                cout << "_varcoeff" << u[k] <<endl;
                // }
            }*/
        }
        break;

        case eTestLinearSphere:
        case eTestNonlinearSphere:
        {
            Array<OneD, NekDouble> u(nq);
            Array<OneD, NekDouble> v(nq);

            Morphogenesis(initialtime, 0, u);
            Morphogenesis(initialtime, 1, v);

            m_fields[0]->SetPhys(u);
            m_fields[1]->SetPhys(v);
        }
        break;

            // case eFHNStandard:
            // case eFHNRogers:
            // case eFHNAlievPanf:
            // {
            //     Array<OneD, NekDouble> Zero(nq, 0.0);
            //     m_fields[0]->SetPhys(PlanePhiWave());
            //     m_fields[1]->SetPhys(Zero);
            // }
            // break;

        default:
        {
            EquationSystem::v_SetInitialConditions(initialtime, false);
        }
        break;
    }

    // forward transform to fill the modal coeffs
    // for (int i = 0; i < m_fields.size(); ++i)
    // {
    //     m_fields[i]->SetPhysState(true);
    //     m_fields[i]->FwdTrans(m_fields[i]->GetPhys(),
    //                           m_fields[i]->UpdateCoeffs());
    // }

    if (dumpInitialConditions)
    {
        std::string outname = m_sessionName + "_initial.chk";
        WriteFld(outname);
    }
}

void MMFNeuralEP::TestPlaneProblem(const NekDouble time,
                                    Array<OneD, NekDouble> &outfield)

{
    int nq = GetTotPoints();

    Array<OneD, NekDouble> x(nq);
    Array<OneD, NekDouble> y(nq);
    Array<OneD, NekDouble> z(nq);

    m_fields[0]->GetCoords(x, y, z);

    outfield = Array<OneD, NekDouble>(nq);
    for (int k = 0; k < nq; k++)
    {
        outfield[k] = exp(-1.0 * m_pi * m_pi * time) * sin(m_pi * x[k]) *
                      cos(m_pi * y[k]);
    }
}

void MMFNeuralEP::TestCubeProblem(const NekDouble time,
                                   Array<OneD, NekDouble> &outfield)

{
    int nq = GetTotPoints();

    Array<OneD, NekDouble> x(nq);
    Array<OneD, NekDouble> y(nq);
    Array<OneD, NekDouble> z(nq);

    m_fields[0]->GetCoords(x, y, z);

    outfield = Array<OneD, NekDouble>(nq);
    for (int k = 0; k < nq; k++)
    {
        outfield[k] = exp(-1.0 * m_pi * m_pi * time) * sin(m_pi * x[k]) *
                      sin(m_pi * y[k]) * sin(m_pi * z[k]);
    }
}

void MMFNeuralEP::Morphogenesis(const NekDouble time, unsigned int field,
                                 Array<OneD, NekDouble> &outfield)
{
    int nq = GetTotPoints();

    int i, m, n, ind;
    NekDouble a_n, d_n, gamma_n;
    NekDouble A_mn, C_mn, theta, phi, radius;

    std::complex<double> Spericharmonic, delta_n, temp;
    std::complex<double> varphi0, varphi1;
    std::complex<double> B_mn, D_mn;

    // Set some parameter values
    int Maxn = 6;
    int Maxm = 2 * Maxn - 1;

    NekDouble A = 2.0;
    NekDouble B = 5.0;

    NekDouble m_mu = 0.001;
    NekDouble m_nu = 0.002;

    NekDouble m_a, m_b, m_c, m_d;

    m_a = B - 1.0;
    m_b = A * A;
    m_c = -1.0 * B;
    m_d = -1.0 * A * A;

    Array<OneD, Array<OneD, NekDouble>> Ainit(Maxn);
    Array<OneD, Array<OneD, NekDouble>> Binit(Maxn);

    for (i = 0; i < Maxn; ++i)
    {
        Ainit[i] = Array<OneD, NekDouble>(Maxm, 0.0);
        Binit[i] = Array<OneD, NekDouble>(Maxm, 0.0);
    }

    Ainit[5][0]  = -0.5839;
    Ainit[5][1]  = -0.8436;
    Ainit[5][2]  = -0.4764;
    Ainit[5][3]  = 0.6475;
    Ainit[5][4]  = 0.1886;
    Ainit[5][5]  = 0.8709;
    Ainit[5][6]  = -0.8338;
    Ainit[5][7]  = 0.1795;
    Ainit[5][8]  = -0.7873;
    Ainit[5][9]  = 0.8842;
    Ainit[5][10] = 0.2943;

    Binit[5][0]  = -0.6263;
    Binit[5][1]  = 0.9803;
    Binit[5][2]  = 0.7222;
    Binit[5][3]  = 0.5945;
    Binit[5][4]  = 0.6026;
    Binit[5][5]  = -0.2076;
    Binit[5][6]  = 0.4556;
    Binit[5][7]  = 0.6024;
    Binit[5][8]  = 0.9695;
    Binit[5][9]  = -0.4936;
    Binit[5][10] = 0.1098;

    Array<OneD, NekDouble> u(nq);
    Array<OneD, NekDouble> v(nq);
    Array<OneD, NekDouble> x(nq);
    Array<OneD, NekDouble> y(nq);
    Array<OneD, NekDouble> z(nq);

    m_fields[0]->GetCoords(x, y, z);
    for (int i = 0; i < nq; ++i)
    {
        radius = sqrt(x[i] * x[i] + y[i] * y[i] + z[i] * z[i]);

        // theta is in [0, pi]
        theta = asin(z[i] / radius) + 0.5 * m_pi;

        // phi is in [0, 2*pi]
        phi = atan2(y[i], x[i]) + m_pi;

        varphi0 = 0.0 * varphi0;
        varphi1 = 0.0 * varphi1;
        for (n = 0; n < Maxn; ++n)
        {
            // Set up parameters
            a_n = m_a - m_mu * (n * (n + 1) / radius / radius);
            d_n = m_d - m_nu * (n * (n + 1) / radius / radius);

            gamma_n = 0.5 * (a_n + d_n);

            temp    = (a_n + d_n) * (a_n + d_n) - 4.0 * (a_n * d_n - m_b * m_c);
            delta_n = 0.5 * sqrt(temp);

            for (m = -n; m <= n; ++m)
            {
                ind  = m + n;
                A_mn = Ainit[n][ind];
                C_mn = Binit[n][ind];

                B_mn = ((a_n - gamma_n) * Ainit[n][ind] + m_b * Binit[n][ind]) /
                       delta_n;
                D_mn = (m_c * Ainit[n][ind] + (d_n - gamma_n) * Binit[n][ind]) /
                       delta_n;

                Spericharmonic =
                    boost::math::spherical_harmonic(n, m, theta, phi);
                varphi0 += exp(gamma_n * time) *
                           (A_mn * cosh(delta_n * time) +
                            B_mn * sinh(delta_n * time)) *
                           Spericharmonic;
                varphi1 += exp(gamma_n * time) *
                           (C_mn * cosh(delta_n * time) +
                            D_mn * sinh(delta_n * time)) *
                           Spericharmonic;
            }
        }

        u[i] = varphi0.real();
        v[i] = varphi1.real();
    }

    switch (field)
    {
        case 0:
        {
            outfield = u;
        }
        break;

        case 1:
        {
            outfield = v;
        }
        break;
    }
}

Array<OneD, NekDouble> MMFNeuralEP::PlanePhiWave()
{
    int nq = GetTotPoints();
    Array<OneD, NekDouble> outarray(nq, 0.0);

    Array<OneD, NekDouble> x(nq);
    Array<OneD, NekDouble> y(nq);
    Array<OneD, NekDouble> z(nq);

    m_fields[0]->GetCoords(x, y, z);

    NekDouble xmin, ymin, xmax;

    xmin = Vmath::Vmin(nq, x, 1);
    xmax = Vmath::Vmax(nq, x, 1);
    ymin = Vmath::Vmin(nq, y, 1);

    NekDouble xp, yp, xp2;
    for (int i = 0; i < nq; i++)
    {
        switch (m_InitWaveType)
        {
            case eLeft:
            {
                NekDouble radiusofinit = 4.0;
                NekDouble frontstiff   = 0.1;

                xp = x[i] - xmin;
                outarray[i] =
                    1.0 / (1.0 + exp((xp - radiusofinit) / frontstiff));
            }
            break;

            case eBothEnds:
            {
                NekDouble radiusofinit = 3.0;
                NekDouble frontstiff   = 0.1;

                xp  = x[i] - xmin;
                xp2 = x[i] - xmax;

                outarray[i] =
                    1.0 / (1.0 +
                           exp((sqrt(xp * xp) - radiusofinit) / frontstiff)) +
                    1.0 / (1.0 +
                           exp((sqrt(xp2 * xp2) - radiusofinit) / frontstiff));
            }
            break;

            case eCenter:
            {
                NekDouble radiusofinit = 6.0;
                NekDouble frontstiff   = 0.1;

                // NekDouble xc = 0.5*(Vmath::Vmax(nq, x, 1) + Vmath::Vmin(nq,
                // x, 1));

                xp = x[i] - xmin;
                outarray[i] =
                    1.0 / (1.0 + exp((xp - radiusofinit) / frontstiff));
            }
            break;

            case eLeftBottomCorner:
            {
                NekDouble radiusofinit = 6.0;
                NekDouble frontstiff   = 0.1;
                NekDouble bs           = 2.0;

                xp = x[i] - xmin;
                yp = y[i] - ymin;
                outarray[i] =
                    1.0 /
                    (1.0 + exp((sqrt(xp * xp + yp * yp) / bs - radiusofinit) /
                               frontstiff));
            }
            break;

            case ePoint:
            {
                NekDouble xloc, yloc, zloc, rad;
                NekDouble radiusofinit = 10.0;

                xloc = x[i] - m_InitPtx;
                yloc = y[i] - m_InitPty;
                zloc = z[i] - m_InitPtz;

                rad = sqrt(xloc * xloc + yloc * yloc + zloc * zloc);

                xloc = xloc / radiusofinit;
                yloc = yloc / radiusofinit;
                zloc = zloc / radiusofinit;

                if (rad < radiusofinit)
                {
                    outarray[i] =
                        exp(-(1.0 / 2.0) *
                            (xloc * xloc + yloc * yloc + zloc * zloc));
                }

                else
                {
                    outarray[i] = 0.0;
                }
            }
            break;

            case eSpiralDock:
            {
                NekDouble radiusofinit = 3.0;
                NekDouble frontstiff   = 0.1;
                xp                     = x[i] - 4.0;
                yp                     = y[i];
                outarray[i] =
                    (1.0 / (1.0 + exp(2.0 * yp))) *
                    (1.0 / (1.0 + exp(-2.0 * xp))) *
                    (1.0 / (1.0 + exp((xp - radiusofinit) / frontstiff)));
            }
            break;

            default:
                break;
        }
    }

    return outarray;
}

void MMFNeuralEP::v_EvaluateExactSolution(unsigned int field,
                                           Array<OneD, NekDouble> &outfield,
                                           const NekDouble time)
{
    switch (m_TestType)
    {
        case eTestPlane:
        {
            TestPlaneProblem(time, outfield);
        }
        break;

        case eTestCube:
        {
            TestCubeProblem(time, outfield);
        }
        break;

        case eTestLinearSphere:
        case eTestNonlinearSphere:
        {
            Morphogenesis(time, field, outfield);
        }
        break;

            // case eFHNStandard:
            // case eFHNRogers:
            // case eFHNAlievPanf:
            // {
            //     int nq   = GetTotPoints();
            //     outfield = Array<OneD, NekDouble>(nq, 0.0);
            // }
            /* Falls through. */
        default:
        {
            EquationSystem::v_EvaluateExactSolution(field, outfield, time);
        }
        break;
    }
}

void MMFNeuralEP::ComputeEuclideanDivMF(
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    Array<OneD, Array<OneD, NekDouble>> &DivMF)
{
    int nq = m_fields[0]->GetNpoints();

    DivMF = Array<OneD, Array<OneD, NekDouble>>(m_expdim);
    for (int j = 0; j < m_expdim; ++j)
    {
        DivMF[j] = Array<OneD, NekDouble>(nq, 0.0);
    }

    Array<OneD, NekDouble> tmp(nq);
    Array<OneD, NekDouble> Dtmp(nq);

    // case eEuclidean:
    for (int j = 0; j < m_expdim; ++j)
    {
        for (int k = 0; k < m_spacedim; ++k)
        {
            Vmath::Vcopy(nq, &movingframes[j][k * nq], 1, &tmp[0], 1);
            m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[k], tmp, Dtmp);
            Vmath::Vadd(nq, &Dtmp[0], 1, &DivMF[j][0], 1, &DivMF[j][0], 1);
        }
    }
}

void MMFNeuralEP::v_GenerateSummary(SolverUtils::SummaryList &s)
{
    MMFSystem::v_GenerateSummary(s);
    SolverUtils::AddSummaryItem(s, "TestType", TestTypeMap[m_TestType]);
    SolverUtils::AddSummaryItem(s, "epsilon0", m_epsilon[0]);
    SolverUtils::AddSummaryItem(s, "epsilon1", m_epsilon[1]);
    SolverUtils::AddSummaryItem(s, "epsilon2", m_epsilon[2]);
    if (m_TestType == eTestLinearSphere)
    {
        SolverUtils::AddSummaryItem(s, "epsilon for u", m_epsu[0]);
        SolverUtils::AddSummaryItem(s, "epsilon for v", m_epsu[1]);
    }
}
 // namespace Nektar
int main(int argc, char *argv[])
{
    LibUtilities::SessionReaderSharedPtr session;
    SpatialDomains::MeshGraphSharedPtr graph;
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