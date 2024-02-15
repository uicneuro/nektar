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
#include <math.h>

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

    // Derive AnisotropyStrength.
    m_AniStrength = Array<OneD, Array<OneD, NekDouble>> (m_expdim);
    m_phieAniStrength = Array<OneD, Array<OneD, NekDouble>> (m_expdim);
    for (int j = 0; j < m_mfdim; ++j)
    {
        m_AniStrength[j] = Array<OneD, NekDouble>(nq, 1.0);
        m_phieAniStrength[j] = Array<OneD, NekDouble>(nq, 1.0);
    }
    
    m_TimeMap = Array<OneD, Array<OneD, NekDouble>>(1);
    for (int i = 0; i < 1; ++i)
    {
        m_TimeMap[i] = Array<OneD, NekDouble>(nq,0.0);
    }

    // Conductance parameters
    m_session->LoadParameter("Chi", m_chi, 28.0);
    m_session->LoadParameter("Cm", m_capMembrane, 0.125);

    // Helmsolver parameter
    m_session->LoadParameter("Helmtau", m_Helmtau, 1.0);
    
    // Resting potential
    m_session->LoadParameter("urest", m_urest, 0.0);

    // NeuralEP paramter on temperature
    m_session->LoadParameter("Temperature", m_Temperature, 24.0);
    m_session->LoadParameter("diameter", m_diameter, 0.001);

    m_session->LoadParameter("FiberLength", m_fiberlen, 0.01);
    m_session->LoadParameter("NodeLength", m_nodelen, 0.01);
    m_session->LoadParameter("MyelinLength", m_myelinlen, 0.2);

    m_session->LoadParameter("Total_Number_Node", m_totNode, 3);
    m_session->LoadParameter("Element_per_Node", m_elemperNode, 4);
    m_session->LoadParameter("Element_per_Myelin", m_elemperMyel, 14);

    m_session->LoadParameter("AnisotropyStrength", m_AnisotropyStrength, 4.0);

    // Relative Extracellular resistance: 1 < \beta < 10
    m_session->LoadParameter("ratio_re_ri", m_ratio_re_ri, 1.0);

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

    // Define FiberType
    if (m_session->DefinesSolverInfo("FIBERTYPE"))
    {
        std::string FiberTypeStr;
        FiberTypeStr = m_session->GetSolverInfo("FIBERTYPE");
        for (int i = 0; i < (int)SIZE_FiberType; ++i)
        {
            if (FiberTypeMap[i] == FiberTypeStr)
            {
                m_fiberType = (FiberType)i;
                break;
            }
        }
    }
    else
    {
        m_fiberType = (FiberType)0;
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
                m_TimeMapScheme = (TimeMapType)i;
                break;
            }
        }
    }
    else
    {
        m_TimeMapScheme = (TimeMapType)0;
    }

    // Either incorporating external current effect or not.
    if (m_session->DefinesSolverInfo("ExtCurrentType"))
    {
        std::string ExtCurrentTypeStr;
        ExtCurrentTypeStr = m_session->GetSolverInfo("ExtCurrentType");
        for (int i = 0; i < (int)SIZE_ExtCurrentType; ++i)
        {
            if (boost::iequals(ExtCurrentTypeMap[i], ExtCurrentTypeStr))
            {
                m_ExtCurrentType = (ExtCurrentType)i;
                break;
            }
        }
    }
    else
    {
        m_ExtCurrentType = (ExtCurrentType)0;
    }

    std::string vNeuronModel;
    m_session->LoadSolverInfo("NEURONMODEL", vNeuronModel,
                                "FrankenHuxley");

    ASSERTL0(vNeuronModel != "", "Neuron Model not specified.");

    m_neuron = GetNeuronModelFactory().CreateInstance(
        vNeuronModel, m_session, m_fields[0]);

    // Rf and Cn are imported
    m_Rf = m_neuron->GetRecistanceValue();
    m_Cm = m_neuron->GetCapacitanceValue(0);
    m_Cn = m_neuron->GetCapacitanceValue(1);

   switch (m_NeuralEPType)
    {
        // case eNeuralEPPT:
        // {
        //     // Ranvier node zone: 0: Myelin, 1: node
        //     m_nfibers  = 1;
        //     m_zoneindex = Array<OneD, Array<OneD, int>>(m_nfibers);
        //     m_zoneindex[0] = Array<OneD, int>(nq, 1);

        //     m_NeuralCm    = Array<OneD, Array<OneD, NekDouble>>(m_nfibers);
        //     m_NeuralCm[0] = ComputeConductivity(m_zoneindex[0]);
        //     break;
        // }

        case eNeuralEP1D:
        {
            // Ranvier node zone: 0: Myelin, 1: node
            m_nfibers  = 1;
            m_zoneindex = Array<OneD, Array<OneD, int>>(m_nfibers);
            m_zoneindex[0] = IndexNodeZone1D(m_fields[0], m_totNode, m_elemperNode, m_elemperMyel); 

              // Get the first and last index of the excitation zone [1,2]
            SetUpDomainZone(m_zoneindex[0], m_excitezone, m_nodezone, m_intrazone, m_extrazone);

            m_NeuralCm    = Array<OneD, Array<OneD, NekDouble>>(m_nfibers);
            m_NeuralCm[0] = ComputeConductivity(m_zoneindex[0]);
            break;
        }

        case eNeuralHelmTest:
        case eNeuralEP2Dmono:
        case eNeuralEP2Dbi:
        {
            // load zoneindexfile
            m_session->LoadSolverInfo("zoneindexfile", m_zoneindexfile, "Null");
            
            m_zoneindex = Array<OneD, Array<OneD, int>>(1);
            m_zoneindex[0]  = Array<OneD, int>(nq, 1); 

            m_npts = m_fields[0]->GetTotPoints(0);
            // m_zoneindex[0] = IndexNodeZone2D(m_fiberlen, m_nodelen, m_myelinlen, m_totNode);
            m_zoneindex[0] = IndexNodeZone2D(m_fiberType);

            // Get the first and last index of the excitation zone [1,2]
            SetUpDomainZone(m_zoneindex[0], m_excitezone, m_nodezone, m_intrazone, m_extrazone);

            if(m_MediumType==eAllNode)
            {
                m_zoneindex[0]  = Array<OneD, int>(nq, 1);
            }

            m_NeuralCm    = Array<OneD, Array<OneD, NekDouble>>(1);
            m_NeuralCm[0] = ComputeConductivity(m_zoneindex[0]);
            break;
        }

        default:
            break;
    }

    // Stimulus
    m_stimulus = NeuralStimulus::LoadStimuli(m_session, m_fields[0]);

    // Derive AnisotropyStrength.
    m_AnisotropyStrength = m_Cn / m_Cm;
    SetUpBiAnisotropy(m_zoneindex[0], m_NeuralCm, m_AniStrength);

    MMFSystem::MMFInitObject(m_AniStrength);

    // Check moving frames and anisotropy
    // CheckOutZoneAni();

    switch (m_NeuralEPType)
    {
        case eNeuralHelmTest:
        case eNeuralEP2Dbi:
        {
            std::string phieMMFdirStr = "TangentX";
            m_session->LoadSolverInfo("phieMMFDir", phieMMFdirStr, "LOCAL");

            Array<OneD, Array<OneD, NekDouble>> phieAniStrength(m_expdim);
            for (int j = 0; j < m_expdim; ++j)
            {
                phieAniStrength[j] = Array<OneD, NekDouble>(nq, 1.0);
            }

            SpatialDomains::GeomMMF phieMMFdir = FindMMFdir(phieMMFdirStr);
            SetUpMovingFrames(phieMMFdir, phieAniStrength, m_unitmovingframes);

            NekDouble sigma_e_ratio = m_AnisotropyStrength / m_ratio_re_ri;

            m_phiemovingframes = Array<OneD, Array<OneD, NekDouble>>(m_expdim);
            for (int j = 0; j < m_expdim; ++j)
            {
                m_phiemovingframes[j] = Array<OneD, NekDouble>(m_spacedim * nq);
            }

            // m_phieMF = \sigma_i + \sigma_e

            for (int j = 0; j < m_expdim; ++j)
                {
                    for (int k = 0; k < m_spacedim; ++k)
                    {
                        for (int i = 0; i < nq; ++i)
                    {
                        m_phiemovingframes[j][k * nq + i] = m_movingframes[j][k*nq+i]; 
                                                               + sigma_e_ratio * m_unitmovingframes[j][k * nq + i];
                    }
                }
            }
            break;
        }

        default:
         break;
    }

    // switch (m_NeuralEPType)
    // {
    //     case eNeuralHelmTest:
    //     case eNeuralEP2Dbi:
    //     {
    //         // Expression of phie moving frames in phieMMFdirStr
    //         // std::string phieMMFdirStr = "LOCAL";
    //         // m_session->LoadSolverInfo("phieMMFDir", phieMMFdirStr, "LOCAL");
    //         // SpatialDomains::GeomMMF phieMMFdir = FindMMFdir(phieMMFdirStr);

    //         // Array<OneD, Array<OneD, NekDouble>> m_phieunitMF(m_mfdim);
    //         // for (int i=0; i<m_mfdim; ++i)
    //         // {
    //         //     m_phieunitMF[i] = Array<OneD, NekDouble>(m_spacedim * nq);
    //         // }

    //         // Array<OneD, Array<OneD, NekDouble>> unitAniStrength(m_mfdim);
    //         // for (int j = 0; j < m_mfdim; ++j)
    //         // {
    //         //     unitAniStrength[j] = Array<OneD, NekDouble>(nq, 1.0);
    //         // }

    //         // std::cout << std::endl;
    //         // std::cout << "Set up moving frames with unitAniStrength ==================== " << std::endl;
    //         // SetUpMovingFrames(phieMMFdir, unitAniStrength, m_phieunitMF);

    //         // // a \vec{e}^{LOC}_1 + b \vec{e}^{LOC}_2 = \vec{e}_1 + \vec{e}_2
    //         // // a = ( \vec{e}_1 \cdot \vec^{LOC}_1 ) + ( \vec{e}_2 \cdot \vec^{LOC}_1 )
    //         // // a = ( \vec{e}_1 \cdot \vec^{LOC}_1 ) + ( \vec{e}_2 \cdot \vec^{LOC}_1 )
    //         // Array<OneD, Array<OneD, NekDouble>> m_phieAniStrength(m_expdim);
    //         // for (int j = 0; j < m_expdim; ++j)
    //         // {
    //         //     m_phieAniStrength[j] = Array<OneD, NekDouble>(nq, 1.0);
    //         // }

    //         // // Expression of phie moving frames in MMFdir
    //         // NekDouble tmp;
    //         // for (int i = 0; i<nq; ++i)
    //         // {
    //         //     for (int j = 0; j < m_expdim; ++j)
    //         //     {
    //         //         m_phieAniStrength[j][i] = ( m_Cn/m_Cm ) / m_ratio_re_ri;
    //         //         for (int k=0; k<m_spacedim; ++k)
    //         //         {
    //         //             tmp = DotproductMF(i, m_movingframes[k], m_phieunitMF[j]);
    //         //             m_phieAniStrength[j][i] += tmp * tmp;
    //         //         }
    //         //     }
    //         // }

    //         // std::cout << "Max phieAnistrength_1  = "
    //         //             << Vmath::Vmax(nq, m_phieAniStrength[0], 1)
    //         //             << ", phieAnistrength_2 = "
    //         //             << Vmath::Vmax(nq, m_phieAniStrength[1], 1)
    //         //             << ", Min phieAnistrength 1 = "
    //         //             << Vmath::Vmin(nq, m_phieAniStrength[0], 1)
    //         //             << ", phieAnistrength 2 = "
    //         //             << Vmath::Vmin(nq, m_phieAniStrength[1], 1) << std::endl;

    //         // std::cout << std::endl;
    //         // std::cout << "Set up moving frames with m_phieAniStrength ==================== " << std::endl;
    //         // SetUpMovingFrames(phieMMFdir, m_phieAniStrength, m_phiemovingframes);

    //         std::string phieMMFdirStr = "TangentX";
    //         m_session->LoadSolverInfo("phieMMFDir", phieMMFdirStr, "LOCAL");

    //         SpatialDomains::GeomMMF phieMMFdir = FindMMFdir(phieMMFdirStr);
    //         SetUpMovingFrames(phieMMFdir, m_phieAniStrength, m_unitmovingframes);

    //         NekDouble sigma_e_ratio = m_AnisotropyStrength / m_ratio_re_ri;

    //         m_phiemovingframes = Array<OneD, Array<OneD, NekDouble>> (m_spacedim);
    //         for (int j = 0; j < m_spacedim; ++j)
    //         {
    //             m_phiemovingframes[j] = Array<OneD, NekDouble>(m_spacedim * nq);
    //         }

    //         // m_phieMF = \sigma_i + \sigma_e
    //         for (int i = 0; i < nq; ++i)
    //         {
    //             for (int j = 0; j < m_shapedim; ++j)
    //             {
    //                 for (int k = 0; k < m_spacedim; ++k)
    //                 {
    //                     m_phiemovingframes[j][k * nq + i] = m_movingframes[j][k*nq+i] 
    //                                                           + sigma_e_ratio * m_unitmovingframes[j][k * nq + i];
    //                 }
    //             }
    //         }
    //         break;
    //     }

    //     default:
    //      break;
    // }

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
            {
                // ComputeVarCoeff1D(m_movingframes, m_varcoeff);
                // m_ode.DefineImplicitSolve(
                //     &MMFNeuralEP::DoImplicitSolveNeuralEP1D, this);
                m_ode.DefineImplicitSolve(&MMFNeuralEP::DoNullSolve, this);
                m_ode.DefineOdeRhs(&MMFNeuralEP::DoOdeRhsNeuralEPPT, this);
                break;
            }

            case eNeuralEP1D:
            {
                // ComputeVarCoeff1D(m_movingframes, m_varcoeff);
                ComputeVarCoeff2D(m_movingframes, m_varcoeff);
                m_ode.DefineImplicitSolve(
                    &MMFNeuralEP::DoImplicitSolveNeuralEP1D, this);
                m_ode.DefineOdeRhs(&MMFNeuralEP::DoOdeRhsNeuralEP1D, this);
                break;
            }

            case eNeuralHelmTest:
            case eNeuralEP2Dmono:
            {
                ComputeVarCoeff2D(m_movingframes, m_varcoeff);
                m_ode.DefineImplicitSolve(
                    &MMFNeuralEP::DoImplicitSolveNeuralEP2Dmono, this);
                m_ode.DefineOdeRhs(&MMFNeuralEP::DoOdeRhsNeuralEP2Dmono, this);
                break;
            }

            case eNeuralEP2Dbi:
            {
                std::cout << "Generating moving frames for the domain with anisotropic fiber" << std::endl;
                ComputeVarCoeff2D(m_movingframes, m_varcoeff);
                std::cout << std::endl;

                std::cout << "Generating moving frames for the distribution of the external potential with point excitation at the Raniver node" << std::endl;
                ComputeVarCoeff2D(m_phiemovingframes, m_phievarcoeff);
                std::cout << std::endl;

                m_ode.DefineImplicitSolve(
                    &MMFNeuralEP::DoImplicitSolveNeuralEP2Dbi, this); 

                m_ode.DefineOdeRhs(&MMFNeuralEP::DoOdeRhsNeuralEP2Dbi, this);
                break;
            }

            default:
                break;
        }
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

        m_zoneindex = Array<OneD, Array<OneD, int>>(1);
        for (int i = 0; i < 1; i++)
        {
            m_zoneindex[i] = Array<OneD, int>(nq, 1);
        }

        std::cout << "Constructing phieunitMF ================================================"
                << std::endl;
        Array<OneD, Array<OneD, NekDouble>> phieAniStrength(m_expdim);
        for (int j = 0; j < m_expdim; ++j)
        {
            phieAniStrength[j] = Array<OneD, NekDouble>(nq, 1.0);
        }

        SetUpMovingFrames(FindMMFdir("LOCAL"), phieAniStrength, m_unitmovingframes);

        std::cout << "Constructing phiemovingframes ================================================"
                << std::endl;

        m_phiemovingframes =
            Array<OneD, Array<OneD, NekDouble>>(m_mfdim);
        NekDouble Helmfactor =
            sqrt((1.0 + m_ratio_re_ri) / m_ratio_re_ri);

        // \sigma_i is zero outside the fiber. 
        for (int j = 0; j < m_expdim; ++j)
        {
            Vmath::Smul(nq, Helmfactor, &phieAniStrength[j][0], 1,
                        &phieAniStrength[j][0], 1);
        }

        SetUpMovingFrames(FindMMFdir("LOCAL"), phieAniStrength,
                            m_phiemovingframes);
        CheckMovingFrames(m_phiemovingframes);

        std::cout << "Constructing m_phievarcoeff ================================================"
                << std::endl;
        ComputeVarCoeff2D(m_phiemovingframes, m_phievarcoeff);

        std::cout << "SolveHelmholtzDiffusion ================================================"
                << std::endl;
        outputarray = Computephie(inputarray[0]);

        Vmath::Smul(nq, m_Rf * m_Cn, outputarray, 1, outputarray, 1);

        Vmath::Vsub(nq, ExactSoln[1], 1, outputarray, 1, outputarray, 1);
        std::cout << "Computephie Error = " << RootMeanSquare(outputarray) << std::endl;

        wait_on_enter();
    }

}

/**
 *
 */
MMFNeuralEP::~MMFNeuralEP()
{
}

NekDouble MMFNeuralEP::DotproductMF(
    const int i, 
    const Array<OneD, const NekDouble> &MF, 
    const Array<OneD, const NekDouble> &MFloc)
{
    int nq   = GetTotPoints();

    NekDouble outarray = 0.0;
    for (int k=0; k<m_spacedim; ++k)
    {
        outarray += MF[k * nq + i] * MFloc[k * nq + i];
    }

    return outarray;
}

void MMFNeuralEP::CheckOutZoneAni()
{
    int nq   = GetTotPoints();

    Array<OneD, NekDouble> x0(nq);
    Array<OneD, NekDouble> x1(nq);
    Array<OneD, NekDouble> x2(nq);

    m_fields[0]->GetCoords(x0, x1, x2);

    for (int i=0; i<nq; ++i)
    {
        std::cout << "i = " << i << ", (x,y) = ( " << x0[i] << " , " << x1[i] << " ), zoneindex = " 
        << m_zoneindex[0][i] << ", Anistrength = ( " << m_AniStrength[0][i] << " , " << m_AniStrength[1][i]
        << " ), e1y = " << m_movingframes[0][i+nq] << ", e2y = " << m_movingframes[1][i] << std::endl;
    }

}

Array<OneD, int> MMFNeuralEP::IndexNodeZone2D(const FiberType fiberType)
{
    int nq   = GetTotPoints();
    
    Array<OneD, int> outarray(nq);

    switch(fiberType)
    {
        case eSingleLinear:
        {
            outarray = IndexNodeSingleFiber(m_fiberlen, m_nodelen, m_myelinlen, m_totNode);
        }
        break;

        case eDoubleLinear:
        {
            // outarray = IndexNodeDoubleFiber(m_fiberlen, m_nodelen, m_myelinlen, m_Node);
        }
        break;

        default:
        break;
    }

    return outarray;
}

Array<OneD, int> MMFNeuralEP::IndexNodeSingleFiber(
    const NekDouble fiberlen, 
    const NekDouble nodelen, 
    const NekDouble myelinlen,
    const int totNode)
{
    int nq   = GetTotPoints();

    Array<OneD, NekDouble> x0(nq);
    Array<OneD, NekDouble> x1(nq);
    Array<OneD, NekDouble> x2(nq);

    m_fields[0]->GetCoords(x0, x1, x2);

    Array<OneD, NekDouble> xcell(nq);
    Array<OneD, NekDouble> ycell(nq);
    Array<OneD, NekDouble> zcell(nq);

    int Nelem = nq/m_npts;

    Array<OneD, NekDouble> xcellavg(Nelem,0.0);
    Array<OneD, NekDouble> ycellavg(Nelem,0.0);
    Array<OneD, NekDouble> zcellavg(Nelem,0.0);

    int index;
    for (int i=0; i<nq; ++i)
    {
        index = i/m_npts;
        xcellavg[index] = xcellavg[index] + x0[i];
        ycellavg[index] = ycellavg[index] + x1[i];
        zcellavg[index] = zcellavg[index] + x2[i];
    }

    Vmath::Smul(Nelem, 1.0/m_npts, xcellavg, 1, xcellavg, 1);
    Vmath::Smul(Nelem, 1.0/m_npts, ycellavg, 1, ycellavg, 1);
    Vmath::Smul(Nelem, 1.0/m_npts, zcellavg, 1, zcellavg, 1);

    for (int i=0; i<nq; ++i)
    {
        index = i/m_npts;

        xcell[i] = xcellavg[index];
        ycell[i] = ycellavg[index];
        zcell[i] = zcellavg[index];
    }

    Array<OneD, int> outarray(nq, -1);

    NekDouble nodestart, nodeend;
    NekDouble fiberleft = -0.5*fiberlen;
    NekDouble fiberright = 0.5*fiberlen;

    for (int i=0; i<nq; ++i)
    {
        if((xcell[i]>fiberleft) && (xcell[i]<fiberright))
        {
            for (int k=0; k<2; ++k)
            {
                if( (ycell[i]>=k*nodelen) && (ycell[i]<=(k+1)*nodelen) )
                {
                    outarray[i] = k;
                }
            }

            for (int k=0; k<totNode; ++k)
            {
                nodestart = 2.0*nodelen + myelinlen + k * (myelinlen + nodelen);
                nodeend = 2.0*nodelen + (k+1) * (myelinlen + nodelen);
                if( (ycell[i]>=nodestart) && (ycell[i]<=nodeend) )
                {
                    outarray[i] = k+2;
                }
            }
        }

        else
        {
            outarray[i] = -2;
        }
    }

    return outarray;
}


// Array<OneD, int> MMFNeuralEP::IndexNodeZone2D(
//     const NekDouble fiberlen, 
//     const NekDouble nodelen, 
//     const NekDouble myelinlen,
//     const int Nnode)
// {
//     int nq   = GetTotPoints();

//     Array<OneD, NekDouble> x0(nq);
//     Array<OneD, NekDouble> x1(nq);
//     Array<OneD, NekDouble> x2(nq);

//     m_fields[0]->GetCoords(x0, x1, x2);

//     Array<OneD, NekDouble> xcell(nq);
//     Array<OneD, NekDouble> ycell(nq);
//     Array<OneD, NekDouble> zcell(nq);

//     int Nelem = nq/m_npts;

//     Array<OneD, NekDouble> xcellavg(Nelem,0.0);
//     Array<OneD, NekDouble> ycellavg(Nelem,0.0);
//     Array<OneD, NekDouble> zcellavg(Nelem,0.0);

//     int index;
//     for (int i=0; i<nq; ++i)
//     {
//         index = i/m_npts;
//         xcellavg[index] = xcellavg[index] + x0[i];
//         ycellavg[index] = ycellavg[index] + x1[i];
//         zcellavg[index] = zcellavg[index] + x2[i];
//     }

//     Vmath::Smul(Nelem, 1.0/m_npts, xcellavg, 1, xcellavg, 1);
//     Vmath::Smul(Nelem, 1.0/m_npts, ycellavg, 1, ycellavg, 1);
//     Vmath::Smul(Nelem, 1.0/m_npts, zcellavg, 1, zcellavg, 1);

//     for (int i=0; i<nq; ++i)
//     {
//         index = i/m_npts;

//         xcell[i] = xcellavg[index];
//         ycell[i] = ycellavg[index];
//         zcell[i] = zcellavg[index];
//     }

//     Array<OneD, int> outarray(nq, -1);

//     NekDouble nodestart, nodeend;
//     NekDouble fiberleft = -0.5*fiberlen;
//     NekDouble fiberright = 0.5*fiberlen;

//     for (int i=0; i<nq; ++i)
//     {
//         if((xcell[i]>fiberleft) && (xcell[i]<fiberright))
//         {
//             for (int k=0; k<2; ++k)
//             {
//                 if( (ycell[i]>=k*nodelen) && (ycell[i]<=(k+1)*nodelen) )
//                 {
//                     outarray[i] = k;
//                 }
//             }

//             for (int k=0; k<Nnode; ++k)
//             {
//                 nodestart = 2.0*nodelen + myelinlen + k * (myelinlen + nodelen);
//                 nodeend = 2.0*nodelen + (k+1) * (myelinlen + nodelen);
//                 if( (ycell[i]>=nodestart) && (ycell[i]<=nodeend) )
//                 {
//                     outarray[i] = k+2;
//                 }
//             }
//         }

//         else
//         {
//             outarray[i] = -2;
//         }
//     }

//     return outarray;
// }

// void MMFNeuralEP::ComputephieMF(
//     const NekDouble ratio_re_ri,
//     const Array<OneD, const int> &zoneindex,
//     Array<OneD, Array<OneD, NekDouble>> &unitfmovingframes, 
//     Array<OneD, Array<OneD, NekDouble>> &phiemovingframes)

// {
//     int nq   = GetTotPoints();

//     std::string phieMMFdirStr = "TangentX";
//     m_session->LoadSolverInfo("phieMMFDir", phieMMFdirStr, "LOCAL");

//     Array<OneD, Array<OneD, NekDouble>> phieAniStrength(m_expdim);
//     for (int j = 0; j < m_expdim; ++j)
//     {
//         phieAniStrength[j] = Array<OneD, NekDouble>(nq, 1.0);
//     }

//     // Create UnitMovingFrames
//     SpatialDomains::GeomMMF phieMMFdir = FindMMFdir(phieMMFdirStr);
//     SetUpMovingFrames(phieMMFdir, phieAniStrength, unitfmovingframes);

//     NekDouble sigma_e = m_AnisotropyStrength / ratio_re_ri;
//     for (int i=0; i<nq; ++i)
//         {
//             for (int j = 0; j < m_expdim; ++j)
//             {
//                 // Node zone
//                 if(zoneindex[i]>=0)
//                 {
//                     phieAniStrength[j][i] = 1.0 + sigma_e;
//                 }

//                 // Myeline zone
//                 else if (zoneindex[i]==-1)
//                 {
//                     phieAniStrength[j][i] = (1.0 + ratio_re_ri) / ratio_re_ri;
//                 }

//                 // Extracellular zone
//                 else
//                 {
//                     phieAniStrength[j][i] = sigma_e;
//                 }
//             }
//     }

//     // Create Phiemovingframes
//     phiemovingframes = Array<OneD, Array<OneD, NekDouble>>(m_mfdim);

//     SetUpMovingFrames(phieMMFdir, phieAniStrength, phiemovingframes);
//     CheckMovingFrames(phiemovingframes);
// }

void MMFNeuralEP::SetUpBiAnisotropy(
            const Array<OneD, const int> &zoneindex,
            const Array<OneD, const Array<OneD, NekDouble>> NeuralCm,
            Array<OneD, Array<OneD, NekDouble>> &AniStrength)
{
    int nq   = GetTotPoints();
    int nvar = m_fields.size();

    for (int j = 0; j < m_expdim; ++j)
    {
        AniStrength[j] = Array<OneD, NekDouble>(nq, 1.0);
    }

    if (m_MediumType == eAnisotropy)
    {
        for (int j = 0; j < m_expdim; ++j)
        {
            Vmath::Smul(nq, m_Cn, &NeuralCm[0][0], 1, &AniStrength[j][0], 1);
            Vmath::Vsqrt(nq, &AniStrength[j][0], 1, &AniStrength[j][0], 1);
        }
    }

    // Let the lenght of moving frames outside the fiber to be zero.
    int cnt = 0;
    for (int i=0; i<nq; ++i)
    {
        if(zoneindex[i] == -2)
        {
            AniStrength[0][i] = 0.0;
            AniStrength[1][i] = 0.0; 
            cnt++;
        }
    }

    std::cout << " =======================================================================" << std::endl;
    std::cout << " Moving frames " << cnt << " / " << nq << " ( " << 100.0*cnt/nq << " % ) are removed" << std::endl;
    std::cout << " =======================================================================" << std::endl;


    if(nvar==2)
    {
        std::cout << "Max Anistrength_1  = "
                    << Vmath::Vmax(nq, AniStrength[0], 1)
                    << ", Anistrength_2 = "
                    << Vmath::Vmax(nq, AniStrength[1], 1)
                    << ", Min Anistrength 1 = "
                    << Vmath::Vmin(nq, AniStrength[0], 1)
                    << ", Anistrength 2 = "
                    << Vmath::Vmin(nq, AniStrength[1], 1) << std::endl;
    }

    else{
        std::cout << "Max Anistrength_1  = "
            << Vmath::Vmax(nq, AniStrength[0], 1)
            << ", Min Anistrength 1 = "
            << Vmath::Vmin(nq, AniStrength[0], 1) << std::endl;
    }
}

void MMFNeuralEP::SetUpDomainZone(
        const Array<OneD, const int> &zoneindex,
        Array<OneD, NekDouble> &excitezone,
        Array<OneD, NekDouble> &nodezone,
        Array<OneD, NekDouble> &intrazone,
        Array<OneD, NekDouble> &extrazone)
{
    int nq   = GetTotPoints();

    int extcnt=0, nodecnt=0, intracnt=0, extracnt=0;

    excitezone = Array<OneD, NekDouble>(nq, 0.0);
    nodezone  = Array<OneD, NekDouble>(nq, 0.0);
    intrazone = Array<OneD, NekDouble>(nq, 0.0);
    extrazone = Array<OneD, NekDouble>(nq, 0.0);
    for (int i=0; i<nq; ++i)
    {
        if( (zoneindex[i]==0) || (zoneindex[i]==1) )
        {
            excitezone[i] = 1.0 / m_Cn;
            extcnt++;
        }

        if(zoneindex[i] >= 0)
        {
            nodezone[i] = 1.0;
            nodecnt++;
        }

        if(zoneindex[i] != -2)
        {
            intrazone[i] = 1.0 ;
            intracnt++;
        }

        if(zoneindex[i] != -1)
        {
            extrazone[i] = 1.0;
            extracnt++;
        }
    }

    int extcntelem = extcnt/m_npts/m_elemperNode;
    int nodecntelem = nodecnt/m_npts/m_elemperNode;
    int myelcntelem = (intracnt-nodecnt)/m_npts/m_elemperMyel;
    int intraelem = intracnt/m_npts;
    int extraelem = extracnt/m_npts;
    int totelem = intraelem + extraelem - nodecntelem*m_elemperNode;
    
    std::cout << "Excite zone = " << 100.0*extcnt/nq << ", Node zone = " << 100.0*nodecntelem/nq <<  ", intra zone = " << 100.0*intracnt/nq << 
    " %, extra zone = " << 100.0*extracnt/nq << " % " << std::endl;

    std::cout << "npts = " << m_npts << ", Excite node elem. = " << extcntelem << ", Node elem. = " << nodecntelem 
    << ", Myelin elem. = " << myelcntelem << ", intra elem. = " << intraelem << ", external elem. = " << extraelem 
    << ", total elem = " << totelem << std::endl;
}

Array<OneD, NekDouble> MMFNeuralEP::ComputeConductivity(
                 const Array<OneD, const int> &zoneindex)
{
    int nq   = GetTotPoints();

    Array<OneD, NekDouble> outarray(nq);

    int cntm = 0, cntn = 0, cnte = 0;
    for (int i = 0; i < nq; ++i)
    {
        // Ranvier node zone
        if (zoneindex[i] >= 0)

        {
            outarray[i] = 1.0 / m_Cn;
            cntm++;
        }

        // Myelin zone
        else if (zoneindex[i] == -1)
        {
            outarray[i] = 1.0 / m_Cm;
            cntn++;
        }

        // Extracellular space: \sigma_i = m_ratio_re_ri * \sigma_e
        else if (zoneindex[i] == -2)
        {
            outarray[i] = 1.0 / m_Cm / m_ratio_re_ri;
            cnte++;
        }
    }

    std::cout << "v_InitObject: Node = " << cntn
    << ", Myelin = " << cntm << ", extracell = " << cnte
    << std::endl;

    return outarray;
}

void MMFNeuralEP::CheckNodeZoneMF(
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    const Array<OneD, const Array<OneD, int>> &NodeZone,
    const Array<OneD, const NekDouble> &inarray)
{
    boost::ignore_unused(NodeZone);

    int nq = GetTotPoints();

    Array<OneD, NekDouble> x0(nq);
    Array<OneD, NekDouble> x1(nq);
    Array<OneD, NekDouble> x2(nq);

    m_fields[0]->GetCoords(x0, x1, x2);

    int i, j, index, npts;;
    NekDouble xp, yp, e1mag, e2mag, inarrayavg;
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

        // dx = x0[m_fields[0]->GetPhys_Offset(i)] - x0[m_fields[0]->GetPhys_Offset(i)+npts-1];
        // dy = x1[m_fields[0]->GetPhys_Offset(i)] - x1[m_fields[0]->GetPhys_Offset(i)+npts-1];
        // dist = sqrt(dx*dx+dy*dy);

        e1mag = sqrt(e1mag / npts);
        e2mag = sqrt(e2mag / npts);
        xp = (xp / npts);
        yp = (yp / npts);
        inarrayavg = inarrayavg/npts;

        // std::cout << "Elemid = " << i << ", Nodeid = " << NodeZone[0][index]
        //         << ", x = " << xp << ", y = " << yp << ", dist = " << dist
        //         << ", e1mag = " << e1mag << ", e2mag = " << e1mag
        //         << ", inarray = " << inarrayavg << std::endl;
    }
}

    Array<OneD, int> MMFNeuralEP::GetInternalBoundaryPoints()
    {
        int nq    = GetNpoints();
        int nTracePts  = GetTraceTotPoints();

        Array<OneD, int> outarray(nq,0);

        Array<OneD, NekDouble> x0(nq);
        Array<OneD, NekDouble> x1(nq);
        Array<OneD, NekDouble> x2(nq);

        m_fields[0]->GetCoords(x0, x1, x2);

        Array<OneD, NekDouble> x0Fwd(nTracePts);
        Array<OneD, NekDouble> x1Fwd(nTracePts);
        Array<OneD, NekDouble> x2Fwd(nTracePts);

        m_fields[0]->ExtractTracePhys(x0, x0Fwd);
        m_fields[0]->ExtractTracePhys(x1, x1Fwd);
        m_fields[0]->ExtractTracePhys(x2, x2Fwd);

        NekDouble xp, yp, distx, disty, dist;
        NekDouble Tol = 0.000001;

        int cnt=0;
        for (int n = 0; n < m_fields[0]->GetBndConditions().size(); ++n)
        {
            if (boost::iequals(m_fields[0]->GetBndConditions()[n]->GetUserDefined(), "Membrane"))
            {
                int id2, index, npts;

                const Array<OneD, const int> &traceBndMap = m_fields[0]->GetTraceBndMap();

                for (int e = 0; e < m_fields[0]->GetBndCondExpansions()[n]->GetExpSize(); ++e)
                {
                    npts = m_fields[0]
                                    ->GetBndCondExpansions()[n]
                                    ->GetExp(e)
                                    ->GetTotPoints();
                    // id1 = m_fields[0]->GetBndCondExpansions()[n]->GetPhys_Offset(e);
                    id2 = m_fields[0]->GetTrace()->GetPhys_Offset(traceBndMap[cnt + e]);

                    for (int i=0;i<npts;++i)
                    {
                        index = id2+i;
                        xp = x0Fwd[index];
                        yp = x1Fwd[index];
                        for (int j=0; j<nq; ++j)
                        {
                            distx = xp - x0[j];
                            disty = yp - x1[j];
                            dist = sqrt(distx*distx + disty*disty);

                            if (dist<Tol)
                            {
                                outarray[j] = 1;
                            }
                        }
                    }
                }
            }
            cnt += m_fields[0]->GetBndCondExpansions()[n]->GetExpSize();
        }

        return outarray;
    }

// Constrcuct Cm vector: 1.0/Cn if node. 1.0/Cm if myeline.
Array<OneD, int> MMFNeuralEP::IndexNodeZone1D(
    const MultiRegions::ExpListSharedPtr &field, const int Nnode, 
    const int NumelemNode, const int NumelemMyel)
{
    int nq = field->GetNpoints();

    Array<OneD, NekDouble> x0(nq);
    Array<OneD, NekDouble> x1(nq);
    Array<OneD, NekDouble> x2(nq);

    field->GetCoords(x0, x1, x2);

    Array<OneD, int> outarray(nq, -1);
    int cntn=0, cntm=0;

    // Node 0
    int indexNode;
    for (int i = 0; i < 2 * NumelemNode; ++i)
    {
        for (int j = 0; j < m_fields[0]->GetTotPoints(i); ++j)
        {
            indexNode = m_fields[0]->GetPhys_Offset(i) + j ;

            // First and last element is all node for easier excitation
            outarray[indexNode] = i / NumelemNode;
            cntn++;
        }
    }

    // Myelin index
    int Myelid;
    int indexMyelin;
    for (int k = 0; k < Nnode; ++k)
    {
        for (int i = 0; i < NumelemMyel; ++i)
        {
            Myelid = (NumelemNode + NumelemMyel)*k + 2 * NumelemNode + i;
            for (int j = 0; j < m_fields[0]->GetTotPoints(Myelid); ++j)
            {
                indexMyelin = m_fields[0]->GetPhys_Offset(Myelid) + j ;

                // First and last element is all node for easier excitation
                outarray[indexMyelin] = -1;
                cntm++;
            }
        }
    }

    // Node index
    int Nodeid;
    for (int k = 0; k < Nnode; ++k)
    {
        for (int i = 0; i < NumelemNode; ++i)
        {
            Nodeid = (NumelemNode + NumelemMyel)*k + 2 * NumelemNode + NumelemMyel + i;
            for (int j = 0; j < m_fields[0]->GetTotPoints(Nodeid); ++j)
            {
                indexNode = m_fields[0]->GetPhys_Offset(Nodeid) + j ;

                // First and last element is all node for easier excitation
                outarray[indexNode] = k + 2;
                cntn++;
            }
        }
    }

    std::cout << "cntn = " << cntn << ", cntm = " << cntm << ", cnte = " << (nq-cntn-cntm) << std::endl;

    // for (int i=0; i<nq; ++i)
    // {
    //     std::cout << "i = " << i << ", y = " << x1[i] << ", zoneindex = " << outarray[i] << std::endl;
    // }

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
        case eMMFZero:
        case eMMFFirst:
        case eTimeMap:
        {
            DoSolveMMFZero();
            break;
        }

        case ePointWise:
        {
            DoSolvePoint();
            break;
        }

        default:
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
            std::cout << "Steps: " << std::setw(8) << std::left << step + 1
                      << " "
                      << "Time: " << std::setw(12) << std::left << m_time
                      << std::endl;

            fulltext.append("\n");
            fulltext.append("Time: " + std::to_string(m_time));
            fulltext.append("\n");

            fulltext.append("CPU Time: " + std::to_string(cpuTime / 60.0) + " min.");
            fulltext.append("\n");

            cpuTime = 0.0;
        }

        // Write out checkpoint files
        if ((m_checksteps && step && !((step + 1) % m_checksteps)) ||
            doCheckTime)
        {
            // phim should be only defined in the intracellular space
            Array<OneD, NekDouble> phi_m(nq);
            Vmath::Vmul(nq, m_intrazone, 1, fields[0], 1, phi_m, 1);

            NekDouble phimMax = Vmath::Vmax(nq, phi_m, 1);
            NekDouble phimMin = Vmath::Vmin(nq, phi_m, 1);

            int phimMaxid = Vmath::Imax(nq, phi_m, 1);
            int phimMinid = Vmath::Imin(nq, phi_m, 1);

            fulltext.append("phi_m, max: " + std::to_string(phimMax) + " at y = " + std::to_string(x1[phimMaxid]) );
            fulltext.append("./ phi_m, min: " + std::to_string(phimMin) + " at y = " + std::to_string(x1[phimMinid]) );

            fulltext.append("\n");

            if( (nvariables==2) && (m_ExtCurrentType == eEphaptic) )
            {
                // phi_e is defined at the node and extracellular space
                Array<OneD, NekDouble> phi_e = Computephie(fields[0]);
                
                NekDouble phieMax = Vmath::Vmax(nq, phi_e, 1);
                NekDouble phieMin = Vmath::Vmin(nq, phi_e, 1);

                int phieMaxid = Vmath::Imax(nq, phi_e, 1);
                int phieMinid = Vmath::Imin(nq, phi_e, 1);

                fulltext.append("phi_e, max: " + std::to_string(phieMax) + " at y = " + std::to_string(x1[phieMaxid]) );
                fulltext.append(", phi_e, min: " + std::to_string(phieMin) + " at y = " + std::to_string(x1[phieMinid]) );

                fulltext.append("\n");

               //  PlotFields(phi_m, phi_e, nchk);
            }

            // if(m_NeuralEPType==eNeuralEP2Dbi)
            // {
            //     DisplayNode2D(fulltext, fields);
            // }

            std::cout << fulltext << "\n" << std::endl;

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
} 
// namespace Nektar

// void MMFNeuralEP::PrintRegionalAvgMax(const Array<OneD, const NekDouble> &field0)
// {
//     int index;
//     int Nodeid, Myelid, Extid;
//     int nq               = GetTotPoints();

//     NekDouble NodeMaxm, MyelineMaxm, ExtMaxm;
//     NekDouble NodeMaxe, MyelineMaxe, ExtMaxe;

//     NekDouble ue, um, elemavgm, elemavge;

//     NekDouble yavgindex, yavg;

//     Array<OneD, NekDouble> x0(nq);
//     Array<OneD, NekDouble> x1(nq);
//     Array<OneD, NekDouble> x2(nq);

//     m_fields[0]->GetCoords(x0, x1, x2);
    
//     std::cout << " ========================================================================================== " << std::endl;

//     // Max um and ue at Node
//     NodeMaxm = 0.0;
//     NodeMaxe = 0.0;

//     Nodeid = 0;
//     yavg = 0.0;
//     for (int i = 0; i < m_ElemNodeEnd; ++i)
//     {
//         elemavgm = 0.0;
//         elemavge = 0.0;
//         yavg = 0.0;
//         for (int j = 0; j < m_fields[0]->GetTotPoints(i); ++j)
//         {
//             index = m_fields[0]->GetPhys_Offset(i) + j;
//             um = field0[index];
//             ue = (m_fields[1]->GetPhys())[index];

//             elemavgm += um;
//             elemavge += ue;
//             yavg += x1[index];
//         }
//         elemavgm = elemavgm/m_fields[0]->GetTotPoints(i);
//         elemavge = elemavge/m_fields[0]->GetTotPoints(i);
//         yavg = yavg/m_fields[0]->GetTotPoints(i);

//         if(elemavgm>NodeMaxm)
//         {
//             NodeMaxm = elemavgm;
//             yavgindex = yavg;
//             Nodeid = i;
//         }

//         if(elemavge>NodeMaxe)
//         {
//             NodeMaxe = elemavge;
//         }
//     }

//     std::cout << "Node id = " << Nodeid << ", : um_max = " << NodeMaxm << " at y = " << yavgindex << ", ue_max = " << NodeMaxe << std::endl;

//     // Max um and ue at Myelin
//     MyelineMaxm = 0.0;
//     MyelineMaxe = 0.0;
    
//     Myelid = m_ElemNodeEnd;
//     yavg = 0.0;
//     for (int i = m_ElemNodeEnd; i < m_ElemMyelenEnd ; ++i)
//     {
//         elemavgm = 0.0;
//         elemavge = 0.0;
//         yavg = 0.0;
//         for (int j = 0; j < m_fields[0]->GetTotPoints(i); ++j)
//         {
//             index = m_fields[0]->GetPhys_Offset(i) + j;
//             um = field0[index];
//             ue = (m_fields[1]->GetPhys())[index];

//             elemavgm += um;
//             elemavge += ue;
//             yavg += x1[index];
//         }
//         elemavgm = elemavgm/m_fields[0]->GetTotPoints(i);
//         elemavge = elemavge/m_fields[0]->GetTotPoints(i);
//         yavg = yavg/m_fields[0]->GetTotPoints(i);

//         if(elemavgm>MyelineMaxm)
//         {
//             MyelineMaxm = elemavgm;
//             yavgindex = yavg;
//             Myelid = i;
//         }

//         if(elemavge>MyelineMaxe)
//         {
//             MyelineMaxe = elemavge;
//         }
//     }

//     std::cout << "Myelid id = " << Myelid << ", : um_max = " << MyelineMaxm << " at y = " << yavgindex << ", ue_max = " << MyelineMaxe << std::endl;

//     // Max um and ue at Exterial space
//     ExtMaxm = 0.0;
//     ExtMaxe = 0.0;

//     Extid = m_ElemMyelenEnd;
//     yavg = 0.0;
//     for (int i = m_ElemMyelenEnd; i < m_fields[0]->GetExpSize() ; ++i)
//     {
//         elemavgm = 0.0;
//         elemavge = 0.0;
//         yavg = 0.0;
//         for (int j = 0; j < m_fields[0]->GetTotPoints(i); ++j)
//         {
//             index = m_fields[0]->GetPhys_Offset(i) + j;
//             um = field0[index];
//             ue = (m_fields[1]->GetPhys())[index];

//             elemavgm += um;
//             elemavge += ue;
//             yavg += x1[index];
//         }
//         elemavgm = elemavgm/m_fields[0]->GetTotPoints(i);
//         elemavge = elemavge/m_fields[0]->GetTotPoints(i);
//         yavg = yavg/m_fields[0]->GetTotPoints(i);

//         if(elemavgm>ExtMaxm)
//         {
//             ExtMaxm = elemavgm;
//             yavgindex = yavg;
//             Extid = i;
//         }

//         if(elemavge>ExtMaxe)
//         {
//             ExtMaxe = elemavge;
//         }
//     }

//     std::cout << "Extid id = " << Extid << ", : um_max = " << ExtMaxm << " at y = " << yavgindex << ", ue_max = " << ExtMaxe << std::endl;

//     std::cout << " ========================================================================================== " << std::endl;
// }

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


void MMFNeuralEP::DisplayNode2D(std::string &fulltext, const Array<OneD, const Array<OneD, NekDouble>> &fields)
{
    int nq               = GetTotPoints();

    Array<OneD, NekDouble> phi_m(nq);
    Vmath::Vmul(nq, m_intrazone, 1, fields[0], 1, phi_m, 1);
    
    Array<OneD, NekDouble> phi_e = Computephie(fields[0]);

    Array<OneD, NekDouble> phimavg(m_totNode+2,0.0);
    Array<OneD, NekDouble> phieavg(m_totNode+2,0.0);

    Array<OneD, int> totnodenumpts(m_totNode+2,0);

    int nodeid;
    for (int i=0; i<nq; ++i)
    {
        if(m_zoneindex[0][i]>=0)
        {
            nodeid = m_zoneindex[0][i];
            phimavg[nodeid] += phi_m[i];
            phieavg[nodeid] += phi_e[i];

            totnodenumpts[nodeid] = totnodenumpts[nodeid] + 1;
        }
    }

    for (int j=0; j<m_totNode+2; ++j)
    {
        if (totnodenumpts[j]>0)
        {
            phimavg[j] = phimavg[j] / totnodenumpts[j];
            phieavg[j] = phieavg[j] / totnodenumpts[j];
        }
    }

    fulltext.append(" \n");
    fulltext.append("(Nodeid,phim,phie): ");
    for (int i = 0; i < m_totNode+1; ++i)
    {
        fulltext.append( "( " + std::to_string(i) + " , " + std::to_string(phimavg[i]) + " , " + std::to_string(phieavg[i]) + " ) ");
    }
    
    fulltext.append(" \n");
}



void MMFNeuralEP::savezoneindex(const Array<OneD, const int> &zoneindex)
{
    int nvar = 1;
    int ncoeffs = m_fields[0]->GetNcoeffs();
    int nq      = m_fields[0]->GetTotPoints();

    std::string outname;
    outname = m_sessionName + "_zoneindex.chk";

    std::vector<Array<OneD, NekDouble>> fieldcoeffs(nvar);
    for (int i = 0; i < nvar; ++i)
    {
        fieldcoeffs[i] = Array<OneD, NekDouble>(ncoeffs);
    }

    std::vector<std::string> variables(nvar);
    variables[0] = "zoneindex";

    Array<OneD, NekDouble> tmp(nq);

    for (int i=0; i<nq; ++i)
    {
        tmp[i] = 1.0 * zoneindex[i];
    }

    m_fields[0]->FwdTrans(tmp, fieldcoeffs[0]);

    WriteFld(outname, m_fields[0], fieldcoeffs, variables);
}


void MMFNeuralEP::PlotFields(const Array<OneD, const NekDouble> &phi_m,
                                   const Array<OneD, const NekDouble> &phi_e,
                                   const int nstep)
{
    int nvar = 4;
    int ncoeffs = m_fields[0]->GetNcoeffs();
    int nq      = m_fields[0]->GetTotPoints();

    std::string outname;
    outname = m_sessionName + "_field_" +
              boost::lexical_cast<std::string>(nstep) + ".chk";

    std::vector<Array<OneD, NekDouble>> fieldcoeffs(nvar);
    for (int i = 0; i < nvar; ++i)
    {
        fieldcoeffs[i] = Array<OneD, NekDouble>(ncoeffs);
    }

    std::vector<std::string> variables(nvar);
    variables[0] = "phi_m";
    variables[1] = "phi_e";
    variables[2] = "phi_i";
    variables[3] = "phieforcing";

    m_fields[0]->FwdTrans(phi_m, fieldcoeffs[0]);
    m_fields[0]->FwdTrans(phi_e, fieldcoeffs[1]);

    // phi_m = phi_i - phi_e
    Array<OneD, NekDouble> phi_i(nq);
    Vmath::Vadd(nq, phi_m, 1, phi_e, 1, phi_i, 1);
    Vmath::Vmul(nq, m_intrazone, 1, phi_i, 1, phi_i, 1);

    m_fields[0]->FwdTrans(phi_i, fieldcoeffs[2]);

    // // Compute \nabla \sigma_i \nabla phi_m and use it as point sources for
    // phi_e. This is equivalently achieved by removing all the point sources in
    // myelinnated fiber region.
    Array<OneD, NekDouble> phieforcing = ComputeMMFDiffusion(m_unitmovingframes, phi_m);

    // Only nonzero for node.
    Vmath::Vmul(nq, m_nodezone, 1, phieforcing, 1, phieforcing, 1);

    m_fields[0]->FwdTrans(phieforcing, fieldcoeffs[3]);

    WriteFld(outname, m_fields[0], fieldcoeffs, variables);
}

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

    StdRegions::ConstFactorMap factors;
    factors[StdRegions::eFactorTau] = m_Helmtau;

    factors[StdRegions::eFactorLambda] = m_Cn * m_Rf / lambda;
    if(nvar==1)
    {
        NekDouble betaratio = (m_ratio_re_ri + 1.0) / m_ratio_re_ri;
        factors[StdRegions::eFactorLambda] = m_Cn * m_Rf * betaratio / lambda;
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
    StdRegions::ConstFactorMap factors;
    factors[StdRegions::eFactorTau] = m_Helmtau;

    // m_ratio_re_ri determines the conductivity only when the variable is one
    NekDouble betaratio = (m_ratio_re_ri + 1.0) / m_ratio_re_ri;
    factors[StdRegions::eFactorLambda] = m_Cn * m_Rf * betaratio / lambda;

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
    StdRegions::ConstFactorMap factors;
    factors[StdRegions::eFactorTau] = m_Helmtau;

    // m_ratio_re_ri determines the conductivity only when the variable is one
    factors[StdRegions::eFactorLambda] = m_Cn * m_Rf / lambda;

    // SetBoundaryConditions(time);
    // SetMembraneBoundaryCondition();

    // Multiply 1.0/timestep
    Vmath::Smul(nq, -factors[StdRegions::eFactorLambda], inarray[0], 1,
                m_fields[0]->UpdatePhys(), 1);

    m_fields[0]->HelmSolve(m_fields[0]->GetPhys(), m_fields[0]->UpdateCoeffs(),
                        factors, m_varcoeff);
    m_fields[0]->BwdTrans(m_fields[0]->GetCoeffs(), outarray[0]);
    m_fields[0]->SetPhysState(true);
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
    m_neuron->TimeIntegrate(m_zoneindex[0], inarray[0], outarray[0], time, m_diameter, m_Temperature);

    for (int i=1; i < nvar; ++i)
    {
        Vmath::Vcopy(nq, m_neuron->GetNeuronSolution(i), 1, m_fields[i]->UpdatePhys(), 1);
    }

    for (unsigned int i = 0; i < m_stimulus.size(); ++i)
    {
        m_stimulus[i]->Update(m_excitezone, outarray, time);
    }
}

void MMFNeuralEP::DoOdeRhsNeuralEP1D(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time)
{
    int nvar = m_fields.size();
    int nq   = m_fields[0]->GetNpoints();

    for (int i=0; i<nvar; ++i)
    {
        outarray[i] = Array<OneD, NekDouble>(nq, 0.0);
    }

    // Compute the reaction function divided by Cm or Cn.
    m_neuron->TimeIntegrate(m_zoneindex[0], inarray[0], outarray[0], time, m_diameter, m_Temperature);

    for (unsigned int i = 0; i < m_stimulus.size(); ++i)
    {
        m_stimulus[i]->Update(m_excitezone, outarray, time);
    }

    // // Add it to the RHS
    // Vmath::Vadd(nq, RHSstimulus[0], 1, outarray[0], 1, outarray[0], 1);

        // switch (m_projectionType)
        // {
        //     case MultiRegions::eDiscontinuous:
        //     {
        //         std::string diffName;

        //         // Do not forwards transform initial condition
        //         m_homoInitialFwd = false;

        //         m_session->LoadSolverInfo("DiffusionType", diffName, "LDG");
        //         m_diffusion = SolverUtils::GetDiffusionFactory().CreateInstance(
        //             diffName, diffName);
        //         m_diffusion->SetFluxVector(&MMFNeuralEP::GetFluxVector, this);
        //         m_diffusion->InitObject(m_session, m_fields);
        //         break;
        //     }

        //     case MultiRegions::eGalerkin:
        //     case MultiRegions::eMixed_CG_Discontinuous:
        //     {
        //         if (m_explicitDiffusion)
        //         {
        //             ASSERTL0(false, "Explicit Galerkin diffusion not set up.");
        //         }
        //     }
        // }

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
            Vmath::Smul(nq, 1.0/m_Rf, &outarrayDiff[i][0], 1, &outarrayDiff[i][0], 1);
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

    for (int i=0; i<nvar; ++i)
    {
        outarray[i] = Array<OneD, NekDouble>(nq, 0.0);
    }

    // Compute the reaction function divided by Cm or Cn.
     m_neuron->TimeIntegrate(m_zoneindex[0], inarray[0], outarray[0], time, m_diameter, m_Temperature);

    // // Add Stimulus
     for (unsigned int j = 0; j < m_stimulus.size(); ++j)
     {
         m_stimulus[j]->Update(m_excitezone, outarray, time);
     }

    if (m_explicitDiffusion)
    {
        // Laplacian only to the first variable
        Array<OneD, NekDouble> Laplacian(nq);
        WeakDGMMFDiffusion(0, inarray[0], Laplacian, time);

        Vmath::Smul(nq, 1.0 / (m_Cn * m_Rf), Laplacian, 1, Laplacian, 1);
        Vmath::Vadd(nq, &Laplacian[0], 1, &outarray[0][0], 1, &outarray[0][0], 1);
    }
}


void MMFNeuralEP::DoOdeRhsNeuralEP2Dbi(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time)
{
    int nvar = m_fields.size();
    int nq   = m_fields[0]->GetNpoints();

    for (int i=0; i<nvar; ++i)
    {
        outarray[i] = Array<OneD, NekDouble>(nq);
    }

    // Compute the reaction function divided by Cm or Cn.
    m_neuron->TimeIntegrate(m_zoneindex[0], inarray[0], outarray[0], time, m_diameter, m_Temperature);

    // Add Stimulus
    for (unsigned int j = 0; j < m_stimulus.size(); ++j)
    {
        m_stimulus[j]->Update(m_excitezone, outarray, time);
    }

    // Compute phi_e to satisfy the following equation
    // \nabla \cdot ( (\signa_e + \sigma_i) \nabla \phi_e) = - \nabla \cdot
    // (\sigma_i \nabla \phi_m)
    Array<OneD, NekDouble> extcurrent(nq,0.0);    
    if(m_ExtCurrentType == eEphaptic)
    {
        Array<OneD, NekDouble> phie = Computephie(inarray[0]);
        
        // Compute (1/C_n/r) * \nabla^2 \phi_e
        extcurrent = ComputeMMFDiffusion(m_movingframes, phie);

        // extra current caused by phi_e only occurs in the intracellular space: / (m_Cn * m_Rf)
        Vmath::Vmul(nq, m_intrazone, 1, extcurrent, 1, extcurrent, 1);
        Vmath::Smul(nq, 1.0/(m_Cn * m_Rf), extcurrent, 1, extcurrent, 1);
    }

    // add divergence of phie to the current
    // Vmath::Svtvp(nq, 1.0 / (m_Cn * m_Rf), &extcurrent[0], 1, &outarray[0][0], 1, &outarray[0][0], 1);
    Vmath::Vadd(nq, &extcurrent[0], 1, &outarray[0][0], 1, &outarray[0][0], 1);

    if (m_explicitDiffusion)
    {
        // Laplacian only to the first variable
        Array<OneD, NekDouble> Laplacian(nq);
        WeakDGMMFDiffusion(0, inarray[0], Laplacian, time);

        Vmath::Smul(nq, 1.0 / (m_Cn * m_Rf), Laplacian, 1, Laplacian, 1);
        Vmath::Vadd(nq, &Laplacian[0], 1, &outarray[0][0], 1, &outarray[0][0], 1);
    }
}

// Input: phi_m
// output: phi_e (m_fields[1]->UpdatePhys()) and outarray (1/C_n/r) * \nabla^2 \phi_e
// Compute phi_e from the given distribution of phi_m
// \nabla \cdot ( (1 + \rho) \mathbf{e}_1 + \mathbf{e}_2 ) ( \nabla \phi_e ))
//                         = - \nabla \cdot \mathbf{e}_1 \nabla \phi_m
Array<OneD, NekDouble> MMFNeuralEP::Computephie(
    const Array<OneD, const NekDouble> &phim)
{
    int nq = m_fields[0]->GetNpoints();

    Array<OneD, NekDouble> outarray(nq);

    // Solve the Poisson equation: \nabla (\sigma_e + \sigma_i ) phi_e = \nabla
    // \sigma_i \nabla phi_m
    StdRegions::ConstFactorMap phiefactors;
    phiefactors[StdRegions::eFactorTau]    = m_Helmtau;
    phiefactors[StdRegions::eFactorLambda] = 0.0;

    // // Compute \nabla \sigma_i \nabla phi_m and use it as point sources for
    // phi_e. This is equivalently achieved by removing all the point sources in
    // myelinnated fiber region.
    Array<OneD, NekDouble> phimLaplacian = ComputeMMFDiffusion(m_movingframes, phim);

    // Only nonzero for node.
    Vmath::Vmul(nq, m_nodezone, 1, phimLaplacian, 1, phimLaplacian, 1);

    // Compute phie distribution
    // SetBoundaryConditions(0.0);

    // Helsolve with pure Neumann boundary condition
    Vmath::Sadd(nq, -1.0 * AvgInt(phimLaplacian), phimLaplacian, 1, phimLaplacian, 1);
    Vmath::Smul(nq, -1.0, phimLaplacian, 1, m_fields[1]->UpdatePhys(), 1);

    // Compute  \nabla \cdot ( (1 + \rho) \mathbf{e}_1 + \mathbf{e}_2 ) ( \nabla \phi_e ))
    //                         = - \nabla \cdot \mathbf{e}_1 \nabla \phi_m
    m_fields[1]->HelmSolve(m_fields[1]->GetPhys(), m_fields[1]->UpdateCoeffs(), phiefactors, m_phievarcoeff);
    m_fields[1]->BwdTrans(m_fields[1]->GetCoeffs(), m_fields[1]->UpdatePhys());
    m_fields[1]->SetPhysState(true);

    outarray = m_fields[1]->GetPhys();

    return outarray;
}

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
                m_stimulus[i]->Update(m_excitezone, tmp, initialtime);
                Vmath::Vmul(nq, m_intrazone, 1, tmp[0], 1, tmp[0], 1);
                m_fields[0]->SetPhys(tmp[0]);
            }
        }

        default:
        {
            EquationSystem::v_SetInitialConditions(initialtime, false);
            break;
        }
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
    for (int i = 0; i < 1; ++i)
    {
        inarray[i] = Array<OneD, NekDouble>(nq);
        Fwd[i]     = Array<OneD, NekDouble>(nTracePts);

        Vmath::Vcopy(nq, &m_fields[i]->GetPhys()[0], 1, &inarray[i][0], 1);
        m_fields[i]->ExtractTracePhys(inarray[i], Fwd[i]);
    }

    Array<OneD, NekDouble> x0(nq);
    Array<OneD, NekDouble> x1(nq);
    Array<OneD, NekDouble> x2(nq);

    m_fields[0]->GetCoords(x0, x1, x2);

    Array<OneD, NekDouble> x0Fwd(nTracePts);
    Array<OneD, NekDouble> x1Fwd(nTracePts);
    Array<OneD, NekDouble> x2Fwd(nTracePts);

    m_fields[0]->ExtractTracePhys(x0, x0Fwd);
    m_fields[0]->ExtractTracePhys(x1, x1Fwd);
    m_fields[0]->ExtractTracePhys(x2, x2Fwd);

    Array<OneD, NekDouble> NodeZoneFwd(nTracePts);

    Array<OneD, NekDouble> NodeZoneDouble(nq);

    for (int i=0; i<nq; ++i)
    {
        NodeZoneDouble[i] = 1.0 * m_zoneindex[0][i];
    }

    m_fields[0]->ExtractTracePhys(NodeZoneDouble, NodeZoneFwd);

    // loop over Boundary Regions
    for (int n = 0; n < m_fields[0]->GetBndConditions().size(); ++n)
    {
        // Wall Boundary Condition
        if (boost::iequals(m_fields[0]->GetBndConditions()[n]->GetUserDefined(), "Membrane"))
        {
            int id1, id2, index, npts;

            const Array<OneD, const int> &traceBndMap = m_fields[0]->GetTraceBndMap();

            for (int e = 0; e < m_fields[0]->GetBndCondExpansions()[n]->GetExpSize(); ++e)
            {
                npts = m_fields[0]
                                ->GetBndCondExpansions()[n]
                                ->GetExp(e)
                                ->GetTotPoints();
                id1 = m_fields[0]->GetBndCondExpansions()[n]->GetPhys_Offset(e);
                id2 = m_fields[0]->GetTrace()->GetPhys_Offset(traceBndMap[cnt + e]);

                for (int i=0;i<npts;++i)
                {
                    index = id2+i;
                    std::cout << "n = " << n << ", e = " << e << ", id2 = " << index << ", Zone = " << NodeZoneFwd[index] 
                    << " at x = " << x0Fwd[index] << ", y = " << x1Fwd[index] << std::endl;
                }

                // Pure Neumann boundary condtiion
                Vmath::Vcopy(npts, &Fwd[0][id2], 1,
                            &(m_fields[0]
                                ->GetBndCondExpansions()[n]
                                ->UpdatePhys())[id1], 1);
            }
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
    boost::ignore_unused(physarray);

    int nq = GetTotPoints();
    //  int nvariables = physarray.size();
    int nTracePts = GetTraceNpoints();

    int id1, id2, npts;
    int eMax = m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExpSize();

    const Array<OneD, const int> &traceBndMap = m_fields[0]->GetTraceBndMap();

    Array<OneD, NekDouble> x0(nq);
    Array<OneD, NekDouble> x1(nq);
    Array<OneD, NekDouble> x2(nq);

    m_fields[0]->GetCoords(x0, x1, x2);

    Array<OneD, NekDouble> x0tmp(nTracePts);
    Array<OneD, NekDouble> x1tmp(nTracePts);
    Array<OneD, NekDouble> x2tmp(nTracePts);

    m_fields[0]->ExtractTracePhys(x0, x0tmp);
    m_fields[0]->ExtractTracePhys(x1, x1tmp);
    m_fields[0]->ExtractTracePhys(x2, x2tmp);

    for (int e = 0; e < eMax; ++e)
    {
        npts = m_fields[0]
                         ->GetBndCondExpansions()[bcRegion]
                         ->GetExp(e)
                         ->GetTotPoints();
        id1 = m_fields[0]->GetBndCondExpansions()[bcRegion]->GetPhys_Offset(e);
        id2 = m_fields[0]->GetTrace()->GetPhys_Offset(traceBndMap[cnt + e]);

        // Pure Neumann boundary condtiion
        Vmath::Vcopy(npts, &Fwd[0][id2], 1,
                    &(m_fields[0]
                        ->GetBndCondExpansions()[bcRegion]
                        ->UpdatePhys())[id1], 1);
    }
}

void MMFNeuralEP::v_EvaluateExactSolution(unsigned int field,
                                          Array<OneD, NekDouble> &outfield,
                                          const NekDouble time)
{
    EquationSystem::v_EvaluateExactSolution(field, outfield, time);
}

void MMFNeuralEP::v_GenerateSummary(SolverUtils::SummaryList &s)
{
    int nvar = m_fields.size();
    int nq = GetTotPoints();

    MMFSystem::v_GenerateSummary(s);

    SolverUtils::AddSummaryItem(s, "NeuralEPType",
                                NeuralEPTypeMap[m_NeuralEPType]);
    SolverUtils::AddSummaryItem(s, "ExtCurrentType", ExtCurrentTypeMap[m_ExtCurrentType]);
    SolverUtils::AddSummaryItem(s, "GlobalSysSoln", m_session->GetSolverInfo("GlobalSysSoln"));
    // SolverUtils::AddSummaryItem(s, "TimeMapScheme", TimeMapTypeMap[m_TimeMapScheme]);
    // SolverUtils::AddSummaryItem(s, "TimeMapStart", m_TimeMapStart);
    // SolverUtils::AddSummaryItem(s, "TimeMapEnd", m_TimeMapEnd);

    SolverUtils::AddSummaryItem(s, "Fiber Type", FiberTypeMap[m_fiberType]);
    SolverUtils::AddSummaryItem(s, "Fiber Length", m_fiberlen);
    SolverUtils::AddSummaryItem(s, "Node Length", m_nodelen);
    SolverUtils::AddSummaryItem(s, "Myelin Length", m_myelinlen);

    SolverUtils::AddSummaryItem(s, "Total_Number_Node", m_totNode);
    SolverUtils::AddSummaryItem(s, "Element_per_Node", m_elemperNode);
    SolverUtils::AddSummaryItem(s, "Element_per_Myelin", m_elemperMyel);

    SolverUtils::AddSummaryItem(s, "Temperature", m_Temperature);
    SolverUtils::AddSummaryItem(s, "diameter", m_diameter);
    SolverUtils::AddSummaryItem(s, "Helmtau", m_Helmtau);
    SolverUtils::AddSummaryItem(s, "nq", nq);
    SolverUtils::AddSummaryItem(s, "npts", m_npts);

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
