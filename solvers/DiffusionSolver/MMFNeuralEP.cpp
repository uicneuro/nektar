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
    for (int j = 0; j < m_expdim; ++j)
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
    
    // Resting potential     NekDouble m_phimrest, m_phimTol, m_dudtTol;
    m_session->LoadParameter("phimrest", m_phimrest, 80.0);
    m_session->LoadParameter("phimTol", m_phimTol, 10.0);
    m_session->LoadParameter("dphimdtTol", m_dphimdtTol, 0.01);

    // NeuralEP paramter on temperature
    m_session->LoadParameter("Temperature", m_Temperature, 24.0);
    m_session->LoadParameter("diameter", m_diameter, 0.01);

    m_session->LoadParameter("g-ratio", m_gratio, 0.8);
    m_session->LoadParameter("relativefiberratio", m_relfiberratio, 0.8);
    m_session->LoadParameter("radiusfiberbundle", m_radiusfiberbundle, 0.01);
    m_session->LoadParameter("radiusaxon", m_radiusaxon, 0.0001);

    m_session->LoadParameter("FiberWidth", m_fiberwidth, 0.01);
    m_session->LoadParameter("FiberGap", m_fibergap, 0.01);
    m_session->LoadParameter("FiberHeightDiff", m_fiberheightdiff, 0.0);

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

        case eTestRanvierNode:
        {
            m_npts = m_fields[0]->GetTotPoints(0);

            m_zoneindex = Array<OneD, Array<OneD, int>>(1);
            m_zoneindex[0] = Array<OneD, int>(nq, 0); 
            
            m_zoneindex[0] = TestRanvierIndex();

            m_NeuralCm    = Array<OneD, Array<OneD, NekDouble>>(1);
            m_NeuralCm[0] = Array<OneD, NekDouble>(nq, 1.0 / m_Cm);
            break;
        }

        case eNeuralHelmTest:
        {
            // load zoneindexfile
            m_npts = m_fields[0]->GetTotPoints(0);

            m_session->LoadParameter("InnerboxEnd", m_InnerboxEnd, 0);

            int index;
            m_zoneindex = Array<OneD, Array<OneD, int>>(1);
            m_zoneindex[0] = Array<OneD, int>(nq, 0);
            for (int i = m_InnerboxEnd; i < m_fields[0]->GetExpSize(); ++i)
            {
                for (int j = 0; j < m_fields[0]->GetTotPoints(i); ++j)
                    {
                        index = m_fields[0]->GetPhys_Offset(i) + j;
                        m_zoneindex[0][index] = -2;
                    }
            }

            // Get the first and last index of the excitation zone [1,2]
            // SetUpDomainZone(m_zoneindex[0], m_excitezone, m_nodezone, m_intrazone, m_extrazone);
            m_NeuralCm    = Array<OneD, Array<OneD, NekDouble>>(1);
            m_NeuralCm[0] = ComputeConductivity(m_zoneindex[0]);
            break;
        }

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

        case eNeuralEP2Dmono:
        case eNeuralEP2Dbi:
        {   
            m_zoneindex = Array<OneD, Array<OneD, int>>(1);
            m_zoneindex[0]  = Array<OneD, int>(nq, 1); 

            m_npts = m_fields[0]->GetTotPoints(0);
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
    
    CheckMovingFrames(m_movingframes);

    // Construct unitmovingframes
    std::string MMFdirStr;
    m_session->LoadSolverInfo("MMFDir", MMFdirStr, "LOCAL");
    m_MMFdir = FindMMFdir(MMFdirStr);

    Array<OneD, Array<OneD, NekDouble>> unitAniStrength(m_expdim);
    for (int j = 0; j < m_expdim; ++j)
    {
        unitAniStrength[j] = Array<OneD, NekDouble>(nq, 1.0);
    }

    for (int i=0; i<nq; ++i)
    {
        for (int j = 0; j < m_expdim; ++j)
        {
            if(m_zoneindex[0][i] == -2)
            {
                unitAniStrength[j][i] = 0.0;
            }
        }
    }

    std::cout << "Unit Moving frames are generated with " << MMFdirStr << " direction ===============" << std::endl;
    SetUpMovingFrames(m_MMFdir, unitAniStrength, m_unitmovingframes);
    
    CheckMovingFrames(m_unitmovingframes);

    // Construct phiemovingframes 
    std::string phieMMFdirStr;
    m_session->LoadSolverInfo("phieMMFDir", phieMMFdirStr, "TangentY");
    SpatialDomains::GeomMMF phieMMFdir = FindMMFdir(phieMMFdirStr);

    Array<OneD, Array<OneD, NekDouble>> phieAniStrength(m_expdim);
    for (int j = 0; j < m_expdim; ++j)
    {
        phieAniStrength[j] = Array<OneD, NekDouble>(nq, 1.0);
    }

    std::cout << "Phie Moving frames are generated with " << phieMMFdirStr 
    << " direction ===============" << std::endl;

    SetUpMovingFrames(phieMMFdir, phieAniStrength, m_phiemovingframes);

    switch (m_NeuralEPType)
    {
        case eTestRanvierNode:
        {
            // Array<OneD, Array<OneD, NekDouble>> sigma_i(m_expdim);
            // Array<OneD, Array<OneD, NekDouble>> sigma_e(m_expdim);
            // for (int j = 0; j < m_expdim; ++j)
            // {
            //     sigma_i[j] = Array<OneD, NekDouble>(nq, 1.0);
            //     sigma_e[j] = Array<OneD, NekDouble>(nq, 1.0);
            // }

            // NekDouble axoncrossA = m_pi*m_radiusaxon*m_radiusaxon;
            // for (int i = 0; i<nq; ++i)
            // {
            //         // Node zone
            //         if(m_zoneindex[0][i]>=0)
            //         {
            //             sigma_e[0][i] = 1.0/m_ratio_re_ri;
            //             sigma_e[1][i] = 1.0/m_ratio_re_ri;
            //         }

            //         // Myelin zone
            //         else if(m_zoneindex[0][i]==-1)
            //         {
            //             sigma_e[0][i] = 1.0/m_ratio_re_ri;
            //             sigma_e[1][i] = 1.0/m_ratio_re_ri;
            //         }

            //         else if(m_zoneindex[0][i]==-2)
            //         {
            //             sigma_e[0][i] = 1.0/axoncrossA;
            //             sigma_e[1][i] = 1.0/axoncrossA;
            //         }
            // }

            // std::cout << "Max sigma_e_1  = "
            //             << Vmath::Vmax(nq, sigma_e[0], 1)
            //             << ", sigma_e_2 = "
            //             << Vmath::Vmax(nq, sigma_e[1], 1)
            //             << ", Min sigma_e_1 = "
            //             << Vmath::Vmin(nq, sigma_e[0], 1)
            //             << ", sigma_e_2 = "
            //             << Vmath::Vmin(nq, sigma_e[1], 1) << std::endl;


            // for (int j = 0; j < m_expdim; ++j)
            // {
            //     Vmath::Vadd(nq, sigma_i[j], 1, sigma_e[j], 1, m_phieAniStrength[j], 1);
            // }

            // // m_phieMF = \sigma_i + \sigma_e
            // for (int i = 0; i < nq; ++i)
            // {
            //     for (int j = 0; j < m_expdim; ++j)
            //     {
            //         for (int k = 0; k < m_spacedim; ++k)
            //         {
            //                 m_phiemovingframes[j][k * nq + i] = sqrt(m_phieAniStrength[j][i]) * m_phiemovingframes[j][k * nq + i];
            //         }
            //     }
            // }
            
            // CheckMovingFrames(m_phiemovingframes);

            // break;            
        }

        case eNeuralEP2Dbi:
        {
            Array<OneD, Array<OneD, NekDouble>> sigma_i(m_expdim);
            Array<OneD, Array<OneD, NekDouble>> sigma_e(m_expdim);
            for (int j = 0; j < m_expdim; ++j)
            {
                sigma_i[j] = Array<OneD, NekDouble>(nq, 1.0);
                sigma_e[j] = Array<OneD, NekDouble>(nq, 1.0);
            }

            // Compute sigma_i
            // Node: 1.0, Myelin: m_Cn / m_Cm, Extraspace: 0.0
            for (int i = 0; i<nq; ++i)
            {
                 for (int j = 0; j < m_expdim; ++j)
                 {
                    if(m_zoneindex[0][i]==-2)
                    {
                        sigma_i[j][i] = 0.0;
                    }

                    else if (m_zoneindex[0][i]==-1)
                    {
                        sigma_i[j][i] = m_AnisotropyStrength;
                    }
                 }
            }

            std::cout << "Max sigma_i_1  = "
                        << Vmath::Vmax(nq, sigma_i[0], 1)
                        << ", sigma_i_2 = "
                        << Vmath::Vmax(nq, sigma_i[1], 1)
                        << ", Min sigma_i_1 = "
                        << Vmath::Vmin(nq, sigma_i[0], 1)
                        << ", sigma_i_2 = "
                        << Vmath::Vmin(nq, sigma_i[1], 1) << std::endl;

            // Compute sigma_e
            NekDouble axoncrossA = m_pi*m_radiusaxon*m_radiusaxon;
            NekDouble PhieMultFactor = m_radiusaxon*m_radiusaxon/(m_relfiberratio*m_gratio*m_gratio*m_radiusfiberbundle*m_radiusfiberbundle);

            // Multiplication factor for fiber bundle of radius R = (a^2/(\rho g^2 R^2))
            std::cout << "Fiber bundle mag factor = " << PhieMultFactor << std::endl;

            for (int i = 0; i<nq; ++i)
            {
                    // Node zone
                    if(m_zoneindex[0][i]>=0)
                    {
                        sigma_e[0][i] = 1.0/m_ratio_re_ri;
                        sigma_e[1][i] = 1.0/m_ratio_re_ri;
                    }

                    // Myelin zone
                    else if(m_zoneindex[0][i]==-1)
                    {
                        sigma_e[0][i] = m_AnisotropyStrength/m_ratio_re_ri;
                        sigma_e[1][i] = m_AnisotropyStrength/m_ratio_re_ri;
                    }

                    else if(m_zoneindex[0][i]==-2)
                    {
                        sigma_e[0][i] = PhieMultFactor / axoncrossA;
                        sigma_e[1][i] = PhieMultFactor / axoncrossA;
                    }
            }

            std::cout << "Max sigma_e_1  = "
                        << Vmath::Vmax(nq, sigma_e[0], 1)
                        << ", sigma_e_2 = "
                        << Vmath::Vmax(nq, sigma_e[1], 1)
                        << ", Min sigma_e_1 = "
                        << Vmath::Vmin(nq, sigma_e[0], 1)
                        << ", sigma_e_2 = "
                        << Vmath::Vmin(nq, sigma_e[1], 1) << std::endl;

            for (int j = 0; j < m_expdim; ++j)
            {
                Vmath::Vadd(nq, sigma_i[j], 1, sigma_e[j], 1, m_phieAniStrength[j], 1);
            }
            
            // m_phieMF = \sigma_i + \sigma_e
            for (int i = 0; i < nq; ++i)
            {
                for (int j = 0; j < m_expdim; ++j)
                {
                    for (int k = 0; k < m_spacedim; ++k)
                    {
                            m_phiemovingframes[j][k * nq + i] = sqrt(m_phieAniStrength[j][i]) * m_phiemovingframes[j][k * nq + i];
                    }
                }
            }

            std::cout << "================================================ " << std::endl;
            std::cout << "Max phieAnistrength_1  = "
                        << Vmath::Vmax(nq, m_phieAniStrength[0], 1)
                        << ", phieAnistrength_2 = "
                        << Vmath::Vmax(nq, m_phieAniStrength[1], 1)
                        << ", Min phieAnistrength 1 = "
                        << Vmath::Vmin(nq, m_phieAniStrength[0], 1)
                        << ", phieAnistrength 2 = "
                        << Vmath::Vmin(nq, m_phieAniStrength[1], 1) << std::endl;
            std::cout << "================================================ " << std::endl;
            
            PlotPhieMF(sigma_i, sigma_e, m_phieAniStrength);

            break;
        }

        // case eNeuralHelmTest:
        // {
        //     // Set up for m_phiemovingframe Poisson solver
        //     std::string phieMMFdirStr;
        //     m_session->LoadSolverInfo("phieMMFDir", phieMMFdirStr, "TangentY");
        //     SpatialDomains::GeomMMF phieMMFdir = FindMMFdir(phieMMFdirStr);

        //     m_phieAniStrength = Array<OneD, Array<OneD, NekDouble>>(m_expdim);

        //     NekDouble rada = 0.01;
        //     NekDouble axoncrossA = m_pi*rada*rada;
        //     m_phieAniStrength[0] = Array<OneD, NekDouble>(nq, 1.0/m_ratio_re_ri);
        //     m_phieAniStrength[1] = Array<OneD, NekDouble>(nq, 1.0/m_ratio_re_ri/axoncrossA);

        //     int index;
        //     for (int i = 0; i < m_InnerboxEnd; ++i)
        //     {
        //         for (int j = 0; j < m_fields[0]->GetTotPoints(i); ++j)
        //             {
        //                 index = m_fields[0]->GetPhys_Offset(i) + j;
        //                 m_phieAniStrength[0][index] += 1.0;
        //                 m_phieAniStrength[1][index] += 1.0;
        //             }
        //     }

        //     // Multiplication factor for fiber bundle of radius R = (a^2/(\rho g^2 R^2))
        //     NekDouble fbundlemag = m_radiusaxon*m_radiusaxon/(m_relfiberratio*m_gratio*m_gratio*m_radiusfiberbundle*m_radiusfiberbundle);

        //     std::cout << "Phie Moving frames are generated with " << phieMMFdirStr << " direction ===============" << std::endl;
        //     std::cout << "Fiber bundle mag factor = " << fbundlemag << std::endl;

        //     SetUpMovingFrames(phieMMFdir, m_phieAniStrength, m_phiemovingframes);

        //     std::cout << "phiAniStr = " << RootMeanSquare(m_phieAniStrength[0]) << ", " << RootMeanSquare(m_phieAniStrength[1]) << std::endl;

        //     // m_phieMF = \sigma_i + \sigma_e
        //     for (int i = 0; i < nq; ++i)
        //     {
        //         for (int j = 0; j < m_expdim; ++j)
        //         {
        //             for (int k = 0; k < m_spacedim; ++k)
        //             {
        //                     m_phiemovingframes[j][k * nq + i] = sqrt(m_phieAniStrength[j][i]) * m_phiemovingframes[j][k * nq + i];
        //             }
        //         }
        //     }

        //     std::cout << "Max phieAnistrength_1  = "
        //                 << Vmath::Vmax(nq, m_phieAniStrength[0], 1)
        //                 << ", phieAnistrength_2 = "
        //                 << Vmath::Vmax(nq, m_phieAniStrength[1], 1)
        //                 << ", Min phieAnistrength 1 = "
        //                 << Vmath::Vmin(nq, m_phieAniStrength[0], 1)
        //                 << ", phieAnistrength 2 = "
        //                 << Vmath::Vmin(nq, m_phieAniStrength[1], 1) << std::endl;

        //     break;
        // }

        default:
         break;
    }

    // Check moving frames
    CheckNodeZoneMF(m_zoneindex, m_movingframes, m_phiemovingframes);

    if (m_explicitDiffusion)
    {
        m_ode.DefineImplicitSolve(&MMFNeuralEP::DoNullSolve, this);
        m_ode.DefineProjection(&MMFNeuralEP::DoOdeProjection, this);
    }

    else
    {
        switch (m_NeuralEPType)
        {
            // case eNeuralEPPT:
            // {
            //     // ComputeVarCoeff1D(m_movingframes, m_varcoeff);
            //     m_ode.DefineImplicitSolve(&MMFNeuralEP::DoNullSolve, this);
            //     m_ode.DefineOdeRhs(&MMFNeuralEP::DoOdeRhsNeuralEPPT, this);
            //     break;
            // }

            case eTestRanvierNode:
            {
                std::cout << std::endl;
                std::cout << "Generating m_varcoeff ================================= " << std::endl;
                ComputeVarCoeff2D(m_movingframes, m_varcoeff);

                std::cout << "Generating m_phievarcoeff ================================= " << std::endl;
                ComputeVarCoeff2D(m_phiemovingframes, m_phievarcoeff);
                std::cout << std::endl;

                int cnt = 0;
                Array<OneD, NekDouble> forcing(nq,0.0);
                for (int i=0;i<nq; ++i)
                {
                    if(m_zoneindex[0][i]>=0)
                    {
                        forcing[i] = 1.0;
                        cnt++;
                    }
                }
                std::cout << "cnt = " << cnt << std::endl;

                StdRegions::ConstFactorMap phiefactors;
                phiefactors[StdRegions::eFactorTau]    = m_Helmtau;
                phiefactors[StdRegions::eFactorLambda] = 0.0;

                Array<OneD, NekDouble> solution(nq,0.0);

                Vmath::Sadd(nq, -1.0 * AvgInt(forcing), forcing, 1, m_fields[1]->UpdatePhys(), 1);
                m_fields[1]->HelmSolve(m_fields[1]->GetPhys(), m_fields[1]->UpdateCoeffs(), phiefactors, m_phievarcoeff);
                m_fields[1]->BwdTrans(m_fields[1]->GetCoeffs(), m_fields[1]->UpdatePhys());
                m_fields[1]->SetPhysState(true);

                solution = m_fields[1]->GetPhys();

                std::cout << "solution, max = " << Vmath::Vmax(nq, solution, 1) << ", min = " << Vmath::Vmin(nq, solution, 1) << std::endl;

                PlotHelmSolve(forcing,solution);

                wait_on_enter();

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

            case eNeuralEP2Dmono:
            {
                ComputeVarCoeff2D(m_movingframes, m_varcoeff);
                m_ode.DefineImplicitSolve(
                    &MMFNeuralEP::DoImplicitSolveNeuralEP2Dmono, this);
                m_ode.DefineOdeRhs(&MMFNeuralEP::DoOdeRhsNeuralEP2Dmono, this);
                break;
            }

            case eNeuralHelmTest:
            case eNeuralEP2Dbi:
            {
                std::cout << std::endl;
                std::cout << "Generating m_varcoeff ================================= " << std::endl;
                ComputeVarCoeff2D(m_movingframes, m_varcoeff);

                std::cout << "Generating m_phievarcoeff ================================= " << std::endl;
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
        Array<OneD, NekDouble> x0(nq);
        Array<OneD, NekDouble> x1(nq);
        Array<OneD, NekDouble> x2(nq);

        m_fields[0]->GetCoords(x0, x1, x2);

        int index;
        Array<OneD, NekDouble> phim(nq,0.0);
        Array<OneD, NekDouble> phieexact(nq,0.0);

        NekDouble coeffa = (1.0+1.0)/(1.0+1.0+1.0/m_ratio_re_ri+1.0/m_ratio_re_ri);
        for (int i = 0; i < m_InnerboxEnd; ++i)
        {
            for (int j = 0; j < m_fields[0]->GetTotPoints(i); ++j)
            {
                index = m_fields[0]->GetPhys_Offset(i) + j;
                phim[index] = cos(m_pi*x0[index])*cos(m_pi*x1[index]);
                phieexact[index] = -coeffa*cos(m_pi*x0[index])*cos(m_pi*x1[index]);
            }
        }

        // NekDouble coeffe;
        for (int i = m_InnerboxEnd; i < m_fields[0]->GetExpSize(); ++i)
        {
            for (int j = 0; j < m_fields[0]->GetTotPoints(i); ++j)
            {
                index = m_fields[0]->GetPhys_Offset(i) + j;

                // coeffe = coeffa;
                // if(x0[index]>1.0)
                // {
                //     coeffe = coeffa / (x0[index]-1) / (x0[index]-1);
                // }

                // else if (x0[index]<-1.0)
                // {
                //     coeffe = coeffa / (x0[index]+1) / (x0[index]+1);
                // }

                phieexact[index] = -coeffa*cos(m_pi*x0[index])*cos(m_pi*x1[index]);
            }
        }
        std::cout << "phim = " << RootMeanSquare(phim) << std::endl;

        StdRegions::ConstFactorMap phiefactors;
        phiefactors[StdRegions::eFactorTau]    = m_Helmtau;
        phiefactors[StdRegions::eFactorLambda] = 0.0;

        // // Compute \nabla \sigma_i \nabla phi_m and use it as point sources for
        // phi_e. This is equivalently achieved by removing all the point sources in
        // myelinnated fiber region.
        Array<OneD, NekDouble> phimLaplacian(nq);
        phimLaplacian = ComputeMMFDiffusion(m_movingframes, phim);
        
        NekDouble extFieldStr = 1.0;
        Vmath::Sadd(nq, -1.0 * AvgInt(phimLaplacian), phimLaplacian, 1, phimLaplacian, 1);
        Vmath::Smul(nq, -1.0 * extFieldStr, phimLaplacian, 1, m_fields[1]->UpdatePhys(), 1);

        m_fields[1]->HelmSolve(m_fields[1]->GetPhys(), m_fields[1]->UpdateCoeffs(), phiefactors, m_phievarcoeff);
        m_fields[1]->BwdTrans(m_fields[1]->GetCoeffs(), m_fields[1]->UpdatePhys());
        m_fields[1]->SetPhysState(true);

        Array<OneD, NekDouble> outarray(nq);
        outarray = m_fields[1]->GetPhys();

        Vmath::Sadd(nq, -1.0 * AvgInt(outarray), outarray, 1, outarray, 1);                    

        std::cout << "phie = " << RootMeanSquare(outarray) << std::endl;

        // Plotting

        int nvar    = 4;
        int ncoeffs = m_fields[0]->GetNcoeffs();

        std::string outname1;
        outname1 = m_sessionName + "_result.chk";

        std::vector<Array<OneD, NekDouble>> fieldcoeffs(nvar);
        for (int i = 0; i < nvar; ++i)
        {
            fieldcoeffs[i] = Array<OneD, NekDouble>(ncoeffs);
        }

        std::vector<std::string> variables(nvar);
        variables[0] = "phim";
        variables[1] = "phie";
        variables[2] = "phieexact";
        variables[3] = "phieerror";

        Array<OneD, NekDouble> phieInerr(nq,0.0);
        Array<OneD, NekDouble> phieOuterr(nq,0.0);
        Array<OneD, NekDouble> phieerr(nq,0.0);

        for (int i = 0; i < m_InnerboxEnd; ++i)
        {
            for (int j = 0; j < m_fields[0]->GetTotPoints(i); ++j)
            {
                index = m_fields[0]->GetPhys_Offset(i) + j;
                phieInerr[index] = outarray[index] - phieexact[index];
            }
        }

        for (int i = m_InnerboxEnd; i < m_fields[0]->GetExpSize(); ++i)
        {
            for (int j = 0; j < m_fields[0]->GetTotPoints(i); ++j)
            {
                index = m_fields[0]->GetPhys_Offset(i) + j;
                phieOuterr[index] = outarray[index] - phieexact[index];
            }
        }


        std::cout << "Error In = " << RootMeanSquare(phieInerr) << ", out = " << RootMeanSquare(phieOuterr) << std::endl;

        Vmath::Vadd(nq, phieInerr, 1, phieOuterr, 1, phieerr, 1);

        m_fields[0]->FwdTransLocalElmt(phim, fieldcoeffs[0]);
        m_fields[0]->FwdTransLocalElmt(outarray, fieldcoeffs[1]);
        m_fields[0]->FwdTransLocalElmt(phieexact, fieldcoeffs[2]);
        m_fields[0]->FwdTransLocalElmt(phieerr, fieldcoeffs[3]);

        WriteFld(outname1, m_fields[0], fieldcoeffs, variables);

        wait_on_enter();
    }

}

/**
 *
 */
MMFNeuralEP::~MMFNeuralEP()
{
}

void MMFNeuralEP::DoOdeProjection(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time)
{
    int i;
    int nvariables = inarray.size();
    SetBoundaryConditions(time);

    switch (m_projectionType)
    {
        case MultiRegions::eDiscontinuous:
        {
            // Just copy over array
            int npoints = GetNpoints();

            for (i = 0; i < nvariables; ++i)
            {
                Vmath::Vcopy(npoints, inarray[i], 1, outarray[i], 1);
            }
            break;
        }
        case MultiRegions::eGalerkin:
        case MultiRegions::eMixed_CG_Discontinuous:
        {
            Array<OneD, NekDouble> coeffs(m_fields[0]->GetNcoeffs());

            for (i = 0; i < nvariables; ++i)
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

Array<OneD, int> MMFNeuralEP::IndexNodeZone2D(
    const FiberType fiber)
    {
        int nq   = GetTotPoints();

        Array<OneD, int> outarray(nq);
        switch(fiber)
        {
            case eSingleLinear:
            {
                outarray = SingleLinearIndex(m_fiberwidth, m_nodelen, m_myelinlen, m_totNode);
                break;
            }

            case eDoubleLinear:
            {
                outarray = DoubleLinearIndex(m_fiberwidth, m_fibergap, m_fiberheightdiff, m_nodelen, m_myelinlen, m_totNode);
                break;
            }

            default:
            break;
        }

        return outarray;
    }

Array<OneD, int> MMFNeuralEP::TestRanvierIndex()
{
    int nq   = GetTotPoints();

    // Compute the center of the element as the coordinate of each grid.
    Array<OneD, NekDouble> xcell(nq);
    Array<OneD, NekDouble> ycell(nq);
    Array<OneD, NekDouble> zcell(nq);

    Getcellavg(xcell,ycell,zcell);

    Array<OneD, int> outarray(nq, -1);

    // NekDouble nodestart, nodeend;
    NekDouble fiberleft = -0.25;
    NekDouble fiberright = 0.25;

    int cnt=0;
    int cmt=0;
    int cet=0;
    for (int i=0; i<nq; ++i)
    {
        if((xcell[i]>fiberleft) && (xcell[i]<fiberright))
        {
            if( (ycell[i]>fiberleft) && (ycell[i]<fiberright) )
            {
                outarray[i] = 0;
                cnt++;
            }

            else
            {
                outarray[i] = -1;
                cmt++;
            }
        }

        else
        {
            outarray[i] = -2;
            cet++;
        }
    }

    std::cout << "cnt = " << cnt << ", cmt = " << cmt << ", cet = " << cet << std::endl;

    return outarray;
}



Array<OneD, int> MMFNeuralEP::SingleLinearIndex(
    const NekDouble fiberwidth, 
    const NekDouble nodelen, 
    const NekDouble myelinlen,
    const int Nnode)
{
    int nq   = GetTotPoints();

    // Compute the center of the element as the coordinate of each grid.

    Array<OneD, NekDouble> xcell(nq);
    Array<OneD, NekDouble> ycell(nq);
    Array<OneD, NekDouble> zcell(nq);

    Getcellavg(xcell,ycell,zcell);

    Array<OneD, int> outarray(nq, -1);

    NekDouble nodestart, nodeend;
    NekDouble fiberleft = -0.5*fiberwidth;
    NekDouble fiberright = 0.5*fiberwidth;

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

            for (int k=0; k<Nnode; ++k)
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

Array<OneD, int> MMFNeuralEP::DoubleLinearIndex(
    const NekDouble fiberwidth, 
    const NekDouble fibergap, 
    const NekDouble fiberheightdiff, 
    const NekDouble nodelen, 
    const NekDouble myelinlen,
    const int Nnode)
{
    int nq   = GetTotPoints();

    // Compute the center of the element as the coordinate of each grid.

    Array<OneD, NekDouble> xcell(nq);
    Array<OneD, NekDouble> ycell(nq);
    Array<OneD, NekDouble> zcell(nq);

    Getcellavg(xcell,ycell,zcell);

    Array<OneD, int> outarray(nq, -1);

    NekDouble fiber1left = -0.5*fiberwidth;
    NekDouble fiber1right = 0.5*fiberwidth;

    NekDouble fiber2left =  0.5*fiberwidth+fibergap;
    NekDouble fiber2right = 1.5*fiberwidth+fibergap;

    NekDouble nodestart, nodeend;
    NekDouble nodebottom, nodetop;

    for (int i=0; i<nq; ++i)
    {
        // first fiber
        if((xcell[i]>fiber1left) && (xcell[i]<fiber1right))
        {
            for (int k=0; k<2; ++k)
            {
                nodebottom = k*nodelen;
                nodetop = (k+1)*nodelen;
                if( (ycell[i]>=nodebottom) && (ycell[i]<=nodetop) )
                {
                    outarray[i] = k;
                }
            }

            for (int k=0; k<Nnode; ++k)
            {
                nodestart = 2.0*nodelen + myelinlen + k * (myelinlen + nodelen);
                nodeend = 2.0*nodelen + (k+1) * (myelinlen + nodelen);
                if( (ycell[i]>=nodestart) && (ycell[i]<=nodeend) )
                {
                    outarray[i] = k+2;
                }
            }
        }

        else if((xcell[i]>fiber2left) && (xcell[i]<fiber2right))
        {
            for (int k=0; k<2; ++k)
            {
                nodebottom = fiberheightdiff + k*nodelen;
                nodetop = fiberheightdiff + (k+1)*nodelen;
                if( (ycell[i]>=nodebottom) && (ycell[i]<=nodetop) )
                {
                    outarray[i] = 100+k;
                }
            }

            for (int k=0; k<Nnode; ++k)
            {
                nodestart = fiberheightdiff + 2.0*nodelen + myelinlen + k * (myelinlen + nodelen);
                nodeend = fiberheightdiff + 2.0*nodelen + (k+1) * (myelinlen + nodelen);
                if( (ycell[i]>=nodestart) && (ycell[i]<=nodeend) )
                {
                    outarray[i] = 100+k+2;
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

void MMFNeuralEP::Getcellavg(
    Array<OneD, NekDouble> &xcell, 
    Array<OneD, NekDouble> &ycell, 
    Array<OneD, NekDouble> &zcell)
{
    int nq   = GetTotPoints();

    Array<OneD, NekDouble> x0(nq);
    Array<OneD, NekDouble> x1(nq);
    Array<OneD, NekDouble> x2(nq);

    m_fields[0]->GetCoords(x0, x1, x2);

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
}


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
        // first node excitation zone
        if( (zoneindex[i]==0) || (zoneindex[i]==1) )
        {
            excitezone[i] = 1.0 / m_Cn;
            extcnt++;
        }

        // second node excitation zone
        if( (zoneindex[i]==100) || (zoneindex[i]==101) )
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
            outarray[i] = 1.0 / m_Cn;
            cnte++;
        }
    }

    std::cout << "v_InitObject: Node = " << cntn
    << ", Myelin = " << cntm << ", extracell = " << cnte
    << std::endl;

    return outarray;
}

void MMFNeuralEP::CheckNodeZoneMF(
    const Array<OneD, const Array<OneD, int>> &NodeZone,
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    const Array<OneD, const Array<OneD, NekDouble>> &phiemovingframes)
{
    boost::ignore_unused(NodeZone);

    int nq = GetTotPoints();

    Array<OneD, NekDouble> x0(nq);
    Array<OneD, NekDouble> x1(nq);
    Array<OneD, NekDouble> x2(nq);

    m_fields[0]->GetCoords(x0, x1, x2);

    int i, j, index, npts;;
    NekDouble xp, yp, e1mag, e2mag;
    NekDouble phie1mag, phie2mag;

    for (i = 0; i < m_fields[0]->GetExpSize(); ++i)
    {
        npts = m_fields[0]->GetTotPoints(i);
        xp = 0.0;
        yp = 0.0;
        e1mag = 0.0;
        e2mag = 0.0;

        phie1mag = 0.0;
        phie2mag = 0.0;
        for (j = 0; j < npts; ++j)
        {
            index = m_fields[0]->GetPhys_Offset(i) + j;

            xp += x0[index];
            yp += x1[index];

            e1mag = e1mag +
                    (movingframes[0][index] * movingframes[0][index] +
                     movingframes[0][nq + index] * movingframes[0][nq + index]);
            e2mag = e2mag +
                    (movingframes[1][index] * movingframes[1][index] +
                     movingframes[1][nq + index] * movingframes[1][nq + index]);

            phie1mag = phie1mag +
                    (phiemovingframes[0][index] * phiemovingframes[0][index] +
                     phiemovingframes[0][nq + index] * phiemovingframes[0][nq + index]);
            phie2mag = phie2mag +
                    (phiemovingframes[1][index] * phiemovingframes[1][index] +
                     phiemovingframes[1][nq + index] * phiemovingframes[1][nq + index]);
        }

        // dx = x0[m_fields[0]->GetPhys_Offset(i)] - x0[m_fields[0]->GetPhys_Offset(i)+npts-1];
        // dy = x1[m_fields[0]->GetPhys_Offset(i)] - x1[m_fields[0]->GetPhys_Offset(i)+npts-1];
        // dist = sqrt(dx*dx+dy*dy);

        e1mag = sqrt(e1mag / npts);
        e2mag = sqrt(e2mag / npts);

        phie1mag = sqrt(phie1mag / npts);
        phie2mag = sqrt(phie2mag / npts);

        xp = (xp / npts);
        yp = (yp / npts);

        std::cout << "Elemid = " << i << ", Nodeid = " << NodeZone[0][index]
                << ", x = " << xp << ", y = " << yp 
                << ", e1mag = " << e1mag << ", e2mag = " << e2mag
                << ", phie1mag = " << phie1mag << ", phie2mag = " << phie2mag << std::endl;
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

    return outarray;
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
    for (i = 0; i < nvariables; ++i)
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

    Array<OneD, NekDouble> velmag(nq, 0.0);
    Array<OneD, NekDouble> velocity(m_spacedim * nq);

    Array<OneD, NekDouble> TimeMap(nq, 0.0);
    Array<OneD, NekDouble> dudtval(nq);
    Array<OneD, NekDouble> dudtvalHistory(nq, 0.0);

    NekDouble Maxphim;
    while (step < m_steps || m_time < m_fintime - NekConstants::kNekZeroTol)
    {
        // Save fields into fieldsold
        Vmath::Vcopy(nq, &fields[0][0], 1, &fields_old[0][0], 1);

        timer.Start();
        fields = m_intScheme->TimeIntegrate(step, m_timestep, m_ode);
        timer.Stop();

        m_time += m_timestep;
        elapsed = timer.TimePerTest(1);
        intTime += elapsed;
        cpuTime += elapsed;

       // Compute TimeMap
       // dudtsign: wavefront = -1.0, waveback = 1.0
        Maxphim = Vmath::Vamax(nq, fields[0], 1);
        Vmath::Vsub(nq, fields[0], 1, fields_old[0], 1, dudtval, 1);
        Vmath::Vmul(nq, m_intrazone, 1, dudtval, 1, dudtval, 1);

        Vmath::Smul(nq, 1.0 / (m_timestep * Maxphim), dudtval, 1, dudtval, 1);

        if ((m_TimeMapStart <= m_time) && (m_TimeMapEnd >= m_time))
        {
            ComputeNeuralTimeMap(m_time, m_phimrest, m_phimTol, m_dphimdtTol, m_zoneindex[0], dudtval, dudtvalHistory, fields[0], TimeMap);
        }

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
            // PrintoutFields(nvariables, fields, fulltext);
            std::cout << fulltext << "\n" << std::endl;

            PlotNeuralTimeMap(fields[0], TimeMap, nchk);
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


void MMFNeuralEP::PrintoutFields(const int nvar, const Array<OneD, const Array<OneD, NekDouble>> &fields, std::string &fulltext)
{
    int nq = GetTotPoints();

    Array<OneD, NekDouble> x0(nq);
    Array<OneD, NekDouble> x1(nq);
    Array<OneD, NekDouble> x2(nq);

    m_fields[0]->GetCoords(x0, x1, x2);

    // phim should be only defined in the intracellular space
    Array<OneD, NekDouble> phi_m; // = m_fields[0]->GetPhys();
    Vmath::Vmul(nq, m_intrazone, 1, fields[0], 1, phi_m, 1);

    NekDouble phimMax = Vmath::Vmax(nq, phi_m, 1);
    NekDouble phimMin = Vmath::Vmin(nq, phi_m, 1);

    int phimMaxid = Vmath::Imax(nq, phi_m, 1);
    int phimMinid = Vmath::Imin(nq, phi_m, 1);

    fulltext.append("phi_m, max: " + std::to_string(phimMax) + " at y = " + std::to_string(x1[phimMaxid]) );
    fulltext.append("./ phi_m, min: " + std::to_string(phimMin) + " at y = " + std::to_string(x1[phimMinid]) );

    fulltext.append("\n");


    Array<OneD, NekDouble> phimLap1 = ComputeMMFDiffusion(m_unitmovingframes, phi_m);
    Array<OneD, NekDouble> phimLap2 = ComputeMMFDiffusion(m_movingframes, phi_m);

    // Only nonzero for node.
    Vmath::Vmul(nq, m_nodezone, 1, phimLap1, 1, phimLap1, 1);
    Vmath::Vmul(nq, m_nodezone, 1, phimLap2, 1, phimLap2, 1);

    Array<OneD, NekDouble> phimLapdiff(nq);
    Vmath::Vsub(nq, phimLap1, 1, phimLap2, 1, phimLapdiff, 1);

   fulltext.append("Lap of phim Diff: Max = " + std::to_string(Vmath::Vmax(nq, phimLapdiff, 1)) + ", Min = " + std::to_string(Vmath::Vmin(nq, phimLapdiff, 1)) );
    fulltext.append("\n");

    if(nvar==2) 
    {
        // phi_e is defined at the node and extracellular space
        Array<OneD, NekDouble> phi_e = m_fields[1]->GetPhys();
        // Vmath::Vmul(nq, m_extrazone, 1, phi_e, 1, phi_e, 1);
        
        NekDouble phieMax = Vmath::Vmax(nq, phi_e, 1);
        NekDouble phieMin = Vmath::Vmin(nq, phi_e, 1);

        int phieMaxid = Vmath::Imax(nq, phi_e, 1);
        int phieMinid = Vmath::Imin(nq, phi_e, 1);

        fulltext.append("phi_e, max: " + std::to_string(phieMax) + " at y = " + std::to_string(x1[phieMaxid]) );
        fulltext.append(", phi_e, min: " + std::to_string(phieMin) + " at y = " + std::to_string(x1[phieMinid]) );
        fulltext.append("\n");

        // Compute (1/C_n/r) * \nabla^2 \phi_e
        Array<OneD, NekDouble> intcurrent = ComputeMMFDiffusion(m_movingframes, phi_m);

        // extra current caused by phi_e only occurs in the intracellular space: / (m_Cn * m_Rf)
        Vmath::Vmul(nq, m_intrazone, 1, intcurrent, 1, intcurrent, 1);
        Vmath::Smul(nq, 1.0/(m_Cn * m_Rf), intcurrent, 1, intcurrent, 1);

        NekDouble intcurrentMax = Vmath::Vmax(nq, intcurrent, 1);
        NekDouble intcurrentMin = Vmath::Vmin(nq, intcurrent, 1);

        Array<OneD, NekDouble> extcurrent = ComputeMMFDiffusion(m_movingframes, phi_e);

        Vmath::Vmul(nq, m_intrazone, 1, extcurrent, 1, extcurrent, 1);
        Vmath::Smul(nq, 1.0/(m_Cn * m_Rf), extcurrent, 1, extcurrent, 1);

        NekDouble extcurrentMax = Vmath::Vmax(nq, extcurrent, 1);
        NekDouble extcurrentMin = Vmath::Vmin(nq, extcurrent, 1);

        Array<OneD, NekDouble> totcurrent(nq);

        Vmath::Vadd(nq, extcurrent, 1, intcurrent, 1, totcurrent, 1);

        NekDouble totcurrentMax = Vmath::Vmax(nq, totcurrent, 1);
        NekDouble totcurrentMin = Vmath::Vmin(nq, totcurrent, 1);

        NekDouble totcurrentMaxid = Vmath::Imax(nq, totcurrent, 1);
        NekDouble totcurrentMinid = Vmath::Imin(nq, totcurrent, 1);

        NekDouble Maxratio = 100.0 * fabs(totcurrentMax) / fabs(intcurrentMax) ;
        NekDouble Minratio = 100.0 * fabs(totcurrentMin) / fabs(intcurrentMin) ;

        fulltext.append("Modified intcurrent, max: " + std::to_string(extcurrentMax) + " ( " + std::to_string(Maxratio) + " %) at y = " + std::to_string(x1[totcurrentMaxid]) );
        fulltext.append("Modified intcurrent, min: " + std::to_string(extcurrentMin) + " ( " + std::to_string(Minratio) + " %) at y = " + std::to_string(x1[totcurrentMinid]) );
        fulltext.append("\n");

        // DisplayNode2D(fulltext, fields);
    }

}

// ComputeNeuralTimeMap(m_time, m_urest, m_uTol, m_dudtTol, m_zoneindex[0], dudtval, dudtvalHistory, fields[0], TimeMap);
void MMFNeuralEP::ComputeNeuralTimeMap(const NekDouble time,
                                    const NekDouble phimrest,
                                    const NekDouble phimTol,
                                    const NekDouble dphimdtTol,
                                    const Array<OneD, const int> &zoneindex,
                                    const Array<OneD, const NekDouble> &dudt,
                                    Array<OneD, NekDouble> &dudtHistory,
                                    const Array<OneD, const NekDouble> &field,
                                    Array<OneD, NekDouble> &TimeMap)
{
    int nq = GetTotPoints();

    NekDouble fnewsum;

    NekDouble phimdiff;
    for (int i = 0; i < nq; ++i)
    {
        if(zoneindex[i]>-2)
        {
            phimdiff = field[i] - phimrest;
            // Only integrate of time if u > Tol, gradu > Tol, du/dt > 0
            if ((phimdiff > phimTol) && (dudt[i] > dphimdtTol))
            {
                // Gradient as the main weight
                fnewsum = dudt[i] + dudtHistory[i];

                if(fabs(fnewsum)>dphimdtTol)
                {
                    TimeMap[i] = (dudt[i] * time + dudtHistory[i] * TimeMap[i]) / fnewsum;
                }

                dudtHistory[i] += dudt[i];
            }
        }

        if ( (zoneindex[i] == 0) || (zoneindex[i] == 1))
        {
            TimeMap[i] = 0.0;
        }
    }
}


void MMFNeuralEP::PlotNeuralTimeMap(
    const Array<OneD, const NekDouble> &phim,
    const Array<OneD, const NekDouble> &TimeMap,
    const int nstep)
{
    int nvar    = 2;
    int nq      = m_fields[0]->GetTotPoints();
    int ncoeffs = m_fields[0]->GetNcoeffs();

    std::string outname1 = m_sessionName + "_timemap_" +
                           boost::lexical_cast<std::string>(nstep) + ".chk";

    std::vector<Array<OneD, NekDouble>> fieldcoeffs(nvar);
    for (int i = 0; i < nvar; ++i)
    {
        fieldcoeffs[i] = Array<OneD, NekDouble>(ncoeffs);
    }

    std::vector<std::string> variables(nvar);
    variables[0] = "phim";
    variables[1] = "TimeMap";

    m_fields[0]->FwdTransLocalElmt(phim, fieldcoeffs[0]);

    // index:0 -> u
    std::cout << "PlotTimeMap: phim: Max = " << Vmath::Vmax(nq, phim, 1)
                << ", Min = " << Vmath::Vmin(nq, phim, 1) << ", Time Map: Max = " << Vmath::Vmax(nq, TimeMap, 1)
                << ", Min = " << Vmath::Vmin(nq, TimeMap, 1) << std::endl;

    m_fields[0]->FwdTransLocalElmt(TimeMap, fieldcoeffs[1]);

    WriteFld(outname1, m_fields[0], fieldcoeffs, variables);
}


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

void MMFNeuralEP::PlotHelmSolve(const Array<OneD, const NekDouble> &forcing,
                                   const Array<OneD, const NekDouble> &solution)
{
    int nvar = 2;
    int ncoeffs = m_fields[0]->GetNcoeffs();
    // int nq      = m_fields[0]->GetTotPoints();

    std::string outname;
    outname = m_sessionName + "_helm.chk";

    std::vector<Array<OneD, NekDouble>> fieldcoeffs(nvar);
    for (int i = 0; i < nvar; ++i)
    {
        fieldcoeffs[i] = Array<OneD, NekDouble>(ncoeffs);
    }

    std::vector<std::string> variables(nvar);
    variables[0] = "forcing";
    variables[1] = "solution";

    m_fields[0]->FwdTransLocalElmt(forcing, fieldcoeffs[0]);
    m_fields[0]->FwdTransLocalElmt(solution, fieldcoeffs[1]);

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
    m_fields[0]->FwdTransLocalElmt(tmp, fieldcoeffs[0]);

    for (int i = 1; i < nvar; ++i)
    {
        tmp = m_neuron->GetNeuronSolution(i);
        m_fields[0]->FwdTransLocalElmt(tmp, fieldcoeffs[i]);
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
    Array<OneD, NekDouble> tmp(nq);    
    if(m_ExtCurrentType == eEphaptic)
    {
        Vmath::Smul(nq, -1.0, inarray[0], 1, tmp, 1);
        Array<OneD, NekDouble> phie = Computephie(tmp);
        
        // Compute (1/C_n/r) * \nabla^2 \phi_e
        extcurrent = ComputeMMFDiffusion(m_movingframes, phie);

        // extra current caused by phi_e only occurs in the intracellular space: / (m_Cn * m_Rf)
        Vmath::Vmul(nq, m_intrazone, 1, extcurrent, 1, extcurrent, 1);
        Vmath::Smul(nq, 1.0/(m_Cn * m_Rf), extcurrent, 1, extcurrent, 1);
    }

    // subtract the current from the divergence of phie
    Vmath::Vadd(nq, &outarray[0][0], 1, &extcurrent[0], 1, &outarray[0][0], 1);

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
    Array<OneD, NekDouble> phimLaplacian = ComputeMMFDiffusion(m_unitmovingframes, phim);

    // Only nonzero for node.
    Vmath::Vmul(nq, m_nodezone, 1, phimLaplacian, 1, phimLaplacian, 1);

    // Compute phie distribution
    // SetBoundaryConditions(0.0);

    // Compute  \nabla \cdot ( (1 + \rho) \mathbf{e}_1 + \mathbf{e}_2 ) ( \nabla \phi_e ))
    //                         = - \nabla \cdot \mathbf{e}_1 \nabla \phi_m
    Vmath::Sadd(nq, -1.0 * AvgInt(phimLaplacian), phimLaplacian, 1, m_fields[1]->UpdatePhys(), 1);
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

void MMFNeuralEP::PlotPhieMF(
    const Array<OneD, const Array<OneD, NekDouble>> &sigma_i,
    const Array<OneD, const Array<OneD, NekDouble>> &sigma_e,
    const Array<OneD, const Array<OneD, NekDouble>> &PhieAniStrength)
{
    int nvar    = 6;
    // int nq      = m_fields[0]->GetTotPoints();
    int ncoeffs = m_fields[0]->GetNcoeffs();

    std::string outname1 = m_sessionName + "_phieMF.chk";

    std::vector<Array<OneD, NekDouble>> fieldcoeffs(nvar);
    for (int i = 0; i < nvar; ++i)
    {
        fieldcoeffs[i] = Array<OneD, NekDouble>(ncoeffs);
    }

    std::vector<std::string> variables(nvar);
    variables[0] = "sigma_i[0]";
    variables[1] = "sigma_i[1]";
    variables[2] = "sigma_e[0]";
    variables[3] = "sigma_e[1]";
    variables[4] = "PhieAniStrength[0]";
    variables[5] = "PhieAniStrength[1]";

    // Compute the gradient of the time map
    m_fields[0]->FwdTransLocalElmt(sigma_i[0], fieldcoeffs[0]);
    m_fields[0]->FwdTransLocalElmt(sigma_i[1], fieldcoeffs[1]);

    m_fields[0]->FwdTransLocalElmt(sigma_e[0], fieldcoeffs[2]);
    m_fields[0]->FwdTransLocalElmt(sigma_e[1], fieldcoeffs[3]);

    m_fields[0]->FwdTransLocalElmt(PhieAniStrength[0], fieldcoeffs[4]);
    m_fields[0]->FwdTransLocalElmt(PhieAniStrength[1], fieldcoeffs[5]);

    WriteFld(outname1, m_fields[0], fieldcoeffs, variables);
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
    SolverUtils::AddSummaryItem(s, "TimeMapScheme", TimeMapTypeMap[m_TimeMapScheme]);
    SolverUtils::AddSummaryItem(s, "TimeMapStart", m_TimeMapStart);
    SolverUtils::AddSummaryItem(s, "TimeMapEnd", m_TimeMapEnd);

    SolverUtils::AddSummaryItem(s, "Fiber Type", FiberTypeMap[m_fiberType]);

    SolverUtils::AddSummaryItem(s, "FiberWidth", m_fiberwidth);
    SolverUtils::AddSummaryItem(s, "FiberGap", m_fibergap);
    SolverUtils::AddSummaryItem(s, "FiberHeightDiff", m_fiberheightdiff);

    SolverUtils::AddSummaryItem(s, "Node Length", m_nodelen);
    SolverUtils::AddSummaryItem(s, "Myelin Length", m_myelinlen);

    SolverUtils::AddSummaryItem(s, "Total_Number_Node", m_totNode);
    SolverUtils::AddSummaryItem(s, "Element_per_Node", m_elemperNode);
    SolverUtils::AddSummaryItem(s, "Element_per_Myelin", m_elemperMyel);

    SolverUtils::AddSummaryItem(s, "phimrest", m_phimrest);
    SolverUtils::AddSummaryItem(s, "phimTol", m_phimTol);
    SolverUtils::AddSummaryItem(s, "dphimdtTol", m_dphimdtTol);

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
