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

    m_pi       = 3.14159265358979323846;

    int nq   = GetTotPoints();

    // Derive AnisotropyStrength.
    m_AniStrength = Array<OneD, Array<OneD, NekDouble>> (m_expdim);
    m_phieAniStrength = Array<OneD, Array<OneD, NekDouble>> (m_expdim);

    for (int j = 0; j < m_expdim; ++j)
    {
        m_AniStrength[j] = Array<OneD, NekDouble>(nq, 1.0);
        m_phieAniStrength[j] = Array<OneD, NekDouble>(nq, 1.0);
    }

    // Conductance parameters
    m_session->LoadParameter("Chi", m_chi, 28.0);
    m_session->LoadParameter("Cm", m_capMembrane, 0.125);

    // Helmsolver parameter
    m_session->LoadParameter("Helmtau", m_Helmtau, 1.0);
    
    // Resting potential     NekDouble m_phimrest, m_phimTol, m_dudtTol;
    m_session->LoadParameter("phimrest", m_phimrest, 80.0);
    // m_session->LoadParameter("phimTol", m_phimTol, 10.0);
    // m_session->LoadParameter("phieTol", m_phieTol, 10.0);
    // m_session->LoadParameter("dphimdtTol", m_dphimdtTol, 1.0);

    // NeuralEP paramter on temperature
    m_session->LoadParameter("Temperature", m_Temperature, 24.0);
    m_session->LoadParameter("axondiameter", m_axondiameter, 0.01);

   // Relative Extracellular resistance: 1 < \beta < 10
    m_session->LoadParameter("ratio_re_ri", m_ratio_re_ri, 1.0);
    m_session->LoadParameter("AnisotropyStrength", m_AnisotropyStrength, 4.0);

    m_session->LoadParameter("FiberAngle", m_fiberangle, 0.0);
    m_session->LoadParameter("FiberWidth", m_fiberwidth, 0.01);
    m_session->LoadParameter("FiberGap", m_fibergap, 0.01);
    m_session->LoadParameter("FiberCurvature", m_fibercurvature, 0.8);

    m_session->LoadParameter("NodeLength", m_nodelen, 0.01);
    m_session->LoadParameter("MyelinLength", m_myelinlen, 0.2);

    m_session->LoadParameter("Total_Number_Node", m_totNode, 3);
    m_session->LoadParameter("nodeinitdown", m_nodeinitdown, 0.01);
    m_session->LoadParameter("nodeinitup", m_nodeinitup, 0.02);

    m_numfiber = 1;

    m_session->LoadParameter("fiber1order", m_fiber1order, 1); // 0: down to up, 1: up to down
    m_session->LoadParameter("fiber1left", m_fiber1left, 0.01);
    m_session->LoadParameter("fiber1right", m_fiber1right, 0.02);

    m_session->LoadParameter("fiber2order", m_fiber2order, 1); // 0: down to up, 1: up to down
    m_session->LoadParameter("fiber2left", m_fiber2left, 0.0);
    m_session->LoadParameter("fiber2right", m_fiber2right, 0.0);

    if( (fabs(m_fiber2left)>0.0) && (fabs(m_fiber2right)>0.0) )
    {
        m_numfiber = 2;
    }

    m_session->LoadParameter("fiber3order", m_fiber3order, 1); // 0: down to up, 1: up to down
    m_session->LoadParameter("fiber3left", m_fiber3left, 0.0);
    m_session->LoadParameter("fiber3right", m_fiber3right, 0.0);

    m_session->LoadParameter("bundleleft", m_bundleleft, 0.0);
    m_session->LoadParameter("bundleright", m_bundleright, 0.05);

    m_session->LoadParameter("g-ratio", m_gratio, 0.8);
    m_session->LoadParameter("relativefiberratio", m_relfiberratio, 0.8);
    m_session->LoadParameter("radiusfiberbundle", m_radiusfiberbundle, 0.01);

    m_session->LoadParameter("Diffext", m_Diffext, 1e-6);

    // m_session->LoadParameter("radiusaxon", m_radiusaxon, 0.01);

    NekDouble axoncrossA = m_pi*m_axondiameter*m_axondiameter;
    NekDouble PhieMultFactor = m_axondiameter*m_axondiameter/(m_relfiberratio*m_gratio*m_gratio*m_radiusfiberbundle*m_radiusfiberbundle);

    // 1.0 /(m_pi * m_relfiberratio*m_gratio*m_gratio*m_radiusfiberbundle*m_radiusfiberbundle)
    m_phiefactor = PhieMultFactor / axoncrossA;

    if( (fabs(m_fiber3left)>0.0) && (fabs(m_fiber3right)>0.0) )
    {
        m_numfiber = 3;
    }

    m_fiberorder = Array<OneD, int>(m_numfiber);
    m_fiberleft = Array<OneD, NekDouble>(m_numfiber);
    m_fiberright = Array<OneD, NekDouble>(m_numfiber);

    m_fiberorder[0] = m_fiber1order;
    m_fiberleft[0] = m_fiber1left;
    m_fiberright[0] = m_fiber1right;

    if(m_numfiber>1)
    {
       m_fiberorder[1] = m_fiber2order;
       m_fiberleft[1] = m_fiber2left;
       m_fiberright[1] = m_fiber2right;
    }

    if(m_numfiber>2)
    {
       m_fiberorder[2] = m_fiber3order;
       m_fiberleft[2] = m_fiber3left;
       m_fiberright[2] = m_fiber3right;
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

    // Either incorporating external current effect or not.
    if (m_session->DefinesSolverInfo("FiberType"))
    {
        std::string FiberTypeStr;
        FiberTypeStr = m_session->GetSolverInfo("FiberType");
        for (int i = 0; i < (int)SIZE_FiberType; ++i)
        {
            if (boost::iequals(FiberTypeMap[i], FiberTypeStr))
            {
                m_FiberType = (FiberType)i;
                break;
            }
        }
    }
    else
    {
        m_FiberType = (FiberType)0;
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

    m_AnisotropyStrength = m_Cn / m_Cm;

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

        case eNeuralHelmSolveSingle:
        case eNeuralHelmSolveDuo:
        {
            m_npts = m_fields[0]->GetTotPoints(0);

            // m_zoneindex = Array<OneD, Array<OneD, int>>(1);
            // m_zoneindex[0] = Array<OneD, int>(nq, 1); 

            m_zoneindexfiber = Array<OneD, Array<OneD, int>>(m_numfiber);
            for (int n=0; n<m_numfiber; ++n)
            {
                m_zoneindexfiber[n]  = Array<OneD, int>(nq, 0); 
            }
            
            if(m_NeuralEPType==eNeuralHelmSolveSingle)
            {
                m_zoneindexfiber[0] = TestRanvierSingleIndex();
            }

            else if(m_NeuralEPType==eNeuralHelmSolveDuo)
            {
                m_zoneindexfiber[0] = TestRanvierDuoIndex();
            }

            // Setup: excitezone, intrazone, extrazone, followed by ploting the zones.
            SetUpDomainZone(m_zoneindexfiber, m_excitezonefiber, 
                            m_intrazonefiber, m_nodezone, m_myelinzone);

            m_extrazone = Array<OneD, NekDouble>(nq) ;
            m_intrazone = Array<OneD, NekDouble>(nq) ;

            Vmath::Vadd(nq, m_nodezone, 1, m_myelinzone, 1, m_intrazone, 1);
            Array<OneD, NekDouble> allone(nq, 1.0);
            Vmath::Vsub(nq, allone, 1, m_myelinzone, 1, m_extrazone, 1);

            // m_NeuralCm    = Array<OneD, Array<OneD, NekDouble>>(1);
            // m_NeuralCm[0] = Array<OneD, NekDouble>(nq, 1.0 / m_Cm);
            m_NeuralCm    = Array<OneD, Array<OneD, NekDouble>>(1);
            m_NeuralCm[0] = Array<OneD, NekDouble>(nq,0.0);

            for (int i = 0; i < nq; ++i)
            {
                // Ranvier node zone
                if (m_zoneindexfiber[0][i] >= 0)

                {
                    m_NeuralCm[0][i] = 1.0 / m_Cn;
                }

                // Myelin zone
                else if (m_zoneindexfiber[0][i] == -1)
                {
                    m_NeuralCm[0][i] = 1.0 / m_Cm;
                }

                // Extracellular space: \sigma_i = m_ratio_re_ri * \sigma_e
                else if (m_zoneindexfiber[0][i] == -2)
                {
                    m_NeuralCm[0][i] = 1.0 / m_Cn;
                }
            }

            break;
        }

        case eNeuralEP1D: 
        {
            // Ranvier node zone: 0: Myelin, 1: node
            m_nfibers  = 1;
            m_zoneindexfiber = Array<OneD, Array<OneD, int>>(m_nfibers);
           // m_zoneindexfiber[0] = IndexNodeZone1D(m_fields[0], m_totNode, m_elemperNode, m_elemperMyel); 

              // Get the first and last index of the excitation zone [1,2]
           // SetUpDomainZone(m_zoneindex[0], m_excitezone, m_nodezone, m_intrazone, m_extrazone);

            m_NeuralCm    = Array<OneD, Array<OneD, NekDouble>>(m_nfibers);
            m_NeuralCm[0] = ComputeConductivity(m_zoneindexfiber[0]);
            break;
        }

        case eNeuralEP2Dmono:
        case eNeuralEP2Dbi:
        {   
            m_npts = m_fields[0]->GetTotPoints(0);

            m_zoneindexfiber = Array<OneD, Array<OneD, int>>(m_numfiber);
            for (int n=0; n<m_numfiber; ++n)
            {
                m_zoneindexfiber[n]  = Array<OneD, int>(nq, 0); 
            }

            IndexNodeZone2D(m_numfiber, m_totNode, m_nodelen, 
                            m_myelinlen, m_nodeinitdown, m_nodeinitup,
                            m_fiberleft, m_fiberright, m_fiberorder, m_zoneindexfiber);

            // If all node zone = index = 1 for myelin. 
            if(m_MediumType==eAllNode)
            {
                for (int n=0; n<m_numfiber; ++n)
                {
                    for (int i=0; i<nq; ++i)
                    {
                         if( m_zoneindexfiber[n][i] == -1)
                         {
                            m_zoneindexfiber[n][i] = 1;
                         }
                    }
                }
            }

            // Construction ZoneIndex for all fibers;
            int tmp;
            m_zoneindex = Array<OneD, int>(nq, -2);
            for (int i=0; i<nq; ++i)
            {
                for (int n=0; n<m_numfiber; ++n)
                {
                    tmp = m_zoneindexfiber[n][i];
                    if(tmp>m_zoneindex[i])
                    {
                        m_zoneindex[i] = tmp;
                    }
                }
            }

            m_excitezonefiber = Array<OneD, Array<OneD, NekDouble>>(m_numfiber);
            m_intrazonefiber = Array<OneD, Array<OneD, NekDouble>>(m_numfiber);
            for (int n=0; n<m_numfiber; ++n)
            {
                m_excitezonefiber[n] = Array<OneD, NekDouble>(nq, 0.0);
                m_intrazonefiber[n] = Array<OneD, NekDouble>(nq, 0.0);
            }

            // Get the first and last index of the excitation zone [1,2]intra
            SetUpDomainZone(m_zoneindexfiber, m_excitezonefiber, 
                            m_intrazonefiber, m_nodezone, m_myelinzone);

            m_extrazone = Array<OneD, NekDouble>(nq) ;
            m_intrazone = Array<OneD, NekDouble>(nq) ;
            m_outerzone = Array<OneD, NekDouble>(nq) ;

            Vmath::Vadd(nq, m_nodezone, 1, m_myelinzone, 1, m_intrazone, 1);
            Array<OneD, NekDouble> allone(nq, 1.0);
            Vmath::Vsub(nq, allone, 1, m_myelinzone, 1, m_extrazone, 1);
            Vmath::Vsub(nq, allone, 1, m_intrazone, 1, m_outerzone, 1);

        // Plotting
            int nvar    = 2;
            int ncoeffs = m_fields[0]->GetNcoeffs();

            std::string outname1 = m_sessionName + "_zone.chk";

            std::vector<Array<OneD, NekDouble>> fieldcoeffs(nvar);
            for (int i = 0; i < nvar; ++i)
            {
                fieldcoeffs[i] = Array<OneD, NekDouble>(ncoeffs);
            }

            std::vector<std::string> variables(nvar);

            variables[0] = "intrazone";
            m_fields[0]->FwdTransLocalElmt(m_intrazone, fieldcoeffs[0]);

            variables[1] = "extrazone";
            m_fields[0]->FwdTransLocalElmt(m_extrazone, fieldcoeffs[1]);

            WriteFld(outname1, m_fields[0], fieldcoeffs, variables);

            m_NeuralCmfiber = ComputeConductivity(m_zoneindexfiber);

            break;
        }

        default:
            break;
    }
    // Stimulus
    m_stimulus = NeuralStimulus::LoadStimuli(m_session, m_fields[0]);

    // Derive AnisotropyStrength.
    // SetUpBiAnisotropy(m_zoneindex[0], m_NeuralCm, m_AniStrength);
    for (int j = 0; j < m_expdim; ++j)
    {
        m_AniStrength[j] = Array<OneD, NekDouble>(nq, 1.0);
        Vmath::Smul(nq, m_Cn, &m_NeuralCmfiber[0], 1, &m_AniStrength[j][0], 1);
    }

    MMFSystem::MMFInitObject(m_AniStrength);
    CheckMovingFrames(m_movingframes);

    Array<OneD, Array<OneD, NekDouble>> phieAniStrength(m_expdim);
    Array<OneD, Array<OneD, NekDouble>> phiediffAniStrength(m_expdim);
    for (int j = 0; j < m_expdim; ++j)
    {
        phieAniStrength[j] = Array<OneD, NekDouble>(nq, 1.0);
        phiediffAniStrength[j] = Array<OneD, NekDouble>(nq, 1.0);
    }

    // Construct phiemovingframes 
    std::string phieMMFdirStr;
    m_session->LoadSolverInfo("phieMMFDir", phieMMFdirStr, "TangentY");
    SpatialDomains::GeomMMF phieMMFdir = FindMMFdir(phieMMFdirStr);

    std::cout << "Phie Moving frames are generated with " << phieMMFdirStr 
    << " direction ===============" << std::endl;

    SetUpMovingFrames(phieMMFdir, m_phieAniStrength, m_phiemovingframes);

    // Construct phiediffmovingframes
    std::string MMFdirStr;
    m_session->LoadSolverInfo("MMFDir", MMFdirStr, "TangentY");
    m_MMFdir = FindMMFdir(MMFdirStr);

    std::cout << std::endl;
    SetUpMovingFrames(m_MMFdir, phiediffAniStrength, m_phiediffmovingframes);

    Array<OneD, Array<OneD, NekDouble>> sigma_i(m_expdim);
    Array<OneD, Array<OneD, NekDouble>> sigma_e(m_expdim);
    Array<OneD, Array<OneD, NekDouble>> sigma_eM(m_expdim);
    for (int j = 0; j < m_expdim; ++j)
    {
        sigma_i[j] = Array<OneD, NekDouble>(nq, 0.0);
        sigma_e[j] = Array<OneD, NekDouble>(nq, 0.0);
        sigma_eM[j] = Array<OneD, NekDouble>(nq, 0.0);
    }

    switch (m_NeuralEPType)
    {
        case eNeuralHelmSolveSingle:
        case eNeuralHelmSolveDuo:
        case eNeuralEP2Dbi:
        {
            ComputeRegionalSigma(m_zoneindex, sigma_i, sigma_e, sigma_eM);
            break;            
        }

        default:
         break;
    }

    for (int j = 0; j < m_expdim; ++j)
    {
        Vmath::Vadd(nq, sigma_i[j], 1, sigma_eM[j], 1, phieAniStrength[j], 1);
    }

    for (int j = 0; j < m_expdim; ++j)
    {
        Vmath::Vcopy(nq, sigma_e[j], 1, phiediffAniStrength[j], 1);
    }
    
    std::cout << "================================================ " << std::endl;
    std::cout << "Max phieAnistrength_1  = "
                << Vmath::Vmax(nq, phieAniStrength[0], 1)
                << ", phieAnistrength_2 = "
                << Vmath::Vmax(nq, phieAniStrength[1], 1)
                << ", Min phieAnistrength 1 = "
                << Vmath::Vmin(nq, phieAniStrength[0], 1)
                << ", phieAnistrength 2 = "
                << Vmath::Vmin(nq, phieAniStrength[1], 1) << std::endl;
    std::cout << "================================================ " << std::endl;

    std::cout << "================================================ " << std::endl;
    std::cout << "Max phiediffAnistrength_1  = "
                << Vmath::Vmax(nq, phiediffAniStrength[0], 1)
                << ", phiediffAnistrength_2 = "
                << Vmath::Vmax(nq, phiediffAniStrength[1], 1)
                << ", Min phiediffAnistrength 1 = "
                << Vmath::Vmin(nq, phiediffAniStrength[0], 1)
                << ", phiediffAnistrength 2 = "
                << Vmath::Vmin(nq, phiediffAniStrength[1], 1) << std::endl;
    std::cout << "================================================ " << std::endl;

    // m_phieMF = \sigma_i + \sigma_e
    for (int i = 0; i < nq; ++i)
    {
        for (int j = 0; j < m_expdim; ++j)
        {
            for (int k = 0; k < m_spacedim; ++k)
            {
                    m_phiediffmovingframes[j][k * nq + i] = sqrt(phiediffAniStrength[j][i]) * m_phiediffmovingframes[j][k * nq + i];
                    m_phiemovingframes[j][k * nq + i] = sqrt(phieAniStrength[j][i]) * m_phiemovingframes[j][k * nq + i];
            }
        }
    }

    std::cout << "Constructing phiemovingframes" << std::endl;
    CheckMovingFrames(m_phiemovingframes);

    std::cout << "Constructing phiediffmovingframes" << std::endl;
    CheckMovingFrames(m_phiediffmovingframes);


    PlotAnisotropy(m_AniStrength, m_phieAniStrength);

    // Check moving frames
    if(m_numfiber==2)
    {
        CheckNodeZoneMF(m_movingframes, m_phiemovingframes);
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
            // case eNeuralEPPT:
            // {
            //     // ComputeVarCoeff1D(m_movingframes, m_varcoeff);
            //     m_ode.DefineImplicitSolve(&MMFNeuralEP::DoNullSolve, this);
            //     m_ode.DefineOdeRhs(&MMFNeuralEP::DoOdeRhsNeuralEPPT, this);
            //     break;
            // }

            case eNeuralHelmSolveSingle:
            case eNeuralHelmSolveDuo:
            {
                std::cout << std::endl;
                std::cout << "Generating m_varcoeff ================================= " << std::endl;
                ComputeVarCoeff2D(m_movingframes, m_varcoeff);

                std::cout << "Generating m_phievarcoeff ================================= " << std::endl;
                ComputeVarCoeff2D(m_phiemovingframes, m_phievarcoeff);
                std::cout << std::endl;

                 Array<OneD, NekDouble> forcing(nq,0.0);
                for (int i=0;i<nq; ++i)
                {
                    if(m_zoneindexfiber[0][i]>=0)
                    {
                        forcing[i] = -100.0;
                     }
                }

                // Compute phim
                StdRegions::ConstFactorMap phimfactors;
                phimfactors[StdRegions::eFactorTau]    = m_Helmtau;

                NekDouble lambda;
                m_session->LoadParameter("Helmlambda", lambda, 0.001);
                phimfactors[StdRegions::eFactorLambda] = m_Cn * m_Rf / lambda;

                Array<OneD, NekDouble> phim(nq,0.0);

                NekDouble intforcing = AvgInt(forcing);

                Vmath::Sadd(nq, -1.0 * intforcing, forcing, 1, m_fields[0]->UpdatePhys(), 1);
                m_fields[0]->HelmSolve(m_fields[0]->GetPhys(), m_fields[0]->UpdateCoeffs(), phimfactors, m_varcoeff);
                m_fields[0]->BwdTrans(m_fields[0]->GetCoeffs(), m_fields[0]->UpdatePhys());
                m_fields[0]->SetPhysState(true);

                phim = m_fields[0]->GetPhys();

                std::cout << "phim, max = " << Vmath::Vmax(nq, phim, 1) << ", min = " << Vmath::Vmin(nq, phim, 1) << std::endl;

                Array<OneD, NekDouble> phimintra(nq);
                Array<OneD, NekDouble> phimextra(nq);

                Array<OneD, NekDouble> x0(nq);
                Array<OneD, NekDouble> x1(nq);
                Array<OneD, NekDouble> x2(nq);

                m_fields[0]->GetCoords(x0, x1, x2);

                NekDouble phim_max = Vmath::Vmax(nq, phim, 1);
                for (int i=0;i<nq;++i)
                {
                    phimextra[i] = (m_zoneindexfiber[0][i]==-2) ? phim[i]/phim_max : 0.0;
                }
                    // if (m_zoneindexfiber[0][i]==-2)
                    // {
                    //     phimextra[i] = phim[i]/Vmath::Vmax(nq, phim, 1);
                    // }

                    // else{
                    //     phimextra[i] = 0.0;
                    // }
                    // }

                std::cout << "phm in ex_zone: L2err = " << RootMeanSquare(phimextra) << ", Linf = " << Vmath::Vamax(nq, phimextra, 1) << std::endl;

                // Compute phie
                StdRegions::ConstFactorMap phiefactors;
                phiefactors[StdRegions::eFactorTau]    = m_Helmtau;
                phiefactors[StdRegions::eFactorLambda] = 0.0;

                Array<OneD, NekDouble> phie(nq,0.0);

                Vmath::Sadd(nq, -1.0 * intforcing, forcing, 1, m_fields[1]->UpdatePhys(), 1);
                m_fields[1]->HelmSolve(m_fields[1]->GetPhys(), m_fields[1]->UpdateCoeffs(), phiefactors, m_phievarcoeff);
                m_fields[1]->BwdTrans(m_fields[1]->GetCoeffs(), m_fields[1]->UpdatePhys());
                m_fields[1]->SetPhysState(true);

                phie = m_fields[1]->GetPhys();

                std::cout << "phie, max = " << Vmath::Vmax(nq, phie, 1) << ", min = " << Vmath::Vmin(nq, phie, 1) << std::endl;

                // Plotting the result
                int nvar = 8;
                int ncoeffs = m_fields[0]->GetNcoeffs();

                std::string outname;
                outname = m_sessionName + "_helm.chk";

                std::vector<Array<OneD, NekDouble>> fieldcoeffs(nvar);
                for (int i = 0; i < nvar; ++i)
                {
                    fieldcoeffs[i] = Array<OneD, NekDouble>(ncoeffs);
                }
                
                Array<OneD, Array<OneD, NekDouble>> mfmag(m_spacedim);
                Array<OneD, Array<OneD, NekDouble>> phiemfmag(m_spacedim);

                for (int j=0; j<m_spacedim; ++j)
                {
                    mfmag[j] = Array<OneD, NekDouble>(nq, 0.0);
                    phiemfmag[j] = Array<OneD, NekDouble>(nq, 0.0);

                    for (int k=0; k<m_spacedim; ++k)
                    {
                        Vmath::Vvtvp(nq, &m_movingframes[j][k*nq], 1, &m_movingframes[j][k*nq], 1, &mfmag[j][0], 1, &mfmag[j][0], 1);
                        Vmath::Vvtvp(nq, &m_phiemovingframes[j][k*nq], 1, &m_phiemovingframes[j][k*nq], 1, &phiemfmag[j][0], 1, &phiemfmag[j][0], 1);
                    }

                    Vmath::Vsqrt(nq, &mfmag[j][0], 1, &mfmag[j][0], 1);
                    Vmath::Vsqrt(nq, &phiemfmag[j][0], 1, &mfmag[j][0], 1);
                }


                std::vector<std::string> variables(nvar);
                variables[0] = "phim";
                variables[1] = "phimextra";
                variables[2] = "phie";
                variables[3] = "forcing";
                variables[4] = "AniX";
                variables[5] = "AniY";
                variables[6] = "phieAniX";
                variables[7] = "phieAniY";

                m_fields[0]->FwdTransLocalElmt(phim, fieldcoeffs[0]);
                m_fields[0]->FwdTransLocalElmt(phimextra, fieldcoeffs[1]);
                m_fields[0]->FwdTransLocalElmt(phie, fieldcoeffs[2]);
                m_fields[0]->FwdTransLocalElmt(forcing, fieldcoeffs[3]);

                m_fields[0]->FwdTransLocalElmt(m_AniStrength[0], fieldcoeffs[4]);
                m_fields[0]->FwdTransLocalElmt(m_AniStrength[1], fieldcoeffs[5]);
                m_fields[0]->FwdTransLocalElmt(m_phieAniStrength[0], fieldcoeffs[6]);
                m_fields[0]->FwdTransLocalElmt(m_phieAniStrength[1], fieldcoeffs[7]);

                WriteFld(outname, m_fields[0], fieldcoeffs, variables);

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
    // if(m_NeuralEPType==eNeuralHelmTest)
    // {
    //     Array<OneD, NekDouble> x0(nq);
    //     Array<OneD, NekDouble> x1(nq);
    //     Array<OneD, NekDouble> x2(nq);

    //     m_fields[0]->GetCoords(x0, x1, x2);

    //     int index;
    //     Array<OneD, NekDouble> phim(nq,0.0);
    //     Array<OneD, NekDouble> phieexact(nq,0.0);

    //     NekDouble coeffa = (1.0+1.0)/(1.0+1.0+1.0/m_ratio_re_ri+1.0/m_ratio_re_ri);
    //     for (int i = 0; i < m_InnerboxEnd; ++i)
    //     {
    //         for (int j = 0; j < m_fields[0]->GetTotPoints(i); ++j)
    //         {
    //             index = m_fields[0]->GetPhys_Offset(i) + j;
    //             phim[index] = cos(m_pi*x0[index])*cos(m_pi*x1[index]);
    //             phieexact[index] = -coeffa*cos(m_pi*x0[index])*cos(m_pi*x1[index]);
    //         }
    //     }

    //     // NekDouble coeffe;
    //     for (int i = m_InnerboxEnd; i < m_fields[0]->GetExpSize(); ++i)
    //     {
    //         for (int j = 0; j < m_fields[0]->GetTotPoints(i); ++j)
    //         {
    //             index = m_fields[0]->GetPhys_Offset(i) + j;

    //             // coeffe = coeffa;
    //             // if(x0[index]>1.0)
    //             // {
    //             //     coeffe = coeffa / (x0[index]-1) / (x0[index]-1);
    //             // }

    //             // else if (x0[index]<-1.0)
    //             // {
    //             //     coeffe = coeffa / (x0[index]+1) / (x0[index]+1);
    //             // }

    //             phieexact[index] = -coeffa*cos(m_pi*x0[index])*cos(m_pi*x1[index]);
    //         }
    //     }
    //     std::cout << "phim = " << RootMeanSquare(phim) << std::endl;

    //     StdRegions::ConstFactorMap phiefactors;
    //     phiefactors[StdRegions::eFactorTau]    = m_Helmtau;
    //     phiefactors[StdRegions::eFactorLambda] = 0.0;

    //     // // Compute \nabla \sigma_i \nabla phi_m and use it as point sources for
    //     // phi_e. This is equivalently achieved by removing all the point sources in
    //     // myelinnated fiber region.
    //     Array<OneD, NekDouble> phimLaplacian(nq);
    //     phimLaplacian = ComputeMMFDiffusion(m_movingframes, phim);
        
    //     NekDouble extFieldStr = 1.0;
    //     Vmath::Sadd(nq, -1.0 * AvgInt(phimLaplacian), phimLaplacian, 1, phimLaplacian, 1);
    //     Vmath::Smul(nq, -1.0 * extFieldStr, phimLaplacian, 1, m_fields[1]->UpdatePhys(), 1);

    //     m_fields[1]->HelmSolve(m_fields[1]->GetPhys(), m_fields[1]->UpdateCoeffs(), phiefactors, m_phievarcoeff);
    //     m_fields[1]->BwdTrans(m_fields[1]->GetCoeffs(), m_fields[1]->UpdatePhys());
    //     m_fields[1]->SetPhysState(true);

    //     Array<OneD, NekDouble> outarray(nq);
    //     outarray = m_fields[1]->GetPhys();

    //     Vmath::Sadd(nq, -1.0 * AvgInt(outarray), outarray, 1, outarray, 1);                    

    //     std::cout << "phie = " << RootMeanSquare(outarray) << std::endl;

    //     // Plotting

    //     int nvar    = 4;
    //     int ncoeffs = m_fields[0]->GetNcoeffs();

    //     std::string outname1;
    //     outname1 = m_sessionName + "_result.chk";

    //     std::vector<Array<OneD, NekDouble>> fieldcoeffs(nvar);
    //     for (int i = 0; i < nvar; ++i)
    //     {
    //         fieldcoeffs[i] = Array<OneD, NekDouble>(ncoeffs);
    //     }

    //     std::vector<std::string> variables(nvar);
    //     variables[0] = "phim";
    //     variables[1] = "phie";
    //     variables[2] = "phieexact";
    //     variables[3] = "phieerror";

    //     Array<OneD, NekDouble> phieInerr(nq,0.0);
    //     Array<OneD, NekDouble> phieOuterr(nq,0.0);
    //     Array<OneD, NekDouble> phieerr(nq,0.0);

    //     for (int i = 0; i < m_InnerboxEnd; ++i)
    //     {
    //         for (int j = 0; j < m_fields[0]->GetTotPoints(i); ++j)
    //         {
    //             index = m_fields[0]->GetPhys_Offset(i) + j;
    //             phieInerr[index] = outarray[index] - phieexact[index];
    //         }
    //     }

    //     for (int i = m_InnerboxEnd; i < m_fields[0]->GetExpSize(); ++i)
    //     {
    //         for (int j = 0; j < m_fields[0]->GetTotPoints(i); ++j)
    //         {
    //             index = m_fields[0]->GetPhys_Offset(i) + j;
    //             phieOuterr[index] = outarray[index] - phieexact[index];
    //         }
    //     }


    //     std::cout << "Error In = " << RootMeanSquare(phieInerr) << ", out = " << RootMeanSquare(phieOuterr) << std::endl;

    //     Vmath::Vadd(nq, phieInerr, 1, phieOuterr, 1, phieerr, 1);

    //     m_fields[0]->FwdTransLocalElmt(phim, fieldcoeffs[0]);
    //     m_fields[0]->FwdTransLocalElmt(outarray, fieldcoeffs[1]);
    //     m_fields[0]->FwdTransLocalElmt(phieexact, fieldcoeffs[2]);
    //     m_fields[0]->FwdTransLocalElmt(phieerr, fieldcoeffs[3]);

    //     WriteFld(outname1, m_fields[0], fieldcoeffs, variables);

    //     wait_on_enter();
    // }

}

MMFNeuralEP::~MMFNeuralEP()
{
}

void MMFNeuralEP::ComputeRegionalSigma(
    const Array<OneD, const int> &zoneindex,
    Array<OneD, Array<OneD, NekDouble>> &sigma_i,
    Array<OneD, Array<OneD, NekDouble>> &sigma_e,
    Array<OneD, Array<OneD, NekDouble>> &sigma_eM)
{
    int nq   = GetTotPoints();

        // Compute sigma_i
        // Node: 1.0, Myelin: m_Cn / m_Cm, Extraspace: 0.0
        int index;
        for (int i = 0; i<nq; ++i)
        {
            index = zoneindex[i];
            for (int j = 0; j < m_expdim; ++j)
            {
                // node for all fibers
                if( index>=0)
                {
                    sigma_i[j][i] = 1.0;
                    sigma_e[j][i] = 1.0/m_ratio_re_ri;
                    sigma_eM[j][i] = 1.0/m_ratio_re_ri;
                }

                // myelin for all fibers
                else if (index==-1)
                {
                    sigma_i[j][i] = m_AnisotropyStrength;
                    sigma_e[j][i] = m_AnisotropyStrength/m_ratio_re_ri;
                    sigma_eM[j][i] = m_AnisotropyStrength/m_ratio_re_ri;
                }

                else if (index==-2)
                {
                    sigma_i[j][i] = 0.0;
                    sigma_e[j][i] = 1.0/m_ratio_re_ri;
                    sigma_eM[j][i] = m_phiefactor/m_ratio_re_ri;
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

        std::cout << "Max sigma_e_1  = "
                    << Vmath::Vmax(nq, sigma_e[0], 1)
                    << ", sigma_e_2 = "
                    << Vmath::Vmax(nq, sigma_e[1], 1)
                    << ", Min sigma_e_1 = "
                    << Vmath::Vmin(nq, sigma_e[0], 1)
                    << ", sigma_e_2 = "
                    << Vmath::Vmin(nq, sigma_e[1], 1) << std::endl;
    }

        // Compute sigma_e
        // for (int i = 0; i<nq; ++i)
        // {
        //     index = zoneindex[i];
        //     // Node zone
        //     for (int j = 0; j < m_expdim; ++j)
        //     {
        //         if(index>=0)
        //         {
        //             sigma_e[j][i] = 1.0/m_ratio_re_ri;
        //         }

        //         // Myelin zone
        //         else if (index==-1)
        //         {
        //             // sigma_eM[j][i] = m_AnisotropyStrength/m_ratio_re_ri;
        //             sigma_e[j][i] = 1.0/m_ratio_re_ri;
        //         }

        //         else if (index==-2)
        //         {
        //             sigma_e[j][i] = 1.0/m_ratio_re_ri;
        //         }
        //     }
        // }

        // Compute sigma_e
        // for (int i = 0; i<nq; ++i)
        // {
        //     index = zoneindex[i];
        //     // Node zone
        //     for (int j = 0; j < m_expdim; ++j)
        //     {
        //         if(index>=0)
        //         {
        //             sigma_eM[j][i] = 1.0/m_ratio_re_ri;
        //         }

        //         // Myelin zone
        //         else if (index==-1)
        //         {
        //             sigma_eM[j][i] = m_AnisotropyStrength/m_ratio_re_ri;
        //         }

        //         else if (index==-2)
        //         {
        //             sigma_eM[j][i] = m_phiefactor/m_ratio_re_ri;
        //         }
        //     }
        // }

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

void MMFNeuralEP::IndexNodeZone2D(        
        const int numfiber, const int totNnode, 
        const NekDouble nodelen, const NekDouble myelinlen,
        const NekDouble nodeinitdown, const NekDouble nodeinitup,
        const Array<OneD, const NekDouble> &fiberleft,
        const Array<OneD, const NekDouble> &fiberright,
        const Array<OneD, const int> &fiberorder,
        Array<OneD, Array<OneD, int>> &outarray)
    {
        int nq   = GetTotPoints();

        // Compute the center of the element as the coordinate of each grid.
        Array<OneD, NekDouble> xcell(nq);
        Array<OneD, NekDouble> ycell(nq);
        Array<OneD, NekDouble> zcell(nq);

        Getcellavg(xcell,ycell,zcell);

        NekDouble xi, yi;
        for (int n=0; n<numfiber; ++n)
         {
            for (int i=0; i<nq; ++i)
            {
                xi = xcell[i];
                yi = ycell[i];

                outarray[n][i] = -2;
                if( (m_FiberType==eLinearAligned) || (m_FiberType==eLinearMisAligned) )
                {
                    if((xi>fiberleft[n]) && (xi<fiberright[n]))
                        {
                            outarray[n][i] = FiberIndex(m_FiberType, n, totNnode, nodelen, myelinlen, 
                                                            nodeinitdown, nodeinitup, fiberorder[n], xi, yi);
                        }
                }

                else{
                    outarray[n][i] = FiberIndex(m_FiberType, n, totNnode, nodelen, myelinlen, 
                                                            nodeinitdown, nodeinitup, fiberorder[n], xi, yi);
                }
            }
        }
    }

Array<OneD, int> MMFNeuralEP::TestRanvierSingleIndex()
{
    int nq   = GetTotPoints();

    // Compute the center of the element as the coordinate of each grid.
    Array<OneD, NekDouble> xcell(nq);
    Array<OneD, NekDouble> ycell(nq);
    Array<OneD, NekDouble> zcell(nq);

    Getcellavg(xcell,ycell,zcell);

    Array<OneD, int> outarray(nq, -1);

    // NekDouble nodestart, nodeend;
    NekDouble fiberleft = -0.2;
    NekDouble fiberright = 0.2;

    for (int i=0; i<nq; ++i)
    {
        if((xcell[i]>fiberleft) && (xcell[i]<fiberright))
        {
            if( (ycell[i]>fiberleft) && (ycell[i]<fiberright) )
            {
                outarray[i] = 0;
            }

            else
            {
                outarray[i] = -1;
            }
        }

        else
        {
            outarray[i] = -2;
        }
    }

    return outarray;
}


Array<OneD, int> MMFNeuralEP::TestRanvierDuoIndex()
{
    int nq   = GetTotPoints();

    // Compute the center of the element as the coordinate of each grid.
    Array<OneD, NekDouble> xcell(nq);
    Array<OneD, NekDouble> ycell(nq);
    Array<OneD, NekDouble> zcell(nq);

    Getcellavg(xcell,ycell,zcell);

    Array<OneD, int> outarray(nq, -1);

    // NekDouble nodestart, nodeend;
    NekDouble fiberx1 = -0.5;
    NekDouble fiberx2 = -0.1;
    NekDouble fiberx3 = 0.1;
    NekDouble fiberx4 = 0.5;

    NekDouble fibery1 = -0.2;
    NekDouble fibery2 = 0.2;

    for (int i=0; i<nq; ++i)
    {
        if((xcell[i]>fiberx1) && (xcell[i]<fiberx2))
        {
            if( (ycell[i]>fibery1) && (ycell[i]<fibery2) )
            {
                outarray[i] = 0;
            }

            else
            {
                outarray[i] = -1;
            }
        }

        else if((xcell[i]>fiberx3) && (xcell[i]<fiberx4))
        {
            if( (ycell[i]>fibery1) && (ycell[i]<fibery2) )
            {
                outarray[i] = 0;
            }

            else
            {
                outarray[i] = -1;
            }
        }

        else
        {
            outarray[i] = -2;
        }
    }

    return outarray;
}


int MMFNeuralEP::FiberIndex(FiberType FiberType, const int fibern,
    const int totNnode, const NekDouble nodelen, const NekDouble myelinlen,
    const NekDouble nodeinitdown, const NekDouble nodeinitup, 
    const int fiberorder, const NekDouble xi, const NekDouble yi)
    {
        int index=0;

        switch(FiberType)
        {
            case eLinearAligned:
            {
                index = LinearAlignedFiberIndex(totNnode, nodelen, myelinlen, nodeinitdown, nodeinitup, fiberorder, yi);                
                break;
            }

            case eLinearMisAligned:
            {
                index = LinearMisAlignedFiberIndex(fibern, totNnode, nodelen, myelinlen, nodeinitdown, nodeinitup, fiberorder, yi);                
                break;
            }

            case eLinearDivergent:
            {
                index = LinearDivergentFiberIndex(fibern, totNnode, nodelen, myelinlen, nodeinitdown, nodeinitup, fiberorder, xi, yi);                
                break;
            }

            case eConstantCurved:
            {
                index = ConstantCurvedFiberIndex(fibern, m_fibercurvature, totNnode, nodelen, myelinlen, xi, yi);                
                break;
            }

            default:
            break;
        }

        return index;
    }


int MMFNeuralEP::ConstantCurvedFiberIndex(const int fibern, 
    const NekDouble fibercurvature, const int totNnode, const NekDouble nodelen, 
    const NekDouble myelinlen, const NekDouble xi, const NekDouble yi)
{
    int output = -2;
    
    NekDouble angle0, nodetheta, myelintheta;
    NekDouble rad, theta, thetatop, thetabottom;
    NekDouble radbottom, radtop;
    NekDouble fiberstart, fiberend;

    angle0 = m_pi/72.0;

    nodetheta = nodelen/fibercurvature;
    myelintheta = myelinlen/fibercurvature;

    fiberstart = angle0 + nodetheta;
    fiberend = angle0 + 3*nodetheta + (totNnode-1)*(nodetheta+myelintheta);

    rad = sqrt(xi*xi + yi*yi);
    theta = atan2(yi/rad,xi/rad);

    if(theta<0)
    {
        theta = theta + 2*m_pi;
    }

    radbottom = fibercurvature + (2*fibern+1)*0.01;
    radtop = radbottom + 0.01;

    if( (rad>radbottom) && (rad<radtop) )
    {
        if ( (theta > fiberstart ) && (theta < fiberend)  )
        {
            output = -1;

            for (int k=0; k<totNnode; ++k)
            {
                thetabottom = angle0 + 2*nodetheta + k*(nodetheta+myelintheta);
                thetatop = thetabottom + nodetheta;

                if ( (theta > thetabottom) && (theta < thetatop) )
                {
                    output = k;
                }
            }
        }
    }
    
    return output;
}


int MMFNeuralEP::LinearAlignedFiberIndex(
    const int totNnode, const NekDouble nodelen, const NekDouble myelinlen,
    const NekDouble nodeinitdown, const NekDouble nodeinitup, 
    const int fiberorder, const NekDouble yi)
{
    int output = -2;
    
    NekDouble  nodestart, nodeend;

    output = -1 ;  // Default of the first fiber = myelin

    // Excitezone
    if( (yi>=nodeinitdown) && (yi<=nodeinitup) )
    {
        if(fiberorder==-1)
        {
            output = totNnode;
        }

        else{
            output = 0;
        }
    }

    for (int k=0; k<totNnode; ++k)
    {
        nodestart = nodeinitup + myelinlen + k * (myelinlen + nodelen);
        nodeend = nodeinitup + (k+1) * (myelinlen + nodelen);
        if( (yi>=nodestart) && (yi<=nodeend) )
        {
            if(fiberorder==-1)
            {
               output = totNnode - 1 - k;
            }

            else{
                output = k + 1;
            }
        }
    }

    if(yi<nodeinitdown)
    {
        output = -2;
    }

    nodeend = nodeinitup + totNnode * (myelinlen + nodelen);
    if(yi>nodeend)
    {
        output = -2;
    }
    
    return output;
}


int MMFNeuralEP::LinearMisAlignedFiberIndex(const int fibern,
    const int totNnode, const NekDouble nodelen, const NekDouble myelinlen,
    const NekDouble nodeinitdown, const NekDouble nodeinitup, 
    const int fiberorder, const NekDouble yi)
{
    boost::ignore_unused(fiberorder);

    int output = -2;
    
    NekDouble  nodestart, nodeend;

    output = -1 ;  // Default of the first fiber = myelin

    // Excitezone
    if( (yi>=nodeinitdown) && (yi<=nodeinitup) )
    {
        output = 0;
    }

    for (int k=0; k<totNnode; ++k)
    {
        if(fibern==0)
        {
            nodestart = nodeinitup + 2*myelinlen + nodelen + k * (2*myelinlen + 2*nodelen);
            nodeend = nodeinitup + (k+1) * (2*myelinlen + 2*nodelen);

            if( (yi>=nodestart) && (yi<=nodeend) )
            {
                output = k + 2;
            }
        }

        else if(fibern==1)
        {
            nodestart = nodeinitup + myelinlen + k * (2*myelinlen + 2*nodelen);
            nodeend = nodeinitup + myelinlen + nodelen + k * (2*myelinlen + 2*nodelen);

           if( (yi>=nodestart) && (yi<=nodeend) )
            {
                output = k + 1;
            }
        }
    }

    if(fibern==0)
    {
        nodestart = nodeinitup + myelinlen ;
        nodeend = nodeinitup + myelinlen + nodelen ;

        if( (yi>=nodestart) && (yi<=nodeend) )
        {
            output = 1;
        }
    }

    if(fibern==1)
    {
        nodestart = nodeinitup + 2*myelinlen + nodelen + (totNnode-1) * (2*myelinlen + 2*nodelen);
        nodeend = nodeinitup + totNnode * (2*myelinlen + 2*nodelen);

        if( (yi>=nodestart) && (yi<=nodeend) )
        {
            output = totNnode + 1;
        }
    }

    if(yi<nodeinitdown)
    {
        output = -2;
    }

    nodeend = nodeinitup + totNnode * (2*myelinlen + 2*nodelen);
    if(yi>nodeend)
    {
        output = -2;
    }
    
    return output;
}

int MMFNeuralEP::LinearDivergentFiberIndex(const int fibern,
    const int totNnode, const NekDouble nodelen, const NekDouble myelinlen,
    const NekDouble nodeinitdown, const NekDouble nodeinitup, 
    const int fiberorder, const NekDouble xi, const NekDouble yi)
{
    boost::ignore_unused(fiberorder);

    int output = -2;
    
    NekDouble  nodestart, nodeend;

    if((xi>m_fiberleft[fibern]) && (xi<m_fiberright[fibern]))
    {
        // Excitezone
        if( (yi>=nodeinitdown) && (yi<=nodeinitup) )
        {
            output = 0;
        }

        for (int k=0; k<2; ++k)
        {

                nodestart = nodeinitup + myelinlen + k * (myelinlen + nodelen);
                nodeend = nodeinitup + (k+1) * (myelinlen + nodelen);
                if( (yi>=nodestart) && (yi<=nodeend) )
                {
                    output = k + 1;
                }
        }
    }

    NekDouble uppeval, loweval;
    NekDouble x1 = 0.01;
    NekDouble x2 = 0.02;
    NekDouble x6 = 0.05;
    NekDouble x7 = 0.06;
    NekDouble y1 = 0.44;

    // First fiber
    if(fibern==0)
    {
        if(yi<=y1)
        {
            if((xi>m_fiberleft[fibern]) && (xi<m_fiberright[fibern]))
            {
                output = -1;
                for (int k=0; k<2; ++k)
                {
                    nodestart = nodeinitup + myelinlen + k * (myelinlen + nodelen);
                    nodeend = nodeinitup + (k+1) * (myelinlen + nodelen);
                    if( (yi>=nodestart) && (yi<=nodeend) )
                    {
                        output = k + 1;
                    }
                }

                if( (yi>=nodeinitdown) && (yi<=nodeinitup) )
                {
                    output = 0;
                }  
            }
        }

        else{
            // upper line
            uppeval = -tan(m_fiberangle)*xi + (y1 + tan(m_fiberangle)*x2);
            loweval = -tan(m_fiberangle)*xi + (y1 + tan(m_fiberangle)*x1);

           if( (yi<uppeval) && (yi>loweval) )
            {
                std::cout << "m_fiberangle = " << m_fiberangle <<", tan(m_fiberangle) = " <<  tan(m_fiberangle) 
                << ", fibern = " << fibern << ", xi = " << xi << ", yi = " << yi 
                << ", uppeval = " << uppeval << ", loweval = " << loweval << std::endl;

                output = -1;

                for (int k=0; k<2; ++k)
                {
                    nodestart = nodeinitup + 2 * (myelinlen + nodelen) - nodelen +(2-k) * (myelinlen + nodelen)*sin(m_fiberangle);
                    nodeend = nodeinitup + 2 * (myelinlen + nodelen) + ((2-k) * (myelinlen + nodelen))*sin(m_fiberangle);

                    if( yi<nodestart)
                    {
                        output = -1;
                    }

                    else if( (yi>=nodestart) && (yi<=nodeend) )
                    {
                        output = totNnode - k - 1;
                    }
                }
             }
        }
    }

    // Second fiber
    if(fibern==1)
    {
        if((xi>m_fiberleft[fibern]) && (xi<m_fiberright[fibern]))
        {
            output = -1;

            for (int k=0; k<totNnode; ++k)
            {
                nodestart = nodeinitup + myelinlen + k * (myelinlen + nodelen);
                nodeend = nodeinitup + (k+1) * (myelinlen + nodelen);
                if( (yi>=nodestart) && (yi<=nodeend) )
                {
                    output = k + 1;
                }
            }

            if( (yi>=nodeinitdown) && (yi<=nodeinitup) )
            {
                output = 0;
            }  
        }
    }

    if(fibern==2)
    {
        if(yi<=y1)
        {
            if((xi>m_fiberleft[fibern]) && (xi<m_fiberright[fibern]))
            {
                output = -1;
                for (int k=0; k<2; ++k)
                {
                    nodestart = nodeinitup + myelinlen + k * (myelinlen + nodelen);
                    nodeend = nodeinitup + (k+1) * (myelinlen + nodelen);
                    if( (yi>=nodestart) && (yi<=nodeend) )
                    {
                        output = k + 1;
                    }
                }

                if( (yi>=nodeinitdown) && (yi<=nodeinitup) )
                {
                    output = 0;
                }  
            }
        }

        else{
            uppeval =  tan(m_fiberangle)*xi + (y1 - tan(m_fiberangle)*x6);
            loweval =  tan(m_fiberangle)*xi + (y1 - tan(m_fiberangle)*x7);

           if( (yi<uppeval) && (yi>loweval) )
            {
                std::cout << "m_fiberangle = " << m_fiberangle <<", tan(m_fiberangle) = " <<  tan(m_fiberangle) 
                << ", fibern = " << fibern << ", xi = " << xi << ", yi = " << yi 
                << ", uppeval = " << uppeval << ", loweval = " << loweval << std::endl;

                output = -1;

                for (int k=0; k<2; ++k)
                {
                    nodestart = nodeinitup + 2 * (myelinlen + nodelen) - nodelen +(2-k) * (myelinlen + nodelen)*sin(m_fiberangle);
                    nodeend = nodeinitup + 2 * (myelinlen + nodelen) + ((2-k) * (myelinlen + nodelen))*sin(m_fiberangle);

                    if( yi<nodestart)
                    {
                        output = -1;
                    }

                    else if( (yi>=nodestart) && (yi<=nodeend) )
                    {
                        output = totNnode - k - 1;
                    }
                }
             }
        }
    }

    if(yi<nodeinitdown)
    {
        output = -2;
    }

    nodeend = nodeinitup + (totNnode-1) * (myelinlen + nodelen);
    if(yi>nodeend)
    {
        output = -2;
    }
    
    return output;
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


// void MMFNeuralEP::SetUpBiAnisotropy(
//             const Array<OneD, const int> &zoneindex,
//             const Array<OneD, const Array<OneD, NekDouble>> NeuralCm,
//             Array<OneD, Array<OneD, NekDouble>> &AniStrength)
// {
//     int nq   = GetTotPoints();
//     int nvar = m_fields.size();

//     for (int j = 0; j < m_expdim; ++j)
//     {
//         AniStrength[j] = Array<OneD, NekDouble>(nq, 1.0);
//     }

//     if (m_MediumType == eAnisotropy)
//     {
//         for (int j = 0; j < m_expdim; ++j)
//         {
//             Vmath::Smul(nq, m_Cn, &NeuralCm[0][0], 1, &AniStrength[j][0], 1);
//         }
//     }

//     // Let the lenght of moving frames outside the fiber to be zero.
//     int cnt = 0;
//     for (int i=0; i<nq; ++i)
//     {
//         if(zoneindex[i] == -2)
//         {
//             AniStrength[0][i] = 0.0;
//             AniStrength[1][i] = 0.0; 
//             cnt++;
//         }
//     }

//     std::cout << " =======================================================================" << std::endl;
//     std::cout << " Moving frames " << cnt << " / " << nq << " ( " << 100.0*cnt/nq << " % ) are removed" << std::endl;
//     std::cout << " =======================================================================" << std::endl;


//     if(nvar==2)
//     {
//         std::cout << "Max Anistrength_1  = "
//                     << Vmath::Vmax(nq, AniStrength[0], 1)
//                     << ", Anistrength_2 = "
//                     << Vmath::Vmax(nq, AniStrength[1], 1)
//                     << ", Min Anistrength 1 = "
//                     << Vmath::Vmin(nq, AniStrength[0], 1)
//                     << ", Anistrength 2 = "
//                     << Vmath::Vmin(nq, AniStrength[1], 1) << std::endl;
//     }

//     else{
//         std::cout << "Max Anistrength_1  = "
//             << Vmath::Vmax(nq, AniStrength[0], 1)
//             << ", Min Anistrength 1 = "
//             << Vmath::Vmin(nq, AniStrength[0], 1) << std::endl;
//     }
// }

// void MMFNeuralEP::SetUpBiAnisotropy(
//             const Array<OneD, const Array<OneD, int>> &zoneindex,
//             const Array<OneD, const NekDouble> NeuralCm,
//             Array<OneD, Array<OneD, NekDouble>> &AniStrength)
// {
//     int nq   = GetTotPoints();
//     int nvar = m_fields.size();

//     for (int j = 0; j < m_expdim; ++j)
//     {
//         AniStrength[j] = Array<OneD, NekDouble>(nq, 1.0);
//     }

//     if (m_MediumType == eAnisotropy)
//     {
//         for (int j = 0; j < m_expdim; ++j)
//         {
//             Vmath::Smul(nq, m_Cn, &NeuralCm[0], 1, &AniStrength[j][0], 1);
//         }
//     }

//     // Let the lenght of moving frames outside the fiber to be zero.
//     // int index, cnt = 0;
//     // for (int i=0; i<nq; ++i)
//     // {
//     //     for (int n=0; n<nq; ++n)
//     //     {
//     //         index = zoneindex[n][i];
//     //         if( index == -2)
//     //         {
//     //             AniStrength[0][i] = 0.0;
//     //             AniStrength[1][i] = 0.0; 
//     //             cnt++;
//     //         }
//     //     }
//     // }

//     // std::cout << " =======================================================================" << std::endl;
//     // std::cout << " Moving frames " << cnt << " / " << nq << " ( " << 100.0*cnt/nq << " % ) are removed" << std::endl;
//     // std::cout << " =======================================================================" << std::endl;


//     if(nvar==2)
//     {
//         std::cout << "Max Anistrength_1  = "
//                     << Vmath::Vmax(nq, AniStrength[0], 1)
//                     << ", Anistrength_2 = "
//                     << Vmath::Vmax(nq, AniStrength[1], 1)
//                     << ", Min Anistrength 1 = "
//                     << Vmath::Vmin(nq, AniStrength[0], 1)
//                     << ", Anistrength 2 = "
//                     << Vmath::Vmin(nq, AniStrength[1], 1) << std::endl;
//     }

//     else{
//         std::cout << "Max Anistrength_1  = "
//             << Vmath::Vmax(nq, AniStrength[0], 1)
//             << ", Min Anistrength 1 = "
//             << Vmath::Vmin(nq, AniStrength[0], 1) << std::endl;
//     }f
// }

void MMFNeuralEP::SetUpDomainZone(
        const Array<OneD, const Array<OneD, int>> &zoneindexfiber,
        Array<OneD, Array<OneD, NekDouble>> &excitezonefiber,
        Array<OneD, Array<OneD, NekDouble>> &intrazonefiber,
        Array<OneD, NekDouble> &nodezone,
        Array<OneD, NekDouble> &myelinzone)
{
    int nq   = GetTotPoints();
    int index;

    // Set up the total node zone and intra zone
    nodezone = Array<OneD, NekDouble>(nq, 0.0);  
    myelinzone = Array<OneD, NekDouble>(nq, 0.0);              
    for (int i=0; i<nq; ++i)
    {
        for (int n=0; n<m_numfiber; ++n)
        {
            index = zoneindexfiber[n][i];

            // excite zone
            if(index == 0)
            {
                excitezonefiber[n][i] = 1.0;
            }

            // intra zone = myelin or nodezone 
            if( index > -2) 
            {
                intrazonefiber[n][i] = 1.0;
            }

            // node zone
            if (index > -1)
            {
                 nodezone[i] = 1.0;
            }

            // myelin zone
            if (index == -1)
            {
                 myelinzone[i] = 1.0;
            }
        }
    }
}

Array<OneD, NekDouble> MMFNeuralEP::ComputeConductivity(
                 const Array<OneD, const int> &zoneindex)
{
    int nq   = GetTotPoints();

    Array<OneD, NekDouble> outarray(nq);

    int cntn = 0, cntm = 0, cnte = 0;
    for (int i = 0; i < nq; ++i)
    {
        // Ranvier node zone
        if ( (zoneindex[i] >= 0) && (zoneindex[i] < 100) )

        {
            outarray[i] = 1.0 / m_Cn;
            cntn++;
        }

        if ( (zoneindex[i] >= 100) && (zoneindex[i] < 200) )

        {
            outarray[i] = 1.0 / m_Cn;
            cntn++;
        }

        // Myelin zone
        if (zoneindex[i] == -1) 
        {
            outarray[i] = 1.0 / m_Cm;
            cntm++;
        }

        // Myelin zone
        if (zoneindex[i] == -101) 
        {
            outarray[i] = 1.0 / m_Cm;
            cntm++;
        }

        // Extracellular space: \sigma_i = m_ratio_re_ri * \sigma_e
        else if (zoneindex[i] == -2)
        {
            outarray[i] = 0.0;
            cnte++;
        }
    }

    cnte = nq - cntn - cntm;

    std::cout << "ComputeConductivity: Node = " << cntn/m_npts << ", Myelinf1 = " 
    << cntm/m_npts << ", extracell = " << cnte/m_npts << std::endl;

    return outarray;
}

Array<OneD, NekDouble> MMFNeuralEP::ComputeConductivity(
                 const Array<OneD, const Array<OneD, int>> &zoneindex)
{
    int nq   = GetTotPoints();

    Array<OneD, NekDouble> outarray(nq, 0.0); // NeuralCm is zero at indexzone == -2

    int index, cntm = 0, cntn = 0, cnte = 0;

    for (int i = 0; i < nq; ++i)
    {
        for (int n=0; n<m_numfiber; ++n)
            {
                index = zoneindex[n][i];

                // Ranvier node zone for all the fibers
                if ( index >= 0)
                {
                    outarray[i] = 1.0 / m_Cn;
                    cntn++;
                }

                // Ranvier node zone for all the fibers
                if (index == -1) 
                {
                    if (m_MediumType == eAnisotropy)
                    {  
                        outarray[i] = 1.0 / m_Cm;
                    }

                    else
                    {
                        outarray[i] = 1.0 / m_Cn;
                    }
                    cntm++;
                }
            }
    }

    cnte = nq - cntn - cntm;

    std::cout << "ComputeConductivity_New: Node = " << cntn/m_npts << ", Myelinf1 = " 
    << cntm/m_npts << ", extracell = " << cnte/m_npts << std::endl;

    return outarray;
}

void MMFNeuralEP::CheckNodeZoneMF(
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    const Array<OneD, const Array<OneD, NekDouble>> &phiemovingframes)
{
    int nq = GetTotPoints();

    Array<OneD, NekDouble> x0(nq);
    Array<OneD, NekDouble> x1(nq);
    Array<OneD, NekDouble> x2(nq);

    m_fields[0]->GetCoords(x0, x1, x2);

    int i, j, index=0, npts;;
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

        e1mag = sqrt(e1mag / npts);
        e2mag = sqrt(e2mag / npts);

        phie1mag = sqrt(phie1mag / npts);
        phie2mag = sqrt(phie2mag / npts);

        xp = (xp / npts);
        yp = (yp / npts);
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
            DoSolveMMF();
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

void MMFNeuralEP::DoSolveMMF()
{
    ASSERTL0(m_intScheme != 0, "No time integration scheme.");

    // int i, nchk = 1;
    const int nq               = GetTotPoints();
    const int nfields          = m_fields.size();
    const int totsteps = (m_steps + 1) / m_checksteps;

    int step = 0, nchk = 1;

    NekDouble intTime = 0.0, cpuTime = 0.0;

    const int nvariables = m_intVariables.empty() ? nfields : m_intVariables.size();
    if (m_intVariables.empty())
    {
        for (int i = 0; i < nfields; ++i)
        {
            m_intVariables.push_back(i);
        }
    }

    // Set up working arrays
    Array<OneD, Array<OneD, NekDouble>> fields(nvariables), fields_old(nvariables);
    Array<OneD, Array<OneD, NekDouble>> dphidt(nvariables), dphidtint(nvariables);
    Array<OneD, Array<OneD, NekDouble>> TimeMap(nvariables);

    for (int i = 0; i < nvariables; ++i)
    {
        fields[i]     = m_fields[m_intVariables[i]]->GetPhys();
        fields_old[i] = Array<OneD, NekDouble>(nq, 0.0);
        dphidt[i]     = Array<OneD, NekDouble>(nq, 0.0);
        dphidtint[i]  = Array<OneD, NekDouble>(nq, 0.0);
        TimeMap[i]    = Array<OneD, NekDouble>(nq, 0.0);

        m_fields[m_intVariables[i]]->SetPhysState(false);
    }

    m_TimeMap = TimeMap;  // Save reference for external access
    m_intScheme->InitializeScheme(m_timestep, fields, m_time, m_ode);

   // Prepare diagnostics arrays
    Array<OneD, NekDouble> timevec(totsteps, 0.0);
    Array<OneD, NekDouble> thredlocf1(totsteps, 0.0), thredlocf2(totsteps, 0.0);
    Array<OneD, int> thredlocf1zone(totsteps, 0), thredlocf2zone(totsteps, 0);
    Array<OneD, NekDouble> phimvecatnode1(totsteps, 0.0), phievecatnode1(totsteps, 0.0), CSDvecatnode1(totsteps, 0.0), rhovecatnode1(totsteps, 0.0);
    Array<OneD, NekDouble> phimvecatnode2(totsteps, 0.0), phievecatnode2(totsteps, 0.0), CSDvecatnode2(totsteps, 0.0), rhovecatnode2(totsteps, 0.0);
    Array<OneD, NekDouble> phimvecatnode3(totsteps, 0.0), phievecatnode3(totsteps, 0.0), CSDvecatnode3(totsteps, 0.0), rhovecatnode3(totsteps, 0.0);
    Array<OneD, NekDouble> phimvecatnode4(totsteps, 0.0), phievecatnode4(totsteps, 0.0), CSDvecatnode4(totsteps, 0.0), rhovecatnode4(totsteps, 0.0);
    Array<OneD, NekDouble> phimvecatnode5(totsteps, 0.0), phievecatnode5(totsteps, 0.0), CSDvecatnode5(totsteps, 0.0), rhovecatnode5(totsteps, 0.0);

    LibUtilities::Timer timer;
    Array<OneD, NekDouble> CSD;
    while (step < m_steps || m_time < m_fintime - NekConstants::kNekZeroTol)
    {
        // Save current solution
        for (int n=0; n<nvariables; ++n)
        {
            Vmath::Vcopy(nq, &fields[n][0], 1, &fields_old[n][0], 1);
        }

        // Time integration
        timer.Start();
        fields = m_intScheme->TimeIntegrate(step, m_timestep, m_ode);
        timer.Stop();

        m_time += m_timestep;
        NekDouble elapsed = timer.TimePerTest(1);
        intTime += elapsed;
        cpuTime += elapsed;

        // Compute normalized time derivatives
        fields[1] = m_fields[1]->GetPhys();
        for (int n = 0; n < std::min(nvariables, 3); ++n)
        {
            NekDouble maxphi = Vmath::Vamax(nq, fields[n], 1);
            Vmath::Vsub(nq, fields[n], 1, fields_old[n], 1, dphidt[n], 1);
            Vmath::Smul(nq, 1.0 / (m_timestep * maxphi), dphidt[n], 1, dphidt[n], 1);
        }

       if (m_time >= m_TimeMapStart && m_time <= m_TimeMapEnd)
        {
            ComputeNeuralTimeMap(m_time, m_zoneindexfiber, fields[0], dphidt[0], dphidtint[0], m_TimeMap[0]);

            CSD = ComputeMMFDiffusion(m_movingframes, fields[1]);
            ComputephieTimeMap(m_time, fields[1], dphidtint[1], m_TimeMap[1]);
            ComputerhoTimeMap(m_time, fields[2], dphidtint[2], m_TimeMap[2]);
        }

        // Info output
        if ((step + 1) % m_infosteps == 0 && m_session->GetComm()->GetRank() == 0)
        {
            std::cout << "Steps: " << std::setw(8) << std::left << step + 1
                      << " Time: " << std::setw(12) << m_time
                      << ", CPU Time = " << cpuTime / 60.0 << " min.\n\n";
            cpuTime = 0.0;
        }

        // Write out checkpoint files
        if ((m_checksteps && step && !((step + 1) % m_checksteps)))
        {
            PlotNeuralEP(fields, m_TimeMap, nchk);
            timevec[nchk] = m_time; 

            Array<OneD, NekDouble> CSD = ComputeMMFDiffusion(m_phiediffmovingframes, fields[1]);
            Vmath::Smul(nq, -m_Diffext, CSD, 1, CSD, 1);

            phimvecatnode1[nchk] = DisplayAtNodes(0, 1, m_zoneindexfiber, fields[0]);
            phievecatnode1[nchk] = DisplayAtNodes(0, 1, m_zoneindexfiber, fields[1]);
            CSDvecatnode1[nchk] = DisplayAtNodes(0, 1, m_zoneindexfiber, CSD);
            rhovecatnode1[nchk] = DisplayAtNodes(0, 1, m_zoneindexfiber, fields[2]);

            phimvecatnode2[nchk] = DisplayAtNodes(0, 2, m_zoneindexfiber, fields[0]);
            phievecatnode2[nchk] = DisplayAtNodes(0, 2, m_zoneindexfiber, fields[1]);
            CSDvecatnode2[nchk] = DisplayAtNodes(0, 2, m_zoneindexfiber, CSD);
            rhovecatnode2[nchk] = DisplayAtNodes(0, 2, m_zoneindexfiber, fields[2]);

            phimvecatnode3[nchk] = DisplayAtNodes(0, 3, m_zoneindexfiber, fields[0]);
            phievecatnode3[nchk] = DisplayAtNodes(0, 3, m_zoneindexfiber, fields[1]);
            CSDvecatnode3[nchk] = DisplayAtNodes(0, 3, m_zoneindexfiber, CSD);
            rhovecatnode3[nchk] = DisplayAtNodes(0, 3, m_zoneindexfiber, fields[2]);

            phimvecatnode4[nchk] = DisplayAtNodes(0, 4, m_zoneindexfiber, fields[0]);
            phievecatnode4[nchk] = DisplayAtNodes(0, 4, m_zoneindexfiber, fields[1]);
            CSDvecatnode4[nchk] = DisplayAtNodes(0, 4, m_zoneindexfiber, CSD);
            rhovecatnode4[nchk] = DisplayAtNodes(0, 4, m_zoneindexfiber, fields[2]);

            phimvecatnode5[nchk] = DisplayAtNodes(0, 5, m_zoneindexfiber, fields[0]);
            phievecatnode5[nchk] = DisplayAtNodes(0, 5, m_zoneindexfiber, fields[1]);
            CSDvecatnode5[nchk] = DisplayAtNodes(0, 5, m_zoneindexfiber, CSD);
            rhovecatnode5[nchk] = DisplayAtNodes(0, 5, m_zoneindexfiber, fields[2]);

            std::cout << "phim at node 1 = " << phimvecatnode1[nchk]  << ", at node 2 = " << phimvecatnode2[nchk] << ", at node 3 = " << phimvecatnode3[nchk] 
           << ", at node 4 = " << phimvecatnode4[nchk] << ", at node 5 = " << phimvecatnode5[nchk] << std::endl;

            std::cout << "phie at node 1 = " << phievecatnode1[nchk]  << ", at node 2 = " << phievecatnode2[nchk] << ", at node 3 = " << phievecatnode3[nchk] 
           << ", at node 4 = " << phievecatnode4[nchk] << ", at node 5 = " << phievecatnode5[nchk] << std::endl;

           std::cout << "CSD at node 1 = " << CSDvecatnode1[nchk]  << ", at node 2 = " << CSDvecatnode2[nchk] << ", at node 3 = " << CSDvecatnode3[nchk] 
           << ", at node 4 = " << CSDvecatnode4[nchk] << ", at node 5 = " << CSDvecatnode5[nchk] << std::endl;

           std::cout << "rho at node 1 = " << rhovecatnode1[nchk]  << ", at node 2 = " << rhovecatnode2[nchk] << ", at node 3 = " << rhovecatnode3[nchk] 
           << ", at node 4 = " << rhovecatnode4[nchk] << ", at node 5 = " << rhovecatnode5[nchk] << std::endl;

            // if (m_numfiber >= 1)
            // {
            //     thredlocf1zone[nchk] = FiberIndex(m_FiberType, 0, m_totNode, m_nodelen, m_myelinlen,
            //                                       m_nodeinitdown, m_nodeinitup, m_fiberorder[0], 0.0, thredlocf1[nchk]);
            // }

            // if (m_numfiber >= 2)
            // {
            //     thredlocf2zone[nchk] = FiberIndex(m_FiberType, 1, m_totNode, m_nodelen, m_myelinlen,
            //                                       m_nodeinitdown, m_nodeinitup, m_fiberorder[1], 0.0, thredlocf2[nchk]);
            // }

            Checkpoint_Output(nchk++);
        }

        ++step;
    } // namespace Nektar

    // Print out summary statistics
    if (m_session->GetComm()->GetRank() == 0)
    {
        std::cout << "Time-integration complete. Total CPU time: " << intTime << "s\n";
    }

    std::cout << " timevec: ";
    for (int i=0; i<totsteps; ++i)
    {
        std::cout << timevec[i] << ", ";
    }
    std::cout << std::endl;

    std::cout << " CSDvecatnode1: ";
    for (int i=0; i<totsteps; ++i)
    {
        std::cout << CSDvecatnode1[i] << ", ";
    }
    std::cout << std::endl;

    std::cout << " rhovecatnode1: ";
    for (int i=0; i<totsteps; ++i)
    {
        std::cout << rhovecatnode1[i] << ", ";
    }
    std::cout << std::endl;

    // std::cout << " thredlocf1: ";
    // for (int i=0; i<totsteps; ++i)
    // {
    //     std::cout << thredlocf1[i] << ", ";
    // }
    // std::cout << std::endl;

    // std::cout << " thredlocf1zone: ";
    // for (int i=0; i<totsteps; ++i)
    // {
    //     std::cout << thredlocf1zone[i] << ", ";
    // }
    // std::cout << std::endl;

    // if(m_numfiber>1)
    // {
    //     std::cout <<  "thredlocf2: ";
    //     for (int i=0; i<totsteps; ++i)
    //     {
    //         std::cout << thredlocf2[i] << ", ";
    //     }

    //     std::cout << std::endl;

    //     std::cout <<  "thredlocf2zone: ";
    //     for (int i=0; i<totsteps; ++i)
    //     {
    //         std::cout << thredlocf2zone[i] << ", ";
    //     }

    //     std::cout << std::endl;
    // }

    for (int i = 0; i < nvariables; ++i)
    {
        m_fields[m_intVariables[i]]->SetPhys(fields[i]);
        m_fields[m_intVariables[i]]->SetPhysState(true);

        m_fields[m_intVariables[i]]->FwdTrans(m_fields[i]->GetPhys(),
                                   m_fields[m_intVariables[i]]->UpdateCoeffs());
    }
} 
// namespace Nektar
NekDouble MMFNeuralEP::DisplayAtNodes(const int fibern, const int nodeindex, 
        const Array<OneD, const Array<OneD, int>> &zoneindexfiber,
        const Array<OneD, const NekDouble> &inarray)
{
    int nq = GetTotPoints();

    NekDouble output=0.0;

    #pragma omp parallel for
    for (int i = 0; i < nq; ++i)
    {
        if( zoneindexfiber[fibern][i] == nodeindex )
        {
            output += inarray[i];
        }
    }

    return output;
}

void MMFNeuralEP::ComputeNeuralTimeMap(const NekDouble time,
                                    const Array<OneD, const Array<OneD, int>> &zoneindex,
                                    const Array<OneD, const NekDouble> &field,
                                    const Array<OneD, const NekDouble> &dphidt,
                                    Array<OneD, NekDouble> &dphidtint,
                                    Array<OneD, NekDouble> &TimeMap)
{
    int nq = GetTotPoints();

    const NekDouble phiTol = 10.0;
    const NekDouble phirest = 80.0;
    const NekDouble dphidtTol = 1.0;

    #pragma omp parallel for
    for (int i = 0; i < nq; ++i)
    {
        bool inZeroZone = false;
        for (int n = 0; n < m_numfiber; ++n)
        {
            if (zoneindex[n][i] == 0)
            {
                TimeMap[i] = 0.0;
                inZeroZone = true;
                break;
            }
        }
        if (inZeroZone)
        {
            continue;
        }

        const NekDouble phidiff = field[i] - phirest;
        const NekDouble dphi = dphidt[i];
        const NekDouble dint = dphidtint[i];

        if (phidiff > phiTol && dphi > dphidtTol)
        {
            const NekDouble fnewsum = dphi + dint;
            if (fnewsum > dphidtTol)
            {
                TimeMap[i] = (dphi * time + dint * TimeMap[i]) / fnewsum;
            }

            dphidtint[i] += dphi;
        }
    }
}

void MMFNeuralEP::ComputephieTimeMap(
    const NekDouble time,
    const Array<OneD, const NekDouble> &field,
    Array<OneD, NekDouble> &fieldint,
    Array<OneD, NekDouble> &TimeMap)
{
    const int nq = GetTotPoints();
    constexpr NekDouble Tol = 0.1;

    #pragma omp parallel for
    for (int i = 0; i < nq; ++i)
    {
        const NekDouble phie = field[i] + 3.0;
        if (phie > Tol)
        {
            const NekDouble fint = fieldint[i];
            const NekDouble fnewsum = phie + fint;

            TimeMap[i] = (phie * time + fint * TimeMap[i]) / fnewsum;
            fieldint[i] += phie;
        }
    }
}

void MMFNeuralEP::ComputerhoTimeMap(
    const NekDouble time,
    const Array<OneD, const NekDouble> &field,
    Array<OneD, NekDouble> &fieldint,
    Array<OneD, NekDouble> &TimeMap)
{
    const int nq = GetTotPoints();
    constexpr NekDouble Tol = 0.001;

    #pragma omp parallel for
    for (int i = 0; i < nq; ++i)
    {
        const NekDouble rho = field[i] + 0.005;
        if (rho > Tol)
        {
            const NekDouble fint = fieldint[i];
            const NekDouble fnewsum = rho + fint;

            TimeMap[i] = (rho * time + fint * TimeMap[i]) / fnewsum;
            fieldint[i] += rho;
        }
    }
}

// PhieCurrent[0] = \int \phi_m dt
// PhieCurrent[1] = \int CSD_e dt
// PhieCurrent[2] = \int \phim CSD_e dt / \int \phi_m dt
// PhieCurrent[3] = \int t CSD_e dt / \int CSD_e dt
void MMFNeuralEP::ComputePhieCurrent(const NekDouble time,
                                    const NekDouble timestep,
                                    const Array<OneD, const Array<OneD, NekDouble>> &field,
                                    Array<OneD, Array<OneD, NekDouble>> &PhieCurrent)
{
    int nq = GetTotPoints();

    NekDouble Tol = 0.1;

    Array<OneD, NekDouble> phim(nq);
    Array<OneD, NekDouble> phie(nq);

    Vmath::Vcopy(nq, &field[0][0], 1, &phim[0], 1);
    Vmath::Vcopy(nq, &field[1][0], 1, &phie[0], 1);

    // CSDe is negative such that it diminishes the total currenty by phi_m
    Array<OneD, NekDouble> CSDe = ComputeMMFDiffusion(m_movingframes, phie);

    Vmath::Vmul(nq, m_intrazone, 1, CSDe, 1, CSDe, 1);
    Vmath::Smul(nq, -1.0/(m_Cn * m_Rf), CSDe, 1, CSDe, 1);

    NekDouble phimsum, phiecurrentsum;
    NekDouble dtphim, dtCSDe;
    for (int i = 0; i < nq; ++i)
    {
        dtphim = timestep * phim[i];
        dtCSDe = timestep * (timestep * CSDe[i]);

        // If phi_m is positive, compute phie_weighted_current
        if (phim[i] > Tol)
        {
            phimsum = dtphim + PhieCurrent[0][i];
            PhieCurrent[2][i] = (dtphim * dtCSDe + PhieCurrent[0][i] * PhieCurrent[1][i]) / phimsum;
            PhieCurrent[0][i] = phimsum;
        }

        // If CSD_e is positive, compute phie_current
        if( CSDe[i] > Tol )
        {
            phiecurrentsum = dtCSDe + PhieCurrent[1][i];
            PhieCurrent[3][i] = (time * dtCSDe + PhieCurrent[3][i] * PhieCurrent[1][i]) / phiecurrentsum;            
            PhieCurrent[1][i] = phiecurrentsum;
        }
    }
}


void MMFNeuralEP::PlotNeuralEP(
    const Array<OneD, const Array<OneD, NekDouble>> &fields,
    const Array<OneD, const Array<OneD, NekDouble>> &TimeMap,
    const int nstep)
{
    const int nvar    = 7;
    const int nq      = m_fields[0]->GetTotPoints();
    const int ncoeffs = m_fields[0]->GetNcoeffs();

    std::string outname1 = m_sessionName + "_field_" +
                           boost::lexical_cast<std::string>(nstep) + ".chk";

    std::vector<Array<OneD, NekDouble>> fieldcoeffs(nvar);
    for (int i = 0; i < nvar; ++i)
    {
        fieldcoeffs[i] = Array<OneD, NekDouble>(ncoeffs);
    }

    std::vector<std::string> variables(nvar);
    variables[0] = "phi_m";
    variables[1] = "phi_e";
    variables[2] = "CSD";
    variables[3] = "rho";
    variables[4] = "TimeMap_phim";
    variables[5] = "TimeMap_phie";
    variables[6] = "TimeMap_rho";

    //     variables[0] = "phi_m";
    Array<OneD, NekDouble> tmp(nq);
    m_fields[0]->FwdTransLocalElmt(fields[0], fieldcoeffs[0]);

    Array<OneD, NekDouble> phim(nq), phie(nq), rho(nq);
    Vmath::Vmul(nq, m_intrazone, 1, fields[0], 1, phim, 1);
    Vmath::Vmul(nq, m_outerzone, 1, fields[1], 1, phie, 1);
    Vmath::Vmul(nq, m_outerzone, 1, fields[2], 1, rho, 1);

    //     variables[1] = "phi_e";
    m_fields[0]->FwdTransLocalElmt(phie, fieldcoeffs[1]);

    //     variables[2] = "CSD";
    Array<OneD, NekDouble> CSD = ComputeMMFDiffusion(m_phiediffmovingframes, fields[1]);
    Vmath::Smul(nq, -m_Diffext, CSD, 1, CSD, 1);
    m_fields[0]->FwdTransLocalElmt(CSD, fieldcoeffs[2]);

    //     variables[3] = "rho";
    m_fields[0]->FwdTransLocalElmt(rho, fieldcoeffs[3]);

    // Max values and indices
    const NekDouble Maxphim  = Vmath::Vmax(nq, phim, 1);
    const int       Maxphimid = Vmath::Imax(nq, phim, 1);
    const NekDouble Maxphie  = Vmath::Vmax(nq, phie, 1);
    const int       Maxphieid = Vmath::Imax(nq, phie, 1);
    const NekDouble MaxCSD   = Vmath::Vmax(nq, CSD, 1);
    const int       MaxCSDid  = Vmath::Imax(nq, CSD, 1);
    const NekDouble Maxrho   = Vmath::Vmax(nq, rho, 1);
    const int       Maxrhoid  = Vmath::Imax(nq, rho, 1);

    // Coordinates
    Array<OneD, NekDouble> x0(nq), x1(nq), x2(nq);
    m_fields[0]->GetCoords(x0, x1, x2);

    std::cout << "phim: Max = " << Maxphim << " at x = " << x0[Maxphimid] << ", y = " << x1[Maxphimid] << '\n';
    std::cout << "phie: Max = " << Maxphie << " at x = " << x0[Maxphieid] << ", y = " << x1[Maxphieid] << '\n';
    std::cout << "CSD: Max = " << MaxCSD << " at x = " << x0[MaxCSDid] << ", y = " << x1[MaxCSDid] << '\n';
    std::cout << "rho : Max = " << Maxrho  << " at x = " << x0[Maxrhoid]  << ", y = " << x1[Maxrhoid]  << "\n\n";

    Vmath::Vmul(nq, m_intrazone, 1, TimeMap[0], 1, phim, 1);
    Vmath::Vmul(nq, m_outerzone, 1, TimeMap[1], 1, phie, 1);
    Vmath::Vmul(nq, m_outerzone, 1, TimeMap[2], 1, rho, 1);
    // TimeMap maxima
    std::cout << "TimeMap: phim = " << Vmath::Vmax(nq, phim, 1)
              << ", phie = " << Vmath::Vmax(nq, phie, 1)
              << ", rho = "  << Vmath::Vmax(nq, rho, 1) << '\n';
    
    // variables[4] = "TimeMap_phim";
    m_fields[0]->FwdTransLocalElmt(TimeMap[0], fieldcoeffs[4]);

    // variables[5] = "TimeMap_phie";
    m_fields[0]->FwdTransLocalElmt(TimeMap[1], fieldcoeffs[5]);

    // variables[6] = "TimeMap_rho";
    m_fields[0]->FwdTransLocalElmt(TimeMap[2], fieldcoeffs[6]);

    WriteFld(outname1, m_fields[0], fieldcoeffs, variables);
}

void MMFNeuralEP::PlotAnisotropy(
    const Array<OneD, const Array<OneD, NekDouble>> &AniStrength, 
    const Array<OneD, const Array<OneD, NekDouble>> &phieAniStrength)
{
    int nvar    = 4;
    int ncoeffs = m_fields[0]->GetNcoeffs();

    std::string outname1 = m_sessionName + "_Anisotropy.chk";

    std::vector<Array<OneD, NekDouble>> fieldcoeffs(nvar);
    for (int i = 0; i < nvar; ++i)
    {
        fieldcoeffs[i] = Array<OneD, NekDouble>(ncoeffs);
    }

    std::vector<std::string> variables(nvar);
    variables[0] = "AniStrength1";
    variables[1] = "AniStrength2";
    variables[2] = "phieAniStrength1";
    variables[3] = "phieAniStrength2";

    // Compute the gradient of the time map
    m_fields[0]->FwdTransLocalElmt(AniStrength[0], fieldcoeffs[0]);
    m_fields[0]->FwdTransLocalElmt(AniStrength[1], fieldcoeffs[1]);

    m_fields[0]->FwdTransLocalElmt(phieAniStrength[0], fieldcoeffs[2]);
    m_fields[0]->FwdTransLocalElmt(phieAniStrength[1], fieldcoeffs[3]);

    WriteFld(outname1, m_fields[0], fieldcoeffs, variables);
}

void MMFNeuralEP::PrintSingleCurrent(const Array<OneD, const NekDouble> &phim,
                                  const Array<OneD, const NekDouble> &dudt,
                                  NekDouble &thredlocf1)
{
    int nq      = m_fields[0]->GetTotPoints();

    Array<OneD, NekDouble> x0(nq);
    Array<OneD, NekDouble> x1(nq);
    Array<OneD, NekDouble> x2(nq);

    m_fields[0]->GetCoords(x0, x1, x2);
    Array<OneD, NekDouble> phie(nq);
    phie = m_fields[1]->GetPhys();
    // Vmath::Vcopy(nq, m_fields[1]->GetPhys(), 1, phie, 1);

    Array<OneD, NekDouble> phieintra1(nq);

    Array<OneD, NekDouble> phieextra(nq);
    Vmath::Vmul(nq, m_intrazonefiber[0], 1, phie, 1, phieintra1, 1);

    Array<OneD, NekDouble> phimintra1(nq);
    Array<OneD, NekDouble> phimextra(nq);

    Vmath::Vmul(nq, m_intrazonefiber[0], 1, phim, 1, phimintra1, 1);
    Vmath::Vmul(nq, m_extrazone, 1, phim, 1, phimextra, 1);

    Array<OneD, NekDouble> dudtintra1(nq);
    Vmath::Vmul(nq, m_intrazonefiber[0], 1, dudt, 1, dudtintra1, 1);

    Vmath::Vmul(nq, m_extrazone, 1, phie, 1, phieextra, 1);
    
    Array<OneD, NekDouble> phimcurrent = ComputeMMFDiffusion(m_movingframes, phim);
    Array<OneD, NekDouble> phiecurrent = ComputeMMFDiffusion(m_movingframes, phie);

    Array<OneD, NekDouble> totcurrent1(nq);

    Array<OneD, NekDouble> phimcurrent1(nq);

    Vmath::Vmul(nq, m_intrazonefiber[0], 1, phimcurrent, 1, phimcurrent1, 1);

    Vmath::Vadd(nq, phimcurrent, 1, phiecurrent, 1, totcurrent1, 1);
    Vmath::Vmul(nq, m_intrazonefiber[0], 1, totcurrent1, 1, totcurrent1, 1);

    // index:0 -> u

    NekDouble Maxphim1 = Vmath::Vmax(nq, phimintra1, 1);
    NekDouble Maxphimextra = Vmath::Vmax(nq, phimextra, 1);

    int Maxphim1index = Vmath::Imax(nq, phimintra1, 1);
    int Maxphiextraindex = Vmath::Imax(nq, phimextra, 1);

    NekDouble phimMaxratio1 = 100.0 * Vmath::Vmax(nq, totcurrent1, 1) / Vmath::Vmax(nq, phimcurrent1, 1);
    NekDouble phimMinratio1 = 100.0 * Vmath::Vmin(nq, totcurrent1, 1) / Vmath::Vmin(nq, phimcurrent1, 1);

    // NekDouble Maxpositf1step = x1[Maxphim1index];
    // NekDouble Maxpositf2step = x1[Maxphim2index];

    std::cout << "fiber1: phim: Max = " << Maxphim1 << " at y = " << x1[Maxphim1index] << ", phimcurret1: Max = " << phimMaxratio1 << "  % , Min = " << phimMinratio1 << " % " << std::endl;
    std::cout << "Extraspace: phim: Max = " << Maxphimextra << " at x = " << x0[Maxphiextraindex] << std::endl;

    NekDouble tmp1;
    NekDouble x1loc=0.0;
    NekDouble dudtf1=0.0;
    for (int i=0;i<nq;++i)
    {
        tmp1 = phimintra1[i];
        
        if( (tmp1>m_phimrest) && (x1loc<x1[i]) )
        {
            x1loc = x1[i];
            dudtf1 = dudt[i];
        }
    }

    thredlocf1 = x1loc;

    std::cout << "fiber1: phim: thredloc at y = " << thredlocf1 << ", dudt = " << dudtf1 << std::endl;
}

void MMFNeuralEP::PrintDuoCurrent(const Array<OneD, const Array<OneD, NekDouble>> &field)
{
    int nq      = m_fields[0]->GetTotPoints();

    Array<OneD, NekDouble> x0(nq);
    Array<OneD, NekDouble> x1(nq);
    Array<OneD, NekDouble> x2(nq);

    m_fields[0]->GetCoords(x0, x1, x2);
    Array<OneD, NekDouble> phim(nq);
    Array<OneD, NekDouble> phie(nq);
    Array<OneD, NekDouble> rho(nq);

    Vmath::Vmul(nq, &m_intrazone[0], 1, &field[0][0], 1, &phim[0], 1);
    Vmath::Vmul(nq, &m_outerzone[0], 1, &field[1][0], 1, &phie[0], 1);
    Vmath::Vmul(nq, &m_outerzone[0], 1, &field[2][0], 1, &rho[0], 1);

    NekDouble Maxphim = Vmath::Vmax(nq, phim, 1);
    int Maxphimid = Vmath::Imax(nq, phim, 1);

    NekDouble Maxphie = Vmath::Vmax(nq, phie, 1);
    int Maxphieid = Vmath::Imax(nq, phie, 1);

    NekDouble Maxrho = Vmath::Vmax(nq, rho, 1);
    int Maxrhoid = Vmath::Imax(nq, rho, 1);

   std::cout << "phim: Max = " << Maxphim << " at x = " << x0[Maxphimid] << ", y = " << x1[Maxphimid] << std::endl;
   std::cout << "phie: Max = " << Maxphie << " at x = " << x0[Maxphieid] << ", y = " << x1[Maxphieid] << std::endl;
   std::cout << "rho: Max = " << Maxrho << " at x = " << x0[Maxrhoid] << ", y = " << x1[Maxrhoid] << std::endl << std::endl;
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


// void MMFNeuralEP::DisplayNode2D(std::string &fulltext, const Array<OneD, const Array<OneD, NekDouble>> &fields)
// {
//     int nq               = GetTotPoints();

//     Array<OneD, NekDouble> phi_m(nq);
//     Vmath::Vmul(nq, m_intrazone, 1, fields[0], 1, phi_m, 1);
    
//     Array<OneD, NekDouble> phi_e = Computephie(m_ExtCurrentType, fields[0]);

//     Array<OneD, NekDouble> phimavg(m_totNode+2,0.0);
//     Array<OneD, NekDouble> phieavg(m_totNode+2,0.0);

//     Array<OneD, int> totnodenumpts(m_totNode+2,0);

//     int nodeid;
//     for (int i=0; i<nq; ++i)
//     {
//         if(m_zoneindex[0][i]>=0)
//         {
//             nodeid = m_zoneindex[0][i];
//             phimavg[nodeid] += phi_m[i];
//             phieavg[nodeid] += phi_e[i];

//             totnodenumpts[nodeid] = totnodenumpts[nodeid] + 1;
//         }
//     }

//     for (int j=0; j<m_totNode+2; ++j)
//     {
//         if (totnodenumpts[j]>0)
//         {
//             phimavg[j] = phimavg[j] / totnodenumpts[j];
//             phieavg[j] = phieavg[j] / totnodenumpts[j];
//         }
//     }

//     fulltext.append(" \n");
//     fulltext.append("(Nodeid,phim,phie): ");
//     for (int i = 0; i < m_totNode+1; ++i)
//     {
//         fulltext.append( "( " + std::to_string(i) + " , " + std::to_string(phimavg[i]) + " , " + std::to_string(phieavg[i]) + " ) ");
//     }
    
//     fulltext.append(" \n");
// }

// void MMFNeuralEP::PlotHelmSolve(const Array<OneD, const NekDouble> &forcing,
//                                    const Array<OneD, const NekDouble> &solution)
// {
//     int nvar = 2;
//     int ncoeffs = m_fields[0]->GetNcoeffs();
//     // int nq      = m_fields[0]->GetTotPoints();

//     std::string outname;
//     outname = m_sessionName + "_helm.chk";

//     std::vector<Array<OneD, NekDouble>> fieldcoeffs(nvar);
//     for (int i = 0; i < nvar; ++i)
//     {
//         fieldcoeffs[i] = Array<OneD, NekDouble>(ncoeffs);
//     }

//     std::vector<std::string> variables(nvar);
//     variables[0] = "forcing";
//     variables[1] = "solution";

//     m_fields[0]->FwdTransLocalElmt(forcing, fieldcoeffs[0]);
//     m_fields[0]->FwdTransLocalElmt(solution, fieldcoeffs[1]);

//     WriteFld(outname, m_fields[0], fieldcoeffs, variables);
// }

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
    const int nq   = m_fields[0]->GetNpoints();

    // Set Helmholtz coefficients
    StdRegions::ConstFactorMap factors;
    factors[StdRegions::eFactorTau] = m_Helmtau;
    factors[StdRegions::eFactorLambda] = m_Cn * m_Rf / lambda;

    // SetBoundaryConditions(time);
    // SetMembraneBoundaryCondition();

    if (outarray[0].size() != nq)
    {
        outarray[0] = Array<OneD, NekDouble>(nq);
    }

    // Multiply 1.0/timestep
    const NekDouble scale = -factors[StdRegions::eFactorLambda];
    Vmath::Smul(nq, scale, inarray[0], 1, m_fields[0]->UpdatePhys(), 1);

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
    m_neuron->TimeIntegrate(m_zoneindexfiber[0], inarray[0], outarray[0], time, m_Temperature);

    for (int i=1; i < nvar; ++i)
    {
        Vmath::Vcopy(nq, m_neuron->GetNeuronSolution(i), 1, m_fields[i]->UpdatePhys(), 1);
    }

    m_stimulus[0]->Update(m_excitezonefiber[0], outarray[0], time);
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
    m_neuron->TimeIntegrate(m_zoneindexfiber[0], inarray[0], outarray[0], time, m_Temperature);

    m_stimulus[0]->Update(m_excitezonefiber[0], outarray[0], time);

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
            Vmath::Vmul(nq, &m_NeuralCmfiber[0], 1, &outarrayDiff[i][0], 1, &outarrayDiff[i][0], 1);

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
     m_neuron->TimeIntegrate(m_zoneindex, inarray[0], outarray[0], time, m_Temperature);

    for (int n=0; n<m_stimulus.size(); ++n)
    {
        m_stimulus[n]->Update(m_excitezonefiber[n], outarray[0], time);
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
    const int nvar = m_fields.size();
    const int nq   = m_fields[0]->GetNpoints();

    // Reuse memory if already allocated
    for (int i = 0; i < nvar; ++i)
    {
        if (outarray[i].size() != nq)
        {
            outarray[i] = Array<OneD, NekDouble>(nq, 0.0);
        }
        else
        {
            Vmath::Zero(nq, outarray[i], 1);
        }
    }

    // 1. Reaction Term (FHN or H-H ion current model)
    m_neuron->TimeIntegrate(m_zoneindex, inarray[0], outarray[0], time, m_Temperature);

    // 2. Apply Stimulus
    for (std::size_t n = 0; n < m_stimulus.size(); ++n)
    {
        m_stimulus[n]->Update(m_excitezonefiber[n], outarray[0], time);
    }

    // 3. Compute phi_e to satisfy bidomain coupling
    // \nabla \cdot ( (\signa_e + \sigma_i) \nabla \phi_e) = - \nabla \cdot
    // (\sigma_i \nabla \phi_m)
    Computephie(inarray[0]);
    const Array<OneD, NekDouble> &phie = m_fields[1]->GetPhys();

    // 4. Compute \nabla \cdot (\sigma_i \nabla \phi_e) and add to membrane current
    static thread_local Array<OneD, NekDouble> phiecurrent;
    if (phiecurrent.size() != nq)
        phiecurrent = Array<OneD, NekDouble>(nq);

    phiecurrent = ComputeMMFDiffusion(m_movingframes, phie);

    // Current caused by extracellular potential affects the total current at the nodes and myelin.
    #pragma omp parallel for
    for (int i = 0; i < nq; ++i)
    {
        if (m_intrazone[i] > 0.0)
        {
            outarray[0][i] += phiecurrent[i] / (m_Cn * m_Rf);
        }
    }

    // 5. Compute \nabla^2 (\nabla^2 \phi_e) = diffusion of extracellular field
    static thread_local Array<OneD, NekDouble> phiediff, phiediff2;
    if (phiediff.size() != nq)
    {
        phiediff = Array<OneD, NekDouble>(nq);
        phiediff2 = Array<OneD, NekDouble>(nq);
    }

    phiediff = ComputeMMFDiffusion(m_phiediffmovingframes, phie);
    phiediff2 = ComputeMMFDiffusion(m_phiediffmovingframes, phiediff);

    // Regular diffusivity value for the extracellular space
    // Physiological Review by Syková & Nicholson (2008)
    const NekDouble Deff = -5e-6; // μm²/ms
    #pragma omp parallel for
    for (int i = 0; i < nq; ++i)
    {
        outarray[2][i] = Deff * phiediff2[i];
    }

    if (m_explicitDiffusion)
    {
        static thread_local Array<OneD, NekDouble> Laplacian;
        if (Laplacian.size() != nq)
        {
            Laplacian = Array<OneD, NekDouble>(nq);
        }

        WeakDGMMFDiffusion(0, inarray[0], Laplacian, time);

        #pragma omp parallel for
        for (int i = 0; i < nq; ++i)
        {
            outarray[0][i] += Laplacian[i] / (m_Cn * m_Rf);
        }
    }
}

// output: phi_e (m_fields[1]->UpdatePhys()) and outarray (1/C_n/r) * \nabla^2 \phi_e
// Compute phi_e from the given distribution of phi_m
// \nabla \cdot ( (1 + \rho) \mathbf{e}_1 + \mathbf{e}_2 ) ( \nabla \phi_e ))
//                         = - \nabla \cdot \mathbf{e}_1 \nabla \phi_m
void MMFNeuralEP::Computephie(const Array<OneD, const NekDouble> &phim)
{
    const int nq = m_fields[0]->GetNpoints();

    // Solve the Poisson equation: \nabla (\sigma_e + \sigma_i ) phi_e = \nabla
    // \sigma_i \nabla phi_m
    StdRegions::ConstFactorMap phiefactors;
    phiefactors[StdRegions::eFactorTau]    = m_Helmtau;
    phiefactors[StdRegions::eFactorLambda] = 0.0;

    // // Compute \nabla \sigma_i \nabla phi_m and use it as point sources for
    // phi_e. This is equivalently achieved by removing all the point sources in
    // myelinnated fiber region.    
    // Allocate only once
    Array<OneD, NekDouble> phimcurrent(nq);
    if (phimcurrent.size() != nq)
    {
        phimcurrent = Array<OneD, NekDouble>(nq, 0.0);
    }
    else
    {
        Vmath::Zero(nq, phimcurrent, 1);
    }

    phimcurrent = ComputeMMFDiffusion(m_movingframes, phim);
    Vmath::Neg(nq, phimcurrent, 1);

    switch(m_ExtCurrentType)
    {
        case eEphaptic:
        {
            Vmath::Vmul(nq, m_nodezone, 1, phimcurrent, 1, phimcurrent, 1);
            break;
        }

        case eNoEphaptic:
        {
            Vmath::Zero(nq, phimcurrent, 1);
            break;
        }

        default:
        break;
    }

    // Compute  \nabla \cdot ( (1 + \rho) \mathbf{e}_1 + \mathbf{e}_2 ) ( \nabla \phi_e ))
    //                         = - \nabla \cdot \mathbf{e}_1 \nabla \phi_m
    const NekDouble avg = AvgInt(phimcurrent);
    #pragma omp parallel for
    for (int i = 0; i < nq; ++i)
    {
        m_fields[1]->UpdatePhys()[i] = phimcurrent[i] - avg;
    }

    m_fields[1]->HelmSolve(m_fields[1]->GetPhys(), m_fields[1]->UpdateCoeffs(), phiefactors, m_phievarcoeff);
    m_fields[1]->BwdTrans(m_fields[1]->GetCoeffs(), m_fields[1]->UpdatePhys());

    // Make it as a value with AvgInt is zero.
    const NekDouble mean = AvgInt(m_fields[1]->GetPhys());
    Vmath::Sadd(nq, -mean, m_fields[1]->GetPhys(), 1, m_fields[1]->UpdatePhys(), 1);

    m_fields[1]->SetPhysState(true);
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

            Array<OneD, NekDouble> initialcondition(nq, 0.0);
            for (int n=0; n<m_numfiber; ++n)
            {
               // Vmath::Svtvp(nq, 100.0, m_excitezone[n], 1, initialcondition, 1, initialcondition, 1);
            }

            // for (int n=0; n<m_numfiber; ++n)
            // {
            //    m_stimulus[n]->Update(m_excitezone[n], initialcondition, initialtime);
            // }

            // Array<OneD, NekDouble> initialcondition(nq, 0.0);
            // for (int n=0; n<m_numfiber; ++n)
            // {
            //     m_stimulus[n]->Update(m_excitezone[n], initialcondition, initialtime);
            // }

            m_fields[0]->SetPhys(initialcondition);
            break;
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


// Change TimeMap to Velocity
    Array<OneD, NekDouble> MMFNeuralEP::ConvertTMtoVel(
        const Array<OneD, const NekDouble> &TimeMap,
        const Array<OneD, const NekDouble> &TmapGrad,
        const Array<OneD, const NekDouble> &TmapGradMag)
{
    int nq = m_fields[0]->GetTotPoints();

    // Array<OneD, NekDouble> outarray(m_spacedim * nq, 0.0);
    Array<OneD, NekDouble> outarray(m_spacedim * nq);

    // Compute VelField \vec{v} = \sum_{i=1}^3 1/(\nabla T \cdot \hat{x}_i)
    // \hat{x}_i
    NekDouble TMgradmag, TM, TMgradrel;
    NekDouble TmapTol = 0.1;
    for (int i = 0; i < nq; i++)
    {
        TMgradmag = TmapGradMag[i];
        TM = TimeMap[i];

        if( ( TMgradmag > TmapTol) && ( TM > TmapTol) )
        {
            TMgradrel = TMgradmag / TM ;

            if (TMgradrel > TmapTol)
            {
                for (int k = 0; k < m_spacedim; ++k)
                {
                    outarray[k*nq + i] = TmapGrad[ k* nq + i] / (TMgradmag * TMgradmag);
                }
            }
        }

        // if ( floor(ValidTimeMap[i]) == 0)
        // {
        //     outarray[i]        = 0.0;
        //     outarray[nq+i]     = 0.0;
        //     outarray[2*nq+i]   = 0.0;
        // }
    }
    
    Array<OneD, NekDouble> tmpx(nq);
    Array<OneD, NekDouble> tmpy(nq);
    Array<OneD, NekDouble> tmpz(nq);

    Vmath::Vcopy(nq, &outarray[0], 1, &tmpx[0], 1);
    Vmath::Vcopy(nq, &outarray[nq], 1, &tmpy[0], 1);
    Vmath::Vcopy(nq, &outarray[2*nq], 1, &tmpz[0], 1);

    // Print out
    std::cout << "(Vx, Vy, Vz) = ( " << RootMeanSquare(tmpx) << " , " 
    << RootMeanSquare(tmpy) << " , " << RootMeanSquare(tmpz) << " ) " << std::endl; 

    return outarray;
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

    SolverUtils::AddSummaryItem(s, "FiberType", FiberTypeMap[m_FiberType]);

    SolverUtils::AddSummaryItem(s, "bundleleft", m_bundleleft);
    SolverUtils::AddSummaryItem(s, "bundleright", m_bundleright);

    SolverUtils::AddSummaryItem(s, "FiberAngle", m_fiberangle);
    SolverUtils::AddSummaryItem(s, "FiberWidth", m_fiberwidth);
    SolverUtils::AddSummaryItem(s, "FiberGap", m_fibergap);
    SolverUtils::AddSummaryItem(s, "Radiusfiberbundle", m_radiusfiberbundle);
    SolverUtils::AddSummaryItem(s, "FiberCurvature", m_fibercurvature);
    SolverUtils::AddSummaryItem(s, "phiefactor", m_phiefactor);

    SolverUtils::AddSummaryItem(s, "Node Length", m_nodelen);
    SolverUtils::AddSummaryItem(s, "Myelin Length", m_myelinlen);

    SolverUtils::AddSummaryItem(s, "Number_Fiber", m_numfiber);
    SolverUtils::AddSummaryItem(s, "Total_Number_Node", m_totNode);

    // SolverUtils::AddSummaryItem(s, "phimrest", m_phimrest);
    // SolverUtils::AddSummaryItem(s, "phimTol", m_phimTol);
    // SolverUtils::AddSummaryItem(s, "dphimdtTol", m_dphimdtTol);

    SolverUtils::AddSummaryItem(s, "Temperature", m_Temperature);
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
