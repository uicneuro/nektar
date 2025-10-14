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

#include <SolverUtils/MMFSystem.h>
#include <MMFSolver/EquationSystems/MMFNeuralEP.h>

#include <CardiacEPSolver/Filters/FilterCellHistoryPoints.h>
#include <CardiacEPSolver/Filters/FilterCheckpointCellModel.h>

#include <SpatialDomains/MeshGraphIO.h>
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

    const int nq = m_fields[0]->GetNpoints();
    const int nvar = m_fields.size();
    switch (nvar)
    {
        // NeuralEP2Dbi model
        case 2:
        {
            m_phimvar = 0;
            m_phievar = 1;
            break;
        }

        // NeuralEP2DbiMulti model for two fibers with a signel phie
        case 3:
        {
            if( (m_NeuralEPType == eNeuralEP2DbiMulti) || (m_NeuralEPType == eNeuralEP2DbiMultiv2) )
            {
                m_phimvar = 1;
            }
            else if(m_NeuralEPType == eNeuralEP2DbiCSD)
            {
                m_phimvar = 0;
            }
            m_phievar = 2;
            break;
        }

        // NeuralEP2DbiMulti model for two fibers with two phies
        case 4:
        {
            m_phievar = 2;
            break;
        }

        default:
        break;
    }

    static constexpr NekDouble PI = 3.14159265358979323846;
    m_pi       = PI;

    // Get the coordinates
    m_x = Array<OneD, NekDouble>(nq);
    m_y = Array<OneD, NekDouble>(nq);
    m_z = Array<OneD, NekDouble>(nq);

    m_fields[0]->GetCoords(m_x, m_y, m_z);

    // Compute the center of the element as the coordinate of each grid.
    m_xcell = Array<OneD, NekDouble>(nq);
    m_ycell = Array<OneD, NekDouble>(nq);
    m_zcell = Array<OneD, NekDouble>(nq);

    GetCellCoordAvg(m_xcell, m_ycell, m_zcell);

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

    m_session->LoadParameter("fiberangle", m_fiberangle, 0.0);
    m_session->LoadParameter("FiberWidth", m_fiberwidth, 0.01);
    m_session->LoadParameter("FiberGap", m_fibergap, 0.01);
    m_session->LoadParameter("FiberCurvature", m_fibercurvature, 0.8);
    m_session->LoadParameter("FiberAngle", m_fiberangle, m_pi/6.0);

    m_session->LoadParameter("NodeLength", m_nodelen, 0.01);
    m_session->LoadParameter("MyelinLength", m_myelinlen, 0.2);

    m_session->LoadParameter("Total_Number_Node", m_totNode, 3);
    m_session->LoadParameter("nodeinitdown", m_nodeinitdown, 0.01);
    m_session->LoadParameter("nodeinitup", m_nodeinitup, 0.02);

    // Total Fiber Length
    m_fiberlength = 3*m_nodelen + m_totNode*(m_nodelen+m_myelinlen);
    std::cout << "\nFiber Length = " << m_fiberlength << std::endl;

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

    m_session->LoadParameter("CSDDiff", m_CSDDiff, 1e-6);

    NekDouble axoncrossA = m_pi*m_axondiameter*m_axondiameter;
    NekDouble PhieMultFactorlower = m_relfiberratio*m_gratio*m_gratio*m_radiusfiberbundle*m_radiusfiberbundle;
    NekDouble PhieMultFactor = m_axondiameter*m_axondiameter/PhieMultFactorlower;

    // 1.0 /(m_pi * m_relfiberratio*m_gratio*m_gratio*m_radiusfiberbundle*m_radiusfiberbundle)
    m_phiefactor = PhieMultFactor / axoncrossA;
    std::cout << "phiefactor = " << m_phiefactor << std::endl;

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

    m_neuron = GetNeuronModelFactory().CreateInstance(
        vNeuronModel, m_session, m_fields[0]);

    // Rf and Cn are imported
    m_Rf = m_neuron->GetRecistanceValue();
    m_Cm = m_neuron->GetCapacitanceValue(0);
    m_Cn = m_neuron->GetCapacitanceValue(1);

    m_AnisotropyStrength = m_Cn / m_Cm;

    if( (m_NeuralEPType==eNeuralEP2DbiMulti) || (m_NeuralEPType==eNeuralEP2DbiMultiv2) )
    {   
        ASSERTL0(m_fields.size()==3, "Number of Variable should be 3");
    }
 
   switch (m_NeuralEPType)
    {
        case eNeuralHelmSolveSingle:
        case eNeuralHelmSolveDuo:
        {
            SetUpNeuralCm();
            break;
        }

        case eNeuralEP2Dmono:
        case eNeuralEP2Dbi:
        case eNeuralEP2DbiMulti:
        case eNeuralEP2DbiMultiv2:
        case eNeuralEP2DbiCSD:
        {   
            // Provide index for nodes and myelins
            IndexNodeZone2D(m_fiberleft, m_fiberright, m_fiberorder, m_zoneindexfiber);

            // Identifying each domain zone accordingly
            SetUpDomainZone();
            break;
        }

        default:
         break;
    }

    // Stimulus
    m_stimulus = NeuralStimulus::LoadStimuli(m_session, m_fields[0]);

   // Compute Anisotropy Strength and conductivity accordingly 
   ComputeAniStrengthfiber(m_zoneindexfiber, m_AniStrengthfiber);
   ComputeGlobalAniStrength(m_AniStrengthfiber, m_AniStrength);

   // Compute phie Anisotropy Strength and conductivity accordingly
   ComputePhieAniStrengthfiber(m_zoneindexfiber, m_phieAniStrengthfiber);
   // ComputeGlobalPhieAniStrength(m_zoneindex, m_phieAniStrength);
   ComputeGlobalPhieAniStrengthfiber(m_zoneindex, m_zoneindexfiber, m_phieAniStrength);


   Array<OneD, Array<OneD, NekDouble>> UnitAnisotroy(m_expdim);
   for (int j = 0; j < m_expdim; ++j)
   {
       UnitAnisotroy[j] = Array<OneD, NekDouble>(nq, 1.0);
   }

   // Create moving frames for phi_m
    MMFSystem::MMFInitObject(UnitAnisotroy);
    CheckMovingFrames(m_movingframes);

    // Setting up moving frames for CSD
    m_CSDmovingframes = Array<OneD, Array<OneD, NekDouble>>(m_spacedim);
    for (int j = 0; j < m_spacedim; ++j)
    {
        m_CSDmovingframes[j] = Array<OneD, NekDouble>(m_spacedim * nq);
        Vmath::Vcopy(m_spacedim * nq, &m_movingframes[j][0], 1, &m_CSDmovingframes[j][0], 1);
    }
    ComputeVarCoeff2D(m_CSDmovingframes, m_CSDvarcoeff);

    // Create moving frames for phi_e
    std::string phieMMFdirStr;
    m_session->LoadSolverInfo("phieMMFDir", phieMMFdirStr, "TangentY");
    SpatialDomains::GeomMMF phieMMFdir = FindMMFdir(phieMMFdirStr);
    SetUpMovingFrames(phieMMFdir, UnitAnisotroy, m_phiemovingframes);

    // Computer moving frames along each fiber  
    std::cout << "\nA: Computemovingframesfiber: numfiber = " << m_numfiber << std::endl;
    Computemovingframesfiber(m_movingframes, m_AniStrengthfiber, m_movingframesfiber);

    // Rescale movingframes
    std::cout << "\nB: Constructing movingframes" << std::endl;
    RescaleMovingFrames(m_AniStrength, m_movingframes);
    CheckMovingFrames(m_movingframes);
    ComputeVarCoeff2D(m_movingframes, m_varcoeff);

    // Rescale movingframes for each fiber
    m_varcoefffiber = Array<OneD, StdRegions::VarCoeffMap>(m_numfiber);
    for (int j = 0; j < m_numfiber; ++j)
    {
        std::cout << "\n C:Constructing m_varcoefffiber, fiber = " << j << std::endl;
        ComputeVarCoeff2D(m_movingframesfiber[j], m_varcoefffiber[j]);
    }

    std::cout << "\nD: Computephiemovingframesfiber = " << m_numfiber << std::endl;
    Computephiemovingframesfiber(m_phiemovingframes, m_phieAniStrengthfiber, m_phiemovingframesfiber);

    RescaleMovingFrames(m_phieAniStrength, m_phiemovingframes);
    ComputeVarCoeff2D(m_phiemovingframes, m_phievarcoeff);
    CheckMovingFrames(m_phiemovingframes);

    // Rescale phiemovingframes for each fiber
    m_phievarcoefffiber = Array<OneD, StdRegions::VarCoeffMap>(m_numfiber);
    for (int j = 0; j < m_numfiber; ++j)
    {
        std::cout << "\n E: Constructing m_phievarcoefffiber, fiber = " << j << std::endl;
        ComputeVarCoeff2D(m_phiemovingframesfiber[j], m_phievarcoefffiber[j]);
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
            case eNeuralHelmSolveSingle:
            case eNeuralHelmSolveDuo:
            {
                TestHelmSolve();
                wait_on_enter();
                break;
            }

            case eNeuralEP2Dmono:
            {
                m_ode.DefineImplicitSolve(&MMFNeuralEP::DoImplicitSolveNeuralEP2Dmono, this);
                m_ode.DefineOdeRhs(&MMFNeuralEP::DoOdeRhsNeuralEP2Dmono, this);
                break;
            }

            case eNeuralEP2Dbi:
            {
                m_ode.DefineImplicitSolve(&MMFNeuralEP::DoImplicitSolveNeuralEP2Dbi, this); 
                m_ode.DefineOdeRhs(&MMFNeuralEP::DoOdeRhsNeuralEP2Dbi, this);
                break;
            }

            case eNeuralEP2DbiMulti:
            {
                m_ode.DefineImplicitSolve(&MMFNeuralEP::DoImplicitSolveNeuralEP2DbiMulti, this); 
                m_ode.DefineOdeRhs(&MMFNeuralEP::DoOdeRhsNeuralEP2DbiMulti, this);
                break;
            }

            case eNeuralEP2DbiMultiv2:
            {
                m_ode.DefineImplicitSolve(&MMFNeuralEP::DoImplicitSolveNeuralEP2DbiMultiv2, this); 
                m_ode.DefineOdeRhs(&MMFNeuralEP::DoOdeRhsNeuralEP2DbiMultiv2, this);
                break;
            }

            case eNeuralEP2DbiCSD:
            {
                m_ode.DefineImplicitSolve(&MMFNeuralEP::DoImplicitSolveNeuralEP2DbiCSD, this); 
                m_ode.DefineOdeRhs(&MMFNeuralEP::DoOdeRhsNeuralEP2DbiCSD, this);
                break;
            }

            default:
                break;
        }
    }
}

MMFNeuralEP::~MMFNeuralEP()
{
}


void MMFNeuralEP::ComputeGlobalPhieAniStrength(
    const Array<OneD, const int> &zoneindex,
    Array<OneD, Array<OneD, NekDouble>> &phieAniStrength)
{
    int nq = GetTotPoints();
    int expdim = m_expdim;

    Array<OneD, Array<OneD, NekDouble>> sigma_i(expdim);
    Array<OneD, Array<OneD, NekDouble>> sigma_e(expdim);
    Array<OneD, Array<OneD, NekDouble>> sigma_eM(expdim);
    for (int j = 0; j < expdim; ++j)
    {
        sigma_i[j] = Array<OneD, NekDouble>(nq, 0.0);
        sigma_e[j] = Array<OneD, NekDouble>(nq, 0.0);
        sigma_eM[j] = Array<OneD, NekDouble>(nq, 0.0);
    }

    ComputeRegionalSigma(zoneindex, sigma_i, sigma_e, sigma_eM);

    phieAniStrength = Array<OneD, Array<OneD, NekDouble>>(expdim);
    for (int j = 0; j < expdim; ++j)
    {
        phieAniStrength[j] = Array<OneD, NekDouble>(nq, 0.0);
        Vmath::Vadd(nq, sigma_i[j], 1, sigma_eM[j], 1, phieAniStrength[j], 1);
    }
}

void MMFNeuralEP::ComputeGlobalPhieAniStrengthfiber(
    const Array<OneD, const int> &zoneindex,
    const Array<OneD, const Array<OneD, int>> &zoneindexfiber,
    Array<OneD, Array<OneD, NekDouble>> &phieAniStrength)
{
    int nq = GetTotPoints();
    int expdim = m_expdim;

    Array<OneD, Array<OneD, NekDouble>> sigma_i(expdim);
    Array<OneD, Array<OneD, NekDouble>> sigma_e(expdim);
    Array<OneD, Array<OneD, NekDouble>> sigma_eM(expdim);
    for (int j = 0; j < expdim; ++j)
    {
        sigma_i[j] = Array<OneD, NekDouble>(nq, 0.0);
        sigma_e[j] = Array<OneD, NekDouble>(nq, 0.0);
        sigma_eM[j] = Array<OneD, NekDouble>(nq, 0.0);
    }

    Array<OneD, Array<OneD, NekDouble>> tmp_i(expdim);
    Array<OneD, Array<OneD, NekDouble>> tmp_e(expdim);
    Array<OneD, Array<OneD, NekDouble>> tmp_eM(expdim);
    for (int j = 0; j < expdim; ++j)
    {
        tmp_i[j] = Array<OneD, NekDouble>(nq, 0.0);
        tmp_e[j] = Array<OneD, NekDouble>(nq, 0.0);
        tmp_eM[j] = Array<OneD, NekDouble>(nq, 0.0);
    }

    // ComputeRegionalSigma(zoneindexfiber, sigma_i, sigma_e, sigma_eM);
    for (int j = 0; j < m_numfiber; ++j)
    {
        ComputeRegionalSigma(zoneindexfiber[j], tmp_i, tmp_e, tmp_eM);
        for (int k = 0; k < m_expdim; ++k)
        {
            Vmath::Vadd(nq, &tmp_i[k][0], 1, &sigma_i[k][0], 1, &sigma_i[k][0], 1);
        }
    }

    ComputeRegionalSigma(zoneindex, tmp_i, sigma_e, sigma_eM);

    phieAniStrength = Array<OneD, Array<OneD, NekDouble>>(expdim);
    for (int j = 0; j < expdim; ++j)
    {
        phieAniStrength[j] = Array<OneD, NekDouble>(nq, 0.0);
        Vmath::Vadd(nq, sigma_i[j], 1, sigma_eM[j], 1, phieAniStrength[j], 1);
    }
}

void MMFNeuralEP::PlotAnisotropyfiber(
        const Array<OneD, const Array<OneD, NekDouble>> &AniStrength, 
        const Array<OneD, const Array<OneD, Array<OneD, NekDouble>>> &AniStrengthfiber, 
        const Array<OneD, const Array<OneD, NekDouble>> &phieAniStrength, 
        const Array<OneD, const Array<OneD, Array<OneD, NekDouble>>> &phieAniStrengthfiber)
{
    int nvar    = m_numfiber * 2 + 2;
    int ncoeffs = m_fields[0]->GetNcoeffs();

    std::string outname1 = m_sessionName + "_Anisotropyfiber.chk";

    std::vector<Array<OneD, NekDouble>> fieldcoeffs(nvar);
    for (int i = 0; i < nvar; ++i)
    {
        fieldcoeffs[i] = Array<OneD, NekDouble>(ncoeffs);
    }

    std::vector<std::string> variables(nvar);
    variables[0] = "AniStrength";
    variables[1] = "AniStrengthfiber1";
    variables[2] = "AniStrengthfiber2";

    variables[3] = "phieAniStrength";
    variables[4] = "phieAniStrengthfiber1";
    variables[5] = "phieAniStrengthfiber2";

    // Compute the gradient of the time map
    m_fields[0]->FwdTransLocalElmt(AniStrength[0], fieldcoeffs[0]);
    m_fields[0]->FwdTransLocalElmt(AniStrengthfiber[0][0], fieldcoeffs[1]);
    m_fields[0]->FwdTransLocalElmt(AniStrengthfiber[1][0], fieldcoeffs[2]);

    m_fields[0]->FwdTransLocalElmt(phieAniStrength[0], fieldcoeffs[3]);
    m_fields[0]->FwdTransLocalElmt(phieAniStrengthfiber[0][0], fieldcoeffs[4]);
    m_fields[0]->FwdTransLocalElmt(phieAniStrengthfiber[1][0], fieldcoeffs[5]);

    WriteFld(outname1, m_fields[0], fieldcoeffs, variables);
}


void MMFNeuralEP::RescaleMovingFrames(
    const Array<OneD, const Array<OneD, NekDouble>> &AniStrength,
    Array<OneD, Array<OneD, NekDouble>> &movingframes)
{
    const int nq       = GetTotPoints();
    const int spacedim = m_spacedim;
    const int expdim   = m_expdim;

    for (int j = 0; j < expdim; ++j)
    {
        for (int k = 0; k < spacedim; ++k)
        {
            for (int i = 0; i < nq; ++i)
            {
                movingframes[j][i + k * nq] =
                    sqrt(AniStrength[j][i]) * movingframes[j][i + k * nq];
            }
        }
    }
}

void MMFNeuralEP::Computemovingframesfiber(
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    const Array<OneD, const Array<OneD, Array<OneD, NekDouble>>> &AniStrengthfiber,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &movingframesfiber)
{
const int numfiber = m_numfiber;
const int spacedim = m_spacedim;
const int nq = GetTotPoints();

const int expdim = m_expdim;

Array<OneD, NekDouble> tmp(nq);
m_movingframesfiber = Array<OneD, Array<OneD, Array<OneD, NekDouble>>>(numfiber);
for (int n = 0; n < numfiber; ++n)
{
    movingframesfiber[n] = Array<OneD, Array<OneD, NekDouble>>(spacedim);
    for (int j = 0; j < expdim; ++j)
    {
        Vmath::Vsqrt(nq, &AniStrengthfiber[n][j][0], 1, &tmp[0], 1);

        movingframesfiber[n][j] = Array<OneD, NekDouble>(spacedim * nq);
        for (int k = 0; k < spacedim; ++k)
        {
            Vmath::Vmul(nq, &tmp[0], 1, &movingframes[j][k*nq], 1, 
                &movingframesfiber[n][j][k*nq], 1);
        }
    }

    // No anisotropy along the surface normal direction
    movingframesfiber[n][expdim] = Array<OneD, NekDouble>(spacedim * nq);
    for (int k = 0; k < spacedim; ++k)
    {
        Vmath::Vcopy(nq, &movingframes[expdim][k*nq], 1, &movingframesfiber[n][expdim][k*nq], 1);
    }

    std::cout << "\nDone: Checking movingframes fiber = " << n << std::endl;
    CheckMovingFrames(movingframesfiber[n]);
}
}


void MMFNeuralEP::Computephiemovingframesfiber(
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    const Array<OneD, const Array<OneD, Array<OneD, NekDouble>>> &AniStrengthfiber,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &phiemovingframesfiber)
{
const int numfiber = m_numfiber;
const int spacedim = m_spacedim;
const int nq = GetTotPoints();

const int expdim = m_expdim;

phiemovingframesfiber = Array<OneD, Array<OneD, Array<OneD, NekDouble>>>(numfiber);
for (int n = 0; n < numfiber; ++n)
{
    phiemovingframesfiber[n] = Array<OneD, Array<OneD, NekDouble>>(spacedim);
    for (int j = 0; j < expdim; ++j)
    {
        phiemovingframesfiber[n][j] = Array<OneD, NekDouble>(spacedim * nq);
    }
}

Array<OneD, NekDouble> tmp(nq, 1.0);
phiemovingframesfiber = Array<OneD, Array<OneD, Array<OneD, NekDouble>>>(numfiber);
for (int n = 0; n < numfiber; ++n)
{
    phiemovingframesfiber[n] = Array<OneD, Array<OneD, NekDouble>>(spacedim);
    for (int j = 0; j < expdim; ++j)
    {
        Vmath::Vsqrt(nq, &AniStrengthfiber[n][j][0], 1, &tmp[0], 1);
 
        phiemovingframesfiber[n][j] = Array<OneD, NekDouble>(spacedim * nq);
        for (int k = 0; k < spacedim; ++k)
        {
            Vmath::Vmul(nq, &tmp[0], 1, &movingframes[j][k*nq], 1, 
                &phiemovingframesfiber[n][j][k*nq], 1);
        }
    }

    // No anisotropy along the surface normal direction
    phiemovingframesfiber[n][expdim] = Array<OneD, NekDouble>(spacedim * nq);
    for (int k = 0; k < spacedim; ++k)
    {
        Vmath::Vcopy(nq, &movingframes[expdim][k*nq], 1, &phiemovingframesfiber[n][expdim][k*nq], 1);
    }

    std::cout << "\nDone: Checking movingframes fiber = " << n << std::endl;
    CheckMovingFrames(phiemovingframesfiber[n]);
}
}


void MMFNeuralEP::ComputePhieAniStrengthfiber(
    const Array<OneD, const Array<OneD, int>> &zoneindexfiber, 
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &phieAniStrengthfiber)
{
    const int numfiber = m_numfiber;
    const int expdim = m_expdim;
    const int nq = GetTotPoints();

    // Compute phieAniStrengthfiber
    Array<OneD, Array<OneD, NekDouble>> sigma_i(numfiber);
    Array<OneD, Array<OneD, NekDouble>> sigma_e(numfiber);
    Array<OneD, Array<OneD, NekDouble>> sigma_eM(numfiber);
    for (int k = 0; k < numfiber; ++k)
    {
        sigma_i[k] = Array<OneD, NekDouble>(nq);
        sigma_e[k] = Array<OneD, NekDouble>(nq);
        sigma_eM[k] = Array<OneD, NekDouble>(nq);
    }

    phieAniStrengthfiber = Array<OneD, Array<OneD, Array<OneD, NekDouble>>>(numfiber);
    for (int j = 0; j < numfiber; ++j)
    {
        phieAniStrengthfiber[j] = Array<OneD, Array<OneD, NekDouble>>(expdim);
        for (int k = 0; k < expdim; ++k)
        {
            phieAniStrengthfiber[j][k] = Array<OneD, NekDouble>(nq, 1.0);
        }

        ComputeRegionalSigma(zoneindexfiber[j], sigma_i, sigma_e, sigma_eM);
        for (int k = 0; k < m_expdim; ++k)
        {
            Vmath::Vadd(nq, &sigma_i[k][0], 1, &sigma_eM[k][0], 1, &phieAniStrengthfiber[j][k][0], 1);
        }
    }
}

// Compute conductivity accordingly 
void MMFNeuralEP::ComputeAniStrengthfiber(
    const Array<OneD, const Array<OneD, int>> &zoneindexfiber, 
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &AniStrengthfiber)
    {
        const int numfiber = m_numfiber;
        const int expdim = m_expdim;
        const int nq = GetTotPoints();    
        const NekDouble Cn = m_Cn;

        Array<OneD, Array<OneD, Array<OneD, NekDouble>>> NeuralCmfiber(numfiber);
        for (int j = 0; j < numfiber; ++j)
        {
            NeuralCmfiber[j] = Array<OneD, Array<OneD, NekDouble>>(m_expdim);
            for (int k = 0; k < m_expdim; ++k)
            {
                NeuralCmfiber[j][k] = Array<OneD, NekDouble>(nq, 1.0);
            }
        }

        ComputeNeuralCmfiber(zoneindexfiber, NeuralCmfiber);

        AniStrengthfiber = Array<OneD, Array<OneD, Array<OneD, NekDouble>>>(numfiber);
        for (int j = 0; j < numfiber; ++j)
        {
            AniStrengthfiber[j] = Array<OneD, Array<OneD, NekDouble>>(m_expdim);
            for (int k = 0; k < expdim; ++k)
            {
                AniStrengthfiber[j][k] = Array<OneD, NekDouble>(nq);
                Vmath::Smul(nq, Cn, &NeuralCmfiber[j][k][0], 1, &AniStrengthfiber[j][k][0], 1);
            }
        }
    }


    void MMFNeuralEP::ComputeNeuralCmfiber(
            const Array<OneD, const Array<OneD, int>> &zoneindexfiber,
            Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &outarrayfiber)
    {
        const int nq       = GetTotPoints();
        const int numfiber = m_numfiber;

        const int npts = m_fields[0]->GetTotPoints(0);

        const NekDouble Cm = m_neuron->GetCapacitanceValue(0);
        const NekDouble Cn = m_neuron->GetCapacitanceValue(1);

        outarrayfiber =
            Array<OneD, Array<OneD, Array<OneD, NekDouble>>>(numfiber);
        for (int n = 0; n < numfiber; ++n)
        {
            outarrayfiber[n] = Array<OneD, Array<OneD, NekDouble>>(m_expdim);
            for (int k = 0; k < m_expdim; ++k)
            {
                outarrayfiber[n][k] = Array<OneD, NekDouble>(nq, 0.0);
            }
        }

        int index, cntm = 0, cntn = 0, cnte = 0;
        for (int n = 0; n < numfiber; ++n)
        {
            cntn = 0;
            cntm = 0;
            cnte = 0;
            for (int k = 0; k < m_expdim; ++k)
            {
                for (int i = 0; i < nq; ++i)
                {
                    index = zoneindexfiber[n][i];

                    // Ranvier node zone for all the fibers
                    if (index >= 0)
                    {
                        outarrayfiber[n][k][i] = 1.0 / Cn;
                        cntn++;
                    }

                    // Ranvier node zone for all the fibers
                    if (index == -1)
                    {
                        if (m_MediumType == eAnisotropy)
                        {
                            outarrayfiber[n][k][i] = 1.0 / Cm;
                        }

                        else
                        {
                            outarrayfiber[n][k][i] = 1.0 / Cn;
                        }
                        cntm++;
                    }
                }
            }

            cnte = nq - cntn - cntm;

            std::cout << "fiber = " << n
                      << ", ComputeConductivity: Node = " << cntn / npts
                      << ", Myelinf1 = " << cntm / npts
                      << ", extracell = " << cnte / npts << std::endl;
        }
    }



Array<OneD, NekDouble> MMFNeuralEP::ComputeNeuralCmfiber(
    const Array<OneD, const int> &zoneindex)
{
    const int nq   = GetTotPoints();
    const int npts = m_fields[0]->GetTotPoints(0);

    const NekDouble Cm = m_neuron->GetCapacitanceValue(0);
    const NekDouble Cn = m_neuron->GetCapacitanceValue(1);

    Array<OneD, NekDouble> outarray(nq);

    int cntn = 0, cntm = 0, cnte = 0;
    for (int i = 0; i < nq; ++i)
    {
        // Ranvier node zone
        if (zoneindex[i] >= 0)

        {
            outarray[i] = 1.0 / Cn;
            cntn++;
        }

        // Myelin zone
        if (zoneindex[i] == -1)
        {
            outarray[i] = 1.0 / Cm;
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

    std::cout << "ComputeNeuralCmfiber: Node = " << cntn / npts
              << ", Myelinf1 = " << cntm / npts
              << ", extracell = " << cnte / npts << std::endl;

    return outarray;
}

// When a node or myeline exist, it is prevalent over the extracellular region
void MMFNeuralEP::ComputeGlobalAniStrength(
    const Array<OneD, const Array<OneD, Array<OneD, NekDouble>>>
        &AniStrengthfiber,
    Array<OneD, Array<OneD, NekDouble>> &outarray)
{
    const int nq       = GetTotPoints();
    const int expdim = m_expdim;
    const int numfiber = m_numfiber;

    outarray = Array<OneD, Array<OneD, NekDouble>>(expdim);
    for (int k = 0; k < m_expdim; ++k)
    {
        outarray[k] = Array<OneD, NekDouble>(nq, -2.0);
    }

    NekDouble tmp;
    for (int n = 0; n < numfiber; ++n)
    {
        for (int k = 0; k < m_expdim; ++k)
        {
            for (int i = 0; i < nq; ++i)
            {
                tmp = AniStrengthfiber[n][k][i];
                if (tmp > outarray[k][i])
                {
                    outarray[k][i] = tmp;
                }
            }
        }
    }
}

void MMFNeuralEP::SetUpNeuralCm()
{
    const int numfiber = m_numfiber;
    const int nq = GetTotPoints();    

    m_zoneindexfiber = Array<OneD, Array<OneD, int>>(numfiber);
    for (int n=0; n<numfiber; ++n)
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

    m_extrazone = Array<OneD, NekDouble>(nq) ;
    m_intrazone = Array<OneD, NekDouble>(nq) ;

    Vmath::Vadd(nq, m_nodezone, 1, m_myelinzone, 1, m_intrazone, 1);
    Array<OneD, NekDouble> allone(nq, 1.0);
    Vmath::Vsub(nq, allone, 1, m_myelinzone, 1, m_extrazone, 1);

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
}


void MMFNeuralEP::ComputeRegionalSigma(
    const Array<OneD, const int> &zoneindex,
    Array<OneD, Array<OneD, NekDouble>> &sigma_i,
    Array<OneD, Array<OneD, NekDouble>> &sigma_e,
    Array<OneD, Array<OneD, NekDouble>> &sigma_eM)
{
    const int nq   = GetTotPoints();
    const NekDouble phiefactor = m_phiefactor;

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
                    sigma_eM[j][i] = phiefactor/m_ratio_re_ri;
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


    void MMFNeuralEP::TestHelmSolve()
    {
        const int nq = GetTotPoints();    
    
        Array<OneD, NekDouble> forcing(nq, 0.0);
        for (int i = 0; i < nq; ++i)
        {
            if (m_zoneindexfiber[0][i] >= 0)
            {
                forcing[i] = -100.0;
            }
        }
    
        // Compute phim
        StdRegions::ConstFactorMap phimfactors;
        phimfactors[StdRegions::eFactorTau] = m_Helmtau;
    
        NekDouble lambda;
        m_session->LoadParameter("Helmlambda", lambda, 0.001);
        phimfactors[StdRegions::eFactorLambda] = m_Cn * m_Rf / lambda;
    
        Array<OneD, NekDouble> phim(nq, 0.0);
    
        NekDouble intforcing = AvgInt(forcing);
    
        Vmath::Sadd(nq, -1.0 * intforcing, forcing, 1,
                    m_fields[0]->UpdatePhys(), 1);
        m_fields[0]->HelmSolve(m_fields[0]->GetPhys(),
                               m_fields[0]->UpdateCoeffs(), phimfactors,
                               m_varcoeff);
        m_fields[0]->BwdTrans(m_fields[0]->GetCoeffs(),
                              m_fields[0]->UpdatePhys());
        m_fields[0]->SetPhysState(true);
    
        phim = m_fields[0]->GetPhys();
    
        std::cout << "phim, max = " << Vmath::Vmax(nq, phim, 1)
                  << ", min = " << Vmath::Vmin(nq, phim, 1) << std::endl;
    
        Array<OneD, NekDouble> phimintra(nq);
        Array<OneD, NekDouble> phimextra(nq);
    
        NekDouble phim_max = Vmath::Vmax(nq, phim, 1);
        for (int i = 0; i < nq; ++i)
        {
            phimextra[i] =
                (m_zoneindexfiber[0][i] == -2) ? phim[i] / phim_max : 0.0;
        }
    
        std::cout << "phm in ex_zone: L2err = " << RootMeanSquare(phimextra)
                  << ", Linf = " << Vmath::Vamax(nq, phimextra, 1) << std::endl;
    
        // Compute phie
        StdRegions::ConstFactorMap phiefactors;
        phiefactors[StdRegions::eFactorTau]    = m_Helmtau;
        phiefactors[StdRegions::eFactorLambda] = 0.0;
    
        Array<OneD, NekDouble> phie(nq, 0.0);
    
        Vmath::Sadd(nq, -1.0 * intforcing, forcing, 1,
                    m_fields[1]->UpdatePhys(), 1);
        m_fields[1]->HelmSolve(m_fields[1]->GetPhys(),
                               m_fields[1]->UpdateCoeffs(), phiefactors,
                               m_phievarcoeff);
        m_fields[1]->BwdTrans(m_fields[1]->GetCoeffs(),
                              m_fields[1]->UpdatePhys());
        m_fields[1]->SetPhysState(true);
    
        phie = m_fields[1]->GetPhys();
    
        std::cout << "phie, max = " << Vmath::Vmax(nq, phie, 1)
                  << ", min = " << Vmath::Vmin(nq, phie, 1) << std::endl;
    
        // Plotting the result
        int nvar    = 8;
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
    
        for (int j = 0; j < m_spacedim; ++j)
        {
            mfmag[j]     = Array<OneD, NekDouble>(nq, 0.0);
            phiemfmag[j] = Array<OneD, NekDouble>(nq, 0.0);
    
            for (int k = 0; k < m_spacedim; ++k)
            {
                Vmath::Vvtvp(nq, &m_movingframes[j][k * nq], 1,
                             &m_movingframes[j][k * nq], 1, &mfmag[j][0], 1,
                             &mfmag[j][0], 1);
                Vmath::Vvtvp(nq, &m_phiemovingframes[j][k * nq], 1,
                             &m_phiemovingframes[j][k * nq], 1,
                             &phiemfmag[j][0], 1, &phiemfmag[j][0], 1);
            }
    
            Vmath::Vsqrt(nq, &mfmag[j][0], 1, &mfmag[j][0], 1);
            Vmath::Vsqrt(nq, &phiemfmag[j][0], 1, &mfmag[j][0], 1);
        }

    }

void MMFNeuralEP::DoOdeProjection(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time)
{
    int nvariables = inarray.size();
    SetBoundaryConditions(time);

    switch (m_projectionType)
    {
        case MultiRegions::eDiscontinuous:
        {
            // Just copy over array
            int npoints = GetNpoints();

            for (int i = 0; i < nvariables; ++i)
            {
                Vmath::Vcopy(npoints, inarray[i], 1, outarray[i], 1);
            }
            break;
        }
        case MultiRegions::eGalerkin:
        case MultiRegions::eMixed_CG_Discontinuous:
        {
            Array<OneD, NekDouble> coeffs(m_fields[0]->GetNcoeffs());

            for (int i = 0; i < nvariables; ++i)
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
        const Array<OneD, const NekDouble> &fiberleft,
        const Array<OneD, const NekDouble> &fiberright,
        const Array<OneD, const int> &fiberorder,
        Array<OneD, Array<OneD, int>> &zoneindexfiber)
    {
        const int nq   = GetTotPoints();
        const int numfiber = m_numfiber;

        const Array<OneD, NekDouble> &xcell = m_xcell;
        const Array<OneD, NekDouble> &ycell = m_ycell;

        int index, cnt;
        NekDouble xi, yi;
        zoneindexfiber = Array<OneD, Array<OneD, int>>(numfiber);
        for (int n=0; n<numfiber; ++n)
         {
            zoneindexfiber[n]  = Array<OneD, int>(nq, -2); 

            cnt = 0;
            for (int i=0; i<nq; ++i)
            {
                xi = xcell[i];
                yi = ycell[i];

                if( (m_FiberType==eLinearAligned) || (m_FiberType==eLinearMisAligned) )
                {
                    if((xi>fiberleft[n]) && (xi<fiberright[n]))
                        {
                            index = FiberIndex(m_FiberType, n, fiberorder[n], xi, yi);
                            if(index >-2)
                            {
                                zoneindexfiber[n][i] = index;
                                cnt++;
                            }
                        }
                }

                else{
                    zoneindexfiber[n][i] = FiberIndex(m_FiberType, n, fiberorder[n], xi, yi);
                }
            }
            std::cout << "for firber n = " << n << ", intracellular space = " << cnt << "/" << nq << std::endl;
        }

        // If all node zone = index = 1 for myelin. 
        if(m_MediumType==eAllNode)
        {
            for (int n=0; n<numfiber; ++n)
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
    }

Array<OneD, int> MMFNeuralEP::TestRanvierSingleIndex()
{
    int nq   = GetTotPoints();

    // Compute the center of the element as the coordinate of each grid.
    const Array<OneD, NekDouble> &xcell = m_xcell;
    const Array<OneD, NekDouble> &ycell = m_ycell;

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
    const Array<OneD, NekDouble> &xcell = m_xcell;
    const Array<OneD, NekDouble> &ycell = m_ycell;
\
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
    const int fiberorder, const NekDouble xi, const NekDouble yi)
    {
        int index=0;

        switch(FiberType)
        {
            case eLinearAligned:
            {
                index = LinearAlignedFiberIndex(fiberorder, yi);                
                break;
            }

            case eLinearMisAligned:
            {
                index = LinearMisAlignedFiberIndex(fibern, yi);                
                break;
            }

            case eLinearDivergent:
            {
                index = LinearDivergentFiberIndex(fibern, xi, yi);                
                break;
            }

            case eLinearCrossing:
            {
                index = LinearCrossingFiberIndex(fibern, fiberorder, xi, yi);                
                break;
            }

            case eConstantCurved:
            {
                index = ConstantCurvedFiberIndex(fibern, m_fibercurvature, xi, yi);                
                break;
            }

            default:
            break;
        }

        return index;
    }


int MMFNeuralEP::ConstantCurvedFiberIndex(const int fibern, 
    const NekDouble fibercurvature, const NekDouble xi, const NekDouble yi)
{
    int output = -2;

  //  const int numfiber = m_numfiber;
    const NekDouble totNnode = m_totNode;
    const NekDouble nodelen = m_nodelen;
    const NekDouble myelinlen = m_myelinlen;
  //  const NekDouble nodeinitdown = m_nodeinitdown;
  //  const NekDouble nodeinitup = m_nodeinitup;

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
    const int fiberorder, const NekDouble yi)
{
    int output;

    const NekDouble totNnode = m_totNode;
    const NekDouble nodelen = m_nodelen;
    const NekDouble myelinlen = m_myelinlen;
    const NekDouble nodeinitdown = m_nodeinitdown;
    const NekDouble nodeinitup = m_nodeinitup;
    
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


int MMFNeuralEP::LinearMisAlignedFiberIndex(const int fibern, const NekDouble yi)
{
    int output = -2;
    
    const NekDouble totNnode = m_totNode;
    const NekDouble nodelen = m_nodelen;
    const NekDouble myelinlen = m_myelinlen;
    const NekDouble nodeinitdown = m_nodeinitdown;
    const NekDouble nodeinitup = m_nodeinitup;

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

int MMFNeuralEP::LinearCrossingFiberIndex(
    const int fibern, const int fiberorder, const NekDouble xi, const NekDouble yi)
{
    int output;

    const NekDouble totNnode = m_totNode;
    const NekDouble nodelen = m_nodelen;
    const NekDouble myelinlen = m_myelinlen;
    const NekDouble nodeinitdown = m_nodeinitdown;
    const NekDouble nodeinitup = m_nodeinitup;
    const NekDouble fiberlength = m_fiberlength;
    
    output = -2 ;  // Default of the first fiber = myelin

    // Excitezone
    if(fibern==0)
    {
        if((xi>m_fiberleft[fibern]) && (xi<m_fiberright[fibern]))
        {
            output = LinearAlignedFiberIndex(fiberorder, yi); 
        }
    }

    else if(fibern==1)
    {
        NekDouble nodestart, nodeend;

        const NekDouble beta = 0.5 * m_pi - m_fiberangle;
        const NekDouble tanb = tan(beta);
        const NekDouble cosb = cos(beta);
        const NekDouble sinb = sin(beta);

        const NekDouble gap = nodelen / cosb;
        const NekDouble upperline = tanb * xi + 0.5 * fiberlength - yi;
        const NekDouble lowerline = tanb * xi + 0.5 * fiberlength - gap - yi;

        const NekDouble x0 = -0.5*fiberlength/sqrt(1 + tanb*tanb);
        const NekDouble y0 = tanb * x0 + 0.5 * fiberlength;
        const NekDouble sp0 = x0 * cosb + y0 * sinb;
        const NekDouble sp = xi * cosb + yi * sinb;
        const NekDouble dist = sp - sp0;

        if( (upperline*lowerline < 0) && ( dist < fiberlength) )
        {
            output = -1;

            // Excitezone
            if( (dist >= nodeinitdown) && (dist <= nodeinitup) )
            {
                output = 0;
            }

            for (int k=0; k<totNnode; ++k)
            {
                nodestart = nodeinitup + myelinlen + k * (myelinlen + nodelen);
                nodeend = nodeinitup + (k+1) * (myelinlen + nodelen);
                if( (dist>=nodestart) && (dist<=nodeend) )
                {
                    output = k + 1;
                }
            }

            // Before first node or after last node
            const NekDouble nodeend = nodeinitup + totNnode * (myelinlen + nodelen);
            if (dist < nodeinitdown || dist > nodeend)
            {
                output = -2;
            }
        }
    }

    return output;
}


int MMFNeuralEP::LinearDivergentFiberIndex(const int fibern, const NekDouble xi, const NekDouble yi)
{
    const NekDouble totNnode = m_totNode;
    const NekDouble nodelen = m_nodelen;
    const NekDouble myelinlen = m_myelinlen;
    const NekDouble nodeinitdown = m_nodeinitdown;
    const NekDouble nodeinitup = m_nodeinitup;

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

void MMFNeuralEP::GetCellCoordAvg(
    Array<OneD, NekDouble> &xcell, 
    Array<OneD, NekDouble> &ycell, 
    Array<OneD, NekDouble> &zcell)
{
    const int nq   = GetTotPoints();
    const int npts = m_fields[0]->GetTotPoints(0);

    const Array<OneD, NekDouble> &x0 = m_x;
    const Array<OneD, NekDouble> &x1 = m_y;
    const Array<OneD, NekDouble> &x2 = m_z;

    int Nelem = nq/npts;

    Array<OneD, NekDouble> xcellavg(Nelem,0.0);
    Array<OneD, NekDouble> ycellavg(Nelem,0.0);
    Array<OneD, NekDouble> zcellavg(Nelem,0.0);

    int index;
    for (int i=0; i<nq; ++i)
    {
        index = i/npts;
        xcellavg[index] = xcellavg[index] + x0[i];
        ycellavg[index] = ycellavg[index] + x1[i];
        zcellavg[index] = zcellavg[index] + x2[i];
    }

    Vmath::Smul(Nelem, 1.0/npts, xcellavg, 1, xcellavg, 1);
    Vmath::Smul(Nelem, 1.0/npts, ycellavg, 1, ycellavg, 1);
    Vmath::Smul(Nelem, 1.0/npts, zcellavg, 1, zcellavg, 1);

    for (int i=0; i<nq; ++i)
    {
        index = i/npts;

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


void MMFNeuralEP::SetUpDomainZone()
{
    const int nq   = GetTotPoints();
    const int numfiber = m_numfiber;

    int index;
    // Construction ZoneIndex for all fibers;
    m_zoneindex = Array<OneD, int>(nq, -2);
    for (int i=0; i<nq; ++i)
    {
        for (int n=0; n<numfiber; ++n)
        {
            index = m_zoneindexfiber[n][i];
            if(index > m_zoneindex[i])
            {
                m_zoneindex[i] = index;
            }
        }
    }

    m_excitezonefiber = Array<OneD, Array<OneD, NekDouble>>(numfiber);
    m_intrazonefiber = Array<OneD, Array<OneD, NekDouble>>(numfiber);
    m_nodezonefiber = Array<OneD, Array<OneD, NekDouble>>(numfiber);
    for (int n=0; n<numfiber; ++n)
    {
        m_excitezonefiber[n] = Array<OneD, NekDouble>(nq, 0.0);
        m_intrazonefiber[n] = Array<OneD, NekDouble>(nq, 0.0);
        m_nodezonefiber[n] = Array<OneD, NekDouble>(nq, 0.0);
    }

    // Set up the total node zone and intra zone
    m_nodezone = Array<OneD, NekDouble>(nq, 0.0);  
    m_myelinzone = Array<OneD, NekDouble>(nq, 0.0);              
    for (int i=0; i<nq; ++i)
    {
        for (int n=0; n<numfiber; ++n)
        {
            index = m_zoneindexfiber[n][i];

            // excite zone
            if(index == 0)
            {
                m_excitezonefiber[n][i] = 1.0;
            }

            // intra zone = myelin or nodezone 
            if( index > -2) 
            {
                m_intrazonefiber[n][i] = 1.0;
            }

            // node zone
            if (index > -1)
            {
                 m_nodezonefiber[n][i] = 1.0;
                 m_nodezone[i] = 1.0;
            }

            // myelin zone
            if (index == -1)
            {
                m_myelinzone[i] = 1.0;
            }
        }
    }

    Array<OneD, NekDouble> allone(nq, 1.0);

    m_extrazone = Array<OneD, NekDouble>(nq) ;
    m_intrazone = Array<OneD, NekDouble>(nq) ;
    m_outerzone = Array<OneD, NekDouble>(nq) ;

    Vmath::Vadd(nq, m_nodezone, 1, m_myelinzone, 1, m_intrazone, 1);
    Vmath::Vsub(nq, allone, 1, m_myelinzone, 1, m_extrazone, 1);
    Vmath::Vsub(nq, allone, 1, m_intrazone, 1, m_outerzone, 1);
}

// void MMFNeuralEP::SetUpDomainZone(
//         const int numfiber,
//         const Array<OneD, const Array<OneD, int>> &zoneindexfiber)
//         Array<OneD, Array<OneD, NekDouble>> &excitezonefiber,
//         Array<OneD, Array<OneD, NekDouble>> &intrazonefiber,
//         Array<OneD, int> &zoneindex,
//         Array<OneD, NekDouble> &nodezone,
//         Array<OneD, NekDouble> &myelinzone,
//         Array<OneD, NekDouble> &intrazone,
//         Array<OneD, NekDouble> &extrazone,
//         Array<OneD, NekDouble> &outerzone)
// {
//     const int nq   = GetTotPoints();
//     int index;

//     // Construction ZoneIndex for all fibers;
//     zoneindex = Array<OneD, int>(nq, -2);
//     for (int i=0; i<nq; ++i)
//     {
//         for (int n=0; n<numfiber; ++n)
//         {
//             index = zoneindexfiber[n][i];
//             if(index > zoneindex[i])
//             {
//                 zoneindex[i] = index;
//             }
//         }
//     }

//     excitezonefiber = Array<OneD, Array<OneD, NekDouble>>(numfiber);
//     intrazonefiber = Array<OneD, Array<OneD, NekDouble>>(numfiber);
//     for (int n=0; n<numfiber; ++n)
//     {
//         excitezonefiber[n] = Array<OneD, NekDouble>(nq, 0.0);
//         intrazonefiber[n] = Array<OneD, NekDouble>(nq, 0.0);
//     }

//     // Set up the total node zone and intra zone
//     nodezone = Array<OneD, NekDouble>(nq, 0.0);  
//     myelinzone = Array<OneD, NekDouble>(nq, 0.0);              
//     for (int i=0; i<nq; ++i)
//     {
//         for (int n=0; n<numfiber; ++n)
//         {
//             index = zoneindexfiber[n][i];

//             // excite zone
//             if(index == 0)
//             {
//                 excitezonefiber[n][i] = 1.0;
//             }

//             // intra zone = myelin or nodezone 
//             if( index > -2) 
//             {
//                 intrazonefiber[n][i] = 1.0;
//             }

//             // node zone
//             if (index > -1)
//             {
//                  nodezone[i] = 1.0;
//             }

//             // myelin zone
//             if (index == -1)
//             {
//                  myelinzone[i] = 1.0;
//             }
//         }
//     }

//     extrazone = Array<OneD, NekDouble>(nq) ;
//     intrazone = Array<OneD, NekDouble>(nq) ;
//     outerzone = Array<OneD, NekDouble>(nq) ;

//     Array<OneD, NekDouble> allone(nq, 1.0);

//     Vmath::Vadd(nq, nodezone, 1, myelinzone, 1, intrazone, 1);
//     Vmath::Vsub(nq, allone, 1, myelinzone, 1, extrazone, 1);
//     Vmath::Vsub(nq, allone, 1, intrazone, 1, outerzone, 1);
// }

void MMFNeuralEP::PlotDomainZone()
{
        switch(m_numfiber)
        {
            case 1:
            {
                PlotDomainZonefib1(m_zoneindex, m_intrazone, m_extrazone, m_outerzone);
                break;
            }

            case 2:
            {
                PlotDomainZonefib2(m_zoneindexfiber, m_zoneindex, m_intrazonefiber, m_intrazone, m_extrazone, m_outerzone);
                break;
            }

            default:
            break;
        }
}

// Plotting Domain Zone
void MMFNeuralEP::PlotDomainZonefib2(
        const Array<OneD, const Array<OneD, int>> zoneindexfiber,
        const Array<OneD, const int> &zoneindex,
        const Array<OneD, const Array<OneD, NekDouble>> &intrazonefiber,
        const Array<OneD, const NekDouble> &intrazone,
        const Array<OneD, const NekDouble> &extrazone,
        const Array<OneD, const NekDouble> &outerzone)
{
    const int nq   = GetTotPoints();
    const int nvar    = 2 * m_numfiber + 4;
    const int ncoeffs = m_fields[0]->GetNcoeffs();

    std::string outname1 = m_sessionName + "_zone.chk";

    std::vector<std::string> variables(nvar);

    variables[0] = "zoneindexfiber1";
    variables[1] = "zoneindexfiber2";
    variables[2] = "zoneindex";
    variables[3] = "intrazonefiber1";
    variables[4] = "intrazonefiber2";
    variables[5] = "intrazone";
    variables[6] = "extrazone";
    variables[7] = "outerzone";

    std::vector<Array<OneD, NekDouble>> fieldcoeffs(nvar);
    for (int i = 0; i < nvar; ++i)
    {
        fieldcoeffs[i] = Array<OneD, NekDouble>(ncoeffs);
    }

    Array<OneD, NekDouble> zoneindextmp(nq);
    for (int nfib=0; nfib<2; ++nfib)
    {
        for (int j=0; j< nq; ++j)
        {
            zoneindextmp[j] = 1.0 * zoneindexfiber[nfib][j];
        }

        m_fields[0]->FwdTransLocalElmt(zoneindextmp, fieldcoeffs[nfib]);
    }

    for (int j=0; j< nq; ++j)
    {
        zoneindextmp[j] = 1.0 * zoneindex[j];
    }

    m_fields[0]->FwdTransLocalElmt(zoneindextmp, fieldcoeffs[2]);

    for (int nfib=0; nfib<2; ++nfib)
    {
        m_fields[0]->FwdTransLocalElmt(intrazonefiber[nfib], fieldcoeffs[3+nfib]);
    }

    m_fields[0]->FwdTransLocalElmt(intrazone, fieldcoeffs[5]);
    m_fields[0]->FwdTransLocalElmt(extrazone, fieldcoeffs[6]);
    m_fields[0]->FwdTransLocalElmt(outerzone, fieldcoeffs[7]);

    WriteFld(outname1, m_fields[0], fieldcoeffs, variables);
}

// Plotting Domain Zone
void MMFNeuralEP::PlotDomainZonefib1(
        const Array<OneD, const int> &zoneindex,
        const Array<OneD, const NekDouble> &intrazone,
        const Array<OneD, const NekDouble> &extrazone,
        const Array<OneD, const NekDouble> &outerzone)
{
    const int nq   = GetTotPoints();
    const int nvar    = 4;
    const int ncoeffs = m_fields[0]->GetNcoeffs();

    std::string outname1 = m_sessionName + "_zone.chk";

    std::vector<std::string> variables(nvar);

    variables[0] = "zoneindex";
    variables[1] = "intrazone";
    variables[2] = "extrazone";
    variables[3] = "outerzone";

    std::vector<Array<OneD, NekDouble>> fieldcoeffs(nvar);
    for (int i = 0; i < nvar; ++i)
    {
        fieldcoeffs[i] = Array<OneD, NekDouble>(ncoeffs);
    }

    Array<OneD, NekDouble> zoneindextmp(nq);
    for (int j=0; j< nq; ++j)
    {
        zoneindextmp[j] = 1.0 * zoneindex[j];
    }

    m_fields[0]->FwdTransLocalElmt(zoneindextmp, fieldcoeffs[0]);
    m_fields[0]->FwdTransLocalElmt(intrazone, fieldcoeffs[1]);
    m_fields[0]->FwdTransLocalElmt(extrazone, fieldcoeffs[2]);
    m_fields[0]->FwdTransLocalElmt(outerzone, fieldcoeffs[3]);

    WriteFld(outname1, m_fields[0], fieldcoeffs, variables);
}

Array<OneD, NekDouble> MMFNeuralEP::ComputeConductivity(
                 const Array<OneD, const int> &zoneindex)
{
    const int nq   = GetTotPoints();
    const int npts = m_fields[0]->GetTotPoints(0);

    const NekDouble Cm = m_neuron->GetCapacitanceValue(0);
    const NekDouble Cn = m_neuron->GetCapacitanceValue(1);

    Array<OneD, NekDouble> outarray(nq);

    int cntn = 0, cntm = 0, cnte = 0;
    for (int i = 0; i < nq; ++i)
    {
        // Ranvier node zone
        if ( (zoneindex[i] >= 0) && (zoneindex[i] < 100) )

        {
            outarray[i] = 1.0 / Cn;
            cntn++;
        }

        if ( (zoneindex[i] >= 100) && (zoneindex[i] < 200) )

        {
            outarray[i] = 1.0 / Cn;
            cntn++;
        }

        // Myelin zone
        if (zoneindex[i] == -1) 
        {
            outarray[i] = 1.0 / Cm;
            cntm++;
        }

        // Myelin zone
        if (zoneindex[i] == -101) 
        {
            outarray[i] = 1.0 / Cm;
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

    std::cout << "ComputeConductivity: Node = " << cntn/npts << ", Myelinf1 = " 
    << cntm/npts << ", extracell = " << cnte/npts << std::endl;

    return outarray;
}

Array<OneD, NekDouble> MMFNeuralEP::ComputeConductivity(
                 const Array<OneD, const Array<OneD, int>> &zoneindexfiber)
{
    const int nq   = GetTotPoints();
    const int numfiber = m_numfiber;

    const int npts = m_fields[0]->GetTotPoints(0);

    const NekDouble Cm = m_neuron->GetCapacitanceValue(0);
    const NekDouble Cn = m_neuron->GetCapacitanceValue(1);

    Array<OneD, NekDouble> outarray(nq, 0.0); // NeuralCm is zero at indexzone == -2

    int index, cntm = 0, cntn = 0, cnte = 0;
    for (int i = 0; i < nq; ++i)
    {
        for (int n = 0; n < numfiber; ++n)
            {
                index = zoneindexfiber[n][i];

                // Ranvier node zone for all the fibers
                if (index >= 0)
                {
                    outarray[i] = 1.0 / Cn;
                    cntn++;
                }

                // Ranvier node zone for all the fibers
                if (index == -1) 
                {
                    if (m_MediumType == eAnisotropy)
                    {  
                        outarray[i] = 1.0 / Cm;
                    }

                    else
                    {
                        outarray[i] = 1.0 / Cn;
                    }
                    cntm++;
                }
            }
    }

    cnte = nq - cntn - cntm;

    std::cout << "ComputeConductivity: Node = " << cntn/npts << ", Myelinf1 = " 
    << cntm/npts << ", extracell = " << cnte/npts << std::endl;

    return outarray;
}

void MMFNeuralEP::CheckNodeZoneMF(
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    const Array<OneD, const Array<OneD, NekDouble>> &phiemovingframes)
{
    int nq = GetTotPoints();

    const Array<OneD, NekDouble> &x0 = m_x;
    const Array<OneD, NekDouble> &x1 = m_y;

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

        const Array<OneD, NekDouble> &x0 = m_x;
        const Array<OneD, NekDouble> &x1 = m_y;
        const Array<OneD, NekDouble> &x2 = m_z;

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

// void MMFNeuralEP::v_DoSolve()
// {
//     switch (m_SolverSchemeType)
//     {
//         case eMMFZero:
//         case eMMFFirst:
//         case eTimeMap:
//         {
//             DoSolveMMF();
//             break;
//         }

//         default:
//          break;
//     }
// }

void MMFNeuralEP::v_DoSolve()
{
    ASSERTL0(m_intScheme != 0, "No time integration scheme.");

    // int i, nchk = 1;
    const int nq      = GetTotPoints();
    const int nvar    = m_fields.size();
    const int phievar = nvar - 1;

    // const int totsteps = (m_steps + 1) / m_checksteps;

    int step = 0, nchk = 1;

    NekDouble intTime = 0.0, cpuTime = 0.0;

    // const int nvariables = m_intVariables.empty() ? nvar : m_intVariables.size();
    if (m_intVariables.empty())
    {
        for (int i = 0; i < nvar; ++i)
        {
            m_intVariables.push_back(i);
        }
    }

    // Set up working arrays
    Array<OneD, Array<OneD, NekDouble>> fields(nvar), fields_old(nvar);
    Array<OneD, Array<OneD, NekDouble>> dphidt(nvar), dphidtint(nvar);
    Array<OneD, Array<OneD, NekDouble>> TimeMap(nvar);

    for (int i = 0; i < nvar; ++i)
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
    // Array<OneD, NekDouble> timevec(totsteps, 0.0);
    // Array<OneD, NekDouble> thredlocf1(totsteps, 0.0), thredlocf2(totsteps, 0.0);
    // Array<OneD, int> thredlocf1zone(totsteps, 0), thredlocf2zone(totsteps, 0);

    LibUtilities::Timer timer;
    while (step < m_steps || m_time < m_fintime - NekConstants::kNekZeroTol)
    {
        // Save current solution
        for (int n=0; n < nvar; ++n)
        {
            Vmath::Vcopy(nq, &fields[n][0], 1, &fields_old[n][0], 1);
        }

        // Time integration
        timer.Start();
        fields = m_intScheme->TimeIntegrate(step, m_timestep);
        timer.Stop();

        m_time += m_timestep;
        NekDouble elapsed = timer.TimePerTest(1);
        intTime += elapsed;
        cpuTime += elapsed;

        // Compute normalized time derivatives
        fields[phievar] = m_fields[phievar]->GetPhys();
        NekDouble factor;
        for (int n = 0; n < nvar; ++n)
        {
            NekDouble maxphi = Vmath::Vamax(nq, fields[n], 1);
            factor = 1.0 / (m_timestep * maxphi);
            Vmath::Vsub(nq, fields[n], 1, fields_old[n], 1, dphidt[n], 1);
            Vmath::Smul(nq, factor, dphidt[n], 1, dphidt[n], 1);
        }

        // Compute neural time map
        ComputephimTimeMap(m_time, fields[0], dphidt[0], dphidtint[0], TimeMap[0]);
        for (int n = 1; n < nvar; ++n)
        {
            ComputephieTimeMap(m_time, fields[n], dphidt[n], dphidtint[n], TimeMap[n]);
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
            // Create .chk files for plotting
            PlotNeuralEP(fields, m_TimeMap, nchk);

            // Print out the values at the nodes
            PrintAtNodes(nvar, m_numfiber, fields);

            // Write out checkpoint files
            Checkpoint_Output(nchk++);
        }

        ++step;
    } // namespace Nektar

    // Print out summary statistics
    if (m_session->GetComm()->GetRank() == 0)
    {
        std::cout << "Time-integration complete. Total CPU time: " << intTime << "s\n";
    }

    for (int i = 0; i < nvar; ++i)
    {
        m_fields[m_intVariables[i]]->SetPhys(fields[i]);
        m_fields[m_intVariables[i]]->SetPhysState(true);

        m_fields[m_intVariables[i]]->FwdTrans(m_fields[i]->GetPhys(),
                                   m_fields[m_intVariables[i]]->UpdateCoeffs());
    }

    // std::cout << " timevec: ";
    // for (int i=0; i<totsteps; ++i)
    // {
    //     std::cout << timevec[i] << ", ";
    // }
    // std::cout << std::endl;

    // std::cout << " CSDvecatnode1: ";
    // for (int i=0; i<totsteps; ++i)
    // {
    //     std::cout << CSDvecatnode1[i] << ", ";
    // }
    // std::cout << std::endl;

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
} 

void MMFNeuralEP::PrintAtNodes(const int nvar, const int numfiber,
                               const Array<OneD, const Array<OneD, NekDouble>> &fields)
{
    const int totNode = m_totNode;

    NekDouble phim1, phim2, phie;
    for (int n = 0; n<totNode; ++n)
    {
        if(numfiber>1)
        {
            if(nvar==2)
            {
                phim1 = DisplayAtNodes(0, n, m_zoneindexfiber, fields[0]);
                phim2 = DisplayAtNodes(1, n, m_zoneindexfiber, fields[0]);
            }

            else if(nvar==3)
            {
                phim1 = DisplayAtNodes(0, n, m_zoneindexfiber, fields[0]);
                phim2 = DisplayAtNodes(1, n, m_zoneindexfiber, fields[1]);
            }

            phie = DisplayAtNodes(0, n, m_zoneindexfiber, fields[1]);

            std::cout << "At node n = " << n << ", phim1 = " << phim1 
            << ", phim2 = " << phim2
            << ", phie = " << phie << std::endl;
        }

        else
        {
            phim1 = DisplayAtNodes(0, n, m_zoneindexfiber, fields[0]);
            phie = DisplayAtNodes(0, n, m_zoneindexfiber, fields[1]);

            std::cout << "At node n = " << n << ", phim = " << phim1 
            << ", phie = " << phie << std::endl;
        }

    }
}

// namespace Nektar
NekDouble MMFNeuralEP::DisplayAtNodes(const int fibern, const int nodeindex, 
        const Array<OneD, const Array<OneD, int>> &zoneindexfiber,
        const Array<OneD, const NekDouble> &inarray)
{
    const int nq = GetTotPoints();
    int cnt=0;

    NekDouble output=0.0;

    for (int i = 0; i < nq; ++i)
    {
        if( zoneindexfiber[fibern][i] == nodeindex )
        {
            output += inarray[i];
            cnt++;
        }
    }

    output = output / cnt;

    return output;
}


void MMFNeuralEP::ComputephimTimeMap(const NekDouble time,
                                    const Array<OneD, const NekDouble> &field,
                                    const Array<OneD, const NekDouble> &dphidt,
                                    Array<OneD, NekDouble> &dphidtint,
                                    Array<OneD, NekDouble> &TimeMap)
{
    const int nq = GetTotPoints();

    const NekDouble phiTol = 10.0;
    const NekDouble phirest = 80.0;
    const NekDouble dphidtTol = 1.0;

    for (int i = 0; i < nq; ++i)
    {
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
        const Array<OneD, const NekDouble> &dphidt,
        Array<OneD, NekDouble> &dphidtint,
        Array<OneD, NekDouble> &TimeMap)
{
    const int nq = GetTotPoints();
    const NekDouble dphidtTol = 0.1;

    NekDouble phie;
    for (int i = 0; i < nq; ++i)
    {
        const NekDouble dphi = dphidt[i];
        phie = field[i];

        // if the field is rising and is largest by now, time mep is the corresponding time.
        if ( (dphi > dphidtTol) && (phie > dphidtint[i]) )
        {
            TimeMap[i] = time;
            dphidtint[i] = phie;
        }
    }
}


// void MMFNeuralEP::ComputephieTimeMap(
//     const NekDouble time,
//     const Array<OneD, const NekDouble> &field,
//     Array<OneD, NekDouble> &fieldint,
//     Array<OneD, NekDouble> &TimeMap)
// {
//     const int nq = GetTotPoints();
//     constexpr NekDouble Tol = 0.01;

//     for (int i = 0; i < nq; ++i)
//     {
//         const NekDouble phie = field[i] + 1.0;
//         if (phie > Tol)
//         {
//             const NekDouble fint = fieldint[i];
//             const NekDouble fnewsum = phie + fint;

//             TimeMap[i] = (phie * time + fint * TimeMap[i]) / fnewsum;
//             fieldint[i] += phie;
//         }
//     }
// }

// void MMFNeuralEP::ComputerhoTimeMap(
//     const NekDouble time,
//     const Array<OneD, const NekDouble> &field,
//     Array<OneD, NekDouble> &fieldint,
//     Array<OneD, NekDouble> &TimeMap)
// {
//     const int nq = GetTotPoints();
//     constexpr NekDouble Tol = 0.000001;

//     #pragma omp parallel for
//     for (int i = 0; i < nq; ++i)
//     {
//         const NekDouble rho = field[i] + 0.002;
//         if (rho > Tol)
//         {
//             const NekDouble fint = fieldint[i];
//             const NekDouble fnewsum = rho + fint;

//             TimeMap[i] = (rho * time + fint * TimeMap[i]) / fnewsum;
//             fieldint[i] += rho;
//         }
//     }
// }

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
        const int nvar = m_fields.size();
        switch(nvar)
        {
            case 2:
            {
                PlotNeuralEPvar2(fields, TimeMap, nstep);
                break;
            }

            case 3:
            {
                if( (m_NeuralEPType == eNeuralEP2DbiMulti) || (m_NeuralEPType == eNeuralEP2DbiMultiv2) )
                {
                    PlotNeuralEPvar3(fields, TimeMap, nstep);
                }
                else if(m_NeuralEPType == eNeuralEP2DbiCSD)
                {
                    PlotNeuralEPvar3CSD(fields, TimeMap, nstep);
                }
                break;
            }

            case 4:
            {
                PlotNeuralEPvar4(fields, TimeMap, nstep);
                break;
            }

            default:
            break;
        }
    }


void MMFNeuralEP::PlotNeuralEPvar2(
    const Array<OneD, const Array<OneD, NekDouble>> &fields,
    const Array<OneD, const Array<OneD, NekDouble>> &TimeMap,
    const int nstep)
{
    const int nvar    = 5;
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
    variables[3] = "TimeMap_phim";
    variables[4] = "TimeMap_phie";

    //     variables[0] = "phi_m";
    Array<OneD, NekDouble> tmp(nq);
    m_fields[0]->FwdTransLocalElmt(fields[0], fieldcoeffs[0]);

    Array<OneD, NekDouble> phim(nq), phie(nq);
    Vmath::Vmul(nq, m_intrazone, 1, fields[0], 1, phim, 1);
    Vmath::Vmul(nq, m_outerzone, 1, fields[1], 1, phie, 1);

    //     variables[1] = "phi_e";
    m_fields[0]->FwdTransLocalElmt(phie, fieldcoeffs[1]);

    //     variables[2] = "CSD";
    Array<OneD, NekDouble> CSD = ComputeMMFDiffusion(m_phiemovingframes, fields[1]);
    m_fields[0]->FwdTransLocalElmt(CSD, fieldcoeffs[2]);

    // Max values and indices
    const NekDouble Maxphim  = Vmath::Vmax(nq, phim, 1);
    const int       Maxphimid = Vmath::Imax(nq, phim, 1);

    const NekDouble Maxphie  = Vmath::Vmax(nq, phie, 1);
    const int       Maxphieid = Vmath::Imax(nq, phie, 1);

    const NekDouble MaxCSD   = Vmath::Vmax(nq, CSD, 1);
    const int       MaxCSDid  = Vmath::Imax(nq, CSD, 1);

    // Coordinates
    const Array<OneD, NekDouble> &x0 = m_x;
    const Array<OneD, NekDouble> &x1 = m_y;

    std::cout << "phim: Max = " << Maxphim << " at x = " << x0[Maxphimid] << ", y = " << x1[Maxphimid] << '\n';
    std::cout << "phie: Max = " << Maxphie << " at x = " << x0[Maxphieid] << ", y = " << x1[Maxphieid] << '\n';
    std::cout << "CSD: Max = " << MaxCSD << " at x = " << x0[MaxCSDid] << ", y = " << x1[MaxCSDid] << '\n';

    Vmath::Vmul(nq, m_intrazone, 1, TimeMap[0], 1, phim, 1);
    Vmath::Vmul(nq, m_outerzone, 1, TimeMap[1], 1, phie, 1);

    // TimeMap maxima
    std::cout << "TimeMap: phim = " << Vmath::Vmax(nq, phim, 1)
              << ", phie = " << Vmath::Vmax(nq, phie, 1) << '\n';
    
    // variables[4] = "TimeMap_phim";
    m_fields[0]->FwdTransLocalElmt(phim, fieldcoeffs[3]);

    // variables[5] = "TimeMap_phie";
    m_fields[0]->FwdTransLocalElmt(phie, fieldcoeffs[4]);

    WriteFld(outname1, m_fields[0], fieldcoeffs, variables);
}

void MMFNeuralEP::PlotNeuralEPvar3(
    const Array<OneD, const Array<OneD, NekDouble>> &fields,
    const Array<OneD, const Array<OneD, NekDouble>> &TimeMap,
    const int nstep)
{
    const int nvar    = 9;
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
    variables[0] = "field";
    variables[1] = "phi_m";
    variables[2] = "phi_e";
    variables[3] = "phi_m1";
    variables[4] = "phi_m2";
    variables[5] = "TimeMap_phim";
    variables[6] = "TimeMap_phim1";
    variables[7] = "TimeMap_phim2";
    variables[8] = "TimeMap_phie";

    Array<OneD, NekDouble> phim(nq);
    Vmath::Vadd(nq, fields[0], 1, fields[1], 1, phim, 1);
    Vmath::Vmul(nq, m_intrazone, 1, phim, 1, phim, 1);

    Array<OneD, NekDouble> phie(nq);
    Vmath::Vmul(nq, m_outerzone, 1, fields[2], 1, phie, 1);

    m_fields[0]->FwdTransLocalElmt(phim, fieldcoeffs[1]);
    m_fields[0]->FwdTransLocalElmt(phie, fieldcoeffs[2]);

    Array<OneD, NekDouble> totfield(nq);
    Vmath::Smul(nq, 100.0 , phie, 1, phie, 1);
    Vmath::Vadd(nq, phim, 1, phie, 1, totfield, 1);
    m_fields[0]->FwdTransLocalElmt(totfield, fieldcoeffs[0]);


    m_fields[0]->FwdTransLocalElmt(fields[0], fieldcoeffs[3]);
    m_fields[0]->FwdTransLocalElmt(fields[1], fieldcoeffs[4]);

    // Max values and indices
    const NekDouble Maxphim1  = Vmath::Vmax(nq, fields[0], 1);
    const int       Maxphim1id = Vmath::Imax(nq, fields[0], 1);

    const NekDouble Maxphim2  = Vmath::Vmax(nq, fields[1], 1);
    const int       Maxphim2id = Vmath::Imax(nq, fields[1], 1);

    const NekDouble Maxphie  = Vmath::Vmax(nq, fields[2], 1);
    const int       Maxphieid = Vmath::Imax(nq, fields[2], 1);

    // Coordinates
    const Array<OneD, NekDouble> &x0 = m_x;
    const Array<OneD, NekDouble> &x1 = m_y;

    std::cout << "phim1: Max = " << Maxphim1 << " at x = " << x0[Maxphim1id] << ", y = " << x1[Maxphim1id] << '\n';
    std::cout << "phim2: Max = " << Maxphim2 << " at x = " << x0[Maxphim2id] << ", y = " << x1[Maxphim2id] << '\n';
    std::cout << "phie: Max = " << Maxphie << " at x = " << x0[Maxphieid] << ", y = " << x1[Maxphieid] << '\n';
 
    // TimeMap maxima
    std::cout << "TimeMap: phim1 = " << Vmath::Vmax(nq, TimeMap[0], 1)
                << ", phim2 = " << Vmath::Vmax(nq, TimeMap[1], 1) << ", phie = " << Vmath::Vmax(nq, TimeMap[2], 1) << std::endl;
    
    // variables[4] = "TimeMap_phim";
    Vmath::Vadd(nq, TimeMap[0], 1, TimeMap[1], 1, totfield, 1);
    m_fields[0]->FwdTransLocalElmt(totfield, fieldcoeffs[5]);

    m_fields[0]->FwdTransLocalElmt(TimeMap[0], fieldcoeffs[6]);
    m_fields[0]->FwdTransLocalElmt(TimeMap[1], fieldcoeffs[7]);
    m_fields[0]->FwdTransLocalElmt(TimeMap[2], fieldcoeffs[8]);

    WriteFld(outname1, m_fields[0], fieldcoeffs, variables);
}

void MMFNeuralEP::PlotNeuralEPvar3CSD(
    const Array<OneD, const Array<OneD, NekDouble>> &fields,
    const Array<OneD, const Array<OneD, NekDouble>> &TimeMap,
    const int nstep)
{
    const int nvar    = 8;
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
    variables[0] = "field";
    variables[1] = "phi_m";
    variables[2] = "rho";
    variables[3] = "phi_e";
    variables[4] = "CSDatnode";
    variables[5] = "TimeMap_phim";
    variables[6] = "TimeMap_rho";
    variables[7] = "TimeMap_phie";

    Array<OneD, NekDouble> phim(nq);
    Vmath::Vmul(nq, m_intrazone, 1, fields[0], 1, phim, 1);

    Array<OneD, NekDouble> rho(nq);
    Vmath::Vcopy(nq, fields[1], 1, rho, 1);

    Array<OneD, NekDouble> phie(nq);
    Vmath::Vmul(nq, m_outerzone, 1, fields[2], 1, phie, 1);

    m_fields[0]->FwdTransLocalElmt(phim, fieldcoeffs[1]);
    m_fields[0]->FwdTransLocalElmt(rho, fieldcoeffs[2]);
    m_fields[0]->FwdTransLocalElmt(phie, fieldcoeffs[3]);

    Array<OneD, NekDouble> totfield(nq);
    Vmath::Smul(nq, 100.0 , phie, 1, phie, 1);
    Vmath::Vadd(nq, phim, 1, phie, 1, totfield, 1);
    m_fields[0]->FwdTransLocalElmt(totfield, fieldcoeffs[0]);

    Array<OneD, NekDouble> CSDatnode(nq);
    CSDatnode = ComputeMMFDiffusion(m_CSDmovingframes, m_fields[2]->GetPhys());
    Vmath::Vmul(nq, m_nodezone, 1, CSDatnode, 1, CSDatnode, 1);
    Vmath::Neg(nq, CSDatnode, 1);

    m_fields[0]->FwdTransLocalElmt(CSDatnode, fieldcoeffs[4]);

    // Max values and indices
    const NekDouble Maxphim  = Vmath::Vmax(nq, fields[0], 1);
    const int       Maxphimid = Vmath::Imax(nq, fields[0], 1);

    const NekDouble Maxrho  = Vmath::Vmax(nq, fields[1], 1);
    const int       Maxrhoid = Vmath::Imax(nq, fields[1], 1);

    const NekDouble Maxphie  = Vmath::Vmax(nq, fields[2], 1);
    const int       Maxphieid = Vmath::Imax(nq, fields[2], 1);

    // Coordinates
    const Array<OneD, NekDouble> &x0 = m_x;
    const Array<OneD, NekDouble> &x1 = m_y;

    std::cout << "phim: Max = " << Maxphim << " at x = " << x0[Maxphimid] << ", y = " << x1[Maxphimid] << '\n';
    std::cout << "rho: Max = " << Maxrho << " at x = " << x0[Maxrhoid] << ", y = " << x1[Maxrhoid] << '\n';
    std::cout << "phie: Max = " << Maxphie << " at x = " << x0[Maxphieid] << ", y = " << x1[Maxphieid] << '\n';
 
    // TimeMap maxima
    std::cout << "TimeMap: phim = " << Vmath::Vmax(nq, TimeMap[0], 1)
                << ", rho = " << Vmath::Vmax(nq, TimeMap[1], 1) << ", phie = " << Vmath::Vmax(nq, TimeMap[2], 1) << std::endl;
    
    // variables[4] = "TimeMap_phim";
    m_fields[0]->FwdTransLocalElmt(TimeMap[0], fieldcoeffs[5]);
    m_fields[0]->FwdTransLocalElmt(TimeMap[1], fieldcoeffs[6]);
    m_fields[0]->FwdTransLocalElmt(TimeMap[2], fieldcoeffs[7]);

    WriteFld(outname1, m_fields[0], fieldcoeffs, variables);
}


void MMFNeuralEP::PlotNeuralEPvar4(
    const Array<OneD, const Array<OneD, NekDouble>> &fields,
    const Array<OneD, const Array<OneD, NekDouble>> &TimeMap,
    const int nstep)
{
    const int nvar    = 11;
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
    variables[0] = "field";
    variables[1] = "phi_m";
    variables[2] = "phi_e";
    variables[3] = "phi_m1";
    variables[4] = "phi_m2";
    variables[5] = "phi_e1";
    variables[6] = "phi_e2";
    variables[7] = "TimeMap_phim1";
    variables[8] = "TimeMap_phim2";
    variables[9] = "TimeMap_phie1";
    variables[10] = "TimeMap_phie2";

    Array<OneD, NekDouble> phim(nq);
    Vmath::Vadd(nq, fields[0], 1, fields[1], 1, phim, 1);
    Vmath::Vmul(nq, m_intrazone, 1, phim, 1, phim, 1);

    Array<OneD, NekDouble> phie(nq);
    Vmath::Vadd(nq, fields[2], 1, fields[3], 1, phie, 1);
    Vmath::Vmul(nq, m_outerzone, 1, phie, 1, phie, 1);

    m_fields[0]->FwdTransLocalElmt(phim, fieldcoeffs[1]);
    m_fields[0]->FwdTransLocalElmt(phie, fieldcoeffs[2]);

    // Multiply 100 to \phi_e for better visualization
    Array<OneD, NekDouble> totfield(nq);
    Vmath::Smul(nq, 100.0 , phie, 1, phie, 1);
    Vmath::Vadd(nq, phim, 1, phie, 1, totfield, 1);
    m_fields[0]->FwdTransLocalElmt(totfield, fieldcoeffs[0]);

    m_fields[0]->FwdTransLocalElmt(fields[0], fieldcoeffs[3]);
    m_fields[0]->FwdTransLocalElmt(fields[1], fieldcoeffs[4]);

    m_fields[0]->FwdTransLocalElmt(fields[2], fieldcoeffs[5]);
    m_fields[0]->FwdTransLocalElmt(fields[3], fieldcoeffs[6]);

    // Max values and indices
    const NekDouble Maxphim1  = Vmath::Vmax(nq, fields[0], 1);
    const int       Maxphim1id = Vmath::Imax(nq, fields[0], 1);

    const NekDouble Maxphim2  = Vmath::Vmax(nq, fields[1], 1);
    const int       Maxphim2id = Vmath::Imax(nq, fields[1], 1);

    const NekDouble Maxphie1  = Vmath::Vmax(nq, fields[2], 1);
    const int       Maxphie1id = Vmath::Imax(nq, fields[2], 1);

    const NekDouble Maxphie2  = Vmath::Vmax(nq, fields[3], 1);
    const int       Maxphie2id = Vmath::Imax(nq, fields[3], 1);

    // Coordinates
    const Array<OneD, NekDouble> &x0 = m_x;
    const Array<OneD, NekDouble> &x1 = m_y;

    std::cout << "phim1: Max = " << Maxphim1 << " at x = " << x0[Maxphim1id] << ", y = " << x1[Maxphim1id] << '\n';
    std::cout << "phim2: Max = " << Maxphim2 << " at x = " << x0[Maxphim2id] << ", y = " << x1[Maxphim2id] << '\n';
    std::cout << "phie1: Max = " << Maxphie1 << " at x = " << x0[Maxphie1id] << ", y = " << x1[Maxphie1id] << '\n';
    std::cout << "phie2: Max = " << Maxphie2 << " at x = " << x0[Maxphie2id] << ", y = " << x1[Maxphie2id] << '\n';

    // TimeMap maxima
    std::cout << "TimeMap: phim1 = " << Vmath::Vmax(nq, TimeMap[0], 1)
                << ", phim2 = " << Vmath::Vmax(nq, TimeMap[1], 1) 
                << ", phie1 = " << Vmath::Vmax(nq, TimeMap[2], 1) 
                << ", phie2 = " << Vmath::Vmax(nq, TimeMap[3], 1)<< std::endl;
    
    // variables[4] = "TimeMap_phim";
    m_fields[0]->FwdTransLocalElmt(TimeMap[0], fieldcoeffs[7]);
    m_fields[0]->FwdTransLocalElmt(TimeMap[1], fieldcoeffs[8]);
    m_fields[0]->FwdTransLocalElmt(TimeMap[2], fieldcoeffs[9]);
    m_fields[0]->FwdTransLocalElmt(TimeMap[3], fieldcoeffs[10]);

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

    const Array<OneD, NekDouble> &x0 = m_x;
    const Array<OneD, NekDouble> &x1 = m_y;

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

    const Array<OneD, NekDouble> &x0 = m_x;
    const Array<OneD, NekDouble> &x1 = m_y;

    Array<OneD, NekDouble> phim(nq);
    Array<OneD, NekDouble> phie(nq);

    Vmath::Vmul(nq, &m_intrazone[0], 1, &field[0][0], 1, &phim[0], 1);
    Vmath::Vmul(nq, &m_outerzone[0], 1, &field[1][0], 1, &phie[0], 1);

    NekDouble Maxphim = Vmath::Vmax(nq, phim, 1);
    int Maxphimid = Vmath::Imax(nq, phim, 1);

    NekDouble Maxphie = Vmath::Vmax(nq, phie, 1);
    int Maxphieid = Vmath::Imax(nq, phie, 1);

   std::cout << "phim: Max = " << Maxphim << " at x = " << x0[Maxphimid] << ", y = " << x1[Maxphimid] << std::endl;
   std::cout << "phie: Max = " << Maxphie << " at x = " << x0[Maxphieid] << ", y = " << x1[Maxphieid] << std::endl;
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

// Implicit solve for NeuralEP 2D solver
void MMFNeuralEP::DoImplicitSolveNeuralEP2Dmono(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time,
    const NekDouble lambda)
{
    (void) time;

    const int nq   = m_fields[0]->GetNpoints();

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
    (void) time;

    const int nq   = m_fields[0]->GetNpoints();

    // Set Helmholtz coefficients
    StdRegions::ConstFactorMap factors;
    factors[StdRegions::eFactorTau] = m_Helmtau;
    factors[StdRegions::eFactorLambda] = m_Cn * m_Rf / lambda;

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

// Implicit solve for NeuralEP 2D solver
void MMFNeuralEP::DoImplicitSolveNeuralEP2DbiMulti(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time,
    const NekDouble lambda)
{
    (void) time;

    const int nq   = m_fields[0]->GetNpoints();
    const int numfiber = m_numfiber;

    // Set Helmholtz coefficients
    StdRegions::ConstFactorMap factors;
    factors[StdRegions::eFactorTau] = m_Helmtau;
    factors[StdRegions::eFactorLambda] = m_Cn * m_Rf / lambda;

    for (int n = 0; n < numfiber; ++n)
    {
        if (outarray[n].size() != nq)
        {
            outarray[n] = Array<OneD, NekDouble>(nq,0.0);
        }
    }

    // Multiply 1.0/timestep
    const NekDouble scale = -factors[StdRegions::eFactorLambda];
    for (int n = 0; n < numfiber; ++n)
    {
        Vmath::Smul(nq, scale, inarray[n], 1, m_fields[n]->UpdatePhys(), 1);
        m_fields[n]->HelmSolve(m_fields[n]->GetPhys(), m_fields[n]->UpdateCoeffs(),
                            factors, m_varcoefffiber[n]);
        m_fields[n]->BwdTrans(m_fields[n]->GetCoeffs(), outarray[n]);
        m_fields[n]->SetPhysState(true);
    }
}

// Implicit solve for NeuralEP 2D solver
void MMFNeuralEP::DoImplicitSolveNeuralEP2DbiMultiv2(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time,
    const NekDouble lambda)
{
    (void) time;

    const int nq   = m_fields[0]->GetNpoints();
    const int numfiber = m_numfiber;

    // Set Helmholtz coefficients
    StdRegions::ConstFactorMap factors;
    factors[StdRegions::eFactorTau] = m_Helmtau;
    factors[StdRegions::eFactorLambda] = m_Cn * m_Rf / lambda;

    for (int n = 0; n < numfiber; ++n)
    {
        if (outarray[n].size() != nq)
        {
            outarray[n] = Array<OneD, NekDouble>(nq,0.0);
        }
    }

    // Multiply 1.0/timestep
    const NekDouble scale = -factors[StdRegions::eFactorLambda];
    for (int n = 0; n < numfiber; ++n)
    {
        Vmath::Smul(nq, scale, inarray[n], 1, m_fields[n]->UpdatePhys(), 1);
        m_fields[n]->HelmSolve(m_fields[n]->GetPhys(), m_fields[n]->UpdateCoeffs(),
                            factors, m_varcoefffiber[n]);
        m_fields[n]->BwdTrans(m_fields[n]->GetCoeffs(), outarray[n]);
        m_fields[n]->SetPhysState(true);
    }
}

// Implicit solve for NeuralEP 2D solver
void MMFNeuralEP::DoImplicitSolveNeuralEP2DbiCSD(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time,
    const NekDouble lambda)
{
    (void) time;

    const int nq   = m_fields[0]->GetNpoints();

    // Set Helmholtz coefficients
    StdRegions::ConstFactorMap factors;
    factors[StdRegions::eFactorTau] = m_Helmtau;
    factors[StdRegions::eFactorLambda] = m_Cn * m_Rf / lambda;

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

    // Set CSD Helmholtz coefficients
    StdRegions::ConstFactorMap CSDfactors;
    CSDfactors[StdRegions::eFactorTau] = m_Helmtau;
    // CSDfactors[StdRegions::eFactorLambda] = m_CSDDiff / lambda;
    CSDfactors[StdRegions::eFactorLambda] = m_Cn * m_Rf / lambda;

    // Multiply 1.0/timestep
    const NekDouble CSDscale = -CSDfactors[StdRegions::eFactorLambda];
    Vmath::Smul(nq, CSDscale, inarray[1], 1, m_fields[1]->UpdatePhys(), 1);

    m_fields[1]->HelmSolve(m_fields[1]->GetPhys(), m_fields[1]->UpdateCoeffs(),
                        CSDfactors, m_CSDvarcoeff);
    m_fields[1]->BwdTrans(m_fields[1]->GetCoeffs(), outarray[1]);
    m_fields[1]->SetPhysState(true);
}

// We Return Y[i] = rhs [i] without no Helomsolver
void MMFNeuralEP::DoNullSolve(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time,
    const NekDouble lambda)
{
    (void) time;
    (void) lambda;

    const int nvariables = inarray.size();
    const int nq         = m_fields[0]->GetNpoints();

    for (int i = 0; i < nvariables; ++i)
    {
        Vmath::Vcopy(nq, &inarray[i][0], 1, &outarray[i][0], 1);
    }
}

void MMFNeuralEP::GetFluxVector(
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    const Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &qfield,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &viscousTensor)
{
    (void) inarray;

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
        m_stimulus[n]->Update(m_zoneindexfiber[n], outarray[0], time);
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
    const int phievar = nvar - 1;
    const NekDouble factor = m_Cn * m_Rf;
    const NekDouble Temp = m_Temperature;

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
    m_neuron->TimeIntegrate(m_zoneindex, inarray[0], outarray[0], time, Temp);

    // 2. Apply Stimulus
    for (std::size_t n = 0; n < m_stimulus.size(); ++n)
    {
      //  m_stimulus[n]->Update(m_zoneindexfiber[n], outarray[0], time);
       m_stimulus[n]->Update(m_zoneindexfiber[n], outarray[0], time);
    }

    // 3. Compute phi_e to satisfy bidomain coupling
    // \nabla \cdot ( (\signa_e + \sigma_i) \nabla \phi_e) = - \nabla \cdot
    // (\sigma_i \nabla \phi_m)
    m_fields[phievar]->UpdatePhys() = ComputeFieldPhie(phievar, inarray[0]);

    // 4. Compute \nabla \cdot (\sigma_i \nabla \phi_e) and add to membrane current
    Array<OneD, NekDouble> phiecurrent(nq);
    phiecurrent = ComputeMMFDiffusion(m_movingframes, m_fields[phievar]->GetPhys());

    // Current caused by extracellular potential affects the total current at the nodes and myelin.
    for (int i = 0; i < nq; ++i)
    {
        if (m_intrazone[i] > 0.0)
        {
            outarray[0][i] += phiecurrent[i] / factor;
        }
    }

    if (m_explicitDiffusion)
    {
        static thread_local Array<OneD, NekDouble> Laplacian;
        if (Laplacian.size() != nq)
        {
            Laplacian = Array<OneD, NekDouble>(nq);
        }

        WeakDGMMFDiffusion(0, inarray[0], Laplacian, time);

        for (int i = 0; i < nq; ++i)
        {
            outarray[0][i] += Laplacian[i] / factor;
        }
    }
}


// var = 0: phim_1 = phim in fiber 1
// var = 1: phim_2 = phim in fiber 2
// var = 2: phie
void MMFNeuralEP::DoOdeRhsNeuralEP2DbiMulti(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time)
{
    const int nvar = m_fields.size();
    const int nq   = m_fields[0]->GetNpoints();
    const int phievar = nvar - 1;
    const int numfiber = m_numfiber;

    const NekDouble factor = m_Cn * m_Rf;
    const NekDouble Temp = m_Temperature;

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
    m_neuron->TimeIntegrateMulti(numfiber, m_zoneindexfiber, inarray, outarray, time, Temp);

    Array<OneD, NekDouble> tmp(nq);
    Array<OneD, NekDouble> phie(nq, 0.0);
    for (int n = 0; n < numfiber; ++n)
    {
        // 2. Apply Stimulus
        m_stimulus[n]->Update(m_zoneindexfiber[n], outarray[n], time);

        // 3. Compute phi_e to satisfy bidomain coupling
        tmp = ComputeFieldPhiefiber(phievar, n, inarray[n]);

        Vmath::Vadd(nq, tmp, 1, phie, 1, phie, 1);
    }

    m_fields[phievar]->UpdatePhys() = phie;

    // 4. Compute \nabla \cdot (\sigma_i \nabla \phi_e) and add to membrane current
    Array<OneD, NekDouble> phiecurrent(nq);
    for (int n=0; n < numfiber; ++n)
    {
        phiecurrent = ComputeMMFDiffusion(m_movingframesfiber[n], m_fields[phievar]->GetPhys());

        // Current caused by extracellular potential affects the total current at the nodes and myelin.
        for (int i = 0; i < nq; ++i)
        {
            if (m_intrazonefiber[n][i] > 0.0)
            {
                outarray[n][i] += phiecurrent[i] / factor;
            }
        }
    }

    if (m_explicitDiffusion)
    {
        static thread_local Array<OneD, NekDouble> Laplacian;
        if (Laplacian.size() != nq)
        {
            Laplacian = Array<OneD, NekDouble>(nq);
        }

        WeakDGMMFDiffusion(0, inarray[0], Laplacian, time);

        for (int i = 0; i < nq; ++i)
        {
            outarray[0][i] += Laplacian[i] / factor;
        }
    }
}

// var = 0: phim_1 = phim in fiber 1
// var = 1: phim_2 = phim in fiber 2
// var = 2: phie
void MMFNeuralEP::DoOdeRhsNeuralEP2DbiMultiv2(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time)
{
    const int nvar = m_fields.size();
    const int nq   = m_fields[0]->GetNpoints();
    const int phievar = nvar - 1;
    const int numfiber = m_numfiber;

    const NekDouble factor = m_Cn * m_Rf;
    const NekDouble Temp = m_Temperature;

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
    m_neuron->TimeIntegrateMulti(numfiber, m_zoneindexfiber, inarray, outarray, time, Temp);

    // 2. Apply Stimulus
    for (int n = 0; n < numfiber; ++n)
    {
        m_stimulus[n]->Update(m_zoneindexfiber[n], outarray[n], time);
    }

    // Add phim for all fibers to produce the total phim
    Array<OneD, NekDouble> phim(nq);
    Vmath::Vadd(nq, inarray[0], 1, inarray[1], 1, phim, 1);
    m_fields[phievar]->UpdatePhys() = ComputeFieldPhiefiberv2(phievar, phim);

    // 4. Compute \nabla \cdot (\sigma_i \nabla \phi_e) and add to membrane current
    Array<OneD, NekDouble> phiecurrent(nq);
    for (int n=0; n < numfiber; ++n)
    {
        phiecurrent = ComputeMMFDiffusion(m_movingframesfiber[n], m_fields[phievar]->GetPhys());

        // Current caused by extracellular potential affects the total current at the nodes and myelin.
        for (int i = 0; i < nq; ++i)
        {
            if (m_intrazonefiber[n][i] > 0.0)
            {
                outarray[n][i] += phiecurrent[i] / factor;
            }
        }
    }

    if (m_explicitDiffusion)
    {
        static thread_local Array<OneD, NekDouble> Laplacian;
        if (Laplacian.size() != nq)
        {
            Laplacian = Array<OneD, NekDouble>(nq);
        }

        WeakDGMMFDiffusion(0, inarray[0], Laplacian, time);

        for (int i = 0; i < nq; ++i)
        {
            outarray[0][i] += Laplacian[i] / factor;
        }
    }
}



void MMFNeuralEP::DoOdeRhsNeuralEP2DbiCSD(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time)
{
    const int nvar = m_fields.size();
    const int phievar = nvar - 1;
    const int nq   = m_fields[0]->GetNpoints();
    const NekDouble factor = m_Cn * m_Rf;
    const NekDouble Temp = m_Temperature;

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
    m_neuron->TimeIntegrate(m_zoneindex, inarray[0], outarray[0], time, Temp);

    // 2. Apply Stimulus
    for (std::size_t n = 0; n < m_stimulus.size(); ++n)
    {
       m_stimulus[n]->Update(m_zoneindexfiber[n], outarray[0], time);
    }

    // 3. Compute phi_e to satisfy bidomain coupling
    // \nabla \cdot ( (\signa_e + \sigma_i) \nabla \phi_e) = - \nabla \cdot
    // (\sigma_i \nabla \phi_m)
    m_fields[phievar]->UpdatePhys() = ComputeFieldPhie(phievar, inarray[0]);

    // 4. Compute \nabla \cdot (\sigma_i \nabla \phi_e) and add to membrane current
    Array<OneD, NekDouble> tmp(nq);
    tmp = ComputeMMFDiffusion(m_movingframes, m_fields[phievar]->GetPhys());

    // Current caused by extracellular potential affects the total current at the nodes and myelin.
    for (int i = 0; i < nq; ++i)
    {
        if (m_intrazone[i] > 0.0)
        {
            outarray[0][i] += tmp[i] / factor;
        }
    }

    // Compute the charge density at the nodes.
    tmp = ComputeMMFDiffusion(m_CSDmovingframes, m_fields[phievar]->GetPhys());
    Vmath::Vmul(nq, m_nodezone, 1, tmp, 1, outarray[1], 1);
    Vmath::Neg(nq, outarray[1], 1);

    if (m_explicitDiffusion)
    {
        static thread_local Array<OneD, NekDouble> Laplacian;
        if (Laplacian.size() != nq)
        {
            Laplacian = Array<OneD, NekDouble>(nq);
        }

        WeakDGMMFDiffusion(0, inarray[0], Laplacian, time);

        for (int i = 0; i < nq; ++i)
        {
            outarray[0][i] += Laplacian[i] / factor;
        }
    }
}


// output: phi_e (m_fields[1]->UpdatePhys()) and outarray (1/C_n/r) * \nabla^2 \phi_e
// Compute phi_e from the given distribution of phi_m
// \nabla \cdot ( (1 + \rho) \mathbf{e}_1 + \mathbf{e}_2 ) ( \nabla \phi_e ))
//                         = - \nabla \cdot \mathbf{e}_1 \nabla \phi_m
Array<OneD, NekDouble> MMFNeuralEP::ComputeFieldPhie(
                                    const int phievar,
                                    const Array<OneD, const NekDouble> &phim)
{
    const int nq = m_fields[0]->GetNpoints();

    Array<OneD, NekDouble> outarray(nq);

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
    for (int i = 0; i < nq; ++i)
    {
        m_fields[phievar]->UpdatePhys()[i] = phimcurrent[i] - avg;
    }

    m_fields[phievar]->HelmSolve(m_fields[phievar]->GetPhys(), m_fields[phievar]->UpdateCoeffs(), 
                                 phiefactors, m_phievarcoeff);

    m_fields[phievar]->BwdTrans(m_fields[phievar]->GetCoeffs(), outarray);

    // Make it as a value with AvgInt is zero.
    const NekDouble mean = AvgInt(outarray);
    Vmath::Sadd(nq, -mean, outarray, 1, outarray, 1);

    return outarray;
}

Array<OneD, NekDouble> MMFNeuralEP::ComputeFieldPhiefiber(
    const int phievar, const int nfiber,
    const Array<OneD, const NekDouble> &phim)
{
    const int nq = m_fields[0]->GetNpoints();

    Array<OneD, NekDouble> outarray(nq);

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

    phimcurrent = ComputeMMFDiffusion(m_movingframesfiber[nfiber], phim);
    Vmath::Neg(nq, phimcurrent, 1);

    switch (m_ExtCurrentType)
    {
        case eEphaptic:
        {
            Vmath::Vmul(nq, m_nodezonefiber[nfiber], 1, phimcurrent, 1,
                        phimcurrent, 1);
        }

        case eNoEphaptic:
        {
            Vmath::Zero(nq, phimcurrent, 1);
            break;
        }

        default:
            break;
    }

    // Compute  \nabla \cdot ( (1 + \rho) \mathbf{e}_1 + \mathbf{e}_2 ) ( \nabla
    // \phi_e ))
    //                         = - \nabla \cdot \mathbf{e}_1 \nabla \phi_m
    const NekDouble avg = AvgInt(phimcurrent);
    for (int i = 0; i < nq; ++i)
    {
        m_fields[phievar]->UpdatePhys()[i] = phimcurrent[i] - avg;
    }

    m_fields[phievar]->HelmSolve(m_fields[phievar]->GetPhys(),
                                 m_fields[phievar]->UpdateCoeffs(), phiefactors,
                                 m_phievarcoefffiber[nfiber]);

    m_fields[phievar]->BwdTrans(m_fields[phievar]->GetCoeffs(), outarray);

    // Make it as a value with AvgInt is zero.
    const NekDouble mean = AvgInt(outarray);
    Vmath::Sadd(nq, -mean, outarray, 1, outarray, 1);

    return outarray;
}

Array<OneD, NekDouble> MMFNeuralEP::ComputeFieldPhiefiberv2(
    const int phievar,
    const Array<OneD, const NekDouble> &phim)
{
    const int nq = m_fields[0]->GetNpoints();
    const int numfiber = m_numfiber;

    Array<OneD, NekDouble> outarray(nq);

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

    Array<OneD, NekDouble> tmp(nq);
    for (int n = 0; n < numfiber; ++n)
    {
        tmp = ComputeMMFDiffusion(m_movingframesfiber[n], phim);
        Vmath::Neg(nq, tmp, 1);

        switch (m_ExtCurrentType)
        {
            case eEphaptic:
            {
                Vmath::Vmul(nq, m_nodezonefiber[n], 1, tmp, 1,
                    tmp, 1);
            }

            case eNoEphaptic:
            {
                Vmath::Zero(nq, tmp, 1);
                break;
            }

            default:
                break;
        }

        Vmath::Vadd(nq, tmp, 1, phimcurrent, 1, phimcurrent, 1);
    }

    // Compute  \nabla \cdot ( (1 + \rho) \mathbf{e}_1 + \mathbf{e}_2 ) ( \nabla
    // \phi_e ))
    //                         = - \nabla \cdot \mathbf{e}_1 \nabla \phi_m
    const NekDouble avg = AvgInt(phimcurrent);
    for (int i = 0; i < nq; ++i)
    {
        m_fields[phievar]->UpdatePhys()[i] = phimcurrent[i] - avg;
    }

    m_fields[phievar]->HelmSolve(m_fields[phievar]->GetPhys(),
                                 m_fields[phievar]->UpdateCoeffs(), phiefactors,
                                 m_phievarcoeff);

    m_fields[phievar]->BwdTrans(m_fields[phievar]->GetCoeffs(), outarray);

    // Make it as a value with AvgInt is zero.
    const NekDouble mean = AvgInt(outarray);
    Vmath::Sadd(nq, -mean, outarray, 1, outarray, 1);

    return outarray;
}


// output: phi_e (m_fields[1]->UpdatePhys()) and outarray (1/C_n/r) * \nabla^2 \phi_e
// Compute phi_e from the given distribution of phi_m
// \nabla \cdot ( (1 + \rho) \mathbf{e}_1 + \mathbf{e}_2 ) ( \nabla \phi_e ))
//                         = - \nabla \cdot \mathbf{e}_1 \nabla \phi_m
// void MMFNeuralEP::ComputePhie(
//                     const int phievar,
//                     const Array<OneD, const NekDouble> &phim)
// {
//     const int nq = m_fields[0]->GetNpoints();

//     Array<OneD, NekDouble> outarray(nq);

//     // Solve the Poisson equation: \nabla (\sigma_e + \sigma_i ) phi_e = \nabla
//     // \sigma_i \nabla phi_m
//     StdRegions::ConstFactorMap phiefactors;
//     phiefactors[StdRegions::eFactorTau]    = m_Helmtau;
//     phiefactors[StdRegions::eFactorLambda] = 0.0;

//     // // Compute \nabla \sigma_i \nabla phi_m and use it as point sources for
//     // phi_e. This is equivalently achieved by removing all the point sources in
//     // myelinnated fiber region.    
//     // Allocate only once
//     Array<OneD, NekDouble> phimcurrent(nq);
//     if (phimcurrent.size() != nq)
//     {
//         phimcurrent = Array<OneD, NekDouble>(nq, 0.0);
//     }
//     else
//     {
//         Vmath::Zero(nq, phimcurrent, 1);
//     }

//     phimcurrent = ComputeMMFDiffusion(m_movingframes, phim);
//     Vmath::Neg(nq, phimcurrent, 1);

//     switch(m_ExtCurrentType)
//     {
//         case eEphaptic:
//         {
//             Vmath::Vmul(nq, m_nodezone, 1, phimcurrent, 1, phimcurrent, 1);
//             break;
//         }

//         case eNoEphaptic:
//         {
//             Vmath::Zero(nq, phimcurrent, 1);
//             break;
//         }

//         default:
//         break;
//     }

//     // Compute  \nabla \cdot ( (1 + \rho) \mathbf{e}_1 + \mathbf{e}_2 ) ( \nabla \phi_e ))
//     //                         = - \nabla \cdot \mathbf{e}_1 \nabla \phi_m
//     const NekDouble avg = AvgInt(phimcurrent);
//     for (int i = 0; i < nq; ++i)
//     {
//         m_fields[phievar]->UpdatePhys()[i] = phimcurrent[i] - avg;
//     }

//     m_fields[phievar]->HelmSolve(m_fields[phievar]->GetPhys(), m_fields[phievar]->UpdateCoeffs(), phiefactors, m_phievarcoeff);
//     m_fields[phievar]->BwdTrans(m_fields[phievar]->GetCoeffs(), m_fields[phievar]->UpdatePhys());

//     // Make it as a value with AvgInt is zero.
//     const NekDouble mean = AvgInt(m_fields[phievar]->GetPhys());
//     Vmath::Sadd(nq, -mean, m_fields[phievar]->GetPhys(), 1, m_fields[phievar]->UpdatePhys(), 1);

//     m_fields[phievar]->SetPhysState(true);
// }


void MMFNeuralEP::v_SetInitialConditions(NekDouble initialtime,
                                         bool dumpInitialConditions,
                                         const int domain)
{
    (void) domain;
    (void) dumpInitialConditions;
    
    const int nq = GetTotPoints();

    switch (m_NeuralEPType)
    {
        case eNeuralEP2Dmono:
        case eNeuralEP2Dbi:
        case eNeuralEP2DbiCSD:
        {
            m_neuron->Initialise();

            // Read initial condition from xml file
            EquationSystem::v_SetInitialConditions(initialtime, false);

            Array<OneD, Array<OneD, NekDouble>> tmp(1);
            tmp[0] = Array<OneD, NekDouble>(nq);

            Array<OneD, NekDouble> initialcondition(nq, 0.0);
            m_fields[0]->SetPhys(initialcondition);
            break;
        }

        case eNeuralEP2DbiMulti:
        case eNeuralEP2DbiMultiv2:
        {
            int numfiber = m_numfiber;

            m_neuron->InitialiseMulti(numfiber);

            // Read initial condition from xml file
            EquationSystem::v_SetInitialConditions(initialtime, false);

            Array<OneD, Array<OneD, NekDouble>> tmp(1);
            tmp[0] = Array<OneD, NekDouble>(nq);

            Array<OneD, NekDouble> initialcondition(nq, 0.0);
            for (int n=0; n < numfiber; ++n)
            {
               m_fields[n]->SetPhys(initialcondition);
            }

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

    SolverUtils::AddSummaryItem(s, "Temperature", m_Temperature);
    SolverUtils::AddSummaryItem(s, "Helmtau", m_Helmtau);
    SolverUtils::AddSummaryItem(s, "nq", nq);

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
    std::string vDriverModule;
    DriverSharedPtr drv;

    try
    {
        // Create session reader.
        session = LibUtilities::SessionReader::CreateInstance(argc, argv);

        // Create MeshGraph
        graph = SpatialDomains::MeshGraphIO::Read(session);

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

// int main(int argc, char *argv[])
// {
//     LibUtilities::SessionReaderSharedPtr session;
//     SpatialDomains::MeshGraphSharedPtr graph;

//     LibUtilities::SessionReaderSharedPtr session1D;
//     SpatialDomains::MeshGraphSharedPtr graph1D;

//     std::string vDriverModule;
//     DriverSharedPtr drv;

//     try
//     {
//         // Create session reader.
//         session = LibUtilities::SessionReader::CreateInstance(argc, argv);

//         // Create MeshGraph
//         graph1D = SpatialDomains::MeshGraphIO::Read(session);

//         // Create driver
//         session->LoadSolverInfo("Driver", vDriverModule, "Standard");
//         drv = GetDriverFactory().CreateInstance(vDriverModule, session, graph);

//         // Execute driver
//         drv->Execute();

//         // Finalise session
//         session->Finalise();
//     }

//     catch (const std::runtime_error &e)
//     {
//         return 1;
//     }
//     catch (const std::string &eStr)
//     {
//         std::cout << "Error: " << eStr << std::endl;
//     }

//     return 0;
// }
