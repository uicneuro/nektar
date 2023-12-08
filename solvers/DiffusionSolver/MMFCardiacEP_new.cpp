///////////////////////////////////////////////////////////////////////////////
//
// File MMFCardiacEP.cpp
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
// Description: MMFCardiacEP.
//
///////////////////////////////////////////////////////////////////////////////
#include <iomanip>
#include <iostream>
#include <stdio.h>
#include <string.h>

#include <boost/algorithm/string.hpp>
#include <boost/math/special_functions/spherical_harmonic.hpp>

#include <DiffusionSolver/EquationSystems/MMFCardiacEP.h>

#include <CardiacEPSolver/Filters/FilterCellHistoryPoints.h>
#include <CardiacEPSolver/Filters/FilterCheckpointCellModel.h>

#include <SolverUtils/Driver.h>
#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>

#include <LibUtilities/BasicUtils/Timer.h>
#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>
#include <LibUtilities/BasicUtils/SessionReader.h>

#ifdef NEKTAR_USE_MPI
#include <LibUtilities/Communication/CommMpi.h>
#endif

using namespace std;
using namespace Nektar::SolverUtils;
using namespace Nektar;

namespace Nektar
{
string MMFCardiacEP::className =
    SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
        "MMFCardiacEP", MMFCardiacEP::create, "MMFCardiacEP equation.");

MMFCardiacEP::MMFCardiacEP(const LibUtilities::SessionReaderSharedPtr &pSession,
                           const SpatialDomains::MeshGraphSharedPtr &pGraph)
    : UnsteadySystem(pSession, pGraph), MMFSystem(pSession, pGraph)
{
}

void MMFCardiacEP::v_InitObject(bool DeclareFields)
{
    UnsteadySystem::v_InitObject(DeclareFields);

    int nq   = GetTotPoints();

    // Derive AnisotropyStrength.
    m_AniStrength = Array<OneD, Array<OneD, NekDouble>> (m_expdim);
    for (int j = 0; j < m_expdim; ++j)
    {
        m_AniStrength[j] = Array<OneD, NekDouble>(nq, 1.0);
    }

    m_TimeMap = Array<OneD, Array<OneD, NekDouble>>(1);
    for (int i = 0; i < 1; ++i)
    {
        m_TimeMap[i] = Array<OneD, NekDouble>(nq,0.0);
    }

    // Conductance parameters
    m_session->LoadParameter("Chi", m_chi, 28.0);
    m_session->LoadParameter("Cm", m_capMembrane, 0.125);

    // Resting potential
    m_session->LoadParameter("urest", m_urest, 0.0);

    // Helmsolver parameter
    m_session->LoadParameter("Helmtau", m_Helmtau, 1.0);

    m_session->LoadParameter("AnisotropyStrength", m_AnisotropyStrength, 4.0);

    m_session->LoadParameter("AniRegionStart", m_AniRegionStart, 0);
    m_session->LoadParameter("AniRegionEnd", m_AniRegionEnd, 1000000);

    m_session->LoadParameter("Gaussiantau", m_Gaussiantau, 10);

    //Aliev-Panfilov Parameter
    m_session->LoadParameter("k", m_k, 0.0);
    m_session->LoadParameter("a", m_a, 0.0);
    m_session->LoadParameter("mu1", m_mu1, 0.0);
    m_session->LoadParameter("mu2", m_mu2, 0.0);
    m_session->LoadParameter("eps", m_eps, 0.0);

    // TimeMap ?
    m_session->LoadParameter("TimeMapScheme", m_TimeMapScheme, 0);
    m_session->LoadParameter("TimeMapStart", m_TimeMapStart, 0.0);
    m_session->LoadParameter("TimeMapEnd", m_TimeMapEnd, 10000.0);

    std::string vCellModel;
    m_session->LoadSolverInfo("CELLMODEL", vCellModel, "FitzHughNagumo");

    ASSERTL0(vCellModel != "", "Cell Model not specified.");

    m_cell = GetCellModelFactory().CreateInstance(vCellModel, m_session,
                                                  m_fields[0]);

    // Stimulus
    m_stimulus = Stimulus::LoadStimuli(m_session, m_fields[0]);

    m_session->LoadParameter("Diffbeta", m_Diffbeta, 0.5);
    m_session->LoadParameter("Diffeta", m_Diffeta, 100.0);
    m_session->LoadParameter("Diffhe", m_Diffhe, 0.5);

    m_session->LoadParameter("PVcond", m_PVcond, 1.0);

    m_session->LoadParameter("ScarSize", m_ScarSize, 0.0);
    m_session->LoadParameter("ScarStr", m_ScarStr, 1.0);
    m_session->LoadParameter("ScarPis", m_ScarPis, 0.2);
    m_session->LoadParameter("ScarLocx", m_ScarLocx, 0.0);
    m_session->LoadParameter("ScarLocy", m_ScarLocy, 0.0);
    m_session->LoadParameter("ScarLocz", m_ScarLocz, 0.0);

    m_session->LoadParameter("RelDivSize", m_RelDivSize, 0.0);
    m_session->LoadParameter("RelDivStr", m_RelDivStr, 1.0);
    m_session->LoadParameter("RelDivPis", m_RelDivPis, 1.0);
    m_session->LoadParameter("RelDivLocx", m_RelDivLocx, 10.0);

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

    // TimeMap and its tangent vector loading
    if(m_SolverSchemeType==eTimeMapMarch)
    {
        std::cout << "======= Start loading TimeMap: eTimeMapMarch  =======" << std::endl;
        m_session->LoadParameter("TimeMapIapp", m_TimeMapIapp, 0.2);
        m_session->LoadParameter("TimeMapnstep", m_TimeMapnstep, 10000);

        // Import TimeMap
        m_session->LoadSolverInfo("TMsessionName", m_TMsessionName, "m_sessionName");
        std::string loadname = m_TMsessionName + "_TimeMap_" +
                            boost::lexical_cast<std::string>(m_TimeMapnstep) + ".chk";

        Array<OneD, NekDouble> m_ValidTM(nq, 1.0);
        LoadTimeMap(loadname, m_ValidTM, m_TimeMap, m_AniStrength, m_TMvelocitymag, m_TMvelocity);
    }

    // Ratio between along the fiber and orthogonal to the fiber
    switch (m_MediumType)
    {
        case eAnisotropy:
        case eHeterogeneousAnisotropy:
        {
            Array<OneD, NekDouble> CardiacFibre;
            LoadCardiacFiber(m_MediumType, m_AnisotropyStrength, m_AniStrength,
                             CardiacFibre);
            MMFSystem::MMFInitObject(m_AniStrength, CardiacFibre);
        }
        break;

        case eHeterogeneousIsotropy:
        case eRegionalHeterogeneous:
        case eRegionalIsotropy:
        case eGaussianHeterogeneous:
        {
            LoadCardiacFiber(m_MediumType, m_AnisotropyStrength, m_AniStrength);
            MMFSystem::MMFInitObject(m_AniStrength);
        }
        break;

        case eIsotropy:
        default:
        {
            Array<OneD, NekDouble> Anitmp(nq, 1.0);

            /*
            // Scar Tissue case
            if (m_ScarSize > 0.0001)
            {
                Anitmp = SmoothCircle(m_ScarStr, m_ScarPis, m_ScarLocx,
                                      m_ScarLocy, m_ScarLocz, m_ScarSize);
            }

            // Relative Divergence Case
            // if(m_RelDivSize>0.001)
            // {
            //     Anitmp = SmoothLine(m_RelDivStr, m_RelDivPis, m_RelDivLocx,
            //     m_RelDivSize);
            // }

            // PV conductance
            // NekDouble diff;
            if ((m_PVcond > 1.0) || (m_PVcond < 1.0))
            {
                Array<OneD, NekDouble> PVAnitmp(nq, 1.0);
                PVAnitmp = SmoothLine(m_PVcond, 1.0, 10.0, 2.0);

                // for (int i=0; i<nq; ++i)
                // {
                //     diff = m_ScarStr - Anitmp[i];
                //     if( fabs(diff) < 0.1 )
                //     {
                //         PVAnitmp[i] = Anitmp[i];
                //     }
                // }

                Vmath::Vcopy(nq, PVAnitmp, 1, Anitmp, 1);
            }

            for (int j = 0; j < m_expdim; ++j)
            {
                Vmath::Vcopy(nq, &Anitmp[0], 1, &m_AniStrength[j][0], 1);
            }
            */

            MMFSystem::MMFInitObject(m_AniStrength);
        }
        break;
    }

    // Setup moving frames for the differentiation of time map
    std::string MMFdirStr = "LOCAL";
    m_session->LoadSolverInfo("TMMMFDir", MMFdirStr, "LOCAL");
    m_TMMMFdir = FindMMFdir(MMFdirStr);
    
    m_TMAniStrength = Array<OneD, Array<OneD, NekDouble>> (1);
    for (int j = 0; j < 1; ++j)
    {
        m_TMAniStrength[j] = Array<OneD, NekDouble>(nq, 1.0);
    }

    std::cout << "================== SetUP Time Map moving frames ==================" << std::endl;
    SetUpMovingFrames(m_TMMMFdir, m_TMAniStrength, m_TMmovingframes);

    if(m_SolverSchemeType==eTimeMapDeform)
    {
        // load original timemap
        std::string TMsessionName_old;
        m_session->LoadSolverInfo("TMsessionName_old", TMsessionName_old, "m_sessionName");

        // consider a new conductivity map
        m_session->LoadSolverInfo("TMsessionName", m_TMsessionName, "m_sessionName");

        ComputeVarCoeff2D(m_movingframes, m_varcoeff);
        Array<OneD, NekDouble> TimeMapDiff = ComputeTimeMapDeform(TMsessionName_old, m_TMsessionName); 

        wait_on_enter();
    }

    if (m_explicitDiffusion)
    {
        m_ode.DefineImplicitSolve(&MMFCardiacEP::DoNullSolve, this);
        m_ode.DefineProjection(&MMFCardiacEP::DoOdeProjection, this);
    }

    else
    {
        // Create varcoeff for Helmsolver
        ComputeVarCoeff2D(m_movingframes, m_varcoeff);
        m_ode.DefineImplicitSolve(&MMFCardiacEP::DoImplicitSolveCardiacEP, this);
    }

    if(m_SolverSchemeType==eTimeMapMarch)
    {
        m_ode.DefineOdeRhs(&MMFCardiacEP::DoOdeRhsCardiacEPTimeMap, this);
    }
    
    else
    {
        m_ode.DefineOdeRhs(&MMFCardiacEP::DoOdeRhsCardiacEP, this);
    }

    // Wait for initialization is done
    // LibUtilities::CommMpi::v_Block();
}

/**
 *
 */
MMFCardiacEP::~MMFCardiacEP()
{
}

Array<OneD, NekDouble> MMFCardiacEP::ComputeTimeMapDeform(const std::string &sessionold, const std::string &sessionnew) 
    {
        int nq      = GetNpoints();

        Array<OneD, Array<OneD, NekDouble>> Velocity_old(m_spacedim);
        Array<OneD, NekDouble> Velocitymag_old(nq);
        Array<OneD, Array<OneD, NekDouble>> AniStrength_old(m_expdim);
        Array<OneD, Array<OneD, NekDouble>> TimeMap_old(1);

        Array<OneD, Array<OneD, NekDouble>> Velocity_new(m_spacedim);
        Array<OneD, NekDouble> Velocitymag_new(nq);
        Array<OneD, Array<OneD, NekDouble>> AniStrength_new(m_expdim);
        Array<OneD, Array<OneD, NekDouble>> TimeMap_new(1);

        Array<OneD, NekDouble> VdiffMag(nq);
        Array<OneD, NekDouble> VdiffDivergence(nq);
        
        /////
        Array<OneD, NekDouble> x0(nq);
        Array<OneD, NekDouble> x1(nq);
        Array<OneD, NekDouble> x2(nq);

        m_fields[0]->GetCoords(x0, x1, x2);

        Array<OneD, NekDouble> x0avg(nq);
        Array<OneD, NekDouble> x1avg(nq);
        Array<OneD, NekDouble> x2avg(nq);

        m_fields[0]->ComputeCellAvg(x0, x1, x2, x0avg, x1avg, x2avg);

        NekDouble x0limitm = 15.0;
        NekDouble x0limitp = 5.0;
        NekDouble x0max = Vmath::Vmax(nq, x0, 1);
        NekDouble x0min = Vmath::Vmin(nq, x0, 1);

        Array<OneD, NekDouble> m_ValidTM(nq, 1.0);
        for (int i=0; i<nq; ++i)
        {
            if( (x0avg[i]<x0min+x0limitm) || (x0avg[i]>x0max-x0limitp) )
            {
                m_ValidTM[i] = 0.0;
            }
        }

        // load old timemap                   
        m_session->LoadParameter("TimeMapnstep", m_TimeMapnstep, 10000);
        
        std::string loadname_old = sessionold + "_TimeMap_" +
                        boost::lexical_cast<std::string>(m_TimeMapnstep) + ".chk";

        LoadTimeMap(loadname_old, m_ValidTM, TimeMap_old, AniStrength_old, Velocitymag_old, Velocity_old);

        Array<OneD, NekDouble> TmapGrad(m_spacedim * nq);
        TmapGrad = ComputeCovGrad(TimeMap_old[0], m_TMmovingframes);

        Array<OneD, NekDouble> TmapGradMag(nq);
        TmapGradMag = ComputeVelocityMag(TmapGrad);


            // Array<OneD, NekDouble> TMVelocity(m_spacedim * nq);
            // TMVelocity = ComputeVelocityTimeMap(m_ValidTimeMap, TimeMap_old);

            // Array<OneD, NekDouble> TMVelMag = ComputeVelocityMag(TMVelocity);


        // for (int i=0; i<nq; ++i)
        // {
        //     // if( (x0[i]<10.0) && (x0[i]>-10.0) )
        //     {
        //         std::cout << "i = " << i << ", TMgradmag = " << TmapGradMag[i] << ", velmag = " << Velocitymag_old[i] << ", x = " << x0[i] << std::endl;
        //     }
        // }

        // NekDouble velmagmax = Vmath::Vmax(nq, Velocitymag_old, 1);
        // int velmagIndex = Vmath::Imax(nq, Velocitymag_old, 1);

        // std::cout << ", velmag max = " << velmagmax << ", x = " << x0[velmagIndex] << std::endl;
        // wait_on_enter();


        // load new timemap                               
        std::string loadname = sessionnew + "_TimeMap_" +
                        boost::lexical_cast<std::string>(m_TimeMapnstep) + ".chk";

        LoadTimeMap(loadname, m_ValidTM, TimeMap_new, AniStrength_new, Velocitymag_new, Velocity_new);

        // Computed the deformed timemap from the old timemap
        Array<OneD, Array<OneD, NekDouble>> Velocity_deformed(m_spacedim);
        ComputeVelocityDeformed(AniStrength_old[0], Velocity_old, AniStrength_new[0], Velocity_deformed);
        
        Array<OneD, Array<OneD, NekDouble>> Vdiff;
        Array<OneD, NekDouble> TimeMapDiff = HelmSolveTimeMapDiff(Velocity_old, Velocity_deformed, Vdiff, VdiffDivergence, VdiffMag);

        PlotDeformedTimeMap(TimeMap_old[0], TimeMap_new[0], TimeMapDiff, Velocity_old, VdiffMag, VdiffDivergence);

        std::cout << "TimeMapExact[0] = " << TimeMap_old[0][0] - TimeMap_new[0][0] << ", TimeMapExact[end] = " << TimeMap_old[0][nq-1] - TimeMap_new[0][nq-1] << std::endl;

        return TimeMapDiff;
    }


    Array<OneD, NekDouble> MMFCardiacEP::HelmSolveTimeMapDiff(
                const Array<OneD, const Array<OneD, NekDouble>> &Velocity_old,
                const Array<OneD, const Array<OneD, NekDouble>> &Velocity_deformed,
                Array<OneD, Array<OneD, NekDouble>> &Vdiff,
                Array<OneD, NekDouble> &VdiffDivergence,
                Array<OneD, NekDouble> &VdiffMag)
    {
        int nq      = GetNpoints();

        Array<OneD, NekDouble> outarray(nq);

        // Compute Vdiff = V_new/Vmag_new^2 - V_old/Vmag_old^2
        Vdiff = Array<OneD, Array<OneD, NekDouble>>(m_spacedim);
        ComputeDeformedDiff(Velocity_old, Velocity_deformed, Vdiff);

        VdiffMag = ComputeVelocityMag(Vdiff);
        std::cout << "Vdiffmag = [ " << Vmath::Vmin(nq, VdiffMag, 1) << " , " << Vmath::Vmax(nq, VdiffMag, 1) << " ] " << std::endl;

        // Compute Euclidean divergence of Vdiff
        VdiffDivergence = ComputeEuclideanDivergence(Vdiff);

        std::cout << "VdiffDivergence = [ " << Vmath::Vmin(nq, VdiffDivergence, 1) 
        << " , " << Vmath::Vmax(nq, VdiffDivergence, 1) << " ] " << std::endl << std::endl;

        // Solve the Poisson equation: \nabla^2 T_{diff} = \nabla \cdot \mathbf{V}_{diff}
        StdRegions::ConstFactorMap factors;
        factors[StdRegions::eFactorTau]    = m_Helmtau;
        factors[StdRegions::eFactorLambda] = 0.0;

        Vmath::Sadd(nq, -1.0 * AvgInt(VdiffDivergence), VdiffDivergence, 1, VdiffDivergence, 1);
        std::cout << "VdiffDiv for pure Neumann = " << AvgInt(VdiffDivergence) << std::endl;
        Vmath::Vcopy(nq, VdiffDivergence, 1, m_fields[0]->UpdatePhys(), 1);

        m_fields[0]->HelmSolve(m_fields[0]->GetPhys(), m_fields[0]->UpdateCoeffs(), factors, m_varcoeff);
        m_fields[0]->BwdTrans(m_fields[0]->GetCoeffs(), m_fields[0]->UpdatePhys());
        m_fields[0]->SetPhysState(true);

        outarray = m_fields[0]->GetPhys();

        std::cout << "outarray[0] = " << outarray[0] << ", outarray[end] = " << outarray[nq-1] << ", outarray_diff = " << outarray[0] - outarray[nq-1] <<  std::endl;

        // Let the Time Map at the left end should be zero
        Vmath::Sadd(nq, -1.0 * outarray[0], outarray, 1, outarray, 1);

        return outarray;
    }


    void MMFCardiacEP::ComputeDeformedDiff(const Array<OneD, const Array<OneD, NekDouble>> &Velocity_old,
                                           const Array<OneD, const Array<OneD, NekDouble>> &Velocity_deformed,
                                           Array<OneD, Array<OneD, NekDouble>> &Vdiff)
    {  
        int nq      = GetNpoints();

        NekDouble Tol=0.001;
        NekDouble Velmag_defm, Velmag_old;

        Array<OneD, NekDouble> Velocitymag_old = ComputeVelocityMag(Velocity_old);
        Array<OneD, NekDouble> Velocitymag_deformed = ComputeVelocityMag(Velocity_deformed);

        for (int k=0;k<m_spacedim;++k)
        {
            Vdiff[k] = Array<OneD, NekDouble>(nq,0.0);
            for (int i=0;i<nq; ++i)
            {
                Velmag_old = Velocitymag_old[i];
                Velmag_defm = Velocitymag_deformed[i];
                if( (Velmag_defm>Tol) && (Velmag_old>Tol))
                {
                     Vdiff[k][i] = Velocity_deformed[k][i]/(Velmag_defm*Velmag_defm) - Velocity_old[k][i]/(Velmag_old*Velmag_old);
                }
            }
        }
    }

void MMFCardiacEP::ComputeVelocityDeformed(const Array<OneD, const NekDouble> &AniStrength_old,
                            const Array<OneD, const Array<OneD, NekDouble>> &Velocity_old,
                            const Array<OneD, const NekDouble> &AniStrength_new,
                            Array<OneD, Array<OneD, NekDouble>> &Velocity_deformed)

    {
            int nq      = GetNpoints();

            Array<OneD, Array<OneD, NekDouble>> AniStrength_deformed(m_expdim);
            for (int j = 0; j < m_expdim; ++j)
            {
                AniStrength_deformed[j] = Array<OneD, NekDouble>(nq, 1.0);
            }

            Array<OneD, Array<OneD, NekDouble>> TimeMap_deformed(1);
            for (int i = 0; i < 1; ++i)
            {
                TimeMap_deformed[i] = Array<OneD, NekDouble>(nq,0.0);
            }

            // convert vectoer coeffcients to vector
            Array<OneD, Array<OneD, NekDouble>> vcoeff_old;
            vector_to_vcoeff(Velocity_old, m_movingframes, vcoeff_old);

            std::cout << "v1 = [ " << Vmath::Vmin(nq, vcoeff_old[0], 1) << " , " << Vmath::Vmax(nq, vcoeff_old[0], 1)
            << " ], v2 = [ " << Vmath::Vmin(nq, vcoeff_old[1], 1) << " , " << Vmath::Vmax(nq, vcoeff_old[1], 1) << " ]" << std::endl;

            Array<OneD, Array<OneD, NekDouble>> vcoeff_deformed(m_expdim);
            for (int j = 0; j < m_expdim; ++j)
            {
                vcoeff_deformed[j] = Array<OneD, NekDouble>(nq);
            }

            Array<OneD, NekDouble> aniratio(nq);
            Vmath::Vdiv(nq, AniStrength_new, 1, AniStrength_old, 1, aniratio, 1);
            Vmath::Vsqrt(nq, aniratio, 1, aniratio, 1);

            std::cout << "aniratio = [ " << Vmath::Vmin(nq, aniratio, 1) << " , " 
            << Vmath::Vmax(nq, aniratio, 1) << " ]" << std::endl;

            // \mathbf{v}^{deformed} = sqrt(d_{new}/d_{old}) v_1 \mathbf{e}^1 + v_2 \mathbf{e}^2
            Vmath::Vmul(nq, aniratio, 1, vcoeff_old[0], 1, vcoeff_deformed[0], 1);
            Vmath::Vcopy(nq, vcoeff_old[1], 1, vcoeff_deformed[1], 1);

            vcoeff_to_vector(vcoeff_deformed, m_movingframes, Velocity_deformed);
    }

void MMFCardiacEP::LoadTimeMap(std::string &loadname,
                                const Array<OneD, const NekDouble> &ValidTM,
                                Array<OneD, Array<OneD, NekDouble>> &TimeMap,
                                Array<OneD, Array<OneD, NekDouble>> &AniStrength,
                                Array<OneD, NekDouble> &Velocitymag,
                                Array<OneD, Array<OneD, NekDouble>> &Velocity)
{
    int nvar    = 6;
    int nq      = GetNpoints();
    int ncoeffs = m_fields[0]->GetNcoeffs();
    std::vector<std::string> variables(nvar);

    variables[0] = "TimeMap";
    variables[1] = "AniStrength";
    variables[2] = "velmag";
    variables[3] = "velx";
    variables[4] = "vely";
    variables[5] = "velz";
        
    Array<OneD, Array<OneD, NekDouble>> tmpc(nvar);
    for (int i = 0; i < nvar; ++i)
    {
       tmpc[i]  = Array<OneD, NekDouble>(ncoeffs,0.0);
    }

    EquationSystem::ImportFld(loadname, variables, tmpc);

    // Time Map
    TimeMap = Array<OneD, Array<OneD, NekDouble>>(1);
    for (int i = 0; i < 1; ++i)
    {
        TimeMap[i] = Array<OneD, NekDouble>(nq,0.0);
    }
    m_fields[0]->BwdTrans(tmpc[0], TimeMap[0]);

    // AniStrength
    AniStrength = Array<OneD, Array<OneD, NekDouble>>(m_expdim);
    for (int j = 0; j < m_expdim; ++j)
    {
        AniStrength[j] = Array<OneD, NekDouble>(nq, 1.0);
    }

    m_fields[0]->BwdTrans(tmpc[1], AniStrength[0]);

    // Time Map Velocity Magnitude 
    m_fields[0]->BwdTrans(tmpc[2], Velocitymag);

    // TMvelocity = Array<OneD, Array<OneD, NekDouble>>(3);
    for (int i=0; i<m_spacedim; ++i)
    {
        Velocity[i] = Array<OneD, NekDouble>(nq);
        m_fields[0]->BwdTrans(tmpc[i+3], Velocity[i]);
    }

    // Realign Time Map by T0
    m_session->LoadParameter("TimeMapDelay", m_TimeMapDelay, 5.0);
    if(m_SolverSchemeType==eTimeMapMarch)
    {
        for (int i = 0; i < nq; ++i)
        {
            if (TimeMap[0][i] > m_TimeMapDelay)
            {
                TimeMap[0][i] = TimeMap[0][i] - m_TimeMapDelay;
            }
        }
    }

    // Apply the Valid region
    Vmath::Vmul(nq, ValidTM, 1, TimeMap[0], 1, TimeMap[0], 1);
    Vmath::Vmul(nq, ValidTM, 1, AniStrength[0], 1, AniStrength[0], 1);
    Vmath::Vmul(nq, ValidTM, 1, Velocitymag, 1, Velocitymag, 1);
    for (int i=0; i<m_spacedim; ++i)
    {
        Vmath::Vmul(nq, ValidTM, 1, Velocity[i], 1, Velocity[i], 1);
    }

    std::cout << "TimeMap Mag = " << RootMeanSquare(TimeMap[0]) << std::endl ;

    std::cout << "TimeMapName = " << loadname << std::endl;
    std::cout << "TimeMap = [ " << Vmath::Vmin(nq, TimeMap[0], 1)  << " , " 
    << Vmath::Vmax(nq, TimeMap[0], 1) << " ] " << std::endl;
    std::cout << "AniStrength = [ " << Vmath::Vmin(nq, AniStrength[0], 1)  << " , " 
    << Vmath::Vmax(nq, AniStrength[0], 1) << " ] " << std::endl;
    std::cout << "Velocitymag = [ " << Vmath::Vmin(nq, Velocitymag, 1)  << " , " 
    << Vmath::Vmax(nq, Velocitymag, 1) << " ] " << std::endl << std::endl;
}

void MMFCardiacEP::LoadCardiacFiber(
    const MediumType CardiacMediumType,
    const NekDouble AnisotropyStrength,
    Array<OneD, Array<OneD, NekDouble>> &AniStrength,
    Array<OneD, NekDouble> &CardiacFibre)
{
    int nq = m_fields[0]->GetNpoints();
    
    switch (CardiacMediumType)
    {
        case eAnisotropy:
        {
            for (int i=0; i<nq; ++i)
            {
                AniStrength[0][i] = sqrt(AnisotropyStrength);            
            }
        }
        break;

        case eAnisotropyFiberMap:
        {
            m_ImportedFiberExist = 1;
            AniStrength[0] = ReadFibermap(AnisotropyStrength, CardiacFibre);
        }
        break;

        case eHeterogeneousIsotropy:
        {
            Array<OneD, NekDouble> tmp = ReadConductivityMap();
            Vmath::Vcopy(nq, &tmp[0], 1, &AniStrength[0][0], 1);
            Vmath::Vcopy(nq, &tmp[0], 1, &AniStrength[1][0], 1);
        }
        break;

        case eHeterogeneousAnisotropy:
        {
            m_ImportedFiberExist = 1;
            AniStrength[0] = ReadFibermap(AnisotropyStrength, CardiacFibre);

            Array<OneD, NekDouble> tmp = ReadConductivityMap();
            Vmath::Vmul(nq, &tmp[0], 1, &AniStrength[0][0], 1,
                        &AniStrength[0][0], 1);
        }
        break;

        case eRegionalHeterogeneous:
        {
            int index;
            for (int i = m_AniRegionStart; i < m_AniRegionEnd; ++i)
                {
                    for (int j = 0; j < m_fields[0]->GetTotPoints(i); ++j)
                        {
                            index = m_fields[0]->GetPhys_Offset(i) + j;
                            AniStrength[0][index] = sqrt(AnisotropyStrength);            
                        }
                }
        }
        break;

        case eRegionalIsotropy:
        {
            int index;
            for (int i = m_AniRegionStart; i < m_AniRegionEnd; ++i)
                {
                    for (int j = 0; j < m_fields[0]->GetTotPoints(i); ++j)
                        {
                            index = m_fields[0]->GetPhys_Offset(i) + j;
                            AniStrength[0][index] = sqrt(AnisotropyStrength);
                            AniStrength[1][index] = sqrt(AnisotropyStrength);                        
                        }
                }
        }
        break;

        case eGaussianHeterogeneous:
        {
            Array<OneD, NekDouble> x0(nq);
            Array<OneD, NekDouble> x1(nq);
            Array<OneD, NekDouble> x2(nq);

            m_fields[0]->GetCoords(x0, x1, x2);

            NekDouble tau;
            for (int i=0; i<nq; ++i)
            {
                tau = x0[i]/m_Gaussiantau;
                AniStrength[0][i] = 4.0 - AnisotropyStrength * exp(-0.5 * tau * tau);
            }
        }
        break;

        default:
            break;
    }

    std::cout << "Anisotropy = [ " << Vmath::Vmin(nq, AniStrength[0], 1) << " , " << Vmath::Vmax(nq, AniStrength[0], 1) << " ] " << std::endl;

    // Plot Cardiac fibre projection map
    // PlotProcessedCardiacFibre(movingframes[0], fcdotk, AniConstruction);
}

// Read the fiber and maintain the same magnitude for all the fiber
Array<OneD, NekDouble> MMFCardiacEP::ReadFibermap(
    const NekDouble AnisotropyStrength, Array<OneD, NekDouble> &CardiacFibre)
{
    int nq = m_fields[0]->GetNpoints();

    cout << "Loading Anisotropic Fibre map ===========" << endl;
    Array<OneD, NekDouble> outarray(nq);

    std::string anisotropy[3] = {"fx", "fy", "fz"};
    CardiacFibre              = Array<OneD, NekDouble>(m_spacedim * nq, 0.0);

    Array<OneD, NekDouble> tmp(nq);
    for (int i = 0; i < m_spacedim; ++i)
    {
        GetFunction("AnisotropicConductivity")->Evaluate(anisotropy[i], tmp);
        Vmath::Vcopy(nq, &tmp[0], 1, &CardiacFibre[i * nq], 1);
    }

    PlotCardiacFibre(CardiacFibre);

    // If there is a fibre, let it be with a strength of m_AnisotropyStrength;
    NekDouble Tol = 1.0e-4;
    NekDouble mag, fx, fy, fz;
    for (int k = 0; k < nq; ++k)
    {
        fx = CardiacFibre[k];
        fy = CardiacFibre[k + nq];
        fz = CardiacFibre[k + 2 * nq];

        mag = sqrt(fx * fx + fy * fy + fz * fz);

        if (mag > Tol)
        {
            outarray[k] = AnisotropyStrength;
        }
    }

    return outarray;
}

Array<OneD, NekDouble> MMFCardiacEP::ReadConductivityMap()
{
    int nq = m_fields[0]->GetNpoints();

    cout << "Loading Isotropic Conductivity map." << endl;

    std::string varName = "intensity";

    Array<OneD, NekDouble> vTemp;
    GetFunction("IsotropicConductivity")->Evaluate(varName, vTemp);

    NekDouble f_min, f_max;
    m_session->LoadParameter("d_min", f_min, -102.0);
    m_session->LoadParameter("d_max", f_max, 32.1);
    const NekDouble scar_min = 0.1;
    const NekDouble scar_max = 1.0;

    // Threshold based on d_min, d_max
    for (int j = 0; j < nq; ++j)
    {
        vTemp[j] = (vTemp[j] < f_min ? f_min : vTemp[j]);
        vTemp[j] = (vTemp[j] > f_max ? f_max : vTemp[j]);
    }

    std::cout << "vTemp: Max = " << Vmath::Vmax(nq, vTemp, 1)
              << ", min =" << Vmath::Vmin(nq, vTemp, 1) << std::endl;

    // Rescale to s \in [0,1] (0 maps to d_max, 1 maps to d_min)
    Vmath::Sadd(nq, -f_min, vTemp, 1, vTemp, 1);
    Vmath::Smul(nq, -1.0 / (f_max - f_min), vTemp, 1, vTemp, 1);
    Vmath::Sadd(nq, 1.0, vTemp, 1, vTemp, 1);
    Vmath::Smul(nq, scar_max - scar_min, vTemp, 1, vTemp, 1);
    Vmath::Sadd(nq, scar_min, vTemp, 1, vTemp, 1);

    return vTemp;
}

void MMFCardiacEP::DoOdeProjection(
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

void MMFCardiacEP::v_DoSolve()
{
    switch (m_SolverSchemeType)
    {
        case eMMFZero:
        {
            DoSolveMMF();
        }
        break;

        case eMMFFirst:
        {
            DoSolveMMFFirst();
        }
        break;

        case eTimeMapMarch:
        {
            DoSolveTimeMap();
        }
        break;

        default:
        break;
    }
}

void MMFCardiacEP::DoSolveTimeMap()
{
    ASSERTL0(m_intScheme != 0, "No time integration scheme.");

    int i, nchk = 1;
    int nvariables = 0;
    int nfields    = m_fields.size();
    int nq         = m_fields[0]->GetNpoints();

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
    for (i = 0; i < nvariables; ++i)
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

    Array<OneD, int> TMcount(nq, 0);
    while (step < m_steps || m_time < m_fintime - NekConstants::kNekZeroTol)
    {
        timer.Start();
        fields = m_intScheme->TimeIntegrate(step, m_timestep, m_ode);
        timer.Stop();

        // Excitation according to TimeMap
        for (int i = 0; i < nq; ++i)
        {
            if ((m_TimeMap[0][i] <= m_time) && (m_TimeMap[0][i] > m_time - m_timestep))
            {
                fields[0][i] += m_TimeMapIapp;
                TMcount[i] += 1;
            }
        }

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

            std::cout << "TMcount = " << Vmath::Vsum(nq, TMcount, 1) << " / "
                      << nq << ", field = " << RootMeanSquare(fields[0]) << std::endl;

            int Iumax = Vmath::Imax(nq, fields[0], 1);
            std::cout << "u_max = " << Vmath::Vmax(nq, fields[0], 1)
                      << " at x = " << x0[Iumax] << ", y = " << x1[Iumax] << ", z = " << x2[Iumax]
                      << std::endl;

            int Iumin = Vmath::Imin(nq, fields[0], 1);
            std::cout << "u_min = " << Vmath::Vmin(nq, fields[0], 1)
                      << " at x = " << x0[Iumin] << ", y = " << x1[Iumin] << ", z = " << x2[Iumin]
                      << std::endl;

            std::stringstream ss;
            ss << cpuTime / 60.0 << " min.";
            std::cout << " CPU Time: " << std::setw(8) << std::left << ss.str()
                      << std::endl;

            cpuTime = 0.0;
        }

        // Transform data into coefficient space
        // for (i = 0; i < nvariables; ++i)
        // {
        //     m_fields[m_intVariables[i]]->SetPhys(fields[i]);
        //     m_fields[m_intVariables[i]]->FwdTrans_IterPerExp(
        //         fields[i], m_fields[m_intVariables[i]]->UpdateCoeffs());
        //     m_fields[m_intVariables[i]]->SetPhysState(false);
        // }

        if ((m_checksteps && step && !((step + 1) % m_checksteps)) ||
            doCheckTime)
        {
            ComputeTimeMapError(nchk, fields);

            Checkpoint_Output(nchk++);
            doCheckTime = false;
        }

        ++step;
    } // namespace Nektar

    if (m_session->GetComm()->GetRank() == 0)
    {
        std::cout << "Time-integration  : " << intTime << "s" << std::endl;
    }
 
    for (i = 0; i < nvariables; ++i)
    {
        m_fields[m_intVariables[i]]->SetPhys(fields[i]);
        m_fields[m_intVariables[i]]->SetPhysState(true);
    }

} // namespace Nektar

void MMFCardiacEP::DoSolveMMFFirst()
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
    Array<OneD, NekDouble> dudtMax(nq, 0.0);
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
    Array<OneD, NekDouble> VelMap(m_spacedim * nq, 0.0);

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
                           dudtvalHistory, dudtMax, IappMap, TimeMap, m_TimeMapScheme);
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
            }

            if ( (RootMeanSquare(TimeMap)>1.0) )
            {
                PlotTimeMap(m_ValidTimeMap, m_AniStrength, TimeMap, nchk);
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
        m_fields[i]->FwdTrans(m_fields[i]->GetPhys(),
                              m_fields[i]->UpdateCoeffs());
    }
} // namespace Nektar

void MMFCardiacEP::DoSolveMMF()
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

    Array<OneD, NekDouble> TimeMap(nq, 0.0);
    Array<OneD, NekDouble> IappMap(nq, 0.0);

    Array<OneD, NekDouble> dudtval(nq);
    Array<OneD, NekDouble> dudtvalHistory(nq, 0.0);
    Array<OneD, NekDouble> dudtMax(nq, 0.0);

    // Aligh Moving Frames along the velocit vector
    Array<OneD, Array<OneD, NekDouble>> MF1st(m_spacedim);
    for (int i = 0; i < m_spacedim; ++i)
    {
        MF1st[i]    = Array<OneD, NekDouble>(m_spacedim * nq);
    }

    Array<OneD, NekDouble> x0(nq);
    Array<OneD, NekDouble> x1(nq);
    Array<OneD, NekDouble> x2(nq);

    m_fields[0]->GetCoords(x0, x1, x2);

    std::cout << "xmax = " << Vmath::Vmax(nq, x0, 1) << ", xmin = " << Vmath::Vmin(nq, x0, 1)
    << ", ymax = " << Vmath::Vmax(nq, x1, 1) << ", ymin = " << Vmath::Vmin(nq, x1, 1) << std::endl;

    while (step < m_steps || m_time < m_fintime - NekConstants::kNekZeroTol)
    {
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

        if ((m_TimeMapStart <= m_time) && (m_TimeMapEnd >= m_time))
        {
            ComputeTimeMap(m_time, m_urest, fields[0], dudtval, m_ValidTimeMap,
                        dudtvalHistory, dudtMax, IappMap, TimeMap, m_TimeMapScheme);
        }

        if (m_session->GetComm()->GetRank() == 0 && !((step + 1) % m_infosteps))
        {
            std::cout << "Steps: " << std::setw(8) << std::left << step + 1
                      << " "
                      << "Time: " << std::setw(12) << std::left << m_time
                      << std::endl;

            std::stringstream ss;
            ss << cpuTime / 60.0 << " min.";
            std::cout << " CPU Time: " << std::setw(8) << std::left << ss.str()
                      << std::endl;

            cpuTime = 0.0;
        }

        // Write out checkpoint files
        if ((m_checksteps && step && !((step + 1) % m_checksteps)) ||
            doCheckTime)
        {
            int Iumax = Vmath::Imax(nq, fields[0], 1);
            std::cout << "u_max = " << Vmath::Vmax(nq, fields[0], 1)
                      << " at x = " << x0[Iumax] << ", y = " << x1[Iumax] << ", z = " << x2[Iumax]
                      << std::endl;

            int Iumin = Vmath::Imin(nq, fields[0], 1);
            std::cout << "u_min = " << Vmath::Vmin(nq, fields[0], 1)
                      << " at x = " << x0[Iumin] << ", y = " << x1[Iumin] << ", z = " << x2[Iumin]
                      << std::endl;
                      
            PlotTimeMap(m_ValidTimeMap, m_AniStrength, TimeMap, nchk);

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
} // namespace Nektar

void MMFCardiacEP::DoImplicitSolveCardiacEP(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time,
    const NekDouble lambda)
{
    boost::ignore_unused(time);

    int nvar = inarray.size();
    int nq   = m_fields[0]->GetNpoints();

    // Set up factors for Helmsolver
    StdRegions::ConstFactorMap factors;
    factors[StdRegions::eFactorTau] = m_Helmtau;

    // 	factors[StdRegions::eFactorLambda] = 1.0 / lambda * m_chi *
    // m_capMembrane;
    factors[StdRegions::eFactorLambda] = 1.0 / lambda;

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

        // m_fields[0]->HelmSolve(m_fields[0]->GetPhys(), m_fields[0]->UpdateCoeffs(),
        //                    factors);
    m_fields[0]->BwdTrans(m_fields[0]->GetCoeffs(), outarray[0]);
    m_fields[0]->SetPhysState(true);

    // No diffusion for the second variable
    for (int i = 1; i < nvar; ++i)
    {
        Vmath::Vcopy(nq, &inarray[i][0], 1, &outarray[i][0], 1);
    }
}

// We Return Y[i] = rhs [i] without no Helomsolver
void MMFCardiacEP::DoNullSolve(
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

void MMFCardiacEP::DoOdeRhsCardiacEPTimeMap(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time)
{
    // Compute the reaction function
    // input: inarray
    // output: outarray
    m_cell->TimeIntegrate(inarray, outarray, time);

    // AlievPanfilovReaction(inarray, outarray, time);

    // std::cout << "inarray0 = " << RootMeanSquare(inarray[0]) << ", outarray = " << RootMeanSquare(outarray[0]) << std::endl;
    // std::cout << "inarray1 = " << RootMeanSquare(inarray[1]) << ", outarray = " << RootMeanSquare(outarray[1]) << std::endl << std::endl;

    // Compute I_stim
    for (unsigned int i = 0; i < m_stimulus.size(); ++i)
    {
        m_stimulus[i]->Update(outarray, time);
    }
 }

void MMFCardiacEP::DoOdeRhsCardiacEP(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time)
{
    int nq = m_fields[0]->GetNpoints();

    // Compute the reaction function
    // input: inarray
    // output: outarray
    m_cell->TimeIntegrate(inarray, outarray, time);

    // Compute I_stim
    for (unsigned int i = 0; i < m_stimulus.size(); ++i)
    {
        m_stimulus[i]->Update(outarray, time);
    }

    if (m_explicitDiffusion)
    {
        // Laplacian only to the first variable
        Array<OneD, NekDouble> Laplacian(nq);
        WeakDGMMFDiffusion(0, inarray[0], Laplacian, time);

        Vmath::Vadd(nq, &Laplacian[0], 1, &outarray[0][0], 1, &outarray[0][0], 1);
    }
}

void MMFCardiacEP::v_SetInitialConditions(NekDouble initialtime,
                                          bool dumpInitialConditions,
                                          const int domain)
{
    boost::ignore_unused(domain, dumpInitialConditions);

    int nvar = m_fields.size();
    int nq = GetTotPoints();
    int ncoeffs = m_fields[0]->GetNcoeffs();

    m_cell->Initialise();

    // Read initial condition from xml file
    EquationSystem::v_SetInitialConditions(initialtime, false);

    Array<OneD, Array<OneD, NekDouble>> tmp(1);
    tmp[0] = Array<OneD, NekDouble>(nq);

    Array<OneD, Array<OneD, NekDouble>> initialcondition(nvar);
    for (int i=0; i<nvar; ++i)
    {
        initialcondition[i] = Array<OneD, NekDouble>(nq,0.0);
    }

    for (unsigned int i = 0; i < m_stimulus.size(); ++i)
    {
        m_stimulus[i]->Update(initialcondition, 0.1);
    }

    Vmath::Vcopy(nq, m_fields[0]->GetPhys(), 1, tmp[0], 1);
    Vmath::Vadd(nq, tmp[0], 1, initialcondition[0], 1, initialcondition[0], 1);

        std::string outname;
        outname = m_sessionName + "_initialcondition.chk";

        std::vector<Array<OneD, NekDouble>> fieldcoeffs(1);
        for (int i = 0; i < 1; ++i)
        {
            fieldcoeffs[i] = Array<OneD, NekDouble>(ncoeffs);
        }

        std::vector<std::string> variables(1);
        variables[0] = "initialcondition";
        
        m_fields[0]->FwdTrans(initialcondition[0], fieldcoeffs[0]);

        WriteFld(outname, m_fields[0], fieldcoeffs, variables);

        m_ValidTimeMap = ComputeTimeMapInitialZone(m_urest, initialcondition[0]);

    // Only the excited regions are considered for m_InitExcitation
    // Vmath::Vsub(nq, tmp[0], 1, initialcondition, 1, initialcondition, 1);
    if (m_SolverSchemeType == eTimeMapMarch)
    {
        // std::cout << "PlotTimeEnergyMap starts ==========================" << std::endl;

        // Array<OneD, NekDouble> VelVector = ComputeVelocityTimeMap(m_ValidTimeMap, m_TimeMap[0]);
 
        // // \nabla \cdot \ell = - \nabla^2 U_Lamb
        // Array<OneD, NekDouble> LambDiv = ComputeLambDiv(m_ValidTimeMap, VelVector);
 
        // // Ion Potential
        // Array<OneD, NekDouble> IonE = HelmsolvePotentialE(m_ValidTimeMap, LambDiv);
 
        // // Plot Energy map
        // PlotTimeEnergyMap(m_TimeMap[0], VelVector, LambDiv, IonE);

        // std::cout << "PlotTimeEnergyMap ends ==========================" << std::endl;
    }

    else
    {
        for (unsigned int i = 0; i < m_stimulus.size(); ++i)
        {
            m_stimulus[i]->Update(tmp, 0.1);
            m_fields[0]->SetPhys(tmp[0]);
        }
    }

    // forward transform to fill the modal coeffs
    for (int i = 0; i < m_fields.size(); ++i)
    {
        m_fields[i]->SetPhysState(true);
        m_fields[i]->FwdTrans(m_fields[i]->GetPhys(),
                              m_fields[i]->UpdateCoeffs());
    }

    std::cout << "Initial: u = " << RootMeanSquare(m_fields[0]->GetPhys()) << ", max u = "
              << Vmath::Vmax(nq, m_fields[0]->GetPhys(), 1) << std::endl;

    if (dumpInitialConditions)
    {
        std::string outname;
        outname = m_sessionName + "_initial.chk";

        WriteFld(outname);
    }
}


Array<OneD, int> MMFCardiacEP::ComputeTimeMapInitialZone(
    const NekDouble urest,
    const Array<OneD, const NekDouble> &inarray)
{
    int nq = GetTotPoints();

    Array<OneD, int> outarray(nq, 1);

    int TMflag, index;
    int cnt             = 0;
    const NekDouble Tol = 0.1;
    NekDouble diff;
    for (int i = 0; i < m_fields[0]->GetExpSize(); ++i)
    {
        TMflag = 0;
        for (int j = 0; j < m_fields[0]->GetTotPoints(i); ++j)
        {
            index = m_fields[0]->GetPhys_Offset(i) + j;
            diff = inarray[index] - urest;
            if ((fabs(diff) > Tol) || (m_MMFActivation[index] == 0))
            {
                TMflag = 1;
            }
        }

        if (TMflag == 1)
        {
            for (int j = 0; j < m_fields[0]->GetTotPoints(i); ++j)
            {
                index           = m_fields[0]->GetPhys_Offset(i) + j;
                outarray[index] = 0;
                cnt++;
            }
        }
    }

    std::cout << "ValidTimeMap = " << cnt << " / " << nq << " ( " << (1.0*cnt/nq)*100.0 << " % ) is De-Acitvated"
              << std::endl;

    return outarray;
}

void MMFCardiacEP::ComputeTimeMapError(const int TMnstep, const Array<OneD, const Array<OneD, NekDouble>> &outfield)
{
    int nvar    = 1;
    int nq      = m_fields[0]->GetNpoints();
    int ncoeffs = m_fields[0]->GetNcoeffs();

    std::vector<std::string> variables(nvar);
    variables[0] = "u";

    std::string loadname = m_TMsessionName + "_" + 
                        boost::lexical_cast<std::string>(TMnstep) + ".chk";

    Array<OneD, Array<OneD, NekDouble>> tmpc(nvar);
    Array<OneD, Array<OneD, NekDouble>> uexact(nvar);
    for (int i = 0; i < nvar; ++i)
    {
        tmpc[i]   = Array<OneD, NekDouble>(ncoeffs);
        uexact[i] = Array<OneD, NekDouble>(nq);
    }

    EquationSystem::ImportFld(loadname, variables, tmpc);
 
    for (int i = 0; i < nvar; ++i)
    {
        m_fields[0]->BwdTrans(tmpc[i], uexact[i]);
    }
 
    Array<OneD, NekDouble> udiff(nq, 0.0);
    Vmath::Vsub(nq, outfield[0], 1, uexact[0], 1, udiff, 1);
    Vmath::Vabs(nq,  udiff, 1,  udiff, 1);

    for (int i=0;i<nq;++i)
    {
        if(m_ValidTimeMap[i]==0)
        {
            udiff[i] = 0.0;
        }
    }

    // Array<OneD, NekDouble> vdiff(nq, 0.0);
    // Array<OneD, NekDouble> tmp = m_cell->GetCellSolution(1);
    // Vmath::Vsub(nq, tmp, 1, uexact[1], 1, vdiff, 1);

    NekDouble L2uerr = RootMeanSquare(udiff) / Vmath::Vamax(nq, uexact[0], 1);

    NekDouble Linfuerr = Vmath::Vamax(nq, udiff, 1) / Vmath::Vamax(nq, uexact[0], 1);
    std::cout << " TimeMap: uExact = " << RootMeanSquare(uexact[0]) << ", u_Error: L2 = " << L2uerr << ", Linf = " << Linfuerr << std::endl;

    PlotTimeMapErr(m_TimeMap, m_ValidTimeMap, outfield[0], uexact[0], udiff, TMnstep);
}

void MMFCardiacEP::PlotTimeMapErr(
                                 const Array<OneD, const Array<OneD, NekDouble>> TimeMap,
                                 const Array<OneD, const int> ValidTimeMap,
                                const Array<OneD, const NekDouble> &field,
                                const Array<OneD, const NekDouble> &uexact,
                                const Array<OneD, const NekDouble> &udiff,
                                const int nstep)
{
    int nvar    = 5;
    int ncoeffs = m_fields[0]->GetNcoeffs();
    int nq = m_fields[0]->GetTotPoints();

    std::string outname1 = m_sessionName + "_TimeMapErr_" +
                           boost::lexical_cast<std::string>(nstep) + ".chk";

    std::vector<Array<OneD, NekDouble>> fieldcoeffs(nvar);
    for (int i = 0; i < nvar; ++i)
    {
        fieldcoeffs[i] = Array<OneD, NekDouble>(ncoeffs);
    }

    std::vector<std::string> variables(nvar);
    variables[0] = "TimeMap";
    variables[1] = "ValidTimeMap";
    variables[2] = "u_TimeMap";
    variables[3] = "u_PDE";
    variables[4] = "u_err";

    Array<OneD, NekDouble> ValidTM(nq);
    for( int i=0;i<nq; i++)
    {
        ValidTM[i] = 1.0*ValidTimeMap[i];
    }

    m_fields[0]->FwdTrans(TimeMap[0], fieldcoeffs[0]);
    m_fields[0]->FwdTrans(ValidTM, fieldcoeffs[1]);
    m_fields[0]->FwdTrans(field, fieldcoeffs[2]);
    m_fields[0]->FwdTrans(uexact, fieldcoeffs[3]);
    m_fields[0]->FwdTrans(udiff, fieldcoeffs[4]);

    WriteFld(outname1, m_fields[0], fieldcoeffs, variables);
}

// Compute Velocity field from Time Map
// flag = 1: Use \vec{v} = \nabla T / \| \nabla T \|^2
// flag = 0:
Array<OneD, NekDouble> MMFCardiacEP::ComputeVelocityTimeMap(
    const Array<OneD, const int> &ValidTimeMap,
    const Array<OneD, const NekDouble> &inarray)
{
    int nq = m_fields[0]->GetTotPoints();

    Array<OneD, NekDouble> outarray(m_spacedim * nq, 0.0);

    Array<OneD, NekDouble> physarray(nq);
    Vmath::Vcopy(nq, inarray, 1, physarray, 1);

    // TmapGrad = \nabla Tmap
    Array<OneD, NekDouble> TmapGrad(m_spacedim * nq);
    TmapGrad = ComputeCovGrad(physarray, m_TMmovingframes);

    Array<OneD, NekDouble> TmapGradMag(nq);
    TmapGradMag = ComputeVelocityMag(TmapGrad);

    // Compute VelField \vec{v} = \sum_{i=1}^3 1/(\nabla T \cdot \hat{x}_i)
    // \hat{x}_i
    NekDouble tmp, TmapGradTol = 0.1;
    for (int i = 0; i < nq; i++)
    {
        tmp = TmapGradMag[i];
        if (tmp > TmapGradTol)
        {
            for (int k = 0; k < m_spacedim; ++k)
            {
                outarray[i + k * nq] = TmapGrad[i + k * nq] / (tmp * tmp);
            }
        }

        if (ValidTimeMap[i] == 0)
        {
            outarray[i]          = 0.0;
            outarray[i + nq]     = 0.0;
            outarray[i + 2 * nq] = 0.0;
        }
    }

    return outarray;
}

// Compute the Lamb vector from the velocity vector
// \boldsymbol{\ell} = ( \nabla \times \mathbf{u} ) \times u
Array<OneD, NekDouble> MMFCardiacEP::ComputeLambDiv(
    const Array<OneD, const int> &ValidTimeMap,
    const Array<OneD, const NekDouble> &inarray, const int PlotIndex)
{

        boost::ignore_unused(PlotIndex);

    int nq = m_fields[0]->GetTotPoints();

    // Compute acceleration along the velocity vector:
    // AccMap[0] = (1/2) \nabla || \vec{v} ||^2 \cdot \vec{v}

    // Compute \nabla \times \vec{u} = V_c \vec{k}
    Array<OneD, NekDouble> CovCurl = ComputeCovCurl(inarray, m_movingframes);

    // Filtering
    for (int i = 0; i < nq; ++i)
    {
        if (ValidTimeMap[i] == 0)
        {
            CovCurl[i] = 0.0;
        }
    }

    std::cout << "CovCurl = " << RootMeanSquare(CovCurl);

    //  \boldsymbol{\ell} = V_c \vec{k} \times \mathbf{u}
    Array<OneD, NekDouble> LambVector(m_spacedim * nq);
    LambVector = VectorCrossProdMF(m_movingframes[2], inarray);

    for (int k = 0; k < m_spacedim; ++k)
    {
        Vmath::Vmul(nq, &CovCurl[0], 1, &LambVector[k * nq], 1,
                    &LambVector[k * nq], 1);
    }

    Array<OneD, NekDouble> LambDiv = ComputeCovDiv(LambVector, m_movingframes);

    // Filtering
    for (int i = 0; i < nq; ++i)
    {
        if (ValidTimeMap[i] == 0)
        {
            LambDiv[i] = 0.0;
        }
    }

    // HelmSolveSmoothing(m_LambDivSmoothL, LambDiv);

    // Find a large LambDiv
    int index, Lambflag;
    NekDouble LambDivTol = 10.0;
    for (int i = 0; i < m_fields[0]->GetExpSize(); ++i)
    {
        Lambflag = 0;
        for (int j = 0; j < m_fields[0]->GetTotPoints(i); ++j)
        {
            index = m_fields[0]->GetPhys_Offset(i) + j;
            if (fabs(LambDiv[index]) > LambDivTol)
            {
                Lambflag = 1;
            }
        }

        for (int j = 0; j < m_fields[0]->GetTotPoints(i); ++j)
        {
            index = m_fields[0]->GetPhys_Offset(i) + j;
            if (Lambflag == 1)
            {
                LambDiv[index] = LambDivTol;
            }
        }
    }

    std::cout << ", Lamb Div vector err: Avg = " << RootMeanSquare(LambDiv)
              << ", max = " << Vmath::Vamax(nq, LambDiv, 1) << std::endl;

    // if (PlotIndex > 0)
    // {
    //     // PlotLambDiv(CovCurl, LambDiv, PlotIndex);
    //     int nvar    = 2;
    //     int ncoeffs = m_fields[0]->GetNcoeffs();

    //     std::string outname1 = m_sessionName + "_LambDiv_" +
    //                            boost::lexical_cast<std::string>(PlotIndex) +
    //                            ".chk";
    //     std::vector<Array<OneD, NekDouble>> fieldcoeffs(nvar);
    //     for (int i = 0; i < nvar; ++i)
    //     {
    //         fieldcoeffs[i] = Array<OneD, NekDouble>(ncoeffs);
    //     }

    //     std::vector<std::string> variables(nvar);
    //     variables[0] = "CovCurl";
    //     variables[1] = "LambDiv";

    //     // Normalized Time Vector
    //     m_fields[0]->FwdTrans(CovCurl, fieldcoeffs[0]);
    //     m_fields[0]->FwdTrans(LambDiv, fieldcoeffs[1]);

    //     WriteFld(outname1, m_fields[0], fieldcoeffs, variables);
    // }

    return LambDiv;
}

Array<OneD, NekDouble> MMFCardiacEP::HelmsolvePotentialE(
    const Array<OneD, const int> &ValidTimeMap,
    const Array<OneD, const NekDouble> &inarray, const int PlotIndex)
{
        boost::ignore_unused(PlotIndex);

    int nq      = m_fields[0]->GetTotPoints();
    int ncoeffs = m_fields[0]->GetNcoeffs();

    Array<OneD, NekDouble> physarray(nq);
    Vmath::Vcopy(nq, inarray, 1, physarray, 1);

    // HelmSolveSmoothing(m_LambDivSmoothL, physarray);

    NekDouble DivAvg = -1.0 * AvgInt(physarray);
    Vmath::Sadd(nq, DivAvg, physarray, 1, physarray, 1);

    Array<OneD, NekDouble> outarray(nq, 0.0);

    StdRegions::ConstFactorMap factors;
    factors[StdRegions::eFactorTau]    = m_Helmtau;
    factors[StdRegions::eFactorLambda] = 0.0;

    Array<OneD, NekDouble> tmpc(ncoeffs);
    m_fields[0]->HelmSolve(physarray, tmpc, factors, m_varcoeff);
    m_fields[0]->BwdTrans(tmpc, outarray);
 
    for (int i = 0; i < nq; ++i)
    {
        if (ValidTimeMap[i] == 0)
        {
            outarray[i] = 0.0;
        }
    }

    // Averaging out
    // DivAvg = -1.0 * AvgInt(outarray);
    // Vmath::Sadd(nq, DivAvg, outarray, 1, outarray, 1);

    // if (PlotIndex > 0)
    // {
    //     // PlotLambDiv(CovCurl, LambDiv, PlotIndex);
    //     int nvar    = 1;
    //     int ncoeffs = m_fields[0]->GetNcoeffs();

    //     std::string outname1 = m_sessionName + "_IonU_" +
    //                            boost::lexical_cast<std::string>(PlotIndex) +
    //                            ".chk";
    //     std::vector<Array<OneD, NekDouble>> fieldcoeffs(nvar);
    //     for (int i = 0; i < nvar; ++i)
    //     {
    //         fieldcoeffs[i] = Array<OneD, NekDouble>(ncoeffs);
    //     }

    //     std::vector<std::string> variables(nvar);
    //     variables[0] = "IonU";

    //     // Normalized Time Vector
    //     m_fields[0]->FwdTrans(outarray, fieldcoeffs[0]);

    //     WriteFld(outname1, m_fields[0], fieldcoeffs, variables);
    // }

    return outarray;
}

void MMFCardiacEP::HelmSolveSmoothing(const NekDouble TimeMapSmoothL,
                                   Array<OneD, NekDouble> &outarray)
{
    int ncoeffs = m_fields[0]->GetNcoeffs();
    int nq      = m_fields[0]->GetTotPoints();

    Array<OneD, NekDouble> tmpc(ncoeffs);

    StdRegions::ConstFactorMap factors;
    factors[StdRegions::eFactorTau] = m_Helmtau;

    // 	factors[StdRegions::eFactorLambda] = 1.0 / lambda * m_chi *
    // m_capMembrane;
    NekDouble pi2L                     = 2.0 * m_pi / TimeMapSmoothL;
    factors[StdRegions::eFactorLambda] = pi2L * pi2L;

    Vmath::Smul(nq, -pi2L * pi2L, outarray, 1, outarray, 1);
    // NekDouble DivAvg = -1.0 * AvgInt(outarray);
    // Vmath::Sadd(nq, DivAvg, outarray, 1, outarray, 1);

    m_fields[0]->HelmSolve(outarray, tmpc, factors, m_varcoeff);
    m_fields[0]->BwdTrans(tmpc, outarray);

    // m_contField->HelmSolve(outarray, tmpc, factors, m_varcoeff);
    // m_contField->BwdTrans(tmpc, outarray);
}

void MMFCardiacEP::PlotTimeEnergyMap(
                              const Array<OneD, const NekDouble> &TimeMap,
                              const Array<OneD, const NekDouble> &VelVector,
                              const Array<OneD, const NekDouble> &LambDiv,
                              const Array<OneD, const NekDouble> &IonE)
{
    int nvar    = 8;
    int ncoeffs = m_fields[0]->GetNcoeffs();
    int nq      = m_fields[0]->GetTotPoints();

    std::vector<Array<OneD, NekDouble>> fieldcoeffs(nvar);
    for (int i = 0; i < nvar; ++i)
    {
        fieldcoeffs[i] = Array<OneD, NekDouble>(ncoeffs);
    }

    std::string outname1 = m_TMsessionName + "_TimeEnergyMap.chk";

    std::vector<std::string> variables(nvar);
    variables[0] = "TimeMap";
    variables[1] = "Velocity_x";
    variables[2] = "Velocity_y";
    variables[3] = "Velocity_z";
    variables[4] = "Kinetic E";
    variables[5] = "Lamb Div";
    variables[6] = "Ion U";
    variables[7] = "Total U";

    // U_kin = 0.5 \| \vec{v} \|^2
    Array<OneD, NekDouble> KineticE(nq, 0.0);
    KineticE = ComputeVelocityMag(VelVector);
    Vmath::Vmul(nq, KineticE, 1, KineticE, 1, KineticE, 1);
    Vmath::Smul(nq, 0.5, KineticE, 1, KineticE, 1);

    m_fields[0]->FwdTrans(TimeMap, fieldcoeffs[0]);

    Array<OneD, NekDouble> tmp(nq);
    for (int k = 0; k < m_spacedim; ++k)
    {
        // Compute the magnitude of the velocity map
        Vmath::Vcopy(nq, &VelVector[k * nq], 1, &tmp[0], 1);
        m_fields[0]->FwdTrans(tmp, fieldcoeffs[1 + k]);
    }

    m_fields[0]->FwdTrans(KineticE, fieldcoeffs[4]);
    m_fields[0]->FwdTrans(LambDiv, fieldcoeffs[5]);
    m_fields[0]->FwdTrans(IonE, fieldcoeffs[6]);

    // Totla U = U_Lamb + U_kin
    Array<OneD, NekDouble> TotalE(nq);
    Vmath::Vadd(nq, KineticE, 1, IonE, 1, TotalE, 1);
    m_fields[0]->FwdTrans(TotalE, fieldcoeffs[7]);

    WriteFld(outname1, m_fields[0], fieldcoeffs, variables);
}


void MMFCardiacEP::ComputeTimeMap(const NekDouble time,
                               const NekDouble urest,
                               const Array<OneD, const NekDouble> &field,
                               const Array<OneD, const NekDouble> &dudt,
                               const Array<OneD, const int> &ValidTimeMap,
                               Array<OneD, NekDouble> &dudtHistory,
                               Array<OneD, NekDouble> &dudtMax,
                               Array<OneD, NekDouble> &IappMap,
                               Array<OneD, NekDouble> &TimeMap,
                               const int TimeMapScheme)
{
    boost::ignore_unused(dudtMax, TimeMapScheme);

    int nq = GetTotPoints();

    NekDouble fnewsum;
    NekDouble uTol = 0.01;
    NekDouble dudtTol = 0.01;

    // Compute WeakDGLaplacian
    Array<OneD, NekDouble> Lapu = ComputeCovariantDiffusion(m_movingframes, field);

    NekDouble udiff;
    for (int i = 0; i < nq; ++i)
    {
        udiff = field[i] - urest;
        // Only integrate of time if u > Tol, gradu > Tol, du/dt > 0
        if ((udiff > uTol) && (dudt[i] > dudtTol))
        {
            // Gradient as the main weight
            fnewsum = dudt[i] + dudtHistory[i];

            if(fabs(fnewsum)>dudtTol)
            {
                TimeMap[i] = (dudt[i] * time + dudtHistory[i] * TimeMap[i]) / fnewsum;
                IappMap[i] = (dudt[i] * Lapu[i] + dudtHistory[i] * IappMap[i]) / fnewsum;
            }

            dudtHistory[i] += dudt[i];
        }
    }

    NekDouble TimeMapMin = Vmath::Vmin(nq, TimeMap, 1);
    for (int i = 0; i < nq; ++i)
    {
        if (ValidTimeMap[i] == 0)
        {
            TimeMap[i] = TimeMapMin;
        }
    }

    // std::cout << "Time Map updated = " << cnt << " / " << nq 
    // << ", TimeMap = [ " << Vmath::Vmin(nq, TimeMap, 1) << " , " << Vmath::Vmax(nq, TimeMap, 1) << " ] " << std::endl;

    // TimeMapforInitZone(ValidTimeMap, dudtHistory, TimeMap);
}

void MMFCardiacEP::PlotTimeMap(
    const Array<OneD, const int> &ValidTimeMap,
    const Array<OneD, const Array<OneD, NekDouble>> &AniStrength,
    const Array<OneD, const NekDouble> &TimeMap,
    const int nstep)
{
    boost::ignore_unused(ValidTimeMap);

    int nvar    = 6;
    int nq      = m_fields[0]->GetTotPoints();
    int ncoeffs = m_fields[0]->GetNcoeffs();

    std::string outname1 = m_sessionName + "_TimeMap_" +
                           boost::lexical_cast<std::string>(nstep) + ".chk";

    std::vector<Array<OneD, NekDouble>> fieldcoeffs(nvar);
    for (int i = 0; i < nvar; ++i)
    {
        fieldcoeffs[i] = Array<OneD, NekDouble>(ncoeffs);
    }

    std::vector<std::string> variables(nvar);
    variables[0] = "TimeMap";
    variables[1] = "AniStrength";
    variables[2] = "velmag";
    variables[3] = "velx";
    variables[4] = "vely";
    variables[5] = "velz";

    // index:0 -> u
    m_fields[0]->FwdTrans(TimeMap, fieldcoeffs[0]);

    Array<OneD, Array<OneD, NekDouble>> TimeMapMF(m_spacedim);
    for (int k=0; k<m_spacedim; ++k)
    {
        TimeMapMF[k] = Array<OneD, NekDouble>(nq, 0.0);
    }

    // Compute the gradient of the time map
    m_fields[0]->FwdTrans(AniStrength[0], fieldcoeffs[1]);

    Array<OneD, NekDouble> TMVelocity(m_spacedim * nq);
    TMVelocity = ComputeVelocityTimeMap(m_ValidTimeMap, TimeMap);

    Array<OneD, NekDouble> TMVelMag = ComputeVelocityMag(TMVelocity);
    m_fields[0]->FwdTrans(TMVelMag, fieldcoeffs[2]);

    Array<OneD, NekDouble> tmp(nq);
    for (int k=0; k<m_spacedim; ++k)
    {
        Vmath::Vcopy(nq, &TMVelocity[k*nq], 1, &tmp[0], 1);
        m_fields[0]->FwdTrans(tmp, fieldcoeffs[k+3]);
    }

    WriteFld(outname1, m_fields[0], fieldcoeffs, variables);

    std::cout << "Time Map: Max = " << Vmath::Vmax(nq, TimeMap, 1)
                << ", Min = " << Vmath::Vmin(nq, TimeMap, 1)
                << ", vel mag max = " << Vmath::Vmax(nq, TMVelMag, 1) << std::endl;
}

    // Velocity = ComputeCovGrad(TimeMap, m_movingframes);

    // m_fields[0]->FwdTrans(AniStrength[0], fieldcoeffs[1]);

    // Array<OneD, NekDouble> VelocityMag(nq);
    // VelocityMag = ComputeVelocityMag(Velocity);

    // m_fields[0]->FwdTrans(VelocityMag, fieldcoeffs[2]);

    // for (int i=0; i<nq; ++i)
    // {
    //     if(ValidTimeMap[i] != 1)
    //     {
    //         Velocity[i] = 0.0;
    //         Velocity[i+nq] = 0.0;
    //         Velocity[i+2*nq] = 0.0;
    //     }
    // }

    // Array<OneD, NekDouble> tmp(nq);
    // for (int k=0; k<m_spacedim; ++k)
    // {
    //     Vmath::Vcopy(nq, &Velocity[k*nq], 1, &tmp[0], 1);
    //     m_fields[0]->FwdTrans(tmp, fieldcoeffs[k+3]);
    // }

void MMFCardiacEP::PlotDeformedTimeMap(
    const Array<OneD, const NekDouble> &TimeMap_old,
    const Array<OneD, const NekDouble> &TimeMap_new,
    const Array<OneD, const NekDouble> &TimeMapDiff,
    const Array<OneD, const Array<OneD, NekDouble>> &Vdiff,
    const Array<OneD, const NekDouble> &VdiffMag,
    const Array<OneD, const NekDouble> &VdiffDivergence)
{
    int nvar    = 10;
    int nq      = m_fields[0]->GetTotPoints();
    int ncoeffs = m_fields[0]->GetNcoeffs();

    std::string outname1 = m_sessionName + "_TimeMapDeform.chk";

    std::vector<Array<OneD, NekDouble>> fieldcoeffs(nvar);
    for (int i = 0; i < nvar; ++i)
    {
        fieldcoeffs[i] = Array<OneD, NekDouble>(ncoeffs);
    }

    Array<OneD, NekDouble> TimeMapDiffExact(nq);
    Vmath::Vsub(nq, TimeMap_new, 1, TimeMap_old, 1, TimeMapDiffExact, 1);
    
    Array<OneD, NekDouble> TimeMapError(nq);
    Vmath::Vsub(nq, TimeMapDiffExact, 1, TimeMapDiff, 1, TimeMapError, 1);

    NekDouble TMerrMax = Vmath::Vamax(nq, TimeMapError, 1) / Vmath::Vamax(nq, TimeMap_old, 1);

    std::cout << "Max. TimeMapError = " << TMerrMax << std::endl;

    std::vector<std::string> variables(nvar);
    variables[0] = "TimeMapDiffExact";
    variables[1] = "TimeMapDiff";
    variables[2] = "Vdiffx";
    variables[3] = "Vdiffy";
    variables[4] = "Vdiffz";
    variables[5] = "VdiffMag";
    variables[6] = "VdiffDivergence";
    variables[7] = "TimeMapError";
    variables[8] = "TimeMap_old";
    variables[9] = "TimeMap_new";

    // index:0 -> u
    m_fields[0]->FwdTrans(TimeMapDiffExact, fieldcoeffs[0]);
    m_fields[0]->FwdTrans(TimeMapDiff, fieldcoeffs[1]);

    m_fields[0]->FwdTrans(Vdiff[0], fieldcoeffs[2]);
    m_fields[0]->FwdTrans(Vdiff[1], fieldcoeffs[3]);
    m_fields[0]->FwdTrans(Vdiff[2], fieldcoeffs[4]);
    m_fields[0]->FwdTrans(VdiffMag, fieldcoeffs[5]);
    m_fields[0]->FwdTrans(VdiffDivergence, fieldcoeffs[6]);

    m_fields[0]->FwdTrans(TimeMapError, fieldcoeffs[7]);
    m_fields[0]->FwdTrans(TimeMap_old, fieldcoeffs[8]);
    m_fields[0]->FwdTrans(TimeMap_new, fieldcoeffs[9]);

    WriteFld(outname1, m_fields[0], fieldcoeffs, variables);
}



void MMFCardiacEP::PlotTimeMapMF(
    const Array<OneD, const NekDouble> &NoboundaryZone,
    const Array<OneD, const NekDouble> &TimeMap,
    const Array<OneD, const Array<OneD, NekDouble>> &TimeMapMF,
    const Array<OneD, const Array<OneD, NekDouble>> &MFFirst,
    const Array<OneD, const Array<OneD, Array<OneD, NekDouble>>>
        &MF1stConnection,
    const Array<OneD, const Array<OneD, NekDouble>> &Relacc, const int nstep)
{
    int nvar    = 9;
    int nq      = m_fields[0]->GetTotPoints();
    int ncoeffs = m_fields[0]->GetNcoeffs();

    std::string outname1 = m_sessionName + "_MFTM_" +
                           boost::lexical_cast<std::string>(nstep) + ".chk";

    std::vector<Array<OneD, NekDouble>> fieldcoeffs(nvar);
    for (int i = 0; i < nvar; ++i)
    {
        fieldcoeffs[i] = Array<OneD, NekDouble>(ncoeffs);
    }

    std::vector<std::string> variables(nvar);
    variables[0] = "TimeMap";
    variables[1] = "ex1";
    variables[2] = "ey1";
    variables[3] = "ez1";
    variables[4] = "AngleDiff";
    variables[5] = "w211";
    variables[6] = "w212";
    variables[7] = "RelAcc";
    variables[8] = "CondBlock";

    // index:0 -> u
    m_fields[0]->FwdTrans(TimeMap, fieldcoeffs[0]);

    // index:[1, 2, 3] -> ex1, ey1, ez1
    Array<OneD, NekDouble> tmp(nq);
    for (int j = 0; j < m_spacedim; ++j)
    {
        Vmath::Vcopy(nq, &TimeMapMF[0][j * nq], 1, &tmp[0], 1);
        m_fields[0]->FwdTrans(tmp, fieldcoeffs[j + 1]);
    }

    // Compute AngleMF: Compute the angle bewteen the original MF and
    // aligned MF
    Array<OneD, NekDouble> AngleDiff(nq, 0.0);
    Array<OneD, NekDouble> MFerr(nq, 0.0);

    NekDouble e1x, e1y, e1z, e1xnew, e1ynew, e1znew;
    NekDouble diffx, diffy, diffz;
    // NekDouble erx, ery, erz, differx, differy, differz;
    for (int i = 0; i < nq; i++)
    {
        // erx = m_polarMF[0][i];
        // ery = m_polarMF[0][i + nq];
        // erz = m_polarMF[0][i + 2 * nq];

        e1x = MFFirst[0][i];
        e1y = MFFirst[0][i + nq];
        e1z = MFFirst[0][i + 2 * nq];

        e1xnew = TimeMapMF[0][i];
        e1ynew = TimeMapMF[0][i + nq];
        e1znew = TimeMapMF[0][i + 2 * nq];

        diffx = e1x - e1xnew;
        diffy = e1y - e1ynew;
        diffz = e1z - e1znew;

        if (NoboundaryZone[i] == 1)
        {
            MFerr[i] = sqrt(diffx * diffx + diffy * diffy + diffz * diffz);
        }
    }

    std::cout << "MFFirst vs. TimeMapMF: MFerr = " << RootMeanSquare(MFerr)
              << std::endl;

    // Angle between MF and fibre
    m_fields[0]->FwdTrans(MFerr, fieldcoeffs[4]);

    Array<OneD, NekDouble> w211(nq);
    Array<OneD, NekDouble> w212(nq);

    Vmath::Vcopy(nq, &MF1stConnection[0][0][0], 1, &w211[0], 1);
    Vmath::Vcopy(nq, &MF1stConnection[0][1][0], 1, &w212[0], 1);

    // Ignore the region where w211 is too big.
    Vmath::Vmul(nq, NoboundaryZone, 1, w211, 1, w211, 1);
    Vmath::Vmul(nq, NoboundaryZone, 1, w212, 1, w212, 1);

    std::cout << "TimeMap: w211 max = " << Vmath::Vamax(nq, w211, 1)
              << std::endl;

    // Connection form w211
    m_fields[0]->FwdTrans(w211, fieldcoeffs[5]);

    // Connection form w212
    m_fields[0]->FwdTrans(w212, fieldcoeffs[6]);

    // Relative Acceleration (I)
    Array<OneD, NekDouble> RelAccetmp(nq);

    Vmath::Vcopy(nq, &Relacc[1][0], 1, &RelAccetmp[0], 1);
    Vmath::Vmul(nq, &NoboundaryZone[0], 1, &RelAccetmp[0], 1, &RelAccetmp[0],
                1);

    m_fields[0]->FwdTrans(RelAccetmp, fieldcoeffs[7]);

    // Compute Conduction block
    Array<OneD, NekDouble> CBlock(nq, 0.0);
    for (int i = 0; i < nq; ++i)
    {
        if ((w212[i] > 0) && (RelAccetmp[i] < 0))
        {
            CBlock[i] = w212[i] - 10.0 * RelAccetmp[i];
        }
    }

    m_fields[0]->FwdTrans(CBlock, fieldcoeffs[8]);

    WriteFld(outname1, m_fields[0], fieldcoeffs, variables);
}

// Compute Velocity field from Time Map
// flag = 1: Use \vec{v} = \nabla T / \| \nabla T \|^2
// flag = 0:
void MMFCardiacEP::ComputeMFTimeMap(const Array<OneD, const int> &ValidTimeMap,
                                 const Array<OneD, const NekDouble> &inarray,
                                 Array<OneD, int> &NewValidTimeMap,
                                 Array<OneD, Array<OneD, NekDouble>> &TMMF)
{
    int nq = m_fields[0]->GetTotPoints();

    Array<OneD, NekDouble> outarray(m_spacedim * nq, 0.0);

    Array<OneD, NekDouble> physarray(nq);
    Vmath::Vcopy(nq, inarray, 1, physarray, 1);

    for (int i = 0; i < nq; ++i)
    {
        NewValidTimeMap[i] = ValidTimeMap[i];
    }

    // TmapGrad = \nabla Tmap
    Array<OneD, NekDouble> TmapGrad(m_spacedim * nq);
    TmapGrad = ComputeCovGrad(physarray, m_movingframes);

    Array<OneD, NekDouble> TmapGradMag(nq);
    TmapGradMag = ComputeVelocityMag(TmapGrad);

    // Compute VelField \vec{v} = \sum_{i=1}^3 1/(\nabla T \cdot \hat{x}_i)
    // \hat{x}_i
    NekDouble tmp, TmapGradTol = 0.01;
    for (int i = 0; i < nq; i++)
    {
        tmp = TmapGradMag[i];
        for (int k = 0; k < m_spacedim; ++k)
        {
            if ((ValidTimeMap[i] == 1) && (tmp > TmapGradTol))
            {
                TMMF[0][i + k * nq] = TmapGrad[i + k * nq] / tmp;
                NewValidTimeMap[i]  = 1;
            }

            else
            {
                TMMF[0][i + k * nq] = m_movingframes[0][i + k * nq];
                NewValidTimeMap[i]  = 0;
            }
        }
    }

    NekDouble MF1x, MF1y, MF1z, MF2x, MF2y, MF2z, MF3x, MF3y, MF3z;
    NekDouble MFmag1, MFmag2;
    for (int i = 0; i < nq; i++)
    {
        MF1x = TMMF[0][i];
        MF1y = TMMF[0][i + nq];
        MF1z = TMMF[0][i + 2 * nq];

        MFmag1 = sqrt(MF1x * MF1x + MF1y * MF1y + MF1z * MF1z);

        TMMF[0][i]          = MF1x / MFmag1;
        TMMF[0][i + nq]     = MF1y / MFmag1;
        TMMF[0][i + 2 * nq] = MF1z / MFmag1;

        MF1x = TMMF[0][i];
        MF1y = TMMF[0][i + nq];
        MF1z = TMMF[0][i + 2 * nq];

        TMMF[2][i]          = m_movingframes[2][i];
        TMMF[2][i + nq]     = m_movingframes[2][i + nq];
        TMMF[2][i + 2 * nq] = m_movingframes[2][i + 2 * nq];

        MF3x = TMMF[2][i];
        MF3y = TMMF[2][i + nq];
        MF3z = TMMF[2][i + 2 * nq];

        MF2x = MF3y * MF1z - MF3z * MF1y;
        MF2y = MF3z * MF1x - MF3x * MF1z;
        MF2z = MF3x * MF1y - MF3y * MF1x;

        MFmag2              = sqrt(MF2x * MF2x + MF2y * MF2y + MF2z * MF2z);
        TMMF[1][i]          = MF2x / MFmag2;
        TMMF[1][i + nq]     = MF2y / MFmag2;
        TMMF[1][i + 2 * nq] = MF2z / MFmag2;
    }

    // Check TimeMap Moving frames
    std::cout << "============= Check Moving frames from TimeMap =============" << std::endl;
}

void MMFCardiacEP::v_EvaluateExactSolution(unsigned int field,
                                           Array<OneD, NekDouble> &outfield,
                                           const NekDouble time)
{
    EquationSystem::v_EvaluateExactSolution(field, outfield, time);
}

void MMFCardiacEP::AlievPanfilovReaction(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time)
{
    boost::ignore_unused(time);
    int nq      = m_fields[0]->GetTotPoints();

    // inarray[0] holds initial physical u values throughout
    // inarray[1] holds initial physical v values throughout

    Array<OneD, NekDouble> m_tmp1(nq);
    Array<OneD, NekDouble> m_tmp2(nq);
    Array<OneD, NekDouble> m_uu(nq);
    Array<OneD, NekDouble> m_uuu(nq);

    // compute u^2: m_u = u*u
    Vmath::Vmul(nq, &inarray[0][0], 1, &inarray[0][0], 1, &m_uu[0], 1);

    // compute u^3: m_u = u*u*u
    Vmath::Vmul(nq, &inarray[0][0], 1, &m_uu[0], 1, &m_uuu[0], 1);
 
    // Ru = au
    Vmath::Smul(nq, m_a, &inarray[0][0], 1, &m_tmp1[0], 1);
    // Ru = (-1-a)u*u + au
    Vmath::Svtvp(nq, (-1.0 - m_a), &m_uu[0], 1, &m_tmp1[0], 1, &m_tmp1[0], 1);
    //        }
    // Ru = u*u*u - (1+a)u*u + au
    Vmath::Vadd(nq, &m_uuu[0], 1, &m_tmp1[0], 1, &m_tmp1[0], 1);
 
    Vmath::Smul(nq, m_k, &m_tmp1[0], 1, &m_tmp1[0], 1);
 
    // Ru = k(u*u*u - (1+a)u*u + au) + I_stim
    Vmath::Vadd(nq, &outarray[0][0], 1, &m_tmp1[0], 1, &outarray[0][0], 1);

    // Ru = k(u*u*u - (1+a)u*u + au) + uv + I_stim
    Vmath::Vvtvp(nq, &inarray[0][0], 1, &inarray[1][0], 1, &m_tmp1[0], 1,
                 &outarray[0][0], 1);
    // Ru = -k(u*u*u - (1+a)u*u + au) - uv - I_stim
    Vmath::Neg(nq, &outarray[0][0], 1);

    // --------------------------------------
    // Compute reaction term g(u,v)
    // --------------------------------------
    // tmp2 = mu2 + u
    Vmath::Sadd(nq, m_mu2, &inarray[0][0], 1, &m_tmp2[0], 1);

    // tmp2 = v/(mu2 + u)
    Vmath::Vdiv(nq, &inarray[1][0], 1, &m_tmp2[0], 1, &m_tmp2[0], 1);

    // tmp2 = mu1*v/(mu2 + u)
    Vmath::Smul(nq, m_mu1, &m_tmp2[0], 1, &m_tmp2[0], 1);

    // tmp1 = Eps + mu1*v/(mu2+u)
    Vmath::Sadd(nq, m_eps, &m_tmp2[0], 1, &m_tmp2[0], 1);

    Vmath::Sadd(nq, (-m_a - 1), &inarray[0][0], 1, &m_tmp1[0], 1);
 
    Vmath::Smul(nq, m_k, &m_tmp1[0], 1, &m_tmp1[0], 1);
 
    // tmp1 = ku(u-a-1) + v
    Vmath::Vvtvp(nq, &inarray[0][0], 1, &m_tmp1[0], 1, &inarray[1][0], 1,
                 &m_tmp1[0], 1);

    // tmp1 = -ku(u-a-1)-v
    Vmath::Neg(nq, &m_tmp1[0], 1);

    // outarray = [Eps + mu1*v/(mu2+u)] * [-ku(u-a-1)-v]
    Vmath::Vmul(nq, &m_tmp1[0], 1, &m_tmp2[0], 1, &outarray[1][0], 1);
}

void MMFCardiacEP::v_GenerateSummary(SolverUtils::SummaryList &s)
{
    MMFSystem::v_GenerateSummary(s);
    AddSummaryItem(s, "SolverSchemeType", SolverSchemeTypeMap[m_SolverSchemeType]);

    SolverUtils::AddSummaryItem(s, "TimeMapScheme", m_TimeMapScheme);
    SolverUtils::AddSummaryItem(s, "TimeMapStart", m_TimeMapStart);
    SolverUtils::AddSummaryItem(s, "TimeMapEnd", m_TimeMapEnd);

    // if(m_SolverSchemeType==eTimeMapMarch)
    // {
    SolverUtils::AddSummaryItem(s, "TimeMapIapp", m_TimeMapIapp);
    SolverUtils::AddSummaryItem(s, "m_TimeMapDelay", m_TimeMapDelay);
    // }

    SolverUtils::AddSummaryItem(s, "AniRegionStart", m_AniRegionStart);
    SolverUtils::AddSummaryItem(s, "AniRegionEnd", m_AniRegionEnd);
    SolverUtils::AddSummaryItem(s, "AnisotropyStrength", m_AnisotropyStrength);
    SolverUtils::AddSummaryItem(s, "urest", m_urest);

    m_cell->GenerateSummary(s);
    
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
