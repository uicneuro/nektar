////////////////////////////////////////////////////////////////////////////
//
// File MMFStaticDiffOp.cpp
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
// Permission is hereby granted, free of charge, to any peron obtaining a
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
// DEALINstd IN THE SOFTWARE.
//
// Description: MMF solve routines
//
///////////////////////////////////////////////////////////////////////////////

// #include <iostream>
// #include <iomanip>
// #include <boost/algorithm/string.hpp>

// #include <MMFSolver/EquationSystems/MMFStaticDiffOp.h>
// #include <MultiRegions/AssemblyMap/AssemblyMapDG.h>

// #include <SolverUtils/Driver.h>
// #include <LibUtilities/BasicUtils/SessionReader.h>

// #include <ADRSolver/EquationSystems/UnsteadyDiffusion.h>

// #include <boost/math/special_functions/spherical_harmonic.hpp>
// using namespace std;
// using namespace Nektar::SolverUtils;
// using namespace Nektar;

#include <MMFSolver/EquationSystems/MMFStaticDiffOp.h>
#include <boost/algorithm/string.hpp>
#include <iomanip>
#include <iostream>

#include <LibUtilities/BasicUtils/Timer.h>

#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>
#include <SolverUtils/MMFSystem.h>

#include <typeinfo>

#include <LibUtilities/BasicUtils/Timer.h>
#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>

namespace Nektar
{

std::string MMFStaticDiffOp::className =
    SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
        "MMFStaticDiffOp", MMFStaticDiffOp::create,
        "MMFStaticDiffOp equation.");

MMFStaticDiffOp::MMFStaticDiffOp(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const SpatialDomains::MeshGraphSharedPtr &pGraph)
    : UnsteadySystem(pSession, pGraph), MMFSystem(pSession, pGraph)
{
}

/**
 * @brief Initialisation object for the unsteady linear advection equation.
 */
void MMFStaticDiffOp::v_InitObject(bool DeclareFields)
{
    // Call to the initialisation object
    UnsteadySystem::v_InitObject(DeclareFields);

    int nq       = m_fields[0]->GetNpoints();
    int shapedim = m_fields[0]->GetShapeDimension();
    Array<OneD, Array<OneD, NekDouble>> Anisotropy(shapedim);
    for (int j = 0; j < shapedim; ++j)
    {
        Anisotropy[j] = Array<OneD, NekDouble>(nq, 1.0);
    }

    MMFSystem::MMFInitObject(Anisotropy);

    std::string Dtype ;
    m_session->LoadSolverInfo("DerivType", Dtype, "Covariant");
    if (Dtype == "Covariant")
        m_DerivType = eCovariant;
    if (Dtype == "Euclidean")
        m_DerivType = eEuclidean;  

 
    // Define TestType
    ASSERTL0(m_session->DefinesSolverInfo("TESTTYPE"),
             "No TESTTYPE defined in session.");
    std::string TestTypeStr = m_session->GetSolverInfo("TESTTYPE");
    for (int i = 0; i < (int)SIZE_TestType; ++i)
    {
        if (TestTypeMap[i] == TestTypeStr)
        {
            m_TestType = (TestType)i;
            break;
        }
    }

    // Read the advection velocities from session file
    m_session->LoadParameter("advx", m_advx, 1.0);
    m_session->LoadParameter("advy", m_advy, 1.0);
    m_session->LoadParameter("advz", m_advz, 1.0);

    // Compute Phi and Theta Direction
    ComputeSphericalTangentVector(m_Phi, m_Theta);

    std::vector<std::string> vel;
    vel.push_back("Vx");
    vel.push_back("Vy");
    vel.push_back("Vz");

    // Resize the advection velocities vector to dimension of the problem
    vel.resize(m_spacedim);

    // Store in the global variable m_velocity the advection velocities
    m_velocity = Array<OneD, Array<OneD, NekDouble>>(m_spacedim);
    for (int k = 0; k < m_spacedim; ++k)
    {
        m_velocity[k] = Array<OneD, NekDouble>(nq,0.0);
    }

    m_velphi = Array<OneD, NekDouble>(nq,0.0);
    m_veltheta = Array<OneD, NekDouble>(nq,0.0);

    NekDouble gms = 9.80616;
    switch (m_TestType)
    {
        case eTestCovariantSphere:
        {
            m_session->LoadParameter("alpha", m_alpha, 0.0);
            m_session->LoadParameter("u0", m_u0, 2.0 * m_pi / 12.0);

            EvaluateCovariantVelocity(m_alpha, m_velphi, m_veltheta, m_velocity);

            Array<OneD, NekDouble> velmag(nq);
            velmag = ComputeVelocityMag(m_velocity);
            PlotFieldVector(velmag, m_velocity, 0);
        }
        break;

        case eTestRiemannCurv:
        {

        }
        break;

        case eTestDivPlane:
        {
            Array<OneD, NekDouble> x0(nq);
            Array<OneD, NekDouble> x1(nq);
            Array<OneD, NekDouble> x2(nq);

            m_fields[0]->GetCoords(x0, x1, x2);

            for (int i = 0; i < nq; i++)
            {
                m_velocity[0][i] = cos(x0[i]);
                m_velocity[1][i] = sin(x1[i]);
                m_velocity[2][i] = 0.0;
            }
        }
        break;

        case eTestDivSphericMetId:
        {
            EvaluateAdvectionVelocity(m_velocity);
        }
        break;

        case eTestDivRossbyMetId:
        case eTestDivRossbyFlow:
        {
            NekDouble SecondToDay = 60.0 * 60.0 * 24.0;
            NekDouble rad_earth   = 6.37122 * 1000000;
            NekDouble Omegams;

            // Nondimensionalized coeffs.
            m_g = (gms * SecondToDay * SecondToDay) / rad_earth;

            m_session->LoadParameter("Omega", Omegams, 7.292 * 0.00001);
            m_Omega = Omegams * SecondToDay;

            m_session->LoadParameter("H0", m_H0, 8000.0);
            m_H0 = m_H0 / rad_earth;

            m_W = 7.848 * 0.000001 * SecondToDay;
            m_K = 7.848 * 0.000001 * SecondToDay;

            EvaluateAdvectionVelocity(m_velocity);
        }
        break;

        case eTestGradPlaneFlow:
        case eTestGradSphericMetID:
        case eTestGradRossbyFlow:
        {
            NekDouble Omegams;
            NekDouble SecondToDay = 60.0 * 60.0 * 24.0;

            m_session->LoadParameter("Omega", Omegams, 7.292 * 0.00001);
            m_Omega = Omegams * SecondToDay;

            m_W = 7.848 * 0.000001 * SecondToDay;
            m_K = 7.848 * 0.000001 * SecondToDay;
        }
        break;

        case eTestLaplacianPlaneFlow:
        case eTestLaplacianMetricID:
        case eTestLaplacianSimpleFlow:
        case eTestLaplacianRossbyFlow:
        {
            NekDouble SecondToDay = 60.0 * 60.0 * 24.0;
            NekDouble Omegams;

            m_session->LoadParameter("Omega", Omegams, 7.292 * 0.00001);
            m_Omega = Omegams * SecondToDay;

            m_session->LoadParameter("W", m_W, 7.848 * 0.000001 * SecondToDay);
            m_session->LoadParameter("K", m_K, 7.848 * 0.000001 * SecondToDay);

            m_session->LoadParameter("Diffbeta", m_Diffbeta, 0.5);
            m_session->LoadParameter("Diffeta", m_Diffeta, 100.0);
            m_session->LoadParameter("Diffhe", m_Diffhe, 0.5);
        }
        break;

        case eTestCurlRossbyMetId:
        case eTestCurlRossbyFlow:
        {
            NekDouble SecondToDay = 60.0 * 60.0 * 24.0;
            NekDouble Omegams;

            m_session->LoadParameter("Omega", Omegams, 7.292 * 0.00001);
            m_Omega = Omegams * SecondToDay;

            m_session->LoadParameter("W", m_W, 7.848 * 0.000001 * SecondToDay);
            m_session->LoadParameter("K", m_K, 7.848 * 0.000001 * SecondToDay);
        }
        break;

        case eTestCurlSphericMetId:
        default:
            break;
    }

    // Compute MF \cdot Theta and MF \cdot Phi: To Get m_ThetacdotMF and
    // m_PhicdotMF;
    // ComputeMFcdotSphericalCoord(m_graddir);

    // Compute m_traceVn = n \cdot v
    // GetNormalVelocity(m_traceVn);
    ComputencdotMF(m_movingframes, m_ncdotMFFwd, m_ncdotMFBwd);

    // Compute the cross producted MF
    ComputeMFtimesMF(m_movingframes, m_CrossProductMF);

    if(m_surfaceType==eSphere)
    {
       std::cout << "Compute Vel vector on the sphere " << std::endl;

        // Compute Velocity vector on the sphere
        Array<OneD, NekDouble> Velvec0(m_spacedim*nq);

        NekDouble sin_varphi, cos_varphi, sin_theta, cos_theta;

        Array<OneD, NekDouble> x0(nq);
        Array<OneD, NekDouble> x1(nq);
        Array<OneD, NekDouble> x2(nq);

        m_fields[0]->GetCoords(x0, x1, x2);

        NekDouble uhat, vhat;

        NekDouble tmp1, tmp2, cos4phi, sin4phi;
        for (int j = 0; j < nq; j++)
        {
            CartesianToNewSpherical(x0[j], x1[j], x2[j], sin_varphi, cos_varphi,
                                    sin_theta, cos_theta);

            tmp1    = 2 * sin_varphi * cos_varphi;
            tmp2    = 2 * cos_varphi * cos_varphi - 1.0;
            sin4phi = 2 * tmp1 * tmp2;

            tmp1    = 2 * cos_varphi * cos_varphi - 1.0;
            cos4phi = 2 * tmp1 * tmp1 - 1.0;

            uhat = m_W * sin_theta +
                    m_K * sin_theta * sin_theta * sin_theta *
                        (4 * cos_theta * cos_theta - sin_theta * sin_theta) *
                        cos4phi;
            vhat = -4.0 * m_K * sin_theta * sin_theta * sin_theta *
                    cos_theta * sin4phi;

            Velvec0[j] = -1.0 * uhat * sin_varphi - vhat * cos_theta * cos_varphi;
            Velvec0[j+nq] = uhat * cos_varphi - vhat * cos_theta * sin_varphi;
            Velvec0[j+2*nq] = vhat * sin_theta;
        }

        Array<OneD, NekDouble> Velveceps(m_spacedim*nq,0.0);
        Array<OneD, NekDouble> V0cdote1(nq, 0.0);
        Array<OneD, NekDouble> V0cdote2(nq, 0.0);
        for (int i=0; i<m_spacedim; ++i)
        {
            Vmath::Vvtvp(nq, &Velvec0[i*nq], 1, &m_movingframes[0][i*nq], 1, &V0cdote1[0], 1, &V0cdote1[0], 1);
            Vmath::Vvtvp(nq, &Velvec0[i*nq], 1, &m_movingframes[1][i*nq], 1, &V0cdote2[0], 1, &V0cdote2[0], 1);
        }

        for (int i=0; i<m_spacedim; ++i)
        {
            Vmath::Vmul(nq, &V0cdote1[0], 1, &m_movingframes[0][i*nq], 1, &Velveceps[i*nq], 1);
            Vmath::Vvtvp(nq, &V0cdote2[0], 1, &m_movingframes[1][i*nq], 1, &Velveceps[i*nq], 1, &Velveceps[i*nq], 1);
        }

        Array<OneD, Array<OneD, NekDouble>> SphericalVector;
        ComputeSphericalVector(SphericalVector);

        Array<OneD, NekDouble> test(nq);
        Array<OneD, NekDouble> diff(m_spacedim*nq);

        // Construct The Moving Frames
        Array<OneD, Array<OneD, NekDouble>> mf_LOCSPH;
        mf_LOCSPH = Array<OneD, Array<OneD, NekDouble>>(m_mfdim);

        for (int j = 0; j < m_mfdim; ++j)
        {
            mf_LOCSPH[j] = Array<OneD, NekDouble>(m_spacedim * nq, 0.0);
        }

        // Get LOCAL moving frames such that \nabla \cdot e^i is the minimum
        GetLOCALMovingframes(mf_LOCSPH);
        ComputeAxisAlignedLOCALMovingframes(m_sphereMF, mf_LOCSPH);
    }


    // Test connections of moving frames
    switch(m_surfaceType)
    {
        case ePlane:
        case ePolar:
        {
            TestPolarConnection1form(m_polarMF, m_MMFActivation);
        }
        break;

        case eSphere:
        {
            TestSphericalConnection1form(m_sphereMF, m_MMFActivation);
        }
        break;

        case ePseudosphere:
        {
            TestPseudosphericalConnection1form(m_pseudosphereMF, m_MMFActivation);
        }
        break;

        default:
        break;
    }

    // Compute Error
    switch (m_TestType)
    {
        case eTestCovariantSphere:
        {
            CheckErrorCovariantDeriv(m_velphi, m_veltheta, m_velocity);
        }
        break;

        case eTestRelAcc:
        {
            CheckErrorRelAcc(m_movingframes);
        }
        break;

        case eTestRiemannCurv:
        {
            CheckErrorRiemCrv(m_movingframes);
        }
        break;

        case eTestDivPlane:
        case eTestDivSphericMetId:
        case eTestDivRossbyMetId:
        case eTestDivRossbyFlow:
        {
            CheckErrorDivergence(SetInitialU());
        }
        break;

        case eTestCurlPlaneFlow:
        case eTestCurlSphericMetId:
        case eTestCurlRossbyMetId:
        case eTestCurlRossbyFlow:
        {
            Array<OneD, NekDouble> uth;
            Array<OneD, NekDouble> uphi;

            SetInitialU(uth, uphi);
            CheckErrorCurl(uth, uphi);
        }
        break;

        case eTestGradPlaneFlow:
        case eTestGradSphericMetID:
        case eTestGradRossbyFlow:
        {
            CheckErrorGradient(SetInitialU());
        }
        break;

        case eTestLaplacianPlaneFlow:
        case eTestLaplacianMetricID:
        case eTestLaplacianSimpleFlow:
        case eTestLaplacianRossbyFlow:
        {
            CheckErrorLaplacian(SetInitialU());
        }
        break;



        default:
            break;
    }

    Wait_On_Enter();

    // If explicit it computes RHS and PROJECTION for the time integration
    if (m_explicitAdvection)
    {
        m_ode.DefineOdeRhs(&MMFStaticDiffOp::DoOdeRhs, this);
        m_ode.DefineProjection(&MMFStaticDiffOp::DoOdeProjection, this);
    }
    // Otherwise it gives an error (no implicit integration)
    else
    {
        ASSERTL0(false, "Implicit unsteady Advection not set up.");
    }

    //  if (!m_explicitDiffusion)
    //  {
    //	  m_ode.DefineImplicitSolve (&MMFStaticDiffOp::DoImplicitSolve,
    // this);
    // }
}

/**
 * @brief Unsteady linear advection equation destructor.
 */
MMFStaticDiffOp::~MMFStaticDiffOp()
{
}

void MMFStaticDiffOp::v_DoSolve()
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
    Array<OneD, Array<OneD, NekDouble>> tmp(nvariables);

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

    int Ntot = m_steps / m_checksteps + 1;
    Array<OneD, NekDouble> dMass(Ntot);

    Array<OneD, NekDouble> zeta(nq);
    Array<OneD, Array<OneD, NekDouble>> fieldsprimitive(nvariables);
    for (int i = 0; i < nvariables; ++i)
    {
        fieldsprimitive[i] = Array<OneD, NekDouble>(nq);
    }

    while (step < m_steps || m_time < m_fintime - NekConstants::kNekZeroTol)
    {
        timer.Start();
        fields = m_intScheme->TimeIntegrate(step, m_timestep, m_ode);
        timer.Stop();

        m_time += m_timestep;
        elapsed = timer.TimePerTest(1);
        intTime += elapsed;
        cpuTime += elapsed;

        // Write out status information
        if (m_session->GetComm()->GetRank() == 0 && !((step + 1) % m_infosteps))
        {
            std::cout << "Steps: " << std::setw(5) << std::left << step + 1
                      << " "
                      << "Time: " << std::setw(6) << std::left << m_time;

            std::stringstream ss;
            ss << cpuTime << "s";
            std::cout << " CPU Time: " << std::setw(8) << std::left << ss.str();

            // Masss = h^*
            // indx = (step + 1) / m_checksteps;
            // dMass[indx] = (m_fields[0]->PhysIntegral(fields[0]) -
            // m_Mass0)/m_Mass0; std::cout << ",   dMass: " << std::setw(8)
            // << std::left << dMass[indx] << std::endl << std::endl;

            cpuTime = 0.0;
        }

        // Transform data into coefficient space
        for (i = 0; i < nvariables; ++i)
        {
            m_fields[m_intVariables[i]]->SetPhys(fields[i]);
            m_fields[m_intVariables[i]]->FwdTransLocalElmt(
                fields[i], m_fields[m_intVariables[i]]->UpdateCoeffs());
            m_fields[m_intVariables[i]]->SetPhysState(false);
        }

        // Write out checkpoint files
        if ((m_checksteps && step && !((step + 1) % m_checksteps)) ||
            doCheckTime)
        {
            Checkpoint_Output(nchk++);
            doCheckTime = false;
        }

        // Step advance
        ++step;
    }

    // Print out summary statistics
    if (m_session->GetComm()->GetRank() == 0)
    {
        if (m_cflSafetyFactor > 0.0)
        {
            std::cout << "CFL safety factor : " << m_cflSafetyFactor
                      << std::endl
                      << "CFL time-step     : " << m_timestep << std::endl;
        }

        if (m_session->GetSolverInfo("Driver") != "SteadyState")
        {
            std::cout << "Time-integration  : " << intTime << "s" << std::endl;
        }

        // Output for error distributions
        Checkpoint_Err(m_time, fields);
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
}

void MMFStaticDiffOp::Checkpoint_Err(
    const NekDouble time,
    const Array<OneD, const Array<OneD, NekDouble>> &fieldphys)
{
    int nvar    = m_fields.size();
    int nq      = m_fields[0]->GetTotPoints();
    int ncoeffs = m_fields[0]->GetNcoeffs();

    std::string outname = m_sessionName + "Err.chk";

    std::vector<Array<OneD, NekDouble>> fieldcoeffs(nvar);
    for (int i = 0; i < nvar; ++i)
    {
        fieldcoeffs[i] = Array<OneD, NekDouble>(ncoeffs);
    }

    std::vector<std::string> variables(nvar);
    variables[0] = "Linf";

    Array<OneD, NekDouble> exactsoln(m_fields[0]->GetNpoints());
    // EvaluateFunction(m_session->GetVariable(0),
    // exactsoln,"ExactSolution", time);
    //  void MMFStaticDiffOp::v_EvaluateExactSolution(unsigned int field,
    //					     Array<OneD, NekDouble>
    //&outfield, 					     const NekDouble
    // time)

    v_EvaluateExactSolution(0, exactsoln, time);

    Vmath::Vsub(nq, exactsoln, 1, fieldphys[0], 1, exactsoln, 1);
    for (int i = 0; i < nq; ++i)
    {
        exactsoln[i] = fabs(exactsoln[i]);
    }

    m_fields[0]->FwdTrans(exactsoln, fieldcoeffs[0]);
    WriteFld(outname, m_fields[0], fieldcoeffs, variables);
}

void MMFStaticDiffOp::PlotCovError(
    const int dir,
    const Array<OneD, const NekDouble> &Errtheta,
    const Array<OneD, const NekDouble> &Errphi,
    const Array<OneD, const NekDouble> &ErrVec)
{
    int nvar    = 3;
    int ncoeffs = m_fields[0]->GetNcoeffs();

    std::vector<std::string> variables(nvar);
    variables[0]  = "ErrTheta";
    variables[1]  = "ErrPhi";
    variables[2]  = "ErrVec";

    std::string outname = m_sessionName + "_dir" + boost::lexical_cast<std::string>(dir)+ "_Err.chk";

    std::vector<Array<OneD, NekDouble>> fieldcoeffs(nvar);
    for (int j = 0; j < nvar; ++j)
    {
        fieldcoeffs[j] = Array<OneD, NekDouble>(ncoeffs);
    }

    m_fields[0]->FwdTrans(Errtheta, fieldcoeffs[0]);
    m_fields[0]->FwdTrans(Errphi, fieldcoeffs[1]);
    m_fields[0]->FwdTrans(ErrVec, fieldcoeffs[2]);

    WriteFld(outname, m_fields[0], fieldcoeffs, variables);
}

/**
 * @brief Compute \int \nabla \cdot \vec{v} to compute errror
 */
void MMFStaticDiffOp::CheckErrorCovariantDeriv(
    const Array<OneD, const NekDouble> &velphi,
    const Array<OneD, const NekDouble> &veltheta,
    const Array<OneD, const Array<OneD, NekDouble>> &velvector)
{
    int ncoeffs    = m_fields[0]->GetNcoeffs();
    int nq         = GetNpoints();

    boost::ignore_unused(velphi, veltheta);

    CheckMeshErr(m_movingframes, m_velocity);

    Array<OneD, NekDouble> tmp(nq);
    Array<OneD, NekDouble> tmpc(ncoeffs);

    // Check Coneection 1-form error
    std::cout << "=========== Starting Coneection 1-form error ===========" << std::endl;

    std::cout << "Counting Activation = "
              << 100.0 * CountActivated(m_MMFActivation) / nq << " %"
              << std::endl;

    Array<OneD, NekDouble> ExactCovDerivTheta(nq);
    Array<OneD, NekDouble> ExactCovDerivPhi(nq);

    Array<OneD, NekDouble> CovDerivTheta(nq);
    Array<OneD, NekDouble> CovDerivPhi(nq);

    Array<OneD, NekDouble> CovDerivSphTheta(nq);
    Array<OneD, NekDouble> CovDerivSphPhi(nq);

    Array<OneD, NekDouble> ErrorCovDerivTheta(nq);
    Array<OneD, NekDouble> ErrorCovDerivPhi(nq);


    Array<OneD, NekDouble> Errtheta(nq);
    Array<OneD, NekDouble> Errphi(nq);

    Array<OneD, NekDouble> ErrorCovDeriv(nq);

    Array<OneD, Array<OneD, NekDouble>> CovDerivLOCAL(m_shapedim);
    Array<OneD, Array<OneD, NekDouble>> CovDerivSph(m_shapedim);
    for (int dir=0; dir<m_shapedim; ++dir)
    {
        CovDerivLOCAL[dir] = Array<OneD, NekDouble>(nq);
        CovDerivSph[dir] = Array<OneD, NekDouble>(nq);
    }

// j=0, i=0:  \nabla_{e_1} \vec{u} \cdot e_1 = \nabla u^1 \cdot e_1 - \omega_{21} (e_1) u^2 
// j=0, i=nq: \nabla_{e_1} \vec{u} \cdot e_2 = \nabla u^2 \cdot e_1 + \omega_{21} (e_1) u^1 
// j=1, i=0:  \nabla_{e_2} \vec{u} \cdot e_1 = \nabla u^1 \cdot e_2 - \omega_{21} (e_2) u^2 
// j=1, i=nq:  \nabla_{e_2} \vec{u} \cdot e_2 = \nabla u^2 \cdot e_2 + \omega_{21} (e_2) u^1 

    // // CovDeriv Test
    // Array<OneD, Array<OneD, NekDouble>> CovDeriv(m_shapedim);
    // for (int j=0; j<m_shapedim; ++j)
    // {
    //     CovDeriv[j] = Array<OneD, NekDouble>(m_shapedim * nq);
    // }

    // Array<OneD, NekDouble> u1(nq);
    // Array<OneD, NekDouble> u2(nq);

    // ComputeCovDeriv(u1, u2, movingframes, CovDeriv);
    
    for (int dir=0; dir<m_shapedim; ++dir)
    {
        std::cout << "alpha = " << m_alpha << ", direction = " << dir << std::endl;

        EvaluateExactCovariantDeriv(dir, ExactCovDerivTheta, ExactCovDerivPhi);
        std::cout << "ExactCov = ( " << RootMeanSquare(ExactCovDerivTheta, m_MMFActivation) << " , " << RootMeanSquare(ExactCovDerivPhi, m_MMFActivation) << " ) " << std::endl << std::endl;

        // spherical
        ComputeCovDeriv(velvector, m_sphereMF[dir], m_sphereMF, CovDerivLOCAL);
        MF_to_Sph(CovDerivLOCAL, m_sphereMF, CovDerivSph);

        ErrorCovDeriv = ComputeVecError(ExactCovDerivTheta, ExactCovDerivPhi, CovDerivSph, Errtheta, Errphi);

        std::cout << "SphCov = ( " << RootMeanSquare(CovDerivSph[0], m_MMFActivation) << " , " << RootMeanSquare(CovDerivSph[1], m_MMFActivation) << " ) " << std::endl;
        std::cout << "Error = " << RootMeanSquare(ErrorCovDeriv, m_MMFActivation) << ", SphinvCov = ( " << RootMeanSquare(Errtheta, m_MMFActivation) << " , " << RootMeanSquare(Errphi, m_MMFActivation) << " ) " << std::endl << std::endl;

        // Local
        ComputeCovDeriv(velvector, m_sphereMF[dir], m_movingframes, CovDerivLOCAL);
        MF_to_Sph(CovDerivLOCAL, m_movingframes, CovDerivSph);
 
        ErrorCovDeriv = ComputeVecError(ExactCovDerivTheta, ExactCovDerivPhi, CovDerivSph, Errtheta, Errphi);

        PlotCovError(dir, Errtheta, Errphi, ErrorCovDeriv);

        std::cout << "LOCALCov_sph = ( " << RootMeanSquare(CovDerivSph[0], m_MMFActivation) << " , " << RootMeanSquare(CovDerivSph[1], m_MMFActivation) << " ) " << std::endl;
        std::cout << "Error = " << RootMeanSquare(ErrorCovDeriv, m_MMFActivation) << ", SphinvCov = ( " << RootMeanSquare(Errtheta, m_MMFActivation) << " , " << RootMeanSquare(Errphi, m_MMFActivation) << " ) " << std::endl << std::endl;
    }

    wait_on_enter();
}

/**
 * @brief Compute \int \nabla \cdot \vec{v} to compute errror
 */
void MMFStaticDiffOp::CheckErrorRelAcc(
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes)
{
    int ncoeffs    = m_fields[0]->GetNcoeffs();
    int nq         = GetNpoints();

    // m_MMFActivation = Array<OneD, int>(nq, 1);

    CheckMeshErr(m_movingframes, m_velocity);

    Array<OneD, NekDouble> tmp(nq);
    Array<OneD, NekDouble> tmpc(ncoeffs);

    // Check Coneection 1-form error
    std::cout << "=========== Starting Coneection 1-form error ===========" << std::endl;

    std::cout << "Counting Activation = "
              << 100.0 * CountActivated(m_MMFActivation) / nq << " %"
              << std::endl;

    // Compute Wavefront Relative Acceleration
    // Array<OneD, NekDouble> RelAcc(nq);

    // RelAcc = ComputeRelAcc(movingframes);

    // ComputeSecondCovDeriv(0, 0, 1, movingframes, RelAcc);
    // ComputeSecondCovDeriv(1, 1, 0, movingframes, SecondCovDeriv110);

    // std::cout << "Rel Acc 001 = ( " << RootMeanSquare(SecondCovDeriv001[0],m_MMFActivation) 
    // << " , " << RootMeanSquare(SecondCovDeriv001[1],m_MMFActivation) << " ) " << std::endl;

    // std::cout << "Rel Acc 110 = ( " << RootMeanSquare(SecondCovDeriv110[0],m_MMFActivation) 
    // << " , " << RootMeanSquare(SecondCovDeriv110[1],m_MMFActivation) << " ) " << std::endl;


    Array<OneD, NekDouble> ExactRelAcc(nq);
    ExactRelAcc = ComputeExactRelAcc();

    std::cout << "ExactRelAcc = " << RootMeanSquare(ExactRelAcc,m_MMFActivation) << std::endl;

    Array<OneD, Array<OneD, NekDouble>> RelAcc001(m_shapedim);
    Array<OneD, Array<OneD, NekDouble>> RelAccbr(m_shapedim);
    Array<OneD, Array<OneD, NekDouble>> RelAccCNT(m_shapedim);
    for (int j=0; j<m_shapedim; ++j)
    {
        RelAcc001[j] = Array<OneD, NekDouble>(nq,0.0);
        RelAccbr[j] = Array<OneD, NekDouble>(nq,0.0);
        RelAccCNT[j] = Array<OneD, NekDouble>(nq,0.0);
    }

    ComputeSecondCovDeriv(0, 0, 1, movingframes, RelAcc001);
    Array<OneD, NekDouble> RelAcc001Err(nq);

    Vmath::Vsub(nq, &RelAcc001[1][0], 1, &ExactRelAcc[0], 1, &RelAcc001Err[0], 1);
    std::cout << "Rel Acc 001 err = " << RootMeanSquare(RelAcc001Err,m_MMFActivation) 
    << " (" <<  RootMeanSquare(RelAcc001[0],m_MMFActivation) << " , " 
    << RootMeanSquare(RelAcc001[1],m_MMFActivation) << ") " << std::endl;

    ComputeRelacc(movingframes, RelAccbr);

    Array<OneD, NekDouble> RelAccErrbr(nq);

    Vmath::Vsub(nq, &RelAccbr[1][0], 1, &ExactRelAcc[0], 1, &RelAccErrbr[0], 1);
    std::cout << "Rel Acc br err = " << RootMeanSquare(RelAccErrbr,m_MMFActivation) 
    << " (" <<  RootMeanSquare(RelAccbr[0],m_MMFActivation) << " , " 
    << RootMeanSquare(RelAccbr[1],m_MMFActivation) << ") " << std::endl;

    ComputeRelaccOmega(movingframes, RelAccCNT);

    Array<OneD, NekDouble> RelAccErrCNT(nq);

    Vmath::Vsub(nq, &RelAccCNT[1][0], 1, &ExactRelAcc[0], 1, &RelAccErrCNT[0], 1);
    std::cout << "Rel Acc CNT err = " << RootMeanSquare(RelAccErrCNT,m_MMFActivation) 
    << " (" <<  RootMeanSquare(RelAccCNT[0],m_MMFActivation) << " , " 
    <<  RootMeanSquare(RelAccCNT[1],m_MMFActivation) <<  ") " << std::endl;


    wait_on_enter();

    // Array<OneD, NekDouble> ExactCovDerivTheta(nq);
    // Array<OneD, NekDouble> ExactCovDerivPhi(nq);

    // Array<OneD, NekDouble> CovDerivTheta(nq);
    // Array<OneD, NekDouble> CovDerivPhi(nq);

    // Array<OneD, NekDouble> CovDerivSphTheta(nq);
    // Array<OneD, NekDouble> CovDerivSphPhi(nq);

    // Array<OneD, NekDouble> ErrorCovDerivTheta(nq);
    // Array<OneD, NekDouble> ErrorCovDerivPhi(nq);

    // Array<OneD, NekDouble> Errtheta(nq);
    // Array<OneD, NekDouble> Errphi(nq);

    // Array<OneD, NekDouble> ErrorCovDeriv(nq);

    // Array<OneD, Array<OneD, NekDouble>> CovDerivLOCAL(m_shapedim);
    // Array<OneD, Array<OneD, NekDouble>> CovDerivSph(m_shapedim);
    // for (int dir=0; dir<m_shapedim; ++dir)
    // {
    //     CovDerivLOCAL[dir] = Array<OneD, NekDouble>(nq);
    //     CovDerivSph[dir] = Array<OneD, NekDouble>(nq);
    // }

    // for (int dir=0; dir<m_shapedim; ++dir)
    // {
    //     std::cout << "alpha = " << m_alpha << ", direction = " << dir << std::endl;

    //     EvaluateExactCovariantDeriv(dir, ExactCovDerivTheta, ExactCovDerivPhi);
    //     std::cout << "ExactCov = ( " << RootMeanSquare(ExactCovDerivTheta, m_MMFActivation) << " , " << RootMeanSquare(ExactCovDerivPhi, m_MMFActivation) << " ) " << std::endl << std::endl;

    //     // spherical

    //     // ComputeCovDeriv(velvector, m_sphereMF[dir], m_sphereMF, CovDerivLOCAL);
    //     MF_to_Sph(CovDerivLOCAL, m_sphereMF, CovDerivSph);

    //     ErrorCovDeriv = ComputeVecError(ExactCovDerivTheta, ExactCovDerivPhi, CovDerivSph, Errtheta, Errphi);

    //     std::cout << "SphCov = ( " << RootMeanSquare(CovDerivSph[0], m_MMFActivation) << " , " << RootMeanSquare(CovDerivSph[1], m_MMFActivation) << " ) " << std::endl;
    //     std::cout << "Error = " << RootMeanSquare(ErrorCovDeriv, m_MMFActivation) << ", SphinvCov = ( " << RootMeanSquare(Errtheta, m_MMFActivation) << " , " << RootMeanSquare(Errphi, m_MMFActivation) << " ) " << std::endl << std::endl;

    //     // Local
    //    // ComputeCovDeriv(velvector, m_sphereMF[dir], m_movingframes, CovDerivLOCAL);
    //     MF_to_Sph(CovDerivLOCAL, m_movingframes, CovDerivSph);
 
    //     ErrorCovDeriv = ComputeVecError(ExactCovDerivTheta, ExactCovDerivPhi, CovDerivSph, Errtheta, Errphi);

    //     PlotCovError(dir, Errtheta, Errphi, ErrorCovDeriv);

    //     std::cout << "LOCALCov_sph = ( " << RootMeanSquare(CovDerivSph[0], m_MMFActivation) << " , " << RootMeanSquare(CovDerivSph[1], m_MMFActivation) << " ) " << std::endl;
    //     std::cout << "Error = " << RootMeanSquare(ErrorCovDeriv, m_MMFActivation) << ", SphinvCov = ( " << RootMeanSquare(Errtheta, m_MMFActivation) << " , " << RootMeanSquare(Errphi, m_MMFActivation) << " ) " << std::endl << std::endl;
    // }

}


/**
 * @brief Compute \int \nabla \cdot \vec{v} to compute errror
 */
void MMFStaticDiffOp::CheckErrorRiemCrv(
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes)
{
    int ncoeffs    = m_fields[0]->GetNcoeffs();
    int nq         = GetNpoints();

    // m_MMFActivation = Array<OneD, int>(nq, 1);

    CheckMeshErr(m_movingframes, m_velocity);

    Array<OneD, NekDouble> tmp(nq);
    Array<OneD, NekDouble> tmpc(ncoeffs);

    // Check Coneection 1-form error
    std::cout << "=========== Starting Coneection 1-form error ===========" << std::endl;

    std::cout << "Counting Activation = "
              << 100.0 * CountActivated(m_MMFActivation) / nq << " %"
              << std::endl;


    // Test Riemann Curvature
    Array<OneD, Array<OneD, NekDouble>> R212(m_shapedim);
    Array<OneD, Array<OneD, NekDouble>> R121(m_shapedim);
    for (int j=0; j<m_shapedim; ++j)
    {
        R212[j] = Array<OneD, NekDouble>(nq,0.0);
        R121[j] = Array<OneD, NekDouble>(nq,0.0);
    }

    // Check LieBracket
    Array<OneD, Array<OneD, NekDouble>> LieBracket121(m_shapedim);
    Array<OneD, Array<OneD, NekDouble>> LieBracket212(m_shapedim);
    for (int j=0; j<m_shapedim; ++j)
    {
        LieBracket121[j] = Array<OneD, NekDouble>(nq,0.0);
        LieBracket212[j] = Array<OneD, NekDouble>(nq,0.0);
    }

    ComputeLieBracket(0,1,0,movingframes,LieBracket121);
    ComputeLieBracket(1,0,1,movingframes,LieBracket212);

    Array<OneD, NekDouble> ExactLB1212(nq);
    Array<OneD, NekDouble> ExactLB2121(nq);
    ComputeExactLieBracket(ExactLB1212, ExactLB2121);

    Array<OneD, NekDouble> LB1212err(nq,0.0); 
    Array<OneD, NekDouble> LB2121err(nq,0.0);

    Vmath::Vsub(nq, &LieBracket121[1][0], 1, &ExactLB2121[0], 1, &LB2121err[0], 1);
    Vmath::Vsub(nq, &LieBracket212[0][0], 1, &ExactLB1212[0], 1, &LB1212err[0], 1);

    std::cout << "LB1212err = " << RootMeanSquare(LB1212err,m_MMFActivation) << " (" << RootMeanSquare(LieBracket121[1],m_MMFActivation);
    std::cout << "), LB2121err = " << RootMeanSquare(LB2121err,m_MMFActivation) << " (" << RootMeanSquare(LieBracket212[0],m_MMFActivation) << ") " << std::endl;

    // Compute Riemann Curvature tensor
    ComputeRiemCrv(0, 1, 0, movingframes, R121);
    ComputeRiemCrv(1, 0, 1, movingframes, R212);

    Array<OneD, NekDouble> ExactR1212(nq);
    Array<OneD, NekDouble> ExactR2121(nq);
    ComputeExactRiemanCvt(ExactR1212, ExactR2121);

    Array<OneD, NekDouble> R1212err(nq,0.0); 
    Array<OneD, NekDouble> R2121err(nq,0.0);

    Vmath::Vsub(nq, &R121[1][0], 1, &ExactR2121[0], 1, &R2121err[0], 1);
    Vmath::Vsub(nq, &R212[0][0], 1, &ExactR1212[0], 1, &R1212err[0], 1);

    std::cout << "R1212err = " << RootMeanSquare(R1212err,m_MMFActivation) << " (" << RootMeanSquare(R212[0],m_MMFActivation);
    std::cout << "), R2121err = " << RootMeanSquare(R2121err,m_MMFActivation) << " (" << RootMeanSquare(R121[1],m_MMFActivation) << ") " << std::endl;

    // // Compute Wavefront Relative Acceleration
    // Array<OneD, Array<OneD, NekDouble>> SecondCovDeriv001(m_shapedim);
    // Array<OneD, Array<OneD, NekDouble>> SecondCovDeriv110(m_shapedim);
    // for (int j=0; j<m_shapedim; ++j)
    // {
    //     SecondCovDeriv001[j] = Array<OneD, NekDouble>(nq,0.0);
    //     SecondCovDeriv110[j] = Array<OneD, NekDouble>(nq,0.0);
    // }

    // ComputeSecondCovDeriv(0, 0, 1, movingframes, SecondCovDeriv001);
    // ComputeSecondCovDeriv(1, 1, 0, movingframes, SecondCovDeriv110);

    // std::cout << "Rel Acc 001 = ( " << RootMeanSquare(SecondCovDeriv001[0],m_MMFActivation) 
    // << " , " << RootMeanSquare(SecondCovDeriv001[1],m_MMFActivation) << " ) " << std::endl;

    // std::cout << "Rel Acc 110 = ( " << RootMeanSquare(SecondCovDeriv110[0],m_MMFActivation) 
    // << " , " << RootMeanSquare(SecondCovDeriv110[1],m_MMFActivation) << " ) " << std::endl;

    // Array<OneD, NekDouble> RelAcc(nq);
    // RelAcc = ComputeRelAcc(movingframes);

    // Array<OneD, NekDouble> ExactRelAcc(nq);
    // ExactRelAcc = ComputeExactRelAcc();

    // Array<OneD, NekDouble> RelAccErr(nq);

    // Vmath::Vsub(nq, &RelAcc[0], 1, &ExactRelAcc[0], 1, &RelAccErr[0], 1);

    // std::cout << "Rel Acc err = " << RootMeanSquare(RelAccErr,m_MMFActivation) << " (" << RootMeanSquare(RelAcc,m_MMFActivation) << ") " << std::endl;

    wait_on_enter();

    // Array<OneD, NekDouble> ExactCovDerivTheta(nq);
    // Array<OneD, NekDouble> ExactCovDerivPhi(nq);

    // Array<OneD, NekDouble> CovDerivTheta(nq);
    // Array<OneD, NekDouble> CovDerivPhi(nq);

    // Array<OneD, NekDouble> CovDerivSphTheta(nq);
    // Array<OneD, NekDouble> CovDerivSphPhi(nq);

    // Array<OneD, NekDouble> ErrorCovDerivTheta(nq);
    // Array<OneD, NekDouble> ErrorCovDerivPhi(nq);

    // Array<OneD, NekDouble> Errtheta(nq);
    // Array<OneD, NekDouble> Errphi(nq);

    // Array<OneD, NekDouble> ErrorCovDeriv(nq);

    // Array<OneD, Array<OneD, NekDouble>> CovDerivLOCAL(m_shapedim);
    // Array<OneD, Array<OneD, NekDouble>> CovDerivSph(m_shapedim);
    // for (int dir=0; dir<m_shapedim; ++dir)
    // {
    //     CovDerivLOCAL[dir] = Array<OneD, NekDouble>(nq);
    //     CovDerivSph[dir] = Array<OneD, NekDouble>(nq);
    // }

    // for (int dir=0; dir<m_shapedim; ++dir)
    // {
    //     std::cout << "alpha = " << m_alpha << ", direction = " << dir << std::endl;

    //     EvaluateExactCovariantDeriv(dir, ExactCovDerivTheta, ExactCovDerivPhi);
    //     std::cout << "ExactCov = ( " << RootMeanSquare(ExactCovDerivTheta, m_MMFActivation) << " , " << RootMeanSquare(ExactCovDerivPhi, m_MMFActivation) << " ) " << std::endl << std::endl;

    //     // spherical

    //     // ComputeCovDeriv(velvector, m_sphereMF[dir], m_sphereMF, CovDerivLOCAL);
    //     MF_to_Sph(CovDerivLOCAL, m_sphereMF, CovDerivSph);

    //     ErrorCovDeriv = ComputeVecError(ExactCovDerivTheta, ExactCovDerivPhi, CovDerivSph, Errtheta, Errphi);

    //     std::cout << "SphCov = ( " << RootMeanSquare(CovDerivSph[0], m_MMFActivation) << " , " << RootMeanSquare(CovDerivSph[1], m_MMFActivation) << " ) " << std::endl;
    //     std::cout << "Error = " << RootMeanSquare(ErrorCovDeriv, m_MMFActivation) << ", SphinvCov = ( " << RootMeanSquare(Errtheta, m_MMFActivation) << " , " << RootMeanSquare(Errphi, m_MMFActivation) << " ) " << std::endl << std::endl;

    //     // Local
    //    // ComputeCovDeriv(velvector, m_sphereMF[dir], m_movingframes, CovDerivLOCAL);
    //     MF_to_Sph(CovDerivLOCAL, m_movingframes, CovDerivSph);
 
    //     ErrorCovDeriv = ComputeVecError(ExactCovDerivTheta, ExactCovDerivPhi, CovDerivSph, Errtheta, Errphi);

    //     PlotCovError(dir, Errtheta, Errphi, ErrorCovDeriv);

    //     std::cout << "LOCALCov_sph = ( " << RootMeanSquare(CovDerivSph[0], m_MMFActivation) << " , " << RootMeanSquare(CovDerivSph[1], m_MMFActivation) << " ) " << std::endl;
    //     std::cout << "Error = " << RootMeanSquare(ErrorCovDeriv, m_MMFActivation) << ", SphinvCov = ( " << RootMeanSquare(Errtheta, m_MMFActivation) << " , " << RootMeanSquare(Errphi, m_MMFActivation) << " ) " << std::endl << std::endl;
    // }

}

Array<OneD, NekDouble> MMFStaticDiffOp::ComputeExactRelAcc()
{
    int nq         = GetNpoints();

    Array<OneD, NekDouble> outarray(nq);

    Array<OneD, NekDouble> x(nq);
    Array<OneD, NekDouble> y(nq);
    Array<OneD, NekDouble> z(nq);

    m_fields[0]->GetCoords(x, y, z);

    NekDouble xp, yp, zp, theta;
    NekDouble sin_varphi, cos_varphi, sin_theta, cos_theta, sech_theta, tanh_theta;

    for (int i=0; i<nq; ++i)
    {
        xp = x[i];
        yp = y[i];
        zp = z[i];

        if( (m_surfaceType==ePlane) || (m_surfaceType==ePolar) )
        {
            CartesianToNewSpherical(xp, yp, zp, sin_varphi, cos_varphi,
                                    sin_theta, cos_theta);

            // outarray[i] = 2.0/rad/rad;
            outarray[i] = 0.0;
        }

        else if(m_surfaceType==eSphere)
        {
            CartesianToNewSpherical(xp, yp, zp, sin_varphi, cos_varphi,
                                    sin_theta, cos_theta);

            outarray[i] = -1.0;
        }

        else if(m_surfaceType==ePseudosphere)
        {
            CartesianToPseudospherical(xp, yp, zp, sin_varphi, cos_varphi,
                            theta, sech_theta, tanh_theta);

            outarray[i] = 1.0;
        }
    }

    return outarray;
}


void MMFStaticDiffOp::ComputeExactRiemanCvt(
    Array<OneD, NekDouble> &ExactR1212, 
    Array<OneD, NekDouble> &ExactR2121)
{
    int nq         = GetNpoints();

    Array<OneD, NekDouble> x(nq);
    Array<OneD, NekDouble> y(nq);
    Array<OneD, NekDouble> z(nq);

    m_fields[0]->GetCoords(x, y, z);

    NekDouble xp, yp, zp, rad;
    NekDouble sin_varphi, cos_varphi, sin_theta, cos_theta;

    for (int i=0; i<nq; ++i)
    {
        xp = x[i];
        yp = y[i];
        zp = z[i];

        rad = sqrt(xp*xp + yp*yp + zp*zp);

        CartesianToNewSpherical(xp, yp, zp, sin_varphi, cos_varphi,
                                sin_theta, cos_theta);

        if(m_surfaceType==ePlane)
        {
            ExactR1212[i] = 2.0/rad/rad;
            ExactR2121[i] = 0.0;
        }

        else if(m_surfaceType==eSphere)
        {
            // ExactR1212[i] = (1.0  )/(sin_theta*sin_theta);

            ExactR1212[i] = 1.0;
        }
    }
}


void MMFStaticDiffOp::ComputeExactLieBracket(
    Array<OneD, NekDouble> &ExactLB1212, 
    Array<OneD, NekDouble> &ExactLB2121)
{
    int nq         = GetNpoints();

    Array<OneD, NekDouble> x(nq);
    Array<OneD, NekDouble> y(nq);
    Array<OneD, NekDouble> z(nq);

    m_fields[0]->GetCoords(x, y, z);

    NekDouble xp, yp, zp, rad;
    NekDouble sin_varphi, cos_varphi, sin_theta, cos_theta;

    for (int i=0; i<nq; ++i)
    {
        xp = x[i];
        yp = y[i];
        zp = z[i];

        rad = sqrt(xp*xp + yp*yp + zp*zp);

        CartesianToNewSpherical(xp, yp, zp, sin_varphi, cos_varphi,
                                sin_theta, cos_theta);

        if(m_surfaceType==ePlane)
        {
            ExactLB1212[i] = -1.0/rad/rad;
            ExactLB2121[i] = -1.0/rad/rad;
        }

        else if(m_surfaceType==eSphere)
        {
            ExactLB1212[i] =  -cos_theta*cos_theta/(sin_theta*sin_theta);
            ExactLB2121[i] =  -cos_theta*cos_theta/(sin_theta*sin_theta);
        }
    }
}



// /**
//  * @brief Compute \int \nabla \cdot \vec{v} to compute errror
//  */
// void MMFStaticDiffOp::CheckErrorRiemCrvSphere(
//     const Array<OneD, const Array<OneD, NekDouble>> &movingframes)
// {
//     int ncoeffs    = m_fields[0]->GetNcoeffs();
//     int nq         = GetNpoints();

//     CheckMeshErr(m_movingframes, m_velocity);

//     Array<OneD, NekDouble> tmp(nq);
//     Array<OneD, NekDouble> tmpc(ncoeffs);

//     // Check Coneection 1-form error
//     std::cout << "=========== Starting Coneection 1-form error ===========" << std::endl;

//     std::cout << "Counting Activation = "
//               << 100.0 * CountActivated(m_MMFActivation) / nq << " %"
//               << std::endl;

//     // Test connections of moving frames
//     TestSphericalConnection1form(m_sphereMF, m_MMFActivation);

//     Array<OneD, NekDouble> ExactCovDerivTheta(nq);
//     Array<OneD, NekDouble> ExactCovDerivPhi(nq);

//     Array<OneD, NekDouble> CovDerivTheta(nq);
//     Array<OneD, NekDouble> CovDerivPhi(nq);

//     Array<OneD, NekDouble> CovDerivSphTheta(nq);
//     Array<OneD, NekDouble> CovDerivSphPhi(nq);

//     Array<OneD, NekDouble> ErrorCovDerivTheta(nq);
//     Array<OneD, NekDouble> ErrorCovDerivPhi(nq);

//     Array<OneD, NekDouble> Errtheta(nq);
//     Array<OneD, NekDouble> Errphi(nq);

//     Array<OneD, NekDouble> ErrorCovDeriv(nq);

//     Array<OneD, Array<OneD, NekDouble>> CovDerivLOCAL(m_shapedim);
//     Array<OneD, Array<OneD, NekDouble>> CovDerivSph(m_shapedim);
//     for (int dir=0; dir<m_shapedim; ++dir)
//     {
//         CovDerivLOCAL[dir] = Array<OneD, NekDouble>(nq);
//         CovDerivSph[dir] = Array<OneD, NekDouble>(nq);
//     }

//     for (int dir=0; dir<m_shapedim; ++dir)
//     {
//         std::cout << "alpha = " << m_alpha << ", direction = " << dir << std::endl;

//         EvaluateExactCovariantDeriv(dir, ExactCovDerivTheta, ExactCovDerivPhi);
//         std::cout << "ExactCov = ( " << RootMeanSquare(ExactCovDerivTheta, m_MMFActivation) << " , " << RootMeanSquare(ExactCovDerivPhi, m_MMFActivation) << " ) " << std::endl << std::endl;

//         // spherical

//         // ComputeCovDeriv(velvector, m_sphereMF[dir], m_sphereMF, CovDerivLOCAL);
//         MF_to_Sph(CovDerivLOCAL, m_sphereMF, CovDerivSph);

//         ErrorCovDeriv = ComputeVecError(ExactCovDerivTheta, ExactCovDerivPhi, CovDerivSph, Errtheta, Errphi);

//         std::cout << "SphCov = ( " << RootMeanSquare(CovDerivSph[0], m_MMFActivation) << " , " << RootMeanSquare(CovDerivSph[1], m_MMFActivation) << " ) " << std::endl;
//         std::cout << "Error = " << RootMeanSquare(ErrorCovDeriv, m_MMFActivation) << ", SphinvCov = ( " << RootMeanSquare(Errtheta, m_MMFActivation) << " , " << RootMeanSquare(Errphi, m_MMFActivation) << " ) " << std::endl << std::endl;

//         // Local
//         // ComputeCovDeriv(velvector, m_sphereMF[dir], m_movingframes, CovDerivLOCAL);
//         MF_to_Sph(CovDerivLOCAL, m_movingframes, CovDerivSph);
 
//         ErrorCovDeriv = ComputeVecError(ExactCovDerivTheta, ExactCovDerivPhi, CovDerivSph, Errtheta, Errphi);

//         PlotCovError(dir, Errtheta, Errphi, ErrorCovDeriv);

//         std::cout << "LOCALCov_sph = ( " << RootMeanSquare(CovDerivSph[0], m_MMFActivation) << " , " << RootMeanSquare(CovDerivSph[1], m_MMFActivation) << " ) " << std::endl;
//         std::cout << "Error = " << RootMeanSquare(ErrorCovDeriv, m_MMFActivation) << ", SphinvCov = ( " << RootMeanSquare(Errtheta, m_MMFActivation) << " , " << RootMeanSquare(Errphi, m_MMFActivation) << " ) " << std::endl << std::endl;
//     }

//     wait_on_enter();
// }



Array<OneD, NekDouble> MMFStaticDiffOp::ComputeVecError(
    const Array<OneD, const NekDouble> &ExactCovDerivTheta, 
    const Array<OneD, const NekDouble> &ExactCovDerivPhi,
    const Array<OneD, const Array<OneD, NekDouble>> &CovDerivSph,
    Array<OneD, NekDouble> &Errortheta,
    Array<OneD, NekDouble> &Errorphi)
{ 
    int nq         = GetNpoints();

    Array<OneD, NekDouble> ErrorCovDeriv(nq);

    Vmath::Vsub(nq, &ExactCovDerivTheta[0], 1, &CovDerivSph[0][0], 1, &Errortheta[0], 1);
    Vmath::Vsub(nq, &ExactCovDerivPhi[0], 1, &CovDerivSph[1][0], 1, &Errorphi[0], 1);
    Vmath::Vmul(nq, Errortheta, 1, Errortheta, 1, ErrorCovDeriv, 1);
    Vmath::Vvtvp(nq, Errorphi, 1, Errorphi, 1, ErrorCovDeriv, 1, ErrorCovDeriv, 1);
    Vmath::Vsqrt(nq, ErrorCovDeriv, 1, ErrorCovDeriv, 1);

    return ErrorCovDeriv;
}

/**
 * @brief Compute \int \nabla \cdot \vec{v} to compute errror
 */
void MMFStaticDiffOp::CheckErrorDivergence(
    const Array<OneD, const NekDouble> &inarray)
{
    int nvariables = m_fields.size();
    int ncoeffs    = m_fields[0]->GetNcoeffs();
    int nq         = GetNpoints();

    CheckMeshErr(m_movingframes, m_velocity);

    Array<OneD, NekDouble> tmp(nq);
    Array<OneD, NekDouble> tmpc(ncoeffs);

    //////////////////////////////////////////////////////////////////////////////////////
    // Derive exact divergence value
    Array<OneD, NekDouble> ExactDivVel(nq);
    ExactDivVel = ComputeExactValue(m_TestType);

    Array<OneD, NekDouble> ExactDivVelProj(nq);
    m_fields[0]->IProductWRTBase(ExactDivVel, tmpc);
    m_fields[0]->BwdTrans(tmpc, ExactDivVelProj);

    std::cout << "ExactDivVel = " << RootMeanSquare(ExactDivVel) << ", ExactDivVelProj = " 
    << RootMeanSquare(ExactDivVelProj) << std::endl;

    //////////////////////////////////////////////////////////////////////////////////////
    // Euclidean Direct computation

    Array<OneD, NekDouble> EuclideanDiv(nq, 0.0);
    EuclideanDiv = ComputeEuclideanDivergence(m_velocity);

    //////////////////////////////////////////////////////////////////////////////////////
    // Direct differentiation along the Spherical coordinate axis
    // DirectDiv = ComputeDivinSphericalCoord(vphi, vth);
    Array<OneD, NekDouble> DirectSpherical(nq, 0.0);
    DirectSpherical = ComputeDivSphericalCoord(m_velocity);

    //////////////////////////////////////////////////////////////////////////////////////
    // Direct differentiation with Connection
    Array<OneD, NekDouble> DirectMFDeriv(nq, 0.0);
    DirectMFDeriv = ComputeCovDiv(m_velocity, m_movingframes);

    // Array<OneD, NekDouble> DirectMFDerivold(nq, 0.0);
    // DirectMFDerivold = ComputeCovariantDivergence(m_velocity, m_movingframes);

    // Vmath::Vsub(nq, DirectMFDeriv, 1, DirectMFDerivold, 1, DirectMFDerivold, 1);

    // std::cout << "New Div error = " << RootMeanSquare(DirectMFDerivold) << std::endl;

    // Compute Div \cdot velocity vector in DG Weak form
    Array<OneD, Array<OneD, NekDouble>> WeakDGDiv(nvariables);
    ComputeWeakDGDivergence(inarray, m_movingframes, m_velocity, WeakDGDiv);

    // Compute Div \cdot velocity vector in DG Weak form
    // Spurious Divergence test

    Array<OneD, NekDouble> SpuriousDiv_pre(nq);
    Array<OneD, NekDouble> SpuriousDiv_k0(nq);
    Array<OneD, NekDouble> SpuriousDiv_D0(nq);
    Array<OneD, NekDouble> SpuriousDiv_D0_k0(nq);

    Array<OneD, NekDouble> velvector(m_spacedim * nq);
    for (int i = 0; i < m_spacedim; ++i)
    {
        Vmath::Vcopy(nq, &m_velocity[i][0], 1, &velvector[i * nq], 1);
    }

    SpuriousDiv_pre = ComputeSurfaceDiv(m_movingframes[2], velvector);

    Array<OneD, Array<OneD, NekDouble>> SphericalVector;
    ComputeSphericalVector(SphericalVector);

    SpuriousDiv_k0 = ComputeSpuriousDivergence(m_movingframes, SphericalVector[2], velvector);

    // Construct The Moving Frames
    Array<OneD, Array<OneD, NekDouble>> mf_LOCSPH;
    mf_LOCSPH = Array<OneD, Array<OneD, NekDouble>>(m_mfdim);

    for (int j = 0; j < m_mfdim; ++j)
    {
        mf_LOCSPH[j] = Array<OneD, NekDouble>(m_spacedim * nq, 0.0);
    }

    GetLOCALMovingframes(mf_LOCSPH);
    ComputeAxisAlignedLOCALMovingframes(m_sphereMF, mf_LOCSPH);
    SpuriousDiv_D0 = ComputeSpuriousDivergence(mf_LOCSPH, m_movingframes[2], velvector);
    Vmath::Neg(nq, SpuriousDiv_D0, 1);
 
    SpuriousDiv_D0_k0 = ComputeSpuriousDivergence(mf_LOCSPH, SphericalVector[2], velvector);
     Vmath::Neg(nq, SpuriousDiv_D0_k0, 1);


    std::cout << "SpuriousDiv_k0 = " << RootMeanSquare(SpuriousDiv_k0) 
    << ", SpuriousDiv_D0 = " << RootMeanSquare(SpuriousDiv_D0) 
    << ", SpuriousDiv_D0_k0 = " << RootMeanSquare(SpuriousDiv_D0_k0) << std::endl;

    Array<OneD, Array<OneD, NekDouble>> WeakDGDiv_SpDiv_pre(nvariables);
    Array<OneD, Array<OneD, NekDouble>> WeakDGDiv_SpDiv_k0(nvariables);
    Array<OneD, Array<OneD, NekDouble>> WeakDGDiv_SpDiv_D0(nvariables);
    Array<OneD, Array<OneD, NekDouble>> WeakDGDiv_SpDiv_D0_k0(nvariables);

    ComputeWeakDGDivergence(inarray, m_movingframes, m_velocity, WeakDGDiv_SpDiv_pre, 1);
    ComputeWeakDGDivergence(inarray, m_movingframes, m_velocity, WeakDGDiv_SpDiv_k0, SpuriousDiv_k0);
    ComputeWeakDGDivergence(inarray, m_movingframes, m_velocity, WeakDGDiv_SpDiv_D0, SpuriousDiv_D0);
    ComputeWeakDGDivergence(inarray, m_movingframes, m_velocity, WeakDGDiv_SpDiv_D0_k0, SpuriousDiv_D0_k0);

    // Printout differences
    std::cout
        << "****************************************************************"
        << std::endl;
    Vmath::Vsub(nq, ExactDivVel, 1, EuclideanDiv, 1, EuclideanDiv, 1);
    std::cout << "Euc Div Error, L2= "
              << RootMeanSquare(EuclideanDiv, m_MMFActivation)
              << ", Linf = " << FindAbsMaximum(EuclideanDiv, m_MMFActivation)
              << std::endl ;

    Vmath::Vsub(nq, ExactDivVel, 1, DirectSpherical, 1, DirectSpherical, 1);
    std::cout << "Direct Spherical Div Error, L2= "
              << RootMeanSquare(DirectSpherical, m_MMFActivation)
              << ", Linf = " << FindAbsMaximum(DirectSpherical, m_MMFActivation)
              << std::endl ;

    Vmath::Vsub(nq, ExactDivVel, 1, DirectMFDeriv, 1, DirectMFDeriv, 1);
    std::cout << "Direct Div Error, L2= "
              << RootMeanSquare(DirectMFDeriv, m_MMFActivation)
              << ", Linf = " << FindAbsMaximum(DirectMFDeriv, m_MMFActivation)
              << std::endl ;

    Vmath::Vsub(nq, ExactDivVelProj, 1, WeakDGDiv[0], 1, WeakDGDiv[0], 1);
    std::cout << "Weak Div Error, L2= "
              << RootMeanSquare(WeakDGDiv[0], m_MMFActivation)
              << ", Linf = " << FindAbsMaximum(WeakDGDiv[0], m_MMFActivation)
              << std::endl ;

    Vmath::Vsub(nq, ExactDivVelProj, 1, WeakDGDiv_SpDiv_pre[0], 1, WeakDGDiv_SpDiv_pre[0], 1);
    std::cout << "Weak Div with SpDiv_pre Error, L2= "
              << RootMeanSquare(WeakDGDiv_SpDiv_pre[0], m_MMFActivation)
              << ", Linf = " << FindAbsMaximum(WeakDGDiv_SpDiv_pre[0], m_MMFActivation)
              << std::endl ;

    Vmath::Vsub(nq, ExactDivVelProj, 1, WeakDGDiv_SpDiv_k0[0], 1, WeakDGDiv_SpDiv_k0[0], 1);
    std::cout << "Weak Div with SpDiv_k0 Error, L2= "
              << RootMeanSquare(WeakDGDiv_SpDiv_k0[0], m_MMFActivation)
              << ", Linf = " << FindAbsMaximum(WeakDGDiv_SpDiv_k0[0], m_MMFActivation)
              << std::endl ;

    Vmath::Vsub(nq, ExactDivVelProj, 1, WeakDGDiv_SpDiv_D0[0], 1, WeakDGDiv_SpDiv_D0[0], 1);
    std::cout << "Weak Div with SpDiv_D0 Error, L2= "
              << RootMeanSquare(WeakDGDiv_SpDiv_D0[0], m_MMFActivation)
              << ", Linf = " << FindAbsMaximum(WeakDGDiv_SpDiv_D0[0], m_MMFActivation)
              << std::endl ; 

    Vmath::Vsub(nq, ExactDivVelProj, 1, WeakDGDiv_SpDiv_D0_k0[0], 1, WeakDGDiv_SpDiv_D0_k0[0], 1);
    std::cout << "Weak Div with SpDiv_D0_k0 Error, L2= "
              << RootMeanSquare(WeakDGDiv_SpDiv_D0_k0[0], m_MMFActivation)
              << ", Linf = " << FindAbsMaximum(WeakDGDiv_SpDiv_D0_k0[0], m_MMFActivation)
              << std::endl ; 

    // NekDouble L2diff = RootMeanSquare(WeakDGDiv[0], m_MMFActivation) - RootMeanSquare(NewWeakDGDiv[0], m_MMFActivation);
    // NekDouble Linfdiff = FindAbsMaximum(WeakDGDiv[0], m_MMFActivation) - FindAbsMaximum(NewWeakDGDiv[0], m_MMFActivation);
    // std::cout << "Difference, L2= "
    //           << L2diff
    //           << ", Linf = " << Linfdiff
    //           << std::endl ;
    std::cout
        << "****************************************************************"
        << std::endl;
}


void MMFStaticDiffOp::CheckErrorCurl(
    const Array<OneD, const NekDouble> &inarrayth,
    const Array<OneD, const NekDouble> &inarrayphi)
{
    int ncoeffs = m_fields[0]->GetNcoeffs();
    int nq      = GetNpoints();

    Array<OneD, NekDouble> tmp(nq);
    Array<OneD, NekDouble> tmpc(ncoeffs);

    CheckMeshErr(m_movingframes, m_velocity);

    // V = inarrayth V_{TH} + inarrayphi V_{phi}
    // V1 = e^1 \cdot V = inarrayth ( e^1 \cdot V_{TH} ) + inarrayphi ( e^1
    // \cdot V_{PHI} ) V2 = e^2 \cdot V = inarrayth ( e^2 \cdot V_{TH} ) +
    // inarrayphi ( e^2 \cdot V_{PHI} )

    // Derive an equivalent vector for the spherical coordinate components
    Array<OneD, Array<OneD, NekDouble>> velocity;
    SphericalToEuclidean(inarrayth, inarrayphi, velocity);

    // Compute Exact Curl
    Array<OneD, NekDouble> ExactCurl(nq);
    ExactCurl = ComputeExactValue(m_TestType);
 
    std::cout << "Exact Curl = " << RootMeanSquare(ExactCurl) << std::endl;

    Array<OneD, NekDouble> ExactCurlProj(nq);
    m_fields[0]->IProductWRTBase(ExactCurl, tmpc);
    m_fields[0]->BwdTrans(tmpc, ExactCurlProj);
 
    // Compute Euclidean Curl
    Array<OneD, NekDouble> EucCurl(nq);
    EucCurl = ComputeEuclideanCurl(velocity, m_movingframes);

    // Compute the direct computation along the Spherical curved axis
    Array<OneD, NekDouble> DirectSpherical(nq);
    DirectSpherical = ComputeCurlSphericalCoord(inarrayphi, inarrayth);
 
    // Compute MMF Covariant Curl
    Array<OneD, NekDouble> MMFCurl(nq);
    MMFCurl = ComputeCovCurl(velocity, m_movingframes);

    // Array<OneD, NekDouble> MMFCurlold(nq);
    // MMFCurlold = ComputeCovariantCurl(velocity, m_movingframes);

    // Vmath::Vsub(nq, MMFCurl, 1, MMFCurlold, 1, MMFCurlold, 1);

    // std::cout << "MMFCurl error = " << RootMeanSquare(MMFCurlold) << std::endl;


     // Compute Weak Curl
    Array<OneD, NekDouble> WeakCurl(nq);
    Array<OneD, Array<OneD, NekDouble>> physfield;
    SphericalToMovingFrames(inarrayth, inarrayphi, physfield);
    WeakDGCurl(physfield, m_movingframes, m_CrossProductMF, WeakCurl, 0);

    // Compute Div \cdot velocity vector in DG Weak form
    // Spurious Divergence test

    Array<OneD, NekDouble> SpuriousDiv_pre(nq);
    Array<OneD, NekDouble> SpuriousDiv_k0(nq);
    Array<OneD, NekDouble> SpuriousDiv_D0(nq);
    Array<OneD, NekDouble> SpuriousDiv_D0_k0(nq);

    // Array<OneD, NekDouble> velvector(m_spacedim * nq);
    // for (int i = 0; i < m_spacedim; ++i)
    // {
    //     Vmath::Vcopy(nq, &m_velocity[i][0], 1, &velvector[i * nq], 1);
    // }

    Array<OneD, NekDouble> velvector(m_spacedim * nq, 0.0);
    for (int i = 0; i < m_spacedim; ++i)
    {
        for (int j = 0; j < m_shapedim; ++j)
        {
            Vmath::Vvtvp(nq, &physfield[j][0], 1,
                            &m_CrossProductMF[j][i * nq], 1, &velvector[i * nq],
                            1, &velvector[i * nq], 1);
        }
    }

    SpuriousDiv_pre = ComputeSurfaceDiv(m_movingframes[2], velvector);

    Array<OneD, Array<OneD, NekDouble>> SphericalVector;
    ComputeSphericalVector(SphericalVector);

    SpuriousDiv_k0 = ComputeSpuriousDivergence(m_movingframes, SphericalVector[2], velvector);

    // Construct The Moving Frames
    Array<OneD, Array<OneD, NekDouble>> mf_LOCSPH;
    mf_LOCSPH = Array<OneD, Array<OneD, NekDouble>>(m_mfdim);

    for (int j = 0; j < m_mfdim; ++j)
    {
        mf_LOCSPH[j] = Array<OneD, NekDouble>(m_spacedim * nq, 0.0);
    }

    GetLOCALMovingframes(mf_LOCSPH);
    ComputeAxisAlignedLOCALMovingframes(m_sphereMF, mf_LOCSPH);
    SpuriousDiv_D0 = ComputeSpuriousDivergence(mf_LOCSPH, m_movingframes[2], velvector);
    Vmath::Neg(nq, SpuriousDiv_D0, 1);


    SpuriousDiv_D0_k0 = ComputeSpuriousDivergence(mf_LOCSPH, SphericalVector[2], velvector);
    Vmath::Neg(nq, SpuriousDiv_D0_k0, 1);


    std::cout << "SpuriousDiv_k0 = " << RootMeanSquare(SpuriousDiv_k0) 
    << ", SpuriousDiv_D0 = " << RootMeanSquare(SpuriousDiv_D0) 
    << ", SpuriousDiv_D0_k0 = " << RootMeanSquare(SpuriousDiv_D0_k0) << std::endl;

    Array<OneD, NekDouble> WeakDGCurl_SpDiv_pre(nq);
    Array<OneD, NekDouble> WeakDGCurl_SpDiv_k0(nq);
    Array<OneD, NekDouble> WeakDGCurl_SpDiv_D0(nq);
    Array<OneD, NekDouble> WeakDGCurl_SpDiv_D0_k0(nq);

    // WeakDGCurl(physfield, m_movingframes, m_CrossProductMF, WeakCurl, 0);
    WeakDGCurl(physfield, m_movingframes, m_CrossProductMF, WeakDGCurl_SpDiv_pre, 1);
    WeakDGCurl(physfield, m_movingframes, m_CrossProductMF, WeakDGCurl_SpDiv_k0, SpuriousDiv_k0);
    WeakDGCurl(physfield, m_movingframes, m_CrossProductMF, WeakDGCurl_SpDiv_D0, SpuriousDiv_D0);
    WeakDGCurl(physfield, m_movingframes, m_CrossProductMF, WeakDGCurl_SpDiv_D0_k0, SpuriousDiv_D0_k0);

    // Compute Weak Curl with projected differentiation
    // Array<OneD, NekDouble> NewWeakCurl(nq);
    // SphericalToMovingFrames(inarrayth, inarrayphi, physfield);
    // WeakDGCurl(physfield, m_movingframes, m_CrossProductMF, NewWeakCurl, 1);

    std::cout
        << "****************************************************************"
        << std::endl;
    Vmath::Vsub(nq, ExactCurl, 1, EucCurl, 1, EucCurl, 1);
    std::cout << "Euc Curl Error, L2= "
              << RootMeanSquare(EucCurl, m_MMFActivation)
              << ", Linf = " << FindAbsMaximum(EucCurl, m_MMFActivation)
              << std::endl ;

    Vmath::Vsub(nq, ExactCurl, 1, DirectSpherical, 1, DirectSpherical, 1);
    std::cout << "Direct Spherical Curl Error, L2= "
              << RootMeanSquare(DirectSpherical, m_MMFActivation)
              << ", Linf = " << FindAbsMaximum(DirectSpherical, m_MMFActivation)
              << std::endl ;

    Vmath::Vsub(nq, ExactCurl, 1, MMFCurl, 1, MMFCurl, 1);
    std::cout << "Direct Curl Error, L2= "
              << RootMeanSquare(MMFCurl, m_MMFActivation)
              << ", Linf = " << FindAbsMaximum(MMFCurl, m_MMFActivation)
              << std::endl ;

    Vmath::Vsub(nq, ExactCurlProj, 1, WeakCurl, 1, WeakCurl, 1);
    std::cout << "Weak Curl Error, L2= "
              << RootMeanSquare(WeakCurl, m_MMFActivation)
              << ", Linf = " << FindAbsMaximum(WeakCurl, m_MMFActivation)
              << std::endl;

    Vmath::Vsub(nq, ExactCurlProj, 1, WeakDGCurl_SpDiv_pre, 1, WeakDGCurl_SpDiv_pre, 1);
    std::cout << "Weak Curl with SpDiv_pre Error, L2= "
              << RootMeanSquare(WeakDGCurl_SpDiv_pre, m_MMFActivation)
              << ", Linf = " << FindAbsMaximum(WeakDGCurl_SpDiv_pre, m_MMFActivation)
              << std::endl ;

    Vmath::Vsub(nq, ExactCurlProj, 1, WeakDGCurl_SpDiv_k0, 1, WeakDGCurl_SpDiv_k0, 1);
    std::cout << "Weak Curl with SpDiv_k0 Error, L2= "
              << RootMeanSquare(WeakDGCurl_SpDiv_k0, m_MMFActivation)
              << ", Linf = " << FindAbsMaximum(WeakDGCurl_SpDiv_k0, m_MMFActivation)
              << std::endl ;

    Vmath::Vsub(nq, ExactCurlProj, 1, WeakDGCurl_SpDiv_D0, 1, WeakDGCurl_SpDiv_D0, 1);
    std::cout << "Weak Curl with SpDiv_D0 Error, L2= "
              << RootMeanSquare(WeakDGCurl_SpDiv_D0, m_MMFActivation)
              << ", Linf = " << FindAbsMaximum(WeakDGCurl_SpDiv_D0, m_MMFActivation)
              << std::endl ; 

    Vmath::Vsub(nq, ExactCurlProj, 1, WeakDGCurl_SpDiv_D0_k0, 1, WeakDGCurl_SpDiv_D0_k0, 1);
    std::cout << "Weak Curl with SpDiv_D0_k0 Error, L2= "
              << RootMeanSquare(WeakDGCurl_SpDiv_D0_k0, m_MMFActivation)
              << ", Linf = " << FindAbsMaximum(WeakDGCurl_SpDiv_D0_k0, m_MMFActivation)
              << std::endl ; 

    // Vmath::Vsub(nq, ExactCurlProj, 1, NewWeakCurl, 1, NewWeakCurl, 1);
    // std::cout << "New Weak MMF Error, L2= "
    //           << RootMeanSquare(NewWeakCurl, m_MMFActivation)
    //           << ", Linf = " << FindAbsMaximum(NewWeakCurl, m_MMFActivation)
    //           << std::endl ;

    // Vmath::Vsub(nq, ExactCurl, 1, WeakCurl, 1, WeakCurl, 1);
    // std::cout << "Weak MMF Error, L2= "
    //           << RootMeanSquare(WeakCurl, m_MMFActivation)
    //           << ", Linf = " << FindAbsMaximum(WeakCurl, m_MMFActivation)
    //           << std::endl;

    // Vmath::Vsub(nq, ExactCurl, 1, NewWeakCurl, 1, NewWeakCurl, 1);
    // std::cout << "New Weak MMF Error, L2= "
    //           << RootMeanSquare(NewWeakCurl, m_MMFActivation)
    //           << ", Linf = " << FindAbsMaximum(NewWeakCurl, m_MMFActivation)
    //           << std::endl ;


    // NekDouble L2diff = RootMeanSquare(WeakCurl, m_MMFActivation) - RootMeanSquare(NewWeakCurl, m_MMFActivation);
    // NekDouble Linfdiff = FindAbsMaximum(WeakCurl, m_MMFActivation) - FindAbsMaximum(NewWeakCurl, m_MMFActivation);
    // std::cout << "Difference, L2= "
    //           << L2diff
    //           << ", Linf = " << Linfdiff
    //           << std::endl ;

    std::cout
        << "****************************************************************"
        << std::endl;
}

/**
 * @brief Compute \int \nabla \cdot \mathbf{e}^1
 * \vec{v} = \vec{e}^1
 */

void MMFStaticDiffOp::CheckErrorGradient(
    const Array<OneD, const NekDouble> &inarray)
{
    int ncoeffs    = m_fields[0]->GetNcoeffs();
    int nq         = GetNpoints();

    Array<OneD, NekDouble> tmp(nq);
    Array<OneD, NekDouble> tmpc(ncoeffs);

    CheckMeshErr(m_movingframes, m_movingframes);

    // Exact solution
    Array<OneD, NekDouble> ExactGrad(m_shapedim * nq);
    ExactGrad = ComputeExactValue(m_TestType);

    Array<OneD, NekDouble> ExactGradProj(m_shapedim*nq);
    for (int j=0; j<m_shapedim; j++)
    {
        Vmath::Vcopy(nq, &ExactGrad[j*nq], 1, &tmp[0], 1);

        m_fields[0]->IProductWRTBase(tmp, tmpc);
        m_fields[0]->BwdTrans(tmpc, tmp);

        Vmath::Vcopy(nq, &tmp[0], 1, &ExactGradProj[j*nq], 1);
    }

    Array<OneD, NekDouble> Phi(m_spacedim * nq);
    Array<OneD, NekDouble> Theta(m_spacedim * nq);

    ComputeSphericalTangentVector(Phi, Theta);

    Array<OneD, NekDouble> ExactGradVec(m_spacedim * nq,0.0);
    Array<OneD, NekDouble> ExactGradProjVec(m_spacedim * nq,0.0);
    for (int k=0; k<m_spacedim; ++k)
    {
        Vmath::Vvtvp(nq, &ExactGrad[0], 1, &Phi[k*nq], 1, &ExactGradVec[k*nq], 1, &ExactGradVec[k*nq], 1);
        Vmath::Vvtvp(nq, &ExactGrad[nq], 1, &Theta[k*nq], 1, &ExactGradVec[k*nq], 1, &ExactGradVec[k*nq], 1);

        Vmath::Vvtvp(nq, &ExactGradProj[0], 1, &Phi[k*nq], 1, &ExactGradProjVec[k*nq], 1, &ExactGradProjVec[k*nq], 1);
        Vmath::Vvtvp(nq, &ExactGradProj[nq], 1, &Theta[k*nq], 1, &ExactGradProjVec[k*nq], 1, &ExactGradProjVec[k*nq], 1);
    }

    // // Exact solution: Physcal space
    // Array<OneD, NekDouble> ExactGrad(m_shapedim * nq);
    // ExactGrad = ComputeExactValue(m_TestType);

    // Array<OneD, NekDouble> Phi(m_spacedim * nq);
    // Array<OneD, NekDouble> Theta(m_spacedim * nq);

    // ComputeSphericalTangentVector(Phi, Theta);

    // Array<OneD, NekDouble> ExactGradVec(m_spacedim * nq,0.0);
    // for (int k=0; k<m_spacedim; ++k)
    // {
    //     Vmath::Vvtvp(nq, &ExactGrad[0], 1, &Phi[k*nq], 1, &ExactGradVec[k*nq], 1, &ExactGradVec[k*nq], 1);
    //     Vmath::Vvtvp(nq, &ExactGrad[nq], 1, &Theta[k*nq], 1, &ExactGradVec[k*nq], 1, &ExactGradVec[k*nq], 1);
    // }

    // // Exact solution: modal space
    // Array<OneD, NekDouble> ExactGradProjVec(m_spacedim * nq,0.0);
    // for (int k=0; k<m_spacedim; k++)
    // {
    //     Vmath::Vcopy(nq, &ExactGradVec[k*nq], 1, &tmp[0], 1);

    //     m_fields[0]->IProductWRTBase(tmp, tmpc);
    //     m_fields[0]->BwdTrans(tmpc, tmp);

    //     Vmath::Vcopy(nq, &tmp[0], 1, &ExactGradProjVec[k*nq], 1);
    // }

    // Array<OneD, NekDouble> ProjErr(nq,0.0);
    // for (int k=0; k<m_spacedim; k++)
    // {
    //     Vmath::Vsub(nq, &ExactGradVec[k*nq], 1, &ExactGradProjVec[k*nq], 1, &tmp[0], 1);
    //     Vmath::Vmul(nq, tmp, 1, tmp, 1, tmp, 1);
    //     Vmath::Vadd(nq, tmp, 1, ProjErr, 1, ProjErr, 1);
    // }
    // Vmath::Vsqrt(nq, ProjErr, 1, ProjErr, 1);

    // std::cout << "ProjError = " << RootMeanSquare(ProjErr, m_MMFActivation) << std::endl;

    // NekDouble sin_varphi, cos_varphi, sin_theta, cos_theta;
    // NekDouble x0j, x1j, x2j, rad;

    // Array<OneD, NekDouble> x0(nq);
    // Array<OneD, NekDouble> x1(nq);
    // Array<OneD, NekDouble> x2(nq);

    // m_fields[0]->GetCoords(x0, x1, x2);
    // for (int i=0; i<nq; ++i)
    // {
    //     x0j = x0[i];
    //     x1j = x1[i];
    //     x2j = x2[i];

    //     CartesianToNewSpherical(x0j, x1j, x2j, rad, sin_varphi, cos_varphi,
    //                             sin_theta, cos_theta);

    //     if(m_MMFActivation[i])
    //     {
    //         ExactGrad[i] = ExactGrad[i]*sin_theta;
    //     }
    // }

    // Euclidean Gradient 
    Array<OneD, NekDouble> EucGrad(m_spacedim * nq, 0.0);
    EucGrad = ComputeEuclideanGradient(inarray);
 
    // Direct Gradient on the spherical coordinate axis
    Array<OneD, NekDouble> DirectSpherical;
    ComputeGradientSphericalCoord(inarray, DirectSpherical);
 
    // Covariant Gradient 
    Array<OneD, NekDouble> MMFGrad(m_spacedim * nq, 0.0);
    MMFGrad = ComputeCovGrad(inarray, m_movingframes);
 
    PlotFieldVector(inarray, MMFGrad);

    // Weak DG Gradient 
    int SurfaceGradient = 0;
    Array<OneD, NekDouble> WeakGrad;
    WeakDGGradientVector(inarray, m_movingframes, WeakGrad, SurfaceGradient);
 
    // New Weak DG Gradient 
    SurfaceGradient = 1;
    Array<OneD, NekDouble> NewWeakGrad;
    WeakDGGradientVector(inarray, m_movingframes, NewWeakGrad, SurfaceGradient);

    std::cout
        << "****************************************************************"
        << std::endl;
    Vmath::Vsub(m_spacedim*nq, ExactGradVec, 1, EucGrad, 1, EucGrad, 1);
    std::cout << "Euc Gradient Error, L2= "
              << RootMeanSquareVector(EucGrad, m_MMFActivation)
              << ", Linf = " << FindAbsMaximumVector(EucGrad, m_MMFActivation)
              << std::endl ;

    Vmath::Vsub(m_spacedim*nq, ExactGradVec, 1, DirectSpherical, 1, DirectSpherical, 1);
    std::cout << "Direct Spherical Error, L2= "
              << RootMeanSquareVector(DirectSpherical, m_MMFActivation)
              << ", Linf = " << FindAbsMaximumVector(DirectSpherical, m_MMFActivation)
              << std::endl ;

    Vmath::Vsub(m_spacedim*nq, ExactGradVec, 1, MMFGrad, 1, MMFGrad, 1);
    std::cout << "Direct MMF Error, L2= "
              << RootMeanSquareVector(MMFGrad, m_MMFActivation)
              << ", Linf = " << FindAbsMaximumVector(MMFGrad, m_MMFActivation)
              << std::endl << std::endl  ;

    Vmath::Vsub(m_spacedim*nq, ExactGradProjVec, 1, WeakGrad, 1, WeakGrad, 1);
    std::cout << "Weak MMF Error, L2= "
              << RootMeanSquareVector(WeakGrad, m_MMFActivation)
              << ", Linf = " << FindAbsMaximumVector(WeakGrad, m_MMFActivation)
              << std::endl;

    Vmath::Vsub(m_spacedim*nq, ExactGradProjVec, 1, NewWeakGrad, 1, NewWeakGrad, 1);
    std::cout << "New Weak MMF Error, L2= "
              << RootMeanSquareVector(NewWeakGrad, m_MMFActivation)
              << ", Linf = " << FindAbsMaximumVector(NewWeakGrad, m_MMFActivation)
              << std::endl;

    std::cout
        << "****************************************************************"
        << std::endl;

} // namespace Nektar


void MMFStaticDiffOp::CheckErrorLaplacian(
    const Array<OneD, const NekDouble> &inarray)
{
    int nq         = GetNpoints();
    int ncoeffs    = m_fields[0]->GetNcoeffs();

    Array<OneD, NekDouble> tmpc(ncoeffs);

    CheckMeshErr(m_movingframes, m_movingframes);
 
    // Check MMFActivaation
    Array<OneD, NekDouble> ExactLap(nq);
    ExactLap = ComputeExactValue(m_TestType);

    Array<OneD, NekDouble> ExactLapProj(nq);
    m_fields[0]->IProductWRTBase(ExactLap, tmpc);
    m_fields[0]->BwdTrans(tmpc, ExactLapProj);

    Array<OneD, NekDouble> Exactdiff(nq);
    Vmath::Vsub(nq, ExactLap, 1, ExactLapProj, 1, Exactdiff, 1);
    std::cout << "ExactLap = " << RootMeanSquare(ExactLap) << "ExactLapProj = " << RootMeanSquare(ExactLapProj) << 
    ", Exact diff = " << RootMeanSquare(Exactdiff, m_MMFActivation) << std::endl;
 
    // Compute Euclidean Diffusion
    Array<OneD, NekDouble> EuclideanLap(nq);
    EuclideanLap = ComputeEuclideanDiffusion(inarray);
 
    // Compute Direct Laplacian in the Spherical coordinate system
    Array<OneD, NekDouble> SphericalLap(nq);
    SphericalLap = ComputeLaplacianSphericalCoord(inarray);
 
    // Compute Covariant Diffusion
    Array<OneD, NekDouble> CovariantLap(nq);
    CovariantLap = ComputeCovariantDiffusion(m_movingframes, inarray, m_DerivType);

     // WeakMMFLaplacian
    // int SurfaceLaplacian;
     
    //SurfaceLaplacian = 0;
    Array<OneD, NekDouble> WeakMMFLap(nq);
    WeakDGMMFLaplacian(0, inarray, WeakMMFLap);
 
    // WeakMMFLaplacian with Surface integration
   // SurfaceLaplacian = 1;
    Array<OneD, NekDouble> NewWeakMMFLap(nq);
    WeakDGMMFLaplacian(0, inarray, NewWeakMMFLap);
 
    std::cout
    << "****************************************************************"
    << std::endl;

    Vmath::Vsub(nq, ExactLap, 1, EuclideanLap, 1, EuclideanLap, 1);
    std::cout << "Euc Lap Error, L2= "
              << RootMeanSquare(EuclideanLap, m_MMFActivation)
              << ", Linf = " << FindAbsMaximum(EuclideanLap, m_MMFActivation)
              << std::endl ;

    Vmath::Vsub(nq, ExactLap, 1, SphericalLap, 1, SphericalLap, 1);
    std::cout << "Direct Spherical Lap Error, L2= "
              << RootMeanSquare(SphericalLap, m_MMFActivation)
              << ", Linf = " << FindAbsMaximum(SphericalLap, m_MMFActivation)
              << std::endl ;

    Vmath::Vsub(nq, ExactLap, 1, CovariantLap, 1, CovariantLap, 1);
    std::cout << "Direct MMF Lap Error, L2= "
              << RootMeanSquare(CovariantLap, m_MMFActivation)
              << ", Linf = " << FindAbsMaximum(CovariantLap, m_MMFActivation)
              << std::endl ;

    Vmath::Vsub(nq, ExactLapProj, 1, WeakMMFLap, 1, WeakMMFLap, 1);
    std::cout << "Weak MMF Lap Error, L2= "
              << RootMeanSquare(WeakMMFLap, m_MMFActivation)
              << ", Linf = " << FindAbsMaximum(WeakMMFLap, m_MMFActivation)
              << std::endl;

    Vmath::Vsub(nq, ExactLapProj, 1, NewWeakMMFLap, 1, NewWeakMMFLap, 1);
    std::cout << "New Weak MMF Error, L2= "
              << RootMeanSquare(NewWeakMMFLap, m_MMFActivation)
              << ", Linf = " << FindAbsMaximum(NewWeakMMFLap, m_MMFActivation)
              << std::endl ;

    std::cout
        << "****************************************************************"
        << std::endl;


}


void MMFStaticDiffOp::CheckErrorHelmSolve(
    const Array<OneD, const NekDouble> &inarray)
{
    int nvariables = m_fields.size();
    int ncoeffs    = m_fields[0]->GetNcoeffs();
    int nq         = GetNpoints();

    // // Compute Div \cdot velocity vector
    Array<OneD, Array<OneD, NekDouble>> Y(nvariables);

    // Get the variables in physical space
    // already in physical space
    for (int i = 0; i < nvariables; ++i)
    {
        Y[i] = Array<OneD, NekDouble>(nq);
    }

    // Check MMFActivaation
    Array<OneD, NekDouble> Yexact(nq);
    Yexact = ComputeExactValue(m_TestType);

    // Compute Weak Eucliean Laplacian
    Array<OneD, Array<OneD, NekDouble>> Yarray(nvariables);
    Array<OneD, Array<OneD, NekDouble>> WeakEuclidean(nvariables);
    Array<OneD, Array<OneD, NekDouble>> DirectEuclidean(nvariables);
    Array<OneD, Array<OneD, NekDouble>> WeakMMF(nvariables);
    for (int i = 0; i < nvariables; ++i)
    {
        Yarray[i]          = Array<OneD, NekDouble>(nq, 0.0);
        WeakEuclidean[i]   = Array<OneD, NekDouble>(nq, 0.0);
        DirectEuclidean[i] = Array<OneD, NekDouble>(nq, 0.0);
        WeakMMF[i]         = Array<OneD, NekDouble>(nq, 0.0);
    }

    Vmath::Vcopy(nq, &inarray[0], 1, &Yarray[0][0], 1);

    ExplicitDiffusion(nvariables, m_fields, Yarray, WeakEuclidean);

    // WeakDGLaplacian(m_movingframes, inarray, eCovariant, WeakMMF[0]);

    std::cout << "End of WeakDGLaplacian" << std::endl;

    // Helmsolve
    NekDouble tau = 1.0;
    Array<OneD, NekDouble> HmMMFsol(nq);
    Array<OneD, NekDouble> HmEucsol(nq);

    Array<OneD, NekDouble> tmpc(ncoeffs);

    StdRegions::VarCoeffMap varcoeff;
    ComputeVarCoeff(m_movingframes, varcoeff);

    StdRegions::ConstFactorMap factors;
    factors[StdRegions::eFactorTau]    = tau;
    factors[StdRegions::eFactorLambda] = 0.0;

    m_fields[0]->HelmSolve(Yexact, tmpc, factors, varcoeff);
    m_fields[0]->BwdTrans(tmpc, HmMMFsol);

    std::cout << "End of HelmSolveMMF" << std::endl;

    Array<OneD, NekDouble> DirEuc(nq);
    DirEuc = ComputeEuclideanDiffusion(inarray);

    Array<OneD, NekDouble> DirCov(nq);
    DirCov = ComputeCovariantDiffusion(m_movingframes, inarray);

    Array<OneD, Array<OneD, NekDouble>> qfieldexact(m_shapedim);
    for (int i = 0; i < m_shapedim; ++i)
    {
        qfieldexact[i] = Array<OneD, NekDouble>(nq);
    }

    Array<OneD, NekDouble> x0(nq);
    Array<OneD, NekDouble> x1(nq);
    Array<OneD, NekDouble> x2(nq);

    m_fields[0]->GetCoords(x0, x1, x2);

    NekDouble xp, yp, zp;
    NekDouble sin_varphi, cos_varphi, sin_theta, cos_theta, cos2phi, sin2phi;
    for (int i = 0; i < nq; ++i)
    {
        xp = x0[i];
        yp = x1[i];
        zp = x2[i];

        CartesianToNewSpherical(xp, yp, zp, sin_varphi, cos_varphi,
                                sin_theta, cos_theta);

        cos2phi = 2.0 * cos_varphi * cos_varphi - 1.0;
        sin2phi = 2.0 * cos_varphi * sin_varphi;

        qfieldexact[0][i] = 2.0 * m_W * sin_theta * cos_theta * cos2phi;
        qfieldexact[1][i] = -2.0 * m_W * sin_theta * sin2phi;
    }
}


/**
 * @brief Compute the right-hand side for the linear advection equation.
 *
 * @param inarray    Given fields.
 * @param outarray   Calculated solution.
 * @param time       Time.
 */
void MMFStaticDiffOp::DoOdeRhs(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time)
{
    int i;
    int nvariables = inarray.size();
    int npoints    = GetNpoints();

    switch (m_projectionType)
    {
        case MultiRegions::eDiscontinuous:
        {
            int ncoeffs = inarray[0].size();
            Array<OneD, Array<OneD, NekDouble>> WeakAdv(nvariables);

            WeakAdv[0] = Array<OneD, NekDouble>(ncoeffs * nvariables);
            for (i = 1; i < nvariables; ++i)
            {
                WeakAdv[i] = WeakAdv[i - 1] + ncoeffs;
            }

            if ((time > 0) && (m_spacedim == 3))
            {
                // Compute \nabla \cdot \vel u according to MMF scheme
                // WeakDGDivergence(inarray, movingframes_VelMF, m_velocity,
                //              WeakAdv);
            }
            else
            {
                // Compute \nabla \cdot \vel u in the Cartesian coordinate
                // WeakDGAdvection(inarray,WeakAdv,false,true);
            }

            for (i = 0; i < nvariables; ++i)
            {
                m_fields[i]->MultiplyByElmtInvMass(WeakAdv[i], WeakAdv[i]);
                m_fields[i]->BwdTrans(WeakAdv[i], outarray[i]);
                Vmath::Neg(npoints, outarray[i], 1);
            }
        }
        break;

        default:
        {
            ASSERTL0(false, "Unknown projection scheme");
        }
        break;
    }
}
/**
 * @brief Compute the projection for the linear advection equation.
 *
 * @param inarray    Given fields.
 * @param outarray   Calculated solution.
 * @param time       Time.
 */
void MMFStaticDiffOp::DoOdeProjection(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time)
{
    // Counter variable
    int i;

    // Number of fields (variables of the problem)
    int nVariables = inarray.size();

    // Set the boundary conditions
    SetBoundaryConditions(time);

    // Switch on the projection type (Discontinuous or Continuous)
    switch (m_projectionType)
    {
        // Discontinuous projection
        case MultiRegions::eDiscontinuous:
        {
            // Number of quadrature points
            int nQuadraturePts = GetNpoints();

            // Just copy over array
            for (i = 0; i < nVariables; ++i)
            {
                Vmath::Vcopy(nQuadraturePts, inarray[i], 1, outarray[i], 1);
            }
            break;
        }

        // Continuous projection
        case MultiRegions::eGalerkin:
        case MultiRegions::eMixed_CG_Discontinuous:
        {
            Array<OneD, NekDouble> coeffs(m_fields[0]->GetNcoeffs(), 0.0);
            for (i = 0; i < nVariables; ++i)
            {
                m_fields[0]->FwdTransLocalElmt(inarray[i], coeffs);
                m_fields[0]->BwdTrans(coeffs, outarray[i]);
            }
            break;
        }

        default:
            ASSERTL0(false, "Unknown projection scheme");
            break;
    }
}

/**
 * @brief Return the flux vector for the linear advection equation.
 *
 * @param i           Component of the flux vector to calculate.
 * @param physfield   Fields.
 * @param flux        Resulting flux.
 */
// void MMFStaticDiffOp::GetFluxVector(
//     const Array<OneD, Array<OneD, NekDouble>> &physfield,
//     Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &flux)
// {
//     ASSERTL1(flux[0].size() == m_velocity.size(),
//              "Dimension of flux array and velocity array do not match");

//     int i, j;
//     int nq = physfield[0].size();

//     for (i = 0; i < flux.size(); ++i)
//     {
//         for (j = 0; j < flux[0].size(); ++j)
//         {
//             Vmath::Vmul(nq, physfield[i], 1, m_velocity[j], 1, flux[i][j],
//             1);
//         }
//     }
// }


void MMFStaticDiffOp::SetBoundaryConditions(NekDouble time)
{
    int i;
    std::string varName;
    int nvariables = m_fields.size();
    int cnt        = 0;
    int nTracePts  = GetTraceTotPoints();
    int nq         = m_fields[0]->GetNpoints();

    for (int n = 0; n < m_fields[0]->GetBndConditions().size(); ++n)
    {

        if (m_fields[0]->GetBndConditions()[n]->IsTimeDependent())
        {
            for (int i = 0; i < nvariables; ++i)
            {
                varName = m_session->GetVariable(i);
                m_fields[i]->EvaluateBoundaryConditions(time, varName);
            }
        }

        else
        {
            Array<OneD, Array<OneD, NekDouble>> Exactsolution(nvariables);
            Exactsolution[0] = SetInitialU();
            for (i = 1; i < nvariables; ++i)
            {
                Exactsolution[i] = Array<OneD, NekDouble>(nq, 0.0);
            }

            Array<OneD, Array<OneD, NekDouble>> Fwd(nvariables);
            for (i = 0; i < nvariables; ++i)
            {
                Fwd[i] = Array<OneD, NekDouble>(nTracePts);
                m_fields[i]->ExtractTracePhys(Exactsolution[i], Fwd[i]);
            }

            int e, id1, id2, npts;

            for (e = 0;
                 e < m_fields[0]->GetBndCondExpansions()[n]->GetExpSize(); ++e)
            {
                npts = m_fields[0]
                           ->GetBndCondExpansions()[n]
                           ->GetExp(e)
                           ->GetNumPoints(0);
                id1 = m_fields[0]->GetBndCondExpansions()[n]->GetPhys_Offset(e);
                id2 = m_fields[0]->GetTrace()->GetPhys_Offset(
                    m_fields[0]->GetTraceMap()->GetBndCondIDToGlobalTraceID(cnt +
                                                                            e));
                                                                            
                // copy boundary adjusted values into the boundary expansion
                for (i = 0; i < nvariables; ++i)
                {
                    Vmath::Vcopy(npts, &Fwd[i][id2], 1,
                                 &(m_fields[i]
                                       ->GetBndCondExpansions()[n]
                                       ->UpdatePhys())[id1],
                                 1);
                }
            }
        }
        cnt += m_fields[0]->GetBndCondExpansions()[n]->GetExpSize();
    }
}

void MMFStaticDiffOp::EvaluateCovariantVelocity(
    const NekDouble alpha,
    Array<OneD, NekDouble> &velphi,
    Array<OneD, NekDouble> &veltheta,
    Array<OneD, Array<OneD, NekDouble>> &velocity)
{
    int nq = m_fields[0]->GetNpoints();

    NekDouble sin_varphi, cos_varphi, sin_theta, cos_theta;
    // NekDouble x0j, x1j, x2j;

    Array<OneD, NekDouble> x0(nq);
    Array<OneD, NekDouble> x1(nq);
    Array<OneD, NekDouble> x2(nq);

    m_fields[0]->GetCoords(x0, x1, x2);

    // Check
    Array<OneD, NekDouble> tmptheta(nq);
    Array<OneD, NekDouble> tmpphi(nq);
    // theta = a*sin(z/r),  phi = a*tan(y/x);
    for (int j = 0; j < nq; j++)
    {
        CartesianToNewSpherical(x0[j], x1[j], x2[j], sin_varphi, cos_varphi, sin_theta, cos_theta);

        if (m_MMFActivation[j])
        {
            // uhat = m_u0 * (sin_theta * cos(m_alpha) + cos_theta * cos_varphi * sin(m_alpha));
            // vhat = -1.0 * m_u0 * sin_varphi * sin(m_alpha);

            veltheta[j]   = -m_u0 * sin_theta * sin_varphi * sin(alpha);
            velphi[j] = m_u0 * sin_theta * ( cos_theta * cos(alpha) + sin_theta * cos_varphi * sin(alpha) );

            Sph_to_Cart(x0[j], x1[j], x2[j], veltheta[j], velphi[j], velocity[0][j], velocity[1][j], velocity[2][j]);
        }
    }

    Cart_to_MF(velocity, m_sphereMF, tmptheta, tmpphi);

    Vmath::Vsub(nq, veltheta, 1, tmptheta, 1, tmptheta, 1);
    Vmath::Vsub(nq, velphi, 1, tmpphi, 1, tmpphi, 1);

    std::cout << "Error: tmptheta = " << RootMeanSquare(tmptheta, m_MMFActivation) 
    << ", tmpphi = " << RootMeanSquare(tmpphi, m_MMFActivation) << std::endl;
}


void MMFStaticDiffOp::EvaluateAdvectionVelocity(
    Array<OneD, Array<OneD, NekDouble>> &velocity)
{
    int nq = m_fields[0]->GetNpoints();

    NekDouble sin_varphi, cos_varphi, sin_theta, cos_theta;
    NekDouble x0j, x1j, x2j;

    Array<OneD, NekDouble> x0(nq);
    Array<OneD, NekDouble> x1(nq);
    Array<OneD, NekDouble> x2(nq);

    m_fields[0]->GetCoords(x0, x1, x2);

    NekDouble uhat, vhat, uth, uphi;
    // theta = a*sin(z/r),  phi = a*tan(y/x);
    for (int j = 0; j < nq; j++)
    {
        x0j = x0[j];
        x1j = x1[j];
        x2j = x2[j];

        CartesianToNewSpherical(x0j, x1j, x2j, sin_varphi, cos_varphi,
                                sin_theta, cos_theta);

        switch (m_TestType)
        {

            case eTestDivSphericMetId:
            {
                if (m_MMFActivation[j])
                {
                    uth  = 1.0 / sin_theta;
                    uphi = 1.0;

                    velocity[0][j] =
                        -1.0 * uphi * sin_varphi - uth * cos_theta * cos_varphi;
                    velocity[1][j] =
                        uphi * cos_varphi - uth * cos_theta * sin_varphi;
                    velocity[2][j] = uth * sin_theta;
                }
            }
            break;

            case eTestDivRossbyMetId:
            case eTestDivRossbyFlow:
            {
                NekDouble tmp1, tmp2, cos4phi, sin4phi;
                tmp1    = 2 * sin_varphi * cos_varphi;
                tmp2    = 2 * cos_varphi * cos_varphi - 1.0;
                sin4phi = 2 * tmp1 * tmp2;

                tmp1    = 2 * cos_varphi * cos_varphi - 1.0;
                cos4phi = 2 * tmp1 * tmp1 - 1.0;

                uhat = m_W * sin_theta +
                       m_K * sin_theta * sin_theta * sin_theta *
                           (4 * cos_theta * cos_theta - sin_theta * sin_theta) *
                           cos4phi;
                vhat = -4.0 * m_K * sin_theta * sin_theta * sin_theta *
                       cos_theta * sin4phi;

                velocity[0][j] =
                    -1.0 * uhat * sin_varphi - vhat * cos_theta * cos_varphi;
                velocity[1][j] =
                    uhat * cos_varphi - vhat * cos_theta * sin_varphi;
                velocity[2][j] = vhat * sin_theta;
            }
            break;

            case eTestGradPlaneFlow:
            case eTestGradSphericMetID:
            case eTestGradRossbyFlow:
            {
                velocity[0][j] = m_movingframes[m_graddir][j];
                velocity[1][j] = m_movingframes[m_graddir][j + nq];
                velocity[2][j] = m_movingframes[m_graddir][j + 2 * nq];
            }
            break;

            default:
                break;
        }
    }
}

// u^{\theta} = - u_0 * \sin \theta * \sin \phi * \sin \alpha
// u^{\phi} = u_0 * \sin \theta * ( cos \theta * \cos \alpha + \sin \theta * \cos \phi * sin \alpha )
// direction = 0: theta
// direction = 1: phi
void MMFStaticDiffOp::EvaluateExactCovariantDeriv(
    const int direction, 
    Array<OneD, NekDouble> &outtheta, 
    Array<OneD, NekDouble> &outphi)
{
    int nq = GetNpoints();

    NekDouble sin_varphi, cos_varphi, sin_theta, cos_theta;
    NekDouble cos2theta, sin2theta;

    outtheta = Array<OneD, NekDouble>(nq,0.0);
    outphi = Array<OneD, NekDouble>(nq,0.0);

    Array<OneD, NekDouble> x0(nq);
    Array<OneD, NekDouble> x1(nq);
    Array<OneD, NekDouble> x2(nq);

    m_fields[0]->GetCoords(x0, x1, x2);

    NekDouble tmp;
    NekDouble Dthuth, Dthuphi, Dphiuth, Dphiuphi;
    for (int i = 0; i < nq; i++)
    {
        CartesianToNewSpherical(x0[i], x1[i], x2[i], sin_varphi, cos_varphi,
                                sin_theta, cos_theta);

        cos2theta = 2.0 * cos_theta * cos_theta - 1.0;
        sin2theta = 2.0 * sin_theta * cos_theta;

        // D_{\theta} uth
        Dthuth = -m_u0 * cos_theta * sin_varphi * sin(m_alpha);
        Dthuphi = m_u0 * ( cos2theta * cos(m_alpha) + sin2theta * cos_varphi * sin(m_alpha) );

        // D_{\phi} uth
        Dphiuth = -m_u0 * cos_varphi * sin(m_alpha);
        Dphiuphi = -m_u0 * sin_theta * sin_varphi * sin(m_alpha);

        tmp = cos_theta * cos(m_alpha) + sin_theta * cos_varphi * sin(m_alpha);

        switch(direction)
        {
            // Covariant derivative for the theta component
            case 0:
            {
                // \nabla_{\theta} u \cdot \theta
                outtheta[i] = Dthuth ;      

                // \nabla_{\theta} u \cdot \phi
                outphi[i] = Dphiuth - m_u0 * cos_theta * tmp;     

            }
            break;

            // Covariant derivative for the phi component
            case 1:
            {
                // \nabla_{\phi} u \cdot \theta
                outtheta[i]   = Dthuphi ; 

                // \nabla_{\phi} u \cdot \phi
                outphi[i]   = Dphiuphi - m_u0 * cos_theta * sin_varphi * sin(m_alpha);
            }
            break;

            default:
             break;
        }
    }
}

Array<OneD, NekDouble> MMFStaticDiffOp::ComputeExactValue(TestType Ttype)
{
    int nq = GetNpoints();

    NekDouble sin_varphi, cos_varphi, sin_theta, cos_theta;
    NekDouble cos2phi, sin2phi, sin4phi, cos4phi, tmp;

    Array<OneD, NekDouble> outarray;
    if ( (Ttype == eTestGradRossbyFlow) || (Ttype == eTestGradSphericMetID) )
    {
        outarray = Array<OneD, NekDouble>(m_shapedim * nq);
    }

    else
    {
        outarray = Array<OneD, NekDouble>(nq);
    }

    Array<OneD, NekDouble> x0(nq);
    Array<OneD, NekDouble> x1(nq);
    Array<OneD, NekDouble> x2(nq);

    m_fields[0]->GetCoords(x0, x1, x2);

    NekDouble xp, yp, zp;

    // theta = a*sin(z/r),  phi = a*tan(y/x);
    for (int i = 0; i < nq; i++)
    {
        xp = x0[i];
        yp = x1[i];
        zp = x2[i];

        CartesianToNewSpherical(xp, yp, zp, sin_varphi, cos_varphi,
                                sin_theta, cos_theta);

        cos2phi = 2.0 * cos_varphi * cos_varphi - 1.0;
        cos4phi = 2.0 * cos2phi * cos2phi - 1.0;

        sin2phi = 2.0 * sin_varphi * cos_varphi;
        sin4phi = 2.0 * sin2phi * cos2phi;

        switch (Ttype)
        {
            case eTestDivPlane:
            {
                outarray[i] = -1 * sin(x0[i]) + cos(x1[i]);
            }
            break;

            case eTestGradPlaneFlow:
            {
                outarray[i] = 0.0;
            }
            break;

            case eTestGradSphericMetID:
            {
                outarray[i] = 0.0;
                outarray[i+nq] = 0.0;
            }
            break;

            case eTestGradRossbyFlow:
            {
                // dfdphi
                outarray[i] = -4.0 * m_K * sin_theta *
                         (4.0 * cos_theta * cos_theta - sin_theta * sin_theta) *
                         sin4phi;

                tmp =
                    3.0 * cos_theta *
                        (4.0 * cos_theta * cos_theta - sin_theta * sin_theta) -
                    10.0 * sin_theta * sin_theta * cos_theta;

                // dfdth
                outarray[i+nq] = m_W * cos_theta + m_K * sin_theta * sin_theta * tmp * cos4phi;
            }
            break;

            case eTestLaplacianPlaneFlow:
            {
                outarray[i] =
                    -2.0 * m_pi * m_pi * sin(m_pi * x0[i]) * cos(m_pi * x1[i]);
            }
            break;

            // case eTestLaplacianSimpleFlow:
            // {
            //     outarray[i] = -2.0 * m_W * cos_theta +
            //                   30.0 * m_K * sin_theta * sin_theta * sin_theta *
            //                       sin_theta * cos_theta * cos4phi;
            // }
            // break;


            case eTestLaplacianSimpleFlow:
            {
                NekDouble dudth, dudth2, dudphi2;

                // (cos_thehta/sin_theta) df/dtheta
                dudth = 2.0 * m_W * cos_theta * cos_theta * cos2phi;

                // d^2f/dtheta^2
                dudth2 = 2.0 * m_W * (2.0 * cos_theta * cos_theta - 1.0) * cos2phi;

                dudphi2 = -4.0 * m_W * cos2phi;

                outarray[i] = dudphi2 + dudth + dudth2;
            }
            break;


            case eTestLaplacianRossbyFlow:
            {
                NekDouble dudth, dudth2, dudphi2, temp;

                dudth = -1.0 * m_W * cos_theta +
                        + m_K * cos_theta * ( 2 * cos_theta * cos_theta - sin_theta * sin_theta ) * cos2phi;

                temp = 2.0 * cos_theta * cos_theta * cos_theta - 4.0 * sin_theta * sin_theta * cos_theta - 3.0 * sin_theta * sin_theta * cos_theta;

                dudth2 = -1.0 * m_W * cos_theta + m_K * temp * cos2phi;

                dudphi2 = -4.0 * m_K * cos_theta * cos2phi;

                outarray[i] = dudphi2 + dudth + dudth2;

            }
            break;


            case eTestCurlRossbyFlow:
            {
                outarray[i] = -2.0 * m_W * cos_theta +
                              30.0 * m_K * sin_theta * sin_theta * sin_theta *
                                  sin_theta * cos_theta * cos4phi;
            }                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
            break;

            case eTestDivSphericMetId:
            case eTestDivRossbyMetId:
            case eTestDivRossbyFlow:
            case eTestCurlPlaneFlow:
            case eTestCurlSphericMetId:
            case eTestCurlRossbyMetId:
            case eTestLaplacianMetricID:
            {
                outarray[i] = 0.0;
            }
            break;

            default:
                break;
        }
    }

    return outarray;
}

// void MMFStaticDiffOp::ComputeMaxlength()
// {
//     int nq = GetNpoints();

//     // Compute distance and arc length between grid points
//     Array<OneD, NekDouble> x0(nq);
//     Array<OneD, NekDouble> x1(nq);
//     Array<OneD, NekDouble> x2(nq);

//     Array<OneD, NekDouble> temp(nq);
//     Array<OneD, NekDouble> dist(nq);

//     m_fields[0]->GetCoords(x0, x1, x2);

//     Vmath::Sadd(nq, -x0[0], x0, 1, temp, 1);
//     Vmath::Vmul(nq, temp, 1, temp, 1, dist, 1);

//     Vmath::Sadd(nq, -x1[0], x1, 1, temp, 1);
//     Vmath::Vvtvp(nq, temp, 1, temp, 1, dist, 1, dist, 1);

//     Vmath::Sadd(nq, -x2[0], x2, 1, temp, 1);
//     Vmath::Vvtvp(nq, temp, 1, temp, 1, dist, 1, dist, 1);

//     Vmath::Vsqrt(nq, dist, 1, dist, 1);

//     int indx;
//     NekDouble Maxdist;
//     Maxdist = Vmath::Vamax(nq, dist, 1);
//     indx    = Vmath::Iamax(nq, dist, 1);

//     NekDouble rad1, rad2, angle;
//     rad1 = sqrt(x0[0] * x0[0] + x1[0] * x1[0] + x2[0] * x2[0]);
//     rad2 =
//         sqrt(x0[indx] * x0[indx] + x1[indx] * x1[indx] + x2[indx] *
//         x2[indx]);

//     angle = acos((x0[0] * x0[indx] + x1[0] * x1[indx] + x2[0] * x2[indx]) /
//                  rad1 / rad2);
//     std::cout << "Max dist = " << Maxdist
//              << ", approx. arc length = " << angle * 1.0 << std::endl;
// }

// NekDouble MMFStaticDiffOp::ComputeCirculatingArclength(const NekDouble
// zlevel,
//                                                        const NekDouble Rhs)
// {

//     NekDouble Tol = 0.0001, Maxiter = 1000, N = 100, tmp;

//     Array<OneD, NekDouble> xp(N + 1);
//     Array<OneD, NekDouble> yp(N + 1);

//     NekDouble intval;
//     switch (m_surfaceType)
//     {
//         case SolverUtils::eSphere:
//         case SolverUtils::eTRSphere:
//         {
//             intval = sqrt(Rhs - zlevel * zlevel);
//         }
//         break;

//         case SolverUtils::eIrregular:
//         {
//             intval = sqrt(0.5 * (Rhs - zlevel * zlevel * zlevel * zlevel -
//                                  zlevel * zlevel));
//         }
//         break;

//         case SolverUtils::eNonconvex:
//         {
//             tmp = 0.5 *
//                   (Rhs - zlevel * zlevel * zlevel * zlevel - zlevel *
//                   zlevel);
//             intval = sqrt(0.5 * (1.0 + sqrt(1.0 + 4.0 * tmp)));
//         }
//         break;

//         default:
//             break;
//     }

//     switch (m_surfaceType)
//     {
//             // Find the half of all the xp and yp on zlevel ....
//         case SolverUtils::eSphere:
//         case SolverUtils::eTRSphere:
//         case SolverUtils::eIrregular:
//         {
//             NekDouble newy, F=0.0, dF=0.0, y0=0.0;
//             for (int j = 0; j < N + 1; ++j)
//             {
//                 xp[j] = j * 2.0 * intval / N - intval;

//                 y0 = 1.0;
//                 for (int i = 0; i < Maxiter; ++i)
//                 {
//                     switch (m_surfaceType)
//                     {
//                             // Find the half of all the xp and yp on zlevel
//                             // ....
//                         case SolverUtils::eSphere:
//                         case SolverUtils::eTRSphere:
//                         {
//                             F = xp[j] * xp[j] + y0 * y0 + zlevel * zlevel -
//                             Rhs; dF = 2.0 * y0;
//                         }
//                         break;

//                         case SolverUtils::eIrregular:
//                         {
//                             F = 2.0 * xp[j] * xp[j] + y0 * y0 * y0 * y0 +
//                                 y0 * y0 + zlevel * zlevel * zlevel * zlevel +
//                                 zlevel * zlevel - Rhs;
//                             dF = 4.0 * y0 * y0 * y0 + 2.0 * y0;
//                         }
//                         break;

//                         default:
//                             break;
//                     }

//                     newy = y0 - F / dF;

//                     if (fabs(F / dF) < Tol)
//                     {
//                         yp[j] = newy;
//                         break;
//                     }

//                     else
//                     {
//                         y0 = newy;
//                     }

//                     ASSERTL0(i < Maxiter,
//                              "Advection Velocity convergence fails");

//                 } // i-loop
//             }
//         }
//         break;

//         case SolverUtils::eNonconvex:
//         {
//             for (int j = 0; j < N + 1; ++j)
//             {
//                 xp[j] = j * 2.0 * intval / N - intval;
//                 tmp   = 0.5 * Rhs -
//                       0.5 * (zlevel * zlevel * zlevel * zlevel +
//                              zlevel * zlevel) -
//                       (xp[j] * xp[j] * xp[j] * xp[j] - xp[j] * xp[j]);
//                 if (tmp < 0)
//                 {
//                     tmp = -1.0 * tmp;
//                 }
//                 yp[j] = sqrt(tmp);
//             } // j-loop
//         }
//         break;

//         default:
//             break;

//     } // switch-loop

//     NekDouble pi        = 3.14159265358979323846;
//     NekDouble arclength = 0.0;
//     for (int j = 0; j < N; ++j)
//     {
//         arclength =
//             arclength + sqrt((yp[j + 1] - yp[j]) * (yp[j + 1] - yp[j]) +
//                              (xp[j + 1] - xp[j]) * (xp[j + 1] - xp[j])) /
//                             pi;
//     }

//     return arclength;
// }

Array<OneD, NekDouble> MMFStaticDiffOp::SetInitialU()
{
    int j;
    int nq = m_fields[0]->GetNpoints();

    NekDouble sin_varphi, cos_varphi, sin_theta, cos_theta;
    NekDouble x0j, x1j, x2j;

    Array<OneD, NekDouble> x0(nq);
    Array<OneD, NekDouble> x1(nq);
    Array<OneD, NekDouble> x2(nq);

    m_fields[0]->GetCoords(x0, x1, x2);

    Array<OneD, NekDouble> initu(nq, 0.0);
    switch (m_TestType)
    {
        case eTestDivPlane:
        case eTestDivSphericMetId:
        case eTestDivRossbyMetId:
        case eTestDivRossbyFlow:
        {
            initu = Array<OneD, NekDouble>(nq, 1.0);
        }
        break;

        case eTestGradPlaneFlow:
        case eTestGradSphericMetID:
        {
            for (j = 0; j < nq; j++)
            {
                initu[j] = 1.0;
            }
        }
        break;

        case eTestGradRossbyFlow:
        {
            NekDouble cos2phi, cos4phi;

            // theta = a*sin(z/r),  phi = a*tan(y/x);
            for (j = 0; j < nq; j++)
            {
                x0j = x0[j];
                x1j = x1[j];
                x2j = x2[j];

                CartesianToNewSpherical(x0j, x1j, x2j, sin_varphi,
                                        cos_varphi, sin_theta, cos_theta);

                cos2phi = 2 * cos_varphi * cos_varphi - 1.0;
                cos4phi = 2 * cos2phi * cos2phi - 1.0;

                initu[j] =
                    m_W * sin_theta +
                    m_K * sin_theta * sin_theta * sin_theta *
                        (4 * cos_theta * cos_theta - sin_theta * sin_theta) *
                        cos4phi;
            }
        }
        break;

        case eTestLaplacianPlaneFlow:
        {
            for (j = 0; j < nq; j++)
            {
                x0j = x0[j];
                x1j = x1[j];
                x2j = x2[j];

                initu[j] = sin(m_pi * x0j) * cos(m_pi * x1j);
            }
        }
        break;

        case eTestLaplacianMetricID:
        {
            for (j = 0; j < nq; j++)
            {
                initu[j] = 1.0;
            }
        }
        break;


        case eTestLaplacianSimpleFlow:
        {
            NekDouble cos2phi; //, cos4phi;
            for (j = 0; j < nq; j++)
            {
                x0j = x0[j];
                x1j = x1[j];
                x2j = x2[j];

                CartesianToNewSpherical(x0j, x1j, x2j, sin_varphi,
                                        cos_varphi, sin_theta, cos_theta);

                cos2phi = 2 * cos_varphi * cos_varphi - 1.0;

                initu[j] = m_W * sin_theta * sin_theta * cos2phi;
            }
        }
        break;

        case eTestLaplacianRossbyFlow:
        {
            NekDouble cos2phi;
            for (j = 0; j < nq; j++)
            {
                x0j = x0[j];
                x1j = x1[j];
                x2j = x2[j];

                CartesianToNewSpherical(x0j, x1j, x2j, sin_varphi,
                                        cos_varphi, sin_theta, cos_theta);

                cos2phi = 2 * cos_varphi * cos_varphi - 1.0;

                initu[j] = m_W * cos_theta + m_K * sin_theta * sin_theta * cos_theta * cos2phi;

                // initu[j] = m_W * sin_theta + m_K * cos_theta * cos_theta * cos_theta *
                //     (4.0 * sin_theta * sin_theta - cos_theta * cos_theta) *
                //     cos4phi;
            }
        }
        break;

        default:
            break;
    }

    return initu;
} // namespace Nektar

void MMFStaticDiffOp::SetInitialU(Array<OneD, NekDouble> &InitialUth,
                                  Array<OneD, NekDouble> &InitialUphi)
{
    int nq = m_fields[0]->GetNpoints();

    NekDouble sin_varphi, cos_varphi, sin_theta, cos_theta;

    Array<OneD, NekDouble> x0(nq);
    Array<OneD, NekDouble> x1(nq);
    Array<OneD, NekDouble> x2(nq);

    m_fields[0]->GetCoords(x0, x1, x2);

    InitialUth  = Array<OneD, NekDouble>(nq);
    InitialUphi = Array<OneD, NekDouble>(nq);
    switch (m_TestType)
    {
        case eTestCurlPlaneFlow:
        {
            for (int j = 0; j < nq; j++)
            {
                CartesianToNewSpherical(x0[j], x1[j], x2[j], sin_varphi,
                                        cos_varphi, sin_theta, cos_theta);

                InitialUth[j]  = x0[j] * x0[j];
                InitialUphi[j] = x1[j] * x1[j];
            }
        }
        break;

        case eTestCurlSphericMetId:
        {
            for (int j = 0; j < nq; j++)
            {
                CartesianToNewSpherical(x0[j], x1[j], x2[j], sin_varphi,
                                        cos_varphi, sin_theta, cos_theta);

                InitialUth[j]  = 1.0;
                InitialUphi[j] = 1.0 / sin_theta;
            }
        }
        break;

        case eTestCurlRossbyMetId:
        {
            NekDouble sin2phi, cos2phi, cos4phi, sin4phi;
            for (int j = 0; j < nq; j++)
            {
                CartesianToNewSpherical(x0[j], x1[j], x2[j], sin_varphi,
                                        cos_varphi, sin_theta, cos_theta);

                sin2phi = 2 * sin_varphi * cos_varphi;
                cos2phi = 2 * cos_varphi * cos_varphi - 1.0;

                sin4phi = 2 * sin2phi * cos2phi;
                cos4phi = 2 * cos2phi * cos2phi - 1.0;
                                
                InitialUth[j] =
                    m_W * sin_theta +
                    m_K * sin_theta * sin_theta * sin_theta *
                        (4.0 * cos_theta * cos_theta - sin_theta * sin_theta) *
                        cos4phi;

                InitialUphi[j] = -4.0 * m_K * sin_theta * sin_theta * sin_theta *
                                cos_theta * sin4phi;
            }
        }
        break;

        case eTestCurlRossbyFlow:
        {
            NekDouble sin2phi, cos2phi, cos4phi, sin4phi;
            for (int j = 0; j < nq; j++)
            {
                CartesianToNewSpherical(x0[j], x1[j], x2[j], sin_varphi,
                                        cos_varphi, sin_theta, cos_theta);

                sin2phi = 2 * sin_varphi * cos_varphi;
                cos2phi = 2 * cos_varphi * cos_varphi - 1.0;

                sin4phi = 2 * sin2phi * cos2phi;
                cos4phi = 2 * cos2phi * cos2phi - 1.0;

                InitialUphi[j] =
                    -m_W * sin_theta -
                    m_K * sin_theta * sin_theta * sin_theta *
                        (4.0 * cos_theta * cos_theta - sin_theta * sin_theta) *
                        cos4phi;

                InitialUth[j] = -4.0 * m_K * sin_theta * sin_theta * sin_theta *
                                cos_theta * sin4phi;
            }
        }
        break;

        default:
            break;
    }
}

void MMFStaticDiffOp::v_SetInitialConditions(const NekDouble initialtime,
                                             bool dumpInitialConditions,
                                             const int domain)
{
    int nq = m_fields[0]->GetNpoints();

    Array<OneD, NekDouble> initu(nq);
    initu = SetInitialU();

    if (initialtime > 0)
    {
    }

    if (domain > 0)
    {
    }

    Array<OneD, NekDouble> zeros(nq, 0.0);
    m_fields[0]->SetPhys(initu);
    m_fields[1]->SetPhys(zeros);
    m_fields[2]->SetPhys(zeros);

    // forward transform to fill the modal coeffs
    for (int i = 0; i < m_fields.size(); ++i)
    {
        m_fields[i]->SetPhysState(true);
        m_fields[i]->FwdTrans(m_fields[i]->GetPhys(),
                              m_fields[i]->UpdateCoeffs());
    }

    if (dumpInitialConditions)
    {
        // dump initial conditions to file
        std::string outname = m_sessionName + "_initial.chk";
        WriteFld(outname);
    }
}

void MMFStaticDiffOp::AdvectionBellPlane(Array<OneD, NekDouble> &outfield)
{
    int nq = m_fields[0]->GetNpoints();

    NekDouble dist, m_radius_limit;

    Array<OneD, NekDouble> x(nq);
    Array<OneD, NekDouble> y(nq);
    Array<OneD, NekDouble> z(nq);

    m_fields[0]->GetCoords(x, y, z);

    // Sets of parameters
    m_radius_limit = 0.5;

    NekDouble x0j, x1j;
    outfield = Array<OneD, NekDouble>(nq, 0.0);
    for (int j = 0; j < nq; ++j)
    {
        x0j = x[j];
        x1j = y[j];

        dist = sqrt(x0j * x0j + x1j * x1j);

        if (dist < m_radius_limit)
        {
            outfield[j] = 0.5 * (1.0 + cos(m_pi * dist / m_radius_limit));
        }
        else
        {
            outfield[j] = 0.0;
        }
    }
}

void MMFStaticDiffOp::AdvectionBellSphere(Array<OneD, NekDouble> &outfield)
{
    int nq = m_fields[0]->GetNpoints();

    NekDouble dist, radius, cosdiff, sin_theta, cos_theta, sin_varphi,
        cos_varphi;
    NekDouble m_theta_c, m_varphi_c, m_radius_limit, m_c0;

    Array<OneD, NekDouble> x(nq);
    Array<OneD, NekDouble> y(nq);
    Array<OneD, NekDouble> z(nq);

    m_fields[0]->GetCoords(x, y, z);

    // Sets of parameters
    m_theta_c      = 0.0;
    m_varphi_c     = 3.0 * m_pi / 2.0;
    m_radius_limit = 7.0 * m_pi / 64.0;
    m_c0           = 0.0;

    NekDouble x0j, x1j, x2j;
    outfield = Array<OneD, NekDouble>(nq, 0.0);
    for (int j = 0; j < nq; ++j)
    {
        x0j = x[j];
        x1j = y[j];
        x2j = z[j];

        radius = sqrt(x0j * x0j + x1j * x1j + x2j * x2j);

        sin_varphi = x1j / sqrt(x0j * x0j + x1j * x1j);
        cos_varphi = x0j / sqrt(x0j * x0j + x1j * x1j);

        sin_theta = x2j / radius;
        cos_theta = sqrt(x0j * x0j + x1j * x1j) / radius;

        cosdiff = cos_varphi * cos(m_varphi_c) + sin_varphi * sin(m_varphi_c);
        dist    = radius * acos(sin(m_theta_c) * sin_theta +
                             cos(m_theta_c) * cos_theta * cosdiff);

        if (dist < m_radius_limit)
        {
            outfield[j] =
                0.5 * (1.0 + cos(m_pi * dist / m_radius_limit)) + m_c0;
        }
        else
        {
            outfield[j] = m_c0;
        }
    }
}

void MMFStaticDiffOp::Test2Dproblem(const NekDouble time,
                                    Array<OneD, NekDouble> &outfield)
{
    int nq = m_fields[0]->GetNpoints();

    Array<OneD, NekDouble> x0(nq);
    Array<OneD, NekDouble> x1(nq);
    Array<OneD, NekDouble> x2(nq);

    m_fields[0]->GetCoords(x0, x1, x2);

    Array<OneD, NekDouble> u(nq);

    for (int i = 0; i < nq; ++i)
    {
        u[i] = cos(m_pi * (x0[i] - m_advx * time)) *
               cos(m_pi * (x1[i] - m_advy * time));
    }

    outfield = u;
}

void MMFStaticDiffOp::Test3Dproblem(const NekDouble time,
                                    Array<OneD, NekDouble> &outfield)
{
    int nq = m_fields[0]->GetNpoints();

    Array<OneD, NekDouble> x0(nq);
    Array<OneD, NekDouble> x1(nq);
    Array<OneD, NekDouble> x2(nq);

    m_fields[0]->GetCoords(x0, x1, x2);

    Array<OneD, NekDouble> u(nq);

    for (int i = 0; i < nq; ++i)
    {
        u[i] = cos(m_pi * (x0[i] - m_advx * time)) *
               cos(m_pi * (x1[i] - m_advy * time)) *
               cos(m_pi * (x2[i] - m_advz * time));
    }

    outfield = u;
}

// Compute the Laplacian operator in the Cartesian coordinates
void MMFStaticDiffOp::ExplicitDiffusion(
    const int nConvectiveFields,
    const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray)
{
    int i, j, k;
    int nDim      = fields[0]->GetCoordim(0);
    int nPts      = fields[0]->GetTotPoints();
    int nCoeffs   = fields[0]->GetNcoeffs();
    int nTracePts = fields[0]->GetTrace()->GetTotPoints();

    Array<OneD, NekDouble> qcoeffs(nCoeffs);

    Array<OneD, Array<OneD, NekDouble>> fluxvector(nDim);

    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> flux(nDim);
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> qfield(nDim);
    for (j = 0; j < nDim; ++j)
    {
        qfield[j] = Array<OneD, Array<OneD, NekDouble>>(nConvectiveFields);
        flux[j]   = Array<OneD, Array<OneD, NekDouble>>(nConvectiveFields);

        for (i = 0; i < nConvectiveFields; ++i)
        {
            qfield[j][i] = Array<OneD, NekDouble>(nPts, 0.0);
            flux[j][i]   = Array<OneD, NekDouble>(nTracePts, 0.0);
        }
    }

    for (k = 0; k < nDim; ++k)
    {
        fluxvector[k] = Array<OneD, NekDouble>(nPts, 0.0);
    }

    // Get flux = u* (nx, ny, yz)
    Getuflux(fields, inarray, flux);

    for (j = 0; j < nDim; ++j)
    {
        for (i = 0; i < nConvectiveFields; ++i)
        {
            // Compute L2 = \int ( \partial phi / \partial x_j) u d x
            fields[i]->IProductWRTDerivBase(j, inarray[i], qcoeffs);

            // Compute -L2
            Vmath::Neg(nCoeffs, qcoeffs, 1);

            // Compute L = -L2 + \int_{\partial} u* n_x dx
            fields[i]->AddTraceIntegral(flux[j][i], qcoeffs);
            // fields[i]->SetPhysState(false);

            // Compute M^{-1} ( -L2 + \int_{\partial} u* n_x dx )
            fields[i]->MultiplyByElmtInvMass(qcoeffs, qcoeffs);
            fields[i]->BwdTrans(qcoeffs, qfield[j][i]);
        }
    }

    // Compute u from q_{\eta} and q_{\xi} to obtain numerical fluxes
    // NumFluxforVector(fields, inarray, qfield, flux[0]);
    Getqflux(fields, inarray, qfield, flux[0]);

    Array<OneD, NekDouble> tmp(nCoeffs);
    Array<OneD, NekDouble> DivSum;
    for (i = 0; i < nConvectiveFields; ++i)
    {
        // Compute tmp = \int \nabla \varphi \cdot \vec{q}
        DivSum = Array<OneD, NekDouble>(nCoeffs, 0.0);
        for (int j = 0; j < nDim; ++j)
        {
            fields[i]->IProductWRTDerivBase(j, qfield[j][i], tmp);
            Vmath::Vadd(nCoeffs, tmp, 1, DivSum, 1, DivSum, 1);
        }

        // Evaulate  <\phi, \hat{F}\cdot n> - outarray[i]
        Vmath::Neg(nCoeffs, DivSum, 1);
        fields[i]->AddTraceIntegral(flux[0][i], DivSum);
        // fields[i]->SetPhysState(false);

        fields[i]->MultiplyByElmtInvMass(DivSum, DivSum);
        fields[i]->BwdTrans(DivSum, outarray[i]);
    }
}

// Compute u* \vec{n} = (u* nx, u* ny, u* nz)
void MMFStaticDiffOp::Getuflux(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
    const Array<OneD, Array<OneD, NekDouble>> &ufield,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &uflux)
{
    int j, k, var;
    int nTracePts  = fields[0]->GetTrace()->GetTotPoints();
    int nvariables = fields.size();
    int nDim       = fields[0]->GetCoordim(0);

    NekDouble uAver, uJump;
    Array<OneD, NekDouble> uFwd(nTracePts), uBwd(nTracePts), uftmp(nTracePts);
    for (var = 0; var < nvariables; ++var)
    {
        fields[var]->GetFwdBwdTracePhys(ufield[var], uFwd, uBwd);

        // Boundary Treatment for Dirichlet and Neumann
        DiffusionBoundaryConditions(euflux, var, uFwd, uBwd);

        // Compute u* = {{ u }} + \beta \vec{n}^+ \cdot [[ u ]]
        for (k = 0; k < nTracePts; ++k)
        {
            uAver    = 0.5 * (uFwd[k] + uBwd[k]);
            uJump    = 0.5 * (uFwd[k] - uBwd[k]);
            uftmp[k] = uAver + uJump;
        }

        // Return u* * (nx, ny, nz)
        for (j = 0; j < nDim; ++j)
        {
            Vmath::Vmul(nTracePts, m_traceNormals[j], 1, uftmp, 1,
                        uflux[j][var], 1);
        }
    }
}

// Compute q* \cdot \vec{n}^+
void MMFStaticDiffOp::Getqflux(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
    const Array<OneD, Array<OneD, NekDouble>> &ufield,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &qfield,
    Array<OneD, Array<OneD, NekDouble>> &qflux)
{
    int j, k, var;
    int nTracePts  = fields[0]->GetTrace()->GetTotPoints();
    int nvariables = fields.size();
    int nDim       = fields[0]->GetCoordim(0);

    Array<OneD, NekDouble> qftmp(nTracePts, 0.0);

    NekDouble qAver, qJump;
    Array<OneD, NekDouble> uFwd(nTracePts);
    Array<OneD, NekDouble> uBwd(nTracePts);

    Array<OneD, Array<OneD, NekDouble>> qFwd(nDim);
    Array<OneD, Array<OneD, NekDouble>> qBwd(nDim);

    for (var = 0; var < nvariables; ++var)
    {
        // Compute q^* \cdot n^+
        fields[0]->GetFwdBwdTracePhys(ufield[var], uFwd, uBwd);

        for (j = 0; j < nDim; ++j)
        {
            qFwd[j] = Array<OneD, NekDouble>(nTracePts, 0.0);
            qBwd[j] = Array<OneD, NekDouble>(nTracePts, 0.0);
            fields[0]->GetFwdBwdTracePhys(qfield[j][var], qFwd[j], qBwd[j]);

            // Boundary Treatment for Dirichlet and Neumann
            DiffusionBoundaryConditions(eqflux, var, qFwd[j], qBwd[j], j);
        }

        for (k = 0; k < nTracePts; ++k)
        {
            qAver = 0.0;
            qJump = 0.0;
            for (int j = 0; j < nDim; ++j)
            {
                // Compute {{ \vec{q} }} \cdot \vec{n}^+
                qAver = qAver +
                        0.5 * (qFwd[j][k] + qBwd[j][k]) * m_traceNormals[j][k];

                // Compute [[ \vec{q} ]]
                qJump = qJump +
                        0.5 * (qFwd[j][k] - qBwd[j][k]) * m_traceNormals[j][k];
            }

            // Compute [[ u ]] \cdot \vec{n}^+
            // uJump = uFwd[k] - uBwd[k];

            qflux[var][k] = qAver - qJump;

            // qflux[var][k] =
            //  qAver - m_Diffbeta * qJump - (m_Diffeta / m_Diffhe) * uJump;
        }
    }
}

// Boundary conditions for uflux
void MMFStaticDiffOp::DiffusionBoundaryConditions(
    FluxType ftype, const int var, const Array<OneD, const NekDouble> &Fwd,
    Array<OneD, NekDouble> &Bwd, int direction)
{
    int i, e, id1, id2;

    // Number of boundary regions
    int nBndEdgePts, nBndEdges;
    int cnt         = 0;
    int nBndRegions = m_fields[var]->GetBndCondExpansions().size();

    for (i = 0; i < nBndRegions; ++i)
    {
        // Number of boundary expansion related to that region
        nBndEdges = m_fields[var]->GetBndCondExpansions()[i]->GetExpSize();

        // Weakly impose boundary conditions by modifying flux values
        for (e = 0; e < nBndEdges; ++e)
        {
            // nBndEdgePts = m_fields[var]
            //                   ->GetBndCondExpansions()[i]
            //                   ->GetExp(e)
            //                   ->GetTotPoints();

            // id1 = m_fields[var]->GetBndCondExpansions()[i]->GetPhys_Offset(e);
            // id2 = m_fields[0]->GetTrace()->GetPhys_Offset(
            //     m_fields[0]->GetTraceMap()->GetBndCondTraceToGlobalTraceMap(
            //         cnt++));

            nBndEdgePts = m_fields[var]
                       ->GetBndCondExpansions()[i]
                       ->GetExp(e)
                       ->GetNumPoints(0);
            id1 = m_fields[var]->GetBndCondExpansions()[i]->GetPhys_Offset(e);
            id2 = m_fields[var]->GetTrace()->GetPhys_Offset(
                m_fields[var]->GetTraceMap()->GetBndCondIDToGlobalTraceID(cnt++));

            if (m_fields[var]
                    ->GetBndConditions()[i]
                    ->GetBoundaryConditionType() == SpatialDomains::eDirichlet)
            {
                if (ftype == euflux)
                {
                    // For Dirichlet boundary condition: uBwd = g_D
                    Vmath::Vcopy(nBndEdgePts,
                                 &(m_fields[var]
                                       ->GetBndCondExpansions()[i]
                                       ->GetPhys())[id1],
                                 1, &Bwd[id2], 1);
                }

                else if (ftype == eqflux)
                {
                    // For Dirichlet boundary condition: qBwd = qFwd
                    Vmath::Vcopy(nBndEdgePts, &Fwd[id2], 1, &Bwd[id2], 1);
                }
            }

            else if ((m_fields[var]->GetBndConditions()[i])
                         ->GetBoundaryConditionType() ==
                     SpatialDomains::eNeumann)
            {
                if (ftype == euflux)
                {
                    // For Neumann boundary condition: uBwd = u+
                    Vmath::Vcopy(nBndEdgePts, &Fwd[id2], 1, &Bwd[id2], 1);
                }
                else if (ftype == eqflux)
                {
                    // For Neumann boundary condition: qBwd = g_N \vec{n}
                    Vmath::Vmul(nBndEdgePts, &m_traceNormals[direction][id2], 1,
                                &(m_fields[var]
                                      ->GetBndCondExpansions()[i]
                                      ->GetPhys())[id1],
                                1, &Bwd[id2], 1);
                }
            }
        }
    }
}

void MMFStaticDiffOp::v_EvaluateExactSolution(unsigned int field,
                                              Array<OneD, NekDouble> &outfield,
                                              const NekDouble time)
{
    int nq = m_fields[0]->GetNpoints();

    if (time > 0.0)
    {
    }

    switch (m_TestType)
    {

        case eTestDivPlane:
        case eTestDivSphericMetId:
        case eTestDivRossbyMetId:
        case eTestDivRossbyFlow:
        case eTestCurlPlaneFlow:
        case eTestCurlSphericMetId:
        case eTestCurlRossbyMetId:
        case eTestCurlRossbyFlow:
        case eTestLaplacianPlaneFlow:
        case eTestLaplacianMetricID:
        case eTestLaplacianSimpleFlow:
        case eTestLaplacianRossbyFlow:
        {
            if (field == 0)
            {
                outfield = Array<OneD, NekDouble>(nq, 1.0);
            }

            else
            {
                outfield = Array<OneD, NekDouble>(nq, 0.0);
            }
        }
        break;

        default:
            break;
    }
}

// void MMFStaticDiffOp::GetFluxVector(
//     const int i, const int j,
//     const Array<OneD, Array<OneD, NekDouble>> &physfield,
//     Array<OneD, Array<OneD, NekDouble>> &derivatives,
//     Array<OneD, Array<OneD, NekDouble>> &flux)
// {
//     for (int k = 0; k < flux.size(); ++k)
//     {
//         Vmath::Zero(GetNpoints(), flux[k], 1);
//     }
//     Vmath::Vcopy(GetNpoints(), physfield[i], 1, flux[j], 1);
// }

void MMFStaticDiffOp::ComputeVarCoeff(
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    StdRegions::VarCoeffMap &varcoeff)
{
    int nq = GetTotPoints();

    StdRegions::VarCoeffType MMFCoeffs[15] = {
        StdRegions::eVarCoeffMF1x,   StdRegions::eVarCoeffMF1y,
        StdRegions::eVarCoeffMF1z,   StdRegions::eVarCoeffMF1Div,
        StdRegions::eVarCoeffMF1Mag, StdRegions::eVarCoeffMF2x,
        StdRegions::eVarCoeffMF2y,   StdRegions::eVarCoeffMF2z,
        StdRegions::eVarCoeffMF2Div, StdRegions::eVarCoeffMF2Mag,
        StdRegions::eVarCoeffMF3x,   StdRegions::eVarCoeffMF3y,
        StdRegions::eVarCoeffMF3z,   StdRegions::eVarCoeffMF3Div,
        StdRegions::eVarCoeffMF3Mag};

    int indx;
    Array<OneD, NekDouble> tmp(nq);
    for (int k = 0; k < m_mfdim; ++k)
    {
        // For Moving Frames
        indx = 5 * k;

        for (int j = 0; j < m_spacedim; ++j)
        {
            varcoeff[MMFCoeffs[indx + j]] = Array<OneD, NekDouble>(nq, 0.0);
            Vmath::Vcopy(nq, &movingframes[k][j * nq], 1,
                         &varcoeff[MMFCoeffs[indx + j]][0], 1);
        }

        // m_DivMF
        varcoeff[MMFCoeffs[indx + 3]] = Array<OneD, NekDouble>(nq, 0.0);

        Array<OneD, Array<OneD, NekDouble>> DivMF;
        ComputeDivMF(m_DerivType, movingframes, DivMF);

        Vmath::Vcopy(nq, &DivMF[k][0], 1, &varcoeff[MMFCoeffs[indx + 3]][0], 1);

        // \| e^k \|
        varcoeff[MMFCoeffs[indx + 4]] = Array<OneD, NekDouble>(nq, 0.0);
        tmp                           = Array<OneD, NekDouble>(nq, 0.0);
        for (int i = 0; i < m_spacedim; ++i)
        {
            Vmath::Vvtvp(nq, &movingframes[k][i * nq], 1,
                         &movingframes[k][i * nq], 1, &tmp[0], 1, &tmp[0], 1);
        }

        Vmath::Vcopy(nq, &tmp[0], 1, &varcoeff[MMFCoeffs[indx + 4]][0], 1);
    }
}

void MMFStaticDiffOp::v_GenerateSummary(SolverUtils::SummaryList &s)
{
    MMFSystem::v_GenerateSummary(s);
    SolverUtils::AddSummaryItem(s, "TestType", TestTypeMap[m_TestType]);
}

void MMFStaticDiffOp::Wait_On_Enter()
{
    std::string dummy;
    std::cout << "Enter to continue..." << std::endl;
    std::getline(std::cin, dummy);
}
} // namespace Nektar
