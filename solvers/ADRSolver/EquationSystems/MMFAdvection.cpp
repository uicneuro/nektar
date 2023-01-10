/////////////////////////////////////////////////////////////////////////////
//
// File MMFAdvection.cpp
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
// Permission is hereby granted, free of charge, to any person obtaining ag
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
// Description: MMF solve routines
//
///////////////////////////////////////////////////////////////////////////////

#include <iomanip>
#include <iostream>

#include <boost/algorithm/string.hpp>
#include <boost/core/ignore_unused.hpp>

#include <ADRSolver/EquationSystems/MMFAdvection.h>
#include <LibUtilities/BasicUtils/Timer.h>
#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>

namespace Nektar
{
std::string MMFAdvection::className =
    SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
        "MMFAdvection", MMFAdvection::create, "MMFAdvection equation.");

MMFAdvection::MMFAdvection(const LibUtilities::SessionReaderSharedPtr &pSession,
                           const SpatialDomains::MeshGraphSharedPtr &pGraph)
    : UnsteadySystem(pSession, pGraph), MMFSystem(pSession, pGraph),
      AdvectionSystem(pSession, pGraph)
{
    m_planeNumber = 0;
}

/**
 * @brief Initialisation object for the unsteady linear advection equation.
 */
void MMFAdvection::v_InitObject(bool DeclareFields)
{
    // Call to the initialisation object
    UnsteadySystem::v_InitObject(DeclareFields);

    int nq       = m_fields[0]->GetNpoints();
    int shapedim = m_fields[0]->GetShapeDimension();
    Array<OneD, Array<OneD, NekDouble>> AniStrength(shapedim);
    for (int j = 0; j < shapedim; ++j)
    {
        AniStrength[j] = Array<OneD, NekDouble>(nq, 1.0);
    }

    MMFSystem::MMFInitObject(AniStrength);

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

    m_session->LoadParameter("Divergence Restore", m_DivergenceRestore, 0);
    m_session->LoadParameter("Angular Frequency", m_waveFreq, m_pi);
    m_session->LoadParameter("Rotational Angle", m_RotAngle, 0.0);

    if (m_DivergenceRestore > 1)
    {
        mf_LOCSPH = Array<OneD, Array<OneD, NekDouble>>(m_mfdim);
        GetLOCALMovingframes(mf_LOCSPH);
        ComputeAxisAlignedLOCALMovingframes(m_sphereMF, mf_LOCSPH);
    }

    m_session->LoadParameter("theta_c", m_theta_c, 0.0);
    m_session->LoadParameter("varphi_c", m_varphi_c, 3.0 * m_pi / 2.0);
    m_session->LoadParameter("radius_limit", m_radius_limit, 7.0 * m_pi / 64.0);

    // Read the advection velocities from session file
    m_session->LoadParameter("advx", m_advx, 1.0);
    m_session->LoadParameter("advy", m_advy, 1.0);
    m_session->LoadParameter("advz", m_advz, 1.0);

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
        m_velocity[k] = Array<OneD, NekDouble>(nq);
    }

    switch (m_surfaceType)
    {
        case SolverUtils::eSphere:
        case SolverUtils::eTRSphere:
        case SolverUtils::eIrregular:
        case SolverUtils::eNonconvex:
        {
            // true = project velocity onto the tangent plane
            EvaluateAdvectionVelocity(m_movingframes, m_velocity);
        }
        break;

        case SolverUtils::ePlane:
        case SolverUtils::eCube:
        {
            GetFunction("AdvectionVelocity")->Evaluate(vel, m_velocity);
        }
        break;

        default:
            break;
    }

    // New Scheme
    ComputeSphericalVector(m_SphericalVector);

    CheckMeshErr(m_movingframes, m_velocity);

    // Compute m_traceVn = n \cdot v
    m_traceVn = GetNormalVelocity(m_velocity);

    // Compute vel \cdot MF
    ComputevelodotMF(m_velocity, m_movingframes);

    // Reflect it into m_ncdotMFFwd and Bwd
    ComputencdotMF(m_movingframes, m_ncdotMFFwd, m_ncdotMFBwd);

    // If explicit it computes RHS and PROJECTION for the time integration
    if (m_explicitAdvection)
    {
        m_ode.DefineOdeRhs(&MMFAdvection::DoOdeRhs, this);
        m_ode.DefineProjection(&MMFAdvection::DoOdeProjection, this);
    }
    // Otherwise it gives an error (no implicit integration)
    else
    {
        ASSERTL0(false, "Implicit unsteady Advection not set up.");
    }
}

/**
 * @brief Unsteady linear advection equation destructor.
 */
MMFAdvection::~MMFAdvection()
{
}

void MMFAdvection::v_DoSolve()
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

    int Ntot, indx;
    // Perform integration in time.
    Ntot = m_steps / m_checksteps + 1;

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
            indx = (step + 1) / m_checksteps - 1;

            dMass[indx] =
                abs((m_fields[0]->PhysIntegral(fields[0]) - m_Mass0) / m_Mass0);
            std::cout << ",   dMass: " << std::setw(8) << std::left
                      << dMass[indx] << std::endl
                      << std::endl;

            NekDouble L2Err, LinfErr;

            Array<OneD, NekDouble> exactsoln(nq);
            v_EvaluateExactSolution(0, exactsoln, 0.0);

            L2Err   = v_L2Error(0, exactsoln, 0);
            LinfErr = v_LinfError(0, exactsoln);

            std::cout << "L2error = " << L2Err << ", LinfErr = " << LinfErr
                      << std::endl;

            cpuTime = 0.0;
        }

        // Transform data into coefficient space
        for (i = 0; i < nvariables; ++i)
        {
            // m_fields[m_intVariables[i]]->SetPhys(fields[i]);
            // m_fields[m_intVariables[i]]->FwdTrans_IterPerExp(
            //     fields[i], m_fields[m_intVariables[i]]->UpdateCoeffs());
            // m_fields[m_intVariables[i]]->SetPhysState(false);
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
    }

    std::cout << "dMass =  ";
    for (int i = 0; i < Ntot; i++)
    {
        std::cout << dMass[i] << " , ";
    }
    std::cout << std::endl << std::endl;

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

void MMFAdvection::Checkpoint_Err(
    const int n, const NekDouble time,
    const Array<OneD, const Array<OneD, NekDouble>> &fieldphys)
{
    int nvar    = m_fields.size();
    int nq      = GetTotPoints();
    int ncoeffs = GetNcoeffs();

    std::string outname =
        m_sessionName + boost::lexical_cast<std::string>(n) + "Err.chk";

    std::vector<Array<OneD, NekDouble>> fieldcoeffs(nvar);
    for (int i = 0; i < nvar; ++i)
    {
        fieldcoeffs[i] = Array<OneD, NekDouble>(ncoeffs);
    }

    std::vector<std::string> variables(nvar);
    variables[0] = "Linf";

    Array<OneD, NekDouble> exactsoln(m_fields[0]->GetNpoints());
    v_EvaluateExactSolution(0, exactsoln, time);

    Vmath::Vsub(nq, exactsoln, 1, fieldphys[0], 1, exactsoln, 1);
    for (int i = 0; i < nq; ++i)
    {
        exactsoln[i] = fabs(exactsoln[i]);
    }

    m_fields[0]->FwdTrans(exactsoln, fieldcoeffs[0]);
    WriteFld(outname, m_fields[0], fieldcoeffs, variables);
}

/**
 * @brief Compute the right-hand side for the linear advection equation.
 *
 * @param inarray    Given fields.
 * @param outarray   Calculated solution.
 * @param time       Time.
 */
void MMFAdvection::DoOdeRhs(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time)
{
    int nvariables = inarray.size();
    int ncoeffs    = GetNcoeffs();
    int nq         = GetNpoints();

    if (time > 0)
    {
    }

    switch (m_projectionType)
    {
        case MultiRegions::eDiscontinuous:
        {
            Array<OneD, Array<OneD, NekDouble>> WeakAdv(nvariables);

            WeakAdv[0] = Array<OneD, NekDouble>(ncoeffs * nvariables);
            for (int i = 1; i < nvariables; ++i)
            {
                WeakAdv[i] = WeakAdv[i - 1] + ncoeffs;
            }

            // Compute \nabla \cdot \vel u according to MMF scheme
            WeakDGDirectionalAdvection(inarray, WeakAdv);

            for (int i = 0; i < nvariables; ++i)
            {
                if (m_DivergenceRestore > 0)
                {
                    Array<OneD, NekDouble> tmpc(ncoeffs);
                    Array<OneD, NekDouble> tmp(nq);

                    Array<OneD, NekDouble> velvector(m_spacedim * nq);
                    for (int k = 0; k < m_spacedim; ++k)
                    {
                        Vmath::Vmul(nq, &inarray[0][0], 1, &m_velocity[k][0], 1,
                                    &velvector[k * nq], 1);
                    }

                    switch (m_DivergenceRestore)
                    {
                        case 1:
                        {
                            tmp = ComputeSpuriousDivergence(
                                m_movingframes, m_SphericalVector[2],
                                velvector);
                            Vmath::Neg(nq, tmp, 1);
                        }
                        break;

                        case 2:
                        {
                            tmp = ComputeSpuriousDivergence(
                                mf_LOCSPH, m_movingframes[2], velvector);
                        }
                        break;

                        case 3:
                        {
                            tmp = ComputeSpuriousDivergence(
                                mf_LOCSPH, m_SphericalVector[2], velvector);
                        }
                        break;

                        default:
                            break;
                    }

                    m_fields[i]->IProductWRTBase(tmp, tmpc);
                    Vmath::Vadd(ncoeffs, tmpc, 1, WeakAdv[i], 1, WeakAdv[i], 1);
                }

                m_fields[i]->MultiplyByElmtInvMass(WeakAdv[i], WeakAdv[i]);
                m_fields[i]->BwdTrans(WeakAdv[i], outarray[i]);

                Vmath::Neg(nq, outarray[i], 1);
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
void MMFAdvection::DoOdeProjection(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time)
{
    // Counter variable
    int i;
    int nQuadraturePts = GetNpoints();
    int nVariables     = inarray.size();

    // Set the boundary conditions
    SetBoundaryConditions(time);

    // Switch on the projection type (Discontinuous or Continuous)
    for (i = 0; i < nVariables; ++i)
    {
        Vmath::Vcopy(nQuadraturePts, inarray[i], 1, outarray[i], 1);
    }
}

/**
 * @brief Return the flux vector for the linear advection equation.
 *
 * @param i           Component of the flux vector to calculate.
 * @param physfield   Fields.
 * @param flux        Resulting flux.
 */
void MMFAdvection::GetFluxVector(
    const Array<OneD, Array<OneD, NekDouble>> &physfield,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &flux)
{
    ASSERTL1(flux[0].size() == m_velocity.size(),
             "Dimension of flux array and velocity array do not match");

    int i, j;
    int nq = physfield[0].size();

    for (i = 0; i < flux.size(); ++i)
    {
        for (j = 0; j < flux[0].size(); ++j)
        {
            Vmath::Vmul(nq, physfield[i], 1, m_velocity[j], 1, flux[i][j], 1);
        }
    }
}

void MMFAdvection::WeakDGDirectionalAdvection(
    const Array<OneD, const Array<OneD, NekDouble>> &InField,
    Array<OneD, Array<OneD, NekDouble>> &OutField)
{
    int i, j;
    int ncoeffs         = GetNcoeffs();
    int nTracePointsTot = GetTraceNpoints();
    int nvariables      = m_fields.size();

    Array<OneD, Array<OneD, NekDouble>> physfield(nvariables);

    // Get the variables in physical space
    // already in physical space
    for (i = 0; i < nvariables; ++i)
    {
        physfield[i] = InField[i];
    }

    Array<OneD, Array<OneD, NekDouble>> WeakDeriv(m_shapedim);
    for (i = 0; i < nvariables; ++i)
    {
        for (j = 0; j < m_shapedim; ++j)
        {
            WeakDeriv[j] = Array<OneD, NekDouble>(ncoeffs, 0.0);

            // Directional derivation with respect to the j'th moving frame
            // tmp[j] = \nabla \physfield[i] \cdot \mathbf{e}^j
            // Implemented at TriExp::v_IProductWRTDirectionalDerivBase_SumFa
            m_fields[0]->IProductWRTDirectionalDerivBase(
                m_movingframes[j], physfield[i], WeakDeriv[j]);
        }

        // Get the numerical flux and add to the modal coeffs
        // if the NumericalFluxs function already includes the
        // normal in the output

        Array<OneD, NekDouble> flux(nTracePointsTot, 0.0);
        Array<OneD, NekDouble> Fwd(nTracePointsTot);
        Array<OneD, NekDouble> Bwd(nTracePointsTot);

        // Evaluate numerical flux in physical space which may in
        // general couple all component of vectors
        m_fields[i]->GetFwdBwdTracePhys(physfield[i], Fwd, Bwd);

        // evaulate upwinded m_fields[i]:
        // m_traceVn = V \cdot n. flux = u+ or u-
        m_fields[i]->GetTrace()->Upwind(m_traceVn, Fwd, Bwd, flux);

        OutField[i] = Array<OneD, NekDouble>(ncoeffs, 0.0);
        for (j = 0; j < m_shapedim; ++j)
        {
            // calculate numflux = (n \cdot MF)*flux
            Vmath::Vmul(nTracePointsTot, &flux[0], 1, &m_ncdotMFFwd[j][0], 1,
                        &Fwd[0], 1);
            Vmath::Vmul(nTracePointsTot, &flux[0], 1, &m_ncdotMFBwd[j][0], 1,
                        &Bwd[0], 1);
            Vmath::Neg(ncoeffs, WeakDeriv[j], 1);

            // FwdBwdtegral because generallize (N \cdot MF)_{FWD} \neq -(N
            // \cdot MF)_{BWD}
            m_fields[i]->AddFwdBwdTraceIntegral(Fwd, Bwd, WeakDeriv[j]);
            m_fields[i]->SetPhysState(false);

            Vmath::Vadd(ncoeffs, &WeakDeriv[j][0], 1, &OutField[i][0], 1,
                        &OutField[i][0], 1);
        }
    }
}

void MMFAdvection::EvaluateAdvectionVelocity(
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    Array<OneD, Array<OneD, NekDouble>> &velocity)
{
    int nq = m_fields[0]->GetNpoints();

    Array<OneD, NekDouble> x0(nq);
    Array<OneD, NekDouble> x1(nq);
    Array<OneD, NekDouble> x2(nq);

    m_fields[0]->GetCoords(x0, x1, x2);

    // theta = a*sin(z/r),  phi = a*tan(y/x);
    NekDouble uhat, vhat, sin_varphi, cos_varphi, sin_theta, cos_theta;
    NekDouble x0j, x1j, x2j;
    for (int j = 0; j < nq; j++)
    {
        x0j = x0[j];
        x1j = x1[j];
        x2j = x2[j];

        CartesianToNewSpherical(x0j, x1j, x2j, sin_varphi, cos_varphi,
                                sin_theta, cos_theta);

        uhat = m_waveFreq * (sin_theta * cos(m_RotAngle) +
                             cos_theta * cos_varphi * sin(m_RotAngle));
        vhat = -1.0 * m_waveFreq * sin_varphi * sin(m_RotAngle);

        velocity[0][j] =
            -1.0 * uhat * sin_varphi - vhat * cos_theta * cos_varphi;
        velocity[1][j] = uhat * cos_varphi - vhat * cos_theta * sin_varphi;
        velocity[2][j] = vhat * sin_theta;
    }

    // Project the veloicty on the tangent plane
    ProjectionOntoMovingFrames(movingframes, velocity, velocity);
}

NekDouble MMFAdvection::ComputeCirculatingArclength(const NekDouble zlevel,
                                                    const NekDouble Rhs)
{

    NekDouble Tol = 0.0001, Maxiter = 1000, N = 100;
    NekDouble newy, F = 0.0, dF = 0.0, y0, tmp;

    Array<OneD, NekDouble> xp(N + 1);
    Array<OneD, NekDouble> yp(N + 1);

    NekDouble intval = 0.0;
    switch (m_surfaceType)
    {
        case SolverUtils::eSphere:
        case SolverUtils::eTRSphere:
        {
            intval = sqrt(Rhs - zlevel * zlevel);
        }
        break;

            //  2*x^2 + y^4 + y + z^4 + z^2 = 1
        case SolverUtils::eIrregular:
        {
            intval = sqrt(0.5 * (Rhs - zlevel * zlevel * zlevel * zlevel -
                                 zlevel * zlevel));
        }
        break;

            //  2 x^2 + 2 (y^4 - y^2 ) + z^4 + z^2 = 2
        case SolverUtils::eNonconvex:
        {
            tmp = 0.5 *
                  (Rhs - zlevel * zlevel * zlevel * zlevel - zlevel * zlevel);
            intval = sqrt(0.5 * (1.0 + sqrt(1.0 + 4.0 * tmp)));
        }
        break;

        default:
            break;
    }

    switch (m_surfaceType)
    {
            // Find the half of all the xp and yp on zlevel ....
        case SolverUtils::eSphere:
        case SolverUtils::eTRSphere:
        case SolverUtils::eIrregular:
        {
            for (int j = 0; j < N + 1; ++j)
            {
                xp[j] = j * 2.0 * intval / N - intval;

                y0 = 1.0;
                for (int i = 0; i < Maxiter; ++i)
                {
                    switch (m_surfaceType)
                    {
                            // Find the half of all the xp and yp on zlevel ....
                        case SolverUtils::eSphere:
                        case SolverUtils::eTRSphere:
                        {
                            F = xp[j] * xp[j] + y0 * y0 + zlevel * zlevel - Rhs;
                            dF = 2.0 * y0;
                        }
                        break;

                        case SolverUtils::eIrregular:
                        {
                            F = 2.0 * xp[j] * xp[j] + y0 * y0 * y0 * y0 +
                                y0 * y0 + zlevel * zlevel * zlevel * zlevel +
                                zlevel * zlevel - Rhs;
                            dF = 4.0 * y0 * y0 * y0 + 2.0 * y0;
                        }
                        break;

                        default:
                            break;
                    }

                    newy = y0 - F / dF;

                    if (fabs(F / dF) < Tol)
                    {
                        yp[j] = newy;
                        break;
                    }

                    else
                    {
                        y0 = newy;
                    }

                    ASSERTL0(i < Maxiter,
                             "Advection Velocity convergence fails");

                } // i-loop
            }
        }
        break;

        case SolverUtils::eNonconvex:
        {
            for (int j = 0; j < N + 1; ++j)
            {
                xp[j] = j * 2.0 * intval / N - intval;
                tmp   = 0.5 * Rhs -
                      0.5 * (zlevel * zlevel * zlevel * zlevel +
                             zlevel * zlevel) -
                      (xp[j] * xp[j] * xp[j] * xp[j] - xp[j] * xp[j]);
                if (tmp < 0)
                {
                    tmp = -1.0 * tmp;
                }
                yp[j] = sqrt(tmp);
            } // j-loop
        }
        break;

        default:
            break;

    } // switch-loop

    NekDouble pi        = 3.14159265358979323846;
    NekDouble arclength = 0.0;
    for (int j = 0; j < N; ++j)
    {
        arclength =
            arclength + sqrt((yp[j + 1] - yp[j]) * (yp[j + 1] - yp[j]) +
                             (xp[j + 1] - xp[j]) * (xp[j + 1] - xp[j])) /
                            pi;
    }

    return arclength;
}

void MMFAdvection::v_SetInitialConditions(const NekDouble initialtime,
                                          bool dumpInitialConditions,
                                          const int domain)
{
    int nq = m_fields[0]->GetNpoints();

    Array<OneD, NekDouble> u(nq);

    if (domain == 0)
    {
    }

    switch (m_TestType)
    {
        case eAdvectionBell:
        {
            AdvectionBellSphere(u);

            m_Mass0 = m_fields[0]->PhysIntegral(u);
            m_fields[0]->SetPhys(u);

            // forward transform to fill the modal coeffs
            for (int i = 0; i < m_fields.size(); ++i)
            {
                m_fields[i]->SetPhysState(true);
                m_fields[i]->FwdTrans(m_fields[i]->GetPhys(),
                                      m_fields[i]->UpdateCoeffs());
            }
        }
        break;

        case eTestPlane:
        {
            Test2Dproblem(initialtime, u);
            m_fields[0]->SetPhys(u);

            m_Mass0 = m_fields[0]->PhysIntegral(u);

            // forward transform to fill the modal coeffs
            for (int i = 0; i < m_fields.size(); ++i)
            {
                m_fields[i]->SetPhysState(true);
                m_fields[i]->FwdTrans(m_fields[i]->GetPhys(),
                                      m_fields[i]->UpdateCoeffs());
            }
        }
        break;

        case eTestPlaneMassConsv:
        {
            AdvectionBellPlane(u);
            m_fields[0]->SetPhys(u);

            m_Mass0 = m_fields[0]->PhysIntegral(u);
            std::cout << "m_Mass0 = " << m_Mass0 << std::endl;

            // forward transform to fill the modal coeffs
            for (int i = 0; i < m_fields.size(); ++i)
            {
                m_fields[i]->SetPhysState(true);
                m_fields[i]->FwdTrans(m_fields[i]->GetPhys(),
                                      m_fields[i]->UpdateCoeffs());
            }
        }
        break;

        case eTestCube:
        {
            Test3Dproblem(initialtime, u);
            m_fields[0]->SetPhys(u);

            // forward transform to fill the modal coeffs
            for (int i = 0; i < m_fields.size(); ++i)
            {
                m_fields[i]->SetPhysState(true);
                m_fields[i]->FwdTrans(m_fields[i]->GetPhys(),
                                      m_fields[i]->UpdateCoeffs());
            }
        }
        break;

        default:
            break;
    }

    if (dumpInitialConditions)
    {
        // dump initial conditions to file
        std::string outname = m_sessionName + "_initial.chk";
        WriteFld(outname);
    }
}

void MMFAdvection::AdvectionBellPlane(Array<OneD, NekDouble> &outfield)
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

        // h = (h0/2)*(1+cos(pi*r/R))
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

void MMFAdvection::AdvectionBellSphere(Array<OneD, NekDouble> &outfield)
{
    int nq = m_fields[0]->GetNpoints();

    NekDouble dist, radius, cosdiff, sin_theta, cos_theta, sin_varphi,
        cos_varphi;

    Array<OneD, NekDouble> x(nq);
    Array<OneD, NekDouble> y(nq);
    Array<OneD, NekDouble> z(nq);

    m_fields[0]->GetCoords(x, y, z);

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
            outfield[j] = 0.5 * (1.0 + cos(m_pi * dist / m_radius_limit));
        }
        else
        {
            outfield[j] = 0.0;
        }
    }

    int ncoeffs = m_fields[0]->GetNcoeffs();

    Array<OneD, NekDouble> tempc(ncoeffs);
    // Smooth cosine bell
    // m_fields[0]->FwdTrans_IterPerExp(outfield, tempc);
    m_fields[0]->FwdTransLocalElmt(outfield, tempc);
    m_fields[0]->BwdTrans(tempc, outfield);
}

void MMFAdvection::Test2Dproblem(const NekDouble time,
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

void MMFAdvection::Test3Dproblem(const NekDouble time,
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

Array<OneD, NekDouble> MMFAdvection::ComputeNablaCdotVelocity()
{
    int nq = m_fields[0]->GetNpoints();

    Array<OneD, NekDouble> velcoeff(nq, 0.0);

    Array<OneD, NekDouble> Dtmp0(nq);
    Array<OneD, NekDouble> Dtmp1(nq);
    Array<OneD, NekDouble> Dtmp2(nq);
    Array<OneD, NekDouble> Drv(nq);

    Array<OneD, NekDouble> vellc(nq, 0.0);

    // m_vellc = \nabla m_vel \cdot tan_i
    Array<OneD, NekDouble> tmp(nq);
    Array<OneD, NekDouble> vessel(nq);

    for (int j = 0; j < m_shapedim; ++j)
    {
        Vmath::Zero(nq, velcoeff, 1);
        for (int k = 0; k < m_spacedim; ++k)
        {
            // a_j = tan_j cdot m_vel
            Vmath::Vvtvp(nq, &m_movingframes[j][k * nq], 1, &m_velocity[k][0],
                         1, &velcoeff[0], 1, &velcoeff[0], 1);
        }

        // d a_j / d x^k
        m_fields[0]->PhysDeriv(velcoeff, Dtmp0, Dtmp1, Dtmp2);

        for (int k = 0; k < m_spacedim; ++k)
        {
            // tan_j^k ( d a_j / d x^k )
            switch (k)
            {
                case (0):
                {
                    Vmath::Vvtvp(nq, &Dtmp0[0], 1, &m_movingframes[j][k * nq],
                                 1, &vellc[0], 1, &vellc[0], 1);
                }
                break;

                case (1):
                {
                    Vmath::Vvtvp(nq, &Dtmp1[0], 1, &m_movingframes[j][k * nq],
                                 1, &vellc[0], 1, &vellc[0], 1);
                }
                break;

                case (2):
                {
                    Vmath::Vvtvp(nq, &Dtmp2[0], 1, &m_movingframes[j][k * nq],
                                 1, &vellc[0], 1, &vellc[0], 1);
                }
                break;
            }
        }
    }
    return vellc;
}

void MMFAdvection::ComputevelodotMF(
    const Array<OneD, const Array<OneD, NekDouble>> &velocity,
    Array<OneD, Array<OneD, NekDouble>> &movingframes)
{
    int nq = m_fields[0]->GetNpoints();

    Array<OneD, NekDouble> veldotMF(nq);
    for (int j = 0; j < m_shapedim; ++j)
    {
        veldotMF = Array<OneD, NekDouble>(nq, 0.0);
        for (int k = 0; k < m_spacedim; ++k)
        {
            Vmath::Vvtvp(nq, &movingframes[j][k * nq], 1, &velocity[k][0], 1,
                         &veldotMF[0], 1, &veldotMF[0], 1);
        }

        // Modify e^i as v^i e^i
        for (int k = 0; k < m_spacedim; k++)
        {
            Vmath::Vmul(nq, &veldotMF[0], 1, &movingframes[j][k * nq], 1,
                        &movingframes[j][k * nq], 1);
        }
    }
}

NekDouble MMFAdvection::v_L2Error(unsigned int field,
                                  const Array<OneD, NekDouble> &exactsoln,
                                  bool Normalised)
{
    if (Normalised)
    {
    }

    NekDouble L2error;

    L2error = m_fields[field]->L2(m_fields[field]->GetPhys(), exactsoln);

    return L2error;
}

NekDouble MMFAdvection::v_LinfError(unsigned int field,
                                    const Array<OneD, NekDouble> &exactsoln)
{
    NekDouble Linferror;

    Linferror = m_fields[field]->Linf(m_fields[field]->GetPhys(), exactsoln);

    return Linferror;
}

// NekDouble MMFAdvection::v_L2Error(unsigned int field,
//                                   const Array<OneD, NekDouble> &exactsoln,
//                                   bool Normalised)
// {
//     int nq = m_fields[0]->GetNpoints();

//     NekDouble dist, radius, cosdiff, sin_theta, cos_theta, sin_varphi,
//         cos_varphi;

//     Array<OneD, NekDouble> x(nq);
//     Array<OneD, NekDouble> y(nq);
//     Array<OneD, NekDouble> z(nq);

//     Array<OneD, NekDouble> uBell(nq);

//     m_fields[0]->GetCoords(x, y, z);

//     NekDouble x0j, x1j, x2j;
//     for (int j = 0; j < nq; ++j)
//     {
//         x0j = x[j];
//         x1j = y[j];
//         x2j = z[j];

//         radius = sqrt(x0j * x0j + x1j * x1j + x2j * x2j);

//         sin_varphi = x1j / sqrt(x0j * x0j + x1j * x1j);
//         cos_varphi = x0j / sqrt(x0j * x0j + x1j * x1j);

//         sin_theta = x2j / radius;
//         cos_theta = sqrt(x0j * x0j + x1j * x1j) / radius;

//         cosdiff = cos_varphi * cos(m_varphi_c) + sin_varphi *
//         sin(m_varphi_c); dist    = radius * acos(sin(m_theta_c) * sin_theta +
//                              cos(m_theta_c) * cos_theta * cosdiff);

//         if (dist < m_radius_limit)
//         {
//             uBell[j] = (m_fields[field]->GetPhys())[j];
//         }

//         else
//         {
//             uBell[j] = (m_fields[field]->GetPhys())[j];
//         }

//         if (Normalised)
//         {
//             uBell[j] = uBell[j];
//         }
//     }

//     NekDouble L2error = m_fields[field]->L2(uBell, exactsoln);
//     return L2error;
// }

// NekDouble MMFAdvection::v_LinfError(unsigned int field,
//                                     const Array<OneD, NekDouble> &exactsoln)
// {
//     int nq = m_fields[0]->GetNpoints();

//     NekDouble dist, radius, cosdiff, sin_theta, cos_theta, sin_varphi,
//         cos_varphi;

//     Array<OneD, NekDouble> x(nq);
//     Array<OneD, NekDouble> y(nq);
//     Array<OneD, NekDouble> z(nq);

//     Array<OneD, NekDouble> uBell(nq);

//     m_fields[0]->GetCoords(x, y, z);

//     NekDouble x0j, x1j, x2j;
//     for (int j = 0; j < nq; ++j)
//     {
//         x0j = x[j];
//         x1j = y[j];
//         x2j = z[j];

//         radius = sqrt(x0j * x0j + x1j * x1j + x2j * x2j);

//         sin_varphi = x1j / sqrt(x0j * x0j + x1j * x1j);
//         cos_varphi = x0j / sqrt(x0j * x0j + x1j * x1j);

//         sin_theta = x2j / radius;
//         cos_theta = sqrt(x0j * x0j + x1j * x1j) / radius;

//         cosdiff = cos_varphi * cos(m_varphi_c) + sin_varphi *
//         sin(m_varphi_c); dist    = radius * acos(sin(m_theta_c) * sin_theta +
//                              cos(m_theta_c) * cos_theta * cosdiff);

//         if (dist < m_radius_limit)
//         {
//             uBell[j] = (m_fields[field]->GetPhys())[j];
//         }

//         else
//         {
//             uBell[j] = (m_fields[field]->GetPhys())[j];
//         }
//     }

//     NekDouble Linferror = m_fields[field]->Linf(uBell, exactsoln);
//     return Linferror;
// }

void MMFAdvection::v_EvaluateExactSolution(unsigned int field,
                                           Array<OneD, NekDouble> &outfield,
                                           const NekDouble time)
{
    // int ncoeffs = GetNcoeffs();

    switch (m_TestType)
    {
        case eAdvectionBell:
        {
            if (field == 0)
            {
                AdvectionBellSphere(outfield);
            }
        }
        break;

        case eTestPlane:
        {
            Test2Dproblem(time, outfield);
        }
        break;

        case eTestPlaneMassConsv:
        {
            AdvectionBellPlane(outfield);
        }
        break;

        case eTestCube:
        {
            Test3Dproblem(time, outfield);
        }
        break;

        default:
            break;
    }

    // Array<OneD, NekDouble> tmpc(ncoeffs);

    // m_fields[0]->IProductWRTBase(outfield, tmpc);
    // m_fields[0]->BwdTrans(tmpc, outfield);
}

void MMFAdvection::v_GenerateSummary(SolverUtils::SummaryList &s)
{
    MMFSystem::v_GenerateSummary(s);
    SolverUtils::AddSummaryItem(s, "TestType", TestTypeMap[m_TestType]);
    SolverUtils::AddSummaryItem(s, "Divergence Restore", m_DivergenceRestore);
    SolverUtils::AddSummaryItem(s, "Rotation Angle", m_RotAngle);

    SolverUtils::AddSummaryItem(s, "theta_c", m_theta_c);
    SolverUtils::AddSummaryItem(s, "varphi_c", m_varphi_c);
    SolverUtils::AddSummaryItem(s, "radius_limit", m_radius_limit);
}

} // namespace Nektar
