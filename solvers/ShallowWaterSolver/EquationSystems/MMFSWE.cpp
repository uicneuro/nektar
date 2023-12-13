/////////////////////////////////////////////////////////////////////////////
//
// File MMFSWE.cpp
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
// Description: MMF solve routines
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/Timer.h>
#include <boost/algorithm/string.hpp>
#include <iomanip>
#include <iostream>

#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>
#include <ShallowWaterSolver/EquationSystems/MMFSWE.h>

#include <LibUtilities/BasicUtils/Timer.h>
#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>

namespace Nektar
{
std::string MMFSWE::className =
    SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
        "MMFSWE", MMFSWE::create, "MMFSWE equation.");

MMFSWE::MMFSWE(const LibUtilities::SessionReaderSharedPtr &pSession,
               const SpatialDomains::MeshGraphSharedPtr &pGraph)
    : UnsteadySystem(pSession, pGraph), MMFSystem(pSession, pGraph)
{
    m_planeNumber = 0;
}

/**
 * @brief Initialisation object for the unsteady linear advection equation.
 */
void MMFSWE::v_InitObject(bool DeclareFields)
{
    // Call to the initialisation object
    UnsteadySystem::v_InitObject(DeclareFields);

    int nq = m_fields[0]->GetNpoints();
    std::cout << "nq = " << nq << std::endl;
    int shapedim = m_fields[0]->GetShapeDimension();
    Array<OneD, Array<OneD, NekDouble>> AniStrength(shapedim);
    for (int j = 0; j < shapedim; ++j)
    {
        AniStrength[j] = Array<OneD, NekDouble>(nq, 1.0);
    }

    MMFSystem::MMFInitObject(AniStrength);

    // ComputeSphericalVector(m_SphericalVector);
    std::cout << "============= Checking Spherical Vector ============="
              << std::endl;
    ComputeTangentUnitVector(m_SphericalVector);
    CheckMovingFrames(m_SphericalVector);
    Array<OneD, NekDouble> GeomError(nq, 0.0);
    Array<OneD, NekDouble> Ones(nq, 1.0);
    for (int k = 0; k < m_spacedim; ++k)
    {
        Vmath::Vvtvp(nq, &m_movingframes[2][k * nq], 1,
                     &m_SphericalVector[2][k * nq], 1, &GeomError[0], 1,
                     &GeomError[0], 1);
    }

    Vmath::Vsub(nq, Ones, 1, GeomError, 1, GeomError, 1);
    std::cout << "e2 cdot k: Error in L2 = " << RootMeanSquare(GeomError)
              << ", in Linf = " << Vmath::Vamax(nq, GeomError, 1) << std::endl;

    m_session->LoadParameter("Divergence Restore", m_DivergenceRestore, 0);

    if (m_DivergenceRestore >= 2)
    {
        mf_LOCSPH = Array<OneD, Array<OneD, NekDouble>>(m_mfdim);
        GetLOCALMovingframes(mf_LOCSPH);
        for (int i = 0; i < m_mfdim; ++i)
        {
            Vmath::Vcopy(3 * nq, &m_movingframes[i][0], 1, &mf_LOCSPH[i][0], 1);
        }
        ComputeAxisAlignedLOCALMovingframes(m_sphereMF, mf_LOCSPH);

        MFDotProd(m_movingframes, mf_LOCSPH, m_LOCAL_cdot_LOCSPH);
        std::cout << "LOCAL cdot LOCSPH = ( "
                  << RootMeanSquare(m_LOCAL_cdot_LOCSPH[0][0]) << " , "
                  << RootMeanSquare(m_LOCAL_cdot_LOCSPH[1][1]) << " , "
                  << RootMeanSquare(m_LOCAL_cdot_LOCSPH[2][2]) << " ) "
                  << std::endl;
    }

    // Load acceleration of gravity
    m_session->LoadParameter("Gravity", m_g, 9.81);

    // Add Coriolois effects
    m_session->LoadParameter("AddCoriolis", m_AddCoriolis, 1);

    // Add Rotation of the sphere along the pole
    m_session->LoadParameter("AddRotation", m_AddRotation, 1);

    // Spurious Diffusion Removal scheme
    m_session->LoadParameter("AddRossbyDisturbance", m_RossbyDisturbance, 0);
    m_session->LoadParameter("PurturbedJet", m_PurturbedJet, 1);

    // Define TestType
    ASSERTL0(m_session->DefinesSolverInfo("TESTTYPE"),
             "No TESTTYPE defined in session.");
    std::string TestTypeStr = m_session->GetSolverInfo("TESTTYPE");
    for (int i = 0; i < (int)SIZE_TestType; ++i)
    {
        if (boost::iequals(TestTypeMap[i], TestTypeStr))
        {
            m_TestType = (TestType)i;
            break;
        }
    }

    // Variable Setting for each test problem
    NekDouble gms = 9.80616;

    switch (m_TestType)
    {
        case eTestSteadyZonal:
        {
            NekDouble SecondToDay = 60.0 * 60.0 * 24.0;
            NekDouble rad_earth   = 6.37122 * 1000000;
            NekDouble Omegams;

            // Nondimensionalized coeffs.
            m_g = (gms * SecondToDay * SecondToDay) / rad_earth;

            m_session->LoadParameter("RotationAngle", m_alpha, 0.0);
            m_session->LoadParameter("u0", m_u0, 2.0 * m_pi / 12.0);
            m_session->LoadParameter("Omega", Omegams, 7.292 * 0.00001);
            m_Omega = Omegams * SecondToDay;

            m_session->LoadParameter("H0", m_H0, 2.94 * 10000);
            m_H0 = m_H0 / (rad_earth * gms);

            m_Hvar = (1.0 / m_g) * (m_Omega * m_u0 + 0.5 * m_u0 * m_u0);
        }
        break;

        case eTestUnsteadyZonal:
        {
            NekDouble SecondToDay = 60.0 * 60.0 * 24.0;
            NekDouble rad_earth   = 6.37122 * 1000000;
            NekDouble Omegams;

            // Nondimensionalized coeffs.
            m_g = (gms * SecondToDay * SecondToDay) / rad_earth;

            m_session->LoadParameter("RotationAngle", m_alpha, 0.0);
            m_session->LoadParameter("u0", m_u0, 2.0 * m_pi / 12.0);
            m_session->LoadParameter("Omega", Omegams, 7.292 * 0.00001);
            m_Omega = Omegams * SecondToDay;

            m_H0 = 133681.0 / (rad_earth * gms); // m^2 / s^2
            m_k2 = 10.0 / (rad_earth * gms);     // m^2 / s^2
        }
        break;

        case eTestRossbyWave:
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

            m_angfreq = 7.848 * 0.000001 * SecondToDay;
            m_K       = 7.848 * 0.000001 * SecondToDay;
        }
        break;

        case eTestIsolatedMountain:
        {
            NekDouble SecondToDay = 60.0 * 60.0 * 24.0;
            NekDouble rad_earth   = 6.37122 * 1000000;
            NekDouble Omegams;

            // Nondimensionalized coeffs.
            m_g = (gms * SecondToDay * SecondToDay) / rad_earth;

            m_session->LoadParameter("RotationAngle", m_alpha, 0.0);

            m_session->LoadParameter("u0", m_u0, 20.0);
            m_u0 = m_u0 * SecondToDay / rad_earth;

            m_session->LoadParameter("Omega", Omegams, 7.292 * 0.00001);
            m_Omega = Omegams * SecondToDay;

            m_H0  = 5960.0 / rad_earth;
            m_hs0 = 2000.0 / rad_earth;
        }
        break;

        case eTestUnstableJet:
        {
            NekDouble SecondToDay = 60.0 * 60.0 * 24.0;
            NekDouble rad_earth   = 6.37122 * 1000000;
            NekDouble Omegams;

            // Nondimensionalized coeffs.
            m_g = (gms * SecondToDay * SecondToDay) / rad_earth;

            m_session->LoadParameter("u0", m_u0, 80.0);
            m_u0 = m_u0 * SecondToDay / rad_earth;

            m_session->LoadParameter("Omega", Omegams, 7.292 * 0.00001);
            m_Omega = Omegams * SecondToDay;

            m_session->LoadParameter("H0", m_H0, 10000.0);
            m_H0 = m_H0 / rad_earth;

            m_uthetamax = 80 * SecondToDay / rad_earth;
            m_theta0    = m_pi / 7.0;
            m_theta1    = m_pi / 2.0 - m_theta0;
            m_en   = exp(-4.0 / (m_theta1 - m_theta0) / (m_theta1 - m_theta0));
            m_hbar = 120.0 / rad_earth;

            std::cout << "m_theta0 = " << m_theta0
                      << ", m_theta1 = " << m_theta1 << ", m_en = " << m_en
                      << ", m_hbar = " << m_hbar << std::endl;
        }
        break;

        default:
            break;
    }

    // TestVorticityComputation
    if (m_surfaceType == SolverUtils::eSphere)
    {
        TestVorticityComputation();
    }

    // If explicit it computes RHS and PROJECTION for the time integration
    if (m_explicitAdvection)
    {
        m_ode.DefineOdeRhs(&MMFSWE::DoOdeRhs, this);
        // m_ode.DefineOdeRhsMMF(&MMFSWE::DoOdeRhsMMF, this);
        m_ode.DefineProjection(&MMFSWE::DoOdeProjection, this);
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
MMFSWE::~MMFSWE()
{
}

void MMFSWE::v_DoSolve()
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
    Array<OneD, Array<OneD, NekDouble>> newfields(nvariables);
    Array<OneD, Array<OneD, NekDouble>> tmp(nvariables);

    // Order storage to list time-integrated fields first.
    for (i = 0; i < nvariables; ++i)
    {
        fields[i] = m_fields[m_intVariables[i]]->GetPhys();
        m_fields[m_intVariables[i]]->SetPhysState(false);

        newfields[i] = Array<OneD, NekDouble>(nq, 0.0);
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
    Array<OneD, NekDouble> dEnergy(Ntot);
    Array<OneD, NekDouble> dVorticity(Ntot);
    Array<OneD, NekDouble> dEnstrophy(Ntot);

    Array<OneD, NekDouble> zeta(nq);
    Array<OneD, Array<OneD, NekDouble>> fieldsprimitive(nvariables);
    for (int i = 0; i < nvariables; ++i)
    {
        fieldsprimitive[i] = Array<OneD, NekDouble>(nq);
    }

    // (h, hu, hv) -> (\eta, u, v)
    ConservativeToPrimitive(fields, fieldsprimitive);

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
            std::cout << "Steps: " << std::setw(8) << std::left << step + 1
                      << " "
                      << "Time: " << std::setw(12) << std::left << m_time;

            std::stringstream ss;
            ss << cpuTime << "s";
            std::cout << " CPU Time: " << std::setw(8) << std::left << ss.str()
                      << std::endl;

            // Printout Mass, Energy, Enstrophy
            ConservativeToPrimitive(fields, fieldsprimitive);

            indx = (step + 1) / m_checksteps - 1;

            // Vorticity zeta
            ComputeVorticity(fieldsprimitive[1], fieldsprimitive[2], zeta);
            dVorticity[indx] =
                std::abs(m_fields[0]->PhysIntegral(zeta) - m_Vorticity0);

            // Masss = h^*
            dMass[indx] =
                abs((ComputeMass(fieldsprimitive[0]) - m_Mass0) / m_Mass0);

            // Energy = 0.5*( h^*(u^2 + v^2) + g ( h^2 - h_s^s ) )
            dEnergy[indx] =
                abs((ComputeEnergy(fieldsprimitive[0], fieldsprimitive[1],
                                   fieldsprimitive[2]) -
                     m_Energy0) /
                    m_Energy0);

            // Enstrophy = 0.5/h^* ( \mathbf{k} \cdot (\nabla \times \mathbf{v}
            // ) + f )^2
            dEnstrophy[indx] =
                ((ComputeEnstrophy(fieldsprimitive[0], fieldsprimitive[1],
                                   fieldsprimitive[2]) -
                  m_Enstrophy0) /
                 m_Enstrophy0);

            std::cout << "dMass = " << std::setw(8) << std::left << dMass[indx]
                      << " "
                      << ", dEnergy = " << std::setw(8) << std::left
                      << dEnergy[indx] << " "
                      << ", dEnstrophy = " << std::setw(8) << std::left
                      << dEnstrophy[indx] << " "
                      << ", dVorticity = " << std::setw(8) << std::left
                      << dVorticity[indx] << std::endl
                      << std::endl;

            cpuTime = 0.0;
        }

        // (h, hu, hv) -> (\eta, u, v)
        ConservativeToPrimitive();

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

            //  (\eta, u, v) -> (h, hu, hv)
            PrimitiveToConservative();
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

    for (i = 0; i < nvariables; ++i)
    {
        m_fields[m_intVariables[i]]->SetPhys(fields[i]);
        m_fields[m_intVariables[i]]->SetPhysState(true);
    }

    // (h, hu, hv) -> (\eta, u, v)
    ConservativeToPrimitive();

    ConservativeToPrimitive(fields, fieldsprimitive);

    // Plot ThetaPhiMap
    CheckPlot_ThetaPhiMap(fieldsprimitive);

    //  Plot error map
    Checkpoint_ErrMap(m_time, fieldsprimitive);

    for (i = 0; i < nvariables; ++i)
    {
        m_fields[i]->FwdTrans(m_fields[i]->GetPhys(),
                              m_fields[i]->UpdateCoeffs());
    }
}

void MMFSWE::CheckPlot_ThetaPhiMap(
    const Array<OneD, const Array<OneD, NekDouble>> &fieldphys)
{
    int nvar    = 5;
    int nq      = m_fields[0]->GetTotPoints();
    int ncoeffs = m_fields[0]->GetNcoeffs();

    NekDouble sin_varphi, cos_varphi, sin_theta, cos_theta;
    NekDouble x0j, x1j, x2j, rad;

    std::string outname = m_sessionName + "_TheatPhiMap.chk";

    std::vector<std::string> var(nvar);
    var[0] = "phi";
    var[1] = "theta";
    var[2] = "eta";
    var[3] = "hstar";
    var[4] = "zeta";

    Array<OneD, NekDouble> x(nq);
    Array<OneD, NekDouble> y(nq);
    Array<OneD, NekDouble> z(nq);

    m_fields[0]->GetCoords(x, y, z);

    // Vorticity zeta
    Array<OneD, NekDouble> zeta(nq);
    ComputeVorticity(fieldphys[1], fieldphys[2], zeta);

    // Array<OneD, Array<OneD, NekDouble> > fieldcoeffs(nvar);
    std::vector<Array<OneD, NekDouble>> fieldcoeffs(nvar);
    for (int i = 0; i < nvar; ++i)
    {
        fieldcoeffs[i] = Array<OneD, NekDouble>(ncoeffs);
    }

    Array<OneD, NekDouble> phi(nq);
    Array<OneD, NekDouble> theta(nq);
    Array<OneD, NekDouble> eta(nq);
    Array<OneD, NekDouble> hstar(nq);
    // Compute phi and theta
    for (int i = 0; i < nq; ++i)
    {
        x0j = x[i];
        x1j = y[i];
        x2j = z[i];

        CartesianToElliptical(x0j, x1j, x2j, rad, sin_varphi, cos_varphi,
                              sin_theta, cos_theta);

        phi[i]   = atan2(sin_varphi, cos_varphi);
        theta[i] = atan2(cos_theta, sin_theta);

        // eta
        eta[i] = fieldphys[0][i];

        // eta
        hstar[i] = fieldphys[0][i] + m_depth[i];
    }

    // m_fields[0]->AdjustThetaAxis(phi);
    // m_fields[0]->AdjustThetaAxis(theta);

    m_fields[0]->FwdTrans(phi, fieldcoeffs[0]);
    m_fields[0]->FwdTrans(theta, fieldcoeffs[1]);
    m_fields[0]->FwdTrans(eta, fieldcoeffs[2]);
    m_fields[0]->FwdTrans(hstar, fieldcoeffs[3]);
    m_fields[0]->FwdTrans(zeta, fieldcoeffs[4]);

    WriteFld(outname, m_fields[0], fieldcoeffs, var);
}

void MMFSWE::DoOdeRhs(const Array<OneD, const Array<OneD, NekDouble>> &inarray,
                      Array<OneD, Array<OneD, NekDouble>> &outarray,
                      const NekDouble time)
{
    boost::ignore_unused(time);

    int nvariables = inarray.size();
    int ncoeffs    = GetNcoeffs();
    int nq         = GetTotPoints();

    // inarray in physical space
    Array<OneD, Array<OneD, NekDouble>> physarray(nvariables);
    Array<OneD, Array<OneD, NekDouble>> modarray(nvariables);

    for (int i = 0; i < nvariables; ++i)
    {
        physarray[i] = Array<OneD, NekDouble>(nq);
        modarray[i]  = Array<OneD, NekDouble>(ncoeffs);
    }

    // (h, hu, hv) -> (\eta, u, v)
    ConservativeToPrimitive(inarray, physarray);

    // Weak Divergence
    if (m_DivergenceRestore < 2)
    {
        WeakDGSWEDivergence(physarray, modarray);
    }

    else
    {
        WeakDGSWEDivergenceCNR(physarray, modarray);
    }

    // Substract 0.5 * g * H * H  / || e^m ||^2 \nalba \cdot e^m
    AddDivForGradient(physarray, modarray);

    for (int i = 0; i < nvariables; ++i)
    {
        Vmath::Neg(ncoeffs, modarray[i], 1);
    }

    // Bottom Elevation Effect
    // Add  g H \nabla m_depth
    AddElevationEffect(physarray, modarray);

    // coriolis forcing
    if (m_AddCoriolis)
    {
        AddCoriolis(physarray, modarray);
    }

    // Add terms concerning the rotation of the moving frame
    if (m_AddRotation)
    {
        AddTimeVariantFrames(physarray, m_movingframes, modarray);
    }

    if (m_DivergenceRestore > 0)
    {
        SpuriousDivRemove(physarray, m_movingframes, modarray);
    }

    for (int i = 0; i < nvariables; ++i)
    {
        m_fields[i]->MultiplyByElmtInvMass(modarray[i], modarray[i]);
        m_fields[i]->BwdTrans(modarray[i], outarray[i]);
    }
}

void MMFSWE::SpuriousDivRemove(
    const Array<OneD, const Array<OneD, NekDouble>> &physarray,
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    Array<OneD, Array<OneD, NekDouble>> &outarray)
{
    int nvariables = m_fields.size();
    int nq         = GetNpoints();
    int ncoeffs    = GetNcoeffs();

    // Create flux vector
    Array<OneD, Array<OneD, NekDouble>> fluxvector(m_shapedim);
    for (int i = 0; i < m_shapedim; ++i)
    {
        fluxvector[i] = Array<OneD, NekDouble>(nq);
    }

    // Get SWE Flux vector
    // 0: [Hu, Hv]
    // 1: [Hu^2 + 0.5 g H^2, Huv]
    // 2: [Huv, Hv^2 + 0.5 g H^2 ]
    int AddGradterm = 0;

    Array<OneD, NekDouble> velvector(m_spacedim * nq);
    for (int i = 0; i < nvariables; ++i)
    {
        // int AddGradterm=0;
        GetSWEFluxVector(i, physarray, fluxvector, AddGradterm);

        // velocity = H \vec{u} = fluxvector[0] e^1 + fluxvector[1] e^2
        velvector = Array<OneD, NekDouble>(m_spacedim * nq, 0.0);
        for (int j = 0; j < m_shapedim; ++j)
        {
            for (int k = 0; k < m_spacedim; ++k)
            {
                Vmath::Vvtvp(nq, &fluxvector[j][0], 1, &movingframes[j][k * nq],
                             1, &velvector[k * nq], 1, &velvector[k * nq], 1);
            }
        }

        Array<OneD, NekDouble> tmpc(ncoeffs);
        Array<OneD, NekDouble> tmp(nq, 0.0);
        switch (m_DivergenceRestore)
        {
            case 1:
            {
                // m_movingframes = LOCAL axis
                tmp = ComputeSpuriousDivergence(
                    movingframes, m_SphericalVector[2], velvector);
            }
            break;

            case 2:
            {
                // tmp = ComputeSpuriousDivergence(movingframes, mf_LOCAL[2],
                // velvector);
                tmp = ComputeSpuriousDivergence(mf_LOCSPH, m_movingframes[2],
                                                velvector);
                // Vmath::Neg(nq, tmp, 1);
            }
            break;

            case 3:
            {
                // m_movingframes = LOCAL axis
                // tmp = ComputeSpuriousDivergence(mf_LOCSPH, movingframes[2],
                // velvector);
                tmp = ComputeSpuriousDivergence(mf_LOCSPH, m_SphericalVector[2],
                                                velvector);
                Vmath::Neg(nq, tmp, 1);
            }
            break;

            default:
                break;
        }
        m_fields[i]->IProductWRTBase(tmp, tmpc);
        Vmath::Vadd(ncoeffs, tmpc, 1, outarray[i], 1, outarray[i], 1);
    }
}

void MMFSWE::WeakDGSWEDivergence(
    const Array<OneD, Array<OneD, NekDouble>> &InField,
    Array<OneD, Array<OneD, NekDouble>> &OutField)
{
    int nvariables      = m_fields.size();
    int nq              = GetNpoints();
    int ncoeffs         = GetNcoeffs();
    int nTracePointsTot = GetTraceNpoints();

    Array<OneD, Array<OneD, NekDouble>> fluxvector(m_shapedim);
    Array<OneD, Array<OneD, NekDouble>> physfield(nvariables);

    for (int j = 0; j < m_shapedim; ++j)
    {
        fluxvector[j] = Array<OneD, NekDouble>(nq);
    }

    // InField is Primitive
    for (int i = 0; i < nvariables; ++i)
    {
        physfield[i] = InField[i];
    }

    // Get the ith component of the  flux vector in (physical space)
    // fluxvector[0] = component for e^1 cdot \nabla \varphi
    // fluxvector[1] = component for e^2 cdot \nabla \varphi
    Array<OneD, NekDouble> tmpc(ncoeffs);
    // Compute Divergence Components
    for (int i = 0; i < nvariables; ++i)
    {
        // 0: [Hu, Hv]
        // 1: [Hu^2 + 0.5gH^2, Huv ]
        // 2: [Huv, Hv^2 + 0.5 g H^2,]
        GetSWEFluxVector(i, physfield, fluxvector, 1);

        OutField[i] = Array<OneD, NekDouble>(ncoeffs, 0.0);
        for (int j = 0; j < m_shapedim; ++j)
        {
            // Directional derivation with respect to the j'th moving frame
            // tmp_j = ( \nabla \phi, fluxvector[j] \mathbf{e}^j )
            m_fields[i]->IProductWRTDirectionalDerivBase(m_movingframes[j],
                                                         fluxvector[j], tmpc);
            Vmath::Vadd(ncoeffs, &tmpc[0], 1, &OutField[i][0], 1,
                        &OutField[i][0], 1);
        }

        // std::cout << "i = " << i << ", OutField = " <<
        // RootMeanSquare(OutField[i]) << std::endl;
    }

    // Numerical Flux
    Array<OneD, Array<OneD, NekDouble>> numfluxFwd(nvariables);
    Array<OneD, Array<OneD, NekDouble>> numfluxBwd(nvariables);

    for (int i = 0; i < nvariables; ++i)
    {
        numfluxFwd[i] = Array<OneD, NekDouble>(nTracePointsTot);
        numfluxBwd[i] = Array<OneD, NekDouble>(nTracePointsTot);
    }

    // NumericalSWEFlux(physfield, numfluxFwd, numfluxBwd);
    NumericalSWEFlux(physfield, numfluxFwd, numfluxBwd);

    // Evaulate  <\phi, \hat{F}\cdot n> - OutField[i]
    for (int i = 0; i < nvariables; ++i)
    {
        Vmath::Neg(ncoeffs, OutField[i], 1);
        m_fields[i]->AddFwdBwdTraceIntegral(numfluxFwd[i], numfluxBwd[i],
                                            OutField[i]);
        m_fields[i]->SetPhysState(false);
    }
}

void MMFSWE::WeakDGSWEDivergenceCNR(
    const Array<OneD, Array<OneD, NekDouble>> &InField,
    Array<OneD, Array<OneD, NekDouble>> &OutField)
{
    int nvariables      = m_fields.size();
    int nq              = GetNpoints();
    int ncoeffs         = GetNcoeffs();
    int nTracePointsTot = GetTraceNpoints();

    Array<OneD, Array<OneD, NekDouble>> Divfluxvector(m_shapedim);
    Array<OneD, Array<OneD, NekDouble>> Gradfluxvector(m_shapedim);

    Array<OneD, Array<OneD, NekDouble>> physfield(nvariables);

    for (int j = 0; j < m_shapedim; ++j)
    {
        Divfluxvector[j]  = Array<OneD, NekDouble>(nq);
        Gradfluxvector[j] = Array<OneD, NekDouble>(nq);
    }

    // InField is Primitive
    for (int i = 0; i < nvariables; ++i)
    {
        physfield[i] = InField[i];
    }

    // Get the ith component of the  flux vector in (physical space)
    // fluxvector[0] = component for e^1 cdot \nabla \varphi
    // fluxvector[1] = component for e^2 cdot \nabla \varphi
    Array<OneD, NekDouble> Divtmpc(ncoeffs);
    Array<OneD, NekDouble> Fluxtmpc(ncoeffs);

    // Compute Divergence Components

    for (int i = 0; i < nvariables; ++i)
    {
        OutField[i] = Array<OneD, NekDouble>(ncoeffs, 0.0);

        // 0: [Hu, Hv]
        // 1: [Hu^2 + 0.5gH^2, Huv ]
        // 2: [Huv, Hv^2 + 0.5 g H^2,]
        GetSWEFluxVector(i, physfield, Divfluxvector);
        GetSWEFluxVector(i + 10, physfield, Gradfluxvector);

        for (int j = 0; j < m_shapedim; ++j)
        {
            // Directional derivation with respect to the j'th moving frame
            // tmp_j = ( \nabla \phi, fluxvector[j] \mathbf{e}^j )
            m_fields[i]->IProductWRTDirectionalDerivBase(
                mf_LOCSPH[j], Divfluxvector[j], Divtmpc);
            Vmath::Vadd(ncoeffs, &Divtmpc[0], 1, &OutField[i][0], 1,
                        &OutField[i][0], 1);

            m_fields[i]->IProductWRTDirectionalDerivBase(
                m_movingframes[j], Gradfluxvector[j], Fluxtmpc);
            Vmath::Vadd(ncoeffs, &Fluxtmpc[0], 1, &OutField[i][0], 1,
                        &OutField[i][0], 1);
        }
    }

    // Numerical Flux
    Array<OneD, Array<OneD, NekDouble>> numfluxFwd(nvariables);
    Array<OneD, Array<OneD, NekDouble>> numfluxBwd(nvariables);

    for (int i = 0; i < nvariables; ++i)
    {
        numfluxFwd[i] = Array<OneD, NekDouble>(nTracePointsTot);
        numfluxBwd[i] = Array<OneD, NekDouble>(nTracePointsTot);
    }

    // NumericalSWEFlux(physfield, numfluxFwd, numfluxBwd);
    NumericalSWEFlux(physfield, numfluxFwd, numfluxBwd);

    // Evaulate  <\phi, \hat{F}\cdot n> - OutField[i]
    for (int i = 0; i < nvariables; ++i)
    {
        Vmath::Neg(ncoeffs, OutField[i], 1);
        m_fields[i]->AddFwdBwdTraceIntegral(numfluxFwd[i], numfluxBwd[i],
                                            OutField[i]);
        m_fields[i]->SetPhysState(false);
    }
}

// Substract 0.5 * g * H * H  / || e^m ||^2 \nalba \cdot e^m
void MMFSWE::AddDivForGradient(Array<OneD, Array<OneD, NekDouble>> &physarray,
                               Array<OneD, Array<OneD, NekDouble>> &outarray)
{
    // routine works for both primitive and conservative formulations
    int ncoeffs = outarray[0].size();
    int nq      = physarray[0].size();

    Array<OneD, NekDouble> h(nq);
    Array<OneD, NekDouble> tmp(nq);
    Array<OneD, Array<OneD, NekDouble>> fluxvector(m_shapedim);
    for (int i = 0; i < m_shapedim; ++i)
    {
        fluxvector[i] = Array<OneD, NekDouble>(nq);
    }

    // Get 0.5 g H*H / || e^m ||^2
    GetSWEFluxVector(3, physarray, fluxvector);

    Array<OneD, NekDouble> tmpc(ncoeffs);
    for (int j = 0; j < m_shapedim; ++j)
    {
        Vmath::Vmul(nq, &fluxvector[0][0], 1, &m_DivMF[j][0], 1, &tmp[0], 1);

        Vmath::Neg(nq, &tmp[0], 1);
        m_fields[0]->IProductWRTBase(tmp, tmpc);

        Vmath::Vadd(ncoeffs, outarray[j + 1], 1, tmpc, 1, outarray[j + 1], 1);
    }
}

// 0: [Hu, Hv]
// 1: [Hu^2 + 0.5gH^2, Huv ]
// 2: [Huv, Hv^2 + 0.5 g H^2,]

void MMFSWE::GetSWEFluxVector(
    const int i, const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &flux, const int AddGradientTerm)
{
    int nvar = 3; // only the dependent variables
    int nq   = m_fields[0]->GetTotPoints();
    Array<OneD, NekDouble> tmp(nq);

    Array<OneD, Array<OneD, NekDouble>> physfield(nvar);
    for (int i = 0; i < nvar; ++i)
    {
        physfield[i] = Array<OneD, NekDouble>(nq, 0.0);
    }

    // if m_DivergenceRestore>=2, we use a sphericl coordinate expansion
    ComputephysfieldinMF(inarray, physfield, m_DivergenceRestore);

    // h
    Vmath::Vadd(nq, physfield[0], 1, m_depth, 1, tmp, 1);

    switch (i)
    {
        // flux function for the h equation = [(\eta + d ) u, (\eta + d) v ]
        case 0:
        {
            // hu in flux 0
            Vmath::Vmul(nq, tmp, 1, physfield[1], 1, flux[0], 1);

            // hv in flux 1
            Vmath::Vmul(nq, tmp, 1, physfield[2], 1, flux[1], 1);
        }
        break;

        // flux function for the hu equation = [Hu^2 + 0.5 g H^2, Huv]
        case 1:
        {
            // hu in flux 1
            Vmath::Vmul(nq, tmp, 1, physfield[1], 1, flux[1], 1);

            // huu in flux 0
            Vmath::Vmul(nq, flux[1], 1, physfield[1], 1, flux[0], 1);

            if (AddGradientTerm)
            {
                //  hh in tmp
                Vmath::Vmul(nq, tmp, 1, tmp, 1, tmp, 1);

                // huu + 0.5 g hh in flux 0
                // Daxpy overwrites flux[0] on exit
                Blas::Daxpy(nq, 0.5 * m_g, tmp, 1, flux[0], 1);
            }

            // huv in flux 1
            Vmath::Vmul(nq, flux[1], 1, physfield[2], 1, flux[1], 1);
        }
        break;

        // flux function for the hv equation = [Huv, Hv^2 + 0.5 g H^2,]
        case 2:
        {
            // hv in flux 0
            Vmath::Vmul(nq, tmp, 1, physfield[2], 1, flux[0], 1);

            // hvv in flux 1
            Vmath::Vmul(nq, flux[0], 1, physfield[2], 1, flux[1], 1);

            // huv in flux 0
            Vmath::Vmul(nq, flux[0], 1, physfield[1], 1, flux[0], 1);

            if (AddGradientTerm)
            {
                //  hh in tmp
                Vmath::Vmul(nq, tmp, 1, tmp, 1, tmp, 1);

                // hvv + 0.5 g hh in flux 1
                Blas::Daxpy(nq, 0.5 * m_g, tmp, 1, flux[1], 1);
            }
        }
        break;

        // SWE Gradient flux function 0.5 g h * h
        case 3:
        {
            Array<OneD, NekDouble> h2(nq);
            flux[0] = Array<OneD, NekDouble>(nq, 0.0);

            //  hh in tmp
            Vmath::Vmul(nq, tmp, 1, tmp, 1, h2, 1);

            // 0.5 g hh in flux 0
            Blas::Daxpy(nq, 0.5 * m_g, h2, 1, flux[0], 1);
        }
        break;

        case 10:
        {
            Vmath::Fill(nq, 0.0, flux[0], 1);
            Vmath::Fill(nq, 0.0, flux[1], 1);
        }
        break;

        // flux function for the hu equation = [Hu^2 + 0.5 g H^2, Huv]
        case 11:
        {
            //  hh in tmp
            Vmath::Vmul(nq, tmp, 1, tmp, 1, tmp, 1);

            // 0.5 g hh in flux 0
            // Daxpy overwrites flux[0] on exit
            Vmath::Smul(nq, 0.5 * m_g, tmp, 1, flux[0], 1);
            Vmath::Fill(nq, 0.0, flux[1], 1);
        }
        break;

        case 12:
        {
            //  hh in tmp
            Vmath::Vmul(nq, tmp, 1, tmp, 1, tmp, 1);

            Vmath::Smul(nq, 0.5 * m_g, tmp, 1, flux[1], 1);
            Vmath::Fill(nq, 0.0, flux[0], 1);
        }
        break;

        default:
            break;
    }
}

void MMFSWE::ComputephysfieldinMF(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &physfield, const int DivergenceRestore)
{
    int nvar = 3; // only the dependent variables
    int nq   = m_fields[0]->GetTotPoints();

    Vmath::Vcopy(nq, &inarray[0][0], 1, &physfield[0][0], 1);
    if (DivergenceRestore >= 2)
    {
        Vmath::Vmul(nq, &inarray[1][0], 1, &m_LOCAL_cdot_LOCSPH[0][0][0], 1,
                    &physfield[1][0], 1);
        Vmath::Vvtvp(nq, &inarray[2][0], 1, &m_LOCAL_cdot_LOCSPH[1][0][0], 1,
                     &physfield[1][0], 1, &physfield[1][0], 1);

        Vmath::Vmul(nq, &inarray[1][0], 1, &m_LOCAL_cdot_LOCSPH[0][1][0], 1,
                    &physfield[2][0], 1);
        Vmath::Vvtvp(nq, &inarray[2][0], 1, &m_LOCAL_cdot_LOCSPH[1][1][0], 1,
                     &physfield[2][0], 1, &physfield[2][0], 1);
    }

    else
    {
        for (int i = 1; i < nvar; ++i)
        {
            Vmath::Vcopy(nq, &inarray[i][0], 1, &physfield[i][0], 1);
        }
    }
}

// newu = u (e^1 \cdot e^1_new) + v (e^2 \cdot e^1_new)
// newv = u (e^1 \cdot e^2_new) + v (e^2 \cdot e^2_new)
// void MMFSWE::Convert_LOCAL_TO_LOCALSPHERE(const Array<OneD, const NekDouble>
// &physu, const Array<OneD, const NekDouble> &physv, Array<OneD, NekDouble>
// &newu, Array<OneD, NekDouble> &newv)
// {
//     int nq = m_fields[0]->GetTotPoints();

//     newu = Array<OneD, NekDouble>(nq, 0.0);
//     newv = Array<OneD, NekDouble>(nq, 0.0);

//     // newu = u (e^1 \cdot e^1_new) + v (e^2 \cdot e^1_new)
//     Vmath::Vvtvp(nq, &physu[0], 1, &m_LOCAL_cdot_CoordAxisDivMF[0][0][0], 1,
//     &newu[0], 1, &newu[0], 1); Vmath::Vvtvp(nq, &physv[0], 1,
//     &m_LOCAL_cdot_CoordAxisDivMF[1][0][0], 1, &newu[0], 1, &newu[0], 1);

//     // newv = u (e^1 \cdot e^2_new) + v (e^2 \cdot e^2_new)
//     Vmath::Vvtvp(nq, &physu[0], 1, &m_LOCAL_cdot_CoordAxisDivMF[0][1][0], 1,
//     &newv[0], 1, &newv[0], 1); Vmath::Vvtvp(nq, &physv[0], 1,
//     &m_LOCAL_cdot_CoordAxisDivMF[1][1][0], 1, &newv[0], 1, &newv[0], 1);
// }

void MMFSWE::NumericalSWEFlux(Array<OneD, Array<OneD, NekDouble>> &physfield,
                              Array<OneD, Array<OneD, NekDouble>> &numfluxFwd,
                              Array<OneD, Array<OneD, NekDouble>> &numfluxBwd)
{
    int i, k;
    int nTraceNumPoints = GetTraceTotPoints();
    int nvariables      = 3; // only the dependent variables

    // get temporary arrays
    Array<OneD, Array<OneD, NekDouble>> Fwd(nvariables);
    Array<OneD, Array<OneD, NekDouble>> Bwd(nvariables);
    Array<OneD, NekDouble> DepthFwd(nTraceNumPoints);
    Array<OneD, NekDouble> DepthBwd(nTraceNumPoints);

    for (i = 0; i < nvariables; ++i)
    {
        Fwd[i] = Array<OneD, NekDouble>(nTraceNumPoints);
        Bwd[i] = Array<OneD, NekDouble>(nTraceNumPoints);
    }

    // get the physical values at the trace
    for (i = 0; i < nvariables; ++i)
    {
        m_fields[i]->GetFwdBwdTracePhys(physfield[i], Fwd[i], Bwd[i]);

        // Copy Fwd to Bwd at boundaries
        CopyBoundaryTrace(Fwd[i], Bwd[i], SolverUtils::eFwdEQBwd);
    }
    m_fields[0]->GetFwdBwdTracePhys(m_depth, DepthFwd, DepthBwd);
    CopyBoundaryTrace(DepthFwd, DepthBwd, SolverUtils::eFwdEQBwd);

    // Compute ncdotMFFwd, ncdotMFBwd, nperpcdotMFFwd, nperpcdotMFBwd
    Array<OneD, Array<OneD, NekDouble>> ncdotMFFwd;
    Array<OneD, Array<OneD, NekDouble>> ncdotMFBwd;

    Array<OneD, Array<OneD, NekDouble>> nperpcdotMFFwd;
    Array<OneD, Array<OneD, NekDouble>> nperpcdotMFBwd;

    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> MFtraceFwd;
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> MFtraceBwd;

    // note that we are using the same depth - i.e. the depth is assumed
    // continuous...
    switch (m_upwindType)
    {
        case SolverUtils::eHLLC:
        {
            // rotate the values to the normal direction
            NekDouble tmpX, tmpY;

            // Fwd[1] = hu^+ = ( h^1^+ (e^1 \cdot n) + h^2^+ ( e^2 \cdot n) )
            // Fwd[2] = hv^+ = ( h^1^+ (e^1 \cdot n^{\perp}) + h^2^+ ( e^2 \cdot
            // n^{\perp}) )
            // Bwd[1] = hu^- = ( h^1^- (e^1 \cdot n) + h^2^- ( e^2 \cdot n) )
            // Bwd[2] = hv^- = ( h^1^- (e^1 \cdot n^{\perp}) + h^2^- ( e^2 \cdot
            // n^{\perp}) )

            for (k = 0; k < nTraceNumPoints; ++k)
            {
                tmpX = Fwd[1][k] * m_ncdotMFFwd[0][k] +
                       Fwd[2][k] * m_ncdotMFFwd[1][k];
                tmpY = Fwd[1][k] * m_nperpcdotMFFwd[0][k] +
                       Fwd[2][k] * m_nperpcdotMFFwd[1][k];
                Fwd[1][k] = tmpX;
                Fwd[2][k] = tmpY;

                tmpX = Bwd[1][k] * m_ncdotMFBwd[0][k] +
                       Bwd[2][k] * m_ncdotMFBwd[1][k];
                tmpY = Bwd[1][k] * m_nperpcdotMFBwd[0][k] +
                       Bwd[2][k] * m_nperpcdotMFBwd[1][k];
                Bwd[1][k] = tmpX;
                Bwd[2][k] = tmpY;
            }

            // Solve the Riemann problem
            NekDouble denomFwd, denomBwd;
            Array<OneD, NekDouble> numfluxF(nvariables);
            Array<OneD, NekDouble> numfluxB(nvariables);

            NekDouble eF1n, eF2n, eB1n, eB2n;
            NekDouble eF1t, eF2t, eB1t, eB2t;
            for (k = 0; k < nTraceNumPoints; ++k)
            {
                RiemannSolverHLLC(k, Fwd[0][k] + DepthFwd[k], Fwd[1][k],
                                  Fwd[2][k], Bwd[0][k] + DepthFwd[k], Bwd[1][k],
                                  Bwd[2][k], numfluxF, numfluxB);

                // uflux = 1/[ ( e^1 \cdot n) ( e^2 \cdot n^{\perp} ) - (e^2
                // \cdot n)(e^1 \cdot n^{\perp} ) ] ( e^2 \cdot n^{\perp} huflux
                // - e^2 \cdot n hvflux
                // vflux = -1/[ ( e^1 \cdot n) ( e^2 \cdot n^{\perp} ) - (e^2
                // \cdot n)(e^1 \cdot n^{\perp} ) ] ( -e^1 \cdot n^{\perp}
                // huflux + e^1 \cdot n hvflux
                eF1n = m_ncdotMFFwd[0][k];
                eF2n = m_ncdotMFFwd[1][k];
                eB1n = m_ncdotMFBwd[0][k];
                eB2n = m_ncdotMFBwd[1][k];

                eF1t = m_nperpcdotMFFwd[0][k];
                eF2t = m_nperpcdotMFFwd[1][k];
                eB1t = m_nperpcdotMFBwd[0][k];
                eB2t = m_nperpcdotMFBwd[1][k];

                denomFwd = eF1n * eF2t - eF2n * eF1t;
                denomBwd = eB1n * eB2t - eB2n * eB1t;

                numfluxFwd[0][k] = numfluxF[0];
                numfluxFwd[1][k] = (1.0 / denomFwd) *
                                   (eF2t * numfluxF[1] - eF2n * numfluxF[2]);
                numfluxFwd[2][k] =
                    (1.0 / denomFwd) *
                    (-1.0 * eF1t * numfluxF[1] + eF1n * numfluxF[2]);

                numfluxBwd[0][k] = 1.0 * numfluxB[0];
                numfluxBwd[1][k] = (1.0 / denomBwd) *
                                   (eB2t * numfluxB[1] - eB2n * numfluxB[2]);
                numfluxBwd[2][k] =
                    (1.0 / denomBwd) *
                    (-1.0 * eB1t * numfluxB[1] + eB1n * numfluxB[2]);
            }
        }
        break;

        case SolverUtils::eAverage:
        case SolverUtils::eLaxFriedrich:
        case SolverUtils::eRusanov:
        {
            Array<OneD, NekDouble> numfluxF(nvariables * (m_shapedim + 1));
            Array<OneD, NekDouble> numfluxB(nvariables * (m_shapedim + 1));

            NekDouble MF1x, MF1y, MF1z, MF2x, MF2y, MF2z;
            NekDouble MB1x, MB1y, MB1z, MB2x, MB2y, MB2z;
            NekDouble MageF1, MageF2, MageB1, MageB2;
            NekDouble eF1_cdot_eB1, eF1_cdot_eB2, eF2_cdot_eB1, eF2_cdot_eB2;
            NekDouble velL, velR;
            for (k = 0; k < nTraceNumPoints; ++k)
            {
                MF1x = m_MFtraceFwd[0][0][k];
                MF1y = m_MFtraceFwd[0][1][k];
                MF1z = m_MFtraceFwd[0][2][k];

                MF2x = m_MFtraceFwd[1][0][k];
                MF2y = m_MFtraceFwd[1][1][k];
                MF2z = m_MFtraceFwd[1][2][k];

                MB1x = m_MFtraceBwd[0][0][k];
                MB1y = m_MFtraceBwd[0][1][k];
                MB1z = m_MFtraceBwd[0][2][k];

                MB2x = m_MFtraceBwd[1][0][k];
                MB2y = m_MFtraceBwd[1][1][k];
                MB2z = m_MFtraceBwd[1][2][k];

                // MFtrace = MFtrace [ j*spacedim + k ], j = shape, k = sapce
                MageF1 = MF1x * MF1x + MF1y * MF1y + MF1z * MF1z;
                MageF2 = MF2x * MF2x + MF2y * MF2y + MF2z * MF2z;
                MageB1 = MB1x * MB1x + MB1y * MB1y + MB1z * MB1z;
                MageB2 = MB2x * MB2x + MB2y * MB2y + MB2z * MB2z;

                eF1_cdot_eB1 = MF1x * MB1x + MF1y * MB1y + MF1z * MB1z;
                eF1_cdot_eB2 = MF1x * MB2x + MF1y * MB2y + MF1z * MB2z;
                eF2_cdot_eB1 = MF2x * MB1x + MF2y * MB1y + MF2z * MB1z;
                eF2_cdot_eB2 = MF2x * MB2x + MF2y * MB2y + MF2z * MB2z;

                if (m_upwindType == SolverUtils::eAverage)
                {
                    AverageFlux(k, Fwd[0][k] + DepthFwd[k], Fwd[1][k],
                                Fwd[2][k], Bwd[0][k] + DepthFwd[k], Bwd[1][k],
                                Bwd[2][k], numfluxF, numfluxB);
                }

                else if (m_upwindType == SolverUtils::eLaxFriedrich)
                {
                    velL = Fwd[1][k] * m_ncdotMFFwd[0][k] +
                           Fwd[2][k] * m_ncdotMFFwd[1][k];
                    velR = -1.0 * (Bwd[1][k] * m_ncdotMFBwd[0][k] +
                                   Bwd[2][k] * m_ncdotMFBwd[1][k]);

                    LaxFriedrichFlux(k, velL, velR, Fwd[0][k] + DepthFwd[k],
                                     Fwd[1][k], Fwd[2][k],
                                     Bwd[0][k] + DepthFwd[k], Bwd[1][k],
                                     Bwd[2][k], MageF1, MageF2, MageB1, MageB2,
                                     eF1_cdot_eB1, eF1_cdot_eB2, eF2_cdot_eB1,
                                     eF2_cdot_eB2, numfluxF, numfluxB);
                }

                else if (m_upwindType == SolverUtils::eRusanov)
                {
                    RusanovFlux(k, Fwd[0][k] + DepthFwd[k], Fwd[1][k],
                                Fwd[2][k], Bwd[0][k] + DepthFwd[k], Bwd[1][k],
                                Bwd[2][k], numfluxF, numfluxB);
                }

                int indx;
                NekDouble tmpF0, tmpF1, tmpB0, tmpB1;
                for (i = 0; i < nvariables; ++i)
                {
                    indx = i * (m_shapedim + 1);

                    tmpF0 = numfluxF[indx] * m_ncdotMFFwd[0][k];
                    tmpF1 = numfluxF[indx + 1] * m_ncdotMFFwd[1][k];

                    tmpB0 = numfluxB[indx] * m_ncdotMFBwd[0][k];
                    tmpB1 = numfluxB[indx + 1] * m_ncdotMFBwd[1][k];

                    numfluxFwd[i][k] = tmpF0 + tmpF1 + numfluxF[indx + 2];
                    numfluxBwd[i][k] = tmpB0 + tmpB1 + numfluxB[indx + 2];
                }
            }
        }
        break;

        default:
        {
            ASSERTL0(false, "populate switch statement for upwind flux");
        }
        break;
    }
}

void MMFSWE::RiemannSolverHLLC(const int index, NekDouble hL, NekDouble uL,
                               NekDouble vL, NekDouble hR, NekDouble uR,
                               NekDouble vR, Array<OneD, NekDouble> &numfluxF,
                               Array<OneD, NekDouble> &numfluxB)
{
    boost::ignore_unused(index);

    NekDouble g = m_g;

    NekDouble cL = sqrt(g * hL);
    NekDouble cR = sqrt(g * hR);

    NekDouble uRF, vRF, uLB, vLB;
    NekDouble hstarF, hstarB;

    // Temporary assignments
    uRF = uR;
    vRF = vR;
    uLB = uL;
    vLB = vL;

    // the two-rarefaction wave assumption
    hstarF = 0.5 * (cL + cR) + 0.25 * (uL - uRF);
    hstarF *= hstarF;
    hstarF *= (1.0 / g);

    hstarB = 0.5 * (cL + cR) + 0.25 * (uLB - uR);
    hstarB *= hstarB;
    hstarB *= (1.0 / g);

    NekDouble hfluxF, hufluxF, hvfluxF;
    NekDouble hfluxB, hufluxB, hvfluxB;
    Computehhuhvflux(hL, uL, vL, hR, uRF, vRF, hstarF, hfluxF, hufluxF,
                     hvfluxF);
    Computehhuhvflux(hL, uLB, vLB, hR, uR, vR, hstarB, hfluxB, hufluxB,
                     hvfluxB);

    numfluxF[0] = hfluxF;
    numfluxF[1] = hufluxF;
    numfluxF[2] = hvfluxF;

    numfluxB[0] = hfluxB;
    numfluxB[1] = hufluxB;
    numfluxB[2] = hvfluxB;
}

void MMFSWE::Computehhuhvflux(NekDouble hL, NekDouble uL, NekDouble vL,
                              NekDouble hR, NekDouble uR, NekDouble vR,
                              NekDouble hstar, NekDouble &hflux,
                              NekDouble &huflux, NekDouble &hvflux)
{
    NekDouble g = m_g;

    NekDouble hC, huC, hvC, SL, SR, Sstar;
    NekDouble cL = sqrt(g * hL);
    NekDouble cR = sqrt(g * hR);

    // Compute SL
    if (hstar > hL)
        SL = uL - cL * sqrt(0.5 * ((hstar * hstar + hstar * hL) / (hL * hL)));
    else
        SL = uL - cL;

    // Compute SR
    if (hstar > hR)
        SR = uR + cR * sqrt(0.5 * ((hstar * hstar + hstar * hR) / (hR * hR)));
    else
        SR = uR + cR;

    if (fabs(hR * (uR - SR) - hL * (uL - SL)) <= 1.0e-15)
        Sstar = 0.0;
    else
        Sstar = (SL * hR * (uR - SR) - SR * hL * (uL - SL)) /
                (hR * (uR - SR) - hL * (uL - SL));

    if (SL >= 0)
    {
        hflux  = hL * uL;
        huflux = uL * uL * hL + 0.5 * g * hL * hL;
        hvflux = hL * uL * vL;
    }
    else if (SR <= 0)
    {
        hflux  = hR * uR;
        huflux = uR * uR * hR + 0.5 * g * hR * hR;
        hvflux = hR * uR * vR;
    }
    else
    {
        if ((SL < 0) && (Sstar >= 0))
        {
            hC  = hL * ((SL - uL) / (SL - Sstar));
            huC = hC * Sstar;
            hvC = hC * vL;

            hflux  = hL * uL + SL * (hC - hL);
            huflux = (uL * uL * hL + 0.5 * g * hL * hL) + SL * (huC - hL * uL);
            hvflux = (uL * vL * hL) + SL * (hvC - hL * vL);
        }
        else
        {
            hC  = hR * ((SR - uR) / (SR - Sstar));
            huC = hC * Sstar;
            hvC = hC * vR;

            hflux  = hR * uR + SR * (hC - hR);
            huflux = (uR * uR * hR + 0.5 * g * hR * hR) + SR * (huC - hR * uR);
            hvflux = (uR * vR * hR) + SR * (hvC - hR * vR);
        }
    }
}

void MMFSWE::AverageFlux(const int index, NekDouble hL, NekDouble uL,
                         NekDouble vL, NekDouble hR, NekDouble uR, NekDouble vR,
                         Array<OneD, NekDouble> &numfluxF,
                         Array<OneD, NekDouble> &numfluxB)
{
    NekDouble MageF1, MageF2, MageB1, MageB2;
    NekDouble eF1_cdot_eB1, eF1_cdot_eB2;
    NekDouble eF2_cdot_eB1, eF2_cdot_eB2;

    NekDouble g = m_g;
    NekDouble uRF, vRF, uLB, vLB;

    ComputeMagAndDot(index, MageF1, MageF2, MageB1, MageB2, eF1_cdot_eB1,
                     eF1_cdot_eB2, eF2_cdot_eB1, eF2_cdot_eB2);

    // uRF = uR component in moving frames e^{Fwd}
    // vRF = vR component in moving frames e^{Fwd}
    uRF = (uR * eF1_cdot_eB1 + vR * eF1_cdot_eB2) / MageF1;
    vRF = (uR * eF2_cdot_eB1 + vR * eF2_cdot_eB2) / MageF2;

    numfluxF[0] = 0.5 * (hL * uL + hR * uRF);
    numfluxF[1] = 0.5 * (hL * vL + hR * vRF);
    numfluxF[2] = 0.0;

    numfluxF[3] =
        0.5 * (hL * uL * uL + hR * uRF * uRF + 0.5 * g * (hL * hL + hR * hR));
    numfluxF[4] = 0.5 * (hL * uL * vL + hR * uRF * vRF);
    numfluxF[5] = 0.0;

    numfluxF[6] = 0.5 * (hL * uL * vL + hR * uRF * vRF);
    numfluxF[7] =
        0.5 * (hL * vL * vL + hR * vRF * vRF + 0.5 * g * (hL * hL + hR * hR));
    numfluxF[8] = 0.0;

    // uLB = uL component in moving frames e^{Bwd}
    // vLB = vL component in moving frames e^{Bwd}
    uLB = (uL * eF1_cdot_eB1 + vL * eF2_cdot_eB1) / MageB1;
    vLB = (uL * eF1_cdot_eB2 + vL * eF2_cdot_eB2) / MageB2;

    numfluxB[0] = 0.5 * (hR * uR + hR * uLB);
    numfluxB[1] = 0.5 * (hR * vR + hR * vLB);
    numfluxB[2] = 0.0;

    numfluxB[3] =
        0.5 * (hR * uR * uR + hR * uLB * uLB + 0.5 * g * (hR * hR + hL * hL));
    numfluxB[4] = 0.5 * (hR * uR * vR + hR * uLB * vLB);
    numfluxB[5] = 0.0;

    numfluxB[6] = 0.5 * (hR * uR * vR + hR * uLB * vLB);
    numfluxB[7] =
        0.5 * (hR * vR * vR + hR * vLB * vLB + 0.5 * g * (hR * hR + hL * hL));
    numfluxB[8] = 0.0;
}

void MMFSWE::LaxFriedrichFlux(
    const int index, const NekDouble velL, const NekDouble velR,
    const NekDouble hL, const NekDouble uL, const NekDouble vL,
    const NekDouble hR, const NekDouble uR, const NekDouble vR,
    const NekDouble MageF1, const NekDouble MageF2, const NekDouble MageB1,
    const NekDouble MageB2, const NekDouble eF1_cdot_eB1,
    const NekDouble eF1_cdot_eB2, const NekDouble eF2_cdot_eB1,
    const NekDouble eF2_cdot_eB2, Array<OneD, NekDouble> &numfluxF,
    Array<OneD, NekDouble> &numfluxB)
{
    boost::ignore_unused(index);

    int nvariables = 3;

    NekDouble g = m_g;
    NekDouble uRF, vRF, uLB, vLB;
    NekDouble lambdaF, lambdaB;

    Array<OneD, NekDouble> EigF(nvariables);
    Array<OneD, NekDouble> EigB(nvariables);

    // Compute Magnitude and Dot product of moving frames for the index
    // ComputeMagAndDot(index, MageF1, MageF2, MageB1, MageB2, eF1_cdot_eB1,
    //                  eF1_cdot_eB2, eF2_cdot_eB1, eF2_cdot_eB2);

    // Get the velocity in the normal to the edge
    // velL = uL * m_ncdotMFFwd[0][index] + vL * m_ncdotMFFwd[1][index];
    // velR = -1.0 * (uR * m_ncdotMFBwd[0][index] + vR *
    // m_ncdotMFBwd[1][index]);

    EigF[0] = velL - sqrt(g * hL);
    EigF[1] = velL;
    EigF[2] = velL + sqrt(g * hL);

    EigB[0] = velR - sqrt(g * hR);
    EigB[1] = velR;
    EigB[2] = velR + sqrt(g * hR);

    lambdaF = Vmath::Vamax(nvariables, EigF, 1);
    lambdaB = Vmath::Vamax(nvariables, EigB, 1);

    // uRF = uR component in moving frames e^{Fwd}
    // vRF = vR component in moving frames e^{Fwd}
    uRF = (uR * eF1_cdot_eB1 + vR * eF1_cdot_eB2) / MageF1;
    vRF = (uR * eF2_cdot_eB1 + vR * eF2_cdot_eB2) / MageF2;

    numfluxF[0] = 0.5 * (hL * uL + hR * uRF);
    numfluxF[1] = 0.5 * (hL * vL + hR * vRF);
    numfluxF[2] = 0.5 * lambdaF * (hL - hR);

    numfluxF[3] = 0.5 * (hL * uL * uL * MageF1 + hR * uRF * uRF * MageB1 +
                         0.5 * g * (hL * hL + hR * hR));
    numfluxF[4] = 0.5 * (hL * uL * vL * MageF1 + hR * uRF * vRF * MageB1);
    numfluxF[5] = 0.5 * lambdaF * (uL * hL - uRF * hR);

    numfluxF[6] = 0.5 * (hL * uL * vL * MageF2 + hR * uRF * vRF * MageB2);
    numfluxF[7] = 0.5 * (hL * vL * vL * MageF2 + hR * vRF * vRF * MageB2 +
                         0.5 * g * (hL * hL + hR * hR));
    numfluxF[8] = 0.5 * lambdaF * (vL * hL - vRF * hR);

    // uLB = uL component in moving frames e^{Bwd}
    // vLB = vL component in moving frames e^{Bwd}
    uLB = (uL * eF1_cdot_eB1 + vL * eF2_cdot_eB1) / MageB1;
    vLB = (uL * eF1_cdot_eB2 + vL * eF2_cdot_eB2) / MageB2;

    numfluxB[0] = 0.5 * (hR * uR + hR * uLB);
    numfluxB[1] = 0.5 * (hR * vR + hR * vLB);
    numfluxB[2] = 0.5 * lambdaB * (hL - hR);

    numfluxB[3] = 0.5 * (hR * uR * uR * MageB1 + hR * uLB * uLB * MageF1 +
                         0.5 * g * (hR * hR + hL * hL));
    numfluxB[4] = 0.5 * (hR * uR * vR * MageB1 + hR * uLB * vLB * MageF1);
    numfluxB[5] = 0.5 * lambdaB * (uLB * hL - uR * hR);

    numfluxB[6] = 0.5 * (hR * uR * vR * MageB2 + hR * uLB * vLB * MageF2);
    numfluxB[7] = 0.5 * (hR * vR * vR * MageB2 + hR * vLB * vLB * MageF2 +
                         0.5 * g * (hR * hR + hL * hL));
    numfluxB[8] = 0.5 * lambdaB * (vLB * hL - vR * hR);
}

void MMFSWE::RusanovFlux(const int index, NekDouble hL, NekDouble uL,
                         NekDouble vL, NekDouble hR, NekDouble uR, NekDouble vR,
                         Array<OneD, NekDouble> &numfluxF,
                         Array<OneD, NekDouble> &numfluxB)
{
    int nvariables = 3;
    NekDouble MageF1, MageF2, MageB1, MageB2;
    NekDouble eF1_cdot_eB1, eF1_cdot_eB2;
    NekDouble eF2_cdot_eB1, eF2_cdot_eB2;

    NekDouble g = m_g;
    NekDouble uRF, vRF, uLB, vLB;
    NekDouble velL, velR;

    Array<OneD, NekDouble> EigF(nvariables);
    Array<OneD, NekDouble> EigB(nvariables);

    // Compute Magnitude and Dot product of moving frames for the index
    ComputeMagAndDot(index, MageF1, MageF2, MageB1, MageB2, eF1_cdot_eB1,
                     eF1_cdot_eB2, eF2_cdot_eB1, eF2_cdot_eB2);

    // Get the velocity in the normal to the edge
    velL = uL * m_ncdotMFFwd[0][index] + vL * m_ncdotMFFwd[1][index];
    velR = -1.0 * (uR * m_ncdotMFBwd[0][index] + vR * m_ncdotMFBwd[1][index]);

    NekDouble SL, SR;
    SL = fabs(velL) + sqrt(g * hL);
    SR = fabs(velR) + sqrt(g * hR);

    NekDouble S;
    if (SL > SR)
        S = SL;
    else
        S = SR;

    // uRF = uR component in moving frames e^{Fwd}
    // vRF = vR component in moving frames e^{Fwd}
    uRF = (uR * eF1_cdot_eB1 + vR * eF1_cdot_eB2) / MageF1;
    vRF = (uR * eF2_cdot_eB1 + vR * eF2_cdot_eB2) / MageF2;

    numfluxF[0] = 0.5 * (hL * uL + hR * uRF);
    numfluxF[1] = 0.5 * (hL * vL + hR * vRF);
    numfluxF[2] = 0.5 * S * (hL - hR);

    numfluxF[3] =
        0.5 * (hL * uL * uL + hR * uRF * uRF + 0.5 * g * (hL * hL + hR * hR));
    numfluxF[4] = 0.5 * (hL * uL * vL + hR * uRF * vRF);
    numfluxF[5] = 0.5 * S * (uL * hL - uRF * hR);

    numfluxF[6] = 0.5 * (hL * uL * vL + hR * uRF * vRF);
    numfluxF[7] =
        0.5 * (hL * vL * vL + hR * vRF * vRF + 0.5 * g * (hL * hL + hR * hR));
    numfluxF[8] = 0.5 * S * (vL * hL - vRF * hR);

    // uLB = uL component in moving frames e^{Bwd}
    // vLB = vL component in moving frames e^{Bwd}
    uLB = (uL * eF1_cdot_eB1 + vL * eF2_cdot_eB1) / MageB1;
    vLB = (uL * eF1_cdot_eB2 + vL * eF2_cdot_eB2) / MageB2;

    numfluxB[0] = 0.5 * (hR * uR + hR * uLB);
    numfluxB[1] = 0.5 * (hR * vR + hR * vLB);
    numfluxB[2] = 0.5 * S * (hL - hR);

    numfluxB[3] =
        0.5 * (hR * uR * uR + hR * uLB * uLB + 0.5 * g * (hR * hR + hL * hL));
    numfluxB[4] = 0.5 * (hR * uR * vR + hR * uLB * vLB);
    numfluxB[5] = 0.5 * S * (uLB * hL - uR * hR);

    numfluxB[6] = 0.5 * (hR * uR * vR + hR * uLB * vLB);
    numfluxB[7] =
        0.5 * (hR * vR * vR + hR * vLB * vLB + 0.5 * g * (hR * hR + hL * hL));
    numfluxB[8] = 0.5 * S * (vLB * hL - vR * hR);
}

void MMFSWE::ComputeMagAndDot(const int index, NekDouble &MageF1,
                              NekDouble &MageF2, NekDouble &MageB1,
                              NekDouble &MageB2, NekDouble &eF1_cdot_eB1,
                              NekDouble &eF1_cdot_eB2, NekDouble &eF2_cdot_eB1,
                              NekDouble &eF2_cdot_eB2)
{
    NekDouble MF1x, MF1y, MF1z, MF2x, MF2y, MF2z;
    NekDouble MB1x, MB1y, MB1z, MB2x, MB2y, MB2z;

    MF1x = m_MFtraceFwd[0][0][index];
    MF1y = m_MFtraceFwd[0][1][index];
    MF1z = m_MFtraceFwd[0][2][index];

    MF2x = m_MFtraceFwd[1][0][index];
    MF2y = m_MFtraceFwd[1][1][index];
    MF2z = m_MFtraceFwd[1][2][index];

    MB1x = m_MFtraceBwd[0][0][index];
    MB1y = m_MFtraceBwd[0][1][index];
    MB1z = m_MFtraceBwd[0][2][index];

    MB2x = m_MFtraceBwd[1][0][index];
    MB2y = m_MFtraceBwd[1][1][index];
    MB2z = m_MFtraceBwd[1][2][index];

    // MFtrace = MFtrace [ j*spacedim + k ], j = shape, k = sapce
    MageF1 = MF1x * MF1x + MF1y * MF1y + MF1z * MF1z;
    MageF2 = MF2x * MF2x + MF2y * MF2y + MF2z * MF2z;
    MageB1 = MB1x * MB1x + MB1y * MB1y + MB1z * MB1z;
    MageB2 = MB2x * MB2x + MB2y * MB2y + MB2z * MB2z;

    eF1_cdot_eB1 = MF1x * MB1x + MF1y * MB1y + MF1z * MB1z;
    eF1_cdot_eB2 = MF1x * MB2x + MF1y * MB2y + MF1z * MB2z;
    eF2_cdot_eB1 = MF2x * MB1x + MF2y * MB1y + MF2z * MB1z;
    eF2_cdot_eB2 = MF2x * MB2x + MF2y * MB2y + MF2z * MB2z;
}

// Add Coriolis factors
void MMFSWE::AddCoriolis(Array<OneD, Array<OneD, NekDouble>> &physarray,
                         Array<OneD, Array<OneD, NekDouble>> &outarray)
{
    int ncoeffs = outarray[0].size();
    int nq      = physarray[0].size();

    Array<OneD, NekDouble> tmp(nq);
    Array<OneD, NekDouble> tmpc(ncoeffs);

    // physarray is primitive
    // conservative formulation compute h
    // h = \eta + d
    Array<OneD, NekDouble> h(nq);
    Vmath::Vadd(nq, physarray[0], 1, m_depth, 1, h, 1);

    // j = 0;
    Vmath::Vmul(nq, m_coriolis, 1, physarray[2], 1, tmp, 1);
    Vmath::Vmul(nq, h, 1, tmp, 1, tmp, 1);
    m_fields[0]->IProductWRTBase(tmp, tmpc);
    Vmath::Vadd(ncoeffs, tmpc, 1, outarray[1], 1, outarray[1], 1);

    // j = 1;
    Vmath::Vmul(nq, m_coriolis, 1, physarray[1], 1, tmp, 1);
    Vmath::Vmul(nq, h, 1, tmp, 1, tmp, 1);
    Vmath::Neg(nq, tmp, 1);
    m_fields[0]->IProductWRTBase(tmp, tmpc);
    Vmath::Vadd(ncoeffs, tmpc, 1, outarray[2], 1, outarray[2], 1);
}

// void MMFSWE::AddCoriolisSDR(Array<OneD, Array<OneD, NekDouble>> &physarray,
//                             Array<OneD, Array<OneD, NekDouble>> &outarray)
// {
//     int ncoeffs = outarray[0].size();
//     int nq = physarray[0].size();

//     Array<OneD, NekDouble> tmp(nq);
//     Array<OneD, NekDouble> tmpc(ncoeffs);

//     // physarray is primitive
//     // conservative formulation compute h
//     // h = \eta + d
//     Array<OneD, NekDouble> h(nq);
//     Vmath::Vadd(nq, physarray[0], 1, m_depth, 1, h, 1);

//     Array<OneD, NekDouble> newu(nq);
//     Array<OneD, NekDouble> newv(nq);

//     Convert_LOCAL_TO_LOCALSPHERE(physarray[1], physarray[2], newu, newv);

//     // newu = u (e^1 \cdot e^1_new) + v (e^2 \cdot e^1_new)
//     // Vmath::Vvtvp(nq, &physu[0], 1, &m_LOCAL_cdot_LOCALSPHERE[0][0][0], 1,
//     &newu[0], 1, &newu[0], 1);
//     // Vmath::Vvtvp(nq, &physv[0], 1, &m_LOCAL_cdot_LOCALSPHERE[1][0][0], 1,
//     &newu[0], 1, &newu[0], 1);

//     // j = 0;
//     Vmath::Vmul(nq, newu, 1, m_LOCAL_cdot_CoordAxisDivMF[0][1], 1, tmp, 1);
//     Vmath::Neg(nq, tmp, 1);
//     Vmath::Vvtvp(nq, newv, 1, m_LOCAL_cdot_CoordAxisDivMF[0][0], 1, tmp, 1,
//     tmp, 1);

//     Vmath::Vmul(nq, m_coriolis, 1, tmp, 1, tmp, 1);
//     Vmath::Vmul(nq, h, 1, tmp, 1, tmp, 1);
//     m_fields[0]->IProductWRTBase(tmp, tmpc);
//     Vmath::Vadd(ncoeffs, tmpc, 1, outarray[1], 1, outarray[1], 1);

//     // j = 1;
//     Vmath::Vmul(nq, newu, 1, m_LOCAL_cdot_CoordAxisDivMF[1][1], 1, tmp, 1);
//     Vmath::Neg(nq, tmp, 1);
//     Vmath::Vvtvp(nq, newv, 1, m_LOCAL_cdot_CoordAxisDivMF[1][0], 1, tmp, 1,
//     tmp, 1);

//     Vmath::Vmul(nq, m_coriolis, 1, physarray[1], 1, tmp, 1);
//     Vmath::Vmul(nq, h, 1, tmp, 1, tmp, 1);
//     Vmath::Neg(nq, tmp, 1);
//     m_fields[0]->IProductWRTBase(tmp, tmpc);
//     Vmath::Vadd(ncoeffs, tmpc, 1, outarray[2], 1, outarray[2], 1);
// }
// int indx = 1;
// for (int j = 0; j < m_shapedim; ++j)
// {
//     if (j == 0)
//     {
//         indx = 2;
//     }

//     else if (j == 1)
//     {
//         indx = 1;
//     }

//     // add to hu equation
//     Vmath::Vmul(nq, m_coriolis, 1, physarray[indx], 1, tmp, 1);
//     Vmath::Vmul(nq, h, 1, tmp, 1, tmp, 1);

//     if (j == 1)
//     {
//         Vmath::Neg(nq, tmp, 1);
//     }

//     // N \cdot (e^1 \times e^2 )
//     // Vmath::Vmul(nq, &m_MF1crossMF2dotSN[0], 1, &tmp[0], 1, &tmp[0], 1);
//     m_fields[0]->IProductWRTBase(tmp, tmpc);
//     Vmath::Vadd(ncoeffs, tmpc, 1, outarray[j + 1], 1, outarray[j + 1], 1);
// }
// }

void MMFSWE::AddElevationEffect(Array<OneD, Array<OneD, NekDouble>> &physarray,
                                Array<OneD, Array<OneD, NekDouble>> &outarray)
{
    int ncoeffs = outarray[0].size();
    int nq      = physarray[0].size();

    Array<OneD, NekDouble> h(nq);
    Array<OneD, NekDouble> tmp(nq);
    Array<OneD, NekDouble> tmpc(ncoeffs);

    // physarray is primitive
    // conservative formulation compute h
    // h = \eta + d
    Vmath::Vadd(nq, physarray[0], 1, m_depth, 1, h, 1);

    for (int j = 0; j < m_shapedim; ++j)
    {
        MMFDirectionalDeriv(m_movingframes[j], m_depth, tmp);
        Vmath::Vmul(nq, h, 1, tmp, 1, tmp, 1);
        Vmath::Smul(nq, m_g, tmp, 1, tmp, 1);

        m_fields[0]->IProductWRTBase(tmp, tmpc);
        Vmath::Vadd(ncoeffs, tmpc, 1, outarray[j + 1], 1, outarray[j + 1], 1);
    }
}

void MMFSWE::AddElevationEffectPhys(
    Array<OneD, Array<OneD, NekDouble>> &physarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray)
{
    int ncoeffs = outarray[0].size();
    int nq      = physarray[0].size();

    Array<OneD, NekDouble> h(nq);
    Array<OneD, NekDouble> tmp(nq);
    Array<OneD, NekDouble> tmpc(ncoeffs);

    // physarray is primitive
    // conservative formulation compute h
    // h = \eta + d
    Vmath::Vadd(nq, physarray[0], 1, m_depth, 1, h, 1);

    for (int j = 0; j < m_shapedim; ++j)
    {
        MMFDirectionalDeriv(m_movingframes[j], m_depth, tmp);
        Vmath::Vmul(nq, h, 1, tmp, 1, tmp, 1);
        Vmath::Smul(nq, m_g, tmp, 1, tmp, 1);

        Vmath::Vadd(nq, tmp, 1, outarray[j + 1], 1, outarray[j + 1], 1);

        // m_fields[0]->IProductWRTBase(tmp, tmpc);
        // Vmath::Vadd(ncoeffs, tmpc, 1, outarray[j + 1], 1, outarray[j + 1],
        // 1);
    }
}

// =================================================
// Add rotational factors
// =================================================
void MMFSWE::AddTimeVariantFrames(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const int DivergenceRestore)
{
    // routine works for both primitive and conservative formulations
    int ncoeffs = outarray[0].size();
    int nq      = m_fields[0]->GetTotPoints();
    int nvar    = 3; // only the dependent variables

    // Compute h
    Array<OneD, NekDouble> h(nq);
    Vmath::Vadd(nq, &inarray[0][0], 1, &m_depth[0], 1, &h[0], 1);

    Array<OneD, NekDouble> tmp(nq);

    Array<OneD, Array<OneD, NekDouble>> physfield(nvar);
    for (int i = 0; i < nvar; ++i)
    {
        physfield[i] = Array<OneD, NekDouble>(nq, 0.0);
    }

    // if m_DivergenceRestore>=2, we use a sphericl coordinate expansion
    ComputephysfieldinMF(inarray, physfield, DivergenceRestore);

    Array<OneD, NekDouble> de0dt_cdot_e0;
    Array<OneD, NekDouble> de0dt_cdot_e1;
    Array<OneD, NekDouble> de1dt_cdot_e0;
    Array<OneD, NekDouble> de1dt_cdot_e1;
    if (m_DivergenceRestore >= 2)
    {
        Compute_demdt_cdot_ek(0, 0, physfield, mf_LOCSPH, de0dt_cdot_e0);
        Compute_demdt_cdot_ek(1, 0, physfield, mf_LOCSPH, de1dt_cdot_e0);
        Compute_demdt_cdot_ek(0, 1, physfield, mf_LOCSPH, de0dt_cdot_e1);
        Compute_demdt_cdot_ek(1, 1, physfield, mf_LOCSPH, de1dt_cdot_e1);
    }

    else
    {
        Compute_demdt_cdot_ek(0, 0, physfield, movingframes, de0dt_cdot_e0);
        Compute_demdt_cdot_ek(1, 0, physfield, movingframes, de1dt_cdot_e0);
        Compute_demdt_cdot_ek(0, 1, physfield, movingframes, de0dt_cdot_e1);
        Compute_demdt_cdot_ek(1, 1, physfield, movingframes, de1dt_cdot_e1);
    }

    // Rott1 = ( u (de0/dt) + v(de1/dt) ) cdot e0
    // Rott2 = ( u (de0/dt) + v(de1/dt) ) cdot e1
    Array<OneD, NekDouble> Rott1(nq);
    Array<OneD, NekDouble> Rott2(nq);
    Vmath::Vmul(nq, inarray[1], 1, de0dt_cdot_e0, 1, Rott1, 1);
    Vmath::Vmul(nq, inarray[1], 1, de0dt_cdot_e1, 1, Rott2, 1);
    Vmath::Vvtvp(nq, inarray[2], 1, de1dt_cdot_e0, 1, Rott1, 1, Rott1, 1);
    Vmath::Vvtvp(nq, inarray[2], 1, de1dt_cdot_e1, 1, Rott2, 1, Rott2, 1);

    // Multiply H and \partial \phi / \partial t which is assumed to be u_{\phi}
    Vmath::Vmul(nq, &h[0], 1, &Rott1[0], 1, &Rott1[0], 1);
    Vmath::Vmul(nq, &h[0], 1, &Rott2[0], 1, &Rott2[0], 1);

    Vmath::Neg(nq, Rott1, 1);
    Vmath::Neg(nq, Rott2, 1);

    Array<OneD, NekDouble> tmpc1(ncoeffs);
    Array<OneD, NekDouble> tmpc2(ncoeffs);
    m_fields[0]->IProductWRTBase(Rott1, tmpc1);
    m_fields[0]->IProductWRTBase(Rott2, tmpc2);

    Vmath::Vadd(ncoeffs, tmpc1, 1, outarray[1], 1, outarray[1], 1);
    Vmath::Vadd(ncoeffs, tmpc2, 1, outarray[2], 1, outarray[2], 1);
}

// =====================================================================
// Compute \frac{d e^{indm}}{dt} \cdot e^{indk}  / \| e^{indk} \|
// Equivalently, [ \frac{d e^{indm}}{d \xi_1} \frac{ d \xi_1}{d \phi} + \frac{d
// e^{indm}}{d \xi_2} \frac{ d \xi_2}{d \phi} ] \frac{d \phi}{dt}
void MMFSWE::Compute_demdt_cdot_ek(
    const int indm, const int indk,
    const Array<OneD, const Array<OneD, NekDouble>> &physarray,
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    Array<OneD, NekDouble> &outarray)
{
    int nq = m_fields[0]->GetNpoints();

    Array<OneD, NekDouble> tmp(nq);

    outarray = Array<OneD, NekDouble>(nq, 0.0);
    for (int j = 0; j < m_shapedim; ++j)
    {
        for (int k = 0; k < m_spacedim; ++k)
        {
            // Compute d e^m / d \xi_1 and d e^m / d \xi_2
            Vmath::Vcopy(nq, &m_movingframes[indm][k * nq], 1, &tmp[0], 1);
            MMFDirectionalDeriv(movingframes[j], tmp, tmp);

            Vmath::Vmul(nq, &physarray[j + 1][0], 1, &tmp[0], 1, &tmp[0], 1);

            Vmath::Vvtvp(nq, &tmp[0], 1, &movingframes[indk][k * nq], 1,
                         &outarray[0], 1, &outarray[0], 1);
        }
    }
}

/**
 * @brief Compute the projection for the linear advection equation.
 *
 * @param inarray    Given fields.
 * @param outarray   Calculated solution.
 * @param time       Time.
 */
void MMFSWE::DoOdeProjection(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time)
{
    switch (m_projectionType)
    {
        case MultiRegions::eDiscontinuous:
        {
            ConservativeToPrimitive(inarray, outarray);
            SetBoundaryConditions(outarray, time);
            PrimitiveToConservative(outarray, outarray);
        }
        break;
        default:
            ASSERTL0(false, "Unknown projection scheme");
            break;
    }
}

void MMFSWE::SetBoundaryConditions(Array<OneD, Array<OneD, NekDouble>> &inarray,
                                   NekDouble time)
{

    int nvariables = m_fields.size();
    int cnt        = 0;

    // loop over Boundary Regions
    for (int n = 0; n < m_fields[0]->GetBndConditions().size(); ++n)
    {

        // Zonal Boundary Condition
        if (m_fields[0]->GetBndConditions()[n]->GetUserDefined() == "eMG")
        {
            if (m_expdim == 1)
            {
                ASSERTL0(false, "Illegal dimension");
            }
            else if (m_expdim == 2)
            {
                // ZonalBoundary2D(n,cnt,inarray);
            }
        }

        // Wall Boundary Condition
        if (m_fields[0]->GetBndConditions()[n]->GetUserDefined() == "eWall")
        {
            if (m_expdim == 1)
            {
                ASSERTL0(false, "Illegal dimension");
            }
            else if (m_expdim == 2)
            {
                WallBoundary2D(n, cnt, inarray);
            }
        }

        // Time Dependent Boundary Condition (specified in meshfile)
        if (m_fields[0]->GetBndConditions()[n]->GetUserDefined() ==
            "eTimeDependent")
        {
            for (int i = 0; i < nvariables; ++i)
            {
                m_fields[i]->EvaluateBoundaryConditions(time);
            }
        }
        cnt += m_fields[0]->GetBndCondExpansions()[n]->GetExpSize();
    }
}

void MMFSWE::WallBoundary2D(int bcRegion, int cnt,
                            Array<OneD, Array<OneD, NekDouble>> &physarray)
{

    int i;
    int nTraceNumPoints = GetTraceTotPoints();
    int nvariables      = physarray.size();

    // get physical values of the forward trace
    Array<OneD, Array<OneD, NekDouble>> Fwd0(nvariables);
    for (i = 0; i < nvariables; ++i)
    {
        Fwd0[i] = Array<OneD, NekDouble>(nTraceNumPoints);
        m_fields[i]->ExtractTracePhys(physarray[i], Fwd0[i]);
    }

    Array<OneD, Array<OneD, NekDouble>> Fwd(nvariables);
    for (i = 0; i < nvariables; ++i)
    {
        Fwd[i] = Array<OneD, NekDouble>(nTraceNumPoints);
        m_fields[i]->ExtractTracePhys(physarray[i], Fwd[i]);
    }

    // Adjust the physical values of the trace to take
    // user defined boundaries into account
    int e, id1, id2, npts;

    for (e = 0; e < m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExpSize();
         ++e)
    {
        npts = m_fields[0]
                   ->GetBndCondExpansions()[bcRegion]
                   ->GetExp(e)
                   ->GetTotPoints();
        id1 = m_fields[0]->GetBndCondExpansions()[bcRegion]->GetPhys_Offset(e);
        id2 = m_fields[0]->GetTrace()->GetPhys_Offset(
            m_fields[0]->GetTraceMap()->GetBndCondIDToGlobalTraceID(cnt + e));

        switch (m_expdim)
        {
            case 1:
            {
                // negate the forward flux
                Vmath::Neg(npts, &Fwd[1][id2], 1);
            }
            break;
            case 2:
            {
                Array<OneD, NekDouble> tmp_n(npts);
                Array<OneD, NekDouble> tmp_t(npts);

                Vmath::Vmul(npts, &Fwd[1][id2], 1, &m_ncdotMFFwd[0][id2], 1,
                            &tmp_n[0], 1);
                Vmath::Vvtvp(npts, &Fwd[2][id2], 1, &m_ncdotMFFwd[1][id2], 1,
                             &tmp_n[0], 1, &tmp_n[0], 1);

                Vmath::Vmul(npts, &Fwd[1][id2], 1, &m_nperpcdotMFFwd[0][id2], 1,
                            &tmp_t[0], 1);
                Vmath::Vvtvp(npts, &Fwd[2][id2], 1, &m_nperpcdotMFFwd[1][id2],
                             1, &tmp_t[0], 1, &tmp_t[0], 1);

                // negate the normal flux
                Vmath::Neg(npts, tmp_n, 1);

                Array<OneD, NekDouble> denom(npts);
                Array<OneD, NekDouble> tmp_u(npts);
                Array<OneD, NekDouble> tmp_v(npts);

                // denom = (e^1 \cdot n ) (e^2 \cdot t) - (e^2 \cdot n ) (e^1
                // \cdot t)
                Vmath::Vmul(npts, &m_ncdotMFFwd[1][id2], 1,
                            &m_nperpcdotMFFwd[0][id2], 1, &denom[0], 1);
                Vmath::Vvtvm(npts, &m_ncdotMFFwd[0][id2], 1,
                             &m_nperpcdotMFFwd[1][id2], 1, &denom[0], 1,
                             &denom[0], 1);

                Vmath::Vmul(npts, &m_ncdotMFFwd[1][id2], 1, &tmp_t[0], 1,
                            &tmp_u[0], 1);
                Vmath::Vvtvm(npts, &m_nperpcdotMFFwd[1][id2], 1, &tmp_n[0], 1,
                             &tmp_u[0], 1, &tmp_u[0], 1);
                Vmath::Vdiv(npts, &tmp_u[0], 1, &denom[0], 1, &tmp_u[0], 1);

                Vmath::Vcopy(npts, &tmp_u[0], 1, &Fwd[1][id2], 1);

                Vmath::Vmul(npts, &m_nperpcdotMFFwd[0][id2], 1, &tmp_n[0], 1,
                            &tmp_v[0], 1);
                Vmath::Vvtvm(npts, &m_ncdotMFFwd[0][id2], 1, &tmp_t[0], 1,
                             &tmp_v[0], 1, &tmp_v[0], 1);
                Vmath::Vdiv(npts, &tmp_v[0], 1, &denom[0], 1, &tmp_v[0], 1);

                Vmath::Vcopy(npts, &tmp_v[0], 1, &Fwd[2][id2], 1);
            }
            break;

            default:
                ASSERTL0(false, "Illegal expansion dimension");
        }

        // copy boundary adjusted values into the boundary expansion
        for (i = 0; i < nvariables; ++i)
        {
            Vmath::Vcopy(npts, &Fwd[i][id2], 1,
                         &(m_fields[i]
                               ->GetBndCondExpansions()[bcRegion]
                               ->UpdatePhys())[id1],
                         1);
        }
    }
}

/**
 *
 */
void MMFSWE::v_DoInitialise()
{
    // Compute m_depth and m_Derivdepth
    EvaluateWaterDepth();
    EvaluateCoriolis();
    SetInitialConditions();
    PrimitiveToConservative();

    // transfer the initial conditions to modal values
    for (int i = 0; i < m_fields.size(); ++i)
    {
        m_fields[i]->SetPhysState(true);
        m_fields[i]->FwdTrans(m_fields[i]->GetPhys(),
                              m_fields[i]->UpdateCoeffs());
    }
}

void MMFSWE::EvaluateWaterDepth(void)
{
    int nq = GetTotPoints();

    m_depth = Array<OneD, NekDouble>(nq, 0.0);

    switch (m_TestType)
    {
        case eTestPlane:
        {
            Vmath::Fill(nq, 1.0, m_depth, 1);
        }
        break;

        case eTestUnsteadyZonal:
        {
            // H_0 = k_1 - k_2 - (1/2/g) (Omega sin \phi)
            Array<OneD, NekDouble> x(nq);
            Array<OneD, NekDouble> y(nq);
            Array<OneD, NekDouble> z(nq);

            m_fields[0]->GetCoords(x, y, z);

            NekDouble x0j, x1j, x2j, rad;
            NekDouble sin_varphi, cos_varphi, sin_theta, cos_theta;

            for (int j = 0; j < nq; ++j)
            {
                x0j = x[j];
                x1j = y[j];
                x2j = z[j];

                CartesianToElliptical(x0j, x1j, x2j, rad, sin_varphi,
                                      cos_varphi, sin_theta, cos_theta);

                m_depth[j] =
                    m_H0 - m_k2 -
                    (0.5 / m_g) * (m_Omega * cos_theta) * (m_Omega * cos_theta);
            }
        }
        break;

        case eTestIsolatedMountain:
        {
            // H_0 = k_1 - k_2 - (1/2/g) (Omega sin \phi)
            Array<OneD, NekDouble> x(nq);
            Array<OneD, NekDouble> y(nq);
            Array<OneD, NekDouble> z(nq);

            m_fields[0]->GetCoords(x, y, z);

            int indx = 0;
            NekDouble x0j, x1j, x2j, rad, dist2;
            NekDouble phi, theta, sin_varphi, cos_varphi, sin_theta, cos_theta;
            NekDouble hRad, phic, thetac;
            NekDouble Tol = 0.000001;

            hRad   = m_pi / 9.0;
            phic   = -m_pi / 2.0;
            thetac = m_pi / 6.0;

            for (int j = 0; j < nq; ++j)
            {
                x0j = x[j];
                x1j = y[j];
                x2j = z[j];

                CartesianToElliptical(x0j, x1j, x2j, rad, sin_varphi,
                                      cos_varphi, sin_theta, cos_theta);

                if ((std::abs(sin(phic) - sin_varphi) +
                     std::abs(sin(thetac) - cos_theta)) < Tol)
                {
                    std::cout << "A point " << j
                              << " is coincient with the singularity "
                              << std::endl;
                    indx = 1;
                }

                phi   = atan2(sin_varphi, cos_varphi);
                theta = atan2(cos_theta, sin_theta);

                // Compute r
                dist2 = (phi - phic) * (phi - phic) +
                        (theta - thetac) * (theta - thetac);

                if (dist2 > hRad * hRad)
                {
                    dist2 = hRad * hRad;
                }

                m_depth[j] = m_H0 - m_hs0 * (1.0 - sqrt(dist2) / hRad);
            }

            if (!indx)
            {
                std::cout << "No point is coincident with the singularity point"
                          << std::endl;
            }
        }
        break;

        case eTestUnstableJet:
        {
            for (int j = 0; j < nq; ++j)
            {
                m_depth[j] = m_H0;
            }
        }
        break;

        case eTestSteadyZonal:
        case eTestRossbyWave:
        {
            Vmath::Zero(nq, m_depth, 1);
        }
        break;

        default:
        {
            Vmath::Zero(nq, m_depth, 1);
        }
        break;
    }

    // Comptue \nabla m_depth \cdot e^i
    m_Derivdepth = Array<OneD, Array<OneD, NekDouble>>(m_shapedim);

    for (int j = 0; j < m_shapedim; j++)
    {
        m_Derivdepth[j] = Array<OneD, NekDouble>(nq);
        MMFDirectionalDeriv(m_movingframes[j], m_depth,
                                          m_Derivdepth[j]);
    }

    std::cout << "Water Depth (m_depth) was generated with mag = "
              << AvgAbsInt(m_depth) << " with max. deriv = ( "
              << Vmath::Vamax(nq, m_Derivdepth[0], 1) << " , "
              << Vmath::Vamax(nq, m_Derivdepth[1], 1) << " ) " << std::endl;
}

void MMFSWE::EvaluateCoriolis(void)
{
    switch (m_TestType)
    {
        case eTestPlane:
        {
            GetFunction("Coriolis")->Evaluate("f", m_coriolis);
        }
        break;

        case eTestSteadyZonal:
        {
            EvaluateCoriolisForZonalFlow(m_coriolis);
        }
        break;

        case eTestUnsteadyZonal:
        case eTestIsolatedMountain:
        case eTestUnstableJet:
        case eTestRossbyWave:
        {
            EvaluateStandardCoriolis(m_coriolis);
        }
        break;

        default:
            break;
    }
}

void MMFSWE::EvaluateCoriolisForZonalFlow(Array<OneD, NekDouble> &outarray)
{
    int nq = GetTotPoints();
    NekDouble sin_theta, cos_theta, sin_varphi, cos_varphi;

    Array<OneD, NekDouble> x(nq);
    Array<OneD, NekDouble> y(nq);
    Array<OneD, NekDouble> z(nq);

    m_fields[0]->GetCoords(x, y, z);

    NekDouble x0j, x1j, x2j, rad;

    outarray = Array<OneD, NekDouble>(nq, 0.0);
    for (int j = 0; j < nq; ++j)
    {
        x0j = x[j];
        x1j = y[j];
        x2j = z[j];

        CartesianToElliptical(x0j, x1j, x2j, rad, sin_varphi, cos_varphi,
                              sin_theta, cos_theta);

        // H = 2 \Omega *(- \cos \phi \cos \theta \sin \alpha + \sin \theta \cos
        // \alpha )
        outarray[j] =
            2.0 * m_Omega *
            (-cos_varphi * sin_theta * sin(m_alpha) + cos_theta * cos(m_alpha));
    }
}

void MMFSWE::EvaluateStandardCoriolis(Array<OneD, NekDouble> &outarray)
{
    int nq = GetTotPoints();

    NekDouble x0j, x1j, x2j, rad;
    NekDouble sin_theta, cos_theta, sin_varphi, cos_varphi;

    Array<OneD, NekDouble> x(nq);
    Array<OneD, NekDouble> y(nq);
    Array<OneD, NekDouble> z(nq);

    m_fields[0]->GetCoords(x, y, z);

    outarray = Array<OneD, NekDouble>(nq, 0.0);
    for (int j = 0; j < nq; ++j)
    {
        x0j = x[j];
        x1j = y[j];
        x2j = z[j];

        CartesianToElliptical(x0j, x1j, x2j, rad, sin_varphi, cos_varphi,
                              sin_theta, cos_theta);

        outarray[j] = 2.0 * m_Omega * cos_theta;
    }
}

void MMFSWE::v_SetInitialConditions(const NekDouble initialtime,
                                    bool dumpInitialConditions,
                                    const int domain)
{
    boost::ignore_unused(domain);

    int nq = GetTotPoints();

    switch (m_TestType)
    {
        case eTestPlane:
        {
            Array<OneD, NekDouble> eta0(nq);
            Array<OneD, NekDouble> u0(nq);
            Array<OneD, NekDouble> v0(nq);

            TestSWE2Dproblem(initialtime, 0, eta0);
            m_fields[0]->SetPhys(eta0);

            TestSWE2Dproblem(initialtime, 1, u0);
            m_fields[1]->SetPhys(u0);

            TestSWE2Dproblem(initialtime, 2, v0);
            m_fields[2]->SetPhys(v0);

            // forward transform to fill the modal coeffs
            for (int i = 0; i < m_fields.size(); ++i)
            {
                m_fields[i]->SetPhysState(true);
                m_fields[i]->FwdTrans(m_fields[i]->GetPhys(),
                                      m_fields[i]->UpdateCoeffs());
            }
        }
        break;

        case eTestSteadyZonal:
        {
            Array<OneD, NekDouble> eta0(nq);
            Array<OneD, NekDouble> u0(nq);
            Array<OneD, NekDouble> v0(nq);
            Array<OneD, NekDouble> zeta0(nq);

            SteadyZonalFlow(0, eta0);
            m_fields[0]->SetPhys(eta0);

            SteadyZonalFlow(1, u0);
            m_fields[1]->SetPhys(u0);

            SteadyZonalFlow(2, v0);
            m_fields[2]->SetPhys(v0);

            ComputeVorticity(u0, v0, zeta0);
            m_Vorticity0 = m_fields[0]->PhysIntegral(zeta0);

            // forward transform to fill the modal coeffs
            for (int i = 0; i < m_fields.size(); ++i)
            {
                m_fields[i]->SetPhysState(true);
                m_fields[i]->FwdTrans(m_fields[i]->GetPhys(),
                                      m_fields[i]->UpdateCoeffs());
            }
        }
        break;

        case eTestUnsteadyZonal:
        {
            Array<OneD, NekDouble> eta0(nq);
            Array<OneD, NekDouble> u0(nq);
            Array<OneD, NekDouble> v0(nq);

            UnsteadyZonalFlow(0, initialtime, eta0);
            m_fields[0]->SetPhys(eta0);

            UnsteadyZonalFlow(1, initialtime, u0);
            m_fields[1]->SetPhys(u0);

            UnsteadyZonalFlow(2, initialtime, v0);
            m_fields[2]->SetPhys(v0);

            // forward transform to fill the modal coeffs
            for (int i = 0; i < m_fields.size(); ++i)
            {
                m_fields[i]->SetPhysState(true);
                m_fields[i]->FwdTrans(m_fields[i]->GetPhys(),
                                      m_fields[i]->UpdateCoeffs());
            }
        }
        break;

        case eTestIsolatedMountain:
        {
            Array<OneD, NekDouble> eta0(nq);
            Array<OneD, NekDouble> u0(nq);
            Array<OneD, NekDouble> v0(nq);

            IsolatedMountainFlow(0, initialtime, eta0);
            m_fields[0]->SetPhys(eta0);

            IsolatedMountainFlow(1, initialtime, u0);
            m_fields[1]->SetPhys(u0);

            IsolatedMountainFlow(2, initialtime, v0);
            m_fields[2]->SetPhys(v0);

            // forward transform to fill the modal coeffs
            for (int i = 0; i < m_fields.size(); ++i)
            {
                m_fields[i]->SetPhysState(true);
                m_fields[i]->FwdTrans(m_fields[i]->GetPhys(),
                                      m_fields[i]->UpdateCoeffs());
            }
        }
        break;

        case eTestUnstableJet:
        {
            Array<OneD, NekDouble> eta0(nq);
            Array<OneD, NekDouble> u0(nq);
            Array<OneD, NekDouble> v0(nq);

            UnstableJetFlow(0, initialtime, eta0);
            m_fields[0]->SetPhys(eta0);

            UnstableJetFlow(1, initialtime, u0);
            m_fields[1]->SetPhys(u0);

            UnstableJetFlow(2, initialtime, v0);
            m_fields[2]->SetPhys(v0);

            // forward transform to fill the modal coeffs
            for (int i = 0; i < m_fields.size(); ++i)
            {
                m_fields[i]->SetPhysState(true);
                m_fields[i]->FwdTrans(m_fields[i]->GetPhys(),
                                      m_fields[i]->UpdateCoeffs());
            }
        }
        break;

        case eTestRossbyWave:
        {
            Array<OneD, NekDouble> eta0(nq);
            Array<OneD, NekDouble> u0(nq);
            Array<OneD, NekDouble> v0(nq);

            RossbyWave(0, eta0);
            m_fields[0]->SetPhys(eta0);

            RossbyWave(1, u0);
            m_fields[1]->SetPhys(u0);

            RossbyWave(2, v0);
            m_fields[2]->SetPhys(v0);

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

    // Projection on the domain
    for (int i = 0; i < m_fields.size(); ++i)
    {
        m_fields[i]->FwdTrans(m_fields[i]->GetPhys(),
                              m_fields[i]->UpdateCoeffs());
        m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(),
                              m_fields[i]->UpdatePhys());
    }

    // Set up initial Mass, Energy, and Enstrophy
    m_Mass0      = ComputeMass(m_fields[0]->GetPhys());
    m_Energy0    = ComputeEnergy(m_fields[0]->GetPhys(), m_fields[1]->GetPhys(),
                                 m_fields[2]->GetPhys());
    m_Enstrophy0 = ComputeEnstrophy(
        m_fields[0]->GetPhys(), m_fields[1]->GetPhys(), m_fields[2]->GetPhys());

    if (dumpInitialConditions)
    {
        // dump initial conditions to file
        std::string outname = m_sessionName + "_initial.chk";
        WriteFld(outname);

        outname = m_sessionName + "_initialCART.chk";
        Checkpoint_Output_Cartesian(outname);

        Array<OneD, Array<OneD, NekDouble>> velocity(m_spacedim);
        for (int i = 0; i < m_spacedim; ++i)
        {
            velocity[i] = Array<OneD, NekDouble>(nq);
        }

        for (int k = 0; k < nq; k++)
        {
            velocity[0][k] =
                (m_fields[1]->GetPhys())[k] * m_movingframes[0][k] +
                (m_fields[2]->GetPhys())[k] * m_movingframes[1][k];
            velocity[1][k] =
                (m_fields[1]->GetPhys())[k] * m_movingframes[0][k + nq] +
                (m_fields[2]->GetPhys())[k] * m_movingframes[1][k + nq];
            velocity[2][k] =
                (m_fields[1]->GetPhys())[k] * m_movingframes[0][k + 2 * nq] +
                (m_fields[2]->GetPhys())[k] * m_movingframes[1][k + 2 * nq];
        }

        CheckMeshErr(m_movingframes, velocity);
    }
}

void MMFSWE::TestSWE2Dproblem(const NekDouble time, unsigned int field,
                              Array<OneD, NekDouble> &outfield)
{

    boost::ignore_unused(time);

    int nq = m_fields[0]->GetNpoints();

    Array<OneD, NekDouble> x0(nq);
    Array<OneD, NekDouble> x1(nq);
    Array<OneD, NekDouble> x2(nq);

    m_fields[0]->GetCoords(x0, x1, x2);

    Array<OneD, NekDouble> eta0(nq);
    Array<OneD, NekDouble> u0(nq);
    Array<OneD, NekDouble> v0(nq);

    Array<OneD, Array<OneD, NekDouble>> uvec(m_spacedim);

    for (int i = 0; i < m_spacedim; ++i)
    {
        uvec[i] = Array<OneD, NekDouble>(nq, 0.0);
    }

    for (int i = 0; i < nq; ++i)
    {
        eta0[i] = (0.771 * 0.395 * 0.395 * (1.0 / cosh(0.395 * x0[i])) *
                   (1.0 / cosh(0.395 * x0[i]))) *
                  (3.0 + 6.0 * x1[i] * x1[i]) / (4.0) *
                  exp(-0.5 * x1[i] * x1[i]);
        uvec[0][i] = (0.771 * 0.395 * 0.395 * (1.0 / cosh(0.395 * x0[i])) *
                      (1.0 / cosh(0.395 * x0[i]))) *
                     (-9.0 + 6.0 * x1[i] * x1[i]) / (4.0) *
                     exp(-0.5 * x1[i] * x1[i]);
        uvec[1][i] = (-2.0 * 0.395 * tanh(0.395 * x0[i])) *
                     (0.771 * 0.395 * 0.395 * (1.0 / cosh(0.395 * x0[i])) *
                      (1.0 / cosh(0.395 * x0[i]))) *
                     (2.0 * x1[i]) * exp(-0.5 * x1[i] * x1[i]);
    }

    u0 = CartesianToMovingframes(m_movingframes, uvec, 0);
    v0 = CartesianToMovingframes(m_movingframes, uvec, 1);

    switch (field)
    {
        case (0):
        {
            outfield = eta0;
        }
        break;

        case (1):
        {
            outfield = u0;
        }
        break;

        case (2):
        {
            outfield = v0;
        }
        break;
    }
}

void MMFSWE::SteadyZonalFlow(unsigned int field,
                             Array<OneD, NekDouble> &outfield)
{
    int nq = GetTotPoints();
    NekDouble uhat, vhat;
    NekDouble sin_theta, cos_theta, sin_varphi, cos_varphi;
    NekDouble x0j, x1j, x2j, rad, tmp;

    Array<OneD, NekDouble> eta(nq, 0.0);
    Array<OneD, NekDouble> u(nq, 0.0);
    Array<OneD, NekDouble> v(nq, 0.0);

    Array<OneD, Array<OneD, NekDouble>> uvec(m_spacedim);
    for (int i = 0; i < m_spacedim; ++i)
    {
        uvec[i] = Array<OneD, NekDouble>(nq, 0.0);
    }

    Array<OneD, NekDouble> x(nq);
    Array<OneD, NekDouble> y(nq);
    Array<OneD, NekDouble> z(nq);

    m_fields[0]->GetCoords(x, y, z);

    for (int j = 0; j < nq; ++j)
    {
        x0j = x[j];
        x1j = y[j];
        x2j = z[j];

        CartesianToElliptical(x0j, x1j, x2j, rad, sin_varphi, cos_varphi,
                              sin_theta, cos_theta);

        // H = H_0 - (1/g)*(a \Omega u_0 + 0.5*u_0^2 )*(- \cos \phi \cos \theta
        // \sin \alpha + \sin \theta \cos \alpha )^2
        tmp = -1.0 * cos_varphi * sin_theta * sin(m_alpha) +
              cos_theta * cos(m_alpha);
        eta[j] = m_H0 - m_Hvar * tmp * tmp;

        // u = (\vec{u} \cdot e^1 )/ || e^1 ||^2 ,   v = (\vec{u} \cdot e^2 )/
        // || e^2 ||^2
        uhat = m_u0 * (sin_theta * cos(m_alpha) +
                       cos_theta * cos_varphi * sin(m_alpha));
        vhat = -1.0 * m_u0 * sin_varphi * sin(m_alpha);

        uvec[0][j] = -1.0 * uhat * sin_varphi - vhat * cos_theta * cos_varphi;
        uvec[1][j] = uhat * cos_varphi - vhat * cos_theta * sin_varphi;
        uvec[2][j] = vhat * sin_theta;
    }

    switch (field)
    {
        case (0):
        {
            outfield = eta;
        }
        break;

        case (1):
        {
            outfield = CartesianToMovingframes(m_movingframes, uvec, 0);
        }
        break;

        case (2):
        {
            outfield = CartesianToMovingframes(m_movingframes, uvec, 1);
        }
        break;
    }
}

void MMFSWE::UnsteadyZonalFlow(unsigned int field, const NekDouble time,
                               Array<OneD, NekDouble> &outfield)
{
    int nq = GetTotPoints();
    NekDouble uhat, vhat;
    NekDouble sin_theta, cos_theta, sin_varphi, cos_varphi;
    NekDouble x0j, x1j, x2j, rad, tmp;

    NekDouble TR, Ttheta;

    Array<OneD, NekDouble> eta(nq, 0.0);
    Array<OneD, NekDouble> u(nq, 0.0);
    Array<OneD, NekDouble> v(nq, 0.0);

    Array<OneD, Array<OneD, NekDouble>> uvec(m_spacedim);
    for (int i = 0; i < m_spacedim; ++i)
    {
        uvec[i] = Array<OneD, NekDouble>(nq, 0.0);
    }

    Array<OneD, NekDouble> x(nq);
    Array<OneD, NekDouble> y(nq);
    Array<OneD, NekDouble> z(nq);

    m_fields[0]->GetCoords(x, y, z);

    for (int j = 0; j < nq; ++j)
    {
        x0j = x[j];
        x1j = y[j];
        x2j = z[j];

        CartesianToElliptical(x0j, x1j, x2j, rad, sin_varphi, cos_varphi,
                              sin_theta, cos_theta);

        // \eta = ( - ( u_0 ( - T_R sin \alpha \cos \theta + \cos \alpha \sin
        // \theta ) + \Omega \sin \theta )^2 + \Omega \sin \theta )^2 + (\Omega
        // \sin \theta )^2
        TR =
            cos_varphi * cos(m_Omega * time) - sin_varphi * sin(m_Omega * time);
        Ttheta =
            sin_varphi * cos(m_Omega * time) + cos_varphi * sin(m_Omega * time);
        tmp = -1.0 * TR * sin(m_alpha) * sin_theta + cos(m_alpha) * cos_theta;

        eta[j] = -1.0 * (m_u0 * tmp + m_Omega * cos_theta) *
                     (m_u0 * tmp + m_Omega * cos_theta) +
                 m_Omega * m_Omega * cos_theta * cos_theta;
        eta[j] = 0.5 * eta[j] / m_g;

        // u = u_0*(TR*\sin \alpha * \sin \theta + \cos \alpha * \cos \theta
        // v = - u_0 Ttheta * \sin \alpha
        uhat =
            m_u0 * (TR * sin(m_alpha) * cos_theta + cos(m_alpha) * sin_theta);
        vhat = -1.0 * m_u0 * Ttheta * sin(m_alpha);

        uvec[0][j] = -1.0 * uhat * sin_varphi - vhat * cos_theta * cos_varphi;
        uvec[1][j] = uhat * cos_varphi - vhat * cos_theta * sin_varphi;
        uvec[2][j] = vhat * sin_theta;
    }

    switch (field)
    {
        case (0):
        {
            outfield = eta;
        }
        break;

        case (1):
        {
            outfield = CartesianToMovingframes(m_movingframes, uvec, 0);
        }
        break;

        case (2):
        {
            outfield = CartesianToMovingframes(m_movingframes, uvec, 1);
        }
        break;
    }
}

NekDouble MMFSWE::ComputeMass(const Array<OneD, const NekDouble> &eta)
{
    int nq = m_fields[0]->GetTotPoints();

    Array<OneD, NekDouble> tmp(nq);
    Vmath::Vadd(nq, eta, 1, m_depth, 1, tmp, 1);

    return m_fields[0]->PhysIntegral(tmp);
}

NekDouble MMFSWE::ComputeEnergy(const Array<OneD, const NekDouble> &eta,
                                const Array<OneD, const NekDouble> &u,
                                const Array<OneD, const NekDouble> &v)
{
    int nq = m_fields[0]->GetTotPoints();

    Array<OneD, NekDouble> tmp(nq);
    Array<OneD, NekDouble> htmp(nq);
    Array<OneD, NekDouble> hstmp(nq);

    Vmath::Vmul(nq, u, 1, u, 1, tmp, 1);
    Vmath::Vvtvp(nq, v, 1, v, 1, tmp, 1, tmp, 1);
    Vmath::Vmul(nq, eta, 1, tmp, 1, tmp, 1);

    Vmath::Sadd(nq, m_H0, eta, 1, htmp, 1);
    Vmath::Vmul(nq, htmp, 1, htmp, 1, htmp, 1);

    Vmath::Sadd(nq, -1.0 * m_H0, m_depth, 1, hstmp, 1);
    Vmath::Vmul(nq, hstmp, 1, hstmp, 1, hstmp, 1);

    Vmath::Vsub(nq, htmp, 1, hstmp, 1, htmp, 1);
    Vmath::Smul(nq, m_g, htmp, 1, htmp, 1);

    Vmath::Vadd(nq, htmp, 1, tmp, 1, tmp, 1);
    Vmath::Smul(nq, 0.5, tmp, 1, tmp, 1);

    return m_fields[0]->PhysIntegral(tmp);
}

NekDouble MMFSWE::ComputeEnstrophy(const Array<OneD, const NekDouble> &eta,
                                   const Array<OneD, const NekDouble> &u,
                                   const Array<OneD, const NekDouble> &v)
{
    int nq = m_fields[0]->GetTotPoints();

    Array<OneD, NekDouble> hstartmp(nq);
    Array<OneD, NekDouble> tmp(nq);

    Vmath::Vadd(nq, eta, 1, m_depth, 1, hstartmp, 1);

    Array<OneD, Array<OneD, NekDouble>> uvec(m_spacedim);
    Array<OneD, Array<OneD, NekDouble>> Curlu(m_spacedim);
    for (int i = 0; i < m_spacedim; ++i)
    {
        uvec[i]  = Array<OneD, NekDouble>(nq);
        Curlu[i] = Array<OneD, NekDouble>(nq);
    }

    ComputeVorticity(u, v, tmp);

    NekDouble Tol = 0.0001;
    for (int i = 0; i < nq; i++)
    {
        if (abs(hstartmp[i]) > Tol)
        {
            tmp[i] = 0.5 * (tmp[i] + m_coriolis[i]) * (tmp[i] + m_coriolis[i]) /
                     hstartmp[i];
        }

        else
        {
            // std::cout << "WARNING: Enstropy is not defined at i= " << i << ",
            // eta = " << eta[i] << ", m_depth = " << m_depth[i] << std::endl;
        }
    }

    // Vmath::Vadd(nq, m_coriolis, 1, tmp, 1, tmp, 1);
    // Vmath::Vmul(nq, tmp, 1, tmp, 1, tmp, 1);
    // Vmath::Vdiv(nq, tmp, 1, hstartmp, 1, tmp, 1);
    // Vmath::Smul(nq, 0.5, tmp, 1, tmp, 1);

    return m_fields[0]->PhysIntegral(tmp);
}

// Vorticity = \nabla v \cdot e^1 + v \nabla \cdot e^1 - ( \nabla u \cdot e^2 +
// u \nabla \cdot e^2 )
void MMFSWE::ComputeVorticity(const Array<OneD, const NekDouble> &u,
                              const Array<OneD, const NekDouble> &v,
                              Array<OneD, NekDouble> &Vorticity)
{

    int nq = m_fields[0]->GetTotPoints();

    Array<OneD, NekDouble> tmp(nq);

    Vorticity = Array<OneD, NekDouble>(nq, 0.0);

    MMFDirectionalDeriv(m_movingframes[0], v, Vorticity);
    Vmath::Vvtvp(nq, &v[0], 1, &m_CurlMF[1][0], 1, &Vorticity[0], 1,
                 &Vorticity[0], 1);

    MMFDirectionalDeriv(m_movingframes[1], u, tmp);
    Vmath::Neg(nq, tmp, 1);
    Vmath::Vvtvp(nq, &u[0], 1, &m_CurlMF[0][0], 1, &tmp[0], 1, &tmp[0], 1);

    Vmath::Vadd(nq, tmp, 1, Vorticity, 1, Vorticity, 1);
}

void MMFSWE::ComputeNablaCdotVelocity(Array<OneD, NekDouble> &vellc)
{
    int nq = m_fields[0]->GetNpoints();

    Array<OneD, NekDouble> velcoeff(nq, 0.0);

    Array<OneD, NekDouble> Dtmp0(nq);
    Array<OneD, NekDouble> Dtmp1(nq);
    Array<OneD, NekDouble> Dtmp2(nq);
    Array<OneD, NekDouble> Drv(nq);

    vellc = Array<OneD, NekDouble>(nq, 0.0);

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
}

void MMFSWE::IsolatedMountainFlow(unsigned int field, const NekDouble time,
                                  Array<OneD, NekDouble> &outfield)
{
    boost::ignore_unused(time);

    int nq = GetTotPoints();

    NekDouble uhat, vhat;
    NekDouble sin_theta, cos_theta, sin_varphi, cos_varphi;
    NekDouble x0j, x1j, x2j, rad;

    Array<OneD, NekDouble> eta(nq, 0.0);
    Array<OneD, NekDouble> u(nq, 0.0);
    Array<OneD, NekDouble> v(nq, 0.0);

    Array<OneD, Array<OneD, NekDouble>> uvec(m_spacedim);
    for (int i = 0; i < m_spacedim; ++i)
    {
        uvec[i] = Array<OneD, NekDouble>(nq, 0.0);
    }

    Array<OneD, NekDouble> x(nq);
    Array<OneD, NekDouble> y(nq);
    Array<OneD, NekDouble> z(nq);

    m_fields[0]->GetCoords(x, y, z);

    for (int j = 0; j < nq; ++j)
    {
        x0j = x[j];
        x1j = y[j];
        x2j = z[j];

        CartesianToElliptical(x0j, x1j, x2j, rad, sin_varphi, cos_varphi,
                              sin_theta, cos_theta);

        // eta = - (1/g) (\Omega u_0 + 0.4 u^2 ) \sin^2 \theta
        eta[j] = (-1.0 / m_g) * (m_Omega * m_u0 + 0.5 * m_u0 * m_u0) *
                 cos_theta * cos_theta;

        uhat = m_u0 * sin_theta;
        vhat = 0.0;

        uvec[0][j] = -1.0 * uhat * sin_varphi - vhat * cos_theta * cos_varphi;
        uvec[1][j] = uhat * cos_varphi - vhat * cos_theta * sin_varphi;
        uvec[2][j] = vhat * sin_theta;
    }

    switch (field)
    {
        case (0):
        {
            outfield = eta;
        }
        break;

        case (1):
        {
            outfield = CartesianToMovingframes(m_movingframes, uvec, 0);
        }
        break;

        case (2):
        {
            outfield = CartesianToMovingframes(m_movingframes, uvec, 1);
        }
        break;
    }
}

void MMFSWE::UnstableJetFlow(unsigned int field, const NekDouble time,
                             Array<OneD, NekDouble> &outfield)
{
    boost::ignore_unused(time);

    int nq = GetTotPoints();

    NekDouble uhat, vhat;
    NekDouble sin_theta, cos_theta, sin_varphi, cos_varphi;
    NekDouble x0j, x1j, x2j, rad;
    NekDouble Ttheta, Tphi;

    Array<OneD, NekDouble> eta(nq, 0.0);
    Array<OneD, NekDouble> u(nq, 0.0);
    Array<OneD, NekDouble> v(nq, 0.0);

    Array<OneD, Array<OneD, NekDouble>> uvec(m_spacedim);
    for (int i = 0; i < m_spacedim; ++i)
    {
        uvec[i] = Array<OneD, NekDouble>(nq, 0.0);
    }

    Array<OneD, NekDouble> x(nq);
    Array<OneD, NekDouble> y(nq);
    Array<OneD, NekDouble> z(nq);

    m_fields[0]->GetCoords(x, y, z);

    int Nint = 1000;
    NekDouble dth, intj;
    for (int j = 0; j < nq; ++j)
    {
        x0j = x[j];
        x1j = y[j];
        x2j = z[j];

        CartesianToElliptical(x0j, x1j, x2j, rad, sin_varphi, cos_varphi,
                              sin_theta, cos_theta);

        Ttheta = atan2(cos_theta, sin_theta);
        Tphi   = atan2(sin_varphi, cos_varphi);

        uhat = ComputeUnstableJetuphi(Ttheta);
        vhat = 0.0;

        // eta = - (1/g) (\Omega u_0 + 0.4 u^2 ) \sin^2 \theta
        dth    = Ttheta / Nint;
        eta[j] = dth * 0.5 *
                 (ComputeUnstableJetEta(0.0) + ComputeUnstableJetEta(Ttheta));
        for (int i = 1; i < Nint - 1; i++)
        {
            intj   = i * dth;
            eta[j] = eta[j] + dth * ComputeUnstableJetEta(intj);
        }
        eta[j] = (-1.0 / m_g) * eta[j];

        // Add perturbation
        if (m_PurturbedJet)
        {
            eta[j] = eta[j] + m_hbar * sin_theta * exp(-9.0 * Tphi * Tphi) *
                                  exp(-225.0 * (m_pi / 4.0 - Ttheta) *
                                      (m_pi / 4.0 - Ttheta));
        }

        uvec[0][j] = -1.0 * uhat * sin_varphi - vhat * cos_theta * cos_varphi;
        uvec[1][j] = uhat * cos_varphi - vhat * cos_theta * sin_varphi;
        uvec[2][j] = vhat * sin_theta;
    }

    // Projection of u onto the tangent plane with conserving the mag. of the
    // velocity.
    Array<OneD, Array<OneD, NekDouble>> uvecproj(m_spacedim);

    for (int i = 0; i < m_spacedim; ++i)
    {
        uvecproj[i] = Array<OneD, NekDouble>(nq, 0.0);
    }

    // u is projected on the tangent plane with preserving its length
    // GramSchumitz(m_surfaceNormal, uvec, uvecproj, true);

    // Change it to the coordinate of moving frames
    // CartesianToMovingframes(0,uvecproj,u);
    // CartesianToMovingframes(1,uvecproj,v);

    u = CartesianToMovingframes(m_movingframes, uvec, 0);
    v = CartesianToMovingframes(m_movingframes, uvec, 1);

    switch (field)
    {
        case (0):
        {
            outfield = eta;
        }
        break;

        case (1):
        {
            outfield = u;
        }
        break;

        case (2):
        {
            outfield = v;
        }
        break;
    }
}

void MMFSWE::RossbyWave(unsigned int field, Array<OneD, NekDouble> &outfield)
{
    int nq = GetTotPoints();
    NekDouble uhat, vhat;
    NekDouble sin_theta, cos_theta, sin_varphi, cos_varphi;
    NekDouble x0j, x1j, x2j, rad;
    NekDouble Ath, Bth, Cth, tmp;

    Array<OneD, NekDouble> eta(nq, 0.0);
    Array<OneD, NekDouble> u(nq, 0.0);
    Array<OneD, NekDouble> v(nq, 0.0);

    Array<OneD, Array<OneD, NekDouble>> uvec(m_spacedim);
    for (int i = 0; i < m_spacedim; ++i)
    {
        uvec[i] = Array<OneD, NekDouble>(nq, 0.0);
    }

    Array<OneD, NekDouble> x(nq);
    Array<OneD, NekDouble> y(nq);
    Array<OneD, NekDouble> z(nq);

    m_fields[0]->GetCoords(x, y, z);

    NekDouble R = 4.0;
    NekDouble sin2theta, sinRtheta, sin6theta, sin2Rtheta, sinRm1theta;
    NekDouble cos2phi, cos4phi, sin4phi, cos8phi;

    // disturbancees of Rossby-Haurwitz Wave
    NekDouble x0d, y0d, z0d, phi0, theta0;

    phi0   = 40.0 * m_pi / 180.0;
    theta0 = 50.0 * m_pi / 180.0;

    x0d = cos(phi0) * cos(theta0);
    y0d = sin(phi0) * cos(theta0);
    z0d = sin(theta0);

    for (int j = 0; j < nq; ++j)
    {
        x0j = x[j];
        x1j = y[j];
        x2j = z[j];

        CartesianToElliptical(x0j, x1j, x2j, rad, sin_varphi, cos_varphi,
                              sin_theta, cos_theta);

        // H = H_0 - (1/g)*(a \Omega u_0 + 0.5*u_0^2 )*(- \cos \phi \cos \theta
        // \sin \alpha + \sin \theta \cos \alpha )^2

        // tmp = cos^{2R} \theta, R = 4;
        sin2theta   = sin_theta * sin_theta;
        sinRm1theta = sin_theta * sin2theta;
        sinRtheta   = sin2theta * sin2theta;
        sin6theta   = sin2theta * sinRtheta;
        sin2Rtheta  = sinRtheta * sinRtheta;

        tmp = (0.5 * m_angfreq) * (2.0 * m_Omega + m_angfreq);
        Ath = tmp * sin2theta +
              0.25 * m_K * m_K * sin6theta *
                  ((R + 1.0) * sinRtheta + (2 * R * R - R - 2.0) * sin2theta -
                   2 * R * R);

        tmp = (2.0 * m_K) * (m_Omega + m_angfreq) / ((R + 1.0) * (R + 2.0));
        Bth = tmp * sinRtheta *
              ((R * R + 2 * R + 2) - (R + 1.0) * (R + 1.0) * sin2theta);

        Cth =
            0.25 * m_K * m_K * sin2Rtheta * ((R + 1.0) * sin2theta - (R + 2.0));

        // cos (2 \phi) = 2 * cos \phi * cos \phi - 1.0
        cos2phi = 2.0 * cos_varphi * cos_varphi - 1.0;
        cos4phi = 2.0 * cos2phi * cos2phi - 1.0;
        cos8phi = 2.0 * cos4phi * cos4phi - 1.0;

        // sin (2 \phi) = 2 * cos \phi * sin \phi
        sin4phi = 4.0 * sin_varphi * cos_varphi * cos2phi;

        eta[j] = m_H0 + (1.0 / m_g) * (Ath + Bth * cos4phi + Cth * cos8phi);

        // disturbances is added
        if (m_RossbyDisturbance)
        {
            eta[j] = eta[j] *
                     (1.0 + (1.0 / 40.0) * (x0j * x0d + x1j * y0d + x2j * z0d));
        }

        // u = (\vec{u} \cdot e^1 )/ || e^1 ||^2 ,   v = (\vec{u} \cdot e^2 )/
        // || e^2 ||^2
        uhat = m_angfreq * sin_theta +
               m_K * sinRm1theta * (R * cos_theta * cos_theta - sin2theta) *
                   cos4phi;
        vhat = -1.0 * m_K * R * sinRm1theta * cos_theta * sin4phi;

        uvec[0][j] = -1.0 * uhat * sin_varphi - vhat * cos_theta * cos_varphi;
        uvec[1][j] = uhat * cos_varphi - vhat * cos_theta * sin_varphi;
        uvec[2][j] = vhat * sin_theta;

        // CartesianToSpherical(x0j, x1j, x2j, sin_varphi, cos_varphi,
        // sin_theta,
        //                      cos_theta);

        // // H = H_0 - (1/g)*(a \Omega u_0 + 0.5*u_0^2 )*(- \cos \phi \cos
        // \theta
        // // \sin \alpha + \sin \theta \cos \alpha )^2

        // // tmp = cos^{2R} \theta, R = 4;
        // cos2theta = cos_theta * cos_theta;
        // cosRm1theta = cos_theta * cos2theta;
        // cosRtheta = cos2theta * cos2theta;
        // cos6theta = cos2theta * cosRtheta;
        // cos2Rtheta = cosRtheta * cosRtheta;

        // tmp = (0.5 * m_angfreq) * (2.0 * m_Omega + m_angfreq);
        // Ath = tmp * cos2theta +
        //       0.25 * m_K * m_K * cos6theta *
        //           ((R + 1.0) * cosRtheta + (2 * R * R - R - 2.0) * cos2theta
        //           -
        //            2 * R * R);

        // tmp = (2.0 * m_K) * (m_Omega + m_angfreq) / ((R + 1.0) * (R + 2.0));
        // Bth = tmp * cosRtheta *
        //       ((R * R + 2 * R + 2) - (R + 1.0) * (R + 1.0) * cos2theta);

        // Cth =
        //     0.25 * m_K * m_K * cos2Rtheta * ((R + 1.0) * cos2theta - (R
        //     + 2.0));

        // // cos (2 \phi) = 2 * cos \phi * cos \phi - 1.0
        // cos2phi = 2.0 * cos_varphi * cos_varphi - 1.0;
        // cos4phi = 2.0 * cos2phi * cos2phi - 1.0;
        // cos8phi = 2.0 * cos4phi * cos4phi - 1.0;

        // // sin (2 \phi) = 2 * cos \phi * sin \phi
        // sin4phi = 4.0 * sin_varphi * cos_varphi * cos2phi;

        // eta[j] = m_H0 + (1.0 / m_g) * (Ath + Bth * cos4phi + Cth * cos8phi);

        // // disturbances is added
        // if (m_RossbyDisturbance)
        // {
        //     eta[j] = eta[j] * (1.0 + (1.0 / 40.0) * (x0j * x0d + x1j * y0d +
        //     x2j * z0d));
        // }

        // // u = (\vec{u} \cdot e^1 )/ || e^1 ||^2 ,   v = (\vec{u} \cdot e^2
        // )/
        // // || e^2 ||^2
        // uhat = m_angfreq * cos_theta +
        //        m_K * cosRm1theta * (R * sin_theta * sin_theta - cos2theta) *
        //            cos4phi;
        // vhat = -1.0 * m_K * R * cosRm1theta * sin_theta * sin4phi;

        // uvec[0][j] = -1.0 * uhat * sin_varphi - vhat * sin_theta *
        // cos_varphi; uvec[1][j] = uhat * cos_varphi - vhat * sin_theta *
        // sin_varphi; uvec[2][j] = vhat * cos_theta;
    }

    // NekDouble etamin, etaaver;
    // etamin = Vmath::Vmin(nq, eta, 1)*rad_earth;
    // etaaver = Vmath::Vsum(nq, eta, 1)/nq*rad_earth;

    // Projection of u onto the tangent plane with conserving the mag. of the
    // velocity.
    Array<OneD, Array<OneD, NekDouble>> uvecproj(m_spacedim);

    for (int i = 0; i < m_spacedim; ++i)
    {
        uvecproj[i] = Array<OneD, NekDouble>(nq, 0.0);
    }

    // u is projected on the tangent plane with preserving its length
    // GramSchumitz(m_surfaceNormal, uvec, uvecproj, true);

    // Change it to the coordinate of moving frames
    // CartesianToMovingframes(0,uvecproj,u);
    // CartesianToMovingframes(1,uvecproj,v);

    u = CartesianToMovingframes(m_movingframes, uvec, 0);
    v = CartesianToMovingframes(m_movingframes, uvec, 1);

    switch (field)
    {
        case (0):
        {
            outfield = eta;
        }
        break;

        case (1):
        {
            outfield = u;
        }
        break;

        case (2):
        {
            outfield = v;
        }
        break;
    }
}

NekDouble MMFSWE::ComputeUnstableJetEta(const NekDouble theta)
{
    NekDouble uphi, f, dh;

    uphi = ComputeUnstableJetuphi(theta);
    f    = 2.0 * m_Omega * sin(theta);

    dh = f * uphi + tan(theta) * uphi * uphi;

    return dh;
}

NekDouble MMFSWE::ComputeUnstableJetuphi(const NekDouble theta)
{
    NekDouble uphi;

    if ((theta > m_theta0) && (theta < m_theta1))
    {
        uphi = (m_uthetamax / m_en) *
               exp(1.0 / (theta - m_theta0) / (theta - m_theta1));
    }

    else
    {
        uphi = 0.0;
    }

    return uphi;
}

void MMFSWE::Checkpoint_Output_Cartesian(std::string outname)
{
    int nq      = m_fields[0]->GetTotPoints();
    int ncoeffs = m_fields[0]->GetNcoeffs();

    NekDouble rad_earth = 6.37122 * 1000000;

    int nvariables = 7;

    // vector u in Cartesian coordinates
    std::vector<std::string> variables(nvariables);

    variables[0] = "eta";
    variables[1] = "hstar";
    variables[2] = "vorticity";
    variables[3] = "ux";
    variables[4] = "uy";
    variables[5] = "uz";
    variables[6] = "null";

    // Obtain \vec{u} in cartesian coordinate
    Array<OneD, Array<OneD, NekDouble>> fieldphys(nvariables);
    std::vector<Array<OneD, NekDouble>> fieldcoeffs(nvariables);
    for (int i = 0; i < nvariables; ++i)
    {
        fieldphys[i]   = Array<OneD, NekDouble>(nq, 0.0);
        fieldcoeffs[i] = Array<OneD, NekDouble>(ncoeffs);
    }

    Vmath::Smul(nq, rad_earth, &(m_fields[0]->GetPhys())[0], 1,
                &fieldphys[0][0], 1);
    // Vmath::Vadd(nq, &(m_fields[0]->GetPhys())[0], 1, &m_depth[0], 1,
    // &fieldphys[1][0], 1);

    Vmath::Vcopy(nq, &m_depth[0], 1, &fieldphys[1][0], 1);
    Vmath::Neg(nq, &fieldphys[1][0], 1);
    Vmath::Sadd(nq, m_H0, &fieldphys[1][0], 1, &fieldphys[1][0], 1);
    Vmath::Smul(nq, rad_earth, &fieldphys[1][0], 1, &fieldphys[1][0], 1);

    Array<OneD, NekDouble> utmp(nq);
    Array<OneD, NekDouble> vtmp(nq);

    Vmath::Vcopy(nq, &(m_fields[1]->GetPhys())[0], 1, &utmp[0], 1);
    Vmath::Vcopy(nq, &(m_fields[2]->GetPhys())[0], 1, &vtmp[0], 1);

    ComputeVorticity(utmp, vtmp, fieldphys[2]);
    // u_x = u e^1_x + v e^2_x
    int indx = 3;
    for (int k = 0; k < m_spacedim; ++k)
    {
        Vmath::Vmul(nq, &utmp[0], 1, &m_movingframes[0][k * nq], 1,
                    &fieldphys[k + indx][0], 1);
        Vmath::Vvtvp(nq, &vtmp[0], 1, &m_movingframes[1][k * nq], 1,
                     &fieldphys[k + indx][0], 1, &fieldphys[k + indx][0], 1);
    }

    for (int i = 0; i < nvariables; ++i)
    {
        m_fields[0]->FwdTrans(fieldphys[i], fieldcoeffs[i]);
    }

    WriteFld(outname, m_fields[0], fieldcoeffs, variables);
}

// (h, hu, hv) -> (\eta, u, v)
void MMFSWE::ConservativeToPrimitive(
    const Array<OneD, const Array<OneD, NekDouble>> &physin,
    Array<OneD, Array<OneD, NekDouble>> &physout)
{
    int nq = GetTotPoints();

    if (physin[0].get() == physout[0].get())
    {
        // copy indata and work with tmp array
        Array<OneD, Array<OneD, NekDouble>> tmp(3);
        for (int i = 0; i < 3; ++i)
        {
            // deep copy
            tmp[i] = Array<OneD, NekDouble>(nq);
            Vmath::Vcopy(nq, physin[i], 1, tmp[i], 1);
        }

        // \eta = h - d
        Vmath::Vsub(nq, tmp[0], 1, m_depth, 1, physout[0], 1);

        // u = hu/h
        Vmath::Vdiv(nq, tmp[1], 1, tmp[0], 1, physout[1], 1);

        // v = hv/h
        Vmath::Vdiv(nq, tmp[2], 1, tmp[0], 1, physout[2], 1);
    }
    else
    {
        // \eta = h - d
        Vmath::Vsub(nq, physin[0], 1, m_depth, 1, physout[0], 1);

        // u = hu/h
        Vmath::Vdiv(nq, physin[1], 1, physin[0], 1, physout[1], 1);

        // v = hv/h
        Vmath::Vdiv(nq, physin[2], 1, physin[0], 1, physout[2], 1);
    }
}

void MMFSWE::PrimitiveToConservative(
    const Array<OneD, const Array<OneD, NekDouble>> &physin,
    Array<OneD, Array<OneD, NekDouble>> &physout)
{

    int nq = GetTotPoints();

    if (physin[0].get() == physout[0].get())
    {
        // copy indata and work with tmp array
        Array<OneD, Array<OneD, NekDouble>> tmp(3);
        for (int i = 0; i < 3; ++i)
        {
            // deep copy
            tmp[i] = Array<OneD, NekDouble>(nq);
            Vmath::Vcopy(nq, physin[i], 1, tmp[i], 1);
        }

        // h = \eta + d
        Vmath::Vadd(nq, tmp[0], 1, m_depth, 1, physout[0], 1);

        // hu = h * u
        Vmath::Vmul(nq, physout[0], 1, tmp[1], 1, physout[1], 1);

        // hv = h * v
        Vmath::Vmul(nq, physout[0], 1, tmp[2], 1, physout[2], 1);
    }
    else
    {
        // h = \eta + d
        Vmath::Vadd(nq, physin[0], 1, m_depth, 1, physout[0], 1);

        // hu = h * u
        Vmath::Vmul(nq, physout[0], 1, physin[1], 1, physout[1], 1);

        // hv = h * v
        Vmath::Vmul(nq, physout[0], 1, physin[2], 1, physout[2], 1);
    }
}

// To obtain [eta, u, v ] from [H, Hu, Hv]
void MMFSWE::ConservativeToPrimitive(void)
{
    int nq = GetTotPoints();

    // u = hu/h
    Vmath::Vdiv(nq, m_fields[1]->GetPhys(), 1, m_fields[0]->GetPhys(), 1,
                m_fields[1]->UpdatePhys(), 1);

    // v = hv/ v
    Vmath::Vdiv(nq, m_fields[2]->GetPhys(), 1, m_fields[0]->GetPhys(), 1,
                m_fields[2]->UpdatePhys(), 1);

    // \eta = h - d
    Vmath::Vsub(nq, m_fields[0]->GetPhys(), 1, m_depth, 1,
                m_fields[0]->UpdatePhys(), 1);
}

void MMFSWE::PrimitiveToConservative(void)
{
    int nq = GetTotPoints();

    // h = \eta + d
    Vmath::Vadd(nq, m_fields[0]->GetPhys(), 1, m_depth, 1,
                m_fields[0]->UpdatePhys(), 1);

    // hu = h * u
    Vmath::Vmul(nq, m_fields[0]->GetPhys(), 1, m_fields[1]->GetPhys(), 1,
                m_fields[1]->UpdatePhys(), 1);

    // hv = h * v
    Vmath::Vmul(nq, m_fields[0]->GetPhys(), 1, m_fields[2]->GetPhys(), 1,
                m_fields[2]->UpdatePhys(), 1);
}

void MMFSWE::TestVorticityComputation(void)
{
    // Construct beta
    int i, k;
    int n = 1, m = 1;
    int nq = m_fields[0]->GetTotPoints();

    NekDouble alpha, beta_theta, beta_phi;

    NekDouble xp, yp, zp, rad, Re;
    NekDouble theta, phi, sin_theta, cos_theta, sin_varphi, cos_varphi;
    NekDouble cosntheta3;

    NekDouble thetax, thetay, thetaz;
    NekDouble phix, phiy, phiz;

    Array<OneD, NekDouble> x(nq);
    Array<OneD, NekDouble> y(nq);
    Array<OneD, NekDouble> z(nq);

    Array<OneD, NekDouble> u(nq);
    Array<OneD, NekDouble> v(nq);

    Array<OneD, NekDouble> vorticitycompt(nq);
    Array<OneD, NekDouble> vorticityexact(nq);

    m_fields[0]->GetCoords(x, y, z);

    Array<OneD, Array<OneD, NekDouble>> uvec(m_spacedim);
    for (i = 0; i < m_spacedim; ++i)
    {
        uvec[i] = Array<OneD, NekDouble>(nq, 0.0);
    }

    // Generate SphericalCoords
    for (k = 0; k < nq; ++k)
    {
        xp = x[k];
        yp = y[k];
        zp = z[k];

        Re = sqrt(xp * xp + yp * yp + zp * zp);

        CartesianToElliptical(xp, yp, zp, rad, sin_varphi, cos_varphi,
                              sin_theta, cos_theta);

        alpha = sin_varphi;

        theta = atan2(cos_theta, sin_theta);
        phi   = atan2(sin_varphi, cos_varphi);

        cosntheta3 = cos(n * theta) * cos(n * theta) * cos(n * theta);
        // cosntheta4 = cosntheta3 * cos(n * theta);

        beta_theta = -4.0 * n * cosntheta3 * cos(m * phi) * sin(n * theta) / Re;
        beta_phi   = -m * cosntheta3 * sin(m * phi) / Re;

        thetax = -1.0 * cos_varphi * sin_theta;
        thetay = -1.0 * sin_varphi * sin_theta;
        thetaz = cos_theta;

        phix = -1.0 * sin_varphi;
        phiy = cos_varphi;
        phiz = 0.0;

        uvec[0][k] = alpha * (beta_theta * thetax + beta_phi * phix);
        uvec[1][k] = alpha * (beta_theta * thetay + beta_phi * phiy);
        uvec[2][k] = alpha * (beta_theta * thetaz + beta_phi * phiz);

        vorticityexact[k] = -4.0 * n / Re / Re * sin_theta * sin_theta *
                            cos_varphi * cos(m * phi) * sin(n * theta);

        // CartesianToSpherical(xp, yp, zp, sin_varphi, cos_varphi, sin_theta,
        //                      cos_theta);

        // alpha = sin_varphi;

        // theta = atan2(sin_theta, cos_theta);
        // phi = atan2(sin_varphi, cos_varphi);

        // cosntheta3 = cos(n * theta) * cos(n * theta) * cos(n * theta);
        // // cosntheta4 = cosntheta3 * cos(n * theta);

        // beta_theta = -4.0 * n * cosntheta3 * cos(m * phi) * sin(n * theta) /
        // Re; beta_phi = -m * cosntheta3 * sin(m * phi) / Re;

        // thetax = -1.0 * cos_varphi * sin_theta;
        // thetay = -1.0 * sin_varphi * sin_theta;
        // thetaz = cos_theta;

        // phix = -1.0 * sin_varphi;
        // phiy = cos_varphi;
        // phiz = 0.0;

        // uvec[0][k] = alpha * (beta_theta * thetax + beta_phi * phix);
        // uvec[1][k] = alpha * (beta_theta * thetay + beta_phi * phiy);
        // uvec[2][k] = alpha * (beta_theta * thetaz + beta_phi * phiz);

        // vorticityexact[k] = -4.0 * n / Re / Re * cos_theta * cos_theta *
        //                     cos_varphi * cos(m * phi) * sin(n * theta);
    }

    u = CartesianToMovingframes(m_movingframes, uvec, 0);
    v = CartesianToMovingframes(m_movingframes, uvec, 1);

    ComputeVorticity(u, v, vorticitycompt);

    Vmath::Vsub(nq, vorticityexact, 1, vorticitycompt, 1, vorticitycompt, 1);

    std::cout << "Vorticity: L2 error = " << AvgAbsInt(vorticitycompt)
              << ", Linf error =  " << Vmath::Vamax(nq, vorticitycompt, 1)
              << std::endl;
}

NekDouble MMFSWE::v_L2Error(unsigned int field,
                            const Array<OneD, NekDouble> &exactsoln,
                            bool Normalised)
{
    boost::ignore_unused(exactsoln);

    int nq            = m_fields[field]->GetNpoints();
    NekDouble L2error = -1.0;

    if (m_NumQuadPointsError == 0)
    {
        if (m_fields[field]->GetPhysState() == false)
        {
            m_fields[field]->BwdTrans(m_fields[field]->GetCoeffs(),
                                      m_fields[field]->UpdatePhys());
        }

        switch (field)
        {
            case (0):
            {
                // I(h) = I (h - h_exact ) / I (h_exact)
                Array<OneD, NekDouble> exactsolution(nq);
                v_EvaluateExactSolution(0, exactsolution, m_time);

                // exactsoln = u - u_T so that L2 compute u_T
                NekDouble L2exact = m_fields[0]->PhysIntegral(exactsolution);

                Vmath::Vsub(nq, &(m_fields[0]->GetPhys())[0], 1,
                            &exactsolution[0], 1, &exactsolution[0], 1);
                Vmath::Vabs(nq, exactsolution, 1, exactsolution, 1);

                L2error = (m_fields[0]->PhysIntegral(exactsolution)) / L2exact;
            }
            break;

            case (1):
            {
                // I2 (u) = I( (u - u_ext)^2 + (v - v_ext)^2 )^{1/2} / I(
                // u_ext^2 + v_ext^2 )^{1/2}
                Array<OneD, NekDouble> exactu(nq);
                Array<OneD, NekDouble> exactv(nq);
                Array<OneD, NekDouble> tmp(nq);

                // L2exact = \int (\sqrt{exactu*exactu+exactv*exactv})
                v_EvaluateExactSolution(1, exactu, m_time);
                v_EvaluateExactSolution(2, exactv, m_time);
                Vmath::Vmul(nq, exactu, 1, exactu, 1, tmp, 1);
                Vmath::Vvtvp(nq, exactv, 1, exactv, 1, tmp, 1, tmp, 1);
                Vmath::Vsqrt(nq, tmp, 1, tmp, 1);

                NekDouble L2exact = m_fields[1]->PhysIntegral(tmp);

                // L2exact = \int
                // (\sqrt{(u-exactu)*(u-exactu)+(v-exactv)*(v-exactv)})
                Vmath::Vsub(nq, &(m_fields[1]->GetPhys())[0], 1, &exactu[0], 1,
                            &exactu[0], 1);
                Vmath::Vsub(nq, &(m_fields[2]->GetPhys())[0], 1, &exactv[0], 1,
                            &exactv[0], 1);
                Vmath::Vmul(nq, exactu, 1, exactu, 1, tmp, 1);
                Vmath::Vvtvp(nq, exactv, 1, exactv, 1, tmp, 1, tmp, 1);
                Vmath::Vsqrt(nq, tmp, 1, tmp, 1);

                L2error = (m_fields[1]->PhysIntegral(tmp)) / L2exact;
            }
            break;

            case (2):
            {
                // L2error = 0.0;
                Array<OneD, NekDouble> exactu(nq);
                Array<OneD, NekDouble> exactv(nq);
                Array<OneD, NekDouble> tmp(nq);
                v_EvaluateExactSolution(1, exactu, m_time);
                v_EvaluateExactSolution(2, exactv, m_time);

                Vmath::Vsub(nq, &(m_fields[1]->GetPhys())[0], 1, &exactu[0], 1,
                            &exactu[0], 1);
                Vmath::Vsub(nq, &(m_fields[2]->GetPhys())[0], 1, &exactv[0], 1,
                            &exactv[0], 1);
                Vmath::Vmul(nq, exactu, 1, exactu, 1, tmp, 1);
                Vmath::Vvtvp(nq, exactv, 1, exactv, 1, tmp, 1, tmp, 1);
                Vmath::Vsqrt(nq, tmp, 1, tmp, 1);

                L2error = RootMeanSquare(tmp);
            }
            break;

            default:
                break;
        }

        if (Normalised == true)
        {
            Array<OneD, NekDouble> one(m_fields[field]->GetNpoints(), 1.0);

            NekDouble Vol = m_fields[field]->PhysIntegral(one);
            m_comm->AllReduce(Vol, LibUtilities::ReduceSum);

            L2error = sqrt(L2error * L2error / Vol);
        }
    }
    else
    {
        Array<OneD, NekDouble> L2INF(2);
        L2INF   = ErrorExtraPoints(field);
        L2error = L2INF[0];
    }

    return L2error;
}

NekDouble MMFSWE::v_LinfError(unsigned int field,
                              const Array<OneD, NekDouble> &exactsoln)
{
    boost::ignore_unused(exactsoln);

    NekDouble LinfError = -1.0;

    if (m_fields[field]->GetPhysState() == false)
    {
        m_fields[field]->BwdTrans(m_fields[field]->GetCoeffs(),
                                  m_fields[field]->UpdatePhys());
    }

    int nq = m_fields[field]->GetNpoints();

    // Obtain \vec{u} in cartesian coordinate
    Array<OneD, NekDouble> x(nq);
    Array<OneD, NekDouble> y(nq);
    Array<OneD, NekDouble> z(nq);

    m_fields[0]->GetCoords(x, y, z);

    switch (field)
    {
        case (0):
        {
            NekDouble Letaint;

            Array<OneD, NekDouble> exactsolution(nq);

            EvaluateExactSolution(field, exactsolution, m_time);
            LinfError = m_fields[field]->Linf(m_fields[field]->GetPhys(),
                                              exactsolution);

            Letaint = Vmath::Vamax(nq, exactsolution, 1);

            Vmath::Vsub(nq, &(m_fields[0]->GetPhys())[0], 1, &exactsolution[0],
                        1, &exactsolution[0], 1);
            // indx = Vmath::Iamax(nq, exactsolution, 1);

            LinfError = fabs(LinfError / Letaint);
        }
        break;

        case (1):
        {
            Array<OneD, NekDouble> exactu(nq);
            Array<OneD, NekDouble> exactv(nq);
            Array<OneD, NekDouble> tmpu(nq);
            Array<OneD, NekDouble> tmpv(nq);
            Array<OneD, NekDouble> Lerr(nq);
            Array<OneD, NekDouble> uT(nq);

            EvaluateExactSolution(1, exactu, m_time);
            EvaluateExactSolution(2, exactv, m_time);

            // Compute max[sqrt{(u-uex)^2 + (v-vex)^2}]
            Vmath::Vcopy(nq, &(m_fields[1]->UpdatePhys())[0], 1, &tmpu[0], 1);
            Vmath::Vcopy(nq, &(m_fields[2]->UpdatePhys())[0], 1, &tmpv[0], 1);

            Vmath::Vsub(nq, &exactu[0], 1, &tmpu[0], 1, &tmpu[0], 1);
            Vmath::Vsub(nq, &exactv[0], 1, &tmpv[0], 1, &tmpv[0], 1);

            Vmath::Vmul(nq, &tmpu[0], 1, &tmpu[0], 1, &tmpu[0], 1);
            Vmath::Vmul(nq, &tmpv[0], 1, &tmpv[0], 1, &tmpv[0], 1);

            Vmath::Vadd(nq, &tmpu[0], 1, &tmpv[0], 1, &Lerr[0], 1);
            Vmath::Vsqrt(nq, &Lerr[0], 1, &Lerr[0], 1);

            // uT = max[sqrt( u_T^2 + v_T^2 ) ]
            Vmath::Vmul(nq, &exactu[0], 1, &exactu[0], 1, &tmpu[0], 1);
            Vmath::Vmul(nq, &exactv[0], 1, &exactv[0], 1, &tmpv[0], 1);
            Vmath::Vadd(nq, &tmpu[0], 1, &tmpv[0], 1, &uT[0], 1);
            Vmath::Vsqrt(nq, &uT[0], 1, &uT[0], 1);

            LinfError = Vmath::Vamax(nq, Lerr, 1) / Vmath::Vamax(nq, uT, 1);
            // indx = Vmath::Iamax(nq, Lerr, 1);
        }
        break;

        case (2):
        {
            LinfError = 0.0;
        }
        break;

        default:
            break;
    }

    return LinfError;
}

void MMFSWE::Checkpoint_ErrMap(
    const NekDouble time,
    const Array<OneD, const Array<OneD, NekDouble>> &fieldphys)
{
    int nvar    = m_fields.size();
    int nq      = m_fields[0]->GetTotPoints();
    int ncoeffs = m_fields[0]->GetNcoeffs();

    std::string outname = m_sessionName + "ErrMap.chk";

    std::vector<Array<OneD, NekDouble>> fieldcoeffs(nvar);
    for (int i = 0; i < nvar; ++i)
    {
        fieldcoeffs[i] = Array<OneD, NekDouble>(ncoeffs);
    }

    std::vector<std::string> variables(nvar);
    variables[0] = "eta";
    variables[1] = "u";
    variables[2] = "v";

    Array<OneD, NekDouble> exactsoln(m_fields[0]->GetNpoints());

    for (int i = 0; i < nvar; ++i)
    {
        v_EvaluateExactSolution(i, exactsoln, time);
        Vmath::Vsub(nq, exactsoln, 1, fieldphys[i], 1, exactsoln, 1);
        Vmath::Vabs(nq, exactsoln, 1, exactsoln, 1);

        m_fields[0]->FwdTrans(exactsoln, fieldcoeffs[i]);
    }
    WriteFld(outname, m_fields[0], fieldcoeffs, variables);
}

void MMFSWE::v_EvaluateExactSolution(unsigned int field,
                                     Array<OneD, NekDouble> &outfield,
                                     const NekDouble time)
{
    switch (m_TestType)
    {
        case eTestSteadyZonal:
        {
            SteadyZonalFlow(field, outfield);
        }
        break;

        case eTestUnsteadyZonal:
        {
            UnsteadyZonalFlow(field, time, outfield);
        }
        break;

        case eTestIsolatedMountain:
        {
            IsolatedMountainFlow(field, time, outfield);
        }
        break;

        case eTestUnstableJet:
        {
            UnstableJetFlow(field, time, outfield);
        }
        break;

        case eTestRossbyWave:
        {
            RossbyWave(field, outfield);
        }
        break;

        case eTestPlane:
        {
            TestSWE2Dproblem(time, field, outfield);
        }
        break;

        default:
            break;
    }
}

void MMFSWE::v_GenerateSummary(SolverUtils::SummaryList &s)
{
    MMFSystem::v_GenerateSummary(s);
    SolverUtils::AddSummaryItem(s, "TestType", TestTypeMap[m_TestType]);
    SolverUtils::AddSummaryItem(s, "Divergence Restore", m_DivergenceRestore);
    SolverUtils::AddSummaryItem(s, "AddCoriolis", m_AddCoriolis);
    SolverUtils::AddSummaryItem(s, "AddRotation", m_AddRotation);

    switch (m_TestType)
    {
        case eTestSteadyZonal:
        {
            SolverUtils::AddSummaryItem(s, "Rotation Angle", m_alpha);
        }
        break;

        case eTestRossbyWave:
        {
            SolverUtils::AddSummaryItem(s, "RossbyDistrubance",
                                        m_RossbyDisturbance);
        }
        break;

        case eTestUnstableJet:
        {
            SolverUtils::AddSummaryItem(s, "PurturbedJet", m_PurturbedJet);
        }
        break;

        default:
            break;
    }
}
} // namespace Nektar
