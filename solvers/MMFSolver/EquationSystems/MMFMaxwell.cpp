/////////////////////////////////////////////////////////////////////////////
//
// File MMFMaxwell.cpp
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
// Description: MMF Maxwell solve routines
//
///////////////////////////////////////////////////////////////////////////////

#include <iomanip>
#include <iostream>
#include <typeinfo>

#include <SolverUtils/MMFSystem.h>

#include <MMFSolver/EquationSystems/MMFMaxwell.h>

#include <boost/core/ignore_unused.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string.hpp>

#include <LibUtilities/BasicUtils/Timer.h>
#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>

#include <LibUtilities/BasicUtils/Timer.h>
#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>

namespace Nektar
{
std::string MMFMaxwell::className =
    SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
        "MMFMaxwell", MMFMaxwell::create, "MMFMaxwell equation.");

MMFMaxwell::MMFMaxwell(const LibUtilities::SessionReaderSharedPtr &pSession,
                       const SpatialDomains::MeshGraphSharedPtr& pGraph)
    : UnsteadySystem(pSession, pGraph), MMFSystem(pSession, pGraph)
{
}

/**
 * @brief Initialisation object for the unsteady linear advection equation.
 */
void MMFMaxwell::v_InitObject(bool DeclareFields)
{
    // Call to the initialisation object
    UnsteadySystem::v_InitObject(DeclareFields);

    int nq       = m_fields[0]->GetNpoints();
    int shapedim = m_fields[0]->GetShapeDimension();

    m_session->LoadParameter("ElemtGroup0", m_ElemtGroup0, 0);
    m_session->LoadParameter("ElemtGroup1", m_ElemtGroup1, 0);
    m_session->LoadParameter("boundaryforSF", m_boundaryforSF, 0);
    m_session->LoadParameter("PrintoutSurfaceCurrent", m_PrintoutSurfaceCurrent,
                             0);
    m_session->LoadParameter("AddRotation", m_AddRotation, 0);


    m_session->LoadParameter("NoInc", m_NoInc, 0);

    // PML parameters
    m_session->LoadParameter("TestPML", m_TestPML, 0);
    m_session->LoadParameter("PMLelement", m_PMLelement, 0);
    m_session->LoadParameter("RecPML", m_RecPML, 0);

    m_session->LoadParameter("AddPML", m_AddPML, 0);
    if (m_AddPML)
    {
        m_RecPML = m_PMLelement;
    }
    m_session->LoadParameter("PMLorder", m_PMLorder, 3);

    m_session->LoadParameter("PMLthickness", m_PMLthickness, 0.0);
    m_session->LoadParameter("PMLstart", m_PMLstart, 0.0);
    m_session->LoadParameter("PMLmaxsigma", m_PMLmaxsigma, 100.0);

    // Point Source parmaters
    m_session->LoadParameter("Psx", m_Psx, 0.0);
    m_session->LoadParameter("Psy", m_Psy, 0.0);
    m_session->LoadParameter("Psz", m_Psz, 0.0);
    m_session->LoadParameter("PSstrength", m_PSstrength, 1.0);
    m_session->LoadParameter("PSduration", m_PSduration, 1.0);

    m_session->LoadParameter("Gaussianradius", m_Gaussianradius, 1.0);

    // Cloaking parameter
    m_session->LoadParameter("CloakNlayer", m_CloakNlayer, 5);
    m_session->LoadParameter("Cloakraddelta", m_Cloakraddelta, 0.0);

    m_varepsilon = Array<OneD, NekDouble>(m_spacedim);
    m_session->LoadParameter("varepsilon1", m_varepsilon[0], 1.0);
    m_session->LoadParameter("varepsilon2", m_varepsilon[1], 1.0);
    m_session->LoadParameter("varepsilon3", m_varepsilon[2], 1.0);

    m_mu = Array<OneD, NekDouble>(m_spacedim);
    m_session->LoadParameter("mu1", m_mu[0], 1.0);
    m_session->LoadParameter("mu2", m_mu[1], 1.0);
    m_session->LoadParameter("mu3", m_mu[2], 1.0);

    Array<OneD, Array<OneD, NekDouble>> AniStrength(shapedim);
    for (int j = 0; j < shapedim; ++j)
    {
        AniStrength[j] = Array<OneD, NekDouble>(nq, 1.0);
    }

    // Add Rectangular PML
    MMFSystem::MMFInitObject(AniStrength);


    m_session->LoadParameter("DivergenceRestore", m_DivergenceRestore, 0);

    if(m_DivergenceRestore>1)
    {
        mf_LOCSPH = Array<OneD, Array<OneD, NekDouble>>(m_mfdim);
        GetLOCALMovingframes(mf_LOCSPH);
        std::cout << "HERE 1";
        ComputeAxisAlignedLOCALMovingframes(m_sphereMF, mf_LOCSPH);
                std::cout << "HERE 2";
    }
    // Compute the cross producted MF
    ComputeMFtimesMF(m_movingframes, m_CrossProductMF);

    m_session->LoadParameter("Frequency", m_freq, sqrt(2.0));

    // Define TestMaxwellType
    if (m_session->DefinesSolverInfo("TESTMAXWELLTYPE"))
    {
        std::string TestMaxwellTypeStr =
            m_session->GetSolverInfo("TESTMAXWELLTYPE");
        for (int i = 0; i < (int)SIZE_TestMaxwellType; ++i)
        {
            if (TestMaxwellTypeMap[i] == TestMaxwellTypeStr)
            {
                m_TestMaxwellType = (TestMaxwellType)i;
                break;
            }
        }
    }

    else
    {
        m_TestMaxwellType = (TestMaxwellType)0;
    }

    // Define Polarization
    if (m_session->DefinesSolverInfo("POLTYPE"))
    {
        std::string PolTypeStr = m_session->GetSolverInfo("POLTYPE");
        for (int i = 0; i < (int)SIZE_PolType; ++i)
        {
            if (PolTypeMap[i] == PolTypeStr)
            {
                m_PolType = (PolType)i;
                break;
            }
        }
    }
    else
    {
        m_PolType = (PolType)0;
    }

    // Define Incident wave Type
    if (m_session->DefinesSolverInfo("INCTYPE"))
    {
        std::string IncTypeStr = m_session->GetSolverInfo("INCTYPE");
        for (int i = 0; i < (int)SIZE_IncType; ++i)
        {
            if (IncTypeMap[i] == IncTypeStr)
            {
                m_IncType = (IncType)i;
                break;
            }
        }
    }
    else
    {
        m_IncType = (IncType)0;
    }

    // Define Source Type
    if (m_session->DefinesSolverInfo("SOURCETYPE"))
    {
        std::string SourceTypeStr = m_session->GetSolverInfo("SOURCETYPE");
        for (int i = 0; i < (int)SIZE_SourceType; ++i)
        {
            if (SourceTypeMap[i] == SourceTypeStr)
            {
                m_SourceType = (SourceType)i;
                break;
            }
        }
    }
    else
    {
        m_SourceType = (SourceType)0;
    }

    Array<OneD, Array<OneD, NekDouble>>m_velocity(m_spacedim);
    for (int k = 0; k < m_spacedim; ++k)
    {
        m_velocity[k] = Array<OneD, NekDouble>(nq);

        Vmath::Vadd(nq, &m_movingframes[0][k*nq], 1, &m_movingframes[1][k*nq], 1, &m_velocity[k][0], 1);
    }

    ComputeSphericalVector(m_SphericalVector);

    CheckMeshErr(m_movingframes, m_velocity);

    // Compute n_timesMFFwd and m_times_timesMFFwd
    ComputentimesMF(m_movingframes, m_ntimesMFFwd, m_ntimesMFBwd,
                    m_ntimes_ntimesMFFwd, m_ntimes_ntimesMFBwd);

    // Compute Covariant Curl of MF
    ComputeCurlMF(SolverUtils::eCovariant, m_movingframes, m_CurlMF);

    // Compute vaepsilon and mu vector (m_epsveci, m_muvec0);
    m_epsvec = Array<OneD, Array<OneD, NekDouble>>(m_spacedim);
    m_muvec  = Array<OneD, Array<OneD, NekDouble>>(m_spacedim);
    for (int k = 0; k < m_spacedim; ++k)
    {
        m_epsvec[k] = Array<OneD, NekDouble>(nq, 1.0);
        m_muvec[k]  = Array<OneD, NekDouble>(nq, 1.0);
    }

    // InitializeCloakVariables();

    // Compute m_Zim and m_Yim
    // ZimFwd = sqrt( muFwd / epsFwd),  ZimBwd = sqrt( muBwd /
    // epsBwd)
    int nTraceNumPoints = GetTraceNpoints();

    m_ZimFwd = Array<OneD, Array<OneD, NekDouble>>(m_shapedim);
    m_ZimBwd = Array<OneD, Array<OneD, NekDouble>>(m_shapedim);
    m_YimFwd = Array<OneD, Array<OneD, NekDouble>>(m_shapedim);
    m_YimBwd = Array<OneD, Array<OneD, NekDouble>>(m_shapedim);

    for (int j = 0; j < m_shapedim; ++j)
    {
        m_ZimFwd[j] = Array<OneD, NekDouble>(nTraceNumPoints, 1.0);
        m_ZimBwd[j] = Array<OneD, NekDouble>(nTraceNumPoints, 1.0);
        m_YimFwd[j] = Array<OneD, NekDouble>(nTraceNumPoints, 1.0);
        m_YimBwd[j] = Array<OneD, NekDouble>(nTraceNumPoints, 1.0);
    }

    ComputeZimYim(m_epsvec, m_muvec, m_ZimFwd, m_ZimBwd, m_YimFwd, m_YimBwd);

    // Compute m_epsvecminus1 and m_muminus1
    m_negepsvecminus1 = Array<OneD, Array<OneD, NekDouble>>(m_spacedim);
    m_negmuvecminus1  = Array<OneD, Array<OneD, NekDouble>>(m_spacedim);
    for (int k = 0; k < m_spacedim; ++k)
    {
        m_negepsvecminus1[k] = Array<OneD, NekDouble>(nq, 0.0);
        m_negmuvecminus1[k]  = Array<OneD, NekDouble>(nq, 0.0);

        if (!m_NoInc)
        {
            Vmath::Sadd(nq, -1.0, m_muvec[k], 1, m_negmuvecminus1[k], 1);
            Vmath::Sadd(nq, -1.0, m_epsvec[k], 1, m_negepsvecminus1[k], 1);

            Vmath::Neg(nq, m_negmuvecminus1[k], 1);
            Vmath::Neg(nq, m_negepsvecminus1[k], 1);
        }
    }

    // Compute de^m/dt \cdot e^k
    if (m_AddRotation)
    {
        m_coriolis = Array<OneD, NekDouble>(nq);
        m_coriolis = EvaluateCoriolis();

       // Computedemdxicdote(m_movingframes, m_dedxi_cdot_e);
    }

    // Generate Sigma Block with thicknes of m_PMLthickness and m_PMLmax
    if (m_AddPML)
    {
        GenerateSigmaPML(m_PMLthickness, m_PMLstart, m_PMLmaxsigma, m_SigmaPML);
    }

    // If explicit it computes RHS and PROJECTION for the time integration
    if (m_explicitAdvection)
    {
        m_ode.DefineOdeRhs(&MMFMaxwell::DoOdeRhs, this);
        m_ode.DefineProjection(&MMFMaxwell::DoOdeProjection, this);
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
MMFMaxwell::~MMFMaxwell()
{
}

void MMFMaxwell::v_DoSolve()
{
    ASSERTL0(m_intScheme != 0, "No time integration scheme.");

    int i, nchk = 1;
    int nq         = GetTotPoints();
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
    Array<OneD, Array<OneD, NekDouble>> tmp(nvariables);

    // Order storage to list time-integrated fields first.
    for (i = 0; i < nvariables; ++i)
    {
        fields[i] = m_fields[m_intVariables[i]]->GetPhys();
        m_fields[m_intVariables[i]]->SetPhysState(false);
    }

    // Compute Initial Energy
    m_Energy0 = ComputeEnergyDensity(fields);

    // Initialise time integration scheme
    m_intScheme->InitializeScheme(m_timestep, fields, m_time, m_ode);

    // Check uniqueness of checkpoint output
    ASSERTL0((m_checktime == 0.0 && m_checksteps == 0) ||
                 (m_checktime > 0.0 && m_checksteps == 0) ||
                 (m_checktime == 0.0 && m_checksteps > 0),
             "Only one of IO_CheckTime and IO_CheckSteps "
             "should be set!");

    int Ntot = m_steps / m_checksteps + 1;

    Array<OneD, NekDouble> TimeSeries(Ntot);
    Array<OneD, NekDouble> Energy(Ntot);

    LibUtilities::Timer timer;
    bool doCheckTime  = false;
    int step          = 0;
    NekDouble intTime = 0.0;
    NekDouble cpuTime = 0.0;
    NekDouble elapsed = 0.0;

    int cntap = 0;
    Array<OneD, NekDouble> Ezantipod;
    int indxantipod = 0;

    switch (m_SourceType)
    {
        case ePointSource:
        {
            Ezantipod = Array<OneD, NekDouble>(m_steps / m_checksteps);

            Array<OneD, NekDouble> x(nq);
            Array<OneD, NekDouble> y(nq);
            Array<OneD, NekDouble> z(nq);

            m_fields[0]->GetCoords(x, y, z);

            NekDouble Tol = 0.000001;
            NekDouble rad;
            for (i = 0; i < nq; ++i)
            {
                rad = sqrt((x[i] + m_Psx) * (x[i] + m_Psx) +
                           (y[i] + m_Psy) * (y[i] + m_Psy) +
                           (z[i] + m_Psz) * (z[i] + m_Psz));
                if (rad < Tol)
                {
                    indxantipod = i;
                    break;
                }
            }
        }
        break;

        case ePlanarSource:
        {
            m_SourceVector = Array<OneD, NekDouble>(nq, 0.0);

            Array<OneD, NekDouble> x(nq);
            Array<OneD, NekDouble> y(nq);
            Array<OneD, NekDouble> z(nq);

            m_fields[0]->GetCoords(x, y, z);

            NekDouble Tol = 0.000001;
            NekDouble rad;
            for (i = 0; i < nq; ++i)
            {
                rad = sqrt((x[i] - m_Psx) * (x[i] - m_Psx));
                if (rad < Tol)
                {
                    m_SourceVector[i] = 1.0;
                }
            }

            std::cout << "*** Area of Planar Source = "
                      << m_fields[0]->PhysIntegral(m_SourceVector) << std::endl;
        }
        break;

        default:
            break;
    }

    int cntpml = 0;
    int P1indx = 0, P2indx = 0, P3indx = 0;
    Array<OneD, NekDouble> P1;
    Array<OneD, NekDouble> P2;
    Array<OneD, NekDouble> P3;
    if (m_TestPML)
    {
        P1 = Array<OneD, NekDouble>(m_steps / m_checksteps);
        P2 = Array<OneD, NekDouble>(m_steps / m_checksteps);
        P3 = Array<OneD, NekDouble>(m_steps / m_checksteps);

        Array<OneD, NekDouble> x(nq);
        Array<OneD, NekDouble> y(nq);
        Array<OneD, NekDouble> z(nq);

        m_fields[0]->GetCoords(x, y, z);

        NekDouble Tol = 0.000001;
        NekDouble rad;
        for (int i = 0; i < nq; ++i)
        {
            rad = sqrt((x[i] + 3.0) * (x[i] + 3.0) + (y[i]) * (y[i]));

            if (rad < Tol)
            {
                P1indx = i;
                break;
            }
        }

        for (int i = 0; i < nq; ++i)
        {
            rad =
                sqrt((x[i] + 3.0) * (x[i] + 3.0) + (y[i] - 1.5) * (y[i] - 1.5));
            if (rad < Tol)
            {
                P2indx = i;
                break;
            }
        }

        for (int i = 0; i < nq; ++i)
        {
            rad =
                sqrt((x[i] + 3.0) * (x[i] + 3.0) + (y[i] - 3.0) * (y[i] - 3.0));
            if (rad < Tol)
            {
                P3indx = i;
                break;
            }
        }
    }

    int indx;
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
            std::cout << "Steps: " << std::setw(8) << std::left << step + 1 << " "
                      << "Time: " << std::setw(12) << std::left << m_time;

            std::stringstream ss;
            ss << cpuTime / 60.0 << " min.";
            std::cout << " CPU Time: " << std::setw(8) << std::left << ss.str()
                      << std::endl;

            cpuTime = 0.0;
        }

        // switch (m_SourceType)
        // {
        //     case ePointSource:
        //     {
        //         if (m_time <= m_PSduration)
        //         {
        //             Array<OneD, NekDouble> Impulse(nq);
        //             Impulse = GaussianPulse(m_time, m_Psx, m_Psy, m_Psz,
        //                                     m_Gaussianradius);
        //             Vmath::Vadd(nq, &Impulse[0], 1,
        //                         &fields[m_intVariables[2]][0], 1,
        //                         &fields[m_intVariables[2]][0], 1);
        //         }
        //     }
        //     break;

        //     case ePlanarSource:
        //     {
        //         Array<OneD, NekDouble> Impulse(nq);
        //         for (int i = 0; i < 3; ++i)
        //         {
        //             Impulse = GetIncidentField(i, m_time);
        //             Vmath::Vmul(nq, m_SourceVector, 1, Impulse, 1, Impulse, 1);
        //             Vmath::Vadd(nq, &Impulse[0], 1,
        //                         &fields[m_intVariables[i]][0], 1,
        //                         &fields[m_intVariables[i]][0], 1);
        //         }
        //     }
        //     break;

        //     default:
        //         break;
        // }

        // Transform data into coefficient space
        for (i = 0; i < nvariables; ++i)
        {
            m_fields[m_intVariables[i]]->SetPhys(fields[i]);
            m_fields[m_intVariables[i]]->FwdTransLocalElmt(
                fields[i], m_fields[m_intVariables[i]]->UpdateCoeffs());
            m_fields[m_intVariables[i]]->SetPhysState(false);
        }
        // for (i = 0; i < nq; ++i)
        //  std::cout << m_fields[0][0][i] <<std::endl;

        // Write out checkpoint files
        if ((m_checksteps && step && !((step + 1) % m_checksteps)) ||
            doCheckTime)
        {
            indx             = (step + 1) / m_checksteps;
            TimeSeries[indx] = m_time;

            // if (m_TestMaxwellType == eScatField2D)
            // {
            //     Checkpoint_TotalFieldOutput(nchk, m_time, fields);
            //     Checkpoint_TotPlotOutput(nchk, m_time, fields);
            // }
            Checkpoint_PlotOutput(nchk, fields);
            // Checkpoint_EDFluxOutput(nchk, m_time, fields);
            // Checkpoint_EnergyOutput(nchk, m_time, fields);
            Checkpoint_Output(nchk++);

            Energy[indx] = ComputeEnergyDensity(fields)-m_Energy0;

            std::cout << "|EHr|: F1 = " << RootMeanSquare(fields[0])
                      << ", F2 = " << RootMeanSquare(fields[1])
                      << ", F3 = " << RootMeanSquare(fields[2])
                      << ", Energy = " << Energy[indx] << std::endl << std::endl;

            if (m_SourceType == ePointSource)
            {
                Ezantipod[cntap++] = fields[2][indxantipod];
            }

            if (m_TestPML)
            {
                P1[cntpml] = fields[2][P1indx];
                P2[cntpml] = fields[2][P2indx];
                P3[cntpml] = fields[2][P3indx];
                cntpml++;
            }
            doCheckTime = false;
        }

        // Step advance
        ++step;
    }

    // Print out summary statistics
    if (m_session->GetComm()->GetRank() == 0)
    {
        std::cout << "Time-integration  : " << intTime << "s" << std::endl;

        std::cout << "TimeSeries = " << std::endl;
        for (int i = 0; i < m_steps / m_checksteps; ++i)
        {
            std::cout << TimeSeries[i] << ", ";
        }
        std::cout << std::endl << std::endl;

        std::cout << "Energy Density = " << std::endl;
        for (int i = 0; i < m_steps / m_checksteps; ++i)
        {
            std::cout << Energy[i] << ", ";
        }
        std::cout << std::endl << std::endl;

        // Print out total Error
        NekDouble L2err, TotalL2err=0.0;
         Array<OneD, NekDouble> exactsoln(nq);
         for (int i=0; i<nvariables; ++i)
        {
            v_EvaluateExactSolution(i, exactsoln, m_time);

            L2err = m_fields[i]->L2(m_fields[i]->GetPhys(), exactsoln);
            TotalL2err += L2err*L2err;
        }
        TotalL2err = sqrt(TotalL2err/3.0);

        std::cout << "Total L2 error = " << TotalL2err << std::endl;

        if (m_PrintoutSurfaceCurrent)
        {
            Printout_SurfaceCurrent(fields, m_time);
        }

        if (m_SourceType == ePointSource)
        {
            std::cout << "Ez at antipod = " << std::endl;
            for (int i = 0; i < m_steps / m_checksteps; ++i)
            {
                std::cout << Ezantipod[i] << ", ";
            }
            std::cout << std::endl << std::endl;
        }
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

/**
 * @brief Compute the right-hand side for the linear advection equation.
 *
 * @param inarray    Given fields.
 * @param outarray   Calculated solution.
 * @param time       Time.
 */
void MMFMaxwell::DoOdeRhs(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time)
{
    int i;
    int nvar    = inarray.size();
    int ncoeffs = GetNcoeffs();
    int nq      = GetTotPoints();

    Array<OneD, Array<OneD, NekDouble>> physarray(nvar);
    Array<OneD, Array<OneD, NekDouble>> modarray(nvar);
    for (i = 0; i < nvar; ++i)
    {
        physarray[i] = Array<OneD, NekDouble>(nq);
        modarray[i]  = Array<OneD, NekDouble>(ncoeffs, 0.0);

        Vmath::Vcopy(nq, &inarray[i][0], 1, &physarray[i][0], 1);
    }

    for (i = 0; i < nvar; i++)
    {
        m_fields[i]->SetPhysState(true);
    }

    // Compute Curl
    switch (m_TestMaxwellType)
    {
        case eScatField2D:
        case eTotField2D:

        {
            // Imaginary part is computed the same as Real part
            Array<OneD, Array<OneD, NekDouble>> tmpin(3);
            Array<OneD, Array<OneD, NekDouble>> tmpout(3);

            for (int i = 0; i < 3; ++i)
            {
                tmpin[i]  = Array<OneD, NekDouble>(nq);
                tmpout[i] = Array<OneD, NekDouble>(ncoeffs, 0.0);

                Vmath::Vcopy(nq, &physarray[i][0], 1, &tmpin[i][0], 1);
            }

            WeakDGMaxwellDirDeriv(physarray, tmpout, time);

            for (int i = 0; i < 3; ++i)
            {
                // For E and H
                Vmath::Vcopy(ncoeffs, &tmpout[i][0], 1, &modarray[i][0], 1);
            }
        }
        break;

        default:
        {
            WeakDGMaxwellDirDeriv(physarray, modarray, time);
        }
        break;
    }

    for (i = 0; i < nvar; ++i)
    {
        m_fields[i]->MultiplyByElmtInvMass(modarray[i], modarray[i]);
        m_fields[i]->BwdTrans(modarray[i], outarray[i]);
    }

    if (m_TestMaxwellType == eMaxwellSphere)
    {
        Array<OneD, NekDouble> F(nq);
        for (int j = 0; j < 2; ++j)
        {
            F = TestMaxwellSphere(time, m_freq, 3 + j);
            Vmath::Vadd(nq, &F[0], 1, &outarray[j][0], 1, &outarray[j][0], 1);
        }
    }

    // Add Absorbing Boundary Conditions
    if (m_AddPML > 0)
    {
        AddPML(physarray, outarray);
    }

    // Add dedt component
    if (m_AddRotation)
    {
        AddCoriolis(physarray, outarray);
        AdddedtMaxwell(physarray, outarray);
    }

    // Divide it by varepsilon or mu
    Array<OneD, NekDouble> dFdt(nq, 0.0);
    switch (m_TestMaxwellType)
    {
        case eMaxwell1D:
        {
            Vmath::Vdiv(nq, outarray[0], 1, m_epsvec[0], 1, outarray[0], 1);
            Vmath::Vdiv(nq, outarray[1], 1, m_muvec[0], 1, outarray[1], 1);
        }
        break;

        // TO BE CHANGED
        case eTestMaxwell2DPECAVGFLUX:
        {
            Array<OneD, NekDouble> Hxdt(nq, 0.0);
            Array<OneD, NekDouble> Hydt(nq, 0.0);

            Hxdt = TestMaxwell2DPEC(time, 10, m_PolType);
            Hydt = TestMaxwell2DPEC(time, 11, m_PolType);

            Array<OneD, NekDouble> x0(nq);
            Array<OneD, NekDouble> x1(nq);
            Array<OneD, NekDouble> x2(nq);

            m_fields[0]->GetCoords(x0, x1, x2);

            NekDouble theta, tmpx, tmpy;
            NekDouble uxx, uxy, uyy, detu, uti, uri;
            Array<OneD, NekDouble> utvec(nq, 1.0);
            Array<OneD, NekDouble> urvec(nq, 1.0);
            Array<OneD, NekDouble> tmpIN(nq);

            // Case I: ut = 4.0, ur = 0.5
            // NekDouble ut=4.0;
            // NekDouble ur=0.5;

            //  m_fields[0]->GenerateElementVector(m_ElemtGroup1, ut, 1.0,
            //  utvec);
            //    m_fields[0]->GenerateElementVector(m_ElemtGroup1, ur, 1.0,
            //    urvec);

            // Case II: ut = 0.5, ur = 1 - 2x^2
            NekDouble ut = 0.5;
            m_fields[0]->GenerateElementVector(m_ElemtGroup1, ut, 1.0, utvec);

            m_fields[0]->GenerateElementVector(m_ElemtGroup1, 1.0, 0.0, tmpIN);

            for (int i = 0; i < nq; i++)
            {
                urvec[i] =
                    tmpIN[i] * (1.0 - 2 * x0[i] * x0[i]) + (1.0 - tmpIN[i]);
            }

            for (int i = 0; i < nq; ++i)
            {
                theta = atan2((x1[i] + 2.0), (x0[i] + 2.0));

                uti = utvec[i];
                uri = urvec[i];

                uxx = uti * cos(theta) * cos(theta) +
                      uri * sin(theta) * sin(theta);
                uyy = uti * sin(theta) * sin(theta) +
                      uri * cos(theta) * cos(theta);
                uxy = (uti - uri) * cos(theta) * sin(theta);

                detu = uxx * uyy - uxy * uxy;

                tmpx = outarray[0][i] + (1.0 - uxx) * Hxdt[i] - uxy * Hydt[i];
                tmpy = outarray[1][i] - uxy * Hxdt[i] + (1.0 - uyy) * Hydt[i];

                outarray[0][i] = (1 / detu) * (uyy * tmpx - uxy * tmpy);
                outarray[1][i] = (1 / detu) * (-uxy * tmpx + uxx * tmpy);
            }
        }
        break;

        case eTestMaxwell2DPEC:
        case eScatField2D:
        case eTotField2D:
        {
            switch (m_PolType)
            {
                case eTransMagnetic:
                {
                    if (m_TestMaxwellType == eTestMaxwell2DPEC)
                    {
                        dFdt = TestMaxwell2DPEC(time, 10, m_PolType);
                        Vmath::Vvtvp(nq, m_negmuvecminus1[0], 1, dFdt, 1,
                                     outarray[0], 1, outarray[0], 1);

                        dFdt = TestMaxwell2DPEC(time, 11, m_PolType);
                        Vmath::Vvtvp(nq, m_negmuvecminus1[1], 1, dFdt, 1,
                                     outarray[1], 1, outarray[1], 1);

                        dFdt = TestMaxwell2DPEC(time, 12, m_PolType);
                        Vmath::Vvtvp(nq, m_negepsvecminus1[2], 1, dFdt, 1,
                                     outarray[2], 1, outarray[2], 1);
                    }

                    if (m_TestMaxwellType == eScatField2D)
                    {
                        dFdt = GetIncidentField(10, time);
                        Vmath::Vvtvp(nq, m_negmuvecminus1[0], 1, dFdt, 1,
                                     outarray[0], 1, outarray[0], 1);

                        dFdt = GetIncidentField(11, time);
                        Vmath::Vvtvp(nq, m_negmuvecminus1[1], 1, dFdt, 1,
                                     outarray[1], 1, outarray[1], 1);

                        dFdt = GetIncidentField(12, time);
                        Vmath::Vvtvp(nq, m_negepsvecminus1[2], 1, dFdt, 1,
                                     outarray[2], 1, outarray[2], 1);
                    }

                    Vmath::Vdiv(nq, outarray[0], 1, m_muvec[0], 1, outarray[0],
                                1);
                    Vmath::Vdiv(nq, outarray[1], 1, m_muvec[1], 1, outarray[1],
                                1);
                    Vmath::Vdiv(nq, outarray[2], 1, m_epsvec[2], 1, outarray[2],
                                1);
                }
                break;

                case eTransElectric:
                {
                    if (m_TestMaxwellType == eTestMaxwell2DPEC)
                    {
                        // (I - \mu^i) d F^{inc} / dt
                        dFdt = TestMaxwell2DPEC(time, 10, m_PolType);
                        Vmath::Vvtvp(nq, m_negepsvecminus1[0], 1, dFdt, 1,
                                     outarray[0], 1, outarray[0], 1);

                        dFdt = TestMaxwell2DPEC(time, 11, m_PolType);
                        Vmath::Vvtvp(nq, m_negepsvecminus1[1], 1, dFdt, 1,
                                     outarray[1], 1, outarray[1], 1);

                        dFdt = TestMaxwell2DPEC(time, 12, m_PolType);
                        Vmath::Vvtvp(nq, m_negmuvecminus1[2], 1, dFdt, 1,
                                     outarray[2], 1, outarray[2], 1);
                    }

                    if (m_TestMaxwellType == eScatField2D)
                    {
                        dFdt = GetIncidentField(10, time);
                        Vmath::Vvtvp(nq, m_negepsvecminus1[0], 1, dFdt, 1,
                                     outarray[0], 1, outarray[0], 1);

                        // Add - wp^2 \int E_2^{inc}
                        dFdt = GetIncidentField(21, time);
                        if (m_DispersiveCloak)
                        {
                            Vmath::Vmul(nq, m_wp2, 1, dFdt, 1, dFdt, 1);
                            Vmath::Vsub(nq, outarray[1], 1, dFdt, 1,
                                        outarray[1], 1);
                        }

                        else
                        {
                            Vmath::Vvtvp(nq, m_negepsvecminus1[1], 1, dFdt, 1,
                                         outarray[1], 1, outarray[1], 1);
                        }

                        dFdt = GetIncidentField(12, time);
                        Vmath::Vvtvp(nq, m_negmuvecminus1[2], 1, dFdt, 1,
                                     outarray[2], 1, outarray[2], 1);
                    }

                    Vmath::Vdiv(nq, outarray[0], 1, m_epsvec[0], 1, outarray[0],
                                1);
                    Vmath::Vdiv(nq, outarray[1], 1, m_epsvec[1], 1, outarray[1],
                                1);
                    Vmath::Vdiv(nq, outarray[2], 1, m_muvec[2], 1, outarray[2],
                                1);
                }
                break;

                default:
                    break;
            }
        }
        break; //	  case SolverUtils::eTestMaxwell2DPEC:

        default:
            break;

    } // switch(m_TestMaxwellType)
}

Array<OneD, NekDouble> MMFMaxwell::GreenDerivCompensate(const int var,
    const Array<OneD, const Array<OneD, NekDouble>> &physarray)
{
    // routine works for both primitive and conservative formulations
    int nq              = GetNpoints();
    int ncoeffs         = GetNcoeffs();

    Array<OneD, NekDouble> tmp(nq);
    Array<OneD, NekDouble> outarray(ncoeffs);

    Array<OneD, Array<OneD, NekDouble>> fluxvector(m_shapedim);
    for (int j = 0; j < m_shapedim; ++j)
    {
        fluxvector[j] = Array<OneD, NekDouble>(nq);
    }

    // TM: (-E3,0), (0,E3), (-H1, H2)
    // TE: (H3,0), (0,-H3), (E1, -E2)
    GetMaxwellFluxVector(var, physarray, fluxvector);

    // CurlMF[0] = e^3 \cdot (\nabla \times e^1)
    // CurlMF[1] = e^3 \cdot (\nabla \times e^2)
    // CurlMF[2] = e^1 \cdot (\nabla \times e^3)
    // CurlMF[3] = e^2 \cdot (\nabla \times e^3)
    switch(var)
    {
        case 0:
        {
            Vmath::Vmul(nq, &fluxvector[0][0], 1, &m_CurlMF[0][0], 1, &tmp[0], 1);
            Vmath::Neg(nq, tmp, 1);
        }
        break;

        case 1:
        {
            Vmath::Vmul(nq, &fluxvector[1][0], 1, &m_CurlMF[1][0], 1, &tmp[0], 1);
            Vmath::Neg(nq, tmp, 1);
        }
        break;

        case 2:
        {
            Vmath::Vmul(nq, &fluxvector[0][0], 1, &m_CurlMF[2][0], 1, &tmp[0], 1);
            Vmath::Neg(nq, tmp, 1);
            Vmath::Vvtvm(nq, &fluxvector[1][0], 1, &m_CurlMF[3][0], 1, &tmp[0], 1,
                        &tmp[0], 1);
        }
        break;

        default:
          break;
    }

    m_fields[0]->IProductWRTBase(tmp, outarray);

    return outarray;
}


Array<OneD, NekDouble> MMFMaxwell::ComputeSDMaxwell(const int var,
    const Array<OneD, const Array<OneD, NekDouble>> &physarray)
{
    int nq              = GetNpoints();
    int ncoeffs         = GetNcoeffs();
 
    Array<OneD, NekDouble> tmpc(ncoeffs);
    Array<OneD, NekDouble> tmp(nq);
    Array<OneD, NekDouble> outarray(ncoeffs,0.0);

    Array<OneD, NekDouble> velvector(m_spacedim * nq);

    Array<OneD, Array<OneD, NekDouble>> fluxvector(m_shapedim);
    for (int j = 0; j < m_shapedim; ++j)
    {
        fluxvector[j] = Array<OneD, NekDouble>(nq);
    }

    GetMaxwellFluxVector(var, physarray, fluxvector);

    for (int k = 0; k < m_spacedim; ++k)
    {
        Vmath::Vmul(nq, &fluxvector[0][0], 1, &m_CrossProductMF[0][k*nq], 1, &velvector[k * nq], 1);
    }

    switch(m_DivergenceRestore)
    {
        case 1:
        {
            tmp = ComputeSpuriousDivergence(m_movingframes, m_SphericalVector[2], velvector);
        }
        break;

        case 2:
        {
            // m_movingframes = LOCALsphere, mf_LOCAL = LOCAL
              tmp = ComputeSpuriousDivergence(mf_LOCSPH, m_movingframes[2], velvector);
              // tmp = ComputeSpuriousDivergence(m_movingframes, mf_LOCAL[2], velvector);
              Vmath::Neg(nq, tmp, 1);
        }
        break;

        case 3:
        {
            // tmp = ComputeSpuriousDivergence(mf_LOCSPH, m_SphericalVector[2], velvector);
            tmp = ComputeSpuriousDivergence(m_movingframes, m_movingframes[2], velvector);
            Vmath::Neg(nq, tmp, 1);
        }
        break;

        default:
        break;
    }

    m_fields[0]->IProductWRTBase(tmp, outarray);

    return outarray;
}
    // switch(var)
    // {
    //     case 0:
    //     {
    //         for (int k = 0; k < m_spacedim; ++k)
    //         {
    //             Vmath::Vmul(nq, &fluxvector[0][0], 1, &m_CrossProductMF[0][k*nq], 1, &velvector[k * nq], 1);
    //         }
    //             tmp = ComputeSurfaceDiv(m_movingframes[2], velvector);
    //             m_fields[0]->IProductWRTBase(tmp, outarray);
    //     }
    //     break;

    //     case 1:
    //     {
    //         for (int k = 0; k < m_spacedim; ++k)
    //         {
    //             Vmath::Vmul(nq, &fluxvector[1][0], 1, &m_CrossProductMF[1][k*nq], 1, &velvector[k * nq], 1);
    //         }
    //             tmp = ComputeSurfaceDiv(m_movingframes[2], velvector);
    //             m_fields[0]->IProductWRTBase(tmp, outarray);
    //     }
    //     break;

    //     case 2:
    //     {
    //         for (int j=0; j < m_shapedim; ++j)
    //         {
    //             for (int k = 0; k < m_spacedim; ++k)
    //             {
    //                 Vmath::Vmul(nq, &fluxvector[j][0], 1, &m_CrossProductMF[j][k*nq], 1, &velvector[k * nq], 1);
    //             }
                
    //             tmp = ComputeSurfaceDiv(m_movingframes[2], velvector);
    //             if(j==0)
    //             {
    //                 Vmath::Neg(nq, tmp, 1);  
    //             }
    //             m_fields[0]->IProductWRTBase(tmp, tmpc);

    //             Vmath::Vadd(ncoeffs, tmpc, 1, outarray, 1, outarray, 1);
    //         }
    //     }
    //     break;

    //     default:
    //     break;    
    // }


/**
 * @brief Calculate weak DG advection in the form \f$ \langle\phi,
 * \hat{F}\cdot n\rangle - (\nabla \phi \cdot F) \f$
 *
 * @param   InField         Fields.
 * @param   OutField        Storage for result.
 * @param   NumericalFluxIncludesNormal     Default: true.
 * @param   InFieldIsPhysSpace              Default: false.
 * @param   nvariables      Number of fields.
 */
void MMFMaxwell::WeakDGMaxwellDirDeriv(
    const Array<OneD, const Array<OneD, NekDouble>> &InField,
    Array<OneD, Array<OneD, NekDouble>> &OutField, const NekDouble time)
{
    int i;
    int nq              = GetNpoints();
    int ncoeffs         = GetNcoeffs();
    int nTracePointsTot = GetTraceNpoints();
    int nvar            = 3;

    Array<OneD, Array<OneD, NekDouble>> fluxvector(m_shapedim);
    for (i = 0; i < m_shapedim; ++i)
    {
        fluxvector[i] = Array<OneD, NekDouble>(nq);
    }

    Array<OneD, Array<OneD, NekDouble>> physfield(nvar);
    for (i = 0; i < nvar; ++i)
    {
        physfield[i] = InField[i];
    }

    Array<OneD, NekDouble> tmpc(ncoeffs);
    for (i = 0; i < nvar; ++i)
    {
        GetMaxwellFluxVector(i, physfield, fluxvector);

        OutField[i] = Array<OneD, NekDouble>(ncoeffs, 0.0);
        for (int j = 0; j < m_shapedim; ++j)
        {
            // Directional derivation with respect to the j'th moving frame
            // tmp_j = ( \nabla \phi, fluxvector[j] \mathbf{e}^j )
            m_fields[i]->IProductWRTDirectionalDerivBase(m_CrossProductMF[j],
                                                         fluxvector[j], tmpc);
            Vmath::Vadd(ncoeffs, &tmpc[0], 1, &OutField[i][0], 1,
                        &OutField[i][0], 1);
        }
    }

    // V the numerical flux and add to the modal coeffs
    // if the NumericalFlux function does not include the
    // normal in the output
    Array<OneD, Array<OneD, NekDouble>> numfluxFwd(nvar);
    Array<OneD, Array<OneD, NekDouble>> numfluxBwd(nvar);
    for (i = 0; i < nvar; ++i)
    {
        numfluxFwd[i] = Array<OneD, NekDouble>(nTracePointsTot, 0.0);
        numfluxBwd[i] = Array<OneD, NekDouble>(nTracePointsTot, 0.0);
    }

    // Evaluate numerical flux in physical space which may in
    // general couple all component of vectors
    NumericalMaxwellFlux(physfield, numfluxFwd, numfluxBwd, time);

    // Evaulate  <\phi, \hat{F}\cdot n> - OutField[i]
    Array<OneD, NekDouble> GreenCompensation(ncoeffs);
    Array<OneD, NekDouble> SpuriousDiv(ncoeffs);
    for (i = 0; i < nvar; ++i)
    {
        Vmath::Neg(ncoeffs, OutField[i], 1);
        m_fields[i]->AddFwdBwdTraceIntegral(numfluxFwd[i], numfluxBwd[i],
                                            OutField[i]);

        GreenCompensation = GreenDerivCompensate(i, physfield);                                       

        // Add GreenCompensation
        Vmath::Vadd(ncoeffs, &GreenCompensation[0], 1, &OutField[i][0], 1, &OutField[i][0], 1);
        m_fields[i]->SetPhysState(false);

        if(m_DivergenceRestore>0)
        {
            SpuriousDiv = ComputeSDMaxwell(i, physfield);
            Vmath::Vadd(ncoeffs, &SpuriousDiv[0], 1, &OutField[i][0], 1, &OutField[i][0], 1);
        }
    }

}


void MMFMaxwell::NumericalMaxwellFlux(
    const Array<OneD, const Array<OneD, NekDouble>> &physfield,
    Array<OneD, Array<OneD, NekDouble>> &numfluxFwd,
    Array<OneD, Array<OneD, NekDouble>> &numfluxBwd, const NekDouble time)
{

    switch (m_TestMaxwellType)
    {
        case eTestMaxwell2DPEC:
        case eTestMaxwell2DPECAVGFLUX:
        case eTestMaxwell2DPMC:
        case eScatField2D:
        case eTotField2D:
        case eMaxwellSphere:
        case eELF2DSurface:
        {
            switch (m_PolType)
            {
                case eTransMagnetic:
                {
                    NumericalMaxwellFluxTM(physfield, numfluxFwd, numfluxBwd, time);
                }
                break;

                case eTransElectric:
                {
                    NumericalMaxwellFluxTE(physfield, numfluxFwd, numfluxBwd, time);
                }
                break;

                default:
                    break;
            }
        }
        break; // eMaxwell2D

        default:
            break;
    } // m_TestMaxwellType
}

void MMFMaxwell::NumericalMaxwellFluxTM(
    const Array<OneD, const Array<OneD, NekDouble>> &physfield,
    Array<OneD, Array<OneD, NekDouble>> &numfluxFwd,
    Array<OneD, Array<OneD, NekDouble>> &numfluxBwd, const NekDouble time)
{
    int nq              = GetNpoints();
    int nTraceNumPoints = GetTraceTotPoints();
    int nvar            = 3;

    // get temporary arrays
    Array<OneD, Array<OneD, NekDouble>> Fwd(nvar);
    Array<OneD, Array<OneD, NekDouble>> Bwd(nvar);
    for (int i = 0; i < nvar; ++i)
    {
        Fwd[i] = Array<OneD, NekDouble>(nTraceNumPoints, 0.0);
        Bwd[i] = Array<OneD, NekDouble>(nTraceNumPoints, 0.0);

        // get the physical values at the trace
        m_fields[i]->GetFwdBwdTracePhys(physfield[i], Fwd[i], Bwd[i]);
    }

    // E^|| = 0 at the PEC boundaries vs. H^|| = 0 at PMC boundaries
    Array<OneD, NekDouble> IncField(nq, 0.0);
    Array<OneD, NekDouble> IncFieldBwd(nTraceNumPoints, 0.0);
    Array<OneD, Array<OneD, NekDouble>> IncFieldFwd(m_spacedim);
    for (int i = 0; i < m_spacedim; ++i)
    {
        IncFieldFwd[i] = Array<OneD, NekDouble>(nTraceNumPoints, 0.0);

        IncField = GetIncidentField(i, time);
        m_fields[0]->GetFwdBwdTracePhys(IncField, IncFieldFwd[i], IncFieldBwd);

        Vmath::Svtvp(nTraceNumPoints, 2.0, &IncFieldFwd[i][0], 1, &Fwd[i][0], 1,
                     &IncFieldFwd[i][0], 1);
    }

    // Total Field Formulation
    // CopyBoundaryTrace(Fwd[0], Bwd[0], SolverUtils::eFwdEQBwd, 0,
    //                   SpatialDomains::ePEC_Forces);
    // CopyBoundaryTrace(Fwd[1], Bwd[1], SolverUtils::eFwdEQBwd, 1,
    //                   SpatialDomains::ePEC_Forces);
    // CopyBoundaryTrace(IncFieldFwd[2], Bwd[2], SolverUtils::eFwdEQNegBwd, 2,
    //                   SpatialDomains::ePEC_Forces);

    // CopyBoundaryTrace(IncFieldFwd[0], Bwd[0], SolverUtils::eFwdEQNegBwd, 0,
    //                   SpatialDomains::ePMC_Forces);
    // CopyBoundaryTrace(IncFieldFwd[1], Bwd[1], SolverUtils::eFwdEQNegBwd, 1,
    //                   SpatialDomains::ePMC_Forces);
    // CopyBoundaryTrace(Fwd[2], Bwd[2], SolverUtils::eFwdEQBwd, 2,
    //                   SpatialDomains::ePMC_Forces);

    // CopyBoundaryTrace(Fwd[0], Bwd[0], SolverUtils::eFwdEQBwd, 0,
    //                   SpatialDomains::ePEC);
    // CopyBoundaryTrace(Fwd[1], Bwd[1], SolverUtils::eFwdEQBwd, 1,
    //                   SpatialDomains::ePEC);
    // CopyBoundaryTrace(Fwd[2], Bwd[2], SolverUtils::eFwdEQNegBwd, 2,
    //                   SpatialDomains::ePEC);

    // CopyBoundaryTrace(Fwd[0], Bwd[0], SolverUtils::eFwdEQNegBwd, 0,
    //                   SpatialDomains::ePMC);
    // CopyBoundaryTrace(Fwd[1], Bwd[1], SolverUtils::eFwdEQNegBwd, 1,
    //                   SpatialDomains::ePMC);
    // CopyBoundaryTrace(Fwd[2], Bwd[2], SolverUtils::eFwdEQBwd, 2,
    //                   SpatialDomains::ePMC);

    Array<OneD, NekDouble> e1Fwd_cdot_ncrossdH(nTraceNumPoints, 0.0);
    Array<OneD, NekDouble> e1Bwd_cdot_ncrossdH(nTraceNumPoints, 0.0);
    Array<OneD, NekDouble> e2Fwd_cdot_ncrossdH(nTraceNumPoints, 0.0);
    Array<OneD, NekDouble> e2Bwd_cdot_ncrossdH(nTraceNumPoints, 0.0);
    Array<OneD, NekDouble> e3Fwd_cdot_dEe3(nTraceNumPoints, 0.0);
    Array<OneD, NekDouble> e3Bwd_cdot_dEe3(nTraceNumPoints, 0.0);

    // Compute  numfluxFwd[dir] = (eFwd^[dir] \cdot n \times e^3) * (YimFwd *
    // EFwd^3 + YimBwd * EBwd^3 )
    ComputeNtimesFz(0, Fwd, Bwd, m_YimFwd[0], m_YimBwd[0], numfluxFwd[0],
                    numfluxBwd[0]);
    ComputeNtimesFz(1, Fwd, Bwd, m_YimFwd[1], m_YimBwd[1], numfluxFwd[1],
                    numfluxBwd[1]);

    // Compute numfluxFwd[2] = eFwd^3 \cdot ( n1e1 \times ( ZimFwd HFwd + ZimBwd
    // HBwd ) ) / 2 {{Z_i}}
    ComputeNtimesF12(Fwd, Bwd, m_ZimFwd[0], m_ZimBwd[0], m_ZimFwd[1],
                     m_ZimBwd[1], numfluxFwd[2], numfluxBwd[2]);

    // Compute e1Fwd_cdot_ncrossdE = eFwd[dir] \cdot \alpha n \times n \times
    // [H] / 2 {{YimFwd}}
    ComputeNtimestimesdFz(0, Fwd, Bwd, m_YimFwd[0], m_YimBwd[0],
                          e1Fwd_cdot_ncrossdH, e1Bwd_cdot_ncrossdH);
    ComputeNtimestimesdFz(1, Fwd, Bwd, m_YimFwd[1], m_YimBwd[1],
                          e2Fwd_cdot_ncrossdH, e2Bwd_cdot_ncrossdH);

    // Compute  \alpha [E3] * ( 1/2{{Zim1}} + 1/2{{Zim2}} )
    ComputeNtimestimesdF12(Fwd, Bwd, m_ZimFwd[0], m_ZimBwd[0], m_ZimFwd[1],
                           m_ZimBwd[1], e3Fwd_cdot_dEe3, e3Bwd_cdot_dEe3);

    Vmath::Svtvp(nTraceNumPoints, -1.0, numfluxFwd[0], 1, e1Fwd_cdot_ncrossdH,
                 1, numfluxFwd[0], 1);
    Vmath::Svtvp(nTraceNumPoints, -1.0, numfluxFwd[1], 1, e2Fwd_cdot_ncrossdH,
                 1, numfluxFwd[1], 1);
    Vmath::Svtvp(nTraceNumPoints, 1.0, numfluxFwd[2], 1, e3Fwd_cdot_dEe3, 1,
                 numfluxFwd[2], 1);

    Vmath::Svtvp(nTraceNumPoints, -1.0, numfluxBwd[0], 1, e1Bwd_cdot_ncrossdH,
                 1, numfluxBwd[0], 1);
    Vmath::Svtvp(nTraceNumPoints, -1.0, numfluxBwd[1], 1, e2Bwd_cdot_ncrossdH,
                 1, numfluxBwd[1], 1);
    Vmath::Svtvp(nTraceNumPoints, 1.0, numfluxBwd[2], 1, e3Bwd_cdot_dEe3, 1,
                 numfluxBwd[2], 1);
}

void MMFMaxwell::NumericalMaxwellFluxTE(
    const Array<OneD, const Array<OneD, NekDouble>> &physfield,
    Array<OneD, Array<OneD, NekDouble>> &numfluxFwd,
    Array<OneD, Array<OneD, NekDouble>> &numfluxBwd, const NekDouble time)

{
    int nq              = GetNpoints();
    int nTraceNumPoints = GetTraceTotPoints();
    int nvar            = 3;

    // Get temporary arrays
    Array<OneD, Array<OneD, NekDouble>> Fwd(nvar);
    Array<OneD, Array<OneD, NekDouble>> Bwd(nvar);
    for (int i = 0; i < nvar; ++i)
    {
        Fwd[i] = Array<OneD, NekDouble>(nTraceNumPoints, 0.0);
        Bwd[i] = Array<OneD, NekDouble>(nTraceNumPoints, 0.0);

        // get the physical values at the trace
        m_fields[i]->GetFwdBwdTracePhys(physfield[i], Fwd[i], Bwd[i]);
    }

    // E = 0 at the PEC boundaries:
    Array<OneD, NekDouble> IncField(nq, 0.0);
    Array<OneD, NekDouble> IncFieldBwd(nTraceNumPoints, 0.0);
    Array<OneD, Array<OneD, NekDouble>> IncFieldFwd(m_spacedim);
    for (int i = 0; i < m_spacedim; ++i)
    {
        IncFieldFwd[i] = Array<OneD, NekDouble>(nTraceNumPoints, 0.0);

        IncField = GetIncidentField(i, time);
        m_fields[0]->GetFwdBwdTracePhys(IncField, IncFieldFwd[i], IncFieldBwd);

        Vmath::Svtvp(nTraceNumPoints, 2.0, &IncFieldFwd[i][0], 1, &Fwd[i][0], 1,
                     &IncFieldFwd[i][0], 1);
    }

    // CopyBoundaryTrace(IncFieldFwd[0], Bwd[0], SolverUtils::eFwdEQNegBwd, 0,
    //                   SpatialDomains::ePEC_Forces);
    // CopyBoundaryTrace(IncFieldFwd[1], Bwd[1], SolverUtils::eFwdEQNegBwd, 1,
    //                   SpatialDomains::ePEC_Forces);
    // CopyBoundaryTrace(Fwd[2], Bwd[2], SolverUtils::eFwdEQBwd, 2,
    //                   SpatialDomains::ePEC_Forces);

    // CopyBoundaryTrace(Fwd[0], Bwd[0], SolverUtils::eFwdEQBwd, 0,
    //                   SpatialDomains::ePMC_Forces);
    // CopyBoundaryTrace(Fwd[1], Bwd[1], SolverUtils::eFwdEQBwd, 1,
    //                   SpatialDomains::ePMC_Forces);
    // CopyBoundaryTrace(IncFieldFwd[2], Bwd[2], SolverUtils::eFwdEQNegBwd, 2,
    //                   SpatialDomains::ePMC_Forces);

    // // Boundary conditions for scattered waves
    // CopyBoundaryTrace(Fwd[0], Bwd[0], SolverUtils::eFwdEQNegBwd, 0,
    //                   SpatialDomains::ePEC);
    // CopyBoundaryTrace(Fwd[1], Bwd[1], SolverUtils::eFwdEQNegBwd, 1,
    //                   SpatialDomains::ePEC);
    // CopyBoundaryTrace(Fwd[2], Bwd[2], SolverUtils::eFwdEQBwd, 2,
    //                   SpatialDomains::ePEC);

    // CopyBoundaryTrace(Fwd[0], Bwd[0], SolverUtils::eFwdEQBwd, 0,
    //                   SpatialDomains::ePMC);
    // CopyBoundaryTrace(Fwd[1], Bwd[1], SolverUtils::eFwdEQBwd, 1,
    //                   SpatialDomains::ePMC);
    // CopyBoundaryTrace(Fwd[2], Bwd[2], SolverUtils::eFwdEQNegBwd, 2,
    //                   SpatialDomains::ePMC);

    Array<OneD, NekDouble> e1Fwd_cdot_ncrossdE(nTraceNumPoints);
    Array<OneD, NekDouble> e1Bwd_cdot_ncrossdE(nTraceNumPoints);
    Array<OneD, NekDouble> e2Fwd_cdot_ncrossdE(nTraceNumPoints);
    Array<OneD, NekDouble> e2Bwd_cdot_ncrossdE(nTraceNumPoints);
    Array<OneD, NekDouble> e3Fwd_cdot_dHe3(nTraceNumPoints);
    Array<OneD, NekDouble> e3Bwd_cdot_dHe3(nTraceNumPoints);

    // Compute  numfluxFwd[dir] = (eFwd^[dir] \cdot n \times e^3) * (ZimFwd *
    // HFwd^3 + ZimBwd * HBwd^3 )
    // Compute  numfluxBwd[dir] = (eBwd^[dir] \cdot n \times e^3) * (ZimFwd *
    // HFwd^3 + ZimBwd * HBwd^3 )
    ComputeNtimesFz(0, Fwd, Bwd, m_ZimFwd[0], m_ZimBwd[0], numfluxFwd[0],
                    numfluxBwd[0]);
    ComputeNtimesFz(1, Fwd, Bwd, m_ZimFwd[1], m_ZimBwd[1], numfluxFwd[1],
                    numfluxBwd[1]);

    // Compute numfluxFwd[2] = eFwd^3 \cdot ( n1e1 \times ( imFwd EFwd + imBwd
    // EBwd ) ) / 2 {{Y_i}}
    ComputeNtimesF12(Fwd, Bwd, m_YimFwd[0], m_YimBwd[0], m_YimFwd[1],
                     m_YimBwd[1], numfluxFwd[2], numfluxBwd[2]);

    // Compute e1Fwd_cdot_ncrossdE = eFwd[dir] \cdot \alpha n \times n \times
    // [E] / 2 {{ZimFwd}}
    ComputeNtimestimesdFz(0, Fwd, Bwd, m_ZimFwd[0], m_ZimBwd[0],
                          e1Fwd_cdot_ncrossdE, e1Bwd_cdot_ncrossdE);
    ComputeNtimestimesdFz(1, Fwd, Bwd, m_ZimFwd[1], m_ZimBwd[1],
                          e2Fwd_cdot_ncrossdE, e2Bwd_cdot_ncrossdE);

    // Compute  - \alpha [H3] * ( 1/2{{Yim1}} + 1/2{{Yim2}} )
    ComputeNtimestimesdF12(Fwd, Bwd, m_YimFwd[0], m_YimBwd[0], m_YimFwd[1],
                           m_YimBwd[1], e3Fwd_cdot_dHe3, e3Bwd_cdot_dHe3);

    Vmath::Svtvp(nTraceNumPoints, 1.0, numfluxFwd[0], 1, e1Fwd_cdot_ncrossdE, 1,
                 numfluxFwd[0], 1);
    Vmath::Svtvp(nTraceNumPoints, 1.0, numfluxFwd[1], 1, e2Fwd_cdot_ncrossdE, 1,
                 numfluxFwd[1], 1);
    Vmath::Svtvp(nTraceNumPoints, -1.0, numfluxFwd[2], 1, e3Fwd_cdot_dHe3, 1,
                 numfluxFwd[2], 1);

    Vmath::Svtvp(nTraceNumPoints, 1.0, numfluxBwd[0], 1, e1Bwd_cdot_ncrossdE, 1,
                 numfluxBwd[0], 1);
    Vmath::Svtvp(nTraceNumPoints, 1.0, numfluxBwd[1], 1, e2Bwd_cdot_ncrossdE, 1,
                 numfluxBwd[1], 1);
    Vmath::Svtvp(nTraceNumPoints, -1.0, numfluxBwd[2], 1, e3Bwd_cdot_dHe3, 1,
                 numfluxBwd[2], 1);
}



// Compute  (e^[dir] \cdot ( n \times e^3)) * (imFwd * Fwd^3 + imBwd * Bwd^3 )
void MMFMaxwell::ComputeNtimesFz(const int dir,
                                const Array<OneD, Array<OneD, NekDouble>> &Fwd,
                                const Array<OneD, Array<OneD, NekDouble>> &Bwd,
                                const Array<OneD, const NekDouble> &imFwd,
                                const Array<OneD, const NekDouble> &imBwd,
                                Array<OneD, NekDouble> &outarrayFwd,
                                Array<OneD, NekDouble> &outarrayBwd)

{
    int nTraceNumPoints = GetTraceTotPoints();

    NekDouble tmpFwd, tmpBwd;
    NekDouble Aver, ntimesz;

    for (int i = 0; i < nTraceNumPoints; ++i)
    {
        Aver = 0.5 * (imFwd[i] + imBwd[i]);

        tmpFwd = 0.0;
        tmpBwd = 0.0;
        for (int k = 0; k < m_spacedim; ++k)
        {
            ntimesz = 0.5 * (imFwd[i] * Fwd[2][i] + imBwd[i] * Bwd[2][i]);

            tmpFwd +=
                m_MFtraceFwd[dir][k][i] * m_ntimesMFFwd[2][k][i] * ntimesz;
            tmpBwd +=
                m_MFtraceBwd[dir][k][i] * m_ntimesMFBwd[2][k][i] * ntimesz;
        }

        outarrayFwd[i] = tmpFwd / Aver;
        outarrayBwd[i] = tmpBwd / Aver;
    }
}

// Compute e^3 \cdot ( n1e1 \times ( imFwd EFwd + imBwd EBwd ) ) / 2{{Y_i}}
void MMFMaxwell::ComputeNtimesF12(const Array<OneD, Array<OneD, NekDouble>> &Fwd,
                                 const Array<OneD, Array<OneD, NekDouble>> &Bwd,
                                 const Array<OneD, const NekDouble> &im1Fwd,
                                 const Array<OneD, const NekDouble> &im1Bwd,
                                 const Array<OneD, const NekDouble> &im2Fwd,
                                 const Array<OneD, const NekDouble> &im2Bwd,
                                 Array<OneD, NekDouble> &outarrayFwd,
                                 Array<OneD, NekDouble> &outarrayBwd)
{
    int nTraceNumPoints = GetTraceTotPoints();

    NekDouble tmpFwd, tmpBwd, Aver1, Aver2, HFwdk, HBwdk;

    Array<OneD, NekDouble> z1HAver(m_spacedim);
    Array<OneD, NekDouble> z2HAver(m_spacedim);
    Array<OneD, NekDouble> n1e1(m_spacedim);
    Array<OneD, NekDouble> n2e2(m_spacedim);

    Array<OneD, NekDouble> n1e1_times_z1HAver(m_spacedim);
    Array<OneD, NekDouble> n2e2_times_z2HAver(m_spacedim);

    for (int i = 0; i < nTraceNumPoints; ++i)
    {
        Aver1 = 0.5 * (im1Fwd[i] + im1Bwd[i]);
        Aver2 = 0.5 * (im2Fwd[i] + im2Bwd[i]);

        for (int k = 0; k < m_spacedim; k++)
        {
            // Compute \vec{HFwd} and \vec{HBwd}
            HFwdk = Fwd[0][i] * m_MFtraceFwd[0][k][i] +
                    Fwd[1][i] * m_MFtraceFwd[1][k][i];
            HBwdk = Bwd[0][i] * m_MFtraceBwd[0][k][i] +
                    Bwd[1][i] * m_MFtraceBwd[1][k][i];

            // Compute z_i {{ \vec{H} }}
            z1HAver[k] = 0.5 * (im1Fwd[i] * HFwdk + im1Bwd[i] * HBwdk);
            z2HAver[k] = 0.5 * (im2Fwd[i] * HFwdk + im2Bwd[i] * HBwdk);

            // Choose e^i for the one in anisotropy region
            n1e1[k] = m_ncdotMFFwd[0][i] * m_MFtraceFwd[0][k][i];
            n2e2[k] = m_ncdotMFFwd[1][i] * m_MFtraceFwd[1][k][i];
        }

        // Compute n1e1 \times z1HAver and n2e2 \times z2HAver
        VectorCrossProd(n1e1, z1HAver, n1e1_times_z1HAver);
        VectorCrossProd(n2e2, z2HAver, n2e2_times_z2HAver);

        // e^3 \cdot ( n1e1 \times z1HAver + n2e2 \times z2HAver)
        tmpFwd = 0.0;
        tmpBwd = 0.0;
        for (int k = 0; k < m_spacedim; k++)
        {
            tmpFwd += m_MFtraceFwd[2][k][i] * (n1e1_times_z1HAver[k] / Aver1 +
                                               n2e2_times_z2HAver[k] / Aver2);
            tmpBwd += m_MFtraceBwd[2][k][i] * (n1e1_times_z1HAver[k] / Aver1 +
                                               n2e2_times_z2HAver[k] / Aver2);
        }

        outarrayFwd[i] = tmpFwd;
        outarrayBwd[i] = tmpBwd;
    }
}

void MMFMaxwell::ComputeNtimestimesdFz(
    const int dir, const Array<OneD, Array<OneD, NekDouble>> &Fwd,
    const Array<OneD, Array<OneD, NekDouble>> &Bwd,
    const Array<OneD, const NekDouble> &imFwd,
    const Array<OneD, const NekDouble> &imBwd,
    Array<OneD, NekDouble> &outarrayFwd, Array<OneD, NekDouble> &outarrayBwd)
{
    int nTraceNumPoints = GetTraceTotPoints();

    Array<OneD, NekDouble> dH(m_spacedim);
    Array<OneD, NekDouble> nFwd(m_spacedim);

    Array<OneD, NekDouble> eiFwd(m_spacedim);
    Array<OneD, NekDouble> eiBwd(m_spacedim);

    Array<OneD, NekDouble> eitimesdHFwd(m_spacedim);
    Array<OneD, NekDouble> eitimesdHBwd(m_spacedim);

    Array<OneD, NekDouble> ntimeseitimesdHFwd(m_spacedim);
    Array<OneD, NekDouble> ntimeseitimesdHBwd(m_spacedim);

    NekDouble Aver, HFwdk, HBwdk, tmpFwd, tmpBwd;
    for (int i = 0; i < nTraceNumPoints; ++i)
    {
        Aver = 0.5 * (imFwd[i] + imBwd[i]);

        // Get [H]
        for (int k = 0; k < m_spacedim; k++)
        {
            HFwdk = Fwd[0][i] * m_MFtraceFwd[0][k][i] +
                    Fwd[1][i] * m_MFtraceFwd[1][k][i];
            HBwdk = Bwd[0][i] * m_MFtraceBwd[0][k][i] +
                    Bwd[1][i] * m_MFtraceBwd[1][k][i];
            dH[k] = HFwdk - HBwdk;

            eiFwd[k] = m_MFtraceFwd[dir][k][i];
            eiBwd[k] = m_MFtraceBwd[dir][k][i];

            nFwd[k] = m_traceNormals[k][i];
        }

        // MFtraceFwd (MFtraceBwd) \times [H]
        // VectorCrossProd(eiFwd, dH, eitimesdHFwd);
        // VectorCrossProd(eiBwd, dH, eitimesdHBwd);
        VectorCrossProd(nFwd, dH, eitimesdHFwd);
        VectorCrossProd(nFwd, dH, eitimesdHBwd);

        // n times eitimesdH
        VectorCrossProd(nFwd, eitimesdHFwd, ntimeseitimesdHFwd);
        VectorCrossProd(nFwd, eitimesdHBwd, ntimeseitimesdHBwd);

        // MFtraceFwd \cdot ntimeseitimesdH
        tmpFwd = 0.0;
        tmpBwd = 0.0;
        for (int k = 0; k < m_spacedim; k++)
        {
            tmpFwd += eiFwd[k] * ntimeseitimesdHFwd[k];
            tmpBwd += eiBwd[k] * ntimeseitimesdHBwd[k];
        }

        outarrayFwd[i] = 0.5 * m_alpha * tmpFwd / Aver;
        outarrayBwd[i] = 0.5 * m_alpha * tmpBwd / Aver;
    }
}

// Compute - \alpha [E3] / ( {{im1}} + {{im2}} )
void MMFMaxwell::ComputeNtimestimesdF12(
    const Array<OneD, Array<OneD, NekDouble>> &Fwd,
    const Array<OneD, Array<OneD, NekDouble>> &Bwd,
    const Array<OneD, const NekDouble> &im1Fwd,
    const Array<OneD, const NekDouble> &im1Bwd,
    const Array<OneD, const NekDouble> &im2Fwd,
    const Array<OneD, const NekDouble> &im2Bwd,
    Array<OneD, NekDouble> &outarrayFwd, Array<OneD, NekDouble> &outarrayBwd)
{
    int nTraceNumPoints = GetTraceTotPoints();

    Array<OneD, NekDouble> directFwd(nTraceNumPoints);
    Array<OneD, NekDouble> directBwd(nTraceNumPoints);

    NekDouble Aver1, Aver2;
    for (int i = 0; i < nTraceNumPoints; ++i)
    {
        Aver1 = im1Fwd[i] + im1Bwd[i];
        Aver2 = im2Fwd[i] + im2Bwd[i];

        outarrayFwd[i] =
            -m_alpha * (Fwd[2][i] - Bwd[2][i]) * (1.0 / Aver1 + 1.0 / Aver2);
        outarrayBwd[i] =
            -m_alpha * (Fwd[2][i] - Bwd[2][i]) * (1.0 / Aver1 + 1.0 / Aver2);
    }
}

/**
 * @brief Compute the projection for the linear advection equation.
 *
 * @param inarray    Given fields.
 * @param outarray   Calculated solution.
 * @param time       Time.
 */
void MMFMaxwell::DoOdeProjection(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time)
{
    boost::ignore_unused(time);

    int var = inarray.size();

    // SetBoundaryConditions(time);

    int nq = GetNpoints();
    for (int i = 0; i < var; ++i)
    {
        Vmath::Vcopy(nq, inarray[i], 1, outarray[i], 1);
    }
}

void MMFMaxwell::v_SetInitialConditions(const NekDouble initialtime,
                                        bool dumpInitialConditions,
                                        const int domain)
{
    boost::ignore_unused(domain);

    int nq   = GetTotPoints();
    int nvar = m_fields.size();

    switch (m_TestMaxwellType)
    {
        case eMaxwell1D:
        {
            m_fields[0]->SetPhys(TestMaxwell1D(initialtime, 0));
            m_fields[1]->SetPhys(TestMaxwell1D(initialtime, 1));
        }
        break;

        case eTestMaxwell2DPEC:
        case eTestMaxwell2DPECAVGFLUX:
        {
            m_fields[0]->SetPhys(TestMaxwell2DPEC(initialtime, 0, m_PolType));
            m_fields[1]->SetPhys(TestMaxwell2DPEC(initialtime, 1, m_PolType));
            m_fields[2]->SetPhys(TestMaxwell2DPEC(initialtime, 2, m_PolType));
        }
        break;

        case eTestMaxwell2DPMC:
        {
            m_fields[0]->SetPhys(TestMaxwell2DPMC(initialtime, 0, m_PolType));
            m_fields[1]->SetPhys(TestMaxwell2DPMC(initialtime, 1, m_PolType));
            m_fields[2]->SetPhys(TestMaxwell2DPMC(initialtime, 2, m_PolType));
        }
        break;

        case eScatField2D:
        case eTotField2D:
        {
            Array<OneD, NekDouble> Zeros(nq, 0.0);

            for (int i = 0; i < nvar; i++)
            {
                m_fields[i]->SetPhys(Zeros);
            }
        }
        break;

        case eMaxwellSphere:
        {
            m_fields[0]->SetPhys(TestMaxwellSphere(initialtime, m_freq, 0));
            m_fields[1]->SetPhys(TestMaxwellSphere(initialtime, m_freq, 1));
            m_fields[2]->SetPhys(TestMaxwellSphere(initialtime, m_freq, 2));
        }
        break;

        case eELF2DSurface:
        {
            m_fields[2]->SetPhys(GaussianPulse(initialtime, m_Psx, m_Psy, m_Psz,
                                               m_Gaussianradius));
        }
        break;

        default:
            break;
    }

    // forward transform to fill the modal coeffs
    for (int i = 0; i < nvar; ++i)
    {
        m_fields[i]->SetPhysState(true);
        m_fields[i]->FwdTrans(m_fields[i]->GetPhys(),
                              m_fields[i]->UpdateCoeffs());
    }

    if (dumpInitialConditions)
    {
        std::string outname = m_sessionName + "_initial.chk";
        WriteFld(outname);

        Array<OneD, Array<OneD, NekDouble>> fields(nvar);
        for (int i = 0; i < nvar; ++i)
        {
            fields[i] = m_fields[i]->GetPhys();
        }

        Checkpoint_PlotOutput(0, fields);
    }
}

void MMFMaxwell::v_EvaluateExactSolution(unsigned int field,
                                         Array<OneD, NekDouble> &outfield,
                                         const NekDouble time)
{
    int nq              = GetNpoints();
    outfield = Array<OneD, NekDouble>(nq);

    switch (m_TestMaxwellType)
    {
        case eMaxwell1D:
        {
            outfield = TestMaxwell1D(time, field);
        }
        break;

        case eTestMaxwell2DPEC:
        case eTestMaxwell2DPECAVGFLUX:
        {
            outfield = TestMaxwell2DPEC(time, field, m_PolType);
        }
        break;

        case eTestMaxwell2DPMC:
        {
            outfield = TestMaxwell2DPMC(time, field, m_PolType);
        }
        break;

        case eMaxwellSphere:
        {
            outfield = TestMaxwellSphere(time, m_freq, field);
        }
        break;

        default:
        {
            outfield = Array<OneD, NekDouble>(nq, 0.0);
        }
        break;
    }
}

Array<OneD, NekDouble> MMFMaxwell::TestMaxwell1D(const NekDouble time,
                                                 unsigned int field)
{
    int nq              = GetNpoints();

    Array<OneD, NekDouble> x0(nq);
    Array<OneD, NekDouble> x1(nq);
    Array<OneD, NekDouble> x2(nq);

    m_fields[0]->GetCoords(x0, x1, x2);

    Array<OneD, NekDouble> E(nq);
    Array<OneD, NekDouble> H(nq);

    // Derive the frequency \omega
    NekDouble omega;
    NekDouble Tol = 0.000000001;
    if (fabs(m_n1 - m_n2) < Tol)
    {
        omega = m_pi / m_n1;
    }

    else
    {
        omega = 2.0 * m_pi / m_n2;

        NekDouble newomega, F, Fprime;
        for (int i = 0; i < 10000; ++i)
        {
            F = m_n1 * tan(m_n2 * omega) + m_n2 * tan(m_n1 * omega);
            Fprime =
                m_n1 * m_n2 * (1.0 / cos(m_n2 * omega) / cos(m_n2 * omega) +
                               1.0 / cos(m_n1 * omega) / cos(m_n1 * omega));

            newomega = omega - F / Fprime;

            if (fabs(newomega - omega) > Tol)
            {
                omega = newomega;
            }

            else
            {
                break;
            }
        }
    }

    // Generate A^k and B^k
    std::complex<double> im = sqrt(std::complex<double>(-1));
    std::complex<double> A1, A2, B1, B2;
    std::complex<double> Ak, Bk, nk;
    std::complex<double> Ec, Hc;

    A1 = m_n2 * cos(m_n2 * omega) / (m_n1 * cos(m_n1 * omega));
    A2 = exp(-1.0 * im * omega * (m_n1 + m_n2));
    B1 = A1 * exp(-2.0 * im * m_n1 * omega);
    B2 = A2 * exp(2.0 * im * m_n2 * omega);

    for (int i = 0; i < nq; ++i)
    {
        if (x0[i] > 0)
        {
            Ak = A2;
            Bk = B2;
            nk = m_n2;
        }

        else
        {
            Ak = A1;
            Bk = B1;
            nk = m_n1;
        }

        Ec = (Ak * exp(im * nk * omega * x0[i]) -
              Bk * exp(-im * nk * omega * x0[i])) *
             exp(im * omega * time);
        Hc = nk * (Ak * exp(im * nk * omega * x0[i]) +
                   Bk * exp(-im * nk * omega * x0[i])) *
             exp(im * omega * time);

        E[i] = Ec.real();
        H[i] = Hc.real();
    }

    Array<OneD, NekDouble> outfield;
    switch (field)
    {
        case (0):
        {
            outfield = E;
        }
        break;

        case (1):
        {
            outfield = H;
        }
        break;
    }

    return outfield;
}

Array<OneD, NekDouble> MMFMaxwell::TestMaxwell2DPEC(
    const NekDouble time, unsigned int field,
    const PolType Polarization)
{
    int nq              = GetNpoints();

    Array<OneD, NekDouble> x0(nq);
    Array<OneD, NekDouble> x1(nq);
    Array<OneD, NekDouble> x2(nq);

    m_fields[0]->GetCoords(x0, x1, x2);

    NekDouble freqm = 1.0, freqn = 1.0;
    NekDouble omega = m_pi * sqrt(freqm * freqm + freqn * freqn);
    NekDouble mpi   = freqm * m_pi;
    NekDouble npi   = freqn * m_pi;

    Array<OneD, NekDouble> F1(nq);
    Array<OneD, NekDouble> F2(nq);
    Array<OneD, NekDouble> Fz(nq);
    Array<OneD, NekDouble> dF1dt(nq);
    Array<OneD, NekDouble> dF2dt(nq);
    Array<OneD, NekDouble> dFzdt(nq);
    NekDouble Fx, Fy, dFxdt, dFydt;

    for (int i = 0; i < nq; ++i)
    {
        switch (Polarization)
        {
            case eTransMagnetic:
            {
                Fx = -1.0 * (npi / omega) * sin(mpi * x0[i]) *
                     cos(npi * x1[i]) * sin(omega * time);
                Fy = (mpi / omega) * cos(mpi * x0[i]) * sin(npi * x1[i]) *
                     sin(omega * time);

                F1[i] =
                    Fx * m_movingframes[0][i] + Fy * m_movingframes[0][i + nq];
                F2[i] =
                    Fx * m_movingframes[1][i] + Fy * m_movingframes[1][i + nq];
                Fz[i] = sin(mpi * x0[i]) * sin(npi * x1[i]) * cos(omega * time);

                dFxdt = (npi)*sin(mpi * x0[i]) * cos(npi * x1[i]) *
                        cos(omega * time);
                dFydt = -1.0 * (mpi)*cos(mpi * x0[i]) * sin(npi * x1[i]) *
                        cos(omega * time);

                dF1dt[i] = dFxdt * m_movingframes[0][i] +
                           dFydt * m_movingframes[0][i + nq];
                dF2dt[i] = dFxdt * m_movingframes[1][i] +
                           dFydt * m_movingframes[1][i + nq];
                dFzdt[i] = omega * sin(mpi * x0[i]) * sin(npi * x1[i]) *
                           sin(omega * time);
            }
            break;

            case eTransElectric:
            {
                Fx = -1.0 * (npi / omega) * cos(mpi * x0[i]) *
                     sin(npi * x1[i]) * sin(omega * time);
                Fy = (mpi / omega) * sin(mpi * x0[i]) * cos(npi * x1[i]) *
                     sin(omega * time);

                F1[i] =
                    Fx * m_movingframes[0][i] + Fy * m_movingframes[0][i + nq];
                F2[i] =
                    Fx * m_movingframes[1][i] + Fy * m_movingframes[1][i + nq];
                Fz[i] = cos(mpi * x0[i]) * cos(npi * x1[i]) * cos(omega * time);

                dFxdt = (npi)*cos(mpi * x0[i]) * sin(npi * x1[i]) *
                        cos(omega * time);
                dFydt = -1.0 * (mpi)*sin(mpi * x0[i]) * cos(npi * x1[i]) *
                        cos(omega * time);

                dF1dt[i] = dFxdt * m_movingframes[0][i] +
                           dFydt * m_movingframes[0][i + nq];
                dF2dt[i] = dFxdt * m_movingframes[1][i] +
                           dFydt * m_movingframes[1][i + nq];
                dFzdt[i] = omega * cos(mpi * x0[i]) * cos(npi * x1[i]) *
                           sin(omega * time);
            }
            break;

            default:
                break;
        }
    }

    Array<OneD, NekDouble> outfield;
    switch (field)
    {
        case (0):
        {
            outfield = F1;
        }
        break;

        case (1):
        {
            outfield = F2;
        }
        break;

        case (2):
        {
            outfield = Fz;
        }
        break;

        case (10):
        {
            outfield = dF1dt;
        }
        break;

        case (11):
        {
            outfield = dF2dt;
        }
        break;

        case (12):
        {
            outfield = dFzdt;
        }
        break;
    }

    return outfield;
}

Array<OneD, NekDouble> MMFMaxwell::TestMaxwell2DPMC(
    const NekDouble time, unsigned int field,
    const PolType Polarization)
{
    int nq              = GetNpoints();

    Array<OneD, NekDouble> x0(nq);
    Array<OneD, NekDouble> x1(nq);
    Array<OneD, NekDouble> x2(nq);

    m_fields[0]->GetCoords(x0, x1, x2);

    NekDouble freqm = 1.0, freqn = 1.0;
    NekDouble omega = m_pi * sqrt(freqm * freqm + freqn * freqn);
    NekDouble mpi   = freqm * m_pi;
    NekDouble npi   = freqn * m_pi;

    Array<OneD, NekDouble> F1(nq);
    Array<OneD, NekDouble> F2(nq);
    Array<OneD, NekDouble> Fz(nq);
    NekDouble Fx, Fy;


    for (int i = 0; i < nq; ++i)
    {
        switch (Polarization)
        {
            case eTransMagnetic:
            {
                Fx = (npi / omega) * cos(mpi * x0[i]) * sin(npi * x1[i]) *
                     sin(omega * time);
                Fy = -(mpi / omega) * sin(mpi * x0[i]) * cos(npi * x1[i]) *
                     sin(omega * time);

                F1[i] =
                    Fx * m_movingframes[0][i] + Fy * m_movingframes[0][i + nq];
                F2[i] =
                    Fx * m_movingframes[1][i] + Fy * m_movingframes[1][i + nq];
                Fz[i] = cos(mpi * x0[i]) * cos(npi * x1[i]) * cos(omega * time);
            }
            break;

            case eTransElectric:
            {
                Fx = (npi / omega) * sin(mpi * x0[i]) * cos(npi * x1[i]) *
                     sin(omega * time);
                Fy = -(mpi / omega) * cos(mpi * x0[i]) * sin(npi * x1[i]) *
                     sin(omega * time);

                F1[i] =
                    Fx * m_movingframes[0][i] + Fy * m_movingframes[0][i + nq];
                F2[i] =
                    Fx * m_movingframes[1][i] + Fy * m_movingframes[1][i + nq];
                Fz[i] = sin(mpi * x0[i]) * sin(npi * x1[i]) * cos(omega * time);
            }
            break;

            default:
                break;
        }
    }

    Array<OneD, NekDouble> outfield;
    switch (field)
    {
        case (0):
        {
            outfield = F1;
        }
        break;

        case (1):
        {
            outfield = F2;
        }
        break;

        case (2):
        {
            outfield = Fz;
        }
        break;
    }

    return outfield;
}

Array<OneD, NekDouble> MMFMaxwell::TestMaxwellSphere(const NekDouble time,
                                                     const NekDouble omega,
                                                     unsigned int field)
{
    int nq = GetNpoints();

    Array<OneD, NekDouble> outfield(nq);

    Array<OneD, NekDouble> x(nq);
    Array<OneD, NekDouble> y(nq);
    Array<OneD, NekDouble> z(nq);

    m_fields[0]->GetCoords(x, y, z);

    Array<OneD, NekDouble> H1(nq);
    Array<OneD, NekDouble> H2(nq);
    Array<OneD, NekDouble> E3(nq);
    Array<OneD, NekDouble> F1(nq);
    Array<OneD, NekDouble> F2(nq);

    Array<OneD, NekDouble> curlv(nq);
    Array<OneD, Array<OneD, NekDouble>> velvec(m_spacedim);
    Array<OneD, Array<OneD, NekDouble>> Fvec(m_spacedim);
    for (int i = 0; i < m_spacedim; ++i)
    {
        velvec[i] = Array<OneD, NekDouble>(nq, 0.0);
        Fvec[i]   = Array<OneD, NekDouble>(nq, 0.0);
    }

    NekDouble xj, yj, zj, sin_varphi, cos_varphi, sin_theta, cos_theta;
    NekDouble vth, vphi, Fth, Fphi;
    for (int i = 0; i < nq; i++)
    {
        xj = x[i];
        yj = y[i];
        zj = z[i];

        // CartesianToSpherical(xj, yj, zj, sin_varphi, cos_varphi, sin_theta,
        //                      cos_theta);

        // vth = -4.0 * sin_varphi * cos_varphi * cos_theta * cos_theta *
        //       cos_theta * sin_theta;
        // vphi = -1.0 * sin_varphi * sin_varphi * cos_theta * cos_theta * cos_theta;

        // velvec[0][i] = -vth * sin_theta * cos_varphi - vphi * sin_varphi;
        // velvec[1][i] = -vth * sin_theta * sin_varphi + vphi * cos_varphi;
        // velvec[2][i] = vth * cos_theta;

        // E3[i] = (-4.0 * cos_theta * cos_theta * sin_theta * cos_varphi *
        //          cos_varphi) *
        //         (1.0 / omega * sin(omega * time));

        // Fth = -omega * vth -
        //       (8.0 / omega) * cos_theta * sin_theta * cos_varphi * sin_varphi;
        // Fphi = -omega * vphi +
        //        (4.0 / omega) * cos_varphi * cos_varphi * cos_theta *
        //            (2.0 * sin_theta * sin_theta - cos_theta * cos_theta);
        // Fvec[0][i] = -Fth * sin_theta * cos_varphi - Fphi * sin_varphi;
        // Fvec[1][i] = -Fth * sin_theta * sin_varphi + Fphi * cos_varphi;
        // Fvec[2][i] = Fth * cos_theta;

        CartesianToNewSpherical(xj, yj, zj, sin_varphi, cos_varphi, sin_theta,
                             cos_theta);

        vth = -4.0 * sin_varphi * cos_varphi * sin_theta * sin_theta *
              sin_theta * cos_theta;
        vphi = -1.0 * sin_varphi * sin_varphi * sin_theta * sin_theta * sin_theta;

        velvec[0][i] = -vth * cos_theta * cos_varphi - vphi * sin_varphi;
        velvec[1][i] = -vth * cos_theta * sin_varphi + vphi * cos_varphi;
        velvec[2][i] = vth * sin_theta;

        E3[i] = (-4.0 * sin_theta * sin_theta * cos_theta * cos_varphi * cos_varphi) *
                (1.0 / omega * sin(omega * time));

        Fth = -omega * vth - (8.0 / omega) * sin_theta * cos_theta * cos_varphi * sin_varphi;

        Fphi = -omega * vphi +
               (4.0 / omega) * cos_varphi * cos_varphi * sin_theta *
                   (2.0 * cos_theta * cos_theta - sin_theta * sin_theta);

        Fvec[0][i] = -Fth * cos_theta * cos_varphi - Fphi * sin_varphi;
        Fvec[1][i] = -Fth * cos_theta * sin_varphi + Fphi * cos_varphi;
        Fvec[2][i] = Fth * sin_theta;
    }



    H1 = CartesianToMovingframes(m_movingframes,velvec, 0);
    H2 = CartesianToMovingframes(m_movingframes,velvec, 1);

    Vmath::Smul(nq, cos(omega * time), H1, 1, H1, 1);
    Vmath::Smul(nq, cos(omega * time), H2, 1, H2, 1);

    F1 = CartesianToMovingframes(m_movingframes,Fvec, 0);
    F2 = CartesianToMovingframes(m_movingframes,Fvec, 1);

    Vmath::Smul(nq, sin(omega * time), F1, 1, F1, 1);
    Vmath::Smul(nq, sin(omega * time), F2, 1, F2, 1);

    switch (field)
    {
        // return H1
        case 0:
        {
            outfield = H1;
        }
        break;

        case 1:
        {
            outfield = H2;
        }
        break;

        case 2:
        {
            outfield = E3;
        }
        break;

        case 3:
        {
            outfield = F1;
        }
        break;

        case 4:
        {
            outfield = F2;
        }
        break;

        default:
            break;
    }

    return outfield;
}

void MMFMaxwell::Printout_SurfaceCurrent(
    Array<OneD, Array<OneD, NekDouble>> &fields, const int time)
{
    int nq              = GetNpoints();
    int nTraceNumPoints = GetTraceTotPoints();

    int totbdryexp =
        m_fields[0]->GetBndCondExpansions()[m_boundaryforSF]->GetExpSize();
    int npts = m_fields[0]
                   ->GetBndCondExpansions()[m_boundaryforSF]
                   ->GetExp(0)
                   ->GetNumPoints(0);
    int totnpts = totbdryexp * npts;

    Array<OneD, NekDouble> Jphi(totnpts);
    Array<OneD, NekDouble> Jrad(totnpts);
    Array<OneD, NekDouble> Jcurrent(totnpts);

    Array<OneD, NekDouble> phiFwd(nTraceNumPoints);
    Array<OneD, NekDouble> radFwd(nTraceNumPoints);

    // Compute phiFwd = acos(-x/r) along the trace
    Array<OneD, NekDouble> x(nq);
    Array<OneD, NekDouble> y(nq);
    Array<OneD, NekDouble> z(nq);

    m_fields[0]->GetCoords(x, y, z);

    Array<OneD, NekDouble> xFwd(nTraceNumPoints);
    Array<OneD, NekDouble> yFwd(nTraceNumPoints);

    m_fields[0]->ExtractTracePhys(x, xFwd);
    m_fields[0]->ExtractTracePhys(y, yFwd);

    for (int i = 0; i < nTraceNumPoints; ++i)
    {
        radFwd[i] = sqrt(xFwd[i] * xFwd[i] / m_MMFfactors[0] / m_MMFfactors[0] +
                         yFwd[i] * yFwd[i] / m_MMFfactors[1] / m_MMFfactors[1]);
        phiFwd[i] = atan2(yFwd[i] / radFwd[i], -1.0 * xFwd[i] / radFwd[i]);
    }

    Array<OneD, NekDouble> ntimesHFwd(nTraceNumPoints);
    ntimesHFwd = ComputeSurfaceCurrent(time, fields);

    // The surface for current should be the first boundary
    // int id2, cnt = 0;
    // for (int e = 0; e < totbdryexp; ++e)
    // {
    //     id2 = m_fields[0]->GetTrace()->GetPhys_Offset(
    //         m_fields[0]->GetTraceMap()->GetBndCondIDToGlobalTraceID(cnt +                      e));

    //     Vmath::Vcopy(npts, &phiFwd[id2], 1, &Jphi[e * npts], 1);
    //     Vmath::Vcopy(npts, &radFwd[id2], 1, &Jrad[e * npts], 1);
    //     Vmath::Vcopy(npts, &ntimesHFwd[id2], 1, &Jcurrent[e * npts], 1);
    // }

    // Vmath::Vmul(totnpts, tmpr, 1, tmpr, 1, Jcurrent, 1);
    // Vmath::Vvtvp(totnpts, tmpi, 1, tmpi, 1, Jcurrent, 1, Jcurrent, 1);
    // Vmath::Vsqrt(totnpts, Jcurrent, 1, Jcurrent, 1);

    // std::cout << "========================================================"
    //          << std::endl;

    // std::cout << "phi = " << std::endl;
    // for (int i = 0; i < totnpts; ++i)
    // {
    //     std::cout << Jphi[i] << ", ";
    // }
    // std::cout << std::endl << std::endl;

    // std::cout << "J = " << std::endl;
    // for (int i = 0; i < totnpts; ++i)
    // {
    //     std::cout << Jcurrent[i] << ", ";
    // }
    // std::cout << std::endl << std::endl;
}

Array<OneD, NekDouble> MMFMaxwell::ComputeSurfaceCurrent(
    const int time, const Array<OneD, const Array<OneD, NekDouble>> &fields)
{
    int nq              = GetNpoints();
    int nTraceNumPoints = GetTraceTotPoints();

    Array<OneD, NekDouble> outfield(nTraceNumPoints, 0.0);

    switch (m_PolType)
    {
        // (n \times H^r)_z = H1r (n \times e^1)_z + H2r (n \times e^2)_z,
        case eTransMagnetic:
        {
            Array<OneD, NekDouble> tmp(nq);
            Array<OneD, NekDouble> tmpFwd(nTraceNumPoints);

            for (int i = 0; i < 2; i++)
            {
                tmp = GetIncidentField(i, time);
                Vmath::Vadd(nq, fields[i], 1, tmp, 1, tmp, 1);

                m_fields[0]->ExtractTracePhys(tmp, tmpFwd);
                Vmath::Vvtvp(nTraceNumPoints, &m_ntimesMFFwd[i][2][0], 1,
                             &tmpFwd[0], 1, &outfield[0], 1, &outfield[0], 1);
            }
        }
        break;

        case eTransElectric:
        {
            Array<OneD, NekDouble> tmp(nq);

            tmp = GetIncidentField(2, time);
            Vmath::Vadd(nq, fields[2], 1, tmp, 1, tmp, 1);
            m_fields[0]->ExtractTracePhys(tmp, outfield);
        }
        break;

        default:
            break;
    }

    return outfield;
}

void MMFMaxwell::GenerateSigmaPML(const NekDouble PMLthickness,
                                  const NekDouble PMLstart,
                                  const NekDouble PMLmaxsigma,
                                  Array<OneD, Array<OneD, NekDouble>> &SigmaPML)
{
    int nq              = GetNpoints();

    // Construct sigmaX and sigmaY for UPML
    Array<OneD, NekDouble> x(nq);
    Array<OneD, NekDouble> y(nq);
    Array<OneD, NekDouble> z(nq);

    m_fields[0]->GetCoords(x, y, z);

    // Construction of SigmaPML
    SigmaPML = Array<OneD, Array<OneD, NekDouble>>(m_shapedim);
    for (int j = 0; j < m_shapedim; j++)
    {
        SigmaPML[j] = Array<OneD, NekDouble>(nq, 0.0);
    }

    // PML region indicator: [ 0 : curvedPML : RecPML]
    // RecPML = [0 : 0 : 1]
    // PMLRegion = [0 : 1 : 1], CurvedPML = PMLRegion - RecPML = [0: 1: 0]
    Array<OneD, NekDouble> RecPML(nq);
    Array<OneD, NekDouble> CurvedPML(nq);
    Array<OneD, NekDouble> PMLRegion(nq);
    m_fields[0]->GenerateElementVector(m_RecPML, 0.0, 1.0, RecPML);
    m_fields[0]->GenerateElementVector(m_PMLelement, 0.0, 1.0, PMLRegion);
    Vmath::Vsub(nq, PMLRegion, 1, RecPML, 1, CurvedPML, 1);

    switch (m_AddPML)
    {
        // RecPML only
        case 1:
        {
            // Rectangular PML
            NekDouble xRlayer, xLlayer, yRlayer, yLlayer;

            xRlayer = Vmath::Vmax(nq, x, 1) - PMLthickness;
            xLlayer = Vmath::Vmin(nq, x, 1) + PMLthickness;
            yRlayer = Vmath::Vmax(nq, y, 1) - PMLthickness;
            yLlayer = Vmath::Vmin(nq, y, 1) + PMLthickness;

            NekDouble xd, yd;
            for (int i = 0; i < nq; i++)
            {
                // SimgaPML along e^1
                if (x[i] >= xRlayer)
                {
                    xd = (x[i] - xRlayer) / PMLthickness;
                }

                else if (x[i] <= xLlayer)
                {
                    xd = (xLlayer - x[i]) / PMLthickness;
                }

                else
                {
                    xd = 0.0;
                }

                SigmaPML[0][i] = RecPML[i] * PMLmaxsigma * (xd * xd * xd);

                // SigmaPML along e^2
                if (y[i] >= yRlayer)
                {
                    yd = (y[i] - yRlayer) / PMLthickness;
                }

                else if (y[i] <= yLlayer)
                {
                    yd = (yLlayer - y[i]) / PMLthickness;
                }

                else
                {
                    yd = 0.0;
                }

                SigmaPML[1][i] = PMLRegion[i] * PMLmaxsigma * (yd * yd * yd);
            }
        }
        break;

        // CurvedPML only
        case 2:
        {
            // Curved PML
            NekDouble relrad, rad;
            for (int i = 0; i < nq; i++)
            {
                rad = sqrt(x[i] * x[i] / m_MMFfactors[0] / m_MMFfactors[0] +
                           y[i] * y[i] / m_MMFfactors[1] / m_MMFfactors[1]);

                if (rad >= PMLstart)
                {
                    relrad = (rad - PMLstart) / PMLthickness;
                    SigmaPML[1][i] =
                        PMLRegion[i] * PMLmaxsigma * pow(relrad, m_PMLorder);
                }
            }
        }
        break;

        // Slanted PML
        case 3:
        {
            NekDouble relrad, radon, radtw, radth, radfo;
            for (int i = 0; i < nq; i++)
            {
                radon = -1.0 * x[i] + y[i] - 7;
                radtw = x[i] + y[i] - 7;
                radth = -x[i] - y[i] - 7;
                radfo = x[i] - y[i] - 7;

                if (radon >= 0.0)
                {
                    relrad = radon / PMLthickness;
                    SigmaPML[1][i] =
                        PMLRegion[i] * PMLmaxsigma * pow(relrad, m_PMLorder);
                }

                if (radtw >= 0.0)
                {
                    relrad = radtw / PMLthickness;
                    SigmaPML[0][i] =
                        PMLRegion[i] * PMLmaxsigma * pow(relrad, m_PMLorder);
                }

                if (radth >= 0.0)
                {
                    relrad = radth / PMLthickness;
                    SigmaPML[0][i] =
                        PMLRegion[i] * PMLmaxsigma * pow(relrad, m_PMLorder);
                }

                if (radfo >= 0.0)
                {
                    relrad = radfo / PMLthickness;
                    SigmaPML[1][i] =
                        PMLRegion[i] * PMLmaxsigma * pow(relrad, m_PMLorder);
                }
            }
        }
        break;
    }

    std::cout << "*** sigma1 = [ " << Vmath::Vmin(nq, &SigmaPML[0][0], 1)
              << " , " << Vmath::Vmax(nq, &SigmaPML[0][0], 1)
              << " ] , sigma2 = [ " << Vmath::Vmin(nq, &SigmaPML[1][0], 1)
              << " , " << Vmath::Vmax(nq, &SigmaPML[1][0], 1) << " ] "
              << std::endl;
}

NekDouble MMFMaxwell::ComputeEnergyDensity(
    Array<OneD, Array<OneD, NekDouble>> &fields)
{
    int nvar            = 3;
    int nq = GetTotPoints();

    NekDouble energy;

    Array<OneD, NekDouble> tmp(nq, 0.0);
    for (int i = 0; i < nvar; ++i)
    {
        Vmath::Vvtvp(nq, &fields[i][0], 1, &fields[i][0], 1, &tmp[0], 1,
                     &tmp[0], 1);
    }

    energy = 0.5 * (m_fields[0]->PhysIntegral(tmp));
    return energy;
}

void MMFMaxwell::ComputeMaterialVector(
    Array<OneD, Array<OneD, NekDouble>> &epsvec,
    Array<OneD, Array<OneD, NekDouble>> &muvec)
{
    switch (m_TestMaxwellType)
    {
        case eMaxwell1D:
        {
            m_fields[0]->GenerateElementVector(m_ElemtGroup1, m_varepsilon[0],
                                               m_varepsilon[1], epsvec[0]);
            m_fields[0]->GenerateElementVector(m_ElemtGroup1, 1.0, 1.0,
                                               muvec[0]);
        }
        break;

        case eTestMaxwell2DPEC:
        case eTestMaxwell2DPECAVGFLUX:
        case eTestMaxwell2DPMC:
        case eScatField2D:
        case eTotField2D:
        {
            switch (m_PolType)
            {
                case eTransMagnetic:
                {
                    m_fields[0]->GenerateElementVector(m_ElemtGroup1, m_mu[0],
                                                       1.0, muvec[0]);
                    m_fields[0]->GenerateElementVector(m_ElemtGroup1, m_mu[1],
                                                       1.0, muvec[1]);
                    m_fields[0]->GenerateElementVector(
                        m_ElemtGroup1, m_varepsilon[2], 1.0, epsvec[2]);

                    // // ONLY FOR VARIABLE ANISOTROPY TEST
                    // int nq  = GetTotPoints();

                    // Array<OneD, NekDouble> tmpIN(nq);

                    // m_fields[0]->GenerateElementVector(m_ElemtGroup1, 1.0,
                    // 0.0, tmpIN);

                    // Array<OneD, NekDouble> x0(nq);
                    // Array<OneD, NekDouble> x1(nq);
                    // Array<OneD, NekDouble> x2(nq);

                    // m_fields[0]->GetCoords(x0,x1,x2);

                    // for (int i=0; i<nq; i++)
                    //   {
                    //     muvec[1][i] = tmpIN[i]*(1.0 - 2*x0[i]*x0[i]) +
                    //     (1.0-tmpIN[i]);
                    //   }
                }
                break;

                case eTransElectric:
                {
                    m_fields[0]->GenerateElementVector(
                        m_ElemtGroup1, m_varepsilon[0], 1.0, epsvec[0]);
                    m_fields[0]->GenerateElementVector(
                        m_ElemtGroup1, m_varepsilon[1], 1.0, epsvec[1]);
                    m_fields[0]->GenerateElementVector(m_ElemtGroup1, m_mu[2],
                                                       1.0, muvec[2]);

                    // // ONLY FOR VARIABLE ANISOTROPY TEST
                    // int nq  = GetTotPoints();

                    // Array<OneD, NekDouble> tmpIN(nq);

                    // m_fields[0]->GenerateElementVector(m_ElemtGroup1, 1.0,
                    // 0.0, tmpIN);

                    // Array<OneD, NekDouble> x0(nq);
                    // Array<OneD, NekDouble> x1(nq);
                    // Array<OneD, NekDouble> x2(nq);

                    // m_fields[0]->GetCoords(x0,x1,x2);

                    // for (int i=0; i<nq; i++)
                    //   {
                    //     epsvec[1][i] = tmpIN[i]*(1.0 - 2*x0[i]*x0[i]) +
                    //     (1.0-tmpIN[i]);
                    //   }
                }
                break;

                default:
                    break; // Pol
            }
        }

        default:
            break; // TestType
    }
}

void MMFMaxwell::ComputeMaterialOpticalCloak(
    const Array<OneD, const NekDouble> &radvec,
    Array<OneD, Array<OneD, NekDouble>> &epsvec,
    Array<OneD, Array<OneD, NekDouble>> &muvec, const bool Dispersion)
{
    boost::ignore_unused(muvec, Dispersion);

    int nq = GetNpoints();

    // Cloaking metamaterial
    // \varepsilon_\theta = (b/(b-a))^2
    // \varepsilon_r = (b/(b-a))^2 ((r-a)/r)^2
    // \mu_z = 1.0m_CloakingOn

    NekDouble m_b = 1.0 / 0.314;
    NekDouble m_a = 1.0;

    NekDouble m_adel = m_a - m_Cloakraddelta;

    NekDouble boveradel = m_b / (m_b - m_adel);

    Array<OneD, NekDouble> Cloakregion(nq, 0.0);

    NekDouble la = m_MMFfactors[0];
    NekDouble lb = m_MMFfactors[1];

    NekDouble ExactCloakArea;
    if (fabs(la * lb - 1.0) < 0.00001)
    {
        ExactCloakArea = m_pi * (m_b * m_b - m_a * m_a);
    }

    else
    {
        ExactCloakArea = m_pi * (3.0 * lb * 3.0 * la - lb * la);
    }

    m_fields[0]->GenerateElementVector(m_ElemtGroup1, 1.0, 0.0, Cloakregion);

    ExactCloakArea = ExactCloakArea - (m_fields[0]->PhysIntegral(Cloakregion));
    std::cout << "*** Error of Cloakregion area = " << ExactCloakArea
              << std::endl;

    NekDouble ratio;

    epsvec[0] = Array<OneD, NekDouble>(nq, 1.0);
    epsvec[1] = Array<OneD, NekDouble>(nq, 1.0);
    m_wp2     = Array<OneD, NekDouble>(nq, 0.0);
    for (int i = 0; i < nq; ++i)
    {
        if (Cloakregion[i] > 0)
        {
            // relrad = m_a +
            // (m_b-m_adel)*(radvec[i]-m_adel)/(Cloakradmax-m_adel);
            ratio = (radvec[i] - m_adel) / radvec[i];

            epsvec[0][i] = boveradel * boveradel;
            if (m_DispersiveCloak)
            {
                epsvec[1][i] = 1.0;
                m_wp2[i]     = m_Incfreq * m_Incfreq *
                               (1.0 - boveradel * boveradel * ratio * ratio) +
                           m_wp2Tol;
            }

            else
            {
                epsvec[1][i] = boveradel * boveradel * (ratio * ratio);
            }
        }
    }
}

void MMFMaxwell::ComputeMaterialMicroWaveCloak(
    const Array<OneD, const NekDouble> &radvec,
    Array<OneD, Array<OneD, NekDouble>> &epsvec,
    Array<OneD, Array<OneD, NekDouble>> &muvec)
{
    int nq = GetNpoints();

    NekDouble m_b = 2.67;
    NekDouble m_a = 1.33;
    NekDouble m_adel;

    m_adel = m_a - m_Cloakraddelta;

    Array<OneD, NekDouble> Cloakregion(nq, 0.0);
    NekDouble ExactCloakArea = m_pi * (m_b * m_b - m_a * m_a);
    m_fields[0]->GenerateElementVector(m_ElemtGroup1, 1.0, 0.0, Cloakregion);

    if (m_ElemtGroup0 > 0)
    {
        Array<OneD, NekDouble> Vacregion(nq, 0.0);
        m_fields[0]->GenerateElementVector(m_ElemtGroup0, 1.0, 0.0, Vacregion);

        Vmath::Vsub(nq, Cloakregion, 1, Vacregion, 1, Cloakregion, 1);
    }

    ExactCloakArea = ExactCloakArea - (m_fields[0]->PhysIntegral(Cloakregion));
    std::cout << "*** Error of Cloakregion area = " << ExactCloakArea
              << std::endl;

    epsvec[0] = Array<OneD, NekDouble>(nq, 1.0);
    epsvec[1] = Array<OneD, NekDouble>(nq, 1.0);

    muvec[0] = Array<OneD, NekDouble>(nq, 1.0);
    muvec[1] = Array<OneD, NekDouble>(nq, 1.0);
    for (int i = 0; i < nq; ++i)
    {
        if (Cloakregion[i] > 0)
        {
            // relrad = m_a +
            // (m_b-m_a)*(radvec[i]-Cloakradmin)/(Cloakradmax-Cloakradmin);
            // ratio = (relrad - m_a + m_Cloakraddelta)/relrad;

            epsvec[0][i] = radvec[i] / (radvec[i] - m_adel);
            epsvec[1][i] = (radvec[i] - m_adel) / radvec[i];
            muvec[2][i]  = (m_b / (m_b - m_adel)) * (m_b / (m_b - m_adel)) *
                          (radvec[i] - m_adel) / radvec[i];

            muvec[0][i]  = epsvec[0][i];
            muvec[1][i]  = epsvec[1][i];
            epsvec[2][i] = muvec[2][i];
        }
    }
}

void MMFMaxwell::AddPML(
    const Array<OneD, const Array<OneD, NekDouble>> &physarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray)
{
    int nq = m_fields[0]->GetTotPoints();
    Array<OneD, NekDouble> tmp(nq);

    Array<OneD, NekDouble> Sigma0plus1Neg(nq);
    Array<OneD, NekDouble> Sigma0minus1(nq);
    Array<OneD, NekDouble> Sigma1minus0(nq);

    Vmath::Vsub(nq, &m_SigmaPML[1][0], 1, &m_SigmaPML[0][0], 1,
                &Sigma1minus0[0], 1);
    Vmath::Vsub(nq, &m_SigmaPML[0][0], 1, &m_SigmaPML[1][0], 1,
                &Sigma0minus1[0], 1);
    Vmath::Vadd(nq, &m_SigmaPML[0][0], 1, &m_SigmaPML[1][0], 1,
                &Sigma0plus1Neg[0], 1);
    Vmath::Neg(nq, Sigma0plus1Neg, 1);

    switch (m_PolType)
    {
        case eTransMagnetic:
        {
            int indxH0 = 0;
            int indxH1 = 1;
            int indxE2 = 2;
            int indxQ0 = 3;
            int indxQ1 = 4;
            int indxP2 = 5;

            // dH0/dt: Add (sigma_0 - \sigma_1) H0 + Q0
            Vmath::Vvtvp(nq, &Sigma0minus1[0], 1, &physarray[indxH0][0], 1,
                         &outarray[indxH0][0], 1, &outarray[indxH0][0], 1);
            Vmath::Vadd(nq, &physarray[indxQ0][0], 1, &outarray[indxH0][0], 1,
                        &outarray[indxH0][0], 1);

            // dH1/dt: Add (sigma_1 - \sigma_0) H1 + Q1
            Vmath::Vvtvp(nq, &Sigma1minus0[0], 1, &physarray[indxH1][0], 1,
                         &outarray[indxH1][0], 1, &outarray[indxH1][0], 1);
            Vmath::Vadd(nq, &physarray[indxQ1][0], 1, &outarray[indxH1][0], 1,
                        &outarray[indxH1][0], 1);

            // dHz/dt: Add -(\sigma_0 + \sigma_1) Ez  + Pz
            Vmath::Vvtvp(nq, &Sigma0plus1Neg[0], 1, &physarray[indxE2][0], 1,
                         &outarray[indxE2][0], 1, &outarray[indxE2][0], 1);
            Vmath::Vadd(nq, &physarray[indxP2][0], 1, &outarray[indxE2][0], 1,
                        &outarray[indxE2][0], 1);

            // dQ0/dt: Assign -\sigma_0 * ( Q0 + (\sigma_0 - \sigma_1) * H0 )
            Vmath::Vvtvp(nq, &Sigma0minus1[0], 1, &physarray[indxH0][0], 1,
                         &physarray[indxQ0][0], 1, &outarray[indxQ0][0], 1);
            Vmath::Vmul(nq, &m_SigmaPML[0][0], 1, &outarray[indxQ0][0], 1,
                        &outarray[indxQ0][0], 1);
            Vmath::Neg(nq, &outarray[indxQ0][0], 1);

            // dQ1/dt: Assign -\sigma_1 * ( Q1 + (\sigma_1 - \sigma_0) * H1 )
            Vmath::Vvtvp(nq, &Sigma1minus0[0], 1, &physarray[indxH1][0], 1,
                         &physarray[indxQ1][0], 1, &outarray[indxQ1][0], 1);
            Vmath::Vmul(nq, &m_SigmaPML[1][0], 1, &outarray[indxQ1][0], 1,
                        &outarray[indxQ1][0], 1);
            Vmath::Neg(nq, &outarray[indxQ1][0], 1);

            if (m_DispersiveCloak)
            {
                Vmath::Vvtvp(nq, &m_wp2[0], 1, &physarray[indxH1][0], 1,
                             &outarray[indxQ1][0], 1, &outarray[indxQ1][0], 1);
            }

            // dP3/dt: Assign - \sigma_1 * \sigma_2 * E_z
            Vmath::Vmul(nq, &m_SigmaPML[0][0], 1, &m_SigmaPML[1][0], 1,
                        &outarray[indxP2][0], 1);
            Vmath::Vmul(nq, &physarray[indxE2][0], 1, &outarray[indxP2][0], 1,
                        &outarray[indxP2][0], 1);
            Vmath::Neg(nq, &outarray[indxP2][0], 1);
        }
        break;

        case eTransElectric:
        {
            int indxE0 = 0;
            int indxE1 = 1;
            int indxH2 = 2;
            int indxQ0 = 3;
            int indxQ1 = 4;
            int indxP2 = 5;

            // dE0/dt: Add (sigma_0 - \sigma_1) E0 - Q0
            Vmath::Vvtvp(nq, &Sigma0minus1[0], 1, &physarray[indxE0][0], 1,
                         &outarray[indxE0][0], 1, &outarray[indxE0][0], 1);
            Vmath::Vsub(nq, &outarray[indxE0][0], 1, &physarray[indxQ0][0], 1,
                        &outarray[indxE0][0], 1);

            // dE1/dt: Add (sigma_1 - \sigma_0) E1 - Q1
            Vmath::Vvtvp(nq, &Sigma1minus0[0], 1, &physarray[indxE1][0], 1,
                         &outarray[indxE1][0], 1, &outarray[indxE1][0], 1);
            Vmath::Vsub(nq, &outarray[indxE1][0], 1, &physarray[indxQ1][0], 1,
                        &outarray[indxE1][0], 1);

            // dHz/dt: Add -(\sigma_0 + \sigma_1) Hz  - Pz
            Vmath::Vvtvp(nq, &Sigma0plus1Neg[0], 1, &physarray[indxH2][0], 1,
                         &outarray[indxH2][0], 1, &outarray[indxH2][0], 1);
            Vmath::Vsub(nq, &outarray[indxH2][0], 1, &physarray[indxP2][0], 1,
                        &outarray[indxH2][0], 1);

            // dQ0/dt: Assign -\sigma_0 * ( Q0 + (\sigma_1 - \sigma_0) * E0 )
            Vmath::Vvtvp(nq, &Sigma1minus0[0], 1, &physarray[indxE0][0], 1,
                         &physarray[indxQ0][0], 1, &outarray[indxQ0][0], 1);
            Vmath::Vmul(nq, &m_SigmaPML[0][0], 1, &outarray[indxQ0][0], 1,
                        &outarray[indxQ0][0], 1);
            Vmath::Neg(nq, &outarray[indxQ0][0], 1);

            // dQ1/dt: Assign -\sigma_1 * ( Q1 + (\sigma_0 - \sigma_1) * E1 )
            Vmath::Vvtvp(nq, &Sigma0minus1[0], 1, &physarray[indxE1][0], 1,
                         &physarray[indxQ1][0], 1, &outarray[indxQ1][0], 1);
            Vmath::Vmul(nq, &m_SigmaPML[1][0], 1, &outarray[indxQ1][0], 1,
                        &outarray[indxQ1][0], 1);
            Vmath::Neg(nq, &outarray[indxQ1][0], 1);

            if (m_DispersiveCloak)
            {
                Vmath::Vvtvp(nq, &m_wp2[0], 1, &physarray[indxE1][0], 1,
                             &outarray[indxQ1][0], 1, &outarray[indxQ1][0], 1);
            }

            // dP3/dt: Assign  \sigma_1 * \sigma_2 * H_z
            Vmath::Vmul(nq, &m_SigmaPML[0][0], 1, &m_SigmaPML[1][0], 1,
                        &outarray[indxP2][0], 1);
            Vmath::Vmul(nq, &physarray[indxH2][0], 1, &outarray[indxP2][0], 1,
                        &outarray[indxP2][0], 1);
        }
        break;

        default:
            break;
    }
}

void MMFMaxwell::Checkpoint_TotalFieldOutput(
    const int n, const NekDouble time,
    const Array<OneD, const Array<OneD, NekDouble>> &fieldphys)
{
    int nvar    = m_fields.size();
    int nq      = m_fields[0]->GetTotPoints();
    int ncoeffs = m_fields[0]->GetNcoeffs();

    std::string outname =
        m_sessionName + "Tot_" + boost::lexical_cast<std::string>(n) + ".chk";

    std::vector<std::string> variables(nvar);
    std::vector<Array<OneD, NekDouble>> fieldcoeffs(nvar);

    for (int i = 0; i < nvar; ++i)
    {
        fieldcoeffs[i] = Array<OneD, NekDouble>(ncoeffs);
    }

    Array<OneD, NekDouble> totfield(nq);
    for (int i = 0; i < nvar; ++i)
    {
        totfield = GetIncidentField(i, time);
        Vmath::Vadd(nq, fieldphys[i], 1, totfield, 1, totfield, 1);

        m_fields[i]->FwdTrans(totfield, fieldcoeffs[i]);
        variables[i] = m_boundaryConditions->GetVariable(i);
    }

    WriteFld(outname, m_fields[0], fieldcoeffs, variables);
}

// Only for nvar = 3
void MMFMaxwell::Checkpoint_PlotOutput(
    const int n, const Array<OneD, const Array<OneD, NekDouble>> &fieldphys)
{
    int nvar    = m_fields.size();
    int nq      = m_fields[0]->GetTotPoints();
    int ncoeffs = m_fields[0]->GetNcoeffs();

    std::string outname =
        m_sessionName + "Plot_" + boost::lexical_cast<std::string>(n) + ".chk";

    std::vector<std::string> variables(nvar);
    variables[0] = "Fx";
    variables[1] = "Fy";
    variables[2] = "Fz";

    std::vector<Array<OneD, NekDouble>> fieldcoeffs(nvar);
    for (int i = 0; i < nvar; ++i)
    {
        fieldcoeffs[i] = Array<OneD, NekDouble>(ncoeffs);
    }

    Array<OneD, NekDouble> tmp(nq);
    for (int k = 0; k < m_spacedim; ++k)
    {
        Vmath::Vmul(nq, &fieldphys[0][0], 1, &m_movingframes[0][k * nq], 1,
                    &tmp[0], 1);
        Vmath::Vvtvp(nq, &fieldphys[1][0], 1, &m_movingframes[1][k * nq], 1,
                     &tmp[0], 1, &tmp[0], 1);

        m_fields[k]->FwdTrans(tmp, fieldcoeffs[k]);
    }

    WriteFld(outname, m_fields[0], fieldcoeffs, variables);
}

void MMFMaxwell::Checkpoint_TotPlotOutput(
    const int n, const NekDouble time,
    const Array<OneD, const Array<OneD, NekDouble>> &fieldphys)
{
    int nvar    = m_fields.size();
    int nq      = m_fields[0]->GetTotPoints();
    int ncoeffs = m_fields[0]->GetNcoeffs();

    std::string outname = m_sessionName + "TotPlot_" +
                          boost::lexical_cast<std::string>(n) + ".chk";

    std::vector<std::string> variables(nvar);
    variables[0] = "Frx";
    variables[1] = "Fry";
    variables[2] = "Frz";
    for (int i = 3; i < nvar; ++i)
    {
        variables[i] = m_boundaryConditions->GetVariable(i);
    }

    std::vector<Array<OneD, NekDouble>> fieldcoeffs(nvar);
    for (int i = 0; i < nvar; ++i)
    {
        fieldcoeffs[i] = Array<OneD, NekDouble>(ncoeffs);
    }

    Array<OneD, NekDouble> tmp(nq);
    Array<OneD, NekDouble> totfield0(nq);
    Array<OneD, NekDouble> totfield1(nq);

    totfield0 = GetIncidentField(0, time);
    Vmath::Vadd(nq, fieldphys[0], 1, totfield0, 1, totfield0, 1);

    totfield1 = GetIncidentField(1, time);
    Vmath::Vadd(nq, fieldphys[1], 1, totfield1, 1, totfield1, 1);

    for (int k = 0; k < m_spacedim; ++k)
    {
        Vmath::Vmul(nq, &totfield0[0], 1, &m_movingframes[0][k * nq], 1,
                    &tmp[0], 1);
        Vmath::Vvtvp(nq, &totfield1[0], 1, &m_movingframes[1][k * nq], 1,
                     &tmp[0], 1, &tmp[0], 1);

        m_fields[k]->FwdTrans(tmp, fieldcoeffs[k]);
    }

    for (int j = 3; j < nvar; ++j)
    {
        m_fields[j]->FwdTrans(fieldphys[j], fieldcoeffs[j]);
    }

    WriteFld(outname, m_fields[0], fieldcoeffs, variables);
}

void MMFMaxwell::Checkpoint_EDFluxOutput(
    const int n, const NekDouble time,
    const Array<OneD, const Array<OneD, NekDouble>> &fieldphys)
{
    boost::ignore_unused(time);

    int nvar    = m_fields.size();
    int nq      = m_fields[0]->GetTotPoints();
    int ncoeffs = m_fields[0]->GetNcoeffs();

    std::string outname = m_sessionName + "EDFlux_" +
                          boost::lexical_cast<std::string>(n) + ".chk";

    std::vector<std::string> variables(nvar);
    variables[0] = "EDFx";
    variables[1] = "EDFy";
    variables[2] = "EDFz";
    for (int i = 3; i < nvar; ++i)
    {
        variables[i] = m_boundaryConditions->GetVariable(i);
    }

    std::vector<Array<OneD, NekDouble>> fieldcoeffs(nvar);
    for (int i = 0; i < nvar; ++i)
    {
        fieldcoeffs[i] = Array<OneD, NekDouble>(ncoeffs);
    }

    Array<OneD, NekDouble> tmp(nq);

    // TE:  H^3 (E^2 e^1 - E^1 e^2 )
    // TM: -E^3 (H^2 e^1 - H^1 e^2 )
    for (int k = 0; k < m_spacedim; ++k)
    {
        Vmath::Vmul(nq, &fieldphys[0][0], 1, &m_movingframes[1][k * nq], 1,
                    &tmp[0], 1);
        Vmath::Vvtvm(nq, &fieldphys[1][0], 1, &m_movingframes[0][k * nq], 1,
                     &tmp[0], 1, &tmp[0], 1);

        Vmath::Vmul(nq, &fieldphys[2][0], 1, &tmp[0], 1, &tmp[0], 1);

        if (m_PolType == eTransMagnetic)
        {
            Vmath::Neg(nq, tmp, 1);
        }

        m_fields[k]->FwdTrans(tmp, fieldcoeffs[k]);
    }

    WriteFld(outname, m_fields[0], fieldcoeffs, variables);
}

void MMFMaxwell::Checkpoint_EnergyOutput(
    const int n, const NekDouble time,
    const Array<OneD, const Array<OneD, NekDouble>> &fieldphys)
{
    boost::ignore_unused(time);

    int nvar    = m_fields.size();
    int nq      = m_fields[0]->GetTotPoints();
    int ncoeffs = m_fields[0]->GetNcoeffs();

    std::string outname = m_sessionName + "Energy_" +
                          boost::lexical_cast<std::string>(n) + ".chk";

    std::vector<std::string> variables(nvar);
    variables[0] = "Energy";
    variables[1] = "EnergyFlux";
    variables[2] = "Zero";
    for (int i = 3; i < nvar; ++i)
    {
        variables[i] = m_boundaryConditions->GetVariable(i);
    }

    std::vector<Array<OneD, NekDouble>> fieldcoeffs(nvar);
    for (int i = 0; i < nvar; ++i)
    {
        fieldcoeffs[i] = Array<OneD, NekDouble>(ncoeffs);
    }

    // Energy = 0.5 * ( E^2 + H^2 )
    Array<OneD, NekDouble> energy(nq, 0.0);
    Array<OneD, NekDouble> totfield(nq);
    for (int k = 0; k < m_spacedim; ++k)
    {
        // totfield = GetIncidentField(k,time);
        // Vmath::Vadd(nq, &fieldphys[k][0], 1, &totfield[0], 1, &totfield[0],
        // 1);

        Vmath::Vvtvp(nq, &fieldphys[k][0], 1, &fieldphys[k][0], 1, &energy[0],
                     1, &energy[0], 1);
    }
    Vmath::Smul(nq, 0.5, energy, 1, energy, 1);
    m_fields[0]->FwdTrans(energy, fieldcoeffs[0]);

    // EnergyFlux =  F3* sqrt( F1^2 + F2^2 )
    Array<OneD, NekDouble> energyflux(nq, 0.0);
    Array<OneD, NekDouble> Zero(nq, 0.0);
    for (int k = 0; k < 2; ++k)
    {
        // totfield = GetIncidentField(k,time);
        // Vmath::Vadd(nq, &fieldphys[k][0], 1, &totfield[0], 1, &totfield[0],
        // 1);

        Vmath::Vvtvp(nq, &fieldphys[k][0], 1, &fieldphys[k][0], 1,
                     &energyflux[0], 1, &energyflux[0], 1);
    }

    Vmath::Vsqrt(nq, energyflux, 1, energyflux, 1);
    Vmath::Vmul(nq, &fieldphys[2][0], 1, &energyflux[0], 1, &energyflux[0], 1);

    m_fields[1]->FwdTrans(energyflux, fieldcoeffs[1]);
    m_fields[2]->FwdTrans(Zero, fieldcoeffs[2]);

    WriteFld(outname, m_fields[0], fieldcoeffs, variables);
}

Array<OneD, NekDouble> MMFMaxwell::GaussianPulse(const NekDouble time,
                                                 const NekDouble Psx,
                                                 const NekDouble Psy,
                                                 const NekDouble Psz,
                                                 const NekDouble Gaussianradius)
{
    int nq      = m_fields[0]->GetTotPoints();
    int ncoeffs = m_fields[0]->GetNcoeffs();

    Array<OneD, NekDouble> x(nq);
    Array<OneD, NekDouble> y(nq);
    Array<OneD, NekDouble> z(nq);

    m_fields[0]->GetCoords(x, y, z);

    Array<OneD, NekDouble> outarray(nq, 0.0);
    Array<OneD, NekDouble> tmpc(ncoeffs);
    NekDouble rad;

    NekDouble SFradius = m_PSduration * 0.1;
    NekDouble SmoothFactor =
        1.0 / (1.0 + exp(-0.5 * (time - m_PSduration) / SFradius));

    for (int j = 0; j < nq; ++j)
    {
        rad = sqrt((x[j] - Psx) * (x[j] - Psx) + (y[j] - Psy) * (y[j] - Psy) +
                   (z[j] - Psz) * (z[j] - Psz));
        outarray[j] = m_PSstrength * SmoothFactor * exp(-1.0 * (rad / Gaussianradius) *
                                         (rad / Gaussianradius));
    }

    // m_fields[0]->FwdTrans_IterPerExp(outarray, tmpc);
    m_fields[0]->FwdTransLocalElmt(outarray, tmpc);
    m_fields[0]->BwdTrans(tmpc, outarray);

    return outarray;
}

Array<OneD, NekDouble> MMFMaxwell::EvaluateCoriolis()
{
    int nq = GetTotPoints();

    NekDouble m_Omega = 1.5486 * 0.000001;

    NekDouble x0j, x1j, x2j;
    NekDouble sin_theta, cos_theta, sin_varphi, cos_varphi;

    Array<OneD, NekDouble> x(nq);
    Array<OneD, NekDouble> y(nq);
    Array<OneD, NekDouble> z(nq);

    m_fields[0]->GetCoords(x, y, z);

    Array<OneD, NekDouble> outarray(nq, 0.0);
    for (int j = 0; j < nq; ++j)
    {
        x0j = x[j];
        x1j = y[j];
        x2j = z[j];

        // CartesianToSpherical(x0j, x1j, x2j, sin_varphi, cos_varphi, sin_theta,
        //                      cos_theta);

        // outarray[j] = 2.0 * m_Omega * sin_theta;

        CartesianToNewSpherical(x0j, x1j, x2j, sin_varphi, cos_varphi, sin_theta, cos_theta);

        outarray[j] = 2.0 * m_Omega * cos_theta;

    }

    return outarray;
}

void MMFMaxwell::AddCoriolis(Array<OneD, Array<OneD, NekDouble>> &physarray,
                             Array<OneD, Array<OneD, NekDouble>> &outarray)
{
    int nq = physarray[0].size();

    Array<OneD, NekDouble> tmp(nq);

    int indx;
    for (int j = 0; j < m_shapedim; ++j)
    {
        if (j == 0)
        {
            indx = 2;
        }

        else if (j == 1)
        {
            indx = 1;
        }

        Vmath::Vmul(nq, m_coriolis, 1, physarray[indx], 1, tmp, 1);

        switch (m_PolType)
        {
            case eTransMagnetic:
            {
                Vmath::Vmul(nq, m_muvec[indx], 1, tmp, 1, tmp, 1);
            }
            break;

            case eTransElectric:
            {
                Vmath::Vmul(nq, m_epsvec[indx], 1, tmp, 1, tmp, 1);
            }
            break;

            default:
                break;
        }

        if (j == 1)
        {
            Vmath::Neg(nq, tmp, 1);
        }

        Vmath::Vadd(nq, tmp, 1, outarray[j], 1, outarray[j], 1);
    }
}

Array<OneD, NekDouble> MMFMaxwell::ComputeRadCloak(const int CloakNlayer)
{
    int nq = GetNpoints();

    NekDouble m_b = 1.0 / 0.314;
    NekDouble m_a = 1.0;

    Array<OneD, int> Layer(CloakNlayer);

    NekDouble la = m_MMFfactors[0];
    NekDouble lb = m_MMFfactors[1];
    if (CloakNlayer == 8)
    {
        // Circular 8 layer
        if (fabs(la * lb - 1.0) < 0.0001)
        {
            Layer[0] = 119;
            Layer[1] = 239;
            Layer[2] = 323;
            Layer[3] = 459;
            Layer[4] = 575;
            Layer[5] = 727;
            Layer[6] = 891;
            Layer[7] = 1059;
        }

        // Ellipsoid 8 layer
        else
        {
            Layer[0] = 115;
            Layer[1] = 219;
            Layer[2] = 335;
            Layer[3] = 459;
            Layer[4] = 607;
            Layer[5] = 767;
            Layer[6] = 951;
            Layer[7] = 1155;
        }
    }

    if (CloakNlayer == 16)
    {
        // Circular 16 layer
        if (fabs(la * lb - 1.0) < 0.0001)
        {
            Layer[0]  = 43;
            Layer[1]  = 87;
            Layer[2]  = 135;
            Layer[3]  = 187;
            Layer[4]  = 239;
            Layer[5]  = 295;
            Layer[6]  = 355;
            Layer[7]  = 415;
            Layer[8]  = 479;
            Layer[9]  = 555;
            Layer[10] = 639;
            Layer[11] = 727;
            Layer[12] = 823;
            Layer[13] = 927;
            Layer[14] = 1039;
            Layer[15] = 1159;
        }

        // Ellipsoid 8 layer
        else
        {
            Layer[0]  = 135;
            Layer[1]  = 259;
            Layer[2]  = 387;
            Layer[3]  = 523;
            Layer[4]  = 667;
            Layer[5]  = 803;
            Layer[6]  = 939;
            Layer[7]  = 1083;
            Layer[8]  = 1235;
            Layer[9]  = 1387;
            Layer[10] = 1539;
            Layer[11] = 1699;
            Layer[12] = 1867;
            Layer[13] = 2035;
            Layer[14] = 2203;
            Layer[15] = 2379;
        }
    }

    Array<OneD, NekDouble> x0(nq);
    Array<OneD, NekDouble> x1(nq);
    Array<OneD, NekDouble> x2(nq);

    m_fields[0]->GetCoords(x0, x1, x2);

    Array<OneD, NekDouble> Formertmp(nq, 0.0);
    Array<OneD, NekDouble> Currenttmp(nq, 0.0);
    Array<OneD, NekDouble> Cloakregion(nq, 0.0);

    m_fields[0]->GenerateElementVector(Layer[0], 1.0, 0.0, Currenttmp);
    Vmath::Vcopy(nq, Currenttmp, 1, Cloakregion, 1);
    for (int i = 1; i < CloakNlayer; ++i)
    {
        m_fields[0]->GenerateElementVector(Layer[i], 1.0, 0.0, Currenttmp);
        m_fields[0]->GenerateElementVector(Layer[i - 1], 1.0, 0.0, Formertmp);

        Vmath::Vsub(nq, Currenttmp, 1, Formertmp, 1, Currenttmp, 1);

        Vmath::Svtvp(nq, 1.0 * (i + 1), Currenttmp, 1, Cloakregion, 1,
                     Cloakregion, 1);
    }

    Array<OneD, NekDouble> radvec(nq);

    for (int i = 0; i < nq; ++i)
    {
        switch (m_MMFdir)
        {
            case SpatialDomains::eTangentCircular:
            {
                if ((Cloakregion[i] > 0) && (CloakNlayer > 0))
                {
                    radvec[i] =
                        1.0 +
                        (m_b - m_a) / CloakNlayer * (Cloakregion[i] - 0.5);
                }

                else
                {
                    radvec[i] =
                        sqrt(x0[i] * x0[i] / la / la + x1[i] * x1[i] / lb / lb);
                }
            }
            break;

            // case SpatialDomains::eTangentIrregular:
            // {
            //     radvec[i] = sqrt(2.0 * x0[i] * x0[i] +
            //                      x1[i] * x1[i] * x1[i] * x1[i] + x1[i] * x1[i]);
            // }
            // break;

            // case SpatialDomains::eTangentNonconvex:
            // {
            //     radvec[i] = sqrt(3.0 * x0[i] * x0[i] +
            //                      x1[i] * x1[i] * x1[i] * x1[i] - x1[i] * x1[i]);
            // }
            // break;

            default:
                break;
        }
    }

    return radvec;
}

void MMFMaxwell::v_GenerateSummary(SolverUtils::SummaryList &s)
{
    MMFSystem::v_GenerateSummary(s);
    SolverUtils::AddSummaryItem(
        s, "TestMaxwellType",
        TestMaxwellTypeMap[m_TestMaxwellType]);
    SolverUtils::AddSummaryItem(s, "PolType",
                                PolTypeMap[m_PolType]);
    SolverUtils::AddSummaryItem(s, "IncType",
                                IncTypeMap[m_IncType]);

    if (m_varepsilon[0] * m_varepsilon[1] * m_varepsilon[2] > 1.0)
    {
        SolverUtils::AddSummaryItem(s, "varepsilon1", m_varepsilon[0]);
        SolverUtils::AddSummaryItem(s, "varepsilon2", m_varepsilon[1]);
        SolverUtils::AddSummaryItem(s, "varepsilon3", m_varepsilon[2]);
    }

    if (m_mu[0] * m_mu[1] * m_mu[2] > 1.0)
    {
        SolverUtils::AddSummaryItem(s, "mu1", m_mu[0]);
        SolverUtils::AddSummaryItem(s, "mu2", m_mu[1]);
        SolverUtils::AddSummaryItem(s, "mu3", m_mu[2]);
    }

    if (m_boundaryforSF > 0)
    {
        SolverUtils::AddSummaryItem(s, "boundarySF", m_boundaryforSF);
    }

    if (m_ElemtGroup1 > 0)
    {
        SolverUtils::AddSummaryItem(s, "CloakNlayer", m_CloakNlayer);
        SolverUtils::AddSummaryItem(s, "ElemtGroup1", m_ElemtGroup1);
    }

    SolverUtils::AddSummaryItem(s, "DivergenceRestore", m_DivergenceRestore);
    SolverUtils::AddSummaryItem(s, "AddRotation", m_AddRotation);

    if (m_AddPML > 0)
    {
        SolverUtils::AddSummaryItem(s, "AddPML", m_AddPML);
        SolverUtils::AddSummaryItem(s, "PMLelement", m_PMLelement);
        SolverUtils::AddSummaryItem(s, "RecPML", m_RecPML);
        SolverUtils::AddSummaryItem(s, "PMLorder", m_PMLorder);
        SolverUtils::AddSummaryItem(s, "PMLthickness", m_PMLthickness);
        SolverUtils::AddSummaryItem(s, "PMLstart", m_PMLstart);
        SolverUtils::AddSummaryItem(s, "PMLmaxsigma", m_PMLmaxsigma);
    }

    if (m_SourceType)
    {
        SolverUtils::AddSummaryItem(s, "SourceType",
                                    SourceTypeMap[m_SourceType]);
        SolverUtils::AddSummaryItem(s, "Psx", m_Psx);
        SolverUtils::AddSummaryItem(s, "Psy", m_Psy);
        SolverUtils::AddSummaryItem(s, "Psz", m_Psz);
        SolverUtils::AddSummaryItem(s, "PSduration", m_PSduration);
        SolverUtils::AddSummaryItem(s, "PSstrength", m_PSstrength);
        SolverUtils::AddSummaryItem(s, "Gaussianradius", m_Gaussianradius);
    }

    if (m_CloakType)
    {
        SolverUtils::AddSummaryItem(s, "CloakType", CloakTypeMap[m_CloakType]);
        SolverUtils::AddSummaryItem(s, "DispersiveCloak", m_DispersiveCloak);
        SolverUtils::AddSummaryItem(s, "CloakNlayer", m_CloakNlayer);
        SolverUtils::AddSummaryItem(s, "Cloakraddelta", m_Cloakraddelta);
    }
}

// dim = dimension, pol = polarization, 0 = TM, 1 = TE.
void MMFMaxwell::ComputeZimYim(
    const Array<OneD, const Array<OneD, NekDouble>> &epsvec,
    const Array<OneD, const Array<OneD, NekDouble>> &muvec,
    Array<OneD, Array<OneD, NekDouble>> &ZimFwd,
    Array<OneD, Array<OneD, NekDouble>> &ZimBwd,
    Array<OneD, Array<OneD, NekDouble>> &YimFwd,
    Array<OneD, Array<OneD, NekDouble>> &YimBwd)

{
    int nTraceNumPoints = GetTraceNpoints();

    switch (m_TestMaxwellType)
    {
        case eMaxwell1D:
        case eScatField1D:
        {
            Array<OneD, NekDouble> Fwdeps(nTraceNumPoints, 1.0);
            Array<OneD, NekDouble> Bwdeps(nTraceNumPoints, 1.0);
            Array<OneD, NekDouble> Fwdmu(nTraceNumPoints, 1.0);
            Array<OneD, NekDouble> Bwdmu(nTraceNumPoints, 1.0);
            m_fields[0]->GetFwdBwdTracePhys(epsvec[0], Fwdeps, Bwdeps);
            m_fields[0]->GetFwdBwdTracePhys(muvec[0], Fwdeps, Bwdeps);

            CopyBoundaryTrace(Fwdeps, Bwdeps, SolverUtils::eFwdEQBwd, 0);
            CopyBoundaryTrace(Fwdmu, Bwdmu, SolverUtils::eFwdEQBwd, 1);

            ZimFwd = Array<OneD, Array<OneD, NekDouble>>(1);
            ZimBwd = Array<OneD, Array<OneD, NekDouble>>(1);
            YimFwd = Array<OneD, Array<OneD, NekDouble>>(1);
            YimBwd = Array<OneD, Array<OneD, NekDouble>>(1);

            ZimFwd[0] = Array<OneD, NekDouble>(nTraceNumPoints, 1.0);
            ZimBwd[0] = Array<OneD, NekDouble>(nTraceNumPoints, 1.0);
            YimFwd[0] = Array<OneD, NekDouble>(nTraceNumPoints, 1.0);
            YimBwd[0] = Array<OneD, NekDouble>(nTraceNumPoints, 1.0);

            // ZimFwd = sqrt( muFwd / epsFwd),  ZimBwd = sqrt( muBwd / epsBwd)
            for (int i = 0; i < nTraceNumPoints; ++i)
            {
                m_ZimFwd[0][i] = sqrt(Fwdmu[i] / Fwdeps[i]);
                m_ZimBwd[0][i] = sqrt(Bwdmu[i] / Bwdeps[i]);

                m_YimFwd[0][i] = 1.0 / m_ZimFwd[0][i];
                m_YimBwd[0][i] = 1.0 / m_ZimBwd[0][i];
            }

            std::cout << "*** ZimFwd = " << RootMeanSquare(m_ZimFwd[0])
                      << ", ZimBwd = " << RootMeanSquare(m_ZimBwd[0])
                      << ", YimFwd = " << RootMeanSquare(m_YimFwd[0])
                      << ", YimBwd = " << RootMeanSquare(m_YimBwd[0])
                      << std::endl;
        }
        break;

        case eTestMaxwell2DPEC:
        case eTestMaxwell2DPECAVGFLUX:
        case eTestMaxwell2DPMC:
        case eScatField2D:
        case eTotField2D:
        case eMaxwellSphere:
        case eELF2DSurface:
        {
            switch (m_PolType)
            {
                case eTransMagnetic:
                {
                    Array<OneD, NekDouble> Fwdmu1(nTraceNumPoints);
                    Array<OneD, NekDouble> Bwdmu1(nTraceNumPoints);
                    Array<OneD, NekDouble> Fwdmu2(nTraceNumPoints);
                    Array<OneD, NekDouble> Bwdmu2(nTraceNumPoints);
                    Array<OneD, NekDouble> Fwdeps3(nTraceNumPoints);
                    Array<OneD, NekDouble> Bwdeps3(nTraceNumPoints);

                    m_fields[0]->GetFwdBwdTracePhys(muvec[0], Fwdmu1, Bwdmu1);
                    m_fields[0]->GetFwdBwdTracePhys(muvec[1], Fwdmu2, Bwdmu2);
                    m_fields[0]->GetFwdBwdTracePhys(epsvec[2], Fwdeps3,
                                                    Bwdeps3);

                    CopyBoundaryTrace(Fwdmu1, Bwdmu1, SolverUtils::eFwdEQBwd,
                                      0);
                    CopyBoundaryTrace(Fwdmu2, Bwdmu2, SolverUtils::eFwdEQBwd,
                                      1);
                    CopyBoundaryTrace(Fwdeps3, Bwdeps3, SolverUtils::eFwdEQBwd,
                                      2);

                    for (int i = 0; i < nTraceNumPoints; ++i)
                    {
                        ZimFwd[0][i] = sqrt(Fwdmu2[i] / Fwdeps3[i]);
                        ZimBwd[0][i] = sqrt(Bwdmu2[i] / Bwdeps3[i]);

                        YimFwd[0][i] = 1.0 / ZimFwd[0][i];
                        YimBwd[0][i] = 1.0 / ZimBwd[0][i];

                        ZimFwd[1][i] = sqrt(Fwdmu1[i] / Fwdeps3[i]);
                        ZimBwd[1][i] = sqrt(Bwdmu1[i] / Bwdeps3[i]);

                        YimFwd[1][i] = 1.0 / ZimFwd[1][i];
                        YimBwd[1][i] = 1.0 / ZimBwd[1][i];
                    }
                }
                break; // eTransMagnetic

                case eTransElectric:
                {
                    Array<OneD, NekDouble> Fwdeps1(nTraceNumPoints);
                    Array<OneD, NekDouble> Bwdeps1(nTraceNumPoints);
                    Array<OneD, NekDouble> Fwdeps2(nTraceNumPoints);
                    Array<OneD, NekDouble> Bwdeps2(nTraceNumPoints);
                    Array<OneD, NekDouble> Fwdmu3(nTraceNumPoints);
                    Array<OneD, NekDouble> Bwdmu3(nTraceNumPoints);

                    m_fields[0]->GetFwdBwdTracePhys(epsvec[0], Fwdeps1,
                                                    Bwdeps1);
                    m_fields[0]->GetFwdBwdTracePhys(epsvec[1], Fwdeps2,
                                                    Bwdeps2);
                    m_fields[0]->GetFwdBwdTracePhys(muvec[2], Fwdmu3, Bwdmu3);

                    CopyBoundaryTrace(Fwdeps1, Bwdeps1, SolverUtils::eFwdEQBwd,
                                      0);
                    CopyBoundaryTrace(Fwdeps2, Bwdeps2, SolverUtils::eFwdEQBwd,
                                      1);
                    CopyBoundaryTrace(Fwdmu3, Bwdmu3, SolverUtils::eFwdEQBwd,
                                      2);

                    for (int i = 0; i < nTraceNumPoints; ++i)
                    {
                        ZimFwd[0][i] = sqrt(Fwdmu3[i] / Fwdeps2[i]);
                        ZimBwd[0][i] = sqrt(Bwdmu3[i] / Bwdeps2[i]);

                        YimFwd[0][i] = 1.0 / ZimFwd[0][i];
                        YimBwd[0][i] = 1.0 / ZimBwd[0][i];

                        ZimFwd[1][i] = sqrt(Fwdmu3[i] / Fwdeps1[i]);
                        ZimBwd[1][i] = sqrt(Bwdmu3[i] / Bwdeps1[i]);

                        YimFwd[1][i] = 1.0 / ZimFwd[1][i];
                        YimBwd[1][i] = 1.0 / ZimBwd[1][i];
                    }
                }
                break; // eTransELectric

                default:
                    break;
            } // PolType

            std::cout << "*** ZimFwd0 = [ "
                      << Vmath::Vmin(nTraceNumPoints, ZimFwd[0], 1) << " , "
                      << Vmath::Vmax(nTraceNumPoints, ZimFwd[0], 1)
                      << " ], ZimBwd0 = [ "
                      << Vmath::Vmin(nTraceNumPoints, ZimBwd[0], 1) << " , "
                      << Vmath::Vmax(nTraceNumPoints, ZimBwd[0], 1) << " ] "
                      << std::endl;
            std::cout << "*** ZimFwd1 = [ "
                      << Vmath::Vmin(nTraceNumPoints, ZimFwd[1], 1) << " , "
                      << Vmath::Vmax(nTraceNumPoints, ZimFwd[1], 1)
                      << " ], ZimBwd1 = [ "
                      << Vmath::Vmin(nTraceNumPoints, ZimBwd[1], 1) << " , "
                      << Vmath::Vmax(nTraceNumPoints, ZimBwd[1], 1) << " ] "
                      << std::endl;
        }
        break; // eMaxwell2D

        default:
            break;
    } // TestMaxwellType
}


void MMFMaxwell::AdddedtMaxwell(
    const Array<OneD, const Array<OneD, NekDouble>> &physarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray)
{
    int i, j;
    int nq = GetTotPoints();

    // m_dedxi_cdot_e[m][j][n][] = de^m / d \xi^j \cdot e^n
    Array<OneD, NekDouble> dedtej(nq);
    Array<OneD, NekDouble> de1dtcdotej(nq);
    Array<OneD, NekDouble> de2dtcdotej(nq);

    Array<OneD, NekDouble> normH(nq);
    Array<OneD, NekDouble> NormedH1(nq, 0.0);
    Array<OneD, NekDouble> NormednegH2(nq, 0.0);

    Vmath::Vmul(nq, &physarray[0][0], 1, &physarray[0][0], 1, &normH[0], 1);
    Vmath::Vvtvp(nq, &physarray[1][0], 1, &physarray[1][0], 1, &normH[0], 1,
                 &normH[0], 1);
    Vmath::Vsqrt(nq, normH, 1, normH, 1);

    NekDouble Tol = 0.001;
    for (i = 0; i < nq; ++i)
    {
        if (normH[i] > Tol)
        {
            NormedH1[i]    = physarray[0][i] / normH[i];
            NormednegH2[i] = -1.0 * physarray[1][i] / normH[i];
        }
    }

    for (j = 0; j < m_shapedim; ++j)
    {
        // Compute de1 / dt \cdot ej = (-H2 de^1/d\xi1 \cdot e^j + H1 de^1/d\xi2
        // \cdot e^j) / sqrt{ H1^2 + H2^2 }
        Vmath::Vmul(nq, &NormednegH2[0], 1, &m_dedxi_cdot_e[0][0][j][0], 1,
                    &de1dtcdotej[0], 1);
        Vmath::Vvtvp(nq, &NormedH1[0], 1, &m_dedxi_cdot_e[0][1][j][0], 1,
                     &de1dtcdotej[0], 1, &de1dtcdotej[0], 1);

        // Compute de2 / dt \cdot ej = (-H2 de2/d\xi1 \cdot e^j + H1 de2/d\xi2
        // \cdot e^j) / sqrt{ H1^2 + H2^2 }
        Vmath::Vmul(nq, &NormednegH2[0], 1, &m_dedxi_cdot_e[1][0][j][0], 1,
                    &de2dtcdotej[0], 1);
        Vmath::Vvtvp(nq, &NormedH1[0], 1, &m_dedxi_cdot_e[1][1][j][0], 1,
                     &de2dtcdotej[0], 1, &de2dtcdotej[0], 1);

        // Add dedt component: (H1 (de1/dt) + H2 (de2/dt) ) \cdot ej
        Vmath::Vmul(nq, &physarray[0][0], 1, &de1dtcdotej[0], 1, &dedtej[0], 1);
        Vmath::Vvtvp(nq, &physarray[1][0], 1, &de2dtcdotej[0], 1, &dedtej[0], 1,
                     &dedtej[0], 1);

        Vmath::Neg(nq, dedtej, 1);

        switch (m_PolType)
        {
            case eTransMagnetic:
            {
                if (j == 0)
                {
                    Vmath::Vmul(nq, m_muvec[0], 1, dedtej, 1, dedtej, 1);
                }

                else if (j == 1)
                {
                    Vmath::Vmul(nq, m_muvec[1], 1, dedtej, 1, dedtej, 1);
                }
            }
            break;

            case eTransElectric:
            {
                if (j == 0)
                {
                    Vmath::Vmul(nq, m_epsvec[0], 1, dedtej, 1, dedtej, 1);
                }

                else if (j == 1)
                {
                    Vmath::Vmul(nq, m_epsvec[1], 1, dedtej, 1, dedtej, 1);
                }
            }
            break;

            default:
                break;
        }

        Vmath::Vadd(nq, &dedtej[0], 1, &outarray[j][0], 1, &outarray[j][0], 1);
    }
}

Array<OneD, NekDouble> MMFMaxwell::GetIncidentField(const int var,
                                                   const NekDouble time)
{
    int nq = m_fields[0]->GetNpoints();

    Array<OneD, NekDouble> x0(nq);
    Array<OneD, NekDouble> x1(nq);
    Array<OneD, NekDouble> x2(nq);

    m_fields[0]->GetCoords(x0, x1, x2);

    // GetSmoothFactor such that wave propages from the left to the object.
    // a = 0.1, ta = 1, f = 1.0./(1.0 + exp( -0.5.*(time-ta)/a ));
    Array<OneD, NekDouble> SmoothFactor(nq, 1.0);
    int m_SmoothFactor = 1;
    switch (m_SmoothFactor)
    {
        case 0:
        {
            for (int i = 0; i < nq; i++)
            {
                SmoothFactor[i] = 1.0 / (1.0 + exp(-1.0 * (time - 1.0) / 0.1));
            }
        }
        break;

        case 1:
        {
            NekDouble xmin = Vmath::Vmin(nq, x0, 1);
            NekDouble xp;
            for (int i = 0; i < nq; i++)
            {
                xp = x0[i] - xmin - time - m_SFinit;
                if (xp > 0.0)
                {
                    SmoothFactor[i] =
                        2.0 / (1.0 + exp(0.5 * (sqrt(xp * xp) - 0.1)));
                }

                else
                {
                    SmoothFactor[i] = 1.0;
                }
            }
        }
        break;

        default:
            break;
    }

    // Generate a factor for smoothly increasing wave
    Array<OneD, NekDouble> F1(nq);
    Array<OneD, NekDouble> F2(nq);
    Array<OneD, NekDouble> F3(nq);

    Array<OneD, NekDouble> dF1dt(nq);
    Array<OneD, NekDouble> dF2dt(nq);
    Array<OneD, NekDouble> dF3dt(nq);

    Array<OneD, NekDouble> F1int(nq);
    Array<OneD, NekDouble> F2int(nq);
    Array<OneD, NekDouble> F3int(nq);

    Array<OneD, NekDouble> outarray(nq);

    switch (m_IncType)
    {
        case ePlaneWave:
        {
            NekDouble cs, sn;
            NekDouble e1y, e2y;
            switch (m_PolType)
            {
                case eTransMagnetic:
                {
                    // H = 0 \hat{x} - e^{ikx} \hat{y}
                    // E = e^{ikx} \hat{z}
                    // outarray1 = Hr1inc
                    // outarray2 = Hr2inc
                    // outarray3 = Ezrinc
                    for (int i = 0; i < nq; i++)
                    {
                        e1y = m_movingframes[0][nq + i];
                        e2y = m_movingframes[1][nq + i];

                        cs = SmoothFactor[i] * cos(m_Incfreq * (x0[i] - time));
                        sn = SmoothFactor[i] * sin(m_Incfreq * (x0[i] - time));

                        F1[i] = -1.0 * cs * e1y;
                        F2[i] = -1.0 * cs * e2y;
                        F3[i] = cs;

                        dF1dt[i] = -1.0 * m_Incfreq * sn * e1y;
                        dF2dt[i] = -1.0 * m_Incfreq * sn * e2y;
                        dF3dt[i] = 1.0 * m_Incfreq * sn;

                        F1int[i] = (1.0 / m_Incfreq) * sn * e1y;
                        F2int[i] = (1.0 / m_Incfreq) * sn * e2y;
                        F3int[i] = (-1.0 / m_Incfreq) * sn;
                    }
                }
                break;

                case eTransElectric:
                {
                    // E = 0 \hat{x} + e^{ikx} \hat{y}
                    // H = e^{ikx} \hat{z}
                    // outarray1 = Er1inc
                    // outarray2 = Er2inc
                    // outarray3 = Hzrinc
                    // outarray1 = Ei1inc
                    // outarray2 = Ei2inc
                    // outarray3 = Hziinc
                    for (int i = 0; i < nq; i++)
                    {
                        e1y = m_movingframes[0][nq + i];
                        e2y = m_movingframes[1][nq + i];

                        cs = SmoothFactor[i] * cos(m_Incfreq * (x0[i] - time));
                        sn = SmoothFactor[i] * sin(m_Incfreq * (x0[i] - time));

                        F1[i] = cs * e1y;
                        F2[i] = cs * e2y;
                        F3[i] = cs;

                        dF1dt[i] = m_Incfreq * sn * e1y;
                        dF2dt[i] = m_Incfreq * sn * e2y;
                        dF3dt[i] = m_Incfreq * sn;

                        F1int[i] = (-1.0 / m_Incfreq) * sn * e1y;
                        F2int[i] = (-1.0 / m_Incfreq) * sn * e2y;
                        F3int[i] = (-1.0 / m_Incfreq) * sn;
                    }
                }
                break;

                default:
                    break;
            }
        }
        break;

        case ePlaneWaveImag:
        {
            NekDouble cs, sn;
            NekDouble e1y, e2y;
            switch (m_PolType)
            {
                case eTransMagnetic:
                {
                    // H = 0 \hat{x} - e^{ikx} \hat{y}
                    // E = e^{ikx} \hat{z}
                    // outarray1 = Hr1inc
                    // outarray2 = Hr2inc
                    // outarray3 = Ezrinc
                    for (int i = 0; i < nq; i++)
                    {
                        e1y = m_movingframes[0][nq + i];
                        e2y = m_movingframes[1][nq + i];

                        cs = SmoothFactor[i] * cos(m_Incfreq * (x0[i] - time));
                        sn = SmoothFactor[i] * sin(m_Incfreq * (x0[i] - time));

                        F1[i] = -1.0 * sn * e1y;
                        F2[i] = -1.0 * sn * e2y;
                        F3[i] = sn;

                        dF1dt[i] = m_Incfreq * cs * e1y;
                        dF2dt[i] = m_Incfreq * cs * e2y;
                        dF3dt[i] = -1.0 * m_Incfreq * cs;

                        F1int[i] = (-1.0 / m_Incfreq) * cs * e1y;
                        F2int[i] = (-1.0 / m_Incfreq) * cs * e2y;
                        F3int[i] = (1.0 / m_Incfreq) * cs;
                    }
                }
                break;

                case eTransElectric:
                {
                    // E = 0 \hat{x} + e^{ikx} \hat{y}
                    // H = e^{ikx} \hat{z}
                    // outarray1 = Er1inc
                    // outarray2 = Er2inc
                    // outarray3 = Hzrinc
                    // outarray1 = Ei1inc
                    // outarray2 = Ei2inc
                    // outarray3 = Hziinc
                    for (int i = 0; i < nq; i++)
                    {
                        e1y = m_movingframes[0][nq + i];
                        e2y = m_movingframes[1][nq + i];

                        cs = SmoothFactor[i] * cos(m_Incfreq * (x0[i] - time));
                        sn = SmoothFactor[i] * sin(m_Incfreq * (x0[i] - time));

                        F1[i] = sn * e1y;
                        F2[i] = sn * e2y;
                        F3[i] = sn;

                        dF1dt[i] = -1.0 * m_Incfreq * cs * e1y;
                        dF2dt[i] = -1.0 * m_Incfreq * cs * e2y;
                        dF3dt[i] = -1.0 * m_Incfreq * cs;

                        F1int[i] = (1.0 / m_Incfreq) * cs * e1y;
                        F2int[i] = (1.0 / m_Incfreq) * cs * e2y;
                        F3int[i] = (1.0 / m_Incfreq) * cs;
                    }
                }
                break;

                default:
                    break;
            }
        }
        break;

        default:
            break;
    }

    switch (var)
    {
        case 0:
        {
            outarray = F1;
        }
        break;

        case 1:
        {
            outarray = F2;
        }
        break;

        case 2:
        {
            outarray = F3;
        }
        break;

        case 10:
        {
            outarray = dF1dt;
        }
        break;

        case 11:
        {
            outarray = dF2dt;
        }
        break;

        case 12:
        {
            outarray = dF3dt;
        }
        break;

        case 20:
        {
            outarray = F1int;
        }
        break;

        case 21:
        {
            outarray = F2int;
        }
        break;

        case 22:
        {
            outarray = F3int;
        }
        break;

        default:
        {
            Vmath::Zero(nq, outarray, 1);
        }
        break;
    }

    return outarray;
}


void MMFMaxwell::GetMaxwellFluxVector(
    const int var, const Array<OneD, const Array<OneD, NekDouble>> &physfield,
    Array<OneD, Array<OneD, NekDouble>> &flux)
{
    switch (m_TestMaxwellType)
    {
        case eMaxwell1D:
        case eScatField1D:
        {
            GetMaxwellFlux1D(var, physfield, flux);
        }
        break;

        case eTestMaxwell2DPEC:
        case eTestMaxwell2DPECAVGFLUX:
        case eTestMaxwell2DPMC:
        case eScatField2D:
        case eTotField2D:
        case eMaxwellSphere:
        case eELF2DSurface:
        {
            GetMaxwellFlux2D(var, physfield, flux);
        }
        break;

        default:
            break;
    }
}

void MMFMaxwell::GetMaxwellFlux1D(
    const int var, const Array<OneD, const Array<OneD, NekDouble>> &physfield,
    Array<OneD, Array<OneD, NekDouble>> &flux)
{
    int nq = m_fields[0]->GetTotPoints();

    switch (var)
    {
        case 0:
        {
            // H in flux 0
            Vmath::Vcopy(nq, physfield[1], 1, flux[0], 1);

            // E in flux 1
            Vmath::Zero(nq, flux[1], 1);
        }
        break;

        case 1:
        {
            // E in flux 0
            Vmath::Vcopy(nq, physfield[0], 1, flux[0], 1);

            // H in flux 1
            Vmath::Zero(nq, flux[1], 1);
        }
        break;
            //----------------------------------------------------

        default:
            break;
    }
}


// TM: (-E3,0), (0,E3), (-H1, H2)
// TE: (H3,0), (0,-H3), (E1, -E2)
void MMFMaxwell::GetMaxwellFlux2D(
    const int var, const Array<OneD, const Array<OneD, NekDouble>> &physfield,
    Array<OneD, Array<OneD, NekDouble>> &flux)
{
    int nq = m_fields[0]->GetTotPoints();

    NekDouble sign = 1.0;
    switch (m_PolType)
    {
            // TransMagnetic
        case eTransMagnetic:
        {
            sign = -1.0;
        }
        break;

            // TransElectric
        case eTransElectric:
        {
            sign = 1.0;
        }
        break;

        default:
            break;
    }

    switch (var)
    {
        case 0:
        {
            // -Ez in flux 1
            Vmath::Smul(nq, sign, physfield[2], 1, flux[0], 1);
            Vmath::Zero(nq, flux[1], 1);
        }
        break;

        case 1:
        {
            // Ez in flux 0
            Vmath::Zero(nq, flux[0], 1);
            Vmath::Smul(nq, -sign, physfield[2], 1, flux[1], 1);
        }
        break;

        case 2:
        {
            Vmath::Smul(nq, sign, physfield[0], 1, flux[0], 1);
            Vmath::Smul(nq, -sign, physfield[1], 1, flux[1], 1);
        }
        break;

        default:
            ASSERTL0(false, "GetFluxVector2D: illegal vector index");
    }
}

void MMFMaxwell::print_MMF(Array<OneD, Array<OneD, NekDouble>> &inarray)
{
    int Ntot = inarray.size();

    NekDouble reval = 0.0;
    for (int i = 0; i < Ntot; ++i)
    {
        std::cout << "[" << i << "] = " << inarray[2][i] << std::endl;
        // reval = reval + inarray[i]*inarray[i];
    }
    reval = sqrt(reval / Ntot);
}
}
