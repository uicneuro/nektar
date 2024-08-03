///////////////////////////////////////////////////////////////////////////////
//
// File MMFBidomain.cpp
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
// Description: MMFBidomain.
//
///////////////////////////////////////////////////////////////////////////////

#include <boost/algorithm/string.hpp>
#include <iomanip>
#include <iostream>

#include <DiffusionSolver/EquationSystems/MMFBidomain.h>
#include <DiffusionSolver/Filters/FilterCellHistoryPoints.h>
#include <DiffusionSolver/Filters/FilterCheckpointCellModel.h>
#include <SolverUtils/Driver.h>

#include <LibUtilities/BasicUtils/Timer.h>

#include <LibUtilities/BasicUtils/SessionReader.h>
#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>

#include <boost/math/special_functions/spherical_harmonic.hpp>
using namespace std;
using namespace Nektar::SolverUtils;
using namespace Nektar;

namespace Nektar
{
string MMFBidomain::className =
    SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
        "MMFBidomain", MMFBidomain::create, "MMFBidomain equation.");

MMFBidomain::MMFBidomain(const LibUtilities::SessionReaderSharedPtr &pSession,
                         const SpatialDomains::MeshGraphSharedPtr &pGraph)
    : UnsteadySystem(pSession, pGraph), MMFSystem(pSession, pGraph)
{
}

void MMFBidomain::v_InitObject()
{
    MMFSystem::v_InitObject();

    m_session->LoadParameter("Helmtau", m_Helmtau, 1.0);

    // Define ProblemType
    if (m_session->DefinesSolverInfo("TESTTYPE"))
    {
        std::string TestTypeStr;
        TestTypeStr = m_session->GetSolverInfo("TESTTYPE");
        for (int i = 0; i < (int)SIZE_TestType; ++i)
        {
            if (boost::iequals(TestTypeMap[i], TestTypeStr))
            {
                m_TestType = (TestType)i;
                break;
            }
        }
    }
    else
    {
        m_TestType = (TestType)0;
    }

    std::string vCellModel;
    m_session->LoadSolverInfo("CELLMODEL", vCellModel, "FitzHughNagumo");

    ASSERTL0(vCellModel != "", "Cell Model not specified.");

    m_cell = GetCellModelFactory().CreateInstance(vCellModel, m_session,
                                                  m_fields[0]);

    // Stimulus
    m_stimulus = Stimulus::LoadStimuli(m_session, m_fields[0]);

    m_Initx = m_stimulus[0]->ReturnStimuliLoc(0);
    m_Inity = m_stimulus[0]->ReturnStimuliLoc(1);
    m_Initz = m_stimulus[0]->ReturnStimuliLoc(2);

    int nq   = m_fields[0]->GetNpoints();
    int nvar = m_fields.num_elements();

    // Diffusivity coefficient for e^j
    m_epsilon = Array<OneD, NekDouble>(m_mfdim);
    m_session->LoadParameter("epsilon0", m_epsilon[0], 1.0);
    m_session->LoadParameter("epsilon1", m_epsilon[1], 1.0);
    m_session->LoadParameter("epsilon2", m_epsilon[2], 1.0);

    m_session->LoadParameter("AnisotropyRegion", m_AnisotropyRegion, 0);
    m_session->LoadParameter("AnisotropyStrength", m_AnisotropyStrength,
                             m_chi * m_capMembrane);

    // Diffusivity coefficient for u^j
    m_epsu = Array<OneD, NekDouble>(nvar + 1);
    m_session->LoadParameter("epsu0", m_epsu[0], 1.0);
    m_session->LoadParameter("epsu1", m_epsu[1], 1.0);

    m_session->LoadParameter("Diffbeta", m_Diffbeta, 0.5);
    m_session->LoadParameter("Diffeta", m_Diffeta, 100.0);
    m_session->LoadParameter("Diffhe", m_Diffhe, 0.5);

    int shapedim = m_fields[0]->GetShapeDimension();

    // Generating Anisotropy Map
    Array<OneD, Array<OneD, NekDouble>> AniStrength(shapedim);
    for (int j = 0; j < shapedim; ++j)
    {
        AniStrength[j] = Array<OneD, NekDouble>(nq, 1.0);
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

    if (m_MediumType == eAnisotropy)
    {
        m_AnisotropyStrength = m_chi * m_capMembrane;
    }

    else
    {
        m_AnisotropyStrength = 1.0;
    }

    // Or, Get it from Anisotropy Map
    Array<OneD, NekDouble> CardiacFibre;
    if (m_session->DefinesFunction("AnisotropicConductivity"))
    {
        cout << "Loading Anisotropic Fibre map." << endl;

        std::string anisotropy[3] = {"fx", "fy", "fz"};
        CardiacFibre = Array<OneD, NekDouble>(m_spacedim * nq, 0.0);

        Array<OneD, NekDouble> tmp(nq);
        for (int i = 0; i < m_spacedim; ++i)
        {
            GetFunction("AnisotropicConductivity")
                ->Evaluate(anisotropy[i], tmp);
            Vmath::Vcopy(nq, &tmp[0], 1, &CardiacFibre[i * nq], 1);
        }

        // If there is a fibre, let it be with a strength of
        // m_AnisotropyStrength;
        NekDouble mag, fx, fy, fz, Tol = 0.00001;
        for (int k = 0; k < nq; ++k)
        {
            fx = CardiacFibre[k];
            fy = CardiacFibre[k + nq];
            fz = CardiacFibre[k + 2 * nq];

            mag = sqrt(fx * fx + fy * fy + fz * fz);

            if (mag > Tol)
            {
                AniStrength[0][k] = m_AnisotropyStrength;
            }
        }

        PlotCardiacFibre(CardiacFibre);
    }

    // Incorporating Anisotropy into Moving Frames
    if (m_AnisotropyRegion)
    {
        Array<OneD, Array<OneD, int>> EWIndex;
        m_fields[0]->GridIndexElementWise(EWIndex);

        int nptsj = EWIndex[0].num_elements();

        for (int i = 0; i < m_AnisotropyRegion; ++i)
        {
            for (int j = 0; j < nptsj; ++j)
            {
                AniStrength[0][EWIndex[i][j]] = m_AnisotropyStrength;
            }
        }
    }

    if (m_session->DefinesFunction("IsotropicConductivity"))
    {
        cout << "Loading Isotropic Conductivity map." << endl;

        std::string varName;
        varName = "intensity";

        Array<OneD, NekDouble> vTemp;
        GetFunction("IsotropicConductivity")->Evaluate(varName, vTemp);

        // If the d_min and d_max parameters are defined, then we need to
        // rescale the isotropic conductivity to convert from the source
        // domain (e.g. late-gad intensity) to conductivity
        // if (m_session->DefinesParameter("d_min") ||
        //     m_session->DefinesParameter("d_max"))
        // {

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
        //}

        Vmath::Smul(nq, m_AnisotropyStrength, &vTemp[0], 1, &AniStrength[0][0],
                    1);
        Vmath::Vcopy(nq, &vTemp[0], 1, &AniStrength[1][0], 1);

        std::cout << "Conductivity: Max = "
                  << Vmath::Vmax(nq, AniStrength[0], 1)
                  << ", min =" << Vmath::Vmin(nq, AniStrength[0], 1)
                  << std::endl;
    }

    if (m_session->DefinesFunction("AnisotropicConductivity"))
    {
        MMFSystem::MMFInitObject(AniStrength, CardiacFibre);
    }

    else
    {
        MMFSystem::MMFInitObject(AniStrength);
    }

    // Plot HHD
    if (m_session->DefinesSolverInfo("GenerateHHDPlot"))
    {
        std::string PlotHHDStr;
        PlotHHDStr = m_session->GetSolverInfo("GenerateHHDPlot");
        if (PlotHHDMap[0] == PlotHHDStr)
        {
            std::cout
                << "PlotHHD initated ======================================="
                << std::endl;
            int nstep;
            m_session->LoadParameter("HHDPlotnstep", nstep, 0);
            GenerateHHDPlot(nstep);
        }
    }

    if (m_session->DefinesSolverInfo("INITWAVETYPE"))
    {
        std::string InitWaveTypeStr;
        InitWaveTypeStr = m_session->GetSolverInfo("INITWAVETYPE");
        for (int i = 0; i < (int)SIZE_TestType; ++i)
        {
            if (boost::iequals(InitWaveTypeMap[i], InitWaveTypeStr))
            {
                m_InitWaveType = (InitWaveType)i;
                break;
            }
        }
    }
    else
    {
        m_InitWaveType = (InitWaveType)0;
    }

    m_ode.DefineOdeRhsMMF(&MMFBidomain::DoOdeRhsBidomain, this);

    if (m_explicitDiffusion)
    {
        m_ode.DefineImplicitSolveMMF(&MMFBidomain::DoNullSolveMMF, this);
        m_ode.DefineProjection(&MMFBidomain::DoOdeProjection, this);
    }

    else
    {
        // Create varcoeff for Helmsolver
        ComputeVarCoeff2D(m_movingframes, m_varcoeff);

        switch (m_TestType)
        {
            case eBidomainCardiac:
            case eBidomainBrain:
            {
                m_ode.DefineImplicitSolveMMF(
                    &MMFBidomain::DoImplicitSolveBidomain, this);
            }
            break;

            default:
                break;
        }
    }
}

/**
 *
 */
MMFBidomain::~MMFBidomain()
{
}

void MMFBidomain::DoNullSolveMMF(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    const Array<OneD, const Array<OneD, NekDouble>> &MF1st,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time,
    const NekDouble lambda)
{
    int nvariables = inarray.num_elements();
    int nq         = m_fields[0]->GetNpoints();

    if (time > 0)
    {
    }

    if (lambda > 0)
    {
    }

    if (RootMeanSquare(MF1st[0]) > 0)
    {
    }

    for (int i = 0; i < nvariables; ++i)
    {
        Vmath::Vcopy(nq, &inarray[i][0], 1, &outarray[i][0], 1);
    }
}

void MMFBidomain::DoOdeProjection(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time)
{
    // Counter variable
    int i;
    int npoints    = GetNpoints();
    int nVariables = inarray.num_elements();

    // Set the boundary conditions
    SetBoundaryConditions(time);

    // Switch on the projection type (Discontinuous or Continuous)
    for (i = 0; i < nVariables; ++i)
    {
        Vmath::Vcopy(npoints, inarray[i], 1, outarray[i], 1);
    }
}

void MMFBidomain::v_DoSolve()
{
    ASSERTL0(m_intScheme != 0, "No time integration scheme.");

    int i, nchk = 1;
    int nq         = GetTotPoints();
    int nvariables = 0;
    int nfields    = m_fields.num_elements();

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
    m_intSoln =
        m_intScheme->InitializeScheme(m_timestep, fields, m_time, m_ode);

    // Check uniqueness of checkpoint output
    ASSERTL0((m_checktime == 0.0 && m_checksteps == 0) ||
                 (m_checktime > 0.0 && m_checksteps == 0) ||
                 (m_checktime == 0.0 && m_checksteps > 0),
             "Only one of IO_CheckTime and IO_CheckSteps "
             "should be set!");

    Array<OneD, Array<OneD, NekDouble>> MF1st(m_spacedim);
    for (int i = 0; i < m_spacedim; ++i)
    {
        MF1st[i] = Array<OneD, NekDouble>(3 * nq);
        Vmath::Smul(3 * nq, 1.0, &m_movingframes[i][0], 1, &MF1st[i][0], 1);
    }

    LibUtilities::Timer timer;
    bool doCheckTime  = false;
    int step          = 0;
    NekDouble intTime = 0.0;
    NekDouble cpuTime = 0.0;
    NekDouble elapsed = 0.0;

    while (step < m_steps || m_time < m_fintime - NekConstants::kNekZeroTol)
    {
        timer.Start();
        fields = m_intScheme->TimeIntegrateMMF(step, m_timestep, m_intSoln,
                                               MF1st, m_ode);
        timer.Stop();

        m_time += m_timestep;
        elapsed = timer.TimePerTest(1);
        intTime += elapsed;
        cpuTime += elapsed;

        if (m_session->GetComm()->GetRank() == 0 && !((step + 1) % m_infosteps))
        {
            std::cout << "max u = "
                      << Vmath::Vamax(nq, m_fields[0]->GetPhys(), 1)
                      << ", max ue = "
                      << Vmath::Vamax(nq, m_fields[1]->GetPhys(), 1)
                      << std::endl;

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

        // Transform data into coefficient space: Only for u. ue is updated at
        // RHS
        m_fields[m_intVariables[0]]->SetPhys(fields[0]);
        m_fields[m_intVariables[0]]->FwdTrans_IterPerExp(
            fields[0], m_fields[m_intVariables[0]]->UpdateCoeffs());
        m_fields[m_intVariables[0]]->SetPhysState(false);

        // Write out checkpoint files
        if ((m_checksteps && step && !((step + 1) % m_checksteps)) ||
            doCheckTime)
        {
            PlotBiDomain(nchk);

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

void MMFBidomain::DoImplicitSolveBidomain(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    const Array<OneD, const Array<OneD, NekDouble>> &MF1st,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time,
    const NekDouble lambda)
{
    int nq = m_fields[0]->GetNpoints();

    // Set up factors for Helmsolver
    StdRegions::ConstFactorMap factors;
    factors[StdRegions::eFactorTau]    = m_Helmtau;
    factors[StdRegions::eFactorLambda] = 1.0 / lambda;

    if (RootMeanSquare(MF1st[0]) > 0)
    {
    }

    if (time > 0)
    {
    }

    // We solve ( \nabla^2 - HHlambda ) Y[i] = rhs [i]
    // inarray = input: \hat{rhs} -> output: \hat{Y}
    // outarray = output: \hat{Y} where \hat = modal coeffs

    // Multiply 1.0/timestep
    // Vmath::Smul(nq, -factors[StdRegions::eFactorLambda], inarray[0], 1,
    //             m_fields[0]->UpdatePhys(), 1);

    // m_fields[0]->HelmSolve(m_fields[0]->GetPhys(),m_fields[0]->UpdateCoeffs(),
    //                        NullFlagList, factors, m_varcoeff);

    // m_fields[0]->BwdTrans(m_fields[0]->GetCoeffs(),
    // m_fields[0]->UpdatePhys()); m_fields[0]->SetPhysState(true);

    // outarray[0] = m_fields[0]->GetPhys();

    // No diffusion for the second variable
    Vmath::Vcopy(nq, &inarray[0][0], 1, &outarray[0][0], 1);
    Vmath::Vcopy(nq, &inarray[1][0], 1, &outarray[1][0], 1);
}

void MMFBidomain::DoOdeRhsBidomain(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    const Array<OneD, const Array<OneD, NekDouble>> &MF1st,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time)
{
    int nq      = m_fields[0]->GetNpoints();
    int ncoeffs = m_fields[0]->GetNcoeffs();

    if (RootMeanSquare(MF1st[0]) > 0)
    {
    }

    // Compute the reaction function
    // input: inarray: inarray[0] and inarray[1]
    // output: outarray: only outarray[0]
    m_cell->TimeIntegrate(inarray, outarray, time);

    // Compute I_stim
    for (unsigned int i = 0; i < m_stimulus.size(); ++i)
    {
        m_stimulus[i]->Update(outarray, time);
    }

    // Compute - nabla \cdot \nabla phi_m
    Array<OneD, NekDouble> tmp(nq);
    tmp = ComputeCovariantDiffusion(m_movingframes, inarray[0], eEuclidean);
    NekDouble DivAvg = -1.0 * AvgInt(tmp);
    Vmath::Sadd(nq, DivAvg, tmp, 1, tmp, 1);
    Vmath::Neg(nq, tmp, 1);

    Vmath::Vcopy(nq, tmp, 1, m_fields[1]->UpdatePhys(), 1);

    // Find \phi_e satisfying \nabla \cdot (\sigma_e + \sigma_i) \nabla \phi_e =
    // - nabla \cdot \sigma_i \nabla phi_m
    StdRegions::ConstFactorMap factorsPoisson;
    factorsPoisson[StdRegions::eFactorTau]    = m_Helmtau;
    factorsPoisson[StdRegions::eFactorLambda] = 0.0;

    m_fields[1]->HelmSolve(m_fields[1]->GetPhys(), m_fields[1]->UpdateCoeffs(),
                           NullFlagList, factorsPoisson, m_varcoeff);
    m_fields[1]->BwdTrans(m_fields[1]->GetCoeffs(), m_fields[1]->UpdatePhys());
    m_fields[1]->SetPhysState(true);

    // tmp = ComputeCovariantDiffusion(m_movingframes, m_fields[1]->GetPhys(),
    // eEuclidean); Vmath::Neg(nq, tmp, 1);

    Array<OneD, NekDouble> tmpc(ncoeffs);
    WeakDGMMFDiffusion(0, m_fields[0]->GetPhys(), tmp);
    Vmath::Neg(nq, tmp, 1);

    Vmath::Vadd(nq, &tmp[0], 1, &outarray[0][0], 1, &outarray[0][0], 1);

    // Array<OneD, NekDouble> ggrad0(nq), ggrad1(nq), ggrad2(nq), ggrad(nq);
    // StdRegions::ConstFactorMap factorsPoisson;
    // factorsPoisson[StdRegions::eFactorTau] = m_Helmtau;
    // factorsPoisson[StdRegions::eFactorLambda] = 0.0;

    // // ----------------------------
    // // Compute \nabla g_i \nabla Vm
    // // ----------------------------
    // m_fields[0]->PhysDeriv(inarray[0], ggrad0, ggrad1, ggrad2);
    // m_fields[0]->PhysDeriv(0, ggrad0, ggrad0);
    // m_fields[0]->PhysDeriv(1, ggrad1, ggrad1);
    // m_fields[0]->PhysDeriv(2, ggrad2, ggrad2);
    // // if (m_session->DefinesFunction("IntracellularAnisotropicConductivity")
    // &&
    // // m_session->DefinesFunction("ExtracellularAnisotropicConductivity"))
    // // {
    // //     Vmath::Vmul(nq, m_vardiffi[StdRegions::eVarCoeffD00], 1, ggrad0,
    // 1,
    // //                 ggrad0, 1);
    // //     Vmath::Vmul(nq, m_vardiffi[StdRegions::eVarCoeffD11], 1, ggrad1,
    // 1,
    // //                 ggrad1, 1);
    // //     Vmath::Vmul(nq, m_vardiffi[StdRegions::eVarCoeffD22], 1, ggrad2,
    // 1,
    // //                 ggrad2, 1);
    // // }
    // // Add partial derivatives together
    // Vmath::Vadd(nq, ggrad0, 1, ggrad1, 1, ggrad, 1);
    // Vmath::Vadd(nq, ggrad2, 1, ggrad, 1, ggrad, 1);

    // Vmath::Smul(nq, -1.0, ggrad, 1, m_fields[1]->UpdatePhys(), 1);

    // std::cout << "time = " << time;
    // std::cout << ", Laplacian phi_m = " <<
    // RootMeanSquare(m_fields[1]->GetPhys());

    // // ----------------------------
    // // Solve Poisson problem for Ve
    // // ----------------------------
    // m_fields[1]->HelmSolve(m_fields[1]->GetPhys(),
    // m_fields[1]->UpdateCoeffs(),
    //                        NullFlagList, factorsPoisson, m_varcoeff);
    // m_fields[1]->BwdTrans(m_fields[1]->GetCoeffs(),
    // m_fields[1]->UpdatePhys()); m_fields[1]->SetPhysState(true);

    // std::cout << ", phi from Poisson = " <<
    // RootMeanSquare(m_fields[1]->GetPhys());

    // // ------------------------------
    // // Compute Laplacian of Ve (forcing term)
    // // ------------------------------
    // m_fields[1]->PhysDeriv(m_fields[1]->GetPhys(), ggrad0, ggrad1, ggrad2);
    // m_fields[1]->PhysDeriv(0, ggrad0, ggrad0);
    // m_fields[1]->PhysDeriv(1, ggrad1, ggrad1);
    // m_fields[1]->PhysDeriv(2, ggrad2, ggrad2);
    // // if (m_session->DefinesFunction("IntracellularAnisotropicConductivity")
    // &&
    // // m_session->DefinesFunction("ExtracellularAnisotropicConductivity"))
    // // {
    // //     Vmath::Vmul(nq, m_vardiffi[StdRegions::eVarCoeffD00], 1, ggrad0,
    // 1,
    // //                 ggrad0, 1);
    // //     Vmath::Vmul(nq, m_vardiffi[StdRegions::eVarCoeffD11], 1, ggrad1,
    // 1,
    // //                 ggrad1, 1);
    // //     Vmath::Vmul(nq, m_vardiffi[StdRegions::eVarCoeffD22], 1, ggrad2,
    // 1,
    // //                 ggrad2, 1);
    // // }

    // // Add partial derivatives together
    // Vmath::Vadd(nq, ggrad0, 1, ggrad1, 1, ggrad, 1);
    // Vmath::Vadd(nq, ggrad2, 1, ggrad, 1, ggrad, 1);

    // Vmath::Vadd(nq, ggrad, 1, outarray[0], 1, outarray[0], 1);

    // std::cout << ", Laplacian phi_e = " << RootMeanSquare(outarray[0]) <<
    // std::endl;
}

void MMFBidomain::v_SetInitialConditions(NekDouble initialtime,
                                         bool dumpInitialConditions,
                                         const int domain)
{
    int nq = GetTotPoints();

    if (domain == 0)
    {
    }

    if (dumpInitialConditions)
    {
    }

    m_cell->Initialise();

    // Read initial condition from xml file
    EquationSystem::v_SetInitialConditions(initialtime, false);

    Array<OneD, Array<OneD, NekDouble>> tmp(1);
    tmp[0] = Array<OneD, NekDouble>(nq);
    Vmath::Vcopy(nq, m_fields[0]->GetPhys(), 1, tmp[0], 1);
    for (unsigned int i = 0; i < m_stimulus.size(); ++i)
    {
        m_stimulus[i]->Update(tmp, initialtime);
        m_fields[0]->SetPhys(tmp[0]);
    }

    // forward transform to fill the modal coeffs
    for (int i = 0; i < m_fields.num_elements(); ++i)
    {
        m_fields[i]->SetPhysState(true);
        m_fields[i]->FwdTrans(m_fields[i]->GetPhys(),
                              m_fields[i]->UpdateCoeffs());
    }

    std::cout << "Initial: max u = "
              << Vmath::Vmax(nq, m_fields[0]->GetPhys(), 1) << std::endl;

    if (dumpInitialConditions)
    {
        std::string outname;
        outname = m_sessionName + "_initial.chk";

        WriteFld(outname);
    }
}

void MMFBidomain::v_GenerateSummary(SolverUtils::SummaryList &s)
{
    MMFSystem::v_GenerateSummary(s);
    SolverUtils::AddSummaryItem(s, "TestType", TestTypeMap[m_TestType]);
    SolverUtils::AddSummaryItem(s, "MediumType", MediumTypeMap[m_MediumType]);
    if (m_AnisotropyRegion)
    {
        SolverUtils::AddSummaryItem(s, "AnisotropyRegion", m_AnisotropyRegion);
    }

    SolverUtils::AddSummaryItem(s, "Helmtau", m_Helmtau);

    m_cell->GenerateSummary(s);
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
