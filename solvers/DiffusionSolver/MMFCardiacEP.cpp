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
    int nvar = m_fields.size();

    // Conductance parameters
    m_session->LoadParameter("Chi", m_chi, 28.0);
    m_session->LoadParameter("Cm", m_capMembrane, 0.125);

    // Helmsolver parameter
    m_session->LoadParameter("Helmtau", m_Helmtau, 1.0);

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

    // TimeMap ?
    if (m_session->DefinesSolverInfo("TimeMapType"))
    {
        std::string TIMEMAPTYPEStr;
        TIMEMAPTYPEStr = m_session->GetSolverInfo("TimeMapType");
        for (int i = 0; i < (int)SIZE_TimeMapType; ++i)
        {
            if (boost::iequals(TimeMapTypeMap[i], TIMEMAPTYPEStr))
            {
                m_TimeMap = (TimeMapType)i;
                break;
            }
        }
    }
    else
    {
        m_TimeMap = (TimeMapType)0;
    }

    std::string vCellModel;
    m_session->LoadSolverInfo("CELLMODEL", vCellModel, "FitzHughNagumo");

    ASSERTL0(vCellModel != "", "Cell Model not specified.");

    m_cell = GetCellModelFactory().CreateInstance(vCellModel, m_session,
                                                  m_fields[0]);

    // Stimulus
    m_stimulus = Stimulus::LoadStimuli(m_session, m_fields[0]);

    // m_Initx = m_stimulus[0]->ReturnStimuliLoc(0);
    // m_Inity = m_stimulus[0]->ReturnStimuliLoc(1);
    // m_Initz = m_stimulus[0]->ReturnStimuliLoc(2);

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

    // Derive AnisotropyStrength.
    Array<OneD, Array<OneD, NekDouble>> AniStrength(m_expdim);
    for (int j = 0; j < m_expdim; ++j)
    {
        AniStrength[j] = Array<OneD, NekDouble>(nq, 1.0);
    }

    // Ratio between along the fiber and orthogonal to the fiber
    switch (m_MediumType)
    {
        case eAnisotropy:
        case eHeterogeneousAnisotropy:
        {
            Array<OneD, NekDouble> CardiacFibre;
            LoadCardiacFiber(m_MediumType, m_AnisotropyStrength, AniStrength,
                             CardiacFibre);
            MMFSystem::MMFInitObject(AniStrength, CardiacFibre);
        }
        break;

        case eHeterogeneousIsotropy:
        case eRegionalHeterogeneous:
        {
            LoadCardiacFiber(m_MediumType, m_AnisotropyStrength, AniStrength);
            MMFSystem::MMFInitObject(AniStrength);
        }
        break;

        case eIsotropy:
        default:
        {
            Array<OneD, NekDouble> Anitmp(nq, 1.0);

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

            Vmath::Vcopy(nq, &Anitmp[0], 1, &AniStrength[0][0], 1);
            Vmath::Vcopy(nq, &Anitmp[0], 1, &AniStrength[1][0], 1);

            MMFSystem::MMFInitObject(AniStrength);
        }
        break;
    }

    // plot Conductivity map
    int ncoeffs = m_fields[0]->GetNcoeffs();

    std::string outname1;
    outname1 = m_sessionName + "_sigmaMap.chk";

    std::vector<Array<OneD, NekDouble>> fieldcoeffs(2);
    for (int i = 0; i < 2; ++i)
    {
        fieldcoeffs[i] = Array<OneD, NekDouble>(ncoeffs);
    }

    std::vector<std::string> variables(nvar);
    variables[0] = "Condx";
    variables[1] = "Condy";

    m_fields[0]->FwdTrans(AniStrength[0], fieldcoeffs[0]);
    m_fields[0]->FwdTrans(AniStrength[1], fieldcoeffs[1]);

    WriteFld(outname1, m_fields[0], fieldcoeffs, variables);

    // Plot HHD
    // if (m_session->DefinesSolverInfo("GenerateHHDPlot"))
    // {
    //     std::string PlotHHDStr;
    //     PlotHHDStr = m_session->GetSolverInfo("GenerateHHDPlot");
    //     if (PlotHHDMap[0] == PlotHHDStr)
    //     {
    //         std::cout << "PlotHHD initated
    //         ======================================="
    //                   << std::endl;
    //         int nstep;
    //         m_session->LoadParameter("HHDPlotnstep", nstep, 0);
    //         GenerateHHDPlot(nstep);
    //     }
    // }

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

    m_ode.DefineOdeRhs(&MMFCardiacEP::DoOdeRhsCardiacEP, this);
}

/**
 *
 */
MMFCardiacEP::~MMFCardiacEP()
{
}

void MMFCardiacEP::LoadCardiacFiber(
    const SolverUtils::MediumType CardiacMediumType,
    const NekDouble AnisotropyStrength,
    Array<OneD, Array<OneD, NekDouble>> &AniStrength,
    Array<OneD, NekDouble> &CardiacFibre)
{
    int nq = m_fields[0]->GetNpoints();

    switch (CardiacMediumType)
    {
        case eAnisotropy:
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
            for (int i = 0; i < m_AnisotropyRegion; ++i)
            {
                for (int j = 0; j < m_fields[0]->GetTotPoints(i); ++j)
                {
                    index                 = m_fields[0]->GetPhys_Offset(i) + j;
                    AniStrength[0][index] = AnisotropyStrength;
                }
            }
        }
        break;

        default:
            break;
    }

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
        case eMMFFirst:
        {
            DoSolveMMFFirst();
        }
        break;

        default:
        {
            DoSolveMMF();
        }
        break;
    }
}

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
    Array<OneD, NekDouble> UnitVelMap(m_spacedim * nq, 0.0);

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
            ComputeTimeMap(m_time, fields[0], dudtval, m_ValidTimeMap,
                           dudtvalHistory, IappMap, TimeMap);
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
            // NekDouble dudtpros, dudtneg;
            // dudtpros = Computedudtpercent(1, dudt);
            // dudtneg  = Computedudtpercent(-1, dudt);

            // NekDouble udiff;
            // udiff = Vmath::Vmax(nq, fields[0], 1) - Vmath::Vmin(nq,
            // fields[0], 1);

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

            // Array<OneD, Array<OneD, NekDouble>> stimulusstrength(nvariables);
            // for (unsigned int i = 0; i < m_stimulus.size(); ++i)
            // {
            //     for (int j=0; j<nvariables; ++j)
            //     {
            //         stimulusstrength[j] = Array<OneD, NekDouble>(nq, 0.0);
            //     }

            //     m_stimulus[0]->Update(stimulusstrength, m_time);

            //     if(Vmath::Vmax(nq, stimulusstrength[0], 1)>0.01)
            //     {
            //         std::cout << " =================================== " <<
            //         std::endl; std::cout << "i = " << i << ", Stimulus: = "
            //         << Vmath::Vmax(nq, stimulusstrength[0], 1) << std::endl;
            //     std::cout << " =================================== " <<
            //     std::endl;
            //     }
            // }

            if (CountActivated(ActivatedHistory) > 0)
            {
                ComputeRelacc(MF1stAligned, Relacc);

                PlotTrajectoryMF(ActivatedHistory, fields[0], MF1stAligned,
                                 MF1stConnection, Relacc, NoBoundaryZone, nchk);

                // Plot relative acceleration and conduction block zone
                // PlotRelacc2D(ActivatedHistory, MF1stConnection,
                // MF1stCurvature, Relacc, RelaccOmega, nchk);
            }

            if (m_TimeMap == eActivated)
            {
                // Compute velocity field
                // VelocityMap = ComputeVelocityField(m_ValidTimeMap, TimeMap);
                ComputeMFTimeMap(m_ValidTimeMap, TimeMap, NewValidTimeMap,
                                 TimeMapMF);

                // Compute2DConnectionCurvature(TimeMapMF, TMMFConnection,
                // TMMFCurvature);

                Compute2DConnection1form(TimeMapMF, TMMFConnection);

                // Compute Relative Acceleration
                ComputeRelacc(TimeMapMF, TMRelacc);

                PlotTimeMapMF(NoBoundaryZone, TimeMap, TimeMapMF, MF1stAligned,
                              TMMFConnection, TMRelacc, nchk);

                std::cout << "Time Map: Max = " << Vmath::Vmax(nq, TimeMap, 1)
                          << ", Min = " << Vmath::Vmin(nq, TimeMap, 1)
                          << std::endl;
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

    Array<OneD, NekDouble> tmpc(ncoeffs);

    Array<OneD, NekDouble> velmag(nq, 0.0);
    Array<OneD, NekDouble> velocity(m_spacedim * nq);

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
        // field time integration
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
            int Iumax = Vmath::Iamax(nq, fields[0], 1);
            std::cout << "u_max = " << Vmath::Vamax(nq, fields[0], 1)
                      << " at x = " << x0[Iumax] << ", y = " << x1[Iumax] << ", z = " << x2[Iumax]
                      << std::endl;

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

    int nq = GetTotPoints();

    m_cell->Initialise();

    std::cout << "SetInitialCond 1" << std::endl;

    // Read initial condition from xml file
    EquationSystem::v_SetInitialConditions(initialtime, false);

    std::cout << "SetInitialCond 2" << std::endl;


    Array<OneD, Array<OneD, NekDouble>> tmp(1);
    tmp[0] = Array<OneD, NekDouble>(nq);

    Array<OneD, NekDouble> initialcondition(nq);
    Vmath::Vcopy(nq, m_fields[0]->GetPhys(), 1, tmp[0], 1);
    Vmath::Vcopy(nq, tmp[0], 1, initialcondition, 1);
    for (unsigned int i = 0; i < m_stimulus.size(); ++i)
    {
        m_stimulus[i]->Update(tmp, initialtime);
        m_fields[0]->SetPhys(tmp[0]);
    }

    std::cout << "SetInitialCond 3" << std::endl;


    // Only the excited regions are considered for m_InitExcitation
    Vmath::Vsub(nq, tmp[0], 1, initialcondition, 1, initialcondition, 1);
    // m_ValidTimeMap = ComputeTimeMapInitialZone(initialcondition);

    if (m_TimeMap == eProcessing)
    {
        TimeMapProcess();
    }
    std::cout << "SetInitialCond 4" << std::endl;

    // forward transform to fill the modal coeffs
    for (int i = 0; i < m_fields.size(); ++i)
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

// Array<OneD, NekDouble> MMFCardiacEP::PlanePhiWave()
// {
//     int nq = GetTotPoints();
//     Array<OneD, NekDouble> outarray(nq, 0.0);

//     Array<OneD, NekDouble> x(nq);
//     Array<OneD, NekDouble> y(nq);
//     Array<OneD, NekDouble> z(nq);

//     m_fields[0]->GetCoords(x, y, z);

//     NekDouble xmin, ymin, xmax;

//     xmin = Vmath::Vmin(nq, x, 1);
//     xmax = Vmath::Vmax(nq, x, 1);
//     ymin = Vmath::Vmin(nq, y, 1);

//     NekDouble xp, yp, xp2;
//     for (int i = 0; i < nq; i++)
//     {
//         switch (m_InitWaveType)
//         {
//             case eLeft:
//             {
//                 NekDouble radiusofinit;
//                 NekDouble frontstiff;

//                 radiusofinit = 3.0;
//                 frontstiff   = 0.5;

//                 xp = x[i] - xmin;
//                 outarray[i] =
//                     1.0 / (1.0 + exp((xp - radiusofinit) / frontstiff));
//             }
//             break;

//             case eBothEnds:
//             {
//                 NekDouble radiusofinit = 3.0;
//                 NekDouble frontstiff   = 0.1;

//                 xp  = x[i] - xmin;
//                 xp2 = x[i] - xmax;

//                 outarray[i] =
//                     1.0 / (1.0 +
//                            exp((sqrt(xp * xp) - radiusofinit) / frontstiff))
//                            +
//                     1.0 / (1.0 +
//                            exp((sqrt(xp2 * xp2) - radiusofinit) /
//                            frontstiff));
//             }
//             break;

//             case eCenter:
//             {
//                 NekDouble radiusofinit = 6.0;
//                 NekDouble frontstiff   = 0.1;

//                 // NekDouble xc = 0.5*(Vmath::Vmax(nq, x, 1) +
//                 // Vmath::Vmin(nq, x, 1));

//                 xp = x[i] - xmin;
//                 outarray[i] =
//                     1.0 / (1.0 + exp((xp - radiusofinit) / frontstiff));
//             }
//             break;

//             case eLeftBottomCorner:
//             {
//                 NekDouble radiusofinit = 6.0;
//                 NekDouble frontstiff   = 0.1;
//                 NekDouble bs           = 2.0;

//                 xp = x[i] - xmin;
//                 yp = y[i] - ymin;
//                 outarray[i] =
//                     1.0 /
//                     (1.0 + exp((sqrt(xp * xp + yp * yp) / bs - radiusofinit)
//                     /
//                                frontstiff));
//             }
//             break;

//             case ePoint:
//             {
//                 NekDouble xloc, yloc, zloc, rad;
//                 NekDouble radiusofinit = 5.0;

//                 xloc = x[i] - m_Initx;
//                 yloc = y[i] - m_Inity;
//                 zloc = z[i] - m_Initz;

//                 rad = sqrt(xloc * xloc + yloc * yloc + zloc * zloc);

//                 xloc = xloc / radiusofinit;
//                 yloc = yloc / radiusofinit;
//                 zloc = zloc / radiusofinit;

//                 if (rad < radiusofinit)
//                 {
//                     outarray[i] =
//                         exp(-(1.0 / 2.0) *
//                             (xloc * xloc + yloc * yloc + zloc * zloc));
//                 }

//                 else
//                 {
//                     outarray[i] = 0.0;
//                 }
//             }
//             break;

//             case eSpiralDock:
//             {
//                 NekDouble radiusofinit = 3.0;
//                 NekDouble frontstiff   = 0.1;
//                 xp                     = x[i] - 4.0;
//                 yp                     = y[i];
//                 outarray[i] =
//                     (1.0 / (1.0 + exp(2.0 * yp))) *
//                     (1.0 / (1.0 + exp(-2.0 * xp))) *
//                     (1.0 / (1.0 + exp((xp - radiusofinit) / frontstiff)));
//             }
//             break;

//             default:
//                 break;
//         } // namespace Nektar
//     }

//     return outarray;
// }

void MMFCardiacEP::v_EvaluateExactSolution(unsigned int field,
                                           Array<OneD, NekDouble> &outfield,
                                           const NekDouble time)
{
    EquationSystem::v_EvaluateExactSolution(field, outfield, time);
}

// Compute \int \nabla u \cdot e^{dir}
// void MMFCardiacEP::WeakDGDirectionalDeriv(
//     const int direction, const Array<OneD, const Array<OneD, NekDouble>>
//     &MF1st, const Array<OneD, const NekDouble> &InField, Array<OneD,
//     NekDouble> &OutField)
// {
//     int ncoeffs         = GetNcoeffs();
//     int nTracePointsTot = GetTraceNpoints();
//     int nq              = GetNpoints();

//     Array<OneD, NekDouble> physfield(nq);

//     // Get the variables in physical space
//     // already in physical space
//     Vmath::Vcopy(nq, InField, 1, physfield, 1);

//     Array<OneD, NekDouble> WeakDeriv(ncoeffs, 0.0);
//     Array<OneD, NekDouble> tmp(nq);

//     // Directional derivation with respect to the j'th moving frame
//     // tmp[j] = \nabla \physfield[i] \cdot \mathbf{e}^j
//     // Implemented at TriExp::v_IProductWRTDirectionalDerivBase_SumFa
//     m_fields[0]->IProductWRTDirectionalDerivBase(MF1st[direction],
//     physfield,
//                                                  WeakDeriv);

//     // if the NumericalFluxs function already includes the normal in the
//     output Array<OneD, NekDouble> Fwd(nTracePointsTot); Array<OneD,
//     NekDouble> Bwd(nTracePointsTot);

//     Array<OneD, NekDouble> flux(nTracePointsTot, 0.0);
//     Array<OneD, NekDouble> fluxFwd(nTracePointsTot);
//     Array<OneD, NekDouble> fluxBwd(nTracePointsTot);

//     // Evaluate numerical flux in physical space which may in
//     // general couple all component of vectors
//     m_fields[0]->GetFwdBwdTracePhys(physfield, Fwd, Bwd);

//     // evaulate upwinded m_fields[i]
//     Array<OneD, NekDouble> traceVn(nTracePointsTot, 0.0);
//     Array<OneD, NekDouble> tmptrace(nTracePointsTot);
//     for (int i = 0; i < m_spacedim; ++i)
//     {
//         Vmath::Vcopy(nq, &MF1st[direction][i * nq], 1, &tmp[0], 1);
//         m_fields[0]->ExtractTracePhys(tmp, tmptrace);
//         Vmath::Vvtvp(nTracePointsTot, m_traceNormals[i], 1, tmptrace, 1,
//                      traceVn, 1, traceVn, 1);
//     }

//     m_fields[0]->GetTrace()->Upwind(traceVn, Fwd, Bwd, flux);

//     Array<OneD, Array<OneD, NekDouble>> ncdotMFFwd;
//     Array<OneD, Array<OneD, NekDouble>> ncdotMFBwd;

//     ComputencdotMF(MF1st, ncdotMFFwd, ncdotMFBwd);

//     OutField = Array<OneD, NekDouble>(ncoeffs, 0.0);
//     // calculate numflux = (n \cdot MF)*flux
//     Vmath::Vmul(nTracePointsTot, &flux[0], 1, &ncdotMFFwd[direction][0],
//     1,
//                 &fluxFwd[0], 1);
//     Vmath::Vmul(nTracePointsTot, &flux[0], 1, &ncdotMFBwd[direction][0],
//     1,
//                 &fluxBwd[0], 1);

//     // FwdBwdtegral because generallize (N \cdot MF)_{FWD} \neq -(N \cdot
//     // MF)_{BWD}
//     Vmath::Neg(ncoeffs, WeakDeriv, 1);
//     m_fields[0]->AddFwdBwdTraceIntegral(fluxFwd, fluxBwd, WeakDeriv);
//     m_fields[0]->SetPhysState(false);

//     Vmath::Vadd(ncoeffs, &WeakDeriv[0], 1, &OutField[0], 1, &OutField[0],
//     1);
// }

void MMFCardiacEP::v_GenerateSummary(SolverUtils::SummaryList &s)
{
    MMFSystem::v_GenerateSummary(s);
    AddSummaryItem(s, "SolverSchemeType", SolverSchemeTypeMap[m_SolverSchemeType]);

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
