///////////////////////////////////////////////////////////////////////////////
//
// File: MMFDiffusion.cpp
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
// Description: MMFDiffusion.
//
///////////////////////////////////////////////////////////////////////////////

#include <iomanip>
#include <iostream>

#include <boost/algorithm/string.hpp>
#include <boost/core/ignore_unused.hpp>

#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>
#include <DiffusionSolver/EquationSystems/MMFDiffusion.h>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>
#include <SolverUtils/Driver.h>

#include <boost/math/special_functions/spherical_harmonic.hpp>
using namespace std;
using namespace Nektar::SolverUtils;
using namespace Nektar;

namespace Nektar
{
string MMFDiffusion::className =
    SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
        "MMFDiffusion", MMFDiffusion::create, "MMFDiffusion equation.");

MMFDiffusion::MMFDiffusion(const LibUtilities::SessionReaderSharedPtr &pSession,
                           const SpatialDomains::MeshGraphSharedPtr &pGraph)
    : UnsteadySystem(pSession, pGraph), MMFSystem(pSession, pGraph)
{
}

void MMFDiffusion::v_InitObject(bool DeclareFields)
{
    UnsteadySystem::v_InitObject(DeclareFields);
    
    int nq    = m_fields[0]->GetNpoints();
    int nvar  = m_fields.size();

    // AniStrength for e^1 and e^2
    m_session->LoadParameter("AniStrength", m_AniStrength, 1.0);
    m_session->LoadParameter("Helmtau", m_Helmtau, 1.0);
    m_session->LoadParameter("EmbededPlane", m_EmbededPlane, 0);

    // Diffusivity coefficient for e^j
    m_epsilon = Array<OneD, NekDouble>(m_spacedim);
    m_session->LoadParameter("epsilon0", m_epsilon[0], 1.0);
    m_session->LoadParameter("epsilon1", m_epsilon[1], 1.0);
    m_session->LoadParameter("epsilon2", m_epsilon[2], 1.0);

    Array<OneD, NekDouble> x0(nq);
    Array<OneD, NekDouble> x1(nq);
    Array<OneD, NekDouble> x2(nq);

    m_fields[0]->GetCoords(x0, x1, x2);

    // m_epsvec = Array<OneD, NekDouble>(nq);
    // for (int i=0; i<nq; ++i)
    // {
    //     if(x1[i]<0)
    //     {
    //         m_epsvec[i] = 1.0;
    //     }

    //     else
    //     {
    //         m_epsvec[i] = 4.0;
    //     }
    // }

    // std::cout << "m_epsvec = " << RootMeanSquare(m_epsvec) << std::endl;

    m_session->LoadParameter("d00", m_d00, 1.0);
    m_session->LoadParameter("d11", m_d11, 1.0);
    m_session->LoadParameter("d22", m_d22, 1.0);

    // Diffusivity coefficient for u^j
    m_epsu = Array<OneD, NekDouble>(nvar + 1);
    m_session->LoadParameter("epsu0", m_epsu[0], 1.0);
    m_session->LoadParameter("epsu1", m_epsu[1], 1.0);

    m_session->LoadParameter("frequency", m_frequency, m_pi);

    m_session->LoadParameter("InitPtx", m_InitPtx, 0.0);
    m_session->LoadParameter("InitPty", m_InitPty, 0.0);
    m_session->LoadParameter("InitPtz", m_InitPtz, 0.0);

    int shapedim = m_fields[0]->GetShapeDimension();
    Array<OneD, Array<OneD, NekDouble>> Anisotropy(shapedim);
    for (int j = 0; j < shapedim; ++j)
    {
        Anisotropy[j] = Array<OneD, NekDouble>(nq, 1.0);
    }

    if (m_session->DefinesParameter("d00"))
    {
        Array<OneD, NekDouble> x0(nq);
        Array<OneD, NekDouble> x1(nq);
        Array<OneD, NekDouble> x2(nq);

        m_fields[0]->GetCoords(x0, x1, x2);

        m_varcoeffXYZ[StdRegions::eVarCoeffD00] = Array<OneD, NekDouble>(nq);
        m_d00vec = Array<OneD, NekDouble>(nq);

        int index;

            for (int i=0; i<nq; ++i)
            {
                m_d00vec[i] = m_d00;
                // m_d00vec[i] = sqrt(m_d00) * ( 2.0 + sin(m_frequency * x0[i]) );

                Anisotropy[0][i] = sqrt(m_d00vec[i]);
                m_varcoeffXYZ[StdRegions::eVarCoeffD00][i] = m_d00vec[i];
            }

        for (int i = 0; i < m_fields[0]->GetExpSize(); ++i)
            {
                for (int j = 0; j < m_fields[0]->GetTotPoints(i); ++j)
                    {
                        index = m_fields[0]->GetPhys_Offset(i) + j;
                    }
                    std::cout << "elemid = " << i << ", x = " << x0[index] << ", Anisotropy[0] = " << m_d00vec[index] << std::endl;
            }


        // for (int i = 0; i < m_fields[0]->GetExpSize(); ++i)
        //     {
        //         xavg=0.0;
        //         for (int j = 0; j < m_fields[0]->GetTotPoints(i); ++j)
        //             {
        //                 index = m_fields[0]->GetPhys_Offset(i) + j;
        //                 xavg = xavg + x0[index];
        //             }

        //         for (int j = 0; j < m_fields[0]->GetTotPoints(i); ++j)
        //             {
        //                 index = m_fields[0]->GetPhys_Offset(i) + j;
        //                 if(xavg>0)
        //                 {
        //                     m_d00vec[index] = sqrt(m_d00);
        //                 }

        //                 else
        //                 {
        //                     m_d00vec[index] = 1.0/sqrt(m_d00);
        //                 }

        //             }
        //             std::cout << "elemid = " << i << ", x = " << x0[index] << ", Anisotropy[0] = " << m_d00vec[index] << std::endl;
        //     }

    }
    if (m_session->DefinesParameter("d11"))
    {
        m_varcoeffXYZ[StdRegions::eVarCoeffD11] = Array<OneD, NekDouble>(nq);
        m_d11vec = Array<OneD, NekDouble>(nq);

        Vmath::Fill(nq, m_d11, &m_d11vec[0], 1);

        Vmath::Vsqrt(nq, m_d11vec, 1, Anisotropy[1], 1);
        m_varcoeffXYZ[StdRegions::eVarCoeffD11] = Array<OneD, NekDouble>(nq, m_d11);
    }
    if (m_session->DefinesParameter("d22"))
    {
        m_varcoeffXYZ[StdRegions::eVarCoeffD22] = Array<OneD, NekDouble>(nq);
        Vmath::Fill(nq, sqrt(m_d22), &Anisotropy[2][0], 1);
    }

    MMFSystem::MMFInitObject(Anisotropy);

    // Define ProblemType
    if (m_session->DefinesSolverInfo("TESTTYPE"))
    {
        std::string TestTypeStr = m_session->GetSolverInfo("TESTTYPE");
        int i;
        for (i = 0; i < (int)SIZE_TestType; ++i)
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

    if (m_session->DefinesSolverInfo("INITWAVETYPE"))
    {
        std::string InitWaveTypeStr = m_session->GetSolverInfo("INITWAVETYPE");
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

    if(m_TestType==eTestPlaneEmbed)
    {
        // Let the moving frames outside the domain be of magnitude zero.
        int index, cnt = 0;
        for (int i = m_EmbededPlane; i < m_fields[0]->GetExpSize(); ++i)
        {
            for (int j = 0; j < m_fields[0]->GetTotPoints(i); ++j)
            {
                index = m_fields[0]->GetPhys_Offset(i) + j;
                for (int k=0; k<m_mfdim; ++k)
                {
                    m_movingframes[k][index] = 0.0;
                    m_movingframes[k][index+nq] = 0.0;
                    m_movingframes[k][index+2*nq] = 0.0;
                }
                cnt++;
            }
        }

        std::cout << "Moving frames " << cnt << " / " << nq << " ( " << 100.0*cnt/nq << " % ) are removed" << std::endl;
    }

    ComputeVarCoeff2D(m_movingframes,m_varcoeff);

    if(m_TestType==eTestPlaneEmbed)
    {
        std::cout << "Compute Varcoeff for phie ====================== " << std::endl;
        m_phiemovingframes = Array<OneD, Array<OneD, NekDouble>>(m_mfdim);
        m_phieMMFdir = FindMMFdir("LOCAL");
        Array<OneD, Array<OneD, NekDouble>> phieAniStrength(m_expdim);
        for (int j = 0; j < m_expdim; ++j)
        {
            phieAniStrength[j] = Array<OneD, NekDouble>(nq, 1.0);
        }

        SetUpMovingFrames(m_phieMMFdir, phieAniStrength, m_phiemovingframes);
                            
        StdRegions::VarCoeffMap m_phievarcoeff;

        ComputeVarCoeff2D(m_phiemovingframes, m_phievarcoeff);

        Array<OneD, NekDouble> ExactSoln(nq);

        TestHelmholtzProblem(0, m_phievarcoeff, m_fields[1]->UpdatePhys());                           
        TestHelmholtzProblem(1, m_phievarcoeff, ExactSoln);                           

        StdRegions::ConstFactorMap phiefactors;
        phiefactors[StdRegions::eFactorTau]    = m_Helmtau;
        phiefactors[StdRegions::eFactorLambda] = 0.0;

        // Compute phie distribution
        SetBoundaryConditions(0.0);

        m_fields[1]->HelmSolve(m_fields[1]->GetPhys(), m_fields[1]->UpdateCoeffs(), phiefactors, m_phievarcoeff);
        m_fields[1]->BwdTrans(m_fields[1]->GetCoeffs(), m_fields[1]->UpdatePhys());
        m_fields[1]->SetPhysState(true);

        Array<OneD, NekDouble> outarray(nq);
        Vmath::Vcopy(nq, m_fields[1]->UpdatePhys(), 1, outarray, 1);

        NekDouble phieavg = -1.0 * AvgInt(outarray);
        Vmath::Sadd(nq, phieavg, outarray, 1, outarray, 1);

        Array<OneD, NekDouble> tmp(nq);
        Vmath::Vsub(nq, outarray, 1, ExactSoln, 1, tmp, 1);
        std::cout << "SolveHelmholtzDiffusion Error = " << RootMeanSquare(tmp) << std::endl;

        Checkpoint_Output_Error(100, m_fields[1]->GetPhys(), ExactSoln);

        wait_on_enter();
    }

    // Test PhysDirectionalDeriv
    if(m_TestType==eTestPlaneAni)
    {
        TestPhysDirectionalDeriv(m_movingframes);
        TestHelmholtzSolver();

        wait_on_enter();
    }

    m_ode.DefineOdeRhs(&MMFDiffusion::DoOdeRhs, this);
    m_ode.DefineProjection(&MMFDiffusion::DoOdeProjection, this);
    m_ode.DefineImplicitSolve(&MMFDiffusion::DoImplicitSolve, this);
}

/**
 *
 */
MMFDiffusion::~MMFDiffusion()
{
}

void MMFDiffusion::DoOdeProjection(
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

/**OdeRhs
 * @param   inarray         Input array.
 * @param   outarray        Output array.
 * @param   time            Current simulation time.
 * @param   lambda          Timestep.
 */
void MMFDiffusion::DoImplicitSolve(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time,
    const NekDouble lambda)
{
    boost::ignore_unused(time);

    int nvariables = inarray.size();
    int nq         = m_fields[0]->GetNpoints();

    StdRegions::ConstFactorMap factors;
    factors[StdRegions::eFactorTau] = m_Helmtau;

    Array<OneD, Array<OneD, NekDouble>> F(nvariables);
    F[0] = Array<OneD, NekDouble>(nq * nvariables);
    for (int n = 1; n < nvariables; ++n)
    {
        F[n] = F[n - 1] + nq;
    }

    // We solve ( \nabla^2 - HHlambda ) Y[i] = rhs [i]
    // inarray = input: \hat{rhs} -> output: \hat{Y}
    // outarray = output: nabla^2 \hat{Y}
    // where \hat = modal coeffs
    // SetBoundaryConditions(time);
    for (int i = 0; i < nvariables; ++i)
    {
        factors[StdRegions::eFactorLambda] = 1.0 / lambda / m_epsu[i];

        // Multiply 1.0/timestep
        Vmath::Smul(nq, -factors[StdRegions::eFactorLambda], inarray[i], 1,
                    F[i], 1);

        m_fields[i]->HelmSolve(F[i], m_fields[i]->UpdateCoeffs(), factors,
                                  m_varcoeff);

        m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(), outarray[i]);
    }
}

/**
 *
 */
void MMFDiffusion::DoOdeRhs(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time)
{
    int nq = GetTotPoints();
    int nvar = m_fields.size();

    switch (m_TestType)
    {
        case eTestLineX:
        {
            Array<OneD, NekDouble> x(nq);
            Array<OneD, NekDouble> y(nq);
            Array<OneD, NekDouble> z(nq);

            m_fields[0]->GetCoords(x, y, z);

            for (int k = 0; k < nq; k++)
            {
                outarray[0][k] = (m_epsvec[k] * m_frequency * m_frequency - 1.0) * exp(-1.0 * time) * cos(m_frequency * x[k]);
            }
        }
        break;

        case eTestLineY:
        {
            Array<OneD, NekDouble> x(nq);
            Array<OneD, NekDouble> y(nq);
            Array<OneD, NekDouble> z(nq);

            m_fields[0]->GetCoords(x, y, z);

            for (int k = 0; k < nq; k++)
            {
                outarray[0][k] = (m_epsvec[k] * m_frequency * m_frequency - 1.0) * exp(-1.0 * time) * cos(m_frequency * y[k]);
            }
        }
        break;

        case eTestPlane:
        {
            Array<OneD, NekDouble> x(nq);
            Array<OneD, NekDouble> y(nq);
            Array<OneD, NekDouble> z(nq);

            m_fields[0]->GetCoords(x, y, z);

            for (int k = 0; k < nq; k++)
            {
                outarray[0][k] = (m_epsilon[0] * m_frequency *
                                 m_frequency + m_epsilon[1] * m_frequency *
                                 m_frequency - m_frequency) * exp(-1.0 * m_frequency * time) *
                                 sin(m_frequency * x[k]) * sin(m_frequency * y[k]);
            }
        }
        break;

        case eTestPlaneAni:
        {
            // StdRegions::VarCoeffType MMFCoeffs[15] = {
            //     StdRegions::eVarCoeffMF1x,   StdRegions::eVarCoeffMF1y,
            //     StdRegions::eVarCoeffMF1z,   StdRegions::eVarCoeffMF1Div,
            //     StdRegions::eVarCoeffMF1Mag, StdRegions::eVarCoeffMF2x,
            //     StdRegions::eVarCoeffMF2y,   StdRegions::eVarCoeffMF2z,
            //     StdRegions::eVarCoeffMF2Div, StdRegions::eVarCoeffMF2Mag,
            //     StdRegions::eVarCoeffMF3x,   StdRegions::eVarCoeffMF3y,
            //     StdRegions::eVarCoeffMF3z,   StdRegions::eVarCoeffMF3Div,
            //     StdRegions::eVarCoeffMF3Mag}; 

            Array<OneD, NekDouble> x(nq);
            Array<OneD, NekDouble> y(nq);
            Array<OneD, NekDouble> z(nq);

            m_fields[0]->GetCoords(x, y, z);

            // Array<OneD, NekDouble> d00(nq);
            // Array<OneD, NekDouble> d11(nq);

            // Vmath::Vcopy(nq, &m_varcoeff[MMFCoeffs[4]][0], 1, &d00[0], 1);
            // Vmath::Vcopy(nq, &m_varcoeff[MMFCoeffs[9]][0], 1, &d11[0], 1);

            for (int k = 0; k < nq; k++)
            {
                outarray[0][k] = (m_d00vec[k] * m_frequency *
                                 m_frequency + m_d11vec[k] * m_frequency *
                                 m_frequency - m_frequency) * exp(-1.0 * m_frequency * time) *
                                 sin(m_frequency * x[k]) * cos(m_frequency * y[k]);
            }
        }
        break;

        case eTestPlaneEmbed:
        {
            Array<OneD, NekDouble> x(nq);
            Array<OneD, NekDouble> y(nq);
            Array<OneD, NekDouble> z(nq);

            m_fields[0]->GetCoords(x, y, z);

            int index;
            outarray[0] = Array<OneD, NekDouble>(nq, 0.0);
            for (int i = 0; i < m_EmbededPlane; ++i)
            {
                for (int j = 0; j < m_fields[0]->GetTotPoints(i); ++j)
                {
                    index = m_fields[0]->GetPhys_Offset(i) + j;
                    outarray[0][index] = (m_epsilon[0] * m_frequency *
                                        m_frequency + m_epsilon[1] * m_frequency *
                                        m_frequency - m_frequency) * exp(-1.0 * m_frequency * time) *
                                        cos(m_frequency * x[index]) * cos(m_frequency * y[index]);        
                }
            }
        }
        break;

        case eTestPlaneNeumann:
        {
            Array<OneD, NekDouble> x(nq);
            Array<OneD, NekDouble> y(nq);
            Array<OneD, NekDouble> z(nq);

            m_fields[0]->GetCoords(x, y, z);

            for (int k = 0; k < nq; k++)
            {
                outarray[0][k] = (m_epsilon[0] * m_frequency *
                                 m_frequency + m_epsilon[1] * m_frequency *
                                 m_frequency - m_frequency) * exp(-1.0 * m_frequency * time) *
                                 cos(m_frequency * x[k]) * cos(m_frequency * y[k]);
            }
        }
        break;


        case eTestCube:
        {

            Array<OneD, NekDouble> x(nq);
            Array<OneD, NekDouble> y(nq);
            Array<OneD, NekDouble> z(nq);

            m_fields[0]->GetCoords(x, y, z);

            for (int k = 0; k < nq; k++)
            {
                outarray[0][k] =
                    (m_epsilon[0] + m_epsilon[1] + m_epsilon[2] - 1.0) * m_frequency *
                    m_frequency * exp(-1.0 * m_frequency * m_frequency * time) * sin(m_frequency * x[k]) *
                    sin(m_frequency * y[k]) * sin(m_frequency * z[k]);
            }
        }
        break;

        case eTestLinearSphere:
        {
            Array<OneD, NekDouble> temp(nq);

            NekDouble A = 2.0;
            NekDouble B = 5.0;

            NekDouble m_a, m_b, m_c, m_d;
            m_a = B - 1.0;
            m_b = A * A;
            m_c = -1.0 * B;
            m_d = -1.0 * A * A;

            temp = Array<OneD, NekDouble>(nq, 0.0);
            Vmath::Svtvp(nq, m_a, &inarray[0][0], 1, &temp[0], 1, &temp[0], 1);
            Vmath::Svtvp(nq, m_b, &inarray[1][0], 1, &temp[0], 1,
                         &outarray[0][0], 1);

            temp = Array<OneD, NekDouble>(nq, 0.0);
            Vmath::Svtvp(nq, m_c, &inarray[0][0], 1, &temp[0], 1, &temp[0], 1);
            Vmath::Svtvp(nq, m_d, &inarray[1][0], 1, &temp[0], 1,
                         &outarray[1][0], 1);
        }
        break;

        case eTestNonlinearSphere:
        {
            NekDouble A = 2.0;
            NekDouble B = 5.0;

            Array<OneD, NekDouble> Aonevec(nq, A);

            // cube = phys0*phys0*phy1
            Array<OneD, NekDouble> cube(nq);
            Vmath::Vmul(nq, &inarray[0][0], 1, &inarray[0][0], 1, &cube[0], 1);
            Vmath::Vmul(nq, &inarray[1][0], 1, &cube[0], 1, &cube[0], 1);

            // outarray[0] = A - B*phy0 + phy0*phy0*phy1 - phy0
            NekDouble coeff = -1.0 * B - 1.0;
            Array<OneD, NekDouble> tmp(nq);
            Vmath::Svtvp(nq, coeff, &inarray[0][0], 1, &cube[0], 1, &tmp[0], 1);
            Vmath::Vadd(nq, &Aonevec[0], 1, &tmp[0], 1, &outarray[0][0], 1);

            // outarray[1] = B*phys0 - phy0*phy0*phy1
            Vmath::Svtvm(nq, B, &inarray[0][0], 1, &cube[0], 1, &outarray[1][0],
                         1);
        }
        break;

        default:
            break;
    }

    if (m_explicitDiffusion)
    {
        switch (m_projectionType)
        {
            case MultiRegions::eDiscontinuous:
            {
                std::string diffName;

                // Do not forwards transform initial condition
                m_homoInitialFwd = false;

                m_session->LoadSolverInfo("DiffusionType", diffName, "LDG");
                m_diffusion = SolverUtils::GetDiffusionFactory().CreateInstance(
                    diffName, diffName);
                m_diffusion->SetFluxVector(&MMFDiffusion::GetFluxVector, this);
                m_diffusion->InitObject(m_session, m_fields);
                break;
            }

            case MultiRegions::eGalerkin:
            case MultiRegions::eMixed_CG_Discontinuous:
            {
                if (m_explicitDiffusion)
                {
                    ASSERTL0(false, "Explicit Galerkin diffusion not set up.");
                }
            }
        }

        Array<OneD, Array<OneD, NekDouble>> outarrayDiff(nvar);
        for (int i = 0; i < nvar; ++i)
        {
            outarrayDiff[i] = Array<OneD, NekDouble>(nq, 0.0);
        }

        m_diffusion->Diffuse(nvar, m_fields, inarray, outarrayDiff);                         

        for (int i = 0; i < nvar; ++i)
        {
            Vmath::Vvtvp(nq, &m_epsvec[0], 1, &outarrayDiff[i][0], 1, &outarray[i][0], 1, &outarray[i][0], 1);
        }
    }

    // if (m_explicitDiffusion)
    // {
    //     int nq = m_fields[0]->GetNpoints();

    //     // Laplacian only to the first variable
    //     Array<OneD, NekDouble> Laplacian(nq);

    //     for (int i=0; i < nvar; ++i)
    //     {
    //         WeakDGMMFDiffusion(i, inarray[i], Laplacian, time);
    //         Vmath::Smul(nq, m_epsu[i], &Laplacian[0], 1, &Laplacian[0], 1);
    //         Vmath::Vadd(nq, &Laplacian[0], 1, &outarray[i][0], 1, &outarray[i][0], 1);
    //     }
    // }
}

void MMFDiffusion::v_SetInitialConditions(NekDouble initialtime,
                                          bool dumpInitialConditions,
                                          const int domain)
{
    boost::ignore_unused(domain);

    int nq = GetTotPoints();

    switch (m_TestType)
    {
        case eTestLineX:
        {
            Array<OneD, NekDouble> u(nq);

            TestLineProblem(0, initialtime, u);
            m_fields[0]->SetPhys(u);
        }
        break;

        case eTestLineY:
        {
            Array<OneD, NekDouble> u(nq);

            TestLineProblem(1, initialtime, u);
            m_fields[0]->SetPhys(u);
        }
        break;

        case eTestPlane:
        {
            Array<OneD, NekDouble> u(nq);

            TestPlaneProblem(initialtime, m_varcoeff, u);
            m_fields[0]->SetPhys(u);
        }
        break;

        case eTestPlaneAni:
        {
            Array<OneD, NekDouble> u(nq);

            TestPlaneAniProblem(initialtime, m_varcoeff, u);
            m_fields[0]->SetPhys(u);
        }
        break;

        case eTestPlaneEmbed:
        {
            Array<OneD, NekDouble> u(nq);

            TestPlaneEmbedProblem(initialtime, m_varcoeff, u);
            m_fields[0]->SetPhys(u);
        }
        break;


        case eTestPlaneNeumann:
        {
            Array<OneD, NekDouble> u(nq);

            TestPlaneNeumannProblem(initialtime, u);
            m_fields[0]->SetPhys(u);
        }
        break;


        case eTestCube:
        {
            Array<OneD, NekDouble> u(nq);

            TestCubeProblem(initialtime, u);
            m_fields[0]->SetPhys(u);
        }
        break;

        case eTestLinearSphere:
        case eTestNonlinearSphere:
        {
            Array<OneD, NekDouble> u(nq);
            Array<OneD, NekDouble> v(nq);

            Morphogenesis(initialtime, 0, u);
            Morphogenesis(initialtime, 1, v);

            m_fields[0]->SetPhys(u);
            m_fields[1]->SetPhys(v);
        }
        break;

        default:
        {
            EquationSystem::v_SetInitialConditions(initialtime, false);
        }
        break;
    }

    std::cout << "Initial: max u = "
              << Vmath::Vmax(nq, m_fields[0]->GetPhys(), 1) << ", u_L2 = " << RootMeanSquare(m_fields[0]->GetPhys()) << std::endl;

    if (dumpInitialConditions)
    {
        std::string outname = m_sessionName + "_initial.chk";
        WriteFld(outname);

        Checkpoint_Output_Error(0,m_fields[0]->GetPhys(),m_fields[0]->GetPhys());
    }
}

void MMFDiffusion::TestPhysDirectionalDeriv(const Array<OneD, const Array<OneD, NekDouble>> &movingframes)

{
    int nq = GetTotPoints();

    std::cout << "Test PhysDirectionalDeriv, ===================================== " << std::endl;

    Array<OneD, NekDouble> tmp(nq);
    
    Array<OneD, NekDouble> Dxtmp(nq);
    Array<OneD, NekDouble> Dytmp(nq);

    Array<OneD, NekDouble> ExtDxtmp(nq);
    Array<OneD, NekDouble> ExtDytmp(nq);

    TestPlaneAniProblem(0.0, m_varcoeff, tmp);

    MMFDirectionalDeriv(movingframes[0], tmp, Dxtmp);
    MMFDirectionalDeriv(movingframes[1], tmp, Dytmp);

    TestPlaneAniDerivProblem(0.0, m_varcoeff, ExtDxtmp, ExtDytmp);

    Vmath::Vsub(nq, Dxtmp, 1, ExtDxtmp, 1, Dxtmp, 1);
    Vmath::Vsub(nq, Dytmp, 1, ExtDytmp, 1, Dytmp, 1);
    std::cout << "Dx error = " << RootMeanSquare(Dxtmp) << ", Dy error = " << RootMeanSquare(Dytmp) << std::endl << std::endl;
}

void MMFDiffusion::TestLineProblem(const int direction, const NekDouble time,
                                    Array<OneD, NekDouble> &outfield)

{
    int nq = GetTotPoints();

    Array<OneD, NekDouble> x(nq);
    Array<OneD, NekDouble> y(nq);
    Array<OneD, NekDouble> z(nq);

    m_fields[0]->GetCoords(x, y, z);

    outfield = Array<OneD, NekDouble>(nq);
    for (int k = 0; k < nq; k++)
    {
        if(direction==0)
        {
            outfield[k] = exp(-1.0 * time) * cos(m_frequency * x[k]);
        }

        else if (direction==1)
        {
            outfield[k] = exp(-1.0 * time) * cos(m_frequency * y[k]);
        }
    }
}

void MMFDiffusion::TestPlaneProblem(const NekDouble time,
                                    StdRegions::VarCoeffMap &varcoeff,
                                    Array<OneD, NekDouble> &outfield)
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

    Array<OneD, NekDouble> x(nq);
    Array<OneD, NekDouble> y(nq);
    Array<OneD, NekDouble> z(nq);

    m_fields[0]->GetCoords(x, y, z);

    Array<OneD, NekDouble> d00(nq);
    Array<OneD, NekDouble> d11(nq);

    Vmath::Vcopy(nq, &varcoeff[MMFCoeffs[4]][0], 1, &d00[0], 1);
    Vmath::Vcopy(nq, &varcoeff[MMFCoeffs[9]][0], 1, &d11[0], 1);

    outfield = Array<OneD, NekDouble>(nq);
    for (int k = 0; k < nq; k++)
    {
        outfield[k] = exp(-1.0 * d00[k] * d11[k] * m_frequency * time) * sin(m_frequency * x[k]) * sin(m_frequency * y[k]);
    }
}

void MMFDiffusion::TestPlaneEmbedProblem(const NekDouble time,
                                    StdRegions::VarCoeffMap &varcoeff,
                                    Array<OneD, NekDouble> &outfield)
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

    Array<OneD, NekDouble> x(nq);
    Array<OneD, NekDouble> y(nq);
    Array<OneD, NekDouble> z(nq);

    m_fields[0]->GetCoords(x, y, z);

    Array<OneD, NekDouble> d00(nq);
    Array<OneD, NekDouble> d11(nq);

    Vmath::Vcopy(nq, &varcoeff[MMFCoeffs[4]][0], 1, &d00[0], 1);
    Vmath::Vcopy(nq, &varcoeff[MMFCoeffs[9]][0], 1, &d11[0], 1);

    int index;
    outfield = Array<OneD, NekDouble>(nq,0.0);
    for (int i = 0; i < m_EmbededPlane; ++i)
    {
        for (int j = 0; j < m_fields[0]->GetTotPoints(i); ++j)
        {
            index = m_fields[0]->GetPhys_Offset(i) + j;
            outfield[index] = exp(-1.0 * m_frequency * time) * cos(m_frequency * x[index]) * cos(m_frequency * y[index]);
        }
    }
}


void MMFDiffusion::TestPlaneNeumannProblem(const NekDouble time,
                                    Array<OneD, NekDouble> &outfield)

{
    int nq = GetTotPoints();

    Array<OneD, NekDouble> x(nq);
    Array<OneD, NekDouble> y(nq);
    Array<OneD, NekDouble> z(nq);

    m_fields[0]->GetCoords(x, y, z);

    outfield = Array<OneD, NekDouble>(nq);
    for (int k = 0; k < nq; k++)
    {
        outfield[k] = exp(-1.0 * m_frequency * time) * cos(m_frequency * x[k]) * cos(m_frequency * y[k]);
    }
}

// \nabla \cdot D \nabla u = f 
// type = 0 (f):  Helmholtz forcing
// type = 1 (u): Exact solution of $u$.
void MMFDiffusion::TestHelmholtzProblem(const int type,
                                    StdRegions::VarCoeffMap &varcoeff,
                                    Array<OneD, NekDouble> &outfield)
{
    boost::ignore_unused(varcoeff);

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

    Array<OneD, NekDouble> x(nq);
    Array<OneD, NekDouble> y(nq);
    Array<OneD, NekDouble> z(nq);

    m_fields[0]->GetCoords(x, y, z);

    Array<OneD, NekDouble> d00(nq);
    Array<OneD, NekDouble> d11(nq);

    Vmath::Vcopy(nq, &varcoeff[MMFCoeffs[4]][0], 1, &d00[0], 1);
    Vmath::Vcopy(nq, &varcoeff[MMFCoeffs[9]][0], 1, &d11[0], 1);

    outfield = Array<OneD, NekDouble>(nq);
    for (int k = 0; k < nq; k++)
    {
        // Helmholtz forcing
        if(type==0)
        {
            outfield[k] = -1.0 * (d00[k] + d11[k]) * m_frequency * m_frequency * sin(m_frequency * x[k]) * cos(m_frequency * y[k]);
        }

        // Helmholtz solution
        else if(type==1)
        {
            outfield[k] = sin(m_frequency * x[k]) * cos(m_frequency * y[k]);
        }
    }
}




void MMFDiffusion::TestPlaneAniProblem(const NekDouble time,
                                    StdRegions::VarCoeffMap &varcoeff,
                                    Array<OneD, NekDouble> &outfield)
{
    boost::ignore_unused(varcoeff);

    int nq = GetTotPoints();

    Array<OneD, NekDouble> x(nq);
    Array<OneD, NekDouble> y(nq);
    Array<OneD, NekDouble> z(nq);

    m_fields[0]->GetCoords(x, y, z);

    outfield = Array<OneD, NekDouble>(nq);
    for (int k = 0; k < nq; k++)
    {
        outfield[k] = exp(-1.0 * m_frequency * time) * sin(m_frequency * x[k]) * cos(m_frequency * y[k]);
    }
}

void MMFDiffusion::TestPlaneAniDerivProblem(const NekDouble time,
                                    StdRegions::VarCoeffMap &varcoeff,
                                    Array<OneD, NekDouble> &Dxoutfield,
                                    Array<OneD, NekDouble> &Dyoutfield)
{
    boost::ignore_unused(varcoeff);

    int nq = GetTotPoints();

    // StdRegions::VarCoeffType MMFCoeffs[3] = {
    //     StdRegions::eVarCoeffD00, StdRegions::eVarCoeffD11, StdRegions::eVarCoeffD22}; 

    //     Vmath::Vcopy(nq, &varcoeff[MMFCoeffs[0]][0], 1, &d00[0], 1);
    // Vmath::Vcopy(nq, &varcoeff[MMFCoeffs[1]][0], 1, &d11[0], 1);

    StdRegions::VarCoeffType MMFCoeffs[15] = {
        StdRegions::eVarCoeffMF1x,   StdRegions::eVarCoeffMF1y,
        StdRegions::eVarCoeffMF1z,   StdRegions::eVarCoeffMF1Div,
        StdRegions::eVarCoeffMF1Mag, StdRegions::eVarCoeffMF2x,
        StdRegions::eVarCoeffMF2y,   StdRegions::eVarCoeffMF2z,
        StdRegions::eVarCoeffMF2Div, StdRegions::eVarCoeffMF2Mag,
        StdRegions::eVarCoeffMF3x,   StdRegions::eVarCoeffMF3y,
        StdRegions::eVarCoeffMF3z,   StdRegions::eVarCoeffMF3Div,
        StdRegions::eVarCoeffMF3Mag}; 

    Array<OneD, NekDouble> x(nq);
    Array<OneD, NekDouble> y(nq);
    Array<OneD, NekDouble> z(nq);

    m_fields[0]->GetCoords(x, y, z);

    Array<OneD, NekDouble> d00(nq);
    Array<OneD, NekDouble> d11(nq);

    Vmath::Vsqrt(nq, &varcoeff[MMFCoeffs[4]][0], 1, &d00[0], 1);
    Vmath::Vsqrt(nq, &varcoeff[MMFCoeffs[9]][0], 1, &d11[0], 1);

    Dxoutfield = Array<OneD, NekDouble>(nq);
    Dyoutfield = Array<OneD, NekDouble>(nq);

    for (int k = 0; k < nq; k++)
    {
        Dxoutfield[k] = m_frequency * d00[k] * exp(-1.0 * m_frequency * time) * cos(m_frequency * x[k]) * cos(m_frequency * y[k]);
        Dyoutfield[k] = -m_frequency * d11[k] * exp(-1.0 * m_frequency * time) * sin(m_frequency * x[k]) * sin(m_frequency * y[k]);
    }
}

void MMFDiffusion::TestCubeProblem(const NekDouble time,
                                   Array<OneD, NekDouble> &outfield)

{
    int nq = GetTotPoints();

    Array<OneD, NekDouble> x(nq);
    Array<OneD, NekDouble> y(nq);
    Array<OneD, NekDouble> z(nq);

    m_fields[0]->GetCoords(x, y, z);

    outfield = Array<OneD, NekDouble>(nq);
    for (int k = 0; k < nq; k++)
    {
        outfield[k] = exp(-1.0 * m_frequency * m_frequency * time) * sin(m_frequency * x[k]) *
                      sin(m_frequency * y[k]) * sin(m_frequency * z[k]);
    }
}

void MMFDiffusion::Morphogenesis(const NekDouble time, unsigned int field,
                                 Array<OneD, NekDouble> &outfield)
{
    int nq = GetTotPoints();

    int i, m, n, ind;
    NekDouble a_n, d_n, gamma_n;
    NekDouble A_mn, C_mn, theta, phi, radius;

    std::complex<double> Spericharmonic, delta_n, temp;
    std::complex<double> varphi0, varphi1;
    std::complex<double> B_mn, D_mn;

    // Set some parameter values
    int Maxn = 6;
    int Maxm = 2 * Maxn - 1;

    NekDouble A = 2.0;
    NekDouble B = 5.0;

    NekDouble m_mu = 0.001;
    NekDouble m_nu = 0.002;

    NekDouble m_a, m_b, m_c, m_d;

    m_a = B - 1.0;
    m_b = A * A;
    m_c = -1.0 * B;
    m_d = -1.0 * A * A;

    Array<OneD, Array<OneD, NekDouble>> Ainit(Maxn);
    Array<OneD, Array<OneD, NekDouble>> Binit(Maxn);

    for (i = 0; i < Maxn; ++i)
    {
        Ainit[i] = Array<OneD, NekDouble>(Maxm, 0.0);
        Binit[i] = Array<OneD, NekDouble>(Maxm, 0.0);
    }

    Ainit[5][0]  = -0.5839;
    Ainit[5][1]  = -0.8436;
    Ainit[5][2]  = -0.4764;
    Ainit[5][3]  = 0.6475;
    Ainit[5][4]  = 0.1886;
    Ainit[5][5]  = 0.8709;
    Ainit[5][6]  = -0.8338;
    Ainit[5][7]  = 0.1795;
    Ainit[5][8]  = -0.7873;
    Ainit[5][9]  = 0.8842;
    Ainit[5][10] = 0.2943;

    Binit[5][0]  = -0.6263;
    Binit[5][1]  = 0.9803;
    Binit[5][2]  = 0.7222;
    Binit[5][3]  = 0.5945;
    Binit[5][4]  = 0.6026;
    Binit[5][5]  = -0.2076;
    Binit[5][6]  = 0.4556;
    Binit[5][7]  = 0.6024;
    Binit[5][8]  = 0.9695;
    Binit[5][9]  = -0.4936;
    Binit[5][10] = 0.1098;

    Array<OneD, NekDouble> u(nq);
    Array<OneD, NekDouble> v(nq);
    Array<OneD, NekDouble> x(nq);
    Array<OneD, NekDouble> y(nq);
    Array<OneD, NekDouble> z(nq);

    m_fields[0]->GetCoords(x, y, z);
    for (int i = 0; i < nq; ++i)
    {
        radius = sqrt(x[i] * x[i] + y[i] * y[i] + z[i] * z[i]);

        // theta is in [0, pi]
        theta = asin(z[i] / radius) + 0.5 * m_pi;

        // phi is in [0, 2*pi]
        phi = atan2(y[i], x[i]) + m_pi;

        varphi0 = 0.0 * varphi0;
        varphi1 = 0.0 * varphi1;
        for (n = 0; n < Maxn; ++n)
        {
            // Set up parameters
            a_n = m_a - m_mu * (n * (n + 1) / radius / radius);
            d_n = m_d - m_nu * (n * (n + 1) / radius / radius);

            gamma_n = 0.5 * (a_n + d_n);

            temp    = (a_n + d_n) * (a_n + d_n) - 4.0 * (a_n * d_n - m_b * m_c);
            delta_n = 0.5 * sqrt(temp);

            for (m = -n; m <= n; ++m)
            {
                ind  = m + n;
                A_mn = Ainit[n][ind];
                C_mn = Binit[n][ind];

                B_mn = ((a_n - gamma_n) * Ainit[n][ind] + m_b * Binit[n][ind]) /
                       delta_n;
                D_mn = (m_c * Ainit[n][ind] + (d_n - gamma_n) * Binit[n][ind]) /
                       delta_n;

                Spericharmonic =
                    boost::math::spherical_harmonic(n, m, theta, phi);
                varphi0 += exp(gamma_n * time) *
                           (A_mn * cosh(delta_n * time) +
                            B_mn * sinh(delta_n * time)) *
                           Spericharmonic;
                varphi1 += exp(gamma_n * time) *
                           (C_mn * cosh(delta_n * time) +
                            D_mn * sinh(delta_n * time)) *
                           Spericharmonic;
            }
        }

        u[i] = varphi0.real();
        v[i] = varphi1.real();
    }

    switch (field)
    {
        case 0:
        {
            outfield = u;
        }
        break;

        case 1:
        {
            outfield = v;
        }
        break;
    }
}

Array<OneD, NekDouble> MMFDiffusion::PlanePhiWave()
{
    int nq = GetTotPoints();
    Array<OneD, NekDouble> outarray(nq, 0.0);

    Array<OneD, NekDouble> x(nq);
    Array<OneD, NekDouble> y(nq);
    Array<OneD, NekDouble> z(nq);

    m_fields[0]->GetCoords(x, y, z);

    NekDouble xmin, ymin, xmax;

    xmin = Vmath::Vmin(nq, x, 1);
    xmax = Vmath::Vmax(nq, x, 1);
    ymin = Vmath::Vmin(nq, y, 1);

    NekDouble xp, yp, xp2;
    for (int i = 0; i < nq; i++)
    {
        switch (m_InitWaveType)
        {
            case eLeft:
            {
                NekDouble radiusofinit = 4.0;
                NekDouble frontstiff   = 0.1;

                xp = x[i] - xmin;
                outarray[i] =
                    1.0 / (1.0 + exp((xp - radiusofinit) / frontstiff));
            }
            break;

            case eBothEnds:
            {
                NekDouble radiusofinit = 3.0;
                NekDouble frontstiff   = 0.1;

                xp  = x[i] - xmin;
                xp2 = x[i] - xmax;

                outarray[i] =
                    1.0 / (1.0 +
                           exp((sqrt(xp * xp) - radiusofinit) / frontstiff)) +
                    1.0 / (1.0 +
                           exp((sqrt(xp2 * xp2) - radiusofinit) / frontstiff));
            }
            break;

            case eCenter:
            {
                NekDouble radiusofinit = 6.0;
                NekDouble frontstiff   = 0.1;

                xp = x[i] - xmin;
                outarray[i] =
                    1.0 / (1.0 + exp((xp - radiusofinit) / frontstiff));
            }
            break;

            case eLeftBottomCorner:
            {
                NekDouble radiusofinit = 6.0;
                NekDouble frontstiff   = 0.1;
                NekDouble bs           = 2.0;

                xp = x[i] - xmin;
                yp = y[i] - ymin;
                outarray[i] =
                    1.0 /
                    (1.0 + exp((sqrt(xp * xp + yp * yp) / bs - radiusofinit) /
                               frontstiff));
            }
            break;

            case ePoint:
            {
                NekDouble xloc, yloc, zloc, rad;
                NekDouble radiusofinit = 10.0;

                xloc = x[i] - m_InitPtx;
                yloc = y[i] - m_InitPty;
                zloc = z[i] - m_InitPtz;

                rad = sqrt(xloc * xloc + yloc * yloc + zloc * zloc);

                xloc = xloc / radiusofinit;
                yloc = yloc / radiusofinit;
                zloc = zloc / radiusofinit;

                if (rad < radiusofinit)
                {
                    outarray[i] =
                        exp(-(1.0 / 2.0) *
                            (xloc * xloc + yloc * yloc + zloc * zloc));
                }

                else
                {
                    outarray[i] = 0.0;
                }
            }
            break;

            case eSpiralDock:
            {
                NekDouble radiusofinit = 3.0;
                NekDouble frontstiff   = 0.1;
                xp                     = x[i] - 4.0;
                yp                     = y[i];
                outarray[i] =
                    (1.0 / (1.0 + exp(2.0 * yp))) *
                    (1.0 / (1.0 + exp(-2.0 * xp))) *
                    (1.0 / (1.0 + exp((xp - radiusofinit) / frontstiff)));
            }
            break;

            default:
                break;
        }
    }

    return outarray;
}

void MMFDiffusion::v_EvaluateExactSolution(unsigned int field,
                                           Array<OneD, NekDouble> &outfield,
                                           const NekDouble time)
{
    switch (m_TestType)
    {
        case eTestLineX:
        {
            TestLineProblem(0, time, outfield);
        }
        break;

        case eTestLineY:
        {
            TestLineProblem(1, time, outfield);
        }
        break;

        case eTestPlane:
        {
            TestPlaneProblem(time, m_varcoeff, outfield);
        }
        break;

        case eTestPlaneEmbed:
        {
            TestPlaneEmbedProblem(time, m_varcoeff, outfield);
        }
        break;

        case eTestPlaneNeumann:
        {
            TestPlaneNeumannProblem(time, outfield);
        }
        break;

        case eTestPlaneAni:
        {
            TestPlaneAniProblem(time, m_varcoeff, outfield);
        }
        break;

        case eTestCube:
        {
            TestCubeProblem(time, outfield);
        }
        break;

        case eTestLinearSphere:
        case eTestNonlinearSphere:
        {
            Morphogenesis(time, field, outfield);
        }
        break;

            // case eFHNStandard:
            // case eFHNRogers:
            // case eFHNAlievPanf:
            // {
            //     int nq   = GetTotPoints();
            //     outfield = Array<OneD, NekDouble>(nq, 0.0);
            // }
            /* Falls through. */
        default:
        {
            EquationSystem::v_EvaluateExactSolution(field, outfield, time);
        }
        break;
    }
}

// void MMFDiffusion::ComputeVarCoeff2D(
//     const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
//     StdRegions::VarCoeffMap &varcoeff)
// {
//     int nq = GetTotPoints();

//     StdRegions::VarCoeffType MMFCoeffs[15] = {
//         StdRegions::eVarCoeffMF1x,   StdRegions::eVarCoeffMF1y,
//         StdRegions::eVarCoeffMF1z,   StdRegions::eVarCoeffMF1Div,
//         StdRegions::eVarCoeffMF1Mag, StdRegions::eVarCoeffMF2x,
//         StdRegions::eVarCoeffMF2y,   StdRegions::eVarCoeffMF2z,
//         StdRegions::eVarCoeffMF2Div, StdRegions::eVarCoeffMF2Mag,
//         StdRegions::eVarCoeffMF3x,   StdRegions::eVarCoeffMF3y,
//         StdRegions::eVarCoeffMF3z,   StdRegions::eVarCoeffMF3Div,
//         StdRegions::eVarCoeffMF3Mag};

//     int indx;
//     Array<OneD, NekDouble> tmp(nq);
//     for (int k = 0; k < m_expdim; ++k)
//     {
//         // For Moving Frames
//         indx = 5 * k;

//         for (int j = 0; j < m_spacedim; ++j)
//         {
//             varcoeff[MMFCoeffs[indx + j]] = Array<OneD, NekDouble>(nq, 0.0);
//             Vmath::Vcopy(nq, &movingframes[k][j * nq], 1,
//                          &varcoeff[MMFCoeffs[indx + j]][0], 1);
//         }

//         // m_DivMF
//         varcoeff[MMFCoeffs[indx + 3]] = Array<OneD, NekDouble>(nq, 0.0);

//         Array<OneD, Array<OneD, NekDouble>> DivMF;
//         // ComputeDivMF(eCovariant, movingframes, DivMF);

//         ComputeEuclideanDivMF(movingframes, DivMF);

//         Vmath::Vcopy(nq, &DivMF[k][0], 1, &varcoeff[MMFCoeffs[indx + 3]][0], 1);
//         // \| e^k \|
//         varcoeff[MMFCoeffs[indx + 4]] = Array<OneD, NekDouble>(nq, 0.0);
//         tmp                           = Array<OneD, NekDouble>(nq, 0.0);
//         for (int i = 0; i < m_spacedim; ++i)
//         {
//             Vmath::Vvtvp(nq, &movingframes[k][i * nq], 1,
//                          &movingframes[k][i * nq], 1, &tmp[0], 1, &tmp[0], 1);
//         }

//         Vmath::Vcopy(nq, &tmp[0], 1, &varcoeff[MMFCoeffs[indx + 4]][0], 1);
//     }

//     std::cout << "m_varcoeff = " << RootMeanSquare(varcoeff[MMFCoeffs[0]])
//               << " , " << RootMeanSquare(varcoeff[MMFCoeffs[1]]) << " , "
//               << RootMeanSquare(varcoeff[MMFCoeffs[2]]) << " , "
//               << RootMeanSquare(varcoeff[MMFCoeffs[3]]) << " , "
//               << RootMeanSquare(varcoeff[MMFCoeffs[4]]) << std::endl;

//     std::cout << " ::::: 2D Varcoeff is Successfully Created ::::: "
//               << std::endl;
// }

void MMFDiffusion::GetFluxVector(
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
            Vmath::Smul(nPts, m_epsilon[j], qfield[j][i], 1, viscousTensor[j][i],
                        1);
        }
    }
}


void MMFDiffusion::ComputeEuclideanDivMF(
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    Array<OneD, Array<OneD, NekDouble>> &DivMF)
{
    int nq = m_fields[0]->GetNpoints();

    DivMF = Array<OneD, Array<OneD, NekDouble>>(m_expdim);
    for (int j = 0; j < m_expdim; ++j)
    {
        DivMF[j] = Array<OneD, NekDouble>(nq, 0.0);
    }

    Array<OneD, NekDouble> tmp(nq);
    Array<OneD, NekDouble> Dtmp(nq);

    // case eEuclidean:
    for (int j = 0; j < m_expdim; ++j)
    {
        for (int k = 0; k < m_spacedim; ++k)
        {
            Vmath::Vcopy(nq, &movingframes[j][k * nq], 1, &tmp[0], 1);
            m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[k], tmp, Dtmp);
            Vmath::Vadd(nq, &Dtmp[0], 1, &DivMF[j][0], 1, &DivMF[j][0], 1);
        }
    }
}

void MMFDiffusion::TestHelmholtzSolver()
{
    int nq               = GetTotPoints();
    int nvar  = m_fields.size();

    int ncoeffs = m_fields[0]->GetNcoeffs();

    Array<OneD, NekDouble> tmpcXYZ(ncoeffs);
    Array<OneD, NekDouble> tmpc(ncoeffs);
    Array<OneD, NekDouble> HelmSolve(nq);
    Array<OneD, NekDouble> HelmMMFSolve(nq);
    Array<OneD, NekDouble> HelmXYZSolve(nq);

    StdRegions::ConstFactorMap factors;
    factors[StdRegions::eFactorTau] = m_Helmtau;
    factors[StdRegions::eFactorLambda] = 0.0;

    std::cout << "Test HelmholtzSolver ===============================================" << std::endl;

    Array<OneD, Array<OneD, NekDouble>> ExactSoln(nvar);

    TestHelmholtzProblem(0,m_varcoeff,m_fields[0]->UpdatePhys());                           
    TestHelmholtzProblem(1,m_varcoeff,ExactSoln[0]);                           

    // Zero field so initial conditions are zero
    Vmath::Zero(m_fields[0]->GetNcoeffs(), m_fields[0]->UpdateCoeffs(), 1);
    TestHelmholtzProblem(0,m_varcoeff,m_fields[0]->UpdatePhys());     

    // SetBoundaryConditions(0.0);                      
    // m_fields[0]->HelmSolve(m_fields[0]->GetPhys(), m_fields[0]->UpdateCoeffs(), factors, m_varcoeffXYZ);
    // m_fields[0]->SetPhysState(false);
    // m_fields[0]->BwdTrans(m_fields[0]->GetCoeffs(), HelmXYZSolve);
    // std::cout << "HelmsolveXYZ done ..................." << std::endl;

    Vmath::Zero(m_fields[0]->GetNcoeffs(), m_fields[0]->UpdateCoeffs(), 1);
    TestHelmholtzProblem(0,m_varcoeff,m_fields[0]->UpdatePhys());                           

    SetBoundaryConditions(0.0);
    m_fields[0]->HelmSolve(m_fields[0]->GetPhys(), m_fields[0]->UpdateCoeffs(), factors, m_varcoeff);
    m_fields[0]->SetPhysState(false);
    m_fields[0]->BwdTrans(m_fields[0]->GetCoeffs(), HelmSolve);
    std::cout << "HelmsolveMMF done ..................." << std::endl;

    Array<OneD, NekDouble> Error(nq,0.0);
    Vmath::Vsub(nq, &ExactSoln[0][0], 1, &HelmSolve[0], 1, &Error[0], 1);

    std::cout << "HelmSolve Error = " << RootMeanSquare(Error) << std::endl;

    Array<OneD, NekDouble> tmp(nq);
    Array<OneD, NekDouble> D2tmp(nq);

    TestHelmholtzProblem(0, m_varcoeff, ExactSoln[0]);   
    TestHelmholtzProblem(1, m_varcoeff, m_fields[0]->UpdatePhys());   
    D2tmp = ComputeCovariantDiffusion(m_movingframes, m_fields[0]->GetPhys());

    Vmath::Vsub(nq, &ExactSoln[0][0], 1, &D2tmp[0], 1, &Error[0], 1);

    std::cout << "Error: CovariantDiffusion = " << RootMeanSquare(Error) << " for Lap_mag = " 
    << RootMeanSquare(D2tmp) << std::endl << std::endl;
}

void MMFDiffusion::v_DoSolve()
{
    ASSERTL0(m_intScheme != 0, "No time integration scheme.");

    int i, nchk = 1;
    int nvariables       = 0;
    int nfields          = m_fields.size();
    std::string fulltext = ""; // initiate fulltext
    int nq = m_fields[0]->GetNpoints();

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

    if(m_TestType==eTestPlaneEmbed)
    {
        // Let the moving frames outside the domain be of magnitude zero.
        int index;
        for (int i = m_EmbededPlane; i < m_fields[0]->GetExpSize(); ++i)
        {
            for (int j = 0; j < m_fields[0]->GetTotPoints(i); ++j)
            {
                index = m_fields[0]->GetPhys_Offset(i) + j;
                fields[0][index] = 0.0; 
            }
        }
    }


    for (i = 0; i < 1; ++i)
    {
        m_fields[m_intVariables[i]]->SetPhys(fields[i]);
        m_fields[m_intVariables[i]]->SetPhysState(true);
    }

    Array<OneD, NekDouble> uexact(nq);
    TestPlaneEmbedProblem(m_time, m_varcoeff, uexact);
    Checkpoint_Output_Error(nchk-1, m_fields[0]->GetPhys(), uexact);
} // namespace Nektar

void MMFDiffusion::v_GenerateSummary(SolverUtils::SummaryList &s)
{
    MMFSystem::v_GenerateSummary(s);
    SolverUtils::AddSummaryItem(s, "TestType", TestTypeMap[m_TestType]);
    SolverUtils::AddSummaryItem(s, "frequency", m_frequency);
    SolverUtils::AddSummaryItem(s, "AniStrength", m_AniStrength);
    SolverUtils::AddSummaryItem(s, "Helmtau", m_Helmtau);
    SolverUtils::AddSummaryItem(s, "EmbededPlane", m_EmbededPlane);

    SolverUtils::AddSummaryItem(s, "epsilon0", m_epsilon[0]);
    SolverUtils::AddSummaryItem(s, "epsilon1", m_epsilon[1]);
    SolverUtils::AddSummaryItem(s, "epsilon2", m_epsilon[2]);
    SolverUtils::AddSummaryItem(s, "d00", m_d00);
    SolverUtils::AddSummaryItem(s, "d11", m_d11);
    SolverUtils::AddSummaryItem(s, "d22", m_d22);
    if (m_TestType == eTestLinearSphere)
    {
        SolverUtils::AddSummaryItem(s, "epsilon for u", m_epsu[0]);
        SolverUtils::AddSummaryItem(s, "epsilon for v", m_epsu[1]);
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
