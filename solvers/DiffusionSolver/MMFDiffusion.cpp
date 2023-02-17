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

// void MMFDiffusion::v_InitObject(bool DeclareFields)
// {
//     UnsteadySystem::v_InitObject(DeclareFields);

//     int nq = GetTotPoints();
//     int nvar = m_fields.size();

//     // Helmsolver parameter
//     m_session->LoadParameter("Helmtau", m_Helmtau, 1.0);

//     // Define ProblemType
//     if (m_session->DefinesSolverInfo("TESTTYPE"))
//     {
//         std::string TestTypeStr;
//         TestTypeStr = m_session->GetSolverInfo("TESTTYPE");
//         for (int i = 0; i < (int)SIZE_TestType; ++i)
//         {
//             if (boost::iequals(TestTypeMap[i], TestTypeStr))
//             {
//                 m_TestType = (TestType)i;
//                 break;
//             }
//         }
//     }
//     else
//     {
//         m_TestType = (TestType)0;
//     }

//     // Define SovlerSchemeType
//     if (m_session->DefinesSolverInfo("SolverSchemeType"))
//     {
//         std::string SolverSchemeTypeStr;
//         SolverSchemeTypeStr = m_session->GetSolverInfo("SolverSchemeType");
//         for (int i = 0; i < (int)SIZE_SolverSchemeType; ++i)
//         {
//             if (boost::iequals(SolverSchemeTypeMap[i], SolverSchemeTypeStr))
//             {
//                 m_SolverSchemeType = (SolverSchemeType)i;
//                 break;
//             }
//         }
//     }
//     else
//     {
//         m_SolverSchemeType = (SolverSchemeType)0;
//     }

//     // TimeMap ?
//     if (m_session->DefinesSolverInfo("TimeMapType"))
//     {
//         std::string TIMEMAPTYPEStr;
//         TIMEMAPTYPEStr = m_session->GetSolverInfo("TimeMapType");
//         for (int i = 0; i < (int)SIZE_TimeMapType; ++i)
//         {
//             if (boost::iequals(TimeMapTypeMap[i], TIMEMAPTYPEStr))
//             {
//                 m_TimeMap = (TimeMapType)i;
//                 break;
//             }
//         }
//     }
//     else
//     {
//         m_TimeMap = (TimeMapType)0;
//     }

//     // Diffusivity coefficient for e^j
//     m_epsilon = Array<OneD, NekDouble>(m_expdim);
//     m_session->LoadParameter("epsilon0", m_epsilon[0], 1.0);
//     m_session->LoadParameter("epsilon1", m_epsilon[1], 1.0);
//     m_session->LoadParameter("epsilon2", m_epsilon[2], 1.0);

//     // Diffusivity coefficient for u^j
//     m_epsu = Array<OneD, NekDouble>(nvar + 1);
//     m_session->LoadParameter("epsu0", m_epsu[0], 1.0);
//     m_session->LoadParameter("epsu1", m_epsu[1], 1.0);

//     m_session->LoadParameter("Diffbeta", m_Diffbeta, 0.5);
//     m_session->LoadParameter("Diffeta", m_Diffeta, 100.0);
//     m_session->LoadParameter("Diffhe", m_Diffhe, 0.5);

//     // Derive AnisotropyStrength.
//     Array<OneD, Array<OneD, NekDouble>> AniStrength(m_expdim);
//     for (int j = 0; j < m_expdim; ++j)
//     {
//         AniStrength[j] = Array<OneD, NekDouble>(nq, 1.0);
//     }

//     MMFSystem::MMFInitObject(AniStrength);

//     // if (m_explicitDiffusion)
//     // {
//     //     m_ode.DefineImplicitSolve(&MMFDiffusion::DoNullSolve, this);
//     //     m_ode.DefineProjection(&MMFDiffusion::DoOdeProjection, this);
//     // }

//     // else
//     // {
//         // Create varcoeff for Helmsolver
//         // m_varcoeff = ComputeVarCoeff2D(m_movingframes);

//     // StdRegions::VarCoeffType MMFCoeffs[15] = {
//     //     StdRegions::eVarCoeffMF1x,   StdRegions::eVarCoeffMF1y,
//     //     StdRegions::eVarCoeffMF1z,   StdRegions::eVarCoeffMF1Div,
//     //     StdRegions::eVarCoeffMF1Mag, StdRegions::eVarCoeffMF2x,
//     //     StdRegions::eVarCoeffMF2y,   StdRegions::eVarCoeffMF2z,
//     //     StdRegions::eVarCoeffMF2Div, StdRegions::eVarCoeffMF2Mag,
//     //     StdRegions::eVarCoeffMF3x,   StdRegions::eVarCoeffMF3y,
//     //     StdRegions::eVarCoeffMF3z,   StdRegions::eVarCoeffMF3Div,
//     //     StdRegions::eVarCoeffMF3Mag};

//     //   int indx;
//     //     for (int k = 0; k < m_mfdim; ++k)
//     //     {
//     //         // For Moving Frames
//     //         indx = 5 * k;

//     //         std::cout << "indx = " << indx << std::endl;
//     //                 std::cout << "MF = ( " <<
//     RootMeanSquare(m_varcoeff[MMFCoeffs[indx + 0]]) << " , "
//     //         << RootMeanSquare(m_varcoeff[MMFCoeffs[indx + 1]]) << " , "
//     //         << RootMeanSquare(m_varcoeff[MMFCoeffs[indx + 2]]) << " ) " <<
//     std::endl;

//     //         std::cout << "DivMF = " <<
//     RootMeanSquare(m_varcoeff[MMFCoeffs[indx + 3]]) << std::endl;

//     //         std::cout << "emag = " <<
//     RootMeanSquare(m_varcoeff[MMFCoeffs[indx + 4]]) << std::endl;
//     //     }

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

//         // std::cout << "indx = " << indx << std::endl;

//         for (int j = 0; j < m_spacedim; ++j)
//         {
//             m_varcoeff[MMFCoeffs[indx + j]] = Array<OneD, NekDouble>(nq,
//             0.0); Vmath::Vcopy(nq, &m_movingframes[k][j * nq], 1,
//                          &m_varcoeff[MMFCoeffs[indx + j]][0], 1);
//         }

//         // std::cout << "MF = ( " <<
//         RootMeanSquare(m_varDiffcoeff[MMFCoeffs[indx + 0]]) << " , "
//         // << RootMeanSquare(m_varDiffcoeff[MMFCoeffs[indx + 1]]) << " , "
//         // << RootMeanSquare(m_varDiffcoeff[MMFCoeffs[indx + 2]]) << " ) " <<
//         std::endl;

//         // m_DivMF
//         m_varcoeff[MMFCoeffs[indx + 3]] = Array<OneD, NekDouble>(nq, 0.0);
//         Vmath::Vcopy(nq, &m_DivMF[k][0], 1, &m_varcoeff[MMFCoeffs[indx +
//         3]][0],
//                      1);

//         // std::cout << "DivMF = " <<
//         RootMeanSquare(m_varDiffcoeff[MMFCoeffs[indx + 3]]) << std::endl;

//         // \| e^k \|
//         m_varcoeff[MMFCoeffs[indx + 4]] = Array<OneD, NekDouble>(nq, 0.0);
//         tmp                             = Array<OneD, NekDouble>(nq, 0.0);
//         for (int i = 0; i < m_spacedim; ++i)
//         {
//             Vmath::Vvtvp(nq, &m_movingframes[k][i * nq], 1,
//                          &m_movingframes[k][i * nq], 1, &tmp[0], 1, &tmp[0],
//                          1);
//         }
//         Vmath::Vcopy(nq, &tmp[0], 1, &m_varcoeff[MMFCoeffs[indx + 4]][0], 1);

//         // std::cout << "emag = " <<
//         RootMeanSquare(m_varDiffcoeff[MMFCoeffs[indx + 4]]) << std::endl <<
//         std::endl;
//     }

//         m_ode.DefineImplicitSolve(&MMFDiffusion::DoImplicitSolve, this);
//             m_ode.DefineOdeRhs(&MMFDiffusion::DoOdeRhs, this);
//     }

void MMFDiffusion::v_InitObject(bool DeclareFields)
{
    UnsteadySystem::v_InitObject(DeclareFields);

    int nq    = m_fields[0]->GetNpoints();
    int nvar  = m_fields.size();
    int MFdim = 3;

    // Diffusivity coefficient for e^j
    m_epsilon = Array<OneD, NekDouble>(MFdim);
    m_session->LoadParameter("epsilon0", m_epsilon[0], 1.0);
    m_session->LoadParameter("epsilon1", m_epsilon[1], 1.0);
    m_session->LoadParameter("epsilon2", m_epsilon[2], 1.0);

    // Diffusivity coefficient for u^j
    m_epsu = Array<OneD, NekDouble>(nvar + 1);
    m_session->LoadParameter("epsu0", m_epsu[0], 1.0);
    m_session->LoadParameter("epsu1", m_epsu[1], 1.0);

    m_session->LoadParameter("InitPtx", m_InitPtx, 0.0);
    m_session->LoadParameter("InitPty", m_InitPty, 0.0);
    m_session->LoadParameter("InitPtz", m_InitPtz, 0.0);

    int shapedim = m_fields[0]->GetShapeDimension();
    Array<OneD, Array<OneD, NekDouble>> Anisotropy(shapedim);
    for (int j = 0; j < shapedim; ++j)
    {
        Anisotropy[j] = Array<OneD, NekDouble>(nq, 1.0);
        Vmath::Fill(nq, sqrt(m_epsilon[j]), &Anisotropy[j][0], 1);
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

    ComputeVarCoeff2D(m_movingframes,m_varcoeff);

    if (!m_explicitDiffusion)
    {
        m_ode.DefineProjection(&MMFDiffusion::DoOdeProjection, this);
        m_ode.DefineOdeRhs(&MMFDiffusion::DoOdeRhs, this);
    }

    else
    {
        m_ode.DefineProjection(&MMFDiffusion::DoOdeProjection, this);
        m_ode.DefineOdeRhs(&MMFDiffusion::DoOdeRhs, this);
        m_ode.DefineImplicitSolve(&MMFDiffusion::DoImplicitSolve, this);
    }

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
    int nvariables = inarray.size();
    int nq         = m_fields[0]->GetNpoints();

    StdRegions::ConstFactorMap factors;
    factors[StdRegions::eFactorTau] = 1.0;

    Array<OneD, Array<OneD, NekDouble>> F(nvariables);
    factors[StdRegions::eFactorLambda] = 1.0 / lambda;
    F[0] = Array<OneD, NekDouble>(nq * nvariables);
    for (int n = 1; n < nvariables; ++n)
    {
        F[n] = F[n - 1] + nq;
    }

    // We solve ( \nabla^2 - HHlambda ) Y[i] = rhs [i]
    // inarray = input: \hat{rhs} -> output: \hat{Y}
    // outarray = output: nabla^2 \hat{Y}
    // where \hat = modal coeffs
    SetBoundaryConditions(time);

    for (int i = 0; i < nvariables; ++i)
    {
        factors[StdRegions::eFactorLambda] = 1.0 / lambda / m_epsu[i];

        // Multiply 1.0/timestep
        Vmath::Smul(nq, -factors[StdRegions::eFactorLambda], inarray[i], 1,
                    F[i], 1);
        m_fields[i]->HelmSolve(F[i], m_fields[i]->UpdateCoeffs(), factors,
                               m_varDiffcoeff);
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

    switch (m_TestType)
    {
        case eTestPlane:
        {
            Array<OneD, NekDouble> x(nq);
            Array<OneD, NekDouble> y(nq);
            Array<OneD, NekDouble> z(nq);

            m_fields[0]->GetCoords(x, y, z);

            for (int k = 0; k < nq; k++)
            {
                outarray[0][k] = (m_epsilon[0] + m_epsilon[1] - 1.0) * m_pi *
                                 m_pi * exp(-1.0 * m_pi * m_pi * time) *
                                 sin(m_pi * x[k]) * cos(m_pi * y[k]);
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
                    (m_epsilon[0] + m_epsilon[1] + m_epsilon[2] - 1.0) * m_pi *
                    m_pi * exp(-1.0 * m_pi * m_pi * time) * sin(m_pi * x[k]) *
                    sin(m_pi * y[k]) * sin(m_pi * z[k]);
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
        int nq = m_fields[0]->GetNpoints();

        // Laplacian only to the first variable
        Array<OneD, NekDouble> Laplacian(nq);
        WeakDGMMFDiffusion(0, inarray[0], Laplacian, time);

        Vmath::Vadd(nq, &Laplacian[0], 1, &outarray[0][0], 1, &outarray[0][0], 1);
    }
}

void MMFDiffusion::v_SetInitialConditions(NekDouble initialtime,
                                          bool dumpInitialConditions,
                                          const int domain)
{
    boost::ignore_unused(domain);

    int nq = GetTotPoints();

    EquationSystem::v_SetInitialConditions(initialtime, false);

    switch (m_TestType)
    {
        case eTestPlane:
        {
            Array<OneD, NekDouble> u(nq);

            TestPlaneProblem(initialtime, u);
            m_fields[0]->SetPhys(u);
        }
        break;

        case eTestCube:
        {
            Array<OneD, NekDouble> u(nq);

            TestCubeProblem(initialtime, u);
            m_fields[0]->SetPhys(u);
            /*for (int k=0; k<nq; ++k)
            {
                //for (int j=0; j<m_spacedim; ++j)
                //{
                cout << "_varcoeff" << u[k] <<endl;
                // }
            }*/
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

            // case eFHNStandard:
            // case eFHNRogers:
            // case eFHNAlievPanf:
            // {
            //     Array<OneD, NekDouble> Zero(nq, 0.0);
            //     m_fields[0]->SetPhys(PlanePhiWave());
            //     m_fields[1]->SetPhys(Zero);
            // }
            // break;

        default:
        {
            EquationSystem::v_SetInitialConditions(initialtime, false);
        }
        break;
    }

    std::cout << "Initial: max um = "
              << Vmath::Vmax(nq, m_fields[0]->GetPhys(), 1) << std::endl;

    if (dumpInitialConditions)
    {
        std::string outname = m_sessionName + "_initial.chk";
        WriteFld(outname);
    }
}

void MMFDiffusion::TestPlaneProblem(const NekDouble time,
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
        outfield[k] = exp(-1.0 * m_pi * m_pi * time) * sin(m_pi * x[k]) *
                      cos(m_pi * y[k]);
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
        outfield[k] = exp(-1.0 * m_pi * m_pi * time) * sin(m_pi * x[k]) *
                      sin(m_pi * y[k]) * sin(m_pi * z[k]);
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

                // NekDouble xc = 0.5*(Vmath::Vmax(nq, x, 1) + Vmath::Vmin(nq,
                // x, 1));

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
        case eTestPlane:
        {
            TestPlaneProblem(time, outfield);
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

void MMFDiffusion::ComputeVarCoeff2D(
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
    for (int k = 0; k < m_expdim; ++k)
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
        // ComputeDivMF(eCovariant, movingframes, DivMF);

        ComputeEuclideanDivMF(movingframes, DivMF);

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

    std::cout << "m_varcoeff = " << RootMeanSquare(varcoeff[MMFCoeffs[0]])
              << " , " << RootMeanSquare(varcoeff[MMFCoeffs[1]]) << " , "
              << RootMeanSquare(varcoeff[MMFCoeffs[2]]) << " , "
              << RootMeanSquare(varcoeff[MMFCoeffs[3]]) << " , "
              << RootMeanSquare(varcoeff[MMFCoeffs[4]]) << std::endl;

    std::cout << " ::::: 2D Varcoeff is Successfully Created ::::: "
              << std::endl;
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

void MMFDiffusion::v_GenerateSummary(SolverUtils::SummaryList &s)
{
    MMFSystem::v_GenerateSummary(s);
    SolverUtils::AddSummaryItem(s, "TestType", TestTypeMap[m_TestType]);
    SolverUtils::AddSummaryItem(s, "epsilon0", m_epsilon[0]);
    SolverUtils::AddSummaryItem(s, "epsilon1", m_epsilon[1]);
    SolverUtils::AddSummaryItem(s, "epsilon2", m_epsilon[2]);
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
