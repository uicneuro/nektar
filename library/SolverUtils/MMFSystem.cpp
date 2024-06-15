///////////////////////////////////////////////////////////////////////////////
//
// File: MMFSystem.cpp
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
// Description: Base class for MMF systems.
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>
#include <SolverUtils/MMFSystem.h>

namespace Nektar
{
namespace SolverUtils
{

MMFSystem::MMFSystem(const LibUtilities::SessionReaderSharedPtr &pSession,
                     const SpatialDomains::MeshGraphSharedPtr &pGraph)
    : UnsteadySystem(pSession, pGraph)
{
}

MMFSystem::~MMFSystem()
{
}

void MMFSystem::MMFInitObject(
    const Array<OneD, const Array<OneD, NekDouble>> &AniStrength,
    const Array<OneD, const NekDouble> &AniDirection)
{
    std::cout << std::endl;
    std::cout << "MMFInitObejct Starts: ==============================================" << std::endl;

    m_pi       = 3.14159265358979323846;
    m_shapedim = m_expdim;
    m_mfdim    = 3;

    ASSERTL0(m_spacedim == 3, "Space Dimension should be 3");

    // Define MMFOrderType
    if (m_session->DefinesSolverInfo("DerivType"))
    {
        std::string DerivTypeStr = m_session->GetSolverInfo("DerivType");
        for (int i = 0; i < (int)SIZE_DerivType; ++i)
        {
            if (DerivTypeMap[i] == DerivTypeStr)
            {
                m_DerivType = (DerivType)i;
                break;
            }
        }
    }
    else
    {
        m_DerivType = (DerivType)0;
    }           

    // Define SurfaceType
    if (m_session->DefinesSolverInfo("SURFACETYPE"))
    {
        std::string SurfaceTypeStr;
        SurfaceTypeStr = m_session->GetSolverInfo("SURFACETYPE");
        for (int i = 0; i < (int)SIZE_SurfaceType; ++i)
        {
            if (SurfaceTypeMap[i] == SurfaceTypeStr)
            {
                m_surfaceType = (SurfaceType)i;
                break;
            }
        }
    }
    else
    {
        m_surfaceType = (SurfaceType)0;
    }

    if (m_surfaceType == SolverUtils::eEllipsoid)
    {
        m_session->LoadParameter("Radx", m_Radx, 1.0);
        m_session->LoadParameter("Rady", m_Rady, 1.0);
        m_session->LoadParameter("Radz", m_Radz, 1.0);
    }

    // Define Gradient Location Type
    if (m_session->DefinesSolverInfo("GRADLOCTYPE"))
    {
        std::string GradLocTypeStr;
        GradLocTypeStr = m_session->GetSolverInfo("GRADLOCTYPE");
        for (int i = 0; i < (int)SIZE_GradLocType; ++i)
        {
            if (boost::iequals(GradLocTypeMap[i], GradLocTypeStr))
            {
                m_GradLocType = (GradLocType)i;
                break;
            }
        }
    }
    else
    {
        m_GradLocType = (GradLocType)0;
    }

    // if discontinuous Galerkin determine numerical flux to use
    for (int i = 0; i < (int)SIZE_UpwindType; ++i)
    {
        bool match;
        m_session->MatchSolverInfo("UPWINDTYPE", UpwindTypeMap[i], match,
                                   false);
        if (match)
        {
            m_upwindType = (UpwindType)i;
            break;
        }
    }

    m_session->LoadParameter("SphereExactRadius", m_SphereExactRadius, 1.0);

    m_session->LoadParameter("Initx", m_Initx, 0.0);
    m_session->LoadParameter("Inity", m_Inity, 0.0);
    m_session->LoadParameter("Initz", m_Initz, m_SphereExactRadius);

    m_session->LoadParameter("ROIx", m_ROIx, m_Initx);
    m_session->LoadParameter("ROIy", m_ROIy, m_Inity);
    m_session->LoadParameter("ROIz", m_ROIz, m_Initz);

    m_session->LoadParameter("GaussianTimeMap", m_GaussianTimeMap, 0);
    m_session->LoadParameter("GaussianRadius", m_GaussianRadius, 1.0);

    m_session->LoadParameter("AdaptNewFramesTol", m_AdaptNewFramesTol, 1.0);
    m_session->LoadParameter("VelActivateTol", m_VelActivationTol, 0.1);
    m_session->LoadParameter("uTol", m_uTol, 0.1);
    m_session->LoadParameter("NoAlignInitRadius", m_NoAlignInitRadius, 5.0);

    m_session->LoadParameter("LDGc11", m_LDGC11, 1.0);

    //  if Velmag < ActivationTol, the new alignment is not activated
    m_session->LoadParameter("c121", m_c121, 0);
    m_session->LoadParameter("c122", m_c122, 0);
    m_session->LoadParameter("c123", m_c123, 0);

    // Factor for Numerical Flux
    // m_session->LoadParameter("alpha", m_alpha, 1.0);

    // Factor for Numerical Flux
    m_session->LoadParameter("Incfreq", m_Incfreq, 1.0);

    // SmoothFactor
    m_session->LoadParameter("SFinit", m_SFinit, 0.0);

    int nq = m_fields[0]->GetNpoints();

    // if 1D, computed trajectory length from the left botom to right top.
    if (m_expdim == 1)
    {
        Array<OneD, NekDouble> x0(nq);
        Array<OneD, NekDouble> x1(nq);
        Array<OneD, NekDouble> x2(nq);

        m_fields[0]->GetCoords(x0, x1, x2);

        NekDouble xdis, ydis, zdis, totlength = 0.0;
        m_seglength = Array<OneD, NekDouble>(nq, 0.0);
        for (int i = 1; i < nq; ++i)
        {
            xdis = x0[i] - x0[i - 1];
            ydis = x1[i] - x1[i - 1];
            zdis = x2[i] - x2[i - 1];
            totlength += sqrt(xdis * xdis + ydis * ydis + zdis * zdis);

            m_seglength[i] = totlength;
        }
    }

    switch (m_surfaceType)
    {
        case SolverUtils::eSphere:
        {
            // Construct Spherical moving frames
            ConstructSphericalMF(m_sphereMF, m_MMFActivation);
        }
        break;

        case SolverUtils::ePseudosphere:
        {
            // Construct Spherical moving frames
            ConstructPseudosphericalMF(m_pseudosphereMF, m_MMFActivation);
        }
        break;

        case SolverUtils::eEllipsoid:
        {
            // Construct Spherical moving frames
            ConstructEllipticalMF(m_Radx, m_Rady, m_Radz, m_sphereMF,
                                  m_MMFActivation);
        }
        break;

        case SolverUtils::ePolar:
        {
            ConstructPolarMF(m_polarMF, m_MMFActivation);
        }
        break;

        case SolverUtils::eTorus:
        case SolverUtils::ePlane:
        default:
        {
            m_MMFActivation = Array<OneD, int>(nq, 1);

            int cnt = Vmath::Vsum(nq, m_MMFActivation, 1);
            std::cout << "MMFActivation = " << cnt
                    << " /  " << nq << " ( " << cnt/nq*100.0 << " % ) " << std::endl;

        }
        break;
    }    

    // SetUpMovingFrames: To generate m_movingframes
    if(AniDirection == NullNekDouble1DArray)
    {
        std::string MMFdirStr;
        m_session->LoadSolverInfo("MMFDir", MMFdirStr, "LOCAL");
        m_MMFdir = FindMMFdir(MMFdirStr);

        SetUpMovingFrames(m_MMFdir, AniStrength, m_movingframes); 
    }

    else
    {
        // GenerateMFbyAniDirection(AniDreiction, m_movingframes);
    }

    switch (m_projectionType)
    {
        case MultiRegions::eDiscontinuous:
        {
            ComputeMFtrace(m_movingframes, m_MFtraceFwd, m_MFtraceBwd);

            // Check Movingframes and surfraceNormal
            // Get: m_ncdotMFFwd,m_ncdotMFBwd,m_nperpcdotMFFwd,m_nperpcdotMFBwd
            if(m_expdim>1)
            {
                ComputencdotMF(m_movingframes, m_ncdotMFFwd, m_ncdotMFBwd, 1);
            }

            else
            {
                int nTracePointsTot = GetTraceNpoints();

                m_ncdotMFFwd = Array<OneD, Array<OneD, NekDouble>>(m_mfdim);
                m_ncdotMFBwd = Array<OneD, Array<OneD, NekDouble>>(m_mfdim);            
                
                m_ncdotMFFwd[0] = Array<OneD, NekDouble>(nTracePointsTot, AniStrength[0]);
                m_ncdotMFBwd[0] = Array<OneD, NekDouble>(nTracePointsTot, AniStrength[0]);
                for (int j = 1; j < m_mfdim; ++j)
                {
                    m_ncdotMFFwd[j] = Array<OneD, NekDouble>(nTracePointsTot, 0.0);
                    m_ncdotMFBwd[j] = Array<OneD, NekDouble>(nTracePointsTot, 0.0);
                }
            }

            ComputenperpcdotMF(m_movingframes, m_nperpcdotMFFwd, m_nperpcdotMFBwd);

            ComputeDivMF(m_DerivType, m_movingframes, m_DivMF);
            ComputeCurlMF(m_DerivType, m_movingframes, m_CurlMF);
            break;
        }

        case MultiRegions::eGalerkin:
        case MultiRegions::eMixed_CG_Discontinuous:
        default:
        break;
    }

    std::cout << "MMFInitObject is done ==============================================" << std::endl;
    std::cout << std::endl;
    // Connection 2-form
    // if(m_expdim>1)
    // {
    //     Compute2DConnection1form(m_movingframes, m_MFConnection);

    //     // Check the Curvature 2-form of the aligned moving frames
    //     Compute2DCurvatureForm(m_movingframes, m_MFConnection, m_MFCurvature);
    // }
}

// Check Connection and Curvature for Spherical coordinate system
void MMFSystem::CheckSphereConnection()
{
    int nq = m_fields[0]->GetNpoints();

    std::cout << "=========== Starting TestSpherical ===========" << std::endl;

    std::cout << "Counting Activation = "
              << 100.0 * CountActivated(m_MMFActivation) / nq << " %"
              << std::endl;

    // Test connections of moving frames
    TestSphericalConnection1form(m_sphereMF, m_MMFActivation);

    Array<OneD, Array<OneD, NekDouble>> MF1st(m_mfdim);
    for (int k = 0; k < m_mfdim; ++k)
    {
        MF1st[k] = Array<OneD, NekDouble>(m_spacedim * nq);
    }

    Array<OneD, NekDouble> x0(nq);
    Array<OneD, NekDouble> x1(nq);
    Array<OneD, NekDouble> x2(nq);

    m_fields[0]->GetCoords(x0, x1, x2);

    NekDouble xp, yp, zp;
    NekDouble sin_varphi, cos_varphi, sin_theta, cos_theta;
    Array<OneD, NekDouble> MFvec(m_spacedim * nq, 0.0);

    // Type 2: Rossby-Harwitz type vector
    MFvec = ComputeRossyHaurtizVelocity();
    for (int i = 0; i < nq; ++i)
    {
        xp = x0[i];
        yp = x1[i];
        zp = x2[i];

        // x = r \sin \theta \cos \varphi
        // y = r \sin \theta \sin \varphi
        // z = r \cos \theta
        CartesianToNewSpherical(xp, yp, zp, sin_varphi, cos_varphi, sin_theta,
                                cos_theta);

        // Type 1: Vector along the spherical coordinate axis of \theta, but
        // factored by \sin \theta MFvec[i]          = sin_theta *
        // m_sphereMF[0][i]; MFvec[i + nq]     = sin_theta * m_sphereMF[0][i +
        // nq]; MFvec[i + 2 * nq] = sin_theta * m_sphereMF[0][i + 2 * nq];

        if (m_MMFActivation[i])
        {
            for (int k = 0; k < m_mfdim; ++k)
            {
                MF1st[k][i]          = m_sphereMF[k][i];
                MF1st[k][i + nq]     = m_sphereMF[k][i + nq];
                MF1st[k][i + 2 * nq] = m_sphereMF[k][i + 2 * nq];
            }
        }

        else
        {
            for (int k = 0; k < m_mfdim; ++k)
            {
                MF1st[k][i]          = m_movingframes[k][i];
                MF1st[k][i + nq]     = m_movingframes[k][i + nq];
                MF1st[k][i + 2 * nq] = m_movingframes[k][i + 2 * nq];
            }
        }
    }

    // Check whether Moving frames are everywhere and orthonormal
    // CheckMovingFrames(MF1st);
    ComputeMFtrace(MF1st, m_MFtraceFwd, m_MFtraceBwd);
    ComputencdotMF(MF1st, m_ncdotMFFwd, m_ncdotMFBwd);

    std::cout << CountActivated(m_MMFActivation) << " / " << nq
              << " spherical MF are constructed " << std::endl;
}

void MMFSystem::CheckPolarConnection()
{
    int nq  = m_fields[0]->GetNpoints();
    int cnt = 0;

    Array<OneD, Array<OneD, NekDouble>> MF1st(m_mfdim);
    for (int k = 0; k < m_mfdim; ++k)
    {
        MF1st[k] = Array<OneD, NekDouble>(m_spacedim * nq);
    }

    Array<OneD, NekDouble> MFvec(m_spacedim * nq, 0.0);
    for (int i = 0; i < nq; ++i)
    {

        if (m_MMFActivation[i])
        {
            for (int k = 0; k < m_mfdim; ++k)
            {
                MF1st[k][i]          = m_polarMF[k][i];
                MF1st[k][i + nq]     = m_polarMF[k][i + nq];
                MF1st[k][i + 2 * nq] = m_polarMF[k][i + 2 * nq];
            }

            MFvec[i]          = m_polarMF[0][i];
            MFvec[i + nq]     = m_polarMF[0][i + nq];
            MFvec[i + 2 * nq] = m_polarMF[0][i + 2 * nq];
            cnt++;
        }

        else
        {
            for (int k = 0; k < m_mfdim; ++k)
            {
                MF1st[k][i]          = m_movingframes[k][i];
                MF1st[k][i + nq]     = m_movingframes[k][i + nq];
                MF1st[k][i + 2 * nq] = m_movingframes[k][i + 2 * nq];
            }
        }
    }

    std::cout << cnt << " / " << nq << " polar MF are constructed "
              << std::endl;

    TestDivMF(MF1st, m_MMFActivation);

    Array<OneD, NekDouble> irrotational;
    Array<OneD, NekDouble> incompressible;
    Array<OneD, NekDouble> harmonic;

    ComputeMFHHD(m_MMFActivation, MFvec, irrotational, incompressible,
                 harmonic);
    PlotHHD(MFvec, irrotational, incompressible, harmonic, 0);

    std::cout << "=========== Ending TestPolar ===========" << std::endl;
}

void MMFSystem::ConstructSphericalMF(
    Array<OneD, Array<OneD, NekDouble>> &SphereMF,
    Array<OneD, int> &SphereMFActivate)
{
    int nq = m_fields[0]->GetNpoints();

    SphereMF = Array<OneD, Array<OneD, NekDouble>>(m_mfdim);
    for (int j = 0; j < m_mfdim; ++j)
    {
        SphereMF[j] = Array<OneD, NekDouble>(m_spacedim * nq, 0.0);
    }

    SphereMFActivate = Array<OneD, int>(nq, 0);

    Array<OneD, NekDouble> x0(nq);
    Array<OneD, NekDouble> x1(nq);
    Array<OneD, NekDouble> x2(nq);

    m_fields[0]->GetCoords(x0, x1, x2);

    NekDouble rad, xp, yp, zp;
    NekDouble sin_varphi, cos_varphi, sin_theta, cos_theta;
    NekDouble Tol = 0.1;
    Array<OneD, NekDouble> radvec;

    int indexj, indexk;
    for (int i = 0; i < m_fields[0]->GetExpSize(); ++i)
    {
        radvec = Array<OneD, NekDouble>(m_fields[0]->GetTotPoints(i));

        // j iteration
        for (int j = 0; j < m_fields[0]->GetTotPoints(i); ++j)
        {
            indexj = m_fields[0]->GetPhys_Offset(i) + j;

            xp = x0[indexj];
            yp = x1[indexj];
            zp = x2[indexj];

            radvec[j] = sqrt((xp - m_Initx) * (xp - m_Initx) +
                             (yp - m_Inity) * (yp - m_Inity) +
                             (zp - m_Initz) * (zp - m_Initz));

            rad = sqrt((xp + m_Initx) * (xp + m_Initx) +
                       (yp + m_Inity) * (yp + m_Inity) +
                       (zp + m_Initz) * (zp + m_Initz));

            if (rad < radvec[j])
            {
                radvec[j] = rad;
            }
        }

        // x = r \sin \theta \cos \varphi
        // y = r \sin \theta \sin \varphi
        // z = r \cos \theta
        // \theta \in [0, \pi], \varphi = [0, 2 \pi]
        // if (fabs(zpavg) < m_SphereExactRadius * (1.0 - Tol))

        // k iteration
        for (int k = 0; k < m_fields[0]->GetTotPoints(i); ++k)
        {
            indexk = m_fields[0]->GetPhys_Offset(i) + k;

            xp = x0[indexk];
            yp = x1[indexk];
            zp = x2[indexk];

            CartesianToNewSpherical(xp, yp, zp, sin_varphi, cos_varphi,
                                    sin_theta, cos_theta);

            if (Vmath::Vmin(m_fields[0]->GetTotPoints(i), radvec, 1) > Tol)
            {
                SphereMFActivate[indexk] = 1;
            }

            // Theta - direction
            SphereMF[0][indexk]          = cos_theta * cos_varphi;
            SphereMF[0][indexk + nq]     = cos_theta * sin_varphi;
            SphereMF[0][indexk + 2 * nq] = -sin_theta;

            // Phi - direction
            SphereMF[1][indexk]          = -1.0 * sin_varphi;
            SphereMF[1][indexk + nq]     = cos_varphi;
            SphereMF[1][indexk + 2 * nq] = 0.0;

            // R - direction
            SphereMF[2][indexk]          = sin_theta * cos_varphi;
            SphereMF[2][indexk + nq]     = sin_theta * sin_varphi;
            SphereMF[2][indexk + 2 * nq] = cos_theta;
        }
    }

    int cnt = Vmath::Vsum(nq, SphereMFActivate, 1);
    std::cout << "MMFActivation = " << cnt
              << " /  " << nq << " ( " << cnt/nq*100.0 << " % ) " << std::endl;
}

void MMFSystem::ConstructPseudosphericalMF(
    Array<OneD, Array<OneD, NekDouble>> &PsedoSphereMF,
    Array<OneD, int> &PseudosphereMFActivate)
{
    int nq = m_fields[0]->GetNpoints();

    PsedoSphereMF = Array<OneD, Array<OneD, NekDouble>>(m_mfdim);
    for (int j = 0; j < m_mfdim; ++j)
    {
        PsedoSphereMF[j] = Array<OneD, NekDouble>(m_spacedim * nq, 0.0);
    }

    PseudosphereMFActivate = Array<OneD, int>(nq, 1);

    Array<OneD, NekDouble> x0(nq);
    Array<OneD, NekDouble> x1(nq);
    Array<OneD, NekDouble> x2(nq);

    m_fields[0]->GetCoords(x0, x1, x2);

    NekDouble xp, yp, zp;
    NekDouble e1x, e1y, e1z;
    NekDouble sin_varphi, cos_varphi, theta, sech_theta, tanh_theta;
    NekDouble Tol = 0.001;
    int flag;

    int indexj, indexk;
    for (int i = 0; i < m_fields[0]->GetExpSize(); ++i)
    {
        flag = 0;
        
        for (int j = 0; j < m_fields[0]->GetTotPoints(i); ++j)
        {
            indexj = m_fields[0]->GetPhys_Offset(i) + j;

            xp = x0[indexj];
            yp = x1[indexj];
            zp = x2[indexj];

            if ((fabs(zp - 0.0) < Tol) || (fabs(zp - 1.0) < Tol))
            {
                flag = 1;
            }
        }

        for (int k = 0; k < m_fields[0]->GetTotPoints(i); ++k)
        {
            indexk = m_fields[0]->GetPhys_Offset(i) + k;

            xp = x0[indexk];
            yp = x1[indexk];
            zp = x2[indexk];

            CartesianToPseudospherical(xp, yp, zp, sin_varphi, cos_varphi,
                                       theta, sech_theta, tanh_theta);

            // Theta - direction
            e1x = -sech_theta * cos_varphi;
            e1y = -sech_theta * sin_varphi;
            e1z = tanh_theta;

            PsedoSphereMF[0][indexk]          = e1x;
            PsedoSphereMF[0][indexk + nq]     = e1y;
            PsedoSphereMF[0][indexk + 2 * nq] = e1z;

            // Phi - direction
            PsedoSphereMF[1][indexk]          = -1.0 * sin_varphi;
            PsedoSphereMF[1][indexk + nq]     = cos_varphi;
            PsedoSphereMF[1][indexk + 2 * nq] = 0.0;

            // R - direction
            PsedoSphereMF[2][indexk]          = sech_theta * cos_varphi;
            PsedoSphereMF[2][indexk + nq]     = sech_theta * sin_varphi;
            PsedoSphereMF[2][indexk + 2 * nq] = theta - tanh_theta;

            if (flag)
            {
                PseudosphereMFActivate[indexk] = 0;
            }
        }
    }

    std::cout << "PseudoSphere MMFActivation = "
              << Vmath::Vsum(nq, PseudosphereMFActivate, 1) << " / " << nq
              << std::endl;
}

void MMFSystem::ConstructEllipticalMF(
    const int Radx, const int Rady, const int Radz,
    Array<OneD, Array<OneD, NekDouble>> &EllipticalMF,
    Array<OneD, int> &EllipticalMFActivate)
{
    int nq = m_fields[0]->GetNpoints();

    EllipticalMF = Array<OneD, Array<OneD, NekDouble>>(m_mfdim);
    for (int j = 0; j < m_mfdim; ++j)
    {
        EllipticalMF[j] = Array<OneD, NekDouble>(m_spacedim * nq, 0.0);
    }

    EllipticalMFActivate = Array<OneD, int>(nq, 0);

    Array<OneD, NekDouble> x0(nq);
    Array<OneD, NekDouble> x1(nq);
    Array<OneD, NekDouble> x2(nq);

    m_fields[0]->GetCoords(x0, x1, x2);

    NekDouble rad, xp, yp, zp;
    NekDouble sin_varphi, cos_varphi, sin_theta, cos_theta;
    NekDouble Tol = 0.1;
    Array<OneD, NekDouble> radvec;

    int indexj, indexk;
    NekDouble tmpx, tmpy, tmpz, mag;
    for (int i = 0; i < m_fields[0]->GetExpSize(); ++i)
    {
        radvec = Array<OneD, NekDouble>(m_fields[0]->GetTotPoints(i));
        for (int j = 0; j < m_fields[0]->GetTotPoints(i); ++j)
        {
            indexj = m_fields[0]->GetPhys_Offset(i) + j;

            xp = x0[indexj];
            yp = x1[indexj];
            zp = x2[indexj];

            radvec[j] = sqrt((xp - m_Initx) * (xp - m_Initx) / Radx / Radx +
                             (yp - m_Inity) * (yp - m_Inity) / Rady / Rady +
                             (zp - m_Initz) * (zp - m_Initz) / Radz / Radz);

            rad = sqrt((xp + m_Initx) * (xp + m_Initx) / Radx / Radx +
                       (yp + m_Inity) * (yp + m_Inity) / Rady / Rady +
                       (zp + m_Initz) * (zp + m_Initz) / Radz / Radz);

            if (rad < radvec[j])
            {
                radvec[j] = rad;
            }
        }

        // x = r \sin \theta \cos \varphi
        // y = r \sin \theta \sin \varphi
        // z = r \cos \theta
        // \theta \in [0, \pi], \varphi = [0, 2 \pi]
        // if (fabs(zpavg) < m_SphereExactRadius * (1.0 - Tol))

        for (int k = 0; k < m_fields[0]->GetTotPoints(i); ++k)
        {
            indexk = m_fields[0]->GetPhys_Offset(i) + k;

            xp = x0[indexk];
            yp = x1[indexk];
            zp = x2[indexk];

            CartesianToElliptical(xp, yp, zp, rad, sin_varphi, cos_varphi,
                                  sin_theta, cos_theta);

            if (Vmath::Vmin(m_fields[0]->GetTotPoints(i), radvec, 1) > Tol)
            {
                EllipticalMFActivate[indexk] = 1;
            }

            // Theta - direction
            tmpx = Radx * cos_theta * cos_varphi;
            tmpy = Rady * cos_theta * sin_varphi;
            tmpz = -Radz * sin_theta;

            mag = sqrt(tmpx * tmpx + tmpy * tmpy + tmpz * tmpz);

            EllipticalMF[0][indexk]          = tmpx / mag;
            EllipticalMF[0][indexk + nq]     = tmpy / mag;
            EllipticalMF[0][indexk + 2 * nq] = tmpz / mag;

            // Phi - direction
            EllipticalMF[1][indexk]          = -1.0 * sin_varphi;
            EllipticalMF[1][indexk + nq]     = cos_varphi;
            EllipticalMF[1][indexk + 2 * nq] = 0.0;

            // // R - direction
            // tmpx = Radx * sin_theta * cos_varphi;
            // tmpy = Rady * sin_theta * sin_varphi;
            // tmpz = Radz * cos_theta;
            // mag = sqrt(tmpx*tmpx + tmpy*tmpy + tmpz*tmpz);
            // EllipticalMF[2][indexk]          = tmpx/mag;
            // EllipticalMF[2][indexk + nq]     = tmpy/mag;
            // EllipticalMF[2][indexk + 2 * nq] = tmpz/mag;
        }
    }

    // e3 = e1 \times e 2
    EllipticalMF[2] = VectorCrossProdMF(EllipticalMF[0], EllipticalMF[1]);

    std::cout << "MMFActivation = " << Vmath::Vsum(nq, EllipticalMFActivate, 1)
              << " / " << nq << std::endl;
}

void MMFSystem::TestSphericalConnection1form(
    const Array<OneD, const Array<OneD, NekDouble>> &SphereMF,
    Array<OneD, int> &SphereMFActivate)
{
    int Nconnection = 4;
    int nq          = m_fields[0]->GetNpoints();

    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> SphereMFConnection(
        m_mfdim);
    for (int i = 0; i < m_mfdim; i++)
    {
        SphereMFConnection[i] =
            Array<OneD, Array<OneD, NekDouble>>(Nconnection);
        for (int j = 0; j < m_mfdim; j++)
        {
            SphereMFConnection[i][j] = Array<OneD, NekDouble>(nq, 0.0);
        }
    }

    // Check the Curvature 1-form of the aligned moving frames
    Compute2DConnection1form(SphereMF, SphereMFConnection);

    // Test moving frames Connection whenever it is possible
    Test2DConnection1formSphere(SphereMFActivate, SphereMF, SphereMFConnection);
}

void MMFSystem::TestPseudosphericalConnection1form(
    const Array<OneD, const Array<OneD, NekDouble>> &PseudosphereMF,
    Array<OneD, int> &PseudosphereMFActivate)
{
    int Nconnection = 4;
    int nq          = m_fields[0]->GetNpoints();

    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> PseudosphereMFConnection(
        m_mfdim);
    for (int i = 0; i < m_mfdim; i++)
    {
        PseudosphereMFConnection[i] =
            Array<OneD, Array<OneD, NekDouble>>(Nconnection);
        for (int j = 0; j < m_mfdim; j++)
        {
            PseudosphereMFConnection[i][j] = Array<OneD, NekDouble>(nq, 0.0);
        }
    }

    // Check the Curvature 1-form of the aligned moving frames
    Compute2DConnection1form(PseudosphereMF, PseudosphereMFConnection);

    // Test moving frames Connection whenever it is possible
    Test2DConnection1formPseudosphere(PseudosphereMFActivate, PseudosphereMF,
                                      PseudosphereMFConnection);
}

void MMFSystem::TestPolarConnection1form(
    const Array<OneD, const Array<OneD, NekDouble>> &PolarMF,
    Array<OneD, int> &PolarMFActivate)
{
    int Nconnection = 4;
    int nq          = m_fields[0]->GetNpoints();

    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> PolarMFConnection(m_mfdim);
    for (int i = 0; i < m_mfdim; i++)
    {
        PolarMFConnection[i] = Array<OneD, Array<OneD, NekDouble>>(Nconnection);
        for (int j = 0; j < m_mfdim; j++)
        {
            PolarMFConnection[i][j] = Array<OneD, NekDouble>(nq, 0.0);
        }
    }

    // Check the Curvature 1-form of the aligned moving frames
    Compute2DConnection1form(PolarMF, PolarMFConnection);

    // Test moving frames Connection whenever it is possible
    Test2DConnection1formPolar(0.0, 0.0, 0.0, PolarMFActivate, PolarMF,
                               PolarMFConnection);
}

void MMFSystem::ConstructPolarMF(Array<OneD, Array<OneD, NekDouble>> &PolarMF,
                                 Array<OneD, int> &PolarMFActivate)
{
    int nq = m_fields[0]->GetNpoints();

    PolarMF = Array<OneD, Array<OneD, NekDouble>>(m_mfdim);
    for (int j = 0; j < m_mfdim; ++j)
    {
        PolarMF[j] = Array<OneD, NekDouble>(m_spacedim * nq, 0.0);
    }

    Array<OneD, NekDouble> x0(nq);
    Array<OneD, NekDouble> x1(nq);
    Array<OneD, NekDouble> x2(nq);

    m_fields[0]->GetCoords(x0, x1, x2);

    NekDouble Tol = 0.1;
    NekDouble xp, yp, rad;
    Array<OneD, NekDouble> radvec;

    PolarMFActivate = Array<OneD, int>(nq, 0);

    int indexj, indexk;
    int cnt = 0;
    for (int i = 0; i < m_fields[0]->GetExpSize(); ++i)
    {
        radvec = Array<OneD, NekDouble>(m_fields[0]->GetTotPoints(i), 0.0);
        for (int j = 0; j < m_fields[0]->GetTotPoints(i); ++j)
        {
            indexj = m_fields[0]->GetPhys_Offset(i) + j;

            xp = x0[indexj];
            yp = x1[indexj];

            radvec[j] = sqrt(xp * xp + yp * yp);
        }

        // x = r \sin \theta \cos \varphi
        // y = r \sin \theta \sin \varphi
        // z = r \cos \theta
        // \theta \in [0, \pi], \varphi = [0, 2 \pi]
        if (Vmath::Vmin(m_fields[0]->GetTotPoints(i), radvec, 1) > Tol)
        {
            for (int k = 0; k < m_fields[0]->GetTotPoints(i); ++k)
            {
                indexk = m_fields[0]->GetPhys_Offset(i) + k;

                xp = x0[indexk];
                yp = x1[indexk];

                rad = sqrt(xp * xp + yp * yp);

                PolarMF[0][indexk]          = xp / rad;
                PolarMF[0][indexk + nq]     = yp / rad;
                PolarMF[0][indexk + 2 * nq] = 0.0;

                PolarMF[1][indexk]          = -yp / rad;
                PolarMF[1][indexk + nq]     = xp / rad;
                PolarMF[1][indexk + 2 * nq] = 0.0;

                PolarMF[2][indexk]          = 0.0;
                PolarMF[2][indexk + nq]     = 0.0;
                PolarMF[2][indexk + 2 * nq] = 1.0;

                PolarMFActivate[indexk] = 1;
            }

            cnt++;
        }
    }

    std::cout << "Polar Activation = " << cnt << " / " << m_fields[0]->GetExpSize() << std::endl;
}

void MMFSystem::ComputeCurl(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray)

{
    int nq = inarray[0].size();

    Array<OneD, NekDouble> tmpx, tmpy, tmpz;
    Array<OneD, NekDouble> Dtmpzdx, Dtmpydx, Dtmpxdy, Dtmpzdy, Dtmpxdz, Dtmpydz;

    tmpx = Array<OneD, NekDouble>(nq);
    tmpy = Array<OneD, NekDouble>(nq);
    tmpz = Array<OneD, NekDouble>(nq);

    Dtmpzdx = Array<OneD, NekDouble>(nq);
    Dtmpydx = Array<OneD, NekDouble>(nq);
    Dtmpxdy = Array<OneD, NekDouble>(nq);
    Dtmpzdy = Array<OneD, NekDouble>(nq);
    Dtmpxdz = Array<OneD, NekDouble>(nq);
    Dtmpydz = Array<OneD, NekDouble>(nq);

    for (int k = 0; k < m_spacedim; ++k)
    {
        Vmath::Vcopy(nq, &inarray[0][0], 1, &tmpx[0], 1);
        Vmath::Vcopy(nq, &inarray[1][0], 1, &tmpy[0], 1);
        Vmath::Vcopy(nq, &inarray[2][0], 1, &tmpz[0], 1);

        m_fields[0]->PhysDeriv(0, tmpz, Dtmpzdx);
        m_fields[0]->PhysDeriv(0, tmpy, Dtmpydx);
        m_fields[0]->PhysDeriv(1, tmpx, Dtmpxdy);
        m_fields[0]->PhysDeriv(1, tmpz, Dtmpzdy);
        m_fields[0]->PhysDeriv(2, tmpx, Dtmpxdz);
        m_fields[0]->PhysDeriv(2, tmpy, Dtmpydz);

        Vmath::Vsub(nq, &Dtmpzdy[0], 1, &Dtmpydz[0], 1, &outarray[0][0], 1);
        Vmath::Vsub(nq, &Dtmpxdz[0], 1, &Dtmpzdx[0], 1, &outarray[1][0], 1);
        Vmath::Vsub(nq, &Dtmpydx[0], 1, &Dtmpxdy[0], 1, &outarray[2][0], 1);
    }
}

// void MMFSystem::SetUpMovingFrames(
//     const Array<OneD, const Array<OneD, NekDouble>> &Anisotropy)
// {
//     int nq = m_fields[0]->GetNpoints();

//     // Construct The Moving Frames
//     m_movingframes = Array<OneD, Array<OneD, NekDouble>>(m_spacedim);
//     for (int j = 0; j < m_spacedim; ++j)
//     {
//         m_movingframes[j] = Array<OneD, NekDouble>(m_spacedim * nq, 0.0);
//     }

//     // Read MMF Geom Info
//     std::string MMFdirStr = "LOCAL";
//     m_session->LoadSolverInfo("MMFDir", MMFdirStr, "LOCAL");
//     m_MMFdir = FindMMFdir(MMFdirStr);

//     // (x-x_0)^2/a^2 + (y-y_0)^2/b^2 = 1
//     // factors[0] = a
//     // factors[1] = b
//     // factors[2] = x_0
//     // factors[3] = y_0
//     m_MMFfactors = Array<OneD, NekDouble>(4);
//     m_session->LoadParameter("MMFCircAxisX", m_MMFfactors[0], 1.0);
//     m_session->LoadParameter("MMFCircAxisY", m_MMFfactors[1], 1.0);
//     m_session->LoadParameter("MMFCircCentreX", m_MMFfactors[2], 0.0);
//     m_session->LoadParameter("MMFCircCentreY", m_MMFfactors[3], 0.0);

//     // Get Tangetn vectors from GeomFactors2D, Orthonormalized = true
//     m_fields[0]->GetMovingFrames(m_MMFdir, m_MMFfactors, m_movingframes);

//     // Multiply Anisotropy to movingframes
//     for (int j = 0; j < m_shapedim; ++j)
//     {
//         for (int k = 0; k < m_spacedim; ++k)
//         {
//             Vmath::Vmul(nq, &Anisotropy[j][0], 1, &m_movingframes[j][k * nq], 1,
//                         &m_movingframes[j][k * nq], 1);
//         }
//     }
//     // Test the moving frames
//     CheckMovingFrames(m_movingframes);
// }

void MMFSystem::SetUpMovingFrames(
    const SpatialDomains::GeomMMF MMFdir,
    const Array<OneD, const Array<OneD, NekDouble>> &Anistrength,
    Array<OneD, Array<OneD, NekDouble>> &movingframes)
{
    int nq = m_fields[0]->GetNpoints();

    // Construct The Moving Frames
    movingframes = Array<OneD, Array<OneD, NekDouble>>(m_mfdim);
    for (int j = 0; j < m_mfdim; ++j)
    {
        movingframes[j] = Array<OneD, NekDouble>(m_spacedim * nq, 0.0);
    }

    // (x-x_0)^2/a^2 + (y-y_0)^2/b^2 = 1
    // factors[0] = a
    // factors[1] = b
    // factors[2] = x_0
    // factors[3] = y_0
    m_MMFfactors = Array<OneD, NekDouble>(4);

    m_session->LoadParameter("MMFCircAxisX", m_MMFfactors[0], 1.0);
    m_session->LoadParameter("MMFCircAxisY", m_MMFfactors[1], 1.0);
    m_session->LoadParameter("MMFCircCentreX", m_MMFfactors[2], 0.0);
    m_session->LoadParameter("MMFCircCentreY", m_MMFfactors[3], 0.0);

    // Get LOCAL moving frames such that \nabla \cdot e^i is the minimum
    // m_fields[0]->GetMovingFrames(m_MMFdir, m_MMFfactors, movingframes);

    // Further alignment of moving frames
    switch (MMFdir)
    {
        case SpatialDomains::eSpherical:
        {
            for (int i = 0; i < nq; ++i)
            {
                if (m_MMFActivation[i])
                {
                    for (int k = 0; k < m_mfdim; ++k)
                    {
                        movingframes[k][i]          = m_sphereMF[k][i];
                        movingframes[k][i + nq]     = m_sphereMF[k][i + nq];
                        movingframes[k][i + 2 * nq] = m_sphereMF[k][i + 2 * nq];
                    }
                }
            }
        }
        break;

        case SpatialDomains::ePolar:
        {
            for (int i = 0; i < nq; ++i)
            {
                if (m_MMFActivation[i])
                {
                    for (int k = 0; k < m_mfdim; ++k)
                    {
                        movingframes[k][i]          = m_polarMF[k][i];
                        movingframes[k][i + nq]     = m_polarMF[k][i + nq];
                        movingframes[k][i + 2 * nq] = m_polarMF[k][i + 2 * nq];
                    }
                }
            }
        }
        break;

        case SpatialDomains::ePseudospherical:
        {
            for (int i = 0; i < nq; ++i)
            {
                if (m_MMFActivation[i])
                {
                    for (int k = 0; k < m_mfdim; ++k)
                    {
                        movingframes[k][i]      = m_pseudosphereMF[k][i];
                        movingframes[k][i + nq] = m_pseudosphereMF[k][i + nq];
                        movingframes[k][i + 2 * nq] =
                            m_pseudosphereMF[k][i + 2 * nq];
                    }
                }
            }
        }
        break;

        case SpatialDomains::eLOCAL:
        {
            // GetLOCALMovingframes(movingframes);
            m_fields[0]->GetMovingFrames(SpatialDomains::eLOCAL1, m_MMFfactors, movingframes);
        }
        break;

        case SpatialDomains::eLOCALSphere:
        case SpatialDomains::eLOCALEllipsoid:
        {
            GetLOCALMovingframes(movingframes);
            ComputeAxisAlignedLOCALMovingframes(m_sphereMF, movingframes);
        }
        break;

        default:
        {
            m_fields[0]->GetMovingFrames(MMFdir, m_MMFfactors, movingframes);
        }
        break;
    }

    // Multiply Anisotropic magnitude to moving frames
    for (int i = 0; i < nq; ++i)
    {
        for (int j = 0; j < m_shapedim; ++j)
        {
            for (int k = 0; k < m_spacedim; ++k)
            {
                movingframes[j][k * nq + i] =
                    sqrt(Anistrength[j][i]) * movingframes[j][k * nq + i];
            }
        }
    }

   // CheckMovingFrames(movingframes);
}

SpatialDomains::GeomMMF MMFSystem::FindMMFdir(std::string MMFdirStr)
{
    SpatialDomains::GeomMMF MMFdir;

    if (MMFdirStr == "TangentX")
        MMFdir = SpatialDomains::eTangentX;
    if (MMFdirStr == "TangentY")
        MMFdir = SpatialDomains::eTangentY;
    if (MMFdirStr == "TangentXY")
        MMFdir = SpatialDomains::eTangentXY;
    if (MMFdirStr == "TangentZ")
        MMFdir = SpatialDomains::eTangentZ;
    if (MMFdirStr == "TangentCircular")
        MMFdir = SpatialDomains::eTangentCircular;
    if (MMFdirStr == "Polar")
        MMFdir = SpatialDomains::ePolar;
    if (MMFdirStr == "Spherical")
        MMFdir = SpatialDomains::eSpherical;
    if (MMFdirStr == "Pseudospherical")
        MMFdir = SpatialDomains::ePseudospherical;
    if (MMFdirStr == "LOCAL")
        MMFdir = SpatialDomains::eLOCAL;
    if (MMFdirStr == "LOCALSphere")
        MMFdir = SpatialDomains::eLOCALSphere;

    return MMFdir;
}

// Simple orientation check if e^3 is away from zero
void MMFSystem::CheckMFOrientation(
    Array<OneD, Array<OneD, NekDouble>> &movingframes,
    Array<OneD, NekDouble> &origin)
{
    int nq = m_fields[0]->GetNpoints();

    if (origin == NullNekDouble1DArray)
    {
        origin = Array<OneD, NekDouble>(3, 0.0);
    }

    Array<OneD, NekDouble> x0(nq);
    Array<OneD, NekDouble> x1(nq);
    Array<OneD, NekDouble> x2(nq);

    m_fields[0]->GetCoords(x0, x1, x2);

    NekDouble vx, vy, vz;
    NekDouble tmpx, tmpy, tmpz, innerdot;
    int cnt = 0;
    for (int i = 0; i < nq; ++i)
    {
        vx = x0[i] - origin[0];
        vy = x1[i] - origin[1];
        vz = x2[i] - origin[2];

        innerdot = vx * movingframes[2][i] + vy * movingframes[2][i + nq] +
                   vz * movingframes[2][i + 2 * nq];

        // if innerdot is negative, e3 points innerward. Thus, swap e1 and e2.
        // e3 = -e3;
        if (innerdot < 0)
        {
            tmpx = movingframes[0][i];
            tmpy = movingframes[0][i + nq];
            tmpz = movingframes[0][i + 2 * nq];

            movingframes[0][i]          = movingframes[1][i];
            movingframes[0][i + nq]     = movingframes[1][i + nq];
            movingframes[0][i + 2 * nq] = movingframes[1][i + 2 * nq];

            movingframes[1][i]          = tmpx;
            movingframes[1][i + nq]     = tmpy;
            movingframes[1][i + 2 * nq] = tmpz;

            movingframes[2][i]          = -1.0 * movingframes[2][i];
            movingframes[2][i + nq]     = -1.0 * movingframes[2][i + nq];
            movingframes[2][i + 2 * nq] = -1.0 * movingframes[2][i + 2 * nq];
            cnt++;
        }
    }

    if (cnt > 0)
    {
        std::cout << "Total " << cnt << " grid points have been flipped"
                  << std::endl;
    }
}

void MMFSystem::GetLOCALMovingframes(
    Array<OneD, Array<OneD, NekDouble>> &movingframes)
{
    switch (m_expdim)
    {
        case 1:
        {
            GetLOCALMovingframes1D(movingframes);
        }
        break;

        case 2:
        case 3:
        {
            GetLOCALMovingframes3D(movingframes);
        }
        break;

        default:
            break;
    }
}

void MMFSystem::GetLOCALFIBERMovingFrames(
    MultiRegions::ExpListSharedPtr &field,
    Array<OneD, Array<OneD, NekDouble>> &movingframes)
{
    int fnq = field->GetNpoints();

    movingframes = Array<OneD, Array<OneD, NekDouble>>(m_mfdim);
    for (int i = 0; i < m_mfdim; ++i)
    {
        movingframes[i] = Array<OneD, NekDouble>(m_spacedim * fnq);
    }

    field->GetMovingFrames(SpatialDomains::eLOCAL1, m_MMFfactors,
                           m_movingframes);
    // Simple orientation check if e^3 is away from zero
    // CheckMFOrientation(movingframes);

    std::cout << "SetUp Fiber MF: "
                 "***********************************************************"
              << std::endl;
}

void MMFSystem::GetLOCALMovingframes1D(
    Array<OneD, Array<OneD, NekDouble>> &movingframes)
{
    int nq = m_fields[0]->GetNpoints();

    movingframes = Array<OneD, Array<OneD, NekDouble>>(m_mfdim);
    for (int i = 0; i < m_mfdim; ++i)
    {
        movingframes[i] = Array<OneD, NekDouble>(m_spacedim * nq);
    }

    m_fields[0]->GetMovingFrames(SpatialDomains::eLOCAL, m_MMFfactors, movingframes);
    // Simple orientation check if e^3 is away from zero
    // CheckMFOrientation(movingframes);

    std::cout << "SetUp 1D MF: "
                 "***********************************************************"
              << std::endl;
}

// Choose moving frames such that \sqrt{Div(e^1)^2 + Div(e^1)^2} is the
// smallest for three coordinates.
void MMFSystem::GetLOCALMovingframes3D(
    Array<OneD, Array<OneD, NekDouble>> &movingframes)
{
    int nq = m_fields[0]->GetNpoints();

    movingframes = Array<OneD, Array<OneD, NekDouble>>(m_mfdim);
    for (int i = 0; i < m_mfdim; ++i)
    {
        movingframes[i] = Array<OneD, NekDouble>(m_spacedim * nq);
    }

    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> MF(m_mfdim);
    for (int i = 0; i < m_mfdim; ++i)
    {
        MF[i] = Array<OneD, Array<OneD, NekDouble>>(m_mfdim);
        for (int j = 0; j < m_mfdim; ++j)
        {
            MF[i][j] = Array<OneD, NekDouble>(m_spacedim * nq);
        }
    }

    // m_fields[0]->GetMovingFrames(m_MMFdir, m_MMFfactors, m_movingframes);
    m_fields[0]->GetMovingFrames(SpatialDomains::eLOCAL1, m_MMFfactors, MF[0]);

    m_fields[0]->GetMovingFrames(SpatialDomains::eLOCAL2, m_MMFfactors, MF[1]);
    if (m_expdim == 2)
    {
        m_fields[0]->GetMovingFrames(SpatialDomains::eLOCAL21, m_MMFfactors,
                                     MF[2]);
    }

    if (m_expdim == 3)
    {
        m_fields[0]->GetMovingFrames(SpatialDomains::eLOCAL3, m_MMFfactors,
                                     MF[2]);
    }

    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> DivMF(m_mfdim);
    for (int i = 0; i < m_mfdim; i++)
    {
        DivMF[i] = Array<OneD, Array<OneD, NekDouble>>(m_mfdim);
        for (int j = 0; j < m_mfdim; j++)
        {
            DivMF[i][j] = Array<OneD, NekDouble>(nq, 0.0);
        }
    }

    // Compute Divergence in two different ways
    for (int i = 0; i < m_mfdim; i++)
    {
        ComputeDivMF(m_DerivType, MF[i], DivMF[i], 0);
    }

    // Elementwise comparison
    Array<OneD, NekDouble> DivMFsum(m_mfdim);
    int MF1cnt = 0, MF2cnt = 0, MF0cnt = 0;
    int MFINDEX = 0, indexj;
    for (int i = 0; i < m_fields[0]->GetExpSize(); ++i)
    {
        DivMFsum = Array<OneD, NekDouble>(m_mfdim, 0.0);
        for (int j = 0; j < m_fields[0]->GetTotPoints(i); ++j)
        {
            indexj = m_fields[0]->GetPhys_Offset(i) + j;

            // For each element, compute \sqrt{Dx^2 + Dy^2}
            for (int k = 0; k < m_mfdim; ++k)
            {
                for (int l = 0; l < m_spacedim; ++l)
                {
                    DivMFsum[k] += DivMF[k][l][indexj] * DivMF[k][l][indexj];
                }

                DivMFsum[k] = sqrt(DivMFsum[k]);
            }
        }

        // Get the index for the smallest DivMFsum
        MFINDEX = Vmath::Imin(m_mfdim, DivMFsum, 1);
        switch (MFINDEX)
        {
            case 0:
            {
                MF0cnt++;
            }
            break;

            case 1:
            {
                MF1cnt++;
            }
            break;

            case 2:
            {
                MF2cnt++;
            }
            break;

            default:
                break;
        }

        // Use such Index for the local frames
        for (int k = 0; k < m_mfdim; ++k)
        {
            for (int j = 0; j < m_fields[0]->GetTotPoints(i); ++j)
            {
                indexj = m_fields[0]->GetPhys_Offset(i) + j;

                movingframes[k][indexj]      = MF[MFINDEX][k][indexj];
                movingframes[k][indexj + nq] = MF[MFINDEX][k][indexj + nq];
                movingframes[k][indexj + 2 * nq] = MF[MFINDEX][k][indexj + 2 * nq];
            }
        }
    }
    
    // Simple orientation check if e^3 is away from zero. origin = (0,0,0)
    CheckMFOrientation(movingframes);

    std::cout << "SetUp " << m_expdim
              << "D MF: "
                 "***********************************************************"
              << std::endl;

    std::cout << "LOCAL construction for " << m_fields[0]->GetExpSize()
              << " elements, adapated from MF1 = " << MF0cnt
              << ", MF2 = " << MF1cnt << ", MF3 = " << MF2cnt << std::endl;
}

void MMFSystem::ConstructAnisotropicFrames(
    const Array<OneD, const NekDouble> &AniDirection,
    Array<OneD, Array<OneD, NekDouble>> &movingframes,
    Array<OneD, NekDouble> &AniConstruction)
{
    int nq = m_fields[0]->GetNpoints();

    int indexj, indexk;
    int fibreElemnt = 0, cntflipped = 0;
    int verticalfibre = 0, crossingfibre = 0;

    Array<OneD, NekDouble> fcdotk(nq, 0.0);
    AniConstruction = Array<OneD, NekDouble>(nq, 0.0);

    NekDouble e1x, e1y, e1z, e2x, e2y, e2z, e3x, e3y, e3z;
    NekDouble fx, fy, fz, mag, fdote3;
    NekDouble fxavg, fyavg, fzavg, avgdirtest;
    NekDouble ejx, ejy, ejz, ekx, eky, ekz;
    NekDouble ejcdotek, ejcdotekmin, magej, magek;
    NekDouble dotTol = 0.0;

    NekDouble magTol    = 0.1;
    NekDouble e1ProjErr = 0.0;
    int FibreConstruction;
    for (int i = 0; i < m_fields[0]->GetExpSize(); ++i)
    {
        // Change moving frames along fibre direction if fibre exist for
        // all points for the element
        FibreConstruction = 1;
        fxavg             = 0.0;
        fyavg             = 0.0;
        fzavg             = 0.0;
        for (int j = 0; j < m_fields[0]->GetTotPoints(i); ++j)
        {
            indexj = m_fields[0]->GetPhys_Offset(i) + j;

            fx = AniDirection[indexj];
            fy = AniDirection[indexj + nq];
            fz = AniDirection[indexj + 2 * nq];

            // \vec{e}^3 is the same as before:
            e3x = movingframes[2][indexj];
            e3y = movingframes[2][indexj + nq];
            e3z = movingframes[2][indexj + 2 * nq];

            // Project f into e^1 \times e^2 plane
            fdote3 = fx * e3x + fy * e3y + fz * e3z;
            fx     = fx - fdote3 * e3x;
            fy     = fy - fdote3 * e3y;
            fz     = fz - fdote3 * e3z;

            fcdotk[indexj] = fdote3;

            fxavg += fx;
            fyavg += fy;
            fzavg += fz;

            mag = sqrt(fx * fx + fy * fy + fz * fz);

            if (mag < magTol)
            {
                FibreConstruction = -1;
            }

            for (int k = j; k < m_fields[0]->GetTotPoints(i); ++k)
            {
                indexk = m_fields[0]->GetPhys_Offset(i) + k;

                ekx = AniDirection[indexk];
                eky = AniDirection[indexk + nq];
                ekz = AniDirection[indexk + 2 * nq];

                e3x = movingframes[2][indexk];
                e3y = movingframes[2][indexk + nq];
                e3z = movingframes[2][indexk + 2 * nq];

                fdote3 = ekx * e3x + eky * e3y + ekz * e3z;
                ekx    = ekx - fdote3 * e3x;
                eky    = eky - fdote3 * e3y;
                ekz    = ekz - fdote3 * e3z;

                magek = sqrt(ekx * ekx + eky * eky + ekz * ekz);

                ejcdotek = (fx * ekx + fy * eky + fz * ekz) / mag / magek;
                if (ejcdotek < 0.0)
                {
                    FibreConstruction = -2;
                }
            }
        } // end of j-loop

        if (FibreConstruction == -1)
        {
            verticalfibre++;
        }

        else if (FibreConstruction == -2)
        {
            crossingfibre++;
        }

        fxavg = fxavg / m_fields[0]->GetTotPoints(i);
        fyavg = fyavg / m_fields[0]->GetTotPoints(i);
        fzavg = fzavg / m_fields[0]->GetTotPoints(i);

        // If fibre's magnitude is too small, i.e.,
        // if it is aligned along the normal direction, then do not
        // construct the fibre
        if (FibreConstruction > 0)
        {
            fibreElemnt++;

            for (int k = 0; k < m_fields[0]->GetTotPoints(i); ++k)
            {
                indexk = m_fields[0]->GetPhys_Offset(i) + k;

                fx = AniDirection[indexk];
                fy = AniDirection[indexk + nq];
                fz = AniDirection[indexk + 2 * nq];

                // \vec{e}^3 is the same as before:
                e3x = movingframes[2][indexk];
                e3y = movingframes[2][indexk + nq];
                e3z = movingframes[2][indexk + 2 * nq];

                // Project f into e^1 \times e^2 plane
                fdote3 = fx * e3x + fy * e3y + fz * e3z;
                fx     = fx - fdote3 * e3x;
                fy     = fy - fdote3 * e3y;
                fz     = fz - fdote3 * e3z;

                mag = sqrt(fx * fx + fy * fy + fz * fz);

                // \vec{e}^1 = \vec{f}
                avgdirtest = fxavg * fx + fyavg * fy + fzavg * fz;
                if (avgdirtest > 0)
                {
                    e1x = fx;
                    e1y = fy;
                    e1z = fz;
                }

                else
                {
                    e1x = -1.0 * fx;
                    e1y = -1.0 * fy;
                    e1z = -1.0 * fz;
                    cntflipped++;
                }

                // e2 is the unit vector, orthogonal to e^1 and e^3
                e2x = (e3y * e1z - e1y * e3z) / mag;
                e2y = (e1x * e3z - e3x * e1z) / mag;
                e2z = (e3x * e1y - e1x * e3y) / mag;

                e1ProjErr += fabs(e1x * e3x + e1y * e3y + e1z * e3z);

                movingframes[0][indexk]          = e1x;
                movingframes[0][indexk + nq]     = e1y;
                movingframes[0][indexk + 2 * nq] = e1z;

                movingframes[1][indexk]          = e2x;
                movingframes[1][indexk + nq]     = e2y;
                movingframes[1][indexk + 2 * nq] = e2z;

                AniConstruction[indexk] = 1.0;
            } // end of k-loop
        }
    } // end of i-loop

    // Test differentiability of fibre

    int Nondifferentiable = 0, indj = 0, indk = 0;
    int cnt = 0;
    for (int i = 0; i < m_fields[0]->GetExpSize(); ++i)
    {
        Nondifferentiable = 0;
        ejcdotekmin       = 1.0;
        for (int j = 0; j < m_fields[0]->GetTotPoints(i); ++j)
        {
            indexj = m_fields[0]->GetPhys_Offset(i) + j;

            ejx = movingframes[0][indexj];
            ejy = movingframes[0][indexj + nq];
            ejz = movingframes[0][indexj + 2 * nq];

            magej = sqrt(ejx * ejx + ejy * ejy + ejz * ejz);
            for (int k = j; k < m_fields[0]->GetTotPoints(i); ++k)
            {
                indexk = m_fields[0]->GetPhys_Offset(i) + k;

                ekx = movingframes[0][indexk];
                eky = movingframes[0][indexk + nq];
                ekz = movingframes[0][indexk + 2 * nq];

                magek = sqrt(ekx * ekx + eky * eky + ekz * ekz);

                ejcdotek = (ejx * ekx + ejy * eky + ejz * ekz) / magej / magek;

                if (ejcdotek < dotTol)
                {
                    Nondifferentiable = 1;
                    if (ejcdotek < ejcdotekmin)
                    {
                        ejcdotekmin = ejcdotek;
                        indj        = j;
                        indk        = k;
                    }
                }
            } // end of k-loop
        }     // end of j-loop

        if (Nondifferentiable)
        {
            cnt++;
            std::cout << "Nel = " << i << ", min. dot = " << ejcdotekmin
                      << " at j = " << indj << ", k = " << indk << std::endl;
            for (int j = 0; j < m_fields[0]->GetTotPoints(i); ++j)
            {
                indexj = m_fields[0]->GetPhys_Offset(i) + j;

                ejx = movingframes[0][indexj];
                ejy = movingframes[0][indexj + nq];
                ejz = movingframes[0][indexj + 2 * nq];

                std::cout << "index = " << indexj << ", " << j << "th, e^1 = ( "
                          << ejx << " , " << ejy << " , " << ejz << " ) "
                          << std::endl;
            }
            std::cout << " " << std::endl;
        }
    }

    std::cout << "Vertical fibre = " << verticalfibre
              << " elements are not activated" << std::endl;
    std::cout << "Crossing fibre = " << crossingfibre
              << " elements are not activated" << std::endl;
    std::cout << cntflipped << " moving frames are flipped " << std::endl
              << std::endl;
    std::cout << "Total " << cnt << " elements are possibly non-differentiable"
              << std::endl;
    std::cout << "AniConstrictopm is activated for "
              << Vmath::Vsum(nq, AniConstruction, 1) << " of " << nq
              << std::endl;

    e1ProjErr = e1ProjErr / nq;

    NekDouble fibreconstr = 1.0 * fibreElemnt / m_fields[0]->GetTotPoints(0) * 100.0;
    std::cout << "Total " << fibreElemnt << " / " << m_fields[0]->GetTotPoints(0) << " ( "
              << fibreconstr << " % ) "
              << " elements are aligned along the fibre direction, "
                 "fiber_Proj_err = "
              << e1ProjErr << std::endl;
}

void MMFSystem::CheckMovingFrames(
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes)
{
    NekDouble t1x, t1y, t1z, t2x, t2y, t2z, t3x, t3y, t3z;
    NekDouble mag1, mag2, dot12 = 0.0, dot23 = 0.0, dot31 = 0.0;
    NekDouble Tol = 0.0001;

    int nq = m_fields[0]->GetNpoints();

    int cntunit1 = 0;
    int cntunit2 = 0;

    int cnt1         = 0;
    int cnt2         = 0;
    NekDouble MFmag1 = 0.0;
    NekDouble MFmag2 = 0.0;

    for (int i = 0; i < nq; ++i)
    {
        t1x = movingframes[0][i];
        t1y = movingframes[0][i + nq];
        t1z = movingframes[0][i + 2 * nq];

        t2x = movingframes[1][i];
        t2y = movingframes[1][i + nq];
        t2z = movingframes[1][i + 2 * nq];

        t3x = movingframes[2][i];
        t3y = movingframes[2][i + nq];
        t3z = movingframes[2][i + 2 * nq];

        mag1 = sqrt(t1x * t1x + t1y * t1y + t1z * t1z);
        if (mag1 > 0.000000001)
        {
            cnt1++;
            MFmag1 += mag1 * mag1;
        }

        mag2 = sqrt(t2x * t2x + t2y * t2y + t2z * t2z);
        if (mag2 > 0.000000001)
        {
            cnt2++;
            MFmag2 += mag2 * mag2;
        }

        if (abs(mag1 - 1.0) < 0.000001)
        {
            cntunit1++;
        }

        if (abs(mag2 - 1.0) < 0.000001)
        {
            cntunit2++;
        }

        dot12 += (t1x * t2x + t1y * t2y + t1z * t2z) *
                 (t1x * t2x + t1y * t2y + t1z * t2z);
        dot23 += (t2x * t3x + t2y * t3y + t2z * t3z) *
                 (t2x * t3x + t2y * t3y + t2z * t3z);
        dot31 += (t3x * t1x + t3y * t1y + t3z * t1z) *
                 (t3x * t1x + t3y * t1y + t3z * t1z);
    }

    dot12 = sqrt(dot12 / nq);
    dot23 = sqrt(dot23 / nq);
    dot31 = sqrt(dot31 / nq);

    MFmag1 = sqrt(MFmag1 / cnt1);
    MFmag2 = sqrt(MFmag2 / cnt1);

    if (cntunit1 == cnt1)
    {
        std::cout << "*** 1st Moving frames are in unit length ***"
                  << std::endl;
    }

    else
    {
        std::cout << "*** 1st Moving frames are NOT in unit length, Avg mag. = "
                  << MFmag1 << " *** " << std::endl;
    }

    if (cntunit2 == cnt2)
    {
        std::cout << "*** 2nd Moving frames are in unit length ***"
                  << std::endl;
    }

    else
    {
        std::cout << "*** 2nd Moving frames are NOT in unit length, Avg mag. = "
                  << MFmag2 << " *** " << std::endl;
    }

    if ((fabs(dot12) + fabs(dot23) + fabs(dot31)) < Tol)
    {
        std::cout << "*** Moving frames are Orthogonal" << std::endl;
    }

    else
    {
        std::cout << "dot12 = " << fabs(dot12) << ", dot23 = " << fabs(dot23)
                  << ", dot31 = " << fabs(dot31) << std::endl;
        std::cout << "*** Moving frames are NOT Orthogonal" << std::endl;
    }

    Array<OneD, Array<OneD, NekDouble>> tmp(m_spacedim);
    for (int j = 0; j < m_spacedim; ++j)
    {
        tmp[j] = Array<OneD, NekDouble>(nq, 0.0);
        for (int k = 0; k < m_spacedim; ++k)
        {
            Vmath::Vvtvp(nq, &movingframes[j][k * nq], 1,
                         &movingframes[j][k * nq], 1, &tmp[j][0], 1, &tmp[j][0],
                         1);
        }
        Vmath::Vsqrt(nq, tmp[j], 1, tmp[j], 1);
    }

    Array<OneD, NekDouble> tmpx(nq), tmpy(nq), tmpz(nq);
    for (int i = 0; i < m_spacedim; ++i)
    {
        Vmath::Vcopy(nq, &movingframes[i][0], 1, &tmpx[0], 1);
        Vmath::Vcopy(nq, &movingframes[i][nq], 1, &tmpy[0], 1);
        Vmath::Vcopy(nq, &movingframes[i][2 * nq], 1, &tmpz[0], 1);

        std::cout << "*** MF " << i << " = ( " << RootMeanSquare(tmpx) << " , "
                  << RootMeanSquare(tmpy) << " , " << RootMeanSquare(tmpz)
                  << " ) " << std::endl;
    }

    // for (int i=0; i<nq; ++i)
    // {
    //     std::cout << "i = " << i << ", MF1 = ( " << movingframes[0][i] << " , " << movingframes[0][i+nq] << " , " << movingframes[0][i+2*nq] 
    //     << " ) , MF2 = ( " << movingframes[1][i] << " , " << movingframes[1][i+nq] << " , " << movingframes[1][i+2*nq]  
    //     << " ) , MF3 = ( " << movingframes[2][i] << " , " << movingframes[2][i+nq] << " , " << movingframes[2][i+2*nq] << " ) " << std::endl; 
    // }
}

// 	RebuildMovingFrames(m_K1, m_K2, m_K3, m_distance);
void MMFSystem::RebuildMovingFrames(
    const int K1, const int K2, const int K3,
    Array<OneD, Array<OneD, NekDouble>> &distance,
    Array<OneD, Array<OneD, NekDouble>> &movingframes)
{
    int nq = m_fields[0]->GetNpoints();

    std::cout << "*** Moving is rebuilt with K1 = " << m_c121
              << ", K2 = " << m_c122 << ", K3 = " << m_c123 << std::endl;

    NekDouble e1x, e1y, e1z, e2x, e2y, e2z, e3x, e3y, e3z;
    NekDouble v1, v2, v3, modi, mage;

    for (int i = 0; i < nq; i++)
    {
        e1x = movingframes[0][i];
        e1y = movingframes[0][i + nq];
        e1z = movingframes[0][i + nq * 2];

        e2x = movingframes[1][i];
        e2y = movingframes[1][i + nq];
        e2z = movingframes[1][i + nq * 2];

        e3x = movingframes[2][i];
        e3y = movingframes[2][i + nq];
        e3z = movingframes[2][i + nq * 2];

        v1 = distance[0][i];
        v2 = distance[1][i];
        v3 = distance[2][i];

        // e1<v> = e1 + (k1*v1 + k2*v2 + K3*v3) e2 + v1 e3
        modi = K1 * v1 + K2 * v2 + K3 * v3;

        e1x = e1x + modi * e2x + v1 * e3x;
        e1y = e1y + modi * e2y + v1 * e3y;
        e1z = e1z + modi * e2z + v1 * e3z;

        // Projection of e1 onto the tangent plane
        modi = e1x * e3x + e1y * e3y + e1z * e3z;
        e1x  = e1x - modi * e3x;
        e1y  = e1y - modi * e3y;
        e1z  = e1z - modi * e3z;

        mage = sqrt(e1x * e1x + e1y * e1y + e1z * e1z);
        e1x  = e1x / mage;
        e1y  = e1y / mage;
        e1z  = e1z / mage;

        // Construct MF2
        e2x = e1y * e3z - e1z * e3y;
        e2y = -e1x * e3z + e1z * e3x;
        e2z = e1x * e3y - e1y * e3x;

        mage = sqrt(e2x * e2x + e2y * e2y + e2z * e2z);
        e2x  = e2x / mage;
        e2y  = e2y / mage;
        e2z  = e2z / mage;

        movingframes[0][i]          = e1x;
        movingframes[0][i + nq]     = e1y;
        movingframes[0][i + nq * 2] = e1z;

        movingframes[1][i]          = e2x;
        movingframes[1][i + nq]     = e2y;
        movingframes[1][i + nq * 2] = e2z;

        movingframes[2][i]          = e3x;
        movingframes[2][i + nq]     = e3y;
        movingframes[2][i + nq * 2] = e3z;
    }
}

void MMFSystem::GetFwdBwdMFTrace(const Array<OneD, const NekDouble> &inarray,
                                 Array<OneD, NekDouble> &Fwd,
                                 Array<OneD, NekDouble> &Bwd)
{
    int nTracePointsTot = GetTraceNpoints();

    Fwd = Array<OneD, NekDouble>(nTracePointsTot);
    Bwd = Array<OneD, NekDouble>(nTracePointsTot);

    m_fields[0]->GetFwdBwdTracePhys(inarray, Fwd, Bwd);

    const Array<OneD, const int> &traceBndMap = m_fields[0]->GetTraceBndMap();

    // Provide the same MF for BWD
    int id2, npts, cnt = 0;
    for (int bcRegion = 0; bcRegion < m_fields[0]->GetBndConditions().size();
         ++bcRegion)
    {
        for (int e = 0;
             e < m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExpSize();
             ++e)
        {
            npts = m_fields[0]
                       ->GetBndCondExpansions()[bcRegion]
                       ->GetExp(e)
                       ->GetNumPoints(0);

            id2 = m_fields[0]->GetTrace()->GetPhys_Offset(traceBndMap[cnt + e]);

            // Not sure why Neumann boundary for diffusion does not work for
            // this substitution.
            if ((m_fields[0]->GetBndConditions()[bcRegion])
                    ->GetBoundaryConditionType() != SpatialDomains::eNeumann)
            {
                Vmath::Vcopy(npts, &Fwd[id2], 1, &Bwd[id2], 1);
            }
        }
        cnt += m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExpSize();
    }
}

void MMFSystem::ComputencdotMF(
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    Array<OneD, Array<OneD, NekDouble>> &ncdotMFFwd,
    Array<OneD, Array<OneD, NekDouble>> &ncdotMFBwd, const int Verbose)

{
    int nq              = m_fields[0]->GetNpoints();
    int nTracePointsTot = GetTraceNpoints();

    // Compute MFjFwd and MFjBwd
    Array<OneD, NekDouble> tmp(nq);

    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> MFtraceFwd(m_expdim);
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> MFtraceBwd(m_expdim);
    for (int j = 0; j < m_expdim; ++j)
    {
        MFtraceFwd[j] = Array<OneD, Array<OneD, NekDouble>>(m_spacedim);
        MFtraceBwd[j] = Array<OneD, Array<OneD, NekDouble>>(m_spacedim);
        for (int k = 0; k < m_spacedim; ++k)
        {
            MFtraceFwd[j][k] = Array<OneD, NekDouble>(nTracePointsTot);
            MFtraceBwd[j][k] = Array<OneD, NekDouble>(nTracePointsTot);

            Vmath::Vcopy(nq, &movingframes[j][k * nq], 1, &tmp[0], 1);

            GetFwdBwdMFTrace(tmp, MFtraceFwd[j][k], MFtraceBwd[j][k]);
        }
    }

    // Compute n \times e^i
    ncdotMFFwd = Array<OneD, Array<OneD, NekDouble>>(m_mfdim);
    ncdotMFBwd = Array<OneD, Array<OneD, NekDouble>>(m_mfdim);
    for (int j = 0; j < m_mfdim; ++j)
    {
        ncdotMFFwd[j] = Array<OneD, NekDouble>(nTracePointsTot, 0.0);
        ncdotMFBwd[j] = Array<OneD, NekDouble>(nTracePointsTot, 0.0);
    }

    for (int j = 0; j < m_expdim; ++j)
    {
        VectorDotProd(m_traceNormals, MFtraceFwd[j], ncdotMFFwd[j]);
        VectorDotProd(m_traceNormals, MFtraceBwd[j], ncdotMFBwd[j]);
    }

    if (Verbose)
    {
        Array<OneD, NekDouble> x0(nq);
        Array<OneD, NekDouble> x1(nq);
        Array<OneD, NekDouble> x2(nq);

        m_fields[0]->GetCoords(x0, x1, x2);

        Array<OneD, NekDouble> pFwd(nTracePointsTot);
        Array<OneD, NekDouble> pBwd(nTracePointsTot);

        m_fields[0]->GetFwdBwdTracePhys(x0, pFwd, pBwd);
    }
}

void MMFSystem::ComputentimesMF(
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &ntimesMFFwd,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &ntimesMFBwd,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &ntimes_ntimesMFFwd,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &ntimes_ntimesMFBwd)
{
    int nTracePointsTot = GetTraceNpoints();

    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> MFtraceFwd;
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> MFtraceBwd;
    ComputeMFtrace(movingframes, MFtraceFwd, MFtraceBwd);

    ntimesMFFwd = Array<OneD, Array<OneD, Array<OneD, NekDouble>>>(m_mfdim);
    ntimesMFBwd = Array<OneD, Array<OneD, Array<OneD, NekDouble>>>(m_mfdim);
    ntimes_ntimesMFFwd =
        Array<OneD, Array<OneD, Array<OneD, NekDouble>>>(m_mfdim);
    ntimes_ntimesMFBwd =
        Array<OneD, Array<OneD, Array<OneD, NekDouble>>>(m_mfdim);
    for (int j = 0; j < m_mfdim; ++j)
    {
        ntimesMFFwd[j]        = Array<OneD, Array<OneD, NekDouble>>(m_spacedim);
        ntimesMFBwd[j]        = Array<OneD, Array<OneD, NekDouble>>(m_spacedim);
        ntimes_ntimesMFFwd[j] = Array<OneD, Array<OneD, NekDouble>>(m_spacedim);
        ntimes_ntimesMFBwd[j] = Array<OneD, Array<OneD, NekDouble>>(m_spacedim);

        for (int k = 0; k < m_spacedim; ++k)
        {
            ntimesMFFwd[j][k]        = Array<OneD, NekDouble>(nTracePointsTot);
            ntimesMFBwd[j][k]        = Array<OneD, NekDouble>(nTracePointsTot);
            ntimes_ntimesMFFwd[j][k] = Array<OneD, NekDouble>(nTracePointsTot);
            ntimes_ntimesMFBwd[j][k] = Array<OneD, NekDouble>(nTracePointsTot);
        }

        VectorCrossProd(m_traceNormals, MFtraceFwd[j], ntimesMFFwd[j]);
        VectorCrossProd(m_traceNormals, MFtraceBwd[j], ntimesMFBwd[j]);
        VectorCrossProd(m_traceNormals, ntimesMFFwd[j], ntimes_ntimesMFFwd[j]);
        VectorCrossProd(m_traceNormals, ntimesMFBwd[j], ntimes_ntimesMFBwd[j]);
    }

    std::cout << "*** m_ntimesMFFwd = ( " << VectorAvgMagnitude(ntimesMFFwd[0])
              << " , " << VectorAvgMagnitude(ntimesMFFwd[1]) << " , "
              << VectorAvgMagnitude(ntimesMFFwd[2]) << " ) " << std::endl;
    std::cout << "*** m_ntimesMFBwd = ( " << VectorAvgMagnitude(ntimesMFBwd[0])
              << " , " << VectorAvgMagnitude(ntimesMFBwd[1]) << " , "
              << VectorAvgMagnitude(ntimesMFBwd[2]) << " ) " << std::endl;
    std::cout << "*** m_ntimes_ntimesMFFwd = ( "
              << VectorAvgMagnitude(ntimes_ntimesMFFwd[0]) << " , "
              << VectorAvgMagnitude(ntimes_ntimesMFFwd[1]) << " , "
              << VectorAvgMagnitude(ntimes_ntimesMFFwd[2]) << " ) "
              << std::endl;
    std::cout << "*** m_ntimes_ntimesMFBwd = ( "
              << VectorAvgMagnitude(ntimes_ntimesMFBwd[0]) << " , "
              << VectorAvgMagnitude(ntimes_ntimesMFBwd[1]) << " , "
              << VectorAvgMagnitude(ntimes_ntimesMFBwd[2]) << " ) "
              << std::endl;
}

void MMFSystem::ComputentimesMF(
    const Array<OneD, const Array<OneD, NekDouble>> &MF1st,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &ntimesMFFwd,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &ntimesMFBwd)
{
    int nTracePointsTot = GetTraceNpoints();

    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> MFtraceFwd;
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> MFtraceBwd;

    ComputeMFtrace(MF1st, MFtraceFwd, MFtraceBwd);

    ntimesMFFwd = Array<OneD, Array<OneD, Array<OneD, NekDouble>>>(m_mfdim);
    ntimesMFBwd = Array<OneD, Array<OneD, Array<OneD, NekDouble>>>(m_mfdim);
    for (int j = 0; j < m_mfdim; ++j)
    {
        ntimesMFFwd[j] = Array<OneD, Array<OneD, NekDouble>>(m_spacedim);
        ntimesMFBwd[j] = Array<OneD, Array<OneD, NekDouble>>(m_spacedim);

        for (int k = 0; k < m_mfdim; ++k)
        {
            ntimesMFFwd[j][k] = Array<OneD, NekDouble>(nTracePointsTot);
            ntimesMFBwd[j][k] = Array<OneD, NekDouble>(nTracePointsTot);
        }

        VectorCrossProd(m_traceNormals, MFtraceFwd[j], ntimesMFFwd[j]);
        VectorCrossProd(m_traceNormals, MFtraceBwd[j], ntimesMFBwd[j]);
    }

    // std::cout << "*** m_ntimesMFFwd = ( " <<
    // VectorAvgMagnitude(ntimesMFFwd[0]) << " , "
    //          << VectorAvgMagnitude(ntimesMFFwd[1]) << " , " <<
    //          VectorAvgMagnitude(ntimesMFFwd[2]) << " ) " <<
    //          std::endl;
    // std::cout << "*** m_ntimesMFBwd = ( " <<
    // VectorAvgMagnitude(ntimesMFBwd[0]) << " , "
    //          << VectorAvgMagnitude(ntimesMFBwd[1]) << " , " <<
    //          VectorAvgMagnitude(ntimesMFBwd[2]) << " ) " <<
    //          std::endl;
}

void MMFSystem::ComputenperpcdotMF(
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    Array<OneD, Array<OneD, NekDouble>> &nperpcdotMFFwd,
    Array<OneD, Array<OneD, NekDouble>> &nperpcdotMFBwd)

{
    int nq              = m_fields[0]->GetNpoints();
    int nTracePointsTot = GetTraceNpoints();

    // Compute MFjFwd and MFjBwd
    Array<OneD, NekDouble> tmp(nq);

    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> MFtraceFwd(m_shapedim);
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> MFtraceBwd(m_shapedim);
    Array<OneD, Array<OneD, NekDouble>> SurfaceNormalFwd;
    Array<OneD, Array<OneD, NekDouble>> SurfaceNormalBwd;

    for (int j = 0; j < m_shapedim; ++j)
    {
        SurfaceNormalFwd = Array<OneD, Array<OneD, NekDouble>>(m_spacedim);
        SurfaceNormalBwd = Array<OneD, Array<OneD, NekDouble>>(m_spacedim);

        for (int k = 0; k < m_spacedim; ++k)
        {
            SurfaceNormalFwd[k] = Array<OneD, NekDouble>(nTracePointsTot);
            SurfaceNormalBwd[k] = Array<OneD, NekDouble>(nTracePointsTot);
        }
    }

    ComputeMFtrace(movingframes, MFtraceFwd, MFtraceBwd);

    VectorCrossProd(MFtraceFwd[0], MFtraceFwd[1], SurfaceNormalFwd);
    VectorCrossProd(MFtraceBwd[0], MFtraceBwd[1], SurfaceNormalBwd);

    // Compute n^{\perp} \times e^i
    Array<OneD, Array<OneD, NekDouble>> TracevectorFwd(m_spacedim);
    Array<OneD, Array<OneD, NekDouble>> TracevectorBwd(m_spacedim);
    for (int k = 0; k < m_spacedim; k++)
    {
        TracevectorFwd[k] = Array<OneD, NekDouble>(nTracePointsTot);
        TracevectorBwd[k] = Array<OneD, NekDouble>(nTracePointsTot);
    }

    VectorCrossProd(m_traceNormals, SurfaceNormalFwd, TracevectorFwd);
    VectorCrossProd(m_traceNormals, SurfaceNormalBwd, TracevectorBwd);

    nperpcdotMFFwd = Array<OneD, Array<OneD, NekDouble>>(m_shapedim);
    nperpcdotMFBwd = Array<OneD, Array<OneD, NekDouble>>(m_shapedim);
    for (int j = 0; j < m_shapedim; ++j)
    {
        nperpcdotMFFwd[j] = Array<OneD, NekDouble>(nTracePointsTot, 0.0);
        nperpcdotMFBwd[j] = Array<OneD, NekDouble>(nTracePointsTot, 0.0);

        VectorDotProd(TracevectorFwd, MFtraceFwd[j], nperpcdotMFFwd[j]);
        VectorDotProd(TracevectorBwd, MFtraceBwd[j], nperpcdotMFBwd[j]);
    }
}

void MMFSystem::ComputeMFcdotSphericalCoord(const int dir)
{
    int nq = m_fields[0]->GetNpoints();

    Array<OneD, Array<OneD, NekDouble>> SphericalVector;
    ComputeSphericalVector(SphericalVector);

    m_ThetacdotMF = Array<OneD, NekDouble>(nq, 0.0);
    m_PhicdotMF   = Array<OneD, NekDouble>(nq, 0.0);

    for (int i = 0; i < m_spacedim; ++i)
    {
        Vmath::Vvtvp(nq, &m_movingframes[dir][i * nq], 1,
                     &SphericalVector[0][i * nq], 1, &m_ThetacdotMF[0], 1,
                     &m_ThetacdotMF[0], 1);
        Vmath::Vvtvp(nq, &m_movingframes[dir][i * nq], 1,
                     &SphericalVector[1][i * nq], 1, &m_PhicdotMF[0], 1,
                     &m_PhicdotMF[0], 1);
    }
}

void MMFSystem::ComputeMFcdotSphericalCoord(
    Array<OneD, Array<OneD, NekDouble>> &ThetacdotMF,
    Array<OneD, Array<OneD, NekDouble>> &PhicdotMF)
{
    int i, j;
    int nq = m_fields[0]->GetNpoints();

    Array<OneD, Array<OneD, NekDouble>> SphericalVector;
    ComputeSphericalVector(SphericalVector);

    for (j = 0; j < m_shapedim; ++j)
    {
        ThetacdotMF[j] = Array<OneD, NekDouble>(nq, 0.0);
        PhicdotMF[j]   = Array<OneD, NekDouble>(nq, 0.0);

        for (i = 0; i < m_spacedim; ++i)
        {
            Vmath::Vvtvp(nq, &m_movingframes[j][i * nq], 1,
                         &SphericalVector[0][i * nq], 1, &ThetacdotMF[j][0], 1,
                         &ThetacdotMF[j][0], 1);
            Vmath::Vvtvp(nq, &m_movingframes[j][i * nq], 1,
                         &SphericalVector[1][i * nq], 1, &PhicdotMF[j][0], 1,
                         &PhicdotMF[j][0], 1);
        }
    }
}

void MMFSystem::ComputeMFcdotSphericalCoord(
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    Array<OneD, Array<OneD, NekDouble>> &ThetacdotMF,
    Array<OneD, Array<OneD, NekDouble>> &PhicdotMF)
{
    int i, j;
    int nq = m_fields[0]->GetNpoints();

    Array<OneD, Array<OneD, NekDouble>> SphericalVector;
    ComputeSphericalVector(SphericalVector);

    for (j = 0; j < m_shapedim; ++j)
    {
        ThetacdotMF[j] = Array<OneD, NekDouble>(nq, 0.0);
        PhicdotMF[j]   = Array<OneD, NekDouble>(nq, 0.0);

        for (i = 0; i < m_spacedim; ++i)
        {
            Vmath::Vvtvp(nq, &movingframes[j][i * nq], 1,
                         &SphericalVector[0][i * nq], 1, &ThetacdotMF[j][0], 1,
                         &ThetacdotMF[j][0], 1);
            Vmath::Vvtvp(nq, &movingframes[j][i * nq], 1,
                         &SphericalVector[1][i * nq], 1, &PhicdotMF[j][0], 1,
                         &PhicdotMF[j][0], 1);
        }
    }
}


// Compute \nabla \cdot \sigma \vec{v} = \nabla^{\sigma} \cdot \vec{v}
Array<OneD, NekDouble> MMFSystem::ComputeMFDivergence(
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    const Array<OneD, const Array<OneD, NekDouble>> &velocity)
{
    int nq = m_fields[0]->GetNpoints();

    Array<OneD, Array<OneD, NekDouble>> inarrayMF(m_mfdim);
    for (int k = 0; k < m_mfdim; ++k)
    {
        inarrayMF[k] = Array<OneD, NekDouble>(nq);
    }

    Cart_to_MF(velocity, movingframes, inarrayMF);

    Array<OneD, NekDouble> outarray(nq, 0.0);
    Array<OneD, NekDouble> Dtmp(nq);
    for (int j = 0; j < m_expdim; ++j)
    {
        MMFDirectionalDeriv(movingframes[j], inarrayMF[j], Dtmp);
        Vmath::Vadd(nq, Dtmp, 1, outarray, 1, outarray, 1);
    }

    return outarray;
}


Array<OneD, NekDouble> MMFSystem::ComputeEuclideanDivergence(
    const Array<OneD, const Array<OneD, NekDouble>> &velocity)
{
    int nq = m_fields[0]->GetNpoints();

    Array<OneD, NekDouble> outarray(nq, 0.0);
    Array<OneD, NekDouble> Dtmp(nq);
    for (int k = 0; k < m_spacedim; ++k)
    {
        m_fields[0]->PhysDeriv(k, velocity[k], Dtmp);
        Vmath::Vadd(nq, Dtmp, 1, outarray, 1, outarray, 1);
    }

    return outarray;
}

// Compute Riemann Curvature: R(X,Y)Z =  \nabla_X \nabla_Y Z - \nabla_Y \nabla_X
// Z - \nabla_[X,Y] Z
void MMFSystem::ComputeRiemCrv(
    const int dirX, const int dirY, const int dirZ,
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    Array<OneD, Array<OneD, NekDouble>> &outarray)
{
    int nq = m_fields[0]->GetNpoints();

    Array<OneD, Array<OneD, NekDouble>> SecCovXYZ(m_shapedim);
    Array<OneD, Array<OneD, NekDouble>> SecCovYXZ(m_shapedim);
    outarray = Array<OneD, Array<OneD, NekDouble>>(m_shapedim);
    for (int j = 0; j < m_shapedim; ++j)
    {
        SecCovXYZ[j] = Array<OneD, NekDouble>(nq, 0.0);
        SecCovYXZ[j] = Array<OneD, NekDouble>(nq, 0.0);
        outarray[j]  = Array<OneD, NekDouble>(nq, 0.0);
    }

    ComputeSecondCovDeriv(dirX, dirY, dirZ, movingframes, SecCovXYZ);
    ComputeSecondCovDeriv(dirY, dirX, dirZ, movingframes, SecCovYXZ);

    for (int j = 0; j < m_shapedim; ++j)
    {
        Vmath::Vsub(nq, &SecCovXYZ[j][0], 1, &SecCovYXZ[j][0], 1,
                    &outarray[j][0], 1);
    }

    // Compute LieBracket term: \nabla_[X,Y] Z
    Array<OneD, Array<OneD, NekDouble>> LieBracket(m_shapedim);
    for (int j = 0; j < m_shapedim; ++j)
    {
        LieBracket[j] = Array<OneD, NekDouble>(nq, 0.0);
    }

    ComputeLieBracket(dirX, dirY, dirZ, movingframes, LieBracket);
    for (int j = 0; j < m_shapedim; ++j)
    {
        Vmath::Vsub(nq, &outarray[j][0], 1, &LieBracket[j][0], 1,
                    &outarray[j][0], 1);
    }
}

void MMFSystem::ComputeRelacc(
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    Array<OneD, Array<OneD, NekDouble>> &outarray)
{
    int nq = m_fields[0]->GetNpoints();

    outarray = Array<OneD, Array<OneD, NekDouble>>(m_shapedim);
    for (int j = 0; j < m_shapedim; ++j)
    {
        outarray[j] = Array<OneD, NekDouble>(nq);
    }

    Array<OneD, Array<OneD, NekDouble>> SecondCovDeriv211(m_shapedim);
    Array<OneD, Array<OneD, NekDouble>> SecondCovDeriv121(m_shapedim);
    for (int j = 0; j < m_shapedim; ++j)
    {
        SecondCovDeriv121[j] = Array<OneD, NekDouble>(nq, 0.0);
        SecondCovDeriv211[j] = Array<OneD, NekDouble>(nq, 0.0);
    }

    ComputeSecondCovDeriv(1, 0, 0, movingframes, SecondCovDeriv211);
    ComputeSecondCovDeriv(0, 1, 0, movingframes, SecondCovDeriv121);

    for (int j = 0; j < m_shapedim; ++j)
    {
        Vmath::Vsub(nq, &SecondCovDeriv211[j][0], 1, &SecondCovDeriv121[j][0],
                    1, &outarray[j][0], 1);
    }
}

void MMFSystem::ComputeRelaccOmega(
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    Array<OneD, Array<OneD, NekDouble>> &outarray)
{
    int nq = m_fields[0]->GetNpoints();

    outarray = Array<OneD, Array<OneD, NekDouble>>(m_shapedim);
    for (int j = 0; j < m_shapedim; ++j)
    {
        outarray[j] = Array<OneD, NekDouble>(nq);
    }

    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> MFConnection(m_mfdim);
    for (int i = 0; i < m_mfdim; i++)
    {
        MFConnection[i] = Array<OneD, Array<OneD, NekDouble>>(m_mfdim);
        for (int j = 0; j < m_mfdim; j++)
        {
            MFConnection[i][j] = Array<OneD, NekDouble>(nq, 0.0);
        }
    }

    // Check the Curvature 1-form of the aligned moving frames
    Compute2DConnection1form(movingframes, MFConnection);

    Array<OneD, NekDouble> w211(nq);
    Array<OneD, NekDouble> w212(nq);

    Vmath::Vcopy(nq, &MFConnection[0][0][0], 1, &w211[0], 1);
    Vmath::Vcopy(nq, &MFConnection[0][1][0], 1, &w212[0], 1);

    std::cout << "w211 = " << RootMeanSquare(w211, m_MMFActivation)
              << ", w212 = " << RootMeanSquare(w212, m_MMFActivation)
              << std::endl;

    Array<OneD, NekDouble> tmp(nq);

    // w_{211} * w_{212} - \nabla_{e_1} w_{211}
    Vmath::Vmul(nq, &w211[0], 1, &w212[0], 1, &outarray[0][0], 1);
    MMFDirectionalDeriv(movingframes[0], w211, tmp);
    Vmath::Vsub(nq, &outarray[0][0], 1, &tmp[0], 1, &outarray[0][0], 1);

    // w_{212} * w_{212} - \nabla_{e_1} w_{212}
    Vmath::Vmul(nq, &w212[0], 1, &w212[0], 1, &outarray[1][0], 1);
    MMFDirectionalDeriv(movingframes[0], w212, tmp);
    Vmath::Vsub(nq, &outarray[1][0], 1, &tmp[0], 1, &outarray[1][0], 1);
}

// Array<OneD, NekDouble> MMFSystem::ComputeRelAcc(
//     const Array<OneD, const Array<OneD, NekDouble>> &movingframes)
// {
//     int nq = m_fields[0]->GetNpoints();

//     Array<OneD, NekDouble> outarray(nq);

//     Array<OneD, Array<OneD, NekDouble>> SecondCovDeriv211(m_shapedim);
//     Array<OneD, Array<OneD, NekDouble>> SecondCovDeriv121(m_shapedim);
//     for (int j=0; j<m_shapedim; ++j)
//     {
//         SecondCovDeriv121[j] = Array<OneD, NekDouble>(nq,0.0);
//         SecondCovDeriv211[j] = Array<OneD, NekDouble>(nq,0.0);
//     }

//     ComputeSecondCovDeriv(0, 0, 1, movingframes, SecondCovDeriv);
//     Vmath::Vcopy(nq, &SecondCovDeriv[1][0], 1, &outarray[0], 1);

//     return outarray;
// }

// Compute Covariant derivative; \nabla_{X} \nabla_{Y} Z
// Compute \nabla_X Y
// \nabla_{e_1} e_2 = - \omega_{21} (e_1) e_1 - \omega_{21} (e_2) e_2
// \nabla_{e_2} e_1 =   \omega_{21} (e_1) e_1 + \omega_{21} (e_2) e_2
void MMFSystem::ComputeSecondCovDeriv(
    const int dirX, const int dirY, const int dirZ,
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    Array<OneD, Array<OneD, NekDouble>> &outarray)
{
    int nq = m_fields[0]->GetNpoints();

    outarray = Array<OneD, Array<OneD, NekDouble>>(m_shapedim);
    for (int j = 0; j < m_shapedim; ++j)
    {
        outarray[j] = Array<OneD, NekDouble>(nq, 0.0);
    }

    Array<OneD, Array<OneD, NekDouble>> CovDeriv(m_shapedim);
    Array<OneD, Array<OneD, NekDouble>> CovDerivX1(m_shapedim);
    Array<OneD, Array<OneD, NekDouble>> CovDerivX2(m_shapedim);
    for (int j = 0; j < m_shapedim; ++j)
    {
        CovDeriv[j]   = Array<OneD, NekDouble>(nq, 0.0);
        CovDerivX1[j] = Array<OneD, NekDouble>(nq, 0.0);
        CovDerivX2[j] = Array<OneD, NekDouble>(nq, 0.0);
    }

    // CovDeriv = Compute \nabla_{Y} Z = CovDeriv[0] e^1 + CovDeriv[1] e^2
    ComputeCovDerivMF(dirY, dirZ, movingframes, CovDeriv);

    // CovDerivX1 = Compute \nabla_{X} e^0
    ComputeCovDerivMF(dirX, 0, movingframes, CovDerivX1);

    // CovDerivX2 = Compute \nabla_{X} e^1
    ComputeCovDerivMF(dirX, 1, movingframes, CovDerivX2);

    // std::cout << "X = " << dirX << ", Y = " << dirY << ", dirZ = " << dirZ <<
    // std::endl; std::cout << "CovDeriv = ( " << RootMeanSquare(CovDeriv[0]) <<
    // " , " << RootMeanSquare(CovDeriv[1]) << " ) " << std::endl;

    Array<OneD, NekDouble> DerivXa(nq);
    Array<OneD, NekDouble> DerivXb(nq);

    // Compute \nabla_{e_X} CovDeriv[0]
    MMFDirectionalDeriv(movingframes[dirX], CovDeriv[0], DerivXa);

    // Compute \nabla_{e_X} CovDeriv[1]
    MMFDirectionalDeriv(movingframes[dirX], CovDeriv[1], DerivXb);

    for (int i = 0; i < nq; ++i)
    {
        outarray[0][i] = DerivXa[i] + CovDeriv[0][i] * CovDerivX1[0][i] +
                         CovDeriv[1][i] * CovDerivX2[0][i];
        outarray[1][i] = DerivXb[i] + CovDeriv[0][i] * CovDerivX1[1][i] +
                         CovDeriv[1][i] * CovDerivX2[1][i];
    }
}

// j=0, i=0:  \nabla_{e_1} \vec{u} \cdot e_1 = \nabla u^1 \cdot e_1 -
// \omega_{21} (e_1) u^2 j=0, i=nq: \nabla_{e_1} \vec{u} \cdot e_2 = \nabla u^2
// \cdot e_1 + \omega_{21} (e_1) u^1 j=1, i=0:  \nabla_{e_2} \vec{u} \cdot e_1 =
// \nabla u^1 \cdot e_2 - \omega_{21} (e_2) u^2 j=1, i=nq:  \nabla_{e_2} \vec{u}
// \cdot e_2 = \nabla u^2 \cdot e_2 + \omega_{21} (e_2) u^1

// j=0, i=0: u1=0, u2=1; \nabla_{e_1} e_1 = -\omega_{21} (e_1) e_1
// j=1, i=0: u1=0, u2=1; \nabla_{e_2} e_1 = -\omega_{21} (e_2) e_1

// j=0, i=nq: u1=1, u2=0; \nabla_{e_1} e_2 = \omega_{21} (e_1) e_2
// j=1, i=nq: u1=1, u2=0; \nabla_{e_2} e_2 = \omega_{21} (e_2) e_2

// void MMFSystem::ComputeCovDerivMF(
//     const Array<OneD, const NekDouble> &u1,
//     const Array<OneD, const NekDouble> &u2,
//     const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
//     Array<OneD, Array<OneD, NekDouble>> &outarray)
void MMFSystem::ComputeLieBracket(
    const int dirX, const int dirY, const int dirZ,
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    Array<OneD, Array<OneD, NekDouble>> &outarray)
{
    int nq = m_fields[0]->GetNpoints();

    Array<OneD, Array<OneD, NekDouble>> coeffXY(m_shapedim);
    Array<OneD, Array<OneD, NekDouble>> coeffYX(m_shapedim);
    Array<OneD, Array<OneD, NekDouble>> coeffs(m_shapedim);
    for (int j = 0; j < m_shapedim; ++j)
    {
        coeffXY[j] = Array<OneD, NekDouble>(nq, 0.0);
        coeffYX[j] = Array<OneD, NekDouble>(nq, 0.0);
        coeffs[j]  = Array<OneD, NekDouble>(nq, 0.0);
    }

    // Compute [X, Y] = \nabla_X Y - \nabla_Y X
    ComputeCovDerivMF(dirX, dirY, movingframes, coeffXY);
    ComputeCovDerivMF(dirY, dirX, movingframes, coeffYX);
    for (int j = 0; j < m_shapedim; ++j)
    {
        Vmath::Vsub(nq, &coeffXY[j][0], 1, &coeffYX[j][0], 1, &coeffs[j][0], 1);
    }

    // std::cout << "[dixX, dirY] = ( " << RootMeanSquare(coeffs[0]) << " , " <<
    // RootMeanSquare(coeffs[1]) << " ) " << std::endl;

    Array<OneD, Array<OneD, NekDouble>> CovDeriv1K(m_shapedim);
    Array<OneD, Array<OneD, NekDouble>> CovDeriv2K(m_shapedim);
    for (int j = 0; j < m_shapedim; ++j)
    {
        CovDeriv1K[j] = Array<OneD, NekDouble>(nq, 0.0);
        CovDeriv2K[j] = Array<OneD, NekDouble>(nq, 0.0);
    }

    // Compute \nabla_1 Z and \nabla_2 Z
    ComputeCovDerivMF(0, dirZ, movingframes, CovDeriv1K);
    ComputeCovDerivMF(1, dirZ, movingframes, CovDeriv2K);

    // std::cout << "CovDeriv1K = ( " << RootMeanSquare(CovDeriv1K[0]) << " , "
    // << RootMeanSquare(CovDeriv1K[1]) << " ) " << std::endl; std::cout <<
    // "CovDeriv2K = ( " << RootMeanSquare(CovDeriv2K[0]) << " , " <<
    // RootMeanSquare(CovDeriv2K[1]) << " ) " << std::endl;

    outarray = Array<OneD, Array<OneD, NekDouble>>(m_shapedim);
    for (int j = 0; j < m_shapedim; ++j)
    {
        outarray[j] = Array<OneD, NekDouble>(nq, 0.0);
    }

    // a \nabla_{e_1} e_K + b \nabla_{e_2} e_K
    for (int j = 0; j < m_shapedim; ++j)
    {
        Vmath::Vvtvp(nq, &coeffs[0][0], 1, &CovDeriv1K[j][0], 1,
                     &outarray[j][0], 1, &outarray[j][0], 1);
        Vmath::Vvtvp(nq, &coeffs[1][0], 1, &CovDeriv2K[j][0], 1,
                     &outarray[j][0], 1, &outarray[j][0], 1);
    }
}

// Compute Covariant derivative; \nabla_{ \nabla_{direction1} direction2
// }\vec{u}
void MMFSystem::ComputeLieBracket(
    const Array<OneD, const Array<OneD, NekDouble>> &velocity,
    const Array<OneD, const NekDouble> &direction1,
    const Array<OneD, const NekDouble> &direction2,
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    Array<OneD, Array<OneD, NekDouble>> &outarray)
{
    int nq = m_fields[0]->GetNpoints();

    outarray = Array<OneD, Array<OneD, NekDouble>>(m_shapedim);
    for (int j = 0; j < m_shapedim; ++j)
    {
        outarray[j] = Array<OneD, NekDouble>(nq, 0.0);
    }

    Array<OneD, Array<OneD, NekDouble>> CovDeriv(m_shapedim);
    for (int j = 0; j < m_shapedim; ++j)
    {
        CovDeriv[j] = Array<OneD, NekDouble>(nq, 0.0);
    }

    // Compute \nabla_{direction1} direction2
    Array<OneD, Array<OneD, NekDouble>> direction2Cart(m_spacedim);
    for (int k = 0; k < m_spacedim; ++k)
    {
        direction2Cart[k] = Array<OneD, NekDouble>(nq, 0.0);
        Vmath::Vcopy(nq, &direction2[k * nq], 1, &direction2Cart[k][0], 1);
    }

    // \nabla_{direction1} direction2 = a e_1 + b e_2
    ComputeCovDeriv(direction2Cart, direction1, movingframes, CovDeriv);

    // \nabla_{a e_1 + b e_2 } \vec{u} = a \nabla_1 \vec{u} + b \nabla_2 \vec{u}
    Array<OneD, NekDouble> u1(nq);
    Array<OneD, NekDouble> u2(nq);

    CartesianToMovingframes(velocity[0], velocity[1], velocity[2], movingframes,
                            u1, u2);

    Array<OneD, Array<OneD, NekDouble>> CovDerivMF(m_shapedim);
    for (int j = 0; j < m_shapedim; ++j)
    {
        CovDerivMF[j] = Array<OneD, NekDouble>(m_shapedim * nq);
    }

    ComputeCovDeriv(u1, u2, movingframes, CovDerivMF);

    // e_1 component
    // a( \nabla_{e_1} \vec{u} \cdot e_1 ) + b ( \nabla_{e_2} \vec{u} \cdot e_1
    // )
    Vmath::Vvtvp(nq, &CovDeriv[0][0], 1, &CovDerivMF[0][0 * nq], 1,
                 &outarray[0][0], 1, &outarray[0][0], 1);
    Vmath::Vvtvp(nq, &CovDeriv[1][0], 1, &CovDerivMF[1][0 * nq], 1,
                 &outarray[0][0], 1, &outarray[0][0], 1);

    // e_2 component
    // a( \nabla_{e_1} \vec{u} \cdot e_2 ) + b ( \nabla_{e_2} \vec{u} \cdot e_2
    // )
    Vmath::Vvtvp(nq, &CovDeriv[0][0], 1, &CovDerivMF[0][1 * nq], 1,
                 &outarray[1][0], 1, &outarray[1][0], 1);
    Vmath::Vvtvp(nq, &CovDeriv[1][0], 1, &CovDerivMF[1][1 * nq], 1,
                 &outarray[1][0], 1, &outarray[1][0], 1);
}

// Compute Covariant derivative; \nabla_{vec{v}} \vec{u}
// j=0, i=0:  \nabla_{e_1} \vec{u} \cdot e_1 = \nabla u^1 \cdot e_1 -
// \omega_{21} (e_1) u^2 j=0, i=nq: \nabla_{e_1} \vec{u} \cdot e_2 = \nabla u^2
// \cdot e_1 + \omega_{21} (e_1) u^1 j=1, i=0:  \nabla_{e_2} \vec{u} \cdot e_1 =
// \nabla u^1 \cdot e_2 - \omega_{21} (e_2) u^2 j=1, i=nq:  \nabla_{e_2} \vec{u}
// \cdot e_2 = \nabla u^2 \cdot e_2 + \omega_{21} (e_2) u^1
void MMFSystem::ComputeCovDeriv(
    const Array<OneD, const Array<OneD, NekDouble>> &velocity,
    const Array<OneD, const NekDouble> &direction,
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    Array<OneD, Array<OneD, NekDouble>> &outarray)
{
    int nq = m_fields[0]->GetNpoints();

    outarray = Array<OneD, Array<OneD, NekDouble>>(m_shapedim);
    for (int j = 0; j < m_shapedim; ++j)
    {
        outarray[j] = Array<OneD, NekDouble>(nq, 0.0);
    }

    Array<OneD, Array<OneD, NekDouble>> dirvector(m_spacedim);
    for (int k = 0; k < m_spacedim; ++k)
    {
        dirvector[k] = Array<OneD, NekDouble>(nq);
        Vmath::Vcopy(nq, &direction[k * nq], 1, &dirvector[k][0], 1);
    }

    Array<OneD, Array<OneD, NekDouble>> dirv(m_shapedim);
    for (int j = 0; j < m_shapedim; ++j)
    {
        dirv[j] = Array<OneD, NekDouble>(nq);
    }

    CartesianToMovingframes(dirvector[0], dirvector[1], dirvector[2],
                            movingframes, dirv[0], dirv[1]);

    Array<OneD, Array<OneD, NekDouble>> CovDeriv(m_shapedim);
    for (int j = 0; j < m_shapedim; ++j)
    {
        CovDeriv[j] = Array<OneD, NekDouble>(m_shapedim * nq);
    }

    Array<OneD, NekDouble> u1(nq);
    Array<OneD, NekDouble> u2(nq);

    CartesianToMovingframes(velocity[0], velocity[1], velocity[2], movingframes,
                            u1, u2);
    ComputeCovDeriv(u1, u2, movingframes, CovDeriv);

    for (int j = 0; j < m_shapedim; ++j)
    {
        Vmath::Vmul(nq, &dirv[0][0], 1, &CovDeriv[0][j * nq], 1,
                    &outarray[j][0], 1);
        Vmath::Vvtvp(nq, &dirv[1][0], 1, &CovDeriv[1][j * nq], 1,
                     &outarray[j][0], 1, &outarray[j][0], 1);
    }
}

// Compute Covariant derivative; \nabla_{vec{v}} \vec{u}
// j=0, i=0:  \nabla_{e_1} \vec{u} \cdot e_1 = \nabla u^1 \cdot e_1 -
// \omega_{21} (e_1) u^2 j=0, i=nq: \nabla_{e_1} \vec{u} \cdot e_2 = \nabla u^2
// \cdot e_1 + \omega_{21} (e_1) u^1 j=1, i=0:  \nabla_{e_2} \vec{u} \cdot e_1 =
// \nabla u^1 \cdot e_2 - \omega_{21} (e_2) u^2 j=1, i=nq:  \nabla_{e_2} \vec{u}
// \cdot e_2 = \nabla u^2 \cdot e_2 + \omega_{21} (e_2) u^1
void MMFSystem::ComputeCovDeriv(
    const Array<OneD, const NekDouble> &u1,
    const Array<OneD, const NekDouble> &u2,
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    Array<OneD, Array<OneD, NekDouble>> &outarray)
{
    int nq = m_fields[0]->GetNpoints();

    outarray = Array<OneD, Array<OneD, NekDouble>>(m_shapedim);
    for (int j = 0; j < m_shapedim; ++j)
    {
        outarray[j] = Array<OneD, NekDouble>(m_shapedim * nq);
    }

    // CovDeriv = ComputeCovariantDerivative(0, u1, u2, movingframes);
    Array<OneD, NekDouble> tmp(nq);
    Array<OneD, NekDouble> duide(nq);

    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> MFConnection;
    Compute2DConnection1form(movingframes, MFConnection);

    // j=0:  \nabla_{e_1} \vec{u} \cdot e_1 = \nabla u^1 \cdot e_1 - \omega_{21}
    // (e_1) u^2
    MMFDirectionalDeriv(movingframes[0], u1, duide);
    Vmath::Vvtvm(nq, &MFConnection[0][0][0], 1, &u2[0], 1, &duide[0], 1,
                 &tmp[0], 1);
    Vmath::Neg(nq, tmp, 1);
    Vmath::Vcopy(nq, &tmp[0], 1, &outarray[0][0], 1);

    // j=0:  \nabla_{e_2} \vec{u} \cdot e_1 = \nabla u^1 \cdot e_2 - \omega_{21}
    // (e_2) u^2
    MMFDirectionalDeriv(movingframes[1], u1, duide);
    Vmath::Vvtvm(nq, &MFConnection[0][1][0], 1, &u2[0], 1, &duide[0], 1,
                 &tmp[0], 1);
    Vmath::Neg(nq, tmp, 1);
    Vmath::Vcopy(nq, &tmp[0], 1, &outarray[0][nq], 1);

    // j=1: \nabla_{e_1} \vec{u} \cdot e_2 = \nabla u^2 \cdot e_1 + \omega_{21}
    // (e_1) u^1
    MMFDirectionalDeriv(movingframes[0], u2, duide);
    Vmath::Vvtvp(nq, &MFConnection[0][0][0], 1, &u1[0], 1, &duide[0], 1,
                 &tmp[0], 1);
    Vmath::Vcopy(nq, &tmp[0], 1, &outarray[1][0], 1);

    // j=1:  \nabla_{e_2} \vec{u} \cdot e_2 = \nabla u^2 \cdot e_2 - \omega_{21}
    // (e_2) u^2
    MMFDirectionalDeriv(movingframes[1], u2, duide);
    Vmath::Vvtvp(nq, &MFConnection[0][1][0], 1, &u1[0], 1, &duide[0], 1,
                 &tmp[0], 1);
    Vmath::Vcopy(nq, &tmp[0], 1, &outarray[1][nq], 1);
}

// Compute Covariant derivative; \nabla_{vec{v}} \vec{u}
// j=0, i=0:  \nabla_{e_1} \vec{u} \cdot e_1 = \nabla u^1 \cdot e_1 -
// \omega_{21} (e_1) u^2 j=0, i=nq: \nabla_{e_1} \vec{u} \cdot e_2 = \nabla u^2
// \cdot e_1 + \omega_{21} (e_1) u^1 j=1, i=0:  \nabla_{e_2} \vec{u} \cdot e_1 =
// \nabla u^1 \cdot e_2 - \omega_{21} (e_2) u^2 j=1, i=nq:  \nabla_{e_2} \vec{u}
// \cdot e_2 = \nabla u^2 \cdot e_2 + \omega_{21} (e_2) u^1

// Compute \nabla_X Y
// \nabla_{e_1} e_2 = - \omega_{21} (e_1) e_1 - \omega_{21} (e_2) e_2
// \nabla_{e_2} e_1 =   \omega_{21} (e_1) e_1 + \omega_{21} (e_2) e_2

void MMFSystem::ComputeCovDerivMF(
    const int dirX, const int dirY,
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    Array<OneD, Array<OneD, NekDouble>> &outarray)
{
    int nq = m_fields[0]->GetNpoints();

    outarray = Array<OneD, Array<OneD, NekDouble>>(m_shapedim);
    for (int j = 0; j < m_shapedim; ++j)
    {
        outarray[j] = Array<OneD, NekDouble>(nq, 0.0);
    }

    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> MFConnection;
    Compute2DConnection1form(movingframes, MFConnection);

    if ((dirX == 0) && (dirY == 1))
    {
        for (int j = 0; j < m_shapedim; ++j)
        {
            Vmath::Smul(nq, -1.0, &MFConnection[0][j][0], 1, &outarray[j][0],
                        1);
        }
    }

    else if ((dirX == 1) && (dirY == 0))
    {
        for (int j = 0; j < m_shapedim; ++j)
        {
            Vmath::Vcopy(nq, &MFConnection[0][j][0], 1, &outarray[j][0], 1);
        }
    }
}

Array<OneD, NekDouble> MMFSystem::ComputeCovDiv(
    const Array<OneD, const NekDouble> &velvector,
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes)
{
    int nq = m_fields[0]->GetNpoints();

    Array<OneD, Array<OneD, NekDouble>> velocity(m_spacedim);
    for (int k = 0; k < m_spacedim; ++k)
    {
        velocity[k] = Array<OneD, NekDouble>(nq);
        Vmath::Vcopy(nq, &velvector[k * nq], 1, &velocity[k][0], 1);
    }

    Array<OneD, NekDouble> outarray(nq);
    outarray = ComputeCovDiv(velocity, movingframes);

    return outarray;
}

Array<OneD, NekDouble> MMFSystem::ComputeCovDiv(
    const Array<OneD, const Array<OneD, NekDouble>> &velocity,
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes)
{
    int nq = m_fields[0]->GetNpoints();

    Array<OneD, NekDouble> outarray(nq, 0.0);

    Array<OneD, NekDouble> u1(nq);
    Array<OneD, NekDouble> u2(nq);

    CartesianToMovingframes(velocity[0], velocity[1], velocity[2], movingframes,
                            u1, u2);

    Array<OneD, Array<OneD, NekDouble>> CovDeriv(m_shapedim);
    for (int j = 0; j < m_shapedim; ++j)
    {
        CovDeriv[j] = Array<OneD, NekDouble>(m_shapedim * nq);
    }

    ComputeCovDeriv(u1, u2, movingframes, CovDeriv);

    Vmath::Vadd(nq, &CovDeriv[0][0], 1, &outarray[0], 1, &outarray[0], 1);
    Vmath::Vadd(nq, &CovDeriv[1][nq], 1, &outarray[0], 1, &outarray[0], 1);

    return outarray;
}

Array<OneD, NekDouble> MMFSystem::ComputeCovCurl(
    const Array<OneD, const NekDouble> &velvector,
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes)
{
    int nq = m_fields[0]->GetNpoints();

    Array<OneD, Array<OneD, NekDouble>> velocity(m_spacedim);
    for (int k = 0; k < m_spacedim; ++k)
    {
        velocity[k] = Array<OneD, NekDouble>(nq);
        Vmath::Vcopy(nq, &velvector[k * nq], 1, &velocity[k][0], 1);
    }

    Array<OneD, NekDouble> outarray(nq);
    outarray = ComputeCovCurl(velocity, movingframes);

    return outarray;
}

Array<OneD, NekDouble> MMFSystem::ComputeCovCurl(
    const Array<OneD, const Array<OneD, NekDouble>> &velocity,
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes)
{
    int nq = m_fields[0]->GetNpoints();

    Array<OneD, NekDouble> outarray(nq, 0.0);

    Array<OneD, NekDouble> u1(nq);
    Array<OneD, NekDouble> u2(nq);

    CartesianToMovingframes(velocity[0], velocity[1], velocity[2], movingframes,
                            u1, u2);

    Array<OneD, Array<OneD, NekDouble>> CovDerivMF(m_shapedim);
    for (int j = 0; j < m_shapedim; ++j)
    {
        CovDerivMF[j] = Array<OneD, NekDouble>(m_shapedim * nq);
    }

    ComputeCovDeriv(u1, u2, movingframes, CovDerivMF);

    // Vmath::Vadd(nq, &CovDerivMF[1][0], 1, &outarray[0], 1, &outarray[0], 1);
    // Vmath::Vsub(nq, &outarray[0], 1, &CovDerivMF[0][nq], 1, &outarray[0], 1);

    Vmath::Vsub(nq, &CovDerivMF[1][0], 1, &CovDerivMF[0][nq], 1, &outarray[0],
                1);

    return outarray;
}

// Array<OneD, NekDouble> MMFSystem::ComputeCovariantDivergence(
//     const Array<OneD, const NekDouble> &vellong,
//     const Array<OneD, const Array<OneD, NekDouble>> &movingframes)
// {
//     int nq = m_fields[0]->GetNpoints();

//     Array<OneD, Array<OneD, NekDouble>> velocity(m_spacedim);
//     for (int k = 0; k < m_spacedim; ++k)
//     {
//         velocity[k] = Array<OneD, NekDouble>(nq);
//         Vmath::Vcopy(nq, &vellong[k * nq], 1, &velocity[k][0], 1);
//     }

//     Array<OneD, NekDouble> outarray(nq);
//     outarray = ComputeCovariantDivergence(velocity, movingframes);

//     return outarray;
// }

// \nabla \cdot \vec{v} = \nabla v_1 \cdot \mathbf{e}^1 - \Gamma^2_{11} v_2
//                          + \nabla v_2 \cdot \mathbf{e}^2 + \Gamma^2_{21} v_1
// Array<OneD, NekDouble> MMFSystem::ComputeCovariantDivergence(
//     const Array<OneD, const Array<OneD, NekDouble>> &velocity,
//     const Array<OneD, const Array<OneD, NekDouble>> &movingframes)
// {
//     int nq = m_fields[0]->GetNpoints();

//     Array<OneD, NekDouble> u1(nq);
//     Array<OneD, NekDouble> u2(nq);

//     CartesianToMovingframes(velocity[0], velocity[1], velocity[2],
//     movingframes,
//                             u1, u2);

//     Array<OneD, NekDouble> du1(nq);
//     Array<OneD, NekDouble> du2(nq);

//     m_fields[0]->PhysDirectionalDeriv(movingframes[0], u1, du1);
//     m_fields[0]->PhysDirectionalDeriv(movingframes[1], u2, du2);

//     Array<OneD, NekDouble> outarray(nq, 0.0);

//     Vmath::Vadd(nq, du1, 1, du2, 1, outarray, 1);

//     // Add Christoffel symbol
//     Array<OneD, Array<OneD, Array<OneD, NekDouble>>> MFConnection;
//     ComputeConnection1form(movingframes, MFConnection);

//     // \nabla \cdot e^1 = \Gamma^2_{12} = w^2_1 <e^2>
//     // Vmath::Vcopy(nq, &MFConnection[0][1][0], 1, &DivMF[0][0], 1);
//     // Vmath::Vvtvp(nq, &MFConnection[0][1][0], 1, &u1[0], 1, &outarray[0],
//     1,
//     //          &outarray[0], 1);
//     Array<OneD, NekDouble> Gammau1(nq);

//     Vmath::Vmul(nq, &MFConnection[0][1][0], 1, &u1[0], 1, &Gammau1[0], 1);
//     Vmath::Vadd(nq, &Gammau1[0], 1, &outarray[0], 1, &outarray[0], 1);

//     // \nabla \cdot e^2 = \Gamma^1_{21} = -w^2_1 <e^1>
//     // Vmath::Vcopy(nq, &MFConnection[0][0][0], 1, &DivMF[1][0], 1);
//     Array<OneD, NekDouble> Gammau2(nq);
//     Vmath::Vmul(nq, &MFConnection[0][0][0], 1, &u2[0], 1, &Gammau2[0], 1);
//     Vmath::Neg(nq, Gammau2, 1);

//     Vmath::Vadd(nq, &Gammau2[0], 1, &outarray[0], 1, &outarray[0], 1);

//     return outarray;
// }

// Compute $ \mathbf{k} \cdot ( \nabla \times \vec{v} )$
// \nabla \times \vec{v} = \frac{\partial v_2}{\partial x_1} +
// \Gamma^2_{11} v^1
//                         - (\frac{\partial v_1}{\partial x_2} +
//                         \Gamma^1_{22} v^2)
// Array<OneD, NekDouble> MMFSystem::ComputeCovariantCurl(
//     const Array<OneD, const NekDouble> &vellong,
//     const Array<OneD, const Array<OneD, NekDouble>> &movingframes)
// {
//     int nq = m_fields[0]->GetNpoints();

//     Array<OneD, Array<OneD, NekDouble>> velocity(m_spacedim);
//     for (int k = 0; k < m_spacedim; ++k)
//     {
//         velocity[k] = Array<OneD, NekDouble>(nq);
//         Vmath::Vcopy(nq, &vellong[k * nq], 1, &velocity[k][0], 1);
//     }

//     Array<OneD, NekDouble> rval(nq);
//     rval = ComputeCovariantCurl(velocity, movingframes);

//     return rval;
// }

// Array<OneD, NekDouble> MMFSystem::ComputeCovariantCurl(
//     const Array<OneD, const Array<OneD, NekDouble>> &velocity,
//     const Array<OneD, const Array<OneD, NekDouble>> &movingframes)
// {
//     int nq = m_fields[0]->GetNpoints();

//     Array<OneD, NekDouble> u1(nq);
//     Array<OneD, NekDouble> u2(nq);

//     Array<OneD, NekDouble> du1dx2(nq);
//     Array<OneD, NekDouble> du2dx1(nq);

//     Array<OneD, NekDouble> outarray(nq, 0.0);

//     CartesianToMovingframes(velocity[0], velocity[1], velocity[2],
//     movingframes,
//                             u1, u2);

//     // Add Christoffel symbol
//     Array<OneD, Array<OneD, Array<OneD, NekDouble>>> MFConnection;
//     ComputeConnection1form(movingframes, MFConnection);

//     Array<OneD, NekDouble> tmp(nq);

//     // du1dx2 = Compute \partial u2 / /partial x^1 + \Gamma^2_{12} u2
//     m_fields[0]->PhysDirectionalDeriv(movingframes[0], u2, du2dx1);

//     Vmath::Vmul(nq, &MFConnection[0][1][0], 1, &u2[0], 1, &tmp[0], 1);
//     Vmath::Vadd(nq, &du2dx1[0], 1, &tmp[0], 1, &du2dx1[0], 1);

//     // du2dx1 = Compute \partial u2 / /partial x^1 - \Gamma^2_{11} u1
//     m_fields[0]->PhysDirectionalDeriv(movingframes[1], u1, du1dx2);

//     Vmath::Vmul(nq, &MFConnection[0][0][0], 1, &u1[0], 1, &tmp[0], 1);
//     Vmath::Vsub(nq, &du1dx2[0], 1, &tmp[0], 1, &du1dx2[0], 1);

//     // \nabla \times \vec{u} = du1dx2 - du2dx1
//     Vmath::Vsub(nq, du2dx1, 1, du1dx2, 1, outarray, 1);

//     return outarray;
// }

// Compute Covariant derivative; \nabla_{vec{v}} \vec{u}
// Array<OneD, NekDouble> MMFSystem::ComputeCovariantDerivative(
//     const Array<OneD, const Array<OneD, NekDouble>> &velocity,
//     const Array<OneD, const NekDouble> &direction,
//     const Array<OneD, const Array<OneD, NekDouble>> &movingframes)
//     {
//         int nq = m_fields[0]->GetNpoints();

//         Array<OneD, NekDouble> outarray(m_shapedim * nq, 0.0);

//         Array<OneD, NekDouble> CovDeriv(m_shapedim * nq, 0.0);

//         Array<OneD, NekDouble> u1(nq);
//         Array<OneD, NekDouble> u2(nq);

//         CartesiantoMovingframes(velocity[0], velocity[1], velocity[2],
//         movingframes, u1, u2);

//         Array<OneD, NekDouble> uth(nq);
//         Array<OneD, NekDouble> uphi(nq);

//         MovingframesToSpherical(u1, u2, movingframes, uth, uphi);

//         std::cout << "Velocity: u1 = ( " << RootMeanSquare(u1) << ", u2 = "
//         << RootMeanSquare(u2)
//         << ", uth = " << RootMeanSquare(uth) << ", uphi = " <<
//         RootMeanSquare(uphi) << std::endl;

//         Array<OneD, Array<OneD, NekDouble>> dirvector(m_spacedim);
//         for (int k=0; k<m_spacedim; ++k)
//         {
//             dirvector[k] = Array<OneD, NekDouble>(nq);
//             Vmath::Vcopy(nq, &direction[k*nq], 1, &dirvector[k][0], 1);
//         }

//         Array<OneD, Array<OneD, NekDouble>> dirv(m_shapedim);
//         for (int j=0; j<m_shapedim; ++j)
//         {
//             dirv[j] = Array<OneD, NekDouble>(nq);
//         }
//         Array<OneD, NekDouble> dirv2(nq);

//         Array<OneD, NekDouble> dirvth(nq);
//         Array<OneD, NekDouble> dirvphi(nq);

//         CartesiantoMovingframes(dirvector[0], dirvector[1], dirvector[2],
//         movingframes, dirv[0], dirv[1]); MovingframesToSpherical(dirv[0],
//         dirv[1], movingframes, dirvth, dirvphi);

//         std::cout << "DirDeriv: vel = ( " << RootMeanSquare(dirv[0]) << " , "
//         << RootMeanSquare(dirv[1])
//         << " ), velSph = ( " << RootMeanSquare(dirvth) << " , " <<
//         RootMeanSquare(dirvphi) << " ) " << std::endl;

//         // \nabla_{e_j} \vec{u}
//         Array<OneD, NekDouble> tmp1(nq);
//         Array<OneD, NekDouble> tmp2(nq);
//         for (int j=0; j<m_shapedim; ++j)
//         {
//             // j=0:
//             // tmp1 = \nabla u^1 \cdot \mathbf{e}_1
//             // tmp2 = \nabla u^2 \cdot \mathbf{e}_1 - \omega_{21}
//             (\mathbf{e}_2) u^2

//             // j = 1:
//             // tmp1 = \nabla u^1 \cdot \mathbf{e}_2 + \omega_{21}
//             (\mathbf{e}_2) u^2
//             // tmp2 = \nabla u^2 \cdot \mathbf{e}_2

//             CovDeriv = ComputeCovariantDerivative(j, u1, u2, movingframes);

//             Vmath::Vcopy(nq, &CovDeriv[0], 1, &tmp1[0], 1);
//             Vmath::Vcopy(nq, &CovDeriv[nq], 1, &tmp2[0], 1);

//             std::cout << "j = " << j << ", ComputeCovDeriv = ( " <<
//             RootMeanSquare(tmp1, m_MMFActivation) << " , " <<
//             RootMeanSquare(tmp2, m_MMFActivation) << " ) " << std::endl;

//         // \nabla_{vec{v}} \vec{u} = v^1 \nabla_{e_1} \vec{u}  + v^2
//         \nabla_{e_2} \vec{u}

//             Vmath::Vvtvp(nq, &dirv[j][0], 1, &CovDeriv[0], 1, &outarray[0],
//             1, &outarray[0], 1); Vmath::Vvtvp(nq, &dirv[j][0], 1,
//             &CovDeriv[nq], 1, &outarray[nq], 1, &outarray[nq], 1);

//           //  Vmath::Vmul(nq, &dirv1[0], 1, &CovDeriv[0][j*nq], 1,
//           &outarray[j*nq], 1);
//           //  Vmath::Vvtvp(nq, &dirv2[0], 1, &CovDeriv[1][j*nq], 1,
//           &outarray[j*nq], 1, &outarray[j*nq], 1);
//         }

//         return outarray;
//     }

// Compute Covariant derivative; \nabla_{e_i} \vec{u}
// input: direction i
// input: velocity: \vec{u} = u^1 \mathbf{e}_1 + u^2 \mathbf{e}_2
// input: movingframes: e_1
// output: outarray (m_shapedim * nq): ( \mathbf{e}_1 component, \mathbf{e}_2
// component ) Array<OneD, NekDouble> MMFSystem::ComputeCovariantDerivative(
//     const int direction,
//     const Array<OneD, const NekDouble> &ui,
//     const Array<OneD, const NekDouble> &u2,
//     const Array<OneD, const Array<OneD, NekDouble>> &movingframes)
// {
//     int nq = m_fields[0]->GetNpoints();

//     // return vector
//     Array<OneD, NekDouble> outarray(m_shapedim * nq);

//     Array<OneD, NekDouble> du1dei(nq);
//     Array<OneD, NekDouble> du2dei(nq);

//     m_fields[0]->PhysDirectionalDeriv(movingframes[0], ui, du1dei);
//     m_fields[0]->PhysDirectionalDeriv(movingframes[1], ui, du2dei);

//     // Array<OneD, Array<OneD, NekDouble>> ThetacdotMF(m_shapedim);
//     // Array<OneD, Array<OneD, NekDouble>> PhicdotMF(m_shapedim);

//     // ComputeMFcdotSphericalCoord(movingframes, ThetacdotMF, PhicdotMF);

//     // std::cout << "Theta * e1 = " << RootMeanSquare(ThetacdotMF[0],
//     m_MMFActivation) << ", Theta * e2 = " << RootMeanSquare(ThetacdotMF[1],
//     m_MMFActivation)
//     //    << ", Phi * e1 = " << RootMeanSquare(PhicdotMF[0], m_MMFActivation)
//     << ", Phi * e2 = " << RootMeanSquare(PhicdotMF[1], m_MMFActivation) <<
//     std::endl;

//     // Array<OneD, NekDouble> du1dth(nq);
//     // Array<OneD, NekDouble> du2dphi(nq);
//     // for (int i=0; i<nq; ++i)
//     // {
//     //     du1dth[i] = du1dei[i]*ThetacdotMF[0][i] +
//     du2dei[i]*ThetacdotMF[1][i];
//     //     du2dphi[i] = du1dei[i]*PhicdotMF[0][i]  +
//     du2dei[i]*PhicdotMF[1][i];
//     // }

//     // std::cout << "du1dei = " << RootMeanSquare(du1dei, m_MMFActivation) <<
//     ", du2dei = " << RootMeanSquare(du2dei, m_MMFActivation)
//     // << ". du1dth = " << RootMeanSquare(du1dth, m_MMFActivation) << ",
//     du2dphi = " << RootMeanSquare(du2dphi, m_MMFActivation) << std::endl;

//     // Add Connection form
//     // w_{21} <e^1> = Connectionform[0][0];
//     // w_{21} <e^2> = Connectionform[0][1];
//     // w_{21} <e^3> = Connectionform[0][2];
//     Array<OneD, Array<OneD, Array<OneD, NekDouble>>> MFConnection;
//     ComputeConnection1form(movingframes, MFConnection);

//     Vmath::Vcopy(nq, &du1dei[0], 1, &outarray[0], 1);
//     Vmath::Vcopy(nq, &du2dei[0], 1, &outarray[nq], 1);

//     Array<OneD, NekDouble> Gamma(nq);
//     Vmath::Vmul(nq, &MFConnection[0][1][0], 1, &u2[0], 1, &Gamma[0], 1);

//     // Vmath::Vmul(nq, &MFConnection[0][0][0], 1, &u1[0], 1, &tmp[0], 1);
//     // Vmath::Vvtvp(nq, &MFConnection[0][1][0], 1, &u2[0], 1, &tmp[0], 1,
//     &tmp[0], 1);

//     switch(direction)
//     {
//         case 0:
//         {
//             // \nabla_{e_1} \vec{u}
//             // = ( \nabla u^1 \cdot e_1, \nabla u^2 \cdot e_1 +
//             \omega_{21}(e_1) u^1 + \omega_{21}(e_2) u^2 ) Vmath::Neg(nq,
//             Gamma, 1); Vmath::Vadd(nq, &Gamma[0], 1, &outarray[nq], 1,
//             &outarray[nq], 1);
//         }
//         break;

//         case 1:
//         {
//             // \nabla_{e_e} \vec{u}
//             // = ( \nabla u^1 \cdot e_1 - \omega_{21}(e_1) u^1 -
//             \omega_{21}(e_2) u^2, \nabla u^2 \cdot e_1 ) Vmath::Vadd(nq,
//             &Gamma[0], 1, &outarray[0], 1, &outarray[0], 1);
//         }
//         break;

//         default:
//         break;
//     }

//     return outarray;
// }

Array<OneD, NekDouble> MMFSystem::ComputeCurlSphericalCoord(
    const Array<OneD, const NekDouble> &inarrayphi,
    const Array<OneD, const NekDouble> &inarrayth)
{
    int nq = GetNpoints();

    Array<OneD, NekDouble> x0(nq);
    Array<OneD, NekDouble> x1(nq);
    Array<OneD, NekDouble> x2(nq);

    m_fields[0]->GetCoords(x0, x1, x2);

    Array<OneD, NekDouble> Phi(m_spacedim * nq);
    Array<OneD, NekDouble> Theta(m_spacedim * nq);

    ComputeSphericalTangentVector(Phi, Theta);

    Array<OneD, NekDouble> d1(nq), d2(nq);
    Array<OneD, NekDouble> DirectCurl(nq, 0.0);

    NekDouble sin_varphi, cos_varphi, sin_theta, cos_theta;
    for (int i = 0; i < nq; i++)
    {
        CartesianToNewSpherical(x0[i], x1[i], x2[i], sin_varphi, cos_varphi,
                                sin_theta, cos_theta);

        d1[i] = sin_theta * inarrayphi[i];
        d2[i] = inarrayth[i];
    }

    MMFDirectionalDeriv(Theta, d1, d1);
    MMFDirectionalDeriv(Phi, d2, d2);

    for (int i = 0; i < nq; i++)
    {
        CartesianToNewSpherical(x0[i], x1[i], x2[i], sin_varphi, cos_varphi,
                                sin_theta, cos_theta);
        if (m_MMFActivation[i])
        {
            DirectCurl[i] = (-1.0 / sin_theta) * (d1[i] + d2[i]);
        }
    }

    return DirectCurl;
}

// vth[i]  = physfield[0][i] * uth * sin_theta;
// vphi[i] = physfield[0][i] * uphi;
Array<OneD, NekDouble> MMFSystem::ComputeDivSphericalCoord(
    const Array<OneD, const Array<OneD, NekDouble>> &velocity)
{
    int nq = GetNpoints();

    Array<OneD, NekDouble> vphi(nq, 0.0);
    Array<OneD, NekDouble> vth(nq, 0.0);

    Array<OneD, NekDouble> DirectDiv(nq);

    Array<OneD, NekDouble> x0(nq);
    Array<OneD, NekDouble> x1(nq);
    Array<OneD, NekDouble> x2(nq);

    m_fields[0]->GetCoords(x0, x1, x2);

    Array<OneD, NekDouble> Phi(m_spacedim * nq);
    Array<OneD, NekDouble> Theta(m_spacedim * nq);

    ComputeSphericalTangentVector(Phi, Theta);

    for (int i = 0; i < m_spacedim; ++i)
    {
        Vmath::Vvtvp(nq, &velocity[i][0], 1, &Phi[i * nq], 1, &vphi[0], 1,
                     &vphi[0], 1);
        Vmath::Vvtvp(nq, &velocity[i][0], 1, &Theta[i * nq], 1, &vth[0], 1,
                     &vth[0], 1);
    }

    NekDouble sin_varphi, cos_varphi, sin_theta, cos_theta;
    for (int i = 0; i < nq; ++i)
    {
        CartesianToNewSpherical(x0[i], x1[i], x2[i], sin_varphi, cos_varphi,
                                sin_theta, cos_theta);

        if (m_MMFActivation[i])
        {
            vphi[i] = vphi[i] / sin_theta;
            vth[i]  = vth[i] * sin_theta;
        }
    }

    DirectDiv = ComputeDivSphericalCoord(vphi, vth);

    return DirectDiv;
}

Array<OneD, NekDouble> MMFSystem::ComputeDivSphericalCoord(
    const Array<OneD, const NekDouble> &vphi,
    const Array<OneD, const NekDouble> &vth)
{
    int nq = GetNpoints();

    Array<OneD, NekDouble> dphi(nq, 0.0);
    Array<OneD, NekDouble> dth(nq, 0.0);

    Array<OneD, NekDouble> DirectDiv(nq, 0.0);

    Array<OneD, NekDouble> x0(nq);
    Array<OneD, NekDouble> x1(nq);
    Array<OneD, NekDouble> x2(nq);

    m_fields[0]->GetCoords(x0, x1, x2);

    Array<OneD, NekDouble> Phi(m_spacedim * nq);
    Array<OneD, NekDouble> Theta(m_spacedim * nq);

    ComputeSphericalTangentVector(Phi, Theta);

    MMFDirectionalDeriv(Phi, vphi, dphi);
    MMFDirectionalDeriv(Theta, vth, dth);

    Vmath::Vadd(nq, dphi, 1, dth, 1, DirectDiv, 1);

    NekDouble sin_varphi, cos_varphi, sin_theta, cos_theta;

    for (int i = 0; i < nq; i++)
    {
        CartesianToNewSpherical(x0[i], x1[i], x2[i], sin_varphi, cos_varphi,
                                sin_theta, cos_theta);

        if (m_MMFActivation[i])
        {
            DirectDiv[i] = DirectDiv[i] / sin_theta;
        }
    }

    return DirectDiv;
}

// Compute \vec{C} \cdot ( ( \vec{A} \cdot \nabla ) \vec{B} )
Array<OneD, NekDouble> MMFSystem::ComputeVecCdotNabla(
    const Array<OneD, const NekDouble> &vecC,
    const Array<OneD, const NekDouble> &vecA,
    const Array<OneD, const NekDouble> &vecB)
{
    int nq = m_fields[0]->GetNpoints();

    Array<OneD, NekDouble> outarray(nq, 0.0);

    Array<OneD, NekDouble> tmp(nq);
    Array<OneD, NekDouble> vectmp(nq);
    for (int i = 0; i < m_spacedim; ++i)
    {
        Vmath::Vcopy(nq, &vecB[i * nq], 1, &vectmp[0], 1);

        MMFDirectionalDeriv(vecA, vectmp, tmp);
        Vmath::Vvtvp(nq, &tmp[0], 1, &vecC[i * nq], 1, &outarray[0], 1,
                     &outarray[0], 1);
    }

    return outarray;
}

void MMFSystem::ComputeDivMF(
    DerivType Dtype,
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    Array<OneD, Array<OneD, NekDouble>> &DivMF, const int Verbose)
{
    int nq = m_fields[0]->GetNpoints();

    switch (Dtype)
    {
        case eEuclidean:
        {
            ComputeEuclideanDivMF(movingframes, DivMF);
        }
        break;

        case eCovariant:
        {
            ComputeCovariantDivMF(movingframes, DivMF);
        }
        break;

        default:
            break;
    }

    if (Verbose)
    {
        Array<OneD, Array<OneD, NekDouble>> DivEucMF(m_mfdim);
        Array<OneD, Array<OneD, NekDouble>> DivCovMF(m_mfdim);
        for (int j = 0; j < m_mfdim; ++j)
        {
            DivEucMF[j] = Array<OneD, NekDouble>(nq, 0.0);
            DivCovMF[j] = Array<OneD, NekDouble>(nq, 0.0);
        }

        ComputeEuclideanDivMF(movingframes, DivEucMF);
        ComputeCovariantDivMF(movingframes, DivCovMF);

        for (int i = 0; i < m_mfdim; ++i)
        {
            Vmath::Vsub(nq, DivEucMF[i], 1, DivCovMF[i], 1, DivCovMF[i], 1);
        }

        std::cout << "Div: Euc - Cov = ( " << RootMeanSquare(DivCovMF[0])
                  << " , " << RootMeanSquare(DivCovMF[1]) << " , "
                  << RootMeanSquare(DivCovMF[2]) << " ) " << std::endl;
    }
}

void MMFSystem::ComputeEuclideanDivMF(
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    Array<OneD, Array<OneD, NekDouble>> &DivMF)
{
    int nq = m_fields[0]->GetNpoints();

    DivMF = Array<OneD, Array<OneD, NekDouble>>(m_mfdim);
    for (int j = 0; j < m_mfdim; ++j)
    {
        DivMF[j] = Array<OneD, NekDouble>(nq, 0.0);
    }

    Array<OneD, NekDouble> tmp(nq);
    Array<OneD, NekDouble> Dtmp(nq);

    // case eEuclidean:
    for (int j = 0; j < m_mfdim; ++j)
    {
        for (int k = 0; k < m_spacedim; ++k)
        {
            Vmath::Vcopy(nq, &movingframes[j][k * nq], 1, &tmp[0], 1);
            m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[k], tmp, Dtmp);
            Vmath::Vadd(nq, &Dtmp[0], 1, &DivMF[j][0], 1, &DivMF[j][0], 1);
        }
    }
}

// w^2_1 <e^1> = Connectionform[0][0];
// w^2_1 <e^2> = Connectionform[0][1];
// w^2_1 <e^3> = Connectionform[0][2];

// w^3_1 <e^1> = Connectionform[1][0];
// w^3_1 <e^2> = Connectionform[1][1];
// w^3_1 <e^3> = Connectionform[1][2];

// w^3_2 <e^1> = Connectionform[2][0];
// w^3_2 <e^2> = Connectionform[2][1];
// w^3_2 <e^3> = Connectionform[2][2];
// case eCovariant:
void MMFSystem::ComputeCovariantDivMF(
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    Array<OneD, Array<OneD, NekDouble>> &DivMF)
{
    int nq = m_fields[0]->GetNpoints();

    DivMF = Array<OneD, Array<OneD, NekDouble>>(m_mfdim);
    for (int j = 0; j < m_mfdim; ++j)
    {
        DivMF[j] = Array<OneD, NekDouble>(nq, 0.0);
    }

    Array<OneD, NekDouble> tmp(nq);
    Array<OneD, NekDouble> Dtmp(nq);

    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> MFConnection;
    Compute2DConnection1form(movingframes, MFConnection);

    // \nabla \cdot e^1 = \Gamma^2_{12} = w^2_1 <e^2>
    // Vmath::Vcopy(nq, &MFConnection[0][1][0], 1, &DivMF[0][0], 1);

    // \nabla \cdot e^2 = \Gamma^1_{21} = -w^2_1 <e^1>
    // Vmath::Vcopy(nq, &MFConnection[0][0][0], 1, &DivMF[1][0], 1);
    // Vmath::Neg(nq, &DivMF[1][0], 1);

    // \nabla \cdot e^1 = \Gamma^2_{12} = w^2_1 <e^2>
    Vmath::Vcopy(nq, &MFConnection[0][1][0], 1, &DivMF[0][0], 1);
    Vmath::Neg(nq, &DivMF[0][0], 1);

    // \nabla \cdot e^2 = \Gamma^1_{21} = -w^2_1 <e^1>
    Vmath::Vcopy(nq, &MFConnection[0][0][0], 1, &DivMF[1][0], 1);

    // \nabla \cdot e^3 = \Gamma^1_{31} + \Gamma^2_{32} = - (
    // w^3_1 <e^1> + w^3_2 <e^2> )
    Vmath::Vadd(nq, &MFConnection[1][0][0], 1, &MFConnection[2][1][0], 1,
                &DivMF[2][0], 1);
    Vmath::Neg(nq, &DivMF[2][0], 1);
}

// void MMFSystem::ComputeDivMF(
//     DerivType Dtype,
//     const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
//     Array<OneD, Array<OneD, NekDouble>> &DivMF, const int Verbose)
// {
//     int nq = m_fields[0]->GetNpoints();

//     Array<OneD, Array<OneD, NekDouble>> DivMFEuc(m_mfdim);
//     Array<OneD, Array<OneD, NekDouble>> DivMFCov(m_mfdim);

//     DivMF = Array<OneD, Array<OneD, NekDouble>>(m_mfdim);
//     DivMFEuc = Array<OneD, Array<OneD, NekDouble>>(m_mfdim);
//     DivMFCov = Array<OneD, Array<OneD, NekDouble>>(m_mfdim);
//     for (int j = 0; j < m_mfdim; ++j)
//     {
//         DivMF[j] = Array<OneD, NekDouble>(nq, 0.0);
//         DivMFEuc[j] = Array<OneD, NekDouble>(nq, 0.0);
//         DivMFCov[j] = Array<OneD, NekDouble>(nq, 0.0);
//     }

//     Array<OneD, NekDouble> tmp(nq);
//     Array<OneD, NekDouble> Dtmp(nq);

//     // case eEuclidean:
//     for (int j = 0; j < m_mfdim; ++j)
//     {
//         for (int k = 0; k < m_spacedim; ++k)
//         {
//             Vmath::Vcopy(nq, &movingframes[j][k * nq], 1, &tmp[0], 1);
//             m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[k], tmp,
//             Dtmp); Vmath::Vadd(nq, &Dtmp[0], 1, &DivMFEuc[j][0], 1,
//             &DivMFEuc[j][0], 1);
//         }
//         std::cout << "j = " << j << ", DivMFEuc = " <<
//         RootMeanSquare(DivMFEuc[j]) << std::endl;
//     }

//     // w^2_1 <e^1> = Connectionform[0][0];
//     // w^2_1 <e^2> = Connectionform[0][1];
//     // w^2_1 <e^3> = Connectionform[0][2];

//     // w^3_1 <e^1> = Connectionform[1][0];
//     // w^3_1 <e^2> = Connectionform[1][1];
//     // w^3_1 <e^3> = Connectionform[1][2];

//     // w^3_2 <e^1> = Connectionform[2][0];
//     // w^3_2 <e^2> = Connectionform[2][1];
//     // w^3_2 <e^3> = Connectionform[2][2];
//     // case eCovariant:

//     Array<OneD, Array<OneD, Array<OneD, NekDouble>>> MFConnection;
//     ComputeConnection1form(movingframes, MFConnection);

//     // \nabla \cdot e^1 = \Gamma^2_{12} = w^2_1 <e^2>
//     // Vmath::Vcopy(nq, &MFConnection[0][1][0], 1, &DivMF[0][0], 1);

//     // \nabla \cdot e^2 = \Gamma^1_{21} = -w^2_1 <e^1>
//     // Vmath::Vcopy(nq, &MFConnection[0][0][0], 1, &DivMF[1][0], 1);
//     // Vmath::Neg(nq, &DivMF[1][0], 1);

//     // \nabla \cdot e^1 = \Gamma^2_{12} = w^2_1 <e^2>
//     Vmath::Vcopy(nq, &MFConnection[0][1][0], 1, &DivMFCov[0][0], 1);

//     // \nabla \cdot e^2 = \Gamma^1_{21} = -w^2_1 <e^1>
//     Vmath::Vcopy(nq, &MFConnection[0][0][0], 1, &DivMFCov[1][0], 1);
//     Vmath::Neg(nq, &DivMFCov[1][0], 1);

//     // \nabla \cdot e^3 = \Gamma^1_{31} + \Gamma^2_{32} = - (
//     // w^3_1 <e^1> + w^3_2 <e^2> )
//     Vmath::Vadd(nq, &MFConnection[1][0][0], 1, &MFConnection[2][1][0],
//                 1, &DivMFCov[2][0], 1);
//     Vmath::Neg(nq, &DivMFCov[2][0], 1);

//     // Assign it according to Euclidean or Covariant
//     if(Dtype==eEuclidean)
//     {
//         for (int i=0; i<m_mfdim; ++i)
//         {
//             Vmath::Vcopy(nq, DivMFEuc[i], 1, DivMF[i], 1);
//         }
//     }

//     else if (Dtype==eCovariant)
//     {
//         for (int i=0; i<m_mfdim; ++i)
//         {
//              Vmath::Vcopy(nq, DivMFCov[i], 1, DivMF[i], 1);
//         }
//     }

//     if(Verbose)
//     {
//         std::cout << "Dtype = " << DerivTypeMap[Dtype] << ", DivMF = ( "
//                 << RootMeanSquare(m_DivMF[0], m_MMFActivation) << " , "
//                 << RootMeanSquare(m_DivMF[1], m_MMFActivation) << " , "
//                 << RootMeanSquare(m_DivMF[2], m_MMFActivation) << " ) "
//                 << std::endl ;

//         for (int i=0; i<m_mfdim; ++i)
//         {
//             Vmath::Vsub(nq, DivMFEuc[i], 1, DivMFCov[i], 1,  DivMFCov[i], 1);
//         }
//         std::cout << "Div: Euc - Cov = ( " << RootMeanSquare(DivMFCov[0]) <<
//         " , "
//         << RootMeanSquare(DivMFCov[1]) << " , " <<
//         RootMeanSquare(DivMFCov[2]) << " ) " << std::endl;
//     }
// }

void MMFSystem::ComputeCurlMF(
    DerivType Dtype,
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    Array<OneD, Array<OneD, NekDouble>> &CurlMF, const int Verbose)
{
    int nCurl = 4;
    int nq    = m_fields[0]->GetNpoints();

    switch (Dtype)
    {
        case eEuclidean:
        {
            ComputeEuclideanCurlMF(movingframes, CurlMF);
        }
        break;

        case eCovariant:
        {
            ComputeCovariantCurlMF(movingframes, CurlMF);
        }
        break;

        default:
            break;
    }

    if (Verbose)
    {
        Array<OneD, Array<OneD, NekDouble>> CurlEucMF(nCurl);
        Array<OneD, Array<OneD, NekDouble>> CurlCovMF(nCurl);
        for (int j = 0; j < nCurl; ++j)
        {
            CurlEucMF[j] = Array<OneD, NekDouble>(nq, 0.0);
            CurlCovMF[j] = Array<OneD, NekDouble>(nq, 0.0);
        }

        ComputeEuclideanCurlMF(movingframes, CurlEucMF);
        ComputeCovariantCurlMF(movingframes, CurlCovMF);

        for (int i = 0; i < nCurl; ++i)
        {
            Vmath::Vsub(nq, CurlEucMF[i], 1, CurlCovMF[i], 1, CurlCovMF[i], 1);
        }

        std::cout << "Curl: Euc - Cov = ( " << RootMeanSquare(CurlCovMF[0])
                  << " , " << RootMeanSquare(CurlCovMF[1]) << " , "
                  << RootMeanSquare(CurlCovMF[2]) << " , "
                  << RootMeanSquare(CurlCovMF[3]) << " ) " << std::endl;
    }
}

// CurlMF[0] = e^3 \cdot (\nabla \times e^1)
// CurlMF[1] = e^3 \cdot (\nabla \times e^2)
// CurlMF[2] = e^1 \cdot (\nabla \times e^3)
// CurlMF[3] = e^2 \cdot (\nabla \times e^3)
void MMFSystem::ComputeCovariantCurlMF(
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    Array<OneD, Array<OneD, NekDouble>> &CurlMF)
{
    int nq    = m_fields[0]->GetNpoints();
    int nCurl = 4;

    CurlMF = Array<OneD, Array<OneD, NekDouble>>(nCurl);
    for (int j = 0; j < nCurl; ++j)
    {
        CurlMF[j] = Array<OneD, NekDouble>(nq, 0.0);
    }

    Array<OneD, NekDouble> tmp(nq);
    Array<OneD, NekDouble> Dtmp(nq);

    // w^2_1 <e^1> = Connectionform[0][0];
    // w^2_1 <e^2> = Connectionform[0][1];
    // w^2_1 <e^3> = Connectionform[0][2];

    // w^3_1 <e^1> = Connectionform[1][0];
    // w^3_1 <e^2> = Connectionform[1][1];
    // w^3_1 <e^3> = Connectionform[1][2];

    // w^3_2 <e^1> = Connectionform[2][0];
    // w^3_2 <e^2> = Connectionform[2][1];
    // w^3_2 <e^3> = Connectionform[2][2];

    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> MFConnection;
    Compute2DConnection1form(movingframes, MFConnection);

    // int index, MFindx, ProjMFindx;
    // CurlMF[0] = e^3 \cdot (\nabla \times e^1)
    // \nabla \times e^1 = - \Gamma^2_{11} = - w^2_1 <e^1>

    // index      = 0;
    // MFindx     = 0;
    // ProjMFindx = 2;

    Vmath::Vcopy(nq, &MFConnection[0][0][0], 1, &CurlMF[0][0], 1);
    // Vmath::Neg(nq, &CurlMF[0][0], 1);

    // CurlMF[1] = e^3 \cdot (\nabla \times e^2)
    // \nabla \times e^2 = \Gamma^2_{12} = w^2_1 <e^2>
    // index      = index + 1;
    // MFindx     = 1;
    // ProjMFindx = 2;

    Vmath::Vcopy(nq, &MFConnection[0][1][0], 1, &CurlMF[1][0], 1);

    // \nabla \cdot e^3 = \Gamma^1_{31} + \Gamma^2_{32} = - (
    // w^3_1 <e^1> + w^3_2 <e^2> )
    // Vmath::Vadd(nq, &MFConnection[1][0][0], 1,
    // &MFConnection[2][1][0],
    //             1, &DivMF[2][0], 1);
    // index      = index + 1;
    // MFindx     = 2;
    // ProjMFindx = 0;

    // CurlMF[2] = e^1 \cdot (\nabla \times e^3) = w^3_2 <e^3>
    // w^3_2 <e^3> = Connectionform[2][2];
    Vmath::Vcopy(nq, &MFConnection[2][2][0], 1, &CurlMF[2][0], 1);

    // index      = index + 1;
    // MFindx     = 2;
    // ProjMFindx = 1;

    // CurlMF[2] = e^2 \cdot (\nabla \times e^3) = - w^3_1 <e^3>
    // w^3_1 <e^3> = Connectionform[1][2];
    Vmath::Vcopy(nq, &MFConnection[1][2][0], 1, &CurlMF[3][0], 1);
    Vmath::Neg(nq, &CurlMF[3][0], 1);
}

// Compute $\nabla \times \mathbf{e}^i$
// CurlMF[0] = e^3 \cdot (\nabla \times e^1)
// CurlMF[1] = e^3 \cdot (\nabla \times e^2)
// CurlMF[2] = e^1 \cdot (\nabla \times e^3)
// CurlMF[3] = e^2 \cdot (\nabla \times e^3)
// Output: CurlMF[i][j][] = $(\nabla \times \mathbf{e}^i) \cdot
// \mathbf{e}^j$
void MMFSystem::ComputeEuclideanCurlMF(
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    Array<OneD, Array<OneD, NekDouble>> &CurlMF)
{
    int nq    = m_fields[0]->GetNpoints();
    int nCurl = 4;

    // Compute Curl of MF: CurlMF[i][j] = (\nabla \times e^i) cdot e^j
    CurlMF = Array<OneD, Array<OneD, NekDouble>>(nCurl);
    for (int i = 0; i < nCurl; ++i)
    {
        CurlMF[i] = Array<OneD, NekDouble>(nq, 0.0);
    }

    Array<OneD, Array<OneD, NekDouble>> CurlMFtmp(m_spacedim);
    for (int i = 0; i < m_spacedim; ++i)
    {
        CurlMFtmp[i] = Array<OneD, NekDouble>(nq, 0.0);
    }

    int index, MFindx, ProjMFindx;
    // CurlMF[0] = e^3 \cdot (\nabla \times e^1)
    index      = 0;
    MFindx     = 0;
    ProjMFindx = 2;
    ComputeCurl(movingframes[MFindx], CurlMFtmp);
    for (int i = 0; i < m_spacedim; ++i)
    {
        Vmath::Vvtvp(nq, &movingframes[ProjMFindx][i * nq], 1, &CurlMFtmp[i][0],
                     1, &CurlMF[index][0], 1, &CurlMF[index][0], 1);
    }

    // CurlMF[1] = e^3 \cdot (\nabla \times e^2)
    index      = index + 1;
    MFindx     = 1;
    ProjMFindx = 2;
    ComputeCurl(movingframes[MFindx], CurlMFtmp);
    for (int i = 0; i < m_spacedim; ++i)
    {
        Vmath::Vvtvp(nq, &movingframes[ProjMFindx][i * nq], 1, &CurlMFtmp[i][0],
                     1, &CurlMF[index][0], 1, &CurlMF[index][0], 1);
    }

    // CurlMF[2] = e^1 \cdot (\nabla \times e^3)
    index      = index + 1;
    MFindx     = 2;
    ProjMFindx = 0;
    ComputeCurl(movingframes[MFindx], CurlMFtmp);
    for (int i = 0; i < m_spacedim; ++i)
    {
        Vmath::Vvtvp(nq, &movingframes[ProjMFindx][i * nq], 1, &CurlMFtmp[i][0],
                     1, &CurlMF[index][0], 1, &CurlMF[index][0], 1);
    }

    // CurlMF[3] = e^2 \cdot (\nabla \times e^3)
    index      = index + 1;
    MFindx     = 2;
    ProjMFindx = 1;
    ComputeCurl(movingframes[MFindx], CurlMFtmp);
    for (int i = 0; i < m_spacedim; ++i)
    {
        Vmath::Vvtvp(nq, &movingframes[ProjMFindx][i * nq], 1, &CurlMFtmp[i][0],
                     1, &CurlMF[index][0], 1, &CurlMF[index][0], 1);
    }
}

void MMFSystem::ComputeExactDivMF(
    const Array<OneD, const int> &Activation,
    Array<OneD, Array<OneD, NekDouble>> &ExactDivMF)
{
    int nq = m_fields[0]->GetNpoints();

    Array<OneD, NekDouble> x0(nq);
    Array<OneD, NekDouble> x1(nq);
    Array<OneD, NekDouble> x2(nq);

    m_fields[0]->GetCoords(x0, x1, x2);

    NekDouble xp, yp, zp, rad;

    ExactDivMF = Array<OneD, Array<OneD, NekDouble>>(m_mfdim);
    for (int j = 0; j < m_mfdim; ++j)
    {
        ExactDivMF[j] = Array<OneD, NekDouble>(nq, 0.0);
    }

    if (m_surfaceType == SolverUtils::ePlane)
    {
        for (int i = 0; i < nq; ++i)
        {
            xp = x0[i];
            yp = x1[i];
            zp = x2[i];

            rad = sqrt(xp * xp + yp * yp);

            if (Activation[i] > 0)
            {
                ExactDivMF[0][i] = 1.0 / rad;
            }
        }
    }

    else if (m_surfaceType == SolverUtils::eSphere)
    {
        for (int i = 0; i < nq; ++i)
        {
            xp = x0[i];
            yp = x1[i];
            zp = x2[i];

            rad = sqrt(xp * xp + yp * yp + zp * zp);

            NekDouble sin_varphi, cos_varphi, sin_theta, cos_theta;
            CartesianToNewSpherical(xp, yp, zp, sin_varphi, cos_varphi,
                                    sin_theta, cos_theta);

            if (Activation[i] > 0)
            {
                ExactDivMF[0][i] = (1.0 / rad) * (cos_theta / sin_theta);
                ExactDivMF[1][i] = 0.0;
                ExactDivMF[2][i] = (2.0 / rad);
            }
        }
    }
}

void MMFSystem::ComputeExactCurlMF(
    const Array<OneD, const int> &Activation,
    Array<OneD, Array<OneD, NekDouble>> &ExactCurlMF)
{
    int nq = m_fields[0]->GetNpoints();

    Array<OneD, NekDouble> x0(nq);
    Array<OneD, NekDouble> x1(nq);
    Array<OneD, NekDouble> x2(nq);

    m_fields[0]->GetCoords(x0, x1, x2);

    NekDouble xp, yp, zp, rad;

    ExactCurlMF = Array<OneD, Array<OneD, NekDouble>>(m_mfdim);
    for (int j = 0; j < m_mfdim; ++j)
    {
        ExactCurlMF[j] = Array<OneD, NekDouble>(nq, 0.0);
    }

    if (m_surfaceType == SolverUtils::ePlane)
    {
        for (int i = 0; i < nq; ++i)
        {
            xp = x0[i];
            yp = x1[i];
            zp = x2[i];

            rad = sqrt(xp * xp + yp * yp);

            if (Activation[i] > 0)
            {
                ExactCurlMF[0][i] = 1.0 / rad;
            }
        }
    }

    else if (m_surfaceType == SolverUtils::eSphere)
    {
        for (int i = 0; i < nq; ++i)
        {
            xp = x0[i];
            yp = x1[i];
            zp = x2[i];

            rad = sqrt(xp * xp + yp * yp + zp * zp);

            NekDouble sin_varphi, cos_varphi, sin_theta, cos_theta;
            CartesianToNewSpherical(xp, yp, zp, sin_varphi, cos_varphi,
                                    sin_theta, cos_theta);

            if (Activation[i] > 0)
            {
                ExactCurlMF[0][i] = 0.0;
                ExactCurlMF[1][i] = (1.0 / rad) * (cos_theta / sin_theta);
            }
        }
    }
}

void MMFSystem::TestDivMF(
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    const Array<OneD, const int> &Activation)
{
    int nq = m_fields[0]->GetNpoints();

    Array<OneD, Array<OneD, NekDouble>> DivMF(m_mfdim);
    Array<OneD, Array<OneD, NekDouble>> DivMFCov(m_mfdim);
    for (int j = 0; j < m_mfdim; ++j)
    {
        DivMF[j]    = Array<OneD, NekDouble>(nq, 0.0);
        DivMFCov[j] = Array<OneD, NekDouble>(nq, 0.0);
    }

    // Compute Divergence in two different ways
    ComputeDivMF(eEuclidean, movingframes, DivMF, 1);
    ComputeDivMF(eCovariant, movingframes, DivMFCov, 1);

    Array<OneD, Array<OneD, NekDouble>> ExactDivMF(m_mfdim);

    Array<OneD, NekDouble> ErrDivEuc1(nq);
    Array<OneD, NekDouble> ErrDivCov1(nq);
    Array<OneD, NekDouble> ErrDivDiff1(nq);

    Array<OneD, NekDouble> ErrDivEuc2(nq);
    Array<OneD, NekDouble> ErrDivCov2(nq);
    Array<OneD, NekDouble> ErrDivDiff2(nq);

    Array<OneD, NekDouble> ErrDivEuc3(nq);
    Array<OneD, NekDouble> ErrDivCov3(nq);
    Array<OneD, NekDouble> ErrDivDiff3(nq);

    ComputeExactDivMF(Activation, ExactDivMF);

    std::cout << "Div Exact: e1 = " << RootMeanSquare(ExactDivMF[0])
              << ", e2 = " << RootMeanSquare(ExactDivMF[1])
              << ", e3 = " << RootMeanSquare(ExactDivMF[2]) << std::endl;

    if (m_surfaceType == SolverUtils::ePlane)
    {

        Vmath::Vsub(nq, &DivMF[0][0], 1, &ExactDivMF[0][0], 1, &ErrDivEuc1[0],
                    1);
        Vmath::Vsub(nq, &DivMFCov[0][0], 1, &ExactDivMF[0][0], 1,
                    &ErrDivCov1[0], 1);
        Vmath::Vsub(nq, &DivMF[0][0], 1, &DivMFCov[0][0], 1, &ErrDivDiff1[0],
                    1);

        std::cout << "Div Err: e1: Euc = " << RootMeanSquare(ErrDivEuc1)
                  << ", Cov = " << RootMeanSquare(ErrDivCov1)
                  << ", Diff = " << RootMeanSquare(ErrDivDiff1) << std::endl;
    }

    else if (m_surfaceType == SolverUtils::eSphere)
    {
        Array<OneD, NekDouble> ErrDivEuc(3, 0.0);
        Array<OneD, NekDouble> ErrDivCov(3, 0.0);
        Array<OneD, NekDouble> ErrDivDiff(3, 0.0);

        NekDouble tmp;
        int cnt = 0;
        for (int i = 0; i < nq; ++i)
        {
            if (m_MMFActivation[i])
            {
                for (int k = 0; k < 3; ++k)
                {
                    tmp = DivMF[k][i] - ExactDivMF[k][i];
                    ErrDivEuc[k] += tmp * tmp;

                    tmp = DivMFCov[k][i] - ExactDivMF[k][i];
                    ErrDivCov[k] += tmp * tmp;

                    tmp = DivMF[k][i] - DivMFCov[k][i];
                    ErrDivDiff[k] += tmp * tmp;
                }
                cnt++;
            }
        }

        for (int k = 0; k < 3; ++k)
        {
            ErrDivEuc[k]  = sqrt(ErrDivEuc[k] / cnt);
            ErrDivCov[k]  = sqrt(ErrDivCov[k] / cnt);
            ErrDivDiff[k] = sqrt(ErrDivDiff[k] / cnt);

            std::cout << "Div Err: e" << k << " : Euc = " << ErrDivEuc[k]
                      << ", Cov = " << ErrDivCov[k]
                      << ", Diff = " << ErrDivDiff[k] << std::endl;
        }
    }
} // namespace SolverUtils

void MMFSystem::TestCurlMF(
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    const Array<OneD, const int> &Activation)
{
    int nq = m_fields[0]->GetNpoints();

    Array<OneD, Array<OneD, NekDouble>> CurlMF(m_mfdim);
    Array<OneD, Array<OneD, NekDouble>> CurlMFCov(m_mfdim);
    for (int j = 0; j < m_mfdim; ++j)
    {
        CurlMF[j]    = Array<OneD, NekDouble>(nq, 0.0);
        CurlMFCov[j] = Array<OneD, NekDouble>(nq, 0.0);
    }

    // Compute Divergence in two different ways
    // CurlMF[0] = e^3 \cdot (\nabla \times e^1)
    // CurlMF[1] = e^3 \cdot (\nabla \times e^2)
    // CurlMF[2] = e^1 \cdot (\nabla \times e^3)
    // CurlMF[3] = e^2 \cdot (\nabla \times e^3)s
    ComputeCurlMF(eEuclidean, movingframes, CurlMF);
    ComputeCurlMF(eCovariant, movingframes, CurlMFCov);

    Array<OneD, Array<OneD, NekDouble>> ExactCurlMF(m_mfdim);

    ComputeExactCurlMF(Activation, ExactCurlMF);
    // if (m_surfaceType == SolverUtils::ePlane)
    // {
    //     Vmath::Vsub(nq, &DivMF[0][0], 1, &ExactDivMF[0][0], 1,
    //     &ErrDivEuc1[0],
    //                 1);
    //     Vmath::Vsub(nq, &DivMFCov[0][0], 1, &ExactDivMF[0][0], 1,
    //                 &ErrDivCov1[0], 1);
    //     Vmath::Vsub(nq, &DivMF[0][0], 1, &DivMFCov[0][0], 1,
    //     &ErrDivDiff1[0],
    //                 1);

    //     std::cout << "Div Err: e1: Euc = " <<
    //     RootMeanSquare(ErrDivEuc1)
    //              << ", Cov = " << RootMeanSquare(ErrDivCov1)
    //              << ", Diff = " << RootMeanSquare(ErrDivDiff1) <<
    //              std::endl;
    // }

    if (m_surfaceType == SolverUtils::eSphere)
    {
        Array<OneD, NekDouble> ErrCurlEuc(2, 0.0);
        Array<OneD, NekDouble> ErrCurlCov(2, 0.0);
        Array<OneD, NekDouble> ErrCurlDiff(2, 0.0);

        NekDouble tmp;
        int cnt = 0;
        for (int i = 0; i < nq; ++i)
        {
            if (m_MMFActivation[i])
            {
                for (int k = 0; k < 2; ++k)
                {
                    tmp = CurlMF[k][i] - ExactCurlMF[k][i];
                    ErrCurlEuc[k] += tmp * tmp;

                    tmp = CurlMFCov[k][i] - ExactCurlMF[k][i];
                    ErrCurlCov[k] += tmp * tmp;

                    tmp = CurlMFCov[k][i] - CurlMF[k][i];
                    ErrCurlDiff[k] += tmp * tmp;
                }
                cnt++;
            }
        }

        for (int k = 0; k < 2; ++k)
        {
            ErrCurlEuc[k]  = sqrt(ErrCurlEuc[k] / cnt);
            ErrCurlCov[k]  = sqrt(ErrCurlCov[k] / cnt);
            ErrCurlDiff[k] = sqrt(ErrCurlDiff[k] / cnt);

            std::cout << "Curl Err: e" << k << " : Euc = " << ErrCurlEuc[k]
                      << ", Cov = " << ErrCurlCov[k]
                      << ", Diff = " << ErrCurlDiff[k] << std::endl;
        }
    }
} // namespace SolverUtils

void MMFSystem::ComputeMFtrace(
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &MFtraceFwd,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &MFtraceBwd)
{
    int nq              = GetTotPoints();
    int nTraceNumPoints = GetTraceTotPoints();

    MFtraceFwd = Array<OneD, Array<OneD, Array<OneD, NekDouble>>>(m_mfdim);
    MFtraceBwd = Array<OneD, Array<OneD, Array<OneD, NekDouble>>>(m_mfdim);

    for (int j = 0; j < m_mfdim; ++j)
    {
        MFtraceFwd[j] = Array<OneD, Array<OneD, NekDouble>>(m_spacedim);
        MFtraceBwd[j] = Array<OneD, Array<OneD, NekDouble>>(m_spacedim);
        for (int i = 0; i < m_spacedim; ++i)
        {
            MFtraceFwd[j][i] = Array<OneD, NekDouble>(nTraceNumPoints, 0.0);
            MFtraceBwd[j][i] = Array<OneD, NekDouble>(nTraceNumPoints, 0.0);
        }
    }

    // m_MFtraceFwd[0] = e^1_{Fwd}, m_MFtraceFwd[1] = e^2_{Fwd}
    Array<OneD, NekDouble> tmp(nq);
    Array<OneD, NekDouble> Fwdtmp(nTraceNumPoints, 0.0);
    Array<OneD, NekDouble> Bwdtmp(nTraceNumPoints, 0.0);
    for (int j = 0; j < m_mfdim; ++j)
    {
        for (int i = 0; i < m_spacedim; ++i)
        {
            Vmath::Vcopy(nq, &movingframes[j][i * nq], 1, &tmp[0], 1);

            GetFwdBwdMFTrace(tmp, Fwdtmp, Bwdtmp);

            Vmath::Vcopy(nTraceNumPoints, &Fwdtmp[0], 1, &MFtraceFwd[j][i][0],
                         1);
            Vmath::Vcopy(nTraceNumPoints, &Bwdtmp[0], 1, &MFtraceBwd[j][i][0],
                         1);
        }
    }
}

//    VectorCrossProd(MF3tmp, MF1tmp, MFtmpCurl[0]);
//    VectorCrossProd(MF2tmp, MF3tmp, MFtmpCurl[1]);
//    VectorCrossProd(MF1tmp, MF2tmp, MFtmpCurl[2]);
void MMFSystem::ComputeMFtimesMF(
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    Array<OneD, Array<OneD, NekDouble>> &outarray)
{
    int nq = GetTotPoints();

    outarray = Array<OneD, Array<OneD, NekDouble>>(m_mfdim);

    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> MFtmpCurl(m_mfdim);
    for (int j = 0; j < m_mfdim; ++j)
    {
        outarray[j]  = Array<OneD, NekDouble>(nq * m_spacedim);
        MFtmpCurl[j] = Array<OneD, Array<OneD, NekDouble>>(m_spacedim);
        for (int k = 0; k < m_spacedim; ++k)
        {
            MFtmpCurl[j][k] = Array<OneD, NekDouble>(nq);
        }
    }

    Array<OneD, Array<OneD, NekDouble>> MF1tmp(m_spacedim);
    Array<OneD, Array<OneD, NekDouble>> MF2tmp(m_spacedim);
    Array<OneD, Array<OneD, NekDouble>> MF3tmp(m_spacedim);
    for (int k = 0; k < m_spacedim; ++k)
    {
        MF1tmp[k] = Array<OneD, NekDouble>(nq);
        MF2tmp[k] = Array<OneD, NekDouble>(nq);
        MF3tmp[k] = Array<OneD, NekDouble>(nq);

        Vmath::Vcopy(nq, &movingframes[0][k * nq], 1, &MF1tmp[k][0], 1);
        Vmath::Vcopy(nq, &movingframes[1][k * nq], 1, &MF2tmp[k][0], 1);
        Vmath::Vcopy(nq, &movingframes[2][k * nq], 1, &MF3tmp[k][0], 1);
    }

    VectorCrossProd(MF3tmp, MF1tmp, MFtmpCurl[0]);
    VectorCrossProd(MF2tmp, MF3tmp, MFtmpCurl[1]);
    VectorCrossProd(MF1tmp, MF2tmp, MFtmpCurl[2]);

    for (int j = 0; j < m_mfdim; ++j)
    {
        for (int k = 0; k < m_spacedim; ++k)
        {
            Vmath::Vcopy(nq, &MFtmpCurl[j][k][0], 1, &outarray[j][k * nq], 1);
        }
    }
}

void MMFSystem::VectorDotProd(
    const Array<OneD, const Array<OneD, NekDouble>> &v1,
    const Array<OneD, const Array<OneD, NekDouble>> &v2,
    Array<OneD, NekDouble> &v3)
{
    int coordim = v1.size();
    int nq      = v1[0].size();

    v3 = Array<OneD, NekDouble>(nq, 0.0);
    for (int i = 0; i < coordim; ++i)
    {
        Vmath::Vvtvp(nq, &v1[i][0], 1, &v2[i][0], 1, &v3[0], 1, &v3[0], 1);
    }
}

NekDouble MMFSystem::VectorDotProd(const Array<OneD, const NekDouble> &v1,
                                   const Array<OneD, const NekDouble> &v2)
{
    NekDouble dotprod;

    dotprod = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];

    return dotprod;
}

Array<OneD, NekDouble> MMFSystem::VectorDiff(
    const Array<OneD, const Array<OneD, NekDouble>> &velocity1,
    const Array<OneD, const Array<OneD, NekDouble>> &velocity2)
{
    int nq = velocity1[0].size();

    Array<OneD, NekDouble> tmp(nq, 0.0);
    Array<OneD, NekDouble> diff(nq, 0.0);
    for (int i = 0; i < m_spacedim; ++i)
    {
        Vmath::Vsub(nq, &velocity1[i][0], 1, &velocity2[i][0], 1, &tmp[0], 1);
        Vmath::Vvtvp(nq, &tmp[0], 1, &tmp[0], 1, &diff[0], 1, &diff[0], 1);
    }

    Vmath::Vsqrt(nq, diff, 1, diff, 1);

    return diff;
}

void MMFSystem::MFDotProd(
    const Array<OneD, const Array<OneD, NekDouble>> &MF1,
    const Array<OneD, const Array<OneD, NekDouble>> &MF2,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &outarray)
{
    int nq = GetNpoints();

    outarray = Array<OneD, Array<OneD, Array<OneD, NekDouble>>>(m_mfdim);
    for (int i = 0; i < m_mfdim; ++i)
    {
        outarray[i] = Array<OneD, Array<OneD, NekDouble>>(m_mfdim);
        for (int j = 0; j < m_mfdim; ++j)
        {
            outarray[i][j] = Array<OneD, NekDouble>(nq, 0.0);
            for (int k = 0; k < m_spacedim; ++k)
            {
                Vmath::Vvtvp(nq, &MF1[i][k * nq], 1, &MF2[j][k * nq], 1,
                             &outarray[i][j][0], 1, &outarray[i][j][0], 1);
            }
        }
    }
}

void MMFSystem::MFDotProduct(
    const Array<OneD, const Array<OneD, NekDouble>> &MF1,
    const Array<OneD, const Array<OneD, NekDouble>> &MF2,
    Array<OneD, Array<OneD, NekDouble>> &outarray)
{
    int nq = GetNpoints();

    outarray = Array<OneD, Array<OneD, NekDouble>>(m_spacedim);
    for (int i = 0; i < m_mfdim; ++i)
    {
        outarray[i] = Array<OneD, NekDouble>(nq, 0.0);
        for (int k = 0; k < m_spacedim; ++k)
        {
            Vmath::Vvtvp(nq, &MF1[i][k * nq], 1, &MF2[i][k * nq], 1,
                         &outarray[i][0], 1, &outarray[i][0], 1);
        }
    }
}

// Vector Cross Product with one dimensional vector
/**
 * Computes the vector cross-product in 3D of \a v1 and \a v2, storing
 * the result in \a v3.
 * @param   v1          First input vector.
 * @param   v2          Second input vector.
 * @param   v3          Output vector computed to be orthogonal to
 *                      both \a v1 and \a v2.
 */
Array<OneD, NekDouble> MMFSystem::VectorCrossProdMF(
    const Array<OneD, const NekDouble> &v1,
    const Array<OneD, const NekDouble> &v2)
{
    int nq = GetNpoints();

    Array<OneD, NekDouble> v3(m_spacedim * nq);
    Array<OneD, NekDouble> temp(nq);

    Vmath::Vmul(nq, &v1[2 * nq], 1, &v2[1 * nq], 1, &temp[0], 1);
    Vmath::Vvtvm(nq, &v1[1 * nq], 1, &v2[2 * nq], 1, &temp[0], 1, &v3[0], 1);

    Vmath::Vmul(nq, &v1[0], 1, &v2[2 * nq], 1, &temp[0], 1);
    Vmath::Vvtvm(nq, &v1[2 * nq], 1, &v2[0], 1, &temp[0], 1, &v3[1 * nq], 1);

    Vmath::Vmul(nq, &v1[1 * nq], 1, &v2[0], 1, &temp[0], 1);
    Vmath::Vvtvm(nq, &v1[0], 1, &v2[1 * nq], 1, &temp[0], 1, &v3[2 * nq], 1);

    return v3;
}

// Vector Cross Product with two dimensional vector
void MMFSystem::VectorCrossProd(
    const Array<OneD, const Array<OneD, NekDouble>> &v1,
    const Array<OneD, const Array<OneD, NekDouble>> &v2,
    Array<OneD, Array<OneD, NekDouble>> &v3)
{
    ASSERTL0(v1.size() == 3, "Input 1 has dimension not equal to 3.");
    ASSERTL0(v2.size() == 3, "Input 2 has dimension not equal to 3.");
    ASSERTL0(v3.size() == 3, "Output vector has dimension not equal to 3.");

    int nq = v1[0].size();
    Array<OneD, NekDouble> temp(nq);

    Vmath::Vmul(nq, v1[2], 1, v2[1], 1, temp, 1);
    Vmath::Vvtvm(nq, v1[1], 1, v2[2], 1, temp, 1, v3[0], 1);

    Vmath::Vmul(nq, v1[0], 1, v2[2], 1, temp, 1);
    Vmath::Vvtvm(nq, v1[2], 1, v2[0], 1, temp, 1, v3[1], 1);

    Vmath::Vmul(nq, v1[1], 1, v2[0], 1, temp, 1);
    Vmath::Vvtvm(nq, v1[0], 1, v2[1], 1, temp, 1, v3[2], 1);
}

// Pointwise Vector Cross Product
void MMFSystem::VectorCrossProd(const Array<OneD, const NekDouble> &v1,
                                const Array<OneD, const NekDouble> &v2,
                                Array<OneD, NekDouble> &outarray)
{
    ASSERTL0(v1.size() == 3, "Input 1 has dimension not equal to 3.");
    ASSERTL0(v2.size() == 3, "Input 2 has dimension not equal to 3.");

    outarray = Array<OneD, NekDouble>(3);

    outarray[0] = v1[1] * v2[2] - v1[2] * v2[1];
    outarray[1] = v1[2] * v2[0] - v1[0] * v2[2];
    outarray[2] = v1[0] * v2[1] - v1[1] * v2[0];
}

void MMFSystem::ComputeCurl(const Array<OneD, const NekDouble> &inarray,
                            Array<OneD, Array<OneD, NekDouble>> &outarray)

{
    int nq = inarray.size() / m_spacedim;

    Array<OneD, Array<OneD, NekDouble>> inarrayMF(m_spacedim);
    for (int i = 0; i < m_spacedim; ++i)
    {
        inarrayMF[i] = Array<OneD, NekDouble>(nq);
    }

    for (int i = 0; i < m_spacedim; ++i)
    {
        Vmath::Vcopy(nq, &inarray[i * nq], 1, &inarrayMF[i][0], 1);
    }

    outarray = Array<OneD, Array<OneD, NekDouble>>(m_spacedim);
    for (int i = 0; i < m_spacedim; ++i)
    {
        outarray[i] = Array<OneD, NekDouble>(nq);
    }

    Array<OneD, NekDouble> tmpx, tmpy, tmpz;
    Array<OneD, NekDouble> Dtmpzdx, Dtmpydx, Dtmpxdy, Dtmpzdy, Dtmpxdz, Dtmpydz;

    tmpx = Array<OneD, NekDouble>(nq);
    tmpy = Array<OneD, NekDouble>(nq);
    tmpz = Array<OneD, NekDouble>(nq);

    Dtmpzdx = Array<OneD, NekDouble>(nq);
    Dtmpydx = Array<OneD, NekDouble>(nq);
    Dtmpxdy = Array<OneD, NekDouble>(nq);
    Dtmpzdy = Array<OneD, NekDouble>(nq);
    Dtmpxdz = Array<OneD, NekDouble>(nq);
    Dtmpydz = Array<OneD, NekDouble>(nq);

    for (int k = 0; k < m_spacedim; ++k)
    {
        Vmath::Vcopy(nq, &inarrayMF[0][0], 1, &tmpx[0], 1);
        Vmath::Vcopy(nq, &inarrayMF[1][0], 1, &tmpy[0], 1);
        Vmath::Vcopy(nq, &inarrayMF[2][0], 1, &tmpz[0], 1);

        // outarray[0] = dz / dy - dy / dz
        // outarray[1] = dx / dz - dz / dx
        // outarray[2] = dy / dx - dx / dy

        m_fields[0]->PhysDeriv(0, tmpz, Dtmpzdx);
        m_fields[0]->PhysDeriv(0, tmpy, Dtmpydx);
        m_fields[0]->PhysDeriv(1, tmpx, Dtmpxdy);
        m_fields[0]->PhysDeriv(1, tmpz, Dtmpzdy);
        m_fields[0]->PhysDeriv(2, tmpx, Dtmpxdz);
        m_fields[0]->PhysDeriv(2, tmpy, Dtmpydz);

        Vmath::Vsub(nq, &Dtmpzdy[0], 1, &Dtmpydz[0], 1, &outarray[0][0], 1);
        Vmath::Vsub(nq, &Dtmpxdz[0], 1, &Dtmpzdx[0], 1, &outarray[1][0], 1);
        Vmath::Vsub(nq, &Dtmpydx[0], 1, &Dtmpxdy[0], 1, &outarray[2][0], 1);
    }
}

Array<OneD, NekDouble> MMFSystem::ComputeEuclideanCurl(
    const Array<OneD, const Array<OneD, NekDouble>> &velocity,
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes)
{
    int nq = GetNpoints();

    Array<OneD, NekDouble> EucCurl(nq, 0.0);

    Array<OneD, NekDouble> Vx(nq), Vy(nq), Vz(nq);
    Vmath::Vcopy(nq, &velocity[0][0], 1, &Vx[0], 1);
    Vmath::Vcopy(nq, &velocity[1][0], 1, &Vy[0], 1);
    Vmath::Vcopy(nq, &velocity[2][0], 1, &Vz[0], 1);

    Array<OneD, NekDouble> D1tmp(nq), D2tmp(nq);
    Array<OneD, NekDouble> Curlx(nq), Curly(nq), Curlz(nq);

    // x component
    m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[1], Vz, D1tmp);
    m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2], Vy, D2tmp);
    Vmath::Vsub(nq, D1tmp, 1, D2tmp, 1, Curlx, 1);

    // y component
    m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2], Vx, D1tmp);
    m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[0], Vz, D2tmp);
    Vmath::Vsub(nq, D1tmp, 1, D2tmp, 1, Curly, 1);

    // z component
    m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[0], Vy, D1tmp);
    m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[1], Vx, D2tmp);
    Vmath::Vsub(nq, D1tmp, 1, D2tmp, 1, Curlz, 1);

    for (int i = 0; i < nq; ++i)
    {
        EucCurl[i] = movingframes[2][i] * Curlx[i] +
                     movingframes[2][i + nq] * Curly[i] +
                     movingframes[2][i + 2 * nq] * Curlz[i];
    }

    return EucCurl;
}

Array<OneD, NekDouble> MMFSystem::ComputeEuclideanGradient(
    const Array<OneD, const NekDouble> &fn)
{
    int nq = GetNpoints();

    Array<OneD, NekDouble> outarray(m_spacedim * nq, 0.0);
    Array<OneD, NekDouble> Dtmp(nq);

    // x component
    for (int i = 0; i < m_spacedim; ++i)
    {
        m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[i], fn, Dtmp);
        Vmath::Vcopy(nq, &Dtmp[0], 1, &outarray[i * nq], 1);
    }

    return outarray;
}

Array<OneD, NekDouble> MMFSystem::ComputeCovGrad(
    const Array<OneD, const NekDouble> &fn,
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes)
{
    int nq = GetNpoints();

    Array<OneD, NekDouble> outarray(3 * nq, 0.0);

    Array<OneD, Array<OneD, NekDouble>> dfdu(m_expdim);
    for (int j=0; j<m_expdim; ++j)
    {
        dfdu[j] = Array<OneD, NekDouble>(nq);
        MMFDirectionalDeriv(movingframes[j], fn, dfdu[j]);
    }

    for (int i = 0; i < nq; ++i)
    {
        for (int k=0; k<m_spacedim; ++k)
        {        
           for (int j=0; j<m_expdim; ++j)
           {
                outarray[i + k * nq] = outarray[i + k * nq] + dfdu[j][i] * movingframes[j][i + k * nq];
           }
        }
    }

    return outarray;
}



    // switch (m_expdim)
    // {
    //     case 1:
    //     {
    //         outarray = ComputeCovGrad1D(fn, movingframes);
    //     }
    //     break;

    //     case 2:
    //     {
    //         outarray = ComputeCovGrad2D(fn, movingframes);
    //     }
    //     break;

    //     case 3:
    //     {
    //         outarray = ComputeCovGrad3D(fn, movingframes);
    //     }
    //     break;

    //     default:
    //         break;
    // }

//     return outarray;
// }

// Array<OneD, NekDouble> MMFSystem::ComputeCovGrad1D(
//     const Array<OneD, const NekDouble> &fn,
//     const Array<OneD, const Array<OneD, NekDouble>> &movingframes)
// {
//     int nq = GetNpoints();

//     Array<OneD, NekDouble> outarray(3 * nq, 0.0);

//     Array<OneD, NekDouble> dfdu1(nq);
//     // m_fields[0]->PhysDirectionalDeriv(movingframes[0], fn, dfdu1);
//     m_fields[0]->PhysDeriv(MultiRegions::eS, fn, dfdu1);

//     // std::cout << "dfdu1 = " << RootMeanSquare(dfdu1) << std::endl;

//     for (int i = 0; i < nq; ++i)
//     {
//         for (int k=0; k<m_spacedim; ++k)
//         {
//             outarray[i+k*nq] = dfdu1[i] * movingframes[0][i+k*nq];
//         }
//     }

//     return outarray;
// }

// // Compute Covariant gradent: \nabla f
// Array<OneD, NekDouble> MMFSystem::ComputeCovGrad2D(
//     const Array<OneD, const NekDouble> &fn,
//     const Array<OneD, const Array<OneD, NekDouble>> &movingframes)
// {
//     int nq = GetNpoints();

//     Array<OneD, NekDouble> outarray(3 * nq, 0.0);

//     Array<OneD, NekDouble> dfdu1(nq);
//     Array<OneD, NekDouble> dfdu2(nq);

//     MMFDirectionalDeriv(movingframes[0], fn, dfdu1);
//     MMFDirectionalDeriv(movingframes[1], fn, dfdu2);

//     NekDouble dfde1 = 0.0, dfde2 = 0.0;
//     for (int i = 0; i < nq; ++i)
//     {
//         dfde1 = dfdu1[i];
//         dfde2 = dfdu2[i];

//         for (int k=0; k<m_spacedim; ++k)
//         {
//             outarray[i+k*nq] = dfde1 * movingframes[0][i+k*nq] + dfde2 * movingframes[1][i+k*nq];
//         }
//     }

//     return outarray;
// }

// // Compute Covariant gradent: \nabla f
// Array<OneD, NekDouble> MMFSystem::ComputeCovGrad3D(
//     const Array<OneD, const NekDouble> &fn,
//     const Array<OneD, const Array<OneD, NekDouble>> &movingframes)
// {
//     int nq = GetNpoints();

//     Array<OneD, NekDouble> outarray(3 * nq, 0.0);

//     Array<OneD, NekDouble> dfdu1(nq);
//     Array<OneD, NekDouble> dfdu2(nq);
//     Array<OneD, NekDouble> dfdu3(nq);

//     MMFDirectionalDeriv(movingframes[0], fn, dfdu1);
//     MMFDirectionalDeriv(movingframes[1], fn, dfdu2);
//     MMFDirectionalDeriv(movingframes[2], fn, dfdu3);

//     NekDouble dfde1 = 0.0, dfde2 = 0.0, dfde3 = 0.0;
//     for (int i = 0; i < nq; ++i)
//     {
//         dfde1 = dfdu1[i];
//         dfde2 = dfdu2[i];
//         dfde3 = dfdu3[i];

//         for (int k=0; k<m_spacedim; ++k)
//         {
//             outarray[i + k * nq] = dfde1 * movingframes[0][i + k * nq] +
//                                 dfde2 * movingframes[1][i + k * nq] +
//                                 dfde3 * movingframes[2][i + k * nq];
//         }
//     }

//     return outarray;
// }

// Curlike gradient
Array<OneD, NekDouble> MMFSystem::ComputeCovJGrad(
    const Array<OneD, const NekDouble> &fn,
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes)
{
    int nq = GetNpoints();

    Array<OneD, NekDouble> outarray(3 * nq, 0.0);

    Array<OneD, NekDouble> dfdu1(nq);
    Array<OneD, NekDouble> dfdu2(nq);

    MMFDirectionalDeriv(movingframes[0], fn, dfdu1);
    MMFDirectionalDeriv(movingframes[1], fn, dfdu2);

    NekDouble dfde1 = 0.0, dfde2 = 0.0;
    for (int i = 0; i < nq; ++i)
    {
        dfde1 = dfdu1[i];
        dfde2 = dfdu2[i];

        outarray[i] = -dfde2 * movingframes[0][i] + dfde1 * movingframes[1][i];

        outarray[i + nq] =
            -dfde2 * movingframes[0][i + nq] + dfde1 * movingframes[1][i + nq];

        outarray[i + 2 * nq] = -dfde2 * movingframes[0][i + 2 * nq] +
                               dfde1 * movingframes[1][i + 2 * nq];
    }

    return outarray;
}

// Compute Direct Euclidean Laplacian
Array<OneD, NekDouble> MMFSystem::ComputeEuclideanDiffusion(
    const Array<OneD, const NekDouble> &inarray)
{
    int nq = m_fields[0]->GetNpoints();

    Array<OneD, NekDouble> tmp(nq);
    Array<OneD, NekDouble> Dtmp(nq);
    Array<OneD, NekDouble> DirEuc(nq, 0.0);
    for (int i = 0; i < m_spacedim; ++i)
    {
        m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[i], inarray, tmp);
        m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[i], tmp, Dtmp);
        Vmath::Vadd(nq, Dtmp, 1, DirEuc, 1, DirEuc, 1);
    }

    return DirEuc;
}

// \nabla^2 u = \nabla u^i \cdot e^i + u^i (\nabla \cdot e^i)
Array<OneD, NekDouble> MMFSystem::ComputeMMFDiffusion(
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    const Array<OneD, const NekDouble> &inarray)
{
    int nq = m_fields[0]->GetNpoints();

    Array<OneD, NekDouble> tmp(nq);
    Array<OneD, NekDouble> Dtmp(nq);

    Array<OneD, NekDouble> outarray(nq, 0.0);
    for (int i = 0; i < m_expdim; ++i)
    {
        // Dtmp = \nabla u \cdot e^i
        // D2tmp = \nabla Dtmp \cdot e^i
        MMFDirectionalDeriv(movingframes[i], inarray, tmp);
        MMFDirectionalDeriv(movingframes[i], tmp, Dtmp);

        Vmath::Vadd(nq, Dtmp, 1, outarray, 1, outarray, 1);
    }

    return outarray;
}

void MMFSystem::MMFDirectionalDeriv(const Array<OneD, const NekDouble> &movingframe, 
                                    const Array<OneD, const NekDouble> &inarray, 
                                    Array<OneD, NekDouble> &outarray)
{
    int nq = m_fields[0]->GetNpoints();

    Array<OneD, NekDouble> tmp(nq);

    outarray = Array<OneD, NekDouble>(nq, 0.0);
    for (int k=0; k<m_spacedim; ++k)
    {
        m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[k], inarray, tmp);
        Vmath::Vvtvp(nq, &movingframe[k*nq], 1, &tmp[0], 1, &outarray[0], 1, &outarray[0], 1);
    }
}


// Array<OneD, NekDouble> MMFSystem::ComputeMMFDiffusion(
//     const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
//     const Array<OneD, const NekDouble> &inarray)
// {
//     int nq = m_fields[0]->GetNpoints();

//     Array<OneD, NekDouble> tmp(nq);
//     Array<OneD, NekDouble> Dtmp(nq);

//     Array<OneD, NekDouble> outarray(nq, 0.0);
//     for (int i = 0; i < m_expdim; ++i)
//     {
//         // Dtmp = \nabla u \cdot e^i
//         // D2tmp = \nabla Dtmp \cdot e^i
//         m_fields[0]->PhysDirectionalDeriv(movingframes[i], inarray, tmp);
//         m_fields[0]->PhysDirectionalDeriv(movingframes[i], tmp, Dtmp);

//         Vmath::Vadd(nq, Dtmp, 1, outarray, 1, outarray, 1);
//     }

//     return outarray;
// }

// \nabla^2 u = \nabla u^i \cdot e^i + u^i (\nabla \cdot e^i)
Array<OneD, NekDouble> MMFSystem::ComputeCovariantDiffusion(
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    const Array<OneD, const NekDouble> &inarray, const DerivType DType)
{
    boost::ignore_unused(DType);

    int nq = m_fields[0]->GetNpoints();

    Array<OneD, NekDouble> Dtmp(nq);
    Array<OneD, NekDouble> D2tmp(nq);

    Array<OneD, Array<OneD, NekDouble>> DivMFCov;
    ComputeDivMF(eCovariant, movingframes, DivMFCov);

    // SetBoundaryConditions(0.0);
    Array<OneD, NekDouble> outarray(nq, 0.0);
    for (int i = 0; i < m_expdim; ++i)
    {
        // Dtmp = \nabla u \cdot e^i
        // D2tmp = \nabla Dtmp \cdot e^i
        // m_fields[0]->PhysDirectionalDeriv(movingframes[i], inarray, Dtmp);
        // m_fields[0]->PhysDirectionalDeriv(movingframes[i], Dtmp, D2tmp);
        MMFDirectionalDeriv(movingframes[i], inarray, Dtmp);
        MMFDirectionalDeriv(movingframes[i], Dtmp, D2tmp);

        Vmath::Vvtvp(nq, Dtmp, 1, DivMFCov[i], 1, D2tmp, 1, D2tmp, 1);
        Vmath::Vadd(nq, D2tmp, 1, outarray, 1, outarray, 1);
    }

    return outarray;
}

Array<OneD, NekDouble> MMFSystem::CartesianToMovingframes(
    const Array<OneD, const Array<OneD, NekDouble>> movingframes,
    const Array<OneD, const Array<OneD, NekDouble>> &uvec, unsigned int field)
{
    int nq = m_fields[0]->GetNpoints();

    Array<OneD, NekDouble> outarray(nq, 0.0);

    // new u0 = ( [u v] \cdot e^1 )/ |e^1|^2
    ////////////////////////////////////////////////////////////////////////////////////
    // Vmath::Vmul(nq, &movingframes[field][0], 1, &uvec[0][0], 1, &outarray[0],
    //             1);
    // Vmath::Vvtvp(nq, &movingframes[field][nq], 1, &uvec[1][0], 1,
    // &outarray[0],
    //              1, &outarray[0], 1);
    // Vmath::Vvtvp(nq, &movingframes[field][2 * nq], 1, &uvec[2][0], 1,
    //              &outarray[0], 1, &outarray[0], 1);
    ////////////////////////////////////////////////////////////////////////////////////
    for (int k = 0; k < m_spacedim; ++k)
    {
        Vmath::Vvtvp(nq, &movingframes[field][k * nq], 1, &uvec[k][0], 1,
                     &outarray[0], 1, &outarray[0], 1);
    }

    return outarray;
}

// v_1 \vec{e}^1 + v_2 \vec{e}^2 = v_x \vec{x} + v_y \vec{y} + v_z
// v_i = v_x e^i_x + v_y e^i_y + v_z e^i_z
void MMFSystem::CartesianToMovingframes(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    const Array<OneD, const Array<OneD, NekDouble>> &MF1st,
    Array<OneD, Array<OneD, NekDouble>> &outarray)
{
    int nq = m_fields[0]->GetNpoints();

    outarray = Array<OneD, Array<OneD, NekDouble>>(m_shapedim);
    for (int j = 0; j < m_shapedim; ++j)
    {
        outarray[j] = Array<OneD, NekDouble>(nq, 0.0);
        for (int i = 0; i < m_spacedim; ++i)
        {
            Vmath::Vvtvp(nq, &inarray[i][0], 1, &MF1st[j][i * nq], 1,
                         &outarray[j][0], 1, &outarray[j][0], 1);
        }
    }
}

void MMFSystem::CartesianToMovingframes(
    const Array<OneD, const NekDouble> &inarrayx,
    const Array<OneD, const NekDouble> &inarrayy,
    const Array<OneD, const NekDouble> &inarrayz,
    const Array<OneD, const Array<OneD, NekDouble>> &MF1st,
    Array<OneD, NekDouble> &outarray1, Array<OneD, NekDouble> &outarray2)
{
    int nq = m_fields[0]->GetNpoints();

    outarray1 = Array<OneD, NekDouble>(nq);
    outarray2 = Array<OneD, NekDouble>(nq);

    // e^1 component
    Vmath::Vmul(nq, &inarrayx[0], 1, &MF1st[0][0], 1, &outarray1[0], 1);
    Vmath::Vvtvp(nq, &inarrayy[0], 1, &MF1st[0][nq], 1, &outarray1[0], 1,
                 &outarray1[0], 1);
    Vmath::Vvtvp(nq, &inarrayz[0], 1, &MF1st[0][2 * nq], 1, &outarray1[0], 1,
                 &outarray1[0], 1);

    // e^2 component
    Vmath::Vmul(nq, &inarrayx[0], 1, &MF1st[1][0], 1, &outarray2[0], 1);
    Vmath::Vvtvp(nq, &inarrayy[0], 1, &MF1st[1][nq], 1, &outarray2[0], 1,
                 &outarray2[0], 1);
    Vmath::Vvtvp(nq, &inarrayz[0], 1, &MF1st[1][2 * nq], 1, &outarray2[0], 1,
                 &outarray2[0], 1);
}

void MMFSystem::Cart_to_MF(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    const Array<OneD, const Array<OneD, NekDouble>> &MF1st,
    Array<OneD, Array<OneD, NekDouble>> &outarray)
{
    Cart_to_MF(inarray, MF1st, outarray[0], outarray[1]);
}

void MMFSystem::Cart_to_MF(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    const Array<OneD, const Array<OneD, NekDouble>> &MF1st,
    Array<OneD, NekDouble> &outarray1, Array<OneD, NekDouble> &outarray2)
{
    int nq = m_fields[0]->GetNpoints();

    outarray1 = Array<OneD, NekDouble>(nq, 0.0);
    outarray2 = Array<OneD, NekDouble>(nq, 0.0);

    // e^1 component
    for (int k = 0; k < m_spacedim; ++k)
    {
        Vmath::Vvtvp(nq, &inarray[k][0], 1, &MF1st[0][k * nq], 1, &outarray1[0],
                     1, &outarray1[0], 1);
        Vmath::Vvtvp(nq, &inarray[k][0], 1, &MF1st[1][k * nq], 1, &outarray2[0],
                     1, &outarray2[0], 1);
    }
}

void MMFSystem::Sph_to_Cart(const NekDouble &xj, const NekDouble &yj,
                            const NekDouble &zj, const NekDouble &uth,
                            const NekDouble &uphi, NekDouble &ux, NekDouble &uy,
                            NekDouble &uz)
{
    NekDouble sin_varphi, cos_varphi, sin_theta, cos_theta;

    CartesianToNewSpherical(xj, yj, zj, sin_varphi, cos_varphi, sin_theta,
                            cos_theta);

    ux = uth * cos_theta * cos_varphi - uphi * sin_varphi;
    uy = uth * cos_theta * sin_varphi + uphi * cos_varphi;
    uz = -uth * sin_theta;
}

// u^th e^th + u^phi e^phi = u^1 e^1 + u^2 e^2
// u^th = u^1 * (e^1 \cdot e^th) + u^2 * (e^2 \cdot e^th)
// u^phi = u^1 * (e^1 \cdot e^phi) + u^2 * (e^2 \cdot e^phi)
void MMFSystem::MF_to_Sph(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    const Array<OneD, const Array<OneD, NekDouble>> &MF1st,
    Array<OneD, Array<OneD, NekDouble>> &outarray)
{
    int nq = m_fields[0]->GetNpoints();

    for (int j = 0; j < m_shapedim; ++j)
    {
        outarray[j] = Array<OneD, NekDouble>(nq, 0.0);
    }

    Array<OneD, Array<OneD, NekDouble>> SphericalVector;
    ComputeSphericalVector(SphericalVector);

    Array<OneD, NekDouble> ThetacdotMF(nq);
    Array<OneD, NekDouble> PhicdotMF(nq);
    for (int j = 0; j < m_shapedim; ++j)
    {
        ThetacdotMF = Array<OneD, NekDouble>(nq, 0.0);
        PhicdotMF   = Array<OneD, NekDouble>(nq, 0.0);

        // Compute e^j \cdot \vec{theta} = ThetacdotMF;
        // Compute e^j \cdot \vec{phi} = PhicdotMF;
        for (int i = 0; i < m_spacedim; ++i)
        {
            Vmath::Vvtvp(nq, &MF1st[j][i * nq], 1, &SphericalVector[0][i * nq],
                         1, &ThetacdotMF[0], 1, &ThetacdotMF[0], 1);
            Vmath::Vvtvp(nq, &MF1st[j][i * nq], 1, &SphericalVector[1][i * nq],
                         1, &PhicdotMF[0], 1, &PhicdotMF[0], 1);
        }

        // j = 0
        Vmath::Vvtvp(nq, &inarray[j][0], 1, &ThetacdotMF[0], 1, &outarray[0][0],
                     1, &outarray[0][0], 1);
        Vmath::Vvtvp(nq, &inarray[j][0], 1, &PhicdotMF[0], 1, &outarray[1][0],
                     1, &outarray[1][0], 1);
    }
}

// \mathbf{v} = v_1 \mathbf{e}^1 + v_2 \mathbf{e}^2
// v_i = \mathbf{v} \cdot \mathbf{e}^i
void MMFSystem::vector_to_vcoeff(
    const Array<OneD, const Array<OneD, NekDouble>> &vector,
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    Array<OneD, Array<OneD, NekDouble>> &vcoeff)
{
    int nq = m_fields[0]->GetNpoints();

    vcoeff = Array<OneD, Array<OneD, NekDouble>>(m_expdim);
    for (int j=0; j<m_expdim; ++j)
    {
        vcoeff[j] = Array<OneD, NekDouble>(nq, 0.0);
        for (int i=0; i<nq; ++i)
        {
            for (int k=0; k<m_spacedim; ++k)
            {
                vcoeff[j][i] = vcoeff[j][i] + vector[k][i] * movingframes[j][i+k*nq];
            }
        }
    }
}

// \mathbf{v} = v_1 \mathbf{e}^1 + v_2 \mathbf{e}^2
void MMFSystem::vcoeff_to_vector(
    const Array<OneD, const Array<OneD, NekDouble>> &vcoeff,
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    Array<OneD, Array<OneD, NekDouble>> &vector)
{
    int nq = m_fields[0]->GetNpoints();

    vector = Array<OneD, Array<OneD, NekDouble>>(m_spacedim);
    for (int k=0; k<m_spacedim; ++k)
    {
        vector[k] = Array<OneD, NekDouble>(nq, 0.0);
        for (int i=0; i<nq; ++i)
        {
            for (int j=0; j<m_expdim; ++j)
            {
                vector[k][i] = vector[k][i] + vcoeff[j][i] * movingframes[j][i+k*nq];
            }
        }
    }
}


// v_1 \vec{e}^1 + v_2 \vec{e}^2 = v_x \vec{x} + v_y \vec{y} + v_z
// \vec{z} = \vec{v} v_x = \vec{v} \cdot \vec{x} = v_1 e^1_x + v_2 e^2_x
// v_y = \vec{v} \cdot \vec{y} = v_1 e^1_y + v_2 e^2_y v_z = \vec{v}
// \cdot \vec{z} = v_1 e^1_z
// + v_2 e^2_z
void MMFSystem::MovingframestoCartesian(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    const Array<OneD, const Array<OneD, NekDouble>> &MF1st,
    Array<OneD, Array<OneD, NekDouble>> &outarray)
{
    int nq = m_fields[0]->GetNpoints();

    outarray = Array<OneD, Array<OneD, NekDouble>>(m_spacedim);
    for (int i = 0; i < m_spacedim; ++i)
    {
        outarray[i] = Array<OneD, NekDouble>(nq, 0.0);
        for (int j = 0; j < m_shapedim; ++j)
        {
            Vmath::Vvtvp(nq, &inarray[j][0], 1, &MF1st[j][i * nq], 1,
                         &outarray[i][0], 1, &outarray[i][0], 1);
        }
    }
}

void MMFSystem::MovingframestoCartesian(
    const Array<OneD, const NekDouble> &inarray1,
    const Array<OneD, const NekDouble> &inarray2,
    const Array<OneD, const Array<OneD, NekDouble>> &MF1st,
    Array<OneD, NekDouble> &outarrayx, Array<OneD, NekDouble> &outarrayy,
    Array<OneD, NekDouble> &outarrayz)
{
    int nq = m_fields[0]->GetNpoints();

    outarrayx = Array<OneD, NekDouble>(nq);
    outarrayy = Array<OneD, NekDouble>(nq);
    outarrayz = Array<OneD, NekDouble>(nq);

    // x component
    Vmath::Vmul(nq, &inarray1[0], 1, &MF1st[0][0], 1, &outarrayx[0], 1);
    Vmath::Vvtvp(nq, &inarray2[0], 1, &MF1st[1][0], 1, &outarrayx[0], 1,
                 &outarrayx[0], 1);

    // y component
    Vmath::Vmul(nq, &inarray1[0], 1, &MF1st[0][nq], 1, &outarrayy[0], 1);
    Vmath::Vvtvp(nq, &inarray2[0], 1, &MF1st[1][nq], 1, &outarrayy[0], 1,
                 &outarrayy[0], 1);

    // z component
    Vmath::Vmul(nq, &inarray1[0], 1, &MF1st[0][2 * nq], 1, &outarrayz[0], 1);
    Vmath::Vvtvp(nq, &inarray2[0], 1, &MF1st[1][2 * nq], 1, &outarrayz[0], 1,
                 &outarrayz[0], 1);
}

void MMFSystem::ProjectionOntoMovingFrames(
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    const Array<OneD, const Array<OneD, NekDouble>> &velocity,
    Array<OneD, Array<OneD, NekDouble>> &outarray)
{
    int nq = m_fields[0]->GetNpoints();

    outarray = Array<OneD, Array<OneD, NekDouble>>(m_spacedim);
    for (int k = 0; k < m_spacedim; ++k)
    {
        outarray[k] = Array<OneD, NekDouble>(nq);
    }

    //  Get u and v from \vec{v} = u e^1 + v e^2
    Array<OneD, NekDouble> uvec;
    Array<OneD, NekDouble> vvec;
    CartesianToMovingframes(velocity[0], velocity[1], velocity[2], movingframes,
                            uvec, vvec);

    // Change it into \vec{v} from u and v
    Array<OneD, NekDouble> velx;
    Array<OneD, NekDouble> vely;
    Array<OneD, NekDouble> velz;
    MovingframestoCartesian(uvec, vvec, movingframes, velx, vely, velz);

    Vmath::Vcopy(nq, &velx[0], 1, &outarray[0][0], 1);
    Vmath::Vcopy(nq, &vely[0], 1, &outarray[1][0], 1);
    Vmath::Vcopy(nq, &velz[0], 1, &outarray[2][0], 1);
}

// void MMFSystem::MovingframesToSpherical(
//     const Array<OneD, const NekDouble> &u1,
//     const Array<OneD, const NekDouble> &u2,
//     const Array<OneD, const Array<OneD, NekDouble>> &MF1st,
//     Array<OneD, NekDouble> &outtheta,
//     Array<OneD, NekDouble> &outphi)
// {
//     int nq = m_fields[0]->GetNpoints();

//     Array<OneD, NekDouble> uvec(m_shapedim * nq);
//     Vmath::Vcopy(nq, &u1[0], 1, &uvec[0], 1);
//     Vmath::Vcopy(nq, &u2[0], 1, &uvec[nq], 1);

//     MovingframesToSpherical(uvec, MF1st, outtheta, outphi);
// }

// void MMFSystem::MovingframesToSpherical(
//     const Array<OneD, const NekDouble> &inarray,
//     const Array<OneD, const Array<OneD, NekDouble>> &MF1st,
//     Array<OneD, NekDouble> &outtheta,
//     Array<OneD, NekDouble> &outphi)
// {
//     int nq = m_fields[0]->GetNpoints();

//     Array<OneD, Array<OneD, NekDouble>> SphericalVector;
//     ComputeSphericalVector(SphericalVector);

//     Array<OneD, NekDouble> ThetacdotMF(nq);
//     Array<OneD, NekDouble> PhicdotMF(nq);

//     outtheta = Array<OneD, NekDouble>(nq, 0.0);
//     outphi = Array<OneD, NekDouble>(nq, 0.0);
//     for (int j=0; j<m_shapedim; ++j)
//     {
//         ThetacdotMF = Array<OneD, NekDouble>(nq, 0.0);
//         PhicdotMF = Array<OneD, NekDouble>(nq, 0.0);

//         // Compute e_j \cdot \theta and e_j \cdot \phi
//         for (int k = 0; k < m_spacedim; ++k)
//         {
//             Vmath::Vvtvp(nq, &MF1st[j][k * nq], 1, &SphericalVector[0][k *
//             nq], 1,
//                         &ThetacdotMF[0], 1, &ThetacdotMF[0], 1);

//             Vmath::Vvtvp(nq, &MF1st[j][k * nq], 1, &SphericalVector[1][k *
//             nq], 1,
//                         &PhicdotMF[0], 1, &PhicdotMF[0], 1);
//         }

//         Vmath::Vvtvp(nq, &inarray[j*nq], 1, &ThetacdotMF[0], 1, &outtheta[0],
//         1, &outtheta[0], 1); Vmath::Vvtvp(nq, &inarray[j*nq], 1,
//         &PhicdotMF[0], 1, &outphi[0], 1, &outphi[0], 1);
//     }
// }

// Find u^{\theta}, u^{\phi}
// u^1 e_1 + u^2 e_2 = u^{theta} \theta + u^{\phi} \phi
// inarray (m_shapedim*nq) = [ u^1 u^2 ]
// outarray: outtheta = u^{\theta}, outphi = u^{\phi}
void MMFSystem::MovingframesToSpherical(
    const Array<OneD, const NekDouble> &u1,
    const Array<OneD, const NekDouble> &u2,
    const Array<OneD, const Array<OneD, NekDouble>> &MF1st,
    Array<OneD, NekDouble> &outtheta, Array<OneD, NekDouble> &outphi)
{
    int nq = m_fields[0]->GetNpoints();

    Array<OneD, Array<OneD, NekDouble>> uvec(m_shapedim);
    for (int j = 0; j < m_shapedim; ++j)
    {
        uvec[j] = Array<OneD, NekDouble>(nq);
    }

    Vmath::Vcopy(nq, &u1[0], 1, &uvec[0][0], 1);
    Vmath::Vcopy(nq, &u2[0], 1, &uvec[1][0], 1);

    MovingframesToSpherical(uvec, MF1st, outtheta, outphi);
}

void MMFSystem::MovingframesToSpherical(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    const Array<OneD, const Array<OneD, NekDouble>> &MF1st,
    Array<OneD, NekDouble> &outtheta, Array<OneD, NekDouble> &outphi)
{
    int nq = m_fields[0]->GetNpoints();

    Array<OneD, Array<OneD, NekDouble>> SphericalVector;
    ComputeSphericalVector(SphericalVector);

    // Array<OneD, NekDouble> ThetacdotMF(nq);
    // Array<OneD, NekDouble> PhicdotMF(nq);

    Array<OneD, Array<OneD, NekDouble>> ThetacdotMF(m_shapedim);
    Array<OneD, Array<OneD, NekDouble>> PhicdotMF(m_shapedim);

    for (int j = 0; j < m_shapedim; ++j)
    {
        ThetacdotMF[j] = Array<OneD, NekDouble>(nq, 0.0);
        PhicdotMF[j]   = Array<OneD, NekDouble>(nq, 0.0);

        for (int i = 0; i < m_spacedim; ++i)
        {
            Vmath::Vvtvp(nq, &MF1st[j][i * nq], 1, &SphericalVector[0][i * nq],
                         1, &ThetacdotMF[j][0], 1, &ThetacdotMF[j][0], 1);
            Vmath::Vvtvp(nq, &MF1st[j][i * nq], 1, &SphericalVector[1][i * nq],
                         1, &PhicdotMF[j][0], 1, &PhicdotMF[j][0], 1);
        }
    }

    outtheta = Array<OneD, NekDouble>(nq, 0.0);
    outphi   = Array<OneD, NekDouble>(nq, 0.0);

    // j = 0
    Vmath::Vmul(nq, &inarray[0][0], 1, &ThetacdotMF[0][0], 1, &outtheta[0], 1);
    Vmath::Vvtvp(nq, &inarray[1][0], 1, &ThetacdotMF[1][0], 1, &outtheta[0], 1,
                 &outtheta[0], 1);

    // j = 1
    Vmath::Vmul(nq, &inarray[0][0], 1, &PhicdotMF[0][0], 1, &outphi[0], 1);
    Vmath::Vvtvp(nq, &inarray[1][0], 1, &PhicdotMF[1][0], 1, &outphi[0], 1,
                 &outphi[0], 1);
}

// x = r \sin \theta \cos \varphi
// y = r \sin \theta \sin \varphi
// z = r \cos \theta
// \theta \in [0, \pi], \varphi = [0, 2 \pi]
// sin \theta <-> cos \theta from CartesianToSpherical
void MMFSystem::CartesianToNewSpherical(
    const NekDouble &x0j, const NekDouble &x1j, const NekDouble &x2j,
    NekDouble &sin_varphi, NekDouble &cos_varphi, NekDouble &sin_theta,
    NekDouble &cos_theta)
{
    NekDouble radius, radxy;
    NekDouble Tol = 0.0000001;

    radius = sqrt(x0j * x0j + x1j * x1j + x2j * x2j);
    radxy  = sqrt(x0j * x0j + x1j * x1j);

    // cos_theta & sin_theta
    sin_theta = radxy / radius;
    cos_theta = x2j / radius;

    if (radxy > Tol)
    {
        sin_varphi = x1j / radxy;
        cos_varphi = x0j / radxy;
    }

    else
    {
        sin_varphi = 0.0;
        if (x2j > 0)
        {
            cos_varphi = 1.0;
        }

        else
        {
            cos_varphi = -1.0;
        }
    }
}

// x = \sech \theta \cos \phi
// y = \sech \theta \sin \phi
// z = \theta - \tamh \theta

// \sech_\theta = \sqrt{ x*x + y*y }
//
void MMFSystem::CartesianToPseudospherical(
    const NekDouble &x0, const NekDouble &x1, const NekDouble &x2,
    NekDouble &sin_varphi, NekDouble &cos_varphi, NekDouble &theta,
    NekDouble &sech_theta, NekDouble &tanh_theta)
{
    NekDouble rad;
    NekDouble Tol = 0.0000000001;
    NekDouble exp_theta;

    rad = sqrt(x0 * x0 + x1 * x1);

    sin_varphi = x1 / rad;
    cos_varphi = x0 / rad;

    sech_theta = rad;
    if (fabs(x2) > Tol)
    {
        exp_theta = (1.0 + sqrt(1.0 - x0 * x0 - x1 * x1)) / rad;

        tanh_theta =
            (exp_theta * exp_theta - 1.0) / (exp_theta * exp_theta + 1.0);
        theta = x2 + tanh_theta;
    }

    else
    {
        tanh_theta = 0.0;
        theta      = 0.0;
    }
}

// x = r \cos \theta \cos \varphi
// y = r \cos \theta \sin \varphi
// z = r \sin \theta
// void MMFSystem::CartesianToSpherical(const NekDouble x0j, const NekDouble
// x1j,
//                                      const NekDouble x2j, NekDouble
//                                      &sin_varphi, NekDouble &cos_varphi,
//                                      NekDouble &sin_theta, NekDouble
//                                      &cos_theta)
// {
//     NekDouble radius;
//     NekDouble radxy;

//     radius = sqrt(x0j * x0j + x1j * x1j + x2j * x2j);
//     radxy  = sqrt(x0j * x0j + x1j * x1j);

//     NekDouble Tol = 0.0000001;
//     if (radxy > Tol)
//     {
//         sin_varphi = x1j / radxy;
//         cos_varphi = x0j / radxy;
//     }

//     else
//     {
//         sin_varphi = 0.0;
//         if (x2j > 0)
//         {
//             cos_varphi = 1.0;
//         }

//         else
//         {
//             cos_varphi = -1.0;
//         }
//     }

//     if (radius > Tol)
//     {
//         sin_theta = x2j / radius;
//         cos_theta = radxy / radius;
//     }

//     else
//     {
//         sin_theta = 0;
//         cos_theta = 1.0;
//     }
// }

// x = a \cos \theta \cos \varphi
// y = b \cos \theta \sin \varphi
// z = c \sin \theta
void MMFSystem::CartesianToElliptical(const NekDouble x0j, const NekDouble x1j,
                                      const NekDouble x2j, NekDouble &rad,
                                      NekDouble &sin_varphi,
                                      NekDouble &cos_varphi,
                                      NekDouble &sin_theta,
                                      NekDouble &cos_theta)
{
    // aoverrad = a/c, boverrad = b/c
    NekDouble Radx = 1.0;
    NekDouble Rady = 1.0;
    NekDouble Radz = 1.0;

    switch (m_surfaceType)
    {
        case SolverUtils::eSphere:
        {
            Radx = 1.0;
            Rady = 1.0;
            Radz = 1.0;
        }
        break;

        case SolverUtils::eEllipsoid:
        {
            Radx = m_Radx;
            Rady = m_Rady;
            Radz = m_Radz;
        }
        break;

        default:
            break;
    }

    // sqrt{ (x/a)^2 + (y/b)^2 + (z/c)^2 } = 1
    rad = sqrt(x0j * x0j / Radx / Radx + x1j * x1j / Rady / Rady +
               x2j * x2j / Radz / Radz);

    // sin_theta = (1/c) \sqrt { (c/a)^2 x^2 + (c/b)^2 y^2 }
    sin_theta = sqrt(x0j * x0j / Radx / Radx + x1j * x1j / Rady / Rady);
    cos_theta = x2j / Radz;

    NekDouble Tol = 0.0000001;
    if (sin_theta > Tol)
    {
        cos_varphi = x0j / (Radx * sin_theta);
        sin_varphi = x1j / (Rady * sin_theta);
    }

    else
    {
        sin_varphi = 0.0;
        if (x2j > 0)
        {
            cos_varphi = 1.0;
        }

        else
        {
            cos_varphi = -1.0;
        }
    }
}

void MMFSystem::ComputeSphericalTangentVector(Array<OneD, NekDouble> &Phi,
                                              Array<OneD, NekDouble> &Theta)
{
    int nq = m_fields[0]->GetNpoints();

    NekDouble sin_varphi, cos_varphi, sin_theta, cos_theta;
    NekDouble x0j, x1j, x2j, rad;

    Array<OneD, NekDouble> x0(nq);
    Array<OneD, NekDouble> x1(nq);
    Array<OneD, NekDouble> x2(nq);

    m_fields[0]->GetCoords(x0, x1, x2);

    Theta = Array<OneD, NekDouble>(m_spacedim * nq);
    Phi   = Array<OneD, NekDouble>(m_spacedim * nq);

    for (int i = 0; i < nq; i++)
    {
        x0j = x0[i];
        x1j = x1[i];
        x2j = x2[i];

        rad = sqrt(x0j * x0j + x1j * x1j + x2j * x2j);

        CartesianToNewSpherical(x0j, x1j, x2j, sin_varphi, cos_varphi,
                                sin_theta, cos_theta);

        Theta[i]          = rad * cos_theta * cos_varphi;
        Theta[i + nq]     = rad * cos_theta * sin_varphi;
        Theta[i + 2 * nq] = -rad * sin_theta;

        Phi[i]          = -rad * sin_varphi * sin_theta;
        Phi[i + nq]     = rad * cos_varphi * sin_theta;
        Phi[i + 2 * nq] = 0.0;
    }
}

void MMFSystem::ComputeTangentUnitVector(
    Array<OneD, Array<OneD, NekDouble>> &TangentUnitVector)
{
    switch (m_surfaceType)
    {
        case SolverUtils::eSphere:
        {
            ComputeSphericalVector(TangentUnitVector);
        }
        break;

        case SolverUtils::eEllipsoid:
        {
            ComputeEllipsoidVector(TangentUnitVector);
        }
        break;

        default:
            break;
    }
}

void MMFSystem::ComputeSphericalVector(
    Array<OneD, Array<OneD, NekDouble>> &SphericalVector)
{
    int nq = m_fields[0]->GetNpoints();

    NekDouble sin_varphi, cos_varphi, sin_theta, cos_theta;
    NekDouble x0j, x1j, x2j;

    Array<OneD, NekDouble> x0(nq);
    Array<OneD, NekDouble> x1(nq);
    Array<OneD, NekDouble> x2(nq);

    m_fields[0]->GetCoords(x0, x1, x2);

    SphericalVector = Array<OneD, Array<OneD, NekDouble>>(m_mfdim);
    for (int i = 0; i < m_mfdim; ++i)
    {
        SphericalVector[i] = Array<OneD, NekDouble>(m_spacedim * nq, 0.0);
    }

    for (int i = 0; i < nq; i++)
    {
        x0j = x0[i];
        x1j = x1[i];
        x2j = x2[i];

        CartesianToNewSpherical(x0j, x1j, x2j, sin_varphi, cos_varphi,
                                sin_theta, cos_theta);

        // Theta - direction
        SphericalVector[0][i]          = cos_theta * cos_varphi;
        SphericalVector[0][i + nq]     = cos_theta * sin_varphi;
        SphericalVector[0][i + 2 * nq] = -1.0 * sin_theta;

        // Phi - direction
        SphericalVector[1][i]          = -1.0 * sin_varphi;
        SphericalVector[1][i + nq]     = cos_varphi;
        SphericalVector[1][i + 2 * nq] = 0.0;

        // r - direction
        SphericalVector[2][i]          = sin_theta * cos_varphi;
        SphericalVector[2][i + nq]     = sin_theta * sin_varphi;
        SphericalVector[2][i + 2 * nq] = cos_theta;
    }
}

void MMFSystem::ComputeEllipsoidVector(
    Array<OneD, Array<OneD, NekDouble>> &EllipsoidVector)
{
    int nq = m_fields[0]->GetNpoints();

    NekDouble sin_varphi, cos_varphi, sin_theta, cos_theta;
    NekDouble x0j, x1j, x2j;

    Array<OneD, NekDouble> x0(nq);
    Array<OneD, NekDouble> x1(nq);
    Array<OneD, NekDouble> x2(nq);

    m_fields[0]->GetCoords(x0, x1, x2);

    EllipsoidVector = Array<OneD, Array<OneD, NekDouble>>(m_mfdim);
    for (int i = 0; i < m_mfdim; ++i)
    {
        EllipsoidVector[i] = Array<OneD, NekDouble>(m_spacedim * nq, 0.0);
    }

    NekDouble rad, Rada, Radb, Radc;

    Rada = m_Radx;
    Radb = m_Rady;
    Radc = m_Radz;

    NekDouble tmpx, tmpy, tmpz, mag;
    for (int i = 0; i < nq; i++)
    {
        x0j = x0[i];
        x1j = x1[i];
        x2j = x2[i];

        CartesianToElliptical(x0j, x1j, x2j, rad, sin_varphi, cos_varphi,
                              sin_theta, cos_theta);

        // Theta - direction
        tmpx = Rada * cos_theta * cos_varphi;
        tmpy = Radb * cos_theta * sin_varphi;
        tmpz = -1.0 * Radc * sin_theta;

        mag = sqrt(tmpx * tmpx + tmpy * tmpy + tmpz * tmpz);

        EllipsoidVector[0][i]          = tmpx / mag;
        EllipsoidVector[0][i + nq]     = tmpy / mag;
        EllipsoidVector[0][i + 2 * nq] = tmpz / mag;

        // Phi - direction
        tmpx = -1.0 * Rada * sin_theta * sin_varphi;
        tmpy = Radb * sin_theta * cos_varphi;
        tmpz = 0.0;

        mag = sqrt(tmpx * tmpx + tmpy * tmpy + tmpz * tmpz);

        EllipsoidVector[1][i]          = tmpx / mag;
        EllipsoidVector[1][i + nq]     = tmpy / mag;
        EllipsoidVector[1][i + 2 * nq] = tmpz / mag;

        // r - direction
        tmpx = Radb * Radc * sin_theta * sin_theta * cos_varphi;
        tmpy = Rada * Radc * sin_theta * sin_theta * sin_varphi;
        tmpz = Rada * Radb * sin_theta * cos_theta;

        mag = sqrt(tmpx * tmpx + tmpy * tmpy + tmpz * tmpz);

        EllipsoidVector[2][i]          = tmpx / mag;
        EllipsoidVector[2][i + nq]     = tmpy / mag;
        EllipsoidVector[2][i + 2 * nq] = tmpz / mag;
    }
}

void MMFSystem::CopyBoundaryTrace(const Array<OneD, const NekDouble> &Fwd,
                                  Array<OneD, NekDouble> &Bwd,
                                  const BoundaryCopyType BDCopyType,
                                  const int var, const std::string BDtype)
{
    int id1, id2, npts, nptselem, cnt = 0;
    Array<OneD, NekDouble> Dirichlet, x0, x1, x2;

    // loop over Boundary Regions
    for (int n = 0; n < m_fields[var]->GetBndConditions().size(); ++n)
    {
        nptselem = m_fields[var]->GetBndCondExpansions()[n]->GetNpoints();

        Dirichlet = Array<OneD, NekDouble>(nptselem);
        x0        = Array<OneD, NekDouble>(nptselem);
        x1        = Array<OneD, NekDouble>(nptselem);
        x2        = Array<OneD, NekDouble>(nptselem);

        if (BDCopyType == eDirichlet)
        {
            m_fields[var]->GetBndCondExpansions()[n]->GetCoords(x0, x1, x2);
            LibUtilities::EquationSharedPtr ifunc =
                m_session->GetFunction("BoundaryConditions", 0);
            ifunc->Evaluate(x0, x1, x2, 0.0, Dirichlet);
        }

        for (int e = 0;
             e < m_fields[var]->GetBndCondExpansions()[n]->GetExpSize(); ++e)
        {
            npts = m_fields[var]
                       ->GetBndCondExpansions()[n]
                       ->GetExp(e)
                       ->GetNumPoints(0);
            id1 = m_fields[var]->GetBndCondExpansions()[n]->GetPhys_Offset(e);
            id2 = m_fields[var]->GetTrace()->GetPhys_Offset(
                m_fields[var]->GetTraceMap()->GetBndCondIDToGlobalTraceID(cnt +
                                                                          e));

            if (m_fields[var]->GetBndConditions()[n]->GetUserDefined() ==
                    BDtype ||
                BDtype == "NoUserDefined")
            {
                switch (BDCopyType)
                {
                    case eDirichlet:
                    {
                        Vmath::Vcopy(npts, &Dirichlet[id1], 1, &Bwd[id2], 1);
                    }
                    break;

                    case eFwdEQBwd:
                    {
                        Vmath::Vcopy(npts, &Fwd[id2], 1, &Bwd[id2], 1);
                    }
                    break;

                    case eFwdEQNegBwd:
                    {
                        Vmath::Vcopy(npts, &Fwd[id2], 1, &Bwd[id2], 1);
                        Vmath::Neg(npts, &Bwd[id2], 1);
                    }
                    break;

                    default:
                        break;
                }
            }
        }

        cnt += m_fields[var]->GetBndCondExpansions()[n]->GetExpSize();
    }
}

// void MMFSystem::CopyBoundaryTrace(
//     const Array<OneD, const NekDouble> &Fwd, Array<OneD, NekDouble> &Bwd,
//     const BoundaryCopyType BDCopyType, const int var,
//     const std::string BDtype)
// {
//     int id1, id2, npts, nptselem, cnt = 0;
//     Array<OneD, NekDouble> Dirichlet, x0, x1, x2;

//     // loop over Boundary Regions
//     int max_bcRegion = m_fields[var]->GetBndConditions().size();
//     for (int bcRegion = 0; bcRegion < max_bcRegion; ++bcRegion)
//     {
//         nptselem =
//             m_fields[var]->GetBndCondExpansions()[bcRegion]->GetNpoints();

//         Dirichlet = Array<OneD, NekDouble>(nptselem);
//         x0        = Array<OneD, NekDouble>(nptselem);
//         x1        = Array<OneD, NekDouble>(nptselem);
//         x2        = Array<OneD, NekDouble>(nptselem);

//         if (BDCopyType == eDirichlet)
//         {
//             m_fields[var]->GetBndCondExpansions()[bcRegion]->GetCoords(x0,
//             x1,
//                                                                        x2);
//             LibUtilities::EquationSharedPtr ifunc =
//                 m_session->GetFunction("BOUNDARYCONDITIONS", 0);
//             ifunc->Evaluate(x0, x1, x2, 0.0, Dirichlet);
//         }

//         int max_e =
//             m_fields[var]->GetBndCondExpansions()[bcRegion]->GetExpSize();
//         for (int e = 0; e < max_e; ++e)
//         {
//             npts = m_fields[var]
//                        ->GetBndCondExpansions()[n]
//                        ->GetExp(e)
//                        ->GetNumPoints(0);
//             id1 =
//             m_fields[var]->GetBndCondExpansions()[n]->GetPhys_Offset(e); id2
//             = m_fields[var]->GetTrace()->GetPhys_Offset(
//                 m_fields[var]->GetTraceMap()->GetBndCondIDToGlobalTraceID(cnt
//                 +
//                                                                           e));

//             if (m_fields[var]->GetBndConditions()[n]->GetUserDefined() ==
//                     BDtype ||
//                 BDtype == "NoUserDefined")
//             {
//                 switch (BDCopyType)
//                 {
//                     case eDirichlet:
//                     {
//                         Vmath::Vcopy(npts, &Dirichlet[id1], 1, &Bwd[id2], 1);
//                     }
//                     break;

//                     case eFwdEQBwd:
//                     {
//                         Vmath::Vcopy(npts, &Fwd[id2], 1, &Bwd[id2], 1);
//                     }
//                     break;

//                     case eFwdEQNegBwd:
//                     {
//                         Vmath::Vcopy(npts, &Fwd[id2], 1, &Bwd[id2], 1);
//                         Vmath::Neg(npts, &Bwd[id2], 1);
//                     }
//                     break;

//                     default:
//                         break;
//                 }
//             }
//         }
//         cnt += m_fields[var]->GetBndCondExpansions()[bcRegion]->GetExpSize();
//     }
// }

// inarray = moving frames
// outarray = connection forms (m_connectionform)
// void MMFSystem::GetConnectionFormVector(
//     const Array<OneD, const Array<OneD, NekDouble>> &vector,
//     Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &Connectionform)
// {
//     int nq = vector[0].size();

//     Array<OneD, Array<OneD, NekDouble>> tmp(m_mfdim);
//     for (int i = 0; i < m_mfdim; ++i)
//     {
//         tmp[i] = Array<OneD, NekDouble>(3 * nq, 0.0);
//     }

//     for (int i = 0; i < m_mfdim; ++i)
//     {
//         Vmath::Vcopy(nq, &vector[i][0], 1, &tmp[0][i * nq], 1);
//         Vmath::Vcopy(nq, &m_movingframes[2][i * nq], 1, &tmp[2][i * nq], 1);
//     }

//     // e^2 = e^3 \times e^1
//     // e2x = e3y * e1z - e3z * e1y
//     // e2y = e3z * e1x - e3x * e1z
//     // e2z = e3x * e1y - e3x * e1y

//     // outarray[0] = v1[1]*v2[2] - v1[2]*v2[1];
//     // outarray[1] = v1[2]*v2[0] - v1[0]*v2[2];
//     // outarray[2] = v1[0]*v2[1] - v1[1]*v2[0];
//     NekDouble e1x, e1y, e1z, e3x, e3y, e3z;
//     for (int j = 0; j < nq; j++)
//     {
//         e1x = tmp[0][j];
//         e1y = tmp[0][j + nq];
//         e1z = tmp[0][j + 2 * nq];

//         e3x = tmp[2][j];
//         e3y = tmp[2][j + nq];
//         e3z = tmp[2][j + 2 * nq];

//         tmp[1][j]          = e3y * e1z - e3z * e1y;
//         tmp[1][j + nq]     = e3z * e1x - e3x * e1z;
//         tmp[1][j + 2 * nq] = e3x * e1y - e3x * e1y;
//     }

//     GetConnectionForm(tmp, Connectionform);
// }

// void MMFSystem::GetConnectionCurvature(
//     const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
//     Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &Connectionform,
//     Array<OneD, Array<OneD, Array<OneD, Array<OneD, NekDouble>>>>
//         &Curvatureform)
// {
//     GetConnectionForm(movingframes, Connectionform);
//     GetCurvatureForm(movingframes, Connectionform, Curvatureform);
// }

// void MMFSystem::Compute2DConnection1form(
//     const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
//     Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &Connection1form)
// {
//     ComputeConnection1form(movingframes, Connection1form);
// }

void MMFSystem::Compute2DConnectionCurvature(
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &Connectionform,
    Array<OneD, Array<OneD, NekDouble>> &Curvatureform)
{
    Compute2DConnection1form(movingframes, Connectionform);
    Compute2DCurvatureForm(movingframes, Connectionform, Curvatureform);
}

// w_{21} <e^1> = Connectionform[0][0];
// w_{21} <e^2> = Connectionform[0][1];
// w_{21} <e^3> = Connectionform[0][2];

// w_{31} <e^1> = Connectionform[1][0];
// w_{31} <e^2> = Connectionform[1][1];
// w_{31} <e^3> = Connectionform[1][2];

// w_{32} <e^1> = Connectionform[2][0];
// w_{32} <e^2> = Connectionform[2][1];
// w_{32} <e^3> = Connectionform[2][2];

// void MMFSystem::ComputeConnection1form(
//     const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
//     Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &Connectionform,
//     const int Verbose)

void MMFSystem::Compute2DConnection1form(
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &Connectionform)
{
    int nq = m_fields[0]->GetNpoints();

    Connectionform = Array<OneD, Array<OneD, Array<OneD, NekDouble>>>(m_mfdim);
    for (int i = 0; i < m_mfdim; i++)
    {
        Connectionform[i] = Array<OneD, Array<OneD, NekDouble>>(m_mfdim);
        for (int j = 0; j < m_mfdim; j++)
        {
            Connectionform[i][j] = Array<OneD, NekDouble>(nq, 0.0);
        }
    }

    // int Ntot = 0;
    // NekDouble MFmag;
    // for (int i = 0; i < nq; ++i)
    // {
    //     MFmag = movingframes[0][i] * movingframes[0][i] +
    //             movingframes[0][i + nq] * movingframes[0][i + nq] +
    //             movingframes[0][i + 2 * nq] * movingframes[0][i + 2 * nq];

    //     MFmag = sqrt(MFmag);
    //     if (MFmag > 0.1)
    //     {
    //         Ntot++;
    //     }
    // }

    // std::cout << " Compute Connection 1form " << std::endl;
    // CheckMovingFrames(movingframes);

    // outarray is 3 x 3 x 3 x nq (MF x connection.row x connection.col
    // x numpts) MFJacobian[i][j][k] = d e^i_j / d x_k
    Array<OneD, Array<OneD, Array<OneD, Array<OneD, NekDouble>>>> MFJacobian;
    ComputeMFJacobian(movingframes, MFJacobian);

    for (int i = 0; i < m_mfdim; i++)
    {
        Connectionform[i] = Array<OneD, Array<OneD, NekDouble>>(m_mfdim);
        for (int j = 0; j < m_mfdim; j++)
        {
            Connectionform[i][j] = Array<OneD, NekDouble>(nq, 0.0);
        }
    }

    // w_{ij} (e_k) = f_j J(f_i) f_k

    // w_{21} <e^1> = f_1 J(f_2) f_1,
    // w_(21) <e^2> = f_1 J(f_2) f_2,
    // w_{21} <e^3> = f_1 J(f_2) f_3,
    // w_{31} <e^1> = f_1 J(f_3) f_1,
    // w_{31} <e^2> = f_1 J(f_3) f_2,
    // w_{31} <e^3> = f_1 J(f_3) f_3,
    // w_{32} <e^1> = f_2 J(f_3) f_1,
    // w_{32} <e^2> = f_2 J(f_3) f_2,
    // w_{32} <e^2> = f_2 J(f_3) f_3
    for (int i = 0; i < nq; i++)
    {
        for (int k = 0; k < m_spacedim; k++)
        {
            // w_{21} <e^k> = f_1 J(f_2) f_k,
            Connectionform[0][k][i] =
                ComputeVMV(1, 0, k, i, movingframes, MFJacobian);

            // w_{31} <e^k> = f_1 J(f_3) f_k,
            Connectionform[1][k][i] =
                ComputeVMV(2, 0, k, i, movingframes, MFJacobian);

            // w_{32} <e^k> = f_2 J(f_3) f_k,
            Connectionform[2][k][i] =
                ComputeVMV(2, 1, k, i, movingframes, MFJacobian);
        }
    }

    // // w^2_1 <e^1> = f_2 J(f_1) f_1,
    // // w^2_1 <e^2> = f_2 J(f_1) f_2,
    // // w^2_1 <e^3> = f_2 J(f_1) f_3,
    // // w^3_1 <e^1> = f_3 J(f_1) f_1,
    // // w^3_1 <e^2> = f_3 J(f_1) f_2,
    // // w^3_1 <e^3> = f_3 J(f_1) f_3,
    // // w^3_2 <e^1> = f_3 J(f_2) f_1,
    // // w^3_2 <e^2> = f_3 J(f_2) f_2,
    // // w^3_2 <e^2> = f_3 J(f_2) f_3
    // for (int k = 0; k < nq; k++)
    // {
    //     // w^2_1 <e^1> = f_2 J(f_1) f_1
    //     Connectionform[0][0][k] =
    //         ComputeVMV(0, 1, 0, k, movingframes, MFJacobian);

    //     // w^2_1 <e^2> = f_2 J(f_1) f_2
    //     Connectionform[0][1][k] =
    //         ComputeVMV(0, 1, 1, k, movingframes, MFJacobian);

    //     // w^2_1 <e^3> = f_2 J(f_1) f_3
    //     Connectionform[0][2][k] =
    //         ComputeVMV(0, 1, 2, k, movingframes, MFJacobian);

    //     // w^3_1 <e^1> = f_3 J(f_1) f_1
    //     Connectionform[1][0][k] =
    //         ComputeVMV(0, 2, 0, k, movingframes, MFJacobian);

    //     // w^3_1 <e^2> = f_3 J(f_1) f_2
    //     Connectionform[1][1][k] =
    //         ComputeVMV(0, 2, 1, k, movingframes, MFJacobian);

    //     // w^3_1 <e^3> = f_3 J(f_1) f_3
    //     Connectionform[1][2][k] =
    //         ComputeVMV(0, 2, 2, k, movingframes, MFJacobian);

    //     // w^3_2 <e^1> = f_3 J(f_2) f_1
    //     Connectionform[2][0][k] =
    //         ComputeVMV(1, 2, 0, k, movingframes, MFJacobian);

    //     // w^3_2 <e^2> = f_3 J(f_2) f_2
    //     Connectionform[2][1][k] =
    //         ComputeVMV(1, 2, 1, k, movingframes, MFJacobian);

    //     // w^3_2 <e^3> = f_3 J(f_2) f_3
    //     Connectionform[2][2][k] =
    //         ComputeVMV(1, 2, 2, k, movingframes, MFJacobian);
    // }

    // Array<OneD, Array<OneD, int>> EWIndex;
    // m_fields[0]->GridIndexElementWise(EWIndex);

    // int Nelemtj = EWIndex.size();
    // int nptsj   = EWIndex[0].size();

    // if (Verbose)
    // {
    //     std::cout << "Connection form for " << Ntot / nptsj << " / " <<
    //     Nelemtj
    //               << " elements" << std::endl;

    //     std::cout << "(w211, w212, w213) = ( "
    //               << RootMeanSquare(Connectionform[0][0], m_MMFActivation) <<
    //               " , "
    //               << RootMeanSquare(Connectionform[0][1], m_MMFActivation) <<
    //               " , "
    //               << RootMeanSquare(Connectionform[0][2], m_MMFActivation) <<
    //               " ) "
    //               << std::endl;
    //     std::cout << "(w311, w312, w313) = ( "
    //               << RootMeanSquare(Connectionform[1][0], m_MMFActivation) <<
    //               " , "
    //               << RootMeanSquare(Connectionform[1][1], m_MMFActivation) <<
    //               " , "
    //               << RootMeanSquare(Connectionform[1][2], m_MMFActivation) <<
    //               " ) "
    //               << std::endl;
    //     std::cout << "(w321, w322, w323) = ( "
    //               << RootMeanSquare(Connectionform[2][0], m_MMFActivation) <<
    //               " , "
    //               << RootMeanSquare(Connectionform[2][1], m_MMFActivation) <<
    //               " , "
    //               << RootMeanSquare(Connectionform[2][2], m_MMFActivation) <<
    //               " ) "
    //               << std::endl
    //               << std::endl;
    // }
}

void MMFSystem::GetCurvatureForm(
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    const Array<OneD, const Array<OneD, Array<OneD, NekDouble>>>
        &Connectionform,
    Array<OneD, Array<OneD, Array<OneD, Array<OneD, NekDouble>>>>
        &Curvatureform)
{
    int nq = m_fields[0]->GetNpoints();

   int Ntot = 0;
    NekDouble MFmag;
    for (int i = 0; i < nq; ++i)
    {
        MFmag = movingframes[0][i] * movingframes[0][i] +
                movingframes[0][i + nq] * movingframes[0][i + nq] +
                movingframes[0][i + 2 * nq] * movingframes[0][i + 2 * nq];
        MFmag = sqrt(MFmag);
        if (MFmag > 0.1)
        {
            Ntot++;
        }
    }

    // w^2_1 <e^1> = Connectionform[0][0];
    // w^2_1 <e^2> = Connectionform[0][1];
    // w^2_1 <e^3> = Connectionform[0][2];

    // w^3_1 <e^1> = Connectionform[1][0];
    // w^3_1 <e^2> = Connectionform[1][1];
    // w^3_1 <e^3> = Connectionform[1][2];

    // w^3_2 <e^1> = Connectionform[2][0];
    // w^3_2 <e^2> = Connectionform[2][1];
    // w^3_2 <e^3> = Connectionform[2][2];

    // R^l_{kij} = \Omega^l_k (e^i, e^j) = d \omega^l_k (e^i) (e^j) +
    // w^l_m \wedge w^m_k (e^i, e^j)
    Array<OneD, NekDouble> wij(nq);
    Array<OneD, NekDouble> dwij(nq);
    for (int i = 0; i < m_mfdim; i++)
    {
        for (int j = 0; j < m_mfdim; j++)
        {
            Vmath::Vcopy(nq, &Connectionform[i][j][0], 1, &wij[0], 1);
            for (int k = 0; k < m_mfdim; ++k)
            {
                MMFDirectionalDeriv(movingframes[k], wij, dwij);
                Vmath::Vadd(nq, &dwij[0], 1, &Curvatureform[i][j][k][0], 1,
                            &Curvatureform[i][j][k][0], 1);
            }
        }
    }

    // Add torsion
    NekDouble w211, w212, w213, w311, w312, w313, w321, w322, w323;
    for (int k = 0; k < nq; ++k)
    {
        w211 = Connectionform[0][0][k];
        w212 = Connectionform[0][1][k];
        w213 = Connectionform[0][2][k];

        w311 = Connectionform[1][0][k];
        w312 = Connectionform[1][1][k];
        w313 = Connectionform[1][2][k];

        w321 = Connectionform[2][0][k];
        w322 = Connectionform[2][1][k];
        w323 = Connectionform[2][2][k];

        // Omega^2_1 (e^1, e^1) = dw^2_1 (e^1,e^1) - w^3_2 \wedge w^3_1
        // (e^1,e^1)
        Curvatureform[0][0][0][k] += 0.0;

        // Omega^2_1 (e^1, e^2) = dw^2_1 (e^1,e^2) - w^3_2 \wedge w^3_1
        // (e^1,e^2)
        Curvatureform[0][0][1][k] += -1.0 * (w322 * w312 - w322 * w311);

        // Omega^2_1 (e^1, e^3) = dw^2_1 (e^1,e^3) - w^3_2 \wedge w^3_1
        // (e^1,e^3)
        Curvatureform[0][0][2][k] += -1.0 * (w321 * w313 - w323 * w311);

        // Omega^2_1 (e^2, e^1) = dw^2_1 (e^2,e^1) - w^3_2 \wedge w^3_1
        // (e^2,e^1)
        Curvatureform[0][1][0][k] = -1.0 * Curvatureform[0][0][1][k];

        // Omega^2_1 (e^2, e^2) = dw^2_1 (e^2,e^2) - w^3_2 \wedge w^3_1
        // (e^2,e^2)
        Curvatureform[0][1][1][k] += 0.0;

        // Omega^2_1 (e^2, e^3) = dw^2_1 (e^2,e^3) - w^3_2 \wedge w^3_1
        // (e^2,e^3)
        Curvatureform[0][1][2][k] += -1.0 * (w322 * w313 - w323 * w312);

        // Omega^2_1 (e^3, e^1) = dw^2_1 (e^3,e^1) - w^3_2 \wedge w^3_1
        // (e^3,e^1)
        Curvatureform[0][2][0][k] = -1.0 * Curvatureform[0][0][2][k];

        // Omega^2_1 (e^3, e^2) = dw^2_1 (e^3,e^2) - w^3_2 \wedge w^3_1
        // (e^3,e^2)
        Curvatureform[0][2][1][k] = -1.0 * Curvatureform[0][1][2][k];

        // Omega^2_1 (e^3, e^3) = dw^2_1 (e^3,e^3) - w^3_2 \wedge w^3_1
        // (e^3,e^3)
        Curvatureform[0][2][2][k] += 0.0;

        // Omega^3_1 (e^1, e^1) = dw^3_1 (e^1,e^1) + w^3_2 \wedge w^2_1
        // (e^1,e^1)
        Curvatureform[1][0][0][k] += 0.0;

        // Omega^3_1 (e^1, e^2) = dw^3_1 (e^1,e^2) + w^3_2 \wedge w^2_1
        // (e^1,e^2)
        Curvatureform[1][0][1][k] += w321 * w212 - w322 * w211;

        // Omega^3_1 (e^1, e^3) = dw^3_1 (e^1,e^3) + w^3_2 \wedge w^2_1
        // (e^1,e^3)
        Curvatureform[1][0][2][k] += w321 * w213 - w323 * w211;

        // Omega^3_1 (e^2, e^1) = dw^3_1 (e^2,e^1) + w^3_2 \wedge w^2_1
        // (e^2,e^1)
        Curvatureform[1][1][0][k] = -1.0 * Curvatureform[1][0][1][k];

        // Omega^3_1 (e^2, e^2) = dw^3_1 (e^2,e^2) + w^3_2 \wedge w^2_1
        // (e^2,e^2)
        Curvatureform[1][1][1][k] += 0.0;

        // Omega^3_1 (e^2, e^3) = dw^3_1 (e^2,e^3) + w^3_2 \wedge w^2_1
        // (e^2,e^3)
        Curvatureform[1][1][2][k] += w322 * w213 - w323 * w212;

        // Omega^3_1 (e^3, e^1) = dw^3_1 (e^3,e^1) + w^3_2 \wedge w^2_1
        // (e^3,e^1)
        Curvatureform[1][2][0][k] = -1.0 * Curvatureform[1][0][2][k];

        // Omega^3_1 (e^3, e^2) = dw^3_1 (e^3,e^2) + w^3_2 \wedge w^2_1
        // (e^3,e^2)
        Curvatureform[1][2][1][k] = -1.0 * Curvatureform[1][1][2][k];

        // Omega^3_1 (e^3, e^3) = dw^3_1 (e^3,e^3) + w^3_2 \wedge w^2_1
        // (e^3,e^3)
        Curvatureform[1][2][2][k] += 0.0;

        // Omega^3_2 (e^1, e^1) = dw^3_2 (e^1,e^1) - w^3_1 \wedge w^2_1
        // (e^1,e^1)
        Curvatureform[2][0][0][k] += 0.0;

        // Omega^3_2 (e^1, e^2) = dw^3_2 (e^1,e^2) - w^3_1 \wedge w^2_1
        // (e^1,e^2)
        Curvatureform[2][0][1][k] += -1.0 * (w311 * w212 - w312 * w211);

        // Omega^3_2 (e^1, e^3) = dw^3_2 (e^1,e^3) - w^3_1 \wedge w^2_1
        // (e^1,e^3)
        Curvatureform[2][0][2][k] += -1.0 * (w311 * w213 - w313 * w211);

        // Omega^3_2 (e^2, e^1) = dw^3_2 (e^2,e^1) - w^3_1 \wedge w^2_1
        // (e^2,e^1)
        Curvatureform[2][1][0][k] = -1.0 * Curvatureform[2][0][1][k];

        // Omega^3_2 (e^2, e^2) = dw^3_2 (e^2,e^2) - w^3_1 \wedge w^2_1
        // (e^2,e^2)
        Curvatureform[2][1][1][k] += 0.0;

        // Omega^3_2 (e^2, e^3) = dw^3_2 (e^2,e^3) - w^3_1 \wedge w^2_1
        // (e^2,e^3)
        Curvatureform[2][1][2][k] += -1.0 * (w312 * w213 - w313 * w212);

        // Omega^3_2 (e^3, e^1) = dw^3_1 (e^3,e^1) - w^3_1 \wedge w^2_1
        // (e^3,e^1)
        Curvatureform[2][2][0][k] = -1.0 * Curvatureform[2][0][2][k];

        // Omega^3_2 (e^3, e^2) = dw^3_1 (e^3,e^2) - w^3_1 \wedge w^2_1
        // (e^3,e^2)
        Curvatureform[2][2][1][k] = -1.0 * Curvatureform[2][1][2][k];

        // Omega^3_2 (e^3, e^3) = dw^3_1 (e^3,e^3) - w^3_1 \wedge w^2_1
        // (e^3,e^3)
        Curvatureform[2][2][2][k] += 0.0;
    }

    std::cout << "(R^2_{111}, R^2_{112}, R^2_{113}) = ( "
              << RootMeanSquare(Curvatureform[0][0][0], Ntot) << " , "
              << RootMeanSquare(Curvatureform[0][0][1], Ntot) << " , "
              << RootMeanSquare(Curvatureform[0][0][2], Ntot) << " ) "
              << std::endl;

    std::cout << "(R^2_{121}, R^2_{122}, R^2_{123}) = ( "
              << RootMeanSquare(Curvatureform[0][1][0], Ntot) << " , "
              << RootMeanSquare(Curvatureform[0][1][1], Ntot) << " , "
              << RootMeanSquare(Curvatureform[0][1][2], Ntot) << " ) "
              << std::endl;

    std::cout << "(R^2_{131}, R^2_{132}, R^2_{133}) = ( "
              << RootMeanSquare(Curvatureform[0][2][0], Ntot) << " , "
              << RootMeanSquare(Curvatureform[0][2][1], Ntot) << " , "
              << RootMeanSquare(Curvatureform[0][2][2], Ntot) << " ) "
              << std::endl;

    std::cout << "(R^3_{111}, R^3_{112}, R^3_{113}) = ( "
              << RootMeanSquare(Curvatureform[1][0][0], Ntot) << " , "
              << RootMeanSquare(Curvatureform[1][0][1], Ntot) << " , "
              << RootMeanSquare(Curvatureform[1][0][2], Ntot) << " ) "
              << std::endl;

    std::cout << "(R^3_{121}, R^3_{122}, R^3_{123}) = ( "
              << RootMeanSquare(Curvatureform[1][1][0], Ntot) << " , "
              << RootMeanSquare(Curvatureform[1][1][1], Ntot) << " , "
              << RootMeanSquare(Curvatureform[1][1][2], Ntot) << " ) "
              << std::endl;

    std::cout << "(R^3_{131}, R^3_{132}, R^3_{133}) = ( "
              << RootMeanSquare(Curvatureform[1][2][0], Ntot) << " , "
              << RootMeanSquare(Curvatureform[1][2][1], Ntot) << " , "
              << RootMeanSquare(Curvatureform[1][2][2], Ntot) << " ) "
              << std::endl;

    std::cout << "(R^3_{211}, R^3_{212}, R^3_{213}) = ( "
              << RootMeanSquare(Curvatureform[2][0][0], Ntot) << " , "
              << RootMeanSquare(Curvatureform[2][0][1], Ntot) << " , "
              << RootMeanSquare(Curvatureform[2][0][2], Ntot) << " ) "
              << std::endl;

    std::cout << "(R^3_{221}, R^3_{222}, R^3_{223}) = ( "
              << RootMeanSquare(Curvatureform[2][1][0], Ntot) << " , "
              << RootMeanSquare(Curvatureform[2][1][1], Ntot) << " , "
              << RootMeanSquare(Curvatureform[2][1][2], Ntot) << " ) "
              << std::endl;

    std::cout << "(R^3_{231}, R^3_{232}, R^3_{233}) = ( "
              << RootMeanSquare(Curvatureform[2][2][0], Ntot) << " , "
              << RootMeanSquare(Curvatureform[2][2][1], Ntot) << " , "
              << RootMeanSquare(Curvatureform[2][2][2], Ntot) << " ) "
              << std::endl
              << std::endl;

    NekDouble R11, R12, R22;
    // R11 = R^2_{121} + R^3_{131}
    R11 = RootMeanSquare(Curvatureform[0][1][0], Ntot) +
          RootMeanSquare(Curvatureform[1][2][0], Ntot);

    // R12 = R^2_{122} + R^3_{132}
    R12 = RootMeanSquare(Curvatureform[0][1][1], Ntot) +
          RootMeanSquare(Curvatureform[1][2][1], Ntot);

    // R22 = R^3_{232}
    R22 = RootMeanSquare(Curvatureform[2][2][1], Ntot);

    std::cout << "Ricci Tensor: R11 = " << R11 << ", R12 = " << R12
              << ", R22 = " << R22 << std::endl;
}

// Compute Riemannian curvature tensor on 2D surfaces
// R^2_{112} = Curvatureform[0]
// R^3_{112} = Curvatureform[1]
// R^3_{212} = Curvatureform[2]
void MMFSystem::Compute2DCurvatureForm(
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    const Array<OneD, const Array<OneD, Array<OneD, NekDouble>>>
        &Connectionform,
    Array<OneD, Array<OneD, NekDouble>> &Curvatureform)
{
    int nq = m_fields[0]->GetNpoints();

    Curvatureform = Array<OneD, Array<OneD, NekDouble>>(m_mfdim);
    for (int i = 0; i < m_mfdim; i++)
    {
        Curvatureform[i] = Array<OneD, NekDouble>(nq, 0.0);
    }

    // NekDouble MFmag;
    // int Ntot = 0;
    // for (int i = 0; i < nq; ++i)
    // {
    //     MFmag = movingframes[0][i] * movingframes[0][i] +
    //             movingframes[0][i + nq] * movingframes[0][i + nq] +
    //             movingframes[0][i + 2 * nq] * movingframes[0][i + 2 * nq];
    //     MFmag = sqrt(MFmag);
    //     if (MFmag > 0.1)
    //     {
    //         Ntot++;
    //     }
    // }

    // w^2_1 <e^1> = Connectionform[0][0];
    // w^2_1 <e^2> = Connectionform[0][1];

    // w^3_1 <e^1> = Connectionform[1][0];
    // w^3_1 <e^2> = Connectionform[1][1];

    // w^3_2 <e^1> = Connectionform[2][0];
    // w^3_2 <e^2> = Connectionform[2][1];

    // R^l_{kij} = \Omega^l_k (e^i, e^j) = d \omega^l_k (e^i) (e^j) +
    // w^l_m \wedge w^m_k (e^i, e^j)
    Array<OneD, NekDouble> wij0(nq);
    Array<OneD, NekDouble> wij1(nq);

    Array<OneD, NekDouble> dwij01(nq);
    Array<OneD, NekDouble> dwij10(nq);

    for (int i = 0; i < m_mfdim; i++)
    {
        Vmath::Vcopy(nq, &Connectionform[i][0][0], 1, &wij0[0], 1);
        Vmath::Vcopy(nq, &Connectionform[i][1][0], 1, &wij1[0], 1);

        MMFDirectionalDeriv(movingframes[1], wij0, dwij01);
        MMFDirectionalDeriv(movingframes[0], wij1, dwij10);

        Vmath::Neg(nq, dwij10, 1);
        Vmath::Vadd(nq, &dwij01[0], 1, &dwij10[0], 1, &Curvatureform[i][0], 1);
    }
}

// fj^T * Jacobian^i * fk
NekDouble MMFSystem::ComputeVMV(
    const int indexi, const int indexj, const int indexk, const int knode,
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    const Array<OneD, const Array<OneD, Array<OneD, Array<OneD, NekDouble>>>>
        &MFJacobian)
{
    int nq = m_fields[0]->GetNpoints();

    // Create fjvec
    Array<OneD, NekDouble> fjvec(m_mfdim);
    fjvec[0] = movingframes[indexj][knode];
    fjvec[1] = movingframes[indexj][nq + knode];
    fjvec[2] = movingframes[indexj][2 * nq + knode];

    // Create fkvec
    Array<OneD, NekDouble> fkvec(m_mfdim);
    fkvec[0] = movingframes[indexk][knode];
    fkvec[1] = movingframes[indexk][nq + knode];
    fkvec[2] = movingframes[indexk][2 * nq + knode];

    // Create Jacobian Matrix
    Array<OneD, Array<OneD, NekDouble>> Jacobian(m_mfdim);
    for (int i = 0; i < m_mfdim; i++)
    {
        Jacobian[i] = Array<OneD, NekDouble>(m_mfdim);
    }

    for (int i = 0; i < m_mfdim; i++)
    {
        for (int j = 0; j < m_mfdim; j++)
        {
            Jacobian[i][j] = MFJacobian[indexi][i][j][knode];
        }
    }

    // fj^T * Jacobian^i * fk
    Array<OneD, NekDouble> locsum(3, 0.0);
    NekDouble rval = 0.0;
    for (int i = 0; i < m_mfdim; i++)
    {
        for (int j = 0; j < m_mfdim; j++)
        {
            locsum[i] = locsum[i] + fjvec[j] * Jacobian[j][i];
        }

        rval = rval + locsum[i] * fkvec[i];
    }

    return rval;
}

void MMFSystem::CheckJacobianError()
{
    int nq = m_fields[0]->GetNpoints();

    Array<OneD, NekDouble> Jacobian(nq);
    m_fields[0]->GetJacobian(Jacobian);

    Array<OneD, NekDouble> JacobianAvg(nq);
    JacobianAvg = ComputeJacobianAvg(Jacobian);

    Vmath::Vsub(nq, Jacobian, 1, JacobianAvg, 1, JacobianAvg, 1);
    Vmath::Vabs(nq, JacobianAvg, 1, JacobianAvg, 1);

    std::cout << "JacErr = " << RootMeanSquare(JacobianAvg)
              << ", max = " << Vmath::Vamax(nq, JacobianAvg, 1) << std::endl;

    Array<OneD, NekDouble> GradientJac(m_spacedim * nq);

    // m_fields[0]->GetTangentVectors(TangentVectors);

    GradientJac = ComputeCovGrad(Jacobian, m_movingframes);

    Array<OneD, NekDouble> GradientJacMag(nq, 0.0);
    Array<OneD, NekDouble> GradJacCdotVel(nq, 0.0);
    for (int i = 0; i < m_spacedim; ++i)
    {
        Vmath::Vvtvp(nq, &GradientJac[i * nq], 1, &GradientJac[i * nq], 1,
                     &GradientJacMag[0], 1, &GradientJacMag[0], 1);
    }
    Vmath::Vsqrt(nq, GradientJacMag, 1, GradientJacMag, 1);

    std::cout << "Gradient of Jac = " << RootMeanSquare(GradientJacMag)
              << ", max = " << Vmath::Vamax(nq, GradientJacMag, 1) << std::endl;

    PlotJacobian(JacobianAvg, GradientJacMag);

    std::cout << std::endl;
}

Array<OneD, NekDouble> MMFSystem::ComputeJacobianAvg(
    const Array<OneD, const NekDouble> &Jacobian)
{
    int nq = m_fields[0]->GetNpoints();

    Array<OneD, NekDouble> outarray(nq, 0.0);

    Array<OneD, NekDouble> x(nq);
    Array<OneD, NekDouble> y(nq);
    Array<OneD, NekDouble> z(nq);

    m_fields[0]->GetCoords(x, y, z);

    int ind;
    NekDouble sum;

    for (int i = 0; i < m_fields[0]->GetExpSize(); ++i)
    {
        sum = 0.0;
        for (int j = 0; j < m_fields[0]->GetTotPoints(i); ++j)
        {
            ind = m_fields[0]->GetPhys_Offset(i) + j;
            sum += Jacobian[ind];
        }
        sum = sum / m_fields[0]->GetTotPoints(i);

        for (int j = 0; j < m_fields[0]->GetTotPoints(i); ++j)
        {
            ind           = m_fields[0]->GetPhys_Offset(i) + j;
            outarray[ind] = sum;
        }
    }

    return outarray;
}

// For an vector of 3 x 1, compute Jacobian of 3 x 3 matrix
void MMFSystem::ComputeMFJacobian(
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    Array<OneD, Array<OneD, Array<OneD, Array<OneD, NekDouble>>>> &outarray)
{
    int nq = m_fields[0]->GetNpoints();

    // outarray is 3 x 3 x 3 x nq (MF x connection.row x connection.col
    // x numpts)
    outarray =
        Array<OneD, Array<OneD, Array<OneD, Array<OneD, NekDouble>>>>(m_mfdim);
    for (int i = 0; i < m_mfdim; i++)
    {
        outarray[i] =
            Array<OneD, Array<OneD, Array<OneD, NekDouble>>>(m_spacedim);
        for (int j = 0; j < m_spacedim; j++)
        {
            outarray[i][j] = Array<OneD, Array<OneD, NekDouble>>(m_spacedim);
            for (int k = 0; k < m_spacedim; k++)
            {
                outarray[i][j][k] = Array<OneD, NekDouble>(nq);
            }
        }
    }

    Array<OneD, NekDouble> dfvec(nq);
    Array<OneD, NekDouble> fvec(nq);
    for (int i = 0; i < m_mfdim; i++)
    {
        for (int j = 0; j < m_spacedim; j++)
        {
            Vmath::Vcopy(nq, &movingframes[i][j * nq], 1, &fvec[0], 1);
            for (int k = 0; k < m_spacedim; k++)
            {
                // Compute d f_{ij} / d x_k
                m_fields[0]->PhysDeriv(k, fvec, dfvec);
                Vmath::Vcopy(nq, &dfvec[0], 1, &outarray[i][j][k][0], 1);
            }
        }
    }
}

NekDouble MMFSystem::AvgInt(const Array<OneD, const NekDouble> &inarray)
{
    int nq = m_fields[0]->GetNpoints();
    Array<OneD, NekDouble> Ones(nq, 1.0);

    if (inarray.size() != nq)
    {
        ASSERTL0(false, "AvgInt Error: Vector size is not correct");
    }

    NekDouble jac = m_fields[0]->PhysIntegral(Ones);

   return (m_fields[0]->PhysIntegral(inarray)) / jac;
}

NekDouble MMFSystem::AvgAbsInt(const Array<OneD, const NekDouble> &inarray)
{
    int nq = m_fields[0]->GetNpoints();
    Array<OneD, NekDouble> Ones(nq, 1.0);
    Array<OneD, NekDouble> tmp(nq);

    if (inarray.size() != nq)
    {
        ASSERTL0(false, "AvgAbsInt Error: Vector size is not correct");
    }

    NekDouble jac = m_fields[0]->PhysIntegral(Ones);

    Vmath::Vabs(nq, inarray, 1, tmp, 1);
    return (m_fields[0]->PhysIntegral(tmp)) / jac;
}

NekDouble MMFSystem::AbsIntegral(const Array<OneD, const NekDouble> &inarray)
{
    int nq = m_fields[0]->GetNpoints();
    Array<OneD, NekDouble> tmp(nq);

    if (inarray.size() != nq)
    {
        ASSERTL0(false, "AbsIntegral Error: Vector size is not correct");
    }

    Vmath::Vabs(nq, inarray, 1, tmp, 1);
    return m_fields[0]->PhysIntegral(tmp);
}

NekDouble MMFSystem::FindAbsMaximum(const Array<OneD, const NekDouble> &inarray,
                                    const Array<OneD, const int> &Activated)
{
    int nqtot = inarray.size();
    // int cn    = 0;

    NekDouble reval = 0.0;
    for (int i = 0; i < nqtot; ++i)
    {
        if (Activated[i] > 0)
        {
            if (fabs(inarray[i]) > reval)
            {
                reval = fabs(inarray[i]);
            }
     //       cn++;
        }
    }

    return reval;
}

NekDouble MMFSystem::FindAbsMaximumVector(
    const Array<OneD, const NekDouble> &inarray,
    const Array<OneD, const int> &Activated)
{
    int nq = inarray.size() / m_spacedim;

    NekDouble vmag, reval = 0.0;
    for (int i = 0; i < nq; ++i)
    {
        if (Activated[i] > 0)
        {
            vmag = sqrt(inarray[i] * inarray[i] +
                        inarray[i + nq] * inarray[i + nq] +
                        inarray[i + 2 * nq] * inarray[i + 2 * nq]);
            if (vmag > reval)
            {
                reval = vmag;
            }
        }
    }

    return reval;
}

NekDouble MMFSystem::RootMeanSquare(const Array<OneD, const NekDouble> &inarray,
                                    const Array<OneD, const int> &Activated)
{
    int nq = inarray.size();
    int cn = 0;

    NekDouble reval = 0.0;
    for (int i = 0; i < nq; ++i)
    {
        if (Activated[i] > 0)
        {
            reval += inarray[i] * inarray[i];
            cn++;
        }
    }
    reval = sqrt(reval / cn);
    return reval;
}

NekDouble MMFSystem::RootMeanSquareVector(
    const Array<OneD, const NekDouble> &inarray,
    const Array<OneD, const int> &Activated)
{
    int nq = inarray.size() / m_spacedim;
    int cn = 0;

    NekDouble reval = 0.0;
    for (int i = 0; i < nq; ++i)
    {
        if (Activated[i])
        {
            reval += inarray[i] * inarray[i] +
                     inarray[i + nq] * inarray[i + nq] +
                     inarray[i + 2 * nq] * inarray[i + 2 * nq];
            cn++;
        }
    }
    reval = sqrt(reval / cn);
    return reval;
}

NekDouble MMFSystem::Average(const Array<OneD, const NekDouble> &inarray)
{
    int nq = inarray.size();

    NekDouble reval = 0.0;
    for (int i = 0; i < nq; ++i)
    {
        reval += inarray[i];
    }

    reval = reval / nq;

    return reval;
}

NekDouble MMFSystem::RootMeanSquare(const Array<OneD, const NekDouble> &inarray,
                                    const int Ntot)
{
    int nq = inarray.size();

    NekDouble reval = 0.0;
    for (int i = 0; i < nq; ++i)
    {
        reval += inarray[i] * inarray[i];
    }

    if (Ntot)
    {
        reval = sqrt(reval / Ntot);
    }

    else
    {
        reval = sqrt(reval / nq);
    }

    return reval;
}

NekDouble MMFSystem::RootMeanSquare(const Array<OneD, const int> &inarray,
                                    const int Ntot)
{
    int nq = inarray.size();

    NekDouble reval = 0.0;
    for (int i = 0; i < nq; ++i)
    {
        reval += 1.0 * inarray[i] * inarray[i];
    }

    if (Ntot)
    {
        reval = sqrt(reval / Ntot);
    }

    else
    {
        reval = sqrt(reval / nq);
    }

    return reval;
}

NekDouble MMFSystem::VectorAvgMagnitude(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray)
{
    int nq = inarray[0].size();

    Array<OneD, NekDouble> tmp(nq, 0.0);
    for (int k = 0; k < m_spacedim; k++)
    {
        Vmath::Vvtvp(nq, &inarray[k][0], 1, &inarray[k][0], 1, &tmp[0], 1,
                     &tmp[0], 1);
    }
    Vmath::Vsqrt(nq, tmp, 1, tmp, 1);

    return RootMeanSquare(tmp);
}

void MMFSystem::BubbleSort(Array<OneD, NekDouble> &refarray,
                           Array<OneD, NekDouble> &sortarray)
{
    int nq = refarray.size();

    bool swapped = true;
    int j        = 0;
    NekDouble tmp;

    while (swapped)
    {
        swapped = false;
        j++;
        for (int i = 0; i < nq - j; i++)
        {
            if (refarray[i] > refarray[i + 1])
            {
                tmp             = refarray[i];
                refarray[i]     = refarray[i + 1];
                refarray[i + 1] = tmp;

                tmp              = sortarray[i];
                sortarray[i]     = sortarray[i + 1];
                sortarray[i + 1] = tmp;

                swapped = true;
            }
        }
    }
}

void MMFSystem::GramSchumitz(
    const Array<OneD, const Array<OneD, NekDouble>> &v1,
    const Array<OneD, const Array<OneD, NekDouble>> &v2,
    Array<OneD, Array<OneD, NekDouble>> &outarray, bool KeepTheMagnitude)
{
    int nq = v1[0].size();
    Array<OneD, NekDouble> tmp(nq, 0.0);
    Array<OneD, NekDouble> mag(nq, 0.0);

    for (int i = 0; i < m_spacedim; ++i)
    {
        // u2 = v2 - < u1 , v2 > ( u1 / < u1, u1 > )
        Vmath::Vvtvp(nq, &v1[i][0], 1, &v2[i][0], 1, &tmp[0], 1, &tmp[0], 1);
        Vmath::Vvtvp(nq, &v1[i][0], 1, &v1[i][0], 1, &mag[0], 1, &mag[0], 1);
    }
    Vmath::Vdiv(nq, &tmp[0], 1, &mag[0], 1, &tmp[0], 1);
    Vmath::Neg(nq, &tmp[0], 1);

    // outarray = Array<OneD, Array<OneD, NekDouble> >(m_spacedim);

    // u2 = v2 - < u1 , v2 > ( u1 / < u1, u1 > )
    for (int i = 0; i < m_spacedim; ++i)
    {
        outarray[i] = Array<OneD, NekDouble>(nq, 0.0);
        Vmath::Vvtvp(nq, &tmp[0], 1, &v1[i][0], 1, &v2[i][0], 1,
                     &outarray[i][0], 1);
    }

    if (KeepTheMagnitude)
    {
        Array<OneD, NekDouble> magorig(nq, 0.0);
        Array<OneD, NekDouble> magnew(nq, 0.0);

        for (int i = 0; i < m_spacedim; ++i)
        {
            Vmath::Vmul(nq, &v2[0][0], 1, &v2[0][0], 1, &magorig[0], 1);
            Vmath::Vvtvp(nq, &v2[1][0], 1, &v2[1][0], 1, &magorig[0], 1,
                         &magorig[0], 1);
            Vmath::Vvtvp(nq, &v2[2][0], 1, &v2[2][0], 1, &magorig[0], 1,
                         &magorig[0], 1);

            Vmath::Vmul(nq, &outarray[0][0], 1, &outarray[0][0], 1, &magnew[0],
                        1);
            Vmath::Vvtvp(nq, &outarray[1][0], 1, &outarray[1][0], 1, &magnew[0],
                         1, &magnew[0], 1);
            Vmath::Vvtvp(nq, &outarray[2][0], 1, &outarray[2][0], 1, &magnew[0],
                         1, &magnew[0], 1);
        }

        for (int i = 0; i < m_spacedim; ++i)
        {
            for (int j = 0; j < nq; ++j)
            {
                if (fabs(magnew[j]) > 0.000000001)
                {
                    outarray[i][j] =
                        outarray[i][j] * sqrt(magorig[j] / magnew[j]);
                }
            }
        }
    }
}

void MMFSystem::PlotMovingFrames(
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    const int nstep)
{
    int nvar    = 12;
    int nq      = m_fields[0]->GetTotPoints();
    int ncoeffs = m_fields[0]->GetNcoeffs();

    std::string outname1;
    outname1 = m_sessionName + "_MF_" +
               boost::lexical_cast<std::string>(nstep) + ".chk";
    std::vector<Array<OneD, NekDouble>> fieldcoeffs(nvar);
    for (int i = 0; i < nvar; ++i)
    {
        fieldcoeffs[i] = Array<OneD, NekDouble>(ncoeffs);
    }

    std::vector<std::string> variables(nvar);
    variables[0]  = "ex1";
    variables[1]  = "ey1";
    variables[2]  = "ez1";
    variables[3]  = "ex2";
    variables[4]  = "ey2";
    variables[5]  = "ez2";
    variables[6]  = "ex3";
    variables[7]  = "ey3";
    variables[8]  = "ez3";
    variables[9]  = "e1mag";
    variables[10] = "e2mag";
    variables[11] = "e3mag";

    Array<OneD, NekDouble> tmp(nq);
    Array<OneD, NekDouble> mag(nq);
    for (int i = 0; i < m_mfdim; ++i)
    {
        mag = Array<OneD, NekDouble>(nq, 0.0);
        for (int j = 0; j < m_spacedim; ++j)
        {
            Vmath::Vcopy(nq, &movingframes[i][j * nq], 1, &tmp[0], 1);
            m_fields[0]->FwdTransLocalElmt(tmp, fieldcoeffs[i * 3 + j]);

            Vmath::Vvtvp(nq, tmp, 1, tmp, 1, mag, 1, mag, 1);
        }

        Vmath::Vsqrt(nq, mag, 1, mag, 1);
        m_fields[0]->FwdTransLocalElmt(mag, fieldcoeffs[9 + i]);
    }

    WriteFld(outname1, m_fields[0], fieldcoeffs, variables);
}

// PlotMFInfo(fields, velmag, MF1stAligned[0],
// ActivationIntensityHistory, ErrorIntensityHistory, ActivatedCount,
// ZoneActivation, nchk);
void MMFSystem::PlotTrajectoryMF(
    const Array<OneD, const int> &ActivatedHistory,
    const Array<OneD, const NekDouble> &fieldu,
    const Array<OneD, const Array<OneD, NekDouble>> &MF,
    const Array<OneD, const Array<OneD, Array<OneD, NekDouble>>>
        &MF1stConnection,
    const Array<OneD, const Array<OneD, NekDouble>> &Relacc,
    Array<OneD, NekDouble> &NoBoundaryZone, const int nstep)
{
    int nvar    = 9;
    int nq      = m_fields[0]->GetTotPoints();
    int ncoeffs = m_fields[0]->GetNcoeffs();

    std::string outname1 = m_sessionName + "_MFTRJ_" +
                           boost::lexical_cast<std::string>(nstep) + ".chk";
    std::vector<Array<OneD, NekDouble>> fieldcoeffs(nvar);
    for (int i = 0; i < nvar; ++i)
    {
        fieldcoeffs[i] = Array<OneD, NekDouble>(ncoeffs);
    }

    std::vector<std::string> variables(nvar);
    variables[0] = "fieldu";
    variables[1] = "ex1";
    variables[2] = "ey1";
    variables[3] = "ez1";
    variables[4] = "Activated";
    variables[5] = "w211";
    variables[6] = "w212";
    variables[7] = "RelAcc";
    variables[8] = "CondBlock";

    // index:0 -> u
    m_fields[0]->FwdTransLocalElmt(fieldu, fieldcoeffs[0]);

    // index:[1, 2, 3] -> ex1, ey1, ez1
    Array<OneD, NekDouble> tmp(nq);
    for (int j = 0; j < m_spacedim; ++j)
    {
        Vmath::Vcopy(nq, &MF[0][j * nq], 1, &tmp[0], 1);
        m_fields[0]->FwdTransLocalElmt(tmp, fieldcoeffs[j + 1]);
    }

    // Angle between MF and fibre
    Array<OneD, NekDouble> Acttmp(nq, 0.0);
    for (int i = 0; i < nq; ++i)
    {
        if (ActivatedHistory[i] == 1)
        {
            Acttmp[i] = 1.0;
        }
    }
    m_fields[0]->FwdTransLocalElmt(Acttmp, fieldcoeffs[4]);

    Array<OneD, NekDouble> w211(nq);
    Array<OneD, NekDouble> w212(nq);

    Vmath::Vcopy(nq, &MF1stConnection[0][0][0], 1, &w211[0], 1);
    Vmath::Vcopy(nq, &MF1stConnection[0][1][0], 1, &w212[0], 1);

    // Elementwise activation if w211 < Tol and ActivationHistory==1
    NoBoundaryZone = Array<OneD, NekDouble>(nq, 1.0);

    // Array<OneD, Array<OneD, int>> EWIndex;
    // m_fields[0]->GridIndexElementWise(EWIndex);

    // int Nelemtj = EWIndex.size();
    // int nptsj   = EWIndex[0].size();
    // int index, Actflag=0;
    // NekDouble Tol=1.0;
    // for (int i = 0; i < Nelemtj; ++i)
    // {
    //     Actflag = 1;
    //     for (int j = 0; j < nptsj; ++j)
    //     {
    //         index = EWIndex[i][j];
    //         if( (fabs(w211[index])>Tol) || (ActivatedHistory[index] == 0) )
    //         {
    //             Actflag = 0;
    //             std::cout << "Opted out: w211 = " << w211[index] <<
    //             std::endl;
    //         }
    //     }

    //     for (int j = 0; j < nptsj; ++j)
    //     {
    //         index = EWIndex[i][j];
    //         if(Actflag == 0)
    //         {
    //             NoBoundaryZone[index] = 0.0;
    //         }
    //     }
    // }

    std::cout << "Trajectory: w211 max = " << Vmath::Vamax(nq, w211, 1)
              << std::endl;

    // Ignore the region where w211 is too big.
    Vmath::Vmul(nq, NoBoundaryZone, 1, w211, 1, w211, 1);
    Vmath::Vmul(nq, NoBoundaryZone, 1, w212, 1, w212, 1);

    // Connection form w211
    m_fields[0]->FwdTransLocalElmt(w211, fieldcoeffs[5]);

    // Connection form w212
    m_fields[0]->FwdTransLocalElmt(w212, fieldcoeffs[6]);

    // Relative Acceleration (I)
    Array<OneD, NekDouble> RelAccetmp(nq);

    Vmath::Vcopy(nq, &Relacc[1][0], 1, &RelAccetmp[0], 1);
    Vmath::Vmul(nq, &NoBoundaryZone[0], 1, &RelAccetmp[0], 1, &RelAccetmp[0],
                1);

    m_fields[0]->FwdTransLocalElmt(RelAccetmp, fieldcoeffs[7]);

    // Compute Conduction block
    // Array<OneD, NekDouble> CBlock(nq,0.0);
    // for (int i=0; i<nq; ++i)
    // {
    //     if( (w212[i]>0) && (RelAccetmp[i]<0) )
    //     {
    //         CBlock[i] = w212[i] - 10.0 * RelAccetmp[i];
    //     }
    // }

    Array<OneD, NekDouble> CBlock(nq, 0.0);
    // NekDouble TolRelAcc = 0.01;
    for (int i = 0; i < nq; ++i)
    {
        if ((w212[i] < 0) && (RelAccetmp[i] > 0))
        {
            CBlock[i] = RelAccetmp[i];
        }
    }

    m_fields[0]->FwdTransLocalElmt(CBlock, fieldcoeffs[8]);

    WriteFld(outname1, m_fields[0], fieldcoeffs, variables);
}

// Riemann Curvature, R2121
// m_fields[0]->FwdTrans(MF1stCurvature[0], fieldcoeffs[7]);

// Divergence of e1
// Array<OneD, Array<OneD, NekDouble>> DivMF;
// ComputeDivMF(eCovariant, MF, DivMF);

// Vmath::Vmul(nq, &Activation[0], 1, &DivMF[0][0], 1, &DivMF[0][0], 1);
// m_fields[0]->FwdTrans(DivMF[0], fieldcoeffs[7]);

// // Curl of e1
// Array<OneD, Array<OneD, NekDouble>> CurlMF(m_mfdim);
// for (int j = 0; j < m_mfdim; ++j)
// {
//     CurlMF[j] = Array<OneD, NekDouble>(nq, 0.0);
// }

// ComputeCurlMF(eCovariant, MF, CurlMF);
// Vmath::Vmul(nq, &Activation[0], 1, &CurlMF[0][0], 1, &CurlMF[0][0], 1);
// m_fields[0]->FwdTrans(CurlMF[0], fieldcoeffs[8]);

// // Mean curvature: H = 0.5 ( w311 + w322)
// Array<OneD, NekDouble> Mean(nq);
// Vmath::Vadd(nq, &MF1stConnection[1][0][0], 1, &MF1stConnection[2][1][0], 1,
//             &Mean[0], 1);
// Vmath::Smul(nq, 0.5, &Mean[0], 1, &Mean[0], 1);
// m_fields[0]->FwdTrans(Mean, fieldcoeffs[10]);

// // Gaussian Curvature: K = w311 * w322 - w312 * w321
// Array<OneD, NekDouble> Gauss(nq);
// Vmath::Vmul(nq, &MF1stConnection[1][0][0], 1, &MF1stConnection[2][1][0], 1,
//             &Gauss[0], 1);
// Vmath::Vmul(nq, &MF1stConnection[1][1][0], 1, &MF1stConnection[2][0][0], 1,
//             &tmp[0], 1);
// Vmath::Neg(nq, tmp, 1);
// Vmath::Vadd(nq, tmp, 1, Gauss, 1, Gauss, 1);

// m_fields[0]->FwdTrans(Gauss, fieldcoeffs[11]);

// void MMFSystem::PlotRelacc2D(
//     const Array<OneD, const int> &ActivatedHistory,
//     const Array<OneD, const Array<OneD, Array<OneD, NekDouble>>>
//     &MFConnection, const Array<OneD, const Array<OneD, NekDouble>>
//     MFCurvature, const Array<OneD, const Array<OneD, NekDouble>> &Relacc,
//     const Array<OneD, const Array<OneD, NekDouble>> &RelaccOmega,
//     const int nstep)
// {
//     int nvar    = 5;
//     int nq      = m_fields[0]->GetTotPoints();
//     int ncoeffs = m_fields[0]->GetNcoeffs();

//     std::string outname1 = m_sessionName + "_RelAcc_" +
//                            boost::lexical_cast<std::string>(nstep) + ".chk";

//     std::vector<Array<OneD, NekDouble>> fieldcoeffs(nvar);
//     for (int i = 0; i < nvar; ++i)
//     {
//         fieldcoeffs[i] = Array<OneD, NekDouble>(ncoeffs);
//     }

//     std::vector<std::string> variables(nvar);
//     variables[0]  = "w212";
//     variables[1]  = "R2121";
//     variables[2]  = "RelAcc";
//     variables[3]  = "RelAccOmega";
//     variables[4]  = "Conduction Block";

//     Array<OneD, NekDouble> w211(nq);
//     Array<OneD, NekDouble> w212(nq);
//     Array<OneD, NekDouble> R2121(nq);

//     Vmath::Vcopy(nq, &MFConnection[0][0][0], 1, &w211[0], 1);
//     Vmath::Vcopy(nq, &MFConnection[0][1][0], 1, &w212[0], 1);

//     Vmath::Vcopy(nq, &MFCurvature[0][0], 1, &R2121[0], 1);

//     Array<OneD, NekDouble> Activation(nq, 1.0);
//     NekDouble Tol=0.1;
//     for (int i=0; i<nq; ++i)
//     {
//         if(fabs(w211[i])>Tol)
//         {
//             Activation[i] = 0.0;
//         }

//         if (ActivatedHistory[i] <= 0)
//         {
//             Activation[i] = 0.0;
//         }
//     }

//     // Ignore the region where w211 is too big.
//     Vmath::Vmul(nq, Activation, 1, w211, 1, w211, 1);
//     Vmath::Vmul(nq, Activation, 1, w212, 1, w212, 1);

//     // Connection form w212
//     m_fields[0]->FwdTrans(w212, fieldcoeffs[0]);

//     // Curvature
//     Vmath::Vmul(nq, Activation, 1, R2121, 1, R2121, 1);
//     m_fields[0]->FwdTrans(R2121, fieldcoeffs[1]);

//    // Relative Acceleration (I)
//     Array<OneD, NekDouble> RelAcce2(nq);

//     Vmath::Vcopy(nq, &Relacc[1][0], 1, &RelAcce2[0], 1);
//     Vmath::Vmul(nq, &Activation[0], 1, &RelAcce2[0], 1, &RelAcce2[0], 1);

//     m_fields[0]->FwdTrans(RelAcce2, fieldcoeffs[2]);

//    // Relative Acceleration (II)
//     Array<OneD, NekDouble> RelAcce2Omega(nq);

//     Vmath::Vcopy(nq, &RelaccOmega[1][0], 1, &RelAcce2Omega[0], 1);
//     Vmath::Vmul(nq, &Activation[0], 1, &RelAcce2Omega[0], 1,
//     &RelAcce2Omega[0], 1);

//     m_fields[0]->FwdTrans(RelAcce2, fieldcoeffs[3]);

//     // Compute Conduction block
//     Array<OneD, NekDouble> Conductionblock(nq,0.0);
//     for (int i=0; i<nq; ++i)
//     {
//         if( (w212[i]>0) && (RelAcce2[i]<0) )
//         {
//             Conductionblock[i] = w212[i];
//         }
//     }
//     m_fields[0]->FwdTrans(Conductionblock, fieldcoeffs[4]);

//     WriteFld(outname1, m_fields[0], fieldcoeffs, variables);
// }


// void MMFSystem::TimeMapforInitZone(
//     const Array<OneD, const int> &ValidTimeMap,
//     const Array<OneD, const NekDouble> &dudtHistory,
//     Array<OneD, NekDouble> &TimeMap)
// {
//     int nq = GetTotPoints();

//     // Substitute the initial excited zone with the minimum or maximum of the
//     // Time Map depending on the distance

//     // Find min and max Timemap for which integration has been computed.
//     NekDouble TimeMapMin = 0.0, TimeMapMax = 0.0;
//     int ITMmin = 0, ITMmax = 0;
//     NekDouble dudtTol = 0.01;
//     for (int i = 0; i < nq; ++i)
//     {
//         if (fabs(dudtHistory[i]) > dudtTol)
//         {
//             if (TimeMapMin > TimeMap[i])
//             {
//                 TimeMapMin = TimeMap[i];
//                 ITMmin     = i;
//             }

//             if (TimeMapMax < TimeMap[i])
//             {
//                 TimeMapMax = TimeMap[i];
//                 ITMmax     = i;
//             }
//         }
//     }

//     Array<OneD, NekDouble> x(nq);
//     Array<OneD, NekDouble> y(nq);
//     Array<OneD, NekDouble> z(nq);

//     m_fields[0]->GetCoords(x, y, z);

//     NekDouble xp1, yp1, zp1, xp2, yp2, zp2;
//     NekDouble distmin, distmax;
//     for (int i = 0; i < nq; ++i)
//     {
//         if (ValidTimeMap[i] == 0)
//         {
//             xp1 = x[i] - x[ITMmin];
//             yp1 = y[i] - y[ITMmin];
//             zp1 = z[i] - z[ITMmin];

//             distmin = sqrt(xp1 * xp1 + yp1 * yp1 + zp1 * zp1);

//             xp2 = x[i] - x[ITMmax];
//             yp2 = y[i] - y[ITMmax];
//             zp2 = z[i] - z[ITMmax];

//             distmax = sqrt(xp2 * xp2 + yp2 * yp2 + zp2 * zp2);

//             if (distmin < distmax)
//             {
//                 TimeMap[i] = TimeMapMin;
//             }

//             else
//             {
//                 TimeMap[i] = TimeMapMax;
//             }
//         }
//     }
// }


Array<OneD, NekDouble> MMFSystem::ComputeSpaceTime(
    const Array<OneD, const NekDouble> &TimeMap)
{
    int nq = m_fields[0]->GetTotPoints();

    NekDouble Tmax = Vmath::Vmax(nq, TimeMap, 1);
    NekDouble Tmin = Vmath::Vmin(nq, TimeMap, 1);

    Array<OneD, NekDouble> x0(nq);
    Array<OneD, NekDouble> x1(nq);
    Array<OneD, NekDouble> x2(nq);

    m_fields[0]->GetCoords(x0, x1, x2);

    int ROIindx = 0;

    Array<OneD, NekDouble> outarray(nq);

    NekDouble dx, dy, dz, rad, xdist, tdist;
    NekDouble maxdist = 0.0, dist = 100.0, cspeed;
    // Compute speed of the wave
    if (m_surfaceType == eSphere)
    {
        maxdist = 2 * m_pi * m_SphereExactRadius;
    }

    else
    {
        for (int i = 0; i < nq; ++i)
        {
            dx = x0[i] - m_ROIx;
            dy = x1[i] - m_ROIy;
            dz = x2[i] - m_ROIz;

            rad = sqrt(dx * dx + dy * dy + dz * dz);

            if (rad < dist)
            {
                ROIindx = i;
            }

            if (rad > maxdist)
            {
                maxdist = rad;
            }
        }
    }

    // Compute the speed of info.
    cspeed = maxdist / Tmax;

    std::cout << "TimeMap: Tmin = " << Tmin << ", Tmax = " << Tmax
              << ", speed = " << cspeed << std::endl;

    for (int i = 0; i < nq; ++i)
    {
        rad = sqrt(x0[i] * x0[i] + x1[i] * x1[i] + x2[i] * x2[i]);

        dx = x0[i] - m_ROIx;
        dy = x1[i] - m_ROIy;
        dz = x2[i] - m_ROIz;

        xdist = sqrt(dx * dx + dy * dy + dz * dz);
        if (m_surfaceType == eSphere)
        {
            xdist = rad * 2 * asin(xdist / m_SphereExactRadius / 2.0);
        }

        tdist       = cspeed * (TimeMap[i] - TimeMap[ROIindx]);
        outarray[i] = -tdist * tdist + xdist * xdist;
    }

    return outarray;
}

Array<OneD, NekDouble> MMFSystem::ComputeDirectionVector(
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    const Array<OneD, const NekDouble> &inarray,
    const Array<OneD, const NekDouble> &dudt)

{
    int nq = m_fields[0]->GetTotPoints();

    Array<OneD, Array<OneD, NekDouble>> qfield(m_expdim);
    for (int j = 0; j < m_expdim; ++j)
    {
        qfield[j] = Array<OneD, NekDouble>(nq, 0.0);
        MMFDirectionalDeriv(movingframes[j], inarray, qfield[j]);
    }

    Array<OneD, NekDouble> dudtapplied(nq, 0.0);
    Vmath::Vcopy(nq, dudt, 1, dudtapplied, 1);
    for (int i = 0; i < nq; i++)
    {
        if ((m_GradLocType == eWaveback) && (dudt[i] > 0))
        {
            dudtapplied[i] = 0.0;
        }

        if ((m_GradLocType == eWavefront) && (dudt[i] < 0))
        {
            dudtapplied[i] = 0.0;
        }
    }

    Array<OneD, NekDouble> outarray(m_spacedim * nq);
    for (int i = 0; i < nq; i++)
    {
        outarray[i] = -1.0 * dudtapplied[i] *
                      (qfield[0][i] * movingframes[0][i] +
                       qfield[1][i] * movingframes[1][i]);

        outarray[i + nq] = -1.0 * dudtapplied[i] *
                           (qfield[0][i] * movingframes[0][i + nq] +
                            qfield[1][i] * movingframes[1][i + nq]);

        outarray[i + 2 * nq] = -1.0 * dudtapplied[i] *
                               (qfield[0][i] * movingframes[0][i + 2 * nq] +
                                qfield[1][i] * movingframes[1][i + 2 * nq]);
    }

    return outarray;
}

// If velmag is larger than Tol, nullify the vector
void MMFSystem::VectorCutOff(const NekDouble Tol,
                             Array<OneD, NekDouble> &vector)
{
    int nq = m_fields[0]->GetTotPoints();

    int index, Velmagflag = 0;
    NekDouble vx, vy, vz, velmag;
    for (int i = 0; i < m_fields[0]->GetExpSize(); ++i)
    {
        Velmagflag = 0;
        for (int j = 0; j < m_fields[0]->GetTotPoints(i); ++j)
        {
            index = m_fields[0]->GetPhys_Offset(i) + j;

            vx = vector[index];
            vy = vector[index + nq];
            vz = vector[index + 2 * nq];

            velmag = sqrt(vx * vx + vy * vy + vz * vz);

            if (velmag > Tol)
            {
                Velmagflag = 1;
            }
        }

        if (Velmagflag == 1)
        {
            for (int j = 0; j < m_fields[0]->GetTotPoints(i); ++j)
            {
                index                  = m_fields[0]->GetPhys_Offset(i) + j;
                vector[index]          = 0.0;
                vector[index + nq]     = 0.0;
                vector[index + 2 * nq] = 0.0;
            }
        }
    }
}

Array<OneD, int> MMFSystem::ComputeZoneActivation(
    const NekDouble uTol, const Array<OneD, const NekDouble> &fieldu,
    const NekDouble NoAlignInitRadius)
{
    int nq = GetNpoints();

    Array<OneD, NekDouble> x0(nq);
    Array<OneD, NekDouble> x1(nq);
    Array<OneD, NekDouble> x2(nq);

    m_fields[0]->GetCoords(x0, x1, x2);

    NekDouble distx, disty, distz, rad;

    Array<OneD, int> outarray(nq, 0);
    for (int i = 0; i < nq; ++i)
    {
        // Only consider the region where u is not negligible
        if (fabs(fieldu[i]) > uTol)
        {
            outarray[i] = 1;
        }

        else
        {
            outarray[i] = 0;
        }

        distx = x0[i] - m_Initx;
        disty = x1[i] - m_Inity;
        distz = x2[i] - m_Initz;

        rad = sqrt(distx * distx + disty * disty + distz * distz);
        if (rad < NoAlignInitRadius)
        {
            // outarray[i] = 0;
        }
    }

    return outarray;
}

// void MMFSystem::PlotTimeMapAnalysis(
//     const Array<OneD, const int> &ActivatedHistory,
//     const Array<OneD, const Array<OneD, NekDouble>> &MF,
//     const Array<OneD, const NekDouble> &TimeMap,
//     const Array<OneD, const int> &APindex, const int nstep)
// {
//     int nvar    = 16;
//     int nq      = m_fields[0]->GetTotPoints();
//     int ncoeffs = m_fields[0]->GetNcoeffs();

//     std::string outname1 = m_sessionName + "_TmapAnalysis_" +
//                            boost::lexical_cast<std::string>(nstep) + ".chk";
//     std::vector<Array<OneD, NekDouble>> fieldcoeffs(nvar);
//     for (int i = 0; i < nvar; ++i)
//     {
//         fieldcoeffs[i] = Array<OneD, NekDouble>(ncoeffs);
//     }

//     std::vector<std::string> variables(nvar);
//     variables[0]  = "ex1";
//     variables[1]  = "ey1";
//     variables[2]  = "ez1";
//     variables[3]  = "TMaprel";
//     variables[4]  = "tx";
//     variables[5]  = "ty";
//     variables[6]  = "tz";
//     variables[7]  = "tvecmag";
//     variables[8]  = "wavevel";
//     variables[9]  = "tcdotMF";
//     variables[10] = "w211_MFvel";
//     variables[11] = "w212_MFvel";
//     variables[12] = "R2121_MFvel";
//     variables[13] = "theta";
//     variables[14] = "rad";
//     variables[15] = "APindex";

//     Array<OneD, NekDouble> tmp(nq);
//     for (int j = 0; j < m_spacedim; ++j)
//     {
//         Vmath::Vcopy(nq, &MF[0][j * nq], 1, &tmp[0], 1);
//         m_fields[0]->FwdTrans(tmp, fieldcoeffs[j]);
//     }

//     // Normalized Time Vector
//     NekDouble Tmax = Vmath::Vamax(nq, TimeMap, 1);
//     Vmath::Smul(nq, 1.0 / Tmax, TimeMap, 1, tmp, 1);
//     m_fields[0]->FwdTrans(tmp, fieldcoeffs[3]);

//     // Time vector is the gradient of time map
//     Array<OneD, NekDouble> tvec(m_spacedim * nq);
//     tvec = ComputeCovGrad(TimeMap, m_movingframes);

//     // Compute the connection of the time vector map
//     for (int k = 0; k < m_spacedim; ++k)
//     {
//         Vmath::Vcopy(nq, &tvec[k * nq], 1, &tmp[0], 1);
//         m_fields[0]->FwdTrans(tmp, fieldcoeffs[k + 4]);
//     }

//     NekDouble tx, ty, tz;
//     NekDouble e1x, e1y, e1z, MFmag, dottmp;

//     Array<OneD, NekDouble> AngleTimeMF(nq);
//     Array<OneD, NekDouble> tvecmag(nq, 0.0);
//     Array<OneD, NekDouble> wavevel(nq, 0.0);
//     int cnt       = 0;
//     NekDouble Tol = 0.0001;
//     for (int i = 0; i < nq; i++)
//     {
//         e1x = MF[0][i];
//         e1y = MF[0][i + nq];
//         e1z = MF[0][i + 2 * nq];

//         tx = tvec[i];
//         ty = tvec[i + nq];
//         tz = tvec[i + 2 * nq];

//         if (ActivatedHistory[i] > 0)
//         {
//             MFmag = sqrt(e1x * e1x + e1y * e1y + e1z * e1z);

//             tvecmag[i] = sqrt(tx * tx + ty * ty + tz * tz);
//             if (tvecmag[i] > Tol)
//             {
//                 wavevel[i] = 1.0 / tvecmag[i];
//             }

//             dottmp         = fabs(e1x * tx + e1y * ty + e1z * tz);
//             AngleTimeMF[i] = dottmp / MFmag / tvecmag[i];
//             cnt++;
//         }
//     }

//     m_fields[0]->FwdTrans(tvecmag, fieldcoeffs[7]);
//     m_fields[0]->FwdTrans(wavevel, fieldcoeffs[8]);
//     m_fields[0]->FwdTrans(AngleTimeMF, fieldcoeffs[9]);

//     // Compute the connection form and Riemann curvature of the MF with
//     moving
//     // frames.
//     Array<OneD, Array<OneD, NekDouble>> MF1stwithVel(m_spacedim);
//     for (int i = 0; i < m_spacedim; ++i)
//     {
//         MF1stwithVel[i] = Array<OneD, NekDouble>(3 * nq, 0.0);

//         Vmath::Smul(3 * nq, 1.0, &MF[i][0], 1, &MF1stwithVel[i][0], 1);

//         // multiply tvec to the first moving frames
//         Vmath::Vmul(nq, &tvecmag[0], 1, &MF1stwithVel[0][i * nq], 1,
//                     &MF1stwithVel[0][i * nq], 1);
//     }

//     Array<OneD, Array<OneD, Array<OneD, NekDouble>>> MF1stVelConnection(
//         m_mfdim);
//     Array<OneD, Array<OneD, NekDouble>> MF1stVelCurvature(m_mfdim);
//     for (int i = 0; i < m_mfdim; i++)
//     {
//         MF1stVelCurvature[i] = Array<OneD, NekDouble>(nq, 0.0);

//         MF1stVelConnection[i] = Array<OneD, Array<OneD, NekDouble>>(m_mfdim);
//         for (int j = 0; j < m_mfdim; j++)
//         {
//             MF1stVelConnection[i][j] = Array<OneD, NekDouble>(nq, 0.0);
//         }
//     }

//     Compute2DConnectionCurvature(MF1stwithVel, MF1stVelConnection,
//                                  MF1stVelCurvature);

//     // Connection form w211
//     m_fields[0]->FwdTrans(MF1stVelConnection[0][0], fieldcoeffs[10]);

//     // Connection form w212
//     m_fields[0]->FwdTrans(MF1stVelConnection[0][1], fieldcoeffs[11]);

//     // Riemann Curvature, R2121
//     m_fields[0]->FwdTrans(MF1stVelCurvature[0], fieldcoeffs[12]);

//     // Compute theta and rad for test
//     Array<OneD, NekDouble> Spheretheta(nq, 0.0);
//     Array<OneD, NekDouble> radius(nq, 0.0);
//     NekDouble xp, yp, zp, rad;
//     NekDouble sin_varphi, cos_varphi, sin_theta, cos_theta;

//     Array<OneD, NekDouble> x0(nq);
//     Array<OneD, NekDouble> x1(nq);
//     Array<OneD, NekDouble> x2(nq);

//     m_fields[0]->GetCoords(x0, x1, x2);
//     for (int i = 0; i < nq; i++)
//     {
//         xp = x0[i];
//         yp = x1[i];
//         zp = x2[i];

//         CartesianToNewSpherical(xp, yp, zp, rad, sin_varphi, cos_varphi,
//                                 sin_theta, cos_theta);

//         if (ActivatedHistory[i] > 0)
//         {
//             Spheretheta[i] = atan2(sin_theta, cos_theta);
//             radius[i]      = sqrt(xp * xp + yp * yp);
//         }
//     }

//     m_fields[0]->FwdTrans(Spheretheta, fieldcoeffs[13]);
//     m_fields[0]->FwdTrans(radius, fieldcoeffs[14]);

//     for (int i = 0; i < nq; i++)
//     {
//         tmp[i] = 1.0 * APindex[i];
//     }
//     m_fields[0]->FwdTrans(tmp, fieldcoeffs[15]);

//     WriteFld(outname1, m_fields[0], fieldcoeffs, variables);
// }

// PlotMFInfo(fields, velmag, MF1stAligned[0],
// ActivationIntensityHistory, ErrorIntensityHistory, ActivatedCount,
// ZoneActivation, nchk);
void MMFSystem::PlotHHD(const Array<OneD, const NekDouble> &MFvec,
                        const Array<OneD, const NekDouble> &irrotational,
                        const Array<OneD, const NekDouble> &incompressible,
                        const Array<OneD, const NekDouble> &harmonic,
                        const int nstep)
{
    int nvar    = 21;
    int nq      = m_fields[0]->GetTotPoints();
    int ncoeffs = m_fields[0]->GetNcoeffs();

    std::string outname1;
    outname1 = m_sessionName + "_HHD_" +
               boost::lexical_cast<std::string>(nstep) + ".chk";

    std::vector<Array<OneD, NekDouble>> fieldcoeffs(nvar);
    for (int i = 0; i < nvar; ++i)
    {
        fieldcoeffs[i] = Array<OneD, NekDouble>(ncoeffs);
    }

    std::vector<std::string> variables(nvar);
    variables[0]  = "irrx";
    variables[1]  = "irry";
    variables[2]  = "irrz";
    variables[3]  = "irrmag";
    variables[4]  = "irrSource";
    variables[5]  = "irrPotential";
    variables[6]  = "incompx";
    variables[7]  = "incompy";
    variables[8]  = "incompz";
    variables[9]  = "incompmag";
    variables[10] = "incompSource";
    variables[11] = "incompPotential";
    variables[12] = "harmox";
    variables[13] = "harmoy";
    variables[14] = "harmoz";
    variables[15] = "harmomag";
    variables[16] = "vLapharmonicmag";
    variables[17] = "harmonPotential";
    variables[18] = "ex";
    variables[19] = "ey";
    variables[20] = "ez";

    Array<OneD, NekDouble> tmp(nq);
    for (int k = 0; k < m_spacedim; ++k)
    {
        Vmath::Vcopy(nq, &MFvec[k * nq], 1, &tmp[0], 1);
        m_fields[0]->FwdTransLocalElmt(tmp, fieldcoeffs[k + 18]);

        Vmath::Vcopy(nq, &irrotational[k * nq], 1, &tmp[0], 1);
        m_fields[0]->FwdTransLocalElmt(tmp, fieldcoeffs[k]);

        Vmath::Vcopy(nq, &incompressible[k * nq], 1, &tmp[0], 1);
        m_fields[0]->FwdTransLocalElmt(tmp, fieldcoeffs[k + 6]);

        Vmath::Vcopy(nq, &harmonic[k * nq], 1, &tmp[0], 1);
        m_fields[0]->FwdTransLocalElmt(tmp, fieldcoeffs[k + 12]);
    }

    // Irrotational
    Vmath::Vcopy(nq, &irrotational[3 * nq], 1, &tmp[0], 1);
    m_fields[0]->FwdTransLocalElmt(tmp, fieldcoeffs[3]);

    Vmath::Vcopy(nq, &irrotational[4 * nq], 1, &tmp[0], 1);
    m_fields[0]->FwdTransLocalElmt(tmp, fieldcoeffs[4]);

    Vmath::Vcopy(nq, &irrotational[5 * nq], 1, &tmp[0], 1);
    m_fields[0]->FwdTransLocalElmt(tmp, fieldcoeffs[5]);

    // Incompressible
    Vmath::Vcopy(nq, &incompressible[3 * nq], 1, &tmp[0], 1);
    m_fields[0]->FwdTransLocalElmt(tmp, fieldcoeffs[9]);

    Vmath::Vcopy(nq, &incompressible[4 * nq], 1, &tmp[0], 1);
    m_fields[0]->FwdTransLocalElmt(tmp, fieldcoeffs[10]);

    Vmath::Vcopy(nq, &incompressible[5 * nq], 1, &tmp[0], 1);
    m_fields[0]->FwdTransLocalElmt(tmp, fieldcoeffs[11]);

    // Harmonic
    Vmath::Vcopy(nq, &harmonic[3 * nq], 1, &tmp[0], 1);
    m_fields[0]->FwdTransLocalElmt(tmp, fieldcoeffs[15]);

    Vmath::Vcopy(nq, &harmonic[4 * nq], 1, &tmp[0], 1);
    m_fields[0]->FwdTransLocalElmt(tmp, fieldcoeffs[16]);

    Vmath::Vcopy(nq, &harmonic[5 * nq], 1, &tmp[0], 1);
    m_fields[0]->FwdTransLocalElmt(tmp, fieldcoeffs[17]);

    WriteFld(outname1, m_fields[0], fieldcoeffs, variables);
}

void MMFSystem::PlotCardiacFibre(const Array<OneD, const NekDouble> &fibre)
{
    int nvar    = 4;
    int nq      = m_fields[0]->GetTotPoints();
    int ncoeffs = m_fields[0]->GetNcoeffs();

    std::string outname1;
    outname1 = m_sessionName + "_CardiacFibre.chk";

    std::vector<Array<OneD, NekDouble>> fieldcoeffs(nvar);
    for (int i = 0; i < nvar; ++i)
    {
        fieldcoeffs[i] = Array<OneD, NekDouble>(ncoeffs);
    }

    std::vector<std::string> variables(nvar);
    variables[0] = "fx";
    variables[1] = "fy";
    variables[2] = "fz";
    variables[3] = "mag";

    Array<OneD, NekDouble> tmp(nq);
    Array<OneD, NekDouble> mag(nq, 0.0);
    for (int j = 0; j < m_spacedim; ++j)
    {
        Vmath::Vcopy(nq, &fibre[j * nq], 1, &tmp[0], 1);
        m_fields[0]->FwdTransLocalElmt(tmp, fieldcoeffs[j]);

        Vmath::Vvtvp(nq, &tmp[0], 1, &tmp[0], 1, &mag[0], 1, &mag[0], 1);
    }
    Vmath::Vsqrt(nq, &mag[0], 1, &mag[0], 1);
    m_fields[0]->FwdTransLocalElmt(mag, fieldcoeffs[3]);

    WriteFld(outname1, m_fields[0], fieldcoeffs, variables);
}

void MMFSystem::PlotProcessedCardiacFibre(
    const Array<OneD, const NekDouble> &fibre,
    const Array<OneD, const NekDouble> &fcdotk,
    const Array<OneD, const NekDouble> &AniConstruction)
{
    int nvar    = 6;
    int nq      = m_fields[0]->GetTotPoints();
    int ncoeffs = m_fields[0]->GetNcoeffs();

    std::string outname1 = m_sessionName + "_ProcessedCardiacFibre.chk";
    std::vector<Array<OneD, NekDouble>> fieldcoeffs(nvar);
    for (int i = 0; i < nvar; ++i)
    {
        fieldcoeffs[i] = Array<OneD, NekDouble>(ncoeffs);
    }

    std::vector<std::string> variables(nvar);
    variables[0] = "fx";
    variables[1] = "fy";
    variables[2] = "fz";
    variables[3] = "mag";
    variables[4] = "fcdotk";
    variables[5] = "AniConstruction";

    Array<OneD, NekDouble> tmp(nq);
    Array<OneD, NekDouble> mag(nq, 0.0);
    for (int j = 0; j < m_spacedim; ++j)
    {
        Vmath::Vcopy(nq, &fibre[j * nq], 1, &tmp[0], 1);
        m_fields[0]->FwdTransLocalElmt(tmp, fieldcoeffs[j]);

        Vmath::Vvtvp(nq, &tmp[0], 1, &tmp[0], 1, &mag[0], 1, &mag[0], 1);
    }
    Vmath::Vsqrt(nq, &mag[0], 1, &mag[0], 1);
    m_fields[0]->FwdTransLocalElmt(mag, fieldcoeffs[3]);

    m_fields[0]->FwdTransLocalElmt(fcdotk, fieldcoeffs[4]);
    m_fields[0]->FwdTransLocalElmt(AniConstruction, fieldcoeffs[5]);

    WriteFld(outname1, m_fields[0], fieldcoeffs, variables);
}

void MMFSystem::PlotMovingFrames(
    const int dir,
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    const int nstep)
{
    int nvar    = 3;
    int nq      = m_fields[0]->GetTotPoints();
    int ncoeffs = m_fields[0]->GetNcoeffs();

    std::string outname1;
    outname1 = m_sessionName + "_MF_" + boost::lexical_cast<std::string>(dir) +
               "_" + boost::lexical_cast<std::string>(nstep) + ".chk";
    std::vector<Array<OneD, NekDouble>> fieldcoeffs(nvar);
    for (int i = 0; i < nvar; ++i)
    {
        fieldcoeffs[i] = Array<OneD, NekDouble>(ncoeffs);
    }

    std::vector<std::string> variables(nvar);
    variables[0] = "ex";
    variables[1] = "ey";
    variables[2] = "ez";

    Array<OneD, NekDouble> tmp(nq);
    for (int j = 0; j < m_spacedim; ++j)
    {
        Vmath::Vcopy(nq, &movingframes[dir][j * nq], 1, &tmp[0], 1);
        m_fields[0]->FwdTransLocalElmt(tmp, fieldcoeffs[j]);
    }

    WriteFld(outname1, m_fields[0], fieldcoeffs, variables);
}

void MMFSystem::PlotFieldVector(const Array<OneD, const NekDouble> &field,
                                const Array<OneD, const NekDouble> &velocity,
                                const int nstep)
{
    int nvar    = 4;
    int nq      = m_fields[0]->GetTotPoints();
    int ncoeffs = m_fields[0]->GetNcoeffs();

    std::string outname;
    outname = m_sessionName + "_Field_" +
              boost::lexical_cast<std::string>(nstep) + ".chk";

    std::vector<Array<OneD, NekDouble>> fieldcoeffs(nvar);
    for (int i = 0; i < nvar; ++i)
    {
        fieldcoeffs[i] = Array<OneD, NekDouble>(ncoeffs);
    }

    std::vector<std::string> variables(nvar);
    variables[0] = "field";
    variables[1] = "vx";
    variables[2] = "vy";
    variables[3] = "vz";

    m_fields[0]->FwdTransLocalElmt(field, fieldcoeffs[0]);

    Array<OneD, NekDouble> tmp(nq);
    for (int i = 0; i < m_spacedim; ++i)
    {
        Vmath::Vcopy(nq, &velocity[i * nq], 1, &tmp[0], 1);
        m_fields[0]->FwdTransLocalElmt(tmp, fieldcoeffs[i + 1]);
    }

    WriteFld(outname, m_fields[0], fieldcoeffs, variables);
}

void MMFSystem::PlotFieldVector(
    const Array<OneD, const NekDouble> &field,
    const Array<OneD, const Array<OneD, NekDouble>> &velocity, const int nstep)
{
    int nvar    = 4;
    int nq      = m_fields[0]->GetTotPoints();
    int ncoeffs = m_fields[0]->GetNcoeffs();

    std::string outname;
    outname = m_sessionName + "_Field_" +
              boost::lexical_cast<std::string>(nstep) + ".chk";

    std::vector<Array<OneD, NekDouble>> fieldcoeffs(nvar);
    for (int i = 0; i < nvar; ++i)
    {
        fieldcoeffs[i] = Array<OneD, NekDouble>(ncoeffs);
    }

    std::vector<std::string> variables(nvar);
    variables[0] = "field";
    variables[1] = "vx";
    variables[2] = "vy";
    variables[3] = "vz";

    m_fields[0]->FwdTransLocalElmt(field, fieldcoeffs[0]);

    Array<OneD, NekDouble> tmp(nq);
    for (int i = 0; i < nvar - 1; ++i)
    {
        Vmath::Vcopy(nq, &velocity[i][0], 1, &tmp[0], 1);
        m_fields[0]->FwdTransLocalElmt(tmp, fieldcoeffs[i + 1]);
    }

    WriteFld(outname, m_fields[0], fieldcoeffs, variables);
}

void MMFSystem::PlotJacobian(const Array<OneD, const NekDouble> &Jacobian,
                             const Array<OneD, const NekDouble> &JacGradMag)
{
    int nvar    = 2;
    int ncoeffs = m_fields[0]->GetNcoeffs();

    std::string outname;
    outname = m_sessionName + "_Jac.chk";

    std::vector<Array<OneD, NekDouble>> fieldcoeffs(nvar);
    for (int i = 0; i < nvar; ++i)
    {
        fieldcoeffs[i] = Array<OneD, NekDouble>(ncoeffs);
    }

    std::vector<std::string> variables(nvar);
    variables[0] = "Jac";
    variables[1] = "JacGradMag";

    m_fields[0]->FwdTransLocalElmt(Jacobian, fieldcoeffs[0]);
    m_fields[0]->FwdTransLocalElmt(JacGradMag, fieldcoeffs[1]);

    WriteFld(outname, m_fields[0], fieldcoeffs, variables);
}

void MMFSystem::PlotAnisotropyFiber(
    const Array<OneD, const NekDouble> &anifibre)
{
    int nvar    = 3;
    int nq      = m_fields[0]->GetTotPoints();
    int ncoeffs = m_fields[0]->GetNcoeffs();

    std::string outname1;
    outname1 = m_sessionName + "_AniFibre.chk";
    std::vector<Array<OneD, NekDouble>> fieldcoeffs(nvar);
    for (int i = 0; i < nvar; ++i)
    {
        fieldcoeffs[i] = Array<OneD, NekDouble>(ncoeffs);
    }

    std::vector<std::string> variables(nvar);
    variables[0] = "fx";
    variables[1] = "fy";
    variables[2] = "fz";

    Array<OneD, NekDouble> tmp(nq);
    for (int i = 0; i < m_spacedim; ++i)
    {
        Vmath::Vcopy(nq, &anifibre[i * nq], 1, &tmp[0], 1);
        m_fields[0]->FwdTransLocalElmt(tmp, fieldcoeffs[i]);
    }

    WriteFld(outname1, m_fields[0], fieldcoeffs, variables);
}

void MMFSystem::PlotDivMF(const Array<OneD, const NekDouble> &DivMF1,
                          const Array<OneD, const NekDouble> &DivMF2,
                          const Array<OneD, const NekDouble> &DivMF3,
                          const int nstep)
{
    int nvar    = 3;
    int nq      = m_fields[0]->GetTotPoints();
    int ncoeffs = m_fields[0]->GetNcoeffs();

    std::string outname;
    outname = m_sessionName + "_DivMF_" +
              boost::lexical_cast<std::string>(nstep) + ".chk";

    std::vector<Array<OneD, NekDouble>> fieldcoeffs(nvar);
    for (int i = 0; i < nvar; ++i)
    {
        fieldcoeffs[i] = Array<OneD, NekDouble>(ncoeffs);
    }

    std::vector<std::string> variables(nvar);
    variables[0] = "DivMF1";
    variables[1] = "DivMF2";
    variables[2] = "DivMF3";

    Array<OneD, NekDouble> tmp(nq);

    Vmath::Vcopy(nq, &DivMF1[0], 1, &tmp[0], 1);
    m_fields[0]->FwdTransLocalElmt(tmp, fieldcoeffs[0]);

    Vmath::Vcopy(nq, &DivMF2[0], 1, &tmp[0], 1);
    m_fields[0]->FwdTransLocalElmt(tmp, fieldcoeffs[1]);

    Vmath::Vcopy(nq, &DivMF3[0], 1, &tmp[0], 1);
    m_fields[0]->FwdTransLocalElmt(tmp, fieldcoeffs[2]);

    WriteFld(outname, m_fields[0], fieldcoeffs, variables);
}

void MMFSystem::PlotConnectionForm(
    const Array<OneD, const Array<OneD, Array<OneD, NekDouble>>>
        &connectionform,
    const int nstep)
{
    int nvar    = 9;
    int nq      = m_fields[0]->GetTotPoints();
    int ncoeffs = m_fields[0]->GetNcoeffs();

    std::string outname;
    outname = m_sessionName + "_Connection_" +
              boost::lexical_cast<std::string>(nstep) + ".chk";

    std::vector<Array<OneD, NekDouble>> fieldcoeffs(nvar);
    for (int i = 0; i < nvar; ++i)
    {
        fieldcoeffs[i] = Array<OneD, NekDouble>(ncoeffs);
    }

    std::vector<std::string> variables(nvar);
    variables[0] = "w121";
    variables[1] = "w122";
    variables[2] = "w123";
    variables[3] = "w131";
    variables[4] = "w132";
    variables[5] = "w133";
    variables[6] = "w231";
    variables[7] = "w232";
    variables[8] = "w233";

    Array<OneD, NekDouble> tmp(nq);
    int i, j, index;
    for (i = 0; i < m_spacedim; ++i)
    {
        for (j = 0; j < m_spacedim; ++j)
        {
            index = m_spacedim * i + j;
            Vmath::Vcopy(nq, &connectionform[i][j][0], 1, &tmp[0], 1);
            m_fields[0]->FwdTransLocalElmt(tmp, fieldcoeffs[index]);
        }
    }

    WriteFld(outname, m_fields[0], fieldcoeffs, variables);
}

void MMFSystem::Plot2DConnectionCurvature(
    const Array<OneD, const Array<OneD, Array<OneD, NekDouble>>>
        &connectionform,
    const Array<OneD, const Array<OneD, NekDouble>> &curvatureform,
    const int nstep)
{
    Plot2DConnectionForm(connectionform, nstep);
    Plot2DMeanGaussCurv(connectionform, nstep);
    Plot2DRiemanTensor(curvatureform, nstep);
}

void MMFSystem::Plot2DConnectionForm(
    const Array<OneD, const Array<OneD, Array<OneD, NekDouble>>>
        &connectionform,
    const int nstep)
{
    int nvar    = 6;
    int nq      = m_fields[0]->GetTotPoints();
    int ncoeffs = m_fields[0]->GetNcoeffs();

    std::string outname;
    outname = m_sessionName + "_Connection_" +
              boost::lexical_cast<std::string>(nstep) + ".chk";

    std::vector<Array<OneD, NekDouble>> fieldcoeffs(nvar);
    for (int i = 0; i < nvar; ++i)
    {
        fieldcoeffs[i] = Array<OneD, NekDouble>(ncoeffs);
    }

    std::vector<std::string> variables(nvar);
    variables[0] = "w211";
    variables[1] = "w212";
    variables[2] = "w311";
    variables[3] = "w312";
    variables[4] = "w321";
    variables[5] = "w322";

    Array<OneD, NekDouble> tmp(nq);
    int i, j, index;
    for (i = 0; i < m_mfdim; ++i)
    {
        for (j = 0; j < m_shapedim; ++j)
        {
            index = m_shapedim * i + j;
            Vmath::Vcopy(nq, &connectionform[i][j][0], 1, &tmp[0], 1);
            m_fields[0]->FwdTransLocalElmt(tmp, fieldcoeffs[index]);
        }
    }

    WriteFld(outname, m_fields[0], fieldcoeffs, variables);
}

void MMFSystem::Plot2DMeanGaussCurv(
    const Array<OneD, const Array<OneD, Array<OneD, NekDouble>>>
        &connectionform,
    const int nstep)
{
    int nvar    = 2;
    int nq      = m_fields[0]->GetTotPoints();
    int ncoeffs = m_fields[0]->GetNcoeffs();

    std::string outname;
    outname = m_sessionName + "_MeanGauss_" +
              boost::lexical_cast<std::string>(nstep) + ".chk";

    std::vector<Array<OneD, NekDouble>> fieldcoeffs(nvar);
    for (int i = 0; i < nvar; ++i)
    {
        fieldcoeffs[i] = Array<OneD, NekDouble>(ncoeffs);
    }

    std::vector<std::string> variables(nvar);
    variables[0] = "Mean Curvature";
    variables[1] = "Gauss Curvature";

    // conection[0][0] = "w211";
    // conection[0][1] = "w212";
    // conection[1][0] = "w311";
    // conection[1][1] = "w312";
    // conection[2][0] = "w321";
    // conection[2][1] = "w322";

    // H = 0.5 ( w311 + w322)
    Array<OneD, NekDouble> Mean(nq);
    Vmath::Vadd(nq, &connectionform[1][0][0], 1, &connectionform[2][1][0], 1,
                &Mean[0], 1);
    Vmath::Smul(nq, 0.5, &Mean[0], 1, &Mean[0], 1);
    m_fields[0]->FwdTransLocalElmt(Mean, fieldcoeffs[0]);

    // K = w311 * w322 - w312 * w321
    Array<OneD, NekDouble> Gauss(nq);
    Array<OneD, NekDouble> tmp(nq);
    Vmath::Vmul(nq, &connectionform[1][0][0], 1, &connectionform[2][1][0], 1,
                &Gauss[0], 1);
    Vmath::Vmul(nq, &connectionform[1][1][0], 1, &connectionform[2][0][0], 1,
                &tmp[0], 1);
    Vmath::Neg(nq, tmp, 1);
    Vmath::Vadd(nq, tmp, 1, Gauss, 1, Gauss, 1);

    m_fields[0]->FwdTransLocalElmt(Gauss, fieldcoeffs[1]);

    WriteFld(outname, m_fields[0], fieldcoeffs, variables);
}

// PlotConnectionError(Exactw121);
void MMFSystem::PlotConnectionError(const Array<OneD, const NekDouble> &Werror,
                                    const int nstep)
{
    int nvar    = 1;
    int ncoeffs = m_fields[0]->GetNcoeffs();

    std::string outname;
    outname = m_sessionName + "_w122_Err_" +
              boost::lexical_cast<std::string>(nstep) + ".chk";

    std::vector<Array<OneD, NekDouble>> fieldcoeffs(nvar);
    for (int i = 0; i < nvar; ++i)
    {
        fieldcoeffs[i] = Array<OneD, NekDouble>(ncoeffs);
    }

    std::vector<std::string> variables(nvar);
    variables[0] = "w122_error";

    m_fields[0]->FwdTransLocalElmt(Werror, fieldcoeffs[0]);

    WriteFld(outname, m_fields[0], fieldcoeffs, variables);
}


void MMFSystem::Plot2DRiemanTensor(
    const Array<OneD, const Array<OneD, NekDouble>> &curvatureform,
    const int nstep)
{
    // int nvar = m_mfdim * m_mfdim * m_mfdim;
    int nvar    = 3;
    int nq      = m_fields[0]->GetTotPoints();
    int ncoeffs = m_fields[0]->GetNcoeffs();

    std::string outname;
    outname = m_sessionName + "_Riemann_" +
              boost::lexical_cast<std::string>(nstep) + ".chk";

    std::vector<Array<OneD, NekDouble>> fieldcoeffs(nvar);
    for (int i = 0; i < nvar; ++i)
    {
        fieldcoeffs[i] = Array<OneD, NekDouble>(ncoeffs);
    }

    std::vector<std::string> variables(nvar);
    variables[0] = "R2112";
    variables[1] = "R3112";
    variables[2] = "R3212";

    Array<OneD, NekDouble> tmp(nq);
    for (int i = 0; i < m_mfdim; ++i)
    {
        Vmath::Vcopy(nq, &curvatureform[i][0], 1, &tmp[0], 1);
        m_fields[0]->FwdTransLocalElmt(tmp, fieldcoeffs[i]);
    }

    WriteFld(outname, m_fields[0], fieldcoeffs, variables);
}
void MMFSystem::PlotCurvatureForm(
    const Array<OneD, const Array<OneD, Array<OneD, Array<OneD, NekDouble>>>>
        &curvatureform,
    const int nstep)
{
    PlotRiemannTensor(curvatureform, nstep);
    PlotRicciTensor(curvatureform, nstep);
}

void MMFSystem::PlotRiemannTensor(
    const Array<OneD, const Array<OneD, Array<OneD, Array<OneD, NekDouble>>>>
        &curvatureform,
    const int nstep)
{
    // int nvar = m_mfdim * m_mfdim * m_mfdim;
    int nvar    = 27;
    int nq      = m_fields[0]->GetTotPoints();
    int ncoeffs = m_fields[0]->GetNcoeffs();

    std::string outname;
    outname = m_sessionName + "_Riemann_" +
              boost::lexical_cast<std::string>(nstep) + ".chk";

    std::vector<Array<OneD, NekDouble>> fieldcoeffs(nvar);
    for (int i = 0; i < nvar; ++i)
    {
        fieldcoeffs[i] = Array<OneD, NekDouble>(ncoeffs);
    }

    std::vector<std::string> variables(nvar);
    variables[0] = "R2111";
    variables[1] = "R2112";
    variables[2] = "R2113";
    variables[3] = "R2121";
    variables[4] = "R2122";
    variables[5] = "R2123";
    variables[6] = "R2131";
    variables[7] = "R2132";
    variables[8] = "R2133";

    variables[9]  = "R3111";
    variables[10] = "R3112";
    variables[11] = "R3113";
    variables[12] = "R3121";
    variables[13] = "R3122";
    variables[14] = "R3123";
    variables[15] = "R3131";
    variables[16] = "R3132";
    variables[17] = "R3133";

    variables[18] = "R3211";
    variables[19] = "R3212";
    variables[20] = "R3213";
    variables[21] = "R3221";
    variables[22] = "R3222";
    variables[23] = "R3223";
    variables[24] = "R3231";
    variables[25] = "R3232";
    variables[26] = "R3233";

    int index;
    Array<OneD, NekDouble> tmp(nq);
    for (int i = 0; i < m_mfdim; ++i)
    {
        for (int j = 0; j < m_mfdim; ++j)
        {
            for (int k = 0; k < m_mfdim; ++k)
            {
                index = m_mfdim * m_mfdim * i + m_mfdim * j + k;
                Vmath::Vcopy(nq, &curvatureform[i][j][k][0], 1, &tmp[0], 1);
                m_fields[0]->FwdTransLocalElmt(tmp, fieldcoeffs[index]);
            }
        }
    }

    WriteFld(outname, m_fields[0], fieldcoeffs, variables);
}

void MMFSystem::PlotRicciTensor(
    const Array<OneD, const Array<OneD, Array<OneD, Array<OneD, NekDouble>>>>
        &curvatureform,
    const int nstep)
{
    int nvar    = 3;
    int nq      = m_fields[0]->GetTotPoints();
    int ncoeffs = m_fields[0]->GetNcoeffs();

    std::string outname;
    outname = m_sessionName + "_Ricci_" +
              boost::lexical_cast<std::string>(nstep) + ".chk";

    std::vector<Array<OneD, NekDouble>> fieldcoeffs(nvar);
    for (int i = 0; i < nvar; ++i)
    {
        fieldcoeffs[i] = Array<OneD, NekDouble>(ncoeffs);
    }

    std::vector<std::string> variables(nvar);
    variables[0] = "R11";
    variables[1] = "R12";
    variables[2] = "R22";

    int index;
    Array<OneD, NekDouble> tmp(nq);

    // R11 = R^2_{121} + R^3_{131}
    index = 0;
    Vmath::Vadd(nq, &curvatureform[0][1][0][0], 1, &curvatureform[1][2][0][0],
                1, &tmp[0], 1);
    m_fields[0]->FwdTransLocalElmt(tmp, fieldcoeffs[index]);

    // R12 = R^2_{122} + R^3_{132}
    // R12 = RootMeanSquare(Curvatureform[0][1][1], Ntot) +
    // RootMeanSquare(Curvatureform[1][2][1], Ntot);
    index = 1;
    Vmath::Vadd(nq, &curvatureform[0][1][1][0], 1, &curvatureform[1][2][1][0],
                1, &tmp[0], 1);
    m_fields[0]->FwdTransLocalElmt(tmp, fieldcoeffs[index]);

    // R22 = R^3_{232}
    // R22 = RootMeanSquare(Curvatureform[2][2][1], Ntot);
    index = 2;
    Vmath::Vcopy(nq, &curvatureform[2][2][1][0], 1, &tmp[0], 1);
    m_fields[0]->FwdTransLocalElmt(tmp, fieldcoeffs[index]);

    WriteFld(outname, m_fields[0], fieldcoeffs, variables);
}

void MMFSystem::ComputeDirectionVector(
    Array<OneD, Array<OneD, NekDouble>> &distance)
{
    int nq = m_fields[0]->GetNpoints();

    distance = Array<OneD, Array<OneD, NekDouble>>(m_mfdim);
    for (int k = 0; k < m_mfdim; ++k)
    {
        distance[k] = Array<OneD, NekDouble>(nq);
    }

    Array<OneD, NekDouble> x0(nq);
    Array<OneD, NekDouble> x1(nq);
    Array<OneD, NekDouble> x2(nq);

    m_fields[0]->GetCoords(x0, x1, x2);

    // Choose the first point as its anchor
    NekDouble xp0, xp1, xp2, rp;
    xp0 = x0[0];
    xp1 = x1[0];
    xp2 = x2[0];
    rp  = sqrt(xp0 * xp0 + xp1 * xp1 + xp2 * xp2);

    // Projectin of distance on the tangent plane: \vec{v} -
    // (\vec{v}\cdot \vec{n}) \vec{n}
    NekDouble vx, vy, vz, rad, length;
    NekDouble nx, ny, nz, magv, vdotn;
    for (int i = 0; i < nq; i++)
    {
        // length is the length of the arc length: (v1 \cdot v2)/||v1||
        // ||v2|
        rad    = sqrt(x0[i] * x0[i] + x1[i] * x1[i] + x2[i] * x2[i]);
        length = acos((x0[i] * xp0 + x1[i] * xp2 + x2[i] * xp2) / rp / rad);

        // Straight distance vector
        vx = x0[i] - xp0;
        vy = x1[i] - xp1;
        vz = x2[i] - xp2;

        // normal vector to the tangent plane
        nx = m_movingframes[2][i];
        ny = m_movingframes[2][i + nq];
        nz = m_movingframes[2][i + 2 * nq];

        // \vec{v} - (\vec{v}\cdot \vec{n}) \vec{n}
        vdotn = vx * nx + vy * ny + vz * nz;
        vx    = vx - vdotn * nx;
        vy    = vy - vdotn * ny;
        vz    = vz - vdotn * nz;

        magv = sqrt(vx * vx + vy * vy + vz * vz);
        if (magv > 0.00000001)
        {
            vx = length * vx / magv;
            vy = length * vy / magv;
            vz = length * vz / magv;
        }

        for (int k = 0; k < m_mfdim; k++)
        {
            distance[k][i] = vx * m_movingframes[k][i] +
                             vy * m_movingframes[k][i + nq] +
                             vz * m_movingframes[k][i + 2 * nq];
        }
    }

    std::cout << "Avg. dist vec = ( " << RootMeanSquare(distance[0]) << " , "
              << RootMeanSquare(distance[1]) << " , "
              << RootMeanSquare(distance[2]) << " ) " << std::endl;
}

// input: velocity ( 3 x nq).
// output: MF1 = velocity, MF2 = orthogonal to MF1, MF3 = same.
// If UseFlatMF = 1, then we do not align MF along the velocity vector
// if its magnitude is zero at least one point.
// void MMFSystem::AlignMFtoVelocity(
//     const NekDouble VelActivationTol,
//     const Array<OneD, const NekDouble> &velocity,
//     const Array<OneD, const NekDouble> &velmag,
//     const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
//     Array<OneD, Array<OneD, NekDouble>> &MF1st, Array<OneD, int> &Activated,
//     Array<OneD, NekDouble> &ActivationIntensity)
// {
//     int nq    = m_fields[0]->GetNpoints();
//     Activated = Array<OneD, int>(nq, 0);

//     // MFmag = \sqrt{ e1x*e1x + e1y*e1y + e1z*e1z}
//     Array<OneD, Array<OneD, NekDouble>> MForigmag(m_mfdim);
//     for (int i = 0; i < m_mfdim; ++i)
//     {
//         MForigmag[i] = Array<OneD, NekDouble>(nq, 0.0);
//         for (int k = 0; k < m_spacedim; ++k)
//         {
//             Vmath::Vvtvp(nq, &movingframes[i][k * nq], 1,
//                          &movingframes[i][k * nq], 1, &MForigmag[i][0], 1,
//                          &MForigmag[i][0], 1);
//         }

//         Vmath::Vsqrt(nq, MForigmag[i], 1, MForigmag[i], 1);
//     }

//     for (int k = 0; k < m_mfdim; ++k)
//     {
//         Vmath::Vcopy(3 * nq, &movingframes[k][0], 1, &MF1st[k][0], 1);
//     }

//     // Let the element is activated when the velocity at all grid points
//     // is not trivial Don't care about non-differentiability here.
//     m_fields[0]->ElementWiseAlignment(velmag, VelActivationTol, Activated);

//     ActivationIntensity = Array<OneD, NekDouble>(nq, 0.0);

//     NekDouble MF1x, MF1y, MF1z, MF2x, MF2y, MF2z, MF3x, MF3y, MF3z;
//     NekDouble MFmag1, MFmag2;
//     for (int i = 0; i < nq; i++)
//     {
//         // if (velmag > TolActivated)
//         if (Activated[i])
//         {
//             ActivationIntensity[i] = velmag[i];

//             // if Activated, align the moving frames along the velocity
//             // vector

//             MF1x = velocity[i];
//             MF1y = velocity[i + nq];
//             MF1z = velocity[i + 2 * nq];

//             MF3x = MF1st[2][i];
//             MF3y = MF1st[2][i + nq];
//             MF3z = MF1st[2][i + 2 * nq];

//             MF2x = MF3y * MF1z - MF3z * MF1y;
//             MF2y = MF3z * MF1x - MF3x * MF1z;
//             MF2z = MF3x * MF1y - MF3y * MF1x;

//             MFmag1           = sqrt(MF1x * MF1x + MF1y * MF1y + MF1z * MF1z);
//             MF1st[0][i]      = MForigmag[0][i] * MF1x / MFmag1;
//             MF1st[0][i + nq] = MForigmag[0][i] * MF1y / MFmag1;
//             MF1st[0][i + 2 * nq] = MForigmag[0][i] * MF1z / MFmag1;

//             MFmag2           = sqrt(MF2x * MF2x + MF2y * MF2y + MF2z * MF2z);
//             MF1st[1][i]      = MForigmag[1][i] * MF2x / MFmag2;
//             MF1st[1][i + nq] = MForigmag[1][i] * MF2y / MFmag2;
//             MF1st[1][i + 2 * nq] = MForigmag[1][i] * MF2z / MFmag2;
//         }
//     }
// }

// input: velocity ( 3 x nq).
// output: MF1 = velocity, MF2 = orthogonal to MF1, MF3 = same.
// If UseFlatMF = 1, then we do not align MF along the velocity vector
// if its magnitude is zero at least one point.
void MMFSystem::AlignMFtoVelocity(
    const Array<OneD, const int> &Activated,
    const Array<OneD, const NekDouble> &velocity,
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    Array<OneD, Array<OneD, NekDouble>> &MF1st)
{
    int nq = m_fields[0]->GetNpoints();

    for (int k = 0; k < m_mfdim; ++k)
    {
        Vmath::Vcopy(3 * nq, &movingframes[k][0], 1, &MF1st[k][0], 1);
    }

    NekDouble MF1x, MF1y, MF1z, MF2x, MF2y, MF2z, MF3x, MF3y, MF3z;
    NekDouble MFmag1, MFmag2;
    for (int i = 0; i < nq; i++)
    {
        // if (velmag > TolActivated)
        if (Activated[i])
        {
            // if Activated, align the moving frames along the velocity
            // vector
            MF1x = velocity[i];
            MF1y = velocity[i + nq];
            MF1z = velocity[i + 2 * nq];

            MF3x = MF1st[2][i];
            MF3y = MF1st[2][i + nq];
            MF3z = MF1st[2][i + 2 * nq];

            MF2x = MF3y * MF1z - MF3z * MF1y;
            MF2y = MF3z * MF1x - MF3x * MF1z;
            MF2z = MF3x * MF1y - MF3y * MF1x;

            MFmag1           = sqrt(MF1x * MF1x + MF1y * MF1y + MF1z * MF1z);
            MF1st[0][i]      = MF1x / MFmag1;
            MF1st[0][i + nq] = MF1y / MFmag1;
            MF1st[0][i + 2 * nq] = MF1z / MFmag1;

            MFmag2           = sqrt(MF2x * MF2x + MF2y * MF2y + MF2z * MF2z);
            MF1st[1][i]      = MF2x / MFmag2;
            MF1st[1][i + nq] = MF2y / MFmag2;
            MF1st[1][i + 2 * nq] = MF2z / MFmag2;
        }
    }
}

void MMFSystem::ActivateRegion(
    const Array<OneD, const int> &ZoneActivation,
    const Array<OneD, const Array<OneD, NekDouble>> &MF1st,
    Array<OneD, Array<OneD, NekDouble>> &MF1stnew, Array<OneD, int> &Activated)
{
    int nq = GetTotPoints();

    // Compare outarray and outarraynew one by one. If the difference is
    // more than tolerance, set ActNow = 0 again
    int indexj, indexk;
    for (int i = 0; i < m_fields[0]->GetExpSize(); ++i)
    {
        for (int j = 0; j < m_fields[0]->GetTotPoints(i); ++j)
        {
            indexj = m_fields[0]->GetPhys_Offset(i) + j;
            // When \| \vec{E}_0 - \vec{E}_{new} \| > Tol, where Tol =
            // ANF*dt to be independent of timestep
            if (ZoneActivation[indexj] == 0)
            {
                for (int k = 0; k < m_fields[0]->GetTotPoints(i); ++k)
                {
                    indexk            = m_fields[0]->GetPhys_Offset(i) + k;
                    Activated[indexk] = 0;

                    // Retrieve the original moving frames
                    for (int m = 0; m < m_mfdim; m++)
                    {
                        MF1stnew[m][indexk]      = MF1st[m][indexk];
                        MF1stnew[m][indexk + nq] = MF1st[m][indexk + nq];
                        MF1stnew[m][indexk + 2 * nq] =
                            MF1st[m][indexk + 2 * nq];
                    }
                }
                break;
            }
        }
    }
}

int MMFSystem::UpdatebyLinearTimeIntegration(
    const Array<OneD, const int> &Activated,
    const Array<OneD, const Array<OneD, NekDouble>> &MF1st,
    Array<OneD, int> &ActivatedHistory,
    Array<OneD, Array<OneD, NekDouble>> &MF1stAlignedInt)
{
    int cnt = 0;
    int nq  = GetNpoints();

    NekDouble ex, ey, ez;
    NekDouble newex, newey, newez;
    for (int i = 0; i < nq; ++i)
    {
        if (Activated[i])
        {
            // Averaged Integration
            NekDouble AC = ActivatedHistory[i];
            for (int j = 0; j < m_mfdim; ++j)
            {
                ex = MF1st[j][i];
                ey = MF1st[j][i + nq];
                ez = MF1st[j][i + 2 * nq];

                newex = MF1stAlignedInt[j][i];
                newey = MF1stAlignedInt[j][i + nq];
                newez = MF1stAlignedInt[j][i + 2 * nq];

                MF1stAlignedInt[j][i]          = (AC * newex + ex) / (AC + 1.0);
                MF1stAlignedInt[j][i + nq]     = (AC * newey + ey) / (AC + 1.0);
                MF1stAlignedInt[j][i + 2 * nq] = (AC * newez + ez) / (AC + 1.0);
            }

            ActivatedHistory[i] += 1;
        }

        if (ActivatedHistory[i] > 0)
        {
            cnt++;
        }
    }

    return cnt;
}

void MMFSystem::UpdateMF1st(
    const Array<OneD, const int> &Activated,
    const Array<OneD, const Array<OneD, NekDouble>> &MF1st,
    const Array<OneD, const NekDouble> &velmag,
    Array<OneD, NekDouble> &VelmagHistory,
    Array<OneD, Array<OneD, NekDouble>> &MF1stAligned,
    Array<OneD, int> &ActivatedHistory)
{
    int nq = GetNpoints();

    Array<OneD, int> outarray(nq, 0);

    NekDouble ex, ey, ez;
    NekDouble newex, newey, newez;
    for (int i = 0; i < nq; ++i)
    {
        if (Activated[i])
        {
            // Weighted Averaged Integration
            NekDouble Wsum   = VelmagHistory[i];
            NekDouble Weight = velmag[i];
            for (int j = 0; j < m_mfdim; ++j)
            {
                ex = MF1st[j][i];
                ey = MF1st[j][i + nq];
                ez = MF1st[j][i + 2 * nq];

                newex = MF1stAligned[j][i];
                newey = MF1stAligned[j][i + nq];
                newez = MF1stAligned[j][i + 2 * nq];

                MF1stAligned[j][i] =
                    (Wsum * newex + Weight * ex) / (Wsum + Weight);
                MF1stAligned[j][i + nq] =
                    (Wsum * newey + Weight * ey) / (Wsum + Weight);
                MF1stAligned[j][i + 2 * nq] =
                    (Wsum * newez + Weight * ez) / (Wsum + Weight);
            }

            VelmagHistory[i] += velmag[i];

            // Activated History update
            ActivatedHistory[i] = 1;
        }
    }
}

// Should be changed as pointwise update
int MMFSystem::NewValueReplacerElementWise(
    const EvalType EType, const Array<OneD, const NekDouble> &Intensity,
    const Array<OneD, const int> &Activated,
    const Array<OneD, const Array<OneD, NekDouble>> &MF1st,
    Array<OneD, NekDouble> &IntensityHistory,
    Array<OneD, int> &ActivatedHistory,
    Array<OneD, Array<OneD, NekDouble>> &MF1stHistory)
{
    int nq = GetTotPoints();

    Array<OneD, int> NewActivated(nq);

    Vmath::Vcopy(nq, Activated, 1, NewActivated, 1);

    int index, Voting;
    bool Voting_Yes;
    for (int i = 0; i < m_fields[0]->GetExpSize(); ++i)
    {
        Voting = 0;
        for (int j = 0; j < m_fields[0]->GetTotPoints(i); ++j)
        {
            index = m_fields[0]->GetPhys_Offset(i) + j;
            if (Activated[index])
            {
                Voting_Yes = false;
                // Stronger Replacer
                if (EType == eStronger)
                {
                    if (Intensity[index] > IntensityHistory[index])
                    {
                        Voting_Yes = true;
                    }
                }

                // Stronger Replacer
                else if (EType == eWeaker)
                {
                    if (ActivatedHistory[i] < 0.000000001)
                    {
                        Voting_Yes = true;
                    }

                    else if (Intensity[index] < IntensityHistory[index])
                    {
                        Voting_Yes = true;
                    }
                }

                if (Voting_Yes)
                {
                    Voting++;
                }
            }
        }

        // Activated if Intensity is larger than Intensity History at
        // all points
        if (Voting < m_fields[0]->GetTotPoints(i))
        {
            for (int j = 0; j < m_fields[0]->GetTotPoints(i); ++j)
            {
                index               = m_fields[0]->GetPhys_Offset(i) + j;
                NewActivated[index] = 0;
            }
        }
    } // end of i-loop

    int cnt = 0;
    // Update MF1st if ActivatedReturn is on
    for (int i = 0; i < nq; ++i)
    {
        // Modify Moving frames
        if (NewActivated[i])
        {
            // Strengthbased Integration
            for (int j = 0; j < m_mfdim; ++j)
            {
                MF1stHistory[j][i]          = MF1st[j][i];
                MF1stHistory[j][i + nq]     = MF1st[j][i + nq];
                MF1stHistory[j][i + 2 * nq] = MF1st[j][i + 2 * nq];
            }

            IntensityHistory[i] = Intensity[i];
            ActivatedHistory[i] = 1.0;
        }

        if (ActivatedHistory[i] > 0)
        {
            cnt++;
        }
    }

    return cnt;
} // namespace SolverUtils

// Should be changed as pointwise update
int MMFSystem::NewValueReplacerPointWise(
    const EvalType EType, const Array<OneD, const NekDouble> &Intensity,
    const Array<OneD, const int> &Activated,
    const Array<OneD, const Array<OneD, NekDouble>> &MF1st,
    Array<OneD, NekDouble> &IntensityHistory,
    Array<OneD, NekDouble> &ActivatedHistory,
    Array<OneD, Array<OneD, NekDouble>> &MF1stHistory)
{
    int nq = GetTotPoints();

    Array<OneD, int> NewActivated(nq);

    Vmath::Vcopy(nq, Activated, 1, NewActivated, 1);

    int index;
    for (int i = 0; i < m_fields[0]->GetExpSize(); ++i)
    {
        for (int j = 0; j < m_fields[0]->GetTotPoints(i); ++j)
        {
            index = m_fields[0]->GetPhys_Offset(i) + j;
            if (Activated[index])
            {
                // Stronger Replacer
                if (EType == eStronger)
                {
                    if (Intensity[index] < IntensityHistory[index])
                    {
                        NewActivated[index] = 0;
                    }
                }

                // Stronger Replacer
                else if (EType == eWeaker)
                {
                    if ((ActivatedHistory[i] > 0.000000001) &&
                        (Intensity[index] > IntensityHistory[index]))
                    {
                        NewActivated[index] = 0;
                    }
                }
            }
        }
    } // end of i-loop

    int cnt = 0;
    // Update MF1st if ActivatedReturn is on
    for (int i = 0; i < nq; ++i)
    {
        // Modify Moving frames
        if (NewActivated[i])
        {
            // Strengthbased Integration
            for (int j = 0; j < m_mfdim; ++j)
            {
                MF1stHistory[j][i]          = MF1st[j][i];
                MF1stHistory[j][i + nq]     = MF1st[j][i + nq];
                MF1stHistory[j][i + 2 * nq] = MF1st[j][i + 2 * nq];
            }

            IntensityHistory[i] = Intensity[i];
            ActivatedHistory[i] = 1.0;
        }

        if (ActivatedHistory[i] > 0)
        {
            cnt++;
        }
    }

    return cnt;
} // namespace SolverUtils

int MMFSystem::WeakerValueReplacer(
    const Array<OneD, const NekDouble> &Intensity,
    const Array<OneD, const int> &Activated,
    const Array<OneD, const Array<OneD, NekDouble>> &MF1st,
    Array<OneD, NekDouble> &IntensityHistory,
    Array<OneD, NekDouble> &ActivatedHistory,
    Array<OneD, Array<OneD, NekDouble>> &MF1stHistory)
{
    int nq = GetTotPoints();

    Array<OneD, int> NewActivated(nq);

    Vmath::Vcopy(nq, Activated, 1, NewActivated, 1);

    int index, Voting;
    for (int i = 0; i < m_fields[0]->GetExpSize(); ++i)
    {
        Voting = 0;
        for (int j = 0; j < m_fields[0]->GetTotPoints(i); ++j)
        {
            index = m_fields[0]->GetPhys_Offset(i) + j;
            if (Activated[index])
            {
                // If it is the first time, then, we count it
                if (ActivatedHistory[i] < 0.000000001)
                {
                    Voting++;
                }

                else if (Intensity[index] < IntensityHistory[index])
                {
                    Voting++;
                }
            }
        }

        // Activated if Intensity is larger than Intensity History at
        // all points
        if (Voting < m_fields[0]->GetTotPoints(i))
        {
            for (int j = 0; j < m_fields[0]->GetTotPoints(i); ++j)
            {
                index               = m_fields[0]->GetPhys_Offset(i) + j;
                NewActivated[index] = 0;
            }
        }
    } // end of i-loop

    int cnt = 0;
    // Update MF1st if ActivatedReturn is on
    for (int i = 0; i < nq; ++i)
    {
        // Modify Moving frames
        if (NewActivated[i])
        {
            // Strengthbased Integration
            for (int j = 0; j < m_mfdim; ++j)
            {
                MF1stHistory[j][i]          = MF1st[j][i];
                MF1stHistory[j][i + nq]     = MF1st[j][i + nq];
                MF1stHistory[j][i + 2 * nq] = MF1st[j][i + 2 * nq];
            }

            IntensityHistory[i] = Intensity[i];
            ActivatedHistory[i] = NewActivated[i];
        }

        if (ActivatedHistory[i] > 0)
        {
            cnt++;
        }
    }

    return cnt;
} // namespace SolverUtils

void MMFSystem::ValidateNewFrames(
    const NekDouble AdaptNewFramesToldt,
    const Array<OneD, const NekDouble> &vecdiff,
    const Array<OneD, const Array<OneD, NekDouble>> &MF1st,
    Array<OneD, Array<OneD, NekDouble>> &MF1stnew, Array<OneD, int> &Activated)
{
    int nq = GetTotPoints();

    // Compare outarray and outarraynew one by one. If the difference is
    // more than tolerance, set ActNow = 0 again
    int indexj, indexk;
    for (int i = 0; i < m_fields[0]->GetExpSize(); ++i)
    {
        for (int j = 0; j < m_fields[0]->GetTotPoints(i); ++j)
        {
            indexj = m_fields[0]->GetPhys_Offset(i) + j;

            // When \| \vec{E}_0 - \vec{E}_{new} \| > Tol, where Tol =
            // ANF*dt to be independent of timestep If vecdiff < 0 is
            // the case when dE/dt or dH/dt is too small
            if (vecdiff[indexj] > AdaptNewFramesToldt)
            {
                for (int k = 0; k < m_fields[0]->GetTotPoints(i); ++k)
                {
                    indexk            = m_fields[0]->GetPhys_Offset(i) + k;
                    Activated[indexk] = 0;

                    // Retrieve the original moving frames
                    for (int m = 0; m < m_mfdim; m++)
                    {
                        MF1stnew[m][indexk]      = MF1st[m][indexk];
                        MF1stnew[m][indexk + nq] = MF1st[m][indexk + nq];
                        MF1stnew[m][indexk + 2 * nq] =
                            MF1st[m][indexk + 2 * nq];
                    }
                }
                break;
            }
        }
    }
}

void MMFSystem::Test2DConnection1formPolar(
    const NekDouble InitPx, const NekDouble InitPy, const NekDouble InitPz,
    const Array<OneD, const int> &Activation,
    const Array<OneD, const Array<OneD, NekDouble>> &MF1st,
    const Array<OneD, const Array<OneD, Array<OneD, NekDouble>>> &Connection)
{
    int cnt = 0;
    int nq  = GetTotPoints();

    Array<OneD, NekDouble> x(nq);
    Array<OneD, NekDouble> y(nq);
    Array<OneD, NekDouble> z(nq);

    m_fields[0]->GetCoords(x, y, z);

    NekDouble Tol = 0.00001;

    // ELF on the plane
    NekDouble diffx, diffy, diffz, rad;
    NekDouble xp, yp, zp;
    NekDouble dx = 0.0, dy = 0.0;

    NekDouble MF1sterr = 0.0;
    NekDouble tmp, w211, w212;
    NekDouble Exactw211 = 0.0, Exactw212 = 0.0;
    NekDouble w211err = 0.0, w212err = 0.0;

    for (int i = 0; i < nq; i++)
    {
        xp = x[i];
        yp = y[i];
        zp = z[i];

        w211 = Connection[0][0][i];
        w212 = Connection[0][1][i];

        diffx = xp - InitPx;
        diffy = yp - InitPy;
        diffz = zp - InitPz;

        rad = sqrt(diffx * diffx + diffy * diffy + diffz * diffz);

        if ((Activation[i] > 0) && (rad > Tol))
        {
            Exactw211 = 0.0;
            Exactw212 = -1.0 / rad;

            dx = MF1st[0][i] - xp / sqrt(xp * xp + yp * yp);
            dy = MF1st[0][i + nq] - yp / sqrt(xp * xp + yp * yp);

            MF1sterr += dx * dx + dy * dy;

            tmp = fabs(w211 - Exactw211);
            w211err += tmp * tmp;

            tmp = fabs(w212 - Exactw212);
            w212err += tmp * tmp;

            cnt++;
        }
    }

    if (cnt)
    {
        MF1sterr = sqrt(MF1sterr / (cnt * 2));

        w211err = sqrt(w211err / cnt);
        w212err = sqrt(w212err / cnt);

        std::cout << "For " << cnt << " / " << nq << " ( " << 100.0 * cnt / nq
                  << " % ) points, Error of MF1st = " << MF1sterr << std::endl;
        std::cout << "w211 err = " << w211err << ", w212 err = " << w212err
                  << std::endl;
    }
}

void MMFSystem::Test2DConnection1formSphere(
    const Array<OneD, const int> &Activation,
    const Array<OneD, const Array<OneD, NekDouble>> &MF1st,
    const Array<OneD, const Array<OneD, Array<OneD, NekDouble>>> &Connection)
{
    int nq = GetTotPoints();

    Array<OneD, NekDouble> x(nq);
    Array<OneD, NekDouble> y(nq);
    Array<OneD, NekDouble> z(nq);

    m_fields[0]->GetCoords(x, y, z);

    NekDouble Tol = 0.00001;

    NekDouble rad;

    int cnt = 0;
    NekDouble xp, yp, zp, sin_varphi, cos_varphi, sin_theta, cos_theta;
    NekDouble dx, dy, dz, theta_x, theta_y, theta_z;

    NekDouble tmp, w212, w311, w322;
    NekDouble Exactw212, Exactw311, Exactw322;
    NekDouble w212err = 0.0, w311err = 0.0, w322err = 0.0;
    NekDouble MF1sterr = 0.0;

    for (int i = 0; i < nq; i++)
    {
        w212 = Connection[0][1][i];
        w311 = Connection[1][0][i];
        w322 = Connection[2][1][i];

        xp = x[i];
        yp = y[i];
        zp = z[i];

        rad = sqrt(xp * xp + yp * yp + zp * zp);

        CartesianToNewSpherical(xp, yp, zp, sin_varphi, cos_varphi, sin_theta,
                                cos_theta);

        theta_x = cos_theta * cos_varphi;
        theta_y = cos_theta * sin_varphi;
        theta_z = -sin_theta;

        // Only for the z=20 initiation case
        if ((Activation[i] > 0) && (fabs(zp) < rad) && fabs(sin_theta) > Tol)
        {
            // omega_{\theta \phi} (e_{\phi}) = \cos \theta / \sin \theta
            Exactw212 = (-1.0 / rad) * (cos_theta / sin_theta);
            Exactw311 = 1.0 / rad;
            Exactw322 = 1.0 / rad;

            tmp = fabs(w212 - Exactw212);
            w212err += tmp * tmp;

            tmp = fabs(w311 - Exactw311);
            w311err += tmp * tmp;

            tmp = fabs(w322 - Exactw322);
            w322err += tmp * tmp;

            dx = MF1st[0][i] - theta_x;
            dy = MF1st[0][i + nq] - theta_y;
            dz = MF1st[0][i + 2 * nq] - theta_z;

            MF1sterr += dx * dx + dy * dy + dz * dz;

            cnt++;
        }
    }

    if (cnt > 0)
    {
        MF1sterr = sqrt(MF1sterr / (cnt * 3));

        w212err = sqrt(w212err / cnt);
        w311err = sqrt(w311err / cnt);
        w322err = sqrt(w322err / cnt);

        std::cout << "Test2DConnection1formSphere ========================="
                  << std::endl;
        std::cout << "For " << cnt << " / " << nq << " ( " << 100.0 * cnt / nq
                  << " % ) points, Error of MF1st = " << MF1sterr << std::endl;
        std::cout << "Error: w212 = " << w212err << ", w311  = " << w311err
                  << ", w322 = " << w322err << std::endl;
    }
}

void MMFSystem::Test2DConnection1formPseudosphere(
    const Array<OneD, const int> &Activation,
    const Array<OneD, const Array<OneD, NekDouble>> &MF1st,
    const Array<OneD, const Array<OneD, Array<OneD, NekDouble>>> &Connection)
{
    int nq = GetTotPoints();

    Array<OneD, NekDouble> x(nq);
    Array<OneD, NekDouble> y(nq);
    Array<OneD, NekDouble> z(nq);

    m_fields[0]->GetCoords(x, y, z);

    // NekDouble Tol = 0.00001;
    // NekDouble rad;

    int cnt = 0;
    NekDouble xp, yp, zp, sin_varphi, cos_varphi, theta, sech_theta, tanh_theta;
    NekDouble dx, dy, dz, theta_x, theta_y, theta_z;

    NekDouble tmp, w211, w212, w311, w322;
    NekDouble Exactw211, Exactw212, Exactw311, Exactw322;
    NekDouble w211err = 0.0, w212err = 0.0, w311err = 0.0, w322err = 0.0;
    NekDouble MF1sterr = 0.0;

    for (int i = 0; i < nq; i++)
    {
        w211 = Connection[0][0][i];
        w212 = Connection[0][1][i];
        w311 = Connection[1][0][i];
        w322 = Connection[2][1][i];

        xp = x[i];
        yp = y[i];
        zp = z[i];

        CartesianToPseudospherical(xp, yp, zp, sin_varphi, cos_varphi, theta,
                                   sech_theta, tanh_theta);

        theta_x = -sech_theta * cos_varphi;
        theta_y = -sech_theta * sin_varphi;
        theta_z = tanh_theta;

        // Only for the z=20 initiation case
        if (Activation[i] > 0)
        {
            Exactw211 = 0.0;
            Exactw212 = -1.0;
            Exactw311 = -1.0 *
                        (-tanh_theta * tanh_theta + sech_theta * sech_theta) *
                        sech_theta / tanh_theta;
            Exactw322 = -1.0;

            tmp = fabs(w211 - Exactw211);
            w211err += tmp * tmp;

            tmp = fabs(w212 - Exactw212);
            w212err += tmp * tmp;

            tmp = fabs(w311 - Exactw311);
            w311err += tmp * tmp;

            tmp = fabs(w322 - Exactw322);
            w322err += tmp * tmp;

            dx = MF1st[0][i] - theta_x;
            dy = MF1st[0][i + nq] - theta_y;
            dz = MF1st[0][i + 2 * nq] - theta_z;

            MF1sterr += dx * dx + dy * dy + dz * dz;

            cnt++;
        }
    }

    if (cnt > 0)
    {
        MF1sterr = sqrt(MF1sterr / (cnt * 3));

        w211err = sqrt(w211err / cnt);
        w212err = sqrt(w212err / cnt);
        w311err = sqrt(w311err / cnt);
        w322err = sqrt(w322err / cnt);

        std::cout << "For " << cnt << " / " << nq << " ( " << 100.0 * cnt / nq
                  << " % ) points, Error of MF1st = " << MF1sterr << std::endl;
        std::cout << "Error: w211 = " << w211err << ", w212 = " << w212err
                  << ", w311  = " << w311err << ", w322 = " << w322err
                  << std::endl;
    }
}

void MMFSystem::Test2DConnectionCurvature(
    const NekDouble InitPx, const NekDouble InitPy, const NekDouble InitPz,
    const Array<OneD, const int> &Activation,
    const Array<OneD, const Array<OneD, NekDouble>> &MF1st,
    const Array<OneD, const Array<OneD, Array<OneD, NekDouble>>> &Connection,
    const Array<OneD, const Array<OneD, NekDouble>> &Curvature)
{
    int nq = GetTotPoints();

    Array<OneD, NekDouble> x(nq);
    Array<OneD, NekDouble> y(nq);
    Array<OneD, NekDouble> z(nq);

    m_fields[0]->GetCoords(x, y, z);

    NekDouble Tol = 0.00001;

    // ELF on the Sphere
    if (m_surfaceType == SolverUtils::eSphere)
    {
        // Array<OneD, NekDouble> Exactw212(nq, 0.0);
        // Array<OneD, NekDouble> Exactw311(nq, 0.0);
        // Array<OneD, NekDouble> Exactw322(nq, 0.0);
        // Array<OneD, NekDouble> MF1sterr(nq, 0.0);
        // Array<OneD, NekDouble> ExactR2112(nq, 0.0);

        NekDouble rad; //, radPole;
        // radPole = sqrt(InitPx * InitPx + InitPy * InitPy + InitPz *
        // InitPz);

        int cnt = 0;
        NekDouble xp, yp, zp, sin_varphi, cos_varphi, sin_theta, cos_theta;
        NekDouble dx, dy, dz, theta_x, theta_y, theta_z;

        NekDouble tmp, w212, w311, w322, R2112;
        NekDouble Exactw212, Exactw311, Exactw322, ExactR2112;
        NekDouble w212err = 0.0, w311err = 0.0, w322err = 0.0;
        NekDouble MF1sterr = 0.0, R2112err = 0.0;

        for (int i = 0; i < nq; i++)
        {
            w212  = Connection[0][1][i];
            w311  = Connection[1][0][i];
            w322  = Connection[2][1][i];
            R2112 = Curvature[0][i];

            xp = x[i];
            yp = y[i];
            zp = z[i];

            rad = sqrt(xp * xp + yp * yp + zp * zp);

            CartesianToNewSpherical(xp, yp, zp, sin_varphi, cos_varphi,
                                    sin_theta, cos_theta);

            // costheta =
            //    (xp * InitPx + yp * InitPy + zp * InitPz) / rad /
            //    radPole;
            // sintheta = sqrt(1.0 - costheta * costheta);

            theta_x = cos_theta * cos_varphi;
            theta_y = cos_theta * sin_varphi;
            theta_z = -sin_theta;

            // Only for the z=20 initiation case
            if ((Activation[i] > 0) && (fabs(zp) < rad) &&
                fabs(sin_theta) > Tol)
            {
                Exactw212 = (1.0 / rad) * (cos_theta / sin_theta);
                Exactw311 = -1.0 / rad;
                Exactw322 = -1.0 / rad;

                ExactR2112 = 1.0 / (rad * rad * sin_theta * sin_theta);

                tmp = fabs(w212 - Exactw212);
                w212err += tmp * tmp;

                tmp = fabs(w311 - Exactw311);
                w311err += tmp * tmp;

                tmp = fabs(w322 - Exactw322);
                w322err += tmp * tmp;

                tmp = fabs(R2112 - ExactR2112);
                R2112err += tmp * tmp;

                dx = MF1st[0][i] - theta_x;
                dy = MF1st[0][i + nq] - theta_y;
                dz = MF1st[0][i + 2 * nq] - theta_z;

                MF1sterr += dx * dx + dy * dy + dz * dz;

                cnt++;
            }
        }

        if (cnt > 0)
        {
            MF1sterr = sqrt(MF1sterr / (cnt * 3));

            w212err  = sqrt(w212err / cnt);
            w311err  = sqrt(w311err / cnt);
            w322err  = sqrt(w322err / cnt);
            R2112err = sqrt(R2112err / cnt);

            std::cout << "For " << cnt << " / " << nq << " ( "
                      << 100.0 * cnt / nq
                      << " % ) points, Error of MF1st = " << MF1sterr
                      << std::endl;
            std::cout << "w212 = " << w212err << ", w311 = " << w311err
                      << ", w322 = " << w322err << ", R2112 = " << R2112err
                      << std::endl;
        }
    }

    else if ((m_surfaceType == SolverUtils::ePlanePoint) ||
             (m_surfaceType == SolverUtils::ePlaneLeft))
    {
        // ELF on the plane
        int cnt = 0;
        NekDouble diffx, diffy, diffz, rad;
        NekDouble xp, yp, zp;
        NekDouble dx = 0.0, dy = 0.0;

        NekDouble MF1sterr = 0.0;
        NekDouble tmp, w211, w212, R2112;
        NekDouble Exactw211 = 0.0, Exactw212 = 0.0, ExactR2112 = 0.0;
        NekDouble w211err = 0.0, w212err = 0.0, R2112err = 0.0;

        for (int i = 0; i < nq; i++)
        {
            xp = x[i];
            yp = y[i];
            zp = z[i];

            w211  = Connection[0][0][i];
            w212  = Connection[0][1][i];
            R2112 = Curvature[0][i];

            diffx = xp - InitPx;
            diffy = yp - InitPy;
            diffz = zp - InitPz;

            rad = sqrt(diffx * diffx + diffy * diffy + diffz * diffz);

            if ((Activation[i] > 0) && (rad > Tol))
            {
                if (m_surfaceType == SolverUtils::ePlanePoint)
                {
                    Exactw211  = 0.0;
                    Exactw212  = 1.0 / rad;
                    ExactR2112 = 1.0 / rad / rad;

                    dx = MF1st[0][i] - xp / sqrt(xp * xp + yp * yp);
                    dy = MF1st[0][i + nq] - yp / sqrt(xp * xp + yp * yp);
                }

                else if (m_surfaceType == SolverUtils::ePlaneLeft)
                {
                    Exactw211  = 0.0;
                    Exactw212  = 0.0;
                    ExactR2112 = 0.0;

                    dx = MF1st[0][i] - 1.0;
                    dy = MF1st[0][i + nq] - 0.0;
                }

                MF1sterr += dx * dx + dy * dy;

                tmp = fabs(w211 - Exactw211);
                w211err += tmp * tmp;

                tmp = fabs(w212 - Exactw212);
                w212err += tmp * tmp;

                tmp = fabs(R2112 - ExactR2112);
                R2112err += tmp * tmp;
                cnt++;
            }
        }
        if (cnt)
        {
            MF1sterr = sqrt(MF1sterr / (cnt * 2));

            w211err  = sqrt(w211err / cnt);
            w212err  = sqrt(w212err / cnt);
            R2112err = sqrt(R2112err / cnt);

            std::cout << "For " << cnt << " / " << nq << " ( "
                      << 100.0 * cnt / nq
                      << " % ) points, Error of MF1st = " << MF1sterr
                      << std::endl;
            std::cout << "w211 = " << w211err << ", w212 = " << w212err
                      << ", R2112 = " << R2112err << std::endl;
        }
    }
} // namespace SolverUtils

/*
void MMFSystem::Test2DConnectionCurvature(
    const NekDouble InitPx, const NekDouble InitPy, const NekDouble
InitPz, const Array<OneD, const int> &Activation, const Array<OneD,
const Array<OneD, NekDouble>> &MF1st, const Array<OneD, const
Array<OneD, Array<OneD, NekDouble>>> &Connection, const Array<OneD,
const Array<OneD, NekDouble>> &Curvature, const int nstep)
{
    int nq = GetTotPoints();

    Array<OneD, NekDouble> x(nq);
    Array<OneD, NekDouble> y(nq);
    Array<OneD, NekDouble> z(nq);

    m_fields[0]->GetCoords(x, y, z);

    NekDouble Tol = 0.00001;

    // ELF on the Sphere
    if (m_surfaceType == SolverUtils::eSphere)
    {
        // Array<OneD, NekDouble> Exactw212(nq, 0.0);
        // Array<OneD, NekDouble> Exactw311(nq, 0.0);
        // Array<OneD, NekDouble> Exactw322(nq, 0.0);
        // Array<OneD, NekDouble> MF1sterr(nq, 0.0);
        // Array<OneD, NekDouble> ExactR2112(nq, 0.0);

        NekDouble rad, radPole;
        radPole = sqrt(InitPx * InitPx + InitPy * InitPy + InitPz *
InitPz);

        int cnt = 0;
        NekDouble xp, yp, zp, sin_varphi, cos_varphi, sin_theta,
cos_theta; NekDouble dx, dy, dz, theta_x, theta_y, theta_z;

        NekDouble tmp, w122, w221, w212, w311, w322;
        NekDouble R1212;
        NekDouble Exactw212, Exactw122, Exactw221, Exactw311, Exactw322,
            ExactR1212;
        NekDouble w212err = 0.0, w122err = 0.0, w221err = 0.0, w311err =
0.0, w322err  = 0.0; NekDouble MF1sterr = 0.0, R1212err = 0.0;

        // w111 = Connectionform[0][0]
        // w112 = Connectionform[0][1]
        // w121 = Connectionform[0][2]
        // w122 = Connectionform[0][3]
        // w211 = Connectionform[1][0]
        // w212 = Connectionform[1][1]
        // w221 = Connectionform[1][2]
        // w222 = Connectionform[1][3]

        for (int i = 0; i < nq; i++)
        {
            w122 = Connection[0][3][i];
            w212 = Connection[1][1][i];
            w221 = Connection[1][2][i];

            w311 = Connection[2][0][i];
            w322 = Connection[2][1][i];

            R1212 = Curvature[0][i];

            xp = x[i];
            yp = y[i];
            zp = z[i];

            CartesianToNewSpherical(xp, yp, zp, rad, sin_varphi,
cos_varphi, sin_theta, cos_theta);

            theta_x = cos_theta * cos_varphi;
            theta_y = cos_theta * sin_varphi;
            theta_z = -sin_theta;

            // Only for the z=20 initiation case
            if ((Activation[i] > 0) && (fabs(zp) < rad) &&
                fabs(sin_theta) > Tol)
            {
                // Exactw212 = (1.0 / rad) * (cos_theta / sin_theta);
                Exactw311 = -1.0 / rad;
                Exactw322 = -1.0 / rad;

                //  ExactR2112 = 1.0 / (rad * rad * sin_theta *
sin_theta);

                Exactw122 = (-1.0 / rad) * (cos_theta / sin_theta);
                Exactw212 = (1.0 / rad) * (cos_theta / sin_theta);
                Exactw221 = (1.0 / rad) * sin_theta * cos_theta;

                ExactR1212 = sin_theta * sin_theta;

                tmp = fabs(w122 - Exactw122);
                w122err += tmp * tmp;

                tmp = fabs(w212 - Exactw212);
                w212err += tmp * tmp;

                tmp = fabs(w221 - Exactw221);
                w221err += tmp * tmp;

                tmp = fabs(w311 - Exactw311);
                w311err += tmp * tmp;

                tmp = fabs(w322 - Exactw322);
                w322err += tmp * tmp;

                tmp = fabs(R1212 - ExactR1212);
                R1212err += tmp * tmp;

                dx = MF1st[0][i] - theta_x;
                dy = MF1st[0][i + nq] - theta_y;
                dz = MF1st[0][i + 2 * nq] - theta_z;

                MF1sterr += dx * dx + dy * dy + dz * dz;

                cnt++;
            }
        }

        if (cnt > 0)
        {
            MF1sterr = sqrt(MF1sterr / (cnt * 3));

            w212err  = sqrt(w212err / cnt);
            w122err  = sqrt(w122err / cnt);
            w221err  = sqrt(w221err / cnt);
            w311err  = sqrt(w311err / cnt);
            w322err  = sqrt(w322err / cnt);
            R1212err = sqrt(R1212err / cnt);

            std::cout << "For " << cnt << " / " << nq << " ( "
                     << 100.0 * cnt / nq
                     << " % ) points, Error of MF1st = " << MF1sterr
                     << std::endl;
            std::cout << "w122 = " << w122err << ", w212 = " << w212err
                     << ", w221 = " << w221err << ", w311 = " << w311err
                     << ", w322 = " << w322err << ", R1212 = " <<
R1212err
                     << std::endl
                     << std::endl;
        }

        // PlotConnectionError(Exactw212, nstep);
    }

else if ((m_surfaceType == SolverUtils::ePlanePoint) ||
         (m_surfaceType == SolverUtils::ePlaneLeft))
{
    // ELF on the plane
    int cnt = 0;
    NekDouble diffx, diffy, diffz, rad;
    NekDouble xp, yp, zp;
    NekDouble dx, dy;

    NekDouble MF1sterr = 0.0;
    NekDouble tmp, w211, w212, R2112;
    NekDouble Exactw211, Exactw212, ExactR2112;
    NekDouble w211err = 0.0, w212err = 0.0, R2112err = 0.0;

    for (int i = 0; i < nq; i++)
    {
        xp = x[i];
        yp = y[i];
        zp = z[i];

        w211  = Connection[0][0][i];
        w212  = Connection[0][1][i];
        R2112 = Curvature[0][i];

        diffx = xp - InitPx;
        diffy = yp - InitPy;
        diffz = zp - InitPz;

        rad = sqrt(diffx * diffx + diffy * diffy + diffz * diffz);

        if ((Activation[i] > 0) && (rad > Tol))
        {
            if (m_surfaceType == SolverUtils::ePlanePoint)
            {
                Exactw211  = 0.0;
                Exactw212  = 1.0 / rad;
                ExactR2112 = 1.0 / rad / rad;

                dx = MF1st[0][i] - xp / sqrt(xp * xp + yp * yp);
                dy = MF1st[0][i + nq] - yp / sqrt(xp * xp + yp * yp);
            }

            else if (m_surfaceType == SolverUtils::ePlaneLeft)
            {
                Exactw211  = 0.0;
                Exactw212  = 0.0;
                ExactR2112 = 0.0;

                dx = MF1st[0][i] - 1.0;
                dy = MF1st[0][i + nq] - 0.0;
            }

            MF1sterr += dx * dx + dy * dy;

            tmp = fabs(w211 - Exactw211);
            w211err += tmp * tmp;

            tmp = fabs(w212 - Exactw212);
            w212err += tmp * tmp;

            tmp = fabs(R2112 - ExactR2112);
            R2112err += tmp * tmp;
            cnt++;
        }
    }
    if (cnt)
    {
        MF1sterr = sqrt(MF1sterr / (cnt * 2));

        w211err  = sqrt(w211err / cnt);
        w212err  = sqrt(w212err / cnt);
        R2112err = sqrt(R2112err / cnt);

        std::cout << "For " << cnt << " / " << nq << " ( " << 100.0 *
cnt / nq
                 << " % ) points, Error of MF1st = " << MF1sterr <<
std::endl; std::cout << "w211 = " << w211err << ", w212 = " << w212err
                 << ", R2112 = " << R2112err << std::endl;
    }
}
} // namespace SolverUtils
*/

NekDouble MMFSystem::CompareVector(
    const Array<OneD, const Array<OneD, NekDouble>> &vector1,
    const Array<OneD, const Array<OneD, NekDouble>> &vector2)
{
    int nq = m_fields[0]->GetNpoints();

    Array<OneD, NekDouble> tmp(nq, 0.0);
    Array<OneD, NekDouble> diff(nq, 0.0);

    Array<OneD, NekDouble> mag(nq, 0.0);
    for (int i = 0; i < m_spacedim; ++i)
    {
        Vmath::Vvtvp(nq, vector1[i], 1, vector1[i], 1, mag, 1, mag, 1);

        Vmath::Vsub(nq, &vector1[i][0], 1, &vector2[i][0], 1, &diff[0], 1);
        Vmath::Vvtvp(nq, &diff[0], 1, &diff[0], 1, &tmp[0], 1, &tmp[0], 1);
    }

    Vmath::Vsqrt(nq, mag, 1, mag, 1);
    Vmath::Vsqrt(nq, tmp, 1, tmp, 1);

    Vmath::Vdiv(nq, tmp, 1, mag, 1, tmp, 1);

    return RootMeanSquare(tmp);
}

void MMFSystem::DeriveVector(const Array<OneD, const NekDouble> &veccoeff1,
                             const Array<OneD, const NekDouble> &veccoeff2,
                             const Array<OneD, const NekDouble> &vecbase1,
                             const Array<OneD, const NekDouble> &vecbase2,
                             Array<OneD, Array<OneD, NekDouble>> &outarray)
{
    int nq = veccoeff1.size();

    outarray = Array<OneD, Array<OneD, NekDouble>>(m_spacedim);
    for (int k = 0; k < m_spacedim; ++k)
    {
        outarray[k] = Array<OneD, NekDouble>(nq);
    }

    for (int i = 0; i < nq; ++i)
    {
        outarray[0][i] =
            veccoeff1[i] * vecbase1[i] + veccoeff2[i] * vecbase2[i];
        outarray[1][i] =
            veccoeff1[i] * vecbase1[i + nq] + veccoeff2[i] * vecbase2[i + nq];
        outarray[2][i] = veccoeff1[i] * vecbase1[i + 2 * nq] +
                         veccoeff2[i] * vecbase2[i + 2 * nq];
    }
}

void MMFSystem::DeriveVector(
    const Array<OneD, const Array<OneD, NekDouble>> &veccoeff,
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    Array<OneD, Array<OneD, NekDouble>> &outarray)
{
    int nq = veccoeff[0].size();

    outarray = Array<OneD, Array<OneD, NekDouble>>(m_spacedim);
    for (int k = 0; k < m_spacedim; ++k)
    {
        outarray[k] = Array<OneD, NekDouble>(nq);
    }

    for (int i = 0; i < nq; ++i)
    {
        outarray[0][i] = veccoeff[0][i] * movingframes[0][i] +
                         veccoeff[1][i] * movingframes[1][i];
        outarray[1][i] = veccoeff[0][i] * movingframes[0][i + nq] +
                         veccoeff[1][i] * movingframes[1][i + nq];
        outarray[2][i] = veccoeff[0][i] * movingframes[0][i + 2 * nq] +
                         veccoeff[1][i] * movingframes[1][i + 2 * nq];
    }
}

// Update fields according to the new frames:
// u e^1 + v e^2 = unew enew^1 + vnew enew^2
// unew = u e^1 \cdot enew^1 + v e^2 \cdot enew^1
// vnew = u e^1 \cdot enew^2 + v e^2 \cdot enew^2
void MMFSystem::DeriveNewFieldComponent(
    const Array<OneD, const Array<OneD, NekDouble>> &MF1,
    const Array<OneD, const NekDouble> &u1,
    const Array<OneD, const NekDouble> &v1,
    const Array<OneD, const Array<OneD, NekDouble>> &MF2,
    Array<OneD, NekDouble> &u2, Array<OneD, NekDouble> &v2)
{
    int nq = m_fields[0]->GetNpoints();

    Array<OneD, NekDouble> e1_cdot_enew1(nq, 0.0);
    Array<OneD, NekDouble> e2_cdot_enew1(nq, 0.0);
    Array<OneD, NekDouble> e1_cdot_enew2(nq, 0.0);
    Array<OneD, NekDouble> e2_cdot_enew2(nq, 0.0);

    for (int k = 0; k < m_spacedim; ++k)
    {
        Vmath::Vvtvp(nq, &MF1[0][k * nq], 1, &MF2[0][k * nq], 1,
                     &e1_cdot_enew1[0], 1, &e1_cdot_enew1[0], 1);
        Vmath::Vvtvp(nq, &MF1[1][k * nq], 1, &MF2[0][k * nq], 1,
                     &e2_cdot_enew1[0], 1, &e2_cdot_enew1[0], 1);
        Vmath::Vvtvp(nq, &MF1[0][k * nq], 1, &MF2[1][k * nq], 1,
                     &e1_cdot_enew2[0], 1, &e1_cdot_enew2[0], 1);
        Vmath::Vvtvp(nq, &MF1[1][k * nq], 1, &MF2[1][k * nq], 1,
                     &e2_cdot_enew2[0], 1, &e2_cdot_enew2[0], 1);
    }

    u2 = Array<OneD, NekDouble>(nq, 0.0);
    Vmath::Vvtvp(nq, e1_cdot_enew1, 1, u1, 1, u2, 1, u2, 1);
    Vmath::Vvtvp(nq, e2_cdot_enew1, 1, v1, 1, u2, 1, u2, 1);

    v2 = Array<OneD, NekDouble>(nq, 0.0);
    Vmath::Vvtvp(nq, e1_cdot_enew2, 1, u1, 1, v2, 1, v2, 1);
    Vmath::Vvtvp(nq, e2_cdot_enew2, 1, v1, 1, v2, 1, v2, 1);
}

// Compute Weak DG Divergence for \vec{v} u
// inarray: u in nodal space
// movingframes & velocity
// outarray: \nabla \cdot \vec{v} u in nodal space
void MMFSystem::ComputeWeakDGDivergence(
    const Array<OneD, const NekDouble> &inarray,
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    const Array<OneD, const Array<OneD, NekDouble>> &velocity,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const int SurfaceDivergence)
{
    int nvariables = m_fields.size();
    int ncoeffs    = m_fields[0]->GetNcoeffs();
    int nq         = GetNpoints();

    Array<OneD, Array<OneD, NekDouble>> physfield(nvariables);

    // Get the variables in physical space
    // already in physical space
    for (int i = 0; i < nvariables; ++i)
    {
        physfield[i] = Array<OneD, NekDouble>(nq, 0.0);
        outarray[i]  = Array<OneD, NekDouble>(nq);
    }

    Vmath::Vcopy(nq, &inarray[0], 1, &physfield[0][0], 1);

    Array<OneD, Array<OneD, NekDouble>> WeakAdv(nvariables);
    WeakAdv[0] = Array<OneD, NekDouble>(ncoeffs * nvariables);
    for (int i = 1; i < nvariables; ++i)
    {
        WeakAdv[i] = WeakAdv[i - 1] + ncoeffs;
    }

    // Compute movingframes with velMF
    Array<OneD, Array<OneD, NekDouble>> veldotMF;
    ComputeveldotMF(movingframes, velocity, veldotMF);

    Array<OneD, Array<OneD, NekDouble>> movingframes_VelMF(m_shapedim);
    for (int j = 0; j < m_shapedim; j++)
    {
        movingframes_VelMF[j] = Array<OneD, NekDouble>(3 * nq);
    }

    // CAUTION: Modify e^i as v^i e^i
    for (int j = 0; j < m_shapedim; j++)
    {
        for (int k = 0; k < m_spacedim; k++)
        {
            Vmath::Vmul(nq, &veldotMF[j][0], 1, &m_movingframes[j][k * nq], 1,
                        &movingframes_VelMF[j][k * nq], 1);
        }
    }

    // Compute \nabla \cdot \vel u according to MMF scheme
    WeakDGDivergence(physfield, movingframes_VelMF, velocity, WeakAdv);

    Array<OneD, NekDouble> tmpc(ncoeffs);
    for (int i = 0; i < nvariables; ++i)
    {
        if (SurfaceDivergence)
        {
            Array<OneD, NekDouble> geometricdiv(nq);
            Array<OneD, NekDouble> velvector(m_spacedim * nq);
            for (int i = 0; i < m_spacedim; ++i)
            {
                Vmath::Vcopy(nq, &velocity[i][0], 1, &velvector[i * nq], 1);
            }

            geometricdiv = ComputeSurfaceDiv(movingframes[2], velvector);

            m_fields[0]->IProductWRTBase(geometricdiv, tmpc);
            Vmath::Vadd(ncoeffs, tmpc, 1, WeakAdv[i], 1, WeakAdv[i], 1);
        }

        //  std::cout << "MultiplybyElmtInvMass" << std::endl;
        //  m_fields[i]->MultiplyByElmtInvMass(WeakAdv[i], WeakAdv[i]);
        m_fields[i]->BwdTrans(WeakAdv[i], outarray[i]);
    }
}

// Compute Weak DG Divergence for \vec{v} u
// inarray: u in nodal space
// movingframes & velocity
// outarray: \nabla \cdot \vec{v} u in nodal space
void MMFSystem::ComputeWeakDGDivergence(
    const Array<OneD, const NekDouble> &inarray,
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    const Array<OneD, const Array<OneD, NekDouble>> &velocity,
    Array<OneD, Array<OneD, NekDouble>> &outarray,
    const Array<OneD, const NekDouble> &SpuriousDiffusion)
{
    int nvariables = m_fields.size();
    int ncoeffs    = m_fields[0]->GetNcoeffs();
    int nq         = GetNpoints();

    Array<OneD, Array<OneD, NekDouble>> physfield(nvariables);

    // Get the variables in physical space
    // already in physical space
    for (int i = 0; i < nvariables; ++i)
    {
        physfield[i] = Array<OneD, NekDouble>(nq, 0.0);
        outarray[i]  = Array<OneD, NekDouble>(nq);
    }

    Vmath::Vcopy(nq, &inarray[0], 1, &physfield[0][0], 1);

    Array<OneD, Array<OneD, NekDouble>> WeakAdv(nvariables);
    WeakAdv[0] = Array<OneD, NekDouble>(ncoeffs * nvariables);
    for (int i = 1; i < nvariables; ++i)
    {
        WeakAdv[i] = WeakAdv[i - 1] + ncoeffs;
    }

    // Compute movingframes with velMF
    Array<OneD, Array<OneD, NekDouble>> veldotMF;
    ComputeveldotMF(movingframes, velocity, veldotMF);

    Array<OneD, Array<OneD, NekDouble>> movingframes_VelMF(m_shapedim);
    for (int j = 0; j < m_shapedim; j++)
    {
        movingframes_VelMF[j] = Array<OneD, NekDouble>(3 * nq);
    }

    // CAUTION: Modify e^i as v^i e^i
    for (int j = 0; j < m_shapedim; j++)
    {
        for (int k = 0; k < m_spacedim; k++)
        {
            Vmath::Vmul(nq, &veldotMF[j][0], 1, &m_movingframes[j][k * nq], 1,
                        &movingframes_VelMF[j][k * nq], 1);
        }
    }

    // Compute \nabla \cdot \vel u according to MMF scheme
    WeakDGDivergence(physfield, movingframes_VelMF, velocity, WeakAdv);

    Array<OneD, NekDouble> tmpc(ncoeffs);
    for (int i = 0; i < nvariables; ++i)
    {
        if (SpuriousDiffusion != NullNekDouble1DArray)
        {
            m_fields[0]->IProductWRTBase(SpuriousDiffusion, tmpc);
            Vmath::Vadd(ncoeffs, tmpc, 1, WeakAdv[i], 1, WeakAdv[i], 1);
        }

        //  std::cout << "MultiplybyElmtInvMass" << std::endl;
        //  m_fields[i]->MultiplyByElmtInvMass(WeakAdv[i], WeakAdv[i]);
        m_fields[i]->BwdTrans(WeakAdv[i], outarray[i]);
    }
}

// Compute Weak DG for (\nabla \times \vec{v}) \cdot e^3
// (\nabla \times \vec{v}) \cdot e^3 =  \nabla \cdot (v2 e^1 - v1 e^2) + v1 e^1
// \cdot (\nabla \times e^3 ) + v2 e^2 \cdot (\nabla \times e^3 ) inarray: u in
// nodal space movingframes & velocity outarray: \nabla \cdot \vec{v} u in nodal
// space
void MMFSystem::ComputeWeakDGCurl(
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    const Array<OneD, const Array<OneD, NekDouble>> &velocity,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const int SurfaceDivergence)
{
    int nvariables = m_fields.size();
    int ncoeffs    = m_fields[0]->GetNcoeffs();
    int nq         = GetNpoints();

    // Get the variables in physical space
    // already in physical space
    Array<OneD, Array<OneD, NekDouble>> physfield(nvariables);
    for (int i = 0; i < nvariables; ++i)
    {
        outarray[i]  = Array<OneD, NekDouble>(nq);
        physfield[i] = Array<OneD, NekDouble>(nq);
    }

    Array<OneD, Array<OneD, NekDouble>> WeakAdv(nvariables);
    WeakAdv[0] = Array<OneD, NekDouble>(ncoeffs * nvariables);
    for (int i = 1; i < nvariables; ++i)
    {
        WeakAdv[i] = WeakAdv[i - 1] + ncoeffs;
    }

    // Compute movingframes with velMF
    Array<OneD, Array<OneD, NekDouble>> veldotMF;
    ComputeveldotMF(movingframes, velocity, veldotMF);

    Array<OneD, Array<OneD, NekDouble>> movingframes_VelMF(m_shapedim);
    for (int j = 0; j < m_shapedim; j++)
    {
        movingframes_VelMF[j] = Array<OneD, NekDouble>(3 * nq);
    }

    // CAUTION: Modify e^i as v^i e^i
    for (int j = 0; j < m_shapedim; j++)
    {
        for (int k = 0; k < m_spacedim; k++)
        {
            Vmath::Vmul(nq, &veldotMF[j][0], 1, &m_movingframes[j][k * nq], 1,
                        &movingframes_VelMF[j][k * nq], 1);
        }
    }

    // Compute \nabla \cdot \vel u according to MMF scheme
    WeakDGDivergence(physfield, movingframes_VelMF, velocity, WeakAdv);

    Array<OneD, NekDouble> tmpc(ncoeffs);
    for (int i = 0; i < nvariables; ++i)
    {
        if (SurfaceDivergence)
        {
            Array<OneD, NekDouble> geometricdiv(nq);
            Array<OneD, NekDouble> velvector(m_spacedim * nq);
            for (int i = 0; i < m_spacedim; ++i)
            {
                Vmath::Vcopy(nq, &velocity[i][0], 1, &velvector[i * nq], 1);
            }

            geometricdiv = ComputeSurfaceDiv(movingframes[2], velvector);

            m_fields[0]->IProductWRTBase(geometricdiv, tmpc);
            Vmath::Vadd(ncoeffs, tmpc, 1, WeakAdv[i], 1, WeakAdv[i], 1);
        }

        //  std::cout << "MultiplybyElmtInvMass" << std::endl;
        //  m_fields[i]->MultiplyByElmtInvMass(WeakAdv[i], WeakAdv[i]);
        m_fields[i]->BwdTrans(WeakAdv[i], outarray[i]);
    }
}

Array<OneD, NekDouble> MMFSystem::ComputeSurfaceDiv(
    const Array<OneD, const NekDouble> &surfaceNormal,
    const Array<OneD, const NekDouble> &velvector)
{
    int nq = m_fields[0]->GetNpoints();

    Array<OneD, NekDouble> k_v_nabla_k(nq);
    Array<OneD, NekDouble> k_k_nabla_v(nq);
    Array<OneD, NekDouble> outarray(nq);

    // Compute \vec{C} \cdot ( ( \vec{A} \cdot \nabla ) \vec{B} )
    k_v_nabla_k = ComputeVecCdotNabla(surfaceNormal, velvector, surfaceNormal);
    k_k_nabla_v = ComputeVecCdotNabla(surfaceNormal, surfaceNormal, velvector);
    Vmath::Neg(nq, k_v_nabla_k, 1);

    Vmath::Vadd(nq, k_v_nabla_k, 1, k_k_nabla_v, 1, outarray, 1);

    return outarray;
}

// Compute SpuriousDivergence
// Input: movingframes, surface normal, velocity vector
// output: scalar
Array<OneD, NekDouble> MMFSystem::ComputeSpuriousDivergence(
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    const Array<OneD, const NekDouble> &surfaceNormal,
    const Array<OneD, const NekDouble> &velvector)
{
    int nq = m_fields[0]->GetNpoints();
    Array<OneD, NekDouble> outarray(nq, 0.0);

    // Compute k_i, v_i
    Array<OneD, Array<OneD, NekDouble>> SNmf(m_spacedim);
    Array<OneD, Array<OneD, NekDouble>> Velmf(m_spacedim);
    for (int i = 0; i < m_spacedim; ++i)
    {
        SNmf[i]  = Array<OneD, NekDouble>(nq, 0.0);
        Velmf[i] = Array<OneD, NekDouble>(nq, 0.0);
        for (int k = 0; k < m_spacedim; ++k)
        {
            Vmath::Vvtvp(nq, &surfaceNormal[k * nq], 1,
                         &movingframes[i][k * nq], 1, &SNmf[i][0], 1,
                         &SNmf[i][0], 1);
            Vmath::Vvtvp(nq, &velvector[k * nq], 1, &movingframes[i][k * nq], 1,
                         &Velmf[i][0], 1, &Velmf[i][0], 1);
        }
    }

    Array<OneD, NekDouble> tmp(nq);
    Array<OneD, NekDouble> Dtmp(nq);
    Array<OneD, NekDouble> locsum(nq);

    Array<OneD, NekDouble> locsum1(nq);
    Array<OneD, NekDouble> locsum2(nq);

    // Compute k_i [ k_1 (\nabla v_i \cdot e^1) + k2 (\nabla v_i \cdot e^2) ]
    for (int i = 0; i < m_spacedim; ++i)
    {
        locsum1 = Array<OneD, NekDouble>(nq, 0.0);

        // Add k_1 (\nabla v_i \cdot e^1) + k2 (\nabla v_i \cdot e^2)
        Vmath::Vcopy(nq, &Velmf[i][0], 1, &tmp[0], 1);
        for (int j = 0; j < m_expdim; ++j)
        {
            MMFDirectionalDeriv(movingframes[j], tmp, Dtmp);
            Vmath::Vvtvp(nq, &SNmf[j][0], 1, &Dtmp[0], 1, &locsum1[0], 1,
                         &locsum1[0], 1);
        }

        // Add -v_1 (\nabla k_i \cdot e^1) - v_2 (\nabla k_i \cdot e^2)
        locsum2 = Array<OneD, NekDouble>(nq, 0.0);

        Vmath::Vcopy(nq, &SNmf[i][0], 1, &tmp[0], 1);
        Vmath::Neg(nq, &tmp[0], 1);
        for (int j = 0; j < m_expdim; ++j)
        {
            MMFDirectionalDeriv(movingframes[j], tmp, Dtmp);
            Vmath::Vvtvp(nq, &Velmf[j][0], 1, &Dtmp[0], 1, &locsum2[0], 1,
                         &locsum2[0], 1);
        }

        Vmath::Vadd(nq, locsum1, 1, locsum2, 1, locsum, 1);
        Vmath::Vvtvp(nq, &SNmf[i][0], 1, &locsum[0], 1, &outarray[0], 1,
                     &outarray[0], 1);
    }

    return outarray;
}

void MMFSystem::ComputeveldotMF(
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    const Array<OneD, const Array<OneD, NekDouble>> &velocity,
    Array<OneD, Array<OneD, NekDouble>> &veldotMF)
{
    int nq = m_fields[0]->GetNpoints();

    veldotMF = Array<OneD, Array<OneD, NekDouble>>(m_shapedim);
    for (int j = 0; j < m_shapedim; ++j)
    {
        veldotMF[j] = Array<OneD, NekDouble>(nq, 0.0);
    }

    for (int j = 0; j < m_shapedim; ++j)
    {
        for (int k = 0; k < m_spacedim; ++k)
        {
            Vmath::Vvtvp(nq, &movingframes[j][k * nq], 1, &velocity[k][0], 1,
                         &veldotMF[j][0], 1, &veldotMF[j][0], 1);
        }
    }
}

// Compute Weak DG Divergence in moving frames
// InField: phyfield u in nodal space
// movingframes: orthonormal bases
// velocity: velocity vector \vec{v}
// OutField: outfield in modal space
void MMFSystem::WeakDGDivergence(
    const Array<OneD, const Array<OneD, NekDouble>> &InField,
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    const Array<OneD, const Array<OneD, NekDouble>> &velocity,
    Array<OneD, Array<OneD, NekDouble>> &OutField)
{
    int nvariables      = m_fields.size();
    int ncoeffs         = GetNcoeffs();
    int nTracePointsTot = GetTraceNpoints();

    Array<OneD, Array<OneD, NekDouble>> physfield(nvariables);
    // Get the variables in physical space already in physical space
    for (int i = 0; i < nvariables; ++i)
    {
        physfield[i] = InField[i];
    }

    Array<OneD, Array<OneD, NekDouble>> WeakDeriv(m_shapedim);
    for (int i = 0; i < nvariables; ++i)
    {
        for (int j = 0; j < m_shapedim; ++j)
        {
            WeakDeriv[j] = Array<OneD, NekDouble>(ncoeffs, 0.0);

            // Directional derivation with respect to the j'th moving
            // frame tmp[j] = \nabla \physfield[i] \cdot \mathbf{e}^j
            // Implemented at
            // TriExp::v_IProductWRTDirectionalDerivBase_SumFa
            m_fields[i]->IProductWRTDirectionalDerivBase(
                movingframes[j], physfield[i], WeakDeriv[j]);
        }

        // if the NumericalFluxs function already includes the normal in
        // the output
        Array<OneD, NekDouble> Fwd(nTracePointsTot);
        Array<OneD, NekDouble> Bwd(nTracePointsTot);

        Array<OneD, NekDouble> flux(nTracePointsTot, 0.0);
        Array<OneD, NekDouble> fluxFwd(nTracePointsTot);
        Array<OneD, NekDouble> fluxBwd(nTracePointsTot);

        // Evaluate numerical flux in physical space which may in
        // general couple all component of vectors
        m_fields[i]->GetFwdBwdTracePhys(physfield[i], Fwd, Bwd);

        // evaulate upwinded m_fields[i]
        Array<OneD, NekDouble> traceVn;
        traceVn = GetNormalVelocity(velocity);
        m_fields[i]->GetTrace()->Upwind(traceVn, Fwd, Bwd, flux);

        Array<OneD, Array<OneD, NekDouble>> ncdotMF_VelMF_Fwd;
        Array<OneD, Array<OneD, NekDouble>> ncdotMF_VelMF_Bwd;
        ComputencdotMF(movingframes, ncdotMF_VelMF_Fwd, ncdotMF_VelMF_Bwd);

        OutField[i] = Array<OneD, NekDouble>(ncoeffs, 0.0);
        for (int j = 0; j < m_shapedim; ++j)
        {
            // calculate numflux = (n \cdot MF)*flux
            Vmath::Vmul(nTracePointsTot, &flux[0], 1, &ncdotMF_VelMF_Fwd[j][0],
                        1, &fluxFwd[0], 1);
            Vmath::Vmul(nTracePointsTot, &flux[0], 1, &ncdotMF_VelMF_Bwd[j][0],
                        1, &fluxBwd[0], 1);
            Vmath::Neg(ncoeffs, WeakDeriv[j], 1);

            // FwdBwdtegral because generallize (N \cdot MF)_{FWD} \neq
            // -(N \cdot MF)_{BWD}
            m_fields[i]->AddFwdBwdTraceIntegral(fluxFwd, fluxBwd, WeakDeriv[j]);
            m_fields[i]->SetPhysState(false);

            Vmath::Vadd(ncoeffs, &WeakDeriv[j][0], 1, &OutField[i][0], 1,
                        &OutField[i][0], 1);
        }
    }
}

// void MMFSystem::WeakDGGradientVector(
//     const Array<OneD, const NekDouble> &inarray,
//     const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
//     Array<OneD, NekDouble> &outarray,
//     const int SurfaceGradient)
// {
//     int nq   = GetNpoints();

//     Array<OneD, NekDouble> physfield(nq, 0.0);
//     Vmath::Vcopy(nq, &inarray[0], 1, &physfield[0], 1);

//     // Compute Directional derivation with respect to the j'th moving frame
//     // tmp[j] = \nabla \physfield[i] \cdot \mathbf{e}^j
//     Array<OneD, NekDouble> WeakAdvtmp(m_spacedim*nq);
//     outarray = Array<OneD, NekDouble>(m_spacedim * nq, 0.0);
//     for (int j = 0; j < m_shapedim; ++j)
//     {
//         WeakDGDirectionalDeriv(j, physfield, movingframes, WeakAdvtmp,
//         SurfaceGradient); Vmath::Vadd(m_spacedim*nq, &WeakAdvtmp[0], 1,
//         &outarray[0], 1, &outarray[0], 1);
//     }
// }

void MMFSystem::WeakDGGradientVector(
    const Array<OneD, const NekDouble> &inarray,
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    Array<OneD, NekDouble> &outarray, const int SurfaceGradient)
{
    int nq = GetNpoints();

    Array<OneD, NekDouble> physfield(nq, 0.0);
    Vmath::Vcopy(nq, &inarray[0], 1, &physfield[0], 1);

    Array<OneD, NekDouble> WeakAdvtmp(nq);
    outarray = Array<OneD, NekDouble>(m_spacedim * nq, 0.0);
    // Compute Directional derivation with respect to the j'th moving frame
    // tmp[j] = \nabla \physfield[i] \cdot \mathbf{e}^j
    for (int j = 0; j < m_shapedim; ++j)
    {
        WeakDGDirectionalDeriv(j, physfield, movingframes, WeakAdvtmp,
                               SurfaceGradient);

        // Three dimensional vector
        for (int k = 0; k < m_spacedim; ++k)
        {
            Vmath::Vvtvp(nq, &WeakAdvtmp[0], 1, &movingframes[j][k * nq], 1,
                         &outarray[k * nq], 1, &outarray[k * nq], 1);
        }
    }
}

// WeakDG without e^i
void MMFSystem::WeakDGDirectionalDeriv(
    const int dir, const Array<OneD, const NekDouble> &InField,
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    Array<OneD, NekDouble> &OutField, const int SurfaceGradient)
{
    int nq              = GetNpoints();
    int ncoeffs         = GetNcoeffs();
    int nTracePointsTot = GetTraceNpoints();

    // Get the variables in physical space
    // already in physical space
    Array<OneD, NekDouble> tmp(nq);
    Array<OneD, NekDouble> tmpc(ncoeffs);

    Array<OneD, NekDouble> physfield(nq);
    Vmath::Vcopy(nq, InField, 1, physfield, 1);

    Array<OneD, Array<OneD, NekDouble>> velocity(m_spacedim);
    for (int k = 0; k < m_spacedim; ++k)
    {
        velocity[k] = Array<OneD, NekDouble>(nq);
    }

    Array<OneD, NekDouble> WeakDeriv(ncoeffs, 0.0);
    m_fields[0]->IProductWRTDirectionalDerivBase(movingframes[dir], physfield,
                                                 WeakDeriv);

    // if the NumericalFluxs function already includes the normal in
    // the output
    Array<OneD, NekDouble> Fwd(nTracePointsTot);
    Array<OneD, NekDouble> Bwd(nTracePointsTot);

    Array<OneD, NekDouble> flux(nTracePointsTot, 0.0);
    Array<OneD, NekDouble> fluxFwd(nTracePointsTot);
    Array<OneD, NekDouble> fluxBwd(nTracePointsTot);

    // Evaluate numerical flux in physical space which may in
    // general couple all component of vectors
    m_fields[0]->GetFwdBwdTracePhys(physfield, Fwd, Bwd);

    // evaulate upwinded m_fields[i]
    Array<OneD, NekDouble> traceVn;
    for (int k = 0; k < m_spacedim; ++k)
    {
        Vmath::Vcopy(nq, &movingframes[dir][k * nq], 1, &velocity[k][0], 1);
    }
    traceVn = GetNormalVelocity(velocity);
    m_fields[0]->GetTrace()->Upwind(traceVn, Fwd, Bwd, flux);

    // calculate numflux = (n \cdot MF)*flux
    Vmath::Vmul(nTracePointsTot, &flux[0], 1, &m_ncdotMFFwd[dir][0], 1,
                &fluxFwd[0], 1);
    Vmath::Vmul(nTracePointsTot, &flux[0], 1, &m_ncdotMFBwd[dir][0], 1,
                &fluxBwd[0], 1);

    // FwdBwdtegral because generallize (N \cdot MF)_{FWD} \neq -(N
    // \cdot MF)_{BWD}
    Vmath::Neg(ncoeffs, WeakDeriv, 1);
    m_fields[0]->AddFwdBwdTraceIntegral(fluxFwd, fluxBwd, WeakDeriv);
    m_fields[0]->SetPhysState(false);

    // Add DivMF
    Vmath::Vmul(nq, &physfield[0], 1, &m_DivMF[dir][0], 1, &tmp[0], 1);
    Vmath::Neg(nq, &tmp[0], 1);
    m_fields[0]->IProductWRTBase(tmp, tmpc);
    Vmath::Vadd(ncoeffs, &tmpc[0], 1, &WeakDeriv[0], 1, &WeakDeriv[0], 1);

    if (SurfaceGradient)
    {
        Array<OneD, NekDouble> velvector(m_spacedim * nq);
        for (int k = 0; k < m_spacedim; ++k)
        {
            Vmath::Vmul(nq, &physfield[0], 1, &movingframes[dir][k * nq], 1,
                        &velvector[k * nq], 1);
        }

        tmp = ComputeSurfaceDiv(movingframes[2], velvector);

        m_fields[0]->IProductWRTBase(tmp, tmpc);
        Vmath::Vadd(ncoeffs, &tmpc[0], 1, &WeakDeriv[0], 1, &WeakDeriv[0], 1);
    }

    m_fields[0]->BwdTrans(WeakDeriv, OutField);
}

// void MMFSystem::WeakDGDirectionalDerivwithMF(
//     const int dir, const Array<OneD, const NekDouble> &InField,
//     const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
//     Array<OneD, NekDouble> &OutField,
//     const int SurfaceGradient)
//     {
//         int nq              = GetNpoints();

//         Array<OneD, NekDouble> outarraytmp(nq);
//         for (int k=0; k<m_spacedim; ++k)
//         {
//             WeakDGDirectionalDerivwithMF(dir, k, InField, movingframes,
//             outarraytmp, SurfaceGradient); Vmath::Vcopy(nq, &outarraytmp[0],
//             1, &OutField[k*nq], 1);
//         }
//     }

// // WeakDG including e^i
// void MMFSystem::WeakDGDirectionalDerivwithMF(
//     const int dir, const int xj, const Array<OneD, const NekDouble> &InField,
//     const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
//     Array<OneD, NekDouble> &OutField,
//     const int SurfaceGradient)
// {
//     int nq              = GetNpoints();
//     int ncoeffs         = GetNcoeffs();
//     int nTracePointsTot = GetTraceNpoints();

//     // Get the variables in physical space
//     // already in physical space
//     Array<OneD, NekDouble> tmp(nq);
//     Array<OneD, NekDouble> tmpc(ncoeffs);

//     Array<OneD, NekDouble> physfield(nq);
//     Vmath::Vcopy(nq, InField, 1, physfield, 1);

//     Array<OneD, NekDouble> MFxj(nq);
//     Vmath::Vcopy(nq, &movingframes[dir][xj * nq], 1, &MFxj[0], 1);

//     Array<OneD, NekDouble> physfieldxj(nq);
//     Vmath::Vmul(nq, &physfield[0], 1, &MFxj[0], 1, &physfieldxj[0], 1);

//     Array<OneD, NekDouble> WeakDeriv(ncoeffs, 0.0);
//     m_fields[0]->IProductWRTDirectionalDerivBase(movingframes[dir],
//     physfield, WeakDeriv);

//     // if the NumericalFluxs function already includes the normal in
//     // the output
//     Array<OneD, NekDouble> Fwd(nTracePointsTot);
//     Array<OneD, NekDouble> Bwd(nTracePointsTot);

//     Array<OneD, NekDouble> flux(nTracePointsTot, 0.0);
//     Array<OneD, NekDouble> fluxFwd(nTracePointsTot);
//     Array<OneD, NekDouble> fluxBwd(nTracePointsTot);

//     // Evaluate numerical flux in physical space which may in
//     // general couple all component of vectors
//     Array<OneD, Array<OneD, NekDouble>> velocity(m_spacedim);
//     for (int k = 0; k < m_spacedim; ++k)
//     {
//         velocity[k] = Array<OneD, NekDouble>(nq);
//     }

//     m_fields[0]->GetFwdBwdTracePhys(physfieldxj, Fwd, Bwd);

//     // evaulate upwinded m_fields[i]
//     Array<OneD, NekDouble> traceVn;
//     for (int k = 0; k < m_spacedim; ++k)
//     {
//         Vmath::Vmul(nq, &MFxj[0], 1, &movingframes[dir][k * nq], 1,
//         &velocity[k][0], 1);
//     }
//     traceVn = GetNormalVelocity(velocity);
//     m_fields[0]->GetTrace()->Upwind(traceVn, Fwd, Bwd, flux);

//     // calculate numflux = (n \cdot MF)*flux
//     Vmath::Vmul(nTracePointsTot, &flux[0], 1, &m_ncdotMFFwd[dir][0], 1,
//     &fluxFwd[0], 1); Vmath::Vmul(nTracePointsTot, &flux[0], 1,
//     &m_ncdotMFBwd[dir][0], 1, &fluxBwd[0], 1);

//     // FwdBwdtegral because generallize (N \cdot MF)_{FWD} \neq -(N
//     // \cdot MF)_{BWD}
//     Vmath::Neg(ncoeffs, WeakDeriv, 1);
//     m_fields[0]->AddFwdBwdTraceIntegral(fluxFwd, fluxBwd, WeakDeriv);
//     m_fields[0]->SetPhysState(false);

//     // Add (\nabla \cdot e^i) * e^i_j * f
//     Vmath::Vmul(nq, &physfieldxj[0], 1, &m_DivMF[dir][0], 1, &tmp[0], 1);
//     Vmath::Neg(nq, &tmp[0], 1);
//     m_fields[0]->IProductWRTBase(tmp, tmpc);
//     // Vmath::Vadd(ncoeffs, &tmpc[0], 1, &WeakDeriv[0], 1, &WeakDeriv[0], 1);

//     // Add (\nabla e^i_j \cdot e^i) * f
//     // Vmath::Vmul(nq, &physfieldxj[0], 1, &m_DivMF[dir][0], 1, &tmp[0], 1);
//     m_fields[0]->PhysDirectionalDeriv(m_movingframes[dir], MFxj, tmp);
//     Vmath::Neg(nq, &tmp[0], 1);
//     Vmath::Vmul(nq, physfield, 1, tmp, 1, tmp, 1);
//     m_fields[0]->IProductWRTBase(tmp, tmpc);
//    // Vmath::Vadd(ncoeffs, &tmpc[0], 1, &WeakDeriv[0], 1, &WeakDeriv[0], 1);

//     if (SurfaceGradient)
//       {
//     //     Array<OneD, NekDouble> velvector(m_spacedim * nq);
//     //     for (int k = 0; k < m_spacedim; ++k)
//     //     {
//     //         Vmath::Vmul(nq, &physfield[0], 1, &movingframes[dir][k * nq],
//     1,
//     //                     &velvector[k * nq], 1);
//     //     }

//     //     tmp = ComputeSurfaceDiv(movingframes[2], velvector);

//     //     m_fields[0]->IProductWRTBase(tmp, tmpc);
//     //     Vmath::Vadd(ncoeffs, &tmpc[0], 1, &WeakDeriv[0], 1, &WeakDeriv[0],
//     1);
//       }

//     m_fields[0]->BwdTrans(WeakDeriv, OutField);
// }

void MMFSystem::WeakDGCurl(
    const Array<OneD, const Array<OneD, NekDouble>> &physfield,
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    const Array<OneD, const Array<OneD, NekDouble>> &CrossProductMF,
    Array<OneD, NekDouble> &outarray, const int SurfaceCurl)
{
    int nq              = GetNpoints();
    int ncoeffs         = GetNcoeffs();
    int nTracePointsTot = GetTraceNpoints();
    int nvar            = 3;

    Array<OneD, Array<OneD, NekDouble>> fluxvector(m_shapedim);
    for (int j = 0; j < m_shapedim; ++j)
    {
        fluxvector[j] = Array<OneD, NekDouble>(nq, 0.0);
    }

    //  Get Flux Vector
    Vmath::Smul(nq, 1.0, physfield[0], 1, fluxvector[0], 1);
    Vmath::Smul(nq, -1.0, physfield[1], 1, fluxvector[1], 1);

    Array<OneD, NekDouble> tmp(nq, 0.0);
    Array<OneD, NekDouble> tmpc(ncoeffs);
    Array<OneD, NekDouble> OutField(ncoeffs, 0.0);
    for (int j = 0; j < m_shapedim; ++j)
    {
        // Directional derivation with respect to the j'th moving frame
        // tmp_j = ( \nabla \phi, fluxvector[j] \mathbf{e}^j )
        m_fields[j]->IProductWRTDirectionalDerivBase(CrossProductMF[j],
                                                     fluxvector[j], tmpc);
        Vmath::Vadd(ncoeffs, &tmpc[0], 1, &OutField[0], 1, &OutField[0], 1);
    }

    // V the numerical flux and add to the modal coeffs
    // if the NumericalFlux function does not include the
    // normal in the outputDi
    Array<OneD, Array<OneD, NekDouble>> numfluxFwd(nvar);
    Array<OneD, Array<OneD, NekDouble>> numfluxBwd(nvar);
    for (int i = 0; i < nvar; ++i)
    {
        numfluxFwd[i] = Array<OneD, NekDouble>(nTracePointsTot, 0.0);
        numfluxBwd[i] = Array<OneD, NekDouble>(nTracePointsTot, 0.0);
    }

    // Evaluate numerical flux in physical space which may in
    // general couple all component of vectors
    NumericalCurlFlux(physfield, numfluxFwd, numfluxBwd);

    // Evaulate  <\phi, \hat{F}\cdot n> - OutField[i]
    m_fields[0]->AddFwdBwdTraceIntegral(numfluxFwd[2], numfluxBwd[2], OutField);

    // Add Green curl compensation
    Vmath::Vmul(nq, &physfield[0][0], 1, &m_CurlMF[2][0], 1, &tmp[0], 1);
    Vmath::Vvtvm(nq, &physfield[1][0], 1, &m_CurlMF[3][0], 1, &tmp[0], 1,
                 &tmp[0], 1);

    Vmath::Neg(nq, tmp, 1);
    m_fields[0]->IProductWRTBase(tmp, tmpc);
    Vmath::Vadd(ncoeffs, tmpc, 1, OutField, 1, OutField, 1);

    if (SurfaceCurl)
    {
        Array<OneD, NekDouble> velvector(m_spacedim * nq, 0.0);
        for (int i = 0; i < m_spacedim; ++i)
        {
            for (int j = 0; j < m_shapedim; ++j)
            {
                Vmath::Vvtvp(nq, &physfield[j][0], 1,
                             &CrossProductMF[j][i * nq], 1, &velvector[i * nq],
                             1, &velvector[i * nq], 1);
            }
        }

        tmp = ComputeSurfaceDiv(movingframes[2], velvector);
        Vmath::Neg(nq, tmp, 1);

        m_fields[0]->IProductWRTBase(tmp, tmpc);
        Vmath::Vadd(ncoeffs, tmpc, 1, OutField, 1, OutField, 1);
    }

    //  m_fields[0]->MultiplyByElmtInvMass(OutField, OutField);
    m_fields[0]->BwdTrans(OutField, outarray);
}

void MMFSystem::WeakDGCurl(
    const Array<OneD, const Array<OneD, NekDouble>> &physfield,
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    const Array<OneD, const Array<OneD, NekDouble>> &CrossProductMF,
    Array<OneD, NekDouble> &outarray,
    const Array<OneD, const NekDouble> &SpuriousDiffusion)
{
    boost::ignore_unused(movingframes);

    int nq              = GetNpoints();
    int ncoeffs         = GetNcoeffs();
    int nTracePointsTot = GetTraceNpoints();
    int nvar            = 3;

    Array<OneD, Array<OneD, NekDouble>> fluxvector(m_shapedim);
    for (int j = 0; j < m_shapedim; ++j)
    {
        fluxvector[j] = Array<OneD, NekDouble>(nq, 0.0);
    }

    //  Get Flux Vector
    Vmath::Smul(nq, 1.0, physfield[0], 1, fluxvector[0], 1);
    Vmath::Smul(nq, -1.0, physfield[1], 1, fluxvector[1], 1);

    Array<OneD, NekDouble> tmp(nq, 0.0);
    Array<OneD, NekDouble> tmpc(ncoeffs);
    Array<OneD, NekDouble> OutField(ncoeffs, 0.0);
    for (int j = 0; j < m_shapedim; ++j)
    {
        // Directional derivation with respect to the j'th moving frame
        // tmp_j = ( \nabla \phi, fluxvector[j] \mathbf{e}^j )
        m_fields[j]->IProductWRTDirectionalDerivBase(CrossProductMF[j],
                                                     fluxvector[j], tmpc);
        Vmath::Vadd(ncoeffs, &tmpc[0], 1, &OutField[0], 1, &OutField[0], 1);
    }

    // V the numerical flux and add to the modal coeffs
    // if the NumericalFlux function does not include the
    // normal in the outputDi
    Array<OneD, Array<OneD, NekDouble>> numfluxFwd(nvar);
    Array<OneD, Array<OneD, NekDouble>> numfluxBwd(nvar);
    for (int i = 0; i < nvar; ++i)
    {
        numfluxFwd[i] = Array<OneD, NekDouble>(nTracePointsTot, 0.0);
        numfluxBwd[i] = Array<OneD, NekDouble>(nTracePointsTot, 0.0);
    }

    // Evaluate numerical flux in physical space which may in
    // general couple all component of vectors
    NumericalCurlFlux(physfield, numfluxFwd, numfluxBwd);

    // Evaulate  <\phi, \hat{F}\cdot n> - OutField[i]
    m_fields[0]->AddFwdBwdTraceIntegral(numfluxFwd[2], numfluxBwd[2], OutField);

    // Add Green curl compensation
    Vmath::Vmul(nq, &physfield[0][0], 1, &m_CurlMF[2][0], 1, &tmp[0], 1);
    Vmath::Vvtvm(nq, &physfield[1][0], 1, &m_CurlMF[3][0], 1, &tmp[0], 1,
                 &tmp[0], 1);

    Vmath::Neg(nq, tmp, 1);
    m_fields[0]->IProductWRTBase(tmp, tmpc);
    Vmath::Vadd(ncoeffs, tmpc, 1, OutField, 1, OutField, 1);

    if (SpuriousDiffusion != NullNekDouble1DArray)
    {
        m_fields[0]->IProductWRTBase(SpuriousDiffusion, tmpc);
        Vmath::Vadd(ncoeffs, tmpc, 1, OutField, 1, OutField, 1);
    }

    //  m_fields[0]->MultiplyByElmtInvMass(OutField, OutField);
    m_fields[0]->BwdTrans(OutField, outarray);
}

void MMFSystem::NumericalCurlFlux(
    const Array<OneD, const Array<OneD, NekDouble>> &physfield,
    Array<OneD, Array<OneD, NekDouble>> &numfluxFwd,
    Array<OneD, Array<OneD, NekDouble>> &numfluxBwd)

{
    int nTraceNumPoints = GetTraceTotPoints();
    int nvar            = physfield.size();

    // get temporary arrays
    Array<OneD, Array<OneD, NekDouble>> Fwd(nvar);
    Array<OneD, Array<OneD, NekDouble>> Bwd(nvar);
    for (int i = 0; i < nvar; ++i)
    {
        Fwd[i] = Array<OneD, NekDouble>(nTraceNumPoints, 0.0);
        Bwd[i] = Array<OneD, NekDouble>(nTraceNumPoints, 0.0);

        numfluxFwd[i] = Array<OneD, NekDouble>(nTraceNumPoints, 0.0);
        numfluxBwd[i] = Array<OneD, NekDouble>(nTraceNumPoints, 0.0);

        // get the physical values at the trace
        m_fields[i]->GetFwdBwdTracePhys(physfield[i], Fwd[i], Bwd[i]);
    }

    // CopyBoundaryTrace(Fwd[0], Bwd[0], SolverUtils::eFwdEQBwd,    0,
    // SpatialDomains::ePEC); CopyBoundaryTrace(Fwd[1], Bwd[1],
    // SolverUtils::eFwdEQBwd,    1, SpatialDomains::ePEC);
    // CopyBoundaryTrace(Fwd[2], Bwd[2], SolverUtils::eFwdEQNegBwd, 2,
    // SpatialDomains::ePEC);

    // CopyBoundaryTrace(Fwd[0], Bwd[0], SolverUtils::eFwdEQNegBwd, 0,
    // SpatialDomains::ePMC); CopyBoundaryTrace(Fwd[1], Bwd[1],
    // SolverUtils::eFwdEQNegBwd, 1, SpatialDomains::ePMC);
    // CopyBoundaryTrace(Fwd[2], Bwd[2], SolverUtils::eFwdEQBwd,    2,
    // SpatialDomains::ePMC);

    CopyBoundaryTrace(Fwd[0], Bwd[0], SolverUtils::eFwdEQBwd, 2);
    CopyBoundaryTrace(Fwd[1], Bwd[1], SolverUtils::eFwdEQBwd, 2);

    // Array<OneD, NekDouble> e1Fwd_cdot_ncrossdH(nTraceNumPoints,0.0);
    // Array<OneD, NekDouble> e1Bwd_cdot_ncrossdH(nTraceNumPoints,0.0);
    // Array<OneD, NekDouble> e2Fwd_cdot_ncrossdH(nTraceNumPoints,0.0);
    // Array<OneD, NekDouble> e2Bwd_cdot_ncrossdH(nTraceNumPoints,0.0);
    Array<OneD, NekDouble> e3Fwd_cdot_dEe3(nTraceNumPoints, 0.0);
    Array<OneD, NekDouble> e3Bwd_cdot_dEe3(nTraceNumPoints, 0.0);

    Array<OneD, NekDouble> ZimFwd(nTraceNumPoints, 1.0);
    Array<OneD, NekDouble> ZimBwd(nTraceNumPoints, 1.0);

    // Compute  numfluxFwd[dir] = (eFwd^[dir] \cdot n \times e^3) *
    // (YimFwd
    // * EFwd^3 + YimBwd * EBwd^3 ) ComputeNtimesFz(0, Fwd, Bwd,
    // m_YimFwd[0], m_YimBwd[0], numfluxFwd[0], numfluxBwd[0]);
    // ComputeNtimesFz(1, Fwd, Bwd, m_YimFwd[1], m_YimBwd[1],
    // numfluxFwd[1], numfluxBwd[1]);

    // Compute numfluxFwd[2] = eFwd^3 \cdot ( n1e1 \times ( ZimFwd HFwd
    // + ZimBwd HBwd ) ) / 2 {{Z_i}}
    ComputeNtimesF12(Fwd, Bwd, ZimFwd, ZimBwd, ZimFwd, ZimBwd, numfluxFwd[2],
                     numfluxBwd[2]);

    // Compute e1Fwd_cdot_ncrossdE = eFwd[dir] \cdot \alpha n \times n
    // \times [H] / 2 {{YimFwd}} ComputeNtimestimesdFz(0, Fwd, Bwd,
    // m_YimFwd[0], m_YimBwd[0], e1Fwd_cdot_ncrossdH,
    // e1Bwd_cdot_ncrossdH); ComputeNtimestimesdFz(1, Fwd, Bwd,
    // m_YimFwd[1], m_YimBwd[1], e2Fwd_cdot_ncrossdH,
    // e2Bwd_cdot_ncrossdH);

    // Compute  \alpha [E3] * ( 1/2{{Zim1}} + 1/2{{Zim2}} )
    // ComputeNtimestimesdF12(Fwd, Bwd, ZimFwd, ZimBwd, ZimFwd, ZimBwd,
    // e3Fwd_cdot_dEe3, e3Bwd_cdot_dEe3);

    // std::cout << "numfluxFwd = " << RootMeanSquare(e3Fwd_cdot_dEe3)
    // << ", numfluxBwd = " << RootMeanSquare(e3Bwd_cdot_dEe3) <<
    // std::endl;

    // Vmath::Svtvp(nTraceNumPoints, -1.0, numfluxFwd[0], 1,
    // e1Fwd_cdot_ncrossdH, 1, numfluxFwd[0], 1);
    // Vmath::Svtvp(nTraceNumPoints, -1.0, numfluxFwd[1], 1,
    // e2Fwd_cdot_ncrossdH, 1, numfluxFwd[1], 1);
    // Vmath::Svtvp(nTraceNumPoints,  1.0, numfluxFwd[2], 1,
    // e3Fwd_cdot_dEe3, 1, numfluxFwd[2], 1);

    // Vmath::Svtvp(nTraceNumPoints, -1.0, numfluxBwd[0], 1,
    // e1Bwd_cdot_ncrossdH, 1, numfluxBwd[0], 1);
    // Vmath::Svtvp(nTraceNumPoints, -1.0, numfluxBwd[1], 1,
    // e2Bwd_cdot_ncrossdH, 1, numfluxBwd[1], 1);
    // Vmath::Svtvp(nTraceNumPoints,  1.0, numfluxBwd[2], 1,
    // e3Bwd_cdot_dEe3, 1, numfluxBwd[2], 1);
}

// Compute e^3 \cdot ( n1e1 \times ( imFwd EFwd + imBwd EBwd ) ) /
// 2{{Y_i}}
void MMFSystem::ComputeNtimesF12(const Array<OneD, Array<OneD, NekDouble>> &Fwd,
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

Array<OneD, NekDouble> MMFSystem::ComputeLaplacianSphericalCoord(
    const Array<OneD, const NekDouble> &inarray)
{
    int nq = GetNpoints();

    Array<OneD, NekDouble> LBDirect(nq, 0.0);

    Array<OneD, NekDouble> Phi(m_spacedim * nq);
    Array<OneD, NekDouble> Theta(m_spacedim * nq);

    ComputeSphericalTangentVector(Phi, Theta);

    Array<OneD, NekDouble> dudphi(nq), du2dphi2(nq), dudth(nq), du2dth2(nq);

    MMFDirectionalDeriv(Phi, inarray, dudphi);
    MMFDirectionalDeriv(Phi, dudphi, du2dphi2);

    MMFDirectionalDeriv(Theta, inarray, dudth);
    MMFDirectionalDeriv(Theta, dudth, du2dth2);

    Array<OneD, NekDouble> x0(nq);
    Array<OneD, NekDouble> x1(nq);
    Array<OneD, NekDouble> x2(nq);

    m_fields[0]->GetCoords(x0, x1, x2);

    NekDouble sin_varphi, cos_varphi, sin_theta, cos_theta;
    NekDouble xp, yp, zp;
    for (int i = 0; i < nq; i++)
    {
        xp = x0[i];
        yp = x1[i];
        zp = x2[i];

        CartesianToNewSpherical(xp, yp, zp, sin_varphi, cos_varphi, sin_theta,
                                cos_theta);

        if (m_MMFActivation[i])
        {
            LBDirect[i] = (1.0 / sin_theta / sin_theta) * du2dphi2[i] +
                          (cos_theta / sin_theta) * dudth[i] + du2dth2[i];
        }
    }

    return LBDirect;
}

// // Explicit Diffusion Compute the Laplacian operator in the MMF coordinates
// void MMFSystem::WeakDGLaplacian(
//     const Array<OneD, const Array<OneD, NekDouble>> &MF1st,
//     const Array<OneD, const NekDouble> &inarray,
//     Array<OneD, NekDouble> &outarray, const DerivType DType,
//     const int SurfaceLaplacian)
// {
//     int nq        = m_fields[0]->GetTotPoints();
//     int ncoeffs   = m_fields[0]->GetNcoeffs();
//     int nTracePts = m_fields[0]->GetTrace()->GetTotPoints();

//     Array<OneD, NekDouble> qcoeffs(ncoeffs);

//     Array<OneD, Array<OneD, NekDouble>> ufluxFwd(m_shapedim);
//     Array<OneD, Array<OneD, NekDouble>> ufluxBwd(m_shapedim);

//     Array<OneD, Array<OneD, NekDouble>> qfield(m_shapedim);

//     for (int j = 0; j < m_shapedim; ++j)
//     {
//         qfield[j]   = Array<OneD, NekDouble>(nq, 0.0);
//         ufluxFwd[j] = Array<OneD, NekDouble>(nTracePts, 0.0);
//         ufluxBwd[j] = Array<OneD, NekDouble>(nTracePts, 0.0);
//     }

//     // Compute Divergence of moving frames
//     Array<OneD, Array<OneD, NekDouble>> DivMF;
//     ComputeDivMF(DType, MF1st, DivMF);

//     // Compute q_{\eta} and q_{\xi} to bbtain numerical fluxes
//     // Get flux = u* ( e^i \cdot \vec{n} )
//     // GetufluxMMF(0, inarray, ufluxFwd, ufluxBwd);

//     Array<OneD, NekDouble> tmp(nq);
//     Array<OneD, NekDouble> tmpc(ncoeffs);
//     for (int j = 0; j < m_shapedim; ++j)
//     {
//         // Compute L2 = \int ( \partial phi / \partial x_j) u d x
//         m_fields[0]->IProductWRTDirectionalDerivBase(MF1st[j], inarray,
//                                                      qcoeffs);
//         // Compute -L2
//         Vmath::Neg(ncoeffs, qcoeffs, 1);

//         // Compute L = -L2 + \int_{\partial} u* n_x dx
//         // m_fields[0]->AddFwdBwdTraceIntegral(ufluxFwd[j], ufluxBwd[j],
//         qcoeffs);

//         // Add -\int ( \nabla \cdot e^i ) u dx
//         Vmath::Vmul(nq, &inarray[0], 1, &DivMF[j][0], 1, &tmp[0], 1);
//         Vmath::Neg(nq, &tmp[0], 1);

//         m_fields[0]->IProductWRTBase(tmp, tmpc);
//         Vmath::Vadd(ncoeffs, tmpc, 1, qcoeffs, 1, qcoeffs, 1);

//         m_fields[0]->SetPhysState(false);

//         // Compute M^{-1} ( -L2 + \int_{\partial} u* n_x dx )
//         m_fields[0]->MultiplyByElmtInvMass(qcoeffs, qcoeffs);
//         m_fields[0]->BwdTrans(qcoeffs, qfield[j]);
//     }

//     // Anisotropy effect
//     Array<OneD, Array<OneD, NekDouble>> Anicoeff(m_shapedim);
//     for (int j = 0; j < m_shapedim; ++j)
//     {
//         Anicoeff[j] = Array<OneD, NekDouble>(nq, 0.0);
//         for (int k = 0; k < m_spacedim; ++k)
//         {
//             Vmath::Vvtvp(nq, &MF1st[j][k * nq], 1, &MF1st[j][k * nq], 1,
//                          &Anicoeff[j][0], 1, &Anicoeff[j][0], 1);
//         }
//         Vmath::Vsqrt(nq, &Anicoeff[j][0], 1, &Anicoeff[j][0], 1);

//         Vmath::Vdiv(nq, &qfield[j][0], 1, &Anicoeff[j][0], 1, &qfield[j][0],
//         1);
//     }

//     outarray = Array<OneD, NekDouble>(ncoeffs, 0.0);
//     // Compute tmp = \int \nabla \varphi \cdot \vec{q}
//     for (int j = 0; j < m_shapedim; ++j)
//     {
//         m_fields[0]->IProductWRTDirectionalDerivBase(MF1st[j], qfield[j],
//         tmpc); Vmath::Vadd(ncoeffs, tmpc, 1, outarray, 1, outarray, 1);

//         if (SurfaceLaplacian)
//         {
//             Array<OneD, NekDouble> geometricdiv(nq);
//             Array<OneD, NekDouble> velvector(m_spacedim * nq);
//             for (int j = 0; j < m_shapedim; ++j)
//             {
//                 for (int i = 0; i < m_spacedim; ++i)
//                 {
//                     Vmath::Vmul(nq, &qfield[j][0], 1, &MF1st[j][i * nq], 1,
//                                 &velvector[i * nq], 1);
//                 }

//                 geometricdiv = ComputeSurfaceDiv(MF1st[2], velvector);

//                 m_fields[0]->IProductWRTBase(geometricdiv, tmpc);
//                 Vmath::Vadd(ncoeffs, tmpc, 1, outarray, 1, outarray, 1);
//             }
//         }
//     }

//     // Compute u from q_{\eta} and q_{\xi} to obtain numerical fluxes
//     Array<OneD, NekDouble> qflux(nTracePts,0.0);
//     // GetqfluxMMF(MF1st, inarray, qfield, qflux);

//     // Evaulate  <\phi, \hat{F}\cdot n> - outarray[i]
//     Vmath::Neg(ncoeffs, outarray, 1);
//     m_fields[0]->AddTraceIntegral(qflux, outarray);
//     m_fields[0]->SetPhysState(false);

//     // m_fields[0]->MultiplyByElmtInvMass(DivSum, DivSum);
//     // m_fields[0]->BwdTrans(DivSum, outarray);
// }

/**
 * @brief Get the normal velocity for the linear advection equation.
 */
Array<OneD, NekDouble> MMFSystem::GetNormalVelocity(
    const Array<OneD, const Array<OneD, NekDouble>> &velocity)
{
    // Number of trace (interface) points
    int nTracePts = GetTraceNpoints();

    // Auxiliary variable to compute the normal velocity
    Array<OneD, NekDouble> tmp(nTracePts);

    // Reset the normal velocity
    Array<OneD, NekDouble> traceVn(nTracePts, 0.0);

    for (int i = 0; i < m_spacedim; ++i)
    {
        m_fields[0]->ExtractTracePhys(velocity[i], tmp);

        Vmath::Vvtvp(nTracePts, m_traceNormals[i], 1, tmp, 1, traceVn, 1,
                     traceVn, 1);
    }

    return traceVn;
}

Array<OneD, int> MMFSystem::Computedudt(
    const NekDouble uTol, const Array<OneD, const NekDouble> &fieldu,
    const Array<OneD, const NekDouble> &fielduold)
{
    int nq = m_fields[0]->GetTotPoints();

    Array<OneD, NekDouble> diff(nq);
    Vmath::Vsub(nq, fieldu, 1, fielduold, 1, diff, 1);

    Array<OneD, int> outarray(nq, 0);

    for (int i = 0; i < nq; ++i)
    {
        if (fieldu[i] > uTol)
        {
            if (diff[i] > 0)
            {
                outarray[i] = 1;
            }

            // dudt < 0: waveback
            if (diff[i] < 0)
            {
                outarray[i] = -1;
            }
        }
    }

    return outarray;
}

void MMFSystem::ComputedudtHistory(const Array<OneD, const int> &dudt,
                                   const Array<OneD, const NekDouble> &fieldu,
                                   Array<OneD, int> &dudtHistory,
                                   Array<OneD, int> &APindex)
{
    int nq = GetTotPoints();

    for (int i = 0; i < nq; i++)
    {
        // At the Wavefront we add index
        if (dudt[i] > 0 && fieldu[i] > m_uTol)
        {
            APindex[i] = APindex[i] + 1;
        }

        dudtHistory[i] = dudt[i];
    }
}

// Compute Velocity magnitude:
// Input: m_spacedim x nq
// output: nq
Array<OneD, NekDouble> MMFSystem::ComputeVelocityMag(
    const Array<OneD, const NekDouble> &velocity)
{
    int nq = m_fields[0]->GetNpoints();

    Array<OneD, NekDouble> outarray(nq, 0.0);

    for (int k = 0; k < m_spacedim; ++k)
    {
        Vmath::Vvtvp(nq, &velocity[k * nq], 1, &velocity[k * nq], 1,
                     &outarray[0], 1, &outarray[0], 1);
    }
    Vmath::Vsqrt(nq, &outarray[0], 1, &outarray[0], 1);

    return outarray;
}

Array<OneD, NekDouble> MMFSystem::ComputeVelocityMag(
    const Array<OneD, const Array<OneD, NekDouble>> &velocity)
{
    int nq = m_fields[0]->GetNpoints();
    Array<OneD, NekDouble> outarray(nq, 0.0);
    for (int k = 0; k < m_spacedim; ++k)
    {
        Vmath::Vvtvp(nq, &velocity[k][0], 1, &velocity[k][0], 1, &outarray[0],
                     1, &outarray[0], 1);
    }
    Vmath::Vsqrt(nq, &outarray[0], 1, &outarray[0], 1);

    return outarray;
}

void MMFSystem::ComputeMFHHD(const Array<OneD, const int> &Activation,
                             const Array<OneD, const NekDouble> &MF0,
                             Array<OneD, NekDouble> &irrotational,
                             Array<OneD, NekDouble> &incompressible,
                             Array<OneD, NekDouble> &harmonic)
{
    int nq = GetTotPoints();

    Array<OneD, Array<OneD, NekDouble>> velocity(m_spacedim);
    for (int k = 0; k < m_spacedim; ++k)
    {
        velocity[k] = Array<OneD, NekDouble>(nq);
        Vmath::Vcopy(nq, &MF0[k * nq], 1, &velocity[k][0], 1);
    }

    ComputeHHD(Activation, velocity, irrotational, incompressible, harmonic);
}

// velocity = irrotational + incompressible + harmonic
void MMFSystem::ComputeHHD(
    const Array<OneD, const int> &Activation,
    const Array<OneD, const Array<OneD, NekDouble>> &velocity,
    Array<OneD, NekDouble> &irrotational,
    Array<OneD, NekDouble> &incompressible, Array<OneD, NekDouble> &harmonic)
{
    int nq      = GetTotPoints();
    int ncoeffs = m_fields[0]->GetNcoeffs();
    int nvar    = 6;

    irrotational   = Array<OneD, NekDouble>(nvar * nq, 0.0);
    incompressible = Array<OneD, NekDouble>(nvar * nq, 0.0);
    harmonic       = Array<OneD, NekDouble>(nvar * nq, 0.0);

    Array<OneD, NekDouble> velmag(nq, 0.0);
    for (int k = 0; k < m_spacedim; ++k)
    {
        Vmath::Vvtvp(nq, &velocity[k][0], 1, &velocity[k][0], 1, &velmag[0], 1,
                     &velmag[0], 1);
    }
    Vmath::Vsqrt(nq, velmag, 1, velmag, 1);

    // We solve ( \nabla^2 - HHlambda ) Y[i] = rhs [i]
    // inarray = input: \hat{rhs} -> output: \hat{Y}
    // outarray = output: \hat{Y} where \hat = modal coeffs
    Array<OneD, NekDouble> Avgtmp(nq);
    Array<OneD, NekDouble> tmp(nq);
    Array<OneD, NekDouble> check(nq);

    Array<OneD, NekDouble> outarray(nq);

    Array<OneD, NekDouble> tmpcoeff(ncoeffs);
    Array<OneD, NekDouble> tmpgrad(m_spacedim * nq, 0.0);

    // Compute varcoeff for Helmsolver
    StdRegions::VarCoeffMap varcoeff;
    ComputeVarCoeff2D(m_movingframes, varcoeff);

    m_session->LoadParameter("HHDtau", m_HHDtau, 0.001);

    StdRegions::ConstFactorMap factors;
    factors[StdRegions::eFactorTau]    = m_HHDtau;
    factors[StdRegions::eFactorLambda] = 0.0;

    // 1) Solve \nabla^2 D = \nabla \cdot velocity
    tmp              = ComputeCovDiv(velocity, m_movingframes);
    NekDouble DivAvg = -1.0 * AvgInt(tmp);
    Vmath::Sadd(nq, DivAvg, tmp, 1, tmp, 1);

    std::cout << "Div. Average of the HHDized vector = " << DivAvg << std::endl;

    SetBoundaryConditions(0.0);
    m_fields[0]->HelmSolve(tmp, tmpcoeff, factors, varcoeff);
    // m_contField->HelmSolve(tmp, tmpcoeff, factors, varcoeff);
    m_fields[0]->BwdTrans(tmpcoeff, outarray);

    tmpgrad = ComputeCovGrad(outarray, m_movingframes);
    // tmpgrad = ComputeEuclideanGradient(outarray);

    Array<OneD, NekDouble> irrmag(nq, 0.0);
    for (int k = 0; k < m_spacedim; ++k)
    {
        Vmath::Vvtvp(nq, &tmpgrad[k * nq], 1, &tmpgrad[k * nq], 1, &irrmag[0],
                     1, &irrmag[0], 1);
        Vmath::Vcopy(nq, &tmpgrad[k * nq], 1, &irrotational[k * nq], 1);
    }
    Vmath::Vsqrt(nq, irrmag, 1, irrmag, 1);

    check = ComputeCovCurl(irrotational, m_movingframes);

    Vmath::Vcopy(nq, &irrmag[0], 1, &irrotational[3 * nq], 1);
    Vmath::Vcopy(nq, &tmp[0], 1, &irrotational[4 * nq], 1);
    Vmath::Vcopy(nq, &outarray[0], 1, &irrotational[5 * nq], 1);

    // 2) Solve \nabla^2 R = - \nabla \cdot J velocity = \vec{k} \cdot \nabla
    // \times velocity
    tmp               = ComputeCovCurl(velocity, m_movingframes);
    NekDouble CurlAvg = -1.0 * AvgInt(tmp);
    Vmath::Sadd(nq, CurlAvg, tmp, 1, tmp, 1);
    std::cout << "Curl. Average of the HHDized vector = " << CurlAvg
              << std::endl;

    SetBoundaryConditions(0.0);
    m_fields[0]->HelmSolve(tmp, tmpcoeff, factors, varcoeff);
    // m_contField->HelmSolve(tmp, tmpcoeff, factors, varcoeff);
    m_fields[0]->BwdTrans(tmpcoeff, outarray);

    // Compute J \nabla R
    tmpgrad = ComputeCovJGrad(outarray, m_movingframes);
    Array<OneD, NekDouble> incompmag(nq, 0.0);
    for (int k = 0; k < m_spacedim; ++k)
    {
        Vmath::Vvtvp(nq, &tmpgrad[k * nq], 1, &tmpgrad[k * nq], 1,
                     &incompmag[0], 1, &incompmag[0], 1);
        Vmath::Vcopy(nq, &tmpgrad[k * nq], 1, &incompressible[k * nq], 1);
    }
    Vmath::Vsqrt(nq, incompmag, 1, incompmag, 1);

    check = ComputeCovDiv(incompressible, m_movingframes);

    Vmath::Vcopy(nq, &incompmag[0], 1, &incompressible[3 * nq], 1);
    Vmath::Vcopy(nq, &tmp[0], 1, &incompressible[4 * nq], 1);

    NekDouble CurlUavg = -1.0 * AvgInt(outarray);
    Vmath::Sadd(nq, CurlUavg, outarray, 1, outarray, 1);
    Vmath::Vcopy(nq, &outarray[0], 1, &incompressible[5 * nq], 1);

    // 3) harmonic = velocity - \int vel dx - irrotational - incompressible
    Array<OneD, NekDouble> harmomag(nq, 0.0);
    for (int k = 0; k < m_spacedim; ++k)
    {
        for (int i = 0; i < nq; ++i)
        {
            if (Activation[i] > 0)
            {
                harmonic[k * nq + i] = velocity[k][i] -
                                       irrotational[k * nq + i] -
                                       incompressible[k * nq + i];
            }
        }
        Vmath::Vvtvp(nq, &harmonic[k * nq], 1, &harmonic[k * nq], 1,
                     &harmomag[0], 1, &harmomag[0], 1);
    }
    Vmath::Vsqrt(nq, harmomag, 1, harmomag, 1);

    // Check vector Laplacian
    // Vector Laplacian of \vec{h} = \nabla (\nabla \cdot \vec{h}) + \nabla
    // \times ( \nabla \times \vec{h} )
    Array<OneD, NekDouble> vecLapharmonic(m_spacedim * nq);
    Array<OneD, NekDouble> vecLapharmonicmag(nq);

    vecLapharmonic = ComputeVectorLaplacian2D(Activation, harmonic);

    NekDouble vLhx, vLhy, vLhz, vecLaphsum = 0.0;
    int cnt = 0;
    for (int i = 0; i < nq; ++i)
    {
        vLhx = vecLapharmonic[i];
        vLhy = vecLapharmonic[i + nq];
        vLhz = vecLapharmonic[i + 2 * nq];

        // if vector Laplacian is too large, then cancel the construction of
        // harmonic field.
        vecLapharmonicmag[i] = sqrt(vLhx * vLhx + vLhy * vLhy + vLhz * vLhz);

        if (Activation[i] > 0)
        {
            vecLaphsum += vecLapharmonicmag[i] * vecLapharmonicmag[i];
            cnt++;
        }
    }
    vecLaphsum = sqrt(vecLaphsum / cnt);

    Vmath::Vcopy(nq, &harmomag[0], 1, &harmonic[3 * nq], 1);
    Vmath::Vcopy(nq, &vecLapharmonicmag[0], 1, &harmonic[4 * nq], 1);

    // Compute the potential for the harmonic vector: \nabla^2 U = \nabla \cdot
    // \vec{h}
    Array<OneD, NekDouble> Divharmon(nq);
    Divharmon = ComputeCovDiv(harmonic, m_movingframes);

    NekDouble DivharmonAvg = -1.0 * AvgInt(Divharmon);
    Vmath::Sadd(nq, DivharmonAvg, Divharmon, 1, Divharmon, 1);

    std::cout << "Div. Average of the harmonic vector = " << DivharmonAvg
              << std::endl;

    SetBoundaryConditions(0.0);

    m_fields[0]->HelmSolve(Divharmon, tmpcoeff, factors, varcoeff);
    // m_contField->HelmSolve(Divharmon, tmpcoeff,   factors,
    // varcoeff);
    m_fields[0]->BwdTrans(tmpcoeff, outarray);

    NekDouble DivhUAvg = -1.0 * AvgInt(outarray);
    Vmath::Sadd(nq, DivhUAvg, outarray, 1, outarray, 1);
    Vmath::Vcopy(nq, &outarray[0], 1, &harmonic[5 * nq], 1);

    Array<OneD, NekDouble> Divvel(nq);
    Array<OneD, NekDouble> Curlvel(nq);
    Array<OneD, NekDouble> Divirrot(nq);
    Array<OneD, NekDouble> Curlirrot(nq);
    Array<OneD, NekDouble> Divincomp(nq);
    Array<OneD, NekDouble> Curlincomp(nq);
    Array<OneD, NekDouble> Curlharmon(nq);

    Divvel  = ComputeCovDiv(velocity, m_movingframes);
    Curlvel = ComputeCovCurl(velocity, m_movingframes);

    Divirrot  = ComputeCovDiv(irrotational, m_movingframes);
    Curlirrot = ComputeCovCurl(irrotational, m_movingframes);

    Divincomp  = ComputeCovDiv(incompressible, m_movingframes);
    Curlincomp = ComputeCovCurl(incompressible, m_movingframes);

    Curlharmon = ComputeCovCurl(harmonic, m_movingframes);

    std::cout << "vel = irrotational + incompressible + harmonic: "
                 "==================, cnt = "
              << cnt << " / " << nq << std::endl;

    std::cout << "Div: " << RootMeanSquare(Divvel, Activation) << " = "
              << RootMeanSquare(Divirrot, Activation) << " + "
              << RootMeanSquare(Divincomp, Activation) << " + "
              << RootMeanSquare(Divharmon, Activation) << std::endl;

    std::cout << "Curl: " << RootMeanSquare(Curlvel, Activation) << " = "
              << RootMeanSquare(Curlirrot, Activation) << " + "
              << RootMeanSquare(Curlincomp, Activation) << " + "
              << RootMeanSquare(Curlharmon, Activation) << std::endl;

    std::cout << "Vel = " << RootMeanSquare(velmag, Activation)
              << ", irrotational = " << RootMeanSquare(irrmag, Activation)
              << ", incompressible = " << RootMeanSquare(incompmag, Activation)
              << ", harmonic = " << vecLaphsum << std::endl
              << std::endl;
}

// Compute \nabla^2 \vec{h} = \nabla (\nabla \cdot \vec{h} ) - \nabla \times (
// \nabla \times \vec{h} )
Array<OneD, NekDouble> MMFSystem::ComputeVectorLaplacian2D(
    const Array<OneD, const int> &Activation,
    const Array<OneD, const NekDouble> &harmonic)
{
    int nq = GetTotPoints();

    Array<OneD, NekDouble> outarray(m_spacedim * nq);

    Array<OneD, NekDouble> Divharmon(nq);
    Array<OneD, NekDouble> Curlharmon(nq);

    Divharmon  = ComputeCovDiv(harmonic, m_movingframes);
    Curlharmon = ComputeCovCurl(harmonic, m_movingframes);

    std::cout << "Vector Laplacian: Div = " << RootMeanSquare(Divharmon)
              << ", Curl = " << RootMeanSquare(Curlharmon) << std::endl;

    // \nabla (\nabla \cdot \vec{h} )
    Array<OneD, NekDouble> tmp1(m_spacedim * nq, 0.0);
    tmp1 = ComputeCovJGrad(outarray, m_movingframes);

    // \nabla \times ( \nabla \times \vec{h} )
    //       = (\nabla h_c \cdot e^2) e^1 - (\nabla h_c \cdot e^1) e^2
    Array<OneD, NekDouble> dhcdx1(nq);
    Array<OneD, NekDouble> dhcdx2(nq);

    MMFDirectionalDeriv(m_movingframes[1], Curlharmon, dhcdx2);
    MMFDirectionalDeriv(m_movingframes[0], Curlharmon, dhcdx1);

    NekDouble tmp2i;
    for (int k = 0; k < m_spacedim; ++k)
    {
        for (int i = 0; i < nq; ++i)
        {
            if (Activation[i] > 0)
            {
                tmp2i = dhcdx2[i] * m_movingframes[0][i + k * nq] -
                        dhcdx1[i] * m_movingframes[1][i + k * nq];
                outarray[i + k * nq] = tmp1[i + k * nq] - tmp2i;
            }
        }
    }

    return outarray;
}

void MMFSystem::ComputeVarCoeff1D(
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    StdRegions::VarCoeffMap &varcoeff)
{
    int nq = GetTotPoints();

    StdRegions::VarCoeffType MMFCoeffs[3] = {StdRegions::eVarCoeffD00,
                                             StdRegions::eVarCoeffD11,
                                             StdRegions::eVarCoeffD22};

    for (int k = 0; k < m_mfdim; ++k)
    {
        varcoeff[MMFCoeffs[k]] = Array<OneD, NekDouble>(nq, 0.0);
    }

    Array<OneD, NekDouble> tmp(nq, 0.0);
    for (int k = 0; k < m_spacedim; ++k)
    {
        Vmath::Vvtvp(nq, &movingframes[0][k * nq], 1, &movingframes[0][k * nq],
                     1, &tmp[0], 1, &tmp[0], 1);
    }

    Vmath::Vsqrt(nq, &tmp[0], 1, &varcoeff[MMFCoeffs[0]][0], 1);

    std::cout << " ::::: 1D Varcoeff is Successfully Created ::::: "
              << std::endl;
}

void MMFSystem::ComputeVarCoeff2D(
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

        // \| e^k \|^2
        varcoeff[MMFCoeffs[indx + 4]] = Array<OneD, NekDouble>(nq, 0.0);
        tmp                           = Array<OneD, NekDouble>(nq, 0.0);
        for (int i = 0; i < m_spacedim; ++i)
        {
            Vmath::Vvtvp(nq, &movingframes[k][i * nq], 1,
                         &movingframes[k][i * nq], 1, &tmp[0], 1, &tmp[0], 1);
        }

        // Vmath::Vsqrt(nq, &tmp[0], 1, &varcoeff[MMFCoeffs[indx + 4]][0], 1);
        Vmath::Vcopy(nq, &tmp[0], 1, &varcoeff[MMFCoeffs[indx + 4]][0], 1);

    std::cout << "k = " << k << ", m_varcoeff = ( " << RootMeanSquare(varcoeff[MMFCoeffs[indx]])
              << " , " << RootMeanSquare(varcoeff[MMFCoeffs[indx+1]]) << " , "
              << RootMeanSquare(varcoeff[MMFCoeffs[indx+2]]) << " , "
              << RootMeanSquare(varcoeff[MMFCoeffs[indx+3]]) << " , "
              << RootMeanSquare(varcoeff[MMFCoeffs[indx+4]]) << " ) " << std::endl;

    }

    std::cout << " ::::: 2D Varcoeff is Successfully Created ::::: "
              << std::endl << std::endl;
}

void MMFSystem::ComputeVarCoeff2DDxDyDz(
    const Array<OneD, const NekDouble> &epsilon,
    StdRegions::VarCoeffMap &varcoeff)
{
    int nq = GetTotPoints();

    StdRegions::VarCoeffType MMFCoeffs[3] = {StdRegions::eVarCoeffD00,
                                            StdRegions::eVarCoeffD11,
                                            StdRegions::eVarCoeffD22};

    for (int j = 0; j < m_spacedim; ++j)
    {
        varcoeff[MMFCoeffs[j]] = Array<OneD, NekDouble>(nq, epsilon[j]);
    }

    std::cout << "m_varcoeff = ( " << RootMeanSquare(varcoeff[MMFCoeffs[0]])
              << " , " << RootMeanSquare(varcoeff[MMFCoeffs[1]]) << " , "
              << RootMeanSquare(varcoeff[MMFCoeffs[2]]) << " ) " << std::endl;

    std::cout << " ::::: 2D [Dx Dy Dz] Varcoeff is Successfully Created ::::: "
              << std::endl << std::endl;
}



// StdRegions::VarCoeffMap MMFSystem::ComputeVarCoeff2D(
//     const Array<OneD, const Array<OneD, NekDouble>> &movingframes)
// {
//     int nq = GetTotPoints();

//     StdRegions::VarCoeffMap varcoeff;

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
//     for (int k = 0; k < m_mfdim; ++k)
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
//         ComputeDivMF(m_DerivType, movingframes, DivMF);

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

//     return varcoeff;

//     std::cout << " ::::: 2D Varcoeff is Successfully Created ::::: "
//               << std::endl;
// }

int MMFSystem::CountActivated(const Array<OneD, const int> &ActivatedHistory)
{
    int nq  = GetTotPoints();
    int cnt = 0;

    for (int i = 0; i < nq; ++i)
    {
        if (ActivatedHistory[i] > 0)
        {
            cnt++;
        }
    }
    return cnt;
}

Array<OneD, NekDouble> MMFSystem::ElementwiseAverage(
    const Array<OneD, const NekDouble> &inarray)
{
    int nq = m_fields[0]->GetNpoints();

    Array<OneD, NekDouble> outarray(nq, 0.0);

    int index;
    NekDouble Avg, Avgjac;
    Array<OneD, NekDouble> inarrayelemt(nq);
    Array<OneD, NekDouble> jacelemt(nq);
    for (int i = 0; i < m_fields[0]->GetTotPoints(i); ++i)
    {
        Vmath::Fill(nq, 0.0, inarrayelemt, 1);
        Vmath::Fill(nq, 0.0, jacelemt, 1);
        for (int j = 0; j < m_fields[0]->GetTotPoints(i); ++j)
        {
            index               = m_fields[0]->GetPhys_Offset(i) + j;
            inarrayelemt[index] = inarray[index];
            jacelemt[index]     = 1.0;
        }
        Avgjac = m_fields[0]->PhysIntegral(jacelemt);
        Avg    = (m_fields[0]->PhysIntegral(inarrayelemt)) / Avgjac;

        for (int j = 0; j < m_fields[0]->GetTotPoints(i); ++j)
        {
            index           = m_fields[0]->GetPhys_Offset(i) + j;
            outarray[index] = Avg;
        }
    }

    return outarray;
}

void MMFSystem::GenerateHHDPlot(const int nstep)
{
    int nq     = m_fields[0]->GetNpoints();
    int ncoeff = m_fields[0]->GetNcoeffs();
    int nvar   = 12;

    std::cout << "Start: Computing and Plotting HHD" << std::endl;

    std::vector<std::string> variables(nvar);
    variables[0]  = "u";
    variables[1]  = "ex1";
    variables[2]  = "ey1";
    variables[3]  = "ez1";
    variables[4]  = "AngleMF";
    variables[5]  = "w211";
    variables[6]  = "w212";
    variables[7]  = "R2121";
    variables[8]  = "Dive1";
    variables[9]  = "Curle1";
    variables[10] = "MeanC";
    variables[11] = "GaussC";

    std::cout << "Import index = " << nstep << std::endl;

    std::string processfile = m_sessionName + "_CNT_" +
                              boost::lexical_cast<std::string>(nstep) + ".chk";

    Array<OneD, Array<OneD, NekDouble>> tmpc(nvar);
    Array<OneD, Array<OneD, NekDouble>> tmp(nvar);
    for (int i = 0; i < nvar; ++i)
    {
        tmpc[i] = Array<OneD, NekDouble>(ncoeff);
        tmp[i]  = Array<OneD, NekDouble>(nq);
    }

    EquationSystem::ImportFld(processfile, variables, tmpc);

    for (int i = 0; i < nvar; ++i)
    {
        m_fields[0]->BwdTrans(tmpc[i], tmp[i]);
    }

    Array<OneD, NekDouble> MFvec(m_spacedim * nq);
    Array<OneD, NekDouble> MFmag(nq, 0.0);
    for (int i = 0; i < m_spacedim; ++i)
    {
        Vmath::Vcopy(nq, &tmp[i + 1][0], 1, &MFvec[i * nq], 1);
    }

    if (m_surfaceType == eSphere)
    {
        // Intermediate process of MFvec....
        std::cout << "Processing of MFvec is intiated" << std::endl;
        Array<OneD, NekDouble> x0(nq);
        Array<OneD, NekDouble> x1(nq);
        Array<OneD, NekDouble> x2(nq);

        m_fields[0]->GetCoords(x0, x1, x2);

        int indexj, Ncnt       = 0;
        NekDouble zpavg, zpTol = 5.0;
        for (int i = 0; i < m_fields[0]->GetExpSize(); ++i)
        {
            zpavg = 0.0;
            for (int j = 0; j < m_fields[0]->GetTotPoints(i); ++j)
            {
                indexj = m_fields[0]->GetPhys_Offset(i) + j;
                zpavg += x2[indexj];
            }
            zpavg = zpavg / m_fields[0]->GetTotPoints(i);

            if (zpavg < zpTol)
            {
                for (int j = 0; j < m_fields[0]->GetTotPoints(i); ++j)
                {
                    indexj                 = m_fields[0]->GetPhys_Offset(i) + j;
                    MFvec[indexj]          = 0.0;
                    MFvec[indexj + nq]     = 0.0;
                    MFvec[indexj + 2 * nq] = 0.0;
                }
                Ncnt++;
            }
        }

        std::cout << "Ncnt = " << Ncnt << " has been cancelled " << std::endl;
    }

    for (int i = 0; i < m_spacedim; ++i)
    {
        Vmath::Vvtvp(nq, &MFvec[i * nq], 1, &MFvec[i * nq], 1, &MFmag[0], 1,
                     &MFmag[0], 1);
    }

    NekDouble MFTol = 0.0001;
    Array<OneD, int> Activation(nq, 0);
    for (int i = 0; i < nq; ++i)
    {
        if (MFmag[i] > MFTol)
        {
            Activation[i] = 1;
        }
    }

    Array<OneD, NekDouble> irrotational;
    Array<OneD, NekDouble> incompressible;
    Array<OneD, NekDouble> harmonic;

    // ComputeMFHHD and PlotHHD
    ComputeMFHHD(Activation, MFvec, irrotational, incompressible, harmonic);
    PlotHHD(MFvec, irrotational, incompressible, harmonic, nstep);

    std::cout << "End: Computing and Plotting HHD" << std::endl;

    ASSERTL0(m_spacedim == 4, "End of programm......");
}

Array<OneD, NekDouble> MMFSystem::ComputeRossyHaurtizVelocity()
{
    int nq = m_fields[0]->GetNpoints();

    Array<OneD, NekDouble> outarray(m_spacedim * nq, 0.0);

    NekDouble SecondToDay = 60.0 * 60.0 * 24.0;
    NekDouble rad_earth   = 6.37122 * 1000000;
    // NekDouble Omegams = 7.292 * 0.00001;

    // Nondimensionalized coeffs.
    // NekDouble gms = 9.80616;
    // NekDouble m_g = (gms * SecondToDay * SecondToDay) / rad_earth;

    // NekDouble m_Omega = Omegams * SecondToDay;

    NekDouble m_H0 = 8000.0;
    m_H0           = m_H0 / rad_earth;

    NekDouble m_W = 7.848 * 0.000001 * SecondToDay;
    NekDouble m_K = 7.848 * 0.000001 * SecondToDay;

    Array<OneD, NekDouble> x0(nq);
    Array<OneD, NekDouble> x1(nq);
    Array<OneD, NekDouble> x2(nq);

    m_fields[0]->GetCoords(x0, x1, x2);

    NekDouble sin_varphi, cos_varphi, sin_theta, cos_theta;
    NekDouble x0j, x1j, x2j;
    NekDouble tmp1, tmp2, cos4phi, sin4phi;
    NekDouble uphi, uth;
    for (int i = 0; i < nq; i++)
    {
        x0j = x0[i];
        x1j = x1[i];
        x2j = x2[i];

        CartesianToNewSpherical(x0j, x1j, x2j, sin_varphi, cos_varphi,
                                sin_theta, cos_theta);

        tmp1    = 2.0 * sin_varphi * cos_varphi;
        tmp2    = 2.0 * cos_varphi * cos_varphi - 1.0;
        sin4phi = 2.0 * tmp1 * tmp2;

        tmp1    = 2.0 * cos_varphi * cos_varphi - 1.0;
        cos4phi = 2.0 * tmp1 * tmp1 - 1.0;

        uphi = m_W * sin_theta +
               m_K * sin_theta * sin_theta * sin_theta *
                   (4 * cos_theta * cos_theta - sin_theta * sin_theta) *
                   cos4phi;
        uth = -4.0 * m_K * sin_theta * sin_theta * sin_theta * cos_theta *
              sin4phi;

        outarray[i] = -1.0 * uphi * sin_varphi - uth * cos_theta * cos_varphi;
        outarray[i + nq]     = uphi * cos_varphi - uth * cos_theta * sin_varphi;
        outarray[i + 2 * nq] = uth * sin_theta;
    }

    return outarray;
}

void MMFSystem::SphericalToMovingFrames(
    const Array<OneD, const NekDouble> &inarrayth,
    const Array<OneD, const NekDouble> &inarrayphi,
    Array<OneD, Array<OneD, NekDouble>> &physfield)
{
    int nq = GetNpoints();

    physfield = Array<OneD, Array<OneD, NekDouble>>(m_shapedim);
    // Get the variables in physical space already in physical space
    for (int i = 0; i < m_shapedim; ++i)
    {
        physfield[i] = Array<OneD, NekDouble>(nq);
    }

    Array<OneD, Array<OneD, NekDouble>> ThetacdotMF(m_shapedim);
    Array<OneD, Array<OneD, NekDouble>> PhicdotMF(m_shapedim);

    ComputeMFcdotSphericalCoord(ThetacdotMF, PhicdotMF);

    Vmath::Vmul(nq, &inarrayth[0], 1, &ThetacdotMF[0][0], 1, &physfield[0][0],
                1);
    Vmath::Vvtvp(nq, &inarrayphi[0], 1, &PhicdotMF[0][0], 1, &physfield[0][0],
                 1, &physfield[0][0], 1);

    Vmath::Vmul(nq, &inarrayth[0], 1, &ThetacdotMF[1][0], 1, &physfield[1][0],
                1);
    Vmath::Vvtvp(nq, &inarrayphi[0], 1, &PhicdotMF[1][0], 1, &physfield[1][0],
                 1, &physfield[1][0], 1);
}

void MMFSystem::SphericalToEuclidean(
    const Array<OneD, const NekDouble> &inarrayth,
    const Array<OneD, const NekDouble> &inarrayphi,
    Array<OneD, Array<OneD, NekDouble>> &velocity)
{
    int nq = GetNpoints();

    Array<OneD, Array<OneD, NekDouble>> physfield(m_shapedim);
    // Get the variables in physical space already in physical space
    for (int i = 0; i < m_shapedim; ++i)
    {
        physfield[i] = Array<OneD, NekDouble>(nq);
    }

    Array<OneD, Array<OneD, NekDouble>> ThetacdotMF(m_shapedim);
    Array<OneD, Array<OneD, NekDouble>> PhicdotMF(m_shapedim);

    ComputeMFcdotSphericalCoord(ThetacdotMF, PhicdotMF);

    Vmath::Vmul(nq, &inarrayth[0], 1, &ThetacdotMF[0][0], 1, &physfield[0][0],
                1);
    Vmath::Vvtvp(nq, &inarrayphi[0], 1, &PhicdotMF[0][0], 1, &physfield[0][0],
                 1, &physfield[0][0], 1);

    Vmath::Vmul(nq, &inarrayth[0], 1, &ThetacdotMF[1][0], 1, &physfield[1][0],
                1);
    Vmath::Vvtvp(nq, &inarrayphi[0], 1, &PhicdotMF[1][0], 1, &physfield[1][0],
                 1, &physfield[1][0], 1);

    velocity = Array<OneD, Array<OneD, NekDouble>>(m_spacedim);
    for (int k = 0; k < m_spacedim; ++k)
    {
        velocity[k] = Array<OneD, NekDouble>(nq);
    }

    for (int i = 0; i < nq; ++i)
    {
        velocity[0][i] = physfield[0][i] * m_movingframes[0][i] +
                         physfield[1][i] * m_movingframes[1][i];
        velocity[1][i] = physfield[0][i] * m_movingframes[0][i + nq] +
                         physfield[1][i] * m_movingframes[1][i + nq];
        velocity[2][i] = physfield[0][i] * m_movingframes[0][i + 2 * nq] +
                         physfield[1][i] * m_movingframes[1][i + 2 * nq];
    }
}

void MMFSystem::WeakDGMMFDiff1D(
    const Array<OneD, const Array<OneD, NekDouble>> &MF1st,
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &qfield,
    Array<OneD, NekDouble> &outarray, const NekDouble time)
{
    boost::ignore_unused(time);

    int nvariables = m_fields.size();
    int nq         = m_fields[0]->GetTotPoints();
    int ncoeffs    = m_fields[0]->GetNcoeffs();
    int nTracePts  = m_fields[0]->GetTrace()->GetTotPoints();

    Array<OneD, NekDouble> qcoeffs(ncoeffs);

    Array<OneD, Array<OneD, NekDouble>> ufluxFwd(m_expdim);
    Array<OneD, Array<OneD, NekDouble>> ufluxBwd(m_expdim);

    qfield = Array<OneD, Array<OneD, NekDouble>>(m_expdim);

    for (int j = 0; j < m_expdim; ++j)
    {
        qfield[j]   = Array<OneD, NekDouble>(nq, 0.0);
        ufluxFwd[j] = Array<OneD, NekDouble>(nTracePts, 0.0);
        ufluxBwd[j] = Array<OneD, NekDouble>(nTracePts, 0.0);
    }

    // Array<OneD, Array<OneD, NekDouble>> DivMF;
    // ComputeDivMF(eCovariant, MF1st, DivMF);

    // Compute q_{\eta} and q_{\xi} to bbtain numerical fluxes
    // Get flux = u* ( e^i \cdot \vec{n} )
    // GetufluxMMF(MF1st, inarray, ufluxFwd, ufluxBwd);
    Array<OneD, Array<OneD, NekDouble>> ncdotMFFwd;
    Array<OneD, Array<OneD, NekDouble>> ncdotMFBwd;

    Array<OneD, NekDouble> uFwd(nTracePts), uBwd(nTracePts);
    for (int var = 0; var < nvariables; ++var)
    {
        m_fields[0]->GetFwdBwdTracePhys(inarray, uFwd, uBwd);

        // Boundary Treatment for Dirichlet and Neumann
        DiffusionScalarBoundary(uFwd, uBwd);

        // Return u* * ( e^i \cdot \vec{n} )
        for (int j = 0; j < m_expdim; ++j)
        {
            // Vmath::Vmul(nTracePts, &uFwd[0], 1, &m_ncdotMFFwd[j][0], 1,
            // &ufluxFwd[j][0], 1); Vmath::Vmul(nTracePts, &uFwd[0], 1,
            // &m_ncdotMFBwd[j][0], 1, &ufluxBwd[j][0], 1);
            Vmath::Vcopy(nTracePts, &uFwd[0], 1, &ufluxFwd[j][0], 1);
            Vmath::Vcopy(nTracePts, &uFwd[0], 1, &ufluxBwd[j][0], 1);
        }
    }

    Array<OneD, NekDouble> tmp(nq);
    Array<OneD, NekDouble> tmpc(ncoeffs);
    for (int j = 0; j < m_expdim; ++j)
    {
        // Compute L2 = \int ( \partial phi / \partial x_j) u d x
        m_fields[0]->IProductWRTDirectionalDerivBase(MF1st[j], inarray,
                                                     qcoeffs);
        // Compute -L2
        Vmath::Neg(ncoeffs, qcoeffs, 1);

        // Compute L = -L2 + \int_{\partial} u* n_x dx
        m_fields[0]->AddFwdBwdTraceIntegral(ufluxFwd[j], ufluxBwd[j], qcoeffs);

        // Compute M^{-1} ( -L2 + \int_{\partial} u* n_x dx )
        m_fields[0]->MultiplyByElmtInvMass(qcoeffs, qcoeffs);
        m_fields[0]->BwdTrans(qcoeffs, qfield[j]);

        // Add -\int ( \nabla \cdot e^i ) u dx
        // Vmath::Vmul(nq, &inarray[0], 1, &m_DivMF[j][0], 1, &tmp[0], 1);
        // Vmath::Neg(nq, &tmp[0], 1);

        Vmath::Vadd(nq, tmp, 1, qfield[j], 1, qfield[j], 1);
    }

    Array<OneD, NekDouble> x0(nq);
    Array<OneD, NekDouble> x1(nq);
    Array<OneD, NekDouble> x2(nq);

    m_fields[0]->GetCoords(x0, x1, x2);

    Array<OneD, NekDouble> qexact(nq);
    Array<OneD, NekDouble> qdiff(nq);
    for (int i = 0; i < nq; ++i)
    {
        qexact[i] = exp(-1.0 * time) * m_pi * cos(m_pi * x0[i]);
    }
    Vmath::Vsub(nq, qfield[0], 1, qexact, 1, qdiff, 1);
    std::cout << "q = " << RootMeanSquare(qfield[0])
              << ", exactq = " << RootMeanSquare(qexact)
              << ", error = " << RootMeanSquare(qdiff) << std::endl;

    // Anisotropy effect
    // Array<OneD, Array<OneD, NekDouble>> Anicoeff(m_expdim);
    // for (int j = 0; j < m_expdim; ++j)
    // {
    //     Anicoeff[j] = Array<OneD, NekDouble>(nq, 0.0);
    //     for (int k = 0; k < m_spacedim; ++k)
    //     {
    //         Vmath::Vvtvp(nq, &MF1st[j][k * nq], 1, &MF1st[j][k * nq], 1,
    //                      &Anicoeff[j][0], 1, &Anicoeff[j][0], 1);
    //     }
    //     Vmath::Vsqrt(nq, &Anicoeff[j][0], 1, &Anicoeff[j][0], 1);

    //     Vmath::Vdiv(nq, &qfield[j][0], 1, &Anicoeff[j][0], 1,
    //     &qfield[j][0], 1);
    // }

    Array<OneD, NekDouble> DivSum(ncoeffs, 0.0);
    // Compute tmp = \int \nabla \varphi \cdot \vec{q}
    for (int j = 0; j < m_expdim; ++j)
    {
        m_fields[0]->IProductWRTDirectionalDerivBase(MF1st[j], qfield[j], tmpc);
        Vmath::Vadd(ncoeffs, tmpc, 1, DivSum, 1, DivSum, 1);
    }

    // Compute u from q_{\eta} and q_{\xi} to obtain numerical fluxes
    Array<OneD, NekDouble> qflux(nTracePts);
    GetqfluxMMF(MF1st, inarray, qfield, qflux);

    // Evaulate  <\phi, \hat{F}\cdot n> - outarray[i]
    Vmath::Neg(ncoeffs, DivSum, 1);
    m_fields[0]->AddTraceIntegral(qflux, DivSum);
    m_fields[0]->SetPhysState(false);

    m_fields[0]->MultiplyByElmtInvMass(DivSum, DivSum);
    m_fields[0]->BwdTrans(DivSum, outarray);
}

void MMFSystem::WeakDGMMFDiff(
    const Array<OneD, const Array<OneD, NekDouble>> &MF1st,
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &qfield,
    Array<OneD, NekDouble> &outarray, const NekDouble time)
{
    boost::ignore_unused(time);

    int nvariables = m_fields.size();
    int nq         = m_fields[0]->GetTotPoints();
    int ncoeffs    = m_fields[0]->GetNcoeffs();
    int nTracePts  = m_fields[0]->GetTrace()->GetTotPoints();

    Array<OneD, NekDouble> qcoeffs(ncoeffs);

    Array<OneD, Array<OneD, NekDouble>> ufluxFwd(m_expdim);
    Array<OneD, Array<OneD, NekDouble>> ufluxBwd(m_expdim);

    qfield = Array<OneD, Array<OneD, NekDouble>>(m_expdim);

    for (int j = 0; j < m_expdim; ++j)
    {
        qfield[j]   = Array<OneD, NekDouble>(nq, 0.0);
        ufluxFwd[j] = Array<OneD, NekDouble>(nTracePts, 0.0);
        ufluxBwd[j] = Array<OneD, NekDouble>(nTracePts, 0.0);
    }

    // Array<OneD, Array<OneD, NekDouble>> DivMF;
    // ComputeDivMF(eCovariant, MF1st, DivMF);
    // lux = u* ( e^i \cdot \vec{n} )
    // GetufluxMMF(MF1st, inarray, ufluxFwd, ufluxBwd);
    Array<OneD, Array<OneD, NekDouble>> ncdotMFFwd;
    Array<OneD, Array<OneD, NekDouble>> ncdotMFBwd;

    Array<OneD, NekDouble> uFwd(nTracePts), uBwd(nTracePts);
    for (int var = 0; var < nvariables; ++var)
    {
        m_fields[0]->GetFwdBwdTracePhys(inarray, uFwd, uBwd);

        // Boundary Treatment for Dirichlet and Neumann
        DiffusionScalarBoundary(uFwd, uBwd);

        // Return u* * ( e^i \cdot \vec{n} )
        for (int j = 0; j < m_expdim; ++j)
        {
            Vmath::Vmul(nTracePts, &uFwd[0], 1, &m_ncdotMFFwd[j][0], 1,
                        &ufluxFwd[j][0], 1);
            Vmath::Vmul(nTracePts, &uFwd[0], 1, &m_ncdotMFBwd[j][0], 1,
                        &ufluxBwd[j][0], 1);
        }
    }

    Array<OneD, NekDouble> tmp(nq);
    Array<OneD, NekDouble> tmpc(ncoeffs);
    for (int j = 0; j < m_expdim; ++j)
    {
        // Compute L2 = \int ( \partial phi / \partial x_j) u d x
        m_fields[0]->IProductWRTDirectionalDerivBase(MF1st[j], inarray,
                                                     qcoeffs);
        // Compute -L2
        Vmath::Neg(ncoeffs, qcoeffs, 1);

        // Compute L = -L2 + \int_{\partial} u* n_x dx
        m_fields[0]->AddFwdBwdTraceIntegral(ufluxFwd[j], ufluxBwd[j], qcoeffs);

        // Compute M^{-1} ( -L2 + \int_{\partial} u* n_x dx )
        m_fields[0]->MultiplyByElmtInvMass(qcoeffs, qcoeffs);
        m_fields[0]->BwdTrans(qcoeffs, qfield[j]);

        // Add -\int ( \nabla \cdot e^i ) u dx
        Vmath::Vmul(nq, &inarray[0], 1, &m_DivMF[j][0], 1, &tmp[0], 1);
        Vmath::Neg(nq, &tmp[0], 1);

        Vmath::Vadd(nq, tmp, 1, qfield[j], 1, qfield[j], 1);
    }

    // Anisotropy effect
    Array<OneD, Array<OneD, NekDouble>> Anicoeff(m_expdim);
    for (int j = 0; j < m_expdim; ++j)
    {
        Anicoeff[j] = Array<OneD, NekDouble>(nq, 0.0);
        for (int k = 0; k < m_spacedim; ++k)
        {
            Vmath::Vvtvp(nq, &MF1st[j][k * nq], 1, &MF1st[j][k * nq], 1,
                         &Anicoeff[j][0], 1, &Anicoeff[j][0], 1);
        }
        Vmath::Vsqrt(nq, &Anicoeff[j][0], 1, &Anicoeff[j][0], 1);

        Vmath::Vdiv(nq, &qfield[j][0], 1, &Anicoeff[j][0], 1, &qfield[j][0], 1);
    }

    Array<OneD, NekDouble> DivSum(ncoeffs, 0.0);
    // Compute tmp = \int \nabla \varphi \cdot \vec{q}
    for (int j = 0; j < m_expdim; ++j)
    {
        m_fields[0]->IProductWRTDirectionalDerivBase(MF1st[j], qfield[j], tmpc);
        Vmath::Vadd(ncoeffs, tmpc, 1, DivSum, 1, DivSum, 1);
    }

    // Compute u from q_{\eta} and q_{\xi} to obtain numerical fluxes
    Array<OneD, NekDouble> qflux(nTracePts);
    GetqfluxMMF(MF1st, inarray, qfield, qflux);
    std::cout << "qflux = " << RootMeanSquare(qflux) << std::endl;

    // Evaulate  <\phi, \hat{F}\cdot n> - outarray[i]
    Vmath::Neg(ncoeffs, DivSum, 1);
    m_fields[0]->AddTraceIntegral(qflux, DivSum);
    m_fields[0]->SetPhysState(false);

    m_fields[0]->MultiplyByElmtInvMass(DivSum, DivSum);
    m_fields[0]->BwdTrans(DivSum, outarray);
}

void MMFSystem::GetufluxMMF(
    const Array<OneD, const Array<OneD, NekDouble>> &MF1st,
    const Array<OneD, const NekDouble> &ufield,
    Array<OneD, Array<OneD, NekDouble>> &ufluxFwd,
    Array<OneD, Array<OneD, NekDouble>> &ufluxBwd)
{
    boost::ignore_unused(MF1st);

    int nTracePts  = m_fields[0]->GetTrace()->GetTotPoints();
    int nvariables = m_fields.size();

    Array<OneD, Array<OneD, NekDouble>> ncdotMFFwd;
    Array<OneD, Array<OneD, NekDouble>> ncdotMFBwd;

    // ComputencdotMF(MF1st, ncdotMFFwd, ncdotMFBwd);

    Array<OneD, NekDouble> uFwd(nTracePts), uBwd(nTracePts);
    for (int var = 0; var < nvariables; ++var)
    {
        m_fields[0]->GetFwdBwdTracePhys(ufield, uFwd, uBwd);

        // Boundary Treatment for Dirichlet and Neumann
        DiffusionScalarBoundary(uFwd, uBwd);

        // Return u* * ( e^i \cdot \vec{n} )
        for (int j = 0; j < m_expdim; ++j)
        {
            Vmath::Vmul(nTracePts, &uFwd[0], 1, &m_ncdotMFFwd[j][0], 1,
                        &ufluxFwd[j][0], 1);
            Vmath::Vmul(nTracePts, &uFwd[0], 1, &m_ncdotMFBwd[j][0], 1,
                        &ufluxBwd[j][0], 1);
        }
    }
}

// Compute q* \cdot \vec{n}^+
void MMFSystem::GetqfluxMMF(
    const Array<OneD, const Array<OneD, NekDouble>> &MF1st,
    const Array<OneD, const NekDouble> &ufield,
    Array<OneD, Array<OneD, NekDouble>> &qfield, Array<OneD, NekDouble> &qflux)
{
    int j, k, var;
    int nq         = GetTotPoints();
    int nTracePts  = m_fields[0]->GetTrace()->GetTotPoints();
    int nvariables = m_fields.size();

    Array<OneD, NekDouble> qftmp(nTracePts, 0.0);

    Array<OneD, NekDouble> uFwd(nTracePts);
    Array<OneD, NekDouble> uBwd(nTracePts);

    Array<OneD, Array<OneD, NekDouble>> qFwd(m_expdim);
    Array<OneD, Array<OneD, NekDouble>> qBwd(m_expdim);
    for (var = 0; var < nvariables; ++var)
    {
        // Compute q^* \cdot n^+
        m_fields[0]->GetFwdBwdTracePhys(ufield, uFwd, uBwd);

        for (j = 0; j < m_expdim; ++j)
        {
            qFwd[j] = Array<OneD, NekDouble>(nTracePts, 0.0);
            qBwd[j] = Array<OneD, NekDouble>(nTracePts, 0.0);
            m_fields[0]->GetFwdBwdTracePhys(qfield[j], qFwd[j], qBwd[j]);
        }

        Array<OneD, Array<OneD, NekDouble>> qVec(m_spacedim);
        Array<OneD, Array<OneD, NekDouble>> qFwdVec(m_spacedim);
        Array<OneD, Array<OneD, NekDouble>> qBwdVec(m_spacedim);
        for (int k = 0; k < m_spacedim; ++k)
        {
            qVec[k]    = Array<OneD, NekDouble>(nq, 0.0);
            qFwdVec[k] = Array<OneD, NekDouble>(nTracePts, 0.0);
            qBwdVec[k] = Array<OneD, NekDouble>(nTracePts, 0.0);
        }

        //   DeriveVector(qfield[0], qfield[1], MF1st[0], MF1st[1], qVec);
        for (int i = 0; i < nq; ++i)
        {
            qVec[0][i] = qfield[0][i] * MF1st[0][i];
            qVec[1][i] = qfield[0][i] * MF1st[0][i + nq];
            qVec[2][i] = qfield[0][i] * MF1st[0][i + 2 * nq];
        }

        for (int i = 0; i < m_spacedim; ++i)
        {
            m_fields[0]->GetFwdBwdTracePhys(qVec[i], qFwdVec[i], qBwdVec[i]);
        }

        // Boundary Treatment for Dirichlet and Neumann
        DiffusionVectorBoundary(MF1st, qFwd, qFwdVec, qBwd, qBwdVec);

        NekDouble qAver, qJump;
        for (k = 0; k < nTracePts; ++k)
        {
            qAver = 0.0;
            qJump = 0.0;
            for (int j = 0; j < m_spacedim; ++j)
            {
                qAver += 0.5 * (qFwdVec[j][k] + qBwdVec[j][k]) *
                         m_traceNormals[j][k];
                qJump += 0.5 * (qFwdVec[j][k] - qBwdVec[j][k]) *
                         m_traceNormals[j][k];
            }
            // Compute [[ u ]] \cdot \vec{n}^+
            qflux[k] = qAver - qJump + 0.5 * (uFwd[k] - uBwd[k]);
        }
    }
}

// Boundary conditions for uflux
void MMFSystem::DiffusionScalarBoundary(const Array<OneD, const NekDouble> &Fwd,
                                        Array<OneD, NekDouble> &Bwd)
{
    int i, e, id1, id2;

    // Number of boundary regions
    int nBndEdgePts, nBndEdges;
    int cnt         = 0;
    int nBndRegions = m_fields[0]->GetBndCondExpansions().size();

    for (i = 0; i < nBndRegions; ++i)
    {
        // Number of boundary expansion related to that region
        nBndEdges = m_fields[0]->GetBndCondExpansions()[i]->GetExpSize();

        // Weakly impose boundary conditions by modifying flux values
        for (e = 0; e < nBndEdges; ++e)
        {
            // nBndEdgePts = m_fields[0]
            //                   ->GetBndCondExpansions()[i]
            //                   ->GetExp(e)
            //                   ->GetTotPoints();

            // id1 = m_fields[0]->GetBndCondExpansions()[i]->GetPhys_Offset(e);
            // id2 = m_fields[0]->GetTrace()->GetPhys_Offset(
            //     m_fields[0]->GetTraceMap()->GetBndCondTraceToGlobalTraceMap(
            //         cnt++));

            nBndEdgePts =
                m_fields[0]->GetBndCondExpansions()[i]->GetExp(e)->GetNumPoints(
                    0);
            id1 = m_fields[0]->GetBndCondExpansions()[i]->GetPhys_Offset(e);
            id2 = m_fields[0]->GetTrace()->GetPhys_Offset(
                m_fields[0]->GetTraceMap()->GetBndCondIDToGlobalTraceID(cnt++));

            if (m_fields[0]
                    ->GetBndConditions()[i]
                    ->GetBoundaryConditionType() == SpatialDomains::eDirichlet)
            {
                // For Dirichlet boundary condition: uBwd = g_D
                Vmath::Vcopy(
                    nBndEdgePts,
                    &(m_fields[0]->GetBndCondExpansions()[i]->GetPhys())[id1],
                    1, &Bwd[id2], 1);
            }

            else if ((m_fields[0]->GetBndConditions()[i])
                         ->GetBoundaryConditionType() ==
                     SpatialDomains::eNeumann)
            {
                // For Neumann boundary condition: uBwd = u+
                Vmath::Vcopy(nBndEdgePts, &Fwd[id2], 1, &Bwd[id2], 1);
            }
        }
    }
}

// Boundary conditions for uflux
void MMFSystem::DiffusionVectorBoundary(
    const Array<OneD, const Array<OneD, NekDouble>> &MF1st,
    const Array<OneD, const Array<OneD, NekDouble>> &qFwd,
    const Array<OneD, const Array<OneD, NekDouble>> &qFwdVec,
    Array<OneD, Array<OneD, NekDouble>> &qBwd,
    Array<OneD, Array<OneD, NekDouble>> &qBwdVec)

{
    boost::ignore_unused(qFwd);

    int id2, cnt = 0;
    int nBndEdgePts, nBndEdges;
    int nTracePts = m_fields[0]->GetTrace()->GetTotPoints();

    // Number of boundary regions
    int nBndRegions = m_fields[0]->GetBndCondExpansions().size();

    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> MFtraceFwd;
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> MFtraceBwd;
    ComputeMFtrace(MF1st, MFtraceFwd, MFtraceBwd);

    for (int i = 0; i < nBndRegions; ++i)
    {
        // Number of boundary expansion related to that region
        nBndEdges = m_fields[0]->GetBndCondExpansions()[i]->GetExpSize();

        // Weakly impose boundary conditions by modifying flux values
        for (int e = 0; e < nBndEdges; ++e)
        {
            // nBndEdgePts = m_fields[0]
            //                   ->GetBndCondExpansions()[i]
            //                   ->GetExp(e)
            //                   ->GetTotPoints();

            // id2 = m_fields[0]->GetTrace()->GetPhys_Offset(
            //     m_fields[0]->GetTraceMap()->GetBndCondTraceToGlobalTraceMap(
            //         cnt++));

            nBndEdgePts =
                m_fields[0]->GetBndCondExpansions()[i]->GetExp(e)->GetNumPoints(
                    0);

            id2 = m_fields[0]->GetTrace()->GetPhys_Offset(
                m_fields[0]->GetTraceMap()->GetBndCondIDToGlobalTraceID(cnt++));

            // For Dirichlet boundary condition: qBwd = qFwd
            if (m_fields[0]
                    ->GetBndConditions()[i]
                    ->GetBoundaryConditionType() == SpatialDomains::eDirichlet)
            {
                Array<OneD, NekDouble> tmp0(nBndEdgePts, 0.0);
                Array<OneD, NekDouble> tmp1(nBndEdgePts, 0.0);
                for (int k = 0; k < m_spacedim; ++k)
                {
                    Vmath::Vvtvp(nBndEdgePts, &qFwdVec[k][id2], 1,
                                 &MFtraceBwd[0][k][id2], 1, &tmp0[0], 1,
                                 &tmp0[0], 1);

                    Vmath::Vvtvp(nBndEdgePts, &qFwdVec[k][id2], 1,
                                 &MFtraceBwd[1][k][id2], 1, &tmp1[0], 1,
                                 &tmp1[0], 1);

                    Vmath::Vcopy(nBndEdgePts, &qFwdVec[k][id2], 1,
                                 &qBwdVec[k][id2], 1);
                }

                Vmath::Vcopy(nBndEdgePts, &tmp0[0], 1, &qBwd[0][id2], 1);
                // Vmath::Vcopy(nBndEdgePts, &tmp1[0], 1, &qBwd[1][id2], 1);
            }

            // For Neumann boundary condition: qBwd = g_N \vec{n}
            // Suppose g_N=0. Create an orthonormal vector to \vec{n}. Then
            // let \vec{q}_Fwd = \vec{n}^{\perp}. Find corresponding q_Bwd
            // components.
            else if ((m_fields[0]->GetBndConditions()[i])
                         ->GetBoundaryConditionType() ==
                     SpatialDomains::eNeumann)
            {
                Array<OneD, Array<OneD, NekDouble>> ZeroNvec(m_spacedim);
                for (int i = 0; i < m_spacedim; ++i)
                {
                    ZeroNvec[i] = Array<OneD, NekDouble>(nTracePts, 0.0);
                }

                // Vector for 0 along the n direction
                for (int j = 0; j < nTracePts; ++j)
                {
                    ZeroNvec[0][j] = 1.0 * m_traceNormals[1][j];
                    ZeroNvec[1][j] = -1.0 * m_traceNormals[0][j];
                }

                Array<OneD, NekDouble> tmp0(nBndEdgePts, 0.0);
                Array<OneD, NekDouble> tmp1(nBndEdgePts, 0.0);
                for (int k = 0; k < m_spacedim; ++k)
                {
                    Vmath::Vvtvp(nBndEdgePts, &ZeroNvec[k][id2], 1,
                                 &MFtraceBwd[0][k][id2], 1, &tmp0[0], 1,
                                 &tmp0[0], 1);

                    Vmath::Vvtvp(nBndEdgePts, &ZeroNvec[k][id2], 1,
                                 &MFtraceBwd[1][k][id2], 1, &tmp1[0], 1,
                                 &tmp1[0], 1);

                    Vmath::Vcopy(nBndEdgePts, &ZeroNvec[k][id2], 1,
                                 &qBwdVec[k][id2], 1);
                }

                Vmath::Vcopy(nBndEdgePts, &tmp0[0], 1, &qBwd[0][id2], 1);
                Vmath::Vcopy(nBndEdgePts, &tmp1[0], 1, &qBwd[1][id2], 1);
            }
        }
    }
}

Array<OneD, NekDouble> MMFSystem::ComputeLaplacianDiff(
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

void MMFSystem::ComputeGradient(
    const Array<OneD, const Array<OneD, NekDouble>> &MF1st,
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &velocity,
    Array<OneD, NekDouble> &velmag)

{
    boost::ignore_unused(inarray);

    int nq = m_fields[0]->GetTotPoints();

    Array<OneD, Array<OneD, NekDouble>> qfield(m_expdim);

    // if (GradCompt == SolverUtils::eWeakGrad)
    // {
    //     ComputeGradientWeak(MF1st, inarray, qfield);
    // }

    // else if (GradCompt == SolverUtils::eDirectGrad)
    // {
    //     ComputeGradientDirect(MF1st, inarray, qfield);
    // }

    for (int i = 0; i < nq; i++)
    {
        velocity[0][i] =
            qfield[0][i] * MF1st[0][i] + qfield[1][i] * MF1st[1][i];
        velocity[1][i] =
            qfield[0][i] * MF1st[0][i + nq] + qfield[1][i] * MF1st[1][i + nq];
        velocity[2][i] = qfield[0][i] * MF1st[0][i + 2 * nq] +
                         qfield[1][i] * MF1st[1][i + 2 * nq];
    }

    // Compute velocity magnitude
    velmag = Array<OneD, NekDouble>(nq, 0.0);
    for (int k = 0; k < m_spacedim; k++)
    {
        Vmath::Vvtvp(nq, velocity[k], 1, velocity[k], 1, velmag, 1, velmag, 1);
    }

    Vmath::Vsqrt(nq, velmag, 1, velmag, 1);
}

void MMFSystem::ComputeGradientWeak(
    const Array<OneD, const Array<OneD, NekDouble>> &MF1st,
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &qfield)

{
    int nq        = m_fields[0]->GetTotPoints();
    int ncoeffs   = m_fields[0]->GetNcoeffs();
    int nTracePts = m_fields[0]->GetTrace()->GetTotPoints();

    Array<OneD, NekDouble> qcoeffs(ncoeffs);

    Array<OneD, NekDouble> ufluxFwd(nTracePts);
    Array<OneD, NekDouble> ufluxBwd(nTracePts);

    for (int j = 0; j < m_expdim; ++j)
    {
        qfield[j] = Array<OneD, NekDouble>(nq, 0.0);
    }

    Array<OneD, Array<OneD, NekDouble>> DivMF;
    ComputeDivMF(eEuclidean, MF1st, DivMF);

    // Compute q_{\eta} and q_{\xi} to bbtain numerical fluxes
    // Get flux = u* ( e^i \cdot \vec{n} )
    // GetufluxMMF(0, inarray, ufluxFwd, ufluxBwd);
    Array<OneD, NekDouble> uflux(nTracePts);
    uflux = ComputeufluxMMF(0, inarray);

    Array<OneD, NekDouble> tmp(nq);
    Array<OneD, NekDouble> tmpc(ncoeffs);
    for (int j = 0; j < m_expdim; ++j)
    {
        // Compute L2 = \int ( \partial phi / \partial x_j) u d x
        m_fields[0]->IProductWRTDirectionalDerivBase(MF1st[j], inarray,
                                                     qcoeffs);

        // Compute -L2
        Vmath::Neg(ncoeffs, qcoeffs, 1);

        // Compute L = -L2 + \int_{\partial} u* n_x dx
        // m_fields[0]->AddFwdBwdTraceIntegral(ufluxFwd[j], ufluxBwd[j],
        // qcoeffs);
        Vmath::Vmul(nTracePts, &uflux[0], 1, &m_ncdotMFFwd[j][0], 1,
                    &ufluxFwd[0], 1);
        Vmath::Vmul(nTracePts, &uflux[0], 1, &m_ncdotMFBwd[j][0], 1,
                    &ufluxBwd[0], 1);

        m_fields[0]->AddFwdBwdTraceIntegral(ufluxFwd, ufluxBwd, qcoeffs);

        // Add -\int ( \nabla \cdot e^i ) u dx
        Vmath::Vmul(nq, &inarray[0], 1, &DivMF[j][0], 1, &tmp[0], 1);
        Vmath::Neg(nq, &tmp[0], 1);

        m_fields[0]->IProductWRTBase(tmp, tmpc);
        Vmath::Vadd(ncoeffs, &tmpc[0], 1, &qcoeffs[0], 1, &qcoeffs[0], 1);
        // m_fields[0]->SetPhysState(false);

        // Compute M^{-1} ( -L2 + \int_{\partial} u* n_x dx )
        m_fields[0]->MultiplyByElmtInvMass(qcoeffs, qcoeffs);
        m_fields[0]->BwdTrans(qcoeffs, qfield[j]);
    }
}

void MMFSystem::ComputeGradientDirect(
    const Array<OneD, const Array<OneD, NekDouble>> &MF1st,
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &qfield)

{
    int nq = GetTotPoints();
    for (int j = 0; j < m_expdim; ++j)
    {
        qfield[j] = Array<OneD, NekDouble>(nq, 0.0);
        MMFDirectionalDeriv(MF1st[j], inarray, qfield[j]);
    }
}

void MMFSystem::ComputeGradientSphericalCoord(
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray)
{
    int nq = GetNpoints();

    outarray = Array<OneD, NekDouble>(m_spacedim * nq);

    Array<OneD, NekDouble> DPhi(nq);
    Array<OneD, NekDouble> DTh(nq);

    Array<OneD, NekDouble> Phi(m_spacedim * nq);
    Array<OneD, NekDouble> Theta(m_spacedim * nq);

    ComputeSphericalTangentVector(Phi, Theta);

    MMFDirectionalDeriv(Phi, inarray, DPhi);
    MMFDirectionalDeriv(Theta, inarray, DTh);

    Array<OneD, NekDouble> x0(nq);
    Array<OneD, NekDouble> x1(nq);
    Array<OneD, NekDouble> x2(nq);

    m_fields[0]->GetCoords(x0, x1, x2);

    NekDouble sin_varphi, cos_varphi, sin_theta, cos_theta;
    NekDouble x0j, x1j, x2j;
    for (int i = 0; i < nq; i++)
    {
        x0j = x0[i];
        x1j = x1[i];
        x2j = x2[i];

        CartesianToNewSpherical(x0j, x1j, x2j, sin_varphi, cos_varphi,
                                sin_theta, cos_theta);

        // DirectGradPhi[i] = -DirectGradPhi[i] / cos_theta;
        // DirectGradTh[i]  = DirectGradTh[i] / cos_theta;
        if (m_MMFActivation[i])
        {
            outarray[i] =
                DPhi[i] / sin_theta / sin_theta * Phi[i] + DTh[i] * Theta[i];

            outarray[i + nq] = DPhi[i] / sin_theta / sin_theta * Phi[i + nq] +
                               DTh[i] * Theta[i + nq];

            outarray[i + 2 * nq] =
                DPhi[i] / sin_theta / sin_theta * Phi[i + 2 * nq] +
                DTh[i] * Theta[i + 2 * nq];
        }
    }
}

void MMFSystem::ComputeAxisAlignedLOCALMovingframes(
    const Array<OneD, const Array<OneD, NekDouble>> &AxisMF,
    Array<OneD, Array<OneD, NekDouble>> &movingframes)
{
    int nq = m_fields[0]->GetNpoints();

    Array<OneD, NekDouble> tmp(nq, 0.0);
    for (int k = 0; k < m_spacedim; ++k)
    {
        Vmath::Vvtvp(nq, &AxisMF[2][k * nq], 1, &AxisMF[2][k * nq], 1, &tmp[0],
                     1, &tmp[0], 1);
    }
    Vmath::Vsqrt(nq, tmp, 1, tmp, 1);

    // e1_sp = e1 - ( \vec{k} \cdot e1) / || \vec{k} ||^2 \vec{k}
    Array<OneD, NekDouble> kedot(nq, 0.0);
    Array<OneD, NekDouble> kmag2(nq, 0.0);
    for (int k = 0; k < m_spacedim; ++k)
    {
        Vmath::Vdiv(nq, &AxisMF[2][k * nq], 1, &tmp[0], 1,
                    &movingframes[2][k * nq], 1);

        Vmath::Vvtvp(nq, &movingframes[2][k * nq], 1, &movingframes[0][k * nq],
                     1, &kedot[0], 1, &kedot[0], 1);
        Vmath::Vvtvp(nq, &movingframes[2][k * nq], 1, &movingframes[2][k * nq],
                     1, &kmag2[0], 1, &kmag2[0], 1);
    }
    Vmath::Vdiv(nq, kedot, 1, kmag2, 1, kedot, 1);
    Vmath::Neg(nq, kedot, 1);

    for (int k = 0; k < m_spacedim; ++k)
    {
        Vmath::Vvtvp(nq, &kedot[0], 1, &movingframes[2][k * nq], 1,
                     &movingframes[0][k * nq], 1, &movingframes[0][k * nq], 1);
    }

    tmp = Array<OneD, NekDouble>(nq, 0.0);
    for (int k = 0; k < m_spacedim; ++k)
    {
        Vmath::Vvtvp(nq, &movingframes[0][k * nq], 1, &movingframes[0][k * nq],
                     1, &tmp[0], 1, &tmp[0], 1);
    }
    Vmath::Vsqrt(nq, tmp, 1, tmp, 1);

    for (int k = 0; k < m_spacedim; ++k)
    {
        Vmath::Vdiv(nq, &movingframes[0][k * nq], 1, &tmp[0], 1,
                    &movingframes[0][k * nq], 1);
    }

    // Compute e3 \times e1 = e2
    movingframes[1] = VectorCrossProdMF(movingframes[2], movingframes[0]);
}

void MMFSystem::CheckMeshErr()
{
    int nq = m_fields[0]->GetTotPoints();

    Array<OneD, NekDouble> x(nq);
    Array<OneD, NekDouble> y(nq);
    Array<OneD, NekDouble> z(nq);

    Array<OneD, NekDouble> rad(nq, 0.0);

    m_fields[0]->GetCoords(x, y, z);

    // Mesh error
    Array<OneD, NekDouble> MeshErr(nq);

    // Surface normal vector error
    Array<OneD, NekDouble> SurfaceNormalError(nq, 1.0);

    NekDouble MErrint = 0.0;
    // x^2 + y^2 + z^2 = R^2
    if (m_surfaceType == eSphere)
    {
        Vmath::Vvtvp(nq, x, 1, x, 1, rad, 1, rad, 1);
        Vmath::Vvtvp(nq, y, 1, y, 1, rad, 1, rad, 1);
        Vmath::Vvtvp(nq, z, 1, z, 1, rad, 1, rad, 1);
        Vmath::Vsqrt(nq, rad, 1, rad, 1);

        Vmath::Sadd(nq, -1.0 * m_SphereExactRadius, rad, 1, MeshErr, 1);
        Vmath::Vabs(nq, MeshErr, 1, MeshErr, 1);

        NekDouble exactint =
            4 * m_pi * m_SphereExactRadius * m_SphereExactRadius;
        Array<OneD, NekDouble> ones(nq, 1.0);
        MErrint = m_fields[0]->PhysIntegral(ones);
        MErrint = fabs(MErrint - exactint);

        // Check SurfaceNormalVector = 1. 0 - e^3 \cdot rad / ||e^3|| / || rad
        // ||
        Array<OneD, NekDouble> tmp(nq, 0.0);
        Vmath::Vvtvp(nq, &x[0], 1, &m_movingframes[2][0], 1, &tmp[0], 1,
                     &tmp[0], 1);
        Vmath::Vvtvp(nq, &y[0], 1, &m_movingframes[2][nq], 1, &tmp[0], 1,
                     &tmp[0], 1);
        Vmath::Vvtvp(nq, &z[0], 1, &m_movingframes[2][2 * nq], 1, &tmp[0], 1,
                     &tmp[0], 1);

        Vmath::Vdiv(nq, tmp, 1, rad, 1, tmp, 1);
        Vmath::Vsub(nq, tmp, 1, SurfaceNormalError, 1, SurfaceNormalError, 1);
        Vmath::Vabs(nq, SurfaceNormalError, 1, SurfaceNormalError, 1);
    }

    // x^2/a^2 + y^2/b^2 + z^2/c^2 = 1
    else if (m_surfaceType == eEllipsoid)
    {
        Array<OneD, NekDouble> tmp(nq);

        Vmath::Vmul(nq, x, 1, x, 1, tmp, 1);
        Vmath::Smul(nq, 1.0 / (m_Radx * m_Radx), tmp, 1, tmp, 1);
        Vmath::Vadd(nq, tmp, 1, rad, 1, rad, 1);

        Vmath::Vmul(nq, y, 1, y, 1, tmp, 1);
        Vmath::Smul(nq, 1.0 / (m_Rady * m_Rady), tmp, 1, tmp, 1);
        Vmath::Vadd(nq, tmp, 1, rad, 1, rad, 1);

        Vmath::Vmul(nq, z, 1, z, 1, tmp, 1);
        Vmath::Smul(nq, 1.0 / (m_Radz * m_Radz), tmp, 1, tmp, 1);
        Vmath::Vadd(nq, tmp, 1, rad, 1, rad, 1);

        Vmath::Sadd(nq, -1.0, rad, 1, MeshErr, 1);
        Vmath::Vabs(nq, MeshErr, 1, MeshErr, 1);
    }

    std::cout << "Mesh Error: "
                 "***********************************************************"
              << std::endl;
    std::cout << "Mesh error (Linf) = " << Vmath::Vamax(nq, MeshErr, 1)
              << std::endl;
    std::cout << "Mesh error (L2) = " << RootMeanSquare(MeshErr) << std::endl;
    std::cout << "Jac error = " << MErrint << std::endl;
    std::cout << "Surface Normal error = " << RootMeanSquare(SurfaceNormalError)
              << ", max = " << Vmath::Vamax(nq, SurfaceNormalError, 1)
              << std::endl;

    // Plot Mesh error
    // int nvar    = 2;
    // int ncoeffs = m_fields[0]->GetNcoeffs();

    // std::string outname = m_sessionName + "_MeshErr.chk";

    // std::vector<Array<OneD, NekDouble>> fieldcoeffs(nvar);
    // for (int i = 0; i < nvar; ++i)
    // {
    //     fieldcoeffs[i] = Array<OneD, NekDouble>(ncoeffs, 0.0);
    // }

    // std::vector<std::string> variables(nvar);
    // variables[0] = "Mesh Error";
    // variables[1] = "SurfaceNormal Error";

    // m_fields[0]->FwdTrans(MeshErr, fieldcoeffs[0]);
    // m_fields[0]->FwdTrans(SurfaceNormalError, fieldcoeffs[1]);

    // WriteFld(outname, m_fields[0], fieldcoeffs, variables);
}

void MMFSystem::CheckMeshErr(
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    const Array<OneD, const Array<OneD, NekDouble>> &velocity)
{
    int nq = m_fields[0]->GetTotPoints();

    Array<OneD, NekDouble> MeshErr(nq);

    Array<OneD, NekDouble> velvector(m_spacedim * nq);
    for (int i = 0; i < m_spacedim; ++i)
    {
        Vmath::Vcopy(nq, &velocity[i][0], 1, &velvector[i * nq], 1);
    }

    // Test dfdu error
    Array<OneD, NekDouble> dfdu3err(nq);
    for (int i = 0; i < m_spacedim; ++i)
    {
        MMFDirectionalDeriv(movingframes[2], velocity[i],
                                          dfdu3err);

        std::cout << "i = " << i << ", dfdu3err = "
                  << RootMeanSquare(dfdu3err, m_MMFActivation) << std::endl;
    }

    Array<OneD, NekDouble> k_F_nabla_k(nq);
    Array<OneD, NekDouble> k_k_nabla_F(nq);

    k_k_nabla_F =
        ComputeVecCdotNabla(movingframes[2], movingframes[2], velvector);
    k_F_nabla_k =
        ComputeVecCdotNabla(movingframes[2], velvector, movingframes[2]);

    std::cout << "G1: k cdot (k cdot nabla) F = "
              << RootMeanSquare(k_k_nabla_F, m_MMFActivation) << std::endl;
    std::cout << "G2: k cdot (F cdot nabla) k = "
              << RootMeanSquare(k_F_nabla_k, m_MMFActivation) << std::endl
              << std::endl;

    // Projection Error
    Array<OneD, Array<OneD, NekDouble>> velocityProj;
    Array<OneD, Array<OneD, NekDouble>> velMF;

    CartesianToMovingframes(velocity, m_movingframes, velMF);
    MovingframestoCartesian(velMF, m_movingframes, velocityProj);

    Array<OneD, NekDouble> tmp(nq);
    Array<OneD, NekDouble> Projerr(nq, 0.0);
    for (int i = 0; i < m_spacedim; ++i)
    {
        Vmath::Vsub(nq, &velocity[i][0], 1, &velocityProj[i][0], 1, &tmp[0], 1);
        Vmath::Vmul(nq, tmp, 1, tmp, 1, tmp, 1);
        Vmath::Vadd(nq, tmp, 1, Projerr, 1, Projerr, 1);
    }
    Vmath::Vsqrt(nq, Projerr, 1, Projerr, 1);
}

Array<OneD, NekDouble> MMFSystem::WeakDGMMFLDG(
    const int var, const Array<OneD, const NekDouble> &inarray,
    const NekDouble time)
{
    int nq = m_fields[0]->GetTotPoints();

    Array<OneD, NekDouble> outarray(nq);
    switch (m_expdim)
    {
        case 1:
        {
            outarray = WeakDGMMFLDG1D(var, inarray, time);
        }
        break;

        case 2:
        {
            outarray = WeakDGMMFLDG2D(var, inarray, time);
        }
        break;

        case 3:
        {
            // outarray = WeakDGMMFLDG3D(var,inarray);
        }
        break;

        default:
            break;
    }
    return outarray;
}

Array<OneD, NekDouble> MMFSystem::WeakDGMMFLDG1D(
    const int var, const Array<OneD, const NekDouble> &inarray,
    const NekDouble time)
{
    boost::ignore_unused(time);

    int nq        = m_fields[var]->GetTotPoints();
    int ncoeffs   = m_fields[var]->GetNcoeffs();
    int nTracePts = m_fields[var]->GetTrace()->GetTotPoints();

    Array<OneD, NekDouble> tmp(nq);
    Array<OneD, NekDouble> tmpc(ncoeffs);
    Array<OneD, NekDouble> qfieldc(ncoeffs);

    Array<OneD, NekDouble> flux(nTracePts);
    Array<OneD, NekDouble> fluxFwd(nTracePts);
    Array<OneD, NekDouble> fluxBwd(nTracePts);

    Array<OneD, Array<OneD, NekDouble>> qfieldMMF(m_shapedim);
    for (int j = 0; j < m_expdim; ++j)
    {
        qfieldMMF[j] = Array<OneD, NekDouble>(nq, 0.0);
    }

    // Compute the numerical flux for the scalar $u$ variable
    flux = ComputeufluxMMF(var, inarray);

    // Compute \vec{q} = \nabla u
    for (int j = 0; j < m_expdim; ++j)
    {
        // Compute L2 = \int ( \partial phi / \partial x_j) u d x
        // m_fields[0]->IProductWRTDirectionalDerivBase(m_movingframes[j], inarray,
        //                                              qfieldc);
        
        m_fields[0]->IProductWRTDerivBase(MultiRegions::eS, inarray, qfieldc);

        // Compute -L2
        Vmath::Neg(ncoeffs, qfieldc, 1);

        // Compute L = -L2 + \int_{\partial} u* n_x dx
        // Vmath::Vmul(nTracePts, &flux[0], 1, &m_ncdotMFFwd[j][0], 1, &fluxFwd[0], 1);
        // Vmath::Vmul(nTracePts, &flux[0], 1, &m_ncdotMFBwd[j][0], 1, &fluxBwd[0], 1);

        // m_fields[var]->AddFwdBwdTraceIntegral(fluxFwd, fluxBwd, qfieldc);
        m_fields[var]->AddTraceIntegral(flux, qfieldc);

        // Add -\int ( \nabla \cdot e^i ) u dx
        // Vmath::Vmul(nq, &inarray[0], 1, &m_DivMF[j][0], 1, &tmp[0], 1);
        // Vmath::Neg(nq, &tmp[0], 1);

        // m_fields[0]->IProductWRTBase(tmp, tmpc);
        // Vmath::Vadd(ncoeffs, tmpc, 1, qfieldc, 1, qfieldc, 1);
        m_fields[0]->SetPhysState(false);

        // Compute M^{-1} ( -L2 + \int_{\partial} u* n_x dx )
        m_fields[0]->MultiplyByElmtInvMass(qfieldc, qfieldc);
        m_fields[0]->BwdTrans(qfieldc, qfieldMMF[j]);
    }

    // Compute u from q_{\eta} and q_{\xi} to obtain numerical fluxes
    fluxFwd =
        ComputeqfluxMMF(var, inarray, m_ncdotMFFwd, m_ncdotMFBwd, qfieldMMF);

    // m_fields[var]->IProductWRTDerivBase(qfiel
    qfieldc = Array<OneD, NekDouble>(ncoeffs, 0.0);
    // for (int j = 0; j < m_expdim; ++j)
    // {
        // m_fields[0]->IProductWRTDirectionalDerivBase(m_movingframes[j],
        //                                              qfieldMMF[j], tmpc);
    m_fields[0]->IProductWRTDerivBase(MultiRegions::eS, qfieldMMF[0], qfieldc);

    //    Vmath::Vadd(ncoeffs, tmpc, 1, qfieldc, 1, qfieldc, 1);
  //  }
    Vmath::Neg(ncoeffs, qfieldc, 1);

    // Evaulate  <\phi, \hat{F}\cdot n> - outarray[i]
    m_fields[var]->AddTraceIntegral(fluxFwd, qfieldc);
    m_fields[var]->SetPhysState(false);

    return qfieldc;
}



Array<OneD, NekDouble> MMFSystem::WeakDGMMFLDG2D(
    const int var, const Array<OneD, const NekDouble> &inarray,
    const NekDouble time)
{
    boost::ignore_unused(time);

    int nq        = m_fields[var]->GetTotPoints();
    int ncoeffs   = m_fields[var]->GetNcoeffs();
    int nTracePts = m_fields[var]->GetTrace()->GetTotPoints();

    Array<OneD, NekDouble> tmp(nq);
    Array<OneD, NekDouble> tmpc(ncoeffs);
    Array<OneD, NekDouble> qfieldc(ncoeffs);

    Array<OneD, NekDouble> flux(nTracePts);
    Array<OneD, NekDouble> fluxFwd(nTracePts);
    Array<OneD, NekDouble> fluxBwd(nTracePts);

    Array<OneD, Array<OneD, NekDouble>> qfieldMMF(m_shapedim);
    for (int j = 0; j < m_expdim; ++j)
    {
        qfieldMMF[j] = Array<OneD, NekDouble>(nq, 0.0);
    }

    // Compute the numerical flux for the scalar $u$ variable
    flux = ComputeufluxMMF(var, inarray);

    // Compute \vec{q} = \nabla u
    for (int j = 0; j < m_expdim; ++j)
    {
        // Compute L2 = \int ( \partial phi / \partial x_j) u d x
        m_fields[0]->IProductWRTDirectionalDerivBase(m_movingframes[j], inarray,
                                                     qfieldc);

        // Compute -L2
        Vmath::Neg(ncoeffs, qfieldc, 1);

        // Compute L = -L2 + \int_{\partial} u* n_x dx
        Vmath::Vmul(nTracePts, &flux[0], 1, &m_ncdotMFFwd[j][0], 1, &fluxFwd[0], 1);
        Vmath::Vmul(nTracePts, &flux[0], 1, &m_ncdotMFBwd[j][0], 1, &fluxBwd[0], 1);

        m_fields[var]->AddFwdBwdTraceIntegral(fluxFwd, fluxBwd, qfieldc);

        // Add -\int ( \nabla \cdot e^i ) u dx
        Vmath::Vmul(nq, &inarray[0], 1, &m_DivMF[j][0], 1, &tmp[0], 1);
        Vmath::Neg(nq, &tmp[0], 1);

        m_fields[0]->IProductWRTBase(tmp, tmpc);
        Vmath::Vadd(ncoeffs, tmpc, 1, qfieldc, 1, qfieldc, 1);
        m_fields[0]->SetPhysState(false);

        // Compute M^{-1} ( -L2 + \int_{\partial} u* n_x dx )
        m_fields[0]->MultiplyByElmtInvMass(qfieldc, qfieldc);
        m_fields[0]->BwdTrans(qfieldc, qfieldMMF[j]);
    }

    // Compute u from q_{\eta} and q_{\xi} to obtain numerical fluxes
    fluxFwd =
        ComputeqfluxMMF(var, inarray, m_ncdotMFFwd, m_ncdotMFBwd, qfieldMMF);

    // m_fields[var]->IProductWRTDerivBase(qfiel
    qfieldc = Array<OneD, NekDouble>(ncoeffs, 0.0);
    for (int j = 0; j < m_expdim; ++j)
    {
        m_fields[0]->IProductWRTDirectionalDerivBase(m_movingframes[j],
                                                     qfieldMMF[j], tmpc);
        Vmath::Vadd(ncoeffs, tmpc, 1, qfieldc, 1, qfieldc, 1);
    }
    Vmath::Neg(ncoeffs, qfieldc, 1);

    // Evaulate  <\phi, \hat{F}\cdot n> - outarray[i]
    m_fields[var]->AddTraceIntegral(fluxFwd, qfieldc);
    m_fields[var]->SetPhysState(false);

    return qfieldc;
}


Array<OneD, NekDouble> MMFSystem::WeakDGMMFLDG2D(
    const int var, const Array<OneD, const NekDouble> &inarray,
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    const NekDouble time)
{
    boost::ignore_unused(time);

    int nq        = m_fields[var]->GetTotPoints();
    int ncoeffs   = m_fields[var]->GetNcoeffs();
    int nTracePts = m_fields[var]->GetTrace()->GetTotPoints();

    Array<OneD, NekDouble> tmp(nq);
    Array<OneD, NekDouble> tmpc(ncoeffs);
    Array<OneD, NekDouble> qfieldc(ncoeffs);

    Array<OneD, NekDouble> flux(nTracePts);
    Array<OneD, NekDouble> fluxFwd(nTracePts);
    Array<OneD, NekDouble> fluxBwd(nTracePts);

    Array<OneD, Array<OneD, NekDouble>> qfieldMMF(m_shapedim);
    for (int j = 0; j < m_expdim; ++j)
    {
        qfieldMMF[j] = Array<OneD, NekDouble>(nq, 0.0);
    }

    // Compute the numerical flux for the scalar $u$ variable
    flux = ComputeufluxMMF(var, inarray);

    Array<OneD, Array<OneD, NekDouble>> ncdotMFFwd;
    Array<OneD, Array<OneD, NekDouble>> ncdotMFBwd;
    ComputencdotMF(movingframes, ncdotMFFwd, ncdotMFBwd, 0);

    Array<OneD, Array<OneD, NekDouble>> DivMF;
    ComputeDivMF(m_DerivType, movingframes, DivMF, 0);

    // Compute \vec{q} = \nabla u
    for (int j = 0; j < m_expdim; ++j)
    {
        // Compute L2 = \int ( \partial phi / \partial x_j) u d x
        m_fields[0]->IProductWRTDirectionalDerivBase(movingframes[j], inarray, qfieldc);

        // Compute -L2
        Vmath::Neg(ncoeffs, qfieldc, 1);

        // Compute L = -L2 + \int_{\partial} u* n_x dx
        Vmath::Vmul(nTracePts, &flux[0], 1, &ncdotMFFwd[j][0], 1, &fluxFwd[0], 1);
        Vmath::Vmul(nTracePts, &flux[0], 1, &ncdotMFBwd[j][0], 1, &fluxBwd[0], 1);

        m_fields[var]->AddFwdBwdTraceIntegral(fluxFwd, fluxBwd, qfieldc);

        // Add -\int ( \nabla \cdot e^i ) u dx
        Vmath::Vmul(nq, &inarray[0], 1, &DivMF[j][0], 1, &tmp[0], 1);
        Vmath::Neg(nq, &tmp[0], 1);

        m_fields[0]->IProductWRTBase(tmp, tmpc);
        Vmath::Vadd(ncoeffs, tmpc, 1, qfieldc, 1, qfieldc, 1);
        m_fields[0]->SetPhysState(false);

        // Compute M^{-1} ( -L2 + \int_{\partial} u* n_x dx )
        m_fields[0]->MultiplyByElmtInvMass(qfieldc, qfieldc);
        m_fields[0]->BwdTrans(qfieldc, qfieldMMF[j]);
    }

    // Compute u from q_{\eta} and q_{\xi} to obtain numerical fluxes
    fluxFwd = ComputeqfluxMMF(var, inarray, ncdotMFFwd, ncdotMFBwd, qfieldMMF);

    // m_fields[var]->IProductWRTDerivBase(qfiel
    qfieldc = Array<OneD, NekDouble>(ncoeffs, 0.0);
    for (int j = 0; j < m_expdim; ++j)
    {
        m_fields[0]->IProductWRTDirectionalDerivBase(movingframes[j],
                                                     qfieldMMF[j], tmpc);
        Vmath::Vadd(ncoeffs, tmpc, 1, qfieldc, 1, qfieldc, 1);
    }
    Vmath::Neg(ncoeffs, qfieldc, 1);

    // Evaulate  <\phi, \hat{F}\cdot n> - outarray[i]
    m_fields[var]->AddTraceIntegral(fluxFwd, qfieldc);
    m_fields[var]->SetPhysState(false);

    return qfieldc;
}

void MMFSystem::WeakDGMMFLaplacian(const int var,
                                   const Array<OneD, const NekDouble> &inarray,
                                   Array<OneD, NekDouble> &outarray)
{
    int ncoeffs = m_fields[var]->GetNcoeffs();

    Array<OneD, NekDouble> qfieldc(ncoeffs);
    qfieldc = WeakDGMMFLDG(var, inarray);

    m_fields[var]->BwdTrans(qfieldc, outarray);
}

void MMFSystem::WeakDGMMFDiffusion(const int var,
                                   const Array<OneD, const NekDouble> &inarray,
                                   Array<OneD, NekDouble> &outarray,
                                   const NekDouble time)
{
    int ncoeffs = m_fields[var]->GetNcoeffs();

    Array<OneD, NekDouble> qfieldc(ncoeffs);
    qfieldc = WeakDGMMFLDG(var, inarray, time);

    m_fields[var]->MultiplyByElmtInvMass(qfieldc, qfieldc);
    m_fields[var]->BwdTrans(qfieldc, outarray);
}

Array<OneD, NekDouble> MMFSystem::WeakDGMMFirstLDG(
    const int var, const Array<OneD, const Array<OneD, NekDouble>> movingframes,
    const Array<OneD, const NekDouble> &inarray)
{
    int nq        = m_fields[var]->GetTotPoints();
    int ncoeffs   = m_fields[var]->GetNcoeffs();
    int nTracePts = m_fields[var]->GetTrace()->GetTotPoints();

    Array<OneD, NekDouble> tmp(nq);
    Array<OneD, NekDouble> tmpc(ncoeffs);
    Array<OneD, NekDouble> qfieldc(ncoeffs);

    Array<OneD, NekDouble> flux(nTracePts);
    Array<OneD, NekDouble> fluxFwd(nTracePts);
    Array<OneD, NekDouble> fluxBwd(nTracePts);

    Array<OneD, Array<OneD, NekDouble>> DivMF;
    ComputeDivMF(m_DerivType, movingframes, DivMF);

    Array<OneD, Array<OneD, NekDouble>> ncdotMFFwd;
    Array<OneD, Array<OneD, NekDouble>> ncdotMFBwd;
    ComputencdotMF(movingframes, ncdotMFFwd, ncdotMFBwd);

    Array<OneD, Array<OneD, NekDouble>> qfieldMMF(m_shapedim);
    for (int j = 0; j < m_shapedim; ++j)
    {
        qfieldMMF[j] = Array<OneD, NekDouble>(nq, 0.0);
    }

    // Compute the numerical flux for the scalar $u$ variable
    flux = ComputeufluxMMF(var, inarray);

    // Compute \vec{q} = \nabla u
    for (int j = 0; j < m_shapedim; ++j)
    {
        // Compute L2 = \int ( \partial phi / \partial x_j) u d x
        m_fields[0]->IProductWRTDirectionalDerivBase(movingframes[j], inarray,
                                                     qfieldc);

        // Compute -L2
        Vmath::Neg(ncoeffs, qfieldc, 1);

        // Compute L = -L2 + \int_{\partial} u* n_x dx
        Vmath::Vmul(nTracePts, &flux[0], 1, &ncdotMFFwd[j][0], 1, &fluxFwd[0],
                    1);
        Vmath::Vmul(nTracePts, &flux[0], 1, &ncdotMFBwd[j][0], 1, &fluxBwd[0],
                    1);

        m_fields[var]->AddFwdBwdTraceIntegral(fluxFwd, fluxBwd, qfieldc);

        // Add -\int ( \nabla \cdot e^i ) u dx
        Vmath::Vmul(nq, &inarray[0], 1, &DivMF[j][0], 1, &tmp[0], 1);
        Vmath::Neg(nq, &tmp[0], 1);

        m_fields[0]->IProductWRTBase(tmp, tmpc);
        Vmath::Vadd(ncoeffs, tmpc, 1, qfieldc, 1, qfieldc, 1);
        m_fields[0]->SetPhysState(false);

        // Compute M^{-1} ( -L2 + \int_{\partial} u* n_x dx )
        m_fields[0]->MultiplyByElmtInvMass(qfieldc, qfieldc);
        m_fields[0]->BwdTrans(qfieldc, qfieldMMF[j]);
    }

    // Compute u from q_{\eta} and q_{\xi} to obtain numerical fluxes
    fluxFwd = ComputeqfluxMMF(var, inarray, ncdotMFFwd, ncdotMFBwd, qfieldMMF);

    // m_fields[var]->IProductWRTDerivBase(qfield, tmpc);
    qfieldc = Array<OneD, NekDouble>(ncoeffs, 0.0);
    for (int j = 0; j < m_shapedim; ++j)
    {
        m_fields[0]->IProductWRTDirectionalDerivBase(movingframes[j],
                                                     qfieldMMF[j], tmpc);
        Vmath::Vadd(ncoeffs, tmpc, 1, qfieldc, 1, qfieldc, 1);
    }
    Vmath::Neg(ncoeffs, qfieldc, 1);

    // Evaulate  <\phi, \hat{F}\cdot n> - outarray[i]
    m_fields[var]->AddTraceIntegral(fluxFwd, qfieldc);
    m_fields[var]->SetPhysState(false);

    return qfieldc;
}

void MMFSystem::WeakDGMMFirstLaplacian(
    const int var, const Array<OneD, const Array<OneD, NekDouble>> movingframes,
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray)
{
    int ncoeffs = m_fields[var]->GetNcoeffs();

    Array<OneD, NekDouble> qfieldc(ncoeffs);
    qfieldc = WeakDGMMFirstLDG(var, movingframes, inarray);

    m_fields[var]->BwdTrans(qfieldc, outarray);
}

Array<OneD, NekDouble> MMFSystem::ComputeufluxMMF(
    const int var, const Array<OneD, const NekDouble> &ufield)
{
    int nTracePts = m_fields[0]->GetTrace()->GetTotPoints();

    Array<OneD, NekDouble> outarray(nTracePts, 0.0);

    // Get the sign of (v \cdot n), v = an arbitrary vector
    // Evaluate upwind flux:
    // uflux = \hat{u} \phi \cdot u = u^{(+,-)} n
    Array<OneD, NekDouble> Fwd(nTracePts);
    Array<OneD, NekDouble> Bwd(nTracePts);
    m_fields[var]->GetFwdBwdTracePhys(ufield, Fwd, Bwd);

    // Upwind
    Vmath::Vcopy(nTracePts, Fwd, 1, outarray, 1);

    // Imposing weak boundary condition with flux
    if (m_fields[0]->GetBndCondExpansions().size())
    {
        ApplyScalarBCs(var, Fwd, outarray);
    }

    return outarray;
}

void MMFSystem::ApplyScalarBCs(const int var,
                               const Array<OneD, const NekDouble> &Fwd,
                               Array<OneD, NekDouble> &penaltyflux)
{
    // Number of boundary regions
    std::size_t nBndRegions = m_fields[var]->GetBndCondExpansions().size();
    std::size_t cnt         = 0;
    for (std::size_t i = 0; i < nBndRegions; ++i)
    {
        if (m_fields[var]->GetBndConditions()[i]->GetBoundaryConditionType() ==
            SpatialDomains::ePeriodic)
        {
            continue;
        }

        // Number of boundary expansion related to that region
        std::size_t nBndEdges =
            m_fields[var]->GetBndCondExpansions()[i]->GetExpSize();

        // Weakly impose boundary conditions by modifying flux values
        for (std::size_t e = 0; e < nBndEdges; ++e)
        {
            // std::size_t nBndEdgePts = m_fields[var]
            //                               ->GetBndCondExpansions()[i]
            //                               ->GetExp(e)
            //                               ->GetTotPoints();
            // std::size_t id1 =
            //     m_fields[var]->GetBndCondExpansions()[i]->GetPhys_Offset(e);
            // std::size_t id2 = m_fields[0]->GetTrace()->GetPhys_Offset(
            //     m_fields[var]->GetTraceMap()->GetBndCondTraceToGlobalTraceMap(
            //         cnt++));

            std::size_t nBndEdgePts = m_fields[var]
                                          ->GetBndCondExpansions()[i]
                                          ->GetExp(e)
                                          ->GetNumPoints(0);
            std::size_t id1 =
                m_fields[var]->GetBndCondExpansions()[i]->GetPhys_Offset(e);
            std::size_t id2 = m_fields[var]->GetTrace()->GetPhys_Offset(
                m_fields[0]->GetTraceMap()->GetBndCondIDToGlobalTraceID(cnt++));

            // AV boundary conditions
            if (boost::iequals(
                    m_fields[var]->GetBndConditions()[i]->GetUserDefined(),
                    "Wall") ||
                boost::iequals(
                    m_fields[var]->GetBndConditions()[i]->GetUserDefined(),
                    "Symmetry") ||
                boost::iequals(
                    m_fields[var]->GetBndConditions()[i]->GetUserDefined(),
                    "WallViscous") ||
                boost::iequals(
                    m_fields[var]->GetBndConditions()[i]->GetUserDefined(),
                    "WallAdiabatic"))
            {
                Vmath::Vcopy(nBndEdgePts, &Fwd[id2], 1, &penaltyflux[id2], 1);
            }

            // For Dirichlet boundary condition: uflux = g_D
            else if (m_fields[var]
                         ->GetBndConditions()[i]
                         ->GetBoundaryConditionType() ==
                     SpatialDomains::eDirichlet)
            {
                Vmath::Vcopy(
                    nBndEdgePts,
                    &(m_fields[var]->GetBndCondExpansions()[i]->GetPhys())[id1],
                    1, &penaltyflux[id2], 1);
            }
            // For Neumann boundary condition: uflux = u+
            else if ((m_fields[var]->GetBndConditions()[i])
                         ->GetBoundaryConditionType() ==
                     SpatialDomains::eNeumann)
            {
                Vmath::Vcopy(nBndEdgePts, &Fwd[id2], 1, &penaltyflux[id2], 1);
            }
        }
    }
}

/**
 * @brief Build the numerical flux for the 2nd order derivatives
 * todo: add variable coeff and h dependence to penalty term
 */
Array<OneD, NekDouble> MMFSystem::ComputeqfluxMMF(
    const int var, const Array<OneD, const NekDouble> &ufield,
    const Array<OneD, const Array<OneD, NekDouble>> &ncdotMFFwd,
    const Array<OneD, const Array<OneD, NekDouble>> &ncdotMFBwd,
    const Array<OneD, const Array<OneD, NekDouble>> &qfieldMMF)
{
    boost::ignore_unused(ncdotMFFwd, ncdotMFBwd);

    int nTracePts = m_fields[0]->GetTrace()->GetTotPoints();

    Array<OneD, NekDouble> uFwd(nTracePts);
    Array<OneD, NekDouble> uBwd(nTracePts);

    Array<OneD, NekDouble> qFwd(nTracePts);
    Array<OneD, NekDouble> qBwd(nTracePts);

    Array<OneD, NekDouble> qflux(nTracePts, 0.0);
    Array<OneD, NekDouble> outarray(nTracePts, 0.0);

    // Evaulate upwind flux:
    // qflux = \hat{q} \cdot u = q \cdot n - C_(11)*(u^+ - u^-)

    // Generate Stability term = - ( u- - u+ )
    Array<OneD, NekDouble> uterm(nTracePts);
    m_fields[var]->GetFwdBwdTracePhys(ufield, uFwd, uBwd);

    Vmath::Vsub(nTracePts, uFwd, 1, uBwd, 1, uterm, 1);
    Vmath::Smul(nTracePts, -m_LDGC11, uterm, 1, uterm, 1);

    for (int j = 0; j < m_expdim; ++j)
    {
        //  Compute Fwd and Bwd value of ufield of jth direction
        m_fields[var]->GetFwdBwdTracePhys(qfieldMMF[j], qFwd, qBwd);

        // Downwind
        Vmath::Vmul(nTracePts, ncdotMFBwd[j], 1, qBwd, 1, qflux, 1);

        // Flux = {Fwd, Bwd} * (nx, ny, nz) + uterm * (nx, ny)
        Vmath::Vadd(nTracePts, uterm, 1, qflux, 1, qflux, 1);

        // Imposing weak boundary condition with flux
        if (m_fields[0]->GetBndCondExpansions().size())
        {
            ApplyVectorBCs(var, j, qFwd, ncdotMFFwd, qflux);
        }

        // q_hat \cdot n = (q_xi \cdot n_xi) or (q_eta \cdot n_eta)
        // n_xi = n_x * tan_xi_x + n_y * tan_xi_y + n_z * tan_xi_z
        // n_xi = n_x * tan_eta_x + n_y * tan_eta_y + n_z*tan_eta_z
        Vmath::Vadd(nTracePts, qflux, 1, outarray, 1, outarray, 1);
    }

    return outarray;
}

/**
 * Diffusion: Imposing weak boundary condition for q with flux
 *  uflux = g_D  on Dirichlet boundary condition
 *  uflux = u_Fwd  on Neumann boundary condition
 */

void MMFSystem::ApplyVectorBCs(
    const std::size_t var, const std::size_t dir,
    const Array<OneD, const NekDouble> &qFwd,
    const Array<OneD, const Array<OneD, NekDouble>> &ncdotMFFwd,
    Array<OneD, NekDouble> &penaltyflux)
{
    std::size_t nBndRegions = m_fields[var]->GetBndCondExpansions().size();
    std::size_t cnt         = 0;

    for (std::size_t i = 0; i < nBndRegions; ++i)
    {
        if (m_fields[var]->GetBndConditions()[i]->GetBoundaryConditionType() ==
            SpatialDomains::ePeriodic)
        {
            continue;
        }

        std::size_t nBndEdges =
            m_fields[var]->GetBndCondExpansions()[i]->GetExpSize();

        // Weakly impose boundary conditions by modifying flux values
        for (std::size_t e = 0; e < nBndEdges; ++e)
        {
            // std::size_t nBndEdgePts = m_fields[var]
            //                               ->GetBndCondExpansions()[i]
            //                               ->GetExp(e)
            //                               ->GetTotPoints();

            // std::size_t id1 =
            //     m_fields[var]->GetBndCondExpansions()[i]->GetPhys_Offset(e);

            // std::size_t id2 = m_fields[0]->GetTrace()->GetPhys_Offset(
            //     m_fields[0]->GetTraceMap()->GetBndCondTraceToGlobalTraceMap(
            //         cnt++));

            std::size_t nBndEdgePts = m_fields[var]
                                          ->GetBndCondExpansions()[i]
                                          ->GetExp(e)
                                          ->GetNumPoints(0);
            std::size_t id1 =
                m_fields[0]->GetBndCondExpansions()[i]->GetPhys_Offset(e);
            std::size_t id2 = m_fields[0]->GetTrace()->GetPhys_Offset(
                m_fields[0]->GetTraceMap()->GetBndCondIDToGlobalTraceID(cnt++));

            // AV boundary conditions
            if (boost::iequals(
                    m_fields[var]->GetBndConditions()[i]->GetUserDefined(),
                    "Wall") ||
                boost::iequals(
                    m_fields[var]->GetBndConditions()[i]->GetUserDefined(),
                    "Symmetry") ||
                boost::iequals(
                    m_fields[var]->GetBndConditions()[i]->GetUserDefined(),
                    "WallViscous") ||
                boost::iequals(
                    m_fields[var]->GetBndConditions()[i]->GetUserDefined(),
                    "WallAdiabatic"))
            {
                Vmath::Zero(nBndEdgePts, &penaltyflux[id2], 1);
            }

            // For Dirichlet boundary condition:
            // qflux = q+ - C_11 (u+ -    g_D) (nx, ny)
            else if (m_fields[var]
                         ->GetBndConditions()[i]
                         ->GetBoundaryConditionType() ==
                     SpatialDomains::eDirichlet)
            {
                Array<OneD, NekDouble> tmp(nBndEdgePts);
                Vmath::Vmul(nBndEdgePts, &ncdotMFFwd[dir][id2], 1, &qFwd[id2],
                            1, &tmp[0], 1);

                Vmath::Vcopy(nBndEdgePts, &tmp[0], 1, &penaltyflux[id2], 1);
            }

            // For Neumann boundary condition: qflux = g_N
            else if ((m_fields[var]->GetBndConditions()[i])
                         ->GetBoundaryConditionType() ==
                     SpatialDomains::eNeumann)
            {
                Vmath::Vmul(
                    nBndEdgePts, &ncdotMFFwd[dir][id1], 1,
                    &(m_fields[var]->GetBndCondExpansions()[i]->GetPhys())[id1],
                    1, &penaltyflux[id2], 1);
            }
        }
    }
}

NekDouble MMFSystem::ComputeDistance(const Array<OneD, const NekDouble> &p1,
                                     const Array<OneD, const NekDouble> &p2)
{
    NekDouble distance;
    NekDouble dx, dy, dz;

    dx = p1[0] - p2[0];
    dy = p1[1] - p2[1];
    dz = p1[2] - p2[2];

    distance = sqrt(dx * dx + dy * dy + dz * dz);

    return distance;
}

void MMFSystem::Checkpoint_Output_1D(
    const int n, const Array<OneD, const NekDouble> &seglength,
    const Array<OneD, const NekDouble> &field,
    const Array<OneD, const NekDouble> &TimeMap)
{
    int nvar    = 3;
    int ncoeffs = m_fields[0]->GetNcoeffs();

    std::string outname1 =
        m_sessionName + "_" + boost::lexical_cast<std::string>(n) + "_1D.chk";
    std::vector<Array<OneD, NekDouble>> fieldcoeffs(nvar);
    for (int i = 0; i < nvar; ++i)
    {
        fieldcoeffs[i] = Array<OneD, NekDouble>(ncoeffs);
    }

    std::vector<std::string> variables(nvar);
    variables[0] = "seglength";
    variables[1] = "field";
    variables[2] = "TimeMap";

    // Normalized Time Vector
    m_fields[0]->FwdTransLocalElmt(seglength, fieldcoeffs[0]);
    m_fields[0]->FwdTransLocalElmt(field, fieldcoeffs[1]);
    m_fields[0]->FwdTransLocalElmt(TimeMap, fieldcoeffs[2]);

    WriteFld(outname1, m_fields[0], fieldcoeffs, variables);
}

void MMFSystem::Checkpoint_Output_Error(
    const int n, const Array<OneD, const NekDouble> &field,
    const Array<OneD, const NekDouble> &exactsoln)
{
    int nvar    = 3;
    int ncoeffs = m_fields[0]->GetNcoeffs();
    int nq = m_fields[0]->GetNpoints();

    std::string outname1 =
        m_sessionName + "_" + boost::lexical_cast<std::string>(n) + "_Error.chk";

    std::vector<Array<OneD, NekDouble>> fieldcoeffs(nvar);
    for (int i = 0; i < nvar; ++i)
    {
        fieldcoeffs[i] = Array<OneD, NekDouble>(ncoeffs);
    }

    std::vector<std::string> variables(nvar);
    variables[0] = "Field";
    variables[1] = "Exactsoln";
    variables[2] = "Error";

    // Normalized Time Vector
    m_fields[0]->FwdTransLocalElmt(field, fieldcoeffs[0]);
    m_fields[0]->FwdTransLocalElmt(exactsoln, fieldcoeffs[1]);

    Array<OneD, NekDouble> Error(nq);
    Vmath::Vsub(nq, field, 1, exactsoln, 1, Error, 1);
    m_fields[0]->FwdTransLocalElmt(Error, fieldcoeffs[2]);

    WriteFld(outname1, m_fields[0], fieldcoeffs, variables);
}

// Gaussian pulse of dudt to make a smoothed dudt in an element.
Array<OneD, NekDouble> MMFSystem::GaussianPulse(
    const Array<OneD, const NekDouble> &dudt,
    const Array<OneD, const NekDouble> &xp,
    const Array<OneD, const NekDouble> &yp,
    const Array<OneD, const NekDouble> &zp, const NekDouble Grad)
{
    int npts = xp.size();

    Array<OneD, NekDouble> outarray(npts, 0.0);

    NekDouble distx, disty, distz;
    NekDouble xi, yi, zi;
    for (int i = 0; i < npts; ++i)
    {
        xi = xp[i];
        yi = yp[i];
        zi = zp[i];
        for (int j = 0; j < npts; ++j)
        {
            distx = xp[j] - xi;
            disty = yp[j] - yi;
            distz = zp[j] - zi;

            if (dudt[i] > 0)
            {
                outarray[j] +=
                    dudt[i] *
                    exp(-1.0 * (distx * distx + disty * disty + distz * distz) /
                        Grad / Grad) /
                    npts;
            }
        }
    }

    return outarray;
}

Array<OneD, NekDouble> MMFSystem::SmoothCircle(
    const NekDouble &v_amp, const NekDouble &m_pis, const NekDouble &m_px,
    const NekDouble &m_py, const NekDouble &m_pz, const NekDouble &m_pr)
{

    int nq = m_fields[0]->GetTotPoints();

    Array<OneD, NekDouble> outarray(nq, 0.0);

    Array<OneD, NekDouble> x0(nq);
    Array<OneD, NekDouble> x1(nq);
    Array<OneD, NekDouble> x2(nq);

    m_fields[0]->GetCoords(x0, x1, x2);

    std::cout << "ScarStr= " << v_amp << ", ScarStrSlope = " << m_pis
              << ", ScarPis = "
              << ", x = " << m_px << " , y = " << m_py << ", z = " << m_pz
              << std::endl;

    NekDouble xlen, ylen, zlen, kernel;
    for (int j = 0; j < nq; ++j)
    {
        xlen = m_pis * x0[j] - m_px;
        ylen = m_pis * x1[j] - m_py;
        zlen = m_pis * x2[j] - m_pz;

        kernel      = 0.5 * (1.0 + tanh(((xlen + m_pr) * (xlen - m_pr) +
                                    (ylen + m_pr) * (ylen - m_pr) +
                                    (zlen + m_pr) * (zlen - m_pr)) /
                                            2.0 +
                                        0.5));
        outarray[j] = (1.0 - v_amp) * kernel + v_amp;
    }

    std::cout << "SmoothCircle = " << RootMeanSquare(outarray) << std::endl;

    return outarray;
}

Array<OneD, NekDouble> MMFSystem::SmoothLine(const NekDouble &v_amp,
                                             const NekDouble &m_pis,
                                             const NekDouble &m_px,
                                             const NekDouble &m_pr)
{

    int nq = m_fields[0]->GetTotPoints();

    Array<OneD, NekDouble> outarray(nq, 0.0);

    Array<OneD, NekDouble> x0(nq);
    Array<OneD, NekDouble> x1(nq);
    Array<OneD, NekDouble> x2(nq);

    m_fields[0]->GetCoords(x0, x1, x2);

    std::cout << "ScarSize = " << v_amp << ", ScarStr = " << v_amp
              << ", ScarPis = " << m_pis << ", x = " << m_px << std::endl;

    NekDouble xlen, kernel;
    for (int j = 0; j < nq; ++j)
    {
        xlen = m_pis * x2[j] - m_px;

        kernel      = 0.5 * (1.0 + tanh(((xlen + m_pr)) / 2.0 + 0.5));
        outarray[j] = (v_amp - 1.0) * kernel + 1.0;
    }

    std::cout << "SmoothLine = " << RootMeanSquare(outarray) << std::endl;

    return outarray;
}

void MMFSystem::ComputeTimeMap(const NekDouble time,
                               const NekDouble urest,
                               const Array<OneD, const NekDouble> &field,
                               const Array<OneD, const NekDouble> &dudt,
                               const Array<OneD, const int> &ValidTimeMap,
                               Array<OneD, NekDouble> &dudtHistory,
                               Array<OneD, NekDouble> &TimeMap)
{
    int nq = GetTotPoints();

    NekDouble fnewsum;
    // NekDouble uTol = 0.01;
    NekDouble dudtTol = 0.1;

    NekDouble udiff;
    for (int i = 0; i < nq; ++i)
    {
        udiff = field[i] - urest;
        // Only integrate of time if u > Tol, gradu > Tol, du/dt > 0
        if ((udiff > m_uTol) && (dudt[i] > dudtTol))
        {
            // Gradient as the main weight
            fnewsum = dudt[i] + dudtHistory[i];

            if(fabs(fnewsum)>dudtTol)
            {
                TimeMap[i] = (dudt[i] * time + dudtHistory[i] * TimeMap[i]) / fnewsum;
            }

            dudtHistory[i] += dudt[i];
        }
    }

    // NekDouble TimeMapMin = Vmath::Vmin(nq, TimeMap, 1);
    for (int i = 0; i < nq; ++i)
    {
        if (ValidTimeMap[i] == 0)
        {
            TimeMap[i] = 0.0;
        }
    }
}

void MMFSystem::PlotTimeMap(
    const Array<OneD, const int> &ValidTimeMap,
    const Array<OneD, const Array<OneD, NekDouble>> &AniStrength,
    const Array<OneD, const NekDouble> &TimeMap,
    const int nstep)
{
    boost::ignore_unused(ValidTimeMap);

    int nvar    = 7;
    int nq      = m_fields[0]->GetTotPoints();
    int ncoeffs = m_fields[0]->GetNcoeffs();

    std::string outname1 = m_sessionName + "_timemap_" +
                           boost::lexical_cast<std::string>(nstep) + ".chk";

    std::vector<Array<OneD, NekDouble>> fieldcoeffs(nvar);
    for (int i = 0; i < nvar; ++i)
    {
        fieldcoeffs[i] = Array<OneD, NekDouble>(ncoeffs);
    }

    std::vector<std::string> variables(nvar);
    variables[0] = "TimeMap";
    variables[1] = "ValidTimeMap";
    variables[2] = "AniStrength[0]";
    variables[3] = "AniStrength[1]";
    variables[4] = "TMgrad_x";
    variables[5] = "TMgrad_y";
    variables[6] = "TMgrad_z";

    // index:0 -> u
    std::cout << "Time Map: Max = " << Vmath::Vmax(nq, TimeMap, 1)
                << ", Min = " << Vmath::Vmin(nq, TimeMap, 1) << std::endl;

    m_fields[0]->FwdTransLocalElmt(TimeMap, fieldcoeffs[0]);

    Array<OneD, NekDouble> ValidTM(nq, 1.0);
    for (int i=0; i<nq; ++i)
    {
        ValidTM[i] = 1.0 * ValidTimeMap[i];
    }
    m_fields[0]->FwdTransLocalElmt(ValidTM, fieldcoeffs[1]);

    Array<OneD, Array<OneD, NekDouble>> TimeMapMF(m_spacedim);
    for (int k=0; k<m_spacedim; ++k)
    {
        TimeMapMF[k] = Array<OneD, NekDouble>(nq, 0.0);
    }

    // Compute the gradient of the time map
    m_fields[0]->FwdTransLocalElmt(AniStrength[0], fieldcoeffs[2]);
    m_fields[0]->FwdTransLocalElmt(AniStrength[1], fieldcoeffs[3]);

    // Compute the gradient of the time map
    // Array<OneD, Array<OneD, NekDouble>> Velocity(m_spacedim);
    // ComputeVelocityTimeMap(ValidTM, TimeMap, Velocity);
    Array<OneD, NekDouble> TmapGrad(m_spacedim * nq);
    TmapGrad = ComputeCovGrad(TimeMap, m_movingframes);

    Array<OneD, NekDouble> tmp(nq);
    for (int k=0; k<m_spacedim; ++k)
    {
        Vmath::Vcopy(nq, &TmapGrad[k * nq], 1, &tmp[0], 1);
        m_fields[0]->FwdTransLocalElmt(tmp, fieldcoeffs[k+4]);
    }

    WriteFld(outname1, m_fields[0], fieldcoeffs, variables);
}


void MMFSystem::v_GenerateSummary(SummaryList &s)
{
    int nq = m_fields[0]->GetNpoints();
    UnsteadySystem::v_GenerateSummary(s);

    AddSummaryItem(s, "Surface", SurfaceTypeMap[m_surfaceType]);
    if (m_expdim == 1)
    {
        AddSummaryItem(s, "Length", m_seglength[nq - 1]);
        AddSummaryItem(s, "Avg edge length", m_seglength[nq - 1] / (nq - 2));
    }
    AddSummaryItem(s, "DerivType", DerivTypeMap[m_DerivType]);

    // AddSummaryItem(s, "AdaptNewFramesTol", m_AdaptNewFramesTol);
    // AddSummaryItem(s, "VelActivateTol", m_VelActivationTol);
    // AddSummaryItem(s, "NoAlignInitRadius", m_NoAlignInitRadius);
    // AddSummaryItem(s, "GradLocType", GradLocTypeMap[m_GradLocType]);

    // AddSummaryItem(s, "GradComptType", GradComptTypeMap[m_GradComptType]);
    // AddSummaryItem(s, "GradIntType", GradIntTypeMap[m_GradIntType]);

    // AddSummaryItem(s, "Total grids", nq);

    // Max and Min EdgeLength
    Array<OneD, NekDouble> EdgeLength;
    Array<OneD, NekDouble> x0(nq);
    Array<OneD, NekDouble> x1(nq);
    Array<OneD, NekDouble> x2(nq);

    m_fields[0]->GetCoords(x0, x1, x2);
    // m_fields[0]->ComputeEdgeLength(x0, x1, x2, EdgeLength);

    // Actual lenght for circular arc
    if (m_surfaceType == eSphere)
    {
        NekDouble rad = sqrt(x0[0] * x0[0] + x1[0] * x1[0] + x2[0] * x2[0]);
        EdgeLength[0] =
            rad * acos(1.0 - 0.5 * EdgeLength[0] * EdgeLength[0] / rad / rad);
        EdgeLength[1] =
            rad * acos(1.0 - 0.5 * EdgeLength[1] * EdgeLength[1] / rad / rad);
    }

    // AddSummaryItem(s, "Max Edge Length", EdgeLength[0]);
    AddSummaryItem(s, "MMFdir", SpatialDomains::GeomMMFMap[m_MMFdir]);
}

} // namespace SolverUtils
} // namespace Nektar
