///////////////////////////////////////////////////////////////////////////////
//
// File MMFSystem.h
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
// Description: MMF system
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_MMFSYSTEM_H
#define NEKTAR_SOLVERUTILS_MMFSYSTEM_H

#include <SolverUtils/UnsteadySystem.h>
#include <boost/core/ignore_unused.hpp>

namespace Nektar
{
namespace SolverUtils
{

enum MMFOrder
{
    eMMFZero,
    eMMFFirst,
    eMMFSecond,
    eMMFRotationTest,
    SIZE_MMFOrderType
};

const char *const MMFOrderMap[] = {
    "MMFZero",
    "MMFFirst",
    "MMFSecond",
    "MMFRotationTest",
};

enum MMFLinearAlign
{
    eMMFDirect,
    eMMFAverage,
    eMMFLinear,
    eMMFCircular,
    SIZE_MMFLinearAlignType
};

const char *const MMFLinearAlignMap[] = {
    "MMFDirect",
    "MMFAverage",
    "MMFLinear",
    "MMFCircular",
};

enum PlotHHD
{
    eYes,
    eNo,
    SIZE_PlotHHDType
};

const char *const PlotHHDMap[] = {
    "Yes",
    "No",
};

enum SurfaceType
{
    eLine,
    ePlane,
    ePlaneLeft,
    ePlanePoint,
    ePlaneEmbed,
    ePolar,
    eSphere,
    eTRSphere,
    ePseudosphere,
    eEllipsoid,
    eTorus,
    eIrregular,
    eNonconvex,
    eAtrium,
    eCube,
    SIZE_SurfaceType
};

const char *const SurfaceTypeMap[] = {
    "Line",      "Plane",     "PlaneLeft",    "PlanePoint", "PlaneEmbed", "Polar",
    "Sphere",    "TRSphere",  "Pseudosphere", "Ellipsoid",  "Torus",
    "Irregular", "Nonconvex", "Atrium",       "Cube",
};

enum BoundaryCopyType
{
    eDirichlet,
    eNeumann,
    eFwdEQBwd,
    eFwdEQNegBwd,
    SIZE_BoundaryCopyType ///< Length of enum list
};

const char *const BoundaryCopyTypeMap[] = {
    "Dirichlet",
    "Neumann",
    "FwdEQBwd",
    "FwdEQNegBwd",
};

enum UpwindType
{
    eNotSet,        ///< flux not defined
    eAverage,       ///< averaged (or centred) flux
    eLaxFriedrich,  ///< Lax-Friedrich flux
    eUpwind,        ///  Upwind
    eRusanov,       ///< Rusanov flux
    eHLL,           ///< Harten-Lax-Leer flux
    eHLLC,          ///< Harten-Lax-Leer Contact wave flux
    SIZE_UpwindType ///< Length of enum list
};

const char *const UpwindTypeMap[] = {
    "NoSet", "Average", "LaxFriedrich", "Upwind", "Rusanov", "HLL", "HLLC",
};

enum DerivType
{
    eEuclidean,
    eCovariant,
    eExact,
    SIZE_DerivType
};

const char *const DerivTypeMap[] = {"Euclidean", "Covariant", "Exact"};

enum EvalType
{
    eStronger,
    eWeaker,
    SIZE_EvalType
};

const char *const EvalTypeMap[] = {
    "Stronger",
    "Weaker",
};

enum GradLocType
{
    eWaveback,
    eWavefront,
    eWholewave,
    SIZE_GradLocType
};

const char *const GradLocTypeMap[] = {
    "Waveback",
    "Wavefront",
    "Wholewave",
};

// enum GradIntType
// {
//     eWint,
//     eInt,
//     eStr,
//     SIZE_GradIntType
// };

// const char *const GradIntTypeMap[] = {
//     "Wint",
//     "Int",
//     "Str",
// };

// enum GradComptType
// {
//     eDirectGrad,
//     eWeakGrad,
//     SIZE_GradComptType
// };

// const char *const GradComptTypeMap[] = {
//     "DirectGrad",
//     "WeakGrad",
// };

enum FaceDirType
{
    eFwd,
    eBwd,
    SIZE_FaceDirType
};

const char *const FaceDirTypeMap[] = {
    "Fwd",
    "Bwd",
};

// enum CoordAxisDivScheme
// {
//     eGlobal,
//     eGlobalWithSDR,
//     eCADNoSDR,
//     eCADWithSDR,
//     SIZE_CoordAxisDivScheme
// };

// const char *const CoordAxisDivSchemeMap[] = {
//     "Global",
//     "GlobalWithSDR",
//     "CADNoSDR",
//     "CADWithSDR",
//  };

/// A base class for PDEs which include an advection component
class MMFSystem : virtual public UnsteadySystem
{
public:
    int m_shapedim;
    int m_mfdim;

    SurfaceType m_surfaceType;
    UpwindType m_upwindType;
    DerivType m_DerivType;

    Array<OneD, NekDouble> m_MMFfactors;

    NekDouble m_pi;

    NekDouble m_SphereExactRadius;
    NekDouble m_NoAlignInitRadius;
    NekDouble m_AdaptNewFramesTol;
    NekDouble m_VelActivationTol;
    NekDouble m_Radx, m_Rady, m_Radz;

    // Gradient computing options
    GradLocType m_GradLocType;
    // GradComptType m_GradComptType;
    // GradIntType m_GradIntType;
    DerivType m_DType;
    FaceDirType m_FaceDirType;
    // CoordAxisDivScheme m_CoordAxisDivScheme;

    // SOLVER_UTILS_EXPORT MMFSystem(const LibUtilities::SessionReaderSharedPtr
    // &pSession);
    SOLVER_UTILS_EXPORT MMFSystem(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const SpatialDomains::MeshGraphSharedPtr &pGraph);

    SOLVER_UTILS_EXPORT virtual ~MMFSystem();

    // SOLVER_UTILS_EXPORT virtual void v_InitObject(bool DeclareFields = true)
    // override;

    SOLVER_UTILS_EXPORT virtual void v_GenerateSummary(SummaryList &s) override;

    SOLVER_UTILS_EXPORT void MMFDirectionalDeriv(const Array<OneD, const NekDouble> &movingframe, 
                          const Array<OneD, const NekDouble> &inarray, 
                          Array<OneD, NekDouble> &outarray);

    SOLVER_UTILS_EXPORT void ComputeDirectionVector(
        Array<OneD, Array<OneD, NekDouble>> &distance);

    SOLVER_UTILS_EXPORT Array<OneD, NekDouble> ComputeDirectionVector(
        const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
        const Array<OneD, const NekDouble> &inarray,
        const Array<OneD, const NekDouble> &dudt);

    SOLVER_UTILS_EXPORT void MMFInitObject(
        const Array<OneD, const Array<OneD, NekDouble>> &AniStrength,
        const Array<OneD, const NekDouble> &AniDirection = NullNekDouble1DArray);

    // SOLVER_UTILS_EXPORT void MMFInitObject(
    // const Array<OneD, const Array<OneD, NekDouble>> &Anisotropy);

    SOLVER_UTILS_EXPORT void CopyBoundaryTrace(
        const Array<OneD, const NekDouble> &Fwd, Array<OneD, NekDouble> &Bwd,
        const BoundaryCopyType BDCopyType, const int var = 0,
        const std::string btype = "NoUserDefined");

    SOLVER_UTILS_EXPORT Array<OneD, NekDouble> ComputeVelocityMag(
        const Array<OneD, const NekDouble> &velocity);

    SOLVER_UTILS_EXPORT Array<OneD, NekDouble> ComputeVelocityMag(
        const Array<OneD, const Array<OneD, NekDouble>> &velocity);

    SOLVER_UTILS_EXPORT Array<OneD, int> Computedudt(
        const NekDouble uTol, const Array<OneD, const NekDouble> &fieldu,
        const Array<OneD, const NekDouble> &fielduold);

    SOLVER_UTILS_EXPORT void ComputedudtHistory(
        const Array<OneD, const int> &dudt,
        const Array<OneD, const NekDouble> &fieldu,
        Array<OneD, int> &dudtHistory, Array<OneD, int> &APindex);

    SOLVER_UTILS_EXPORT Array<OneD, int> ComputeZoneActivation(
        const NekDouble uTol, const Array<OneD, const NekDouble> &fieldu,
        const NekDouble NoAlignInitRadius);

    SOLVER_UTILS_EXPORT void GenerateHHDPlot(const int nstep);

    SOLVER_UTILS_EXPORT Array<OneD, NekDouble> ComputeVectorLaplacian2D(
        const Array<OneD, const int> &Activation,
        const Array<OneD, const NekDouble> &harmonic);

protected:
    // True or False to Rebuild MF with the Maurer-Cartan matrix
    int m_ImportedFiberExist = 0;    

    NekDouble m_Helmtau, m_HHDtau;
    NekDouble m_ZoneVarIndex;
    NekDouble m_uTol;
    NekDouble m_c121, m_c122, m_c123;
    NekDouble m_Initx, m_Inity, m_Initz;
    NekDouble m_ROIx, m_ROIy, m_ROIz;
    NekDouble m_alpha;
    NekDouble m_Incfreq;
    NekDouble m_SFinit;

    // Divergence Restore scheme variables
    int m_DivergenceRestore;
    Array<OneD, Array<OneD, NekDouble>> m_SphericalVector;
    Array<OneD, Array<OneD, NekDouble>> mf_LOCSPH;

    // Variable for LDG scheme
    NekDouble m_LDGC11;
    
    // Spherical coordinate vectors
    Array<OneD, int> m_MMFActivation;

    // Relative divergence related variables
    NekDouble m_LambDivSmoothL;

    Array<OneD, int> m_ValidTimeMap;

    // Gaussian Time Map variables
    int m_GaussianTimeMap;
    NekDouble m_GaussianRadius;

    Array<OneD, NekDouble> m_seglength;
    Array<OneD, Array<OneD, NekDouble>> m_sphereMF;
    Array<OneD, Array<OneD, NekDouble>> m_pseudosphereMF;

    Array<OneD, Array<OneD, NekDouble>> m_ellipticMF;
    Array<OneD, Array<OneD, NekDouble>> m_polarMF;

    // Moving frames
    Array<OneD, Array<OneD, NekDouble>> m_movingframes;
    Array<OneD, Array<OneD, NekDouble>> m_surfaceNormal;

    Array<OneD, Array<OneD, NekDouble>> m_ncdotMFFwd;
    Array<OneD, Array<OneD, NekDouble>> m_ncdotMFBwd;

    Array<OneD, Array<OneD, NekDouble>> m_nperpcdotMFFwd;
    Array<OneD, Array<OneD, NekDouble>> m_nperpcdotMFBwd;

    Array<OneD, NekDouble> m_ThetacdotMF;
    Array<OneD, NekDouble> m_PhicdotMF;

    Array<OneD, Array<OneD, NekDouble>> m_DivMF;
    Array<OneD, Array<OneD, NekDouble>> m_CurlMF;
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> m_CurlMFold;

    Array<OneD, Array<OneD, NekDouble>> m_LOCALDivMF;
    Array<OneD, Array<OneD, NekDouble>> m_LOCALCurlMF;

    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> m_MFConnection;
    Array<OneD, Array<OneD, NekDouble>> m_MFCurvature;

    // MFdim \times spacedim \times npts
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> m_MFtraceFwd;
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> m_MFtraceBwd;

    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> m_ntimesMFFwd;
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> m_ntimesMFBwd;
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> m_ntimes_ntimesMFFwd;
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> m_ntimes_ntimesMFBwd;

    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> m_connectionform;
    Array<OneD, Array<OneD, Array<OneD, Array<OneD, NekDouble>>>>
        m_curvatureform;

    SpatialDomains::GeomMMF m_MMFdir;
    SpatialDomains::GeomMMF m_DivMMFdir;

    void GetLOCALMovingframes(
        Array<OneD, Array<OneD, NekDouble>> &movingframes);

    void GetLOCALFIBERMovingFrames(
        MultiRegions::ExpListSharedPtr &field,
        Array<OneD, Array<OneD, NekDouble>> &movingframes);

    void GetLOCALMovingframes1D(
        Array<OneD, Array<OneD, NekDouble>> &movingframes);

    void GetLOCALMovingframes3D(
        Array<OneD, Array<OneD, NekDouble>> &movingframes);

    void ConstructAnisotropicFrames(
        const Array<OneD, const NekDouble> &AniDirection,
        Array<OneD, Array<OneD, NekDouble>> &movingframes,
        Array<OneD, NekDouble> &AniIndex);

    void RebuildMovingFrames(const int K1, const int K2, const int K3,
                             Array<OneD, Array<OneD, NekDouble>> &distance,
                             Array<OneD, Array<OneD, NekDouble>> &movingframes);

    void CheckMovingFrames(
        const Array<OneD, const Array<OneD, NekDouble>> &movingframes);

    void ComputeCurl(const Array<OneD, const Array<OneD, NekDouble>> &inarray,
                     Array<OneD, Array<OneD, NekDouble>> &outarray);

    void ComputeExactDivMF(const Array<OneD, const int> &Activation,
                           Array<OneD, Array<OneD, NekDouble>> &ExactDivMF);

    void ComputeExactCurlMF(const Array<OneD, const int> &Activation,
                            Array<OneD, Array<OneD, NekDouble>> &ExactCurlMF);

    void TestDivMF(
        const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
        const Array<OneD, const int> &Activation);

    void TestCurlMF(
        const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
        const Array<OneD, const int> &Activation);

    void ConstructSphericalMF(Array<OneD, Array<OneD, NekDouble>> &SphereMF,
                              Array<OneD, int> &SphereMFActivate);

    void ConstructPseudosphericalMF(
        Array<OneD, Array<OneD, NekDouble>> &PsedoSphereMF,
        Array<OneD, int> &PseudosphereMFActivate);

    void ConstructEllipticalMF(
        const int Radx, const int Rady, const int Radz,
        Array<OneD, Array<OneD, NekDouble>> &EllipticalMF,
        Array<OneD, int> &EllipticalMFActivate);

    inline void wait_on_enter()
    {
        std::string dummy;
        std::cout << "Enter to continue..." << std::endl;
        std::getline(std::cin, dummy);
    }

    SOLVER_UTILS_EXPORT void vector_to_vcoeff(
        const Array<OneD, const Array<OneD, NekDouble>> &vector,
        const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
        Array<OneD, Array<OneD, NekDouble>> &vcoeff);

    SOLVER_UTILS_EXPORT void vcoeff_to_vector(
    const Array<OneD, const Array<OneD, NekDouble>> &vcoeff,
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    Array<OneD, Array<OneD, NekDouble>> &vector);

    SOLVER_UTILS_EXPORT SpatialDomains::GeomMMF FindMMFdir(
        std::string MMFdirStr);

    SOLVER_UTILS_EXPORT void SetUpMovingFrames(
        const SpatialDomains::GeomMMF MMFdir,
        const Array<OneD, const Array<OneD, NekDouble>> &Anisotropy,
        Array<OneD, Array<OneD, NekDouble>> &movingframes);

    SOLVER_UTILS_EXPORT void ComputeAxisAlignedLOCALMovingframes(
        const Array<OneD, const Array<OneD, NekDouble>> &AxisMF,
        Array<OneD, Array<OneD, NekDouble>> &movingframes);

    SOLVER_UTILS_EXPORT void TestPseudosphericalConnection1form(
        const Array<OneD, const Array<OneD, NekDouble>> &PseudosphereMF,
        Array<OneD, int> &PseudosphereMFActivate);

    SOLVER_UTILS_EXPORT void TestSphericalConnection1form(
        const Array<OneD, const Array<OneD, NekDouble>> &SphereMF,
        Array<OneD, int> &SphereMFActivate);

    SOLVER_UTILS_EXPORT void TestPolarConnection1form(
        const Array<OneD, const Array<OneD, NekDouble>> &PolarMF,
        Array<OneD, int> &PolarMFActivate);

    SOLVER_UTILS_EXPORT void ConstructPolarMF(
        Array<OneD, Array<OneD, NekDouble>> &PolarMF,
        Array<OneD, int> &PolarMFActivate);

    // SOLVER_UTILS_EXPORT void TestPolarConnectionForm(
    //     Array<OneD, Array<OneD, NekDouble>> &PolarMF,
    //     Array<OneD, int> &PolarMFActivate);

    // SOLVER_UTILS_EXPORT void AlignMFtoVelocity(
    //     const NekDouble VelActivationTol,
    //     const Array<OneD, const NekDouble> &velocity,
    //     const Array<OneD, const NekDouble> &velmag,
    //     const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    //     Array<OneD, Array<OneD, NekDouble>> &MF1st, Array<OneD, int>
    //     &Activated, Array<OneD, NekDouble> &ActivationIntensity);

    SOLVER_UTILS_EXPORT void AlignMFtoVelocity(
        const Array<OneD, const int> &Activated,
        const Array<OneD, const NekDouble> &velocity,
        const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
        Array<OneD, Array<OneD, NekDouble>> &MF1st);

    SOLVER_UTILS_EXPORT void ActivateRegion(
        const Array<OneD, const int> &ZoneActivation,
        const Array<OneD, const Array<OneD, NekDouble>> &MF1st,
        Array<OneD, Array<OneD, NekDouble>> &MF1stnew,
        Array<OneD, int> &Activated);

    SOLVER_UTILS_EXPORT int UpdatebyLinearTimeIntegration(
        const Array<OneD, const int> &Activated,
        const Array<OneD, const Array<OneD, NekDouble>> &MF1st,
        Array<OneD, int> &ActivatedHistory,
        Array<OneD, Array<OneD, NekDouble>> &MF1stAlignedInt);

    SOLVER_UTILS_EXPORT void UpdateMF1st(
        const Array<OneD, const int> &Activated,
        const Array<OneD, const Array<OneD, NekDouble>> &MF1st,
        const Array<OneD, const NekDouble> &velmag,
        Array<OneD, NekDouble> &VelmagHistory,
        Array<OneD, Array<OneD, NekDouble>> &MF1stAligned,
        Array<OneD, int> &ActivatedHistory);

    SOLVER_UTILS_EXPORT int NewValueReplacerElementWise(
        const EvalType EType, const Array<OneD, const NekDouble> &Intensity,
        const Array<OneD, const int> &Activated,
        const Array<OneD, const Array<OneD, NekDouble>> &MF1st,
        Array<OneD, NekDouble> &IntensityHistory,
        Array<OneD, int> &ActivatedHistory,
        Array<OneD, Array<OneD, NekDouble>> &MF1stHistory);

    SOLVER_UTILS_EXPORT int NewValueReplacerPointWise(
        const EvalType EType, const Array<OneD, const NekDouble> &Intensity,
        const Array<OneD, const int> &Activated,
        const Array<OneD, const Array<OneD, NekDouble>> &MF1st,
        Array<OneD, NekDouble> &IntensityHistory,
        Array<OneD, NekDouble> &ActivatedHistory,
        Array<OneD, Array<OneD, NekDouble>> &MF1stHistory);

    SOLVER_UTILS_EXPORT int WeakerValueReplacer(
        const Array<OneD, const NekDouble> &Intensity,
        const Array<OneD, const int> &Activated,
        const Array<OneD, const Array<OneD, NekDouble>> &MF1st,
        Array<OneD, NekDouble> &IntensityHistory,
        Array<OneD, NekDouble> &ActivatedHistory,
        Array<OneD, Array<OneD, NekDouble>> &MF1stHistory);

    SOLVER_UTILS_EXPORT int StrongerValueReplacer(
        const Array<OneD, const NekDouble> &Intensity,
        const Array<OneD, const int> &Activated,
        const Array<OneD, const Array<OneD, NekDouble>> &MF1st,
        Array<OneD, NekDouble> &IntensityHistory,
        Array<OneD, NekDouble> &ActivatedHistory,
        Array<OneD, Array<OneD, NekDouble>> &MF1stHistory);

    SOLVER_UTILS_EXPORT void ValidateNewFrames(
        const NekDouble AdaptNewFramesToldt,
        const Array<OneD, const NekDouble> &vecdiff,
        const Array<OneD, const Array<OneD, NekDouble>> &MF1st,
        Array<OneD, Array<OneD, NekDouble>> &MF1stnew,
        Array<OneD, int> &Activated);

    SOLVER_UTILS_EXPORT void ComputencdotMF(
        const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
        Array<OneD, Array<OneD, NekDouble>> &ncdotMFFwd,
        Array<OneD, Array<OneD, NekDouble>> &ncdotMFBwd, const int Verbose = 0);

    SOLVER_UTILS_EXPORT void ComputentimesMF(
        const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
        Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &ntimesMFFwd,
        Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &ntimesMFBwd,
        Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &ntimes_ntimesMFFwd,
        Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &ntimes_ntimesMFBwd);

    SOLVER_UTILS_EXPORT void ComputentimesMF(
        const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
        Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &ntimesMFFwd,
        Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &ntimesMFBwd);

    SOLVER_UTILS_EXPORT void ComputenperpcdotMF(
        const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
        Array<OneD, Array<OneD, NekDouble>> &nperpcdotMFFwd,
        Array<OneD, Array<OneD, NekDouble>> &nperpcdotMFBwd);

    SOLVER_UTILS_EXPORT void ComputeRiemCrv(
        const int dirX, const int dirY, const int dirZ,
        const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
        Array<OneD, Array<OneD, NekDouble>> &outarray);

    SOLVER_UTILS_EXPORT void ComputeRelacc(
        const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
        Array<OneD, Array<OneD, NekDouble>> &outarray);

    SOLVER_UTILS_EXPORT void ComputeRelaccOmega(
        const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
        Array<OneD, Array<OneD, NekDouble>> &outarray);

    SOLVER_UTILS_EXPORT void ComputeSecondCovDeriv(
        const int dirX, const int dirY, const int dirZ,
        const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
        Array<OneD, Array<OneD, NekDouble>> &outarray);

    SOLVER_UTILS_EXPORT void ComputeLieBracket(
        const int dirX, const int dirY, const int dirZ,
        const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
        Array<OneD, Array<OneD, NekDouble>> &outarray);

    SOLVER_UTILS_EXPORT void ComputeLieBracket(
        const Array<OneD, const Array<OneD, NekDouble>> &velocity,
        const Array<OneD, const NekDouble> &direction1,
        const Array<OneD, const NekDouble> &direction2,
        const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
        Array<OneD, Array<OneD, NekDouble>> &outarray);

    SOLVER_UTILS_EXPORT void ComputeCovDeriv(
        const Array<OneD, const Array<OneD, NekDouble>> &velocity,
        const Array<OneD, const NekDouble> &direction,
        const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
        Array<OneD, Array<OneD, NekDouble>> &outarray);

    SOLVER_UTILS_EXPORT void ComputeCovDeriv(
        const Array<OneD, const NekDouble> &u1,
        const Array<OneD, const NekDouble> &u2,
        const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
        Array<OneD, Array<OneD, NekDouble>> &outarray);

    SOLVER_UTILS_EXPORT void ComputeCovDerivMF(
        const int dirX, const int dirY,
        const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
        Array<OneD, Array<OneD, NekDouble>> &outarray);

    SOLVER_UTILS_EXPORT Array<OneD, NekDouble> ComputeCovDiv(
        const Array<OneD, const NekDouble> &velvector,
        const Array<OneD, const Array<OneD, NekDouble>> &movingframes);

    SOLVER_UTILS_EXPORT Array<OneD, NekDouble> ComputeCovDiv(
        const Array<OneD, const Array<OneD, NekDouble>> &velocity,
        const Array<OneD, const Array<OneD, NekDouble>> &movingframes);

    SOLVER_UTILS_EXPORT Array<OneD, NekDouble> ComputeCovCurl(
        const Array<OneD, const NekDouble> &velvector,
        const Array<OneD, const Array<OneD, NekDouble>> &movingframes);

    SOLVER_UTILS_EXPORT Array<OneD, NekDouble> ComputeCovCurl(
        const Array<OneD, const Array<OneD, NekDouble>> &velocity,
        const Array<OneD, const Array<OneD, NekDouble>> &movingframes);

    // SOLVER_UTILS_EXPORT Array<OneD, NekDouble> ComputeCovariantDerivative(
    //     const int direction,
    //     const Array<OneD, const NekDouble> &ui,
    //     const Array<OneD, const NekDouble> &u2,
    //     const Array<OneD, const Array<OneD, NekDouble>> &movingframes);

    // SOLVER_UTILS_EXPORT Array<OneD, NekDouble> ComputeCovariantDivergence(
    //     const Array<OneD, const NekDouble> &vellong,
    //     const Array<OneD, const Array<OneD, NekDouble>> &movingframes);

    // SOLVER_UTILS_EXPORT Array<OneD, NekDouble> ComputeCovariantDivergence(
    //     const Array<OneD, const Array<OneD, NekDouble>> &velocity,
    //     const Array<OneD, const Array<OneD, NekDouble>> &movingframes);

    // SOLVER_UTILS_EXPORT Array<OneD, NekDouble> ComputeCovariantCurl(
    //     const Array<OneD, const NekDouble> &vellong,
    //     const Array<OneD, const Array<OneD, NekDouble>> &movingframes);

    // SOLVER_UTILS_EXPORT Array<OneD, NekDouble> ComputeCovariantCurl(
    //     const Array<OneD, const Array<OneD, NekDouble>> &velocity,
    //     const Array<OneD, const Array<OneD, NekDouble>> &movingframes);

    SOLVER_UTILS_EXPORT Array<OneD, NekDouble> ComputeDivSphericalCoord(
        const Array<OneD, const NekDouble> &vphi,
        const Array<OneD, const NekDouble> &vth);

    SOLVER_UTILS_EXPORT Array<OneD, NekDouble> ComputeDivSphericalCoord(
        const Array<OneD, const Array<OneD, NekDouble>> &velocity);

    SOLVER_UTILS_EXPORT Array<OneD, NekDouble> ComputeVecCdotNabla(
        const Array<OneD, const NekDouble> &vecC,
        const Array<OneD, const NekDouble> &vecA,
        const Array<OneD, const NekDouble> &vecB);

    SOLVER_UTILS_EXPORT Array<OneD, NekDouble> ComputeMFDivergence(
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    const Array<OneD, const Array<OneD, NekDouble>> &velocity);

    SOLVER_UTILS_EXPORT Array<OneD, NekDouble> ComputeEuclideanDivergence(
        const Array<OneD, const Array<OneD, NekDouble>> &velocity);

    SOLVER_UTILS_EXPORT Array<OneD, NekDouble> ComputeCurlSphericalCoord(
        const Array<OneD, const NekDouble> &inarrayphi,
        const Array<OneD, const NekDouble> &inarrayth);

    SOLVER_UTILS_EXPORT Array<OneD, NekDouble> ComputeEuclideanGradient(
        const Array<OneD, const NekDouble> &fn);

    SOLVER_UTILS_EXPORT Array<OneD, NekDouble> ComputeCovGrad(
        const Array<OneD, const NekDouble> &fn,
        const Array<OneD, const Array<OneD, NekDouble>> &movingframes);

    // SOLVER_UTILS_EXPORT Array<OneD, NekDouble> ComputeCovGrad1D(
    //     const Array<OneD, const NekDouble> &fn,
    //     const Array<OneD, const Array<OneD, NekDouble>> &movingframes);

    // SOLVER_UTILS_EXPORT Array<OneD, NekDouble> ComputeCovGrad2D(
    //     const Array<OneD, const NekDouble> &fn,
    //     const Array<OneD, const Array<OneD, NekDouble>> &movingframes);

    // SOLVER_UTILS_EXPORT Array<OneD, NekDouble> ComputeCovGrad3D(
    //     const Array<OneD, const NekDouble> &fn,
    //     const Array<OneD, const Array<OneD, NekDouble>> &movingframes);

    SOLVER_UTILS_EXPORT Array<OneD, NekDouble> ComputeCovJGrad(
        const Array<OneD, const NekDouble> &fn,
        const Array<OneD, const Array<OneD, NekDouble>> &movingframes);

    SOLVER_UTILS_EXPORT Array<OneD, NekDouble> ComputeEuclideanDiffusion(
        const Array<OneD, const NekDouble> &inarray);

    SOLVER_UTILS_EXPORT Array<OneD, NekDouble> ComputeMMFDiffusion(
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    const Array<OneD, const NekDouble> &inarray);

    SOLVER_UTILS_EXPORT Array<OneD, NekDouble> ComputeCovariantDiffusion(
        const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
        const Array<OneD, const NekDouble> &inarray,
        const DerivType DType = eCovariant);

    SOLVER_UTILS_EXPORT void ComputeDivMF(
        DerivType Dtype,
        const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
        Array<OneD, Array<OneD, NekDouble>> &DivMF, const int Verbose = 0);

    SOLVER_UTILS_EXPORT void ComputeEuclideanDivMF(
        const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
        Array<OneD, Array<OneD, NekDouble>> &DivMF);

    SOLVER_UTILS_EXPORT void ComputeCovariantDivMF(
        const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
        Array<OneD, Array<OneD, NekDouble>> &DivMF);

    SOLVER_UTILS_EXPORT void ComputeCurlMF(
        DerivType Dtype,
        const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
        Array<OneD, Array<OneD, NekDouble>> &CurlMF, const int Verbose = 0);

    SOLVER_UTILS_EXPORT void ComputeEuclideanCurlMF(
        const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
        Array<OneD, Array<OneD, NekDouble>> &CurlMF);

    SOLVER_UTILS_EXPORT void ComputeCovariantCurlMF(
        const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
        Array<OneD, Array<OneD, NekDouble>> &CurlMF);

    SOLVER_UTILS_EXPORT void ComputeCurl(
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray);

    SOLVER_UTILS_EXPORT Array<OneD, NekDouble> ComputeEuclideanCurl(
        const Array<OneD, const Array<OneD, NekDouble>> &physfield,
        const Array<OneD, const Array<OneD, NekDouble>> &movingframes);

    SOLVER_UTILS_EXPORT void ComputeMFtrace(
        const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
        Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &MFtraceFwd,
        Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &MFtraceBwd);

    SOLVER_UTILS_EXPORT void ComputeMFcdotSphericalCoord(const int dir);
    SOLVER_UTILS_EXPORT void ComputeMFcdotSphericalCoord(
        Array<OneD, Array<OneD, NekDouble>> &ThetacdotMF,
        Array<OneD, Array<OneD, NekDouble>> &PhicdotMF);

    SOLVER_UTILS_EXPORT void ComputeMFcdotSphericalCoord(
        const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
        Array<OneD, Array<OneD, NekDouble>> &ThetacdotMF,
        Array<OneD, Array<OneD, NekDouble>> &PhicdotMF);

    SOLVER_UTILS_EXPORT void VectorDotProd(
        const Array<OneD, const Array<OneD, NekDouble>> &v1,
        const Array<OneD, const Array<OneD, NekDouble>> &v2,
        Array<OneD, NekDouble> &v3);

    SOLVER_UTILS_EXPORT NekDouble
    VectorDotProd(const Array<OneD, const NekDouble> &v1,
                  const Array<OneD, const NekDouble> &v2);

    SOLVER_UTILS_EXPORT Array<OneD, NekDouble> VectorDiff(
        const Array<OneD, const Array<OneD, NekDouble>> &velocity1,
        const Array<OneD, const Array<OneD, NekDouble>> &velocity2);

    SOLVER_UTILS_EXPORT void MFDotProd(
        const Array<OneD, const Array<OneD, NekDouble>> &MF1,
        const Array<OneD, const Array<OneD, NekDouble>> &MF2,
        Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &outarray);

    SOLVER_UTILS_EXPORT void MFDotProduct(
        const Array<OneD, const Array<OneD, NekDouble>> &MF1,
        const Array<OneD, const Array<OneD, NekDouble>> &MF2,
        Array<OneD, Array<OneD, NekDouble>> &outarray);

    SOLVER_UTILS_EXPORT void VectorCrossProd(
        const Array<OneD, const Array<OneD, NekDouble>> &v1,
        const Array<OneD, const Array<OneD, NekDouble>> &v2,
        Array<OneD, Array<OneD, NekDouble>> &v3);

    Array<OneD, NekDouble> VectorCrossProdMF(
        const Array<OneD, const NekDouble> &v1,
        const Array<OneD, const NekDouble> &v2);

    SOLVER_UTILS_EXPORT void VectorCrossProd(
        const Array<OneD, const NekDouble> &v1,
        const Array<OneD, const NekDouble> &v2,
        Array<OneD, NekDouble> &outarray);

    SOLVER_UTILS_EXPORT void Cart_to_MF(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        const Array<OneD, const Array<OneD, NekDouble>> &MF1st,
        Array<OneD, NekDouble> &outarray1, Array<OneD, NekDouble> &outarray2);

    SOLVER_UTILS_EXPORT void Cart_to_MF(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        const Array<OneD, const Array<OneD, NekDouble>> &MF1st,
        Array<OneD, Array<OneD, NekDouble>> &outarray);

    SOLVER_UTILS_EXPORT void Sph_to_Cart(const NekDouble &x, const NekDouble &y,
                                         const NekDouble &z,
                                         const NekDouble &uth,
                                         const NekDouble &uphi, NekDouble &ux,
                                         NekDouble &uy, NekDouble &uz);

    SOLVER_UTILS_EXPORT void MF_to_Sph(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        const Array<OneD, const Array<OneD, NekDouble>> &MF1st,
        Array<OneD, Array<OneD, NekDouble>> &outarray);

    SOLVER_UTILS_EXPORT Array<OneD, NekDouble> CartesianToMovingframes(
        const Array<OneD, const Array<OneD, NekDouble>> movingframes,
        const Array<OneD, const Array<OneD, NekDouble>> &uvec,
        unsigned int field);

    SOLVER_UTILS_EXPORT void CartesianToMovingframes(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        const Array<OneD, const Array<OneD, NekDouble>> &MF1st,
        Array<OneD, Array<OneD, NekDouble>> &outarray);

    SOLVER_UTILS_EXPORT void CartesianToMovingframes(
        const Array<OneD, const NekDouble> &inarrayx,
        const Array<OneD, const NekDouble> &inarrayy,
        const Array<OneD, const NekDouble> &inarrayz,
        const Array<OneD, const Array<OneD, NekDouble>> &MF1st,
        Array<OneD, NekDouble> &outarray1, Array<OneD, NekDouble> &outarray2);

    SOLVER_UTILS_EXPORT void MovingframestoCartesian(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        const Array<OneD, const Array<OneD, NekDouble>> &MF1st,
        Array<OneD, Array<OneD, NekDouble>> &outarray);

    SOLVER_UTILS_EXPORT void MovingframestoCartesian(
        const Array<OneD, const NekDouble> &inarray1,
        const Array<OneD, const NekDouble> &inarray2,
        const Array<OneD, const Array<OneD, NekDouble>> &MF1st,
        Array<OneD, NekDouble> &outarrayx, Array<OneD, NekDouble> &outarrayy,
        Array<OneD, NekDouble> &outarrayz);

    // SOLVER_UTILS_EXPORT void MovingframesToSpherical(
    // const Array<OneD, const NekDouble> &inarray,
    // const Array<OneD, const Array<OneD, NekDouble>> &MF1st,
    // Array<OneD, NekDouble> &outtheta,
    // Array<OneD, NekDouble> &outphi);

    SOLVER_UTILS_EXPORT void MovingframesToSpherical(
        const Array<OneD, const NekDouble> &u1,
        const Array<OneD, const NekDouble> &u2,
        const Array<OneD, const Array<OneD, NekDouble>> &MF1st,
        Array<OneD, NekDouble> &outtheta, Array<OneD, NekDouble> &outphi);

    SOLVER_UTILS_EXPORT void MovingframesToSpherical(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        const Array<OneD, const Array<OneD, NekDouble>> &MF1st,
        Array<OneD, NekDouble> &outtheta, Array<OneD, NekDouble> &outphi);

    SOLVER_UTILS_EXPORT void ProjectionOntoMovingFrames(
        const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
        const Array<OneD, const Array<OneD, NekDouble>> &velocity,
        Array<OneD, Array<OneD, NekDouble>> &outarray);

    SOLVER_UTILS_EXPORT void ComputeMFtimesMF(
        const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
        Array<OneD, Array<OneD, NekDouble>> &outarray);

    SOLVER_UTILS_EXPORT void CartesianToNewSpherical(
        const NekDouble &x0j, const NekDouble &x1j, const NekDouble &x2j,
        NekDouble &sin_varphi, NekDouble &cos_varphi, NekDouble &sin_theta,
        NekDouble &cos_theta);

    SOLVER_UTILS_EXPORT void CartesianToPseudospherical(
        const NekDouble &x0j, const NekDouble &x1j, const NekDouble &x2j,
        NekDouble &sin_varphi, NekDouble &cos_varphi, NekDouble &theta,
        NekDouble &sech_theta, NekDouble &tanh_theta);

    // SOLVER_UTILS_EXPORT void CartesianToNewSpherical(
    //     const NekDouble x0j, const NekDouble x1j, const NekDouble x2j,
    //     NekDouble &sin_varphi, NekDouble &cos_varphi,
    //     NekDouble &sin_theta, NekDouble &cos_theta);

    // SOLVER_UTILS_EXPORT void CartesianToSpherical(
    //     const NekDouble x0j, const NekDouble x1j, const NekDouble x2j,
    //     NekDouble &sin_varphi, NekDouble &cos_varphi, NekDouble &sin_theta,
    //     NekDouble &cos_theta);

    SOLVER_UTILS_EXPORT void CartesianToElliptical(
        const NekDouble x0j, const NekDouble x1j, const NekDouble x2j,
        NekDouble &rad, NekDouble &sin_varphi, NekDouble &cos_varphi,
        NekDouble &sin_theta, NekDouble &cos_theta);

    SOLVER_UTILS_EXPORT void ComputeTangentUnitVector(
        Array<OneD, Array<OneD, NekDouble>> &TangentUnitVector);

    SOLVER_UTILS_EXPORT void ComputeSphericalVector(
        Array<OneD, Array<OneD, NekDouble>> &SphericalVector);

    SOLVER_UTILS_EXPORT void ComputeEllipsoidVector(
        Array<OneD, Array<OneD, NekDouble>> &EllipsoidVector);

    SOLVER_UTILS_EXPORT void ComputeSphericalTangentVector(
        Array<OneD, NekDouble> &Phi, Array<OneD, NekDouble> &Theta);

    SOLVER_UTILS_EXPORT void Compute2DConnection1form(
        const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
        Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &Connection1form);

    SOLVER_UTILS_EXPORT void Compute2DConnectionCurvature(
        const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
        Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &Connectionform,
        Array<OneD, Array<OneD, NekDouble>> &Curvatureform);

    // SOLVER_UTILS_EXPORT void ComputeConnection1form(
    //     const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    //     Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &Connectionform,
    //     const int Verbose=0);

    SOLVER_UTILS_EXPORT void Compute2DCurvatureForm(
        const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
        const Array<OneD, const Array<OneD, Array<OneD, NekDouble>>>
            &Connectionform,
        Array<OneD, Array<OneD, NekDouble>> &Curvatureform);

    SOLVER_UTILS_EXPORT void CheckJacobianError();

    Array<OneD, NekDouble> ComputeJacobianAvg(
        const Array<OneD, const NekDouble> &Jacobian);

    SOLVER_UTILS_EXPORT void ComputeMFJacobian(
        const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
        Array<OneD, Array<OneD, Array<OneD, Array<OneD, NekDouble>>>>
            &outarray);

    SOLVER_UTILS_EXPORT void GetCurvatureForm(
        const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
        const Array<OneD, const Array<OneD, Array<OneD, NekDouble>>>
            &Connectionform,
        Array<OneD, Array<OneD, Array<OneD, Array<OneD, NekDouble>>>>
            &Curvatureform);

    SOLVER_UTILS_EXPORT void Test2DConnection1formPolar(
        const NekDouble InitPx, const NekDouble InitPy, const NekDouble InitPz,
        const Array<OneD, const int> &Activation,
        const Array<OneD, const Array<OneD, NekDouble>> &MF1st,
        const Array<OneD, const Array<OneD, Array<OneD, NekDouble>>>
            &Connection);

    SOLVER_UTILS_EXPORT void Test2DConnection1formSphere(
        const Array<OneD, const int> &Activation,
        const Array<OneD, const Array<OneD, NekDouble>> &MF1st,
        const Array<OneD, const Array<OneD, Array<OneD, NekDouble>>>
            &Connection);

    SOLVER_UTILS_EXPORT void Test2DConnection1formPseudosphere(
        const Array<OneD, const int> &Activation,
        const Array<OneD, const Array<OneD, NekDouble>> &MF1st,
        const Array<OneD, const Array<OneD, Array<OneD, NekDouble>>>
            &Connection);

    void Test2DConnectionCurvature(
        const NekDouble InitPx, const NekDouble InitPy, const NekDouble InitPz,
        const Array<OneD, const int> &Activation,
        const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
        const Array<OneD, const Array<OneD, Array<OneD, NekDouble>>>
            &Connection,
        const Array<OneD, const Array<OneD, NekDouble>> &Curvature);

    NekDouble ComputeVMV(
        const int indexi, const int indexj, const int indexk, const int knode,
        const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
        const Array<OneD,
                    const Array<OneD, Array<OneD, Array<OneD, NekDouble>>>>
            &MFJacobian);

    SOLVER_UTILS_EXPORT NekDouble
    AvgInt(const Array<OneD, const NekDouble> &inarray);
    SOLVER_UTILS_EXPORT NekDouble
    AvgAbsInt(const Array<OneD, const NekDouble> &inarray);
    SOLVER_UTILS_EXPORT NekDouble
    AbsIntegral(const Array<OneD, const NekDouble> &inarray);

    SOLVER_UTILS_EXPORT NekDouble
    FindAbsMaximum(const Array<OneD, const NekDouble> &inarray,
                   const Array<OneD, const int> &Activated);

    SOLVER_UTILS_EXPORT NekDouble
    FindAbsMaximumVector(const Array<OneD, const NekDouble> &inarray,
                         const Array<OneD, const int> &Activated);

    SOLVER_UTILS_EXPORT NekDouble Average(const Array<OneD, const NekDouble> &inarray);

    SOLVER_UTILS_EXPORT NekDouble
    RootMeanSquare(const Array<OneD, const NekDouble> &inarray,
                   const Array<OneD, const int> &Activated);

    SOLVER_UTILS_EXPORT NekDouble
    RootMeanSquareVector(const Array<OneD, const NekDouble> &inarray,
                         const Array<OneD, const int> &Activated);

    SOLVER_UTILS_EXPORT NekDouble RootMeanSquare(
        const Array<OneD, const NekDouble> &inarray, const int Ntot = 0);

    SOLVER_UTILS_EXPORT NekDouble
    RootMeanSquare(const Array<OneD, const int> &inarray, const int Ntot = 0);

    SOLVER_UTILS_EXPORT NekDouble VectorAvgMagnitude(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray);

    SOLVER_UTILS_EXPORT NekDouble
    CompareVector(const Array<OneD, const Array<OneD, NekDouble>> &vector1,
                  const Array<OneD, const Array<OneD, NekDouble>> &vector2);

    SOLVER_UTILS_EXPORT void DeriveVector(
        const Array<OneD, const NekDouble> &veccoeff1,
        const Array<OneD, const NekDouble> &veccoeff2,
        const Array<OneD, const NekDouble> &vecbase1,
        const Array<OneD, const NekDouble> &vecbase2,
        Array<OneD, Array<OneD, NekDouble>> &outarray);

    SOLVER_UTILS_EXPORT void DeriveVector(
        const Array<OneD, const Array<OneD, NekDouble>> &veccoeff,
        const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
        Array<OneD, Array<OneD, NekDouble>> &outarray);

    SOLVER_UTILS_EXPORT void DeriveNewFieldComponent(
        const Array<OneD, const Array<OneD, NekDouble>> &MF1,
        const Array<OneD, const NekDouble> &u1,
        const Array<OneD, const NekDouble> &v1,
        const Array<OneD, const Array<OneD, NekDouble>> &MF2,
        Array<OneD, NekDouble> &u2, Array<OneD, NekDouble> &v2);

    SOLVER_UTILS_EXPORT void GramSchumitz(
        const Array<OneD, const Array<OneD, NekDouble>> &v1,
        const Array<OneD, const Array<OneD, NekDouble>> &v2,
        Array<OneD, Array<OneD, NekDouble>> &outarray,
        bool KeepTheMagnitude = true);

    SOLVER_UTILS_EXPORT void BubbleSort(Array<OneD, NekDouble> &refarray,
                                        Array<OneD, NekDouble> &sortarray);

    SOLVER_UTILS_EXPORT void PlotCardiacFibre(
        const Array<OneD, const NekDouble> &fibre);

    SOLVER_UTILS_EXPORT void PlotProcessedCardiacFibre(
        const Array<OneD, const NekDouble> &fibre,
        const Array<OneD, const NekDouble> &fcdotk,
        const Array<OneD, const NekDouble> &AniConstruction);

    SOLVER_UTILS_EXPORT void PlotMovingFrames(
        const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
        const Array<OneD, const Array<OneD, NekDouble>> &velocity,
        const int nstep);

    SOLVER_UTILS_EXPORT void PlotTrajectoryMF(
        const Array<OneD, const int> &ActivatedHistory,
        const Array<OneD, const NekDouble> &fieldu,
        const Array<OneD, const Array<OneD, NekDouble>> &MF,
        const Array<OneD, const Array<OneD, Array<OneD, NekDouble>>>
            &MF1stConnection,
        const Array<OneD, const Array<OneD, NekDouble>> &Relacc,
        Array<OneD, NekDouble> &NoBoundaryZone, const int nstep);

    SOLVER_UTILS_EXPORT Array<OneD, NekDouble> ComputeLaplacianDiff(
        const Array<OneD, const NekDouble> &Laplacian,
        const Array<OneD, const NekDouble> &LaplacianNew);

    SOLVER_UTILS_EXPORT void TimeMapProcess();

    SOLVER_UTILS_EXPORT void ComputeMFTimeMap(
        const Array<OneD, const int> &ValidTimeMap,
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, int> &NewValidTimeMap,
        Array<OneD, Array<OneD, NekDouble>> &TMMF);

    SOLVER_UTILS_EXPORT Array<OneD, NekDouble> ComputeSpaceTime(
        const Array<OneD, const NekDouble> &TimeMap);

    SOLVER_UTILS_EXPORT void VectorCutOff(const NekDouble Tol,
                                          Array<OneD, NekDouble> &vector);

    SOLVER_UTILS_EXPORT void HelmSolveSmoothing(
        const NekDouble TimeMapSmoothL, Array<OneD, NekDouble> &outarray);

    SOLVER_UTILS_EXPORT void PlotHHD(
        const Array<OneD, const NekDouble> &MFvec,
        const Array<OneD, const NekDouble> &irrotational,
        const Array<OneD, const NekDouble> &incompressible,
        const Array<OneD, const NekDouble> &harmonic, const int nstep);

    SOLVER_UTILS_EXPORT void PlotMovingFrames(
        const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
        const int nstep = 0);

    SOLVER_UTILS_EXPORT void PlotMovingFrames(
        const int dir,
        const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
        const int nstep = 0);

    SOLVER_UTILS_EXPORT void PlotAnisotropyFiber(
        const Array<OneD, const NekDouble> &anifibre);

    SOLVER_UTILS_EXPORT void PlotFieldVector(
        const Array<OneD, const NekDouble> &field,
        const Array<OneD, const Array<OneD, NekDouble>> &velocity,
        const int nstep = 0);

    SOLVER_UTILS_EXPORT void PlotFieldVector(
        const Array<OneD, const NekDouble> &field,
        const Array<OneD, const NekDouble> &velocity, const int nstep = 0);

    SOLVER_UTILS_EXPORT void PlotJacobian(
        const Array<OneD, const NekDouble> &Jacobian,
        const Array<OneD, const NekDouble> &JacGradMag);

    SOLVER_UTILS_EXPORT
    void PlotDivMF(const Array<OneD, const NekDouble> &DivMF1,
                   const Array<OneD, const NekDouble> &DivMF2,
                   const Array<OneD, const NekDouble> &DivMF3,
                   const int nstep = 0);

    SOLVER_UTILS_EXPORT void PlotConnectionForm(
        const Array<OneD, const Array<OneD, Array<OneD, NekDouble>>>
            &connectionform,
        const int nstep = 0);

    SOLVER_UTILS_EXPORT void Plot2DConnectionCurvature(
        const Array<OneD, const Array<OneD, Array<OneD, NekDouble>>>
            &connectionform,
        const Array<OneD, const Array<OneD, NekDouble>> &curvatureform,
        const int nstep);

    SOLVER_UTILS_EXPORT void Plot2DConnectionForm(
        const Array<OneD, const Array<OneD, Array<OneD, NekDouble>>>
            &connectionform,
        const int nstep = 0);

    SOLVER_UTILS_EXPORT void Plot2DMeanGaussCurv(
        const Array<OneD, const Array<OneD, Array<OneD, NekDouble>>>
            &connectionform,
        const int nstep = 0);

    SOLVER_UTILS_EXPORT
    void PlotConnectionError(const Array<OneD, const NekDouble> &Werror,
                             const int nstep = 0);

    SOLVER_UTILS_EXPORT void Plot2DRiemanTensor(
        const Array<OneD, const Array<OneD, NekDouble>> &curvatureform,
        const int nstep = 0);

    SOLVER_UTILS_EXPORT void PlotCurvatureForm(
        const Array<OneD,
                    const Array<OneD, Array<OneD, Array<OneD, NekDouble>>>>
            &curvatureform,
        const int nstep = 0);

    SOLVER_UTILS_EXPORT void PlotRiemannTensor(
        const Array<OneD,
                    const Array<OneD, Array<OneD, Array<OneD, NekDouble>>>>
            &curvatureform,
        const int nstep = 0);

    SOLVER_UTILS_EXPORT void PlotRicciTensor(
        const Array<OneD,
                    const Array<OneD, Array<OneD, Array<OneD, NekDouble>>>>
            &curvatureform,
        const int nstep = 0);

    SOLVER_UTILS_EXPORT void ComputeWeakDGDivergence(
        const Array<OneD, const NekDouble> &inarray,
        const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
        const Array<OneD, const Array<OneD, NekDouble>> &velocity,
        Array<OneD, Array<OneD, NekDouble>> &outarray,
        const int SurfaceDivergence = 0);

    SOLVER_UTILS_EXPORT void ComputeWeakDGDivergence(
        const Array<OneD, const NekDouble> &inarray,
        const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
        const Array<OneD, const Array<OneD, NekDouble>> &velocity,
        Array<OneD, Array<OneD, NekDouble>> &outarray,
        const Array<OneD, const NekDouble> &SpuriousDiffusion);

    SOLVER_UTILS_EXPORT void ComputeWeakDGCurl(
        const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
        const Array<OneD, const Array<OneD, NekDouble>> &velocity,
        Array<OneD, Array<OneD, NekDouble>> &outarray,
        const int SurfaceDivergence = 0);

    SOLVER_UTILS_EXPORT Array<OneD, NekDouble> ComputeSurfaceDiv(
        const Array<OneD, const NekDouble> &surfaceNormal,
        const Array<OneD, const NekDouble> &velvector);

    SOLVER_UTILS_EXPORT Array<OneD, NekDouble> ComputeSpuriousDivergence(
        const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
        const Array<OneD, const NekDouble> &surfaceNormal,
        const Array<OneD, const NekDouble> &velvector);

    SOLVER_UTILS_EXPORT void ComputeveldotMF(
        const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
        const Array<OneD, const Array<OneD, NekDouble>> &velocity,
        Array<OneD, Array<OneD, NekDouble>> &veldotMF);

    /// Variable diffusivity
    NekDouble m_chi;
    NekDouble m_capMembrane;
    NekDouble m_conductivity;

    SOLVER_UTILS_EXPORT void WeakDGMMFDiff1D(
        const Array<OneD, const Array<OneD, NekDouble>> &MF1st,
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &qfield,
        Array<OneD, NekDouble> &outarray, const NekDouble time);

    SOLVER_UTILS_EXPORT void WeakDGMMFDiff(
        const Array<OneD, const Array<OneD, NekDouble>> &MF1st,
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &qfield,
        Array<OneD, NekDouble> &outarray, const NekDouble time);

    SOLVER_UTILS_EXPORT void GetufluxMMF(
        const Array<OneD, const Array<OneD, NekDouble>> &MF1st,
        const Array<OneD, const NekDouble> &ufield,
        Array<OneD, Array<OneD, NekDouble>> &ufluxFwd,
        Array<OneD, Array<OneD, NekDouble>> &ufluxBwd);

    SOLVER_UTILS_EXPORT void GetqfluxMMF(
        const Array<OneD, const Array<OneD, NekDouble>> &MF1st,
        const Array<OneD, const NekDouble> &ufield,
        Array<OneD, Array<OneD, NekDouble>> &qfield,
        Array<OneD, NekDouble> &qflux);

    SOLVER_UTILS_EXPORT void DiffusionScalarBoundary(
        const Array<OneD, const NekDouble> &Fwd, Array<OneD, NekDouble> &Bwd);

    SOLVER_UTILS_EXPORT void DiffusionVectorBoundary(
        const Array<OneD, const Array<OneD, NekDouble>> &MF1st,
        const Array<OneD, const Array<OneD, NekDouble>> &qFwd,
        const Array<OneD, const Array<OneD, NekDouble>> &qFwdVec,
        Array<OneD, Array<OneD, NekDouble>> &qBwd,
        Array<OneD, Array<OneD, NekDouble>> &qBwdVec);

    SOLVER_UTILS_EXPORT void WeakDGDivergence(
        const Array<OneD, const Array<OneD, NekDouble>> &InField,
        const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
        const Array<OneD, const Array<OneD, NekDouble>> &velocity,
        Array<OneD, Array<OneD, NekDouble>> &OutField);

    SOLVER_UTILS_EXPORT void WeakDGGradientVector(
        const Array<OneD, const NekDouble> &inarray,
        const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
        Array<OneD, NekDouble> &outarray, const int SurfaceGradient = 0);

    SOLVER_UTILS_EXPORT void WeakDGDirectionalDeriv(
        const int dir, const Array<OneD, const NekDouble> &InField,
        const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
        Array<OneD, NekDouble> &OutField, const int SurfaceGradient);

    SOLVER_UTILS_EXPORT void WeakDGDirectionalDerivwithMF(
        const int dir, const Array<OneD, const NekDouble> &InField,
        const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
        Array<OneD, NekDouble> &OutField, const int SurfaceGradient);

    SOLVER_UTILS_EXPORT void WeakDGDirectionalDerivwithMF(
        const int dir, const int xj,
        const Array<OneD, const NekDouble> &InField,
        const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
        Array<OneD, NekDouble> &OutField, const int SurfaceGradient = 0);

    SOLVER_UTILS_EXPORT void WeakDGCurl(
        const Array<OneD, const Array<OneD, NekDouble>> &physfield,
        const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
        const Array<OneD, const Array<OneD, NekDouble>> &CrossProductMF,
        Array<OneD, NekDouble> &OutField, const int SurfaceCurl = 0);

    SOLVER_UTILS_EXPORT void WeakDGCurl(
        const Array<OneD, const Array<OneD, NekDouble>> &physfield,
        const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
        const Array<OneD, const Array<OneD, NekDouble>> &CrossProductMF,
        Array<OneD, NekDouble> &OutField,
        const Array<OneD, const NekDouble> &SpuriousDiffusion);

    SOLVER_UTILS_EXPORT Array<OneD, NekDouble> ComputeLaplacianSphericalCoord(
        const Array<OneD, const NekDouble> &inarray);

    SOLVER_UTILS_EXPORT void ComputeHHD(
        const Array<OneD, const int> &Activation,
        const Array<OneD, const Array<OneD, NekDouble>> &velocity,
        Array<OneD, NekDouble> &irrotational,
        Array<OneD, NekDouble> &incompressible,
        Array<OneD, NekDouble> &harmonic);

    SOLVER_UTILS_EXPORT void ComputeMFHHD(
        const Array<OneD, const int> &Activation,
        const Array<OneD, const NekDouble> &MF0,
        Array<OneD, NekDouble> &irrotational,
        Array<OneD, NekDouble> &incompressible,
        Array<OneD, NekDouble> &harmonic);

    SOLVER_UTILS_EXPORT void ComputeVarCoeff1D(
        const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
        StdRegions::VarCoeffMap &varcoeff);

    SOLVER_UTILS_EXPORT void ComputeVarCoeff2D(
        const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
        StdRegions::VarCoeffMap &varcoeff);

    SOLVER_UTILS_EXPORT void ComputeVarCoeff2DDxDyDz(
    const Array<OneD, const NekDouble> &epsilon,
    StdRegions::VarCoeffMap &varcoeff);

    SOLVER_UTILS_EXPORT int CountActivated(
        const Array<OneD, const int> &ActivatedHistory);

    SOLVER_UTILS_EXPORT Array<OneD, NekDouble> ElementwiseAverage(
        const Array<OneD, const NekDouble> &inarray);

    SOLVER_UTILS_EXPORT void SphericalToMovingFrames(
        const Array<OneD, const NekDouble> &inarrayth,
        const Array<OneD, const NekDouble> &inarrayphi,
        Array<OneD, Array<OneD, NekDouble>> &physfield);

    SOLVER_UTILS_EXPORT void SphericalToEuclidean(
        const Array<OneD, const NekDouble> &inarrayth,
        const Array<OneD, const NekDouble> &inarrayphi,
        Array<OneD, Array<OneD, NekDouble>> &velocity);

    SOLVER_UTILS_EXPORT void ComputeGradient(
        const Array<OneD, const Array<OneD, NekDouble>> &MF1st,
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &velocity,
        Array<OneD, NekDouble> &velmag);

    SOLVER_UTILS_EXPORT void ComputeGradient(
        const Array<OneD, const Array<OneD, NekDouble>> &MF1st,
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &velocity,
        Array<OneD, NekDouble> &velmag);

    SOLVER_UTILS_EXPORT void ComputeGradientDirect(
        const Array<OneD, const Array<OneD, NekDouble>> &MF1st,
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &qfield);

    SOLVER_UTILS_EXPORT void ComputeGradientWeak(
        const Array<OneD, const Array<OneD, NekDouble>> &MF1st,
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &qfield);

    SOLVER_UTILS_EXPORT void ComputeGradientSphericalCoord(
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray);

    SOLVER_UTILS_EXPORT void CheckMeshErr(
        const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
        const Array<OneD, const Array<OneD, NekDouble>> &velocity);

    SOLVER_UTILS_EXPORT void CheckMeshErr();

    SOLVER_UTILS_EXPORT void CheckMFOrientation(
        Array<OneD, Array<OneD, NekDouble>> &movingframes,
        Array<OneD, NekDouble> &origin = NullNekDouble1DArray);

    SOLVER_UTILS_EXPORT Array<OneD, NekDouble> WeakDGMMFLDG(
        const int var, const Array<OneD, const NekDouble> &inarray,
        const NekDouble time = 0.0);

    SOLVER_UTILS_EXPORT Array<OneD, NekDouble> WeakDGMMFLDG1D(
        const int var, const Array<OneD, const NekDouble> &inarray,
        const NekDouble time = 0.0);

    SOLVER_UTILS_EXPORT Array<OneD, NekDouble> WeakDGMMFLDG2D(
        const int var, const Array<OneD, const NekDouble> &inarray,
        const NekDouble time = 0.0);

    SOLVER_UTILS_EXPORT Array<OneD, NekDouble> WeakDGMMFLDG2D(
        const int var, const Array<OneD, const NekDouble> &inarray,
        const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
        const NekDouble time);

    SOLVER_UTILS_EXPORT Array<OneD, NekDouble> WeakDGMMFLDG3D(
        const int var, const Array<OneD, const NekDouble> &inarray);

    SOLVER_UTILS_EXPORT Array<OneD, NekDouble> WeakDGMMFLDGFiber(
        MultiRegions::ExpListSharedPtr &fiberfield,
        const Array<OneD, const NekDouble> &inarray,
        const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
        const Array<OneD, const Array<OneD, NekDouble>> &ncdotMFFwd,
        const Array<OneD, const Array<OneD, NekDouble>> &ncdotMFBwd,
        const Array<OneD, const Array<OneD, NekDouble>> &DivMF,
        const NekDouble time);

    Array<OneD, NekDouble> WeakDGMMFirstLDG(
        const int var,
        const Array<OneD, const Array<OneD, NekDouble>> movingframes,
        const Array<OneD, const NekDouble> &inarray);

    SOLVER_UTILS_EXPORT void WeakDGMMFLaplacian(
        const int var, const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray);

    SOLVER_UTILS_EXPORT void WeakDGMMFirstLaplacian(
        const int var,
        const Array<OneD, const Array<OneD, NekDouble>> movingframes,
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray);

    SOLVER_UTILS_EXPORT void WeakDGMMFDiffusion(
        const int var, const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray, const NekDouble time = 0.0);

    Array<OneD, NekDouble> ComputeufluxMMF(
        const int var, const Array<OneD, const NekDouble> &ufield);

    Array<OneD, NekDouble> ComputeqfluxMMF(
        const int var, const Array<OneD, const NekDouble> &ufield,
        const Array<OneD, const Array<OneD, NekDouble>> &ncdotMFFwd,
        const Array<OneD, const Array<OneD, NekDouble>> &ncdotMFBwd,
        const Array<OneD, const Array<OneD, NekDouble>> &qfield);

    SOLVER_UTILS_EXPORT void ApplyScalarBCs(
        const int var, const Array<OneD, const NekDouble> &Fwd,
        Array<OneD, NekDouble> &penaltyflux);

    SOLVER_UTILS_EXPORT void ApplyVectorBCs(
        const std::size_t var, const std::size_t dir,
        const Array<OneD, const NekDouble> &qFwd,
        const Array<OneD, const Array<OneD, NekDouble>> &ncdotMFFwd,
        Array<OneD, NekDouble> &penaltyflux);

    void NumericalCurlFlux(
        const Array<OneD, const Array<OneD, NekDouble>> &physfield,
        Array<OneD, Array<OneD, NekDouble>> &numfluxFwd,
        Array<OneD, Array<OneD, NekDouble>> &numfluxBwd);

    void ComputeNtimesF12(const Array<OneD, Array<OneD, NekDouble>> &Fwd,
                          const Array<OneD, Array<OneD, NekDouble>> &Bwd,
                          const Array<OneD, const NekDouble> &im1Fwd,
                          const Array<OneD, const NekDouble> &im1Bwd,
                          const Array<OneD, const NekDouble> &im2Fwd,
                          const Array<OneD, const NekDouble> &im2Bwd,
                          Array<OneD, NekDouble> &outarrayFwd,
                          Array<OneD, NekDouble> &outarrayBwd);

    Array<OneD, NekDouble> GetNormalVelocity(
        const Array<OneD, const Array<OneD, NekDouble>> &velocity);

    // MultiRegions::ContField2DSharedPtr m_contField;

    SOLVER_UTILS_EXPORT void CheckSphereConnection();
    SOLVER_UTILS_EXPORT void CheckPolarConnection();
    SOLVER_UTILS_EXPORT Array<OneD, NekDouble> ComputeRossyHaurtizVelocity();

    SOLVER_UTILS_EXPORT NekDouble
    ComputeDistance(const Array<OneD, const NekDouble> &p1,
                    const Array<OneD, const NekDouble> &p2);

    SOLVER_UTILS_EXPORT void Checkpoint_Output_1D(
        const int n, const Array<OneD, const NekDouble> &trjlength,
        const Array<OneD, const NekDouble> &field,
        const Array<OneD, const NekDouble> &TimeMap);

    SOLVER_UTILS_EXPORT Array<OneD, NekDouble> GaussianPulse(
        const Array<OneD, const NekDouble> &dudt,
        const Array<OneD, const NekDouble> &xp,
        const Array<OneD, const NekDouble> &yp,
        const Array<OneD, const NekDouble> &zp, const NekDouble Grad);

    SOLVER_UTILS_EXPORT Array<OneD, NekDouble> SmoothCircle(
        const NekDouble &v_amp, const NekDouble &m_pis, const NekDouble &m_px,
        const NekDouble &m_py, const NekDouble &m_pz, const NekDouble &m_pr);

    SOLVER_UTILS_EXPORT Array<OneD, NekDouble> SmoothLine(
        const NekDouble &v_amp, const NekDouble &m_pis, const NekDouble &m_px,
        const NekDouble &m_pr);

    SOLVER_UTILS_EXPORT void ComputeTimeMap(const NekDouble time,
                               const NekDouble urest,
                               const Array<OneD, const NekDouble> &field,
                               const Array<OneD, const NekDouble> &dudt,
                               const Array<OneD, const int> &ValidTimeMap,
                               Array<OneD, NekDouble> &dudtHistory,
                               Array<OneD, NekDouble> &TimeMap);

    SOLVER_UTILS_EXPORT void PlotTimeMap(
                            const Array<OneD, const int> &ValidTimeMap,
                            const Array<OneD, const Array<OneD, NekDouble>> &AniStrength,
                            const Array<OneD, const NekDouble> &TimeMap,
                            const int nstep);

    SOLVER_UTILS_EXPORT void Checkpoint_Output_Error(
        const int n, const Array<OneD, const NekDouble> &field,
        const Array<OneD, const NekDouble> &exactsoln);

    SOLVER_UTILS_EXPORT void GetFwdBwdMFTrace(
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &Fwd, Array<OneD, NekDouble> &Bwd);
};

// Shared pointer to an MMFSystem class
typedef boost::shared_ptr<MMFSystem> MMFSystemSharedPtr;

} // namespace SolverUtils
} // namespace Nektar

#endif
