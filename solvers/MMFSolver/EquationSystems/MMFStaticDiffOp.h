///////////////////////////////////////////////////////////////////////////////
//
// File MMFStaticDiffOp.h
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
// Description: MMF advection solve routines
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_MMFSOLVER_EQUATIONSYSTEMS_MMFSTATICDIFFOP_H
#define NEKTAR_SOLVERS_MMFSOLVER_EQUATIONSYSTEMS_MMFSTATICDIFFOP_H

#include <SolverUtils/AdvectionSystem.h>
#include <SolverUtils/MMFSystem.h>
#include <SolverUtils/RiemannSolvers/RiemannSolver.h>
#include <SolverUtils/UnsteadySystem.h>

#include <SolverUtils/Diffusion/Diffusion.h>

using namespace Nektar::SolverUtils;

enum FluxType
{
    euflux,
    eqflux,
};

const char *const FluxTypeMap[] = {
    "uflux",
    "qflux",
};

enum TestType
{
    eTestCovariantPlane,
    eTestCovariantSphere,
    eTestRelAcc,
    eTestRiemannCurv,
    eTestDivPlane,
    eTestDivSphericMetId,
    eTestDivRossbyMetId,
    eTestDivRossbyFlow,
    eTestCurlPlaneFlow,
    eTestCurlSphericMetId,
    eTestCurlRossbyMetId,
    eTestCurlRossbyFlow,
    eTestGradPlaneFlow,
    eTestGradSphericMetID,
    eTestGradRossbyFlow,
    eTestLaplacianPlaneFlow,
    eTestLaplacianMetricID,
    eTestLaplacianSimpleFlow,
    eTestLaplacianRossbyFlow,
    SIZE_TestType ///< Length of enum list
};

const char *const TestTypeMap[] = {"TestCovariantPlane",
                                   "TestCovariantSphere",
                                   "TestRelAcc",
                                   "TestRiemannCurv",
                                   "TestDivPlane",
                                   "TestDivSphericMetId",
                                   "TestDivRossbyMetId",
                                   "TestDivRossbyFlow",
                                   "TestCurlPlaneFlow",
                                   "TestCurlSphericMetId",
                                   "TestCurlRossbyMetId",
                                   "TestCurlRossbyFlow",
                                   "TestGradPlaneFlow",
                                   "TestGradSphericMetId",
                                   "TestGradRossbyFlow",
                                   "TestLaplacianPlaneFlow",
                                   "TestLaplacianMetricID",
                                   "TestLaplacianSimpleFlow",
                                   "TestLaplacianRossbyFlow"
};

namespace Nektar
{
class MMFStaticDiffOp : public SolverUtils::MMFSystem
{
public:
    friend class MemoryManager<MMFStaticDiffOp>;

    static SolverUtils::EquationSystemSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const SpatialDomains::MeshGraphSharedPtr &pGraph)
    {
        SolverUtils::EquationSystemSharedPtr p =
            MemoryManager<MMFStaticDiffOp>::AllocateSharedPtr(pSession, pGraph);
        p->InitObject();
        return p;
    }

    /// Name of class
    static std::string className;

    TestType m_TestType;

    /// Initialise the object
    virtual void v_InitObject(bool DeclareFields = true) override;

    virtual void v_DoSolve() override;

    /// Destructor
    virtual ~MMFStaticDiffOp();

protected:
    NekDouble m_advx, m_advy, m_advz;
    NekDouble m_waveFreq, m_RotAngle;
    NekDouble m_alpha, m_u0, m_Omega, m_H0, m_Hvar; // TestSteadyZonal
    NekDouble m_g, m_K, m_W;                        // Gravity

    NekDouble m_Diffbeta; //  \beta_e for LDG
    NekDouble m_Diffeta;  // \eta for LDG
    NekDouble m_Diffhe;   // h_e for LDG

    NekDouble m_Mass0;
    int m_graddir;

    Array<OneD, NekDouble> m_Phi;
    Array<OneD, NekDouble> m_Theta;

    SolverUtils::DiffusionSharedPtr m_diffusion;
    SolverUtils::RiemannSolverSharedPtr m_riemannSolver;

    /// Advection velocity
    Array<OneD, Array<OneD, NekDouble>> m_velocity;
    Array<OneD, NekDouble> m_velphi;
    Array<OneD, NekDouble> m_veltheta;

    Array<OneD, NekDouble> m_traceVn;

    // Array<OneD, Array<OneD, NekDouble>> m_veldotMF;
    // Array<OneD, Array<OneD, NekDouble>> m_movingframes_VelMF;

    Array<OneD, Array<OneD, NekDouble>> m_CrossProductMF;

    // Plane (used only for Discontinous projection
    //        with 3DHomogenoeus1D expansion)
    int m_planeNumber;

    /// Session reader
    MMFStaticDiffOp(const LibUtilities::SessionReaderSharedPtr &pSession,
                    const SpatialDomains::MeshGraphSharedPtr &pGraph);

    Array<OneD, NekDouble> SetInitialU();
    void SetInitialU(Array<OneD, NekDouble> &InitialUth,
                     Array<OneD, NekDouble> &InitialUphi);

    // void ComputeMaxlength();
    
    /// Evaluate the flux at each solution point
    // void GetFluxVector(const Array<OneD, Array<OneD, NekDouble>> &physfield,
    //                    Array<OneD, Array<OneD, Array<OneD, NekDouble>>>
    //                    &flux);

    /// Evaluate the flux at each solution point using dealiasing
    void GetFluxVectorDeAlias(
        const Array<OneD, Array<OneD, NekDouble>> &physfield,
        Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &flux);

    /// Compute the RHS
    void DoOdeRhs(const Array<OneD, const Array<OneD, NekDouble>> &inarray,
                  Array<OneD, Array<OneD, NekDouble>> &outarray,
                  const NekDouble time);

    /// Compute the projection
    void DoOdeProjection(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time);

    void DoImplicitSolve(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray, NekDouble time,
        NekDouble lambda);

    // Compute the velocity vector
    void EvaluateCovariantVelocity(
        const NekDouble alpha,
        Array<OneD, NekDouble> &velphi,
        Array<OneD, NekDouble> &veltheta,
        Array<OneD, Array<OneD, NekDouble>> &velocity);

    void EvaluateAdvectionVelocity(
        Array<OneD, Array<OneD, NekDouble>> &velocity);

    Array<OneD, NekDouble> ComputeExactRelAcc();

    void ComputeExactRiemanCvt(
        Array<OneD, NekDouble> &ExactR1212, 
        Array<OneD, NekDouble> &ExactR2121);

    void ComputeExactLieBracket(
    Array<OneD, NekDouble> &LB1212, 
    Array<OneD, NekDouble> &LB2121);

    // Compute Covariant Derivative of Velocity vector
    void CheckErrorCovariantDeriv(
    const Array<OneD, const NekDouble> &velphi,
    const Array<OneD, const NekDouble> &veltheta,
    const Array<OneD, const Array<OneD, NekDouble>> &velvector);

    void CheckErrorRelAcc(
        const Array<OneD, const Array<OneD, NekDouble>> &movingframes);

    void CheckErrorRiemCrv(
        const Array<OneD, const Array<OneD, NekDouble>> &movingframes);

    Array<OneD, NekDouble> ComputeVecError(
        const Array<OneD, const NekDouble> &ExactCovDerivTheta, 
        const Array<OneD, const NekDouble> &ExactCovDerivPhi,
        const Array<OneD, const Array<OneD, NekDouble>> &CovDerivSph,
        Array<OneD, NekDouble> &Errortheta,
        Array<OneD, NekDouble> &Errorphi);

    void PlotCovError(
        const int dir,
        const Array<OneD, const NekDouble> &Errtheta,
        const Array<OneD, const NekDouble> &Errphi,
        const Array<OneD, const NekDouble> &ErrVec);

    // Compute Divergence of Velocity vector
    void CheckErrorDivergence(const Array<OneD, const NekDouble> &inarray);

    // Compute the curl of a vector (vth, vphi)
    void CheckErrorCurl(const Array<OneD, const NekDouble> &inarrayth,
                        const Array<OneD, const NekDouble> &inarraypphi);

    // Compute Gradient of a function u
    void CheckErrorGradient(const Array<OneD, const NekDouble> &inarray);

    // Compute the Laplacition of a function u
    void CheckErrorLaplacian(const Array<OneD, const NekDouble> &inarray);

    // Compute the Helmsolver of a function u
    void CheckErrorHelmSolve(const Array<OneD, const NekDouble> &inarray);

    Array<OneD, NekDouble> ComputeExactValue(TestType Ttype);

    // Return the exact covariant derivative for the given test velocity vector
    // outarray[0]: component in the theta-direction
    // outarray[1]: component in the phi-direction
    
    void EvaluateExactCovariantDeriv(
        const int direction, 
        Array<OneD, NekDouble> &outtheta, 
        Array<OneD, NekDouble> &outphi);

    // Compute arclenght of the surface at zlebel
    // NekDouble ComputeCirculatingArclength(const NekDouble zlevel,
    //                                       const NekDouble Rhs);

    /// Get the normal velocity
    // void GetNormalVelocity(Array<OneD, NekDouble> &traceVn);

    void ComputeNablaCdotVelocity(Array<OneD, NekDouble> &vellc);

    void ComputeNablaVmcdotMF(
        const Array<OneD, const Array<OneD, NekDouble>> &veldotMF,
        Array<OneD, NekDouble> &nablavncdotMF);

    void AdvectionBellPlane(Array<OneD, NekDouble> &outfield);
    void AdvectionBellSphere(Array<OneD, NekDouble> &outfield);

    void Test2Dproblem(const NekDouble time, Array<OneD, NekDouble> &outfield);
    void Test3Dproblem(const NekDouble time, Array<OneD, NekDouble> &outfield);

    void Checkpoint_Err(
        const NekDouble time,
        const Array<OneD, const Array<OneD, NekDouble>> &fieldphys);

    void PlotVelocityVector();
    void Wait_On_Enter();

    void SetBoundaryConditions(NekDouble time);

    void ExplicitDiffusion(
        const int nConvectiveFields,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
        const Array<OneD, Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray);

    void Getuflux(const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                  const Array<OneD, Array<OneD, NekDouble>> &ufield,
                  Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &uflux);

    void Getqflux(const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                  const Array<OneD, Array<OneD, NekDouble>> &ufield,
                  Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &qfield,
                  Array<OneD, Array<OneD, NekDouble>> &qflux);

    void DiffusionBoundaryConditions(FluxType ftype, const int var,
                                     const Array<OneD, const NekDouble> &Fwd,
                                     Array<OneD, NekDouble> &Bwd,
                                     int direction = 0);

    // void GetFluxVector(const int i, const int j,
    //                    const Array<OneD, Array<OneD, NekDouble>> &physfield,
    //                    Array<OneD, Array<OneD, NekDouble>> &derivatives,
    //                    Array<OneD, Array<OneD, NekDouble>> &flux);

    void ComputeVarCoeff(
        const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
        StdRegions::VarCoeffMap &varcoeff);

    /// Print Summary
    virtual void v_GenerateSummary(SolverUtils::SummaryList &s) override;

    virtual void v_SetInitialConditions(const NekDouble initialtime,
                                        bool dumpInitialConditions,
                                        const int domain) override;

    virtual void v_EvaluateExactSolution(unsigned int field,
                                         Array<OneD, NekDouble> &outfield,
                                         const NekDouble time) override;

private:
};
} // namespace Nektar

#endif
