///////////////////////////////////////////////////////////////////////////////
//
// File: MMFDiffusion.h
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
// Description: MMFDiffusion
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_MMFDIFFUSION_H
#define NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_MMFDIFFUSION_H

#include <SolverUtils/MMFSystem.h>
#include <SolverUtils/Diffusion/Diffusion.h>
#include <SolverUtils/UnsteadySystem.h>

using namespace Nektar::SolverUtils;

namespace Nektar
{

enum TestType
{
    eTestLineX,
    eTestLineY,
    eTestPlaneAni,
    eTestPlane,
    eTestPlaneEmbed,
    eTestPlaneNeumann,
    eTestCube,
    eTestLinearSphere,
    eTestNonlinearSphere,
    SIZE_TestType ///< Length of enum list
};

const char *const TestTypeMap[] = {
    "TestLineX", "TestLineY", "TestPlaneAni", "TestPlane", "TestPlaneEmbed", "TestPlaneNeumann", "TestRanvierNode",
    "TestCube", "TestLinearSphere", "TestNonlinearSphere",
};

enum SolverSchemeType
{
    eDefault,
    eMMFFirst,
    eTimeMap,
    SIZE_SolverSchemeType,
};

const char *const SolverSchemeTypeMap[] = {
    "Default",
    "MMFFirst",
    "TimeMap",
};

enum InitWaveType
{
    ePoint,
    eLeft,
    eBothEnds,
    eCenter,
    eLeftBottomCorner,
    eSpiralDock,
    SIZE_InitWaveType ///< Length of enum list
};

const char *const InitWaveTypeMap[] = {
    "Point", "Left", "BothEnd", "Center", "LeftBottomCorner", "SpiralDock",
};

enum FluxType
{
    euflux,
    eqflux,
};

const char *const FluxTypeMap[] = {
    "qflux",
    "uflux",
};

enum TimeMapType
{
    eActivated,
    eDeActivated,
    eProcessing,
    SIZE_TimeMapType ///< Length of enum list
};

const char *const TimeMapTypeMap[] = {
    "Activated",
    "DeActivated",
    "Processing",
};

/// A model for cardiac conduction.
class MMFDiffusion : public SolverUtils::MMFSystem
{
public:
    friend class MemoryManager<MMFDiffusion>;

    /// Creates an instance of this class
    static SolverUtils::EquationSystemSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const SpatialDomains::MeshGraphSharedPtr &pGraph)
    {
        SolverUtils::EquationSystemSharedPtr p =
            MemoryManager<MMFDiffusion>::AllocateSharedPtr(pSession, pGraph);
        p->InitObject();
        return p;
    }

    /// Name of class
    static std::string className;

    TestType m_TestType;
    SolverSchemeType m_SolverSchemeType;

    virtual void v_DoSolve() override;

    /// Desctructor
    virtual ~MMFDiffusion();

protected:
    bool m_useSpecVanVisc;
    NekDouble m_frequency;
    NekDouble m_d00, m_d11, m_d22;
    StdRegions::VarCoeffMap m_varcoeff;

    Array<OneD, Array<OneD, NekDouble>> m_phiemovingframes;
    SpatialDomains::GeomMMF m_phieMMFdir;

    int m_EmbededPlane;

    Array<OneD, NekDouble> m_d00vec;
    Array<OneD, NekDouble> m_d11vec;
    Array<OneD, NekDouble> m_d22vec;

    StdRegions::VarCoeffMap m_varcoeffXYZ;

    SolverUtils::DiffusionSharedPtr m_diffusion;
    SolverUtils::RiemannSolverSharedPtr m_riemannSolver;

    int m_Convectiven;
    TimeMapType m_TimeMap;

    NekDouble m_AniStrength;

    // Temperature parameter
    NekDouble m_TimeMapStart;
    NekDouble m_TimeMapEnd;

    NekDouble m_Helmtau;

    NekDouble m_Diffbeta, m_Diffeta, m_Diffhe; // h_e for LDG

    // Neural EP:
    NekDouble m_beta; // Relative Extracellular resistance: 1 < \beta < 10

    /// Constructor
    MMFDiffusion(const LibUtilities::SessionReaderSharedPtr &pSession,
                 const SpatialDomains::MeshGraphSharedPtr &pGraph);

    InitWaveType m_InitWaveType;

    virtual void v_InitObject(bool DeclareField = true) override;
    
    void TestHelmholtzSolver();

    void DoOdeProjection(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time);

    /// Solve for the diffusion term.
    void DoImplicitSolve(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray, NekDouble time,
        NekDouble lambda);

    /// Computes the reaction terms \f$f(u,v)\f$ and \f$g(u,v)\f$.

    void DoOdeRhs(const Array<OneD, const Array<OneD, NekDouble>> &inarray,
                  Array<OneD, Array<OneD, NekDouble>> &outarray,
                  const NekDouble time);
                  
    void TestPhysDirectionalDeriv(const Array<OneD, const Array<OneD, NekDouble>> &movingframes);

    void TestLineProblem(const int direction, const NekDouble time,
                                    Array<OneD, NekDouble> &outfield);
                                    
    void TestPlaneProblem(const NekDouble time,
                          StdRegions::VarCoeffMap &varcoeff,
                          Array<OneD, NekDouble> &outfield);

    void TestPlaneEmbedProblem(const NekDouble time,
                                    StdRegions::VarCoeffMap &varcoeff,
                                    Array<OneD, NekDouble> &outfield);                   

    void TestPlaneNeumannProblem(const NekDouble time,
                          Array<OneD, NekDouble> &outfield);

    void TestPlaneAniProblem(const NekDouble time,
                        StdRegions::VarCoeffMap &varcoeff,
                        Array<OneD, NekDouble> &outfield);

   void TestPlaneAniDerivProblem(const NekDouble time,
                                    StdRegions::VarCoeffMap &varcoeff,
                                    Array<OneD, NekDouble> &Dxoutfield,
                                    Array<OneD, NekDouble> &Dyoutfield);                     
                                    
                                    
    void TestHelmholtzProblem(const int type,
                            StdRegions::VarCoeffMap &varcoeff,
                            Array<OneD, NekDouble> &outfield);

    void TestCubeProblem(const NekDouble time,
                         Array<OneD, NekDouble> &outfield);

    void Morphogenesis(const NekDouble time, unsigned int field,
                       Array<OneD, NekDouble> &outfield);

    void ComputeEuclideanDivMF(
        const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
        Array<OneD, Array<OneD, NekDouble>> &DivMF);

    Array<OneD, NekDouble> PlanePhiWave();

    void GetFluxVector(
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    const Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &qfield,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &viscousTensor);

    /// Sets a custom initial condition.
    virtual void v_SetInitialConditions(NekDouble initialtime,
                                        bool dumpInitialConditions,
                                        const int domain) override;

    /// Prints a summary of the model parameters.
    virtual void v_GenerateSummary(SolverUtils::SummaryList &s) override;

    virtual void v_EvaluateExactSolution(unsigned int field,
                                         Array<OneD, NekDouble> &outfield,
                                         const NekDouble time) override;

    NekDouble m_InitPtx, m_InitPty, m_InitPtz;

private:
    /// Variable diffusivity
    Array<OneD, NekDouble> m_epsvec;
    Array<OneD, NekDouble> m_epsilon;
    Array<OneD, NekDouble> m_epsu;
};

} // namespace Nektar

#endif
