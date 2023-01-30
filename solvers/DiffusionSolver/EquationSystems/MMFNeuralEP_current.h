///////////////////////////////////////////////////////////////////////////////
//
// File: MMFNeuralEP.h
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
// Description: MMFNeuralEP
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_MMFNeuralEP_H
#define NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_MMFNeuralEP_H

#include <SolverUtils/MMFSystem.h>
#include <SolverUtils/UnsteadySystem.h>

#include <CardiacEPSolver/CellModels/CellModel.h>
#include <CardiacEPSolver/Stimuli/Stimulus.h>

#include <DiffusionSolver/NeuronModels/NeuronModel.h>

// using namespace Nektar::SolverUtils;

namespace Nektar
{

enum TestType
{
    eTestLine,
    eTestLineY,
    eTestPlane,
    eTestCube,
    eTestLinearSphere,
    eTestNonlinearSphere,
    SIZE_TestType ///< Length of enum list
};

const char *const TestTypeMap[] = {
    "TestLine", "TestLineY",        "TestPlane",
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

enum NeuralEPType
{
    eNeuralTest,
    eNeuralEP1D,
    eNeuralEP2p1D,
    eNeuralEP2D,
    eNeuralEP2DEmbed,
    SIZE_NeuralEPType ///< Length of enum list
};

const char *const NeuralEPTypeMap[] = {
    "NeuralTest",
    "NeuralEP1D",
    "NeuralEP1D",
    "NeuralEP2p1D",
    "NeuralEP2D",
    "NeuralEP2DEmbed",
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
    eDeActivated,
    eActivated,
    eProcessing,
    SIZE_TimeMapType ///< Length of enum list
};

const char *const TimeMapTypeMap[] = {
    "DeActivated",
    "Activated",
    "Processing",
};


/// A model for cardiac conduction.
class MMFNeuralEP : public SolverUtils::MMFSystem
{
public:
    friend class MemoryManager<MMFNeuralEP>;

    /// Creates an instance of this class
    static SolverUtils::EquationSystemSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const SpatialDomains::MeshGraphSharedPtr &pGraph)
    {
        SolverUtils::EquationSystemSharedPtr p =
            MemoryManager<MMFNeuralEP>::AllocateSharedPtr(pSession, pGraph);
        p->InitObject();
        return p;
    }

    /// Name of class
    static std::string className;

    TestType m_TestType;
    NeuralEPType m_NeuralEPType;
    SolverSchemeType m_SolverSchemeType;

    StdRegions::VarCoeffMap m_varDiffcoeff;

    Array<OneD, Array<OneD, int>> m_NodeZone;
    Array<OneD, Array<OneD, NekDouble>> m_NeuralCm;
    Array<OneD, Array<OneD, NekDouble>> m_phieNeuralCm;

    Array<OneD, NekDouble> m_NeuralCmRf;
    Array<OneD, int> m_NodeElement;
    /// Desctructor
    virtual ~MMFNeuralEP();

protected:
    int m_nfibers, m_ElemNodeEnd, m_ElemMyelenEnd;
    int m_numelemperNode, m_Convectiven;
    TimeMapType m_TimeMap;

    NeuronModelSharedPtr m_neuron;

    NekDouble m_AnisotropyStrength;

    // Temperature parameter
    NekDouble m_Temperature;
    NekDouble m_TimeMapStart;
    NekDouble m_TimeMapEnd;

    NekDouble m_Helmtau;

    SpatialDomains::GeomMMF m_phieMMFdir;

    // NeuralEP1D: m_beta_e = r_e/r_i
    NekDouble m_ratio_re_ri;

    // NeuralEP: Capacitance vectors for myeline or Ranvier node.
    int m_Rnodelength, m_Rnodegap;

    NekDouble m_Diffbeta, m_Diffeta, m_Diffhe; // h_e for LDG

    // Neural EP:
    NekDouble m_beta; // Relative Extracellular resistance: 1 < \beta < 10

    /// Constructor
    MMFNeuralEP(const LibUtilities::SessionReaderSharedPtr &pSession,
                 const SpatialDomains::MeshGraphSharedPtr &pGraph);

    InitWaveType m_InitWaveType;

    virtual void v_InitObject(bool DeclareField = true) override;

    /// Solve for the diffusion term.
    void DoImplicitSolve(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray, NekDouble time,
        NekDouble lambda);

    /// Computes the reaction terms \f$f(u,v)\f$ and \f$g(u,v)\f$.
    void DoOdeRhs(const Array<OneD, const Array<OneD, NekDouble>> &inarray,
                  Array<OneD, Array<OneD, NekDouble>> &outarray,
                  const NekDouble time);

    void TestPlaneProblem(const NekDouble time,
                          Array<OneD, NekDouble> &outfield);

    void TestCubeProblem(const NekDouble time,
                         Array<OneD, NekDouble> &outfield);

    void Morphogenesis(const NekDouble time, unsigned int field,
                       Array<OneD, NekDouble> &outfield);

    void ComputeVarCoeff2D(
        const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
        StdRegions::VarCoeffMap &varcoeff);

    void ComputeEuclideanDivMF(
        const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
        Array<OneD, Array<OneD, NekDouble>> &DivMF);

    Array<OneD, NekDouble> PlanePhiWave();

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

    Array<OneD, NekDouble> m_epsilon;
    Array<OneD, NekDouble> m_epsu;
};

} // namespace Nektar

#endif
