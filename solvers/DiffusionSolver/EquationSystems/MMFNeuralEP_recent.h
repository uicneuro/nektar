///////////////////////////////////////////////////////////////////////////////
//
// File MMFNeuralEP.h
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
// Description: MMFNeuralEP
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_MMFNEURALEP_H
#define NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_MMFNEURALEP_H

#include <SolverUtils/MMFSystem.h>
#include <SolverUtils/UnsteadySystem.h>
#include <SolverUtils/Diffusion/Diffusion.h>

#include <CardiacEPSolver/CellModels/CellModel.h>
#include <DiffusionSolver/NeuronModels/NeuralStimuli/NeuralStimulus.h>

#include <DiffusionSolver/NeuronModels/NeuronModel.h>

using namespace Nektar::SolverUtils;


namespace Nektar
{

enum NeuralEPType
{
    eNeuralHelmSolveSingle,
    eNeuralHelmSolveDuo,
    eNeuralEPPT,
    eNeuralEP1D,
    eNeuralEP2Dmono,
    eNeuralEP2Dbi,
    SIZE_NeuralEPType ///< Length of enum list
};

const char *const NeuralEPTypeMap[] = {
    "NeuralHelmSolveSingle",
    "NeuralHelmSolveDuo",
    "NeuralEPPT",
    "NeuralEP1D",
    "NeuralEP2Dmono",
    "NeuralEP2Dbi",
};

enum SolverSchemeType
{
    eMMFZero,
    eMMFFirst,
    ePointWise,
    eTimeMap,
    SIZE_SolverSchemeType,
};

const char *const SolverSchemeTypeMap[] = {
    "MMFZero",
    "MMFFirst",
    "PointWise",
    "TimeMap",
};

enum MediumType
{
    eIsotropy,
    eAnisotropy,
    eHeterogeneousIsotropy,
    eHeterogeneousAnisotropy,
    eRegionalHeterogeneous,
    eAllNode,
    SIZE_MediumType
};

const char *const MediumTypeMap[] = {
    "Isotropy",
    "Anisotropy",
    "HeterogeneousIsotropy",
    "HeterogeneousAnisotropy",
    "RegionalHeterogeneous",
    "AllNode",
};

enum FiberType
{
    eSingleLinear,
    eDoubleLinear,
    SIZE_FiberType
};

const char *const FiberTypeMap[] = {
    "SingleLinear",
    "DoubleLinear",
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
    SIZE_TimeMapType ///< Length of enum list
};

const char *const TimeMapTypeMap[] = {
    "Activated",
    "DeActivated",
};

enum ExtCurrentType
{
    eEphapticNode,
    eEphapticIntra,
    eIsolated,
    SIZE_ExtCurrentType ///< Length of enum list
};

const char *const ExtCurrentTypeMap[] = {
    "EphapticNode",
    "EphapticIntra",
    "Isolated",
};

enum NodeIndexType
{
    eSequential,
    eNodefirst,
    SIZE_NodeIndexType ///< Length of enum list
};

const char *const NodeIndexTypeMap[] = {
    "Sequential",
    "Nodefirst",
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

    NeuralEPType m_NeuralEPType;
    SolverSchemeType m_SolverSchemeType;

    ExtCurrentType m_ExtCurrentType;

    NodeIndexType m_NodeIndexType;

    NekDouble d_max, d_min;

    virtual void v_InitObject(bool DeclareField = true) override;
    virtual void v_DoSolve() override;

    /// Desctructor
    virtual ~MMFNeuralEP();
protected:

    SolverUtils::DiffusionSharedPtr m_diffusion;
    SolverUtils::RiemannSolverSharedPtr m_riemannSolver;

    MediumType m_MediumType;
    FiberType m_fiberType;

    StdRegions::VarCoeffMap m_varcoeff;
    StdRegions::VarCoeffMap m_phievarcoeff;

    int m_npts; // Number of points for each element
    int m_nfibers, m_ElemNodeEnd, m_ElemMyelenEnd, m_ElemExtEnd;
    int m_Convectiven;
    int m_totNode, m_elemperNode, m_elemperMyel;
    int m_zonestart, m_zoneend;

    NekDouble m_gratio, m_relfiberratio, m_radiusfiberbundle, m_radiusaxon;

    NekDouble m_fiber1left, m_fiber1right, m_fiber2left, m_fiber2right, m_nodeinitdown, m_nodeinitup;

    NekDouble m_fiberwidth, m_fibergap, m_fiberheightdiff;
    NekDouble m_nodelen, m_myelinlen;
    NekDouble m_InitPtx, m_InitPty, m_InitPtz;
    NekDouble m_Rf, m_Cn, m_Cm;
    NekDouble m_phimrest, m_phimTol, m_dphimdtTol;

    std::string m_zoneindexfile;

    Array<OneD, NekDouble> m_StimAtNode;

    TimeMapType m_TimeMapScheme;
    
    Array<OneD, NekDouble> ComputeConductivity(
                 const Array<OneD, const int> &zoneindex);

    void Generatephiemovingframes(
    const NekDouble ratio_re_ri,
    Array<OneD, Array<OneD, NekDouble>> &helmfmovingframes, 
    Array<OneD, Array<OneD, NekDouble>> &phiemovingframes);

    // other moving frames neede for NeuralEP
    Array<OneD, Array<OneD, NekDouble>> m_unitmovingframes;
    Array<OneD, Array<OneD, NekDouble>> m_phiemovingframes;

    // Elements for fiber 2D: Start and End index
    int m_fiber2DElemStart, m_fiber2DElemEnd;

    // Temperature parameter
    NekDouble m_Temperature;
    NekDouble m_diameter;

    NekDouble m_TimeMapStart;
    NekDouble m_TimeMapEnd;

    SpatialDomains::GeomMMF m_phieMMFdir;

    // NeuralEP1D: m_beta_e = r_e/r_i
    NekDouble m_ratio_re_ri;

    NekDouble m_Diffbeta, m_Diffeta, m_Diffhe; // h_e for LDG

    // Scar tisseu related variables
    NekDouble m_PVcond;
    NekDouble m_ScarSize, m_ScarStr, m_ScarPis, m_ScarLocx, m_ScarLocy, m_ScarLocz;
    NekDouble m_RelDivSize, m_RelDivStr, m_RelDivPis, m_RelDivLocx;

    // Relative divergence related variables
    NekDouble m_LambDivSmoothL;

    // Neural EP:
    NekDouble m_beta; // Relative Extracellular resistance: 1 < \beta < 10

    // NeuralEP: Capacitance vectors for myeline or Ranvier node.
    StdRegions::VarCoeffMap m_varDiffcoeff;

    Array<OneD, int> m_ValidTimeMap;

    Array<OneD, Array<OneD, int>> m_zoneindex;

    Array<OneD, NekDouble> m_excitezone;
    Array<OneD, NekDouble> m_excitezone1;
    Array<OneD, NekDouble> m_excitezone2;

    Array<OneD, NekDouble> m_nodezone;
    Array<OneD, NekDouble> m_nodezone1;
    Array<OneD, NekDouble> m_nodezone2;

    Array<OneD, NekDouble> m_intrazone;
    Array<OneD, NekDouble> m_intrazone1;
    Array<OneD, NekDouble> m_intrazone2;

    Array<OneD, NekDouble> m_myelinzone1;
    Array<OneD, NekDouble> m_myelinzone2;

    Array<OneD, NekDouble> m_extrazone;

    Array<OneD, Array<OneD, NekDouble>> m_AniStrength;
    Array<OneD, Array<OneD, NekDouble>> m_phieAniStrength;

    Array<OneD, Array<OneD, NekDouble>> m_NeuralCm;
    Array<OneD, Array<OneD, NekDouble>> m_phieNeuralCm;

    Array<OneD, Array<OneD, NekDouble>> m_TimeMap;

    Array<OneD, int> m_InternalBoundary;

    Array<OneD, NekDouble> m_NeuralCmRf;
    Array<OneD, int> m_NodeElement;

    // NeuralEP2D variables
    // Array<OneD, LibUtilities::TimeIntegrationWrapperSharedPtr>
    // m_fiberintScheme; Array<OneD,
    // LibUtilities::TimeIntegrationSchemeOperators> m_fiberode; Array<OneD,
    // LibUtilities::TimeIntegrationSolutionSharedPtr> m_fiberintSoln;

    // Moving frames
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> m_fibermovingframes;

    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> m_ncdotfiberMFFwd;
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> m_ncdotfiberMFBwd;
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> m_DivfiberMF;

    Array<OneD, MultiRegions::ExpListSharedPtr> m_fiberfields;
    Array<OneD, Array<OneD, Array<OneD, int>>> m_fiberindex;

    Array<OneD, StdRegions::VarCoeffMap> m_fibervarcoeff;

    /// Constructor
    MMFNeuralEP(const LibUtilities::SessionReaderSharedPtr &pSession,
                const SpatialDomains::MeshGraphSharedPtr &pGraph);

    LibUtilities::SessionReaderSharedPtr m_session1D;
    SpatialDomains::MeshGraphSharedPtr m_graph1D;

    Array<OneD, LibUtilities::SessionReaderSharedPtr> m_fibersession;
    Array<OneD, SpatialDomains::MeshGraphSharedPtr> m_fibergraph;

    Array<OneD, NekDouble> ExtractFiberValue(
        const int nfib, const Array<OneD, const NekDouble> &inarray);

    void UpdateFibertoField(const int nfib,
                            const Array<OneD, const NekDouble> &fiberinarray,
                            Array<OneD, NekDouble> &outarray);

    InitWaveType m_InitWaveType;

    // void DoSolveMMFFirst();
    void DoSolveMMF();
    void DoSolvePoint();
    
    // Coefficients for Anisotropy
    int m_AnisotropyRegion;
    int m_InnerboxEnd;
    NekDouble m_AnisotropyStrength;
    
    void CheckOutZoneAni();

    Array<OneD, int> GetInternalBoundaryPoints();

    void PlotAnisotropyFiber(const Array<OneD, const NekDouble> &anifibre);
    
    void PlotAnisotropy(
    const Array<OneD, const Array<OneD, NekDouble>> &AniStrength, 
    const Array<OneD, const Array<OneD, NekDouble>> &phieAniStrength);

    void PlotPhieMF(
    const Array<OneD, const Array<OneD, NekDouble>> &sigma_i,
    const Array<OneD, const Array<OneD, NekDouble>> &sigma_e,
    const Array<OneD, const Array<OneD, NekDouble>> &PhieAniStrength);

    // void DisplayNode1D(std::string &fulltext, const Array<OneD, const Array<OneD, NekDouble>> &fields);
    void DisplayNode2D(std::string &fulltext, const Array<OneD, const Array<OneD, NekDouble>> &fields);

    void CheckNodeZoneMF(
    const Array<OneD, const Array<OneD, int>> &NodeZone,
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    const Array<OneD, const Array<OneD, NekDouble>> &phiemovingframes);

    void GetFluxVector(
        const Array<OneD, Array<OneD, NekDouble>> &inarray,
        const Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &qfield,
        Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &viscousTensor);

    // void ImportFiberXml(
    //     const int nfibers,
    //     Array<OneD, MultiRegions::ExpListSharedPtr> &fiberfields,
    //     Array<OneD, Array<OneD, Array<OneD, int>>> &fiberindex,
    //     Array<OneD, LibUtilities::SessionReaderSharedPtr> &fibersession,
    //     Array<OneD, SpatialDomains::MeshGraphSharedPtr> &fibergraph);

    /// Solve for the diffusion term.
    void DoImplicitSolve(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time,
        const NekDouble lambda);

    void DoImplicitSolveNeuralEP1D(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time,
        const NekDouble lambda);

    void DoImplicitSolveNeuralEP2Dmono(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time,
        const NekDouble lambda);

    void DoImplicitSolveNeuralEP2Dbi(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time,
        const NekDouble lambda);

    void DoImplicitSolveNeuralEP2DEmbbi(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time,
        const NekDouble lambda);

    // void DoImplicitSolveNeuralEP2p1Dfiber(
    //     const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    //     Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time,
    //     const NekDouble lambda);

    void DoImplicitSolveNeuralEP2Dfiber(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time,
        const NekDouble lambda);

    void DoNullSolve(const Array<OneD, const Array<OneD, NekDouble>> &inarray,
                     Array<OneD, Array<OneD, NekDouble>> &outarray,
                     const NekDouble time, const NekDouble lambda);

    /// Computes the reaction terms \f$f(u,v)\f$ and \f$g(u,v)\f$.
    void DoOdeRhs(const Array<OneD, const Array<OneD, NekDouble>> &inarray,
                  Array<OneD, Array<OneD, NekDouble>> &outarray,
                  const NekDouble time);

    void DoOdeRhsNeuralEPPT(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time);

    void DoOdeRhsNeuralEP1D(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time);

    void DoOdeRhsNeuralEP2p1D(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time);

    void DoOdeRhsNeuralEP2Dmono(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time);

    void DoOdeRhsNeuralEP2Dbi(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time);

    void DoOdeRhsNeuralEP2DEmbbi(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time);

    Array<OneD, NekDouble> Computephie(
        ExtCurrentType extcurrent,
        const Array<OneD, const NekDouble> &phim);

    void MembraneBoundary2D(int bcRegion, int cnt,
                            Array<OneD, Array<OneD, NekDouble>> &Fwd,
                            Array<OneD, Array<OneD, NekDouble>> &physarray);

    void DoOdeRhsNeuralEP2p1Dfiber(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        const Array<OneD, const Array<OneD, NekDouble>> &MF1st,
        Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time);

    void DoOdeRhsNeuralEP2Dfiber(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        const Array<OneD, const Array<OneD, NekDouble>> &MF1st,
        Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time);

    void DoOdeProjection(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time);

    // void SetMembraneBoundaryCondition(const NekDouble time = 0.0);

    void PlotFHIonCurrent(const Array<OneD, const NekDouble> &inarray,
                          const int nstep);

    void SetAxonWallBoundaryConditions(
        Array<OneD, Array<OneD, NekDouble>> &inarray);

    void AxonWallBoundary2D(int bcRegion, int cnt,
                            Array<OneD, Array<OneD, NekDouble>> &Fwd,
                            Array<OneD, Array<OneD, NekDouble>> &physarray);

    void TestPlaneProblem(const NekDouble time,
                          Array<OneD, NekDouble> &outfield);

    void TestCubeProblem(const NekDouble time,
                         Array<OneD, NekDouble> &outfield);

    void Morphogenesis(const NekDouble time, unsigned int field,
                       Array<OneD, NekDouble> &outfield);

    void ComputeEuclideanDivMF(
        const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
        Array<OneD, Array<OneD, NekDouble>> &DivMF);
        
        
    void PrintRegionalAvgMax(const Array<OneD, const NekDouble> &field0);

    Array<OneD, NekDouble> PlanePhiWave();
    
    Array<OneD, int> TestRanvierSingleIndex();
    Array<OneD, int> TestRanvierDuoIndex();

    Array<OneD, int> IndexNodeZone1D(
        const MultiRegions::ExpListSharedPtr &field, const int Nnode, 
        const int NumelemNode, const int NumelemMyel);

    Array<OneD, int> IndexNodeZone2D(
        const FiberType fiber);

    Array<OneD, int> SingleLinearIndex(
        const NekDouble fiberwidth, 
        const NekDouble nodelen, 
        const NekDouble myelinlen,
        const int Nnode);

int DoubleLinearIndex(
    const int totNnode,
    const NekDouble nodelen, const NekDouble myelinlen,
    const NekDouble fiber1left, const NekDouble fiber1right, 
    const NekDouble fiber2left, const NekDouble fiber2right,
    const NekDouble nodeinitdown, const NekDouble nodeinitup,
    const NekDouble xi, const NekDouble yi);

    // int DoubleLinearIndex(
    //     const NekDouble nodelen, 
    //     const NekDouble myelinlen,
    //     const int Nnode,
    //     const NekDouble xi, 
    //     const NekDouble yi);

    void Getcellavg(
        Array<OneD, NekDouble> &xcell, 
        Array<OneD, NekDouble> &ycell, 
        Array<OneD, NekDouble> &zcell);

    void SetUpBiAnisotropy(
            const Array<OneD, const int> &zoneindex,
            const Array<OneD, const Array<OneD, NekDouble>> NeuralCm,
            Array<OneD, Array<OneD, NekDouble>> &AniStrength);

    void SetUpDomainZone(
        const Array<OneD, const int> &zoneindex,
        Array<OneD, NekDouble> &excitezone,
        Array<OneD, NekDouble> &nodezone,
        Array<OneD, NekDouble> &intrazone,
        Array<OneD, NekDouble> &extrazone);

    void SetUpDomainDuoZone(
        const Array<OneD, const int> &zoneindex,
        Array<OneD, NekDouble> &excitezone1,
        Array<OneD, NekDouble> &excitezone2,
        Array<OneD, NekDouble> &nodezone1,
        Array<OneD, NekDouble> &nodezone2,
        Array<OneD, NekDouble> &intrazone1,
        Array<OneD, NekDouble> &intrazone2,
        Array<OneD, NekDouble> &extrazone);


    void ComputeNeuralTimeMap(const NekDouble time,
                            const Array<OneD, const int> &zoneindex,
                            const Array<OneD, const NekDouble> &field,
                            const Array<OneD, const NekDouble> &dudt,
                            Array<OneD, NekDouble> &dudtHistory,
                            Array<OneD, NekDouble> &TimeMap);
                                
    void PlotNeuralTimeMap(
        const Array<OneD, const NekDouble> &phim,
        const Array<OneD, const NekDouble> &TimeMap,
        const int nstep);
        
    void PlotZone(const Array<OneD, const int> &zoneindex,
Array<OneD, NekDouble> &excitezone1, 
Array<OneD, NekDouble> &excitezone2, 
Array<OneD, NekDouble> &nodezone1, 
Array<OneD, NekDouble> &nodezone2, 
Array<OneD, NekDouble> &intrazone1, 
Array<OneD, NekDouble> &intrazone2, 
Array<OneD, NekDouble> &extrazone);

    void PrintSingleCurrent(const Array<OneD, const NekDouble> &phim);
    void PrintDuoCurrent(const Array<OneD, const NekDouble> &phim,
                        const Array<OneD, const NekDouble> &dudt,
                        NekDouble &thredlocf1, NekDouble &thredlocf2);

    Array<OneD, NekDouble> ConvertTMtoVel(
        const Array<OneD, const NekDouble> &TimeMap,
        const Array<OneD, const NekDouble> &TmapGrad,
        const Array<OneD, const NekDouble> &TmapGradMag);

    /// Sets a custom initial condition.
    virtual void v_SetInitialConditions(NekDouble initialtime,
                                        bool dumpInitialConditions,
                                        const int domain) override;

    /// Prints a summary of the model parameters.
    virtual void v_GenerateSummary(SolverUtils::SummaryList &s) override;

    virtual void v_EvaluateExactSolution(unsigned int field,
                                         Array<OneD, NekDouble> &outfield,
                                         const NekDouble time) override;

private:
    /// Variable diffusivity
    NekDouble m_chi;
    NekDouble m_capMembrane;
    NekDouble m_conductivity;

    CellModelSharedPtr m_cell;

    NeuronModelSharedPtr m_neuron;
    Array<OneD, NeuronModelSharedPtr> m_fiberneurons;

    std::vector<NeuralStimulusSharedPtr> m_stimulus;
    std::vector<NeuralStimulusSharedPtr> m_fiberstimulus;

    Array<OneD, NekDouble> ComputeLaplacianDiff(
        const Array<OneD, const NekDouble> &Laplacian,
        const Array<OneD, const NekDouble> &LaplacianNew);

    /// Stimulus current
    NekDouble m_stimDuration;

    void LoadStimuli();

    int CountActivation(Array<OneD, int> &Activated);
};

} // namespace Nektar

#endif
