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
#include <MMFSolver/EquationSystems/NeuronModels/NeuralStimuli/NeuralStimulus.h>

#include <MMFSolver/EquationSystems/NeuronModels/NeuronModel.h>

using namespace Nektar::SolverUtils;


namespace Nektar
{

enum NeuralEPType
{
    eNeuralHelmSolveSingle,
    eNeuralHelmSolveDuo,
    eNeuralEP2Dmono,
    eNeuralEP2Dbi,
    eNeuralEP2DbiMulti,
    eNeuralEP2DbiMultiv2,
    eNeuralEP2DbiCSD,
    SIZE_NeuralEPType ///< Length of enum list
};

const char *const NeuralEPTypeMap[] = {
    "NeuralHelmSolveSingle",
    "NeuralHelmSolveDuo",
    "NeuralEP2Dmono",
    "NeuralEP2Dbi",
    "NeuralEP2DbiMulti",
    "NeuralEP2DbiMultiv2",
    "NeuralEP2DbiCSD",
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

enum ExtCurrentType
{
    eEphaptic,
    eNoEphaptic,
    SIZE_ExtCurrentType ///< Length of enum list
};

const char *const ExtCurrentTypeMap[] = {
    "Ephaptic",
    "NoEphaptic",
};

enum FiberType
{
    eLinearAligned,
    eLinearMisAligned,
    eLinearDivergent,
    eLinearCrossing,
    eConstantCurved,
    SIZE_FiberType ///< Length of enum list
};

const char *const FiberTypeMap[] = {
    "LinearAligned",
    "LinearMisAligned",
    "LinearDivergent",
    "LinearCrossing",
    "ConstantCurved",
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

    FiberType m_FiberType;
    ExtCurrentType m_ExtCurrentType;

    NekDouble d_max, d_min;

    virtual void v_InitObject(bool DeclareField = true) override;
    virtual void v_DoSolve() override;

    /// Desctructor
    virtual ~MMFNeuralEP();
protected:

    SolverUtils::DiffusionSharedPtr m_diffusion;
    SolverUtils::RiemannSolverSharedPtr m_riemannSolver;

    MediumType m_MediumType;

    StdRegions::VarCoeffMap m_varcoeff;
    StdRegions::VarCoeffMap m_phievarcoeff;
    StdRegions::VarCoeffMap m_CSDvarcoeff;

    Array<OneD, StdRegions::VarCoeffMap> m_varcoefffiber;
    Array<OneD, StdRegions::VarCoeffMap> m_phievarcoefffiber;

    int m_phimvar, m_phievar;
    int m_npts, m_nfibers, m_ElemNodeEnd, m_ElemMyelenEnd, m_ElemExtEnd, m_Convectiven;
    int m_numfiber, m_totNode, m_zonestart, m_zoneend, m_myeline, m_node, m_external;
    int m_fiber2DElemStart, m_fiber2DElemEnd, m_fiber1order, m_fiber2order, m_fiber3order;
    int m_AnisotropyRegion, m_InnerboxEnd;

    NekDouble m_pi;
    
    NekDouble m_gratio, m_relfiberratio, m_radiusfiberbundle, m_fibercurvature, m_phiefactor;
    NekDouble m_axondiameter, m_fiberangle, m_fiberlength, m_nodeinitdown, m_nodeinitup, m_fiberwidth, m_fibergap;
    NekDouble m_fiber1left, m_fiber1right, m_fiber2left, m_fiber2right, m_fiber3left, m_fiber3right, m_bundleleft, m_bundleright;

    NekDouble m_nodelen, m_myelinlen, m_InitPtx, m_InitPty, m_InitPtz;
    NekDouble m_Rf, m_Cn, m_Cm, m_phimrest, m_Diffext, m_AnisotropyStrength;

    // NeuralEP1D: m_beta_e = r_e/r_i
    NekDouble m_ratio_re_ri, m_Diffbeta, m_Diffeta, m_Diffhe; // h_e for LDG

    // Scar tisseu related variables
    NekDouble m_PVcond, m_ScarSize, m_ScarStr, m_ScarPis, m_ScarLocx, m_ScarLocy, m_ScarLocz;
    NekDouble m_RelDivSize, m_RelDivStr, m_RelDivPis, m_RelDivLocx;

    // Temperature parameter
    NekDouble m_Temperature, m_TimeMapStart, m_TimeMapEnd;

    Array<OneD, NekDouble> m_x, m_y, m_z;
    Array<OneD, NekDouble> m_xcell, m_ycell, m_zcell;
    
    Array<OneD, int> m_fiberorder;

    Array<OneD, NekDouble> m_fiberleft;
    Array<OneD, NekDouble> m_fiberright;
    
    Array<OneD, NekDouble> ComputeConductivity(
                 const Array<OneD, const int> &zoneindex);

    Array<OneD, NekDouble> ComputeConductivity(
                 const Array<OneD, const Array<OneD, int>> &zoneindex);

    // other moving frames neede for NeuralEP
    NekDouble m_CSDDiff;
    Array<OneD, Array<OneD, NekDouble>> m_CSDmovingframes;
    Array<OneD, Array<OneD, NekDouble>> m_phiemovingframes;

    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> m_phiemovingframesfiber;

    SpatialDomains::GeomMMF m_phieMMFdir;

    // NeuralEP: Capacitance vectors for myeline or Ranvier node.
    Array<OneD, int> m_ValidTimeMap;
    Array<OneD, int> m_zoneindex;
    Array<OneD, Array<OneD, int>> m_zoneindexfiber;

    Array<OneD, Array<OneD, NekDouble>> m_excitezonefiber;
    Array<OneD, Array<OneD, NekDouble>> m_intrazonefiber;
    Array<OneD, Array<OneD, NekDouble>> m_nodezonefiber;

    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> m_movingframesfiber;
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> m_AniStrengthfiber;
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> m_phieAniStrengthfiber;
    Array<OneD, Array<OneD, NekDouble>> m_phiediffAniStrength;

    Array<OneD, Array<OneD, NekDouble>> m_AniStrengthglobal;
    Array<OneD, Array<OneD, NekDouble>> m_phieAniStrengthglobal;

    Array<OneD, NekDouble> m_nodezone;
    Array<OneD, NekDouble> m_myelinzone;
    Array<OneD, NekDouble> m_intrazone;
    Array<OneD, NekDouble> m_extrazone;
    Array<OneD, NekDouble> m_outerzone;

    Array<OneD, Array<OneD, NekDouble>> m_AniStrength;
    Array<OneD, Array<OneD, NekDouble>> m_phieAniStrength;
    Array<OneD, Array<OneD, NekDouble>> m_phieAniStrengthv2;

    void ComputeRegionalSigma(
        const Array<OneD, const int> &zoneindex,
        Array<OneD, Array<OneD, NekDouble>> &sigma_i,
        Array<OneD, Array<OneD, NekDouble>> &sigma_e,
        Array<OneD, Array<OneD, NekDouble>> &sigma_eM);

    // Array<OneD, NekDouble> m_NeuralCm;
    Array<OneD, Array<OneD, NekDouble>> m_NeuralCm;
    Array<OneD, NekDouble> m_NeuralCmfiber;

    Array<OneD, Array<OneD, NekDouble>> m_TimeMap;
    Array<OneD, Array<OneD, NekDouble>> m_PhieCurrent;

    Array<OneD, int> m_InternalBoundary;
    Array<OneD, int> m_NodeElement;

    /// Constructor
    MMFNeuralEP(const LibUtilities::SessionReaderSharedPtr &pSession,
                const SpatialDomains::MeshGraphSharedPtr &pGraph);

    Array<OneD, int> GetInternalBoundaryPoints();
    
    void SetUpNeuralCm();

    void ComputeAniStrengthfiber(
        const Array<OneD, const Array<OneD, int>> &zoneindexfiber, 
        Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &AniStrengthfiber);
    
    Array<OneD, NekDouble> ComputeNeuralCmfiber(
        const Array<OneD, const int> &zoneindex);
    
    void ComputeNeuralCmfiber(
            const Array<OneD, const Array<OneD, int>> &zoneindexfiber,
            Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &outarrayfiber);

    void ComputeGlobalAniStrength(
    const Array<OneD, const Array<OneD, Array<OneD, NekDouble>>>
        &AniStrengthfiber,
    Array<OneD, Array<OneD, NekDouble>> &outarray);

    void ComputeGlobalPhieAniStrength(
        const Array<OneD, const int> &zoneindex,
        Array<OneD, Array<OneD, NekDouble>> &phieAniStrength);

    void ComputeGlobalPhieAniStrengthfiber(
        const Array<OneD, const int> &zoneindex,
        const Array<OneD, const Array<OneD, int>> &zoneindexfiber,
        Array<OneD, Array<OneD, NekDouble>> &phieAniStrength);

    void PlotAnisotropyfiber(
        const Array<OneD, const Array<OneD, NekDouble>> &AniStrength, 
        const Array<OneD, const Array<OneD, Array<OneD, NekDouble>>> &AniStrengthfiber, 
        const Array<OneD, const Array<OneD, NekDouble>> &phieAniStrength, 
        const Array<OneD, const Array<OneD, Array<OneD, NekDouble>>> &phieAniStrengthfiber);
        
    void PlotAnisotropy(
    const Array<OneD, const Array<OneD, NekDouble>> &AniStrength, 
    const Array<OneD, const Array<OneD, NekDouble>> &phieAniStrength);

    void PlotPhieMF(
    const Array<OneD, const Array<OneD, NekDouble>> &sigma_i,
    const Array<OneD, const Array<OneD, NekDouble>> &sigma_e,
    const Array<OneD, const Array<OneD, NekDouble>> &PhieAniStrength);

    NekDouble DisplayAtNodes(const int fibern, const int nodeindex, 
                                const Array<OneD, const Array<OneD, int>> &zoneindexfiber,
                                const Array<OneD, const NekDouble> &inarray);

    void CheckNodeZoneMF(
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    const Array<OneD, const Array<OneD, NekDouble>> &phiemovingframes);

    void GetFluxVector(
        const Array<OneD, Array<OneD, NekDouble>> &inarray,
        const Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &qfield,
        Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &viscousTensor);

    /// Solve for the diffusion term.
    void DoImplicitSolveNeuralEP2Dmono(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time,
        const NekDouble lambda);

    void DoImplicitSolveNeuralEP2Dbi(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time,
        const NekDouble lambda);

        void DoImplicitSolveNeuralEP2DbiMulti(
            const Array<OneD, const Array<OneD, NekDouble>> &inarray,
            Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time,
            const NekDouble lambda);

            void DoImplicitSolveNeuralEP2DbiMultiv2(
                const Array<OneD, const Array<OneD, NekDouble>> &inarray,
                Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time,
                const NekDouble lambda);

    void DoImplicitSolveNeuralEP2DbiCSD(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time,
        const NekDouble lambda);

    void DoNullSolve(const Array<OneD, const Array<OneD, NekDouble>> &inarray,
                     Array<OneD, Array<OneD, NekDouble>> &outarray,
                     const NekDouble time, const NekDouble lambda);

    /// Computes the reaction terms \f$f(u,v)\f$ and \f$g(u,v)\f$.
    void DoOdeRhsNeuralEP2Dmono(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time);

    void DoOdeRhsNeuralEP2Dbi(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time);

    void DoOdeRhsNeuralEP2DbiMulti(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time);

        void DoOdeRhsNeuralEP2DbiMultiv2(
            const Array<OneD, const Array<OneD, NekDouble>> &inarray,
            Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time);

    void DoOdeRhsNeuralEP2DbiCSD(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time);

    Array<OneD, NekDouble> ComputeFieldPhie(
        const int phievar,
        const Array<OneD, const NekDouble> &phim);

    Array<OneD, NekDouble> ComputeFieldPhiefiber(
            const Array<OneD, const Array<OneD, NekDouble>> &inarray);

            Array<OneD, NekDouble> ComputeFieldPhiefiberv2(
                const int phievar,
                const Array<OneD, const Array<OneD, NekDouble>> &inarray);

    void DoOdeProjection(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time);

    void PlotFHIonCurrent(const Array<OneD, const NekDouble> &inarray,
                          const int nstep);
    
    Array<OneD, int> TestRanvierSingleIndex();
    Array<OneD, int> TestRanvierDuoIndex();

    void IndexNodeZone2D(       
            const Array<OneD, const NekDouble> &fiberleft,
            const Array<OneD, const NekDouble> &fiberright,
            const Array<OneD, const int> &fiberorder,
            Array<OneD, Array<OneD, int>> &zoneindexfiber);

    int FiberIndex(FiberType FiberType, const int fibern,
        const int fiberorder, const NekDouble xi, const NekDouble yi);

    int LinearAlignedFiberIndex(const int fiberorder, const NekDouble yi);

    int LinearMisAlignedFiberIndex(const int fibern, const NekDouble yi);

    int LinearDivergentFiberIndex(const int fibern, const NekDouble xi, const NekDouble yi);

    int LinearCrossingFiberIndex(const int fibern, const int fiberorder, const NekDouble xi, const NekDouble yi);

    int ConstantCurvedFiberIndex(const int fibern,  
        const NekDouble fibercurvature, const NekDouble xi, const NekDouble yi);

    void GetCellCoordAvg(
        Array<OneD, NekDouble> &xcell, 
        Array<OneD, NekDouble> &ycell, 
        Array<OneD, NekDouble> &zcell);

    void SetUpDomainZone();

void PlotDomainZone();

void TestHelmSolve();

void PlotDomainZonefib2(
        const Array<OneD, const Array<OneD, int>> zoneindexfiber,
        const Array<OneD, const int> &zoneindex,
        const Array<OneD, const Array<OneD, NekDouble>> &intrazonefiber,
        const Array<OneD, const NekDouble> &intrazone,
        const Array<OneD, const NekDouble> &extrazone,
        const Array<OneD, const NekDouble> &outerzone);

void PlotDomainZonefib1(
        const Array<OneD, const int> &zoneindex,
        const Array<OneD, const NekDouble> &intrazone,
        const Array<OneD, const NekDouble> &extrazone,
        const Array<OneD, const NekDouble> &outerzone);

void ComputeNeuralTimeMap(
    const NekDouble time, 
    const Array<OneD, const Array<OneD, NekDouble>> &fields,
    const Array<OneD, const Array<OneD, NekDouble>> &dphidts,
    Array<OneD, Array<OneD, NekDouble>> &dphidtints, 
    Array<OneD, Array<OneD, NekDouble>> &TimeMaps);

void ComputephimTimeMap(const NekDouble time,
                        const Array<OneD, const NekDouble> &field,
                        const Array<OneD, const NekDouble> &dphidt,
                        Array<OneD, NekDouble> &dphidtint,
                        Array<OneD, NekDouble> &TimeMap);

void ComputephieTimeMap(
    const NekDouble dt,
    const Array<OneD, const NekDouble> &field,
    const Array<OneD, const NekDouble> &dphidt,
    Array<OneD, NekDouble> &dphidtint,
    Array<OneD, NekDouble> &TimeMap);

    void ComputerhoTimeMap(
        const NekDouble time,
        const Array<OneD, const NekDouble> &field,
        const Array<OneD, const NekDouble> &rhovec,
        Array<OneD, NekDouble> &rhovecint,
        Array<OneD, NekDouble> &TimeMap);

void RescaleMovingFrames(
    const Array<OneD, const Array<OneD, NekDouble>> &AniStrength,
    Array<OneD, Array<OneD, NekDouble>> &movingframes);

void ComputerhoTimeMap(const NekDouble time,
                        const Array<OneD, const NekDouble> &field,
                        Array<OneD, NekDouble> &fieldint,
                        Array<OneD, NekDouble> &TimeMap);                        

void Computemovingframesfiber(
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    const Array<OneD, const Array<OneD, Array<OneD, NekDouble>>> &AniStrengthfiber,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &movingframesfiber);

void Computephiemovingframesfiber(
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes,
    const Array<OneD, const Array<OneD, Array<OneD, NekDouble>>> &AniStrengthfiber,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &movingframesfiber);

void ComputePhieAniStrengthfiber(
    const Array<OneD, const Array<OneD, int>> &zoneindexfiber, 
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &phieAniStrengthfiber);

void ComputePhieCurrent(const NekDouble time,
                        const NekDouble timestep,
                        const Array<OneD, const Array<OneD, NekDouble>> &field,
                        Array<OneD, Array<OneD, NekDouble>> &PhieCurrent);

void PlotNeuralEP(
    const Array<OneD, const Array<OneD, NekDouble>> &fields,
    const Array<OneD, const Array<OneD, NekDouble>> &TimeMap,
    const int nstep);

    void PlotNeuralEPvar2(
    const Array<OneD, const Array<OneD, NekDouble>> &fields,
    const Array<OneD, const Array<OneD, NekDouble>> &TimeMap,
    const int nstep);

    void PlotNeuralEPvar3(
    const Array<OneD, const Array<OneD, NekDouble>> &fields,
    const Array<OneD, const Array<OneD, NekDouble>> &TimeMap,
    const int nstep);

    void PlotNeuralEPvar3CSD(
        const Array<OneD, const Array<OneD, NekDouble>> &fields,
        const Array<OneD, const Array<OneD, NekDouble>> &TimeMap,
        const int nstep);

    void PlotNeuralEPvar4(
    const Array<OneD, const Array<OneD, NekDouble>> &fields,
    const Array<OneD, const Array<OneD, NekDouble>> &TimeMap,
    const int nstep);

void PrintAtNodes(const int nfields, const int numfiber,
                const Array<OneD, const Array<OneD, NekDouble>> &fields);

void PrintSingleCurrent(const Array<OneD, const NekDouble> &phim,
                                  const Array<OneD, const NekDouble> &dudt,
                                  NekDouble &thredlocf1);
                                  
void PrintDuoCurrent(const Array<OneD, const Array<OneD, NekDouble>> &field);

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
};

} // namespace Nektar

#endif
