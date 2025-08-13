///////////////////////////////////////////////////////////////////////////////
//
// File MMFBidomain.h
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
// Description: MMFBidomain
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_MMFBidomain_H
#define NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_MMFBidomain_H

#include <MultiRegions/ContField2D.h>
#include <SolverUtils/MMFSystem.h>
#include <SolverUtils/UnsteadySystem.h>

#include <DiffusionSolver/CellModels/CellModel.h>
#include <DiffusionSolver/Stimuli/Stimulus.h>

// using namespace Nektar::SolverUtils;

namespace Nektar
{

enum TestType
{
    eBidomainCardiac,
    eBidomainBrain,
    SIZE_TestType ///< Length of enum list
};

const char *const TestTypeMap[] = {
    "BidomainCardiac",
    "BidomainBrain",
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

enum MediumType
{
    eIsotropy,
    eAnisotropy,
    SIZE_MediumType
};

const char *const MediumTypeMap[] = {
    "Isotropy",
    "Anisotropy",
};

/// A model for cardiac conduction.
class MMFBidomain : public SolverUtils::MMFSystem
{
public:
    friend class MemoryManager<MMFBidomain>;

    /// Creates an instance of this class
    static SolverUtils::EquationSystemSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const SpatialDomains::MeshGraphSharedPtr &pGraph)
    {
        SolverUtils::EquationSystemSharedPtr p =
            MemoryManager<MMFBidomain>::AllocateSharedPtr(pSession, pGraph);
        p->InitObject();
        return p;
    }

    /// Name of class
    static std::string className;

    TestType m_TestType;

    MediumType m_MediumType;

    NekDouble d_max, d_min;

    /// Desctructor
    virtual ~MMFBidomain();

protected:
    int m_Convectiven;
    NekDouble m_TMT0;

    NekDouble m_Diffbeta; //  \beta_e for LDG
    NekDouble m_Diffeta;  // \eta for LDG
    NekDouble m_Diffhe;   // h_e for LDG

    /// Constructor
    MMFBidomain(const LibUtilities::SessionReaderSharedPtr &pSession,
                const SpatialDomains::MeshGraphSharedPtr &pGraph);

    InitWaveType m_InitWaveType;

    virtual void v_InitObject();

    virtual void v_DoSolve();

    // When the MMF1st alignment is initiated
    NekDouble m_TimeToActivateAlignment;

    // Coefficients for Anisotropy
    int m_AnisotropyRegion;
    NekDouble m_AnisotropyStrength;
    void PlotAnisotropyFiber(const Array<OneD, const NekDouble> &anifibre);

    void DoNullSolveMMF(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        const Array<OneD, const Array<OneD, NekDouble>> &MF1st,
        Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time,
        const NekDouble lambda);

    void DoImplicitSolveBidomain(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        const Array<OneD, const Array<OneD, NekDouble>> &MF1st,
        Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time,
        const NekDouble lambda);

    void DoImplicitSolveBidomainRoth(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        const Array<OneD, const Array<OneD, NekDouble>> &MF1st,
        Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time,
        const NekDouble lambda);

    void DoOdeRhsBidomain(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        const Array<OneD, const Array<OneD, NekDouble>> &MF1st,
        Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time);

    void DoOdeProjection(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time);

    /// Sets a custom initial condition.
    virtual void v_SetInitialConditions(NekDouble initialtime,
                                        bool dumpInitialConditions,
                                        const int domain);

    /// Prints a summary of the model parameters.
    virtual void v_GenerateSummary(SolverUtils::SummaryList &s);

private:
    /// Variable diffusivity
    NekDouble m_chi;
    NekDouble m_capMembrane;
    NekDouble m_conductivity;

    CellModelSharedPtr m_cell;

    std::vector<StimulusSharedPtr> m_stimulus;

    Array<OneD, NekDouble> ComputeLaplacianDiff(
        const Array<OneD, const NekDouble> &Laplacian,
        const Array<OneD, const NekDouble> &LaplacianNew);

    Array<OneD, NekDouble> m_epsilon;
    Array<OneD, NekDouble> m_epsu;

    /// Stimulus current
    NekDouble m_stimDuration;

    void LoadStimuli();

    int CountActivation(Array<OneD, int> &Activated);
};

} // namespace Nektar

#endif
