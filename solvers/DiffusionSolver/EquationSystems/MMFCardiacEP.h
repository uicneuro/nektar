///////////////////////////////////////////////////////////////////////////////
//
// File MMFCardiacEP.h
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
// Description: MMFCardiacEP
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_MMFCARDIACEP_H
#define NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_MMFCARDIACEP_H

#include <SolverUtils/MMFSystem.h>
#include <SolverUtils/UnsteadySystem.h>

#include <CardiacEPSolver/CellModels/CellModel.h>
#include <CardiacEPSolver/Stimuli/Stimulus.h>

namespace Nektar
{

enum CardiacEPType
{
    eCardiacEP1D,
    eCardiacEP2D,
    eCardiacEPSpiral,
    SIZE_CardiacEPType ///< Length of enum list
};

const char *const CardiacEPTypeMap[] = {
    "CardiacEP1D", "CardiacEP2D", "CardiacEPSpiral",
};

enum SolverSchemeType
{
    eDefault,
    eMMFFirst,
    eTimeMap,
    SIZE_SolverSchemeType,
};

const char *const SolverSchemeTypeMap[] = {
    "Default", "MMFFirst", "TimeMap",
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
class MMFCardiacEP : public SolverUtils::MMFSystem
{
public:
    friend class MemoryManager<MMFCardiacEP>;

    /// Creates an instance of this class
    static SolverUtils::EquationSystemSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const SpatialDomains::MeshGraphSharedPtr &pGraph)
    {
        SolverUtils::EquationSystemSharedPtr p =
            MemoryManager<MMFCardiacEP>::AllocateSharedPtr(pSession, pGraph);
        p->InitObject();
        return p;
    }

    /// Name of class
    static std::string className;

    CardiacEPType m_CardiacEPType;
    SolverSchemeType m_SolverSchemeType;

    virtual void v_InitObject(bool DeclareField = true) override;
    virtual void v_DoSolve() override;

    /// Desctructor
    virtual ~MMFCardiacEP();

protected:
    int m_Convectiven;

    TimeMapType m_TimeMap;
    NekDouble m_TimeMapStart;
    NekDouble m_TimeMapEnd;

    NekDouble m_Diffbeta, m_Diffeta, m_Diffhe;   // h_e for LDG

    // Scar tisseu related variables
    NekDouble m_PVcond;
    NekDouble m_ScarSize, m_ScarStr, m_ScarPis, m_ScarLocx, m_ScarLocy, m_ScarLocz;
    NekDouble m_RelDivSize, m_RelDivStr, m_RelDivPis, m_RelDivLocx;

    // Relative divergence related variables
    NekDouble m_LambDivSmoothL;

    Array<OneD, int> m_ValidTimeMap;
    
    /// Constructor
    MMFCardiacEP(const LibUtilities::SessionReaderSharedPtr &pSession,
                 const SpatialDomains::MeshGraphSharedPtr &pGraph);

    InitWaveType m_InitWaveType;

    void DoSolveMMF();
    void DoSolveMMFFirst();

    Array<OneD, NekDouble> ReadFibermap(const NekDouble AnisotropyStrength, Array<OneD, NekDouble> &CardiacFibre);
    Array<OneD, NekDouble> ReadConductivityMap();

    void LoadCardiacFiber(
                const SolverUtils::MediumType CardiacMediumType,
                const NekDouble AnisotropyStrength, 
                Array<OneD, Array<OneD, NekDouble>> &AniStrength,
                Array<OneD, NekDouble> &CardiacFibre = NullNekDouble1DArray);

    // Coefficients for Anisotropy
    int m_AnisotropyRegion;
    NekDouble m_AnisotropyStrength;

    void DoImplicitSolveCardiacEP(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time,
        const NekDouble lambda);

    void DoNullSolve(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time,
        const NekDouble lambda);

    void DoOdeRhsCardiacEP(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time);

    void DoOdeProjection(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time);

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

private:
    /// Variable diffusivity
    NekDouble m_chi;
    NekDouble m_capMembrane;
    NekDouble m_conductivity;

    CellModelSharedPtr m_cell;

    std::vector<StimulusSharedPtr> m_stimulus;
    std::vector<StimulusSharedPtr> m_fiberstimulus;

    /// Stimulus current
    NekDouble m_stimDuration;
};

} // namespace Nektar

#endif
