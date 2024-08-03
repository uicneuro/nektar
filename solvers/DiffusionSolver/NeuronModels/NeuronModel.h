///////////////////////////////////////////////////////////////////////////////
//
// File NeuronModel.h
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
// Description: Neuron model base class.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_ADRSOLVER_NeuronMODELS_NeuronMODEL
#define NEKTAR_SOLVERS_ADRSOLVER_NeuronMODELS_NeuronMODEL

#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
//#include <SpatialDomains/SpatialData.h>
#include <MultiRegions/ExpList.h>
#include <SolverUtils/Core/Misc.h>
#include <SolverUtils/MMFSystem.h>
#include <StdRegions/StdNodalTetExp.h>
#include <StdRegions/StdNodalTriExp.h>

#include <SolverUtils/UnsteadySystem.h>

namespace Nektar
{
// Forward declaration
class NeuronModel;

typedef std::vector<std::pair<std::string, std::string>> SummaryList;

/// A shared pointer to an EquationSystem object
typedef std::shared_ptr<NeuronModel> NeuronModelSharedPtr;
/// Datatype of the NekFactory used to instantiate classes derived from
/// the EquationSystem class.
typedef LibUtilities::NekFactory<std::string, NeuronModel,
                                 const LibUtilities::SessionReaderSharedPtr &,
                                 const MultiRegions::ExpListSharedPtr &>
    NeuronModelFactory;
NeuronModelFactory &GetNeuronModelFactory();

/// Neuron model base class.
class NeuronModel
{
public:
    NeuronModel(const LibUtilities::SessionReaderSharedPtr &pSession,
              const MultiRegions::ExpListSharedPtr &pField);

    virtual ~NeuronModel()
    {
    }

    /// Initialise the Neuron model storage and set initial conditions
    void Initialise();

    void TimeIntegrate( const Array<OneD, const int> &zoneindex,
                        const Array<OneD, const NekDouble> &inarray,
                        Array<OneD, NekDouble> &outarray, 
                        const NekDouble time,
                        const NekDouble diameter = 0.001,
                        const NekDouble Tc = 24.0);

    void Update(const Array<OneD, const int> &zoneindex,
                const Array<OneD, const Array<OneD, NekDouble>> &inarray,
                Array<OneD, Array<OneD, NekDouble>> &outarray,
                const NekDouble time, const NekDouble diameter, const NekDouble Tc)
    {
        v_Update(zoneindex, inarray, outarray, time, diameter, Tc);
    }

    /// Print a summary of the Neuron model
    void GenerateSummary(SummaryList &s)
    {
        v_GenerateSummary(s);
    }

    unsigned int GetNumNeuronVariables()
    {
        return m_nvar;
    }

    // index = 0 : myelin
    // index = 1 : nodal, demyelinated 
    NekDouble GetCapacitanceValue(const int index)
    {
        NekDouble capval = 1.0;
        if (index==0)
        {
            capval = m_Cm;
        }

        else if (index==1)
        {
            capval = m_Cn;
        }

        else
        {
            ASSERTL0(false, "Wrong index for capacitance");
        }

        return capval;
    }

    NekDouble GetRecistanceValue()
    {
        return m_ra;
    }

    std::string GetNeuronVarName(unsigned int idx)
    {
        return v_GetNeuronVarName(idx);
    }

    Array<OneD, NekDouble> GetNeuronSolutionCoeffs(unsigned int idx);

    Array<OneD, NekDouble> GetNeuronSolution(unsigned int idx);

protected:

    NekDouble m_pi;
    /// Session
    LibUtilities::SessionReaderSharedPtr m_session;
    /// Transmembrane potential field from PDE system
    MultiRegions::ExpListSharedPtr m_field;
    /// Number of physical points.
    int m_nq;
    /// Number of variables in Neuron model (inc. transmembrane voltage)
    int m_nvar;

    // Myelin capacitance
    NekDouble m_Cm;

    // Nodal and demyelinated capacitance
    NekDouble m_Cn;

    // axoplasmic resistance
    NekDouble m_ra;

    /// Timestep for pde model
    NekDouble m_lastTime;
    /// Number of substeps to take
    int m_substeps;

    /// Neuron model solution variables
    Array<OneD, Array<OneD, NekDouble>> m_NeuronSol;
    
    /// Neuron model integration workspace
    Array<OneD, Array<OneD, NekDouble>> m_wsp;

    /// Flag indicating whether nodal projection in use
    bool m_useNodal;
    /// StdNodalTri for Neuron model calculations
    StdRegions::StdNodalTriExpSharedPtr m_nodalTri;
    StdRegions::StdNodalTetExpSharedPtr m_nodalTet;
    /// Temporary array for nodal projection
    Array<OneD, Array<OneD, NekDouble>> m_nodalTmp;

    /// Indices of Neuron model variables which are concentrations
    std::vector<int> m_concentrations;
    /// Indices of Neuron model variables which are gates
    std::vector<int> m_gates;
    /// Storage for gate tau values
    Array<OneD, Array<OneD, NekDouble>> m_gates_tau;

    virtual void v_Update(
        const Array<OneD, const int> &zoneindex,
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray,
        const NekDouble time, const NekDouble diameter, const NekDouble Tc) = 0;

    virtual void v_GenerateSummary(SummaryList &s) = 0;

    virtual std::string v_GetNeuronVarName(unsigned int idx)
    {
        return "Var" + boost::lexical_cast<std::string>(idx);
    }

    virtual void v_SetInitialConditions() = 0;

    void LoadNeuronModel();
};

} // namespace Nektar

#endif /* NeuronMODEL_H_ */
