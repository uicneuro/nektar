///////////////////////////////////////////////////////////////////////////////
//
// File NeuronModel.cpp
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

#include <LibUtilities/BasicUtils/VmathArray.hpp>
#include <DiffusionSolver/NeuronModels/NeuronModel.h>
#include <StdRegions/StdNodalTriExp.h>
#include <SolverUtils/MMFSystem.h>

#include <DiffusionSolver/EquationSystems/MMFDiffusion.h>

//#include <LibUtilities/LinearAlgebra/Blas.hpp>

using namespace std;

namespace Nektar
{
NeuronModelFactory &GetNeuronModelFactory()
{
    static NeuronModelFactory instance;
    return instance;
}

/**
 * @class NeuronModel
 *
 * The NeuronModel class and derived classes implement a range of Neuron model
 * ODE systems. A Neuron model comprises a system of ion concentration
 * variables and zero or more gating variables. Gating variables are
 * time-integrated using the Rush-Larsen method and for each variable y,
 * the corresponding y_inf and tau_y value is computed by Update(). The tau
 * values are stored in separate storage to inarray/outarray, #m_gates_tau.
 */

/**
 * Neuron model base class constructor.
 */
NeuronModel::NeuronModel(const LibUtilities::SessionReaderSharedPtr &pSession,
                     const MultiRegions::ExpListSharedPtr &pField)
{
    m_session  = pSession;
    m_field    = pField;
    m_lastTime = 0.0;
    m_useNodal = false;
    m_pi       = 3.14159265358979323846;

    m_session->LoadParameter("Substeps", m_substeps, 1);

    // Number of points in nodal space is the number of coefficients
    // in modified basis
    std::set<enum LibUtilities::ShapeType> s;
    for (unsigned int i = 0; i < m_field->GetNumElmts(); ++i)
    {
        s.insert(m_field->GetExp(i)->DetShapeType());
    }

    m_nq = pField->GetTotPoints();
}

/**
 * Initialise the Neuron model. Allocate workspace and variable storage.
 */
void NeuronModel::Initialise()
{
    ASSERTL1(m_nvar > 0, "Neuron model must have at least 1 variable.");

    m_NeuronSol = Array<OneD, Array<OneD, NekDouble>>(m_nvar);
    m_wsp     = Array<OneD, Array<OneD, NekDouble>>(m_nvar);
    for (unsigned int i = 0; i < m_nvar; ++i)
    {
        m_NeuronSol[i] = Array<OneD, NekDouble>(m_nq);
        m_wsp[i]     = Array<OneD, NekDouble>(m_nq);
    }
    m_gates_tau = Array<OneD, Array<OneD, NekDouble>>(m_gates.size());
    for (unsigned int i = 0; i < m_gates.size(); ++i)
    {
        m_gates_tau[i] = Array<OneD, NekDouble>(m_nq);
    }

    if (m_session->DefinesFunction("NeuronModelInitialConditions"))
    {
      //  LoadNeuronModel();
    }
    else
    {
       v_SetInitialConditions();
    }
}

/**
 * Integrates the Neuron model for one PDE time-step. Neuron model is
 * sub-stepped.
 *
 * Ion concentrations and membrane potential are integrated using forward
 * Euler, while gating variables are integrated using the Rush-Larsen
 * scheme.
 */
void NeuronModel::TimeIntegrate(
    const Array<OneD, const int> &zoneindex,
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray, 
    const NekDouble time,
    const NekDouble diameter,
    const NekDouble Tc)
{
    // int phys_offset = 0;
    // int coef_offset = 0;
    Array<OneD, NekDouble> tmp;

    NekDouble delta_t = (time - m_lastTime) / m_substeps;

    // Copy new transmembrane potential into Neuron model
    Vmath::Vcopy(m_nq, inarray, 1, m_NeuronSol[0], 1);

    // Perform final Neuron model step : m_wsp is the Reaction function from
    // m_NeuronSol[0] of membrane potential.
    Update(zoneindex, m_NeuronSol, m_wsp, time, diameter, Tc);

    // Output dV/dt from last step but integrate remaining Neuron model vars
    // Transform Neuron model I_total from nodal to modal space
    Vmath::Vcopy(m_nq, m_wsp[0], 1, outarray, 1);

    // Ion concentrations
    for (unsigned int j = 0; j < m_concentrations.size(); ++j)
    {
        Vmath::Svtvp(m_nq, delta_t, m_wsp[m_concentrations[j]], 1,
                     m_NeuronSol[m_concentrations[j]], 1,
                     m_NeuronSol[m_concentrations[j]], 1);
    }

    // Gating variables: Rush-Larsen scheme:
    // y_i = y_i^{infty} - ( y_i^{\infty} - y_i (0) ) * exp (-dt / tau_i ) 
    // m_wsp = _inf,  m_NeuronSol = y_i (0)
    for (unsigned int j = 0; j < m_gates.size(); ++j)
    {
        Vmath::Sdiv(m_nq, -delta_t, m_gates_tau[j], 1, m_gates_tau[j], 1);
        Vmath::Vexp(m_nq, m_gates_tau[j], 1, m_gates_tau[j], 1);
        Vmath::Vsub(m_nq, m_NeuronSol[m_gates[j]], 1, m_wsp[m_gates[j]], 1, m_NeuronSol[m_gates[j]], 1);
        Vmath::Vvtvp(m_nq, m_NeuronSol[m_gates[j]], 1, m_gates_tau[j], 1, m_wsp[m_gates[j]], 1, m_NeuronSol[m_gates[j]], 1);
    }

    m_lastTime = time;
}


Array<OneD, NekDouble> NeuronModel::GetNeuronSolutionCoeffs(unsigned int idx)
{
    ASSERTL0(idx < m_nvar, "Index out of range for Neuron model.");

    Array<OneD, NekDouble> outarray(m_field->GetNcoeffs());

    m_field->FwdTrans(m_NeuronSol[idx], outarray);

    return outarray;
}

Array<OneD, NekDouble> NeuronModel::GetNeuronSolution(unsigned int idx)
{
    return m_NeuronSol[idx];
}

void NeuronModel::LoadNeuronModel()
{
    const bool root           = (m_session->GetComm()->GetRank() == 0);
    const std::string fncName = "NeuronModelInitialConditions";
    const int nvar            = m_NeuronSol[0].size();
    std::string varName;
    Array<OneD, NekDouble> coeffs(m_field->GetNcoeffs());
    Array<OneD, NekDouble> tmp;
    int j = 0;

    SpatialDomains::MeshGraphSharedPtr vGraph = m_field->GetGraph();

    if (root)
    {
        cout << "Neuron model initial conditions: " << endl;
    }

    // First determine all the files we need to load
    std::set<std::string> filelist;
    for (j = 1; j < nvar; ++j)
    {
        // Get the name of the jth variable
        varName = GetNeuronVarName(j);

        if (m_session->GetFunctionType(fncName, varName) ==
            LibUtilities::eFunctionTypeFile)
        {
            filelist.insert(m_session->GetFunctionFilename(fncName, varName));
        }
    }

    // Read files
    typedef std::vector<LibUtilities::FieldDefinitionsSharedPtr> FDef;
    typedef std::vector<std::vector<NekDouble>> FData;
    std::map<std::string, FDef> FieldDef;
    std::map<std::string, FData> FieldData;
    LibUtilities::FieldMetaDataMap fieldMetaDataMap;

    for (auto &setIt : filelist)
    {
        if (root)
        {
            cout << "  - Reading file: " << setIt << endl;
        }
        FieldDef[setIt]  = FDef(0);
        FieldData[setIt] = FData(0);
        LibUtilities::FieldIOSharedPtr fld =
            LibUtilities::FieldIO::CreateForFile(m_session, setIt);
        fld->Import(setIt, FieldDef[setIt], FieldData[setIt], fieldMetaDataMap);
    }

    // Get time of checkpoint from file if available
    auto iter = fieldMetaDataMap.find("Time");
    if (iter != fieldMetaDataMap.end())
    {
        m_lastTime = boost::lexical_cast<NekDouble>(iter->second);
    }

    // Load each Neuron model variable
    // j=0 and j=1 are for transmembrane or intra/extra-Neuronular volt.
    Vmath::Zero(m_nq, m_NeuronSol[0], 1);
    for (j = 1; j < m_NeuronSol.size(); ++j)
    {
        // Get the name of the jth variable
        varName = GetNeuronVarName(j);

        // Check if this variable is defined in a file or analytically
        if (m_session->GetFunctionType(fncName, varName) ==
            LibUtilities::eFunctionTypeFile)
        {
            const std::string file =
                m_session->GetFunctionFilename(fncName, varName);

            if (root)
            {
                cout << "  - Field " << varName << ": from file " << file
                     << endl;
            }

            // Extract the data into the modal coefficients
            for (int i = 0; i < FieldDef[file].size(); ++i)
            {
                m_field->ExtractDataToCoeffs(
                    FieldDef[file][i], FieldData[file][i], varName, coeffs);
            }

            // If using nodal Neuron model then we do a modal->nodal transform
            // otherwise we do a backward transform onto physical points.
            if (m_useNodal)
            {
                for (unsigned int i = 0; i < m_field->GetNumElmts(); ++i)
                {
                    int coef_offset = m_field->GetCoeff_Offset(i);
                    if (m_field->GetExp(0)->DetShapeType() ==
                        LibUtilities::eTriangle)
                    {
                        m_nodalTri->ModalToNodal(coeffs + coef_offset,
                                                 tmp = m_NeuronSol[j] +
                                                       coef_offset);
                    }
                    else
                    {
                        m_nodalTet->ModalToNodal(coeffs + coef_offset,
                                                 tmp = m_NeuronSol[j] +
                                                       coef_offset);
                    }
                }
            }
            else
            {
                m_field->BwdTrans(coeffs, m_NeuronSol[j]);
            }
        }
        else if (m_session->GetFunctionType(fncName, varName) ==
                 LibUtilities::eFunctionTypeExpression)
        {
            LibUtilities::EquationSharedPtr equ =
                m_session->GetFunction(fncName, varName);

            if (root)
            {
                cout << "  - Field " << varName << ": " << equ->GetExpression()
                     << endl;
            }

            const unsigned int nphys = m_field->GetNpoints();
            Array<OneD, NekDouble> x0(nphys);
            Array<OneD, NekDouble> x1(nphys);
            Array<OneD, NekDouble> x2(nphys);
            m_field->GetCoords(x0, x1, x2);

            if (m_useNodal)
            {
                Array<OneD, NekDouble> phys(nphys);
                Array<OneD, NekDouble> tmpCoeffs(
                    max(m_nodalTri->GetNcoeffs(), m_nodalTet->GetNcoeffs()));

                equ->Evaluate(x0, x1, x2, phys);
                for (unsigned int i = 0; i < m_field->GetNumElmts(); ++i)
                {
                    int phys_offset = m_field->GetPhys_Offset(i);
                    int coef_offset = m_field->GetCoeff_Offset(i);
                    if (m_field->GetExp(0)->DetShapeType() ==
                        LibUtilities::eTriangle)
                    {
                        m_field->GetExp(0)->FwdTrans(phys + phys_offset,
                                                     tmpCoeffs);
                        m_nodalTri->ModalToNodal(tmpCoeffs, tmp = m_NeuronSol[j] +
                                                                  coef_offset);
                    }
                    else
                    {
                        m_field->GetExp(0)->FwdTrans(phys + phys_offset,
                                                     tmpCoeffs);
                        m_nodalTet->ModalToNodal(tmpCoeffs, tmp = m_NeuronSol[j] +
                                                                  coef_offset);
                    }
                }
            }
            else
            {
                equ->Evaluate(x0, x1, x2, m_NeuronSol[j]);
            }
        }
    }
}

} // namespace Nektar
