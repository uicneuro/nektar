///////////////////////////////////////////////////////////////////////////////
//
// File MMFWMTract.h
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
// Description: MMF Maxwell solve routines
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_MMFWMTract_H
#define NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_MMFWMTract_H

#include <SolverUtils/MMFSystem.h>
#include <SolverUtils/UnsteadySystem.h>

enum WMTractType
{
    eDefault,
    SIZE_WMTractType ///< Length of enum list
};

const char *const WMTractTypeMap[] = {
    "Default", 
};

namespace Nektar
{
class MMFWMTract : public SolverUtils::MMFSystem
{
public:
    friend class MemoryManager<MMFWMTract>;

    /// Creates an instance of this class
    static SolverUtils::EquationSystemSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const SpatialDomains::MeshGraphSharedPtr& pGraph)
    {
        SolverUtils::EquationSystemSharedPtr p =
            MemoryManager<MMFWMTract>::AllocateSharedPtr(pSession, pGraph);
        p->InitObject();
        return p;
    }
    /// Name of class
    static std::string className;

    /// Initialise the object
    virtual void v_InitObject(bool DeclareFields = true) override;
    virtual void v_DoSolve() override;

    /// Destructor
    virtual ~MMFWMTract();

protected:

    WMTractType m_WMTractType;

    /// Session reader
    MMFWMTract(const LibUtilities::SessionReaderSharedPtr &pSession,
               const SpatialDomains::MeshGraphSharedPtr& pGraph);

    void PlotDTIVector(const Array<OneD, const NekDouble> &DTIVector);

    /// Compute the RHS
    void DoOdeRhs(const Array<OneD, const Array<OneD, NekDouble>> &inarray,
                  Array<OneD, Array<OneD, NekDouble>> &outarray,
                  const NekDouble time);

    /// Compute the projection
    void DoOdeProjection(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time);

    void read_csv(std::string fname, Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &csv);

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
}

#endif
