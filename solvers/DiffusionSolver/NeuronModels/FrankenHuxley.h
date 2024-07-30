///////////////////////////////////////////////////////////////////////////////
//
// File FrankenHuxley.h
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
// Description: Luo-Rudy 1991 Neuron model
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_FRANKENHUXLEY_H
#define NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_FRANKENHUXLEY_H

#include <DiffusionSolver/NeuronModels/NeuronModel.h>
namespace Nektar
{
    class FrankenHuxley : public NeuronModel
    {

    public:
        /// Creates an instance of this class
        static NeuronModelSharedPtr create(
                const LibUtilities::SessionReaderSharedPtr& pSession,
                const MultiRegions::ExpListSharedPtr& pField)
        {
            return MemoryManager<FrankenHuxley>::AllocateSharedPtr(pSession, pField);
        }

        /// Name of class
        static std::string className;

        /// Constructor
        FrankenHuxley(const LibUtilities::SessionReaderSharedPtr& pSession,
                  const MultiRegions::ExpListSharedPtr& pField);

        /// Destructor
        virtual ~FrankenHuxley() {}

    protected:

        NekDouble ComputeexpM1(const NekDouble x, const NekDouble y);
        NekDouble ComputeIon(const NekDouble E, const NekDouble ci, const NekDouble co, const NekDouble Tc);
        NekDouble efun(const NekDouble z);

        /// Computes the reaction terms $f(u,v)$ and $g(u,v)$.
        virtual void v_Update(
                const Array<OneD, const int> &zoneindex,
                const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
                      Array<OneD,        Array<OneD, NekDouble> >&outarray,
                const NekDouble time, 
                const NekDouble var_membrane__d = 0.001,
                const NekDouble var_membrane__Tc = 24.0);

        /// Prints a summary of the model parameters.
        virtual void v_GenerateSummary(SummaryList& s);

        /// Set initial conditions for the Neuron model
        virtual void v_SetInitialConditions();
    };
}

#endif
