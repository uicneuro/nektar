///////////////////////////////////////////////////////////////////////////////
//
// File HodgkinHuxley.h
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

#ifndef NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_HODGKINHUXLEY_H
#define NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_HODGKINHUXLEY_H

#include <DiffusionSolver/NeuronModels/NeuronModel.h>
namespace Nektar
{
    class HodgkinHuxley : public NeuronModel
    {

    public:
        /// Creates an instance of this class
        static NeuronModelSharedPtr create(
                const LibUtilities::SessionReaderSharedPtr& pSession,
                const MultiRegions::ExpListSharedPtr& pField)
        {
            return MemoryManager<HodgkinHuxley>::AllocateSharedPtr(pSession, pField);
        }

        /// Name of class
        static std::string className;

        /// Constructor
        HodgkinHuxley(const LibUtilities::SessionReaderSharedPtr& pSession,
                  const MultiRegions::ExpListSharedPtr& pField);

        /// Destructor
        virtual ~HodgkinHuxley() {}

    protected:

        NekDouble Compute_alpha(const NekDouble membrane_V, const NekDouble alphaA, const NekDouble alphaB, const NekDouble alphaC);
        NekDouble Compute_beta(const NekDouble membrane_V, const NekDouble betaA, const NekDouble betaB, const NekDouble betaC);
        NekDouble Compute_beta_h(const NekDouble membrane_V, const NekDouble betaA, const NekDouble betaB, const NekDouble betaC);
        NekDouble ComputeZfunction(const NekDouble membrane_V, const NekDouble Yexp, const NekDouble Yint, const NekDouble Tc);
        NekDouble RootMeanSquare(const Array<OneD, const NekDouble> &inarray, const int Ntot = 1);

        /// Computes the reaction terms $f(u,v)$ and $g(u,v)$.
        virtual void v_Update(
                const Array<OneD, const int> &RvNodeZone,
                const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
                      Array<OneD,        Array<OneD, NekDouble> >&outarray,
                const NekDouble time, const NekDouble diameter, const NekDouble Tc);

        /// Prints a summary of the model parameters.
        virtual void v_GenerateSummary(SummaryList& s);

        /// Set initial conditions for the Neuron model
        virtual void v_SetInitialConditions();
    };
}

#endif
