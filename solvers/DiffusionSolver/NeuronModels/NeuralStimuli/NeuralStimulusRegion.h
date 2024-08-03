///////////////////////////////////////////////////////////////////////////////
//
// File: NeuralStimulusRegion.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
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
// Description: Rectangular stimulus header file
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_NEURALSOLVER_STIMULI_NEURALSTIMULUSREGION
#define NEKTAR_SOLVERS_NEURALSOLVER_STIMULI_NEURALSTIMULUSREGION

#include <DiffusionSolver/NeuronModels/NeuralStimuli/NeuralStimulus.h>
#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <MultiRegions/ExpList.h>

namespace Nektar
{

/// Protocol base class.
class NeuralStimulusRegion : public NeuralStimulus
{
public:
    /// Creates an instance of this class
    static NeuralStimulusSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const MultiRegions::ExpListSharedPtr &pField, const TiXmlElement *pXml)
    {
        return MemoryManager<NeuralStimulusRegion>::AllocateSharedPtr(pSession, pField,
                                                              pXml);
    }

    /// Name of class
    static std::string className;

    friend class MemoryManager<NeuralStimulusRegion>;

    virtual ~NeuralStimulusRegion()
    {
    }

    /// Initialise the stimulus storage and set initial conditions
    void Initialise();

protected:
    NekDouble m_px1;
    NekDouble m_py1;
    NekDouble m_pz1;
    NekDouble m_px2;
    NekDouble m_py2;
    NekDouble m_pz2;
    NekDouble m_pis;
    NekDouble m_strength;
    NekDouble m_chiCapMembrane;

    virtual void v_Update(const Array<OneD, const NekDouble> &excitezone,
                          Array<OneD, NekDouble> &outarray,
                                const NekDouble time) override;

    virtual void v_GenerateSummary(SolverUtils::SummaryList &s) override;

private:
    NeuralStimulusRegion(const LibUtilities::SessionReaderSharedPtr &pSession,
                 const MultiRegions::ExpListSharedPtr &pField,
                 const TiXmlElement *pXml);
};
} // namespace Nektar

#endif
