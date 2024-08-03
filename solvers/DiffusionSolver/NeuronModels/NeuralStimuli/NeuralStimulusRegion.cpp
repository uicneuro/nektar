///////////////////////////////////////////////////////////////////////////////
//
// File: NeuralStimulusRegion.cpp
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
// Description: Rectangular stimulus class
//
///////////////////////////////////////////////////////////////////////////////
#include <LibUtilities/BasicUtils/VmathArray.hpp>
#include <iostream>
#include <tinyxml.h>

#include <DiffusionSolver/NeuronModels/NeuralStimuli/NeuralStimulusRegion.h>

namespace Nektar
{
std::string NeuralStimulusRegion::className =
    GetStimulusFactory().RegisterCreatorFunction(
        "NeuralStimulusRegion", NeuralStimulusRegion::create, "Regional stimulus.");

/**
 * @class NeuralStimulusRegion
 *
 * The Stimulus class and derived classes implement a range of stimuli.
 * The stimulus contains input stimuli that can be applied throughout the
 * domain, on specified regions determined by the derived classes of
 * Stimulus, at specified frequencies determined by the derived classes of
 * Protocol.
 */

/**
 * Stimulus base class constructor.
 */
NeuralStimulusRegion::NeuralStimulusRegion(const LibUtilities::SessionReaderSharedPtr &pSession,
                           const MultiRegions::ExpListSharedPtr &pField,
                           const TiXmlElement *pXml)
    : NeuralStimulus(pSession, pField, pXml)
{
    m_session = pSession;
    m_field   = pField;
    m_nq      = pField->GetTotPoints();

    NekDouble chi, Cm;
    m_session->LoadParameter("Chi", chi, 28.0);
    m_session->LoadParameter("Cm", Cm, 0.125);

    m_chiCapMembrane = chi * Cm;

    if (!pXml)
    {
        return;
    }

    const TiXmlElement *pXmlparameter;

    // Get the dimension of the expansion
    pXmlparameter = pXml->FirstChildElement("p_x1");
    m_px1         = atof(pXmlparameter->GetText());

    pXmlparameter = pXml->FirstChildElement("p_y1");
    m_py1         = atof(pXmlparameter->GetText());

    pXmlparameter = pXml->FirstChildElement("p_z1");
    m_pz1         = atof(pXmlparameter->GetText());

    pXmlparameter = pXml->FirstChildElement("p_x2");
    m_px2         = atof(pXmlparameter->GetText());

    pXmlparameter = pXml->FirstChildElement("p_y2");
    m_py2         = atof(pXmlparameter->GetText());

    pXmlparameter = pXml->FirstChildElement("p_z2");
    m_pz2         = atof(pXmlparameter->GetText());

    pXmlparameter = pXml->FirstChildElement("p_is");
    m_pis         = atof(pXmlparameter->GetText());

    pXmlparameter = pXml->FirstChildElement("p_strength");
    m_strength    = atof(pXmlparameter->GetText());
}

/**
 * Initialise the stimulus. Allocate workspace and variable storage.
 */
void NeuralStimulusRegion::Initialise()
{
}

/**
 *
 */
void NeuralStimulusRegion::v_Update(const Array<OneD, const NekDouble> &excitezone,
                                Array<OneD, NekDouble> &outarray,
                                const NekDouble time)
{
    if (m_field->GetNumElmts() == 0)
    {
        return;
    }

    // int dim = m_field->GetShapeDimension();
    int dim;
    NekDouble Tol=1.0e-7;
    if(fabs(m_pz1-m_pz2)<Tol)
    {
        dim = 2;
        if(fabs(m_py1-m_py2)<Tol)
        {
            dim = 1;
        }
    }

    // const NekDouble var_membrane__cml = 1.87e-11; 
    const NekDouble var_membrane__cnd = 3.14e-9;  // (F/cm) Waxman 1978, Same

    // Retrieve coordinates of quadrature points
    int nq = m_field->GetNpoints();
    
    Array<OneD, NekDouble> x0(nq);
    Array<OneD, NekDouble> x1(nq);
    Array<OneD, NekDouble> x2(nq);
    m_field->GetCoords(x0, x1, x2);

    // Get the protocol amplitude
    NekDouble v_amp = m_Protocol->GetAmplitude(time) * m_strength / ( m_chiCapMembrane * var_membrane__cnd );

    Array<OneD, NekDouble> vampzone(nq);
    Vmath::Smul(nq, v_amp, excitezone, 1, vampzone, 1);

    switch (dim)
    {
        case 1:
            for (int j = 0; j < nq; j++)
            {
                outarray[j] += vampzone[j] * ((tanh(m_pis * (x0[j] - m_px1)) -
                                            tanh(m_pis * (x0[j] - m_px2))) /
                                           2.0);
            }
            break;

        case 2:
            for (int j = 0; j < nq; j++)
            {
                outarray[j] += vampzone[j] * (((tanh(m_pis * (x0[j] - m_px1)) -
                                             tanh(m_pis * (x0[j] - m_px2))) *
                                            (tanh(m_pis * (x1[j] - m_py1)) -
                                             tanh(m_pis * (x1[j] - m_py2)))) /
                                           2.0);
            }
            break;

        case 3:
            for (int j = 0; j < nq; j++)
            {
                outarray[j] += vampzone[j] * (((tanh(m_pis * (x0[j] - m_px1)) -
                                             tanh(m_pis * (x0[j] - m_px2))) *
                                            (tanh(m_pis * (x1[j] - m_py1)) -
                                             tanh(m_pis * (x1[j] - m_py2))) *
                                            (tanh(m_pis * (x2[j] - m_pz1)) -
                                             tanh(m_pis * (x2[j] - m_pz2)))) /
                                           2.0);
            }
            break;
    }
}

/**
 *
 */
void NeuralStimulusRegion::v_GenerateSummary(SolverUtils::SummaryList &s)
{
    boost::ignore_unused(s);
}
} // namespace Nektar
