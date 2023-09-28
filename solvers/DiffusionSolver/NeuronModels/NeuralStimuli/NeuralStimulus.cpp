///////////////////////////////////////////////////////////////////////////////
//
// File: NeuralStimulus.cpp
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
// Description: NeuralStimulus base class.
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/VmathArray.hpp>
#include <tinyxml.h>

#include <DiffusionSolver/NeuronModels/NeuralStimuli/NeuralStimulus.h>

using namespace std;

namespace Nektar
{
NeuralStimulusFactory &GetStimulusFactory()
{
    static NeuralStimulusFactory instance;
    return instance;
}

/**
 * @class NeuralStimulus
 *
 * The NeuralStimulus class and derived classes implement a range of stimuli.
 * The stimulus contains input stimuli that can be applied throughout the
 * domain, on specified regions determined by the derived classes of
 * NeuralStimulus, at specified frequencies determined by the derived classes of
 * Protocol.
 *
 */

/**
 * Stimulus base class constructor.
 */
NeuralStimulus::NeuralStimulus(const LibUtilities::SessionReaderSharedPtr &pSession,
                   const MultiRegions::ExpListSharedPtr &pField,
                   const TiXmlElement *pXml)
{
    m_session = pSession;
    m_field   = pField;
    m_nq      = pField->GetTotPoints();

    const TiXmlElement *vProtocol = pXml->FirstChildElement("PROTOCOL");
    string vTypeP                 = vProtocol->Attribute("TYPE");

    m_Protocol =
        GetProtocolFactory().CreateInstance(vTypeP, pSession, vProtocol);

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
void NeuralStimulus::Initialise()
{
}

/**
 *
 */
vector<NeuralStimulusSharedPtr> NeuralStimulus::LoadStimuli(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const MultiRegions::ExpListSharedPtr &pField)
{
    vector<NeuralStimulusSharedPtr> vStimList;

    TiXmlElement *vStimuli = pSession->GetElement("Nektar/Stimuli");
    if (vStimuli)
    {
        TiXmlElement *vStimulus = vStimuli->FirstChildElement("NEURALSTIMULUS");
        while (vStimulus)
        {
            string vType = vStimulus->Attribute("TYPE");

            vStimList.push_back(GetStimulusFactory().CreateInstance(
                vType, pSession, pField, vStimulus));
            vStimulus = vStimulus->NextSiblingElement("NEURALSTIMULUS");
        }
    }
    return vStimList;
}

void NeuralStimulus::RegionalUpdate(const int zonestart,
                                    const int zoneend,
                                    const Array<OneD, const NekDouble> &excitezone,
                                    Array<OneD, Array<OneD, NekDouble>> &outarray,
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

    // Retrieve coordinates of quadrature points
    int nq = m_field->GetNpoints();
    Array<OneD, NekDouble> x0(nq);
    Array<OneD, NekDouble> x1(nq);
    Array<OneD, NekDouble> x2(nq);
    m_field->GetCoords(x0, x1, x2);

    // Get the protocol amplitude
    NekDouble v_amp =
        m_Protocol->GetAmplitude(time) * m_strength / m_chiCapMembrane;

    Array<OneD, NekDouble> vampzone(nq);
    Vmath::Smul(nq, v_amp, excitezone, 1, vampzone, 1);

    switch (dim)
    {
        case 1:
            for (int j = zonestart; j < zoneend; j++)
            {
                outarray[0][j] += vampzone[j] * ((tanh(m_pis * (x0[j] - m_px1)) -
                                            tanh(m_pis * (x0[j] - m_px2))) /
                                           2.0);
            }
            break;
        case 2:
            for (int j = zonestart; j < zoneend; j++)
            {
                outarray[0][j] += vampzone[j] * (((tanh(m_pis * (x0[j] - m_px1)) -
                                             tanh(m_pis * (x0[j] - m_px2))) *
                                            (tanh(m_pis * (x1[j] - m_py1)) -
                                             tanh(m_pis * (x1[j] - m_py2)))) /
                                           2.0);
            }
            break;
        case 3:
            for (int j = zonestart; j < zoneend; j++)
            {
                outarray[0][j] += vampzone[j] * (((tanh(m_pis * (x0[j] - m_px1)) -
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


} // namespace Nektar
