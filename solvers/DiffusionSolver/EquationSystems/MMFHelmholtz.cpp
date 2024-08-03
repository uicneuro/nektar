///////////////////////////////////////////////////////////////////////////////
//
// File: MMFHelmholtz.cpp
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
// Description: MMFHelmholtz solve routines
//
///////////////////////////////////////////////////////////////////////////////

#include <DiffusionSolver/EquationSystems/MMFHelmholtz.h>

using namespace std;

namespace Nektar
{
string MMFHelmholtz::className1 =
    GetEquationSystemFactory().RegisterCreatorFunction("MMFHelmholtz",
                                                       MMFHelmholtz::create);
string MMFHelmholtz::className2 =
    GetEquationSystemFactory().RegisterCreatorFunction(
        "SteadyDiffusionReaction", MMFHelmholtz::create);

MMFHelmholtz::MMFHelmholtz(const LibUtilities::SessionReaderSharedPtr &pSession,
                     const SpatialDomains::MeshGraphSharedPtr &pGraph)
    : UnsteadySystem(pSession, pGraph), MMFPoisson(pSession, pGraph)
{
    if (pSession->DefinesParameter("Lambda"))
    {
        m_factors[StdRegions::eFactorLambda] =
            m_session->GetParameter("Lambda");
    }
}

void MMFHelmholtz::v_InitObject(bool DeclareFields)
{
    MMFPoisson::v_InitObject(DeclareFields);
}

MMFHelmholtz::~MMFHelmholtz()
{
}

void MMFHelmholtz::v_GenerateSummary(SolverUtils::SummaryList &s)
{
    MMFPoisson::v_GenerateSummary(s);
}

Array<OneD, bool> MMFHelmholtz::v_GetSystemSingularChecks()
{
    if (m_factors[StdRegions::eFactorLambda] == 0)
    {
        return Array<OneD, bool>(m_session->GetVariables().size(), true);
    }
    else
    {
        return Array<OneD, bool>(m_session->GetVariables().size(), false);
    }
}
} // namespace Nektar
