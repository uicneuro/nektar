///////////////////////////////////////////////////////////////////////////////
//
// File: MMFPoisson.cpp
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
// Description: MMFPoisson solve routines
//
///////////////////////////////////////////////////////////////////////////////

#include <DiffusionSolver/EquationSystems/MMFPoisson.h>

using namespace std;

namespace Nektar
{
string MMFPoisson::className1 = GetEquationSystemFactory().RegisterCreatorFunction(
    "MMFPoisson", MMFPoisson::create);
string MMFPoisson::className2 = GetEquationSystemFactory().RegisterCreatorFunction(
    "SteadyDiffusion", MMFPoisson::create);

MMFPoisson::MMFPoisson(const LibUtilities::SessionReaderSharedPtr &pSession,
                 const SpatialDomains::MeshGraphSharedPtr &pGraph)
    : UnsteadySystem(pSession, pGraph), MMFLaplace(pSession, pGraph)
{
    m_factors[StdRegions::eFactorLambda] = 0.0;
    m_factors[StdRegions::eFactorTau]    = 1.0;
}

void MMFPoisson::v_InitObject(bool DeclareFields)
{

    MMFLaplace::v_InitObject(DeclareFields);

    GetFunction("Forcing")->Evaluate(m_session->GetVariables(), m_fields);
    std::cout << "v_InitObject: i = " << 0 << ", phys = " << RootMeanSquare(m_fields[0]->GetPhys()) << std::endl;

}

MMFPoisson::~MMFPoisson()
{
}

void MMFPoisson::v_GenerateSummary(SolverUtils::SummaryList &s)
{
    MMFLaplace::v_GenerateSummary(s);
    
    for (int i = 0; i < m_fields.size(); ++i)
    {
        stringstream name;
        name << "Forcing func [" << i << "]";
        SolverUtils::AddSummaryItem(
            s, name.str(),
            GetFunction("Forcing")->Describe(m_session->GetVariable(0)));
    }
}

Array<OneD, bool> MMFPoisson::v_GetSystemSingularChecks()
{
    return Array<OneD, bool>(m_session->GetVariables().size(), true);
}
} // namespace Nektar
