///////////////////////////////////////////////////////////////////////////////
//
// File: MMFLaplace.cpp
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
// Description: MMFLaplace solve routines
//
///////////////////////////////////////////////////////////////////////////////

#include <iomanip>
#include <iostream>

#include <boost/algorithm/string.hpp>
#include <boost/core/ignore_unused.hpp>

#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>
#include <DiffusionSolver/EquationSystems/MMFDiffusion.h>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>
#include <SolverUtils/Driver.h>
#include <DiffusionSolver/EquationSystems/MMFLaplace.h>

using namespace std;
using namespace Nektar::SolverUtils;
using namespace Nektar;

namespace Nektar
{
string MMFLaplace::className = GetEquationSystemFactory().RegisterCreatorFunction(
    "MMFLaplace", MMFLaplace::create);

MMFLaplace::MMFLaplace(const LibUtilities::SessionReaderSharedPtr &pSession,
                 const SpatialDomains::MeshGraphSharedPtr &pGraph)
    : UnsteadySystem(pSession, pGraph), MMFSystem(pSession, pGraph)
{
    m_factors[StdRegions::eFactorLambda] = 0.0;
    m_factors[StdRegions::eFactorTau]    = 1.0;
}

void MMFLaplace::v_InitObject(bool DeclareFields)
{
    EquationSystem::v_InitObject(DeclareFields);

    int nq    = m_fields[0]->GetNpoints();
    int nvar  = m_fields.size();

    // AniStrength for e^1 and e^2
    m_session->LoadParameter("AniStrength", m_AniStrength, 1.0);
    m_session->LoadParameter("HelmMMF", m_HelmMMF, 0);

    // Diffusivity coefficient for e^j
    m_epsilon = Array<OneD, NekDouble>(m_spacedim);
    m_session->LoadParameter("epsilon0", m_epsilon[0], 1.0);
    m_session->LoadParameter("epsilon1", m_epsilon[1], 1.0);
    m_session->LoadParameter("epsilon2", m_epsilon[2], 1.0);

    // Diffusivity coefficient for u^j
    m_epsu = Array<OneD, NekDouble>(nvar + 1);
    m_session->LoadParameter("epsu0", m_epsu[0], 1.0);
    m_session->LoadParameter("epsu1", m_epsu[1], 1.0);

    int shapedim = m_fields[0]->GetShapeDimension();
    Array<OneD, Array<OneD, NekDouble>> Anisotropy(shapedim);
    for (int j = 0; j < shapedim; ++j)
    {
        Anisotropy[j] = Array<OneD, NekDouble>(nq, 1.0);
        Vmath::Fill(nq, m_AniStrength, &Anisotropy[j][0], 1);
    }

    MMFSystem::MMFInitObject(Anisotropy);

    ComputeVarCoeff2D(m_movingframes,m_varcoeff);

}

MMFLaplace::~MMFLaplace()
{
}

void MMFLaplace::v_GenerateSummary(SolverUtils::SummaryList &s)
{
    EquationSystem::SessionSummary(s);
    SolverUtils::AddSummaryItem(s, "Lambda",
                                m_factors[StdRegions::eFactorLambda]);

    SolverUtils::AddSummaryItem(s, "Lambda",
                                m_factors[StdRegions::eFactorLambda]);

}


void MMFLaplace::v_DoSolve()
{

        GetFunction("Forcing")->Evaluate(m_session->GetVariables(), m_fields);

    for (int i = 0; i < m_fields.size(); ++i)
    {
        std::cout << "v_DoSolve: i = " << i << ", phys = " << RootMeanSquare(m_fields[i]->GetPhys()) << std::endl;
        // Zero field so initial conditions are zero
        Vmath::Zero(m_fields[i]->GetNcoeffs(), m_fields[i]->UpdateCoeffs(), 1);
        if(m_HelmMMF)
        {
            std::cout << "Helmsolver with MMF" << std::endl;            
            m_fields[i]->HelmSolve(m_fields[i]->GetPhys(),
                            m_fields[i]->UpdateCoeffs(), m_factors, m_varcoeff);
        }

        else
        {
            std::cout << "Regular Helmsolver" << std::endl;            
            m_fields[i]->HelmSolve(m_fields[i]->GetPhys(),
                                m_fields[i]->UpdateCoeffs(), m_factors);
        }

        m_fields[i]->SetPhysState(false);
    }
}

Array<OneD, bool> MMFLaplace::v_GetSystemSingularChecks()
{
    return Array<OneD, bool>(m_session->GetVariables().size(), true);
}
} // namespace Nektar
