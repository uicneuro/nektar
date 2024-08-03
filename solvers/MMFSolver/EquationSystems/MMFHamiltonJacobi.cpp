/////////////////////////////////////////////////////////////////////////////
//
// File MMFHamiltonJacobi.cpp
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

#include <MMFSolver/EquationSystems/MMFHamiltonJacobi.h>

#include <boost/core/ignore_unused.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string.hpp>
#include <iomanip>
#include <iostream>

#include <LibUtilities/BasicUtils/Timer.h>
#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>
#include <SolverUtils/MMFSystem.h>

#include <typeinfo>

#include <LibUtilities/BasicUtils/Timer.h>
#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>

namespace Nektar
{
    std::string MMFHamiltonJacobi::className =
        SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
            "MMFHamiltonJacobi", MMFHamiltonJacobi::create, "MMFHamiltonJacobi equation.");

    MMFHamiltonJacobi::MMFHamiltonJacobi(const LibUtilities::SessionReaderSharedPtr &pSession,
                        const SpatialDomains::MeshGraphSharedPtr& pGraph)
        : UnsteadySystem(pSession, pGraph), MMFSystem(pSession, pGraph)
    {
    }

    /**
     * @brief Initialisation object for the unsteady linear advection equation.
     */
    void MMFHamiltonJacobi::v_InitObject(bool DeclareFields)
    {

        int nq       = m_fields[0]->GetNpoints();

        // Call to the initialisation object
        UnsteadySystem::v_InitObject(DeclareFields);

        // Add Rectangular PML
        Array<OneD, Array<OneD, NekDouble>> AniStrength(m_expdim);
        for (int j = 0; j < m_expdim; ++j)
        {
            AniStrength[j] = Array<OneD, NekDouble>(nq, 1.0);
        }
        MMFSystem::MMFInitObject(AniStrength);

        // SetUp alpha, beta, gamma for Lax-Friedrich constants
        m_LFfactors = Array<OneD, NekDouble>(m_expdim);
        m_session->LoadParameter("LFalpha", m_LFfactors[0], 0.0);

        if(m_expdim>1)
        {
            m_session->LoadParameter("LFbeta", m_LFfactors[1], 0.0);
        }

        if(m_expdim>2)
        {
            m_session->LoadParameter("LFgamma", m_LFfactors[2], 0.0);
        }

        // Define TestMaxwellType
        if (m_session->DefinesSolverInfo("HAMILTONJACOBYTYPE"))
        {
            std::string HamiltonJacobiTypeStr =
                m_session->GetSolverInfo("HAMILTONJACOBYTYPE");
            for (int i = 0; i < (int)SIZE_HamiltonJacobiType; ++i)
            {
                if (HamiltonJacobiTypeMap[i] == HamiltonJacobiTypeStr)
                {
                    m_HJtype = (HamiltonJacobiType)i;
                    break;
                }
            }
        }

        else
        {
            m_HJtype = (HamiltonJacobiType)0;
        }

        // If explicit it computes RHS and PROJECTION for the time integration
        m_ode.DefineOdeRhs(&MMFHamiltonJacobi::DoOdeRhs, this);
        m_ode.DefineProjection(&MMFHamiltonJacobi::DoOdeProjection, this);
        }

    /**
     * @brief Unsteady linear advection equation destructor.
     */
    MMFHamiltonJacobi::~MMFHamiltonJacobi()
    {
    }

    void MMFHamiltonJacobi::v_DoSolve()
    {
        ASSERTL0(m_intScheme != 0, "No time integration scheme.");

        int i, nchk = 1;
        int nvariables = 0;
        int nfields    = m_fields.size();
        int nq         = m_fields[0]->GetNpoints();

        if (m_intVariables.empty())
        {
            for (i = 0; i < nfields; ++i)
            {
                m_intVariables.push_back(i);
            }
            nvariables = nfields;
        }
        else
        {
            nvariables = m_intVariables.size();
        }

        // Set up wrapper to fields data storage.
        Array<OneD, Array<OneD, NekDouble>> fields(nvariables);
        Array<OneD, Array<OneD, NekDouble>> tmp(nvariables);

        // Order storage to list time-integrated fields first.
        for (i = 0; i < nvariables; ++i)
        {
            fields[i] = m_fields[m_intVariables[i]]->GetPhys();
            m_fields[m_intVariables[i]]->SetPhysState(false);
        }

        // Initialise time integration scheme
        m_intScheme->InitializeScheme(m_timestep, fields, m_time, m_ode);

        // Check uniqueness of checkpoint output
        ASSERTL0((m_checktime == 0.0 && m_checksteps == 0) ||
                    (m_checktime > 0.0 && m_checksteps == 0) ||
                    (m_checktime == 0.0 && m_checksteps > 0),
                "Only one of IO_CheckTime and IO_CheckSteps "
                "should be set!");

        LibUtilities::Timer timer;
        bool doCheckTime  = false;
        int step          = 0;
        NekDouble intTime = 0.0;
        NekDouble cpuTime = 0.0;
        NekDouble elapsed = 0.0;

        int Ntot = m_steps / m_checksteps + 1;
        Array<OneD, NekDouble> dMass(Ntot);

        Array<OneD, NekDouble> zeta(nq);
        Array<OneD, Array<OneD, NekDouble>> fieldsprimitive(nvariables);
        for (int i = 0; i < nvariables; ++i)
        {
            fieldsprimitive[i] = Array<OneD, NekDouble>(nq);
        }

        while (step < m_steps || m_time < m_fintime - NekConstants::kNekZeroTol)
        {
            timer.Start();
            fields = m_intScheme->TimeIntegrate(step, m_timestep, m_ode);
            timer.Stop();

            m_time += m_timestep;
            elapsed = timer.TimePerTest(1);
            intTime += elapsed;
            cpuTime += elapsed;

            // Write out status information
            if (m_session->GetComm()->GetRank() == 0 && !((step + 1) % m_infosteps))
            {
                std::cout << "Steps: " << std::setw(5) << std::left << step + 1
                        << " "
                        << "Time: " << std::setw(6) << std::left << m_time;

                std::stringstream ss;
                ss << cpuTime << "s";
                std::cout << " CPU Time: " << std::setw(8) << std::left << ss.str()<< std::endl;

                std::cout << "u = " << RootMeanSquare(fields[0]) << ", max u = " << Vmath::Vamax(nq, fields[0], 1) 
                << ", endval = " << fields[0][0] - fields[0][nq-1] << std::endl;

                cpuTime = 0.0;
            }

            // Transform data into coefficient space
            for (i = 0; i < nvariables; ++i)
            {
                m_fields[m_intVariables[i]]->SetPhys(fields[i]);
                m_fields[m_intVariables[i]]->FwdTransLocalElmt(
                    fields[i], m_fields[m_intVariables[i]]->UpdateCoeffs());
                m_fields[m_intVariables[i]]->SetPhysState(false);
            }

            // Write out checkpoint files
            if ((m_checksteps && step && !((step + 1) % m_checksteps)) ||
                doCheckTime)
            {
                Array<OneD, NekDouble> exactsoln;
                v_EvaluateExactSolution(0, exactsoln, m_time);                   
                Checkpoint_Output_Error(nchk++, fields[0], exactsoln);
                doCheckTime = false;
            }

            // Step advance
            ++step;
        }

        // Print out summary statistics
        if (m_session->GetComm()->GetRank() == 0)
        {
            if (m_cflSafetyFactor > 0.0)
            {
                std::cout << "CFL safety factor : " << m_cflSafetyFactor
                        << std::endl
                        << "CFL time-step     : " << m_timestep << std::endl;
            }

            if (m_session->GetSolverInfo("Driver") != "SteadyState")
            {
                std::cout << "Time-integration  : " << intTime << "s" << std::endl;
            }
        }

        for (i = 0; i < nvariables; ++i)
        {
            m_fields[m_intVariables[i]]->SetPhys(fields[i]);
            m_fields[m_intVariables[i]]->SetPhysState(true);
        }

        for (i = 0; i < nvariables; ++i)
        {
            m_fields[i]->FwdTrans(m_fields[i]->GetPhys(), m_fields[i]->UpdateCoeffs());
        }
    }


    /**
     * @brief Compute the projection for the linear advection equation.
     *
     * @param inarray    Given fields.
     * @param outarray   Calculated solution.
     * @param time       Time.
     */
    void MMFHamiltonJacobi::DoOdeProjection(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time)
    {
        boost::ignore_unused(time);

        int var = inarray.size();

        SetBoundaryConditions(time);

        int nq = GetNpoints();
        for (int i = 0; i < var; ++i)
        {
            Vmath::Vcopy(nq, inarray[i], 1, outarray[i], 1);
        }
    }

    /**
     * @brief Compute the right-hand side for the linear advection equation.
     *
     * @param inarray    Given fields.
     * @param outarray   Calculated solution.
     * @param time       Time.
     */
    void MMFHamiltonJacobi::DoOdeRhs(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time)
    {
        boost::ignore_unused(time);

        int nvar    = inarray.size();
        int ncoeffs = GetNcoeffs();
        int nq      = GetTotPoints();

        Array<OneD, Array<OneD, NekDouble>> physarray(nvar);
        Array<OneD, Array<OneD, NekDouble>> modarray(nvar);
        for (int i = 0; i < nvar; ++i)
        {
            physarray[i] = Array<OneD, NekDouble>(nq);
            modarray[i]  = Array<OneD, NekDouble>(ncoeffs, 0.0);

            Vmath::Vcopy(nq, &inarray[i][0], 1, &physarray[i][0], 1);
        }
        // Compute HamiltonJacobian
       // WeakDGHamiltonJacobi(m_HJalpha, m_HJbeta, physarray, outarray, time);
        WeakDGDirectionalHamiltonJacobi(m_LFfactors, physarray, outarray, time);
        for (int i = 0; i < nvar; ++i)
        {
             Vmath::Neg(nq, outarray[i], 1);
        }

        // std::cout << "outarray[0] - outarray[end] = " << outarray[0][0] - outarray[0][nq-1] << std::endl;

        // WeakDGHamiltonJacobi(m_HJalpha, m_HJbeta, physarray, modarray, time);

        // for (int i = 0; i < nvar; ++i)
        // {
        //     Vmath::Neg(ncoeffs, modarray[i], 1);
        //     m_fields[i]->MultiplyByElmtInvMass(modarray[i],modarray[i]);
        //     m_fields[i]->BwdTrans(modarray[i], outarray[i]);
        // }
    }


    void MMFHamiltonJacobi::WeakDGDirectionalHamiltonJacobi(
        const Array<OneD, const NekDouble> &LFfactors, 
        const Array<OneD, const Array<OneD, NekDouble>> &InField,
        Array<OneD, Array<OneD, NekDouble>> &OutField, const NekDouble time)
    {
        switch(m_expdim)
        {
            case 1:
            {
                WeakDGDirectionalHamiltonJacobi1D(LFfactors, InField, OutField, time);
            }
            break;

            case 2:
            {
                // WeakDGDirectionalHamiltonJacobi2D(alpha, beta, InField, OutField, time);
            }
            break;

            case 3:
            {
                // WeakDGHamiltonJacobi3D(alpha, beta, InField, OutField, time);
            }
            break;

            default:
            break;
        }
    }

    void MMFHamiltonJacobi::WeakDGHamiltonJacobi(
        const Array<OneD, const NekDouble> &LFfactors, 
        const Array<OneD, const Array<OneD, NekDouble>> &InField,
        Array<OneD, Array<OneD, NekDouble>> &OutField, const NekDouble time)
    {

        switch(m_expdim)
        {
            case 1:
            {
                // WeakDGHamiltonJacobi1D(LFfactors, InField, OutField, time);
            }
            break;

            case 2:
            {
                WeakDGHamiltonJacobi2D(LFfactors, InField, OutField, time);
            }
            break;

            case 3:
            {
                // WeakDGHamiltonJacobi3D(alpha, beta, InField, OutField, time);
            }
            break;

            default:
            break;
        }
    }

    /**
     * @brief Calculate weak DG advection in the form \f$ \langle\phi,
     * \hat{F}\cdot n\rangle - (\nabla \phi \cdot F) \f$
     *
     * @param   InField         Fields.
     * @param   OutField        Storage for result.
     * @param   NumericalFluxIncludesNormal     Default: true.
     * @param   InFieldIsPhysSpace              Default: false.
     * @param   nvariables      Number of fields.
    //  */
    // void MMFHamiltonJacobi::WeakDGHamiltonJacobi1D(
    //     const Array<OneD, const NekDouble> &LFfactors, const Array<OneD, const Array<OneD, NekDouble>> &InField,
    //     Array<OneD, Array<OneD, NekDouble>> &OutField, const NekDouble time)
    // {
    //     boost::ignore_unused(time);

    //     int nq              = GetNpoints();
    //     int ncoeffs         = GetNcoeffs();
    //     int nTracePointsTot = GetTraceNpoints();
    //     int nvar    = InField.size();
    //     int m_pmdim = 2;

    //     Array<OneD, Array<OneD, NekDouble>> physfield(nvar);
    //     for (int i = 0; i < nvar; ++i)
    //     {
    //         physfield[i] = InField[i];
    //     }

    //     Array<OneD, Array<OneD, NekDouble>> HJinput(m_expdim);
    //     Array<OneD, Array<OneD, NekDouble>> HJflux(m_expdim*m_pmdim);
    //     Array<OneD, Array<OneD, NekDouble>> HJfluxc(m_expdim*m_pmdim);
    //     for (int j = 0; j < m_expdim; ++j)
    //     {
    //         for (int pm = 0; pm < m_pmdim; ++pm)
    //         {
    //             HJflux[2*j+pm] = Array<OneD, NekDouble>(nq);
    //             HJfluxc[2*j+pm] = Array<OneD, NekDouble>(ncoeffs);

    //         }
    //     }

    //     Array<OneD, NekDouble> Fwd(nTracePointsTot);
    //     Array<OneD, NekDouble> Bwd(nTracePointsTot);
    //     Array<OneD, NekDouble> Fwdtmp(nTracePointsTot);
    //     Array<OneD, NekDouble> Bwdtmp(nTracePointsTot);
    //     Array<OneD, NekDouble> pqflux(nTracePointsTot);

    //     Array<OneD, NekDouble> tmp(nq);
    //     Array<OneD, NekDouble> tmpc(ncoeffs);

    //     Array<OneD, Array<OneD, NekDouble>> ncdotMFFwd(m_mfdim);
    //     Array<OneD, Array<OneD, NekDouble>> ncdotMFBwd(m_mfdim);
    //     for (int j = 0; j < m_mfdim; ++j)
    //     {
    //         ncdotMFFwd[j] = Array<OneD, NekDouble>(nTracePointsTot, 0.0);
    //         ncdotMFBwd[j] = Array<OneD, NekDouble>(nTracePointsTot, 0.0);
    //     }

    //     Vmath::Fill(nTracePointsTot, 1.0, &ncdotMFFwd[0][0], 1);
    //     Vmath::Fill(nTracePointsTot, -1.0, &ncdotMFBwd[0][0], 1);

    //     // p1 = \int \nabla \phi \cdot \mathbf{e}_1
    //     // p2 = \int \nabla \phi \cdot \mathbf{e}_1
    //     // q1 = \int \nabla \phi \cdot \mathbf{e}_2
    //     // q2 = \int \nabla \phi \cdot \mathbf{e}_2
    //     Array<OneD, NekDouble> diffFB(nTracePointsTot);
    //     for (int i=0; i<nvar; ++i)
    //     {
    //         // Fwd = u^+, Bwd = u^-
    //         m_fields[i]->GetFwdBwdTracePhys(physfield[i], Fwd, Bwd);
    //         for (int j = 0; j < m_expdim; ++j)
    //         {
    //             HJinput[j] = Array<OneD, NekDouble>(nq, 0.0);
    //             for (int pm = 0; pm < m_pmdim; ++pm)
    //             {
    //                 // Directional derivation with respect to the j'th moving frame
    //                 // tmp_j = ( \nabla \phi, fluxvector[j] \mathbf{e}^j )
    //                 // m_fields[i]->IProductWRTDirectionalDerivBase(m_movingframes[j], physfield[i], tmpc);
    //                 m_fields[i]->IProductWRTDerivBase(0, physfield[i], tmpc);

    //                 if (pm==1)
    //                 {
    //                     m_fields[i]->GetTrace()->Upwind(ncdotMFFwd[j], Fwd, Bwd, pqflux);
    //                 }

    //                 else if(pm==0)
    //                 {
    //                     m_fields[i]->GetTrace()->Upwind(ncdotMFBwd[j], Bwd, Fwd, pqflux);
    //                 }

    //                 Vmath::Vadd(nTracePointsTot, Fwd, 1, Bwd, 1, pqflux, 1);
    //                 Vmath::Smul(nTracePointsTot, 0.5, pqflux, 1, pqflux, 1);

    //                 Vmath::Vmul(nTracePointsTot, &pqflux[0], 1, &ncdotMFFwd[j][0], 1, &Fwdtmp[0], 1);
    //                 // Vmath::Vmul(nTracePointsTot, &pqflux[0], 1, &ncdotMFBwd[j][0], 1, &Bwdtmp[0], 1);
    //                 Vmath::Neg(ncoeffs, tmpc, 1);
    //                 m_fields[i]->AddTraceIntegral(Fwdtmp, tmpc);     

    //                 m_fields[i]->MultiplyByElmtInvMass(tmpc, tmpc);
    //                 Vmath::Vcopy(ncoeffs, tmpc, 1, HJfluxc[2*j+pm], 1);
    //                 m_fields[i]->BwdTrans(tmpc, HJflux[2*j+pm]);

    //                 Vmath::Vadd(nq, &HJflux[2*j+pm][0], 1, &HJinput[j][0], 1, &HJinput[j][0], 1);
    //             }

    //             Vmath::Smul(nq, 0.5, &HJinput[j][0], 1, &HJinput[j][0], 1);
    //         }

    //         OutField[i] = HamiltonianFunction(HJinput[0], HJinput[1]);

    //         // Add -0.5 \alpha (p1 - p2)
    //         Vmath::Vsub(nq, &HJflux[0][0], 1, &HJflux[1][0], 1, &tmp[0], 1);
    //         Vmath::Smul(nq, -0.5*LFfactors[0], &tmp[0], 1, &tmp[0], 1);
    //         Vmath::Vadd(nq, tmp, 1, OutField[i], 1, OutField[i], 1);
    //     }
    // }


    void MMFHamiltonJacobi::WeakDGDirectionalHamiltonJacobi1D(
        const Array<OneD, const NekDouble> &LFfactors, 
        const Array<OneD, const Array<OneD, NekDouble>> &InField,
        Array<OneD, Array<OneD, NekDouble>> &OutField, 
        const NekDouble time)
    {
        // boost::ignore_unused(time);

        int nq              = GetNpoints();
        int ncoeffs         = GetNcoeffs();
        int nTracePointsTot = GetTraceNpoints();
        int nvar    = InField.size();
        int m_pmdim = 2;
        
        Array<OneD, Array<OneD, NekDouble>> physfield(nvar);
        for (int i = 0; i < nvar; ++i)
        {
            physfield[i] = InField[i];
        }

        Array<OneD, Array<OneD, NekDouble>> HJinput(m_expdim);
        Array<OneD, Array<OneD, NekDouble>> HJflux(m_expdim);
        for (int j = 0; j < m_expdim; ++j)
        {
            HJflux[j] = Array<OneD, NekDouble>(nq);
        }

        Array<OneD, NekDouble> Fwd(nTracePointsTot);
        Array<OneD, NekDouble> Bwd(nTracePointsTot);
        Array<OneD, NekDouble> Fwdtmp(nTracePointsTot);
        Array<OneD, NekDouble> Bwdtmp(nTracePointsTot);
        Array<OneD, NekDouble> pqflux(nTracePointsTot);

        Array<OneD, NekDouble> tmp(nq);
        Array<OneD, NekDouble> tmpc(ncoeffs);

        Array<OneD, NekDouble> flux(nq);
        Array<OneD, NekDouble> fluxc(ncoeffs);

        // p1 = \int \nabla \phi \cdot \mathbf{e}_1
        // p2 = \int \nabla \phi \cdot \mathbf{e}_1
        // q1 = \int \nabla \phi \cdot \mathbf{e}_2
        // q2 = \int \nabla \phi \cdot \mathbf{e}_2
        NekDouble csign;
        for (int i=0; i<nvar; ++i)
        {
            // Fwd = u^+, Bwd = u^-
            m_fields[i]->GetFwdBwdTracePhys(physfield[i], Fwd, Bwd);

            for (int j = 0; j < m_expdim; ++j)
            {
                // Directional derivation with respect to the j'th moving frame
                // tmp_j = ( \nabla \phi, fluxvector[j] \mathbf{e}^j )
                m_fields[i]->IProductWRTDirectionalDerivBase(m_movingframes[j], physfield[i], tmpc);
                // m_fields[i]->IProductWRTDerivBase(1, physfield[i], tmpc);
                Vmath::Neg(ncoeffs, tmpc, 1);
                                    
                m_fields[i]->MultiplyByElmtInvMass(tmpc,fluxc);
                m_fields[i]->BwdTrans(fluxc,flux);

                HJinput[j] = Array<OneD, NekDouble>(nq, 0.0);
                HJflux[j] = Array<OneD, NekDouble>(nq, 0.0);
                for (int pm = 0; pm < m_pmdim; ++pm)
                {
                    Vmath::Vcopy(ncoeffs, tmpc, 1, fluxc, 1);

                    if(pm==0)
                    {
                       csign = 1.0;
                        m_fields[i]->AddTraceIntegral(Fwd, fluxc);     
                      // Vmath::Vmul(nTracePointsTot, &Fwd[0], 1, &m_ncdotMFFwd[j][0], 1, &Fwdtmp[0], 1);
                     //  Vmath::Vmul(nTracePointsTot, &Bwd[0], 1, &m_ncdotMFBwd[j][0], 1, &Bwdtmp[0], 1);
                    }

                    else if(pm==1)
                    {
                       csign=-1.0;
                        m_fields[i]->AddTraceIntegral(Bwd, fluxc);     
                     //  Vmath::Vmul(nTracePointsTot, &Bwd[0], 1, &m_ncdotMFFwd[j][0], 1, &Fwdtmp[0], 1);
                     //  Vmath::Vmul(nTracePointsTot, &Fwd[0], 1, &m_ncdotMFBwd[j][0], 1, &Bwdtmp[0], 1);
                    }

                  //   Vmath::Vadd(nTracePointsTot, Fwd, 1, Bwd, 1, pqflux, 1);
                  //  Vmath::Smul(nTracePointsTot, 0.5, pqflux, 1, pqflux, 1);

                    //m_fields[i]->AddFwdBwdTraceIntegral(Fwdtmp, Bwdtmp, fluxc);
                 //   m_fields[i]->AddTraceIntegral(Fwdtmp, fluxc);     
                    m_fields[i]->SetPhysState(false);
                    m_fields[i]->MultiplyByElmtInvMass(fluxc,fluxc);
                    m_fields[i]->BwdTrans(fluxc,flux);

                    // Compute (p1+p2)/2
                    Vmath::Svtvp(nq, 0.5, &flux[0], 1, &HJinput[j][0], 1, &HJinput[j][0], 1);

                    // Compute (p1-p2)
                    Vmath::Svtvp(nq, csign, &flux[0], 1, &HJflux[j][0], 1, &HJflux[j][0], 1);
                }
            }
            Array<OneD, NekDouble> exactp(nq);
            Array<OneD, NekDouble> exactperr(nq);
            exactp = Evaluatepq(time);

            // Forcible boudary condtions
            // HJinput[0][0] = exactp[0];
            // HJinput[0][nq-1] = exactp[nq-1];

            OutField[i] = HamiltonianFunction(exactp, exactp);
            // OutField[i] = HamiltonianFunction(exactp, HJinput[1]);

            // Add -0.5 \alpha (p1 - p2)
            for (int j=0; j<m_expdim; ++j)
            {
                Vmath::Svtvp(nq, LFfactors[j], &HJflux[j][0], 1, &OutField[i][0], 1, &OutField[i][0], 1);
            }
        }
    }


    /**
     * @brief Calculate weak DG advection in the form \f$ \langle\phi,
     * \hat{F}\cdot n\rangle - (\nabla \phi \cdot F) \f$
     *
     * @param   InField         Fields.
     * @param   OutField        Storage for result.
     * @param   NumericalFluxIncludesNormal     Default: true.
     * @param   InFieldIsPhysSpace              Default: false.
     * @param   nvariables      Number of fields.
     */
    void MMFHamiltonJacobi::WeakDGHamiltonJacobi2D(
        const Array<OneD, const NekDouble> &LFfactors, 
        const Array<OneD, const Array<OneD, NekDouble>> &InField,
        Array<OneD, Array<OneD, NekDouble>> &OutField, const NekDouble time)
    {
        boost::ignore_unused(time);

        int nq              = GetNpoints();
        int ncoeffs         = GetNcoeffs();
        int nTracePointsTot = GetTraceNpoints();
        int nvar    = InField.size();
        int m_pmdim = 2;
        
        Array<OneD, Array<OneD, NekDouble>> physfield(nvar);
        for (int i = 0; i < nvar; ++i)
        {
            physfield[i] = InField[i];
        }

        Array<OneD, Array<OneD, NekDouble>> HJinput(m_expdim);
        Array<OneD, Array<OneD, NekDouble>> HJflux(m_expdim*m_pmdim);
        for (int j = 0; j < m_expdim; ++j)
        {
            HJinput[j] = Array<OneD, NekDouble>(nq, 0.0);
            for (int pm = 0; pm < m_pmdim; ++pm)
            {
                HJflux[2*j+pm] = Array<OneD, NekDouble>(nq);
            }
        }

        Array<OneD, NekDouble> Fwd(nTracePointsTot);
        Array<OneD, NekDouble> Bwd(nTracePointsTot);
        Array<OneD, NekDouble> Fwdtmp(nTracePointsTot);
        Array<OneD, NekDouble> Bwdtmp(nTracePointsTot);
        Array<OneD, NekDouble> pqflux(nTracePointsTot);

        Array<OneD, NekDouble> tmp(nq);
        Array<OneD, NekDouble> tmpc(ncoeffs);

        // p1 = \int \nabla \phi \cdot \mathbf{e}_1
        // p2 = \int \nabla \phi \cdot \mathbf{e}_1
        // q1 = \int \nabla \phi \cdot \mathbf{e}_2
        // q2 = \int \nabla \phi \cdot \mathbf{e}_2
        for (int i=0; i<nvar; ++i)
        {
            // Fwd = u^+, Bwd = u^-
            m_fields[i]->GetFwdBwdTracePhys(physfield[i], Fwd, Bwd);

            for (int j = 0; j < m_expdim; ++j)
            {
                for (int pm = 0; pm < m_pmdim; ++pm)
                {
                    // Directional derivation with respect to the j'th moving frame
                    // tmp_j = ( \nabla \phi, fluxvector[j] \mathbf{e}^j )
                    m_fields[i]->IProductWRTDirectionalDerivBase(m_movingframes[j], physfield[i], tmpc);

                    if (pm==0)
                    {
                        m_fields[i]->GetTrace()->Upwind(m_ncdotMFFwd[j], Fwd, Bwd, pqflux);
                    }

                    else if(pm==1)
                    {
                        m_fields[i]->GetTrace()->Upwind(m_ncdotMFBwd[j], Bwd, Fwd, pqflux);
                    }

                    Vmath::Vmul(nTracePointsTot, &pqflux[0], 1, &m_ncdotMFFwd[j][0], 1, &Fwdtmp[0], 1);
                    Vmath::Vmul(nTracePointsTot, &pqflux[0], 1, &m_ncdotMFBwd[j][0], 1, &Bwdtmp[0], 1);

                    Vmath::Neg(ncoeffs, tmpc, 1);
                    m_fields[i]->AddFwdBwdTraceIntegral(Fwdtmp, Bwdtmp, tmpc);

                    m_fields[i]->MultiplyByElmtInvMass(tmpc, tmpc);
                    m_fields[i]->BwdTrans(tmpc, tmp);

                    Vmath::Vcopy(nq, &tmp[0], 1, &HJflux[2*j+pm][0], 1);
                    Vmath::Vadd(nq, &tmp[0], 1, &HJinput[j][0], 1, &HJinput[j][0], 1);
                }

                Vmath::Smul(nq, 0.5, &HJinput[j][0], 1, &HJinput[j][0], 1);
            }

            Array<OneD, NekDouble> x0(nq);
            Array<OneD, NekDouble> x1(nq);
            Array<OneD, NekDouble> x2(nq);

            m_fields[0]->GetCoords(x0, x1, x2);

            Array<OneD, NekDouble> Exactp(nq);
            for (int i=0; i<nq; ++i)
            {
                // tmp = exp(-1.0*time)*tan(x0[i]/2.0);
                // Exactp[i] = -1.0*cos(2.0*atan(tmp));
                Exactp[i] = cos(x0[i]);
            }
            Vmath::Vsub(nq, HJinput[0], 1, Exactp, 1, Exactp, 1);

            // Compute \int \hat{H} (p_1,p_2,q_1,q_2) dx
            tmp = HamiltonianFunction(HJinput[0], HJinput[1]);
            Vmath::Vcopy(nq, &tmp[0], 1, &OutField[i][0], 1);

            // Add -0.5 \alpha (p1 - p2)
            Vmath::Vsub(nq, &HJflux[0][0], 1, &HJflux[1][0], 1, &tmp[0], 1);
            Vmath::Smul(nq, -0.5*LFfactors[0], &tmp[0], 1, &tmp[0], 1);
            Vmath::Vadd(nq, tmp, 1, OutField[i], 1, OutField[i], 1);
            
            // Add -0.5 \beta (q1 - q2)
            Vmath::Vsub(nq, &HJflux[2][0], 1, &HJflux[3][0], 1, &tmp[0], 1);
            Vmath::Smul(nq, -0.5*LFfactors[1], &tmp[0], 1, &tmp[0], 1);
            Vmath::Vadd(nq, tmp, 1, OutField[i], 1, OutField[i], 1);
        }
    }


    Array<OneD, NekDouble> MMFHamiltonJacobi::HamiltonianFunction(
        const Array<OneD, const NekDouble> &parray, 
        const Array<OneD, const NekDouble> &qarray)
    {
        int nq              = GetNpoints();

        Array<OneD, NekDouble> outarray(nq); 

        Array<OneD, NekDouble> x0(nq);
        Array<OneD, NekDouble> x1(nq);
        Array<OneD, NekDouble> x2(nq);

        m_fields[0]->GetCoords(x0, x1, x2);

        Array<OneD, NekDouble> tmp(nq);
        switch(m_HJtype)
        {
            case eLinear1Dx:
            {
                for (int i=0; i<nq; ++i)
                {
                    outarray[i] = sin(x0[i])*parray[i];
                }
            }
            break;

            case eLinear1Dy:
            {
                for (int i=0; i<nq; ++i)
                {
                    outarray[i] = sin(x1[i])*parray[i];
                }
            }
            break;

            case eLinearNS1Dy:
            {
                int cnt=0;
                NekDouble signcos;
                NekDouble Tol=0.000001;
                for (int i=0; i<nq; ++i)
                {
                   if( cos(x1[i])>Tol )
                   {
                       signcos = 1.0;
                   }

                   else if ( cos(x1[i])<-1.0*Tol )
                   {
                       signcos = -1.0;
                   }

                   else 
                   {
                       signcos = 0.0;
                       cnt++;
                   }

                    outarray[i] = signcos*parray[i];
                }
                //std::cout << "zero signcos = " << cnt << std::endl;
            }
            break;

            case eBurgers1D:
            case eBurgersNS1D:
            {
                Vmath::Vmul(nq, parray, 1, parray, 1, tmp, 1);
                Vmath::Smul(nq, 0.5, tmp, 1, outarray, 1);
            }
            break;

            case eEikonal1D:
            {
                Vmath::Vabs(nq, parray, 1, outarray, 1);
            }
            break;

            case eBurgers:
            {
                Array<OneD, NekDouble> ones(nq, 1.0);

                Vmath::Vadd(nq, parray, 1, qarray, 1, tmp, 1);
                Vmath::Vadd(nq, ones, 1, tmp, 1, tmp, 1);
                Vmath::Vmul(nq, tmp, 1, tmp, 1, tmp, 1);
                Vmath::Smul(nq, 0.5, tmp, 1, outarray, 1);
            }
            break;

            case eNonconvex:
            {
                Array<OneD, NekDouble> ones(nq, 1.0);

                Vmath::Vadd(nq, parray, 1, qarray, 1, tmp, 1);
                Vmath::Vadd(nq, ones, 1, tmp, 1, tmp, 1);
                for (int i=0; i<nq; ++i)
                {
                    outarray[i] = -1.0*cos(tmp[i]);
                }
            }
            break;

            case eRiemann:
            {
                Vmath::Vadd(nq, parray, 1, qarray, 1, tmp, 1);
                for (int i=0; i<nq; ++i)
                {
                    outarray[i] = sin(tmp[i]);
                }
            }
            break;

            default:
            break;
        }

        return outarray;
    }

    void MMFHamiltonJacobi::v_SetInitialConditions(const NekDouble initialtime,
                                            bool dumpInitialConditions,
                                            const int domain)
    {
        int nq              = GetNpoints();

        boost::ignore_unused(domain);
        boost::ignore_unused(initialtime);

        // int nq   = GetTotPoints();
        int nvar = m_fields.size();

        Array<OneD, NekDouble> x0(nq);
        Array<OneD, NekDouble> x1(nq);
        Array<OneD, NekDouble> x2(nq);

        m_fields[0]->GetCoords(x0, x1, x2);

        Array<OneD, NekDouble> inittmp(nq);
        switch(m_HJtype)
        {
            case eLinear1Dx:
            {
                for (int i=0; i<nq; ++i)
                {
                   inittmp[i] = sin(x0[i]);
                }
            }
            break;

            case eLinear1Dy:
            case eLinearNS1Dy:
            case eEikonal1D:
            case eBurgers1D:
            {
                for (int i=0; i<nq; ++i)
                {
                   inittmp[i] = sin(x1[i]);
                }
            }
            break;

            case eBurgersNS1D:
            {
                Array<OneD, NekDouble> x0(nq);
                Array<OneD, NekDouble> x1(nq);
                Array<OneD, NekDouble> x2(nq);

                m_fields[0]->GetCoords(x0, x1, x2);

                NekDouble xp;
                for (int i=0; i<nq; ++i)
                {
                   xp = x1[i];
                   if(xp<=m_pi)
                   {
                       inittmp[i] = m_pi - xp;
                   }

                   else
                   {
                        inittmp[i] = xp - m_pi;
                   }
                }
            }
            break;

            default:
             break;
        }

        m_fields[0]->SetPhys(inittmp);
        
        // forward transform to fill the modal coeffs
        for (int i = 0; i < nvar; ++i)
        {
            m_fields[i]->SetPhysState(true);
            m_fields[i]->FwdTrans(m_fields[i]->GetPhys(), m_fields[i]->UpdateCoeffs());
        }

        if (dumpInitialConditions)
        {
            std::cout << "inittmp = " << RootMeanSquare(inittmp) << std::endl;
            std::string outname = m_sessionName + "_initial.chk";
            WriteFld(outname);
        }
    }

    void MMFHamiltonJacobi::v_EvaluateExactSolution(unsigned int field,
                                            Array<OneD, NekDouble> &outfield,
                                            const NekDouble time)
    {
        boost::ignore_unused(field);

        int nq              = GetNpoints();
        outfield = Array<OneD, NekDouble>(nq,0.0);
        
        Array<OneD, NekDouble> x0(nq);
        Array<OneD, NekDouble> x1(nq);
        Array<OneD, NekDouble> x2(nq);

        m_fields[0]->GetCoords(x0, x1, x2);

        NekDouble tmp;
        outfield = Array<OneD, NekDouble>(nq, 0.0);
        switch(m_HJtype)
        {
            case eLinear1Dx:
            {
                for (int i=0; i<nq; ++i)
                {
                    tmp = exp(-1.0*time)*tan(x0[i]/2.0);
                    outfield[i] = sin(2.0*atan(tmp));
                }
            }
            break;

            case eLinear1Dy:
            {
                for (int i=0; i<nq; ++i)
                {
                    tmp = exp(-1.0*time)*tan(x1[i]/2.0);
                    outfield[i] = sin(2.0*atan(tmp));
                }
            }
            break;

            case eLinearNS1Dy:
            case eEikonal1D:
            {
                Array<OneD, NekDouble> x0(nq);
                Array<OneD, NekDouble> x1(nq);
                Array<OneD, NekDouble> x2(nq);

                m_fields[0]->GetCoords(x0, x1, x2);

                NekDouble xp;
                for (int i=0; i<nq; ++i)
                {
                    xp = x1[i];
                    if(time<=0.5*m_pi)
                    {
                        if(xp<=0.5*m_pi)
                        {
                           outfield[i] = sin(xp-time);
                        }

                        else if (xp<=1.5*m_pi-time)
                        {
                           outfield[i] = sin(xp+time);
                        }

                        else if (xp<=1.5*m_pi+time)
                        {
                           outfield[i] = -1.0;
                        }

                        else
                        {
                            outfield[i] = sin(xp-time);
                        }
                    }

                    else if(time<=m_pi)
                    {
                        if(xp<=time-0.5*m_pi)
                        {
                           outfield[i] = -1.0;
                        }

                        else if (xp<=0.5*m_pi)
                        {
                           outfield[i] = sin(xp-time);
                        }

                        else if (xp<=1.5*m_pi-time)
                        {
                           outfield[i] = sin(xp+time);
                        }

                        else
                        {
                            outfield[i] = -1.0;
                        }
                    }

                    else
                    {
                        outfield[i] = -1.0;
                    }
                }
            }
            break;

            case eBurgers1D:
            {
                Array<OneD, NekDouble> x0(nq);
                Array<OneD, NekDouble> x1(nq);
                Array<OneD, NekDouble> x2(nq);

                m_fields[0]->GetCoords(x0, x1, x2);

                NekDouble xp, xpts;
                for (int i=0; i<nq; ++i)
                {
                    xp = x1[i];
                    xpts = BurgersFindx0(xp, time);
                    outfield[i] = sin(xpts) + 0.5*cos(xpts)*cos(xpts)*time;
                }
            }
            break;

            default:
             break;
        }
    }

    // Compute the exact p and q
   Array<OneD, NekDouble> MMFHamiltonJacobi::Evaluatepq(const NekDouble time)
    {
        int nq              = GetNpoints();
        
        Array<OneD, NekDouble> x0(nq);
        Array<OneD, NekDouble> x1(nq);
        Array<OneD, NekDouble> x2(nq);

        m_fields[0]->GetCoords(x0, x1, x2);

        NekDouble tmp;

        Array<OneD, NekDouble> outfield(nq);
        switch(m_HJtype)
        {
            case eLinear1Dx:
            {
                for (int i=0; i<nq; ++i)
                {
                    tmp = exp(-1.0*time)*tan(x0[i]/2.0);
                    outfield[i] = cos(2.0*atan(tmp));
                }
            }
            break;

            case eLinear1Dy:
            {
                for (int i=0; i<nq; ++i)
                {
                    tmp = exp(-1.0*time)*tan(x1[i]/2.0);
                    outfield[i] = cos(2.0*atan(tmp));
                }
            }
            break;

            case eLinearNS1Dy:
            case eEikonal1D:
            {
                Array<OneD, NekDouble> x0(nq);
                Array<OneD, NekDouble> x1(nq);
                Array<OneD, NekDouble> x2(nq);

                m_fields[0]->GetCoords(x0, x1, x2);

                NekDouble xp;
                for (int i=0; i<nq; ++i)
                {
                    xp = x1[i];
                    if(time<=0.5*m_pi)
                    {
                        if(xp<=0.5*m_pi)
                        {
                           outfield[i] = cos(xp-time);
                        }

                        else if (xp<=1.5*m_pi-time)
                        {
                           outfield[i] = cos(xp+time);
                        }

                        else if (xp<=1.5*m_pi+time)
                        {
                           outfield[i] = 0.0;
                        }

                        else
                        {
                            outfield[i] = cos(xp-time);
                        }
                    }

                    else if(time<=m_pi)
                    {
                        if(xp<=time-0.5*m_pi)
                        {
                           outfield[i] = 0.0;
                        }

                        else if (xp<=0.5*m_pi)
                        {
                           outfield[i] = cos(xp-time);
                        }

                        else if (xp<=1.5*m_pi-time)
                        {
                           outfield[i] = cos(xp+time);
                        }

                        else
                        {
                            outfield[i] = 0.0;
                        }
                    }

                    else
                    {
                        outfield[i] = 0.0;
                    }
                }
            }
            break;

            case eBurgers1D:
            {
                Array<OneD, NekDouble> x0(nq);
                Array<OneD, NekDouble> x1(nq);
                Array<OneD, NekDouble> x2(nq);

                m_fields[0]->GetCoords(x0, x1, x2);

                NekDouble xp, xpts;
                for (int i=0; i<nq; ++i)
                {
                    xp = x1[i];
                    xpts = BurgersFindx0(xp, time);
                    outfield[i] = cos(xpts) - cos(xpts)*sin(xpts)*time;
                }
            }
            break;

            default:
             break;
        }

        return outfield;
    }

    // F = x - x0 - cos(x0) * time
    NekDouble MMFHamiltonJacobi::BurgersFindx0(const NekDouble xp, const NekDouble time)
    {
        int i;
        int Maxiter = 100000;
        NekDouble Tol = 0.0001;

        NekDouble xpts, xptsnew;
        xpts = 0.0;
        for ( i=0;i<Maxiter; ++i)
        {
            xptsnew = xpts - (xp - xpts - cos(xpts)*time)/( -1.0 + sin(xpts)*time );

            if(fabs(xptsnew-xpts)<Tol)
            {
                break;
            }
            
            xpts = xptsnew;
        }

        return xpts;
    }

    void MMFHamiltonJacobi::v_GenerateSummary(SolverUtils::SummaryList &s)
    {
        MMFSystem::v_GenerateSummary(s);
        SolverUtils::AddSummaryItem(s, "HamiltonJacobiType", HamiltonJacobiTypeMap[m_HJtype]);
        SolverUtils::AddSummaryItem(s, "LF alpha", m_LFfactors[0]);
        SolverUtils::AddSummaryItem(s, "LF beta", m_LFfactors[1]);
    }
}
