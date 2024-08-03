/////////////////////////////////////////////////////////////////////////////
//
// File MMFWMTract.cpp
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
// Description: MMF White Matter Tractography routines
//
///////////////////////////////////////////////////////////////////////////////

#include <MMFSolver/EquationSystems/MMFWMTract.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <stdlib.h>

#include <boost/core/ignore_unused.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string.hpp>
#include <iomanip>

#include <LibUtilities/BasicUtils/Timer.h>
#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>
#include <SolverUtils/MMFSystem.h>

#include <typeinfo>

#include <LibUtilities/BasicUtils/Timer.h>
#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>


namespace Nektar
{
    std::string MMFWMTract::className =
        SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
            "MMFWMTract", MMFWMTract::create, "MMFWMTract equation.");

    MMFWMTract::MMFWMTract(const LibUtilities::SessionReaderSharedPtr &pSession,
                        const SpatialDomains::MeshGraphSharedPtr& pGraph)
        : UnsteadySystem(pSession, pGraph), MMFSystem(pSession, pGraph)
    {
    }

    /**
     * @brief Initialisation object for the unsteady linear advection equation.
     */
    void MMFWMTract::v_InitObject(bool DeclareFields)
    {
        // Call to the initialisation object
        UnsteadySystem::v_InitObject(DeclareFields);

        int nq = GetTotPoints();

        // Add Rectangular PML
        Array<OneD, Array<OneD, NekDouble>> AniStrength(m_expdim);
        for (int j = 0; j < m_expdim; ++j)
        {
            AniStrength[j] = Array<OneD, NekDouble>(nq, 1.0);
        }
        MMFSystem::MMFInitObject(AniStrength);

        // Define TestMaxwellType
        if (m_session->DefinesSolverInfo("WMTractType"))
        {
            std::string WMTractTypeStr =
                m_session->GetSolverInfo("WMTractType");
            for (int i = 0; i < (int)SIZE_WMTractType; ++i)
            {
                if (WMTractTypeMap[i] == WMTractTypeStr)
                {
                    m_WMTractType = (WMTractType)i;
                    break;
                }
            }
        }

        else
        {
            m_WMTractType = (WMTractType)0;
        }
        
        // Read csv file =================================================================
        // csv file name should be the same as xml file.
        Array<OneD, Array<OneD, Array<OneD, NekDouble>>> DTIcsv;
        std::string csvfile = m_sessionName + ".csv";

        read_csv(csvfile, DTIcsv);
        
        // Assign a vector or value to each element =================================================================
        int Nelemtj = m_fields[0]->GetExpSize();
        int nptsj   = m_fields[0]->GetTotPoints(0);

        std::cout << "Nelemtj = " << Nelemtj << ", nptsj = " << nptsj << ", nq = " << nq << std::endl;
        
        Array<OneD, NekDouble> x0(nq);
        Array<OneD, NekDouble> x1(nq);
        Array<OneD, NekDouble> x2(nq);
        
        m_fields[0]->GetCoords(x0, x1, x2);

        int Ns = cbrt(Nelemtj);
        std::cout << "Num elemt = " << Ns << std::endl;
        for (int i = 0; i < Ns; ++i)
        {
            for (int j = 0; j < Ns; ++j)
            {
                for (int k = 0; k < Ns; ++k)
                {
                    std::cout << "DTIVec: ( " << i << " , " << j << " , " << k << " ) = " << DTIcsv[i][j][k] << std::endl;
                }
            }
        }

        // Each element assigned by (i,j,k) index 
        Array<OneD, NekDouble> DTIVector(nq);

        int index, elemid;
        for (int k = 0; k< Ns; ++k)
        {
            for (int j = 0; j < Ns; ++j)
            {
                for (int i = 0; i < Ns; ++i)
                {
                    // Option 1
                    elemid = i + Ns * j + Ns * Ns * k;

                    // Option 2
                    // elemid = j + Ns * i + Ns * Ns * k;

                    std::cout << "elemid = " << elemid << std::endl;

                    for (int pt = 0; pt < nptsj; ++pt)
                    {
                        index = m_fields[0]->GetPhys_Offset(elemid) + j;

                        DTIVector[index] = DTIcsv[i][j][k];
                    }
                }
            }
        }

        // for (int i=0; i<nq; ++i)
        // {
        //     std::cout << "DTI: i = " << i << ", val = " << DTIVector[i] << std::endl;
        // }

        PlotDTIVector(DTIVector);

        wait_on_enter();

        // If explicit it computes RHS and PROJECTION for the time integration
        m_ode.DefineOdeRhs(&MMFWMTract::DoOdeRhs, this);
        m_ode.DefineProjection(&MMFWMTract::DoOdeProjection, this);
    }

    /**
     * @brief Unsteady linear advection equation destructor.
     */
    MMFWMTract::~MMFWMTract()
    {
    }

    void MMFWMTract::PlotDTIVector(const Array<OneD, const NekDouble> &DTIVector)
    {
        int nvar    = 3;
        int ncoeffs = m_fields[0]->GetNcoeffs();

        std::string outname1;
        outname1 = m_sessionName + "_DTIVec.chk";

        std::vector<Array<OneD, NekDouble>> fieldcoeffs(nvar);
        for (int i = 0; i < nvar; ++i)
        {
            fieldcoeffs[i] = Array<OneD, NekDouble>(ncoeffs,0.0);
        }

        std::vector<std::string> variables(nvar);
        variables[0] = "DTIvecx";
        variables[1] = "DTIvecy";
        variables[2] = "DTIvecz";

        m_fields[0]->FwdTrans(DTIVector, fieldcoeffs[0]);

        WriteFld(outname1, m_fields[0], fieldcoeffs, variables);
    }

    void MMFWMTract::read_csv(std::string fname, Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &csv)
    {        
        int dim[3] = { 0, 0, 0 }; // size of the array
        int i = 0;
        std::fstream file(fname, std::ios::in);
        std::string line, word;
        std::vector<std::vector<int> > csv_indices;
        std::vector<int> csv_idx;
        std::vector<double> csv_values;

        if (file.is_open())
        {
            while(getline(file, line))
            {
                // For each line (row)
                i = 0;
                std::stringstream str(line);
                csv_idx.clear();
                while (getline(str, word, ','))
                {
                    // For each column
                    if (i < 3)
                    {
                        int d = std::stoi(word.c_str());
                        if (d + 1 > dim[i])
                            dim[i] = d + 1;
                        csv_idx.push_back(d);
                    } else if (i == 3) {
                        csv_values.push_back(std::strtod(word.c_str(), NULL));
                    }
                    i++;
                }
                csv_indices.push_back(csv_idx);
            }

            // Allocate memory for csv
            // double*** csv = (double***)malloc(sizeof(double**) * dim[0]);
            csv = Array<OneD, Array<OneD, Array<OneD, NekDouble>>>(dim[0]);
            // for (int x = 0; x < dim[0]; x++)
            // {
            //     csv[x] = (double**)malloc(sizeof(double**) * dim[1]);
            //     for (int y = 0; y < dim[1]; y++)
            //     {
            //         csv[x][y] = (double*)malloc(sizeof(double) * dim[2]);
            //     }
            // }
            for (int x = 0; x < dim[0]; x++)
            {
                // csv[x] = (double**)malloc(sizeof(double**) * dim[1]);
                csv[x] = Array<OneD, Array<OneD, NekDouble>>(dim[1]);
                for (int y = 0; y < dim[1]; y++)
                {
                    // csv[x][y] = (double*)malloc(sizeof(double) * dim[2]);
                    csv[x][y] = Array<OneD, NekDouble>(dim[2]);
                }
            }


            // Move vector to array
            int x, y, z;
            for (i = 0; i < csv_values.size(); i++)
            {
                x = csv_indices[i][0];
                y = csv_indices[i][1];
                z = csv_indices[i][2];
                csv[x][y][z] = csv_values[i];
            }


            // TEST: print array
            // for (int x = 0; x < dim[0]; x++)
            // {
            //     for (int y = 0; y < dim[1]; y++)
            //     {
            //         for (int z = 0; z < dim[2]; z++)
            //         {
            //             cout << x << " " << y << " " << z << " ";
            //             cout << csv[x][y][z] << "\n";
            //         }
            //     }
            // }
        }
    }


    void MMFWMTract::v_DoSolve()
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
    void MMFWMTract::DoOdeProjection(
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
    void MMFWMTract::DoOdeRhs(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time)
    {
        boost::ignore_unused(time);

        int nvar    = inarray.size();
        int nq      = GetTotPoints();

        for (int i = 0; i < nvar; ++i)
        {
            outarray[i] = Array<OneD, NekDouble>(nq, 0.0);
        }
    }

    void MMFWMTract::v_SetInitialConditions(const NekDouble initialtime,
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

        // forward transform to fill the modal coeffs
        for (int i = 0; i < nvar; ++i)
        {
            m_fields[i]->SetPhysState(true);
            m_fields[i]->FwdTrans(m_fields[i]->GetPhys(), m_fields[i]->UpdateCoeffs());
        }

        if (dumpInitialConditions)
        {
            std::string outname = m_sessionName + "_initial.chk";
            WriteFld(outname);
        }
    }

    void MMFWMTract::v_EvaluateExactSolution(unsigned int field,
                                            Array<OneD, NekDouble> &outfield,
                                            const NekDouble time)
    {
        boost::ignore_unused(time,field);

        int nq              = GetNpoints();
        outfield = Array<OneD, NekDouble>(nq,0.0);
    }


    void MMFWMTract::v_GenerateSummary(SolverUtils::SummaryList &s)
    {
        MMFSystem::v_GenerateSummary(s);
    }
}
