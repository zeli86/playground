/* * atus-pro testing - atus-pro testing playgroung
 * Copyright (C) 2017 Želimir Marojević <zelimir.marojevic@gmail.com>
 *
 * This file is part of atus-pro testing.
 *
 * atus-pro testing is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * atus-pro testing is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with atus-pro testing.  If not, see <http://www.gnu.org/licenses/>.
 */

template <int dim, int no_time_steps>
  void MySolver<dim,no_time_steps>::rt_propagtion_forward()
  {
    m_computing_timer.enter_section(__func__);
    
    for( int i=1; i<no_time_steps; i++ )
    {
      m_timeindex = i-1;

      m_workspace = m_Psi;
      MPI::MyComplexTools::AssembleSystem_LIN_Step( dof_handler, fe, constraints, m_workspace, 0.5*m_dt, system_matrix, system_rhs );
      solve_eq1();
      m_workspace = m_sol;
      MPI::MyComplexTools::AssembleSystem_NL_Step( dof_handler, fe, constraints, m_workspace, m_potential, m_dt, m_gs,  system_matrix, system_rhs );
      solve_eq1();
      m_workspace = m_sol;
      MPI::MyComplexTools::AssembleSystem_LIN_Step( dof_handler, fe, constraints, m_workspace, 0.5*m_dt, system_matrix, system_rhs );
      solve_eq1();
      m_Psi = m_sol;

      m_all_Psi[i] = m_Psi;

      if(m_root) printf( "f: %g \n", double(i)*m_dt );
    }
   
    m_computing_timer.exit_section();
  }
      
  template <int dim, int no_time_steps>
  void MySolver<dim,no_time_steps>::solve_eq1 ()
  {
    m_computing_timer.enter_section(__func__);

    m_sol=0;
    PETScWrappers::PreconditionSOR preconditioner;
    PETScWrappers::PreconditionSOR::AdditionalData data;
    preconditioner.initialize(system_matrix, data);

    SolverControl solver_control (dof_handler.n_dofs(), 1e-15);
    PETScWrappers::SolverGMRES solver(solver_control, mpi_communicator);
    
    solver.solve(system_matrix, m_sol, system_rhs, preconditioner);
    constraints.distribute (m_sol);

/*
    SolverControl solver_control;
    PETScWrappers::SparseDirectMUMPS solver(solver_control, mpi_communicator);
    solver.set_symmetric_mode(false);
    solver.solve(system_matrix, m_sol, system_rhs);
    constraints.distribute (m_sol);
*/
    m_computing_timer.exit_section();
  }   