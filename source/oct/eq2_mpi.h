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
  void MySolver<dim,no_time_steps>::rt_propagtion_backward( const int ex )
  {
    m_computing_timer.enter_section(__func__);

    m_workspace_ng = m_Psi_d;
    m_workspace_ng -= m_all_Psi[no_time_steps-1];

    constraints.distribute(m_workspace_ng);
    m_Psi = m_workspace_ng;

    MyComplexTools::MPI::AssembleSystem_mulvz( dof_handler, fe, constraints, m_Psi, std::complex<double>(0,-1), system_matrix, system_rhs );
    m_sol=0;
    solve_eq1();
    m_all_p[no_time_steps-1] = m_sol;
    double N = MyComplexTools::MPI::Particle_Number( mpi_communicator, dof_handler, fe, m_Psi );
    
    for( int i=no_time_steps-2; i>0; i-- )
    {
      MyComplexTools::MPI::AssembleSystem_LIN_Step( dof_handler, fe, constraints, m_Psi, -m_dt, system_matrix, system_rhs );
      solve_eq1();
      m_potential.set_time( double(i)*m_dt );
      MyComplexTools::MPI::AssembleSystem_NL_Step( dof_handler, fe, constraints, m_Psi, m_potential, -m_dt, 2*m_gs,  system_matrix, system_rhs );
      solve_cg();
      
      m_all_p[i] = m_sol;      

//      output_vec( "p_" + to_string(i) + ".vtu", m_all_p[i] );
//      double N = MyComplexTools::MPI::Particle_Number( mpi_communicator, dof_handler, fe, m_Psi );
      if(m_root) printf( "b: %g\n", double(i)*m_dt );
    }
    m_computing_timer.exit_section();
  }
