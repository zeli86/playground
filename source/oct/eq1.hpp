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
  void MySolver<dim,no_time_steps>::rt_propagtion_forward( const int ex )
  {
    m_Psi = m_all_Psi[0];
    for( int i=1; i<no_time_steps; i++ )
    {
      MyComplexTools::AssembleSystem_LIN_Step( dof_handler, fe, m_Psi, m_dt, system_matrix, system_rhs );
      solve();
      m_potential.set_time( double(i)*m_dt );
      MyComplexTools::AssembleSystem_NL_Step( dof_handler, fe, m_Psi, m_potential, m_dt, m_gs,  system_matrix, system_rhs );
      solve();

      //output_vec( "Psi_" + to_string(ex) + "_" + to_string(i) + ".vtu", m_sol );

      m_all_Psi[i] = m_Psi;

      //double N = MyComplexTools::Particle_Number( dof_handler, fe, m_Psi );
      //printf( "f: %g %g\n", double(i)*m_dt, N );
    }
  }
  
    
