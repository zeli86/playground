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
  void MySolver<dim,no_time_steps>::compute_correction()
  {
    m_computing_timer.enter_section(__func__);
    
    const QGauss<dim>  quadrature_formula(fe.degree+1);
    FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_quadrature_points|update_JxW_values);

    const unsigned n_q_points = quadrature_formula.size();
    vector<Vector<double>> Psi(n_q_points,Vector<double>(2));
    vector<Vector<double>> p(n_q_points,Vector<double>(2));
    
    double retval, tmp1;

    /*
    // loop over all lambdas
    for( int s=0; s<no_lam; s++ )
    {
      // loop over all time steps
      for( int ti=1; ti<no_time_steps; ti++ )
      {
	      Potential.m_sel = 1+s;
	      Potential.m_timeindex = ti;

	      tmp1=0.0;
	
        m_workspace = m_all_Psi[ti];
        m_workspace_2 = m_all_p[ti];
	
        typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
        for (; cell!=endc; ++cell)
        {
          if( cell->is_locally_owned() )
          {
            fe_values.reinit (cell);
            fe_values.get_function_values( m_workspace, Psi );
            fe_values.get_function_values( m_workspace_2, p );
            
            for ( unsigned qp=0; qp<n_q_points; qp++ )
            {
              tmp1 += fe_values.JxW(qp) * Potential.value(fe_values.quadrature_point(qp)) * (Psi[qp][0]*p[qp][0]+Psi[qp][1]*p[qp][1]);
            }
          }
        }
        MPI_Allreduce( &tmp1, &retval, 1, MPI_DOUBLE, MPI_SUM, mpi_communicator);
        //m_all_lambdas_grad[s][ti] += retval;
        m_all_lambdas_grad[s][ti] = retval;

        printf( "%d %d %g\n", s, ti, retval );
      }
    }    
    */

    /*
    LAPACKFullMatrix<double> lap(no_time_steps,no_time_steps);
    for( int i=0; i<no_time_steps; i++ ) lap(i,i) = 2;
    for( int i=1; i<no_time_steps-1; i++ ) lap(i,i+1) = -1; // rechte Nebendiagonale 
    for( int i=1; i<no_time_steps-1; i++ ) lap(i,i-1) = -1; // linke Nebendiagonale
    */
    
//     ofstream bla( "lap.txt" ); 
//     lap.print_formatted( bla );

//    lap.compute_lu_factorization();
    

    m_computing_timer.exit_section();
  }

