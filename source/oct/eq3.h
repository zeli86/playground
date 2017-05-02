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

  template <int no_time_steps, int no_lam>
  void MySolver<no_time_steps,no_lam>::compute_all_lambdas_tt()
  {
    const double fak = 1/(m_dt*m_dt);
    
    for( int s=0; s<no_lam; s++ )
      for( int i=1; i<no_time_steps-1; i++ )
        m_all_lambdas_grad[s][i] = (m_all_lambdas[s][i-1]+m_all_lambdas[s][i+1]-2*m_all_lambdas[s][i])*fak;
      
    for( int s=0; s<no_lam; s++ )
    {
      m_all_lambdas_grad[s][0] = (m_all_lambdas[s][0]-2*m_all_lambdas[s][1]+m_all_lambdas[s][2])*fak;
      m_all_lambdas_grad[s][no_time_steps-1] = (m_all_lambdas[s][no_time_steps-1]-2*m_all_lambdas[s][no_time_steps-2]+m_all_lambdas[s][no_time_steps-3])*fak;
    } 
  }
  
  template <int no_time_steps, int no_lam>
  void MySolver<no_time_steps,no_lam>::compute_correction()
  {
    compute_all_lambdas_tt();
    
    const QGauss<1> quadrature_formula(fe.degree+1);
    FEValues<1> fe_values (fe, quadrature_formula, update_values|update_quadrature_points|update_JxW_values);

    const unsigned int n_q_points = quadrature_formula.size();
    vector<Vector<double>> Psi(n_q_points,Vector<double>(2));
    vector<Vector<double>> p(n_q_points,Vector<double>(2));
    
    double retval;

    CPotential<no_time_steps,no_lam> Potential ( m_all_lambdas, 0 );
    
    // loop over all lambdas
    for( int s=0; s<no_lam; s++ )
    {
      // loop over all time steps
      for( int ti=1; ti<no_time_steps; ti++ )
      {
        Potential.m_sel = 1+s;
        Potential.m_timeindex = ti;

        retval=0;
    
        m_workspace = m_all_Psi[ti];
        m_workspace_2 = m_all_p[ti];

        typename DoFHandler<1>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
        for (; cell!=endc; ++cell)
        {
          if( cell->is_locally_owned() )
          {
            fe_values.reinit (cell);
            fe_values.get_function_values( m_workspace, Psi );
            fe_values.get_function_values( m_workspace_2, p );
            
            for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
            {
              retval += fe_values.JxW(q_point)*Potential.value(fe_values.quadrature_point(q_point))*(Psi[q_point][0]*p[q_point][0]+Psi[q_point][1]*p[q_point][1]);
            }
          }
        }
        //m_all_lambdas_grad[s][ti] += retval;
        m_all_lambdas_grad[s][ti] = retval;
      }
    }    

    for( int s=0; s<no_lam; s++ )
    {
      m_all_lambdas_grad[s][0] = 0;
      m_all_lambdas_grad[s][no_time_steps-1] = 0;
    }
    
    output_lambdas_grad( "grad_l2.txt" );
    
    LAPACKFullMatrix<double> lap(no_time_steps,no_time_steps);
    for( int i=0; i<no_time_steps; i++ ) lap(i,i) = 2;
    for( int i=1; i<no_time_steps-1; i++ ) lap(i,i+1) = -1; // rechte Nebendiagonale 
    for( int i=1; i<no_time_steps-1; i++ ) lap(i,i-1) = -1; // linke Nebendiagonale
    
    lap.compute_lu_factorization();
    
    Vector<double> tmp_vec(no_time_steps);
    
    for( int s=0; s<no_lam; s++ )
    {
      for( int ti=1; ti<no_time_steps-1; ti++ )
      {
        tmp_vec[ti] = m_dt*m_dt*m_all_lambdas_grad[s][ti];
      }

      lap.apply_lu_factorization(tmp_vec,false);

      for( int ti=0; ti<no_time_steps; ti++ )
      {
        m_all_lambdas_grad[s][ti] = tmp_vec[ti];
      }
    }
    
    output_lambdas_grad( "grad_sob.txt" );
    
    for( int s=0; s<no_lam; s++ )
    {
      for( int ti=1; ti<no_time_steps-1; ti++ )
      {
        m_all_lambdas[s][ti] -= 0.1*m_all_lambdas_grad[s][ti];
      }
    }
    
    for( int s=0; s<no_lam; s++ )
    {
      double tmp=0;
      for( int ti=1; ti<no_time_steps-1; ti++ )
      {
        tmp += m_all_lambdas_grad[s][ti]*m_all_lambdas_grad[s][ti];
      }
      printf( "l2norm sob grad = %g\n", sqrt(fabs(m_dt)*tmp));
    }    
  }

