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
  void MySolver<dim,no_time_steps>::rt_propagtion_backward( const int )
  {
    m_workspace = m_Psi_d;
    m_workspace -= m_all_Psi[no_time_steps-1];

    double cost = MyComplexTools::Particle_Number( dof_handler, fe, m_workspace );
    printf( "cost = %g\n", cost );   

    MyComplexTools::AssembleSystem_mulvz( dof_handler, fe, m_workspace, std::complex<double>(0,-1), system_matrix, system_rhs );
    solve();
    m_all_p[no_time_steps-1] = m_Psi;
    //output_vec( "p_" +  to_string(ex) + "_" + to_string(no_time_steps-1) + ".vtu", m_all_p[no_time_steps-1] );

    for( int i=no_time_steps-2; i>0; i-- )
    {
      MyComplexTools::AssembleSystem_LIN_Step( dof_handler, fe, m_Psi, -m_dt, system_matrix, system_rhs );
      solve();
      assemble_system(i);
      solve();
      m_all_p[i] = m_Psi;

      //output_vec( "p_" + to_string(i) + ".vtu", m_all_p[i] );
      
      //double N = MyComplexTools::Particle_Number( dof_handler, fe, m_Psi );
      //printf( "b: %g %g\n", double(i)*m_dt, N );
    }    
  }

  template <int dim, int no_time_steps>
  void MySolver<dim,no_time_steps>::assemble_system ( const int idx )
  {
    m_workspace = m_all_Psi[idx];

    const QGauss<dim> quadrature_formula(fe.degree+1);

    m_potential.set_time( double(idx)*m_dt );

    //pcout << double(idx)*m_dt << "\t" <<  double(idx+1)*m_dt << endl;
    
    const FEValuesExtractors::Scalar rt (0);
    const FEValuesExtractors::Scalar it (1);
        
    system_matrix = 0;
    system_rhs = 0;

    FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_gradients|update_quadrature_points|update_JxW_values);

    const unsigned dofs_per_cell = fe.dofs_per_cell;
    const unsigned n_q_points = quadrature_formula.size();

    FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);
    Vector<double> cell_rhs (dofs_per_cell);

    vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    vector<Vector<double>> Psi(n_q_points,Vector<double>(2));
    vector<Vector<double>> p(n_q_points,Vector<double>(2));
    vector<vector<Tensor<1,dim>>> p_grad(n_q_points, vector<Tensor<1,dim>>(2));
 
    const double dt = -m_dt; //reversed time direction
    
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
      cell_matrix = 0;
      cell_rhs = 0;

      fe_values.reinit (cell);
      fe_values.get_function_values(m_Psi, p);
      fe_values.get_function_gradients(m_Psi, p_grad);
      fe_values.get_function_values(m_workspace, Psi );

      for( unsigned qp=0; qp<n_q_points; qp++ )
      {
        double JxW = fe_values.JxW(qp);
        double ar = m_gs*(Psi[qp][0]*Psi[qp][0]-Psi[qp][1]*Psi[qp][1]);
        double ai = 2*m_gs*Psi[qp][0]*Psi[qp][1];
        double A = m_potential.value(fe_values.quadrature_point(qp)) + 2*m_gs*(Psi[qp][0]*Psi[qp][0]+Psi[qp][1]*Psi[qp][1]);
        double beta = sqrt(ar*ar+ai*ai+A*A);
        double f1 = cos(dt*beta);
        double f2 = sin(dt*beta)/beta;
        double M00 = f1 + f2*ai;
        double M01 = f2*(A-ar);
        double M10 = f2*(-A-ar);
        double M11 = f1 - f2*ai;
        
        for (unsigned i=0; i<dofs_per_cell; i++ )
        {
          for (unsigned j=0; j<dofs_per_cell; j++ )
          {
            cell_matrix(i,j) += JxW*(fe_values[rt].value(i,qp)*fe_values[rt].value(j,qp) + fe_values[it].value(i,qp)*fe_values[it].value(j,qp));                        
          }
          cell_rhs(i) += JxW*((M00*p[qp][0]+M01*p[qp][1])*fe_values[rt].value(i,qp) + 
                              (M10*p[qp][0]+M11*p[qp][1])*fe_values[it].value(i,qp));
        }
      }
      cell->get_dof_indices (local_dof_indices);
      for ( unsigned i = 0; i < dofs_per_cell; i++ )
      {
        system_rhs(local_dof_indices[i]) += cell_rhs(i);
        for ( unsigned j = 0; j < dofs_per_cell; j++ )
          system_matrix.add (local_dof_indices[i], local_dof_indices[j], cell_matrix(i, j));    
      } 
    }
  }    
