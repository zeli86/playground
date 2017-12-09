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

//    std::complex<double> z = MyComplexTools::MPI::L2_dot_product(mpi_communicator, dof_handler, fe, m_Psi_d, m_Psi );
//    z *= std::complex<double>(0,1);
//    if(m_root) printf( "z == %g + i %g \n", real(z), imag(z) );
//    save( "Psi_T.bin" );
    
/*
    m_workspace_ng = m_Psi_d;
    m_workspace_ng -= m_all_Psi[no_time_steps-1];
    constraints.distribute(m_workspace_ng);
    m_Psi = m_workspace_ng;
*/

    MyComplexTools::MPI::AssembleSystem_mulvz( dof_handler, fe, constraints, m_Psi, std::complex<double>(0,1), system_matrix, system_rhs );
    m_sol=0;
    solve_eq1();
    m_all_p[no_time_steps-1] = m_sol;
    m_N_pT = MyComplexTools::MPI::Particle_Number( mpi_communicator, dof_handler, fe, m_Psi );
    if(m_root) printf( "m_N_pT == %g\n", m_N_pT );

    output_vec( "p_" + to_string(ex) + ".vtu", m_all_p[no_time_steps-1] );
    for( int i=no_time_steps-2; i>0; i-- )
    {
      MyComplexTools::MPI::AssembleSystem_LIN_Step( dof_handler, fe, constraints, m_Psi, -m_dt, system_matrix, system_rhs );
      solve_eq1();
      assemble_system(i);
      solve_cg();
            
      m_all_p[i] = m_sol;      

      //output_vec( "p_" + to_string(i) + ".vtu", m_all_p[i] );
      //double N = MyComplexTools::MPI::Particle_Number( mpi_communicator, dof_handler, fe, m_Psi );
      //if(m_root) printf( "b: %g %g\n", double(i)*m_dt, N );
    }
    m_computing_timer.exit_section();
  }

  template <int dim, int no_time_steps>
  void MySolver<dim,no_time_steps>::assemble_system ( const int idx )
  {
    m_computing_timer.enter_section(__func__);

    m_potential.set_time( double(idx)*m_dt );

    constraints.distribute(m_all_Psi[idx]);
    m_workspace = m_all_Psi[idx];

    const QGauss<dim> quadrature_formula(fe.degree+1);

    const FEValuesExtractors::Scalar rt (0);
    const FEValuesExtractors::Scalar it (1);
        
    system_matrix = 0;
    system_rhs = 0;

    FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_quadrature_points|update_JxW_values);

    const unsigned dofs_per_cell = fe.dofs_per_cell;
    const unsigned n_q_points = quadrature_formula.size();

    FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);
    Vector<double> cell_rhs (dofs_per_cell);

    vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    vector<Vector<double>> Psi(n_q_points,Vector<double>(2));
    vector<Vector<double>> p(n_q_points,Vector<double>(2));
 
    const double dt = -m_dt; //reversed time direction
    
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
      if( cell->is_locally_owned() )
      {
        cell_matrix = 0;
        cell_rhs = 0;

        fe_values.reinit (cell);
        fe_values.get_function_values(m_Psi, p);
        fe_values.get_function_values(m_workspace, Psi );

        for( unsigned qp=0; qp<n_q_points; qp++ )
        {
          double JxW = fe_values.JxW(qp);
          double a = p[qp][0];
          double b = p[qp][1];
          double c = Psi[qp][0];
          double d = Psi[qp][1];
          double pot = m_potential.value(fe_values.quadrature_point(qp)) + 2*m_gs*(c*c+d*d);

          double beta = sqrt(c*c+d*d+pot*pot);
          double f1 = cos(dt*beta);
          double f2 = sin(dt*beta)/beta;
          double M00 = f1 + f2*d;
          double M01 = f2*(pot-c);
          double M10 = f2*(-pot-c);
          double M11 = f1 - f2*d;
          
          for (unsigned i=0; i<dofs_per_cell; i++ )
          {
            for (unsigned j=0; j<dofs_per_cell; j++ )
            {
              cell_matrix(i,j) += JxW*(fe_values[rt].value(i,qp)*fe_values[rt].value(j,qp) + fe_values[it].value(i,qp)*fe_values[it].value(j,qp));                        
            }
            cell_rhs(i) += JxW*((M00*a-M01*b)*fe_values[rt].value(i,qp) + (M01*a+M11*b)*fe_values[it].value(i,qp));
          }
        }
        cell->get_dof_indices (local_dof_indices);
        constraints.distribute_local_to_global(cell_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs);
      }
    }
    system_matrix.compress(VectorOperation::add);
    system_rhs.compress(VectorOperation::add);
    m_computing_timer.exit_section();
  }
