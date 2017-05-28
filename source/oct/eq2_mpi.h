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

  template <int dim, int no_time_steps, int no_lam>
  void MySolver<dim,no_time_steps,no_lam>::rt_propagtion_backward()
  {
    m_computing_timer.enter_section(__func__);
    compute_initial_p();
    
    for( int i=no_time_steps-2; i>0; i-- )
    {
      m_timeindex = i;
      m_workspace = m_all_Psi[i];
      m_workspace_2 = m_all_Psi[i+1];

      assemble_system_2();
      solve();
      m_all_p[i] = m_Psi;
      
      double N = Particle_Number(m_Psi);
      if(m_root)printf( "b: %g %g\n", double(i)*m_dt, N );
    }
    m_computing_timer.exit_section();
  }

  template <int dim, int no_time_steps, int no_lam>
  void MySolver<dim,no_time_steps,no_lam>::compute_initial_p()
  {
    m_computing_timer.enter_section(__func__);

    double re, im;
    Scalar_product( m_Psi_d, m_Psi, re, im );
    m_overlap = sqrt(re*re+im*im);
    printf( "overlap = %g\n", m_overlap );

    assemble_system_4_initial_p( im, -re );
    solve_eq2();
    m_all_p[no_time_steps-1] = m_Psi;
//     output_vec( "init_p.vtu", m_Psi );
    
    m_computing_timer.exit_section();
  }

  template <int dim, int no_time_steps, int no_lam>
  void MySolver<dim,no_time_steps,no_lam>::assemble_system_2 ()
  {
    m_computing_timer.enter_section(__func__);
    const QGauss<dim> quadrature_formula(fe.degree+1);

    CPotential<dim,no_time_steps,no_lam> Potential0 ( m_all_lambdas, m_omega, m_timeindex );
    CPotential<dim,no_time_steps,no_lam> Potential1 ( m_all_lambdas, m_omega, m_timeindex+1 );
    
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
    vector<Vector<double>> Psi0(n_q_points,Vector<double>(2));
    vector<Vector<double>> Psi1(n_q_points,Vector<double>(2));
    vector<Vector<double>> p(n_q_points,Vector<double>(2));
    vector<vector<Tensor<1,dim>>> p_grad(n_q_points, vector<Tensor<1,dim>>(2));
 
    double JxW, pot, tmp01, tmp02;
    
    const double dt = -m_dt; //reversed time direction
    const double dth = 0.5*dt;
    const double gamdt = m_gs*dt;
    
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
      if( cell->is_locally_owned() )
      {
        cell_matrix = 0;
        cell_rhs = 0;

        fe_values.reinit (cell);
        fe_values.get_function_values(m_Psi, p);
        fe_values.get_function_gradients(m_Psi, p_grad);

        fe_values.get_function_values(m_workspace, Psi1 );
        fe_values.get_function_values(m_workspace_2, Psi0 );

        for( unsigned int qp=0; qp<n_q_points; qp++ )
        {
          JxW = fe_values.JxW(qp);
          pot = 0.5*(Potential0.value(fe_values.quadrature_point(qp)) + Potential1.value(fe_values.quadrature_point(qp)));
          pot += m_gs*(Psi0[qp][0]*Psi0[qp][0] + Psi0[qp][1]*Psi0[qp][1] + Psi1[qp][0]*Psi1[qp][0] + Psi1[qp][1]*Psi1[qp][1]);
  
          tmp01 = 0.5*gamdt*(Psi1[qp][0]*Psi1[qp][1]+Psi0[qp][0]*Psi0[qp][1]);
          tmp02 = 0.5*m_gs*(Psi0[qp][0]*Psi0[qp][0]+Psi1[qp][0]*Psi1[qp][0]-Psi0[qp][1]*Psi0[qp][1]-Psi1[qp][1]*Psi1[qp][1]);
  
          for (unsigned i=0; i<dofs_per_cell; i++ )
          {
            for (unsigned j=0; j<dofs_per_cell; j++ )
            {
              cell_matrix(i,j) += JxW*((1.0-tmp01)*fe_values[rt].value(i,qp)*fe_values[rt].value(j,qp) - dth*(fe_values[rt].gradient(i,qp)*fe_values[it].gradient(j,qp) + (pot-tmp02)*fe_values[rt].value(i,qp)*fe_values[it].value(j,qp)) + 
                                       (1.0+tmp01)*fe_values[it].value(i,qp)*fe_values[it].value(j,qp) + dth*(fe_values[it].gradient(i,qp)*fe_values[rt].gradient(j,qp) + (pot+tmp02)*fe_values[it].value(i,qp)*fe_values[rt].value(j,qp)));
            }
            cell_rhs(i) += JxW*((1.0+tmp01)*p[qp][0]*fe_values[rt].value(i,qp) + dth*(p_grad[qp][1]*fe_values[rt].gradient(i,qp) + (pot-tmp02)*p[qp][1]*fe_values[rt].value(i,qp)) + 
                                (1.0-tmp01)*p[qp][1]*fe_values[it].value(i,qp) - dth*(p_grad[qp][0]*fe_values[it].gradient(i,qp) + (pot+tmp02)*p[qp][0]*fe_values[it].value(i,qp)));
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

  template <int dim, int no_time_steps, int no_lam>
  void MySolver<dim,no_time_steps,no_lam>::assemble_system_4_initial_p ( const double a, const double b )
  {
    m_computing_timer.enter_section(__func__);
    const QGauss<dim> quadrature_formula(fe.degree+1);

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
 
    m_Psi = m_Psi_d;
    
    double JxW;
    
    double a1 = a/(a*a+b*b);
    double b1 = b/(a*a+b*b);
    
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
      if( cell->is_locally_owned() )
      {
        cell_matrix = 0;
        cell_rhs = 0;

        fe_values.reinit (cell);
        fe_values.get_function_values(m_Psi, Psi);

        for( unsigned qp=0; qp<n_q_points; qp++ )
        {
          JxW = fe_values.JxW(qp);
          for ( unsigned i=0; i<dofs_per_cell; i++ )
          {
            for ( unsigned j=0; j<dofs_per_cell; j++ )
            {
               cell_matrix(i,j) += JxW*(a1*fe_values[rt].value(i,qp)*fe_values[rt].value(j,qp) + b1*fe_values[rt].value(i,qp)*fe_values[it].value(j,qp) + a1*fe_values[it].value(i,qp)*fe_values[it].value(j,qp) - b1*fe_values[it].value(i,qp)*fe_values[rt].value(j,qp)); 
            }
            cell_rhs(i) += JxW* (Psi[qp][0]*fe_values[rt].value(i,qp) + Psi[qp][1]*fe_values[it].value(i,qp));
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

  template <int dim, int no_time_steps, int no_lam>
  void MySolver<dim,no_time_steps,no_lam>::solve_eq2 ()
  {
    m_computing_timer.enter_section(__func__);

/*  
    SolverControl solver_control (m_sol.size(), 1e-15);
    
    m_sol=0;
    PETScWrappers::SolverGMRES solver (solver_control, mpi_communicator);
    PETScWrappers::PreconditionSOR preconditioner(system_matrix);
    solver.solve (system_matrix, m_sol, system_rhs, preconditioner);
    
    constraints.distribute (m_sol);
    m_Psi = m_sol;
*/
    SolverControl solver_control;
    PETScWrappers::SparseDirectMUMPS solver(solver_control, mpi_communicator);
    solver.set_symmetric_mode(false);
    solver.solve(system_matrix, m_sol, system_rhs);
    constraints.distribute (m_sol);
    m_Psi = m_sol;

    m_computing_timer.exit_section();
  }  