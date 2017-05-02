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
  void MySolver<dim,no_time_steps,no_lam>::rt_propagtion_forward()
  {
    m_computing_timer.enter_section(__func__);
    
    for( int i=1; i<no_time_steps; i++ )
    {
      m_timeindex = i-1;
      DoIter();
      m_all_Psi[i] = m_Psi;
      
      if(m_root) printf( "f: %g \n", double(i)*m_dt );
    }
   
    m_computing_timer.exit_section();
  }
  
  template <int dim, int no_time_steps, int no_lam>
  void MySolver<dim,no_time_steps,no_lam>::assemble_system ()
  {
    m_computing_timer.enter_section(__func__);
    const QGauss<dim> quadrature_formula(fe.degree+1);

    CPotential<dim,no_time_steps,no_lam> Potential0 ( m_all_lambdas, m_timeindex );
    CPotential<dim,no_time_steps,no_lam> Potential1 ( m_all_lambdas, m_timeindex+1 );

    const FEValuesExtractors::Scalar rt (0);
    const FEValuesExtractors::Scalar it (1);

    system_matrix=0;

    FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_gradients|update_quadrature_points|update_JxW_values);

    const unsigned dofs_per_cell = fe.dofs_per_cell;
    const unsigned n_q_points = quadrature_formula.size();

    FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);

    vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    vector<Vector<double>> Psi(n_q_points,Vector<double>(2));
    vector<vector<Tensor<1,dim>>> Psi_grad(n_q_points, vector<Tensor<1,dim>>(2));
    vector<Vector<double>> Psi_t(n_q_points,Vector<double>(2));
    vector<vector<Tensor<1,dim>>> Psi_t_grad(n_q_points, vector<Tensor<1,dim>>(2));
 
    double JxW, pot=0, tmp1a, tmp1b, tmp2, sum_re, sum_req, sum_im, sum_imq;

    const double fak2 = 0.5*m_dt;
    const double fak4 = 0.25*m_gs*m_dt;
    const double fak8 = 0.125*m_gs*m_dt;

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
      if( cell->is_locally_owned() )
      {
        cell_matrix = 0;

        fe_values.reinit (cell);
        fe_values.get_function_values(m_Psi, Psi);
        fe_values.get_function_gradients(m_Psi, Psi_grad);
        fe_values.get_function_values(m_Psi_t, Psi_t);
        fe_values.get_function_gradients(m_Psi_t, Psi_t_grad);

        for( unsigned qp=0; qp<n_q_points; qp++ )
        {
          JxW = fe_values.JxW(qp);
          pot = 0.5*(Potential0.value(fe_values.quadrature_point(qp))+Potential1.value(fe_values.quadrature_point(qp))); 
  
          sum_re = Psi[qp][0]+Psi_t[qp][0];
          sum_im = Psi[qp][1]+Psi_t[qp][1];
          sum_req = sum_re*sum_re;
          sum_imq = sum_im*sum_im;
          tmp1a = fak8*(sum_req + 3*sum_imq);
          tmp1b = fak8*(sum_imq + 3*sum_req);
          tmp2 = fak4*sum_re*sum_im;

          for (unsigned int i=0; i<dofs_per_cell; i++ )
          {
            for (unsigned int j=0; j<dofs_per_cell; j++ )
            {
              cell_matrix(i,j) += JxW*(1.0-tmp2)*fe_values[rt].value(i,qp)*fe_values[rt].value(j,qp);
              cell_matrix(i,j) += JxW*(1.0+tmp2)*fe_values[it].value(i,qp)*fe_values[it].value(j,qp);
              cell_matrix(i,j) -= JxW*tmp1a*fe_values[rt].value(i,qp)*fe_values[it].value(j,qp);
              cell_matrix(i,j) -= JxW*fak2*(fe_values[rt].gradient(i,qp)*fe_values[it].gradient(j,qp) + pot*fe_values[rt].value(i,qp)*fe_values[it].value(j,qp));
              cell_matrix(i,j) += JxW*tmp1b*fe_values[it].value(i,qp)*fe_values[rt].value(j,qp);
              cell_matrix(i,j) += JxW*fak2*(fe_values[it].gradient(i,qp)*fe_values[rt].gradient(j,qp) + pot*fe_values[it].value(i,qp)*fe_values[rt].value(j,qp));
            }
          }
        }
        cell->get_dof_indices (local_dof_indices);
        constraints.distribute_local_to_global(cell_matrix, local_dof_indices, system_matrix );
      }
    }
    system_matrix.compress(VectorOperation::add);
    m_computing_timer.exit_section();
  }
  
  template <int dim, int no_time_steps, int no_lam>
  void MySolver<dim,no_time_steps,no_lam>::assemble_rhs ()
  {
    m_computing_timer.enter_section(__func__);
    const QGauss<dim> quadrature_formula(fe.degree+1);

    CPotential<dim,no_time_steps,no_lam> Potential0 ( m_all_lambdas, m_timeindex );
    CPotential<dim,no_time_steps,no_lam> Potential1 ( m_all_lambdas, m_timeindex+1 );

    const FEValuesExtractors::Scalar rt (0);
    const FEValuesExtractors::Scalar it (1);

    system_rhs=0;

    FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_gradients|update_quadrature_points|update_JxW_values);

    const unsigned dofs_per_cell = fe.dofs_per_cell;
    const unsigned n_q_points = quadrature_formula.size();

    Vector<double> cell_rhs (dofs_per_cell);

    vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    vector<Vector<double>> Psi(n_q_points,Vector<double>(2));
    vector<vector<Tensor<1,dim>>> Psi_grad(n_q_points, vector<Tensor<1,dim>>(2));
    vector<Vector<double>> Psi_t(n_q_points,Vector<double>(2));
    vector<vector<Tensor<1,dim>>> Psi_t_grad(n_q_points, vector<Tensor<1,dim>>(2));
 
    double JxW, pot=0, tmp1, sum_re, sum_im;

    const double fak2 = 0.5*m_dt;
    const double fak8 = 0.125*m_gs*m_dt;

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
      if( cell->is_locally_owned() )
      {
        cell_rhs = 0;

        fe_values.reinit (cell);
        fe_values.get_function_values(m_Psi, Psi);
        fe_values.get_function_gradients(m_Psi, Psi_grad);
        fe_values.get_function_values(m_Psi_t, Psi_t);
        fe_values.get_function_gradients(m_Psi_t, Psi_t_grad);

        for( unsigned int qp=0; qp<n_q_points; qp++ )
        {
          JxW = fe_values.JxW(qp);
          pot = 0.5*(Potential0.value(fe_values.quadrature_point(qp))+Potential1.value(fe_values.quadrature_point(qp))); 
          sum_re = Psi[qp][0]+Psi_t[qp][0];
          sum_im = Psi[qp][1]+Psi_t[qp][1];  
          tmp1 = fak8*(sum_re*sum_re + sum_im*sum_im);
  
          for ( unsigned i=0; i<dofs_per_cell; i++ )
          {
            cell_rhs(i) -= JxW*fak2*((Psi_grad[qp][1]+Psi_t_grad[qp][1])*fe_values[rt].gradient(i,qp) + pot*sum_im*fe_values[rt].value(i,qp));
            cell_rhs(i) += JxW*fak2*((Psi_grad[qp][0]+Psi_t_grad[qp][0])*fe_values[it].gradient(i,qp) + pot*sum_re*fe_values[it].value(i,qp));
            cell_rhs(i) += JxW*((Psi_t[qp][0]-Psi[qp][0])*fe_values[rt].value(i,qp) - tmp1*sum_im*fe_values[rt].value(i,qp) + 
                                (Psi_t[qp][1]-Psi[qp][1])*fe_values[it].value(i,qp) + tmp1*sum_re*fe_values[it].value(i,qp));
          }
        }
        cell->get_dof_indices (local_dof_indices);
        constraints.distribute_local_to_global( cell_rhs, local_dof_indices, system_rhs );
      }
    }
    system_rhs.compress(VectorOperation::add);
    m_res = system_rhs.l2_norm();

    m_computing_timer.exit_section();
  }  
  
  template <int dim, int no_time_steps, int no_lam>
  void MySolver<dim,no_time_steps,no_lam>::DoIter()
  {
    m_computing_timer.enter_section(__func__);

    m_Psi_t = 0;
    m_res = 0;
    assemble_rhs();
    do
    {
      assemble_system();
      solve_eq1();
      
      system_rhs=m_Psi_t;
      system_rhs.add( -1.0, m_sol );
      constraints.distribute(system_rhs);
      m_Psi_t=system_rhs;

      assemble_rhs();
    }
    while( m_res > 1e-15 ); 

    m_Psi = m_Psi_t;
    m_computing_timer.exit_section();
  }
    
  template <int dim, int no_time_steps, int no_lam>
  void MySolver<dim,no_time_steps,no_lam>::solve_eq1 ()
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

    m_computing_timer.exit_section();
  }    
