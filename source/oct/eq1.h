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
  void MySolver<no_time_steps,no_lam>::rt_propagtion_forward()
  {
    for( int i=1; i<no_time_steps; i++ )
    {
      m_timeindex = i-1;
      DoIter();
      m_all_Psi[i] = m_Psi;

//       double N = Particle_Number(m_Psi);
//       printf( "f: %g %g\n", double(i)*m_dt, N );

//       CPotential<no_time_steps,no_lam> Potential ( &m_all_lambdas, m_current_time_index );
//       VectorTools::interpolate(dof_handler, Potential, pot_vec );;
//       output_vec( "pot-" + to_string(m_current_time_index) + ".vtu", pot_vec );
//       output_vec( "fw-" + to_string(i) + ".gnuplot", m_Psi );
    }
    output_vec( "fw-60.gnuplot", m_Psi );
    double N = Particle_Number(m_Psi);
    printf( "N = %g\n", N );
  }

  template <int no_time_steps, int no_lam>
  void MySolver<no_time_steps,no_lam>::assemble_system ()
  {
    const QGauss<1> quadrature_formula(fe.degree+1);

    CPotential<no_time_steps,no_lam> Potential0 ( m_all_lambdas, m_timeindex );
    CPotential<no_time_steps,no_lam> Potential1 ( m_all_lambdas, m_timeindex+1 );
    
    const FEValuesExtractors::Scalar rt (0);
    const FEValuesExtractors::Scalar it (1);

    system_matrix=0;

    FEValues<1> fe_values (fe, quadrature_formula, update_values|update_gradients|update_quadrature_points|update_JxW_values);

    const unsigned dofs_per_cell = fe.dofs_per_cell;
    const unsigned n_q_points    = quadrature_formula.size();

    FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);

    vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    vector<Vector<double>> Psi(n_q_points,Vector<double>(2));
    vector<vector<Tensor<1,1>>> Psi_grad(n_q_points, vector<Tensor<1,1>>(2));
    vector<Vector<double>> Psi_t(n_q_points,Vector<double>(2));
    vector<vector<Tensor<1,1>>> Psi_t_grad(n_q_points, vector<Tensor<1,1>>(2));
 
    double JxW, pot=0, tmp1a, tmp1b, tmp2, sum_re, sum_im, sum_req, sum_imq;

    const double fak2 = 0.5*m_dt;
    const double fak4 = 0.25*m_gs*m_dt;
    const double fak8 = 0.125*m_gs*m_dt;

    DoFHandler<1>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
      cell_matrix=0;

      cell->get_dof_indices (local_dof_indices);

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
	        
	for (unsigned i=0; i<dofs_per_cell; i++ )
        {
          for (unsigned j=0; j<dofs_per_cell; j++ )
          {
            cell_matrix(i,j) += JxW*((1.0-tmp2)*fe_values[rt].value(i,qp)*fe_values[rt].value(j,qp) + (1.0+tmp2)*fe_values[it].value(i,qp)*fe_values[it].value(j,qp)
                                     - tmp1a*fe_values[rt].value(i,qp)*fe_values[it].value(j,qp) - fak2*(fe_values[rt].gradient(i,qp)*fe_values[it].gradient(j,qp) + pot*fe_values[rt].value(i,qp)*fe_values[it].value(j,qp))
                                     + tmp1b*fe_values[it].value(i,qp)*fe_values[rt].value(j,qp) + fak2*(fe_values[it].gradient(i,qp)*fe_values[rt].gradient(j,qp) + pot*fe_values[it].value(i,qp)*fe_values[rt].value(j,qp)));
	  }
	}
	
        for (unsigned int i=0; i<dofs_per_cell; ++i)
        {
          for (unsigned int j=0; j<dofs_per_cell; ++j) 
            system_matrix.add (local_dof_indices[i],local_dof_indices[j], cell_matrix(i,j));
        }
      }
    }
    map<types::global_dof_index,double> boundary_values;
    VectorTools::interpolate_boundary_values (dof_handler, 0, ZeroFunction<1>(2), boundary_values);
    MatrixTools::apply_boundary_values (boundary_values, system_matrix, sol, system_rhs);
  }
  
  template <int no_time_steps, int no_lam>
  void MySolver<no_time_steps,no_lam>::assemble_rhs ()
  {
    const QGauss<1> quadrature_formula(fe.degree+1);

    CPotential<no_time_steps,no_lam> Potential0 ( m_all_lambdas, m_timeindex );
    CPotential<no_time_steps,no_lam> Potential1 ( m_all_lambdas, m_timeindex+1 );
    
    const FEValuesExtractors::Scalar rt (0);
    const FEValuesExtractors::Scalar it (1);

    system_rhs=0;

    FEValues<1> fe_values (fe, quadrature_formula, update_values|update_gradients|update_quadrature_points|update_JxW_values);

    const unsigned dofs_per_cell = fe.dofs_per_cell;
    const unsigned n_q_points    = quadrature_formula.size();

    vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    vector<Vector<double>> Psi(n_q_points,Vector<double>(2));
    vector<vector<Tensor<1,1>>> Psi_grad(n_q_points, vector<Tensor<1,1>>(2));
    vector<Vector<double>> Psi_t(n_q_points,Vector<double>(2));
    vector<vector<Tensor<1,1>>> Psi_t_grad(n_q_points, vector<Tensor<1,1>>(2));
 
    double JxW, pot=0, tmp1, tmp2, sum_re, sum_im;

    const double fak2 = 0.5*m_dt;
    const double fak8 = 0.125*m_gs*m_dt;

    DoFHandler<1>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {

      fe_values.reinit (cell);
      fe_values.get_function_values(m_Psi, Psi);
      fe_values.get_function_gradients(m_Psi, Psi_grad);
      fe_values.get_function_values(m_Psi_t, Psi_t);
      fe_values.get_function_gradients(m_Psi_t, Psi_t_grad);

      cell->get_dof_indices (local_dof_indices);
      for( unsigned qp=0; qp<n_q_points; qp++ )
      {
        JxW = fe_values.JxW(qp);
        pot = 0.5*(Potential0.value(fe_values.quadrature_point(qp))+Potential1.value(fe_values.quadrature_point(qp))); 

	sum_re = Psi[qp][0]+Psi_t[qp][0];
	sum_im = Psi[qp][1]+Psi_t[qp][1];
        tmp1 = fak8*(sum_re*sum_re + sum_im*sum_im);
	
        for (unsigned i=0; i<dofs_per_cell; i++ )
        {
          system_rhs(local_dof_indices[i]) += JxW*(-fak2*((Psi_grad[qp][1]+Psi_t_grad[qp][1])*fe_values[rt].gradient(i,qp) + pot*sum_im*fe_values[rt].value(i,qp))
					      +fak2*((Psi_grad[qp][0]+Psi_t_grad[qp][0])*fe_values[it].gradient(i,qp) + pot*sum_re*fe_values[it].value(i,qp))
					      +(Psi_t[qp][0]-Psi[qp][0])*fe_values[rt].value(i,qp) - tmp1*sum_im*fe_values[rt].value(i,qp)  
					      +(Psi_t[qp][1]-Psi[qp][1])*fe_values[it].value(i,qp) + tmp1*sum_re*fe_values[it].value(i,qp));
        }
      }
    }
    
    map<types::global_dof_index,double> boundary_values;
    VectorTools::interpolate_boundary_values (dof_handler, 0, ZeroFunction<1>(2), boundary_values);
    MatrixTools::apply_boundary_values (boundary_values, system_matrix, sol, system_rhs);

    m_workspace = system_rhs;
    VectorTools::integrate_difference ( dof_handler,  m_workspace, ZeroFunction<1>(2), m_error_per_cell,  QGauss<1>(fe.degree+1), VectorTools::L2_norm);
    m_res = m_error_per_cell.l2_norm();    
  }
  
  template <int no_time_steps, int no_lam>
  void MySolver<no_time_steps,no_lam>::solve_eq1 ()
  {
    SparseDirectUMFPACK A_direct;
    A_direct.initialize(system_matrix);
    A_direct.vmult (sol, system_rhs);
    
//     SolverControl solver_control (sol.size(), 1e-15);
//     
//     sol=0;
//     PETScWrappers::SolverGMRES solver (solver_control, mpi_communicator);
//     PETScWrappers::PreconditionSOR preconditioner(system_matrix);
//     solver.solve (system_matrix, sol, system_rhs, preconditioner);
//     
//     constraints.distribute (sol);
//     m_Psi = sol;
  }
  
  template <int no_time_steps, int no_lam>
  void MySolver<no_time_steps,no_lam>::DoIter()
  {
    m_Psi_t = 0;
    assemble_rhs();
    //cout << "m_res = " << m_res << endl;       
    do
    {	
      assemble_system();
      solve_eq1();
      m_Psi_t.add( -1, sol );
      assemble_rhs();
      //cout << "m_res = " << m_res << endl;       
    }
    while( m_res >  1e-14 ); 

    m_Psi = m_Psi_t;
  }
    
