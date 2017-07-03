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
  void MySolver<dim,no_time_steps>::rt_propagtion_backward()
  {
    m_computing_timer.enter_section(__func__);

    m_workspace_ng = m_Psi_d;
    m_workspace_ng -= m_all_Psi[no_time_steps-1];

    constraints.distribute(m_workspace_ng);
    m_Psi = m_workspace_ng;

    MPI::MyComplexTools::AssembleSystem_mulvz( dof_handler, fe, constraints, m_Psi, std::complex<double>(0,-1), system_matrix, system_rhs );
    solve_eq2();
    m_all_p[no_time_steps-1] = m_Psi;
    output_vec( "p_" + to_string(no_time_steps-1) + ".vtu", m_all_p[no_time_steps-1] );

    for( int i=no_time_steps-2; i>0; i-- )
    {
      assemble_system(i);
      solve();
      m_all_p[i] = m_sol;

      //output_vec( "p_" + to_string(i) + ".vtu", m_all_p[i] );
      
      double N = MPI::MyComplexTools::Particle_Number( mpi_communicator, dof_handler, fe, m_Psi );
      if(m_root)printf( "b: %g %g\n", double(i)*m_dt, N );
    }
    m_computing_timer.exit_section();
  }

  template <int dim, int no_time_steps>
  void MySolver<dim,no_time_steps>::assemble_system ( const int idx )
  {
    m_computing_timer.enter_section(__func__);

    constraints.distribute(m_all_Psi[idx]);
    constraints.distribute(m_all_Psi[idx+1]);
    m_workspace = m_all_Psi[idx];
    m_workspace_2 = m_all_Psi[idx+1];

    const QGauss<dim> quadrature_formula(fe.degree+1);

    CPotential<dim,no_time_steps> potential2 = m_potential;
    m_potential.set_time( double(idx)*m_dt );
    potential2.set_time( double(idx+1)*m_dt );

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
    vector<Vector<double>> Psi0(n_q_points,Vector<double>(2));
    vector<Vector<double>> Psi1(n_q_points,Vector<double>(2));
    vector<Vector<double>> p(n_q_points,Vector<double>(2));
    vector<vector<Tensor<1,dim>>> p_grad(n_q_points, vector<Tensor<1,dim>>(2));
 
    const double dt = -m_dt; //reversed time direction
    const double dth = 0.5*dt;
    const double gamdth = m_gs*dth;
    
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

        for( unsigned qp=0; qp<n_q_points; qp++ )
        {
          double JxW = fe_values.JxW(qp);
          double alpha = 0.5 * m_gs * ((Psi0[qp][0]*Psi0[qp][0] - Psi0[qp][1]*Psi0[qp][1]) + (Psi1[qp][0]*Psi1[qp][0] - Psi1[qp][1]*Psi1[qp][1])); 
          double beta = gamdth * (Psi0[qp][0]*Psi0[qp][1] + Psi1[qp][0]*Psi1[qp][1]);

          double pot = 0.5*(potential2.value(fe_values.quadrature_point(qp)) + m_potential.value(fe_values.quadrature_point(qp))) 
                     + m_gs*(Psi0[qp][0]*Psi0[qp][0] + Psi0[qp][1]*Psi0[qp][1] + Psi1[qp][0]*Psi1[qp][0] + Psi1[qp][1]*Psi1[qp][1]) + alpha;
  
          for (unsigned i=0; i<dofs_per_cell; i++ )
          {
            for (unsigned j=0; j<dofs_per_cell; j++ )
            {
              cell_matrix(i,j) += JxW*((1.0-beta) * fe_values[rt].value(i,qp)*fe_values[rt].value(j,qp) - dth * (fe_values[rt].gradient(i,qp)*fe_values[it].gradient(j,qp) + pot*fe_values[rt].value(i,qp)*fe_values[it].value(j,qp)) 
                                      +(1.0-beta) * fe_values[it].value(i,qp)*fe_values[it].value(j,qp) + dth * (fe_values[it].gradient(i,qp)*fe_values[rt].gradient(j,qp) + pot*fe_values[it].value(i,qp)*fe_values[rt].value(j,qp)));
            }
            cell_rhs(i) += JxW*((1.0+beta)*p[qp][0]*fe_values[rt].value(i,qp) + dth * (p_grad[qp][1]*fe_values[it].gradient(i,qp) + pot*p[qp][1]*fe_values[it].value(i,qp)) 
                               +(1.0+beta)*p[qp][1]*fe_values[it].value(i,qp) - dth * (p_grad[qp][0]*fe_values[rt].gradient(i,qp) + pot*p[qp][0]*fe_values[rt].value(i,qp)));
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

  template <int dim, int no_time_steps>
  void MySolver<dim,no_time_steps>::solve_eq2 ()
  {
    m_computing_timer.enter_section(__func__);

    
    m_sol=0;

    SolverControl solver_control (m_sol.size(), 1e-15);
    PETScWrappers::SolverGMRES solver (solver_control, mpi_communicator);
    PETScWrappers::PreconditionSOR preconditioner(system_matrix);
    solver.solve (system_matrix, m_sol, system_rhs, preconditioner);
    
    constraints.distribute (m_sol);
    m_Psi = m_sol;

/*
    SolverControl solver_control;
    PETScWrappers::SparseDirectMUMPS solver(solver_control, mpi_communicator);
    solver.set_symmetric_mode(false);
    solver.solve(system_matrix, m_sol, system_rhs);
    constraints.distribute (m_sol);
    m_Psi = m_sol;
*/
    m_computing_timer.exit_section();
  }  