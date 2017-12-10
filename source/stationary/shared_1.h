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

/** Želimir Marojević
 */

#ifndef __shared_1_h__
#define __shared_1_h__
#include <stdio.h>

  template <int dim>
  void MySolver<dim>::compute_E_lin( LA::MPI::Vector& vec, double& T, double& N, double& W )
  {
    m_computing_timer.enter_section(__func__);
    
    constraints.distribute(vec);
    m_workspace_1 = vec;
    
    CPotential<dim> Potential( m_omega );
    const QGauss<dim>  quadrature_formula(fe.degree+1);
    FEValues<dim> fe_values (fe, quadrature_formula, update_gradients|update_values|update_JxW_values|update_quadrature_points);

    const unsigned int n_q_points = quadrature_formula.size();
    vector<double> vec_vals(n_q_points);
    vector<Tensor<1, dim> > vec_grad(n_q_points);

    double JxW, vec_val_q;
    
    double tmp1[3]={}, res[3]={};
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for ( ; cell!=endc; ++cell )
    {
      if( cell->is_locally_owned() )
      {
        fe_values.reinit (cell);
        fe_values.get_function_values(m_workspace_1, vec_vals);
        fe_values.get_function_gradients(m_workspace_1, vec_grad);
        for ( unsigned qp=0; qp<n_q_points; qp++ )
        {
          JxW = fe_values.JxW(qp);
          vec_val_q = vec_vals[qp]*vec_vals[qp];
          tmp1[0] += JxW*( vec_grad[qp]*vec_grad[qp] + Potential.value(fe_values.quadrature_point(qp))*vec_val_q );
          tmp1[1] += JxW*vec_val_q;
          tmp1[2] += JxW*vec_val_q*vec_val_q;
        }
      }
    }
    MPI_Allreduce( tmp1, res, 3, MPI_DOUBLE, MPI_SUM, mpi_communicator);
    T=res[0]; N=res[1]; W=res[2];
    m_computing_timer.exit_section();
  }

  template <int dim>
  void MySolver<dim>::estimate_error ( double& err )
  {
    m_computing_timer.enter_section(__func__);
    
    CPotential<dim> Potential( m_omega );
    const QGauss<dim> quadrature_formula(fe.degree+1);
   
    constraints.distribute(m_Psi_ref);
    m_workspace_1=m_Psi_ref;
    
    m_system_rhs=0;
    m_system_matrix=0;

    FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_gradients|update_JxW_values|update_quadrature_points);

    const unsigned dofs_per_cell = fe.dofs_per_cell;
    const unsigned n_q_points = quadrature_formula.size();

    Vector<double> cell_rhs (dofs_per_cell);
    FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);

    vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    vector<double> vals(n_q_points);
    vector<Tensor<1,dim>> grads(n_q_points);

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for ( ; cell!=endc; ++cell )
    {
      if( cell->is_locally_owned() )
      {
        cell_rhs=0;
        cell_matrix=0;

        fe_values.reinit (cell);
        fe_values.get_function_values(m_workspace_1, vals);
        fe_values.get_function_gradients(m_workspace_1, grads);

        for ( unsigned qp=0; qp<n_q_points; qp++ )
        {
          double JxW = fe_values.JxW(qp);
          double Q1 = Potential.value(fe_values.quadrature_point(qp)) - m_mu + m_gs*(vals[qp]*vals[qp]);

          for ( unsigned i=0; i<dofs_per_cell; i++ )
          {
            cell_rhs(i) += JxW*(grads[qp]*fe_values.shape_grad(i,qp) + Q1*vals[qp]*fe_values.shape_value(i,qp));
            for ( unsigned j=0; j<dofs_per_cell; j++ )
              cell_matrix(i,j) += JxW*(fe_values.shape_value(i,qp)*fe_values.shape_value(j,qp));
          }
        }
        cell->get_dof_indices (local_dof_indices);
        constraints.distribute_local_to_global(cell_matrix, cell_rhs, local_dof_indices, m_system_matrix, m_system_rhs);
      }
    }
    m_system_rhs.compress(VectorOperation::add);   
    m_system_matrix.compress(VectorOperation::add);   

    solve();
    
    m_workspace_1=m_newton_update;
    VectorTools::integrate_difference ( dof_handler, m_workspace_1, ZeroFunction<dim>(2), m_error_per_cell, QGauss<dim>(fe.degree+2), VectorTools::L2_norm);    
    const double total_local_error = m_error_per_cell.l2_norm();
    err = std::sqrt (Utilities::MPI::sum (total_local_error * total_local_error, MPI_COMM_WORLD));     
  
    m_computing_timer.exit_section();
  }
 
  template <int dim>
  bool MySolver<dim>::solve ()
  {
    m_computing_timer.enter_section(__func__);
    pcout << "Solving..." << endl;
/*        
    SolverControl solver_control;
    PETScWrappers::SparseDirectMUMPS solver(solver_control, mpi_communicator);
    solver.set_symmetric_mode(false);
    solver.solve(m_system_matrix, m_newton_update, m_system_rhs);
    constraints.distribute (m_newton_update);
*/    

    m_newton_update = 0;
    SolverControl solver_control (m_newton_update.size(), m_res*1e-4);
    //PETScWrappers::SolverGMRES solver (solver_control, mpi_communicator);
    PETScWrappers::SolverBicgstab solver (solver_control, mpi_communicator);
    
    //PETScWrappers::PreconditionBlockJacobi::AdditionalData adata;
    //PETScWrappers::PreconditionBlockJacobi preconditioner(m_system_matrix,adata);    
    
    PETScWrappers::PreconditionParaSails::AdditionalData adata;
    PETScWrappers::PreconditionParaSails preconditioner(m_system_matrix,adata);    

    try 
    {
      solver.solve(m_system_matrix, m_newton_update, m_system_rhs, preconditioner);
    }
    catch( ExceptionBase& e )
    {
      pcout << e.what() << endl;
      //pcout << "Possible singular matrix!" << endl;
      return false;
    }
    constraints.distribute (m_newton_update);

    m_computing_timer.exit_section();
    return true;
  }
#endif
