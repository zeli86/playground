
#pragma once

#include <boost/log/trivial.hpp>
#include <stdio.h>

template <int dim>
void MySolver<dim>::estimate_error(double& err)
{
  TimerOutput::Scope timing_section(m_computing_timer, "");

  CPotential<dim> Potential(m_omega);
  const QGauss<dim> quadrature_formula(this->m_FE.degree + 1);

  this->m_constraints.distribute(this->m_Psi_Ref);
  this->m_Workspace[0] = this->m_Psi_Ref;

  this->m_System_RHS = 0;
  this->m_System_Matrix = 0;

  FEValues<dim> fe_values(this->m_FE, quadrature_formula, update_values | update_gradients | update_JxW_values | update_quadrature_points);

  const unsigned dofs_per_cell = this->m_FE.dofs_per_cell;
  const unsigned n_q_points = quadrature_formula.size();

  Vector<double> cell_rhs(dofs_per_cell);
  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);

  vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  vector<double> vals(n_q_points);
  vector<Tensor<1, dim>> grads(n_q_points);

  typename DoFHandler<dim>::active_cell_iterator cell = this->m_DOF_Handler.begin_active(), endc = this->m_DOF_Handler.end();
  for (; cell != endc; ++cell)
  {
    if (cell->is_locally_owned())
    {
      cell_rhs = 0;
      cell_matrix = 0;

      fe_values.reinit(cell);
      fe_values.get_function_values(this->m_Workspace[0], vals);
      fe_values.get_function_gradients(this->m_Workspace[0], grads);

      for (unsigned qp = 0; qp < n_q_points; qp++)
      {
        const double JxW = fe_values.JxW(qp);
        const double Q1 = Potential.value(fe_values.quadrature_point(qp)) - m_rMu + m_gs * (vals[qp] * vals[qp]);

        for (unsigned i = 0; i < dofs_per_cell; i++)
        {
          cell_rhs(i) += JxW * (grads[qp] * fe_values.shape_grad(i, qp) + Q1 * vals[qp] * fe_values.shape_value(i, qp));
          for (unsigned j = 0; j < dofs_per_cell; j++)
          {
            cell_matrix(i, j) += JxW * (fe_values.shape_value(i, qp) * fe_values.shape_value(j, qp));
          }
        }
      }
      cell->get_dof_indices(local_dof_indices);
      this->m_constraints.distribute_local_to_global(cell_matrix, cell_rhs, local_dof_indices, this->m_System_Matrix, this->m_System_RHS);
    }
  }
  this->m_System_RHS.compress(VectorOperation::add);
  this->m_System_Matrix.compress(VectorOperation::add);

  solve();

  this->m_Workspace[0] = this->m_Search_Direction;
  VectorTools::integrate_difference(this->m_DOF_Handler,  this->m_Workspace[0], ZeroFunction<dim>(2), this->m_error_per_cell, QGauss<dim>(this->m_FE.degree + 2), VectorTools::L2_norm);
  const double total_local_error = this->m_error_per_cell.l2_norm();
  err = std::sqrt(Utilities::MPI::sum(total_local_error * total_local_error, MPI_COMM_WORLD));
}

