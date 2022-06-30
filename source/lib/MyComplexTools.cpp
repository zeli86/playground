
#include <cassert>
#include <cmath>
#include <complex>
#include <vector>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/generic_linear_algebra.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/vector_tools.h>
#include "MyComplexTools.h"

namespace MyComplexTools
{
  using namespace std;
  using namespace dealii;

  template <int dim>
  void Interpolate_R_to_C(const DoFHandler<dim>& dof_handler, const FE_Q<dim>& fe,
                          const Vector<double>& vec,
                          const DoFHandler<dim>& dof_handler_2,
                          const FESystem<dim>& fe_2,
                          const dealii::AffineConstraints<double>& constraints,
                          Vector<double>& ret)
  {
    // assert( vec.has_ghost_elements() == true );
    // assert( ret.has_ghost_elements() == true );

    const QGauss<dim> quadrature_formula(fe.degree + 1);
    const FEValuesExtractors::Scalar rt(0);
    const FEValuesExtractors::Scalar it(1);

    Vector<double> rhs(dof_handler_2.n_dofs());
    Vector<double> sol(dof_handler_2.n_dofs());

    SparseMatrix<double> matrix;
    SparsityPattern sparsity_pattern;
    DynamicSparsityPattern dsp(dof_handler_2.n_dofs(), dof_handler_2.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler_2, dsp);
    sparsity_pattern.copy_from(dsp);
    matrix.reinit(sparsity_pattern);

    rhs = 0;
    matrix = 0;

    FEValues<dim> fe_values(fe, quadrature_formula, update_values | update_JxW_values);
    FEValues<dim> fe_values_2(fe_2, quadrature_formula, update_values | update_JxW_values);

    const unsigned dofs_per_cell = fe_2.dofs_per_cell;
    const unsigned n_q_points = quadrature_formula.size();

    vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    vector<double> vals(n_q_points);
    Vector<double> cell_rhs(dofs_per_cell);
    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                   endc = dof_handler.end();
    typename DoFHandler<dim>::active_cell_iterator cell_2 = dof_handler_2.begin_active();
    for (; cell != endc; ++cell, ++cell_2)
    {
      cell_rhs = 0;
      cell_matrix = 0;

      fe_values.reinit(cell);
      fe_values_2.reinit(cell_2);
      fe_values.get_function_values(vec, vals);

      for (unsigned qp = 0; qp < n_q_points; ++qp)
      {
        const double JxW = fe_values_2.JxW(qp);
        const double tmp1 = vals[qp];

        for (unsigned i = 0; i < dofs_per_cell; ++i)
        {
          cell_rhs(i) += JxW * tmp1 * fe_values_2[rt].value(i, qp);
          for (unsigned j = 0; j < dofs_per_cell; ++j)
          {
            cell_matrix(i, j) += JxW * (fe_values_2[rt].value(i, qp) * fe_values_2[rt].value(j, qp) +
                                        fe_values_2[it].value(i, qp) * fe_values_2[it].value(j, qp));
          }
        }
      }
      cell_2->get_dof_indices(local_dof_indices);
      for (unsigned i = 0; i < dofs_per_cell; ++i)
      {
        rhs(local_dof_indices[i]) += cell_rhs(i);
        for (unsigned j = 0; j < dofs_per_cell; ++j)
        {
          matrix.add(local_dof_indices[i], local_dof_indices[j], cell_matrix(i, j));
        }
      }
    }

    map<types::global_dof_index, double> boundary_values;
    VectorTools::interpolate_boundary_values(dof_handler, 0, ZeroFunction<dim>(), boundary_values);
    MatrixTools::apply_boundary_values(boundary_values, matrix, sol, rhs);

    SparseDirectUMFPACK A_direct;
    A_direct.initialize(matrix);
    A_direct.vmult(sol, rhs);

    ret = sol;
  }
} // end of namespace

namespace MyComplexTools
{
  namespace MPI
  {
    namespace LA
    {
      using namespace dealii::LinearAlgebraPETSc;
    }
    using namespace dealii;


    template <int dim>
    std::complex<double> L2_dot_product(MPI_Comm mpi_communicator,
                                        const DoFHandler<dim>& dof_handler,
                                        const FESystem<dim>& fe, const LA::MPI::Vector& vec1,
                                        const LA::MPI::Vector& vec2)
    {
      assert(vec1.has_ghost_elements() == true);
      assert(vec2.has_ghost_elements() == true);

      double tmp[] = {0, 0};
      double retval[] = {0, 0};

      const QGauss<dim> quadrature_formula(fe.degree + 1);
      FEValues<dim> fe_values(fe, quadrature_formula, update_values | update_JxW_values);

      const unsigned n_q_points = quadrature_formula.size();
      vector<Vector<double>> vec_vals_1(n_q_points, Vector<double>(2));
      vector<Vector<double>> vec_vals_2(n_q_points, Vector<double>(2));

      typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                     endc = dof_handler.end();
      for (; cell != endc; ++cell)
      {
        if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);
          fe_values.get_function_values(vec1, vec_vals_1);
          fe_values.get_function_values(vec2, vec_vals_2);
          for (unsigned qp = 0; qp < n_q_points; ++qp)
          {
            const double JxW = fe_values.JxW(qp);
            tmp[0] += JxW * (vec_vals_1[qp][0] * vec_vals_2[qp][0] + vec_vals_1[qp][1] * vec_vals_2[qp][1]);
            tmp[1] += JxW * (vec_vals_1[qp][0] * vec_vals_2[qp][1] - vec_vals_1[qp][1] * vec_vals_2[qp][0]);
          }
        }
      }
      MPI_Allreduce(tmp, retval, 2, MPI_DOUBLE, MPI_SUM, mpi_communicator);
      return std::complex<double>(retval[0], retval[1]);
    }


    template <int dim>
    void Interpolate_R_to_C(MPI_Comm mpi_communicator,
                            const DoFHandler<dim>& dof_handler, const FE_Q<dim>& fe,
                            const LA::MPI::Vector& vec,
                            const DoFHandler<dim>& dof_handler_2,
                            const FESystem<dim>& fe_2,
                            const dealii::AffineConstraints<double>& constraints,
                            LA::MPI::Vector& ret)
    {
      assert(vec.has_ghost_elements() == true);
      assert(ret.has_ghost_elements() == true);

      const QGauss<dim> quadrature_formula(fe.degree + 1);
      const FEValuesExtractors::Scalar rt(0);
      const FEValuesExtractors::Scalar it(1);

      IndexSet locally_owned_dofs, locally_relevant_dofs;

      locally_owned_dofs = dof_handler_2.locally_owned_dofs();
      DoFTools::extract_locally_relevant_dofs(dof_handler_2, locally_relevant_dofs);

      LA::MPI::Vector rhs(locally_owned_dofs, mpi_communicator);
      LA::MPI::Vector sol(locally_owned_dofs, mpi_communicator);

      DynamicSparsityPattern dsp(locally_relevant_dofs);
      DoFTools::make_sparsity_pattern(dof_handler_2, dsp, constraints, false);
      SparsityTools::distribute_sparsity_pattern(dsp, dof_handler_2.n_locally_owned_dofs_per_processor(), mpi_communicator, locally_relevant_dofs);
      LA::MPI::SparseMatrix matrix;
      matrix.reinit(locally_owned_dofs, locally_owned_dofs, dsp, mpi_communicator);

      rhs = 0;
      sol = 0;
      matrix = 0;

      FEValues<dim> fe_values(fe, quadrature_formula, update_values | update_JxW_values);
      FEValues<dim> fe_values_2(fe_2, quadrature_formula, update_values | update_JxW_values);

      const unsigned dofs_per_cell = fe_2.dofs_per_cell;
      const unsigned n_q_points = quadrature_formula.size();

      vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

      vector<double> vals(n_q_points);
      Vector<double> cell_rhs(dofs_per_cell);
      FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);

      typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                     endc = dof_handler.end();
      typename DoFHandler<dim>::active_cell_iterator cell_2 = dof_handler_2.begin_active();
      for (; cell != endc; ++cell, ++cell_2)
      {
        if (cell->is_locally_owned())
        {
          cell_rhs = 0;
          cell_matrix = 0;

          fe_values.reinit(cell);
          fe_values_2.reinit(cell_2);
          fe_values.get_function_values(vec, vals);

          for (unsigned qp = 0; qp < n_q_points; qp++)
          {
            const double JxW = fe_values_2.JxW(qp);
            const double tmp1 = vals[qp];

            for (unsigned i = 0; i < dofs_per_cell; ++i)
            {
              cell_rhs(i) += JxW * tmp1 * fe_values_2[rt].value(i, qp);
              for (unsigned j = 0; j < dofs_per_cell; ++j)
              {
                cell_matrix(i, j) += JxW * (fe_values_2[rt].value(i, qp) * fe_values_2[rt].value(j, qp) +
                                            fe_values_2[it].value(i, qp) * fe_values_2[it].value(j, qp));
              }
            }
          }
          cell_2->get_dof_indices(local_dof_indices);
          constraints.distribute_local_to_global(cell_matrix, cell_rhs, local_dof_indices, matrix, rhs);
        }
      }
      rhs.compress(VectorOperation::add);
      matrix.compress(VectorOperation::add);

      /*
      SolverControl solver_control;
      PETScWrappers::SparseDirectMUMPS solver(solver_control, mpi_communicator);
      solver.set_symmetric_mode(false);
      solver.solve(matrix, sol, rhs);
      constraints.distribute (sol);
      */

      SolverControl solver_control(sol.size(), 1e-15);
      PETScWrappers::SolverBicgstab solver(solver_control, mpi_communicator);
      PETScWrappers::PreconditionBlockJacobi::AdditionalData adata;
      PETScWrappers::PreconditionBlockJacobi preconditioner(matrix, adata);
      solver.solve(matrix, sol, rhs, preconditioner);
      constraints.distribute(sol);

      ret = sol;
    }
  }
}
