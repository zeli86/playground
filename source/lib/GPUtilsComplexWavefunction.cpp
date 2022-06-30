
#include "GPUtilsComplexWavefunction.hpp"
#include "mpi.h"

#include <deal.II/base/quadrature_lib.h>

namespace utils
{
  namespace complex_wavefunction
  {

    template <int iDim>
    void assemble_mulvz
    (
      IComplexWavefunction<iDim>* pSolver,
      const Vector<double>& vWavefunction,
      const std::complex<double> z,
      SparseMatrix<double>& oMatrix,
      Vector<double>& vRhs
    )
    {
      oMatrix = 0;
      vRhs = 0;

      const auto& fe = pSolver->get_fe();
      const auto& dof_handler = pSolver->get_dof_handler();

      const QGauss<iDim> quadrature_formula(fe.degree + 1);
      const FEValuesExtractors::Scalar rt(0);
      const FEValuesExtractors::Scalar it(1);

      FEValues<iDim> fe_values(fe, quadrature_formula, update_values | update_quadrature_points | update_JxW_values);

      const unsigned dofs_per_cell = fe.dofs_per_cell;
      const unsigned n_q_points = quadrature_formula.size();

      FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
      Vector<double> cell_rhs(dofs_per_cell);

      std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
      std::vector<Vector<double>> vals(n_q_points, Vector<double>(2));

      const double a = std::real(z);
      const double b = std::imag(z);

      typename DoFHandler<iDim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
      for (; cell != endc; ++cell)
      {
        cell_matrix = 0;
        cell_rhs = 0;

        fe_values.reinit(cell);
        fe_values.get_function_values(vWavefunction, vals);

        for (unsigned qp = 0; qp < n_q_points; qp++)
        {
          const double JxW = fe_values.JxW(qp);

          for (unsigned i = 0; i < dofs_per_cell; i++)
          {
            for (unsigned j = 0; j < dofs_per_cell; j++)
            {
              cell_matrix(i, j) += JxW * (fe_values[rt].value(i, qp) * fe_values[rt].value(j, qp) + fe_values[it].value(i, qp) * fe_values[it].value(j, qp));
            }
            const double c = vals[qp][0];
            const double d = vals[qp][1];
            cell_rhs(i) += JxW * ((a * c - b * d) * fe_values[rt].value(i, qp) + (b * c + a * d) * fe_values[it].value(i, qp));
          }
        }
        cell->get_dof_indices(local_dof_indices);
        for (unsigned i = 0; i < dofs_per_cell; i++)
        {
          vRhs(local_dof_indices[i]) += cell_rhs(i);
          for (unsigned j = 0; j < dofs_per_cell; j++)
          {
            oMatrix.add(local_dof_indices[i], local_dof_indices[j], cell_matrix(i, j));
          }
        }
      }
    }
    template void assemble_mulvz<1>(IComplexWavefunction<1>*, const Vector<double>&, const std::complex<double>, SparseMatrix<double>&, Vector<double>&);
    template void assemble_mulvz<2>(IComplexWavefunction<2>*, const Vector<double>&, const std::complex<double>, SparseMatrix<double>&, Vector<double>&);


    template <int iDim>
    void assemble_mulvz
    (
      IComplexWavefunction<iDim>* pSolver,
      const LA::MPI::Vector& vWavefunction,
      const std::complex<double> z,
      LA::MPI::SparseMatrix& oMatrix,
      LA::MPI::Vector& vRhs
    )
    {
      assert(vWavefunction.has_ghost_elements() == true);
      assert(vRhs.has_ghost_elements() == false);

      const auto& fe = pSolver->get_fe();
      const auto& dof_handler = pSolver->get_dof_handler();

      oMatrix = 0;
      vRhs = 0;

      const QGauss<iDim> quadrature_formula(fe.degree + 1);
      const FEValuesExtractors::Scalar rt(0);
      const FEValuesExtractors::Scalar it(1);

      FEValues<iDim> fe_values(fe, quadrature_formula, update_values | update_JxW_values);

      const unsigned dofs_per_cell = fe.dofs_per_cell;
      const unsigned n_q_points = quadrature_formula.size();

      FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
      Vector<double> cell_rhs(dofs_per_cell);
      std::vector<Vector<double>> vals(n_q_points, Vector<double>(2));

      std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

      const double a = std::real(z);
      const double b = std::imag(z);

      typename DoFHandler<iDim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
      for (; cell != endc; ++cell)
      {
        if (cell->is_locally_owned())
        {
          cell_rhs = 0;
          cell_matrix = 0;

          fe_values.reinit(cell);
          fe_values.get_function_values(vWavefunction, vals);

          for (unsigned qp = 0; qp < n_q_points; ++qp)
          {
            const double JxW = fe_values.JxW(qp);

            for (unsigned i = 0; i < dofs_per_cell; ++i)
            {
              for (unsigned j = 0; j < dofs_per_cell; ++j)
              {
                cell_matrix(i, j) += JxW * (fe_values[rt].value(i, qp) * fe_values[rt].value(j, qp) + fe_values[it].value(i, qp) * fe_values[it].value(j, qp));
              }
              double c = vals[qp][0];
              double d = vals[qp][1];
              cell_rhs(i) += JxW * ((a * c - b * d) * fe_values[rt].value(i, qp) + (b * c + a * d) * fe_values[it].value(i, qp));
            }
          }
          cell->get_dof_indices(local_dof_indices);
          pSolver->get_constraints().distribute_local_to_global(cell_matrix, cell_rhs, local_dof_indices, oMatrix, vRhs);
        }
      }
      vRhs.compress(VectorOperation::add);
      oMatrix.compress(VectorOperation::add);
    }
    template void assemble_mulvz<2>(IComplexWavefunction<2>*, const LA::MPI::Vector&, const std::complex<double> z, LA::MPI::SparseMatrix&, LA::MPI::Vector&);
    template void assemble_mulvz<3>(IComplexWavefunction<3>*, const LA::MPI::Vector&, const std::complex<double> z, LA::MPI::SparseMatrix&, LA::MPI::Vector&);


    template <int iDim>
    void assemble_lin_step
    (
      IComplexWavefunction<iDim>* pSolver,
      const Vector<double>& vWavefunction,
      const double dt,
      SparseMatrix<double>& oMatrix,
      Vector<double>& vRhs
    )
    {
      const auto& fe = pSolver->get_fe();
      const auto& dof_handler = pSolver->get_dof_handler();

      oMatrix = 0;
      vRhs = 0;

      const QGauss<iDim> quadrature_formula(fe.degree + 1);
      const FEValuesExtractors::Scalar rt(0);
      const FEValuesExtractors::Scalar it(1);
      const double dth = 0.5 * dt;

      FEValues<iDim> fe_values(fe, quadrature_formula,  update_values | update_gradients | update_quadrature_points | update_JxW_values);

      const unsigned dofs_per_cell = fe.dofs_per_cell;
      const unsigned n_q_points = quadrature_formula.size();

      FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
      Vector<double> cell_rhs(dofs_per_cell);

      std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
      std::vector<Vector<double>> vals(n_q_points, Vector<double>(2));
      std::vector<std::vector<Tensor<1, iDim>>> vals_grad(n_q_points, std::vector<Tensor<1, iDim>>(2));

      typename DoFHandler<iDim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
      for (; cell != endc; ++cell)
      {
        cell_matrix = 0;
        cell_rhs = 0;

        fe_values.reinit(cell);
        fe_values.get_function_values(vWavefunction, vals);
        fe_values.get_function_gradients(vWavefunction, vals_grad);

        for (unsigned qp = 0; qp < n_q_points; ++qp)
        {
          const double JxW = fe_values.JxW(qp);

          for (unsigned i = 0; i < dofs_per_cell; ++i)
          {
            for (unsigned j = 0; j < dofs_per_cell; ++j)
            {
              cell_matrix(i, j) += JxW * (fe_values[rt].value(i, qp) * fe_values[rt].value(j, qp) -
                                          dth * fe_values[rt].gradient(i, qp) * fe_values[it].gradient(j, qp) +
                                          fe_values[it].value(i, qp) * fe_values[it].value(j, qp) +
                                          dth * fe_values[it].gradient(i, qp) * fe_values[rt].gradient(j, qp));
            }
            cell_rhs(i) += JxW * (vals[qp][0] * fe_values[rt].value(i, qp) + dth * vals_grad[qp][1] * fe_values[rt].gradient(i, qp) +
                                  vals[qp][1] * fe_values[it].value(i, qp) - dth * vals_grad[qp][0] * fe_values[it].gradient(i, qp));
          }
        }
        cell->get_dof_indices(local_dof_indices);
        for (unsigned i = 0; i < dofs_per_cell; ++i)
        {
          vRhs(local_dof_indices[i]) += cell_rhs(i);
          for (unsigned j = 0; j < dofs_per_cell; ++j)
          {
            oMatrix.add(local_dof_indices[i], local_dof_indices[j], cell_matrix(i, j));
          }
        }
      }
    }
    template void assemble_lin_step<1>(IComplexWavefunction<1>*, const Vector<double>& vWavefunction, const double dt, SparseMatrix<double>&, Vector<double>&);
    template void assemble_lin_step<2>(IComplexWavefunction<2>*, const Vector<double>& vWavefunction, const double dt, SparseMatrix<double>&, Vector<double>&);

    template <int iDim>
    void assemble_lin_step
    (
      IComplexWavefunction<iDim>* pSolver,
      const LA::MPI::Vector& vWavefunction,
      const double dt,
      LA::MPI::SparseMatrix& oMatrix,
      LA::MPI::Vector& oRhs
    )
    {
      assert(vWavefunction.has_ghost_elements() == true);
      assert(oRhs.has_ghost_elements() == false);

      const auto& fe = pSolver->get_fe();
      const auto& dof_handler = pSolver->get_dof_handler();

      oMatrix = 0;
      oRhs = 0;

      const QGauss<iDim> quadrature_formula(fe.degree + 1);
      const FEValuesExtractors::Scalar rt(0);
      const FEValuesExtractors::Scalar it(1);
      const double dth = 0.5 * dt;

      FEValues<iDim> fe_values(fe, quadrature_formula, update_values | update_gradients | update_quadrature_points | update_JxW_values);

      const unsigned dofs_per_cell = fe.dofs_per_cell;
      const unsigned n_q_points = quadrature_formula.size();

      FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
      Vector<double> cell_rhs(dofs_per_cell);

      std::vector<Vector<double>> vals(n_q_points, Vector<double>(2));
      std::vector<std::vector<Tensor<1, iDim>>> vals_grad(n_q_points, std::vector<Tensor<1, iDim>>(2));
      std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

      typename DoFHandler<iDim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
      for (; cell != endc; ++cell)
      {
        if (cell->is_locally_owned())
        {
          cell_matrix = 0;
          cell_rhs = 0;

          fe_values.reinit(cell);
          fe_values.get_function_values(vWavefunction, vals);
          fe_values.get_function_gradients(vWavefunction, vals_grad);

          for (unsigned qp = 0; qp < n_q_points; ++qp)
          {
            const double JxW = fe_values.JxW(qp);

            for (unsigned i = 0; i < dofs_per_cell; ++i)
            {
              for (unsigned j = 0; j < dofs_per_cell; ++j)
              {
                cell_matrix(i, j) += JxW * (fe_values[rt].value(i, qp) * fe_values[rt].value(j, qp) -
                                            dth * fe_values[rt].gradient(i, qp) *
                                            fe_values[it].gradient(j, qp) + fe_values[it].value(i, qp) * fe_values[it].value(j, qp) +
                                            dth * fe_values[it].gradient(i, qp) * fe_values[rt].gradient(j, qp));
              }
              cell_rhs(i) += JxW * (vals[qp][0] * fe_values[rt].value(i, qp) + dth * vals_grad[qp][1] * fe_values[rt].gradient(i, qp) +
                                    vals[qp][1] * fe_values[it].value(i, qp) - dth * vals_grad[qp][0] * fe_values[it].gradient(i, qp));
            }
          }
          cell->get_dof_indices(local_dof_indices);
          pSolver->get_constraints().distribute_local_to_global(cell_matrix, cell_rhs, local_dof_indices, oMatrix, oRhs);
        }
      }
      oRhs.compress(VectorOperation::add);
      oMatrix.compress(VectorOperation::add);
    }
    template void assemble_lin_step<2>(IComplexWavefunction<2>*, const LA::MPI::Vector&, const double dt, LA::MPI::SparseMatrix&, LA::MPI::Vector&);
    template void assemble_lin_step<3>(IComplexWavefunction<3>*, const LA::MPI::Vector&, const double dt, LA::MPI::SparseMatrix&, LA::MPI::Vector&);

    template <int iDim>
    void assemble_lin_step
    (
      IComplexWavefunction<iDim>* pSolver,
      const Vector<double>& vWavefunction,
      const Function<iDim>& oPotential,
      const double dt,
      SparseMatrix<double>& oMatrix,
      Vector<double>& vRhs
    )
    {
      const auto& fe = pSolver->get_fe();
      const auto& dof_handler = pSolver->get_dof_handler();

      oMatrix = 0;
      vRhs = 0;

      const QGauss<iDim> quadrature_formula(fe.degree + 1);
      const FEValuesExtractors::Scalar rt(0);
      const FEValuesExtractors::Scalar it(1);
      const double dth = 0.5 * dt;

      FEValues<iDim> fe_values(fe, quadrature_formula, update_values | update_gradients | update_quadrature_points | update_JxW_values);

      const unsigned dofs_per_cell = fe.dofs_per_cell;
      const unsigned n_q_points = quadrature_formula.size();

      FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
      Vector<double> cell_rhs(dofs_per_cell);

      std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
      std::vector<Vector<double>> vals(n_q_points, Vector<double>(2));
      std::vector<std::vector<Tensor<1, iDim>>> vals_grad(n_q_points, std::vector<Tensor<1, iDim>>(2));

      typename DoFHandler<iDim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
      for (; cell != endc; ++cell)
      {
        cell_matrix = 0;
        cell_rhs = 0;

        fe_values.reinit(cell);
        fe_values.get_function_values(vWavefunction, vals);
        fe_values.get_function_gradients(vWavefunction, vals_grad);

        for (unsigned qp = 0; qp < n_q_points; ++qp)
        {
          const double JxW = fe_values.JxW(qp);
          const double pot = oPotential.value(fe_values.quadrature_point(qp));

          for (unsigned i = 0; i < dofs_per_cell; ++i)
          {
            for (unsigned j = 0; j < dofs_per_cell; ++j)
            {
              cell_matrix(i, j) += JxW * (fe_values[rt].value(i, qp) * fe_values[rt].value(j, qp) - dth * (fe_values[rt].gradient(i, qp) * fe_values[it].gradient(j, qp) + pot * fe_values[rt].value(i, qp) * fe_values[it].value(j, qp)) + fe_values[it].value(i, qp) * fe_values[it].value(j, qp) + dth * (fe_values[it].gradient(i, qp) * fe_values[rt].gradient(j, qp) + pot * fe_values[it].value(i, qp) * fe_values[rt].value(j, qp)));
            }
            cell_rhs(i) += JxW * (vals[qp][0] * fe_values[rt].value(i, qp) + dth * (vals_grad[qp][1] * fe_values[rt].gradient(i, qp) + pot * vals[qp][1] * fe_values[rt].value(i, qp)) +
                                  vals[qp][1] * fe_values[it].value(i, qp) - dth * (vals_grad[qp][0] * fe_values[it].gradient(i, qp) + pot * vals[qp][0] * fe_values[it].value(i, qp)));
          }
        }
        cell->get_dof_indices(local_dof_indices);
        for (unsigned i = 0; i < dofs_per_cell; ++i)
        {
          vRhs(local_dof_indices[i]) += cell_rhs(i);
          for (unsigned j = 0; j < dofs_per_cell; ++j)
          {
            oMatrix.add(local_dof_indices[i], local_dof_indices[j], cell_matrix(i, j));
          }
        }
      }
    }
    template void assemble_lin_step<1>(IComplexWavefunction<1>*, const Vector<double>&, const Function<1>&, const double, SparseMatrix<double>&, Vector<double>&);
    template void assemble_lin_step<2>(IComplexWavefunction<2>*, const Vector<double>&, const Function<2>&, const double, SparseMatrix<double>&, Vector<double>&);

    template <int iDim>
    void assemble_lin_step
    (
      IComplexWavefunction<iDim>* pSolver,
      const LA::MPI::Vector& vWavefunction,
      const Function<iDim>& oPotential,
      const double dt,
      LA::MPI::SparseMatrix& oMatrix,
      LA::MPI::Vector& vRhs
    )
    {
      assert(vWavefunction.has_ghost_elements() == true);
      assert(vRhs.has_ghost_elements() == false);

      const auto& fe = pSolver->get_fe();
      const auto& dof_handler = pSolver->get_dof_handler();

      oMatrix = 0;
      vRhs = 0;

      const QGauss<iDim> quadrature_formula(fe.degree + 1);
      const FEValuesExtractors::Scalar rt(0);
      const FEValuesExtractors::Scalar it(1);
      const double dth = 0.5 * dt;

      FEValues<iDim> fe_values(fe, quadrature_formula, update_values | update_gradients | update_quadrature_points | update_JxW_values);

      const unsigned dofs_per_cell = fe.dofs_per_cell;
      const unsigned n_q_points = quadrature_formula.size();

      FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
      Vector<double> cell_rhs(dofs_per_cell);
      std::vector<Vector<double>> vals(n_q_points, Vector<double>(2));
      std::vector<std::vector<Tensor<1, iDim>>> vals_grad(n_q_points, std::vector<Tensor<1, iDim>>(2));
      std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

      typename DoFHandler<iDim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
      for (; cell != endc; ++cell)
      {
        if (cell->is_locally_owned())
        {
          cell_matrix = 0;
          cell_rhs = 0;

          fe_values.reinit(cell);
          fe_values.get_function_values(vWavefunction, vals);
          fe_values.get_function_gradients(vWavefunction, vals_grad);

          for (unsigned qp = 0; qp < n_q_points; ++qp)
          {
            const double JxW = fe_values.JxW(qp);
            const double pot = oPotential.value(fe_values.quadrature_point(qp));

            for (unsigned i = 0; i < dofs_per_cell; ++i)
            {
              for (unsigned j = 0; j < dofs_per_cell; ++j)
              {
                cell_matrix(i, j) += JxW * (fe_values[rt].value(i, qp) * fe_values[rt].value(j, qp) -
                                            dth * (fe_values[rt].gradient(i, qp) * fe_values[it].gradient(j, qp) + pot * fe_values[rt].value(i, qp) * fe_values[it].value(j, qp)) +
                                            fe_values[it].value(i, qp) * fe_values[it].value(j, qp) +
                                            dth * (fe_values[it].gradient(i, qp) * fe_values[rt].gradient(j, qp) + pot * fe_values[it].value(i, qp) * fe_values[rt].value(j, qp)));
              }
              cell_rhs(i) += JxW * (vals[qp][0] * fe_values[rt].value(i, qp) + dth * (vals_grad[qp][1] * fe_values[rt].gradient(i, qp) + pot * vals[qp][1] * fe_values[rt].value(i, qp)) +
                                    vals[qp][1] * fe_values[it].value(i, qp) - dth * (vals_grad[qp][0] * fe_values[it].gradient(i, qp) + pot * vals[qp][0] * fe_values[it].value(i, qp)));
            }
          }
          cell->get_dof_indices(local_dof_indices);
          pSolver->get_constraints().distribute_local_to_global(cell_matrix, cell_rhs, local_dof_indices, oMatrix, vRhs);
        }
      }
      vRhs.compress(VectorOperation::add);
      oMatrix.compress(VectorOperation::add);
    }
    template void assemble_lin_step<2>(IComplexWavefunction<2>*, const LA::MPI::Vector&, const Function<2>&, const double, LA::MPI::SparseMatrix&, LA::MPI::Vector&);
    template void assemble_lin_step<3>(IComplexWavefunction<3>*, const LA::MPI::Vector&, const Function<3>&, const double, LA::MPI::SparseMatrix&, LA::MPI::Vector&);

    template <int iDim>
    void assemble_nl_step
    (
      IComplexWavefunction<iDim>* pSolver,
      const Vector<double>& vWavefunction,
      const Function<iDim>& oPotential,
      const double dt,
      const double gam,
      SparseMatrix<double>& oMatrix,
      Vector<double>& vRhs
    )
    {
      const auto& fe = pSolver->get_fe();
      const auto& dof_handler = pSolver->get_dof_handler();

      oMatrix = 0;
      vRhs = 0;

      const QGauss<iDim> quadrature_formula(fe.degree + 1);
      const FEValuesExtractors::Scalar rt(0);
      const FEValuesExtractors::Scalar it(1);
      const double dth = 0.5 * dt;

      FEValues<iDim> fe_values(fe, quadrature_formula, update_values | update_quadrature_points | update_JxW_values);

      const unsigned dofs_per_cell = fe.dofs_per_cell;
      const unsigned n_q_points = quadrature_formula.size();

      FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
      Vector<double> cell_rhs(dofs_per_cell);

      std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
      std::vector<Vector<double>> vals(n_q_points, Vector<double>(2));
      std::vector<std::vector<Tensor<1, iDim>>> vals_grad(n_q_points, std::vector<Tensor<1, iDim>>(2));

      typename DoFHandler<iDim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
      for (; cell != endc; ++cell)
      {
        cell_matrix = 0;
        cell_rhs = 0;

        fe_values.reinit(cell);
        fe_values.get_function_values(vWavefunction, vals);
        fe_values.get_function_gradients(vWavefunction, vals_grad);

        for (unsigned qp = 0; qp < n_q_points; ++qp)
        {
          const double JxW = fe_values.JxW(qp);
          const double pot = oPotential.value(fe_values.quadrature_point(qp));
          const double c = vals[qp][0];
          const double d = vals[qp][1];
          const double phi = -dt * (gam * (vals[qp][0] * vals[qp][0] + vals[qp][1] * vals[qp][1]) + oPotential.value(fe_values.quadrature_point(qp)));
          double a, b;
          sincos(phi, &b, &a);

          for (unsigned i = 0; i < dofs_per_cell; ++i)
          {
            for (unsigned j = 0; j < dofs_per_cell; ++j)
            {
              cell_matrix(i, j) += JxW * (fe_values[rt].value(i, qp) * fe_values[rt].value(j, qp) +
                                          fe_values[it].value(i, qp) * fe_values[it].value(j, qp));
            }
            cell_rhs(i) += JxW * ((a * c - b * d) * fe_values[rt].value(i, qp) + (b * c + a * d) * fe_values[it].value(i, qp));
          }
        }
        cell->get_dof_indices(local_dof_indices);
        for (unsigned i = 0; i < dofs_per_cell; ++i)
        {
          vRhs(local_dof_indices[i]) += cell_rhs(i);
          for (unsigned j = 0; j < dofs_per_cell; ++j)
          {
            oMatrix.add(local_dof_indices[i], local_dof_indices[j], cell_matrix(i, j));
          }
        }
      }
    }
    template void assemble_nl_step<1>(IComplexWavefunction<1>*, const Vector<double>&, const Function<1>&, const double, const double, SparseMatrix<double>&, Vector<double>&);
    template void assemble_nl_step<2>(IComplexWavefunction<2>*, const Vector<double>&, const Function<2>&, const double, const double, SparseMatrix<double>&, Vector<double>&);

    template <int iDim>
    void assemble_nl_step
    (
      IComplexWavefunction<iDim>* pSolver,
      const LA::MPI::Vector& vWavefunction,
      const Function<iDim>& oPotential,
      const double dt,
      const double gam,
      LA::MPI::SparseMatrix& oMatrix,
      LA::MPI::Vector& vRhs
    )
    {
      assert(vWavefunction.has_ghost_elements() == true);
      assert(vRhs.has_ghost_elements() == false);

      const auto& fe = pSolver->get_fe();
      const auto& dof_handler = pSolver->get_dof_handler();

      oMatrix = 0;
      vRhs = 0;

      const QGauss<iDim> quadrature_formula(fe.degree + 1);
      const FEValuesExtractors::Scalar rt(0);
      const FEValuesExtractors::Scalar it(1);

      FEValues<iDim> fe_values(fe, quadrature_formula, update_values | update_quadrature_points | update_JxW_values);

      const unsigned dofs_per_cell = fe.dofs_per_cell;
      const unsigned n_q_points = quadrature_formula.size();

      FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
      Vector<double> cell_rhs(dofs_per_cell);
      std::vector<Vector<double>> vals(n_q_points, Vector<double>(2));

      std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

      typename DoFHandler<iDim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
      for (; cell != endc; ++cell)
      {
        if (cell->is_locally_owned())
        {
          cell_rhs = 0;
          cell_matrix = 0;

          fe_values.reinit(cell);
          fe_values.get_function_values(vWavefunction, vals);

          for (unsigned qp = 0; qp < n_q_points; ++qp)
          {
            const double JxW = fe_values.JxW(qp);
            const double c = vals[qp][0];
            const double d = vals[qp][1];
            const double phi = -dt * (gam * (vals[qp][0] * vals[qp][0] + vals[qp][1] * vals[qp][1]) + oPotential.value(fe_values.quadrature_point(qp)));
            double a, b;
            sincos(phi, &b, &a);

            for (unsigned i = 0; i < dofs_per_cell; ++i)
            {
              for (unsigned j = 0; j < dofs_per_cell; ++j)
              {
                cell_matrix(i, j) += JxW * (fe_values[rt].value(i, qp) * fe_values[rt].value(j, qp) +
                                            fe_values[it].value(i, qp) * fe_values[it].value(j, qp));
              }
              cell_rhs(i) += JxW * ((a * c - b * d) * fe_values[rt].value(i, qp) + (b * c + a * d) * fe_values[it].value(i, qp));
            }
          }
          cell->get_dof_indices(local_dof_indices);
          pSolver->get_constraints().distribute_local_to_global(cell_matrix, cell_rhs, local_dof_indices, oMatrix, vRhs);
        }
      }
      vRhs.compress(VectorOperation::add);
      oMatrix.compress(VectorOperation::add);
    }
    template void assemble_nl_step<2>(IComplexWavefunction<2>*, const LA::MPI::Vector&, const Function<2>&, const double, const double, LA::MPI::SparseMatrix&, LA::MPI::Vector&);
    template void assemble_nl_step<3>(IComplexWavefunction<3>*, const LA::MPI::Vector&, const Function<3>&, const double, const double, LA::MPI::SparseMatrix&, LA::MPI::Vector&);

    template <int iDim>
    void assemble_nl_step
    (
      IComplexWavefunction<iDim>* pSolver,
      const LA::MPI::Vector& vWavefunction,
      const double gamdt,
      LA::MPI::SparseMatrix& oMatrix,
      LA::MPI::Vector& vRhs
    )
    {
      assert(vWavefunction.has_ghost_elements() == true);
      assert(vRhs.has_ghost_elements() == false);

      const auto& fe = pSolver->get_fe();
      const auto& dof_handler = pSolver->get_dof_handler();

      oMatrix = 0;
      vRhs = 0;

      const QGauss<iDim> quadrature_formula(fe.degree + 1);
      const FEValuesExtractors::Scalar rt(0);
      const FEValuesExtractors::Scalar it(1);

      FEValues<iDim> fe_values(fe, quadrature_formula, update_values | update_quadrature_points | update_JxW_values);

      const unsigned dofs_per_cell = fe.dofs_per_cell;
      const unsigned n_q_points = quadrature_formula.size();

      FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
      Vector<double> cell_rhs(dofs_per_cell);
      std::vector<Vector<double>> vals(n_q_points, Vector<double>(2));
      std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

      typename DoFHandler<iDim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
      for (; cell != endc; ++cell)
      {
        if (cell->is_locally_owned())
        {
          cell_rhs = 0;
          cell_matrix = 0;

          fe_values.reinit(cell);
          fe_values.get_function_values(vWavefunction, vals);

          for (unsigned qp = 0; qp < n_q_points; ++qp)
          {
            const double JxW = fe_values.JxW(qp);
            const double c = vals[qp][0];
            const double d = vals[qp][1];
            const double phi = -gamdt * (vals[qp][0] * vals[qp][0] + vals[qp][1] * vals[qp][1]);
            double a, b;
            sincos(phi, &b, &a);

            for (unsigned i = 0; i < dofs_per_cell; ++i)
            {
              for (unsigned j = 0; j < dofs_per_cell; ++j)
              {
                cell_matrix(i, j) += JxW * (fe_values[rt].value(i, qp) * fe_values[rt].value(j, qp) +
                                            fe_values[it].value(i, qp) * fe_values[it].value(j, qp));
              }
              cell_rhs(i) += JxW * ((a * c - b * d) * fe_values[rt].value(i, qp) + (b * c + a * d) * fe_values[it].value(i, qp));
            }
          }
          cell->get_dof_indices(local_dof_indices);
          pSolver->get_constraints().distribute_local_to_global(cell_matrix, cell_rhs, local_dof_indices, oMatrix, vRhs);
        }
      }
      vRhs.compress(VectorOperation::add);
      oMatrix.compress(VectorOperation::add);
    }
    template void assemble_nl_step<2>(IComplexWavefunction<2>*, const LA::MPI::Vector&, const double, LA::MPI::SparseMatrix&, LA::MPI::Vector&);
    template void assemble_nl_step<3>(IComplexWavefunction<3>*, const LA::MPI::Vector&, const double, LA::MPI::SparseMatrix&, LA::MPI::Vector&);
  }
}