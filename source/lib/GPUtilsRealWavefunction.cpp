
#include "GPUtilsRealWavefunction.hpp"
#include "deal.II/base/mpi.h"
#include "deal.II/base/quadrature_lib.h"
#include "deal.II/lac/generic_linear_algebra.h"
#include <tuple>
#include <vector>

namespace LA
{
  using namespace dealii::LinearAlgebraPETSc;
}

namespace utils
{
  namespace real_wavefunction
  {
    /**
     * @brief
     *
     * @tparam iDim
     * @param pBase
     * @param vWavefunction
     * @param oPotential
     * @return std::tuple<double, double, double>
     */
    template <int iDim>
    std::tuple<double, double, double>
    GP
    (
      IRealWavefunction<iDim>* pBase,
      const Vector<double>& vWavefunction,
      const Function<iDim>& oPotential
    )
    {
      const auto& fe = pBase->get_fe();
      const auto& dof_handler = pBase->get_dof_handler();

      const QGauss<iDim>  quadrature_formula(fe.degree + 1);
      FEValues<iDim> fe_values(fe, quadrature_formula, update_gradients | update_values | update_JxW_values | update_quadrature_points);

      const unsigned int n_q_points = quadrature_formula.size();
      std::vector<double> vec_vals(n_q_points);
      std::vector<Tensor<1, iDim>> vec_grad(n_q_points);

      Point<3> tmp1{0, 0, 0};
      typename DoFHandler<iDim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
      for (; cell != endc; ++cell)
      {
        fe_values.reinit(cell);
        fe_values.get_function_values(vWavefunction, vec_vals);
        fe_values.get_function_gradients(vWavefunction, vec_grad);
        for (unsigned qp = 0; qp < n_q_points; qp++)
        {
          const double JxW = fe_values.JxW(qp);
          const double vec_val_q = vec_vals[qp] * vec_vals[qp];
          tmp1[0] += JxW * (vec_grad[qp] * vec_grad[qp] + oPotential.value(fe_values.quadrature_point(qp)) * vec_val_q);
          tmp1[1] += JxW * vec_val_q;
          tmp1[2] += JxW * vec_val_q * vec_val_q;
        }
      }
      return std::make_tuple(tmp1[0], tmp1[1], tmp1[2]);
    }

    template std::tuple<double, double, double> GP<1>(IRealWavefunction<1>*, const Vector<double>&, const Function<1>&);
    template std::tuple<double, double, double> GP<2>(IRealWavefunction<2>*, const Vector<double>&, const Function<2>&);

    /**
     * @brief
     *
     * @tparam iDim
     * @param pBase
     * @param vWavefunction
     * @param oPotential
     * @param oMpiCommunicator
     * @return std::tuple<double, double, double>
     */
    template <int iDim>
    std::tuple<double, double, double>
    GP
    (
      IRealWavefunction<iDim>* pBase,
      const LA::MPI::Vector& vWavefunction,
      const Function<iDim>& oPotential,
      MPI_Comm oMpiCommunicator
    )
    {
      assert(vWavefunction.has_ghost_elements());

      const auto& fe = pBase->get_fe();
      const auto& dof_handler = pBase->get_dof_handler();

      const QGauss<iDim>  quadrature_formula(fe.degree + 1);
      FEValues<iDim> fe_values(fe, quadrature_formula, update_gradients | update_values | update_JxW_values | update_quadrature_points);

      const unsigned int n_q_points = quadrature_formula.size();
      std::vector<double> vec_vals(n_q_points);
      std::vector<Tensor<1, iDim>> vec_grad(n_q_points);

      Point<3> tmp1{0, 0, 0};
      typename DoFHandler<iDim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
      for (; cell != endc; ++cell)
      {
        if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);
          fe_values.get_function_values(vWavefunction, vec_vals);
          fe_values.get_function_gradients(vWavefunction, vec_grad);
          for (unsigned qp = 0; qp < n_q_points; qp++)
          {
            const double JxW = fe_values.JxW(qp);
            const double vec_val_q = vec_vals[qp] * vec_vals[qp];
            tmp1[0] += JxW * (vec_grad[qp] * vec_grad[qp] + oPotential.value(fe_values.quadrature_point(qp)) * vec_val_q);
            tmp1[1] += JxW * vec_val_q;
            tmp1[2] += JxW * vec_val_q * vec_val_q;
          }
        }
      }

      MPI_Allreduce(tmp1.begin_raw(), tmp1.begin_raw(), iDim, MPI_DOUBLE, MPI_SUM, oMpiCommunicator);
      return std::make_tuple(tmp1[0], tmp1[1], tmp1[2]);
    }

    template std::tuple<double, double, double> GP<2>(IRealWavefunction<2>*, const LA::MPI::Vector&, const Function<2>&, MPI_Comm);
    template std::tuple<double, double, double> GP<3>(IRealWavefunction<3>*, const LA::MPI::Vector&, const Function<3>&, MPI_Comm);

    /**
     * @brief
     *
     * @tparam iDim
     * @param pBase
     * @param vWavefunction
     * @param oPotential
     * @return double
     */
    template <int iDim>
    double
    MU
    (
      IRealWavefunction<iDim>* pBase,
      const Vector<double>& vWavefunction,
      const Function<iDim>& oPotential
    )
    {
      const auto& fe = pBase->get_fe();
      const auto& dof_handler = pBase->get_dof_handler();

      const QGauss<iDim>  quadrature_formula(fe.degree + 1);
      FEValues<iDim> fe_values(fe, quadrature_formula, update_gradients | update_values | update_JxW_values | update_quadrature_points);

      const unsigned int n_q_points = quadrature_formula.size();
      std::vector<double> vec_vals(n_q_points);
      std::vector<Tensor<1, iDim>> vec_grad(n_q_points);

      Point<2> tmp{0, 0};
      typename DoFHandler<iDim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
      for (; cell != endc; ++cell)
      {
        fe_values.reinit(cell);
        fe_values.get_function_values(vWavefunction, vec_vals);
        fe_values.get_function_gradients(vWavefunction, vec_grad);
        for (unsigned qp = 0; qp < n_q_points; qp++)
        {
          const double JxW = fe_values.JxW(qp);
          const double valq = vec_vals[qp] * vec_vals[qp];
          tmp[0] += JxW * (vec_grad[qp] * vec_grad[qp] + oPotential.value(fe_values.quadrature_point(qp)) * valq);
          tmp[1] += JxW * valq;
        }
      }
      return (tmp[0] / tmp[1]);
    }

    template double MU<1>(IRealWavefunction<1>*, const Vector<double>&, const Function<1>&);
    template double MU<2>(IRealWavefunction<2>*, const Vector<double>&, const Function<2>&);


    /**
     * @brief
     *
     * @tparam iDim
     * @param pBase
     * @param vWavefunction
     * @param oPotential
     * @param oMpiCommunicator
     * @return double
     */
    template <int iDim>
    double
    MU
    (
      IRealWavefunction<iDim>* pBase,
      const LA::MPI::Vector& vWavefunction,
      const Function<iDim>& oPotential,
      MPI_Comm oMpiCommunicator
    )
    {
      assert(vWavefunction.has_ghost_elements());

      const auto& fe = pBase->get_fe();
      const auto& dof_handler = pBase->get_dof_handler();

      const QGauss<iDim>  quadrature_formula(fe.degree + 1);
      FEValues<iDim> fe_values(fe, quadrature_formula, update_gradients | update_values | update_JxW_values | update_quadrature_points);

      const unsigned int n_q_points = quadrature_formula.size();
      std::vector<double> vec_vals(n_q_points);
      std::vector<Tensor<1, iDim>> vec_grad(n_q_points);

      Point<2> tmp{0, 0};
      typename DoFHandler<iDim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
      for (; cell != endc; ++cell)
      {
        if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);
          fe_values.get_function_values(vWavefunction, vec_vals);
          fe_values.get_function_gradients(vWavefunction, vec_grad);
          for (unsigned qp = 0; qp < n_q_points; qp++)
          {
            const double JxW = fe_values.JxW(qp);
            const double valq = vec_vals[qp] * vec_vals[qp];
            tmp[0] += JxW * (vec_grad[qp] * vec_grad[qp] + oPotential.value(fe_values.quadrature_point(qp)) * valq);
            tmp[1] += JxW * valq;
          }
        }
      }

      const auto rNom =  Utilities::MPI::sum(tmp[0], oMpiCommunicator);
      const auto rDen =  Utilities::MPI::sum(tmp[1], oMpiCommunicator);
      return (rNom / rDen);
    }

    template double MU<2>(IRealWavefunction<2>*, const LA::MPI::Vector&, const Function<2>&, MPI_Comm);
    template double MU<3>(IRealWavefunction<3>*, const LA::MPI::Vector&, const Function<3>&, MPI_Comm);

    /**
     * @brief
     *
     * @tparam dim
     * @param pBase
     * @param vWavefunction
     * @param oPotential
     * @param rMu
     * @param rG
     * @return SparseMatrix<double>
     */
    template <int dim>
    SparseMatrix<double> assemble_jacobian
    (
      IRealWavefunction<dim>* pBase,
      const Vector<double>& vWavefunction,
      const Function<dim>& oPotential,
      const double rMu,
      const double rG
    )
    {
      const auto& fe = pBase->get_fe();
      const auto& dof_handler = pBase->get_dof_handler();

      const QGauss<dim> quadrature_formula(fe.degree + 1);

      SparseMatrix<double> oMatrix{};

      FEValues<dim> fe_values(fe, quadrature_formula, update_values | update_gradients | update_JxW_values | update_quadrature_points);

      const unsigned dofs_per_cell = fe.dofs_per_cell;
      const unsigned n_q_points = quadrature_formula.size();

      FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);

      std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
      std::vector<double> vals(n_q_points);

      typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
      for (; cell != endc; ++cell)
      {
        cell_matrix = 0;

        fe_values.reinit(cell);
        fe_values.get_function_values(vWavefunction, vals);

        for (unsigned qp = 0; qp < n_q_points; ++qp)
        {
          const double JxW = fe_values.JxW(qp);
          const double Q2 = oPotential.value(fe_values.quadrature_point(qp)) - rMu + 3.0 * rG * vals[qp] * vals[qp];

          for (unsigned i = 0; i < dofs_per_cell; ++i)
            for (unsigned j = 0; j < dofs_per_cell; ++j)
            {
              cell_matrix(i, j) += JxW * (fe_values.shape_grad(i, qp) * fe_values.shape_grad(j, qp) + Q2 * fe_values.shape_value(i, qp) * fe_values.shape_value(j, qp));
            }
        }
        cell->get_dof_indices(local_dof_indices);
        for (unsigned i = 0; i < dofs_per_cell; ++i)
          for (unsigned j = 0; j < dofs_per_cell; ++j)
          {
            oMatrix.add(local_dof_indices[i], local_dof_indices[j], cell_matrix(i, j));
          }
      }
      return oMatrix;
    }

    template SparseMatrix<double> assemble_jacobian<1>(IRealWavefunction<1>*, const Vector<double>&, const Function<1>&, const double, const double);
    template SparseMatrix<double> assemble_jacobian<2>(IRealWavefunction<2>*, const Vector<double>&, const Function<2>&, const double, const double);

    template <int dim>
    void assemble_jacobian
    (
      IRealWavefunction<dim>* pBase,
      const LA::MPI::Vector& vWavefunction,
      const Function<dim>& oPotential,
      LA::MPI::SparseMatrix& oMatrix,
      const double mu,
      const double gs
    )
    {
      assert(vWavefunction.has_ghost_elements() == true);

      const auto& fe = pBase->get_fe();
      const auto& dof_handler = pBase->get_dof_handler();

      const QGauss<dim> quadrature_formula(fe.degree + 1);

      oMatrix = 0;

      FEValues<dim> fe_values(fe, quadrature_formula, update_values | update_gradients | update_JxW_values | update_quadrature_points);

      const unsigned dofs_per_cell = fe.dofs_per_cell;
      const unsigned n_q_points = quadrature_formula.size();

      FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);

      std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
      std::vector<double> vals(n_q_points);

      typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
      for (; cell != endc; ++cell)
      {
        if (cell->is_locally_owned())
        {
          cell_matrix = 0;

          fe_values.reinit(cell);
          fe_values.get_function_values(vWavefunction, vals);

          for (unsigned qp = 0; qp < n_q_points; ++qp)
          {
            double JxW = fe_values.JxW(qp);
            double Q2 = oPotential.value(fe_values.quadrature_point(qp)) - mu + 3.0 * gs * vals[qp] * vals[qp];

            for (unsigned i = 0; i < dofs_per_cell; ++i)
              for (unsigned j = 0; j < dofs_per_cell; ++j)
              {
                cell_matrix(i, j) += JxW * (fe_values.shape_grad(i, qp) * fe_values.shape_grad(j, qp) + Q2 * fe_values.shape_value(i, qp) * fe_values.shape_value(j, qp));
              }
          }
          cell->get_dof_indices(local_dof_indices);
          pBase->get_constraints().distribute_local_to_global(cell_matrix, local_dof_indices, oMatrix);
        }
      }
      oMatrix.compress(VectorOperation::add);
    }

    template void assemble_jacobian<2>(IRealWavefunction<2>*, const LA::MPI::Vector&, const Function<2>&, LA::MPI::SparseMatrix&, const double, const double);
    template void assemble_jacobian<3>(IRealWavefunction<3>*, const LA::MPI::Vector&, const Function<3>&, LA::MPI::SparseMatrix&, const double, const double);

    /**
     * @brief
     *
     * @tparam dim
     * @param pBase
     * @param vWavefunction
     * @param oPotential
     * @param rMu
     * @param rG
     * @param vGradient
     * @return double
     */
    template <int dim>
    double assemble_L2gradient
    (
      IRealWavefunction<dim>* pBase,
      const Vector<double>& vWavefunction,
      const Function<dim>& oPotential,
      const double rMu,
      const double rG,
      Vector<double>& vGradient
    )
    {
      const auto& fe = pBase->get_fe();
      const auto& dof_handler = pBase->get_dof_handler();

      const QGauss<dim> quadrature_formula(fe.degree + 1);

      std::fill(vGradient.begin(), vGradient.end(), 0);

      FEValues<dim> fe_values(fe, quadrature_formula, update_values | update_gradients | update_JxW_values | update_quadrature_points);

      const unsigned dofs_per_cell = fe.dofs_per_cell;
      const unsigned n_q_points = quadrature_formula.size();

      Vector<double> cell_rhs(dofs_per_cell);

      std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
      std::vector<Tensor<1, dim>> grad_vals(n_q_points);
      std::vector<double> vals(n_q_points);

      typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
      for (; cell != endc; ++cell)
      {
        cell_rhs = 0;
        fe_values.reinit(cell);
        fe_values.get_function_values(vWavefunction, vals);
        fe_values.get_function_gradients(vWavefunction, grad_vals);

        for (unsigned qp = 0; qp < n_q_points; ++qp)
        {
          const double JxW = fe_values.JxW(qp);
          const double Q1 = oPotential.value(fe_values.quadrature_point(qp)) - rMu + rG * vals[qp] * vals[qp];

          for (unsigned i = 0; i < dofs_per_cell; ++i)
          {
            cell_rhs(i) += JxW * (grad_vals[qp] * fe_values.shape_grad(i, qp) + Q1 * vals[qp] * fe_values.shape_value(i, qp));
          }
        }
        cell->get_dof_indices(local_dof_indices);
        for (unsigned i = 0; i < dofs_per_cell; ++i)
        {
          vGradient(local_dof_indices[i]) += cell_rhs(i);
        }
      }
      return vGradient.l2_norm();
    }

    template double assemble_L2gradient<1>(IRealWavefunction<1>*, const Vector<double>&, const Function<1>&, const double, const double, Vector<double>&);
    template double assemble_L2gradient<2>(IRealWavefunction<2>*, const Vector<double>&, const Function<2>&, const double, const double, Vector<double>&);

    /**
     * @brief
     *
     * @tparam dim
     * @param pBase
     * @param vWavefunction
     * @param oPotential
     * @param rMu
     * @param rG
     * @param vGradient
     * @return double
     */
    template <int dim>
    double assemble_L2gradient
    (
      IRealWavefunction<dim>* pBase,
      const LA::MPI::Vector& vWavefunction,
      const Function<dim>& oPotential,
      const double rMu,
      const double rG,
      LA::MPI::Vector& vGradient
    )
    {
      assert(vWavefunction.has_ghost_elements() == true);

      const auto& fe = pBase->get_fe();
      const auto& dof_handler = pBase->get_dof_handler();

      const QGauss<dim> quadrature_formula(fe.degree + 1);

      vGradient = 0;

      FEValues<dim> fe_values(fe, quadrature_formula, update_values | update_gradients | update_JxW_values | update_quadrature_points);

      const unsigned dofs_per_cell = fe.dofs_per_cell;
      const unsigned n_q_points = quadrature_formula.size();

      Vector<double> cell_rhs(dofs_per_cell);

      std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
      std::vector<Tensor<1, dim>> grad_vals(n_q_points);
      std::vector<double> vals(n_q_points);

      typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
      for (; cell != endc; ++cell)
      {
        if (cell->is_locally_owned())
        {
          cell_rhs = 0;

          fe_values.reinit(cell);
          fe_values.get_function_values(vWavefunction, vals);
          fe_values.get_function_gradients(vWavefunction, grad_vals);

          for (unsigned qp = 0; qp < n_q_points; ++qp)
          {
            const double JxW = fe_values.JxW(qp);
            const double Q1 = oPotential.value(fe_values.quadrature_point(qp)) - rMu + rG * vals[qp] * vals[qp];

            for (unsigned i = 0; i < dofs_per_cell; ++i)
            {
              cell_rhs(i) += JxW * (grad_vals[qp] * fe_values.shape_grad(i, qp) + Q1 * vals[qp] * fe_values.shape_value(i, qp));
            }
          }
          cell->get_dof_indices(local_dof_indices);
          pBase->get_constraints().distribute_local_to_global(cell_rhs, local_dof_indices, vGradient);
        }
      }
      vGradient.compress(VectorOperation::add);
      return vGradient.l2_norm();
    }

    template double assemble_L2gradient<2>(IRealWavefunction<2>*, const LA::MPI::Vector&, const Function<2>&, const double, const double, LA::MPI::Vector&);
    template double assemble_L2gradient<3>(IRealWavefunction<3>*, const LA::MPI::Vector&, const Function<3>&, const double, const double, LA::MPI::Vector&);
  }
}