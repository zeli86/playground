
#include "UtilsComplexWavefunction.hpp"
#include "mpi.h"

#include <deal.II/base/quadrature_lib.h>

namespace utils
{
  namespace complex_wavefunction
  {
    using namespace dealii;

    template <int iDim>
    double
    particle_number
    (
      IComplexWavefunction<iDim>* pBase,
      const Vector<double>& vec
    )
    {
      const auto& fe = pBase->get_fe();
      const auto& dof_handler = pBase->get_dof_handler();

      double retval{0};
      const QGauss<iDim> quadrature_formula(fe.degree + 1);

      FEValues<iDim> fe_values(fe, quadrature_formula, update_values | update_quadrature_points | update_JxW_values);

      const unsigned n_q_points = quadrature_formula.size();

      std::vector<Vector<double>> vals(n_q_points, Vector<double>(2));

      typename DoFHandler<iDim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
      for (; cell != endc; ++cell)
      {
        fe_values.reinit(cell);
        fe_values.get_function_values(vec, vals);

        for (unsigned qp = 0; qp < n_q_points; qp++)
        {
          retval += fe_values.JxW(qp) * (vals[qp][0] * vals[qp][0] + vals[qp][1] * vals[qp][1]);
        }
      }
      return retval;
    }
    template double particle_number(IComplexWavefunction<1>*, const Vector<double>&);
    template double particle_number(IComplexWavefunction<2>*, const Vector<double>&);

    template <int iDim>
    double
    particle_number
    (
      IComplexWavefunction<iDim>* pBase,
      const LA::MPI::Vector& vec,
      MPI_Comm mpi_communicator
    )
    {
      assert(vec.has_ghost_elements() == true);
      const auto& fe = pBase->get_fe();
      const auto& dof_handler = pBase->get_dof_handler();
      double tmp = 0;

      const QGauss<iDim> quadrature_formula(fe.degree + 1);
      FEValues<iDim> fe_values(fe, quadrature_formula, update_values | update_JxW_values);

      const unsigned n_q_points = quadrature_formula.size();
      std::vector<Vector<double>> vec_vals(n_q_points, Vector<double>(2));

      typename DoFHandler<iDim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
      for (; cell != endc; ++cell)
      {
        if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);
          fe_values.get_function_values(vec, vec_vals);
          for (unsigned qp = 0; qp < n_q_points; ++qp)
          {
            tmp += fe_values.JxW(qp) * (vec_vals[qp][0] * vec_vals[qp][0] + vec_vals[qp][1] * vec_vals[qp][1]);
          }
        }
      }
      return Utilities::MPI::sum(tmp, mpi_communicator);
    }
    template double particle_number<2>(IComplexWavefunction<2>*, const LA::MPI::Vector&, MPI_Comm);
    template double particle_number<3>(IComplexWavefunction<3>*, const LA::MPI::Vector&, MPI_Comm);

    template <int iDim>
    std::array<double, iDim>
    expectation_value_position
    (
      IComplexWavefunction<iDim>* pBase,
      const Vector<double>& vec
    )
    {
      const auto& fe = pBase->get_fe();
      const auto& dof_handler = pBase->get_dof_handler();

      const QGauss<iDim> quadrature_formula(fe.degree + 1);

      FEValues<iDim> fe_values(fe, quadrature_formula, update_values | update_quadrature_points | update_JxW_values);

      const unsigned dofs_per_cell = fe.dofs_per_cell;
      const unsigned n_q_points = quadrature_formula.size();

      std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
      std::vector<Vector<double>> vals(n_q_points, Vector<double>(2));

      Point<iDim> oTotalIntegral = {};

      typename DoFHandler<iDim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
      for (; cell != endc; ++cell)
      {
        fe_values.reinit(cell);
        fe_values.get_function_values(vec, vals);

        for (unsigned qp = 0; qp < n_q_points; ++qp)
        {
          oTotalIntegral += fe_values.JxW(qp) * (vals[qp][0] * vals[qp][0] + vals[qp][1] * vals[qp][1]) * fe_values.quadrature_point(qp);
        }
      }

      std::array<double, iDim> retval;
      for (int i = 0; i < iDim; ++i)
      {
        retval[i] = oTotalIntegral[i];
      }
      return retval;
    }
    template std::array<double, 1> expectation_value_position<1>(IComplexWavefunction<1>*, const Vector<double>&);
    template std::array<double, 2> expectation_value_position<2>(IComplexWavefunction<2>*, const Vector<double>&);

    template <int iDim>
    std::array<double, iDim>
    expectation_value_position
    (
      IComplexWavefunction<iDim>* pBase,
      const LA::MPI::Vector& vec,
      MPI_Comm mpi_communicator
    )
    {
      assert(vec.has_ghost_elements() == true);
      const auto& fe = pBase->get_fe();
      const auto& dof_handler = pBase->get_dof_handler();

      Point<iDim> oLocalIntegral = {};

      const QGauss<iDim> quadrature_formula(fe.degree + 1);
      FEValues<iDim> fe_values(fe, quadrature_formula, update_values | update_quadrature_points | update_JxW_values);

      const unsigned n_q_points = quadrature_formula.size();
      std::vector<Vector<double>> vec_vals(n_q_points, Vector<double>(2));

      typename DoFHandler<iDim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
      for (; cell != endc; ++cell)
      {
        if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);
          fe_values.get_function_values(vec, vec_vals);
          for (unsigned qp = 0; qp < n_q_points; ++qp)
          {
            oLocalIntegral += fe_values.JxW(qp) * (vec_vals[qp][0] * vec_vals[qp][0] + vec_vals[qp][1] * vec_vals[qp][1]) * fe_values.quadrature_point(qp);
          }
        }
      }
      std::array<double, iDim> oTotalIntegral = {};
      MPI_Allreduce(oLocalIntegral.begin_raw(), oTotalIntegral.data(), iDim, MPI_DOUBLE, MPI_SUM, mpi_communicator);
      return oTotalIntegral;
    }
    template std::array<double, 2> expectation_value_position<2>(IComplexWavefunction<2>*, const LA::MPI::Vector&, MPI_Comm);
    template std::array<double, 3> expectation_value_position<3>(IComplexWavefunction<3>*, const LA::MPI::Vector&, MPI_Comm);

    template <int iDim>
    std::array<double, iDim>
    expectation_value_momentum(IComplexWavefunction<iDim>* pBase, const Vector<double>& vec)
    {
      const auto& fe = pBase->get_fe();
      const auto& dof_handler = pBase->get_dof_handler();

      Point<iDim> oTotalIntegral = {};

      const QGauss<iDim> quadrature_formula(fe.degree + 1);
      FEValues<iDim> fe_values(fe, quadrature_formula, update_values | update_gradients | update_JxW_values);

      const unsigned n_q_points = quadrature_formula.size();
      std::vector<Vector<double>> vec_vals(n_q_points, Vector<double>(2));
      std::vector<std::vector<Tensor<1, iDim>>> vec_grads(n_q_points, std::vector<Tensor<1, iDim>>(2));

      typename DoFHandler<iDim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
      for (; cell != endc; ++cell)
      {
        fe_values.reinit(cell);
        fe_values.get_function_values(vec, vec_vals);
        for (unsigned qp = 0; qp < n_q_points; ++qp)
        {
          const double JxW = fe_values.JxW(qp);
          for (unsigned i = 0; i < iDim; ++i)
          {
            oTotalIntegral[i] += JxW * (vec_vals[qp][0] * vec_grads[qp][1][i] - vec_vals[qp][1] * vec_grads[qp][0][i]);
          }
        }
      }
      std::array<double, iDim> retval;
      for (int i = 0; i < iDim; ++i)
      {
        retval[i] = oTotalIntegral[i];
      }
      return retval;
    }
    template std::array<double, 1> expectation_value_momentum<1>(IComplexWavefunction<1>*, const Vector<double>&);
    template std::array<double, 2> expectation_value_momentum<2>(IComplexWavefunction<2>*, const Vector<double>&);

    template <int iDim>
    std::array<double, iDim>
    expectation_value_momentum
    (
      IComplexWavefunction<iDim>* pBase,
      const LA::MPI::Vector& vec,
      MPI_Comm mpi_communicator
    )
    {
      const auto& fe = pBase->get_fe();
      const auto& dof_handler = pBase->get_dof_handler();

      assert(vec.has_ghost_elements() == true);

      Point<iDim> oLocalIntegral = {};

      const QGauss<iDim> quadrature_formula(fe.degree + 1);
      FEValues<iDim> fe_values(fe, quadrature_formula, update_values | update_gradients | update_JxW_values);

      const unsigned n_q_points = quadrature_formula.size();
      std::vector<Vector<double>> vec_vals(n_q_points, Vector<double>(2));
      std::vector<std::vector<Tensor<1, iDim>>> vec_grads(n_q_points, std::vector<Tensor<1, iDim>>(2));

      typename DoFHandler<iDim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
      for (; cell != endc; ++cell)
      {
        if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);
          fe_values.get_function_values(vec, vec_vals);
          for (unsigned qp = 0; qp < n_q_points; qp++)
          {
            const double JxW = fe_values.JxW(qp);
            for (unsigned i = 0; i < iDim; ++i)
            {
              oLocalIntegral[i] += JxW * (vec_vals[qp][0] * vec_grads[qp][1][i] - vec_vals[qp][1] * vec_grads[qp][0][i]);
            }
          }
        }
      }
      std::array<double, iDim> oTotalIntegral = {};
      MPI_Allreduce(oLocalIntegral.begin_raw(), oTotalIntegral.data(), iDim, MPI_DOUBLE, MPI_SUM, mpi_communicator);
      return oTotalIntegral;
    }
    template std::array<double, 2> expectation_value_momentum<2>(IComplexWavefunction<2>*, const LA::MPI::Vector&, MPI_Comm);
    template std::array<double, 3> expectation_value_momentum<3>(IComplexWavefunction<3>*, const LA::MPI::Vector&, MPI_Comm);

    template <int iDim>
    std::array<double, iDim>
    expectation_value_width
    (
      IComplexWavefunction<iDim>* pBase,
      const Vector<double>& vec,
      const Point<iDim>& pos
    )
    {
      const auto& fe = pBase->get_fe();
      const auto& dof_handler = pBase->get_dof_handler();

      assert(vec.has_ghost_elements() == true);

      Point<iDim> oTotalIntegral = {};

      const QGauss<iDim> quadrature_formula(fe.degree + 1);
      FEValues<iDim> fe_values(fe, quadrature_formula, update_values | update_quadrature_points | update_JxW_values);

      const unsigned n_q_points = quadrature_formula.size();
      std::vector<Vector<double>> vec_vals(n_q_points, Vector<double>(2));

      typename DoFHandler<iDim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
      for (; cell != endc; ++cell)
      {
        fe_values.reinit(cell);
        fe_values.get_function_values(vec, vec_vals);
        for (unsigned qp = 0; qp < n_q_points; ++qp)
        {
          const double JxWxn = fe_values.JxW(qp) * (vec_vals[qp][0] * vec_vals[qp][0] + vec_vals[qp][1] * vec_vals[qp][1]);
          Point<iDim> spacept = fe_values.quadrature_point(qp);
          for (unsigned i = 0; i < iDim; i++)
          {
            oTotalIntegral[i] += JxWxn * (spacept[i] - pos[i]) * (spacept[i] - pos[i]);
          }
        }
      }
      std::array<double, iDim> retval;
      for (int i = 0; i < iDim; ++i)
      {
        retval[i] = oTotalIntegral[i];
      }
      return retval;
    }
    template std::array<double, 1> expectation_value_width<1>(IComplexWavefunction<1>*, const Vector<double>&, const Point<1>&);
    template std::array<double, 2> expectation_value_width<2>(IComplexWavefunction<2>*, const Vector<double>&, const Point<2>&);

    template <int iDim>
    std::array<double, iDim>
    expectation_value_width(IComplexWavefunction<iDim>* pBase, const LA::MPI::Vector& vec, const Point<iDim>& pos, MPI_Comm mpi_communicator)
    {
      const auto& fe = pBase->get_fe();
      const auto& dof_handler = pBase->get_dof_handler();

      assert(vec.has_ghost_elements() == true);

      Point<iDim> oLocalIntegral = {};

      const QGauss<iDim> quadrature_formula(fe.degree + 1);
      FEValues<iDim> fe_values(fe, quadrature_formula, update_values | update_quadrature_points | update_JxW_values);

      const unsigned n_q_points = quadrature_formula.size();
      std::vector<Vector<double>> vec_vals(n_q_points, Vector<double>(2));

      typename DoFHandler<iDim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
      for (; cell != endc; ++cell)
      {
        if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);
          fe_values.get_function_values(vec, vec_vals);
          for (unsigned qp = 0; qp < n_q_points; ++qp)
          {
            const double JxWxn = fe_values.JxW(qp) * (vec_vals[qp][0] * vec_vals[qp][0] + vec_vals[qp][1] * vec_vals[qp][1]);
            Point<iDim> spacept = fe_values.quadrature_point(qp);
            for (unsigned i = 0; i < iDim; i++)
            {
              oLocalIntegral[i] += JxWxn * (spacept[i] - pos[i]) * (spacept[i] - pos[i]);
            }
          }
        }
      }
      std::array<double, iDim> oTotalIntegral = {};
      MPI_Allreduce(oLocalIntegral.begin_raw(), oTotalIntegral.data(), iDim, MPI_DOUBLE, MPI_SUM, mpi_communicator);
      return oTotalIntegral;
    }
    template std::array<double, 2> expectation_value_width<2>(IComplexWavefunction<2>*, const LA::MPI::Vector&, const Point<2>&, MPI_Comm);
    template std::array<double, 3> expectation_value_width<3>(IComplexWavefunction<3>*, const LA::MPI::Vector&, const Point<3>&, MPI_Comm);
  }
}