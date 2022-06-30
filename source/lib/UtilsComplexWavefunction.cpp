
#include "UtilsComplexWavefunction.hpp"
#include "mpi.h"

#include <deal.II/base/quadrature_lib.h>

namespace utils
{
  namespace complex_wavefunction
  {
    template <int dim>
    double
    particle_number
    (
      IComplexWavefunction<dim>* pBase,
      const Vector<double>& vec
    )
    {
      const auto& fe = pBase->get_fe();
      const auto& dof_handler = pBase->get_dof_handler();

      double retval{0};
      const QGauss<dim> quadrature_formula(fe.degree + 1);

      FEValues<dim> fe_values(fe, quadrature_formula, update_values | update_quadrature_points | update_JxW_values);

      const unsigned n_q_points = quadrature_formula.size();

      std::vector<Vector<double>> vals(n_q_points, Vector<double>(2));

      typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
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

    template <int dim>
    double
    particle_number
    (
      IComplexWavefunction<dim>* pBase,
      const LA::MPI::Vector& vec,
      MPI_Comm mpi_communicator
    )
    {
      assert(vec.has_ghost_elements() == true);
      const auto& fe = pBase->get_fe();
      const auto& dof_handler = pBase->get_dof_handler();
      double tmp = 0;

      const QGauss<dim> quadrature_formula(fe.degree + 1);
      FEValues<dim> fe_values(fe, quadrature_formula, update_values | update_JxW_values);

      const unsigned n_q_points = quadrature_formula.size();
      std::vector<Vector<double>> vec_vals(n_q_points, Vector<double>(2));

      typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
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
      return Utilities::MPI::sum(tmp, mpi_communicator);;
    }
    template double particle_number<2>(IComplexWavefunction<2>*, const LA::MPI::Vector&, MPI_Comm);
    template double particle_number<3>(IComplexWavefunction<3>*, const LA::MPI::Vector&, MPI_Comm);

    template <int dim>
    Point<dim>
    expectation_value_position
    (
      IComplexWavefunction<dim>* pBase,
      const Vector<double>& vec
    )
    {
      const auto& fe = pBase->get_fe();
      const auto& dof_handler = pBase->get_dof_handler();

      const QGauss<dim> quadrature_formula(fe.degree + 1);

      FEValues<dim> fe_values(fe, quadrature_formula, update_values | update_quadrature_points | update_JxW_values);

      const unsigned dofs_per_cell = fe.dofs_per_cell;
      const unsigned n_q_points = quadrature_formula.size();

      std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
      std::vector<Vector<double>> vals(n_q_points, Vector<double>(2));

      Point<dim> retval = {};

      typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
      for (; cell != endc; ++cell)
      {
        fe_values.reinit(cell);
        fe_values.get_function_values(vec, vals);

        for (unsigned qp = 0; qp < n_q_points; ++qp)
        {
          double JxWxn = fe_values.JxW(qp) * (vals[qp][0] * vals[qp][0] + vals[qp][1] * vals[qp][1]);
          Point<dim> spacept = fe_values.quadrature_point(qp);
          retval += JxWxn * spacept;
        }
      }
      return retval;
    }
    template Point<1> expectation_value_position(IComplexWavefunction<1>*, const Vector<double>&);
    template Point<2> expectation_value_position(IComplexWavefunction<2>*, const Vector<double>&);

    template <int dim>
    Point<dim>
    expectation_value_position
    (
      IComplexWavefunction<dim>* pBase,
      const LA::MPI::Vector& vec,
      MPI_Comm mpi_communicator
    )
    {
      assert(vec.has_ghost_elements() == true);
      const auto& fe = pBase->get_fe();
      const auto& dof_handler = pBase->get_dof_handler();

      Point<dim> tmp = {};

      const QGauss<dim> quadrature_formula(fe.degree + 1);
      FEValues<dim> fe_values(fe, quadrature_formula, update_values | update_quadrature_points | update_JxW_values);

      const unsigned n_q_points = quadrature_formula.size();
      std::vector<Vector<double>> vec_vals(n_q_points, Vector<double>(2));

      typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
      for (; cell != endc; ++cell)
      {
        if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);
          fe_values.get_function_values(vec, vec_vals);
          for (unsigned qp = 0; qp < n_q_points; ++qp)
          {
            const double JxWxn = fe_values.JxW(qp) * (vec_vals[qp][0] * vec_vals[qp][0] + vec_vals[qp][1] * vec_vals[qp][1]);
            Point<dim> spacept = fe_values.quadrature_point(qp);
            tmp += JxWxn * spacept;
          }
        }
      }
      return Utilities::MPI::sum(tmp, mpi_communicator);
    }
    template Point<2> expectation_value_position<2>(IComplexWavefunction<2>*, const LA::MPI::Vector&, MPI_Comm);
    template Point<3> expectation_value_position<3>(IComplexWavefunction<3>*, const LA::MPI::Vector&, MPI_Comm);

    template <int dim>
    Point<dim>
    expectation_value_momentum(IComplexWavefunction<dim>* pBase, const Vector<double>& vec)
    {
      const auto& fe = pBase->get_fe();
      const auto& dof_handler = pBase->get_dof_handler();

      Point<dim >retval = {};

      Point<dim> tmp;

      const QGauss<dim> quadrature_formula(fe.degree + 1);
      FEValues<dim> fe_values(fe, quadrature_formula, update_values | update_gradients | update_JxW_values);

      const unsigned n_q_points = quadrature_formula.size();
      std::vector<Vector<double>> vec_vals(n_q_points, Vector<double>(2));
      std::vector<Vector<Tensor<1, dim>>> vec_grads(n_q_points, Vector<Tensor<1, dim>>(2));

      typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
      for (; cell != endc; ++cell)
      {
        fe_values.reinit(cell);
        fe_values.get_function_values(vec, vec_vals);
        for (unsigned qp = 0; qp < n_q_points; ++qp)
        {
          double JxW = fe_values.JxW(qp);
          for (unsigned i = 0; i < dim; ++i)
          {
            retval[i] += JxW * (vec_vals[qp][0] * vec_grads[qp][1][i] - vec_vals[qp][1] * vec_grads[qp][0][i]);
          }
        }
      }
      return retval;
    }
    template Point<1> expectation_value_momentum(IComplexWavefunction<1>*, const Vector<double>&);
    template Point<2> expectation_value_momentum(IComplexWavefunction<2>*, const Vector<double>&);

    template <int dim>
    Point<dim>
    expectation_value_momentum
    (
      IComplexWavefunction<dim>* pBase,
      const LA::MPI::Vector& vec,
      MPI_Comm mpi_communicator
    )
    {
      const auto& fe = pBase->get_fe();
      const auto& dof_handler = pBase->get_dof_handler();

      assert(vec.has_ghost_elements() == true);

      Point<dim> tmp = {};

      const QGauss<dim> quadrature_formula(fe.degree + 1);
      FEValues<dim> fe_values(fe, quadrature_formula, update_values | update_gradients | update_JxW_values);

      const unsigned n_q_points = quadrature_formula.size();
      std::vector<Vector<double>> vec_vals(n_q_points, Vector<double>(2));
      std::vector<Vector<Tensor<1, dim>>> vec_grads(n_q_points, Vector<Tensor<1, dim>>(2));

      typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
      for (; cell != endc; ++cell)
      {
        if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);
          fe_values.get_function_values(vec, vec_vals);
          for (unsigned qp = 0; qp < n_q_points; qp++)
          {
            const double JxW = fe_values.JxW(qp);
            for (unsigned i = 0; i < dim; ++i)
            {
              tmp[i] += JxW * (vec_vals[qp][0] * vec_grads[qp][1][i] - vec_vals[qp][1] * vec_grads[qp][0][i]);
            }
          }
        }
      }
      return Utilities::MPI::sum(tmp, mpi_communicator);
    }
    template Point<2> expectation_value_momentum<2>(IComplexWavefunction<2>*, const LA::MPI::Vector&, MPI_Comm);
    template Point<3> expectation_value_momentum<3>(IComplexWavefunction<3>*, const LA::MPI::Vector&, MPI_Comm);

    template <int dim>
    Point<dim> expectation_value_width
    (
      IComplexWavefunction<dim>* pBase,
      const Vector<double>& vec,
      const Point<dim>& pos
    )
    {
      const auto& fe = pBase->get_fe();
      const auto& dof_handler = pBase->get_dof_handler();

      assert(vec.has_ghost_elements() == true);

      Point<dim> tmp = {};

      const QGauss<dim> quadrature_formula(fe.degree + 1);
      FEValues<dim> fe_values(fe, quadrature_formula, update_values | update_quadrature_points | update_JxW_values);

      const unsigned n_q_points = quadrature_formula.size();
      std::vector<Vector<double>> vec_vals(n_q_points, Vector<double>(2));

      typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
      for (; cell != endc; ++cell)
      {
        fe_values.reinit(cell);
        fe_values.get_function_values(vec, vec_vals);
        for (unsigned qp = 0; qp < n_q_points; ++qp)
        {
          const double JxWxn = fe_values.JxW(qp) * (vec_vals[qp][0] * vec_vals[qp][0] + vec_vals[qp][1] * vec_vals[qp][1]);
          Point<dim> spacept = fe_values.quadrature_point(qp);
          for (unsigned i = 0; i < dim; i++)
          {
            tmp[i] += JxWxn * (spacept[i] - pos[i]) * (spacept[i] - pos[i]);
          }
        }
      }
      return tmp;
    }
    template Point<1> expectation_value_width<1>(IComplexWavefunction<1>*, const Vector<double>&, const Point<1>&);
    template Point<2> expectation_value_width<2>(IComplexWavefunction<2>*, const Vector<double>&, const Point<2>&);

    template <int dim>
    Point<dim>
    expectation_value_width(IComplexWavefunction<dim>* pBase, const LA::MPI::Vector& vec, const Point<dim>& pos, MPI_Comm mpi_communicator)
    {
      const auto& fe = pBase->get_fe();
      const auto& dof_handler = pBase->get_dof_handler();

      assert(vec.has_ghost_elements() == true);

      Point<dim> tmp = {};

      const QGauss<dim> quadrature_formula(fe.degree + 1);
      FEValues<dim> fe_values(fe, quadrature_formula, update_values | update_quadrature_points | update_JxW_values);

      const unsigned n_q_points = quadrature_formula.size();
      std::vector<Vector<double>> vec_vals(n_q_points, Vector<double>(2));

      typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
      for (; cell != endc; ++cell)
      {
        if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);
          fe_values.get_function_values(vec, vec_vals);
          for (unsigned qp = 0; qp < n_q_points; ++qp)
          {
            const double JxWxn = fe_values.JxW(qp) * (vec_vals[qp][0] * vec_vals[qp][0] + vec_vals[qp][1] * vec_vals[qp][1]);
            Point<dim> spacept = fe_values.quadrature_point(qp);
            for (unsigned i = 0; i < dim; i++)
            {
              tmp[i] += JxWxn * (spacept[i] - pos[i]) * (spacept[i] - pos[i]);
            }
          }
        }
      }
      return Utilities::MPI::sum(tmp, mpi_communicator);
    }
    template Point<2> expectation_value_width<2>(IComplexWavefunction<2>*, const LA::MPI::Vector&, const Point<2>&, MPI_Comm);
    template Point<3> expectation_value_width<3>(IComplexWavefunction<3>*, const LA::MPI::Vector&, const Point<3>&, MPI_Comm);
  }
}