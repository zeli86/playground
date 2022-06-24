
#include "UtilsRealWavefunction.hpp"
#include "deal.II/base/point.h"
#include "deal.II/base/quadrature_lib.h"
#include "deal.II/lac/generic_linear_algebra.h"
#include "mpi.h"
#include <array>
#include <vector>


namespace utils
{
  namespace real_wavefunction
  {
    template<int iDim>
    double particle_number
    (
      IBaseRealWavefunction<iDim>* pBase,
      const std::vector<double>& vWavefunction
    )
    {
      double retval{0};

      const auto& fe = pBase->get_fe();
      const auto& dof_handler = pBase->get_dof_handler();

      const dealii::QGauss<iDim> quadrature_formula(fe.degree + 1);

      FEValues<iDim> fe_values(fe, quadrature_formula, update_values | update_quadrature_points | update_JxW_values);

      const unsigned n_q_points = quadrature_formula.size();

      std::vector<double> vals(n_q_points);

      typename DoFHandler<iDim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
      for (; cell != endc; ++cell)
      {
        fe_values.reinit(cell);
        fe_values.get_function_values(vWavefunction, vals);

        for (unsigned qp = 0; qp < n_q_points; qp++)
        {
          retval += fe_values.JxW(qp) * vals[qp] * vals[qp];
        }
      }
      return retval;
    }

    template double particle_number<1>(IBase<DoFHandler<1>, FE_Q<1>, AffineConstraints<double>>*, const std::vector<double>&);
    template double particle_number<2>(IBase<DoFHandler<2>, FE_Q<2>, AffineConstraints<double>>*, const std::vector<double>&);

    template<int iDim>
    double particle_number
    (
      IBaseRealWavefunction<iDim>* pBase,
      const LA::MPI::Vector& vWavefunction,
      MPI_Comm mpi_communicator
    )
    {
      assert(vWavefunction.has_ghost_elements() == true);

      const auto& fe = pBase->get_fe();
      const auto& dof_handler = pBase->get_dof_handler();

      double tmp = 0;

      const QGauss<iDim>  quadrature_formula(fe.degree + 1);
      FEValues<iDim> fe_values(fe, quadrature_formula, update_values | update_JxW_values);

      const unsigned n_q_points = quadrature_formula.size();
      std::vector<double> vec_vals(n_q_points);

      typename DoFHandler<iDim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
      for (; cell != endc; ++cell)
      {
        if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);
          fe_values.get_function_values(vWavefunction, vec_vals);
          for (unsigned qp = 0; qp < n_q_points; ++qp)
          {
            tmp += fe_values.JxW(qp) * vec_vals[qp] * vec_vals[qp];
          }
        }
      }
      return Utilities::MPI::sum(tmp, mpi_communicator);
    }

    template double particle_number<2>(IBase<DoFHandler<2>, FE_Q<2>, AffineConstraints<double>>*, const LA::MPI::Vector&, MPI_Comm);
    template double particle_number<3>(IBase<DoFHandler<3>, FE_Q<3>, AffineConstraints<double>>*, const LA::MPI::Vector&, MPI_Comm);

    template<int iDim>
    Point<iDim> expectation_value_position
    (
      IBaseRealWavefunction<iDim>* pBase,
      const std::vector<double>& vWavefunction
    )
    {
      Point<iDim> retval = {};

      const auto& fe = pBase->get_fe();
      const auto& dof_handler = pBase->get_dof_handler();

      const dealii::QGauss<iDim> quadrature_formula(fe.degree + 1);

      FEValues<iDim> fe_values(fe, quadrature_formula, update_values | update_quadrature_points | update_JxW_values);

      const unsigned n_q_points = quadrature_formula.size();

      std::vector<double> vals(n_q_points);

      typename DoFHandler<iDim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
      for (; cell != endc; ++cell)
      {
        fe_values.reinit(cell);
        fe_values.get_function_values(vWavefunction, vals);

        for (unsigned qp = 0; qp < n_q_points; qp++)
        {
          retval += fe_values.JxW(qp) * vals[qp] * fe_values.quadrature_point(qp) * vals[qp];
        }
      }
      return retval;
    }

    template Point<1> expectation_value_position<1>(IBase<DoFHandler<1>, FE_Q<1>, AffineConstraints<double>>*, const std::vector<double>&);
    template Point<2> expectation_value_position<2>(IBase<DoFHandler<2>, FE_Q<2>, AffineConstraints<double>>*, const std::vector<double>&);

    template<int iDim>
    Point<iDim> expectation_value_position
    (
      IBaseRealWavefunction<iDim>* pBase,
      const LA::MPI::Vector& vec,
      MPI_Comm mpi_communicator
    )
    {
      assert(vec.has_ghost_elements() == true);

      const auto& fe = pBase->get_fe();
      const auto& dof_handler = pBase->get_dof_handler();

      Point<iDim> retval = {};

      const QGauss<iDim>  quadrature_formula(fe.degree + 1);
      FEValues<iDim> fe_values(fe, quadrature_formula, update_values | update_quadrature_points | update_JxW_values);

      const unsigned n_q_points = quadrature_formula.size();
      std::vector<double> vec_vals(n_q_points);

      typename DoFHandler<iDim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
      for (; cell != endc; ++cell)
      {
        if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);
          fe_values.get_function_values(vec, vec_vals);
          for (unsigned qp = 0; qp < n_q_points; ++qp)
          {
            double JxWxn = fe_values.JxW(qp) * vec_vals[qp] * vec_vals[qp];
            retval += JxWxn * fe_values.quadrature_point(qp);
          }
        }
      }
      return Utilities::MPI::sum(retval, mpi_communicator);
    }

    template Point<2> expectation_value_position<2>(IBaseRealWavefunction<2>*, const LA::MPI::Vector&, MPI_Comm);
    template Point<3> expectation_value_position<3>(IBaseRealWavefunction<3>*, const LA::MPI::Vector&, MPI_Comm);

    template<int iDim>
    Point<iDim> expectation_value_width
    (
      IBaseRealWavefunction<iDim>* pBase,
      const std::vector<double>& vWavefunction,
      const Point<iDim>& pos
    )
    {
      Point<iDim> retval = {};

      const auto& fe = pBase->get_fe();
      const auto& dof_handler = pBase->get_dof_handler();

      const dealii::QGauss<iDim> quadrature_formula(fe.degree + 1);

      FEValues<iDim> fe_values(fe, quadrature_formula, update_values | update_quadrature_points | update_JxW_values);

      const unsigned n_q_points = quadrature_formula.size();

      std::vector<double> vals(n_q_points);

      typename DoFHandler<iDim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
      for (; cell != endc; ++cell)
      {
        fe_values.reinit(cell);
        fe_values.get_function_values(vWavefunction, vals);
        for (unsigned qp = 0; qp < n_q_points; qp++)
        {
          const double JxWxn = fe_values.JxW(qp) * vals[qp] * vals[qp];
          auto spacept = fe_values.quadrature_point(qp);
          for (unsigned i = 0; i < iDim; i++)
          {
            retval[i] += JxWxn * ((spacept[i] - pos[i]) * (spacept[i] - pos[i]));
          }
        }
      }
      return retval;
    }

    template Point<1> expectation_value_width<1>(IBaseRealWavefunction<1>*, const std::vector<double>&, const Point<1>&);
    template Point<2> expectation_value_width<2>(IBaseRealWavefunction<2>*, const std::vector<double>&, const Point<2>&);

    template<int iDim>
    Point<iDim> expectation_value_width
    (
      IBaseRealWavefunction<iDim>* pBase,
      const LA::MPI::Vector& vec,
      const Point<iDim>& pos,
      MPI_Comm mpi_communicator
    )
    {
      assert(vec.has_ghost_elements() == true);

      const auto& fe = pBase->get_fe();
      const auto& dof_handler = pBase->get_dof_handler();

      Point<iDim> tmp = {};

      const QGauss<iDim>  quadrature_formula(fe.degree + 1);
      FEValues<iDim> fe_values(fe, quadrature_formula, update_values | update_quadrature_points | update_JxW_values);

      const unsigned n_q_points = quadrature_formula.size();
      std::vector<double> vec_vals(n_q_points);

      typename DoFHandler<iDim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
      for (; cell != endc; ++cell)
      {
        if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);
          fe_values.get_function_values(vec, vec_vals);
          for (unsigned qp = 0; qp < n_q_points; qp++)
          {
            const double JxWxn = fe_values.JxW(qp) * vec_vals[qp] * vec_vals[qp];
            auto spacept = fe_values.quadrature_point(qp);
            for (unsigned i = 0; i < iDim; i++)
            {
              tmp[i] += JxWxn * (spacept[i] - pos[i]) * (spacept[i] - pos[i]);
            }
          }
        }
      }
      return Utilities::MPI::sum(tmp, mpi_communicator);
    }

    template Point<2> expectation_value_width<2>(IBaseRealWavefunction<2>*, const LA::MPI::Vector&, const Point<2>&, MPI_Comm);
    template Point<3> expectation_value_width<3>(IBaseRealWavefunction<3>*, const LA::MPI::Vector&, const Point<3>&, MPI_Comm);
  }
}