
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
    using namespace dealii;

    template<int iDim>
    double particle_number
    (
      IRealWavefunction<iDim>* pBase,
      const Vector<double>& vWavefunction
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

    template double particle_number<1>(IBase<DoFHandler<1>, FE_Q<1>, AffineConstraints<double>>*, const Vector<double>&);
    template double particle_number<2>(IBase<DoFHandler<2>, FE_Q<2>, AffineConstraints<double>>*, const Vector<double>&);

    template<int iDim>
    double particle_number
    (
      IRealWavefunction<iDim>* pBase,
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
    std::array<double, iDim> expectation_value_position
    (
      IRealWavefunction<iDim>* pBase,
      const Vector<double>& vWavefunction,
      const std::optional<double> oParticleNumber
    )
    {
      Point<iDim> oTotalIntegral = {};

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
          oTotalIntegral += fe_values.JxW(qp) * vals[qp] * fe_values.quadrature_point(qp) * vals[qp];
        }
      }

      double rN = 1;
      if (oParticleNumber)
      {
        rN = oParticleNumber.value();
      }
      else
      {
        rN = particle_number(pBase, vWavefunction);
      }

      std::array<double, iDim> retval;
      for (int i = 0; i < iDim; ++i)
      {
        retval[i] = oTotalIntegral[i] / rN;
      }
      return retval;
    }

    template std::array<double, 1> expectation_value_position<1>(IBase<DoFHandler<1>, FE_Q<1>, AffineConstraints<double>>*, const Vector<double>&, const std::optional<double>);
    template std::array<double, 2> expectation_value_position<2>(IBase<DoFHandler<2>, FE_Q<2>, AffineConstraints<double>>*, const Vector<double>&, const std::optional<double>);

    template<int iDim>
    std::array<double, iDim> expectation_value_position
    (
      IRealWavefunction<iDim>* pBase,
      const Vector<double>& vWavefunction
    )
    {
      return expectation_value_position(pBase, vWavefunction, std::nullopt);
    }
    template std::array<double, 1> expectation_value_position<1>(IBase<DoFHandler<1>, FE_Q<1>, AffineConstraints<double>>*, const Vector<double>&);
    template std::array<double, 2> expectation_value_position<2>(IBase<DoFHandler<2>, FE_Q<2>, AffineConstraints<double>>*, const Vector<double>&);

    template<int iDim>
    std::array<double, iDim> expectation_value_position
    (
      IRealWavefunction<iDim>* pBase,
      const LA::MPI::Vector& vec,
      MPI_Comm mpi_communicator,
      std::optional<double> oParticleNumber
    )
    {
      assert(vec.has_ghost_elements() == true);

      const auto& fe = pBase->get_fe();
      const auto& dof_handler = pBase->get_dof_handler();

      Point<iDim> oLocalIntegral = {};

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
            oLocalIntegral += JxWxn * fe_values.quadrature_point(qp);
          }
        }
      }

      double rN = 1;
      if (oParticleNumber)
      {
        rN = oParticleNumber.value();
      }
      else
      {
        rN = particle_number(pBase, vec, mpi_communicator);
      }
      std::array<double, iDim> oTotalIntegral;
      MPI_Allreduce(oLocalIntegral.begin_raw(), oTotalIntegral.data(), iDim, MPI_DOUBLE, MPI_SUM, mpi_communicator);

      for (auto& oVal : oTotalIntegral)
      {
        oVal /= rN;
      }
      return oTotalIntegral;
    }

    template std::array<double, 2> expectation_value_position<2>(IRealWavefunction<2>*, const LA::MPI::Vector&, MPI_Comm, const std::optional<double>);
    template std::array<double, 3> expectation_value_position<3>(IRealWavefunction<3>*, const LA::MPI::Vector&, MPI_Comm, const std::optional<double>);

    template<int iDim>
    std::array<double, iDim> expectation_value_position
    (
      IRealWavefunction<iDim>* pBase,
      const LA::MPI::Vector& vec,
      MPI_Comm mpi_communicator
    )
    {
      return expectation_value_position(pBase, vec, mpi_communicator, std::nullopt);
    }
    template std::array<double, 2> expectation_value_position<2>(IRealWavefunction<2>*, const LA::MPI::Vector&, MPI_Comm);
    template std::array<double, 3> expectation_value_position<3>(IRealWavefunction<3>*, const LA::MPI::Vector&, MPI_Comm);

    template<int iDim>
    std::array<double, iDim> expectation_value_width
    (
      IRealWavefunction<iDim>* pBase,
      const Vector<double>& vWavefunction,
      const Point<iDim>& pos,
      std::optional<double> oParticleNumber
    )
    {
      Point<iDim> oResult = {};

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
            oResult[i] += JxWxn * ((spacept[i] - pos[i]) * (spacept[i] - pos[i]));
          }
        }
      }

      double rN = 1;
      if (oParticleNumber)
      {
        rN = oParticleNumber.value();
      }
      else
      {
        rN = particle_number(pBase, vWavefunction);
      }

      std::array<double, iDim> retval;
      for (int i = 0; i < iDim; ++i)
      {
        retval[i] = oResult[i] / rN;
      }
      return retval;
    }

    template std::array<double, 1> expectation_value_width<1>(IRealWavefunction<1>*, const Vector<double>&, const Point<1>&, const std::optional<double>);
    template std::array<double, 2> expectation_value_width<2>(IRealWavefunction<2>*, const Vector<double>&, const Point<2>&, const std::optional<double>);

    template<int iDim>
    std::array<double, iDim> expectation_value_width
    (
      IRealWavefunction<iDim>* pBase,
      const Vector<double>& vWavefunction,
      const Point<iDim>& pos
    )
    {
      return expectation_value_width(pBase, vWavefunction, pos, std::nullopt);
    }

    template std::array<double, 1> expectation_value_width<1>(IRealWavefunction<1>*, const Vector<double>&, const Point<1>&);
    template std::array<double, 2> expectation_value_width<2>(IRealWavefunction<2>*, const Vector<double>&, const Point<2>&);

    template<int iDim>
    std::array<double, iDim> expectation_value_width
    (
      IRealWavefunction<iDim>* pBase,
      const LA::MPI::Vector& vec,
      const Point<iDim>& pos,
      MPI_Comm mpi_communicator,
      std::optional<double> oParticleNumber
    )
    {
      assert(vec.has_ghost_elements() == true);

      const auto& fe = pBase->get_fe();
      const auto& dof_handler = pBase->get_dof_handler();

      Point<iDim> oLocalIntegral = {};

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
              oLocalIntegral[i] += JxWxn * (spacept[i] - pos[i]) * (spacept[i] - pos[i]);
            }
          }
        }
      }

      double rN = 1;
      if (oParticleNumber)
      {
        rN = oParticleNumber.value();
      }
      else
      {
        rN = particle_number(pBase, vec, mpi_communicator);
      }

      std::array<double, iDim> oTotalIntegral = {};
      MPI_Allreduce(oLocalIntegral.begin_raw(), oTotalIntegral.data(), iDim, MPI_DOUBLE, MPI_SUM, mpi_communicator);

      for (auto& oVal : oTotalIntegral)
      {
        oVal /= rN;
      }
      return oTotalIntegral;
    }

    template std::array<double, 2> expectation_value_width<2>(IRealWavefunction<2>*, const LA::MPI::Vector&, const Point<2>&, MPI_Comm, const std::optional<double>);
    template std::array<double, 3> expectation_value_width<3>(IRealWavefunction<3>*, const LA::MPI::Vector&, const Point<3>&, MPI_Comm, const std::optional<double>);

    template<int iDim>
    std::array<double, iDim> expectation_value_width
    (
      IRealWavefunction<iDim>* pBase,
      const LA::MPI::Vector& vec,
      const Point<iDim>& pos,
      MPI_Comm mpi_communicator
    )
    {
      return expectation_value_width(pBase, vec, pos, mpi_communicator);
    }

    template std::array<double, 2> expectation_value_width<2>(IRealWavefunction<2>*, const LA::MPI::Vector&, const Point<2>&, MPI_Comm);
    template std::array<double, 3> expectation_value_width<3>(IRealWavefunction<3>*, const LA::MPI::Vector&, const Point<3>&, MPI_Comm);
  }
}