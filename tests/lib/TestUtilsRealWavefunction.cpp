//
// atus-pro testing - atus-pro testing playgroung
// Copyright (C) 2022 Želimir Marojević <zelimir.marojevic@gmail.com>
//
// This file is part of atus-pro testing.
//
// atus-pro testing is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// atus-pro testing is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with atus-pro testing.  If not, see <http://www.gnu.org/licenses/>.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "UtilsRealWavefunction.hpp"

#include <iomanip>
#include <cmath>

template<typename T, int iDim>
class CMock : public utils::real_wavefunction::IRealWavefunction<iDim>
{
public:
  template<typename... Args>
  CMock(Args... args)
  : m_oTriangulation(std::forward<Args>(args)...)
  , m_oDofHandler(m_oTriangulation)
  {
    const dealii::Point<iDim, double> oGridCornerOne{ -10, -10, -10 };
    const dealii::Point<iDim, double> oGridCornerTwo{ 10, 10, 10 };
    dealii::GridGenerator::hyper_rectangle(m_oTriangulation, oGridCornerTwo, oGridCornerOne);
    m_oTriangulation.refine_global(10);
    m_oDofHandler.distribute_dofs(m_oFe);

    dealii::IndexSet oLocallyRelevantDofs;
    dealii::DoFTools::extract_locally_relevant_dofs(m_oDofHandler, oLocallyRelevantDofs);

    m_oConstraints.clear();
    m_oConstraints.reinit(oLocallyRelevantDofs);
    dealii::DoFTools::make_hanging_node_constraints(m_oDofHandler, m_oConstraints);
    dealii::VectorTools::interpolate_boundary_values(
      m_oDofHandler, 0, dealii::ZeroFunction<iDim>(), m_oConstraints);
    m_oConstraints.close();
  }
  dealii::DoFHandler<iDim>& get_dof_handler() { return m_oDofHandler; }

  dealii::FE_Q<iDim>& get_fe() { return m_oFe; }

  dealii::AffineConstraints<double>& get_constraints() { return m_oConstraints; }

protected:
  T                                 m_oTriangulation;
  dealii::DoFHandler<iDim>          m_oDofHandler;
  dealii::FE_Q<iDim>                m_oFe{ 2 };
  dealii::AffineConstraints<double> m_oConstraints;
};

TEST_CASE("particle_number", "[widget]")
{
  using namespace utils::real_wavefunction;

  SECTION("exp(-x^2)")
  {
    class CRealFunction : public dealii::Function<1>
    {
    public:
      virtual double value(const dealii::Point<1>& p, const unsigned int = 0) const
      {
        return exp(-p(0) * p(0));
      }
    } oFunction;

    CMock<dealii::Triangulation<1>, 1> oMock(dealii::Triangulation<1>::MeshSmoothing::none);
    dealii::Vector<double>             vPhi;
    vPhi.reinit(oMock.get_dof_handler().n_dofs());
    dealii::VectorTools::interpolate(oMock.get_dof_handler(), oFunction, vPhi);
    const auto rN = particle_number(dynamic_cast<IRealWavefunction<1>*>(&oMock), vPhi);

    INFO("rN = " << std::setprecision(15) << rN);
    CHECK(Catch::Matchers::WithinAbs(1.253314137, 1e-8).match(rN));
  }

  SECTION("exp(-x^2/2)/pi^{1/4}")
  {
    class CRealFunction : public dealii::Function<1>
    {
    public:
      virtual double value(const dealii::Point<1>& p, const unsigned int = 0) const
      {
        return 0.751125544464943 * exp(-0.5 * p(0) * p(0));
      }
    } oFunction;

    CMock<dealii::Triangulation<1>, 1> oMock(dealii::Triangulation<1>::MeshSmoothing::none);
    dealii::Vector<double>             vPhi;
    vPhi.reinit(oMock.get_dof_handler().n_dofs());
    dealii::VectorTools::interpolate(oMock.get_dof_handler(), oFunction, vPhi);
    const auto rN = particle_number(dynamic_cast<IRealWavefunction<1>*>(&oMock), vPhi);

    INFO("rN = " << std::setprecision(15) << rN);
    CHECK(Catch::Matchers::WithinAbs(1, 1e-8).match(rN));
  }

  SECTION("y*exp(-x^2-y^2)")
  {
    class CRealFunction : public dealii::Function<2>
    {
    public:
      virtual double value(const dealii::Point<2>& p, const unsigned int = 0) const
      {
        return p(1) * exp(-p(0) * p(0) - p(1) * p(1));
      }
    } oFunction;

    CMock<dealii::Triangulation<2>, 2> oMock(dealii::Triangulation<2>::MeshSmoothing::none);
    dealii::Vector<double>             vPhi;
    vPhi.reinit(oMock.get_dof_handler().n_dofs());
    dealii::VectorTools::interpolate(oMock.get_dof_handler(), oFunction, vPhi);
    const auto rN = particle_number(dynamic_cast<IRealWavefunction<2>*>(&oMock), vPhi);

    INFO("rN = " << std::setprecision(15) << rN);
    CHECK(Catch::Matchers::WithinAbs(0.392699081698727, 1e-8).match(rN));
  }
}

TEST_CASE("expectation_value_position", "[widget]")
{
  using namespace utils::real_wavefunction;
  namespace dpd = dealii::parallel::distributed;

  SECTION("exp(-x^2)")
  {
    class CRealFunction : public dealii::Function<1>
    {
    public:
      virtual double value(const dealii::Point<1>& p, const unsigned int = 0) const
      {
        return exp(-p(0) * p(0));
      }
    } oFunction;

    CMock<dealii::Triangulation<1>, 1> oMock(dealii::Triangulation<1>::MeshSmoothing::none);
    dealii::Vector<double>             vPhi;
    vPhi.reinit(oMock.get_dof_handler().n_dofs());
    dealii::VectorTools::interpolate(oMock.get_dof_handler(), oFunction, vPhi);
    const auto aPosition
      = expectation_value_position(dynamic_cast<IRealWavefunction<1>*>(&oMock), vPhi);

    INFO("<x> = " << std::setprecision(15) << aPosition[0]);
    CHECK(Catch::Matchers::WithinAbs(0, 1e-8).match(aPosition[0]));
  }

  SECTION("exp(-(x-1)^2-(y-2)^2)")
  {
    class CRealFunction : public dealii::Function<2>
    {
    public:
      virtual double value(const dealii::Point<2>& p, const unsigned int = 0) const
      {
        return exp(-std::pow(p(0) - 1, 2) - std::pow(p(1) - 2, 2));
      }
    } oFunction;

    CMock<dealii::Triangulation<2>, 2> oMock(dealii::Triangulation<2>::MeshSmoothing::none);
    dealii::Vector<double>             vPhi;
    vPhi.reinit(oMock.get_dof_handler().n_dofs());
    dealii::VectorTools::interpolate(oMock.get_dof_handler(), oFunction, vPhi);
    const auto aPosition
      = expectation_value_position(dynamic_cast<IRealWavefunction<2>*>(&oMock), vPhi);

    INFO("<x> = " << std::setprecision(15) << aPosition[0]);
    CHECK(Catch::Matchers::WithinAbs(1, 1e-8).match(aPosition[0]));
    INFO("<y> = " << std::setprecision(15) << aPosition[1]);
    CHECK(Catch::Matchers::WithinAbs(2, 1e-8).match(aPosition[1]));
  }

  SECTION("exp(-(x-1)^2-(y-2)^2)")
  {
    class CRealFunction : public dealii::Function<2>
    {
    public:
      virtual double value(const dealii::Point<2>& p, const unsigned int = 0) const
      {
        return exp(-std::pow(p(0) - 1, 2) - std::pow(p(1) - 2, 2));
      }
    } oFunction;

    MPI_Comm oComm{ MPI_COMM_WORLD };

    CMock<dpd::Triangulation<2>, 2> oMock(
      oComm, dpd::Triangulation<2>::MeshSmoothing(dpd::Triangulation<2>::smoothing_on_refinement));
    dealii::Vector<double> vPhi;
    vPhi.reinit(oMock.get_dof_handler().n_dofs());
    dealii::VectorTools::interpolate(oMock.get_dof_handler(), oFunction, vPhi);
    const auto aPosition
      = expectation_value_position(dynamic_cast<IRealWavefunction<2>*>(&oMock), vPhi);

    INFO("<x> = " << std::setprecision(15) << aPosition[0]);
    CHECK(Catch::Matchers::WithinAbs(1, 1e-8).match(aPosition[0]));
    INFO("<y> = " << std::setprecision(15) << aPosition[1]);
    CHECK(Catch::Matchers::WithinAbs(2, 1e-8).match(aPosition[1]));
  }
}

TEST_CASE("expectation_value_width", "[widget]")
{
  using namespace utils::real_wavefunction;

  SECTION("exp(-(x-1)^2)")
  {
    class CRealFunction : public dealii::Function<1>
    {
    public:
      virtual double value(const dealii::Point<1>& p, const unsigned int = 0) const
      {
        return exp(-std::pow(p(0) - 1, 2));
      }
    } oFunction;

    CMock<dealii::Triangulation<1>, 1> oMock(dealii::Triangulation<1>::MeshSmoothing::none);
    dealii::Vector<double>             vPhi;
    vPhi.reinit(oMock.get_dof_handler().n_dofs());
    dealii::VectorTools::interpolate(oMock.get_dof_handler(), oFunction, vPhi);
    const auto aPosition = expectation_value_width(
      dynamic_cast<IRealWavefunction<1>*>(&oMock), vPhi, dealii::Point<1>{ 1 });

    INFO("<x-<x>> = " << std::setprecision(15) << aPosition[0]);
    CHECK(Catch::Matchers::WithinAbs(0.25, 1e-8).match(aPosition[0]));
  }

  SECTION("exp(-(x-1)^2-(y-2)^2)")
  {
    class CRealFunction : public dealii::Function<2>
    {
    public:
      virtual double value(const dealii::Point<2>& p, const unsigned int = 0) const
      {
        return exp(-std::pow(p(0) - 1, 2) - 0.5 * std::pow(p(1), 2));
      }
    } oFunction;

    CMock<dealii::Triangulation<2>, 2> oMock(dealii::Triangulation<2>::MeshSmoothing::none);
    dealii::Vector<double>             vPhi;
    vPhi.reinit(oMock.get_dof_handler().n_dofs());
    dealii::VectorTools::interpolate(oMock.get_dof_handler(), oFunction, vPhi);
    const auto aPosition = expectation_value_width(
      dynamic_cast<IRealWavefunction<2>*>(&oMock), vPhi, dealii::Point<2>{ 1, 0 });

    INFO("<x-<x>> = " << std::setprecision(15) << aPosition[0]);
    CHECK(Catch::Matchers::WithinAbs(0.25, 1e-8).match(aPosition[0]));
    INFO("<y-<y>> = " << std::setprecision(15) << aPosition[1]);
    CHECK(Catch::Matchers::WithinAbs(0.5, 1e-8).match(aPosition[1]));
  }
}
