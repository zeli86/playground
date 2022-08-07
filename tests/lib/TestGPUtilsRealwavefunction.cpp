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
#include "GPUtilsRealWavefunction.hpp"

#include <iomanip>
#include <cmath>

template<int iDim>
class CMock : public utils::real_wavefunction::IRealWavefunction<iDim>
{
public:
  CMock() : m_oTriangulation(), m_oDofHandler(m_oTriangulation)
  {
    const dealii::Point<iDim, double> oGridCornerOne{-10, -10, -10};
    const dealii::Point<iDim, double> oGridCornerTwo{10, 10, 10};
    dealii::GridGenerator::hyper_rectangle(m_oTriangulation, oGridCornerTwo, oGridCornerOne);
    m_oTriangulation.refine_global(10);
    m_oDofHandler.distribute_dofs(m_oFe);

    dealii::IndexSet oLocallyRelevantDofs;
    dealii::DoFTools::extract_locally_relevant_dofs(m_oDofHandler, oLocallyRelevantDofs);

    m_oConstraints.clear();
    m_oConstraints.reinit(oLocallyRelevantDofs);
    dealii::DoFTools::make_hanging_node_constraints(m_oDofHandler, m_oConstraints);
    dealii::VectorTools::interpolate_boundary_values(m_oDofHandler, 0, dealii::ZeroFunction<iDim>(), m_oConstraints);
    m_oConstraints.close();
  }

  dealii::DoFHandler<iDim>& get_dof_handler()
  {
    return m_oDofHandler;
  }

  dealii::FE_Q<iDim>& get_fe()
  {
    return m_oFe;
  }

  dealii::AffineConstraints<double>& get_constraints()
  {
    return m_oConstraints;
  }
protected:
  dealii::Triangulation<iDim> m_oTriangulation;
  dealii::DoFHandler<iDim> m_oDofHandler;
  dealii::FE_Q<iDim> m_oFe{2};
  dealii::AffineConstraints<double> m_oConstraints;
};

class CPotential : public dealii::Function<1>
{
public:
  virtual double value(const dealii::Point<1>& p, const unsigned int = 0) const
  {
    return 0.25 * p(0) * p(0);
  }
} oPotential1d;

// class CPotential : public dealii::Function<2>
// {
// public:
//   virtual double value(const dealii::Point<2>& p, const unsigned int = 0) const
//   {
//     return 0.25 * (p(0) * p(0) + p(1) * p(1));
//   }
// } oPotential2d;


TEST_CASE("GP", "[widget]")
{
  using namespace  utils::real_wavefunction;

  SECTION("pi^{-1/4}*exp(-0.5*x^2)")
  {
    class CTestFunction : public dealii::Function<1>
    {
    public:
      virtual double value(const dealii::Point<1>& p, const unsigned int = 0) const
      {
        return 0.751125544464943 * exp(-0.5 * p(0) * p(0));
      }
    } oFunction;

    CMock<1> oMock;
    dealii::Vector<double> vPhi;
    vPhi.reinit(oMock.get_dof_handler().n_dofs());
    dealii::VectorTools::interpolate(oMock.get_dof_handler(), oFunction, vPhi);
    const auto tpContributions = GP(dynamic_cast<IRealWavefunction<1>*>(&oMock), vPhi, oPotential1d);

    INFO("Ekin = " << std::setprecision(15) << std::get<0>(tpContributions));
    INFO("Epot = " << std::setprecision(15) << std::get<1>(tpContributions));
    INFO("Eint = " << std::setprecision(15) << std::get<2>(tpContributions));
    CHECK(Catch::Matchers::WithinAbs(0.625, 1e-8).match(std::get<0>(tpContributions)));
    CHECK(Catch::Matchers::WithinAbs(1, 1e-8).match(std::get<1>(tpContributions)));
    CHECK(Catch::Matchers::WithinAbs(0.398942280401433, 1e-8).match(std::get<2>(tpContributions)));
  }
}

TEST_CASE("mu", "[widget]")
{
  using namespace  utils::real_wavefunction;

  SECTION("pi^{-1/4}*exp(-0.5*x^2)")
  {
    class CTestFunction : public dealii::Function<1>
    {
    public:
      virtual double value(const dealii::Point<1>& p, const unsigned int = 0) const
      {
        return 0.751125544464943 * exp(-0.5 * p(0) * p(0));
      }
    } oFunction;

    CMock<1> oMock;
    dealii::Vector<double> vPhi;
    vPhi.reinit(oMock.get_dof_handler().n_dofs());
    dealii::VectorTools::interpolate(oMock.get_dof_handler(), oFunction, vPhi);
    const auto rMu = MU(dynamic_cast<IRealWavefunction<1>*>(&oMock), vPhi, oPotential1d);

    INFO("mu = " << std::setprecision(15) << rMu);
    CHECK(Catch::Matchers::WithinAbs(1, 1e-8).match(rMu));
  }
}