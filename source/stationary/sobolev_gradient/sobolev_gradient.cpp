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
//

/*
    Želimir Marojević, Ertan Göklü, Claus Lämmerzahl
    (ZARM - Center of Applied Space Technology and Microgravity, Germany, http://www.zarm.uni-bremen.de/)

    Public use and modification of this code are allowed provided that the following papers are cited:

    Marojević, Želimir, Ertan Göklü, und Claus Lämmerzahl. "ATUS-PRO: A FEM-based solver for the time-dependent and stationary Gross–Pitaevskii equation",
    Computer Physics Communications, Vol 202, 2016, p. 216--232. doi:10.1016/j.cpc.2015.12.004.

    W. Bangerth and D. Davydov and T. Heister and L. Heltai and G. Kanschat and M. Kronbichler and M. Maier and B. Turcksin and D. Wells.
    "The \texttt{deal.II} Library, Version 8.4", Journal of Numerical Mathematics, vol 24, 2016.

    The authors would be grateful for all information and/or comments regarding the use of the code.
*/

// #include <deal.II/lac/dynamic_sparsity_pattern.h>
// #include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_direct.h>
// #include <deal.II/lac/sparsity_tools.h>
// #include <deal.II/lac/vector.h>

// #include <deal.II/base/function.h>
// #include <deal.II/base/function_parser.h>
// #include <deal.II/base/quadrature_lib.h>
// #include <deal.II/base/timer.h>
// #include <deal.II/base/utilities.h>

// #include <deal.II/dofs/dof_accessor.h>
// #include <deal.II/dofs/dof_handler.h>
// #include <deal.II/dofs/dof_tools.h>

// #include <deal.II/fe/fe_q.h>
// #include <deal.II/fe/fe_values.h>

// #include <deal.II/grid/grid_generator.h>
// #include <deal.II/grid/tria_accessor.h>
// #include <deal.II/grid/tria_iterator.h>

#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>
// #include <deal.II/numerics/data_out.h>

#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <cmath>
#include <limits>

#include "functions.h"
#include "sobolev_gradient.hpp"

namespace solver::stationary
{
  using namespace std;
  using namespace dealii;

  enum Status { SUCCESS, FAILED, ZERO_SOL, SLOW_CONV };

  template <int dim>
  CSobolevGradient<dim>::CSobolevGradient(const std::string& xmlfilename) :
    m_oParameters(xmlfilename),
    triangulation(),
    dof_handler(triangulation)
  {
  }

  template <int dim>
  CSobolevGradient<dim>::~CSobolevGradient()
  {
    dof_handler.clear();
  }

  template <int dim>
  void CSobolevGradient<dim>::make_grid()
  {
    Point<dim, double> pt1;
    Point<dim, double> pt2;

    std::vector<double> oGridCornerOne{0, 0, 0};
    std::vector<double> oGridCornerTwo{0, 0, 0};

    m_oParameters.get("grid.grid_corner_one", oGridCornerOne);
    m_oParameters.get("grid.grid_corber_two", oGridCornerTwo);

    for (int i = 0; i < dim; i++)
    {
      pt1(i) = oGridCornerOne[i];
      pt2(i) = oGridCornerTwo[i];
    }

    GridGenerator::hyper_rectangle(triangulation, pt2, pt1);
    triangulation.refine_global(m_iGlobalRefinement);
  }

  template <int dim>
  void CSobolevGradient<dim>::setup_system()
  {
    dof_handler.distribute_dofs(fe);

    m_vPhi.reinit(dof_handler.n_dofs());
    m_vPhiSobolev.reinit(dof_handler.n_dofs());
    m_vL2Gradient.reinit(dof_handler.n_dofs());
    m_solution.reinit(dof_handler.n_dofs());
    m_vSobolevGradient.reinit(dof_handler.n_dofs());
    m_error_per_cell.reinit(triangulation.n_active_cells());

    m_oConstraints.clear();
    m_oConstraints.reinit(m_oLocallyRelevantDofs);
    dealii::DoFTools::make_hanging_node_constraints(dof_handler, m_oConstraints);
    dealii::VectorTools::interpolate_boundary_values(dof_handler, 0, ZeroFunction<dim>(), m_oConstraints);
    m_oConstraints.close();

    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, dsp);
    m_sparsity_pattern.copy_from(dsp);
    m_system_matrix.reinit(m_sparsity_pattern);
  }

  template <int dim>
  void CSobolevGradient<dim>::project_gradient()
  {
    const QGauss<dim> quadrature_formula(fe.degree + 1);
    FEValues<dim> fe_values(fe, quadrature_formula, update_values | update_JxW_values | update_quadrature_points);

    const unsigned n_q_points = quadrature_formula.size();
    vector<double> vals_1(n_q_points);
    vector<double> vals_2(n_q_points);
    vector<double> vals_3(n_q_points);

    double sum[] = {0, 0};

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell != endc; ++cell)
    {
      fe_values.reinit(cell);
      fe_values.get_function_values(m_vPhiSobolev, vals_1);
      fe_values.get_function_values(m_vSobolevGradient, vals_2);
      fe_values.get_function_values(m_vPhi, vals_3);

      for (unsigned qp = 0; qp < n_q_points; ++qp)
      {
        const double JxW = fe_values.JxW(qp);
        sum[0] += JxW * (vals_3[qp] * vals_2[qp]);
        sum[1] += JxW * (vals_3[qp] * vals_1[qp]);
      }
    }
    m_vSobolevGradient.add(-sum[0] / sum[1], m_vPhiSobolev);
  }


  template <int dim>
  void CSobolevGradient<dim>::assemble_system()
  {
    const QGauss<dim> quadrature_formula(fe.degree + 1);

    m_system_matrix = 0;

    FEValues<dim> fe_values(fe, quadrature_formula, update_values | update_gradients | update_JxW_values | update_quadrature_points);

    const auto dofs_per_cell = fe.dofs_per_cell;
    const auto n_q_points = quadrature_formula.size();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell != endc; ++cell)
    {
      cell_matrix = 0;
      fe_values.reinit(cell);

      for (unsigned qp = 0; qp < n_q_points; ++qp)
      {
        const double JxW = fe_values.JxW(qp);

        for (unsigned i = 0; i < dofs_per_cell; ++i)
          for (unsigned j = 0; j < dofs_per_cell; ++j)
          {
            cell_matrix(i, j) += JxW * (fe_values.shape_grad(i, qp) * fe_values.shape_grad(j, qp) + fe_values.shape_value(i, qp) * fe_values.shape_value(j, qp));
          }
      }
      cell->get_dof_indices(local_dof_indices);
      for (unsigned i = 0; i < dofs_per_cell; ++i)
        for (unsigned j = 0; j < dofs_per_cell; ++j)
        {
          m_system_matrix.add(local_dof_indices[i], local_dof_indices[j], cell_matrix(i, j));
        }
    }
  }

  template <int dim>
  void CSobolevGradient<dim>::compute_Psi_sob()
  {
    m_vL2Gradient = 0;

    CPotential<dim> Potential(m_omega);
    const QGauss<dim> quadrature_formula(fe.degree + 1);

    FEValues<dim> fe_values(fe, quadrature_formula, update_values | update_gradients | update_JxW_values);

    const unsigned dofs_per_cell = fe.dofs_per_cell;
    const unsigned n_q_points = quadrature_formula.size();

    vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
    vector<double> vals(n_q_points);
    Vector<double> cell_rhs(dofs_per_cell);
    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell != endc; ++cell)
    {
      cell_rhs = 0;
      cell_matrix = 0;
      fe_values.reinit(cell);
      fe_values.get_function_values(m_vPhi, vals);
      cell->get_dof_indices(local_dof_indices);

      for (unsigned qp = 0; qp < n_q_points; ++qp)
      {
        const double JxW = fe_values.JxW(qp);

        for (unsigned i = 0; i < dofs_per_cell; ++i)
        {
          cell_rhs(i) += JxW * vals[qp] * fe_values.shape_value(i, qp);
          for (unsigned j = 0; j < dofs_per_cell; ++j)
          {
            cell_matrix(i, j) += JxW * (fe_values.shape_grad(i, qp) * fe_values.shape_grad(j, qp) + fe_values.shape_value(i, qp) * fe_values.shape_value(j, qp));
          }
        }
      }
      cell->get_dof_indices(local_dof_indices);
      for (unsigned i = 0; i < dofs_per_cell; ++i)
      {
        m_vL2Gradient(local_dof_indices[i]) += cell_rhs(i);
        for (unsigned j = 0; j < dofs_per_cell; ++j)
        {
          m_system_matrix.add(local_dof_indices[i], local_dof_indices[j], cell_matrix(i, j));
        }
      }
    }
    solve();
    m_vPhiSobolev = m_solution;
  }

  template <int dim>
  void CSobolevGradient<dim>::solve()
  {
    map<types::global_dof_index, double> boundary_values;
    VectorTools::interpolate_boundary_values(dof_handler, 0, ZeroFunction<dim>(), boundary_values);
    MatrixTools::apply_boundary_values(boundary_values, m_system_matrix, m_solution, m_vL2Gradient);

    SparseDirectUMFPACK A_direct;
    A_direct.initialize(m_system_matrix);
    A_direct.vmult(m_solution, m_vL2Gradient);
  }


  template <int dim>
  void CSobolevGradient<dim>::output_results(string filename)
  {
    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(m_vPhi, "Psi");
    data_out.add_data_vector(m_error_per_cell, "error_per_cell");
    data_out.build_patches();

    ofstream output(filename.c_str());
    data_out.write_gnuplot(output);
  }

  template <int dim>
  int CSobolevGradient<dim>::DoIter(string path)
  {
    using namespace utils::real_wavefunction;
    int retval = Status::SUCCESS;

    CPotential<dim> oPotential(m_omega);
    assemble_L2gradient(dynamic_cast<IRealWavefunction<dim>*>(this), m_vPhi, oPotential, 0, m_rG, m_vL2Gradient);

    m_res = 0;
    do
    {
      assemble_system();
      solve();
      m_vSobolevGradient = m_solution;

      compute_Psi_sob();
      project_gradient();
      m_res = m_vSobolevGradient.l2_norm();

      m_vPhi.add(-0.1, m_vSobolevGradient);

      // compute_mu(m_rMu);
      m_N = particle_number(dynamic_cast<IRealWavefunction<dim>*>(this), m_vPhi);

      if (fabs(m_N - 1) > 1e-5)
      {
        m_vPhi *= 1 / sqrt(m_N);
      }

      assemble_L2gradient(dynamic_cast<IRealWavefunction<dim>*>(this), m_vPhi, oPotential, 0, m_rG, m_vL2Gradient);

      if (m_res < m_rEpsilon)
      {
        retval = Status::SUCCESS;
        break;
      }

    }
    while (true);

    if (m_N < 1e-5)
    {
      retval = Status::ZERO_SOL;
    }

    return retval;
  }

  template <int dim>
  void CSobolevGradient<dim>::run()
  {
    using namespace utils::real_wavefunction;

    int status;

    make_grid();
    setup_system();

    CPotential<dim> oPotential(m_omega);

    unsigned QN[] = {0, 0, 0};
    CEigenfunctions<dim> Ef1(QN, m_omega);
    VectorTools::interpolate(dof_handler, Ef1, m_vPhi);

    auto tpContributions = GP(dynamic_cast<IRealWavefunction<dim>*>(this), m_vPhi, oPotential);

    m_vPhi *= sqrt(1 / std::get<1>(tpContributions));

    std::cout << setprecision(9);
    std::cout << "T = " << std::get<0>(tpContributions) << std::endl;
    std::cout << "N = " << std::get<1>(tpContributions) << std::endl;
    std::cout << "W = " << std::get<2>(tpContributions) << std::endl;

    status = DoIter("");

    if (status == Status::SUCCESS)
    {
      output_results("final_gs_one.gnuplot");
      save("final_gs_one.bin");
    }
  }

  template <int dim>
  void CSobolevGradient<dim>::save(string filename)
  {
    ofstream ofs(filename);
    m_vPhi.block_write(ofs);
  }

  template class CSobolevGradient<1>;
  template class CSobolevGradient<2>;
} // end of namespace
