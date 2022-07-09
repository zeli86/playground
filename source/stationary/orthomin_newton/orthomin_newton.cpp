/*
    atus-pro is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    atus-pro is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with atus-pro.  If not, see <http://www.gnu.org/licenses/>.
*/
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

// #include <deal.II/base/logstream.h>
// #include <deal.II/base/utilities.h>

// #include <deal.II/lac/vector.h>
// #include <deal.II/lac/full_matrix.h>
// #include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparse_direct.h>
// #include <deal.II/lac/dynamic_sparsity_pattern.h>

// #include <deal.II/grid/tria.h>
// #include <deal.II/grid/grid_generator.h>
// #include <deal.II/grid/tria_accessor.h>
// #include <deal.II/grid/tria_iterator.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>

#include <fstream>
#include <iostream>
#include <limits>
#include <stdlib.h>
#include <stdio.h>

#include "global.h"
#include "functions.h"
#include "MyParameterHandler.h"
#include "UtilsRealWavefunction.hpp"
#include "muParser.h"
#include "orthomin_newton.hpp"

/// namespace BreedSolver_1
namespace BreedSolver_1
{
  using namespace std;
  using namespace dealii;

  //#include "cBase.h"

  template <int dim>
  MySolver<dim>::MySolver(const std::string& xml_filename)
  {
  }

  template <int dim>
  MySolver<dim>::~MySolver()
  {
  }

  template <int dim>
  void MySolver<dim>::setup_system()
  {
    m_DOF_Handler.distribute_dofs(m_FE);

    const auto ndofs = m_DOF_Handler.n_dofs();

    m_Psi_Ref.reinit(ndofs);
    m_Search_Direction.reinit(ndofs);
    m_System_RHS.reinit(ndofs);
    m_error_per_cell.reinit(m_Triangulation.n_active_cells());

    for (int i = 0; i < 2; ++i)
    {
      m_Psi[i].reinit(ndofs);
      m_Workspace[i].reinit(ndofs);
    }

    // m_constraints.clear();
    // m_constraints.reinit(ndofs);
    // DoFTools::make_hanging_node_constraints(m_DOF_Handler, m_constraints);
    // VectorTools::interpolate_boundary_values(m_DOF_Handler, 0, ZeroFunction<dim>(), m_constraints);
    // m_constraints.close();

    DynamicSparsityPattern dsp(ndofs);
    DoFTools::make_sparsity_pattern(m_DOF_Handler, dsp);
    SparsityPattern sp;
    sp.copy_from(dsp);
    m_System_Matrix.reinit(sp);
  }

  template <int dim>
  void MySolver<dim>::make_grid()
  {
  }

  template <int dim>
  void MySolver<dim>::compute_contributions()
  {
    CPotential<dim> Potential(m_omega);
    const QGauss<dim>  quadrature_formula(this->m_FE.degree + 1);
    FEValues<dim> fe_values(this->m_FE, quadrature_formula, update_gradients | update_values | update_JxW_values | update_quadrature_points);

    const unsigned n_q_points = quadrature_formula.size();

    vector<double> Psi_1(n_q_points);
    vector<double> Psi_2(n_q_points);
    vector<Tensor<1, dim>> Psi_1_grad(n_q_points);
    vector<Tensor<1, dim>> Psi_2_grad(n_q_points);

    std::array<double, 8> total_contributions;
    total_contributions.fill(0);

    typename DoFHandler<dim>::active_cell_iterator cell = this->m_DOF_Handler.begin_active(), endc = this->m_DOF_Handler.end();
    for (; cell != endc; ++cell)
    {
      fe_values.reinit(cell);
      fe_values.get_function_values(this->m_Psi[0], Psi_1);
      fe_values.get_function_values(this->m_Psi[1], Psi_2);
      fe_values.get_function_gradients(this->m_Psi[0], Psi_1_grad);
      fe_values.get_function_gradients(this->m_Psi[1], Psi_2_grad);

      for (unsigned qp = 0; qp < n_q_points; ++qp)
      {
        const double JxW = fe_values.JxW(qp);
        const double p12 = Psi_1[qp] * Psi_2[qp];
        const double p1q = Psi_1[qp] * Psi_1[qp];
        const double p2q = Psi_2[qp] * Psi_2[qp];
        const double Q = Potential.value(fe_values.quadrature_point(qp)) - m_rMu;

        total_contributions[0] += JxW * p1q * p1q;
        total_contributions[1] += JxW * p2q * p2q;
        total_contributions[2] += JxW * p1q * p2q;
        total_contributions[3] += JxW * p1q * p12;
        total_contributions[4] += JxW * p2q * p12;
        total_contributions[5] += JxW * (Psi_1_grad[qp] * Psi_1_grad[qp] + Q * p1q);
        total_contributions[6] += JxW * (Psi_2_grad[qp] * Psi_2_grad[qp] + Q * p2q);
        total_contributions[7] += JxW * (Psi_1_grad[qp] * Psi_2_grad[qp] + Q * p12);
      }
    }

    m_coeffs["t0^2"]      = 0.5 * total_contributions[5];
    m_coeffs["t0_t1"]     = total_contributions[7];
    m_coeffs["t1^2"]      = 0.5 * total_contributions[6];
    m_coeffs["t0^4"]      = 0.25 * m_gs * total_contributions[0];
    m_coeffs["t0^1_t1^3"] = 0.25 * m_gs * 4 * total_contributions[4];
    m_coeffs["t0^2_t1^2"] = 0.25 * m_gs * 6 * total_contributions[2];
    m_coeffs["t0^3_t1^1"] = 0.25 * m_gs * 4 * total_contributions[3];
    m_coeffs["t1^4"]      = 0.25 * m_gs * total_contributions[1];
  }


  template <int dim>
  void MySolver<dim>::assemble_system()
  {
    CPotential<dim> Potential(m_omega);
    const QGauss<dim> quadrature_formula(this->m_FE.degree + 1);

    this->m_System_Matrix = 0;

    FEValues<dim> fe_values(this->m_FE, quadrature_formula, update_values | update_gradients | update_JxW_values | update_quadrature_points);

    const unsigned dofs_per_cell = this->m_FE.dofs_per_cell;
    const unsigned n_q_points = quadrature_formula.size();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
    vector<Tensor<1, dim>> Psi_ref_grad(n_q_points);
    vector<double> Psi_ref(n_q_points);

    typename DoFHandler<dim>::active_cell_iterator cell = this->m_DOF_Handler.begin_active(), endc = this->m_DOF_Handler.end();
    for (; cell != endc; ++cell)
    {
      cell_matrix = 0;

      fe_values.reinit(cell);
      fe_values.get_function_values(this->m_Psi_Ref, Psi_ref);
      fe_values.get_function_gradients(this->m_Psi_Ref, Psi_ref_grad);

      for (unsigned qp = 0; qp < n_q_points; ++qp)
      {
        const double JxW = fe_values.JxW(qp);
        const double Q2 = Potential.value(fe_values.quadrature_point(qp)) - m_rMu + 3.0 * m_gs * Psi_ref[qp] * Psi_ref[qp];

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
          this->m_System_Matrix.add(local_dof_indices[i], local_dof_indices[j], cell_matrix(i, j));
        }
    }
  }


  template <int dim>
  void MySolver<dim>::assemble_rhs()
  {
    /*vector<bool> boundary_dofs (this->m_DOF_Handler.n_dofs());
    DoFTools::extract_boundary_dofs (this->m_DOF_Handler, ComponentMask(), boundary_dofs);
    for (unsigned int i = 0; i < this->m_DOF_Handler.n_dofs(); ++i)
      if (boundary_dofs[i] == true) this->m_Psi_Ref(i) = 0.0;*/

    CPotential<dim> Potential(m_omega);
    const QGauss<dim> quadrature_formula(this->m_FE.degree + 1);

    this->m_System_RHS = 0;

    FEValues<dim> fe_values(this->m_FE, quadrature_formula, update_values | update_gradients | update_JxW_values | update_quadrature_points);

    const unsigned dofs_per_cell = this->m_FE.dofs_per_cell;
    const unsigned n_q_points = quadrature_formula.size();

    Vector<double> cell_rhs(dofs_per_cell);
    vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
    vector<Tensor<1, dim>> Psi_ref_grad(n_q_points);
    vector<double> Psi_ref(n_q_points);

    typename DoFHandler<dim>::active_cell_iterator cell = this->m_DOF_Handler.begin_active(), endc = this->m_DOF_Handler.end();
    for (; cell != endc; ++cell)
    {
      cell_rhs = 0;
      fe_values.reinit(cell);
      fe_values.get_function_values(this->m_Psi_Ref, Psi_ref);
      fe_values.get_function_gradients(this->m_Psi_Ref, Psi_ref_grad);

      for (unsigned qp = 0; qp < n_q_points; ++qp)
      {
        const double JxW = fe_values.JxW(qp);
        const double Q1 = Potential.value(fe_values.quadrature_point(qp)) - m_rMu + m_gs * Psi_ref[qp] * Psi_ref[qp];

        for (unsigned i = 0; i < dofs_per_cell; ++i)
        {
          cell_rhs(i) += JxW * (Psi_ref_grad[qp] * fe_values.shape_grad(i, qp) + Q1 * Psi_ref[qp] * fe_values.shape_value(i, qp));
        }
      }
      cell->get_dof_indices(local_dof_indices);
      for (unsigned i = 0; i < dofs_per_cell; ++i)
      {
        m_System_RHS(local_dof_indices[i]) += cell_rhs(i);
      }
    }
  }


  template <int dim>
  void MySolver<dim>::solve()
  {
    map<types::global_dof_index, double> boundary_values;
    VectorTools::interpolate_boundary_values(this->m_DOF_Handler, 0, ZeroFunction<dim>(), boundary_values);
    MatrixTools::apply_boundary_values(boundary_values, this->m_System_Matrix, this->m_Search_Direction, this->m_System_RHS);

    SparseDirectUMFPACK A_direct;
    A_direct.solve(this->m_System_Matrix, this->m_System_RHS);
    this->m_Search_Direction = this->m_System_RHS;
    m_res = this->m_Search_Direction.l2_norm();
  }


  template <int dim>
  int MySolver<dim>::DoIter(string path)
  {
    CPotential<dim> Potential(m_omega);

    int retval = Status::CONTINUE;

    m_t[0] = m_ti;
    m_t[1] = m_ti;

    m_res = std::numeric_limits<double>::max();
    m_counter = 0;
    this->m_Psi_Ref = 0;
    this->m_Psi_Ref.add(m_t[0], this->m_Psi[0], m_t[1], this->m_Psi[1]);
    assemble_rhs();
    //m_res_old = m_res;
    do
    {
      cout << "--------------------------------------------------------------------------------" << endl;
      cout << "- " << m_counter << " - " << path << endl;

      assemble_system();
      solve();

      double tau;
      //MyRealTools::compute_stepsize(this->m_DOF_Handler, this->m_FE, Potential, this->m_Psi_Ref, this->m_Search_Direction, m_rMu, m_gs, tau);
      this->m_Psi[1].add(tau * m_t[1] / fabs(m_t[1]), this->m_Search_Direction);

      if (m_counter % m_NA == 0)
      {
        output_results(path);
      }

      this->find_ortho_min();

      this->m_Psi_Ref = 0;
      this->m_Psi_Ref.add(m_t[0], this->m_Psi[0], m_t[1], this->m_Psi[1]);
      assemble_rhs();

      m_counter++;

      if (this->l2norm_t() < 1e-4)
      {
        retval = Status::ZERO_SOL;
      }
      if (m_res < m_epsilon[0])
      {
        retval = Status::SUCCESS;
      }
    }
    while (retval == Status::CONTINUE);

    this->m_Psi_Ref = 0;
    this->m_Psi_Ref.add(m_t[0], this->m_Psi[0], m_t[1], this->m_Psi[1]);

    // Standard Newton
    do
    {
      cout << "--------------------------------------------------------------------------------" << endl;
      cout << "-- " << path << " - " << m_counter << endl;

      assemble_system();
      solve();

      this->m_Psi_Ref.add(-1, this->m_Search_Direction);
      assemble_rhs();

      if (m_counter % m_NA == 0)
      {
        output_results(path);
      }

      if (m_res < m_epsilon[1])
      {
        retval = Status::SUCCESS;
        break;
      }

      m_counter++;
    }
    while (true);

    //m_N = MyRealTools::Particle_Number(this->m_DOF_Handler, this->m_FE, this->m_Psi_Ref);

    if (m_N < 1e-5)
    {
      retval = Status::ZERO_SOL;
    }
    return retval;
  }


  template <int dim>
  void MySolver<dim>::run()
  {
    string path;
    char shellcmd[255];
    double T, N, W;

    make_grid();
    this->setup_system();

    CEigenfunctions<dim> Ef1(m_QN1, m_omega);
    VectorTools::interpolate(m_DOF_Handler, Ef1, m_Psi[0]);

    //m_Psi[0] *= 1.0 / sqrt(particle_number<dim>(...));
    m_Psi[1] = 0;

    compute_E_lin(this->m_Psi[0], T, N, W);
    double m_rMu_0 = T / N;
    m_rMu = ceil(10.0 * m_rMu_0) / 10.0 + m_gs / fabs(m_gs) * m_dmu;

    cout << "T = " << T << endl;
    cout << "N = " << N << endl;
    cout << "W = " << W << endl;
    cout << "m_rMu = " << m_rMu << endl;
    cout << "m_gs = " << m_gs << endl;

    output_guess();
    for (unsigned i = 0; i < m_Ndmu; ++i)
    {
      sprintf(shellcmd, "mkdir %.4d/", i);
      system(shellcmd);
      sprintf(shellcmd, "%.4d/", i);
      path = shellcmd;

      // Nehari Reset
      // sqrt((m_rMu*N-T)/(m_gs*W)); if this->m_Psi[1] == 0
      // sqrt((m_rMu*N-T)/(4.0*m_gs*W)); if this->m_Psi[1] == this->m_Psi[0]
#ifdef NEHARI
      m_ti = sqrt((m_rMu * N - T) / (m_gs * W));
#endif

      const auto status = DoIter(path);

      if (status == Status::SUCCESS)
      {
        // output_results(path, "final");
        // this->save(path + "final.bin");
        // this->save_one(path + "final_one.bin");

        // vector<double> newgs = {m_gs * m_N};
        // m_ph.Set_Physics("gs_1", newgs);
        // m_ph.SaveXMLFile(path + "params_one.xml");
        // newgs[0] = m_gs;
        // m_ph.Set_Physics("gs_1", newgs);
        // m_ph.SaveXMLFile(path + "params.xml");

        this->m_Psi[0] = this->m_Psi_Ref;
        this->m_Psi[1] = 0;
      }
      else
      {
        break;
      }
      m_rMu += m_gs / fabs(m_gs) * m_dmu;
      compute_E_lin(this->m_Psi_Ref, T, N, W);
    }
  }


  template <int dim>
  void MySolver<dim>::output_guess()
  {
    string filename = "guess.gnuplot";

    CPotential<dim> Pot(m_omega);
    VectorTools::interpolate(this->m_DOF_Handler, Pot, this->m_Workspace[0]);

    DataOut<dim> data_out;
    data_out.attach_dof_handler(this->m_DOF_Handler);
    data_out.add_data_vector(this->m_Psi[0], "Psi_0");
    data_out.add_data_vector(this->m_Psi[1], "Psi_1");
    data_out.add_data_vector(this->m_Workspace[0], "Potential");
    data_out.build_patches();

    ofstream output(filename.c_str());
    data_out.write_gnuplot(output);
  }


  template <int dim>
  void MySolver<dim>::output_results(string path, string prefix)
  {
    string filename = path + prefix + "-" + Utilities::int_to_string(m_counter, 5) + ".gnuplot";

    DataOut<dim> data_out;
    data_out.attach_dof_handler(this->m_DOF_Handler);
    data_out.add_data_vector(this->m_Psi_Ref, "Psi_ref");
    data_out.add_data_vector(this->m_Psi[0], "Psi_1");
    data_out.add_data_vector(this->m_Psi[1], "Psi_2");
    data_out.add_data_vector(this->m_Search_Direction, "this->m_Search_Direction");
    data_out.add_data_vector(this->m_error_per_cell, "error_per_cell");
    data_out.build_patches();

    ofstream output(filename.c_str());
    data_out.write_gnuplot(output);
  }
} // end of namespace

