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

#include <deal.II/lac/sparse_direct.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>

#include <fstream>
#include <iostream>
#include <limits>
#include <stdlib.h>
#include <stdio.h>
#include <nlopt.hpp>

#include "functions.h"
#include "MyParameterHandler.h"
#include "UtilsRealWavefunction.hpp"
#include "orthomin_newton.hpp"
//#include "muParser.h"


namespace BreedSolver_1
{
  using namespace std;
  using namespace dealii;

  template <int dim>
  double myfunc(const std::vector<double>& t, std::vector<double>& grad, void* pSolver)
  {
    const auto& oSolver = *reinterpret_cast<MySolver<dim>*>(pSolver);
    const auto& mC = oSolver.m_coeffs;

    const double t0 = t.at(0);
    const double t1 = t.at(1);

    const double retval = mC.at("t0^2") * t0 * t0
                          + mC.at("t0_t1") * t0 * t1
                          + mC.at("t1^2") * t1 * t1
                          + mC.at("t0^4") * t0 * t0 * t0 * t0
                          + mC.at("t0^1_t1^3") * t0 * t1 * t1 * t1
                          + mC.at("t0^2_t1^2") * t0 * t0 * t1 * t1
                          + mC.at("t0^3_t1^1") * t0 * t0 * t0 * t1
                          + mC.at("t1^4") * t1 * t1 * t1 * t1;

    if (!grad.empty())
    {
      grad.at(0) = 2 * mC.at("t0^2") * t0
                   + mC.at("t0_t1") * t0
                   + 4 * mC.at("t0^4") * t0 * t0 * t0
                   + mC.at("t0^1_t1^3") * t1 * t1 * t1
                   + 2 * mC.at("t0^2_t1^2") * t0 * t1 * t1
                   + 3 * mC.at("t0^3_t1^1") * t0 * t0 * t1;

      grad.at(1) = mC.at("t0_t1") * t1
                   + 2 * mC.at("t1^2") * t1
                   + 3 * mC.at("t0^1_t1^3") * t0 * t1 * t1
                   + 2 * mC.at("t0^2_t1^2") * t0 * t0 * t1
                   + mC.at("t0^3_t1^1") * t0 * t0 * t0
                   + 4 * mC.at("t1^4") * t1 * t1 * t1;
    }
    return retval;
  }


  template <int dim>
  MySolver<dim>::MySolver(const std::string& xml_filename)
    :
    m_oParameters(xml_filename),
    m_oTriangulation(),
    m_oDofHandler(m_oTriangulation)
  {
  }

  template <int dim>
  MySolver<dim>::~MySolver()
  {
    m_oDofHandler.clear();
  }

  template <int dim>
  void MySolver<dim>::setup_system()
  {
    m_oDofHandler.distribute_dofs(m_FE);

    const auto ndofs = m_oDofHandler.n_dofs();

    dealii::DoFTools::extract_locally_relevant_dofs(m_oDofHandler, m_locally_relevant_dofs);

    m_Psi_Ref.reinit(ndofs);
    m_Search_Direction.reinit(ndofs);
    m_System_RHS.reinit(ndofs);
    m_error_per_cell.reinit(m_oTriangulation.n_active_cells());

    for (int i = 0; i < 2; ++i)
    {
      m_Psi[i].reinit(ndofs);
      m_Workspace[i].reinit(ndofs);
    }

    m_constraints.clear();
    m_constraints.reinit(m_locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints(m_oDofHandler, m_constraints);
    VectorTools::interpolate_boundary_values(m_oDofHandler, 0, ZeroFunction<dim>(), m_constraints);
    m_constraints.close();

    DynamicSparsityPattern dsp(ndofs);
    DoFTools::make_sparsity_pattern(m_oDofHandler, dsp);
    SparsityPattern sp;
    sp.copy_from(dsp);
    m_System_Matrix.reinit(sp);
  }

  template <int dim>
  void MySolver<dim>::make_grid()
  {
    Point<dim, double> pt1;
    Point<dim, double> pt2;

    std::vector<double> oGridCornerOne{0, 0, 0};
    std::vector<double> oGridCornerTwo{0, 0, 0};

    m_oParameters.get("grid.grid_corner_one", oGridCornerOne);
    m_oParameters.get("grid.grid_corner_two", oGridCornerTwo);

    for (int i = 0; i < dim; i++)
    {
      pt1(i) = oGridCornerOne[i];
      pt2(i) = oGridCornerTwo[i];
    }

    GridGenerator::hyper_rectangle(m_oTriangulation, pt2, pt1);
    m_oTriangulation.refine_global(m_iGlobalRefinement);
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

    typename DoFHandler<dim>::active_cell_iterator cell = this->m_oDofHandler.begin_active(), endc = this->m_oDofHandler.end();
    for (; cell != endc; ++cell)
    {
      fe_values.reinit(cell);
      fe_values.get_function_values(m_Psi[0], Psi_1);
      fe_values.get_function_values(m_Psi[1], Psi_2);
      fe_values.get_function_gradients(m_Psi[0], Psi_1_grad);
      fe_values.get_function_gradients(m_Psi[1], Psi_2_grad);

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
  void MySolver<dim>::find_ortho_min()
  {
    compute_contributions();

    nlopt::opt opt(nlopt::GN_DIRECT_L, 2);

    static const vector<double> ub{25, 25};
    static const vector<double> lb{0, -25};

    opt.set_upper_bounds(ub);
    opt.set_lower_bounds(lb);

    opt.set_min_objective(myfunc<dim>, this);

    //opt.set_xtol_rel(1e-10);
    opt.set_ftol_rel(1e-10);
    std::vector<double> t(2);
    t[0] = m_t[0];
    t[1] = m_t[1];
    double minf = std::numeric_limits<double>::max();

    try
    {
      nlopt::result result = opt.optimize(t, minf);
      //BOOST_LOG_TRIVIAL(debug) << "found minimum at f(" << t[0] << "," << t[1] << ") = " << std::setprecision(10) << minf;
      // std::cout << "result " << result << std::endl;
    }
    catch (std::exception& e)
    {
      std::cout << "nlopt failed: " << e.what() << std::endl;
    }

    m_t[0] = t[0];
    m_t[1] = t[1];
  }

  template <int dim>
  void MySolver<dim>::assemble_system()
  {
    CPotential<dim> Potential(m_omega);
    const QGauss<dim> quadrature_formula(m_FE.degree + 1);

    m_System_Matrix = 0;

    FEValues<dim> fe_values(m_FE, quadrature_formula, update_values | update_gradients | update_JxW_values | update_quadrature_points);

    const unsigned dofs_per_cell = m_FE.dofs_per_cell;
    const unsigned n_q_points = quadrature_formula.size();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
    vector<Tensor<1, dim>> Psi_ref_grad(n_q_points);
    vector<double> Psi_ref(n_q_points);

    typename DoFHandler<dim>::active_cell_iterator cell = m_oDofHandler.begin_active(), endc = m_oDofHandler.end();
    for (; cell != endc; ++cell)
    {
      cell_matrix = 0;

      fe_values.reinit(cell);
      fe_values.get_function_values(m_Psi_Ref, Psi_ref);
      fe_values.get_function_gradients(m_Psi_Ref, Psi_ref_grad);

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
    /*vector<bool> boundary_dofs (this->m_oDofHandler.n_dofs());
    DoFTools::extract_boundary_dofs (this->m_oDofHandler, ComponentMask(), boundary_dofs);
    for (unsigned int i = 0; i < this->m_oDofHandler.n_dofs(); ++i)
      if (boundary_dofs[i] == true) this->m_Psi_Ref(i) = 0.0;*/

    CPotential<dim> Potential(m_omega);
    const QGauss<dim> quadrature_formula(m_FE.degree + 1);

    m_System_RHS = 0;

    FEValues<dim> fe_values(m_FE, quadrature_formula, update_values | update_gradients | update_JxW_values | update_quadrature_points);

    const unsigned dofs_per_cell = m_FE.dofs_per_cell;
    const unsigned n_q_points = quadrature_formula.size();

    Vector<double> cell_rhs(dofs_per_cell);
    vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
    vector<Tensor<1, dim>> Psi_ref_grad(n_q_points);
    vector<double> Psi_ref(n_q_points);

    typename DoFHandler<dim>::active_cell_iterator cell = m_oDofHandler.begin_active(), endc = m_oDofHandler.end();
    for (; cell != endc; ++cell)
    {
      cell_rhs = 0;
      fe_values.reinit(cell);
      fe_values.get_function_values(m_Psi_Ref, Psi_ref);
      fe_values.get_function_gradients(m_Psi_Ref, Psi_ref_grad);

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
    VectorTools::interpolate_boundary_values(m_oDofHandler, 0, ZeroFunction<dim>(), boundary_values);
    MatrixTools::apply_boundary_values(boundary_values, m_System_Matrix, m_Search_Direction, m_System_RHS);

    SparseDirectUMFPACK A_direct;
    A_direct.solve(m_System_Matrix, m_System_RHS);
    m_Search_Direction = m_System_RHS;
    m_res = m_Search_Direction.l2_norm();
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
    m_Psi_Ref = 0;
    m_Psi_Ref.add(m_t[0], m_Psi[0], m_t[1], m_Psi[1]);
    assemble_rhs();
    //m_res_old = m_res;
    do
    {
      cout << "--------------------------------------------------------------------------------" << endl;
      cout << "- " << m_counter << " - " << path << endl;

      assemble_system();
      solve();

      double tau;
      //MyRealTools::compute_stepsize(this->m_oDofHandler, this->m_FE, Potential, this->m_Psi_Ref, this->m_Search_Direction, m_rMu, m_gs, tau);
      m_Psi[1].add(tau * m_t[1] / fabs(m_t[1]), m_Search_Direction);

      if (m_counter % m_NA == 0)
      {
        output_results(path);
      }

      find_ortho_min();

      m_Psi_Ref = 0;
      m_Psi_Ref.add(m_t[0], m_Psi[0], m_t[1], m_Psi[1]);
      assemble_rhs();

      m_counter++;

      // if (this->l2norm_t() < 1e-4)
      // {
      //   retval = Status::ZERO_SOL;
      // }
      if (m_res < m_epsilon[0])
      {
        retval = Status::SUCCESS;
      }
    }
    while (retval == Status::CONTINUE);

    m_Psi_Ref = 0;
    m_Psi_Ref.add(m_t[0], m_Psi[0], m_t[1], m_Psi[1]);

    // Standard Newton
    // do
    // {
    // }
    // while (true);

    //m_N = MyRealTools::Particle_Number(this->m_oDofHandler, this->m_FE, this->m_Psi_Ref);

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
    VectorTools::interpolate(m_oDofHandler, Ef1, m_Psi[0]);

    //m_Psi[0] *= 1.0 / sqrt(particle_number<dim>(...));
    m_Psi[1] = 0;

    //compute_E_lin(this->m_Psi[0], T, N, W);
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

        m_Psi[0] = m_Psi_Ref;
        m_Psi[1] = 0;
      }
      else
      {
        break;
      }
      m_rMu += m_gs / fabs(m_gs) * m_dmu;
      //compute_E_lin(m_Psi_Ref, T, N, W);
    }
  }


  template <int dim>
  void MySolver<dim>::output_guess()
  {
    string filename = "guess.gnuplot";

    CPotential<dim> Pot(m_omega);
    VectorTools::interpolate(this->m_oDofHandler, Pot, this->m_Workspace[0]);

    DataOut<dim> data_out;
    data_out.attach_dof_handler(this->m_oDofHandler);
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
    data_out.attach_dof_handler(this->m_oDofHandler);
    data_out.add_data_vector(this->m_Psi_Ref, "Psi_ref");
    // data_out.add_data_vector(this->m_Psi[0], "Psi_1");
    // data_out.add_data_vector(this->m_Psi[1], "Psi_2");
    // data_out.add_data_vector(this->m_Search_Direction, "this->m_Search_Direction");
    // data_out.add_data_vector(this->m_error_per_cell, "error_per_cell");
    data_out.build_patches();

    ofstream output(filename.c_str());
    data_out.write_gnuplot(output);
  }

  template class MySolver<1>;
  template class MySolver<2>;
}
