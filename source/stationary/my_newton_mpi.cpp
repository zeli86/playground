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

#include "default_includes.h"

#include <az_aztec.h>
#include <boost/log/trivial.hpp>
#include <deal.II/base/mpi.h>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <locale>
#include <limits>
#include <cmath>
#include <array>
#include <algorithm>
#include <vector>

#include "global.h"
#include "mpi.h"
#include "functions.h"
//#include "ref_pt_list.h"
#include "MyParameterHandler.h"
#include "UtilsRealWavefunction.hpp"
#include "GPUtilsRealWavefunction.hpp"
#include "muParser.h"

namespace BreedSolver
{
  using namespace std;
  using namespace dealii;

  template <int dim>
  class MySolver : public utils::real_wavefunction::IRealWavefunction<dim>
  {
  public:
    explicit MySolver(const std::string);

    void run2b();

    double m_I[8];
  protected:

    MPI_Comm mpi_communicator;

    int DoIter(string = "");

    void make_grid();
    void setup_system();

    //void estimate_error(double&);

    bool solve();
    void compute_contributions();

    void output_vector(LA::MPI::Vector&, string);

    double m_t[2] = {0, 0};
    double m_ti{1};
    double m_N;
    std::vector<double> m_omega{0.25, 0.25, 0.25};
    unsigned m_counter;
    //dealii::TimerOutput m_computing_timer;
    double m_final_error;
    unsigned m_NA{100};
    double m_Ndmu;
    double m_dmu;
    unsigned m_QN1[3] = {1, 1, 1};
    double m_res;
    double m_resp;
    double m_res_old;
    double m_epsilon{1e-10};
    int m_iMaxIter{100};
    int m_iTotalNoCells{0};
    int m_iTotalNoActiveCells{0};
    int m_iGlobalRefinement{8};

    std::map<std::string, double> m_coeffs;

  private:
    dealii::DoFHandler<dim>& get_dof_handler()
    {
      return m_oDofHandler;
    }

    dealii::FE_Q<dim>& get_fe()
    {
      return m_oFe;
    }

    dealii::AffineConstraints<double>& get_constraints()
    {
      return m_constraints;
    }

    MyParameterHandler m_oParameters;
    bool m_root{false};
    tMpiTriangulation<dim> m_oTriangulation;
    tDoFHandler<dim> m_oDofHandler;
    tConstraints m_oConstraints;
    dealii::FE_Q<dim> m_oFe{2};

    tIndexSet m_locally_owned_dofs;
    tIndexSet m_locally_relevant_dofs;

    tMpiVector m_vPsiRef;
    double m_rG{0};
    double m_rMu{0};

    tMpiVector m_Psi_Ref; // non ghosted
    tMpiVector m_System_RHS; // non ghosted
    tMpiVector m_Search_Direction; // non ghosted
    tMpiVector m_Workspace_NG;
    dealii::Vector<double> m_error_per_cell;

    std::array<tMpiVector, 2> m_Psi; // non ghosted
    std::array<tMpiVector, 2> m_Workspace; // ghosted

    tMpiMatrix m_System_Matrix;

    dealii::AffineConstraints<double> m_constraints;
  };

  /**
   * Constructor
   */
  template <int dim>
  MySolver<dim>::MySolver(const std::string sConfigFile)
    :
    m_oParameters(sConfigFile),
    m_root(dealii::Utilities::MPI::this_mpi_process(mpi_communicator) == 0),
    m_oTriangulation(mpi_communicator, typename Triangulation<dim>::MeshSmoothing(Triangulation<dim>::limit_level_difference_at_vertices | Triangulation<dim>::eliminate_refined_inner_islands | Triangulation<dim>::smoothing_on_refinement | Triangulation<dim>::smoothing_on_coarsening)),
    m_oDofHandler(m_oTriangulation)
  {
  }

  template <int dim>
  void MySolver<dim>::make_grid()
  {
    // TimerOutput::Scope timing_section(m_computing_timer, "");

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

    GridGenerator::hyper_rectangle(m_oTriangulation, pt2, pt1);
    m_oTriangulation.refine_global(m_iGlobalRefinement);

    //m_iTotalNoCells = Utilities::MPI::sum(tmp[0], mpi_communicator);
    //m_iTotalNoActiveCells = Utilities::MPI::sum(tmp[1], mpi_communicator);
  }

  template <int dim>
  void MySolver<dim>::setup_system()
  {
    //dealii::TimerOutput::Scope timing_section(m_computing_timer, "");

    m_oDofHandler.distribute_dofs(m_oFe);

    m_locally_owned_dofs = m_oDofHandler.locally_owned_dofs();

    dealii::DoFTools::extract_locally_relevant_dofs(m_oDofHandler, m_locally_relevant_dofs);

    m_Psi_Ref.reinit(m_locally_owned_dofs, m_locally_relevant_dofs, mpi_communicator);
    m_Search_Direction.reinit(m_locally_owned_dofs, mpi_communicator);
    m_System_RHS.reinit(m_locally_owned_dofs, mpi_communicator);
    m_Workspace_NG.reinit(m_locally_owned_dofs, mpi_communicator);
    m_error_per_cell.reinit(m_oTriangulation.n_active_cells());

    m_Psi[0].reinit(m_locally_owned_dofs, mpi_communicator);
    m_Workspace[0].reinit(m_locally_owned_dofs, m_locally_relevant_dofs, mpi_communicator);
    m_Psi[1].reinit(m_locally_owned_dofs, mpi_communicator);
    m_Workspace[1].reinit(m_locally_owned_dofs, m_locally_relevant_dofs, mpi_communicator);

    m_oConstraints.clear();
    m_oConstraints.reinit(m_locally_relevant_dofs);
    dealii::DoFTools::make_hanging_node_constraints(m_oDofHandler, m_oConstraints);
    dealii::VectorTools::interpolate_boundary_values(m_oDofHandler, 0, ZeroFunction<dim>(), m_oConstraints);
    m_oConstraints.close();

    dealii::DynamicSparsityPattern csp(m_locally_relevant_dofs);
    dealii::DoFTools::make_sparsity_pattern(m_oDofHandler, csp, m_oConstraints, false);
    dealii::SparsityTools::distribute_sparsity_pattern(csp, m_oDofHandler.locally_owned_dofs(), mpi_communicator, m_locally_relevant_dofs);
    m_System_Matrix.reinit(m_locally_owned_dofs, m_locally_owned_dofs, csp, mpi_communicator);
  }

  template <int dim>
  void MySolver<dim>::compute_contributions()
  {
    using namespace utils::real_wavefunction;
    //TimerOutput::Scope timing_section(m_computing_timer, "");

    this->update_workspace();

    std::cout << "N0 = " << particle_number(dynamic_cast<IRealWavefunction<dim>*>(this), this->m_Workspace[0], mpi_communicator) << "\n";
    std::cout << "N1 = " << particle_number(dynamic_cast<IRealWavefunction<dim>*>(this), this->m_Workspace[1], mpi_communicator) << "\n";

    CPotential<dim> Potential_fct(m_omega);
    const QGauss<dim> quadrature_formula(m_oFe.degree + 1);
    FEValues<dim> fe_values(m_oFe, quadrature_formula, update_gradients | update_values | update_JxW_values | update_quadrature_points);

    const unsigned n_q_points = quadrature_formula.size();

    vector<double> Psi_0(n_q_points);
    vector<double> Psi_1(n_q_points);
    vector<Tensor<1, dim>> Psi_0_grad(n_q_points);
    vector<Tensor<1, dim>> Psi_1_grad(n_q_points);

    std::array<double, 8> local_contributions;
    std::array<double, 8> total_contributions;
    local_contributions.fill(0);
    total_contributions.fill(0);

    typename DoFHandler<dim>::active_cell_iterator cell = m_oDofHandler.begin_active(), endc = m_oDofHandler.end();
    for (; cell != endc; ++cell)
    {
      if (cell->is_locally_owned())
      {
        fe_values.reinit(cell);
        fe_values.get_function_values(this->m_Workspace[0], Psi_0);
        fe_values.get_function_values(this->m_Workspace[1], Psi_1);
        fe_values.get_function_gradients(this->m_Workspace[0], Psi_0_grad);
        fe_values.get_function_gradients(this->m_Workspace[1], Psi_1_grad);

        for (unsigned qp = 0; qp < n_q_points; ++qp)
        {
          const double JxW = fe_values.JxW(qp);
          const double p01 = Psi_0[qp] * Psi_1[qp];
          const double p0q = Psi_0[qp] * Psi_0[qp];
          const double p1q = Psi_1[qp] * Psi_1[qp];
          const double Q = Potential_fct.value(fe_values.quadrature_point(qp)) - m_rMu;

          local_contributions[0] += JxW * p0q * p0q;
          local_contributions[1] += JxW * p1q * p1q;
          local_contributions[2] += JxW * p0q * p1q;
          local_contributions[3] += JxW * p0q * p01;
          local_contributions[4] += JxW * p1q * p01;
          local_contributions[5] += JxW * (Psi_0_grad[qp] * Psi_0_grad[qp] + Q * p0q);
          local_contributions[6] += JxW * (Psi_1_grad[qp] * Psi_1_grad[qp] + Q * p1q);
          local_contributions[7] += JxW * (Psi_0_grad[qp] * Psi_1_grad[qp] + Q * p01);
        }
      }
    }

    MPI_Allreduce(local_contributions.data(), total_contributions.data(), total_contributions.size(), MPI_DOUBLE, MPI_SUM, mpi_communicator);

    this->m_coeffs["t0^2"]      = 0.5 * total_contributions[5];
    this->m_coeffs["t0_t1"]     = total_contributions[7];
    this->m_coeffs["t1^2"]      = 0.5 * total_contributions[6];
    this->m_coeffs["t0^4"]      = 0.25 * m_rG * total_contributions[0];
    this->m_coeffs["t0^1_t1^3"] = 0.25 * m_rG * 4 * total_contributions[4];
    this->m_coeffs["t0^2_t1^2"] = 0.25 * m_rG * 6 * total_contributions[2];
    this->m_coeffs["t0^3_t1^1"] = 0.25 * m_rG * 4 * total_contributions[3];
    this->m_coeffs["t1^4"]      = 0.25 * m_rG * total_contributions[1];
  }


  template <int dim>
  void MySolver<dim>::output_vector(LA::MPI::Vector& vec, string filename)
  {
    // TimerOutput::Scope timing_section(m_computing_timer, "");

    // this->m_constraints.distribute(vec);
    // this->m_Workspace[0] = vec;

    // DataOut<dim> data_out;
    // data_out.attach_m_oDofHandler(m_oDofHandler);
    // data_out.add_data_vector(this->m_Workspace[0], "vec");
    // data_out.build_patches(gl_subdivisions);
    // data_out.write_vtu_in_parallel(filename.c_str(), mpi_communicator);
  }

  template <int dim>
  int MySolver<dim>::DoIter(string path)
  {
    using namespace utils::real_wavefunction;

    int retval = Status::SUCCESS;

    m_t[0] = m_ti;
    m_t[1] = m_ti;

    CPotential<dim> Potential(m_omega);

    this->do_linear_superposition();
    this->m_Workspace[0] = this->m_Psi_Ref;
    assemble_L2gradient<dim>(dynamic_cast<IRealWavefunction<dim>*>(this), this->m_Workspace[0], Potential, m_rMu, m_rG, this->m_System_RHS);
    m_res_old = m_res;
    m_counter = 0;
    bool bempty;

    do
    {
      BOOST_LOG_TRIVIAL(info) << std::string('-', 80);
      BOOST_LOG_TRIVIAL(info) << "- " << path << " - " << m_counter;

      this->m_Workspace[0] = this->m_Psi_Ref;
      assemble_jacobian<dim>(dynamic_cast<IRealWavefunction<dim>*>(this), this->m_Workspace[0], Potential, m_rMu, m_rG, this->m_System_Matrix);
      bool bsucc = solve();
      if (!bsucc)
      {
        return Status::SINGULAR;
      }

      /*
            m_Workspace[0] = m_Psi_1;
            m_Workspace[1] = m_Search_Direction;
            orthonormalize(...);
            m_Search_Direction = m_Workspace_NG;
      */
      double tau = -0.5;
      this->m_Workspace[1] = this->m_Search_Direction;
      //compute_stepsize(...);
      std::cout << "NL2 = " << particle_number(dynamic_cast<IRealWavefunction<dim>*>(this), this->m_Workspace[1], mpi_communicator) << "\n";

      m_Psi[1].add(tau * m_t[1] / fabs(m_t[1]), m_Search_Direction);
      m_constraints.distribute(m_Psi[1]);

      bempty = this->find_ortho_min();

      if (bempty)
      {
        return Status::FAILED;
      }

      this->do_linear_superposition();
      this->m_Workspace[0] = this->m_Psi_Ref;
      assemble_L2gradient<dim>(dynamic_cast<IRealWavefunction<dim>*>(this), this->m_Workspace[0], Potential, m_rMu, m_rG, this->m_System_RHS);

      m_resp = m_res_old - m_res;
      m_res_old = m_res;

      if (m_counter % m_NA == 0)
      {
        this->output_results(path);
      }

      // columns& cols = m_table.new_line();
      // m_table.insert(cols, MyTable::COUNTER, double(m_counter));
      // m_table.insert(cols, MyTable::RES, m_res);
      // m_table.insert(cols, MyTable::RESP, m_resp);
      // m_table.insert(cols, MyTable::STEPSIZE, tau);
      // m_table.insert(cols, MyTable::MU, m_rMu);
      // m_table.insert(cols, MyTable::GS, m_gs);
      // m_table.insert(cols, MyTable::t1, m_t[0]);
      // m_table.insert(cols, MyTable::t2, m_t[1]);
      // m_table.insert(cols, MyTable::l2norm_t, this->l2norm_t());

      m_counter++;
      // if (m_root)
      // {
      //   cout << m_table;
      // }
      // if (m_res < m_epsilon[0])
      // {
      //   retval = Status::SUCCESS;
      //   break;
      // }
      // if (this->l2norm_t() < 1e-4)
      // {
      //   retval = Status::ZERO_SOL;
      //   break;
      // }
      // if (m_counter == m_maxiter)
      // {
      //   retval = Status::MAXITER;
      //   break;
      // }
      // if (isnan(m_res))
      // {
      //   retval = Status::FAILED;
      //   break;
      // }
    }
    while (true);

    // particle_number(m_Psi_Ref);

    // Standard Newton
    // do
    // {
    //   BOOST_LOG_TRIVIAL(info) << std::string('-', 80);
    //   BOOST_LOG_TRIVIAL(info) << "-- " << path << " - " << m_counter;

    //   m_Workspace[0] = m_Psi_Ref;
    //   assemble_jacobian<dim>(...);
    //   bool bsucc = solve();
    //   if (!bsucc)
    //   {
    //     return Status::SINGULAR;
    //   }

    //   double tau;
    //   m_Workspace[1] = m_Search_Direction;
    //   compute_stepsize(...);

    //   m_Psi_Ref.add(tau * m_t[1] / fabs(m_t[1]), m_Search_Direction);
    //   m_constraints.distribute(m_Psi_Ref);

    //   m_Workspace[0] = m_Psi_Ref;
    //   assembleRHS_L2gradient<dim>(...);

    //   m_resp = m_res_old - m_res;
    //   m_res_old = m_res;

    // if (m_counter % m_NA == 0)
    // {
    //   this->output_results(path);
    // }

    // columns& cols = m_table.new_line();
    // m_table.insert(cols, MyTable::COUNTER, double(m_counter));
    // m_table.insert(cols, MyTable::RES, m_res);
    // m_table.insert(cols, MyTable::RESP, m_resp);
    // m_table.insert(cols, MyTable::STEPSIZE, tau);
    // m_table.insert(cols, MyTable::MU, m_rMu);
    // m_table.insert(cols, MyTable::GS, m_rG);
    // m_table.insert(cols, MyTable::t1, m_t[0]);
    // m_table.insert(cols, MyTable::t2, m_t[1]);
    // m_table.insert(cols, MyTable::l2norm_t, this->l2norm_t());
    // m_table.insert(cols, MyTable::PARTICLE_NUMBER, m_N);

    // m_counter++;

    // if (m_root)
    // {
    //   cout << m_table;
    // }
    // if (m_res < m_epsilon[1])
    // {
    //   retval = Status::SUCCESS;
    //   break;
    // }
    // if (m_resp < 0)
    // {
    //   retval = Status::SUCCESS;
    //   break;
    // }
    // if (m_counter == m_maxiter)
    // {
    //   retval = Status::MAXITER;
    //   break;
    // }
    // }
    // while (true);

    m_Workspace[0] = m_Psi_Ref;
    m_N = particle_number(dynamic_cast<IRealWavefunction<dim>*>(this), this->m_Workspace[0], mpi_communicator);

    if (m_N < 1e-5)
    {
      retval = Status::ZERO_SOL;
    }
    return retval;
  }


  template <int dim>
  void MySolver<dim>::run2b()
  {
    using namespace utils::real_wavefunction;

    string path;
    char shellcmd[255];
    double T, N, W;
    int status;

    make_grid();
    setup_system();

    CEigenfunctions<dim> Ef1(m_QN1, m_omega);
    //CEigenfunctions<dim> Ef2( m_QN2, m_omega );
    CPotential<dim> Potential(m_omega);

    VectorTools::interpolate(m_oDofHandler, Ef1, this->m_Psi[0]);

    this->m_Workspace[0] = this->m_Psi[0];
    this->m_Psi[0] *= 1.0 / sqrt(particle_number(dynamic_cast<IRealWavefunction<dim>*>(this), this->m_Workspace[0], mpi_communicator));
    this->m_Psi[1] = 0;

    this->m_Workspace[0] = this->m_Psi[0];
    auto tp = GP(dynamic_cast<IRealWavefunction<dim>*>(this), Potential, this->m_Workspace[0], mpi_communicator);

    double m_rMu_0 = T / N;
    m_rMu = ceil(10.0 * m_rMu_0) / 10.0 + m_rG / fabs(m_rG) * m_dmu;

    //output_guess();
    for (int i = 0; i < m_Ndmu; ++i)
    {
      sprintf(shellcmd, "mkdir %.4d/", i);
      if (m_root)
      {
        system(shellcmd);
      }
      sprintf(shellcmd, "%.4d/", i);
      path = shellcmd;

      // nehari
      // sqrt((m_rMu*N-T)/(m_gs*W)); if this->m_Psi[1] == 0
      // sqrt((m_rMu*N-T)/(4.0*m_gs*W)); if this->m_Psi[1] == m_Psi_1
      // m_ti = sqrt((m_rMu * N - T) / (m_rG * W));

      // BOOST_LOG_TRIVIAL(info) << "T = " << T;
      // BOOST_LOG_TRIVIAL(info) << "N = " << N;
      // BOOST_LOG_TRIVIAL(info) << "W = " << W;
      // BOOST_LOG_TRIVIAL(info) << "m_rMu = " << m_rMu;
      // BOOST_LOG_TRIVIAL(info) << "m_ti = " << m_ti;

      status = DoIter(path);

      // columns& cols = m_results.new_line();
      // m_results.insert(cols, MyTable::MU, m_rMu);
      // m_results.insert(cols, MyTable::GS, m_gs);
      // m_results.insert(cols, MyTable::PARTICLE_NUMBER, m_N);
      // m_results.insert(cols, MyTable::COUNTER, double(m_counter));
      // m_results.insert(cols, MyTable::STATUS, double(status));

      // if (status == Status::SUCCESS)
      // {
      //   estimate_error(m_final_error);

      //   this->save(path + "final.bin");

      //   m_ph.SaveXMLFile(path + "params.xml");

      //   this->output_results(path, "final");
      //   this->dump_info_xml(path);
      //   this->m_Psi[0]  = this->m_Psi_Ref;
      //   this->m_Psi[1] = 0;
      // }
      // else if (status == Status::MAXITER || status == Status::SINGULAR)
      // {
      //   this->m_Psi_Ref = this->m_Psi[0];
      //   this->m_Psi[1] = 0;
      // }
      // else
      // {
      //   //this->m_Psi[1] = 0;
      //   break;
      // }
      // m_rMu += m_rG / fabs(m_rG) * m_dmu;
      // m_Workspace[0] = m_Psi[0];
      // Compute_E_lin(mpi_communicator, m_oDofHandler, m_oFe, this->m_constraints, Potential, this->m_Workspace[0], T, N, W);

      // BOOST_LOG_TRIVIAL(info) << "T = " << T << endl;
      // BOOST_LOG_TRIVIAL(info) << "N = " << N << endl;
      // BOOST_LOG_TRIVIAL(info) << "W = " << W << endl;
    }
  }
} // end of namespace

int main(int argc, char* argv[])
{
  using namespace dealii;
  deallog.depth_console(0);

  MyParameterHandler params("params.xml");
  int dim = 0;

  try
  {
    params.get("parameter.spatial_dimension", dim);
  }
  catch (mu::Parser::exception_type& e)
  {
    cout << "Message:  " << e.GetMsg() << "\n";
    cout << "Formula:  " << e.GetExpr() << "\n";
    cout << "Token:    " << e.GetToken() << "\n";
    cout << "Position: " << e.GetPos() << "\n";
    cout << "Errc:     " << e.GetCode() << "\n";
  }

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);
  {
    switch (dim)
    {
    case 2:
    {
      BreedSolver::MySolver<2> solver("params.xml");
      solver.run2b();
      break;
    }
    case 3:
    {
      BreedSolver::MySolver<3> solver("params.xml");
      solver.run2b();
      break;
    }
    default:
      cout << "You have found a new dimension!" << endl;
    }
  }
  return EXIT_SUCCESS;
}
