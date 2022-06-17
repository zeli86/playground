//
// atus-pro testing - atus-pro testing playgroung
// Copyright (C) 2020 Želimir Marojević <zelimir.marojevic@gmail.com>
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

#include <boost/log/trivial.hpp>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <locale>
#include <limits>
#include <cmath>
#include <array>
#include <algorithm>

#include "global.h"
#include "mpi.h"
#include "functions.h"
//#include "ref_pt_list.h"
#include "MyParameterHandler.h"
#include "MyRealTools.h"
#include "muParser.h"

namespace BreedSolver
{
  using namespace std;
  using namespace dealii;

#include "cBaseMPI.h"

  template <int dim>
  class MySolver : public cBaseMPI<dim, 2, dealii::FE_Q<dim>, false>
  {
  public:
    explicit MySolver(const std::string);

    void run2b();

    double m_I[8];
  protected:
    int DoIter(string = "");

    void make_grid();
    void make_grid_custom();

    void estimate_error(double&);

    bool solve();
    void compute_contributions();

    void output_vector(LA::MPI::Vector&, string);

    dealii::FE_Q<dim> m_FEQ{gl_degree_fe};

    double m_mu_punkt;

    using cBaseMPI<dim, 2, dealii::FE_Q<dim>, false>::mpi_communicator;
    using cBaseMPI<dim, 2, dealii::FE_Q<dim>, false>::m_root;
    using cBaseMPI<dim, 2, dealii::FE_Q<dim>, false>::m_t;
    using cBaseMPI<dim, 2, dealii::FE_Q<dim>, false>::m_ti;
    using cBaseMPI<dim, 2, dealii::FE_Q<dim>, false>::m_N;
    using cBaseMPI<dim, 2, dealii::FE_Q<dim>, false>::m_omega;
    using cBaseMPI<dim, 2, dealii::FE_Q<dim>, false>::m_mu;
    using cBaseMPI<dim, 2, dealii::FE_Q<dim>, false>::m_gs;
    using cBaseMPI<dim, 2, dealii::FE_Q<dim>, false>::m_counter;
    using cBaseMPI<dim, 2, dealii::FE_Q<dim>, false>::m_computing_timer;
    using cBaseMPI<dim, 2, dealii::FE_Q<dim>, false>::m_ph;
    using cBaseMPI<dim, 2, dealii::FE_Q<dim>, false>::m_final_error;
    using cBaseMPI<dim, 2, dealii::FE_Q<dim>, false>::m_NA;
    using cBaseMPI<dim, 2, dealii::FE_Q<dim>, false>::m_Ndmu;
    using cBaseMPI<dim, 2, dealii::FE_Q<dim>, false>::m_dmu;
    using cBaseMPI<dim, 2, dealii::FE_Q<dim>, false>::m_QN1;
    using cBaseMPI<dim, 2, dealii::FE_Q<dim>, false>::m_res;
    using cBaseMPI<dim, 2, dealii::FE_Q<dim>, false>::m_resp;
    using cBaseMPI<dim, 2, dealii::FE_Q<dim>, false>::m_res_old;
    using cBaseMPI<dim, 2, dealii::FE_Q<dim>, false>::m_epsilon;
    using cBaseMPI<dim, 2, dealii::FE_Q<dim>, false>::m_maxiter;
    using cBaseMPI<dim, 2, dealii::FE_Q<dim>, false>::m_total_no_cells;
    using cBaseMPI<dim, 2, dealii::FE_Q<dim>, false>::m_total_no_active_cells;
    using cBaseMPI<dim, 2, dealii::FE_Q<dim>, false>::m_global_refinement;
  };

  /**
   * Constructor
   */
  template <int dim>
  MySolver<dim>::MySolver(const std::string xmlfilename) : cBaseMPI<dim, 2, dealii::FE_Q<dim>, false>(xmlfilename, m_FEQ)
  {
  }

#include "shared_1.h"
#include "grid.h"


  template <int dim>
  void MySolver<dim>::compute_contributions()
  {
    TimerOutput::Scope timing_section(m_computing_timer, "");

    this->update_workspace();

    std::cout << "N0 = " << MyRealTools::MPI::Particle_Number(mpi_communicator, this->m_DOF_Handler, this->m_FE, this->m_Workspace[0]) << "\n";
    std::cout << "N1 = " << MyRealTools::MPI::Particle_Number(mpi_communicator, this->m_DOF_Handler, this->m_FE, this->m_Workspace[1]) << "\n";

    CPotential<dim> Potential_fct(m_omega);
    const QGauss<dim> quadrature_formula(this->m_FE.degree + 1);
    FEValues<dim> fe_values(this->m_FE, quadrature_formula, update_gradients | update_values | update_JxW_values | update_quadrature_points);

    const unsigned n_q_points = quadrature_formula.size();

    vector<double> Psi_0(n_q_points);
    vector<double> Psi_1(n_q_points);
    vector<Tensor<1, dim>> Psi_0_grad(n_q_points);
    vector<Tensor<1, dim>> Psi_1_grad(n_q_points);

    std::array<double, 8> local_contributions;
    std::array<double, 8> total_contributions;
    local_contributions.fill(0);
    total_contributions.fill(0);

    typename DoFHandler<dim>::active_cell_iterator cell = this->m_DOF_Handler.begin_active(), endc = this->m_DOF_Handler.end();
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
          const double Q = Potential_fct.value(fe_values.quadrature_point(qp)) - m_mu;

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
    this->m_coeffs["t0^4"]      = 0.25 * m_gs.at(0) * total_contributions[0];
    this->m_coeffs["t0^1_t1^3"] = 0.25 * m_gs.at(0) * 4 * total_contributions[4];
    this->m_coeffs["t0^2_t1^2"] = 0.25 * m_gs.at(0) * 6 * total_contributions[2];
    this->m_coeffs["t0^3_t1^1"] = 0.25 * m_gs.at(0) * 4 * total_contributions[3];
    this->m_coeffs["t1^4"]      = 0.25 * m_gs.at(0) * total_contributions[1];


  }


  template <int dim>
  void MySolver<dim>::output_vector(LA::MPI::Vector& vec, string filename)
  {
    TimerOutput::Scope timing_section(m_computing_timer, "");

    this->m_constraints.distribute(vec);
    this->m_Workspace[0] = vec;

    DataOut<dim> data_out;
    data_out.attach_this->m_DOF_Handler(this->m_DOF_Handler);
    data_out.add_data_vector(this->m_Workspace[0], "vec");
    data_out.build_patches(gl_subdivisions);
    data_out.write_vtu_in_parallel(filename.c_str(), mpi_communicator);
  }

  template <int dim>
  int MySolver<dim>::DoIter(string path)
  {
    int retval = Status::SUCCESS;

    m_t[0] = m_ti;
    m_t[1] = m_ti;

    CPotential<dim> Potential(m_omega);

    this->do_linear_superposition();
    this->m_Workspace[0] = this->m_Psi_Ref;
    MyRealTools::MPI::AssembleRHS_L2gradient<dim>(this->m_DOF_Handler, this->m_FE, this->m_constraints, this->m_Workspace[0], Potential, m_mu, m_gs.at(0), m_res, this->m_System_RHS);
    m_res_old = m_res;
    m_counter = 0;
    bool bempty;

    do
    {
      BOOST_LOG_TRIVIAL(info) << std::string('-', 80);
      BOOST_LOG_TRIVIAL(info) << "- " << path << " - " << m_counter;

      this->m_Workspace[0] = this->m_Psi_Ref;
      MyRealTools::MPI::AssembleSystem_Jacobian<dim>(this->m_DOF_Handler, this->m_FE, this->m_constraints, this->m_Workspace[0], Potential, m_mu, m_gs.at(0), this->m_System_Matrix);
      bool bsucc = solve();
      if (!bsucc)
      {
        return Status::SINGULAR;
      }

      /*
            this->m_Workspace[0] = m_Psi_1;
            this->m_Workspace[1] = this->m_Search_Direction;
            MyRealTools::MPI::orthonormalize(mpi_communicator, this->m_DOF_Handler, fe, this->m_constraints, this->m_Workspace[1], this->m_Workspace[0], this->m_Workspace_NG );
            this->m_Search_Direction = this->m_Workspace_NG;
      */
      double tau = -0.5;
      this->m_Workspace[1] = this->m_Search_Direction;
      MyRealTools::MPI::compute_stepsize(mpi_communicator, this->m_DOF_Handler, this->m_FE, Potential, this->m_Workspace[0], this->m_Workspace[1], m_mu, m_gs.at(0), tau);
      std::cout << "NL2 = " << MyRealTools::MPI::Particle_Number(mpi_communicator, this->m_DOF_Handler, this->m_FE, this->m_Workspace[1]) << "\n";


      this->m_Psi[1].add(tau * m_t[1] / fabs(m_t[1]), this->m_Search_Direction);
      this->m_constraints.distribute(this->m_Psi[1]);

      bempty = this->find_ortho_min();

      if (bempty)
      {
        return Status::FAILED;
      }

      this->do_linear_superposition();
      this->m_Workspace[0] = this->m_Psi_Ref;
      MyRealTools::MPI::AssembleRHS_L2gradient<dim>(this->m_DOF_Handler, this->m_FE, this->m_constraints, this->m_Workspace[0], Potential, m_mu, m_gs.at(0), m_res, this->m_System_RHS);

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
      // m_table.insert(cols, MyTable::MU, m_mu);
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

    this->particle_number(this->m_Psi_Ref);

    // Standard Newton
    do
    {
      BOOST_LOG_TRIVIAL(info) << std::string('-', 80);
      BOOST_LOG_TRIVIAL(info) << "-- " << path << " - " << m_counter;

      this->m_Workspace[0] = this->m_Psi_Ref;
      MyRealTools::MPI::AssembleSystem_Jacobian<dim>(this->m_DOF_Handler, this->m_FE, this->m_constraints, this->m_Workspace[0], Potential, m_mu, m_gs.at(0), this->m_System_Matrix);
      bool bsucc = solve();
      if (!bsucc)
      {
        return Status::SINGULAR;
      }

      double tau;
      this->m_Workspace[1] = this->m_Search_Direction;
      MyRealTools::MPI::compute_stepsize(mpi_communicator, this->m_DOF_Handler, this->m_FE, Potential, this->m_Workspace[0], this->m_Workspace[1], m_mu, m_gs.at(0), tau);

      this->m_Psi_Ref.add(tau * m_t[1] / fabs(m_t[1]), this->m_Search_Direction);
      this->m_constraints.distribute(this->m_Psi_Ref);

      this->m_Workspace[0] = this->m_Psi_Ref;
      MyRealTools::MPI::AssembleRHS_L2gradient<dim>(this->m_DOF_Handler, this->m_FE, this->m_constraints, this->m_Workspace[0], Potential, m_mu, m_gs.at(0), m_res, this->m_System_RHS);

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
      // m_table.insert(cols, MyTable::MU, m_mu);
      // m_table.insert(cols, MyTable::GS, m_gs.at(0));
      // m_table.insert(cols, MyTable::t1, m_t[0]);
      // m_table.insert(cols, MyTable::t2, m_t[1]);
      // m_table.insert(cols, MyTable::l2norm_t, this->l2norm_t());
      // m_table.insert(cols, MyTable::PARTICLE_NUMBER, m_N);

      m_counter++;

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
    }
    while (true);

    this->m_Workspace[0] = this->m_Psi_Ref;
    m_N = MyRealTools::MPI::Particle_Number(mpi_communicator, this->m_DOF_Handler, this->m_FE, this->m_Workspace[0]);

    if (m_N < 1e-5)
    {
      retval = Status::ZERO_SOL;
    }
    return retval;
  }


  template <int dim>
  void MySolver<dim>::run2b()
  {
    string path;
    char shellcmd[255];
    double T, N, W;
    int status;

    make_grid_custom();
    this->setup_system();

    CEigenfunctions<dim> Ef1(m_QN1, m_omega);
    //CEigenfunctions<dim> Ef2( m_QN2, m_omega );
    CPotential<dim> Potential(m_omega);

    VectorTools::interpolate(this->m_DOF_Handler, Ef1, this->m_Psi[0]);

    this->m_Workspace[0] = this->m_Psi[0];
    this->m_Psi[0] *= 1.0 / sqrt(MyRealTools::MPI::Particle_Number(mpi_communicator, this->m_DOF_Handler, this->m_FE, this->m_Workspace[0]));
    this->m_Psi[1] = 0;

    this->m_Workspace[0] = this->m_Psi[0];
    MyRealTools::MPI::Compute_E_lin(mpi_communicator, this->m_DOF_Handler, this->m_FE, this->m_constraints, Potential, this->m_Workspace[0], T, N, W);
    double m_mu_0 = T / N;
    m_mu = ceil(10.0 * m_mu_0) / 10.0 + m_gs.at(0) / fabs(m_gs.at(0)) * m_dmu;

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
      // sqrt((m_mu*N-T)/(m_gs*W)); if this->m_Psi[1] == 0
      // sqrt((m_mu*N-T)/(4.0*m_gs*W)); if this->m_Psi[1] == m_Psi_1
      m_ti = sqrt((m_mu * N - T) / (m_gs.at(0) * W));

      BOOST_LOG_TRIVIAL(info) << "T = " << T;
      BOOST_LOG_TRIVIAL(info) << "N = " << N;
      BOOST_LOG_TRIVIAL(info) << "W = " << W;
      BOOST_LOG_TRIVIAL(info) << "m_mu = " << m_mu;
      BOOST_LOG_TRIVIAL(info) << "m_ti = " << m_ti;

      status = DoIter(path);

      // columns& cols = m_results.new_line();
      // m_results.insert(cols, MyTable::MU, m_mu);
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
      m_mu += m_gs.at(0) / fabs(m_gs.at(0)) * m_dmu;
      this->m_Workspace[0] = this->m_Psi[0];
      MyRealTools::MPI::Compute_E_lin(mpi_communicator, this->m_DOF_Handler, this->m_FE, this->m_constraints, Potential, this->m_Workspace[0], T, N, W);

      BOOST_LOG_TRIVIAL(info) << "T = " << T << endl;
      BOOST_LOG_TRIVIAL(info) << "N = " << N << endl;
      BOOST_LOG_TRIVIAL(info) << "W = " << W << endl;
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
    params.GetParameter("parameter.spatial_dimension", dim);
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
