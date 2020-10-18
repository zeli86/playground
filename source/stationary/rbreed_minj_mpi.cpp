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

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>

#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <locale>
#include <limits>
#include <cmath>
#include <array>
#include <unistd.h>

#include <nlopt.hpp>

#include "global.h"
#include "mpi.h"
#include "functions.h"
#include "MyParameterHandler.h"
#include "my_table.h"
#include "MyRealTools.h"

#define STR1(x) #x
#define STR2(x) STR1(x)

namespace BreedSolver
{
  using namespace std;
  using namespace dealii;

#include "cBaseMPI.h"

  template <int dim>
  class MySolver; // class forward declaration

  template <int dim>
  class MySolver : public cBaseMPI<dim, 2, dealii::FE_Q<dim>>
  {
  public:
    explicit MySolver(const std::string);

    void run2b();

    void set_t (const double a, const double b)
    {
      m_t[0] = a;
      m_t[1] = b;
    };

    double m_I[8]; // remove me

  protected:

    dealii::FE_Q<dim> m_FEQ;

    using cBaseMPI<dim, 2, dealii::FE_Q<dim>>::mpi_communicator;
    using cBaseMPI<dim, 2, dealii::FE_Q<dim>>::m_root;
    using cBaseMPI<dim, 2, dealii::FE_Q<dim>>::m_rank;
    using cBaseMPI<dim, 2, dealii::FE_Q<dim>>::m_t;
    using cBaseMPI<dim, 2, dealii::FE_Q<dim>>::m_ti;
    using cBaseMPI<dim, 2, dealii::FE_Q<dim>>::m_N;
    using cBaseMPI<dim, 2, dealii::FE_Q<dim>>::m_omega;
    using cBaseMPI<dim, 2, dealii::FE_Q<dim>>::m_mu;
    using cBaseMPI<dim, 2, dealii::FE_Q<dim>>::m_gs;
    using cBaseMPI<dim, 2, dealii::FE_Q<dim>>::m_counter;
    using cBaseMPI<dim, 2, dealii::FE_Q<dim>>::pcout;
    using cBaseMPI<dim, 2, dealii::FE_Q<dim>>::m_computing_timer;
    using cBaseMPI<dim, 2, dealii::FE_Q<dim>>::m_ph;
    using cBaseMPI<dim, 2, dealii::FE_Q<dim>>::m_final_error;
    using cBaseMPI<dim, 2, dealii::FE_Q<dim>>::m_NA;
    using cBaseMPI<dim, 2, dealii::FE_Q<dim>>::m_Ndmu;
    using cBaseMPI<dim, 2, dealii::FE_Q<dim>>::m_dmu;
    using cBaseMPI<dim, 2, dealii::FE_Q<dim>>::m_QN1;
    using cBaseMPI<dim, 2, dealii::FE_Q<dim>>::m_res;
    using cBaseMPI<dim, 2, dealii::FE_Q<dim>>::m_resp;
    using cBaseMPI<dim, 2, dealii::FE_Q<dim>>::m_res_old;
    using cBaseMPI<dim, 2, dealii::FE_Q<dim>>::m_epsilon;
    using cBaseMPI<dim, 2, dealii::FE_Q<dim>>::m_maxiter;
    using cBaseMPI<dim, 2, dealii::FE_Q<dim>>::m_total_no_cells;
    using cBaseMPI<dim, 2, dealii::FE_Q<dim>>::m_total_no_active_cells;
    using cBaseMPI<dim, 2, dealii::FE_Q<dim>>::m_global_refinement;

    void estimate_error(double&);

    int DoIter(string = "");

    void make_grid();
    void make_grid_custom();

    bool solve();
    void compute_E_lin(LA::MPI::Vector&, double&, double&, double&);
    void find_min_J();

    void compute_contributions();

    string m_guess_str;
    MyTable m_table;
    MyTable m_results;
  };


#include "shared_1.h"
#include "grid.h"


  /**
   * Constructor
   */
  template <int dim>
  MySolver<dim>::MySolver(const string xmlfilename) : m_FEQ(gl_degree_fe), cBaseMPI<dim, 2, dealii::FE_Q<dim>> (xmlfilename, m_FEQ)
  {
    try
    {
      m_omega = m_ph.Get_Physics("omega");
      m_gs = m_ph.Get_Physics("gs_1", 0);
      m_QN1[0] = int (m_ph.Get_Physics("QN1", 0));
      m_QN1[1] = int (m_ph.Get_Physics("QN1", 1));
      m_QN1[2] = int (m_ph.Get_Physics("QN1", 2));

      m_ti = m_ph.Get_Algorithm("ti", 0);
      m_t[0] = m_ti;
      m_t[1] = m_ti;
      m_NA = int (m_ph.Get_Algorithm("NA", 0));
      m_Ndmu = m_ph.Get_Algorithm("Ndmu", 0);
      m_dmu = m_ph.Get_Algorithm("dmu", 0);
      m_epsilon = m_ph.Get_Algorithm("epsilon");
    }
    catch (const std::string& info)
    {
      std::cerr << info << endl;
      MPI_Abort(mpi_communicator, 0);
    }

    MPI_Comm_rank(mpi_communicator, &m_rank);
  }

  template <int dim>
  void MySolver<dim>::compute_contributions()
  {
    TimerOutput::Scope timing_section(m_computing_timer, "compute contributions");

    this->update_workspace();

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
    this->m_coeffs["t0^4"]      = 0.25 * m_gs * total_contributions[0];
    this->m_coeffs["t0^1_t1^3"] = 0.25 * m_gs * 4 * total_contributions[4];
    this->m_coeffs["t0^2_t1^2"] = 0.25 * m_gs * 6 * total_contributions[2];
    this->m_coeffs["t0^3_t1^1"] = 0.25 * m_gs * 4 * total_contributions[3];
    this->m_coeffs["t1^4"]      = 0.25 * m_gs * total_contributions[1];
  }


  template<int dim>
  void MySolver<dim>::find_min_J()
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
      // std::cout << "found minimum at f(" << t[0] << "," << t[1] << ") = " << std::setprecision(10) << minf << std::endl;
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
  int MySolver<dim>::DoIter(string path)
  {
    int retval = Status::SUCCESS;

    m_table.clear();

    m_t[0] = m_ti;
    m_t[1] = m_ti;

    CPotential<dim> Potential(m_omega);

    this->do_linear_superposition();
    this->m_Workspace[0] = this->m_Psi_Ref;
    MyRealTools::MPI::AssembleRHS_L2gradient<dim> (this->m_DOF_Handler, this->m_FE, this->m_constraints, this->m_Workspace[0], Potential, m_mu, m_gs, m_res, this->m_System_RHS);
    m_res_old = m_res;
    m_counter = 0;
    bool bempty;

    do
    {
      pcout << "--------------------------------------------------------------------------------" << endl;
      pcout << "- " << path << " - " << m_counter << endl;

      this->m_Workspace[0] = this->m_Psi_Ref;
      MyRealTools::MPI::AssembleSystem_Jacobian<dim> (this->m_DOF_Handler, this->m_FE, this->m_constraints, this->m_Workspace[0], Potential, m_mu, m_gs, this->m_System_Matrix);
      bool bsucc = solve();
      if (!bsucc)
      {
        return Status::SINGULAR;
      }

      double tau = -0.5;
      this->m_Workspace[1] = this->m_Search_Direction;
      //MyRealTools::MPI::compute_stepsize(mpi_communicator, this->m_DOF_Handler, this->m_FE, Potential, this->m_Workspace[0], this->m_Workspace[1], m_mu, m_gs, tau);

      this->m_Psi[1].add(tau * m_t[1] / fabs(m_t[1]), this->m_Search_Direction);
      this->m_constraints.distribute(this->m_Psi[1]);

      find_min_J();

      // if (bempty)
      // {
      //   return Status::FAILED;
      // }

      this->do_linear_superposition();
      this->m_Workspace[0] = this->m_Psi_Ref;
      MyRealTools::MPI::AssembleRHS_L2gradient<dim> (this->m_DOF_Handler, this->m_FE, this->m_constraints, this->m_Workspace[0], Potential, m_mu, m_gs, m_res, this->m_System_RHS);

      m_resp = m_res_old - m_res;
      m_res_old = m_res;

      // if (m_counter % m_NA == 0)
      // {
      this->output_results(path);
      // }

      columns& cols = m_table.new_line();
      m_table.insert(cols, MyTable::COUNTER, double (m_counter));
      m_table.insert(cols, MyTable::RES, m_res);
      m_table.insert(cols, MyTable::RESP, m_resp);
      m_table.insert(cols, MyTable::STEPSIZE, tau);
      m_table.insert(cols, MyTable::MU, m_mu);
      m_table.insert(cols, MyTable::GS, m_gs);
      m_table.insert(cols, MyTable::t1, m_t[0]);
      m_table.insert(cols, MyTable::t2, m_t[1]);
      m_table.insert(cols, MyTable::l2norm_t, this->l2norm_t());

      m_counter++;
      if (m_root)
      {
        cout << m_table;
      }
      if (m_res < m_epsilon[0])
      {
        retval = Status::SUCCESS;
        break;
      }
      if (this->l2norm_t() < 1e-4)
      {
        retval = Status::ZERO_SOL;
        break;
      }
      if (m_counter == m_maxiter)
      {
        retval = Status::MAXITER;
        break;
      }
      if (isnan(m_res))
      {
        retval = Status::FAILED;
        break;
      }

    }
    while (true);

    // Standard Newton
    do
    {
      pcout << "--------------------------------------------------------------------------------" << endl;
      pcout << "-- " << path << " - " << m_counter << endl;

      this->m_Workspace[0] = this->m_Psi_Ref;
      MyRealTools::MPI::AssembleSystem_Jacobian<dim> (this->m_DOF_Handler, this->m_FE, this->m_constraints, this->m_Workspace[0], Potential, m_mu, m_gs, this->m_System_Matrix);
      bool bsucc = solve();
      if (!bsucc)
      {
        return Status::SINGULAR;
      }

      double tau = 0.1;
      this->m_Workspace[1] = this->m_Search_Direction;
      MyRealTools::MPI::compute_stepsize(mpi_communicator, this->m_DOF_Handler, this->m_FE, Potential, this->m_Workspace[0], this->m_Workspace[1], m_mu, m_gs, tau);

      this->m_Psi_Ref.add(tau * m_t[1] / fabs(m_t[1]), this->m_Search_Direction);
      this->m_Psi_Ref.add(tau, this->m_Search_Direction);
      this->m_constraints.distribute(this->m_Psi_Ref);

      this->m_Workspace[0] = this->m_Psi_Ref;
      MyRealTools::MPI::AssembleRHS_L2gradient<dim> (this->m_DOF_Handler, this->m_FE, this->m_constraints, this->m_Workspace[0], Potential, m_mu, m_gs, m_res, this->m_System_RHS);

      m_resp = m_res_old - m_res;
      m_res_old = m_res;

      if (m_counter % m_NA == 0)
      {
        this->output_results(path);
      }

      columns& cols = m_table.new_line();
      m_table.insert(cols, MyTable::COUNTER, double (m_counter));
      m_table.insert(cols, MyTable::RES, m_res);
      m_table.insert(cols, MyTable::RESP, m_resp);
      m_table.insert(cols, MyTable::STEPSIZE, tau);
      m_table.insert(cols, MyTable::MU, m_mu);
      m_table.insert(cols, MyTable::GS, m_gs);
      m_table.insert(cols, MyTable::t1, m_t[0]);
      m_table.insert(cols, MyTable::t2, m_t[1]);
      m_table.insert(cols, MyTable::l2norm_t, this->l2norm_t());
      m_table.insert(cols, MyTable::PARTICLE_NUMBER, m_N);

      m_counter++;

      if (m_root)
      {
        cout << m_table;
      }
      if (m_res < m_epsilon[1])
      {
        retval = Status::SUCCESS;
        break;
      }
      if (m_resp < 0)
      {
        retval = Status::SUCCESS;
        break;
      }
      if (m_counter == m_maxiter)
      {
        retval = Status::MAXITER;
        break;
      }
    }
    while (true);

    this->m_Workspace[0] = this->m_Psi_Ref;
    m_N = MyRealTools::MPI::Particle_Number(mpi_communicator, this->m_DOF_Handler, this->m_FE, this->m_Workspace[0]);

    if (m_N < 1e-5)
    {
      retval = Status::ZERO_SOL;
    }

    string filename = path + "log.csv";
    if (m_root)
    {
      m_table.dump_2_file(filename);
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
    VectorTools::interpolate(this->m_DOF_Handler, Ef1, this->m_Psi[0]);

    this->m_Workspace[0] = this->m_Psi[0];
    this->m_Psi[0] *= 1.0 / sqrt(MyRealTools::MPI::Particle_Number(mpi_communicator, this->m_DOF_Handler, this->m_FE, this->m_Workspace[0]));
    this->m_Psi[1] = 0;

    compute_E_lin(this->m_Psi[0], T, N, W);
    double m_mu_0 = T / N;
    m_mu = ceil(10.0 * m_mu_0) / 10.0 + m_gs / fabs(m_gs) * m_dmu;

    this->output_guess();
    m_results.clear();
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
      //m_ti = sqrt((m_mu*N-T)/(4.0*m_gs*W)); // if m_Psi_2 == m_Psi_1
      m_ti = sqrt((m_mu * N - T) / (m_gs * W));

      pcout << "T = " << T << endl;
      pcout << "N = " << N << endl;
      pcout << "W = " << W << endl;
      pcout << "m_mu = " << m_mu << endl;
      pcout << "m_ti = " << m_ti << endl;

      status = DoIter(path);

      columns& cols = m_results.new_line();
      m_results.insert(cols, MyTable::MU, m_mu);
      m_results.insert(cols, MyTable::GS, m_gs);
      m_results.insert(cols, MyTable::PARTICLE_NUMBER, m_N);
      m_results.insert(cols, MyTable::COUNTER, double (m_counter));
      m_results.insert(cols, MyTable::STATUS, double (status));

      if (status == Status::SUCCESS)
      {
        this->output_results(path, "final");
        string filename = path + "final.bin";
        this->save(filename);
        this->dump_info_xml(path);
        this->m_Psi[0] = this->m_Psi_Ref;
        this->m_Psi[1] = 0;
      }
      else if (status == Status::SLOW_CONV)
      {
        this->m_Psi[1] = 0;
      }
      else
      {
        break;
      }
      compute_E_lin(this->m_Psi_Ref, T, N, W);    // TODO: kommentier mich aus, falls ich kein nehari reset habe
      m_mu += m_gs / fabs(m_gs) * m_dmu;
    }
    if (m_root)
    {
      m_results.dump_2_file("results.csv");
    }
  }
} // end of namespace


int main(int argc, char* argv[])
{
  using namespace dealii;
  deallog.depth_console(0);

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);
  {
    BreedSolver::MySolver<2> solver("params.xml");
    solver.run2b();
  }
  return EXIT_SUCCESS;
}