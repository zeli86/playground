
#include "default_includes.h"

#include <boost/log/trivial.hpp>
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
#include "MyRealTools.h"
#include "MyLogging.h"

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

    double m_I[8]; // remove me

  protected:

    dealii::FE_Q<dim> m_FEQ{gl_degree_fe};

    void estimate_error(double&);

    int DoIter(string = "");

    void make_grid();
    void make_grid_custom();

    bool solve();
    void find_min_J();

    void compute_contributions();

    string m_guess_str;
  };

#include "shared_1.h"
#include "grid.h"


  /**
   * Constructor
   */
  template <int dim>
  MySolver<dim>::MySolver(const string xmlfilename) : cBaseMPI<dim, 2, dealii::FE_Q<dim>> (xmlfilename, m_FEQ)
  {
    if (m_root)
    {
      disableLogging();
    }
  }

  template <int dim>
  void MySolver<dim>::compute_contributions()
  {

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
      BOOST_LOG_TRIVIAL(debug) << "found minimum at f(" << t[0] << "," << t[1] << ") = " << std::setprecision(10) << minf;
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

    m_t[0] = m_ti;
    m_t[1] = m_ti;

    CPotential<dim> Potential(m_omega);

    this->do_linear_superposition();
    this->m_Workspace[0] = this->m_Psi_Ref;
    MyRealTools::MPI::AssembleRHS_L2gradient<dim> (this->m_DOF_Handler, this->m_FE, this->m_constraints, this->m_Workspace[0], Potential, m_rMu, m_rG, m_res, this->m_System_RHS);
    m_res_old = m_res;
    m_counter = 0;
    bool bempty;

    do
    {
      BOOST_LOG_TRIVIAL(info) << std::string('-', 80);
      BOOST_LOG_TRIVIAL(info) << "- " << path << " - " << m_counter;

      this->m_Workspace[0] = this->m_Psi_Ref;
      MyRealTools::MPI::AssembleSystem_Jacobian<dim> (this->m_DOF_Handler, this->m_FE, this->m_constraints, this->m_Workspace[0], Potential, m_rMu, m_rG, this->m_System_Matrix);
      bool bsucc = solve();
      if (!bsucc)
      {
        return Status::SINGULAR;
      }

      double tau = -0.5;
      this->m_Workspace[1] = this->m_Search_Direction;
      MyRealTools::MPI::compute_stepsize(mpi_communicator, this->m_DOF_Handler, this->m_FE, Potential, this->m_Workspace[0], this->m_Workspace[1], m_rMu, m_rG, tau);

      this->m_Psi[1].add(tau * m_t[1] / fabs(m_t[1]), this->m_Search_Direction);
      this->m_constraints.distribute(this->m_Psi[1]);

      find_min_J();

      // if (bempty)
      // {
      //   return Status::FAILED;
      // }

      this->do_linear_superposition();
      this->m_Workspace[0] = this->m_Psi_Ref;
      MyRealTools::MPI::AssembleRHS_L2gradient<dim> (this->m_DOF_Handler, this->m_FE, this->m_constraints, this->m_Workspace[0], Potential, m_rMu, m_rG, m_res, this->m_System_RHS);

      m_resp = m_res_old - m_res;
      m_res_old = m_res;

      // if (m_counter % m_NA == 0)
      // {
      //this->output_results(path);
      // }

      m_counter++;
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
      BOOST_LOG_TRIVIAL(info) << std::string('-', 80);
      BOOST_LOG_TRIVIAL(info) << "-- " << path << " - " << m_counter << endl;

      this->m_Workspace[0] = this->m_Psi_Ref;
      MyRealTools::MPI::AssembleSystem_Jacobian<dim> (this->m_DOF_Handler, this->m_FE, this->m_constraints, this->m_Workspace[0], Potential, m_rMu, m_rG, this->m_System_Matrix);
      bool bsucc = solve();
      if (!bsucc)
      {
        return Status::SINGULAR;
      }

      double tau = 0.1;
      this->m_Workspace[1] = this->m_Search_Direction;
      MyRealTools::MPI::compute_stepsize(mpi_communicator, this->m_DOF_Handler, this->m_FE, Potential, this->m_Workspace[0], this->m_Workspace[1], m_rMu, m_rG, tau);

      this->m_Psi_Ref.add(tau * m_t[1] / fabs(m_t[1]), this->m_Search_Direction);
      this->m_Psi_Ref.add(tau, this->m_Search_Direction);
      this->m_constraints.distribute(this->m_Psi_Ref);

      this->m_Workspace[0] = this->m_Psi_Ref;
      MyRealTools::MPI::AssembleRHS_L2gradient<dim> (this->m_DOF_Handler, this->m_FE, this->m_constraints, this->m_Workspace[0], Potential, m_rMu, m_rG, m_res, this->m_System_RHS);

      m_resp = m_res_old - m_res;
      m_res_old = m_res;

      if (m_counter % m_NA == 0)
      {
        this->output_results(path);
      }

      m_counter++;

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

    CPotential<dim> Potential(m_omega);

    this->m_Workspace[0] = this->m_Psi[0];
    this->m_Psi[0] *= 1.0 / sqrt(MyRealTools::MPI::Particle_Number(mpi_communicator, this->m_DOF_Handler, this->m_FE, this->m_Workspace[0]));
    this->m_Psi[1] = 0;

    this->m_Workspace[0] = this->m_Psi[0];
    MyRealTools::MPI::Compute_E_lin(mpi_communicator, this->m_DOF_Handler, this->m_FE, this->m_constraints, Potential, this->m_Workspace[0], T, N, W);
    double m_rMu_0 = T / N;
    m_rMu = ceil(10.0 * m_rMu_0) / 10.0 + m_rG / fabs(m_rG) * m_dmu;

    this->output_guess();
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
      m_ti = sqrt((m_rMu * N - T) / (m_rG * W));

      BOOST_LOG_TRIVIAL(info) << "T = " << T;
      BOOST_LOG_TRIVIAL(info) << "N = " << N;
      BOOST_LOG_TRIVIAL(info) << "W = " << W;
      BOOST_LOG_TRIVIAL(info) << "m_rMu = " << m_rMu;
      BOOST_LOG_TRIVIAL(info) << "m_ti = " << m_ti;

      status = DoIter(path);

      if (status == Status::SUCCESS)
      {
        this->output_results(path, "final");
        string filename = path + "final.bin";
        this->save(filename);
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
      this->m_Workspace[0] = this->m_Psi[0];
      MyRealTools::MPI::Compute_E_lin(mpi_communicator, this->m_DOF_Handler, this->m_FE, this->m_constraints, Potential, this->m_Workspace[0], T, N, W);
      m_rMu += m_rG / fabs(m_rG) * m_dmu;
    }
  }
} // end of namespace


int main(int argc, char* argv[])
{
  using namespace dealii;
  deallog.depth_console(0);

  initLogging("log.txt");

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);
  {
    BreedSolver::MySolver<2> solver("params.xml");
    solver.run2b();
  }
  return EXIT_SUCCESS;
}