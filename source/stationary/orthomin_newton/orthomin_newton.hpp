
#pragma once

//#include "global.h"
#include "functions.h"
//#include "ref_pt_list.h"
#include "MyParameterHandler.h"
#include "UtilsRealWavefunction.hpp"
#include "GPUtilsRealWavefunction.hpp"
#include "muParser.h"
#include <limits>

namespace BreedSolver_1
{

  enum Status { SUCCESS, FAILED, ZERO_SOL, SLOW_CONV, MAXITER, SINGULAR, CONTINUE };

  template <int dim>
  class MySolver : public utils::real_wavefunction::IRealWavefunction<dim>
  {
  public:
    explicit MySolver(const std::string&);
    ~MySolver();

    void run();

  private:

    dealii::DoFHandler<dim>& get_dof_handler()
    {
      return m_DOF_Handler;
    }

    dealii::FE_Q<dim>& get_fe()
    {
      return m_FE;
    }

    dealii::AffineConstraints<double>& get_constraints()
    {
      return m_constraints;
    }

    void make_grid();
    void setup_system();
    void assemble_system();
    void assemble_rhs();
    int DoIter(string = "");

    void solve();
    void compute_contributions();
    void compute_E_lin(Vector<double>&, double&, double&, double&);
    void output_results(string, string = "step");
    void output_guess();

    typedef typename dealii::Vector<double> aVector;
    typedef typename dealii::SparseMatrix<double> aSparseMatrix;
    typedef typename dealii::Triangulation<dim> aTriangulation;

    aVector m_Psi_Ref;
    aVector m_System_RHS;
    aVector m_Search_Direction;
    dealii::Vector<double> m_error_per_cell;

    std::array<aVector, 2> m_Psi;
    std::array<aVector, 2> m_Workspace;

    aSparseMatrix m_System_Matrix;

    AffineConstraints<double> m_constraints;

    double m_t[2];
    double m_t_guess[2];

    double m_res{std::numeric_limits<double>::max()};
    double m_ti;
    double m_final_error;
    double m_N{0};
    double m_rMu{0};
    double m_dmu{0.1};
    double m_gs{1};
    vector<double> m_omega;
    vector<double> m_epsilon;

    int m_rank = -1;

    unsigned m_counter;
    unsigned m_maxiter{500};
    unsigned m_global_refinement{9};
    unsigned m_total_no_cells;
    unsigned m_total_no_active_cells;
    unsigned m_NA{10};
    unsigned m_Ndmu{10};
    unsigned m_QN1[3];

    MyParameterHandler m_ph;

    aTriangulation m_Triangulation;
    FE_Q<dim> m_FE;
    DoFHandler<dim> m_DOF_Handler;

    std::map<std::string, double> m_coeffs;
  };

}