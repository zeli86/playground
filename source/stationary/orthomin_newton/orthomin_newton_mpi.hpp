//#include "default_includes.h"

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

#include <deal.II/dofs/dof_tools.h>
#include <deal.II/grid/grid_generator.h>

//#include "global.h"
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
    virtual ~MySolver();

    void run2b();

    double m_I[8];

    std::map<std::string, double> m_coeffs;

  protected:

    MPI_Comm mpi_communicator;

    int DoIter(string = "");

    void make_grid();
    void setup_system();

    //void estimate_error(double&);

    bool solve();
    void compute_contributions();
    int find_ortho_min();

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
}