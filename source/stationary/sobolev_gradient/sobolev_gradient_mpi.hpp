
#include <deal.II/lac/generic_linear_algebra.h>

namespace LA
{
  using namespace dealii::LinearAlgebraPETSc;
}

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>

#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_precondition.h>
#include <deal.II/lac/sparsity_tools.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_out.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_q.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>
//#include <deal.II/numerics/error_estimator.h>
//#include <deal.II/numerics/derivative_approximation.h>

#include <deal.II/base/utilities.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/function_parser.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/timer.h>

#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/solution_transfer.h>

#include "MyParameterHandler.h"
#include "UtilsRealWavefunction.hpp"
#include "GPUtilsRealWavefunction.hpp"

namespace BreedSolver
{
  using namespace std;
  using namespace dealii;

  enum Status { SUCCESS, FAILED, ZERO_SOL, SLOW_CONV };

  template <int dim>
  class MySolver
  {
  public:
    explicit MySolver(const std::string&);
    virtual ~MySolver();

    void run();

    MPI_Comm mpi_communicator;
  protected:

    double m_res, m_res_old, m_resp;
    double m_final_error = 0.0;
    double m_N = 0.0;

    double m_rMu = 0;
    double m_rG = 1;
    std::vector<double> m_gs{1, 1};
    std::vector<double> m_omega{3, 1};
    std::vector<double> m_epsilon{2, 1e-10};

    bool m_root = false;
    int m_rank = -1;

    unsigned m_counter = 0;
    unsigned m_total_no_cells;
    unsigned m_total_no_active_cells;
    unsigned m_NA;

    int DoIter(string = "");

    void make_grid();
    void setup_system();
    void assemble_rhs();
    void assemble_system();
    void compute_Psi_sob();
    void save(string);
    void Project_gradient();

    void solve();
    void estimate_error(double&);

    void output_results(string, string = "step");
    void output_guess();

    MyParameterHandler m_oParameters;
    parallel::distributed::Triangulation<dim> m_oTriangulation;
    FE_Q<dim> fe;
    DoFHandler<dim> dof_handler;
    IndexSet locally_owned_dofs;
    IndexSet locally_relevant_dofs;
    AffineConstraints<double> constraints;

    LA::MPI::SparseMatrix m_system_matrix;
    LA::MPI::Vector m_system_rhs;
    LA::MPI::Vector m_Psi;
    LA::MPI::Vector m_Psi_C_ghosted;
    LA::MPI::Vector m_sob_grad;
    LA::MPI::Vector m_Psi_sob;
    LA::MPI::Vector m_workspace_1;
    LA::MPI::Vector m_workspace_2;
    LA::MPI::Vector m_workspace_3;
    Vector<double> m_error_per_cell;

    int m_iGlobalRefinement = 7;

    string m_guess_str;
  };
}