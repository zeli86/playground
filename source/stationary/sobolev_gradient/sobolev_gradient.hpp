
#include "MyParameterHandler.h"
#include "UtilsRealWavefunction.hpp"
#include "GPUtilsRealWavefunction.hpp"
#include <limits>

namespace solver::stationary
{
  using namespace std;
  using namespace dealii;

  template <int dim>
  class CSobolevGradient : public utils::real_wavefunction::IRealWavefunction<dim>
  {
  public:
    explicit CSobolevGradient(const std::string&);
    virtual ~CSobolevGradient();

    void run();

  protected:
    dealii::DoFHandler<dim>& get_dof_handler()
    {
      return dof_handler;
    }

    dealii::FE_Q<dim>& get_fe()
    {
      return fe;
    }

    dealii::AffineConstraints<double>& get_constraints()
    {
      return m_oConstraints;
    }
    void make_grid();
    void setup_system();
    void assemble_system();
    void compute_Psi_sob();
    void save(string);
    void project_gradient();

    void solve();
    int DoIter(string = "");

    void output_results(string);

  protected:
    double m_res, m_res_old, m_resp;
    double m_N = 0.0;
    double m_rG{1.0};
    vector<double> m_omega;
    double m_rEpsilon{1e-10};
    unsigned m_NA{100};

    MyParameterHandler m_oParameters;
    Triangulation<dim> triangulation;
    FE_Q<dim> fe{2};
    DoFHandler<dim> dof_handler;

    AffineConstraints<double> m_oConstraints;
    IndexSet m_oLocallyOwnedDofs;
    IndexSet m_oLocallyRelevantDofs;

    SparsityPattern m_sparsity_pattern;
    SparseMatrix<double> m_system_matrix;
    Vector<double> m_vPhi;
    Vector<double> m_vPhiSobolev;
    Vector<double> m_vL2Gradient;
    Vector<double> m_vSobolevGradient;
    Vector<double> m_solution;
    Vector<double> m_error_per_cell;
    int m_iGlobalRefinement{7};
  };
}
