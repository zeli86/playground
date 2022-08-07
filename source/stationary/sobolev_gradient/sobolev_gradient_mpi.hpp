

#include "MyParameterHandler.h"
#include "UtilsRealWavefunction.hpp"
#include "GPUtilsRealWavefunction.hpp"
#include <limits>

namespace solver::mpi::stationary
{
  using namespace std;
  using namespace dealii;

  enum Status { SUCCESS, FAILED, ZERO_SOL, SLOW_CONV };

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
      return m_oDofHandler;
    }

    dealii::FE_Q<dim>& get_fe()
    {
      return m_oFe;
    }

    dealii::AffineConstraints<double>& get_constraints()
    {
      return m_oConstraints;
    }

    int DoIter(string);

    void make_grid();
    void setup_system();
    void assemble_system();
    void compute_psi_sobolev();
    void project_gradient();
    void solve();
    void save(string);
    void output_results(string, string);

  private:

    std::vector<double> m_omega{0.25, 0.25, 0.25};
    double m_rEpsilon{1e-10};
    double m_res{std::numeric_limits<double>::max()};
    double m_N{0.0};
    double m_rMu{0};
    double m_rG{1};
    bool m_root{false};
    int m_rank{-1};
    int m_NA{10};
    int m_iMaxIter{100};

    MyParameterHandler m_oParameters;
    MPI_Comm mpi_communicator;
    parallel::distributed::Triangulation<dim> m_oTriangulation;
    FE_Q<dim> m_oFe;
    DoFHandler<dim> m_oDofHandler;
    AffineConstraints<double> m_oConstraints;
    IndexSet m_oLocallyOwnedDofs;
    IndexSet m_oLocallyRelevantDofs;

    LA::MPI::SparseMatrix m_oMatrix;
    LA::MPI::Vector m_vRhs;
    LA::MPI::Vector m_vPhi;
    LA::MPI::Vector m_vSobGradient;
    LA::MPI::Vector m_vPhiSobolev;
    LA::MPI::Vector m_vWorkspace1;
    LA::MPI::Vector m_vWorkspace2;
    LA::MPI::Vector m_vWorkspace3;

    int m_iGlobalRefinement{7};

    string m_guess_str;
  };
}