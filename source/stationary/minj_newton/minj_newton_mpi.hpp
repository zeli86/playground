
#include "MyParameterHandler.h"
#include "UtilsComplexWavefunction.hpp"
#include "GPUtilsComplexWavefunction.hpp"

template <int dim>
class MySolver : public utils::complex_wavefunction::IComplexWavefunction<dim>
{
public:
  explicit MySolver(const std::string);
  virtual ~MySolver() = default;

  void run2b();

protected:

  void setup_system();


  void compute_contributions();

  int DoIter(string = "");

  void make_grid_custom();
  void assemble_rhs();
  void assemble_system();
  void save_one(string);

  void solve();
  void compute_E_lin(LA::MPI::Vector&, double&, double&, double&);
  void estimate_error(double&);

  string m_guess_str;
};