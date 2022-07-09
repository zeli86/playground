

// #define STR1(x) #x
// #define STR2(x) STR1(x)

// #include "ref_pt_list.h"
// #include "nlopt.h"
// #include "pugixml.hpp"

// template <int dim>
// double myfunc(unsigned n, const double* t, double* grad, void* my_func_data)
// {
//   MySolver<dim>* sol = reinterpret_cast<MySolver<dim>*>(my_func_data);
//   if (grad)
//   {
//     grad[0] = t[0] * sol->m_I[5] + t[1] * sol->m_I[7] + t[0] * t[0] * t[0] * sol->m_I[0] + 3.0 * t[0] * t[1] * t[1] * sol->m_I[2] + 3.0 * t[0] * t[0] * t[1] * sol->m_I[3] + t[1] * t[1] * t[1] * sol->m_I[4];
//     grad[1] = t[1] * sol->m_I[6] + t[0] * sol->m_I[7] + t[1] * t[1] * t[1] * sol->m_I[1] + 3.0 * t[0] * t[1] * t[1] * sol->m_I[4] + 3.0 * t[0] * t[0] * t[1] * sol->m_I[2] + t[0] * t[0] * t[0] * sol->m_I[3];
//   }
//   return (0.25 * sol->m_I[0] * pow(t[0], 4) + 0.5 * sol->m_I[5] * pow(t[0], 2) + 0.25 * sol->m_I[1] * pow(t[1], 4) + t[0] * sol->m_I[4] * pow(t[1], 3) + 1.5 * sol->m_I[2] * pow(t[0], 2) * pow(t[1], 2) + 0.5 * sol->m_I[6] * pow(t[1], 2) + pow(t[0], 3) * sol->m_I[3] * t[1] + t[0] * sol->m_I[7] * t[1]);
// }

template <int dim, int N>
class CBase
{
public:
  explicit CBase(const std::string&);
  virtual ~CBase() {};

  void find_ortho_min();
  void dump_info_xml(const string = "");

  double l2norm_t();

  virtual void compute_contributions() = 0;
protected:
  void screening();

  double m_t[N];
  double m_t_guess[N];

  double m_res;
  double m_res_old;
  double m_resp;
  double m_ti;
  double m_final_error;
  double m_rMu;
  double m_dmu;
  double m_gs;
  double m_N;
  vector<double> m_epsilon;
  vector<double> m_omega;

  unsigned m_counter;
  unsigned m_global_refinement;
  unsigned m_total_no_cells;
  unsigned m_total_no_active_cells;
  unsigned m_NA;
  unsigned m_Ndmu;

  unsigned m_QN1[3];

  MyParameterHandler m_ph;
};

template <int dim, int N>
CBase<dim, N>::CBase(const std::string& xmlfilename)
  :
  m_ph(xmlfilename)
{
  m_QN1[0] = unsigned(m_ph.Get_Physics("QN1", 0));
  m_omega = m_ph.Get_Physics("omega");
  m_gs = m_ph.Get_Physics("gs_1", 0);

  m_global_refinement = unsigned(m_ph.Get_Mesh("global_refinements", 0));

  m_ti = m_ph.Get_Algorithm("ti", 0);
  m_epsilon = m_ph.Get_Algorithm("epsilon");
  m_NA = int(m_ph.Get_Algorithm("NA", 0));
  m_Ndmu = m_ph.Get_Algorithm("Ndmu", 0);
  m_dmu = m_ph.Get_Algorithm("dmu", 0);

  m_counter = 0;
}

template <int dim, int N>
void CBase<dim, N>::screening()
{
  m_ref_pt_list_tmp.reset(10, 20);

  for (auto& it : m_ref_pt_list_tmp.m_list)
  {
    nlopt_opt opt;
    opt = nlopt_create(NLOPT_LD_MMA, N);
    nlopt_set_xtol_rel(opt, 1e-10);
    nlopt_set_min_objective(opt, myfunc<dim>, this);

    double* x = new double[N];
    double minf; /* the minimum objective value, upon return */

    /* some initial guess */
    for (int i = 0; i < N; i++)
    {
      x[i] = it.ti[i];
    }

    int status = nlopt_optimize(opt, x, &minf);
    /*
        if ( status < 0) {
            printf("nlopt failed!\n");
        }
        else {
            printf("found minimum at f(%g,%g) = %0.10g\n", x[0], x[1], minf);
        }
    */
    it.status = status;
    it.failed = (status < 0);
    if (!isfinite(minf))
    {
      continue;
    }
    it.f = minf;
    for (int i = 0; i < N; i++)
    {
      it.t[i] = x[i];
      //it.df[i] =
    }
    //it.DumpItem( std::cout );

    delete [] x;
    nlopt_destroy(opt);
  }
  m_ref_pt_list_tmp.condense();
}

#endif
