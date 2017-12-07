/* * atus-pro testing - atus-pro testing playgroung
 * Copyright (C) 2017 Želimir Marojević <zelimir.marojevic@gmail.com>
 *
 * This file is part of atus-pro testing.
 *
 * atus-pro testing is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * atus-pro testing is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with atus-pro testing.  If not, see <http://www.gnu.org/licenses/>.
 */


#ifndef __CBASE_H__
#define __CBASE_H__

#define STR1(x) #x
#define STR2(x) STR1(x)

#include "ref_pt_list.h"
#include "nlopt.h"

template <int dim>
class MySolver;

template <int dim>
double myfunc(unsigned n, const double * t, double * grad, void * my_func_data)
{
  MySolver<dim> * sol = reinterpret_cast<MySolver<dim>*>(my_func_data); 
  if (grad) 
  {
    grad[0] = t[0]*sol->m_I[5] + t[1]*sol->m_I[7] + t[0]*t[0]*t[0]*sol->m_I[0] + 3.0*t[0]*t[1]*t[1]*sol->m_I[2] + 3.0*t[0]*t[0]*t[1]*sol->m_I[3] + t[1]*t[1]*t[1]*sol->m_I[4];
    grad[1] = t[1]*sol->m_I[6] + t[0]*sol->m_I[7] + t[1]*t[1]*t[1]*sol->m_I[1] + 3.0*t[0]*t[1]*t[1]*sol->m_I[4] + 3.0*t[0]*t[0]*t[1]*sol->m_I[2] + t[0]*t[0]*t[0]*sol->m_I[3];
  }
  return (0.25*sol->m_I[0]*pow(t[0],4) + 0.5*sol->m_I[5]*pow(t[0],2) + 0.25*sol->m_I[1]*pow(t[1],4) + t[0]*sol->m_I[4]*pow(t[1],3) +1.5*sol->m_I[2]*pow(t[0],2)*pow(t[1],2) +0.5*sol->m_I[6]*pow(t[1],2) +pow(t[0],3)*sol->m_I[3]*t[1] +t[0]*sol->m_I[7]*t[1]);
}

enum Status { SUCCESS, FAILED, ZERO_SOL, SLOW_CONV, MAXITER, SINGULAR };

template <int dim> 
class CBase
{
  public:
    CBase( const std::string );
    virtual ~CBase() {};
    
    int find_ortho_min( bool=true );
    void dump_info_xml( const string="" );

    const double l2norm_t() const;
    
    virtual void compute_contributions()=0;
 
    MPI_Comm mpi_communicator;
  protected:
    void screening();

    double m_t[dim];
    double m_t_guess[dim];

    double m_xmin, m_xmax;
    double m_ymin, m_ymax;
    double m_zmin, m_zmax;

    double m_res[2];
    double m_res_inf[2];
    double m_res_old[2];
    double m_resp[2];
    double m_res_over_resp[2];
    double m_ti;    
    double m_final_error;
    
    double m_mu[2];
    double m_dmu;
    vector<double> m_epsilon;
    vector<double> m_gs;
    vector<double> m_gs2;
    vector<double> m_omega;

    double m_N[2];

    bool m_root;
    int m_rank;

    unsigned m_counter;
    unsigned m_maxiter;
    unsigned m_global_refinement;
    unsigned m_total_no_cells;
    unsigned m_total_no_active_cells;    
    unsigned m_NA;
    unsigned m_Ndmu;
    
    unsigned m_QN1[3];
    unsigned m_QN2[3];
    unsigned m_QN3[3];
    unsigned m_QN4[3];
    
    ofstream m_computing_timer_log;
    TimerOutput m_computing_timer;    
    MyParameterHandler m_ph;
    ConditionalOStream pcout;

    MyUtils::ref_pt_list<dim> m_ref_pt_list;
    MyUtils::ref_pt_list<dim> m_ref_pt_list_tmp;
    string m_guess_str;
};

template <int dim>
CBase<dim>::CBase( const std::string xmlfilename  ) 
  : 
  mpi_communicator(MPI_COMM_WORLD), 
  m_computing_timer_log("benchmark.txt"),
  m_computing_timer(mpi_communicator, m_computing_timer_log, TimerOutput::summary, TimerOutput:: cpu_and_wall_times ), 
  m_ph(xmlfilename),
  pcout (cout, (Utilities::MPI::this_mpi_process(mpi_communicator) == 0))
{
  try 
  {
    m_omega = m_ph.Get_Physics("omega");
    m_gs = m_ph.Get_Physics("gs_1");
    m_gs2 = m_ph.Get_Physics("gs_1");
    m_QN1[0] = int(m_ph.Get_Physics("QN1",0));
    m_QN1[1] = int(m_ph.Get_Physics("QN1",1));
    m_QN1[2] = int(m_ph.Get_Physics("QN1",2));

    m_xmin = m_ph.Get_Mesh("xrange",0);
    m_xmax = m_ph.Get_Mesh("xrange",1);
    m_ymin = m_ph.Get_Mesh("yrange",0);
    m_ymax = m_ph.Get_Mesh("yrange",1);
    m_zmin = m_ph.Get_Mesh("zrange",0);
    m_zmax = m_ph.Get_Mesh("zrange",1);
    m_global_refinement = unsigned(m_ph.Get_Mesh("global_refinements",0));

    m_ti = m_ph.Get_Algorithm("ti",0);
    m_epsilon = m_ph.Get_Algorithm("epsilon");
    m_t[0] = m_ti;
    m_t[1] = m_ti;
    m_t_guess[0] = m_ti;
    m_t_guess[1] = m_ti;

    m_NA = int(m_ph.Get_Algorithm("NA",0));
    m_Ndmu = m_ph.Get_Algorithm("Ndmu",0); 
    m_dmu = m_ph.Get_Algorithm("dmu",0);
    m_maxiter = 1000;
  }
  catch( const std::string info )
  {
    std::cerr << info << endl;
    MPI_Abort( mpi_communicator, 0 );
  }

  /*
  m_mu[0] = m_ph.get_double("mu");
  m_mu[1] = m_ph.get_double("mu2");
  m_guess_str = m_ph.get( "guess_fct" );
  */
  
  m_root = (Utilities::MPI::this_mpi_process(mpi_communicator) == 0);
  MPI_Comm_rank(mpi_communicator, &m_rank);
  
  m_counter=0;
  m_final_error=0;
}

template <int dim>
const double CBase<dim>::l2norm_t() const
{
  double retval=0;
  for( int i=0; i<dim; i++ )
    retval += m_t[i]*m_t[i];
return sqrt(retval);
}

template <int dim>
void CBase<dim>::screening()
{
  m_ref_pt_list_tmp.reset( 5, 20 );

  for( auto& it : m_ref_pt_list_tmp.m_list )
  {
    nlopt_opt opt;
    opt = nlopt_create(NLOPT_LD_MMA, dim);
    nlopt_set_xtol_rel(opt, 1e-10);
    nlopt_set_min_objective(opt, myfunc<dim>, this);

    double * x = new double[dim];  
    double minf; /* the minimum objective value, upon return */
   
    /* some initial guess */
    for( int i=0; i<dim; i++ )
      x[i] = it.ti[i]; 
    
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
    if( !isfinite(minf) ) continue;
    it.f = minf;
    for( int i=0; i<dim; i++ )
    {
      it.t[i] = x[i];
      //it.df[i] = 
    }
    //it.DumpItem( std::cout );

    delete [] x;
    nlopt_destroy(opt);
  }

  m_ref_pt_list_tmp.remove_zero_ref_points();  
}

template <int dim>
int CBase<dim>::find_ortho_min( bool bdel_f_non_zero )
{
  compute_contributions();

  m_computing_timer.enter_section(__func__);

  if( m_root )
  {
    double l2_norm_t_old = 0;
    double min = std::numeric_limits<double>::max();
    for( int i=0; i<dim; i++ )
      l2_norm_t_old += m_t[i]*m_t[i];

    l2_norm_t_old = sqrt(l2_norm_t_old);

    CBase<dim>::screening();

    m_ref_pt_list_tmp.remove_zero_ref_points();
    m_ref_pt_list_tmp.remove_duplicates();
    //if(bdel_f_non_zero) m_ref_pt_list_tmp.remove_f_non_zero();
    ///m_ref_pt_list_tmp.Dump( std::cout );
    //m_ref_pt_list_tmp.Dump( std::cout );

    
    for( auto it : m_ref_pt_list_tmp.m_list )
    {
      double tmp1 = fabs(it.l2norm_t() - l2_norm_t_old);
      if( tmp1 < min )
      {
        min = tmp1;
        for( int i=0; i<dim; i++ )
            m_t[i] = it.t[i];
      }
    }
  }
  int retval = m_ref_pt_list_tmp.m_list.empty();
  MPI_Bcast( m_t, dim, MPI_DOUBLE, 0, mpi_communicator );
  MPI_Bcast( &retval, dim, MPI_INT, 0, mpi_communicator );
  m_computing_timer.exit_section();
return retval;
}

template <int dim>
void CBase<dim>::dump_info_xml ( const string path )
{
  //TODO: change to pugixml
  string filename = path + "info.xml";

  wofstream fs2;
  fs2.open(filename);

  locale utf8_locale("en_US.UTF8");
  fs2.imbue(utf8_locale); 
  fs2 << L"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<INFO>\n";
  fs2 << L"<MU>" << m_mu[0] << L"</MU>\n";
  fs2 << L"<MU2>" << m_mu[1] << L"</MU2>\n";
  fs2 << L"<GS>" << m_gs[0] << L"</GS>\n";
  fs2 << L"<GS2>" << m_gs[1] << L"</GS2>\n";
  fs2 << L"<N>" << m_N[0] << "L</N>\n";
  fs2 << L"<N2>" << m_N[1] << "L</N2>\n";
  fs2 << L"<XMIN>" << m_xmin << L"</XMIN>\n";
  fs2 << L"<XMAX>" << m_xmax << L"</XMAX>\n";
  fs2 << L"<YMIN>" << m_ymin << L"</YMIN>\n";
  fs2 << L"<YMAX>" << m_ymax << L"</YMAX>\n";
  fs2 << L"<ZMIN>" << m_zmin << L"</ZMIN>\n";
  fs2 << L"<ZMAX>" << m_zmax << L"</ZMAX>\n";
  fs2 << L"<FINAL_ERROR>" << m_final_error << L"</FINAL_ERROR>\n";
  fs2 << L"<REVISION>" << STR2(GIT_SHA1) << L"</REVISION>\n";
  fs2 << L"</INFO>\n";
  fs2.close();
}

#endif