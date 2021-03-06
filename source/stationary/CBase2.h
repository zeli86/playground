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


#include "MyParameterHandler.h"
#include <algorithm>
#include <cassert>
#include <limits>
#include <random>
#include <vector>
#include <array>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multiroots.h>

#ifndef __CBASE2_H__
#define __CBASE2_H__

#define STR1(x) #x
#define STR2(x) STR1(x)

enum Status { SUCCESS, FAILED, ZERO_SOL, SLOW_CONV };

double operator*( const vector<double>& rhs, const vector<double>& lhs )
{
  assert( rhs.size() == lhs.size() );
    
  double retval=0;
  const unsigned dim=rhs.size();
  for( unsigned i=0; i<dim; i++ )
    retval += rhs[i]*lhs[i];

return retval;
}

bool my_fun0( const vector<double>& rhs, const vector<double>& lhs )
{
  assert( rhs.size() == lhs.size() );
    
  double n1=0, n2=0;
  const unsigned dim=rhs.size();
  for( unsigned i=0; i<dim; i++ )
  {
    n1 += rhs[i]*rhs[i];
    n2 += lhs[i]*lhs[i];
  }    
  return (n1 < n2);
}

bool my_fun1( const vector<double>& rhs, const vector<double>& lhs )
{
  assert( rhs.size() == lhs.size() );
    
  double n1=0, n2=0;
  const unsigned dim=rhs.size();
  for( unsigned i=0; i<dim; i++ )
  {
    n1 += rhs[i]*rhs[i];
    n2 += lhs[i]*lhs[i];
  }  
  return (fabs(n1-n2)<1e-5);
}

class CBase2
{
  public:
    CBase2( const std::string );
    virtual ~CBase2() {};
    
    void dump_info_xml( const string="" );
    double l2norm_t() { return m_t[0]*m_t[0]+m_t[1]*m_t[1]; };
  
    virtual void compute_contributions()=0;
 
    MPI_Comm mpi_communicator;
    gsl_multiroot_function_fdf m_fun;

    void find_new_t();
  protected:
    double m_T[2];
    double m_W[5];
    double m_V2[9];
    double m_I12; 
    double m_x[2];
    double m_f[2];
    double m_f_df[4];

    void generate_initial_points();
    void screening();
    void select_t();
    void compute_f_df();
    
    double m_t[2];
    double m_t_old[2];
    double m_xmin, m_xmax;
    double m_ymin, m_ymax;
    double m_zmin, m_zmax;
    double m_res;
    double m_res_old;
    double m_resp;
    double m_res_over_resp;
    double m_ti;    
    double m_final_error;
    double m_mu;
    double m_dmu;
    double m_N;
    double m_L_halbe;

    vector<double> m_epsilon;
    vector<double> m_gs;
    vector<double> m_omega;
    vector<vector<double>> m_found_t;
    vector<vector<double>> m_t_guess;

    bool m_root;
    int m_rank;
    int m_max_iter;

    unsigned m_counter;
    unsigned m_global_refinement;
    unsigned m_total_no_cells;
    unsigned m_total_no_active_cells;    
    unsigned m_NA;
    unsigned m_Ndmu;
    unsigned m_no_initial_points;
    unsigned m_QN1[3];
    //unsigned m_QN2[3];
    
    ofstream m_computing_timer_log;
    TimerOutput m_computing_timer;    
    MyParameterHandler m_ph;
    ConditionalOStream pcout;

    string m_guess_str;
};

CBase2::CBase2( const std::string xmlfilename  ) 
  : 
  mpi_communicator(MPI_COMM_WORLD), 
  m_computing_timer_log("benchmark.txt"),
  m_computing_timer(mpi_communicator, m_computing_timer_log, TimerOutput::summary, TimerOutput:: cpu_and_wall_times ), 
  m_ph(xmlfilename),
  pcout (cout, (Utilities::MPI::this_mpi_process(mpi_communicator) == 0))
{
  m_root = (Utilities::MPI::this_mpi_process(mpi_communicator) == 0);
  MPI_Comm_rank(mpi_communicator, &m_rank);
 
  try 
  {
    m_omega = m_ph.Get_Physics("omega");
    m_gs = m_ph.Get_Physics("gs_1");
    m_QN1[0] = unsigned(m_ph.Get_Physics("QN1",0));
    m_QN1[1] = unsigned(m_ph.Get_Physics("QN1",1));
    m_QN1[2] = unsigned(m_ph.Get_Physics("QN1",2));

    m_xmin = m_ph.Get_Mesh("xrange",0);
    m_xmax = m_ph.Get_Mesh("xrange",1);
    m_ymin = m_ph.Get_Mesh("yrange",0);
    m_ymax = m_ph.Get_Mesh("yrange",1);
    m_zmin = m_ph.Get_Mesh("zrange",0);
    m_zmax = m_ph.Get_Mesh("zrange",1);
    m_global_refinement = unsigned(m_ph.Get_Mesh("global_refinements",0));

    m_ti = m_ph.Get_Algorithm("ti",0);
    m_t[0] = m_ti;
    m_t[1] = m_ti;
    m_NA = int(m_ph.Get_Algorithm("NA",0));
    m_Ndmu = m_ph.Get_Algorithm("Ndmu",0); 
    m_dmu = m_ph.Get_Algorithm("dmu",0);
    m_epsilon = m_ph.Get_Algorithm("epsilon");
  }
  catch( const std::string info )
  {
    std::cerr << info << endl;
    MPI_Abort( mpi_communicator, 0 );
  }

  m_counter=0;
  m_max_iter=200;
  m_final_error=0;
  m_L_halbe = 20;
  m_no_initial_points = 100;
  generate_initial_points();
  memset( m_V2,0,9*sizeof(double) );
}

void CBase2::generate_initial_points()
{
  for( auto& el : m_t_guess )
    el.clear();

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(-m_L_halbe, m_L_halbe);
  
  vector<double> tmp;
  for( unsigned n=0; n<m_no_initial_points; n++ ) 
  {
    tmp.push_back(dis(gen));
    tmp.push_back(dis(gen));
    m_t_guess.push_back(tmp);
    tmp.clear();
  }
}

void CBase2::screening()
{
  m_found_t.clear();
  compute_contributions();

  int iter;
  vector<double> vec(2);

  for( unsigned n=0; n<m_no_initial_points; n++ )
  {
    m_x[0] = m_t_guess[n][0];
    m_x[1] = m_t_guess[n][1];

    double f_norm, x_norm;
    iter=0;
    do    
    {
      compute_f_df();

      double fak = 1.0/(m_f_df[0]*m_f_df[3]-m_f_df[2]*m_f_df[1]);

      double a2 = fak*m_f_df[3];
      double b2 = -fak*m_f_df[2];
      double c2 = -fak*m_f_df[1];
      double d2 = fak*m_f_df[0];

      m_x[0] -= (a2*m_f[0]+b2*m_f[1]);
      m_x[1] -= (c2*m_f[0]+d2*m_f[1]);      

      f_norm = m_f[0]*m_f[0]+m_f[1]*m_f[1];
      x_norm = m_x[0]*m_x[0]+m_x[1]*m_x[1];
      iter++;
      //pcout << n << ", " << iter << " , x_norm = " << x_norm << ", f_norm = " << f_norm << endl;
    }
    while ( iter < m_max_iter && f_norm > 1e-10 );

    //pcout << iter << " , x_norm = " << x_norm << ", f_norm = " << f_norm << endl;

    if( iter < m_max_iter && f_norm < 1e-10 && x_norm > 0.1 )
    {
      vec[0] = m_x[0];
      vec[1] = m_x[1];
      m_found_t.push_back(vec);
    }
  }
  
 for ( vector<vector<double>>::iterator it=m_found_t.begin(); it!=m_found_t.end(); /*it++*/ ) 
 {
   if( fabs((*it)[0]) < 1e-5 || fabs((*it)[1]) < 1e-5 ) 
     it = m_found_t.erase(it);
   else 
     ++it;
  }

  std::sort( m_found_t.begin(), m_found_t.end(), my_fun0 );
  m_found_t.erase( std::unique(m_found_t.begin(), m_found_t.end(), my_fun1), m_found_t.end() );  
}

void CBase2::select_t ()
{
  int sel=-1;
  double old_norm = m_t_old[0]*m_t_old[0]+m_t_old[1]*m_t_old[1]; 
  double diff = std::numeric_limits<double>::max();
  for( int i=0; i<m_found_t.size(); i++ )
  {
    //double tmp = fabs( m_found_t[i][0]*m_found_t[i][0] - old_norm);
    double tmp = fabs( m_found_t[i][0]*m_found_t[i][0] + m_found_t[i][1]*m_found_t[i][1] - old_norm);
    if( tmp < diff )
    {
      sel=i;
      diff=tmp;
    }
    //pcout << "i = " << i << ", " <<  m_found_t[i][0] << ", " << m_found_t[i][1] << endl;
  }
  m_t[0] = m_found_t[sel][0];
  m_t[1] = m_found_t[sel][1];
  
  //pcout << "sel = " << sel << endl;
  //pcout << "selected " <<  m_t[0] << ", " <<  m_t[1] << ", " << fabs( m_t[0]*m_t[0] + m_t[1]*m_t[1] - old_norm) << endl;
}

void CBase2::find_new_t ()
{
  m_t_old[0] = m_t[0];
  m_t_old[1] = m_t[1];

  screening();
/*  
  for( auto el : m_found_t )
  {
    pcout << el*el;
    for( auto el2 : el )
      pcout << "\t" << el2;
    pcout << endl;
  }
*/  
  if( m_root ) select_t();
  MPI_Bcast( m_t, 2, MPI_DOUBLE, 0, mpi_communicator);
}

void CBase2::compute_f_df()
{
  memset(m_f,0,2*sizeof(double));
  memset(m_f_df,0,4*sizeof(double));

  const double t0 = m_x[0];
  const double t1 = m_x[1];
  const double t0q = t0*t0;
  const double t1q = t1*t1;
  const double t0k = t0q*t0;
  const double t1k = t1q*t1;

  m_f[0] = t0*m_T[0] + t1*m_I12;
  m_f[0] += t0k*m_W[0] + 3*t0*t1q*m_W[2] + 3*t0q*t1*m_W[3] + t1k*m_W[4];
  m_f[0] += t0k*m_V2[0] + t0q*t1*(m_V2[1]+m_V2[6]) + t0*t1q*(m_V2[2]+m_V2[7]) + t1k*m_V2[8];
  m_f[1] = t1*m_T[1] + t0*m_I12; 
  m_f[1] += t1k*m_W[1] + 3*t0*t1q*m_W[4] + 3*t0q*t1*m_W[2] + t0k*m_W[3];
  m_f[1] += t0k*m_V2[6] + t0q*t1*(m_V2[7]+m_V2[3]) + t0*t1q*(m_V2[8]+m_V2[4]) + t1k*m_V2[5];

  // df0/dt0
  m_f_df[0] = m_T[0] + 3*t0q*m_W[0] + 3*t1q*m_W[2] + 6*t0*t1*m_W[3];
  m_f_df[0] += 3*t0q*m_V2[0] + 2*t0*t1*(m_V2[1]+m_V2[6]) + t1q*(m_V2[2]+m_V2[7]);

  // df0/dt1
  m_f_df[2] = m_I12 + 6*t0*t1*m_W[2] + 3*t0q*m_W[3] + 3*t1q*m_W[4];
  m_f_df[2] += t0q*(m_V2[1]+m_V2[6]) + 2*t0*t1*(m_V2[2]+m_V2[7]) + 3*t1q*m_V2[8];

  // df1/dt0
  m_f_df[1] = m_I12 + 3*t1q*m_W[4] + 6*t0*t1*m_W[2] + 3*t0q*m_W[3];
  m_f_df[1] += 3*t0q*m_V2[6] + 2*t0*t1*(m_V2[7]+m_V2[3]) + t1q*(m_V2[8]+m_V2[4]);

  // df1/dt1
  m_f_df[3] = m_T[1] + 3*t1q*m_W[1] + 6*t0*t1*m_W[4] + 3*t0q*m_W[2];
  m_f_df[3] += t0q*(m_V2[7]+m_V2[3]) + 2*t0*t1*(m_V2[8]+m_V2[4]) + 3*t1q*m_V2[5];
}

void CBase2::dump_info_xml ( const string path )
{
  string filename = path + "info.xml";

  wofstream fs2;
  fs2.open(filename);

  locale utf8_locale("en_US.UTF8");
  fs2.imbue(utf8_locale); 
  fs2 << L"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<INFO>\n";
  fs2 << L"<MU>" << m_mu << L"</MU>\n";
  fs2 << L"<GS>" << m_gs[0] << L"</GS>\n";
  fs2 << L"<N>" << m_N << "L</N>\n";
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
