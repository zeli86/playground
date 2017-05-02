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


#ifndef __CBASE_NO_MPI_H__
#define __CBASE_NO_MPI_H__

#include "ref_pt_list.h"

enum Status { SUCCESS, FAILED, ZERO_SOL, SLOW_CONV, CONTINUE };

template <int dim> 
class CBase
{
  public:
    CBase( const std::string& );
    virtual ~CBase() {};
    
    void find_ortho_min( bool=true );

    const double l2norm_t() const;

    virtual void compute_contributions()=0;
  protected:
    void screening();

    double m_t[dim];
    double m_t_guess[dim];

    double m_xmin;
    double m_xmax;

    double m_res;
    double m_res_inf;
    double m_res_old;
    double m_resp;
    double m_res_over_resp;
    double m_ti;
    vector<double> m_epsilon;
    
    double m_mu;
    double m_dmu;
    double m_gs;
    vector<double> m_omega;

    double m_N;

    unsigned m_counter;
    unsigned m_global_refinement;
    unsigned m_total_no_cells;
    unsigned m_total_no_active_cells;    
    unsigned m_NA;
    unsigned m_Ndmu;
    
    unsigned m_QN1[3];
    
    ofstream m_computing_timer_log;
    TimerOutput m_computing_timer;    
    MyParameterHandler m_ph;

    MyUtils::ref_pt_list<dim> m_ref_pt_list;
    MyUtils::ref_pt_list<dim> m_ref_pt_list_tmp;
};

template <int dim>
CBase<dim>::CBase( const std::string& xmlfilename ) 
  : 
  m_computing_timer_log("benchmark.txt"),
  m_computing_timer( m_computing_timer_log, TimerOutput::summary, TimerOutput:: cpu_and_wall_times ), 
  m_ph(xmlfilename)
{
  m_QN1[0] = unsigned(m_ph.Get_Physics("QN1",0));
  m_omega = m_ph.Get_Physics("omega");
  m_gs = m_ph.Get_Physics("gs_1",0);
  
  m_global_refinement = unsigned(m_ph.Get_Mesh("global_refinements",0));
  m_xmin = m_ph.Get_Mesh("xMin",0);
  m_xmax = m_ph.Get_Mesh("xMax",1);

  m_ti = m_ph.Get_Algorithm("ti",0);
  m_epsilon = m_ph.Get_Algorithm("epsilon");
  m_t[0] = m_ti;
  m_t[1] = m_ti;
  m_t_guess[0] = m_ti;
  m_t_guess[1] = m_ti;

  m_NA = int(m_ph.Get_Algorithm("NA",0));
  m_Ndmu = m_ph.Get_Algorithm("Ndmu",0); 
  m_dmu = m_ph.Get_Algorithm("dmu",0);

  m_counter=0;
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
    size_t iter = 0;
    int status;
    double size;

    gsl_vector* x = gsl_vector_alloc(dim); 
    gsl_vector* ss = gsl_vector_alloc(dim); 
    gsl_vector_set_all (ss, 1.0);        
        
    gsl_multimin_function my_func;

    my_func.n = dim;
    my_func.f = fun<DIMENSION>;
    my_func.params = this;

    for( int i=0; i<dim; i++ )
      gsl_vector_set (x, i, it.ti[i]); 

    const gsl_multimin_fminimizer_type* T = gsl_multimin_fminimizer_nmsimplex2;
    gsl_multimin_fminimizer* s = gsl_multimin_fminimizer_alloc (T, dim);

    gsl_multimin_fminimizer_set (s, &my_func, x, ss );
    do
    {
      iter++;
      status = gsl_multimin_fminimizer_iterate (s);
      if (status) break;

      size = gsl_multimin_fminimizer_size (s);
      status = gsl_multimin_test_size (size, 1e-10);
    }
    while( status == GSL_CONTINUE && iter < 1000 );
    //printf ("status = %s\n", gsl_strerror (status));

    it.f = s->fval;
    for( int i=0; i<dim; i++ )
    {
      it.t[i] = gsl_vector_get (s->x, i);
    }
    //it.DumpItem( std::cout );

    gsl_vector_free (x);
    gsl_vector_free (ss);
    gsl_multimin_fminimizer_free (s);
  }
  m_ref_pt_list_tmp.remove_zero_ref_points();
}

template <int dim>
void CBase<dim>::find_ortho_min( bool bdel_f_non_zero )
{
  m_computing_timer.enter_section(__func__);
  compute_contributions();

  double l2_norm_t_old = 0;
  double min = std::numeric_limits<double>::max();
  for( int i=0; i<dim; i++ )
    l2_norm_t_old += m_t[i]*m_t[i];

  l2_norm_t_old = sqrt(l2_norm_t_old);

  CBase<dim>::screening();

  m_ref_pt_list_tmp.remove_zero_ref_points();
  m_ref_pt_list_tmp.remove_duplicates();
  if(bdel_f_non_zero) m_ref_pt_list_tmp.remove_f_non_zero();
  
  for( auto it : m_ref_pt_list_tmp.m_list )
  {
    if( fabs(it.l2norm_t() - l2_norm_t_old) < min )
    {
      min = fabs(it.l2norm_t() - l2_norm_t_old);
      for( int i=0; i<dim; i++ )
          m_t[i] = it.t[i];
    }
  }
  m_computing_timer.exit_section();  
}

#endif