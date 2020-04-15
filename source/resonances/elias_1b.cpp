//
// atus-pro testing - atus-pro testing playgroung
// Copyright (C) 2017 Želimir Marojević <zelimir.marojevic@gmail.com>
//
// This file is part of atus-pro testing.
//
// atus-pro testing is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// atus-pro testing is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with atus-pro testing.  If not, see <http://www.gnu.org/licenses/>.
//

/** Želimir Marojević
 */

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/function_parser.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/precondition.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_out.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_q.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_sf_erf.h>

#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <locale>
#include <limits>
#include <cmath>
#include <array>
#include <unistd.h>

#include "global.h"
#include "MyParameterHandler.h"
#include "my_table.h"

#define STR1(x) #x
#define STR2(x) STR1(x)

namespace BreedSolver
{
  using namespace std;
  using namespace dealii;

  enum Status { SUCCESS, FAILED, ZERO_SOL, SLOW_CONV };

  class MySolver
  {
  public:
    MySolver( const std::string );
    virtual ~MySolver();

    void run2b ();
        
  protected:
    FunctionParser<1> m_Potential;
    FunctionParser<1> m_F;
    FunctionParser<1> m_fctp_re;
    FunctionParser<1> m_fctp_im;

    double m_xmin, m_xmax;
    double m_res;
    double m_res_old;
    double m_resp;
    double m_N;
    
    double m_mu; // real part of chemical potential
    double m_E; // energy
    double m_Gamma; // imaginary part of chemical potential (life time)
    double m_gs; // self interaction parameter
    double m_sigma; // mass of field
    double m_q; // GM/c^2
    double m_potmin; 
    double m_potmax; 
    double m_s; // computed shift of the initial guess
    double m_alpha; // computed witdh of the initial guess
    double m_bcval; // Psi(xmin)
    vector<double> m_epsilon;

    unsigned m_counter;
    unsigned m_NA;
    unsigned m_Ndmu;
    unsigned m_global_refinements;
        
    int DoIter( string="" );

    void make_grid_custom();
    void setup_system();
    void assemble_rhs();
    void assemble_system( Vector<double>&, Vector<double>& );
    void assemble_system_2( Vector<double>&, Vector<double>& );
    void save( string );
    void save_one( string );
    void compute_potmax();
    void compute_potmin();
    void compute_guess_parameter();
    void Project_gradient();
    void compute_mu_and_E();
    void dump_info_xml( string="" );
        
    double Particle_Number();
    
    void solve( Vector<double>&, Vector<double>& );
    void output_results ( string, string = "step" );
    void output_guess ();

    MyParameterHandler m_ph;
    Triangulation<1> triangulation;
    FE_Q<1> fe;
    SparsityPattern sparsity_pattern;
    DoFHandler<1> dof_handler;

    SparseMatrix<double> system_matrix;
    Vector<double> m_system_rhs_re;
    Vector<double> m_system_rhs_im;
    Vector<double> m_system_rhs_2_re;
    Vector<double> m_system_rhs_2_im;
    Vector<double> m_Psi_re;
    Vector<double> m_Psi_im;
    Vector<double> m_Psi_sob_re;
    Vector<double> m_Psi_sob_im;
    Vector<double> m_sob_grad_re;
    Vector<double> m_sob_grad_im;    
    Vector<double> m_workspace;
    Vector<double> m_error_per_cell;

    string m_Psi_1_str;
    string m_Psi_2_str;
    string m_potential_str;
    string m_F_str;
    MyTable m_table;
    MyTable m_results;
  };

/*************************************************************************************************/
/**
 * Constructor
 */
  MySolver::MySolver ( const string xmlfilename ) 
    : 
    m_ph(xmlfilename),
    triangulation (Triangulation<1>::MeshSmoothing(Triangulation<1>::limit_level_difference_at_vertices|Triangulation<1>::eliminate_refined_inner_islands|Triangulation<1>::smoothing_on_refinement|Triangulation<1>::smoothing_on_coarsening)),
    fe (gl_degree_fe),
    dof_handler (triangulation)
  {
    try
    {
      m_gs = m_ph.Get_Physics("gs_1",0);
      m_Gamma = m_ph.Get_Physics("Gamma",0); // has to negative
      m_sigma = m_ph.Get_Physics("sigma",0); 
      m_q = m_ph.Get_Physics("q",0); 
      m_bcval = m_ph.Get_Physics("bcval",0); 

      m_xmin = m_ph.Get_Mesh("xrange",0);
      m_xmax = m_ph.Get_Mesh("xrange",1);
      m_global_refinements = unsigned(m_ph.Get_Mesh("global_refinements",0));

      m_NA = int(m_ph.Get_Algorithm("NA",0));
      m_epsilon = m_ph.Get_Algorithm("epsilon");
      
      m_F_str = m_ph.Get_Parameter("F");
      m_potential_str = m_ph.Get_Parameter("Potential");
    }
    catch( const std::string info )
    {
      std::cerr << info << endl;
      exit(0);
    }    
    m_counter=0;
   
    map<string,double> constants;
    constants["PI"] = numbers::PI;
    constants["sigma"] = m_sigma;
    constants["q"] = m_q;
   
    m_Potential.initialize( "r", m_potential_str, constants );
    m_F.initialize( "r", m_F_str, constants );
  }

  MySolver::~MySolver ()
  {
    dof_handler.clear ();
  }
  
  struct helper
  {
    double sigma;
    double q;  
    double sign;
  };

  struct helper2
  {
    double x0;
    double w;  
  };
  
  double fn1 (double r, void * params)
  {
    helper * p = reinterpret_cast<helper*>(params);
    double sigma = p->sigma;
    double q = p->q;
    return p->sign*((r*r*r + 3*q - 3*r) * (3*r*r*r*sigma*sigma - 2*r*r*r + 3*q) / 9 / (r*r*r*r));
  }
  
  int fn2 ( const gsl_vector * x, void * params, gsl_vector * f )
  {
    helper2 * p = reinterpret_cast<helper2*>(params);
      
    const double x0 = p->x0;
    const double w = p->w;

    const double alp = gsl_vector_get (x, 0);
    const double alpq = alp*alp;
    const double s = gsl_vector_get (x, 1);

    const double a = 0.25*alpq*(x0+s)*exp(-2*(s-x0)*(s-x0)/alpq);
    const double b = 0.0625*alp*sqrt(2*M_PI)*(alpq+4*s*s)*(1+gsl_sf_erf(sqrt(2)*(s-x0)/alp));
    const double c = exp(-(x0-s)*(x0-s)/alpq);

    gsl_vector_set (f, 0, a+b-1);
    gsl_vector_set (f, 1, c-w);

    return GSL_SUCCESS;
  }  
  
  void MySolver::compute_potmax()
  {
    helper params;
    params.sigma = m_sigma;
    params.q = m_q;
    params.sign = 1;
    
    int status;
    int iter=0, max_iter=1000;
    
    double a=0.05, b=2.0, m=0.15;
    
    gsl_function F;
    F.function = &fn1;
    F.params = &params;

    const gsl_min_fminimizer_type * T = gsl_min_fminimizer_brent;
    gsl_min_fminimizer * s = gsl_min_fminimizer_alloc (T);
    
    gsl_min_fminimizer_set (s, &F, m, a, b);

    do
    {
      iter++;
      status = gsl_min_fminimizer_iterate (s);

      m = gsl_min_fminimizer_x_minimum (s);
      a = gsl_min_fminimizer_x_lower (s);
      b = gsl_min_fminimizer_x_upper (s);

      status = gsl_min_test_interval (a, b, 1e-8, 0.0);

      //if (status == GSL_SUCCESS) printf ("Converged:\n");
      //printf ("%5d %.7g %.7g\n", iter, m, b - a);
    }
    while (status == GSL_CONTINUE && iter < max_iter);

    m_potmax = m;
    
    gsl_min_fminimizer_free (s);
  }
  
  void MySolver::compute_potmin()
  {
    helper params;
    params.sigma = m_sigma;
    params.q = m_q;
    params.sign = -1;
    
    int status;
    int iter=0, max_iter=1000;
    
    double a=m_potmax, b=5.0, m=1;
    
    gsl_function F;
    F.function = &fn1;
    F.params = &params;
    
    //cout << "f(a)" << fn1(a,&params) << endl;
    //cout << "f(m)" << fn1(m,&params) << endl;
    //cout << "f(b)" << fn1(b,&params) << endl;

    const gsl_min_fminimizer_type * T = gsl_min_fminimizer_brent;
    gsl_min_fminimizer * s = gsl_min_fminimizer_alloc (T);
    
    gsl_min_fminimizer_set (s, &F, m, a, b);

    do
    {
      iter++;
      status = gsl_min_fminimizer_iterate (s);

      m = gsl_min_fminimizer_x_minimum (s);
      a = gsl_min_fminimizer_x_lower (s);
      b = gsl_min_fminimizer_x_upper (s);

      status = gsl_min_test_interval (a, b, 1e-8, 0.0);

      //if (status == GSL_SUCCESS) printf ("Converged:\n");
      //printf ("%5d %.7g %.7g\n", iter, m, b - a);
    }
    while (status == GSL_CONTINUE && iter < max_iter);

    m_potmin = m;

    gsl_min_fminimizer_free (s);
  } 
  
int
print_state (size_t iter, gsl_multiroot_fsolver * s)
{
  printf ("iter = %3u x = % .3f % .3f "
          "f(x) = % .3e % .3e\n",
          iter,
          gsl_vector_get (s->x, 0), 
          gsl_vector_get (s->x, 1),
          gsl_vector_get (s->f, 0), 
          gsl_vector_get (s->f, 1));
}  
  
  void MySolver::compute_guess_parameter()
  {
    const gsl_multiroot_fsolver_type *T;
    gsl_multiroot_fsolver *s;

    int status;
    size_t iter=0;

    struct helper2 p = {m_xmin, m_bcval};
    gsl_multiroot_function f = {&fn2, 2, &p};

    gsl_vector * x = gsl_vector_alloc (2);
    gsl_vector_set (x, 0, 1);
    gsl_vector_set (x, 1, m_potmin);

    T = gsl_multiroot_fsolver_hybrids;
    s = gsl_multiroot_fsolver_alloc (T, 2);
    gsl_multiroot_fsolver_set (s, &f, x);

    print_state (iter, s);
    do
    {
      iter++;
      status = gsl_multiroot_fsolver_iterate (s);

      print_state (iter, s);

      if (status) break;

      status = gsl_multiroot_test_residual (s->f, 1e-7);
    }
    while (status == GSL_CONTINUE && iter < 1000);

    m_alpha = gsl_vector_get (s->x, 0);
    m_s = gsl_vector_get (s->x, 1);
    
    cout << "m_alpha = " << m_alpha << endl;
    cout << "m_s = " << m_s << endl;

    printf ("status = %s\n", gsl_strerror (status));

    gsl_multiroot_fsolver_free (s);
    gsl_vector_free (x);
  }

  double MySolver::Particle_Number()
  {
    double retval=0;
   
    const QGauss<1>  quadrature_formula(fe.degree+1);
    FEValues<1> fe_values (fe, quadrature_formula, update_values|update_JxW_values|update_quadrature_points);

    const unsigned n_q_points = quadrature_formula.size();
    vector<double> vals_re(n_q_points);
    vector<double> vals_im(n_q_points);

    DoFHandler<1>::active_cell_iterator cell=dof_handler.begin_active(), endc=dof_handler.end();
    for( ; cell!=endc; ++cell )
    {
      fe_values.reinit (cell);
      fe_values.get_function_values( m_Psi_re, vals_re );      
      fe_values.get_function_values( m_Psi_im, vals_im );      

      for( unsigned qp=0; qp<n_q_points; ++qp )
        retval += fe_values.JxW(qp)*fe_values.quadrature_point(qp)[0]*fe_values.quadrature_point(qp)[0]*(vals_re[qp]*vals_re[qp]+vals_im[qp]*vals_im[qp]);
    }
  return retval;
  }
  
  void MySolver::Project_gradient()
  {
    const QGauss<1>  quadrature_formula(fe.degree+1);
    FEValues<1> fe_values (fe, quadrature_formula, update_values|update_JxW_values|update_quadrature_points);

    double tmp[4] = {};

    const unsigned n_q_points = quadrature_formula.size();
    vector<double> vals_Psi_re(n_q_points);
    vector<double> vals_Psi_im(n_q_points);
    vector<double> vals_Psi_sob_re(n_q_points);
    vector<double> vals_Psi_sob_im(n_q_points);
    vector<double> vals_sob_grad_re(n_q_points);
    vector<double> vals_sob_grad_im(n_q_points);
 
    DoFHandler<1>::active_cell_iterator cell=dof_handler.begin_active(), endc=dof_handler.end();
    for( ; cell!=endc; ++cell )
    {
      fe_values.reinit (cell);
      fe_values.get_function_values( m_Psi_re, vals_Psi_re );      
      fe_values.get_function_values( m_Psi_im, vals_Psi_im );      
      fe_values.get_function_values( m_Psi_sob_re, vals_Psi_sob_re );      
      fe_values.get_function_values( m_Psi_sob_im, vals_Psi_sob_im );      
      fe_values.get_function_values( m_sob_grad_re, vals_sob_grad_re );      
      fe_values.get_function_values( m_sob_grad_im, vals_sob_grad_im );

      for( unsigned qp=0; qp<n_q_points; ++qp )
      {
        double JxW = fe_values.JxW(qp)*fe_values.quadrature_point(qp)[0]*fe_values.quadrature_point(qp)[0];
        tmp[0] += JxW*(vals_Psi_re[qp]*vals_sob_grad_re[qp] + vals_Psi_im[qp]*vals_sob_grad_im[qp]);
        tmp[1] += JxW*(vals_Psi_re[qp]*vals_sob_grad_im[qp] - vals_Psi_im[qp]*vals_sob_grad_re[qp]);
        tmp[2] += JxW*(vals_Psi_re[qp]*vals_Psi_sob_re[qp] + vals_Psi_im[qp]*vals_Psi_sob_im[qp]);
        tmp[3] += JxW*(vals_Psi_re[qp]*vals_Psi_sob_im[qp] - vals_Psi_im[qp]*vals_Psi_sob_re[qp]);
      }
    }
    
    const double fak = (tmp[0]*tmp[2]+tmp[1]*tmp[3])/(tmp[2]*tmp[2]+tmp[3]*tmp[3]);
    
    m_sob_grad_re.add( -fak, m_Psi_sob_re );
    m_sob_grad_im.add( -fak, m_Psi_sob_im );
  }
  
  void MySolver::assemble_rhs ()
  {
    m_system_rhs_re=0;
    m_system_rhs_im=0;
    m_system_rhs_2_re=0;
    m_system_rhs_2_im=0;
    
    const QGauss<1> quadrature_formula(fe.degree+1);

    FEValues<1> fe_values (fe, quadrature_formula, update_values|update_gradients|update_JxW_values|update_quadrature_points);

    const unsigned dofs_per_cell = fe.dofs_per_cell;
    const unsigned n_q_points = quadrature_formula.size();

    vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    
    vector<double> vals_re(n_q_points);
    vector<double> vals_im(n_q_points);
    vector<Tensor<1,1>> grads_re(n_q_points);
    vector<Tensor<1,1>> grads_im(n_q_points);

    double JxW, Pot, req, imq, tmp, F;
    DoFHandler<1>::active_cell_iterator cell=dof_handler.begin_active(), endc=dof_handler.end();
    for( ; cell!=endc; ++cell )
    {
      fe_values.reinit (cell);
      fe_values.get_function_values(m_Psi_re, vals_re);
      fe_values.get_function_gradients(m_Psi_re, grads_re);
      fe_values.get_function_values(m_Psi_im, vals_im);
      fe_values.get_function_gradients(m_Psi_im, grads_im);
     
      cell->get_dof_indices (local_dof_indices);
        
      for( unsigned qp=0; qp<n_q_points; ++qp )
      {
        JxW = fe_values.JxW(qp)*fe_values.quadrature_point(qp)[0]*fe_values.quadrature_point(qp)[0];
        Pot = m_Potential.value(fe_values.quadrature_point(qp));
        req = vals_re[qp]*vals_re[qp];
        imq = vals_im[qp]*vals_im[qp];
        F = m_F.value(fe_values.quadrature_point(qp));

        for( unsigned i=0; i<dofs_per_cell; ++i )
        {
          m_system_rhs_re(local_dof_indices[i]) += JxW*(grads_re[qp]*fe_values.shape_grad(i,qp) + (Pot+m_gs*F*(req-3*imq))*vals_re[qp]*fe_values.shape_value(i,qp) + m_Gamma*vals_im[qp]*fe_values.shape_value(i,qp));
          m_system_rhs_im(local_dof_indices[i]) += JxW*(grads_im[qp]*fe_values.shape_grad(i,qp) + (Pot+m_gs*F*(3*req-imq))*vals_im[qp]*fe_values.shape_value(i,qp) - m_Gamma*vals_re[qp]*fe_values.shape_value(i,qp));
          m_system_rhs_2_re(local_dof_indices[i]) += JxW*vals_re[qp]*fe_values.shape_value(i,qp);
          m_system_rhs_2_im(local_dof_indices[i]) += JxW*vals_im[qp]*fe_values.shape_value(i,qp);
        }
      }
    }
  }

  void MySolver::assemble_system ( Vector<double>& rhs, Vector<double>& sol )
  {
    system_matrix=0;
    const QGauss<1> quadrature_formula(fe.degree+1);

    FEValues<1> fe_values (fe, quadrature_formula, update_values|update_gradients|update_JxW_values|update_quadrature_points);

    const unsigned dofs_per_cell = fe.dofs_per_cell;
    const unsigned n_q_points = quadrature_formula.size();

    FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);

    vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    double JxW, Pot, fak;
    DoFHandler<1>::active_cell_iterator cell=dof_handler.begin_active(), endc=dof_handler.end();
    for( ; cell!=endc; ++cell )
    {
      cell_matrix=0;

      fe_values.reinit (cell);
      cell->get_dof_indices (local_dof_indices);

      for( unsigned qp=0; qp<n_q_points; ++qp )
      {
        JxW = fe_values.JxW(qp)*fe_values.quadrature_point(qp)[0]*fe_values.quadrature_point(qp)[0];
        Pot = m_Potential.value(fe_values.quadrature_point(qp));
        

        for( unsigned i=0; i<dofs_per_cell; ++i )
          for( unsigned j=0; j<dofs_per_cell; ++j )
            cell_matrix(i,j) += JxW*( (1+Pot)*fe_values.shape_value(i,qp)*fe_values.shape_value(j,qp) + fe_values.shape_grad(i,qp)*fe_values.shape_grad(j,qp));
      }
      
      for( unsigned i=0; i<dofs_per_cell; ++i )
      {
        for( unsigned j=0; j<dofs_per_cell; ++j ) 
          system_matrix.add (local_dof_indices[i],local_dof_indices[j], cell_matrix(i,j));
      }      
    }
    
    map<types::global_dof_index,double> boundary_values;
    VectorTools::interpolate_boundary_values (dof_handler, 0, ZeroFunction<1>(), boundary_values);
    MatrixTools::apply_boundary_values (boundary_values, system_matrix, sol, rhs);
  }
  
  void MySolver::assemble_system_2 ( Vector<double>& rhs, Vector<double>& sol )
  {
    system_matrix=0;
    const QGauss<1> quadrature_formula(fe.degree+1);

    FEValues<1> fe_values (fe, quadrature_formula, update_values|update_JxW_values|update_quadrature_points);

    const unsigned dofs_per_cell = fe.dofs_per_cell;
    const unsigned n_q_points = quadrature_formula.size();

    FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);

    vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    double JxW, Pot, fak;
    DoFHandler<1>::active_cell_iterator cell=dof_handler.begin_active(), endc=dof_handler.end();
    for( ; cell!=endc; ++cell )
    {
      cell_matrix=0;

      fe_values.reinit (cell);
      cell->get_dof_indices (local_dof_indices);

      for( unsigned qp=0; qp<n_q_points; ++qp )
      {
        JxW = fe_values.JxW(qp)*fe_values.quadrature_point(qp)[0]*fe_values.quadrature_point(qp)[0];

        for( unsigned i=0; i<dofs_per_cell; ++i )
          for( unsigned j=0; j<dofs_per_cell; ++j )
            cell_matrix(i,j) += JxW*fe_values.shape_value(i,qp)*fe_values.shape_value(j,qp);
      }
      
      for( unsigned i=0; i<dofs_per_cell; ++i )
      {
        for( unsigned j=0; j<dofs_per_cell; ++j ) 
          system_matrix.add (local_dof_indices[i],local_dof_indices[j], cell_matrix(i,j));
      }
    }
    
    map<types::global_dof_index,double> boundary_values;
    VectorTools::interpolate_boundary_values (dof_handler, 0, ZeroFunction<1>(), boundary_values);
    MatrixTools::apply_boundary_values (boundary_values, system_matrix, sol, rhs);
  }  
  
  void MySolver::compute_mu_and_E()
  {
    m_mu=0;
    m_E=0;
    
    const QGauss<1>  quadrature_formula(fe.degree+1);
    FEValues<1> fe_values (fe, quadrature_formula, update_gradients|update_values|update_JxW_values|update_quadrature_points);

    const unsigned n_q_points = quadrature_formula.size();
    
    vector<double> vals_re(n_q_points);
    vector<double> vals_im(n_q_points);
    vector<Tensor<1,1>> grads_re(n_q_points);
    vector<Tensor<1,1>> grads_im(n_q_points);

    double z1[] = {0,0};
    double z2[] = {0,0};
    double JxW, req, imq, F, Pot;
    DoFHandler<1>::active_cell_iterator cell=dof_handler.begin_active(), endc=dof_handler.end();
    for( ; cell!=endc; ++cell )
    {
      fe_values.reinit (cell);
      
      fe_values.get_function_values(m_Psi_re, vals_re);
      fe_values.get_function_gradients(m_Psi_re, grads_re);
      fe_values.get_function_values(m_Psi_im, vals_im);
      fe_values.get_function_gradients(m_Psi_im, grads_im);

      for( unsigned qp=0; qp<n_q_points; ++qp )
      {
        JxW = fe_values.JxW(qp)*fe_values.quadrature_point(qp)[0]*fe_values.quadrature_point(qp)[0];
        req = vals_re[qp]*vals_re[qp];
        imq = vals_im[qp]*vals_im[qp];
        F = m_F.value(fe_values.quadrature_point(qp));
        Pot = m_Potential.value(fe_values.quadrature_point(qp));
        
        z2[0] += JxW*(req-imq);
        z2[1] += 2*JxW*vals_re[qp]*vals_im[qp];
        z1[0] += JxW*( grads_re[qp]*grads_re[qp] - grads_im[qp]*grads_im[qp] + Pot*(req-imq) + F*m_gs*(req*req+imq*imq-6*req*imq) ); 
        z1[1] += JxW*( 2*grads_re[qp]*grads_im[qp] + 2*Pot*vals_re[qp]*vals_im[qp] + 4*F*m_gs*(req*vals_re[qp]*vals_im[qp]-vals_re[qp]*vals_im[qp]*imq)); 
        m_E += JxW*( grads_re[qp]*grads_re[qp] + grads_im[qp]*grads_im[qp] + (Pot+0.5*m_gs*F*(req+imq))*(req+imq) );
      }
    }
    
    double z12[]={0,0};
    double fak = 1/(z2[0]*z2[0]+z2[1]*z2[1]);
    z12[0] = fak*(z1[0]*z2[0]+z1[1]*z2[1]);
    z12[1] = fak*(z1[1]*z2[0]-z1[0]*z2[1]);

    m_mu = z12[0];
  }
 
  void MySolver::solve ( Vector<double>& rhs, Vector<double>& sol )
  {
    SparseDirectUMFPACK A_direct;
    A_direct.initialize(system_matrix);
    A_direct.vmult (sol, rhs);
  }
  
  void MySolver::make_grid_custom ()
  {
    cout << "m_xmin = " << m_xmin << endl;
    Point<1,double> pt1(m_xmin);
    Point<1,double> pt2(m_xmax);
    
    GridGenerator::hyper_rectangle(triangulation, pt2, pt1);
    triangulation.refine_global(m_global_refinements);
  }
  
  void MySolver::setup_system()
  {
    dof_handler.distribute_dofs (fe);

    m_system_rhs_re.reinit(dof_handler.n_dofs());
    m_system_rhs_im.reinit(dof_handler.n_dofs());
    m_system_rhs_2_re.reinit(dof_handler.n_dofs());
    m_system_rhs_2_im.reinit(dof_handler.n_dofs());
    m_Psi_re.reinit(dof_handler.n_dofs());
    m_Psi_im.reinit(dof_handler.n_dofs());
    m_Psi_sob_re.reinit(dof_handler.n_dofs());
    m_Psi_sob_im.reinit(dof_handler.n_dofs());
    m_sob_grad_re.reinit(dof_handler.n_dofs());
    m_sob_grad_im.reinit(dof_handler.n_dofs());
    m_workspace.reinit(dof_handler.n_dofs());
    m_error_per_cell.reinit(triangulation.n_active_cells());
    
    DynamicSparsityPattern dsp(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern (dof_handler, dsp);

    sparsity_pattern.copy_from(dsp);
    system_matrix.reinit (sparsity_pattern);
  }

  void MySolver::output_guess ()
  {
    VectorTools::interpolate ( dof_handler, m_Potential, m_workspace );
    VectorTools::interpolate ( dof_handler, m_F, m_system_rhs_re );

    DataOut<1> data_out;
    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (m_Psi_re, "Psi re");
    data_out.add_data_vector (m_Psi_im, "Psi im");
    data_out.add_data_vector (m_workspace, "V(r)");
    data_out.add_data_vector (m_system_rhs_re, "F(r)");
    data_out.build_patches ();

    ofstream output ("guess.gnuplot");
    data_out.write_gnuplot (output);
  }

  void MySolver::output_results ( string path, string prefix )
  {
    string filename = path + prefix + "-" + Utilities::int_to_string (m_counter,5) + ".gnuplot";

    DataOut<1> data_out;
    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (m_Psi_re, "Psi re");
    data_out.add_data_vector (m_Psi_im, "Psi im");
    data_out.add_data_vector (m_sob_grad_re, "sob grad re");
    data_out.add_data_vector (m_sob_grad_im, "sob grad im");
    data_out.build_patches ();

    ofstream output (filename.c_str());
    data_out.write_gnuplot (output);
  }

  int MySolver::DoIter ( string path )
  {
    int retval = Status::SUCCESS;
    
    m_table.clear();
    
    assemble_rhs();
    
    m_res=0;
    m_res_old=m_res;
    m_counter=0;
    
    do
    {
      cout << "--------------------------------------------------------------------------------" << endl;
      cout << "- " << path << " - " << m_counter << endl;

      assemble_system( m_system_rhs_re, m_sob_grad_re );
      solve( m_system_rhs_re, m_sob_grad_re );
      
      assemble_system( m_system_rhs_im, m_sob_grad_im );
      solve( m_system_rhs_im, m_sob_grad_im );
     
      assemble_system( m_system_rhs_2_re, m_Psi_sob_re );
      solve( m_system_rhs_2_re, m_Psi_sob_re );
      
      assemble_system( m_system_rhs_2_im, m_Psi_sob_im );
      solve( m_system_rhs_2_im, m_Psi_sob_im );

      Project_gradient();

      double tmp1 = m_sob_grad_re.l2_norm();
      double tmp2 = m_sob_grad_im.l2_norm();
      m_res = sqrt(tmp1*tmp1+tmp2*tmp2);

      m_Psi_re.add( -0.1, m_sob_grad_re ); 
      m_Psi_im.add( -0.1, m_sob_grad_im ); 
      
      m_N = Particle_Number();
      
      if( fabs(m_N-1) > 1e-4 )
      {
        m_Psi_re *= 1/sqrt(m_N);
        m_Psi_im *= 1/sqrt(m_N);
      }
      
      m_resp = m_res_old-m_res;
      m_res_old = m_res;

      if( m_counter % m_NA == 0 ) output_results(path);

      columns& cols = m_table.new_line();
      m_table.insert( cols, MyTable::COUNTER, double(m_counter) );
      m_table.insert( cols, MyTable::RES, m_res );
      m_table.insert( cols, MyTable::RESP, m_resp );
      m_table.insert( cols, MyTable::GS, m_gs );
      m_table.insert( cols, MyTable::PARTICLE_NUMBER, m_N );

      cout << m_table;
      if( m_res < m_epsilon[0] ) { break; }

      assemble_rhs();
      m_counter++;
    }while( true );
    
    if( m_N < 1e-5 ) retval = Status::ZERO_SOL;
    
    compute_mu_and_E();
    
    string filename = path + "log.csv";
    m_table.dump_2_file(filename);
    
    return retval;
  }

  void MySolver::run2b ()
  {
    string path;
    char shellcmd[255];
    double T, N, W;
    int status;

    compute_potmax();
    compute_potmin();
    
    cout << "m_potmin = " << m_potmin << endl;

    m_xmin = m_potmax;

    make_grid_custom();
    setup_system();
    
    compute_guess_parameter();

    map<string,double> constants;
    constants["PI"] = numbers::PI;

    string tmp1 = "(r-" + to_string(m_s) + ")";
    m_Psi_1_str = "exp(-" + tmp1 + "^2/" + to_string(m_alpha*m_alpha) + ")";
    cout << "m_Psi_1_str = " << m_Psi_1_str << endl;
    m_fctp_re.initialize( "r", m_Psi_1_str, constants );
    m_fctp_im.initialize( "r", "0*r", constants );

    VectorTools::interpolate (dof_handler, m_fctp_re, m_Psi_re );
    VectorTools::interpolate (dof_handler, m_fctp_im, m_Psi_im );

//    double fak = 1/sqrt(Particle_Number());
//    m_Psi_re *= fak;
//    m_Psi_im *= fak;

    output_guess();
    m_results.clear();
    for( int i=0; i<1; ++i )
    {
      sprintf( shellcmd, "mkdir %.4d/", i );
      system(shellcmd);
      sprintf( shellcmd, "%.4d/", i );
      path = shellcmd;
      
     
      status = DoIter(path);

      columns& cols = m_results.new_line();
      m_results.insert( cols, MyTable::MU, m_mu );
      m_results.insert( cols, MyTable::GS, m_gs );
      m_results.insert( cols, MyTable::PARTICLE_NUMBER, m_N );
      m_results.insert( cols, MyTable::COUNTER, double(m_counter) );
      m_results.insert( cols, MyTable::STATUS, double(status) );

      if( status == Status::SUCCESS )
      {
        m_counter=0;
        output_results(path,"Cfinal");
        dump_info_xml(path);
      }
      //compute_E_lin( m_Psi_ref, T, N, W ); 
    }
    m_results.dump_2_file( "results.csv" );
  }
  
  void MySolver::save( string filename )
  {
  }

  void MySolver::save_one( string filename )
  {
  }
  
  void MySolver::dump_info_xml ( string path )
  {
    string filename = path + "info.xml";

    wofstream fs2;
    fs2.open(filename);

    locale utf8_locale("en_US.UTF8");
    fs2.imbue(utf8_locale); 
    fs2 << L"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<INFO>\n";
    fs2 << L"<MU>" << m_mu << L"</MU>\n";
    fs2 << L"<GS>" << m_gs << L"</GS>\n";
    fs2 << L"<N>" << m_N << "L</N>\n";
    fs2 << L"<XMIN>" << m_xmin << L"</XMIN>\n";
    fs2 << L"<XMAX>" << m_xmax << L"</XMAX>\n";
    fs2 << L"<REVISION>" << STR2(GIT_SHA1) << L"</REVISION>\n";
    fs2 << L"</INFO>\n";
    fs2.close();
  }  
} // end of namespace 

int main ( int argc, char *argv[] )
{
  using namespace dealii;
  deallog.depth_console (0);

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv );
  {
    BreedSolver::MySolver solver("params.xml");
    solver.run2b ();
  }
return EXIT_SUCCESS;
}