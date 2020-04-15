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

    void run ();
    void run2b ();

    void compute_J( double& );
    void set_t( const double a, const double b ) { m_t[0]=a; m_t[1]=b; };
    
  protected:
    FunctionParser<1> m_Potential;

    double m_t[2];
    double m_xmin, m_xmax;
    double m_res;
    double m_res_old;
    double m_resp;
    double m_ti;
    double m_final_error;
    double m_N;
    
    double m_mu;
    double m_Gamma;
    double m_dmu;
    double m_gs;
    vector<double> m_epsilon;
    vector<double> m_omega;

    unsigned m_counter;
    unsigned m_NA;
    unsigned m_Ndmu;
    unsigned m_global_refinements;
        
    double l2norm_t() { return sqrt(m_t[0]*m_t[0]+m_t[1]*m_t[1]); };
    
    int DoIter( string="" );

    void make_grid_custom();
    void setup_system( const bool );
    void assemble_rhs();
    void assemble_system();
    void do_superposition();
    void save( string );
    void save_one( string );
    void dump_info_xml( string="" );
    
    double Particle_Number( Vector<double>& );
    
    void solve();
    void compute_E_lin( Vector<double>&, double&, double&, double& );
    void estimate_error( double& );
    void find_min_J();

    void output_results ( string, string = "step" );
    void output_guess ();

    MyParameterHandler m_ph;
    Triangulation<1> triangulation;
    FESystem<1> fe;
    SparsityPattern sparsity_pattern;
    DoFHandler<1> dof_handler;

    SparseMatrix<double> system_matrix;
    Vector<double> newton_update;
    Vector<double> system_rhs;
    Vector<double> m_Psi_ref;
    Vector<double> m_Psi_1;
    Vector<double> m_Psi_2;
    Vector<double> m_workspace;
    Vector<double> m_error_per_cell;

    string m_guess_str;
    string m_potential_str;
    MyTable m_table;
    MyTable m_results;
  };

  double fun_neu( const gsl_vector* x, void *params )
  {
    MySolver* sol = reinterpret_cast<MySolver*>(params); 
    double retval=0;
    
    sol->set_t( gsl_vector_get(x,0), gsl_vector_get(x,1) );
    sol->compute_J(retval);
  return retval;
  }

/*************************************************************************************************/
/**
 * Constructor
 */
  MySolver::MySolver ( const string xmlfilename ) 
    : 
    m_ph(xmlfilename),
    triangulation (Triangulation<1>::MeshSmoothing(Triangulation<1>::limit_level_difference_at_vertices|Triangulation<1>::eliminate_refined_inner_islands|Triangulation<1>::smoothing_on_refinement|Triangulation<1>::smoothing_on_coarsening)),
    fe ( FE_Q<1>(gl_degree_fe), 2),
    dof_handler (triangulation)
  {
    try
    {
      m_omega = m_ph.Get_Physics("omega");
      m_gs = m_ph.Get_Physics("gs_1",0);

      m_xmin = m_ph.Get_Mesh("xrange",0);
      m_xmax = m_ph.Get_Mesh("xrange",1);
      m_global_refinements = unsigned(m_ph.Get_Mesh("global_refinements",0));

      m_ti = m_ph.Get_Algorithm("ti",0);
      m_t[0] = m_ti;
      m_t[1] = m_ti;
      m_NA = int(m_ph.Get_Algorithm("NA",0));
      m_Ndmu = m_ph.Get_Algorithm("Ndmu",0); 
      m_dmu = m_ph.Get_Algorithm("dmu",0);
      m_epsilon = m_ph.Get_Algorithm("epsilon");
      m_Gamma=-0.0001;
    }
    catch( const std::string info )
    {
      std::cerr << info << endl;
      exit(0);
    }    
    m_counter=0;
    m_final_error=0;    
    
    m_guess_str = "exp(-0.5*x*x)";
    m_potential_str = "0.5*x*x";
    
    map<string,double> constants;
    constants["PI"] = numbers::PI;
    constants["omega"] = m_omega[0];
    m_Potential.initialize( FunctionParser<1>::default_variable_names(), m_potential_str, constants );
  }

  MySolver::~MySolver ()
  {
    dof_handler.clear ();
  }
  
  void MySolver::compute_E_lin( Vector<double>& vec, double& T, double& N, double& W )
  {
    const QGauss<1>  quadrature_formula(fe.degree+1);
    FEValues<1> fe_values (fe, quadrature_formula, update_gradients|update_values|update_JxW_values|update_quadrature_points);

    const unsigned n_q_points = quadrature_formula.size();
    vector<Vector<double>> vals(n_q_points,Vector<double>(2));    
    vector<vector<Tensor<1,1>>> grad(n_q_points, vector<Tensor<1,1>>(2));
    
    double JxW, vec_val_q;
    
    T=0; N=0; W=0;
    typename DoFHandler<1>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for( ; cell!=endc; ++cell )
    {
      fe_values.reinit (cell);
      fe_values.get_function_values(vec, vals);
      fe_values.get_function_gradients(vec, grad);
      for( unsigned qp=0; qp<n_q_points; ++qp )
      {
        vec_val_q = vals[qp]*vals[qp];
        JxW = fe_values.JxW(qp);
        T += JxW*( grad[qp][0]*grad[qp][0] + grad[qp][1]*grad[qp][1] + m_Potential.value(fe_values.quadrature_point(qp))*vec_val_q );
        N += JxW*vec_val_q;
        W += JxW*vec_val_q*vec_val_q;
      }
    }
  }

  double MySolver::Particle_Number( Vector<double>& vec )
  {
    double retval=0;
   
    const QGauss<1>  quadrature_formula(fe.degree+1);
    FEValues<1> fe_values (fe, quadrature_formula, update_values|update_JxW_values|update_quadrature_points);

    const unsigned n_q_points = quadrature_formula.size();
    vector<Vector<double>> vals(n_q_points,Vector<double>(2));

    DoFHandler<1>::active_cell_iterator cell=dof_handler.begin_active(), endc=dof_handler.end();
    for( ; cell!=endc; ++cell )
    {
      fe_values.reinit (cell);
      fe_values.get_function_values( vec, vals );      

      for( unsigned qp=0; qp<n_q_points; ++qp )
        retval += fe_values.JxW(qp)*(vals[qp]*vals[qp]);
    }
  return retval;
  }
  
  void MySolver::estimate_error ( double& err )
  {
    /*
    const FEValuesExtractors::Scalar rt (0);
    const FEValuesExtractors::Scalar it (1);    
    
    CPotential<dim> Potential( m_omega );
    const QGauss<dim> quadrature_formula(fe.degree+1);
    //const QGauss<dim-1> face_quadrature_formula(fe.degree+1);
   
    system_rhs=0;
    system_matrix=0;

    FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_gradients|update_JxW_values|update_quadrature_points);

    const unsigned dofs_per_cell = fe.dofs_per_cell;
    const unsigned n_q_points    = quadrature_formula.size();

    Vector<double> cell_rhs (dofs_per_cell);
    FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);

    vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    vector<Vector<double>> vals(n_q_points,Vector<double>(2));
    vector<vector<Tensor<1,dim>>> grads(n_q_points, vector<Tensor<1,dim>>(2));

    double JxW, Q1, pq, tmp;
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for ( ; cell!=endc; ++cell )
    {
        cell_rhs=0;
        cell_matrix=0;

        fe_values.reinit (cell);
        fe_values.get_function_values(m_workspace_1, vals);
        fe_values.get_function_gradients(m_workspace_1, grads);

        for( unsigned qp=0; qp<n_q_points; ++qp )
        {
          JxW = fe_values.JxW(qp);
          Q1 = Potential.value(fe_values.quadrature_point(qp)) - m_mu + m_gs*(vals[qp]*vals[qp]);

          for( unsigned i=0; i<dofs_per_cell; ++i )
          {
            cell_rhs(i) += JxW*(grads[qp][0]*fe_values[rt].gradient(i,qp) + Q1*vals[qp][0]*fe_values[rt].value(i,qp) + grads[qp][1]*fe_values[it].gradient(i,qp) + Q1*vals[qp][1]*fe_values[it].value(i,qp));
            for( unsigned j=0; j<dofs_per_cell; ++j )
              cell_matrix(i,j) += JxW*(fe_values[rt].value(i,qp)*fe_values[rt].value(j,qp)+fe_values[it].value(i,qp)*fe_values[it].value(j,qp));
          }
        }
    }

    solve();
    
    m_workspace_1=newton_update;
    VectorTools::integrate_difference ( dof_handler, m_workspace_1, ZeroFunction<dim>(2), m_error_per_cell, QGauss<dim>(fe.degree+2), VectorTools::L2_norm);    
    const double total_local_error = m_error_per_cell.l2_norm();
    err = std::sqrt (Utilities::MPI::sum (total_local_error * total_local_error, MPI_COMM_WORLD)); 
    */
  }
  
  void MySolver::assemble_rhs ()
  {
    /*
    vector<bool> boundary_dofs (dof_handler.n_dofs());
    DoFTools::extract_boundary_dofs(dof_handler, ComponentMask(), boundary_dofs);
    for( unsigned i=0; i<dof_handler.n_dofs(); ++i )
      if( boundary_dofs[i] ) m_Psi_ref(i)=0;        
    */

    system_rhs=0;
    
    const FEValuesExtractors::Scalar rt (0);
    const FEValuesExtractors::Scalar it (1);
    const QGauss<1> quadrature_formula(fe.degree+1);

    FEValues<1> fe_values (fe, quadrature_formula, update_values|update_gradients|update_JxW_values|update_quadrature_points);

    const unsigned dofs_per_cell = fe.dofs_per_cell;
    const unsigned n_q_points = quadrature_formula.size();

    vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    
    vector<Vector<double>> vals(n_q_points,Vector<double>(2));
    vector<vector<Tensor<1,1>>> grads(n_q_points, vector<Tensor<1,1>>(2));
    
    double JxW, Pot, req, imq, tmp;
    DoFHandler<1>::active_cell_iterator cell=dof_handler.begin_active(), endc=dof_handler.end();
    for( ; cell!=endc; ++cell )
    {
      fe_values.reinit (cell);
      fe_values.get_function_values(m_Psi_ref, vals);
      fe_values.get_function_gradients(m_Psi_ref, grads);
      cell->get_dof_indices (local_dof_indices);
        
      for( unsigned qp=0; qp<n_q_points; ++qp )
      {
        JxW = fe_values.JxW(qp);
        Pot = m_Potential.value(fe_values.quadrature_point(qp))-m_mu;
        req = vals[qp][0]*vals[qp][0];
        imq = vals[qp][1]*vals[qp][1];

        for( unsigned i=0; i<dofs_per_cell; ++i )
            system_rhs(local_dof_indices[i]) += JxW*(grads[qp][0]*fe_values[rt].gradient(i,qp) + (Pot+m_gs*(req-3*imq))*vals[qp][0]*fe_values[rt].value(i,qp) + m_Gamma*vals[qp][1]*fe_values[rt].value(i,qp) 
                                                    +grads[qp][1]*fe_values[it].gradient(i,qp) + (Pot+m_gs*(3*req-imq))*vals[qp][1]*fe_values[it].value(i,qp) - m_Gamma*vals[qp][0]*fe_values[it].value(i,qp) );
      }
    }
    m_res = system_rhs.l2_norm();
  }

  void MySolver::assemble_system ()
  {
    /*
    vector<bool> boundary_dofs (dof_handler.n_dofs());
    DoFTools::extract_boundary_dofs(dof_handler, ComponentMask(), boundary_dofs);
    for( unsigned i=0; i<dof_handler.n_dofs(); ++i )
      if( boundary_dofs[i] ) m_Psi_ref(i)=0;        
    */
   
    system_matrix=0;
    
    const FEValuesExtractors::Scalar rt (0);
    const FEValuesExtractors::Scalar it (1);
    const QGauss<1> quadrature_formula(fe.degree+1);

    FEValues<1> fe_values (fe, quadrature_formula, update_values|update_gradients|update_JxW_values|update_quadrature_points);

    const unsigned dofs_per_cell = fe.dofs_per_cell;
    const unsigned n_q_points = quadrature_formula.size();

    FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);

    vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    vector<vector<Tensor<1,1>>> Psi_ref_grad(n_q_points, vector<Tensor<1,1>>(2));
    vector<Vector<double>> Psi_ref(n_q_points,Vector<double>(2));

    double JxW, Pot, fak;
    DoFHandler<1>::active_cell_iterator cell=dof_handler.begin_active(), endc=dof_handler.end();
    for( ; cell!=endc; ++cell )
    {
      cell_matrix=0;

      fe_values.reinit (cell);
      fe_values.get_function_values(m_Psi_ref, Psi_ref);
      fe_values.get_function_gradients(m_Psi_ref, Psi_ref_grad);
      cell->get_dof_indices (local_dof_indices);

      for( unsigned qp=0; qp<n_q_points; ++qp )
      {
        JxW = fe_values.JxW(qp);
        Pot = m_Potential.value(fe_values.quadrature_point(qp)) - m_mu + 3*m_gs*(Psi_ref[qp][0]*Psi_ref[qp][0] - Psi_ref[qp][1]*Psi_ref[qp][1]);
        fak = 6*m_gs*Psi_ref[qp][0]*Psi_ref[qp][1];

        for( unsigned i=0; i<dofs_per_cell; ++i )
          for( unsigned j=0; j<dofs_per_cell; ++j )
            cell_matrix(i,j) += JxW*(fe_values[rt].gradient(i,qp)*fe_values[rt].gradient(j,qp) + Pot*fe_values[rt].value(i,qp)*fe_values[rt].value(j,qp) 
                                    +(m_Gamma-fak)*fe_values[rt].value(i,qp)*fe_values[it].value(j,qp) 
                                    +(fak-m_Gamma)*fe_values[it].value(i,qp)*fe_values[rt].value(j,qp)
                                    +fe_values[it].gradient(i,qp)*fe_values[it].gradient(j,qp) + Pot*fe_values[it].value(i,qp)*fe_values[it].value(j,qp));
      }
      
      for( unsigned i=0; i<dofs_per_cell; ++i )
      {
        for( unsigned j=0; j<dofs_per_cell; ++j ) 
          system_matrix.add (local_dof_indices[i],local_dof_indices[j], cell_matrix(i,j));
      }      
    }
    
    map<types::global_dof_index,double> boundary_values;
    VectorTools::interpolate_boundary_values (dof_handler, 0, ZeroFunction<1>(2), boundary_values);
    MatrixTools::apply_boundary_values (boundary_values, system_matrix, newton_update, system_rhs);
  }
  
  void MySolver::compute_J( double& retval )
  {
    do_superposition();
    
    const QGauss<1>  quadrature_formula(fe.degree+1);
    FEValues<1> fe_values (fe, quadrature_formula, update_gradients|update_values|update_JxW_values|update_quadrature_points);

    const unsigned n_q_points = quadrature_formula.size();
    vector<Vector<double>> vals(n_q_points,Vector<double>(2));
    vector<vector<Tensor<1,1>>> grads(n_q_points, vector<Tensor<1,1>>(2));

    double psum=0, req, imq;
    DoFHandler<1>::active_cell_iterator cell=dof_handler.begin_active(), endc=dof_handler.end();
    for( ; cell!=endc; ++cell )
    {
      fe_values.reinit (cell);
      fe_values.get_function_values(m_Psi_ref, vals);
      fe_values.get_function_gradients(m_Psi_ref, grads);

      for( unsigned qp=0; qp<n_q_points; ++qp )
      {
        req = vals[qp][0]*vals[qp][0];
        imq = vals[qp][1]*vals[qp][1];
        psum += fe_values.JxW(qp)*( grads[qp][0]*grads[qp][0] - grads[qp][1]*grads[qp][1] + (m_Potential.value(fe_values.quadrature_point(qp))-m_mu)*(req-imq) + 0.5*m_gs*(req*req+imq*imq-6*req*imq) - m_Gamma*vals[qp][0]*vals[qp][1] ); 
      }
    }
    retval = 0.5*psum;
  }
  
  void MySolver::find_min_J()
  {
    size_t iter=0;
    int status;
    double size;

    gsl_multimin_function my_func;
    my_func.n = 2;
    my_func.f = fun_neu;
    my_func.params = this;

    gsl_vector * ss = gsl_vector_alloc (2);
    gsl_vector * x = gsl_vector_alloc (2);
    gsl_vector_set_all (ss, 1.0);
    gsl_vector_set (x, 0, m_t[0]);
    gsl_vector_set (x, 1, m_t[1]);

    const gsl_multimin_fminimizer_type * T = gsl_multimin_fminimizer_nmsimplex2;
    gsl_multimin_fminimizer * s = gsl_multimin_fminimizer_alloc (T, 2);

    gsl_multimin_fminimizer_set (s, &my_func, x, ss);

    do
    {
      iter++;
      status = gsl_multimin_fminimizer_iterate (s);

      if (status) break;

      size = gsl_multimin_fminimizer_size (s);
      
      //printf ("%5d %10.3e %10.3e f() = %7.3f size = %.3f\n", iter, gsl_vector_get (s->x, 0), gsl_vector_get (s->x, 1), s->fval, size);

      status = gsl_multimin_test_size(size,1e-6);
    }
    while (status == GSL_CONTINUE && iter < 1000);

    m_t[0] = gsl_vector_get (s->x, 0);
    m_t[1] = gsl_vector_get (s->x, 1);

    gsl_multimin_fminimizer_free (s);
    gsl_vector_free(ss);
    gsl_vector_free (x);
  }    
  
  void MySolver::solve ()
  {
    SparseDirectUMFPACK A_direct;
    A_direct.initialize(system_matrix);
    A_direct.vmult (newton_update, system_rhs);
  }
  
  void MySolver::make_grid_custom ()
  {
    Point<1,double> pt1(m_xmin); 
    Point<1,double> pt2(m_xmax);
    
    GridGenerator::hyper_rectangle(triangulation, pt2, pt1);
    triangulation.refine_global(m_global_refinements);
/*    
    double isovalues[] = {15,13,11};

    for( unsigned int step=0; step<sizeof(isovalues)/sizeof(double); step++ )
    {
      typename parallel::distributed::Triangulation::active_cell_iterator cell = triangulation.begin_active(), endc = triangulation.end();
      for( ; cell!=endc; ++cell )
        for( unsigned v=0; v<GeometryInfo<dim>::vertices_per_cell; v++ )
        {
          Point<dim> p = cell->vertex(v);
          if( fabs(Potential_fct.value(p))  < isovalues[step] )
          {
            cell->set_refine_flag ();
            break;
          }
        }
      triangulation.execute_coarsening_and_refinement ();
    }
*/
/*
    GridOutFlags::Msh opt(true, true);
    string filename = "grid-" + Utilities::int_to_string(triangulation.locally_owned_subdomain(), 4);
    ofstream output ((filename + ".msh").c_str());
    GridOut grid_out;
    grid_out.set_flags(opt);
    grid_out.write_msh (triangulation,output);
*/
  }
  
  void MySolver::do_superposition()
  {
    m_Psi_ref=0;
    m_Psi_ref.add(m_t[0],m_Psi_1,m_t[1],m_Psi_2);
    
    vector<bool> boundary_dofs (dof_handler.n_dofs());
    DoFTools::extract_boundary_dofs(dof_handler, ComponentMask(), boundary_dofs);
    for( unsigned i=0; i<dof_handler.n_dofs(); ++i )
      if( boundary_dofs[i] ) m_Psi_ref(i)=0;          
  }

  void MySolver::setup_system( const bool initial_step )
  {
    if( initial_step )
    {
      dof_handler.distribute_dofs (fe);

      m_Psi_1.reinit (dof_handler.n_dofs());
      m_Psi_2.reinit (dof_handler.n_dofs());
    }

    m_Psi_ref.reinit (dof_handler.n_dofs());
    newton_update.reinit (dof_handler.n_dofs());
    system_rhs.reinit (dof_handler.n_dofs());
    m_workspace.reinit (dof_handler.n_dofs());
    m_error_per_cell.reinit(triangulation.n_active_cells());
    
    DynamicSparsityPattern dsp(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern (dof_handler, dsp);

    sparsity_pattern.copy_from(dsp);
    system_matrix.reinit (sparsity_pattern);
  }

  void MySolver::output_guess ()
  {
  }

  void MySolver::output_results ( string path, string prefix )
  {
    string filename = path + prefix + "-" + Utilities::int_to_string (m_counter,5) + ".gnuplot";

    DataOut<1> data_out;
    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (m_Psi_ref, "Psi_ref");
    //data_out.add_data_vector (m_error_per_cell, "error_per_cell");
    data_out.build_patches ();

    ofstream output (filename.c_str());
    data_out.write_gnuplot (output);
  }

  int MySolver::DoIter ( string path )
  {
    int retval = Status::SUCCESS;
    
    m_table.clear();
    
    m_t[0] = m_ti;
    m_t[1] = m_ti;
    
    do_superposition();
    assemble_rhs();
    
    m_res=0;
    m_res_old=m_res;
    m_counter=0;
    do
    {
      cout << "--------------------------------------------------------------------------------" << endl;
      cout << "- " << path << " - " << m_counter << endl;

      assemble_system();
      solve();

      m_Psi_2.add( -0.1*m_t[1]/fabs(m_t[1]), newton_update); 

      find_min_J();
      
      do_superposition();
      m_N = Particle_Number(m_Psi_ref);
      assemble_rhs();
      
      m_resp = m_res_old-m_res;
      m_res_old = m_res;

      if( m_counter % m_NA == 0 ) output_results(path);

      columns& cols = m_table.new_line();
      m_table.insert( cols, MyTable::COUNTER, double(m_counter) );
      m_table.insert( cols, MyTable::RES, m_res );
      m_table.insert( cols, MyTable::RESP, m_resp );
      m_table.insert( cols, MyTable::MU, m_mu );
      m_table.insert( cols, MyTable::GS, m_gs );
      m_table.insert( cols, MyTable::t1, m_t[0] );
      m_table.insert( cols, MyTable::t2, m_t[1] );
      m_table.insert( cols, MyTable::l2norm_t, l2norm_t() );
      m_table.insert( cols, MyTable::PARTICLE_NUMBER, m_N );
      //m_table.insert( cols, MyTable::total_no_cells, double(m_total_no_cells) );
      //m_table.insert( cols, MyTable::total_no_active_cells, double(m_total_no_active_cells) );

      cout << m_table;
      if( m_res < m_epsilon[0] ) { retval=Status::SUCCESS; break; }
      if( l2norm_t() < 1e-4 ) { retval=Status::ZERO_SOL; break; }

      m_counter++;
    }while( true );
    
    // Standard Newton 
    do
    {
      cout << "--------------------------------------------------------------------------------" << endl;
      cout << "-- " << path << " - " << m_counter << endl;

      assemble_system();
      solve();

      m_Psi_ref.add( -0.1, newton_update); 

      m_N = Particle_Number(m_Psi_ref);
      assemble_rhs();
      
      m_resp = m_res_old-m_res;
      m_res_old = m_res;

      if( m_counter % m_NA == 0 ) output_results(path);

      columns& cols = m_table.new_line();
      m_table.insert( cols, MyTable::COUNTER, double(m_counter) );
      m_table.insert( cols, MyTable::RES, m_res );
      m_table.insert( cols, MyTable::RESP, m_resp );
      m_table.insert( cols, MyTable::MU, m_mu );
      m_table.insert( cols, MyTable::GS, m_gs );
      m_table.insert( cols, MyTable::t1, m_t[0] );
      m_table.insert( cols, MyTable::t2, m_t[1] );
      m_table.insert( cols, MyTable::l2norm_t, l2norm_t() );
      m_table.insert( cols, MyTable::PARTICLE_NUMBER, m_N );
      //m_table.insert( cols, MyTable::total_no_cells, double(m_total_no_cells) );
      //m_table.insert( cols, MyTable::total_no_active_cells, double(m_total_no_active_cells) );

      cout << m_table;
      if( m_res < m_epsilon[1] ) { retval=Status::SUCCESS; break; }

      m_counter++;
    }while( true );    

    do_superposition();
    m_N = Particle_Number(m_Psi_ref);
    
    if( m_N < 1e-5 ) retval = Status::ZERO_SOL;
    
    string filename = path + "log.csv";
    m_table.dump_2_file(filename);
    
    return retval;
  }

  void MySolver::run()
  {
    int status;
    string path;
    char shellcmd[255];
    double T, N, W;

    make_grid_custom();
    setup_system(true);

    map<string,double> constants;
    constants["PI"] = numbers::PI;
    constants["omega"] = m_omega[0];
    FunctionParser<1> guess_fct;
    guess_fct.initialize( FunctionParser<1>::default_variable_names(), m_guess_str, constants );

    VectorTools::interpolate (dof_handler, guess_fct, m_Psi_1);
    m_Psi_1 *= 1.0/sqrt(Particle_Number(m_Psi_1));
    m_Psi_2 = 0;
      
    compute_E_lin( m_Psi_1, T, N, W );
    output_guess();

    m_mu = T/N+m_gs/fabs(m_gs)*m_dmu;

    cout << "T = " << T << endl;
    cout << "N = " << N << endl;
    cout << "W = " << W << endl;
    cout << "m_mu = " << m_mu << endl;

    status = DoIter();
    
    cout << "L2_norm of m_Psi_ref: " << m_N << endl;

    if( status == Status::SUCCESS )
    {
      output_results("","final");
      dump_info_xml();
    }

    ofstream ofs("log.txt");
    ofs << m_table;
  }

  void MySolver::run2b ()
  {
    string path;
    char shellcmd[255];
    double T, N, W;
    int status;

    make_grid_custom();
    setup_system(true);

    map<string,double> constants;
    constants["PI"] = numbers::PI;
    constants["omega"] = m_omega[0];
    FunctionParser<1> guess_fct;
    guess_fct.initialize( FunctionParser<1>::default_variable_names(), m_guess_str, constants );
    VectorTools::interpolate (dof_handler, guess_fct, m_Psi_1 );

    m_Psi_1 *= 1.0/sqrt(Particle_Number(m_Psi_1));
    m_Psi_2 = 0; 

    compute_E_lin( m_Psi_1, T, N, W );
    double m_mu_0 = T/N;
    m_mu = ceil(10.0*m_mu_0)/10.0 + m_gs/fabs(m_gs)*m_dmu;
    
    output_guess();
    m_results.clear();
    for( int i=0; i<m_Ndmu; ++i )
    {
      sprintf( shellcmd, "mkdir %.4d/", i );
      system(shellcmd);
      sprintf( shellcmd, "%.4d/", i );
      path = shellcmd;
      
      // nehari
      //m_ti = sqrt((m_mu*N-T)/(4.0*m_gs*W)); // if m_Psi_2 == m_Psi_1
      m_ti = sqrt((m_mu*N-T)/(m_gs*W));

      cout << "T = " << T << endl;
      cout << "N = " << N << endl;
      cout << "W = " << W << endl;
      cout << "m_mu = " << m_mu << endl;      
      cout << "m_ti = " << m_ti << endl;      
      
      status = DoIter(path);

      columns& cols = m_results.new_line();
      m_results.insert( cols, MyTable::MU, m_mu );
      m_results.insert( cols, MyTable::GS, m_gs );
      m_results.insert( cols, MyTable::PARTICLE_NUMBER, m_N );
      m_results.insert( cols, MyTable::COUNTER, double(m_counter) );
      m_results.insert( cols, MyTable::STATUS, double(status) );

      if( status == Status::SUCCESS )
      {
        estimate_error(m_final_error);
        output_results(path,"Cfinal");
        string filename = path + "Cfinal.bin";
        save( filename );
        filename = path + "Cfinal-1.bin";
        save_one( filename );
        dump_info_xml(path);
        m_Psi_1 = m_Psi_ref;
        m_Psi_2 = 0;
      }
      else if( status == Status::SLOW_CONV )
      {
        m_Psi_2 = 0; 
      }
      else
      {
        break;
      }
      compute_E_lin( m_Psi_ref, T, N, W ); // TODO: kommentier mich aus, falls ich kein nehari reset habe
      m_mu += m_gs/fabs(m_gs)*m_dmu;
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
    fs2 << L"<FINAL_ERROR>" << m_final_error << L"</FINAL_ERROR>\n";
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