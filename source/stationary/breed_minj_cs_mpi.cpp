//
// atus-pro testing - atus-pro testing playgroung
// Copyright (C) 2020 Želimir Marojević <zelimir.marojevic@gmail.com>
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

#include "default_includes.h"

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
#include "mpi.h"
#include "functions_cs.h"
#include "MyParameterHandler.h"
#include "my_table.h"

#define STR1(x) #x
#define STR2(x) STR1(x)
#define _FAKTOR_ 6.28318530718
//#define _FAKTOR_ 1

namespace BreedSolver
{
  using namespace std;
  using namespace dealii;

  template <int dim>
  class MySolver; // class forward declaration
  
  enum Status { SUCCESS, FAILED, ZERO_SOL, SLOW_CONV };

  template<int dim>
  double fun_neu( const gsl_vector* x, void *params )
  {
    MySolver<dim>* sol = reinterpret_cast<MySolver<dim>*>(params); 
    double retval=0;
    
    sol->set_t( gsl_vector_get(x,0), gsl_vector_get(x,1) );
    sol->compute_J(retval);
  return retval;
  }
  
  template <int dim>
  class MySolver
  {
  public:
    explicit MySolver( const std::string );
    virtual ~MySolver();

    void run ();
    void run2b ();

    void compute_J( double& );
    void set_t( const double a, const double b ) { m_t[0]=a; m_t[1]=b; };
    
    MPI_Comm mpi_communicator;
  protected:
    double m_t[2];
    double m_xmin, m_xmax;
    double m_ymin, m_ymax;
    double m_zmin, m_zmax;
    double m_res;
    double m_res_old;
    double m_resp;
    double m_ti;
    double m_final_error;
    double m_N;
    
    double m_mu;
    double m_dmu;
    double m_gs;
    vector<double> m_epsilon;
    vector<double> m_omega;

    bool m_root;
    int m_myrank;
    
    unsigned m_counter;
    unsigned m_global_refinement;
    unsigned m_total_no_cells;
    unsigned m_total_no_active_cells;    
    unsigned m_NA;
    unsigned m_Ndmu;
    unsigned m_QN1[3];
        
    double l2norm_t() { return sqrt(m_t[0]*m_t[0]+m_t[1]*m_t[1]); };
    
    int DoIter( string="" );

    void make_grid_custom();
    void setup_system( const bool );
    void assemble_rhs();
    void assemble_system();
    void do_superposition();
    void save( string );
    void save_one( string );
    void dump_info_xml( const string );
    
    double Particle_Number( LA::MPI::Vector& );
    
    void solve();
    void compute_E_lin( LA::MPI::Vector&, double&, double&, double& );
    void estimate_error( double& );
    void find_min_J();

    void output_results ( string, string = "step" );
    void output_guess ();

    ofstream m_computing_timer_log;
    TimerOutput m_computing_timer;    
    MyParameterHandler m_ph;
    ConditionalOStream pcout;
    parallel::distributed::Triangulation<dim> triangulation;
    FESystem<dim> fe;
    DoFHandler<dim> dof_handler;
    IndexSet locally_owned_dofs;
    IndexSet locally_relevant_dofs;
    AffineConstraints<double> constraints;

    LA::MPI::SparseMatrix system_matrix;
    LA::MPI::Vector newton_update;
    LA::MPI::Vector system_rhs;
    LA::MPI::Vector m_Psi_ref;
    LA::MPI::Vector m_Psi_1;
    LA::MPI::Vector m_Psi_2;
    LA::MPI::Vector m_workspace_1;
    LA::MPI::Vector m_workspace_2;
    LA::MPI::Vector m_workspace_ng;
    Vector<double> m_error_per_cell;

    string m_guess_str;
    MyTable m_table;
    MyTable m_results;
  };

/*************************************************************************************************/
/**
 * Constructor
 */
  template <int dim>
  MySolver<dim>::MySolver ( const std::string xmlfilename ) 
    : 
    mpi_communicator(MPI_COMM_WORLD), 
    m_computing_timer_log("benchmark.txt"),
    m_computing_timer(mpi_communicator, m_computing_timer_log, TimerOutput::summary, TimerOutput::cpu_and_wall_times ), 
    m_ph(xmlfilename),
    pcout (cout, (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)),
    triangulation (mpi_communicator, typename Triangulation<dim>::MeshSmoothing(Triangulation<dim>::limit_level_difference_at_vertices|Triangulation<dim>::eliminate_refined_inner_islands|Triangulation<dim>::smoothing_on_refinement|Triangulation<dim>::smoothing_on_coarsening)),
    fe (FE_Q<dim>(gl_degree_fe), 2),
    dof_handler (triangulation)
  {
    try
    {
      m_omega = m_ph.Get_Physics("omega");
      m_gs = m_ph.Get_Physics("gs_1",0);
      m_QN1[0] = int(m_ph.Get_Physics("QN1",0));
      m_QN1[1] = int(m_ph.Get_Physics("QN1",1));
      m_QN1[2] = int(m_ph.Get_Physics("QN1",2));

      m_xmin = m_ph.Get_Mesh("xrange",0);
      m_xmax = m_ph.Get_Mesh("xrange",1);
      m_ymin = m_ph.Get_Mesh("yrange",0);
      m_ymax = m_ph.Get_Mesh("yrange",1);
      m_zmin = m_ph.Get_Mesh("zrange",0);
      m_zmax = m_ph.Get_Mesh("zrange",1);
      
      m_ti = m_ph.Get_Algorithm("ti",0);
      m_t[0] = m_ti;
      m_t[1] = m_ti;

      m_NA = int(m_ph.Get_Algorithm("NA",0));
      m_Ndmu = m_ph.Get_Algorithm("Ndmu",0); 
      m_dmu = m_ph.Get_Algorithm("dmu",0);
    }
    catch( const std::string& info )
    {
      std::cerr << info << endl;
      MPI_Abort( mpi_communicator, 0 );
    }    

    m_root = (Utilities::MPI::this_mpi_process(mpi_communicator) == 0);
    MPI_Comm_rank(mpi_communicator, &m_myrank);
    m_counter = 0;
    m_final_error=0;    
  }

  template <int dim>
  MySolver<dim>::~MySolver ()
  {
    dof_handler.clear ();
  }
  
  template <int dim>
  void MySolver<dim>::compute_E_lin( LA::MPI::Vector& vec, double& T, double& N, double& W )
  {
    TimerOutput::Scope timing_section(m_computing_timer, "");
    
    constraints.distribute(vec);
    m_workspace_1 = vec;
    
    CPotential Potential( m_omega, m_QN1[2] );
    const QGauss<dim>  quadrature_formula(fe.degree+1);
    FEValues<dim> fe_values (fe, quadrature_formula, update_gradients|update_values|update_JxW_values|update_quadrature_points);

    const unsigned int n_q_points = quadrature_formula.size();
    vector<Vector<double>> vals(n_q_points,Vector<double>(2));    
    vector<vector<Tensor<1,dim>>> grad(n_q_points, vector<Tensor<1,dim>>(2));
    
    double tmp[]={0,0,0}, res[]={0,0,0};
    
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for ( ; cell!=endc; ++cell)
    {
      if( cell->is_locally_owned() )
      {
        fe_values.reinit (cell);
        fe_values.get_function_values(m_workspace_1, vals);
        fe_values.get_function_gradients(m_workspace_1, grad);
        for ( unsigned qp=0; qp<n_q_points; ++qp )
        {
          const double vec_val_q = vals[qp]*vals[qp];
          const double JxW = fabs(fe_values.quadrature_point(qp)[1])*fe_values.JxW(qp);
          tmp[0] += JxW*( grad[qp][0]*grad[qp][0] + grad[qp][1]*grad[qp][1] + Potential.value(fe_values.quadrature_point(qp))*vec_val_q );
          tmp[1] += JxW*vec_val_q;
          tmp[2] += JxW*vec_val_q*vec_val_q;
        }
      }
    }

    MPI_Allreduce( tmp, res, 3, MPI_DOUBLE, MPI_SUM, mpi_communicator);
    T=res[0]; 
    N=res[1];
    W=res[2];
    
    
  }

  template<int dim>
  double MySolver<dim>::Particle_Number( LA::MPI::Vector& vec )
  {
    TimerOutput::Scope timing_section(m_computing_timer, "");
    double tmp1=0, retval=0;
   
    constraints.distribute(vec);
    m_workspace_1 = vec;

    const QGauss<dim>  quadrature_formula(fe.degree+1);
    FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_JxW_values|update_quadrature_points);

    const unsigned int n_q_points = quadrature_formula.size();
    vector<Vector<double>> vals(n_q_points,Vector<double>(2));

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
      if( cell->is_locally_owned() )
      {
        fe_values.reinit (cell);
        fe_values.get_function_values( m_workspace_1, vals );      

        for ( unsigned qp=0; qp<n_q_points; ++qp )
          tmp1 += fabs(fe_values.quadrature_point(qp)[1])*fe_values.JxW(qp)*(vals[qp]*vals[qp]);
      }
    }
    tmp1 *= _FAKTOR_;
    
    MPI_Allreduce( &tmp1, &retval, 1, MPI_DOUBLE, MPI_SUM, mpi_communicator);
    
  return retval;
  }
  
  template <int dim>
  void MySolver<dim>::estimate_error ( double& err )
  {
    TimerOutput::Scope timing_section(m_computing_timer, "");

    const FEValuesExtractors::Scalar rt (0);
    const FEValuesExtractors::Scalar it (1);    
    
    CPotential Potential( m_omega, m_QN1[2] );
    const QGauss<dim> quadrature_formula(fe.degree+1);
    //const QGauss<dim-1> face_quadrature_formula(fe.degree+1);
   
    constraints.distribute(m_Psi_ref);
    m_workspace_1=m_Psi_ref;
    
    system_rhs=0;
    system_matrix=0;

    FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_gradients|update_JxW_values|update_quadrature_points);
    //FEFaceValues<dim> fe_face_values ( fe, face_quadrature_formula, update_gradients|update_values|update_quadrature_points|update_normal_vectors|update_JxW_values);

    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    const unsigned int n_q_points    = quadrature_formula.size();
    //const unsigned int n_face_q_points = face_quadrature_formula.size();

    Vector<double> cell_rhs (dofs_per_cell);
    FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);

    vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    vector<Vector<double>> vals(n_q_points,Vector<double>(2));
    vector<vector<Tensor<1,dim>>> grads(n_q_points, vector<Tensor<1,dim>>(2));

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for ( ; cell!=endc; ++cell )
    {
      if( cell->is_locally_owned() )
      {
        cell_rhs=0;
        cell_matrix=0;

        fe_values.reinit (cell);
        fe_values.get_function_values(m_workspace_1, vals);
        fe_values.get_function_gradients(m_workspace_1, grads);

        for ( unsigned qp=0; qp<n_q_points; ++qp )
        {
          const double JxW = fe_values.JxW(qp)*fabs(fe_values.quadrature_point(qp)[1]);
          const double Q1 = Potential.value(fe_values.quadrature_point(qp)) - m_mu + m_gs*(vals[qp]*vals[qp]);

          for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            cell_rhs(i) += JxW*(grads[qp][0]*fe_values[rt].gradient(i,qp) + Q1*vals[qp][0]*fe_values[rt].value(i,qp) + grads[qp][1]*fe_values[it].gradient(i,qp) + Q1*vals[qp][1]*fe_values[it].value(i,qp));
	          for (unsigned int j=0; j<dofs_per_cell; ++j)
	            cell_matrix(i,j) += JxW*(fe_values[rt].value(i,qp)*fe_values[rt].value(j,qp)+fe_values[it].value(i,qp)*fe_values[it].value(j,qp));
          }
        }
        cell->get_dof_indices (local_dof_indices);
        constraints.distribute_local_to_global(cell_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs);
      }
    }
    system_rhs.compress(VectorOperation::add);   
    system_matrix.compress(VectorOperation::add);   

    solve();
    
    m_workspace_1=newton_update;
    VectorTools::integrate_difference ( dof_handler, m_workspace_1, ZeroFunction<dim>(2), m_error_per_cell, QGauss<dim>(fe.degree+2), VectorTools::L2_norm);    
    const double total_local_error = m_error_per_cell.l2_norm();
    err = std::sqrt (Utilities::MPI::sum (total_local_error * total_local_error, MPI_COMM_WORLD)); 
    
  }
  
  template <int dim>
  void MySolver<dim>::assemble_rhs ()
  {
    TimerOutput::Scope timing_section(m_computing_timer, "");

    const FEValuesExtractors::Scalar rt (0);
    const FEValuesExtractors::Scalar it (1);
    
    CPotential Potential( m_omega, m_QN1[2] );
    const QGauss<dim> quadrature_formula(fe.degree+1);
   
    system_rhs=0;
    
    constraints.distribute(m_Psi_ref);
    m_workspace_1=m_Psi_ref;

    FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_gradients|update_JxW_values|update_quadrature_points);

    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    const unsigned int n_q_points    = quadrature_formula.size();

    vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    
    Vector<double> cell_rhs (dofs_per_cell);
    vector<Vector<double>> vals(n_q_points,Vector<double>(2));
    vector<vector<Tensor<1,dim>>> grads(n_q_points, vector<Tensor<1,dim>>(2));
    
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for ( ; cell!=endc; ++cell )
    {
      if( cell->is_locally_owned() )
      {
        cell_rhs=0;

        fe_values.reinit (cell);
        fe_values.get_function_values(m_workspace_1, vals);
        fe_values.get_function_gradients(m_workspace_1, grads);

        for ( unsigned qp=0; qp<n_q_points; ++qp )
        {
          const double JxW = fe_values.JxW(qp)*fabs(fe_values.quadrature_point(qp)[1]);
          const double Q1 = Potential.value(fe_values.quadrature_point(qp)) - m_mu + m_gs*(vals[qp]*vals[qp]);

          for (unsigned int i=0; i<dofs_per_cell; ++i)
            cell_rhs(i) += JxW*(grads[qp][0]*fe_values[rt].gradient(i,qp) + Q1*vals[qp][0]*fe_values[rt].value(i,qp) 
	                       +grads[qp][1]*fe_values[it].gradient(i,qp) + Q1*vals[qp][1]*fe_values[it].value(i,qp));
        }
        cell->get_dof_indices (local_dof_indices);
        constraints.distribute_local_to_global(cell_rhs, local_dof_indices, system_rhs);
      }
    }
    system_rhs.compress(VectorOperation::add);   
    m_res = system_rhs.l2_norm();
    
  }

  template <int dim>
  void MySolver<dim>::assemble_system ()
  {
    TimerOutput::Scope timing_section(m_computing_timer, "");

    const FEValuesExtractors::Scalar rt (0);
    const FEValuesExtractors::Scalar it (1);
    
    CPotential Potential( m_omega, m_QN1[2] );
    const QGauss<dim> quadrature_formula(fe.degree+1);
    //const QGauss<dim-1> face_quadrature_formula(fe.degree+1);
       
    constraints.distribute(m_Psi_ref);
    m_workspace_1=m_Psi_ref;
    system_matrix=0;

    FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_gradients|update_JxW_values|update_quadrature_points);
    //FEFaceValues<dim> fe_face_values ( fe, face_quadrature_formula, update_gradients|update_values|update_quadrature_points|update_normal_vectors|update_JxW_values);

    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    const unsigned int n_q_points    = quadrature_formula.size();
    //const unsigned int n_face_q_points = face_quadrature_formula.size();

    FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);

    vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    vector<vector<Tensor<1,dim>>> Psi_ref_grad(n_q_points, vector<Tensor<1,dim>>(2));
    vector<Vector<double>> Psi_ref(n_q_points,Vector<double>(2));

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for ( ; cell!=endc; ++cell )
    {
      if( cell->is_locally_owned() )
      {
        cell_matrix=0;

        fe_values.reinit (cell);
        fe_values.get_function_values(m_workspace_1, Psi_ref);
        fe_values.get_function_gradients(m_workspace_1, Psi_ref_grad);

        for ( unsigned qp=0; qp<n_q_points; ++qp )
        {
          const double JxW = fe_values.JxW(qp)*fabs(fe_values.quadrature_point(qp)[1]);
	        const double fak = m_gs*Psi_ref[qp][0]*Psi_ref[qp][1];
	        const double Pot = Potential.value(fe_values.quadrature_point(qp)) - m_mu;
	        const double req = Psi_ref[qp][0]*Psi_ref[qp][0];
	        const double imq = Psi_ref[qp][1]*Psi_ref[qp][1];

          for (unsigned int i=0; i<dofs_per_cell; ++i)
            for (unsigned int j=0; j<dofs_per_cell; ++j)
              cell_matrix(i,j) += JxW*(fe_values[rt].gradient(i,qp)*fe_values[rt].gradient(j,qp) + (Pot+m_gs*(3*req+imq))*fe_values[rt].value(i,qp)*fe_values[rt].value(j,qp) 
	                              +fak*fe_values[rt].value(i,qp)*fe_values[it].value(j,qp) 
	                              +fak*fe_values[it].value(i,qp)*fe_values[rt].value(j,qp)
				      +fe_values[it].gradient(i,qp)*fe_values[it].gradient(j,qp) + (Pot+m_gs*(3*imq+req))*fe_values[it].value(i,qp)*fe_values[it].value(j,qp));
        }
        cell->get_dof_indices (local_dof_indices);
        constraints.distribute_local_to_global(cell_matrix, local_dof_indices, system_matrix);
      }
    }
    system_matrix.compress(VectorOperation::add);
    
  }
  
  template <int dim>
  void MySolver<dim>::compute_J( double& retval )
  {
    TimerOutput::Scope timing_section(m_computing_timer, "");

    do_superposition();
    constraints.distribute(m_Psi_ref);
    m_workspace_1 = m_Psi_ref;
    
    CPotential Potential( m_omega, m_QN1[2] );
    const QGauss<dim>  quadrature_formula(fe.degree+1);
    FEValues<dim> fe_values (fe, quadrature_formula, update_gradients|update_values|update_JxW_values|update_quadrature_points);

    const unsigned int n_q_points = quadrature_formula.size();
    vector<Vector<double>> vals(n_q_points,Vector<double>(2));
    vector<vector<Tensor<1,dim>>> grads(n_q_points, vector<Tensor<1,dim>>(2));

    double psum=0;
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for ( ; cell!=endc; ++cell )
    {
      if( cell->is_locally_owned() )
      {
        fe_values.reinit (cell);
        fe_values.get_function_values(m_workspace_1, vals);
        fe_values.get_function_gradients(m_workspace_1, grads);

	      for ( unsigned qp=0; qp<n_q_points; ++qp )
        {
          const double uq = vals[qp]*vals[qp];
	        psum += fe_values.JxW(qp)*fabs(fe_values.quadrature_point(qp)[1])*( grads[qp][0]*grads[qp][0] + grads[qp][1]*grads[qp][1] + (Potential.value(fe_values.quadrature_point(qp))-m_mu)*uq + 0.5*m_gs*uq*uq ); 
        }
      }
    }
    psum *= 0.5;
    MPI_Allreduce( &psum, &retval, 1, MPI_DOUBLE, MPI_SUM, mpi_communicator);
    
  }
  
  template<int dim>
  void MySolver<dim>::find_min_J()
  {
    size_t iter=0;
    int status;
    double size;

    gsl_multimin_function my_func;
    my_func.n = 2;
    my_func.f = fun_neu<dim>;
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
      status = gsl_multimin_test_size(size,1e-6);
    }
    while (status == GSL_CONTINUE && iter < 1000);

    m_t[0] = gsl_vector_get (s->x, 0);
    m_t[1] = gsl_vector_get (s->x, 1);

    gsl_multimin_fminimizer_free (s);
    gsl_vector_free(ss);
    gsl_vector_free (x);
  }    
  
  template <int dim>
  void MySolver<dim>::solve ()
  {
    TimerOutput::Scope timing_section(m_computing_timer, "");
    pcout << "Solving..." << endl;
    
    SolverControl solver_control;
    
    PETScWrappers::SparseDirectMUMPS solver(solver_control, mpi_communicator);
    solver.set_symmetric_mode(false);
    solver.solve(system_matrix, newton_update, system_rhs);
    constraints.distribute (newton_update);
    
  }
  
  template <int dim>
  void MySolver<dim>::make_grid_custom ()
  {
    TimerOutput::Scope timing_section(m_computing_timer, "");
    Point<dim,double> pt1(0,0); 
    Point<dim,double> pt2(m_xmax,m_ymax);
    
    CPotential Potential_fct ( m_omega, m_QN1[2] );
    
    GridGenerator::hyper_rectangle(triangulation, pt2, pt1);
    triangulation.refine_global(6);
    
    //double isovalues[] = {60,50,40};
    double isovalues[] = {60,50};

    for( unsigned int step=0; step<sizeof(isovalues)/sizeof(double); step++ )
    {
      typename parallel::distributed::Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active(), endc = triangulation.end();
      for( ; cell!=endc; ++cell )
        for (unsigned int v=0; v < GeometryInfo<dim>::vertices_per_cell; ++v )
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
    
    typename parallel::distributed::Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active(), endc = triangulation.end();
    for(; cell!=endc; ++cell)
    {
      for (unsigned int f=0; f < GeometryInfo<dim>::faces_per_cell; ++f)
      {
      	const Point<dim> face_center = cell->face(f)->center();
      	if (cell->face(f)->at_boundary() && !(face_center[1]==0) )    
      	{
	        cell->face(f)->set_all_boundary_ids(1);
	      }
      }
    }

    unsigned int tmp1[2], tmp2[2];
    tmp1[0] = triangulation.n_cells();
    tmp1[1] = triangulation.n_active_cells();

    MPI_Allreduce( tmp1, tmp2, 2, MPI_UNSIGNED, MPI_SUM, mpi_communicator);

/*
    GridOutFlags::Msh opt(true, true);
    string filename = "grid-" + Utilities::int_to_string(triangulation.locally_owned_subdomain(), 4);
    ofstream output ((filename + ".msh").c_str());
    GridOut grid_out;
    grid_out.set_flags(opt);
    grid_out.write_msh (triangulation,output);
*/
    m_total_no_cells = tmp2[0];
    m_total_no_active_cells = tmp2[1];
    
  }
  
  template<int dim>
  void MySolver<dim>::do_superposition()
  {
    m_Psi_ref=0;
    m_Psi_ref.add(m_t[0],m_Psi_1,m_t[1],m_Psi_2);
    constraints.distribute (m_Psi_ref);
  }

  template <int dim>
  void MySolver<dim>::setup_system( const bool initial_step )
  {
    TimerOutput::Scope timing_section(m_computing_timer, "");

    if( initial_step )
    {
      dof_handler.distribute_dofs (fe);
      
      locally_owned_dofs = dof_handler.locally_owned_dofs ();
      DoFTools::extract_locally_relevant_dofs (dof_handler, locally_relevant_dofs);

      m_Psi_1.reinit (locally_owned_dofs, mpi_communicator);
      m_Psi_2.reinit (locally_owned_dofs, mpi_communicator);
    }

    m_Psi_ref.reinit (locally_owned_dofs, mpi_communicator);
    newton_update.reinit (locally_owned_dofs, mpi_communicator);
    system_rhs.reinit(locally_owned_dofs, mpi_communicator);
    m_workspace_ng.reinit (locally_owned_dofs, mpi_communicator);
    m_workspace_1.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    m_workspace_2.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    m_error_per_cell.reinit(triangulation.n_active_cells());
    
    int myrank;
    MPI_Comm_rank( mpi_communicator, &myrank );    
    cout << "(" << myrank << ") locally_owned_dofs = " << system_rhs.local_size()  << endl;
    
    system_rhs = 0;

    vector<bool> mask (dof_handler.get_fe().n_components(), true );
    
    constraints.clear ();
    constraints.reinit (locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints (dof_handler, constraints);
    VectorTools::interpolate_boundary_values (dof_handler, 1, ZeroFunction<dim>(), constraints, ComponentMask(mask));
    if( m_QN1[2] > 0 )
      VectorTools::interpolate_boundary_values (dof_handler, 0, ZeroFunction<dim>(), constraints, ComponentMask(mask));
    constraints.close ();

    DynamicSparsityPattern dsp (locally_relevant_dofs);

    DoFTools::make_sparsity_pattern (dof_handler, dsp, constraints, false);
    SparsityTools::distribute_sparsity_pattern (dsp, dof_handler.n_locally_owned_dofs_per_processor(), mpi_communicator, locally_relevant_dofs);

    system_matrix.reinit (locally_owned_dofs, locally_owned_dofs, dsp, mpi_communicator);    
    
  }

  template <int dim>
  void MySolver<dim>::output_guess ()
  {
    TimerOutput::Scope timing_section(m_computing_timer, "");
    
    constraints.distribute(m_Psi_1);
    constraints.distribute(m_Psi_2);
    m_workspace_1=m_Psi_1;
    m_workspace_2=m_Psi_2;
    
    CPotential Potential_fct ( m_omega, m_QN1[2] );
    VectorTools::interpolate (dof_handler, Potential_fct, m_workspace_ng );
    constraints.distribute(m_workspace_ng);
    m_Psi_ref=m_workspace_ng;
    
    DataOut<dim> data_out;
    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (m_workspace_1, "m_Psi_1");
    data_out.add_data_vector (m_workspace_2, "m_Psi_2");
    data_out.add_data_vector (m_Psi_ref, "m_Potential");
    data_out.build_patches ();
    data_out.write_vtu_in_parallel ("guess.vtu",mpi_communicator);

    
  }

  template <int dim>
  void MySolver<dim>::output_results ( string path, string prefix )
  {
    TimerOutput::Scope timing_section(m_computing_timer, "");

    Vector<float> subdomain (triangulation.n_active_cells());
    for (unsigned int i=0; i<subdomain.size(); ++i)
      subdomain(i) = triangulation.locally_owned_subdomain();

    m_workspace_1=m_Psi_ref;
    
    DataOut<dim> data_out;
    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (m_workspace_1, "Psi_sol");
    data_out.add_data_vector (m_error_per_cell, "error per cell");
    data_out.add_data_vector (subdomain, "subdomain");
    data_out.build_patches ();

    std::string filename = path + prefix + "-" + Utilities::int_to_string (m_counter,5) + ".vtu";
    data_out.write_vtu_in_parallel (filename.c_str(), mpi_communicator);

        
  }

  template <int dim>
  int MySolver<dim>::DoIter ( string path )
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
      pcout << "--------------------------------------------------------------------------------" << endl;
      pcout << "- " << path << " - " << m_counter << endl;

      assemble_system();
      solve();

      m_Psi_2.add( -1e-3*m_t[1]/fabs(m_t[1]), newton_update); 
      constraints.distribute(m_Psi_2);

      find_min_J();
      
      do_superposition();
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

      if( m_root ) cout << m_table;
      if( m_res < m_epsilon[0] ) { retval=Status::SUCCESS; break; }
      if( l2norm_t() < 1e-4 ) { retval=Status::ZERO_SOL; break; }

      m_counter++;
    }while( true );

    
    // Standard Newton 
    do
    {
      pcout << "--------------------------------------------------------------------------------" << endl;
      pcout << "-- " << path << " - " << m_counter << endl;

      assemble_system();
      solve();

      m_Psi_ref.add( -0.1, newton_update); 
      constraints.distribute(m_Psi_ref);

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

      if( m_root ) cout << m_table;
      if( m_res < m_epsilon[1] ) { retval=Status::SUCCESS; break; }

      m_counter++;
    }while( true );
    
    m_N = Particle_Number(m_Psi_ref);
    
    if( m_N < 1e-5 ) retval = Status::ZERO_SOL;
    
    string filename = path + "log.csv";
    if( m_root ) m_table.dump_2_file(filename);
    
    return retval;
  }

  template <int dim>
  void MySolver<dim>::run()
  {
    int status;
    string path;
    char shellcmd[255];
    double T, N, W;

    make_grid_custom();
    setup_system(true);

    CEigenfunctions Ef1( m_QN1, m_omega );

    VectorTools::interpolate (dof_handler, Ef1, m_Psi_1 );
    m_Psi_1 *= 1.0/sqrt(Particle_Number(m_Psi_1));
    m_Psi_2 = 0;
      
    compute_E_lin( m_Psi_1, T, N, W );
    output_guess();

    m_mu = T/N+m_gs/fabs(m_gs)*m_dmu;

    pcout << "T = " << T << endl;
    pcout << "N = " << N << endl;
    pcout << "W = " << W << endl;
    pcout << "m_mu = " << m_mu << endl;

    status = DoIter();
    
    pcout << "L2_norm of m_Psi_ref: " << m_N << endl;

    if( status == Status::SUCCESS )
    {
      output_results("","final");
      dump_info_xml();
    }

    if( m_root )
    {
      ofstream ofs("log.txt");
      ofs << m_table;
    }
  }

  template <int dim>
  void MySolver<dim>::run2b ()
  {
    string path;
    char shellcmd[255];
    double T, N, W;
    int status;

    make_grid_custom();
    setup_system(true);

    CEigenfunctions Ef1( m_QN1, m_omega );
    VectorTools::interpolate (dof_handler, Ef1, m_Psi_1 );

    m_Psi_1 *= 1.0/sqrt(Particle_Number(m_Psi_1));
    //m_Psi_2 = m_Psi_1;
    m_Psi_2 = 0; 

    compute_E_lin( m_Psi_1, T, N, W );
    double m_mu_0 = T/N;
    m_mu = ceil(10.0*m_mu_0)/10.0 + m_gs/fabs(m_gs)*m_dmu;
    
    output_guess();
    m_results.clear();
    for( int i=0; i<m_Ndmu; ++i )
    {
      sprintf( shellcmd, "mkdir %.4d/", i );
      if( m_root ) system(shellcmd);
      sprintf( shellcmd, "%.4d/", i );
      path = shellcmd;
      
      // nehari
      //m_ti = sqrt((m_mu*N-T)/(4.0*m_gs*W)); // if m_Psi_2 == m_Psi_1
      m_ti = sqrt(2*(m_mu*N-T)/(m_gs*W));

      pcout << "T = " << T << endl;
      pcout << "N = " << N << endl;
      pcout << "W = " << W << endl;
      pcout << "m_mu = " << m_mu << endl;      
      pcout << "m_ti = " << m_ti << endl;      
      
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
    if( m_root ) m_results.dump_2_file( "results.csv" );
  }
  
  template<int dim>
  void MySolver<dim>::save( string filename )
  {
    m_workspace_1=m_Psi_ref;
    parallel::distributed::SolutionTransfer<dim,LA::MPI::Vector> solution_transfer(dof_handler);
    solution_transfer.prepare_for_serialization(m_workspace_1);
    triangulation.save( filename.c_str() );
  }
  
  template<int dim>
  void MySolver<dim>::save_one( string filename )
  {
    double tmp = Particle_Number(m_Psi_ref);
    m_workspace_1=m_Psi_ref;
    m_workspace_1*=sqrt(1/tmp);
    parallel::distributed::SolutionTransfer<dim,LA::MPI::Vector> solution_transfer(dof_handler);
    solution_transfer.prepare_for_serialization(m_workspace_1);
    triangulation.save( filename.c_str() );
  }
  
  template <int dim>
  void MySolver<dim>::dump_info_xml ( const string path )
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
    fs2 << L"<YMIN>" << m_ymin << L"</YMIN>\n";
    fs2 << L"<YMAX>" << m_ymax << L"</YMAX>\n";
    fs2 << L"<ZMIN>" << m_zmin << L"</ZMIN>\n";
    fs2 << L"<ZMAX>" << m_zmax << L"</ZMAX>\n";
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
    BreedSolver::MySolver<2> solver("params.xml");
    solver.run2b ();
  }
return EXIT_SUCCESS;
}
