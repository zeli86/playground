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

/** Želimir Marojević
 */

#include <deal.II/lac/generic_linear_algebra.h>

namespace LA
{
  using namespace dealii::LinearAlgebraPETSc;
}

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/petsc_parallel_sparse_matrix.h>
#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_precondition.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/function_parser.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/solution_transfer.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>

#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <locale>
#include <limits>
#include <cmath>
#include <array>
#include <unistd.h>

#include "global.h"
#include "mpi.h"
#include "functions.h"
#include "MyParameterHandler.h"
#include "MyRealTools.h"
#include "my_table.h"

#define STR1(x) #x
#define STR2(x) STR1(x)

namespace BreedSolver
{
  using namespace std;
  using namespace dealii;

  enum Status { SUCCESS, FAILED, ZERO_SOL, SLOW_CONV };

  template <int dim>
  class MySolver
  {
  public:
    MySolver( const std::string );
    virtual ~MySolver();

    void run ( const std::string& );

    MPI_Comm mpi_communicator;
  protected:
 
    double m_xmin, m_xmax;
    double m_ymin, m_ymax;
    double m_zmin, m_zmax;
    double m_res, m_res_old, m_resp;
    double m_final_error;
    double m_N;
    double m_mu_0;
    double m_mu_0_punkt;
    double m_mu;
    double m_delta_mu;    
    double m_delta_s;
    double m_gs;
    vector<double> m_omega;
    vector<double> m_epsilon;

    bool m_root;
    int m_rank;
    
    unsigned m_counter;
    unsigned m_global_refinement;
    unsigned m_total_no_cells;
    unsigned m_total_no_active_cells;    
    unsigned m_NA;
    
    int DoIter( string="" );
    
    void make_grid();
    void make_grid_custom();
    void setup_system();
    void assemble_rhs();
    void assemble_system();
    void compute_mu();
    void save( string );
    void load( const std::string& );
    void dump_info_xml( const string );
    void compute_tangent();
    void compute_newton_direction_natural();
    void compute_newton_direction_arclength();
    
    double Particle_Number( LA::MPI::Vector& );
    
    void solve();
    void compute_E_lin( LA::MPI::Vector&, double&, double&, double& );
    void estimate_error( double& );

    void output_results ( string, string = "step" );
    void output_guess ();

    ofstream m_computing_timer_log;
    TimerOutput m_computing_timer;    
    MyParameterHandler m_oPropertyHandler;
    ConditionalOStream pcout;
    parallel::distributed::Triangulation<dim> triangulation;
    FE_Q<dim> fe;
    DoFHandler<dim> dof_handler;
    FESystem<dim> fe_2;
    DoFHandler<dim> dof_handler_2;
    IndexSet locally_owned_dofs;
    IndexSet locally_relevant_dofs;
    IndexSet locally_owned_dofs_2;
    IndexSet locally_relevant_dofs_2;
    ConstraintMatrix constraints;
    ConstraintMatrix constraints_2;

    LA::MPI::SparseMatrix m_system_matrix, m_system_matrix_2;
    LA::MPI::Vector m_system_rhs, m_system_rhs_2;
    LA::MPI::Vector m_Psi;
    LA::MPI::Vector m_Psi_C_ghosted;
    LA::MPI::Vector m_workspace_1;
    LA::MPI::Vector m_workspace_ng;
    LA::MPI::Vector m_y;
    LA::MPI::Vector m_z;
    LA::MPI::Vector m_Psi_0;
    LA::MPI::Vector m_Psi_0_punkt;
    LA::MPI::Vector m_newton_direction;
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
    m_computing_timer(mpi_communicator, m_computing_timer_log, TimerOutput::summary, TimerOutput:: cpu_and_wall_times ), 
    m_ph(xmlfilename),
    pcout (cout, (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)),
    triangulation (mpi_communicator, typename Triangulation<dim>::MeshSmoothing(Triangulation<dim>::limit_level_difference_at_vertices|Triangulation<dim>::eliminate_refined_inner_islands|Triangulation<dim>::smoothing_on_refinement|Triangulation<dim>::smoothing_on_coarsening)),
    fe (gl_degree_fe),
    dof_handler (triangulation),
    fe_2 (FE_Q<dim>(gl_degree_fe), 2),
    dof_handler_2 (triangulation)    
  {
    try
    {
      m_omega = m_ph.Get_Physics("omega");
      m_gs = m_ph.Get_Physics("gs_1",0);

      m_xmin = m_ph.Get_Mesh("xrange",0);
      m_xmax = m_ph.Get_Mesh("xrange",1);
      m_ymin = m_ph.Get_Mesh("yrange",0);
      m_ymax = m_ph.Get_Mesh("yrange",1);
      m_zmin = m_ph.Get_Mesh("zrange",0);
      m_zmax = m_ph.Get_Mesh("zrange",1);    

      m_NA = int(m_ph.Get_Algorithm("NA",0));
      m_epsilon = m_ph.Get_Algorithm("epsilon"); 
      m_guess_str = m_ph.Get_Parameter( "guess_fct" );
    }
    catch( const std::string info )
    {
      std::cerr << info << endl;
      MPI_Abort( mpi_communicator, 0 );
    }    
    
    m_root = (Utilities::MPI::this_mpi_process(mpi_communicator) == 0);
    MPI_Comm_rank(mpi_communicator, &m_rank);

    m_counter = 0;
    m_final_error = 0;

    m_delta_s = 0;    
  }

  template <int dim>
  MySolver<dim>::~MySolver ()
  {
    dof_handler.clear();
    dof_handler_2.clear();
  }

  template <int dim>
  void MySolver<dim>::compute_tangent ()
  {
    TimerOutput::Scope timing_section(m_computing_timer, "");

    constraints.distribute(m_Psi);
    m_workspace_1=m_Psi;
    
    CPotential<dim> Potential( m_omega );
    MyRealTools::MPI::AssembleSystem_tangent( dof_handler, fe, constraints, m_workspace_1, Potential, m_mu, m_gs, m_system_matrix, m_system_rhs );

    pcout << "Solving..." << endl;
    SolverControl solver_control;
    PETScWrappers::SparseDirectMUMPS solver(solver_control, mpi_communicator);
    solver.set_symmetric_mode(false);
    solver.solve(m_system_matrix, m_Psi_0_punkt, m_system_rhs);
    constraints.distribute (m_Psi_0_punkt);    

    double N = Particle_Number( );
    m_mu_0_punkt =  1/sqrt(1+N); // +/- 1 / sqrt(1+N) 

    m_Psi_0_punkt *= m_mu_0_punkt;
    
    
  }

  template <int dim>
  void MySolver<dim>::compute_newton_direction_natural ()
  {
    TimerOutput::Scope timing_section(m_computing_timer, "");

    CPotential<dim> Potential( m_omega );
    const QGauss<dim> quadrature_formula(fe.degree+1);
    
    constraints.distribute(m_Psi);
    m_workspace_1=m_Psi;
    
    m_system_matrix=0;
    m_system_rhs=0;

    FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_gradients|update_JxW_values|update_quadrature_points);

    const unsigned dofs_per_cell = fe.dofs_per_cell;
    const unsigned n_q_points = quadrature_formula.size();

    Vector<double> cell_rhs (dofs_per_cell);
    FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);

    vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    vector<Tensor<1, dim> > Psi_grad(n_q_points);
    vector<double> Psi(n_q_points);

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for ( ; cell!=endc; ++cell )
    {
      if( cell->is_locally_owned() )
      {
        cell_rhs = 0;
        cell_matrix = 0;

        fe_values.reinit (cell);
        fe_values.get_function_values(m_workspace_1, Psi);
        fe_values.get_function_gradients(m_workspace_1, Psi_grad);

        for ( unsigned qp=0; qp<n_q_points; ++qp )
        {
          double JxW = fe_values.JxW(qp);
          double pq = m_gs*Psi[qp]*Psi[qp];
          double pot = Potential.value(fe_values.quadrature_point(qp)); 
          double Q1 = pot + pq;
          double Q2 = pot + 3.0*pq;

          for ( unsigned i=0; i<dofs_per_cell; ++i )
          {
            cell_rhs(i) -= JxW*(Psi_grad[qp]*fe_values.shape_grad(i,qp) + Q1*Psi[qp]*fe_values.shape_value(i,qp));
            for ( unsigned j=0; j<dofs_per_cell; ++j )
              cell_matrix(i,j) += JxW*(fe_values.shape_grad(i,qp)*fe_values.shape_grad(j,qp) + Q2*fe_values.shape_value(i,qp)*fe_values.shape_value(j,qp));
          }
        }
        cell->get_dof_indices (local_dof_indices);
        constraints.distribute_local_to_global(cell_matrix, cell_rhs, local_dof_indices, m_system_matrix, m_system_rhs);
      }
    }
    m_system_rhs.compress(VectorOperation::add);
    m_system_matrix.compress(VectorOperation::add);

    pcout << "Solving..." << endl;
    SolverControl solver_control;
    PETScWrappers::SparseDirectMUMPS solver(solver_control, mpi_communicator);
    solver.set_symmetric_mode(false);
    solver.solve(m_system_matrix, m_newton_direction, m_system_rhs);
    constraints.distribute (m_newton_direction);    
    
    
  }  

  template <int dim>
  void MySolver<dim>::compute_newton_direction_arclength ()
  {
    TimerOutput::Scope timing_section(m_computing_timer, "");

    CPotential<dim> Potential( m_omega );
    const QGauss<dim> quadrature_formula(fe.degree+1);
    
    constraints.distribute(m_Psi);
    m_workspace_1=m_Psi;
    
    m_system_matrix=0;
    m_system_rhs=0;

    FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_gradients|update_JxW_values|update_quadrature_points);

    const unsigned dofs_per_cell = fe.dofs_per_cell;
    const unsigned n_q_points = quadrature_formula.size();

    Vector<double> cell_rhs (dofs_per_cell);
    Vector<double> cell_rhs_2 (dofs_per_cell);
    FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);

    vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    vector<Tensor<1, dim> > Psi_grad(n_q_points);
    vector<double> Psi(n_q_points);

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for ( ; cell!=endc; ++cell )
    {
      if( cell->is_locally_owned() )
      {
        cell_rhs = 0;
        cell_rhs_2 = 0;
        cell_matrix = 0;

        fe_values.reinit (cell);
        fe_values.get_function_values(m_workspace_1, Psi);
        fe_values.get_function_gradients(m_workspace_1, Psi_grad);

        for ( unsigned qp=0; qp<n_q_points; ++qp )
        {
          double JxW = fe_values.JxW(qp);
          double pq = m_gs*Psi[qp]*Psi[qp];
          double pot = Potential.value(fe_values.quadrature_point(qp)); 
          double Q1 = pot + pq;
          double Q2 = pot + 3.0*pq;

          for ( unsigned i=0; i<dofs_per_cell; ++i )
          {
            cell_rhs(i) -= JxW*m_mu*Psi[qp]*fe_values.shape_value(i);
            cell_rhs_2(i) -= JxW*(Psi_grad[qp]*fe_values.shape_grad(i,qp) + Q1*Psi[qp]*fe_values.shape_value(i,qp));
            for ( unsigned j=0; j<dofs_per_cell; ++j )
              cell_matrix(i,j) += JxW*(fe_values.shape_grad(i,qp)*fe_values.shape_grad(j,qp) + Q2*fe_values.shape_value(i,qp)*fe_values.shape_value(j,qp));
          }
        }
        cell->get_dof_indices (local_dof_indices);
        constraints.distribute_local_to_global(cell_matrix, cell_rhs, local_dof_indices, m_system_matrix, m_system_rhs);
        constraints.distribute_local_to_global(cell_rhs_2, local_dof_indices, m_system_rhs);
      }
    }
    m_system_rhs.compress(VectorOperation::add);
    m_system_rhs_2.compress(VectorOperation::add);
    m_system_matrix.compress(VectorOperation::add);

    m_system_matrix.copy_from( m_system_matrix );

    pcout << "Solving m_y ..." << endl;
    SolverControl solver_control;
    PETScWrappers::SparseDirectMUMPS solver(solver_control, mpi_communicator);
    solver.set_symmetric_mode(false);
    solver.solve(m_system_matrix, m_y, m_system_rhs);
    constraints.distribute (m_y);    

    pcout << "Solving m_z ..." << endl;
    solver.solve(m_system_matrix_2, m_z, m_system_rhs_2);
    constraints.distribute (m_z);

    m_workspace_ng = m_Psi;
    m_workspace_ng -= m_Psi_0;
    double N = m_Psi_0_punkt * m_workspace_ng + m_mu_0_punkt * (m_mu - m_mu_0) + m_delta_s;
    m_delta_mu = (-N-m_Psi_0_punkt*m_z)/(m_mu_0_punkt-m_Psi_0*m_y);

    m_y *= m_delta_mu;
    m_newton_direction = m_z;
    m_newton_direction -= m_y;
    
    
  }    

  template <int dim>
  void MySolver<dim>::compute_E_lin( LA::MPI::Vector& vec, double& T, double& N, double& W )
  {
    TimerOutput::Scope timing_section(m_computing_timer, "");
    
    constraints.distribute(vec);
    m_workspace_1 = vec;
    
    CPotential<dim> Potential( m_omega );
    const QGauss<dim>  quadrature_formula(fe.degree+1);
    FEValues<dim> fe_values (fe, quadrature_formula, update_gradients|update_values|update_JxW_values|update_quadrature_points);

    const unsigned int n_q_points = quadrature_formula.size();
    vector<double> vals(n_q_points);    
    vector<Tensor<1,dim>> grad(n_q_points);
    
    double JxW, vec_val_q;
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
          vec_val_q = vals[qp]*vals[qp];
          JxW = fe_values.JxW(qp);
          tmp[0] += JxW*( grad[qp]*grad[qp] + Potential.value(fe_values.quadrature_point(qp))*vec_val_q );
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

    const unsigned n_q_points = quadrature_formula.size();
    vector<double> vals(n_q_points);

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
      if( cell->is_locally_owned() )
      {
        fe_values.reinit (cell);
        fe_values.get_function_values( m_workspace_1, vals );      

        for ( unsigned qp=0; qp<n_q_points; ++qp )
          tmp1 += fe_values.JxW(qp)*(vals[qp]*vals[qp]);
      }
    }
    MPI_Allreduce( &tmp1, &retval, 1, MPI_DOUBLE, MPI_SUM, mpi_communicator);
    
  return retval;
  }
  
 
  template <int dim>
  void MySolver<dim>::estimate_error ( double& err )
  {
    TimerOutput::Scope timing_section(m_computing_timer, "");
    
    compute_mu();
    
    CPotential<dim> Potential( m_omega );
    const QGauss<dim> quadrature_formula(fe.degree+1);
   
    constraints.distribute(m_Psi);
    m_workspace_1=m_Psi;
    
    m_system_rhs=0;
    m_system_matrix=0;

    FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_gradients|update_JxW_values|update_quadrature_points);

    const unsigned dofs_per_cell = fe.dofs_per_cell;
    const unsigned n_q_points = quadrature_formula.size();

    Vector<double> cell_rhs (dofs_per_cell);
    FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);

    vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    vector<double> vals(n_q_points);
    vector<Tensor<1,dim>> grads(n_q_points);

    double JxW, Q1, pq, tmp;
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
          JxW = fe_values.JxW(qp);
          Q1 = Potential.value(fe_values.quadrature_point(qp)) - m_mu + m_gs*(vals[qp]*vals[qp]);

          for ( unsigned i=0; i<dofs_per_cell; ++i )
          {
            cell_rhs(i) += JxW*(grads[qp]*fe_values.shape_grad(i,qp) + Q1*vals[qp]*fe_values.shape_value(i,qp));
            for ( unsigned j=0; j<dofs_per_cell; ++j )
              cell_matrix(i,j) += JxW*(fe_values.shape_value(i,qp)*fe_values.shape_value(j,qp));
          }
        }
        cell->get_dof_indices (local_dof_indices);
        constraints.distribute_local_to_global(cell_matrix, cell_rhs, local_dof_indices, m_system_matrix, m_system_rhs);
      }
    }
    m_system_rhs.compress(VectorOperation::add);   
    m_system_matrix.compress(VectorOperation::add);   

    solve();
    
    m_workspace_1=m_newton_direction;
    VectorTools::integrate_difference ( dof_handler, m_workspace_1, ZeroFunction<dim>(2), m_error_per_cell, QGauss<dim>(fe.degree+2), VectorTools::L2_norm);    
    const double total_local_error = m_error_per_cell.l2_norm();
    err = std::sqrt (Utilities::MPI::sum (total_local_error * total_local_error, MPI_COMM_WORLD)); 
    
  }
  
  template <int dim>
  void MySolver<dim>::assemble_rhs ()
  {
    TimerOutput::Scope timing_section(m_computing_timer, "");

    CPotential<dim> Potential( m_omega );
    const QGauss<dim> quadrature_formula(fe.degree+1);

    constraints.distribute(m_Psi);
    m_workspace_1=m_Psi;
    
    m_system_rhs=0;

    FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_gradients|update_JxW_values|update_quadrature_points);

    const unsigned dofs_per_cell = fe.dofs_per_cell;
    const unsigned n_q_points = quadrature_formula.size();

    vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    
    Vector<double> cell_rhs (dofs_per_cell);
    vector<double> vals(n_q_points);
    vector<Tensor<1,dim>> grads(n_q_points);
    
    double JxW, Q1;
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
          JxW = fe_values.JxW(qp);
          Q1 = Potential.value(fe_values.quadrature_point(qp)) + m_gs*(vals[qp]*vals[qp]); 

          for (unsigned int i=0; i<dofs_per_cell; ++i)
            cell_rhs(i) += JxW*(grads[qp]*fe_values.shape_grad(i,qp) + Q1*vals[qp]*fe_values.shape_value(i,qp));
        }
        cell->get_dof_indices (local_dof_indices);
        constraints.distribute_local_to_global(cell_rhs, local_dof_indices, m_system_rhs);
      }
    }
    m_system_rhs.compress(VectorOperation::add);   
    
  }

  
  template<int dim>
  void MySolver<dim>::load( const std::string& filename )
  {
    make_grid();
    triangulation.load( filename.c_str() );

    setup_system();

    parallel::distributed::SolutionTransfer<dim,LA::MPI::Vector> solution_transfer(dof_handler);
    solution_transfer.deserialize(m_Psi_0);
    constraints.distribute(m_Psi_0);
  }    

  template <int dim>
  void MySolver<dim>::assemble_system ()
  {
    TimerOutput::Scope timing_section(m_computing_timer, "");

    CPotential<dim> Potential( m_omega );
    const QGauss<dim> quadrature_formula(fe.degree+1);
    
    constraints.distribute(m_Psi);
    m_workspace_1=m_Psi;    
    
    m_system_matrix=0;

    FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_gradients|update_JxW_values|update_quadrature_points);

    const unsigned dofs_per_cell = fe.dofs_per_cell;
    const unsigned n_q_points = quadrature_formula.size();

    FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);
    vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    
    double JxW, Q1;
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for ( ; cell!=endc; ++cell )
    {
      if( cell->is_locally_owned() )
      {
        cell_matrix=0;
        fe_values.reinit (cell);

        for ( unsigned qp=0; qp<n_q_points; ++qp )
        {
          JxW = fe_values.JxW(qp);
	        //Q1 = Potential.value(fe_values.quadrature_point(qp)); 

          for (unsigned int i=0; i<dofs_per_cell; ++i)
            for (unsigned int j=0; j<dofs_per_cell; ++j)
              cell_matrix(i,j) += JxW*(fe_values.shape_grad(i,qp)*fe_values.shape_grad(j,qp) + fe_values.shape_value(i,qp)*fe_values.shape_value(j,qp));
        }
        cell->get_dof_indices (local_dof_indices);
        constraints.distribute_local_to_global(cell_matrix, local_dof_indices, m_system_matrix);
      }
    }
    m_system_matrix.compress(VectorOperation::add);
    
  }
   
  
  template <int dim>
  void MySolver<dim>::compute_mu()
  {
    TimerOutput::Scope timing_section(m_computing_timer, "");

    CPotential<dim> Potential( m_omega );
    const QGauss<dim>  quadrature_formula(fe.degree+1);
    FEValues<dim> fe_values (fe, quadrature_formula, update_gradients|update_values|update_JxW_values|update_quadrature_points);

    const unsigned n_q_points = quadrature_formula.size();
    vector<double> vals(n_q_points);
    vector<Tensor<1,dim>> grads(n_q_points);
    
    constraints.distribute(m_Psi);
    m_workspace_1=m_Psi;   
    
    double JxW, psum=0;
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
          double uq = vals[qp]*vals[qp];
          psum += fe_values.JxW(qp)*(grads[qp]*grads[qp] + (Potential.value(fe_values.quadrature_point(qp)) + m_gs*uq)*uq); 
        }
      }
    }    
    MPI_Allreduce( &psum, &m_mu, 1, MPI_DOUBLE, MPI_SUM, mpi_communicator);
    
  }
  
  
  template <int dim>
  void MySolver<dim>::solve ()
  {
    TimerOutput::Scope timing_section(m_computing_timer, "");
    pcout << "Solving..." << endl;
    
    SolverControl solver_control;
    
    PETScWrappers::SparseDirectMUMPS solver(solver_control, mpi_communicator);
    solver.set_symmetric_mode(false);
    solver.solve(m_system_matrix, m_newton_direction, m_system_rhs);
    constraints.distribute (m_newton_direction);
    
  }
  
  /*
  template <int dim>
  void MySolver<dim>::solve ()
  {
    m_sob_grad=0;    
    SolverControl solver_control ( m_sob_grad.size(), 1e-15 );
    PETScWrappers::SolverCG cg (solver_control, mpi_communicator);

    PETScWrappers::PreconditionSOR preconditioner;
    PETScWrappers::PreconditionSOR::AdditionalData data;
    preconditioner.initialize(m_system_matrix, data);    
    cg.solve (m_system_matrix, m_sob_grad, m_system_rhs, preconditioner);      
    constraints.distribute (m_sob_grad);
  }
  */
  
  template <int dim>
  void MySolver<dim>::make_grid_custom ()
  {
    TimerOutput::Scope timing_section(m_computing_timer, "");
    Point<dim,double> pt1(0,m_ymin); 
    Point<dim,double> pt2(m_xmax,m_ymax);
    
    CPotential<dim> Potential_fct ( m_omega );
    
    GridGenerator::hyper_rectangle(triangulation, pt2, pt1);
    triangulation.refine_global(5);
    
    double isovalues[] = {20,15,10};

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

  template <int dim>
  void MySolver<dim>::make_grid ()
  {
    TimerOutput::Scope timing_section(m_computing_timer, "");
    Point<dim,double> pt1(0,m_ymin); 
    Point<dim,double> pt2(m_xmax,m_ymax);
    
    GridGenerator::hyper_rectangle(triangulation, pt2, pt1);
    triangulation.refine_global(1);
    
    
  }
  
  template <int dim>
  void MySolver<dim>::setup_system ()
  {
    TimerOutput::Scope timing_section(m_computing_timer, "");

    dof_handler.distribute_dofs (fe);
      
    locally_owned_dofs = dof_handler.locally_owned_dofs ();
    DoFTools::extract_locally_relevant_dofs (dof_handler, locally_relevant_dofs);
    m_Psi.reinit (locally_owned_dofs, mpi_communicator);
    m_Psi_0.reinit (locally_owned_dofs, mpi_communicator);
    m_Psi_0_punkt.reinit (locally_owned_dofs, mpi_communicator);
    m_z.reinit (locally_owned_dofs, mpi_communicator);
    m_y.reinit (locally_owned_dofs, mpi_communicator);
    m_system_rhs.reinit(locally_owned_dofs, mpi_communicator);
    m_system_rhs_2.reinit(locally_owned_dofs, mpi_communicator);
    m_workspace_1.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    m_workspace_ng.reinit (locally_owned_dofs, mpi_communicator);
    m_error_per_cell.reinit(triangulation.n_active_cells());
    
    cout << "(" << m_rank << ") locally_owned_dofs = " << m_Psi.local_size()  << endl;
    
    vector<bool> mask (dof_handler.get_fe().n_components(), true );
    
    constraints.clear ();
    constraints.reinit (locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints (dof_handler, constraints);
    VectorTools::interpolate_boundary_values (dof_handler, 0, ZeroFunction<dim>(), constraints, ComponentMask(mask));
    constraints.close ();

    DynamicSparsityPattern dsp (locally_relevant_dofs);
    DoFTools::make_sparsity_pattern (dof_handler, dsp, constraints, false);
    SparsityTools::distribute_sparsity_pattern (dsp, dof_handler.n_locally_owned_dofs_per_processor(), mpi_communicator, locally_relevant_dofs);
    m_system_matrix.reinit (locally_owned_dofs, locally_owned_dofs, dsp, mpi_communicator);    

    // stuff for the second dof handler
    dof_handler_2.distribute_dofs (fe_2);

    locally_owned_dofs_2 = dof_handler_2.locally_owned_dofs ();
    DoFTools::extract_locally_relevant_dofs (dof_handler_2, locally_relevant_dofs_2);

    m_Psi_C_ghosted.reinit (locally_owned_dofs_2, locally_relevant_dofs_2, mpi_communicator);

    constraints_2.clear ();
    constraints_2.reinit (locally_relevant_dofs_2);
    DoFTools::make_hanging_node_constraints (dof_handler_2, constraints_2);
    VectorTools::interpolate_boundary_values (dof_handler_2, 0, ZeroFunction<dim>(2), constraints_2, ComponentMask(mask));
    constraints_2.close ();

    
  }

  template <int dim>
  void MySolver<dim>::output_guess ()
  {
    TimerOutput::Scope timing_section(m_computing_timer, "");
    
    constraints.distribute(m_Psi);
    m_workspace_1=m_Psi;
    
    DataOut<dim> data_out;
    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (m_workspace_1, "m_Psi");
    data_out.build_patches ();
    data_out.write_vtu_in_parallel ("guess.vtu",mpi_communicator);

    
  }

  template <int dim>
  void MySolver<dim>::output_results ( string path, string prefix )
  {
    TimerOutput::Scope timing_section(m_computing_timer, "");

    string filename;

    Vector<float> subdomain (triangulation.n_active_cells());
    for (unsigned int i=0; i<subdomain.size(); ++i)
      subdomain(i) = triangulation.locally_owned_subdomain();

    constraints.distribute(m_Psi);
    m_workspace_1=m_Psi;

    DataOut<dim> data_out;
    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (m_workspace_1, "Psi");
    data_out.add_data_vector (m_error_per_cell, "error per cell");
    data_out.add_data_vector (subdomain, "subdomain");
    data_out.build_patches ();

    filename = path + prefix + "-" + Utilities::int_to_string (m_counter,5) + ".vtu";
    data_out.write_vtu_in_parallel (filename.c_str(), mpi_communicator);

        
  }

  template <int dim>
  int MySolver<dim>::DoIter ( string path )
  {
    int retval = Status::SUCCESS;
    
    m_table.clear();

    pugi::xml_document doc;
    if (!doc.load_file("info.xml")) throw;
    
    m_mu_0 = stod(doc.child( "INFO" ).child( "MU" ).child_value());
   
    compute_tangent();
/*
    compute_newton_direction_natural();

    compute_newton_direction_arclength();

    m_mu += m_delta_mu;

    m_Psi.add( m_stepsize, m_newton_update); 
    constraints.distribute(m_Psi);
*/
   
    return retval;
  }

  template <int dim>
  void MySolver<dim>::run( const std::string& filename )
  {
    load( filename );
   
  }
    
  template<int dim>
  void MySolver<dim>::save( string filename )
  {
    m_workspace_1=m_Psi;
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

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);
  {
    BreedSolver::MySolver<2> solver("params.xml");
    solver.run ( "final.bin" );
  }
return EXIT_SUCCESS;
}