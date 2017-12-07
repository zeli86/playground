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
#include <gsl/gsl_sf.h>

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
#include "my_table.h"
#include "MyComplexTools.h"

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

    void run ();

    MPI_Comm mpi_communicator;
  protected:
 
    double m_xmin, m_xmax;
    double m_ymin, m_ymax;
    double m_zmin, m_zmax;
    double m_res, m_res_old, m_resp;
    double m_final_error;
    double m_N;
    
    double m_mu;
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
    void setup_system();
    void assemble_rhs();
    void assemble_system();
    void compute_Psi_sob();
    void compute_mu();
    void save( string );
    void save_C( string );
    void load( string );
    void dump_info_xml( const string );
    void Project_gradient();
    
    double Particle_Number( LA::MPI::Vector& );
    
    void solve();
    void compute_E_lin( LA::MPI::Vector&, double&, double&, double& );
    void estimate_error( double& );

    void output_results ( string, string = "step" );
    void output_guess ();

    ofstream m_computing_timer_log;
    TimerOutput m_computing_timer;    
    MyParameterHandler m_ph;
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

    LA::MPI::SparseMatrix m_system_matrix;
    LA::MPI::Vector m_system_rhs;
    LA::MPI::Vector m_Psi;
    LA::MPI::Vector m_Psi_final;
    LA::MPI::Vector m_sob_grad;
    LA::MPI::Vector m_Psi_sob;
    LA::MPI::Vector m_workspace_1;
    LA::MPI::Vector m_workspace_2;
    LA::MPI::Vector m_workspace_3;

    LA::MPI::Vector m_Psi_C_initial_ghosted;
    LA::MPI::Vector m_Psi_C_final_ghosted;

    Vector<double> m_error_per_cell;

    //string m_guess_str;
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

      m_NA = m_ph.Get_Algorithm("NA",0);
      m_epsilon = m_ph.Get_Algorithm("epsilon");
    }
    catch( const std::string info )
    {
      std::cerr << info << endl;
      MPI_Abort( mpi_communicator, 0 );
    }    

    m_root = (Utilities::MPI::this_mpi_process(mpi_communicator) == 0);
    MPI_Comm_rank(mpi_communicator, &m_rank);

    m_counter=0;
    m_final_error=0;    
  }

  template <int dim>
  MySolver<dim>::~MySolver ()
  {
    dof_handler.clear();
    dof_handler_2.clear();
  }
 
  template <int dim>
  void MySolver<dim>::compute_E_lin( LA::MPI::Vector& vec, double& T, double& N, double& W )
  {
    m_computing_timer.enter_section(__func__);
    
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
        for ( unsigned qp=0; qp<n_q_points; qp++ )
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
    
    m_computing_timer.exit_section();
  }

  template<int dim>
  double MySolver<dim>::Particle_Number( LA::MPI::Vector& vec )
  {
    m_computing_timer.enter_section(__func__);
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

        for ( unsigned qp=0; qp<n_q_points; qp++ )
          tmp1 += fe_values.JxW(qp)*(vals[qp]*vals[qp]);
      }
    }
    MPI_Allreduce( &tmp1, &retval, 1, MPI_DOUBLE, MPI_SUM, mpi_communicator);
    m_computing_timer.exit_section();
  return retval;
  }
  
  template<int dim>
  void MySolver<dim>::Project_gradient()
  {
    m_computing_timer.enter_section(__func__);
    double tmp1[]={0,0}, sum[]={0,0};
   
    constraints.distribute(m_Psi_sob);
    m_workspace_1 = m_Psi_sob;
    
    constraints.distribute(m_sob_grad);
    m_workspace_2 = m_sob_grad;
    
    constraints.distribute(m_Psi);
    m_workspace_3 = m_Psi;

    const QGauss<dim>  quadrature_formula(fe.degree+1);
    FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_JxW_values|update_quadrature_points);

    const unsigned n_q_points = quadrature_formula.size();
    vector<double> vals_1(n_q_points);
    vector<double> vals_2(n_q_points);
    vector<double> vals_3(n_q_points);
    
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
      if( cell->is_locally_owned() )
      {
        fe_values.reinit (cell);
        fe_values.get_function_values( m_workspace_1, vals_1 );      
        fe_values.get_function_values( m_workspace_2, vals_2 );      
        fe_values.get_function_values( m_workspace_3, vals_3 );      

        for ( unsigned qp=0; qp<n_q_points; qp++ )
        {
          double JxW = fe_values.JxW(qp);
          tmp1[0] += JxW*(vals_3[qp]*vals_2[qp]);
          tmp1[1] += JxW*(vals_3[qp]*vals_1[qp]);
        }
      }
    }
    MPI_Allreduce( tmp1, sum, 2, MPI_DOUBLE, MPI_SUM, mpi_communicator);
    
    m_sob_grad.add( -sum[0]/sum[1], m_Psi_sob );
    
    m_computing_timer.exit_section();
  }
  
  template <int dim>
  void MySolver<dim>::estimate_error ( double& err )
  {
    m_computing_timer.enter_section(__func__);
    
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

        for ( unsigned qp=0; qp<n_q_points; qp++ )
        {
          JxW = fe_values.JxW(qp);
          Q1 = Potential.value(fe_values.quadrature_point(qp)) - m_mu + m_gs*(vals[qp]*vals[qp]);

          for ( unsigned i=0; i<dofs_per_cell; i++ )
          {
            cell_rhs(i) += JxW*(grads[qp]*fe_values.shape_grad(i,qp) + Q1*vals[qp]*fe_values.shape_value(i,qp));
            for ( unsigned j=0; j<dofs_per_cell; j++ )
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
    
    m_workspace_1=m_sob_grad;
    VectorTools::integrate_difference ( dof_handler, m_workspace_1, ZeroFunction<dim>(2), m_error_per_cell, QGauss<dim>(fe.degree+2), VectorTools::L2_norm);    
    const double total_local_error = m_error_per_cell.l2_norm();
    err = std::sqrt (Utilities::MPI::sum (total_local_error * total_local_error, MPI_COMM_WORLD)); 
    m_computing_timer.exit_section();
  }
  
  template <int dim>
  void MySolver<dim>::assemble_rhs ()
  {
    m_computing_timer.enter_section(__func__);

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
    
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for ( ; cell!=endc; ++cell )
    {
      if( cell->is_locally_owned() )
      {
        cell_rhs=0;

        fe_values.reinit (cell);
        fe_values.get_function_values(m_workspace_1, vals);
        fe_values.get_function_gradients(m_workspace_1, grads);

        for ( unsigned qp=0; qp<n_q_points; qp++ )
        {
          double JxW = fe_values.JxW(qp);
          double Q1 = Potential.value(fe_values.quadrature_point(qp)) + m_gs*(vals[qp]*vals[qp]); 

          for (unsigned i=0; i<dofs_per_cell; ++i)
            cell_rhs(i) += JxW*(grads[qp]*fe_values.shape_grad(i,qp) + Q1*vals[qp]*fe_values.shape_value(i,qp));
        }
        cell->get_dof_indices (local_dof_indices);
        constraints.distribute_local_to_global(cell_rhs, local_dof_indices, m_system_rhs);
      }
    }
    m_system_rhs.compress(VectorOperation::add);   
    m_computing_timer.exit_section();
  }

  template<int dim>
  void MySolver<dim>::load( string filename )
  {
    make_grid();
    triangulation.load( filename.c_str() );

    setup_system();

    parallel::distributed::SolutionTransfer<dim,LA::MPI::Vector> solution_transfer(dof_handler);
    solution_transfer.deserialize(m_Psi_final);
    constraints.distribute(m_Psi_final);
  }    

  template <int dim>
  void MySolver<dim>::assemble_system ()
  {
    m_computing_timer.enter_section(__func__);

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

        for ( unsigned qp=0; qp<n_q_points; qp++ )
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
    m_computing_timer.exit_section();
  }

  template <int dim>
  void MySolver<dim>::compute_Psi_sob ()
  {
    m_computing_timer.enter_section(__func__);

    const QGauss<dim> quadrature_formula(fe.degree+1);
    
    m_system_matrix=0;
    m_system_rhs=0;
    
    constraints.distribute(m_Psi);
    m_workspace_1=m_Psi;    

    FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_gradients|update_JxW_values|update_quadrature_points);

    const unsigned dofs_per_cell = fe.dofs_per_cell;
    const unsigned n_q_points = quadrature_formula.size();

    Vector<double> cell_rhs (dofs_per_cell);
    vector<double> vals(n_q_points);
    FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);
    vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    double JxW, val;
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for ( ; cell!=endc; ++cell )
    {
      if( cell->is_locally_owned() )
      {
        cell_matrix=0;
        cell_rhs=0;
        
        fe_values.reinit (cell);
        fe_values.get_function_values(m_workspace_1, vals);

        for ( unsigned qp=0; qp<n_q_points; qp++ )
        {
          JxW = fe_values.JxW(qp);
          val = vals[qp];

          for ( unsigned i=0; i<dofs_per_cell; i++ )
          {
            cell_rhs(i)=JxW*val*fe_values.shape_value(i,qp);
            for ( unsigned j=0; j<dofs_per_cell; j++ )
              cell_matrix(i,j) += JxW*(fe_values.shape_grad(i,qp)*fe_values.shape_grad(j,qp) + fe_values.shape_value(i,qp)*fe_values.shape_value(j,qp));
          
          }
        }
        cell->get_dof_indices (local_dof_indices);
        constraints.distribute_local_to_global(cell_matrix, cell_rhs, local_dof_indices, m_system_matrix, m_system_rhs);
      }
    }
    m_system_matrix.compress(VectorOperation::add);
    m_system_rhs.compress(VectorOperation::add);
    
    pcout << "Solving..." << endl;
    SolverControl solver_control;
    PETScWrappers::SparseDirectMUMPS solver(solver_control, mpi_communicator);
    solver.set_symmetric_mode(false);
    solver.solve(m_system_matrix, m_Psi_sob, m_system_rhs);
    constraints.distribute (m_Psi_sob);    
    
/*
    m_Psi_sob=0;    
    SolverControl solver_control ( m_Psi_sob.size(), 1e-15 );
    PETScWrappers::SolverCG cg (solver_control, mpi_communicator);

    PETScWrappers::PreconditionSOR preconditioner;
    PETScWrappers::PreconditionSOR::AdditionalData data;
    preconditioner.initialize(m_system_matrix, data);    

    cg.solve (m_system_matrix, m_Psi_sob, m_system_rhs, preconditioner);    
*/    
    m_computing_timer.exit_section();
  }
  
  template <int dim>
  void MySolver<dim>::compute_mu()
  {
    m_computing_timer.enter_section(__func__);

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

        for ( unsigned qp=0; qp<n_q_points; qp++ )
        {
          double uq = vals[qp]*vals[qp];
          psum += fe_values.JxW(qp)*(grads[qp]*grads[qp] + (Potential.value(fe_values.quadrature_point(qp)) + m_gs*uq)*uq); 
        }
      }
    }    
    MPI_Allreduce( &psum, &m_mu, 1, MPI_DOUBLE, MPI_SUM, mpi_communicator);
    m_computing_timer.exit_section();
  }
  
  
  template <int dim>
  void MySolver<dim>::solve ()
  {
    m_computing_timer.enter_section(__func__);
    pcout << "Solving..." << endl;
    
    SolverControl solver_control;
    
    PETScWrappers::SparseDirectMUMPS solver(solver_control, mpi_communicator);
    solver.set_symmetric_mode(false);
    solver.solve(m_system_matrix, m_sob_grad, m_system_rhs);
    constraints.distribute (m_sob_grad);
    m_computing_timer.exit_section();
    /*
    m_Psi_sob=0;    
    SolverControl solver_control ( m_Psi_sob.size(), 1e-15 );
    PETScWrappers::SolverCG cg (solver_control, mpi_communicator);

    PETScWrappers::PreconditionSOR preconditioner;
    PETScWrappers::PreconditionSOR::AdditionalData data;
    preconditioner.initialize(m_system_matrix, data);    

    cg.solve (m_system_matrix, m_Psi_sob, m_system_rhs, preconditioner);
*/   
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
  void MySolver<dim>::make_grid ()
  {
    m_computing_timer.enter_section(__func__);
    Point<dim,double> pt1;
    Point<dim,double> pt2;

    double min[] = {m_xmin, m_ymin, m_zmin};
    double max[] = {m_xmax, m_ymax, m_zmax};

    for( int i=0; i<dim; i++ )
    {
      pt1(i) = min[i];
      pt2(i) = max[i];
    }
    
    GridGenerator::hyper_rectangle(triangulation, pt2, pt1);
    triangulation.refine_global(1);
    
    m_computing_timer.exit_section();
  }
  
  template <int dim>
  void MySolver<dim>::setup_system ()
  {
    m_computing_timer.enter_section(__func__);

    dof_handler.distribute_dofs (fe);
      
    locally_owned_dofs = dof_handler.locally_owned_dofs ();
    DoFTools::extract_locally_relevant_dofs (dof_handler, locally_relevant_dofs);
    m_Psi.reinit (locally_owned_dofs, mpi_communicator);
    m_Psi_final.reinit (locally_owned_dofs, mpi_communicator);
    m_Psi_sob.reinit (locally_owned_dofs, mpi_communicator);
    m_system_rhs.reinit(locally_owned_dofs, mpi_communicator);
    m_sob_grad.reinit (locally_owned_dofs, mpi_communicator);
    m_workspace_1.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    m_workspace_2.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    m_workspace_3.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    m_error_per_cell.reinit(triangulation.n_active_cells());
    
    cout << "(" << m_rank << ") locally_owned_dofs = " << m_Psi.local_size()  << endl;
    
    vector<bool> mask (dof_handler.get_fe().n_components(), true );
    
    constraints.clear ();
    constraints.reinit (locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints (dof_handler, constraints);
    VectorTools::interpolate_boundary_values (dof_handler, 0, ZeroFunction<dim>(), constraints, ComponentMask(mask));
    constraints.close ();

    DynamicSparsityPattern csp (locally_relevant_dofs);
    DoFTools::make_sparsity_pattern (dof_handler, csp, constraints, false);
    SparsityTools::distribute_sparsity_pattern (csp, dof_handler.n_locally_owned_dofs_per_processor(), mpi_communicator, locally_relevant_dofs);
    m_system_matrix.reinit (locally_owned_dofs, locally_owned_dofs, csp, mpi_communicator);    

    // stuff for the second dof handler
    dof_handler_2.distribute_dofs (fe_2);

    locally_owned_dofs_2 = dof_handler_2.locally_owned_dofs ();
    DoFTools::extract_locally_relevant_dofs (dof_handler_2, locally_relevant_dofs_2);

    m_Psi_C_final_ghosted.reinit (locally_owned_dofs_2, locally_relevant_dofs_2, mpi_communicator);
    m_Psi_C_initial_ghosted.reinit (locally_owned_dofs_2, locally_relevant_dofs_2, mpi_communicator);

    constraints_2.clear ();
    constraints_2.reinit (locally_relevant_dofs_2);
    DoFTools::make_hanging_node_constraints (dof_handler_2, constraints_2);
    VectorTools::interpolate_boundary_values (dof_handler_2, 0, ZeroFunction<dim>(2), constraints_2, ComponentMask(mask));
    constraints_2.close ();

    m_computing_timer.exit_section();
  }

  template <int dim>
  void MySolver<dim>::output_guess ()
  {
    m_computing_timer.enter_section(__func__);
    
    constraints.distribute(m_Psi);
    m_workspace_1=m_Psi;
    
    DataOut<dim> data_out;
    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (m_workspace_1, "m_Psi");
    data_out.build_patches ();
    data_out.write_vtu_in_parallel ("guess.vtu",mpi_communicator);

    m_computing_timer.exit_section();
  }

  template <int dim>
  void MySolver<dim>::output_results ( string path, string prefix )
  {
    m_computing_timer.enter_section(__func__);

    string filename;

    constraints.distribute(m_Psi);
    m_workspace_1=m_Psi;

    DataOut<dim> data_out;
    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (m_workspace_1, "Psi");
    //data_out.add_data_vector (m_error_per_cell, "error per cell");
    data_out.build_patches ();

    filename = path + prefix + "-" + Utilities::int_to_string (m_counter,5) + ".vtu";
    data_out.write_vtu_in_parallel (filename.c_str(), mpi_communicator);

    m_computing_timer.exit_section();    
  }

  template <int dim>
  int MySolver<dim>::DoIter ( string path )
  {
    int retval = Status::SUCCESS;
    
    double err;
    
    m_table.clear();
    
    assemble_rhs();
    
    m_res=0;
    m_counter=0;
    do
    {
      pcout << "--------------------------------------------------------------------------------" << endl;
      pcout << "- " << m_counter << endl;
     
      assemble_system();
      solve();
      
      compute_Psi_sob();
      Project_gradient();
      m_res = m_sob_grad.l2_norm();
      
      m_Psi.add( -1e-2, m_sob_grad);

      //compute_mu(m_mu);
      m_N=Particle_Number(m_Psi);
      
      if( fabs(m_N-1) > 1e-5 ) m_Psi *= 1/sqrt(m_N);
      
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
      m_table.insert( cols, MyTable::PARTICLE_NUMBER, m_N );
      //m_table.insert( cols, MyTable::total_no_cells, double(m_total_no_cells) );
      //m_table.insert( cols, MyTable::total_no_active_cells, double(m_total_no_active_cells) );

      if( m_root ) cout << m_table;
      if( m_res < m_epsilon[1] ) { retval=Status::SUCCESS; break; }

      m_counter++;
    }while( true );

    m_N = Particle_Number(m_Psi);
    
    if( m_N < 1e-5 ) retval = Status::ZERO_SOL;
    
    string filename = path + "log.csv";
    if( m_root ) m_table.dump_2_file(filename);
    
    return retval;
  }

  template <int dim>
  void MySolver<dim>::run()
  {
    int status;
    double T, N, W;

    load("final_one.bin");
    
    unsigned QN[]={0,0,0};
    CEigenfunctions<dim> Ef1( QN, m_omega );
    VectorTools::interpolate (dof_handler, Ef1, m_Psi );
   
    compute_E_lin( m_Psi, T, N, W );
    m_Psi *= sqrt(1/N);

    output_guess();
    
    pcout << setprecision(9);
    pcout << "T = " << T << endl;
    pcout << "N = " << N << endl;
    pcout << "W = " << W << endl;
    //pcout << "m_mu = " << m_mu << endl;

    status = DoIter("");

    constraints.distribute(m_Psi);
    m_workspace_1 = m_Psi;
    MyComplexTools::MPI::Interpolate_R_to_C( mpi_communicator, dof_handler, fe, m_workspace_1, dof_handler_2, fe_2, constraints_2, m_Psi_C_initial_ghosted );

    constraints.distribute(m_Psi_final);
    m_workspace_1 = m_Psi_final;
    MyComplexTools::MPI::Interpolate_R_to_C( mpi_communicator, dof_handler, fe, m_workspace_1, dof_handler_2, fe_2, constraints_2, m_Psi_C_final_ghosted );

    save_C( "oct_0.bin" );

    if( status == Status::SUCCESS )
    {
      //output_results("","final");
      //output_results("","final");
      //dump_info_xml("");
    }

    if( m_root )
    {
      ofstream ofs("log.txt");
      ofs << m_table;
    }
  }
    
  template<int dim>
  void MySolver<dim>::save( string filename )
  {
    m_workspace_1=m_Psi;
    parallel::distributed::SolutionTransfer<dim,LA::MPI::Vector> solution_transfer(dof_handler);
    solution_transfer.prepare_serialization(m_workspace_1);
    triangulation.save( filename.c_str() );
  }

  template<int dim>
  void MySolver<dim>::save_C( string filename )
  {
    std::vector<const LA::MPI::Vector*> x_system (2);
    x_system[0] = &m_Psi_C_initial_ghosted;
    x_system[1] = &m_Psi_C_final_ghosted;
    
    parallel::distributed::SolutionTransfer<dim,LA::MPI::Vector> solution_transfer(dof_handler_2);
    solution_transfer.prepare_serialization(x_system);
    triangulation.save( filename.c_str() );

/*    
    DataOut<dim> data_out;
    data_out.attach_dof_handler (dof_handler_2);
    data_out.add_data_vector (m_Psi_C_ghosted, "Psi");
    data_out.build_patches ();
    string filename = "C_final.vtu";
    data_out.write_vtu_in_parallel (filename.c_str(), mpi_communicator);    
*/   

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
    BreedSolver::MySolver<2> solver("params_one.xml");
    solver.run ();
  }
return EXIT_SUCCESS;
}