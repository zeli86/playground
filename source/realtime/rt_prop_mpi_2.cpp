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
 * Purpose: Real time propagation for the Gross-Pitaevskii equation (cartesian coordinates)
 * Method: Operator splitting - the linear step is solved with Crank-Nicolson 
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

#include <deal.II/lac/petsc_parallel_sparse_matrix.h>
#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_precondition.h>
#include <deal.II/lac/sparsity_tools.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_out.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_q.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/fe_field_function.h>

#include <deal.II/base/utilities.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/function_parser.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/timer.h>

#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/solution_transfer.h>

#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <iomanip>

#include "global.h"
#include "mpi.h"
#include "MyParameterHandler.h"
#include "my_table.h"
#include "functions.h"
#include "MyComplexTools.h"

namespace realtime_propagation
{
  enum Status { SUCCESS, FAILED };

  using namespace std;
  using namespace dealii;

  template <int dim>
  class MySolver
  {
  public:
    MySolver( const std::string& );
    ~MySolver();

    void run ();
    double Particle_Number( LA::MPI::Vector& );
    void Expectation_value_momentum( LA::MPI::Vector&, double* );
    void Expectation_value_position( LA::MPI::Vector&, double* );
    void Expectation_value_variance( LA::MPI::Vector&, double*, double* );

  protected:
    void make_grid();
    void setup_system();
    
    void Do_NL_Step();
    void Do_Lin_Step( const double );
    void output_results ( string );
    void load( string );
    void save( string );

    MyParameterHandler m_ph;
    MPI_Comm mpi_communicator;
    parallel::distributed::Triangulation<dim> triangulation;
    FESystem<dim> fe;
    DoFHandler<dim> dof_handler;
    IndexSet locally_owned_dofs;
    IndexSet locally_relevant_dofs;
    ConstraintMatrix constraints;

    LA::MPI::SparseMatrix system_matrix;
    LA::MPI::Vector system_rhs;
    LA::MPI::Vector sol;
    LA::MPI::Vector m_Psi; 

    ConditionalOStream pcout;
    ofstream m_computing_timer_log;
    TimerOutput m_computing_timer;    

    bool m_root;

    double m_gs;
    double m_t;
    double m_dt;
    double m_dth; // 0.5*m_dt
    double m_res;
    vector<double> m_omega;

    double m_xmin, m_xmax;
    double m_ymin, m_ymax;
    double m_zmin, m_zmax;

    int m_rank;
    unsigned m_NA;
    unsigned m_NK;
    
    MyTable m_table;    
  };

  template <int dim>
  MySolver<dim>::MySolver ( const std::string& xmlfilename ) 
    : 
    m_ph(xmlfilename),
    mpi_communicator (MPI_COMM_WORLD),
    triangulation (mpi_communicator, typename Triangulation<dim>::MeshSmoothing(Triangulation<dim>::smoothing_on_refinement|Triangulation<dim>::smoothing_on_coarsening)),
    fe (FE_Q<dim>(gl_degree_fe), 2),
    dof_handler (triangulation),
    pcout (cout, (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)),
    m_computing_timer_log("benchmark.txt"),
    m_computing_timer(mpi_communicator, m_computing_timer_log, TimerOutput::summary, TimerOutput:: cpu_and_wall_times )
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
      m_NK = m_ph.Get_Algorithm("NK",0);
      m_dt = m_ph.Get_Algorithm("dt",0);    
    }
    catch( const std::string info )
    {
      std::cerr << info << endl;
      MPI_Abort( mpi_communicator, 0 );
    }    
    m_dth = 0.5*m_dt;
    m_t = 0;
    m_root = (Utilities::MPI::this_mpi_process(mpi_communicator) == 0);
    MPI_Comm_rank(mpi_communicator, &m_rank);
  }

  template <int dim>
  MySolver<dim>::~MySolver ()
  {
    dof_handler.clear ();
  }

  template<int dim>
  double MySolver<dim>::Particle_Number( LA::MPI::Vector& vec )
  {
    m_computing_timer.enter_section(__func__);
    double tmp1=0.0;
    
    const QGauss<dim> quadrature_formula(fe.degree+1);
    FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_JxW_values);

    const unsigned n_q_points = quadrature_formula.size();
    vector<Vector<double>> vec_vals(n_q_points,Vector<double>(2));

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for( ; cell!=endc; cell++ )
    {
      if( cell->is_locally_owned() )
      {
        fe_values.reinit (cell);
        fe_values.get_function_values( vec, vec_vals );
        for( unsigned qp=0; qp<n_q_points; qp++ )
         tmp1 += fe_values.JxW(qp)*(vec_vals[qp][0]*vec_vals[qp][0]+vec_vals[qp][1]*vec_vals[qp][1]);
      }
    }

    double retval;
    MPI_Allreduce( &tmp1, &retval, 1, MPI_DOUBLE, MPI_SUM, mpi_communicator);
    m_computing_timer.exit_section();
  return retval;
  }  
  
  template <int dim>
  void MySolver<dim>::make_grid ()
  {
    m_computing_timer.enter_section(__func__);
#if DIMENSION==2
    Point<dim,double> pt1( m_xmin, m_ymin );
    Point<dim,double> pt2( m_xmax, m_ymax );
#endif
#if DIMENSION==3
    Point<dim,double> pt1( m_xmin, m_ymin, m_zmin );
    Point<dim,double> pt2( m_xmax, m_ymax, m_zmax );
#endif

    GridGenerator::hyper_rectangle(triangulation, pt2, pt1);
    triangulation.refine_global(1);
    m_computing_timer.exit_section();
  }
  
  template <int dim>
  void MySolver<dim>::setup_system()
  {
    m_computing_timer.enter_section(__func__);

    dof_handler.distribute_dofs (fe);
   
    locally_owned_dofs = dof_handler.locally_owned_dofs ();
    DoFTools::extract_locally_relevant_dofs (dof_handler, locally_relevant_dofs);

    m_Psi.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    sol.reinit(locally_owned_dofs, mpi_communicator); // no ghosts
    system_rhs.reinit(locally_owned_dofs, mpi_communicator); // no ghosts

    constraints.clear ();
    constraints.reinit (locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints (dof_handler, constraints);
    VectorTools::interpolate_boundary_values (dof_handler, 0, ZeroFunction<dim>(2), constraints);
    constraints.close ();

    DynamicSparsityPattern dsp (locally_relevant_dofs);
    DoFTools::make_sparsity_pattern (dof_handler, dsp, constraints, false);
    SparsityTools::distribute_sparsity_pattern (dsp, dof_handler.n_locally_owned_dofs_per_processor(), mpi_communicator, locally_relevant_dofs);
    system_matrix.reinit (locally_owned_dofs, locally_owned_dofs, dsp, mpi_communicator);

    m_computing_timer.exit_section();
  }  
  
  template<int dim>
  void MySolver<dim>::Expectation_value_position( LA::MPI::Vector& vec, double* retval )
  {
    m_computing_timer.enter_section(__func__);
    double tmp[] = {0,0,0}, JxWxn;
    
    const QGauss<dim>  quadrature_formula(fe.degree+1);
    FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_quadrature_points|update_JxW_values);

    const unsigned int n_q_points = quadrature_formula.size();
    vector<Vector<double>> vec_vals(n_q_points,Vector<double>(2));
    Point<dim> spacept;

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for( ; cell!=endc; cell++ )
    {
      if( cell->is_locally_owned() )
      {
        fe_values.reinit (cell);
        fe_values.get_function_values( vec, vec_vals );
        for( unsigned qp=0; qp<n_q_points; qp++ )
        {
          JxWxn = fe_values.JxW(qp)*(vec_vals[qp][0]*vec_vals[qp][0]+vec_vals[qp][1]*vec_vals[qp][1]);
          spacept = fe_values.quadrature_point(qp);
          tmp[0] += spacept[0]*JxWxn;
          tmp[1] += spacept[1]*JxWxn;
#if dim == 3
          tmp[2] += spacept[2]*JxWxn;
#endif
        }
      }
    }
    MPI_Allreduce( tmp, retval, 3, MPI_DOUBLE, MPI_SUM, mpi_communicator);
    m_computing_timer.exit_section();
  }  

  template<int dim>
  void MySolver<dim>::Expectation_value_variance( LA::MPI::Vector& vec, double* retval, double* x1 )
  {
    m_computing_timer.enter_section(__func__);
    double tmp[] = {0,0,0}, JxWxn;
    
    const QGauss<dim>  quadrature_formula(fe.degree+1);
    FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_quadrature_points|update_JxW_values);

    const unsigned n_q_points = quadrature_formula.size();
    vector<Vector<double>> vec_vals(n_q_points,Vector<double>(2));
    Point<dim> spacept;

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for( ; cell!=endc; cell++ )
    {
      if( cell->is_locally_owned() )
      {
        fe_values.reinit (cell);
        fe_values.get_function_values( vec, vec_vals );
        for( unsigned qp=0; qp<n_q_points; qp++ )
        {
          JxWxn = fe_values.JxW(qp)*(vec_vals[qp][0]*vec_vals[qp][0]+vec_vals[qp][1]*vec_vals[qp][1]);
          spacept = fe_values.quadrature_point(qp);
          tmp[0] += (spacept[0]-x1[0])*(spacept[0]-x1[0])*JxWxn;
          tmp[1] += (spacept[1]-x1[1])*(spacept[1]-x1[1])*JxWxn;
#if dim == 3
          tmp[2] += (spacept[2]-x1[2])*(spacept[2]-x1[2])*JxWxn;
#endif
        }
      }
    }
    MPI_Allreduce( tmp, retval, 3, MPI_DOUBLE, MPI_SUM, mpi_communicator);
    m_computing_timer.exit_section();
  }  
  
  template<int dim>
  void MySolver<dim>::Expectation_value_momentum( LA::MPI::Vector& vec, double* retval )
  {
    m_computing_timer.enter_section(__func__);
    double tmp[] = {0,0,0}, JxW;
    
    const QGauss<dim>  quadrature_formula(fe.degree+1);
    FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_gradients|update_JxW_values);

    const unsigned int n_q_points = quadrature_formula.size();
    vector<Vector<double>> vec_vals(n_q_points,Vector<double>(2));
    vector<vector<Tensor<1,dim>>> vec_grads(n_q_points, vector<Tensor<1,dim>>(2));

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for( ; cell!=endc; cell++ )
    {
      if( cell->is_locally_owned() )
      {
        fe_values.reinit (cell);
        fe_values.get_function_values( vec, vec_vals );
        fe_values.get_function_gradients( vec, vec_grads );

        for( unsigned qp=0; qp<n_q_points; qp++ )
        {
          JxW = fe_values.JxW(qp);
          tmp[0] += JxW*(vec_vals[qp][0]*vec_grads[qp][1][0] - vec_vals[qp][1]*vec_grads[qp][0][0]);
          tmp[1] += JxW*(vec_vals[qp][0]*vec_grads[qp][1][1] - vec_vals[qp][1]*vec_grads[qp][0][1]);
#if dim == 3
          tmp[2] += JxW*(vec_vals[qp][0]*vec_grads[qp][1][2] - vec_vals[qp][1]*vec_grads[qp][0][2]); 
#endif
        }
      }
    }
    MPI_Allreduce( tmp, retval, 3, MPI_DOUBLE, MPI_SUM, mpi_communicator);
    m_computing_timer.exit_section();
  }  

  template <int dim>
  void MySolver<dim>::output_results ( string path ) 
  {
    m_computing_timer.enter_section(__func__);
    string filename;
    
    vector<std::string> solution_names;
    solution_names.push_back ("Re_Psi");
    solution_names.push_back ("Im_Psi");    

    DataOut<dim> data_out;
    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector ( m_Psi, solution_names );
    data_out.build_patches ();
    
    filename = path + "solution-" + to_string(m_t) + ".vtu";
    data_out.write_vtu_in_parallel ( filename.c_str(), mpi_communicator );

    m_computing_timer.exit_section();
  }   

  template<int dim>
  void MySolver<dim>::save( string filename )
  {
    parallel::distributed::SolutionTransfer<dim,LA::MPI::Vector> solution_transfer(dof_handler);
    solution_transfer.prepare_serialization(m_Psi);

    triangulation.save( filename.c_str() );
  }

  template<int dim>
  void MySolver<dim>::load( string filename )
  {
    make_grid();
    triangulation.load( filename.c_str() );

    setup_system();

    parallel::distributed::SolutionTransfer<dim,LA::MPI::Vector> solution_transfer(dof_handler);
    solution_transfer.deserialize(system_rhs);

    m_Psi = system_rhs;
  }    

  template<int dim> 
  void MySolver<dim>::Do_NL_Step()
  {
    m_computing_timer.enter_section(__func__);

    pcout << "Computing nonlinear step, t = " << m_dt << endl;

    MPI::MyComplexTools::AssembleSystem_NL_Step<dim>( dof_handler, fe, constraints, m_Psi, m_gs*m_dt, system_matrix, system_rhs );

    SolverControl solver_control (sol.size(), 1e-15);
    PETScWrappers::SolverCG solver (solver_control, mpi_communicator);
    PETScWrappers::PreconditionSSOR preconditioner(system_matrix);
    solver.solve (system_matrix, sol, system_rhs, preconditioner);
    
    constraints.distribute (sol);
    m_Psi = sol;

    m_computing_timer.exit_section();   
  }
  
  template<int dim> 
  void MySolver<dim>::Do_Lin_Step( const double dt )
  {
    m_computing_timer.enter_section(__func__);
    pcout << "Computing linear step (" << dt << "), t = " << m_t << endl;

    CPotential<dim> Potential ( m_omega );
    MPI::MyComplexTools::AssembleSystem_LIN_Step<dim>( dof_handler, fe, constraints, m_Psi, Potential, dt, system_matrix, system_rhs );

    SolverControl solver_control (sol.size(), 1e-15);
    PETScWrappers::SolverCG solver (solver_control, mpi_communicator);
    PETScWrappers::PreconditionSSOR preconditioner(system_matrix);
    solver.solve (system_matrix, sol, system_rhs, preconditioner);
    
    constraints.distribute (sol);
    m_Psi = sol;    

    m_t += dt;
    m_computing_timer.exit_section();   
  }
  
/* 
  template <int dim>
  void MySolver<dim>::solve ()
  {
    m_computing_timer.enter_section(__func__);
    SolverControl solver_control;
    
    sol=0;
    PETScWrappers::SparseDirectMUMPS solver(solver_control, mpi_communicator);
    solver.set_symmetric_mode(false);
    solver.solve(system_matrix, sol, system_rhs);

    constraints.distribute (sol);
    m_Psi = sol;
    m_computing_timer.exit_section();
  }
*/
 
  template <int dim>
  void MySolver<dim>::run()
  {
    double T, N, W;

    load( "final.bin" );

    double p[] = {0,0,0};
    double pos[] = {0,0,0};
    double var[] = {0,0,0};

    output_results("");

    N = Particle_Number(m_Psi);
    pcout << "N == " << N << endl;
    
    Expectation_value_position( m_Psi, pos );
    Expectation_value_momentum( m_Psi, p );
    Expectation_value_variance( m_Psi, var, pos );

    pcout << "N == " << Particle_Number(m_Psi) << endl;
    pcout << "p == " << p[0]/N << ", " << p[1]/N << ", " << p[2]/N << endl;
    pcout << "pos == " << pos[0]/N << ", " << pos[1]/N << ", " << pos[2]/N << endl;
    pcout << "var == " << var[0]/N << ", " << var[1]/N << ", " << var[2]/N << endl;

    for( unsigned i=1; i<=m_NA; i++ )
    {
      Do_Lin_Step( m_dth );
      for( unsigned j=2; j<=m_NK; j++ )
      {
	      Do_NL_Step();
        Do_Lin_Step( m_dt );
      }
      Do_NL_Step();
      Do_Lin_Step( m_dth );

      Expectation_value_position( m_Psi, pos );
      Expectation_value_momentum( m_Psi, p );
      Expectation_value_variance( m_Psi, var, pos );

      pcout << "N == " << Particle_Number(m_Psi) << endl;
      pcout << "p == " << p[0]/N << ", " << p[1]/N << ", " << p[2]/N << endl;
      pcout << "pos == " << pos[0]/N << ", " << pos[1]/N << ", " << pos[2]/N << endl;
      pcout << "var == " << var[0]/N << ", " << var[1]/N << ", " << var[2]/N << endl;

//       save( "last.bin" );
      output_results("");
    }
  }
} // end of namespace 

int main ( int argc, char *argv[] )
{
  using namespace dealii;
  deallog.depth_console (0);

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv );
  {
    realtime_propagation::MySolver<DIMENSION> solver("params.xml");
    solver.run();
  }  
return EXIT_SUCCESS;
}