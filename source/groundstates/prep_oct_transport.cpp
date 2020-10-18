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

#include <deal.II/lac/petsc_parallel_sparse_matrix.h>
#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_precondition.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/derivative_approximation.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/function_parser.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/solution_transfer.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>

#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <locale>
#include <limits>
#include <cmath>
#include <array>

#include "global.h"
#include "mpi.h"
#include "functions.h"
#include "my_table.h"
#include "ref_pt_list.h"
#include "MyParameterHandler.h"
#include "MyRealTools.h"
#include "MyComplexTools.h"

namespace BreedSolver
{
  template <int dim>
  class MySolver;

  using namespace std;
  using namespace dealii;

  #include "my_solver_mpi_ortho_funcs.h"

#ifdef __variant_1__
  using namespace variant_1;
#endif
#ifdef __variant_2__
  using namespace variant_2;
#endif

  #include "CBase.h"

  template <int dim>
  class MySolver : public CBase<2>
  {
  public:
    MySolver( const std::string );
    virtual ~MySolver();

    void run ();

  protected:
    int DoIter();

    void save( string );
    void save_one( string );
    void make_grid();
    void make_grid_custom();
    void setup_system();
    void do_superposition();
    void estimate_error( double& );
    void Interpolate_R_to_C( LA::MPI::Vector& );

    void solve();
    void compute_contributions();

    void output_results ( string, string = "step" );
    void output_vector ( LA::MPI::Vector&, string );
    void output_guess ();

    parallel::distributed::Triangulation<dim> triangulation;
    FE_Q<dim> fe;
    DoFHandler<dim> dof_handler;
    FESystem<dim> fe_2;
    DoFHandler<dim> dof_handler_2;
    IndexSet locally_owned_dofs, locally_owned_dofs_2;
    IndexSet locally_relevant_dofs, locally_relevant_dofs_2;
    ConstraintMatrix constraints, constraints_2;

    LA::MPI::SparseMatrix m_system_matrix;
    LA::MPI::Vector m_system_rhs;
    LA::MPI::Vector m_newton_update;
    LA::MPI::Vector m_Psi_C_ghosted_i;
    LA::MPI::Vector m_Psi_C_ghosted_f;
    LA::MPI::Vector m_Psi;
    LA::MPI::Vector m_workspace_1;
    LA::MPI::Vector m_workspace_ng;
    Vector<double> m_error_per_cell;

    MyTable m_table;

    FunctionParser<dim> m_Potential;
  };

/*************************************************************************************************/
/**
 * Constructor
 */

  template <int dim>
  MySolver<dim>::MySolver ( const std::string xmlfilename ) 
    : 
    CBase<2>(xmlfilename),
    triangulation (mpi_communicator, typename Triangulation<dim>::MeshSmoothing(Triangulation<dim>::limit_level_difference_at_vertices|Triangulation<dim>::eliminate_refined_inner_islands|Triangulation<dim>::smoothing_on_refinement|Triangulation<dim>::smoothing_on_coarsening)),
    fe (gl_degree_fe),
    dof_handler (triangulation),
    fe_2 (FE_Q<dim>(gl_degree_fe),2),
    dof_handler_2 (triangulation)
  {
  }

  template <int dim>
  MySolver<dim>::~MySolver ()
  {
    dof_handler.clear ();
    dof_handler_2.clear ();
  }

  template <int dim>
  void MySolver<dim>::compute_contributions()
  {
  }

  template <int dim>
  void MySolver<dim>::make_grid_custom ()
  {
    TimerOutput::Scope timing_section(m_computing_timer, "");

    Point<dim,double> pt1;
    Point<dim,double> pt2;

    double min[] = {m_xmin, m_ymin, m_zmin};
    double max[] = {m_xmax, m_ymax, m_zmax};

    for( int i=0; i<dim; ++i )
    {
      pt1(i) = min[i];
      pt2(i) = max[i];
    }
    
    GridGenerator::hyper_rectangle(triangulation, pt2, pt1);
    
    double vals[] = {10,12,14};
    
    for( unsigned step=0; step<sizeof(vals)/sizeof(double); step++ )
    {
      typename parallel::distributed::Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active(), endc = triangulation.end();
      for( ; cell!=endc; ++cell )
        for (unsigned v=0; v < GeometryInfo<dim>::vertices_per_cell; ++v )
        {
          Point<dim> p = cell->vertex(v);
          if( fabs(p(0)) < vals[step] && fabs(p(1)) < vals[step] )
          {
            cell->set_refine_flag ();
            break;
          }
        }
      triangulation.execute_coarsening_and_refinement ();
    }
    
  }

  template <int dim>
  void MySolver<dim>::setup_system()
  {
    TimerOutput::Scope timing_section(m_computing_timer, "");

    dof_handler.distribute_dofs (fe);

    locally_owned_dofs = dof_handler.locally_owned_dofs ();
    DoFTools::extract_locally_relevant_dofs (dof_handler, locally_relevant_dofs);
    
    m_Psi.reinit (locally_owned_dofs, mpi_communicator);
    m_newton_update.reinit (locally_owned_dofs, mpi_communicator);
    m_system_rhs.reinit(locally_owned_dofs, mpi_communicator);
    m_workspace_ng.reinit (locally_owned_dofs, mpi_communicator);
    m_workspace_1.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    m_error_per_cell.reinit(triangulation.n_active_cells());
    
    cout << "(" << m_rank << ") locally_owned_dofs = " << m_Psi.local_size()  << endl;
     
    m_system_rhs=0;

    constraints.clear ();
    constraints.reinit (locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints (dof_handler, constraints);
    VectorTools::interpolate_boundary_values (dof_handler, 0, ZeroFunction<dim>(), constraints);
    constraints.close ();

    DynamicSparsityPattern csp (locally_relevant_dofs);
    DoFTools::make_sparsity_pattern (dof_handler, csp, constraints, false);
    SparsityTools::distribute_sparsity_pattern (csp, dof_handler.n_locally_owned_dofs_per_processor(), mpi_communicator, locally_relevant_dofs);
    m_system_matrix.reinit (locally_owned_dofs, locally_owned_dofs, csp, mpi_communicator);
    
    // stuff for the second dof handler
    dof_handler_2.distribute_dofs (fe_2);

    locally_owned_dofs_2 = dof_handler_2.locally_owned_dofs ();
    DoFTools::extract_locally_relevant_dofs (dof_handler_2, locally_relevant_dofs_2);

    m_Psi_C_ghosted_i.reinit (locally_owned_dofs_2, locally_relevant_dofs_2, mpi_communicator);
    m_Psi_C_ghosted_f.reinit (locally_owned_dofs_2, locally_relevant_dofs_2, mpi_communicator);

    constraints_2.clear ();
    constraints_2.reinit (locally_relevant_dofs_2);
    DoFTools::make_hanging_node_constraints (dof_handler_2, constraints_2);
    VectorTools::interpolate_boundary_values (dof_handler_2, 0, ZeroFunction<dim>(2), constraints_2);
    constraints_2.close ();

    
  }

  template <int dim>
  void MySolver<dim>::output_results ( string path, string prefix )
  {
    TimerOutput::Scope timing_section(m_computing_timer, "");

    constraints.distribute(m_Psi);
    m_workspace_1=m_Psi;
    
    string filename;

    Vector<float> subdomain (triangulation.n_active_cells());
    for (unsigned i=0; i<subdomain.size(); ++i)
      subdomain(i) = triangulation.locally_owned_subdomain();

    DataOut<dim> data_out;
    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (m_workspace_1, "Psi_sol");
    data_out.add_data_vector (m_error_per_cell, "error per cell");
    data_out.add_data_vector (subdomain, "subdomain");
    data_out.build_patches ();

    filename = path + prefix + "-" + Utilities::int_to_string (m_counter,5) + ".vtu";
    data_out.write_vtu_in_parallel (filename.c_str(), mpi_communicator);

        
  }

  template <int dim>
  void MySolver<dim>::output_vector ( LA::MPI::Vector& vec, string filename )
  {
    TimerOutput::Scope timing_section(m_computing_timer, "");

    constraints.distribute(vec);
    m_workspace_1=vec;

    DataOut<dim> data_out;
    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (m_workspace_1, "vec");
    data_out.build_patches ();
    data_out.write_vtu_in_parallel (filename.c_str(), mpi_communicator);

        
  }  

  template <int dim>
  int MySolver<dim>::DoIter ()
  {
    int retval = Status::SUCCESS;
    
    m_table.clear();
    
    // Standard Newton 
    do
    {
      pcout << "--------------------------------------------------------------------------------" << endl;
      pcout << "--- " << m_counter << endl;

      m_workspace_1 = m_Psi;
      MyRealTools::MPI::AssembleSystem_Jacobian<dim>( dof_handler, fe, constraints, m_workspace_1, m_Potential, m_mu[0], m_gs[0], m_system_matrix);
      solve();

      m_Psi.add( -0.1, m_newton_update); 
      constraints.distribute(m_Psi);

      m_workspace_1 = m_Psi;
      MyRealTools::MPI::AssembleRHS_L2gradient<dim>( dof_handler, fe, constraints, m_workspace_1, m_Potential, m_mu[0], m_gs[0], m_res[0], m_system_rhs );
      
      m_resp[0] = m_res_old[0]-m_res[0];
      m_res_old[0] = m_res[0];

      //if( m_counter % m_NA == 0 ) output_results(path);

      columns& cols = m_table.new_line();
      m_table.insert( cols, MyTable::COUNTER, double(m_counter) );
      m_table.insert( cols, MyTable::RES, m_res[0] );
      m_table.insert( cols, MyTable::RESP, m_resp[0] );
      m_table.insert( cols, MyTable::MU, m_mu[0] );
      m_table.insert( cols, MyTable::l2norm_t, l2norm_t() );
      m_table.insert( cols, MyTable::PARTICLE_NUMBER, m_N[0] );

      if( m_root ) cout << m_table;
      if( m_res[0] < m_epsilon[1] ) { retval=Status::SUCCESS; break; }

      m_counter++;
    }while( true );
    
    m_workspace_1 = m_Psi;
    m_N[0] = MyRealTools::MPI::Particle_Number( mpi_communicator, dof_handler, fe, m_workspace_1 );
    
    if( m_N[0] < 1e-5 ) retval = Status::ZERO_SOL;
   
    return retval;
  }
  
  template <int dim>
  void MySolver<dim>::run ()
  {
    int status;

    make_grid_custom();
    setup_system();

    
    map<string,double> constants;
    constants["pi"] = numbers::PI;

    // \Psi_{initial}
    string potential_str = "0.5*(x+3)^2+0.5*y^2"; // change me 
    string guess_str = "exp(-0.5*((x+3)^2+y^2))"; // change me

    FunctionParser<dim> guess_fct;
    guess_fct.initialize( FunctionParser<dim>::default_variable_names(), guess_str , constants );
    VectorTools::interpolate (dof_handler, guess_fct, m_Psi );

    m_workspace_1 = m_Psi;
    m_Psi *= 1.0/sqrt(MyRealTools::MPI::Particle_Number( mpi_communicator, dof_handler, fe, m_workspace_1 ));
 
    MyRealTools::MPI::Compute_mu( mpi_communicator, dof_handler, fe, m_Potential, m_workspace_1, m_mu[0] );

    m_mu[0] += m_gs[0]/fabs(m_gs[0])*m_dmu;
    pcout << "m_mu = " << m_mu[0] << endl;

    m_Potential.initialize( FunctionParser<dim>::default_variable_names(), potential_str , constants );    
      
    status = DoIter();

    constraints.distribute(m_Psi);
    m_workspace_1 = m_Psi;
    m_N[0] = MyRealTools::MPI::Particle_Number( mpi_communicator, dof_handler, fe, m_workspace_1 );
    m_Psi *= 1.0/sqrt(m_N[0]);

    Interpolate_R_to_C( m_Psi_C_ghosted_i );

    // \Psi_{final}
    potential_str = "0.5*(x-3)^2+0.5*y^2"; // change me
    guess_str = "exp(-0.5*((x-3)^2+y^2))"; // change me

    guess_fct.initialize( FunctionParser<dim>::default_variable_names(), guess_str , constants );
    VectorTools::interpolate (dof_handler, guess_fct, m_Psi );

    m_workspace_1 = m_Psi;
    m_Psi *= 1.0/sqrt(MyRealTools::MPI::Particle_Number( mpi_communicator, dof_handler, fe, m_workspace_1 ));

    m_Potential.initialize( FunctionParser<dim>::default_variable_names(), potential_str , constants );    
      
    status = DoIter();

    m_Psi *= 1.0/sqrt(m_N[0]);

    Interpolate_R_to_C( m_Psi_C_ghosted_f );

    // output
    vector<double> newgs = {m_gs[0]*m_N[0]}; 
    m_ph.Set_Physics( "gs_1", newgs );
    m_ph.SaveXMLFile( "params_one.xml" );

    std::vector<const LA::MPI::Vector*> x_system (2);
    x_system[0] = &m_Psi_C_ghosted_i;
    x_system[1] = &m_Psi_C_ghosted_f;
    
    parallel::distributed::SolutionTransfer<dim,LA::MPI::Vector> solution_transfer(dof_handler_2);
    solution_transfer.prepare_for_serialization(x_system);
    triangulation.save( "oct_0.bin" );    
  }

  template <int dim>
  void MySolver<dim>::solve ()
  {
    TimerOutput::Scope timing_section(m_computing_timer, "");
    pcout << "Solving..." << endl;
    
    SolverControl solver_control;
    
    PETScWrappers::SparseDirectMUMPS solver(solver_control, mpi_communicator);
    solver.set_symmetric_mode(false);
    solver.solve(m_system_matrix, m_newton_update, m_system_rhs);
    constraints.distribute (m_newton_update);

    
  }  
  
  template<int dim>
  void MySolver<dim>::Interpolate_R_to_C( LA::MPI::Vector& retval )
  {
    TimerOutput::Scope timing_section(m_computing_timer, "");

    constraints.distribute(m_Psi);
    m_workspace_1=m_Psi;

    MyComplexTools::MPI::Interpolate_R_to_C( mpi_communicator, dof_handler, fe, m_workspace_1, dof_handler_2, fe_2, constraints_2, retval );

        
  }  

  template<int dim>
  void MySolver<dim>::save( string filename )
  {
    constraints.distribute(m_Psi);
    m_workspace_1 = m_Psi;
    parallel::distributed::SolutionTransfer<dim,LA::MPI::Vector> solution_transfer(dof_handler);
    solution_transfer.prepare_for_serialization(m_workspace_1);
    triangulation.save( filename.c_str() );
  }
  
  template<int dim>
  void MySolver<dim>::save_one( string filename )
  {
    m_workspace_1 = m_Psi;
    double tmp = MyRealTools::MPI::Particle_Number( mpi_communicator, dof_handler, fe, m_workspace_1 );
    m_workspace_ng=m_Psi;
    m_workspace_ng*=sqrt(1/tmp);
    m_workspace_1 = m_workspace_ng; 
    parallel::distributed::SolutionTransfer<dim,LA::MPI::Vector> solution_transfer(dof_handler);
    solution_transfer.prepare_for_serialization(m_workspace_1);
    triangulation.save( filename.c_str() );
  }
} // end of namespace 

int main ( int argc, char *argv[] )
{
  using namespace dealii;
  deallog.depth_console (0);

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv );
  {
    BreedSolver::MySolver<DIMENSION> solver("params.xml");
    solver.run ();
  }
return EXIT_SUCCESS;
}
