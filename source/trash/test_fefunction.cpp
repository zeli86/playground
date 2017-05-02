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


#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/parameter_handler.h>

#include <deal.II/lac/generic_linear_algebra.h>

namespace LA
{
  using namespace dealii::LinearAlgebraPETSc;
}

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/compressed_simple_sparsity_pattern.h>

#include <deal.II/lac/petsc_parallel_sparse_matrix.h>
#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_precondition.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
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
#include <deal.II/numerics/fe_field_function.h>

#include <deal.II/base/utilities.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/function_parser.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/solution_transfer.h>


#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <iomanip>

#include "mpi.h"
#include "functions.h"

namespace realtime_propagation
{
  using namespace std;
  using namespace dealii;

  template <int dim>
  class MySolver
  {
  public:
    MySolver();
    ~MySolver();

    void run ();
    
  protected:
    void make_grid();
    void setup_system();
    void output_results ( string );
    void assemble_system();
    void solve();

    LA::MPI::SparseMatrix system_matrix;
    LA::MPI::Vector system_rhs;
    MPI_Comm mpi_communicator;
    parallel::distributed::Triangulation<dim> triangulation;
    FESystem<dim> fe;
    DoFHandler<dim> dof_handler;
    IndexSet locally_owned_dofs;
    IndexSet locally_relevant_dofs;
    ConstraintMatrix constraints;

    LA::MPI::Vector m_Psi_1; 
    LA::MPI::Vector m_Psi_2; 

    ConditionalOStream pcout;
    ofstream m_computing_timer_log;
    TimerOutput m_computing_timer;    

    bool m_root;

    double m_xmin;
    double m_xmax;
    double m_ymin;
    double m_ymax;
    double m_zmin;
    double m_zmax;

    unsigned int m_NA;
    unsigned int m_NK;
    unsigned int m_global_refinement;
  };

  template <int dim>
  MySolver<dim>::MySolver () 
    : 
    mpi_communicator (MPI_COMM_WORLD),
    triangulation (mpi_communicator, typename Triangulation<dim>::MeshSmoothing(Triangulation<dim>::smoothing_on_refinement|Triangulation<dim>::smoothing_on_coarsening)),
    fe (FE_Q<dim>(1), 2),
    dof_handler (triangulation),
    pcout (cout, (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)),
    m_computing_timer_log("benchmark.txt"),
    m_computing_timer(mpi_communicator, m_computing_timer_log, TimerOutput::summary, TimerOutput:: cpu_and_wall_times )
  {
    m_global_refinement = 8;
    m_xmin = -5.0;
    m_xmax = -m_xmin;
    m_ymin = -5.0;
    m_ymax = -m_ymin;
    m_zmin = -5.0;
    m_zmax = -m_zmin;

    m_root = (Utilities::MPI::this_mpi_process(mpi_communicator) == 0);
  }

  template <int dim>
  MySolver<dim>::~MySolver ()
  {
    dof_handler.clear ();
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
    triangulation.refine_global(m_global_refinement);

//     unsigned int tmp1[2], tmp2[2];
//     tmp1[0] = triangulation.n_cells();
//     tmp1[1] = triangulation.n_active_cells();
//     MPI_Allreduce( tmp1, tmp2, 2, MPI_UNSIGNED, MPI_SUM, mpi_communicator);

    m_computing_timer.exit_section();
  }
  
  template <int dim>
  void MySolver<dim>::setup_system()
  {
    m_computing_timer.enter_section(__func__);
    dof_handler.distribute_dofs (fe);

    locally_owned_dofs = dof_handler.locally_owned_dofs ();
    DoFTools::extract_locally_relevant_dofs (dof_handler, locally_relevant_dofs);

    system_rhs.reinit(locally_owned_dofs, mpi_communicator);
    m_Psi_1.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    m_Psi_2.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);

    vector<bool> mask (dof_handler.get_fe().n_components(), true );    
    
    constraints.clear ();
    constraints.reinit (locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints (dof_handler, constraints);
    VectorTools::interpolate_boundary_values (dof_handler, 0, ZeroFunction<dim>(2), constraints, ComponentMask(mask));
    constraints.close ();
    
    CompressedSimpleSparsityPattern csp (locally_relevant_dofs);

    DoFTools::make_sparsity_pattern (dof_handler, csp, constraints, false);
    SparsityTools::distribute_sparsity_pattern (csp, dof_handler.n_locally_owned_dofs_per_processor(), mpi_communicator, locally_relevant_dofs);

    system_matrix.reinit (locally_owned_dofs, locally_owned_dofs, csp, mpi_communicator);

    m_computing_timer.exit_section();
  }  
  
  template <int dim>
  void MySolver<dim>::output_results ( string path ) 
  {
    m_computing_timer.enter_section(__func__);
    string filename;
    
    vector<std::string> solution_names;

    Vector<float> subdomain (triangulation.n_active_cells());
    for (unsigned int i=0; i<subdomain.size(); ++i)
      subdomain(i) = triangulation.locally_owned_subdomain();

    m_Psi_1.update_ghost_values();
    m_Psi_2.update_ghost_values();
    
    DataOut<dim> data_out;
    data_out.attach_dof_handler (dof_handler);
    solution_names.push_back ("Re Psi_1");
    solution_names.push_back ("Im Psi_1");    
    data_out.add_data_vector (m_Psi_1, solution_names);
    solution_names.clear();
    solution_names.push_back ("Re Psi_2");
    solution_names.push_back ("Im Psi_2");    
    data_out.add_data_vector (m_Psi_2, solution_names);
    data_out.add_data_vector (subdomain, "subdomain");
    data_out.build_patches ();
    
    filename = path + "out.vtu";
    data_out.write_vtu_in_parallel ( filename.c_str(), mpi_communicator );

    m_computing_timer.exit_section();
  }
  
  template <int dim>
  void MySolver<dim>::assemble_system ()
  {
    m_computing_timer.enter_section(__func__);
    const QGauss<dim> quadrature_formula(fe.degree+1);

    const FEValuesExtractors::Scalar rt (0);
    const FEValuesExtractors::Scalar it (1);
        
    
    const double a = 1;
    const double b = 1;
    
    const double a1 = a/(a*a+b*b);
    const double b1 = b/(a*a+b*b);
    
    system_matrix = 0;
    system_rhs = 0;

    FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_gradients|update_quadrature_points|update_JxW_values);

    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    const unsigned int n_q_points    = quadrature_formula.size();

    FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);
    Vector<double> cell_rhs (dofs_per_cell);

    vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    vector<Vector<double>> Psi(n_q_points,Vector<double>(2));
 
    double JxW;
    
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
      if( cell->is_locally_owned() )
      {
        cell_matrix = 0;
        cell_rhs = 0;

        fe_values.reinit (cell);
        fe_values.get_function_values(m_Psi_1, Psi);
	
        for( unsigned int qp=0; qp<n_q_points; qp++ )
        {
	  JxW = fe_values.JxW(qp);
	  for (unsigned int i=0; i<dofs_per_cell; i++ )
	  {
	    for (unsigned int j=0; j<dofs_per_cell; j++ )
	    {
               cell_matrix(i,j) += JxW*(-b1*fe_values[rt].value(i,qp)*fe_values[rt].value(j,qp) + a1*fe_values[rt].value(i,qp)*fe_values[it].value(j,qp) - b1*fe_values[it].value(i,qp)*fe_values[it].value(j,qp) - a1*fe_values[it].value(i,qp)*fe_values[rt].value(j,qp)); 
               //cell_matrix(i,j) += JxW*(a1*fe_values[rt].value(i,qp)*fe_values[rt].value(j,qp) + b1*fe_values[rt].value(i,qp)*fe_values[it].value(j,qp) + a1*fe_values[it].value(i,qp)*fe_values[it].value(j,qp) - b1*fe_values[it].value(i,qp)*fe_values[rt].value(j,qp)); 
            }
            cell_rhs(i) += JxW* (Psi[qp][0]*fe_values[rt].value(i,qp) + Psi[qp][1]*fe_values[it].value(i,qp));
          }
        }
        cell->get_dof_indices (local_dof_indices);
        constraints.distribute_local_to_global(cell_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs);
      }
    }
    system_matrix.compress(VectorOperation::add);
    system_rhs.compress(VectorOperation::add);
    m_computing_timer.exit_section();
  }

  template <int dim>
  void MySolver<dim>::solve ()
  {
    m_computing_timer.enter_section(__func__);
    LA::MPI::Vector tmp_vec(locally_owned_dofs, mpi_communicator);
    SolverControl solver_control;
    
    PETScWrappers::SparseDirectMUMPS solver(solver_control, mpi_communicator);
    solver.set_symmetric_mode(false);
    solver.solve(system_matrix, tmp_vec, system_rhs);

    constraints.distribute (tmp_vec);
    m_Psi_2 = tmp_vec;
    m_computing_timer.exit_section();
  }
  
  template <int dim>
  void MySolver<dim>::run ()
  {
    double T, N, W;

    make_grid();
    setup_system();
    
    std::map<std::string,double> constants;
    constants["pi"] = numbers::PI;    
    // Define the variables that will be used inside the expressions
    std::string variables = "x,y";
    // Define the expressions of the individual components of a
    // vector valued function with two components:
    std::vector<std::string> expressions(2);
    expressions[0] = "exp(-0.5*(x^2+y^2))";
    expressions[1] = "0"; //"x*exp(-0.5*(x^2+y^2))";
    // function parser with 3 variables and 2 components
    FunctionParser<dim> vector_function(2);
    // And populate it with the newly created objects.
    vector_function.initialize(variables, expressions, constants );
    
    VectorTools::interpolate (dof_handler, vector_function, m_Psi_1 );        
    constraints.distribute(m_Psi_1);
    m_Psi_1.update_ghost_values();

    
    assemble_system();
    solve();
    
    output_results("");
  }
} // end of namespace 

int main ( int argc, char *argv[] )
{
  int myrank;

  using namespace dealii;
  deallog.depth_console (0);

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv );
  {
    MPI_Comm_rank( MPI_COMM_WORLD, &myrank );
    realtime_propagation::MySolver<DIMENSION> solver;
    solver.run();
  }
  
return EXIT_SUCCESS;
}
