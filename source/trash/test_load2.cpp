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
#include <deal.II/numerics/derivative_approximation.h>
#include <deal.II/numerics/fe_field_function.h>
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
#include <iomanip>
#include <stdlib.h>
#include <locale>
#include <limits>
#include <cmath>
#include <array>
#include <unistd.h>

#include "mpi.h"
#include "functions.h"
#include "SinCos.h"
#include "ParameterReader.h"
#include "my_table.h"
#include "ref_pt_list.h"

#define STR1(x) #x
#define STR2(x) STR1(x)

namespace BreedSolver
{
  using namespace std;
  using namespace dealii;

  template <int dim>
  class MySolver
  {
  public:
    MySolver ();
    virtual ~MySolver();

    void run ();
  protected:
    MPI_Comm mpi_communicator;

    void make_grid ();
    void setup_system ();
    void save ();
    void output_results ();

    ConditionalOStream pcout;
    parallel::distributed::Triangulation<dim> triangulation;
    FE_Q<dim> fe;
    DoFHandler<dim> dof_handler;
    IndexSet locally_owned_dofs;
    IndexSet locally_relevant_dofs;
    ConstraintMatrix constraints;

    LA::MPI::Vector m_vec;
    LA::MPI::Vector m_vec_ng; // no ghosts
  };

  template <int dim>
  MySolver<dim>::MySolver () 
    : 
    mpi_communicator(MPI_COMM_WORLD), 
    pcout (cout, (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)),
    triangulation (mpi_communicator, typename Triangulation<dim>::MeshSmoothing(Triangulation<dim>::limit_level_difference_at_vertices|Triangulation<dim>::eliminate_refined_inner_islands|Triangulation<dim>::smoothing_on_refinement|Triangulation<dim>::smoothing_on_coarsening)),
    fe (2),
    dof_handler (triangulation)
  {
  }

  template <int dim>
  MySolver<dim>::~MySolver ()
  {
    dof_handler.clear ();
  }

  /*
  template<int dim>
  void MySolver<dim>::load ()
  {
    make_grid();
    triangulation.load( filename.c_str() );

    setup_system(true);

    parallel::distributed::SolutionTransfer<dim,LA::MPI::Vector> solution_transfer(dof_handler);
    solution_transfer.deserialize(m_Psi_vgl);
    constraints.distribute(m_Psi_vgl);
  }    
  */

  template <int dim>
  void MySolver<dim>::make_grid ()
  {
    Point<dim,double> pt1(-3,-3); 
    Point<dim,double> pt2(3,3);
    
    GridGenerator::hyper_rectangle(triangulation, pt2, pt1);
    triangulation.refine_global(7);
  }
  
  template <int dim>
  void MySolver<dim>::setup_system ()
  {
    dof_handler.distribute_dofs (fe);
      
    locally_owned_dofs = dof_handler.locally_owned_dofs ();
    DoFTools::extract_locally_relevant_dofs (dof_handler, locally_relevant_dofs);

    m_vec.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    m_vec_ng.reinit (locally_owned_dofs, mpi_communicator);
    
    int myrank;
    MPI_Comm_rank( mpi_communicator, &myrank );    
    cout << "(" << myrank << ") locally_owned_dofs = " << m_vec_ng.local_size()  << endl;
    
    vector<bool> mask (dof_handler.get_fe().n_components(), true );
    constraints.clear ();
    constraints.reinit (locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints (dof_handler, constraints);
    VectorTools::interpolate_boundary_values (dof_handler, 0, ZeroFunction<dim>(), constraints, ComponentMask(mask));
    constraints.close ();
  }

  template <int dim>
  void MySolver<dim>::output_results ()
  {
    Vector<float> subdomain (triangulation.n_active_cells());
    for (unsigned int i=0; i<subdomain.size(); ++i)
      subdomain(i) = triangulation.locally_owned_subdomain();

    DataOut<dim> data_out;
    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector ( m_vec, "Psi");
    data_out.add_data_vector ( subdomain, "subdomain");    
    data_out.build_patches ();
    data_out.write_vtu_in_parallel ( "out-1.vtu", mpi_communicator);
  }

  template <int dim>
  void MySolver<dim>::run()
  {
    make_grid();
    setup_system();

    string fct_str = "x*exp(-x*x-y*y)";
    map<string,double> constants;
    constants["pi"] = numbers::PI;
    FunctionParser<dim> guess_fct;
    guess_fct.initialize( FunctionParser<dim>::default_variable_names(), fct_str, constants );
    VectorTools::interpolate (dof_handler, guess_fct, m_vec );
    constraints.distribute(m_vec);
    
    output_results();
    save();
  }
  
  template<int dim>
  void MySolver<dim>::save ()
  {
    parallel::distributed::SolutionTransfer<dim,LA::MPI::Vector> solution_transfer(dof_handler);
    solution_transfer.prepare_serialization(m_vec);
    triangulation.save( "out-1.bin" );
  }
} // end of namespace 

int main ( int argc, char *argv[] )
{
  using namespace dealii;

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv );
  {
    BreedSolver::MySolver<2> solver;
    solver.run ();
  }
return EXIT_SUCCESS;
}