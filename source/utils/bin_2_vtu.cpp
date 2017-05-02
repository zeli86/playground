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
 * Purpose: convert deal.ii native file format files into vtu files
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
#include <deal.II/lac/dynamic_sparsity_pattern.h>

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
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/solution_transfer.h>

#include <fstream>
#include <iostream>
#include <regex>
#include <stdlib.h>

#include "mpi.h"
#include "ParameterReader.h"
#include "global.h"

namespace HelperPrograms
{
  enum Status { SUCCESS, FAILED };

  using namespace std;
  using namespace dealii;

  #include "utils_complex.h"

  template <int dim>
  class MySolver
  {
  public:
    MySolver( ParameterHandler & );
    ~MySolver();

    void run ( string );

  protected:
    void make_grid();
    void setup_system();

    void output_results ( string );
    
    void load( string );
    void load_2( string );
    void save( string );

    ParameterHandler &m_ph;
    MPI_Comm mpi_communicator;
    parallel::distributed::Triangulation<dim> triangulation;
    FESystem<dim> fe;
    DoFHandler<dim> dof_handler;
    IndexSet locally_owned_dofs;
    IndexSet locally_relevant_dofs;
    ConstraintMatrix constraints;

    LA::MPI::Vector m_Psi;
    LA::MPI::Vector m_Psi_0;

    ConditionalOStream pcout;

    bool m_root;

    double m_xmin;
    double m_xmax;
    double m_ymin;
    double m_ymax;
    double m_zmin;
    double m_zmax;

    int m_rank;
    unsigned int m_NA;
    unsigned int m_NK;
    unsigned int m_global_refinement;
    unsigned int m_QN1[3];
  };

  template <int dim>
  MySolver<dim>::MySolver ( ParameterHandler &ph ) 
    : 
    m_ph(ph),
    mpi_communicator (MPI_COMM_WORLD),
    triangulation (mpi_communicator, typename Triangulation<dim>::MeshSmoothing(Triangulation<dim>::smoothing_on_refinement|Triangulation<dim>::smoothing_on_coarsening)),
    fe (FE_Q<dim>(gl_degree_fe), 2),
    dof_handler (triangulation),
    pcout (cout, (Utilities::MPI::this_mpi_process(mpi_communicator) == 0))
  {
    m_ph.enter_subsection("Mesh & geometry parameters");    
    m_global_refinement = m_ph.get_integer("global_refinements");
    m_xmin = m_ph.get_double("xMin");
    m_xmax = m_ph.get_double("xMax");
    m_ymin = m_ph.get_double("yMin");
    m_ymax = m_ph.get_double("yMax");
    m_zmin = m_ph.get_double("zMin");
    m_zmax = m_ph.get_double("zMax");
    m_ph.leave_subsection();

    m_root = (Utilities::MPI::this_mpi_process(mpi_communicator) == 0);
    MPI_Comm_rank(mpi_communicator, &m_rank);
  }

  template <int dim>
  MySolver<dim>::~MySolver ()
  {
    dof_handler.clear ();
  }

  template <int dim>
  void MySolver<dim>::make_grid ()
  {
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
  }
  
  template <int dim>
  void MySolver<dim>::setup_system()
  {
    dof_handler.distribute_dofs (fe);

    locally_owned_dofs = dof_handler.locally_owned_dofs ();
    DoFTools::extract_locally_relevant_dofs (dof_handler, locally_relevant_dofs);

    m_Psi.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    m_Psi_0.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
  }

  template <int dim>
  void MySolver<dim>::output_results ( string filename ) 
  {
    ComputeIntensity<dim> intensities("Intensity Psi");
    ComputeIntensity<dim> intensities2("Intensity Psi_0");
    //ComputePhase<dim> phase;

    vector<std::string> solution_names;

    DataOut<dim> data_out;
    data_out.attach_dof_handler (dof_handler);
    solution_names.push_back ("Re Psi");
    solution_names.push_back ("Im Psi");    
    data_out.add_data_vector (m_Psi, solution_names);
    data_out.add_data_vector (m_Psi, intensities);
    
    solution_names.clear();
    solution_names.push_back ("Re Psi_0");
    solution_names.push_back ("Im Psi_0");    
    data_out.add_data_vector (m_Psi_0, solution_names);
    data_out.add_data_vector (m_Psi_0, intensities2);
    data_out.build_patches ();
    data_out.write_vtu_in_parallel ( filename.c_str(), mpi_communicator );
  }
  
  template <int dim>
  void MySolver<dim>::run( string filename )
  {
    load_2( filename );
    string filename2 = filename;
    filename2.erase( filename.length()-3, 3 );
    filename2 += "vtu";
    output_results( filename2 );
  }

  template<int dim>
  void MySolver<dim>::load_2( string filename )
  {
    make_grid();
    triangulation.load( filename.c_str() );

    setup_system();

    std::vector<LA::MPI::Vector*> x_system (2);
    x_system[0] = &m_Psi;
    x_system[1] = &m_Psi_0;

    parallel::distributed::SolutionTransfer<dim,LA::MPI::Vector> solution_transfer(dof_handler);
    solution_transfer.deserialize(x_system);
    
    LA::MPI::Vector tmp;
    tmp.reinit (locally_owned_dofs, mpi_communicator);
    
    tmp=m_Psi;
    constraints.distribute(tmp);
    m_Psi=tmp;

    tmp=m_Psi_0;
    constraints.distribute(tmp);
    m_Psi_0=tmp;
  }
} // end of namespace 

int main ( int argc, char *argv[] )
{
  int myrank;

  using namespace dealii;
  deallog.depth_console (0);

  ParameterHandler  prm;
  ParameterReader   param(prm);

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv );
  {
    MPI_Comm_rank( MPI_COMM_WORLD, &myrank );
    if( !param.read_parameters("params.prm") )
    {
      MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
      return EXIT_FAILURE;
    }

    std::string filename = argv[1];
    HelperPrograms::MySolver<DIMENSION> solver(prm);
    solver.run(filename);
  }  
return EXIT_SUCCESS;
}
