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

#include <deal.II/lac/generic_linear_algebra.h>

namespace LA
{
  using namespace dealii::LinearAlgebraPETSc;
}

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/constraint_matrix.h>


#include <deal.II/lac/petsc_parallel_sparse_matrix.h>
#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_precondition.h>
#include <deal.II/lac/sparsity_tools.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/intergrid_map.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_q.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>

#include <deal.II/base/utilities.h>
#include <deal.II/base/aligned_vector.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/function_parser.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/data_out_base.h>
#include <deal.II/base/table.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/solution_transfer.h>

#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <algorithm>

#include "global.h"
#include "mpi.h"
#include "MyParameterHandler.h"
#include "my_table.h"
#include "atus2.h"
#include "anyoption.h"

namespace TestPrograms
{
  using namespace std;
  using namespace dealii;

  template <int dim>
  class MySolver
  {
  public:
    MySolver( const std::string& );
    ~MySolver();

    void run ( string, string );

  protected:
    void make_grid();
    void setup_system();
    void setup_boundary_ids();

    MyParameterHandler m_ph;
    MPI_Comm mpi_communicator;
    parallel::distributed::Triangulation<dim> triangulation;
    FE_Q<dim> fe;
    DoFHandler<dim> dof_handler;
    IndexSet locally_owned_dofs;
    IndexSet locally_relevant_dofs;
    ConstraintMatrix constraints;

    LA::MPI::Vector m_Psi_1;
    LA::MPI::Vector m_Psi_2;

    ConditionalOStream pcout;

    bool m_root;

    double m_xmin, m_xmax;
    double m_ymin, m_ymax;
    double m_zmin, m_zmax;
    double m_mu;
    double m_gs;
    vector<double> m_omega;
    
    int m_rank;
  };

  template <int dim>
  MySolver<dim>::MySolver ( const std::string& xmlfilename ) 
    : 
    m_ph(xmlfilename),
    mpi_communicator (MPI_COMM_WORLD),
    triangulation (mpi_communicator, typename Triangulation<dim>::MeshSmoothing(Triangulation<dim>::smoothing_on_refinement|Triangulation<dim>::smoothing_on_coarsening)),
    fe (gl_degree_fe),
    dof_handler (triangulation),
    m_root(Utilities::MPI::this_mpi_process(mpi_communicator) == 0),
    pcout (cout, m_root)
  {
    m_xmin = m_ph.Get_Mesh("xrange",0);
    m_xmax = m_ph.Get_Mesh("xrange",1);
    m_ymin = m_ph.Get_Mesh("yrange",0);
    m_ymax = m_ph.Get_Mesh("yrange",1);
    m_zmin = m_ph.Get_Mesh("zrange",0);
    m_zmax = m_ph.Get_Mesh("zrange",1);
  
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
  }
  
  template <int dim>
  void MySolver<dim>::setup_system()
  {
    dof_handler.distribute_dofs (fe);

    locally_owned_dofs = dof_handler.locally_owned_dofs ();
    DoFTools::extract_locally_relevant_dofs (dof_handler, locally_relevant_dofs);

    m_Psi_1.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    m_Psi_2.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);

    constraints.clear ();
    constraints.reinit (locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints (dof_handler, constraints);
    constraints.close ();
  }

  template <int dim>
  void MySolver<dim>::run( string filename, string filename2 )
  {
    make_grid();
    triangulation.load( filename.c_str() );
    
    setup_system();

    LA::MPI::Vector tmp, tmp2;
    tmp.reinit (locally_owned_dofs, mpi_communicator);

    parallel::distributed::SolutionTransfer<dim,LA::MPI::Vector> solution_transfer(dof_handler);
    solution_transfer.deserialize(tmp);
    constraints.distribute(tmp);
    m_Psi_1=tmp;

    /// 

    cout << "HERE1" << endl;
    parallel::distributed::Triangulation<dim> triangulation_2(mpi_communicator, typename Triangulation<dim>::MeshSmoothing(Triangulation<dim>::smoothing_on_refinement|Triangulation<dim>::smoothing_on_coarsening));

    Point<dim,double> pt1;
    Point<dim,double> pt2;

    double min[] = {m_xmin, m_ymin, m_zmin};
    double max[] = {m_xmax, m_ymax, m_zmax};

    for( int i=0; i<dim; i++ )
    {
      pt1(i) = min[i];
      pt2(i) = max[i];
    }

    GridGenerator::hyper_rectangle(triangulation_2, pt2, pt1);

    triangulation_2.load( filename2.c_str() );
    DoFHandler<dim> dof_handler_2(triangulation_2);
    dof_handler_2.distribute_dofs (fe);

    IndexSet locally_owned_dofs_2;
    IndexSet locally_relevant_dofs_2;
    ConstraintMatrix constraints_2;

    locally_owned_dofs_2 = dof_handler_2.locally_owned_dofs ();
    DoFTools::extract_locally_relevant_dofs (dof_handler_2, locally_relevant_dofs_2);

    tmp.reinit (locally_owned_dofs_2, mpi_communicator);
    tmp2.reinit (locally_owned_dofs_2, locally_relevant_dofs_2, mpi_communicator);

    constraints_2.clear ();
    constraints_2.reinit (locally_relevant_dofs_2);
    DoFTools::make_hanging_node_constraints (dof_handler_2, constraints_2);
    constraints_2.close ();

    cout << "HERE1c" << endl;

    parallel::distributed::SolutionTransfer<dim,LA::MPI::Vector> solution_transfer_2(dof_handler_2);
    solution_transfer_2.deserialize(tmp);
    constraints_2.distribute(tmp);
    tmp2=tmp;

    ///
    cout << "HERE2" << endl;
    InterGridMap<DoFHandler<dim> > grid_2_to_1_map;
    grid_2_to_1_map.make_mapping (dof_handler_2,dof_handler);    

    VectorTools::interpolate_to_different_mesh( grid_2_to_1_map, tmp2, constraints, m_Psi_2 );

    DataOut<dim> data_out;
    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (m_Psi_1, "Psi_1" );
    data_out.add_data_vector (m_Psi_2, "Psi_2" );
    data_out.build_patches ();
    data_out.write_vtu_in_parallel ( "delme.vtu", mpi_communicator );    
  }
} // end of namespace 

int main ( int argc, char *argv[] )
{
  using namespace dealii;
  deallog.depth_console (0);

  std::string filename, filename2;

  AnyOption * opt = new AnyOption();
  int dim=0;

  opt->addUsage( "" );
  opt->addUsage( "Usage: binR_to_atus2 [options] filename" );
  opt->addUsage( "" );
  opt->addUsage( " --help -h   Prints this help " );
  opt->addUsage( " --dim       2 or 3" );
  opt->addUsage( "" );
  opt->setFlag(  "help", 'h' );   
  opt->setOption( "dim" );   

  opt->processCommandArgs( argc, argv );

  if( opt->getFlag( "help" ) || opt->getFlag( 'h' ) ) opt->printUsage();

  if( opt->getValue("dim") == nullptr ) 
  {
    opt->printUsage();
    delete opt;     
    return EXIT_FAILURE;
  }

  dim = atoi(opt->getValue("dim"));
  
  if( !(dim == 2 || dim == 3) ) 
  {
    opt->printUsage();
    delete opt;     
    return EXIT_FAILURE;
  }
 
  if( opt->getArgc() > 0 ) {
    filename = opt->getArgv(0);
    filename2 = opt->getArgv(1);}
  else opt->printUsage();
  delete opt; 

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);
  {
    if( dim == 2 )
    {
      TestPrograms::MySolver<2> solver("params.xml");
      solver.run(filename,filename2);
    }
    if( dim == 3 )
    {
      TestPrograms::MySolver<3> solver("params.xml");
      solver.run(filename,filename2);
    }
  }  
return EXIT_SUCCESS;
}