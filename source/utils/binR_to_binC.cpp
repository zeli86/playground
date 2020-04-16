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

#include <deal.II/lac/generic_linear_algebra.h>

namespace LA
{
  using namespace dealii::LinearAlgebraPETSc;
}

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>

#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>
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
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/index_set.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/solution_transfer.h>

#include <iostream>

#include "global.h"
#include "mpi.h"
#include "MyParameterHandler.h"
#include "MyComplexTools.h"
#include "muParser.h"
#include "cxxopts.hpp"

namespace HelperPrograms
{
  using namespace std;
  using namespace dealii;

  template <int dim>
  class binR_to_binC
  {
  public:
    binR_to_binC( MyParameterHandler &, const string );
    ~binR_to_binC();

    void run ();

  protected:
    void make_grid();
    void setup_system();
    void Interpolate_R_to_C( string );
    void load( const string& );

    MyParameterHandler & m_ph;
    MPI_Comm mpi_communicator;    
    parallel::distributed::Triangulation<dim> triangulation;
    FE_Q<dim> fe;
    DoFHandler<dim> dof_handler;
    FESystem<dim> fe_2;
    DoFHandler<dim> dof_handler_2;
    IndexSet locally_owned_dofs, locally_owned_dofs_2;
    IndexSet locally_relevant_dofs, locally_relevant_dofs_2;
    AffineConstraints<double> constraints, constraints_2;

    string m_bin_filename;

    LA::MPI::Vector m_Psi_C_ghosted;
    LA::MPI::Vector m_Psi_R_ghosted;

    bool m_root;

    double m_xmin, m_xmax;
    double m_ymin, m_ymax;
    double m_zmin, m_zmax;   
  };

/*************************************************************************************************/

/**
 * Constructor
 */
  template <int dim>
  binR_to_binC<dim>::binR_to_binC ( MyParameterHandler & ph, const string fn ) 
    : 
    m_ph(ph),
    mpi_communicator (MPI_COMM_WORLD),
    triangulation (mpi_communicator, typename Triangulation<dim>::MeshSmoothing(Triangulation<dim>::limit_level_difference_at_vertices|Triangulation<dim>::eliminate_refined_inner_islands|Triangulation<dim>::smoothing_on_refinement|Triangulation<dim>::smoothing_on_coarsening)),
    fe (gl_degree_fe),
    dof_handler (triangulation),
    fe_2 (FE_Q<dim>(gl_degree_fe),2),
    dof_handler_2 (triangulation),
    m_bin_filename(fn)
  {
     m_xmin = m_ph.Get_Mesh("xrange",0);
     m_xmax = m_ph.Get_Mesh("xrange",1);
     m_ymin = m_ph.Get_Mesh("yrange",0);
     m_ymax = m_ph.Get_Mesh("yrange",1);
     m_zmin = m_ph.Get_Mesh("zrange",0);
     m_zmax = m_ph.Get_Mesh("zrange",1);   
  }

  template <int dim>
  binR_to_binC<dim>::~binR_to_binC ()
  {
    dof_handler.clear ();
    dof_handler_2.clear ();
  }
  
  template <int dim>
  void binR_to_binC<dim>::make_grid ()
  {
    Point<dim,double> pt1, pt2;

    double min[] = {m_xmin, m_ymin, m_zmin};
    double max[] = {m_xmax, m_ymax, m_zmax};

    for( int i=0; i<dim; ++i )
    {
      pt1(i) = min[i];
      pt2(i) = max[i];
    }

    GridGenerator::hyper_rectangle(triangulation, pt2, pt1);
  }

  template <int dim>
  void binR_to_binC<dim>::setup_system()
  {
    // stuff for the first dof handler
    dof_handler.distribute_dofs (fe);

    locally_owned_dofs = dof_handler.locally_owned_dofs ();
    DoFTools::extract_locally_relevant_dofs (dof_handler, locally_relevant_dofs);
    
    m_Psi_R_ghosted.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
   
    constraints.clear ();
    constraints.reinit (locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints (dof_handler, constraints);
    VectorTools::interpolate_boundary_values (dof_handler, 0, ZeroFunction<dim>(), constraints);
    constraints.close ();

    // stuff for the second dof handler
    dof_handler_2.distribute_dofs (fe_2);

    locally_owned_dofs_2 = dof_handler_2.locally_owned_dofs ();
    DoFTools::extract_locally_relevant_dofs (dof_handler_2, locally_relevant_dofs_2);

    m_Psi_C_ghosted.reinit (locally_owned_dofs_2, locally_relevant_dofs_2, mpi_communicator);

    constraints_2.clear ();
    constraints_2.reinit (locally_relevant_dofs_2);
    DoFTools::make_hanging_node_constraints (dof_handler_2, constraints_2);
    VectorTools::interpolate_boundary_values (dof_handler_2, 0, ZeroFunction<dim>(2), constraints_2);
    constraints_2.close ();
  }

  template <int dim>
  void binR_to_binC<dim>::run()
  {
    load(m_bin_filename);
    setup_system();
    Interpolate_R_to_C("C" + m_bin_filename);
  }
  
  template<int dim>
  void binR_to_binC<dim>::Interpolate_R_to_C( string filename )
  {
    MyComplexTools::MPI::Interpolate_R_to_C( mpi_communicator, dof_handler, fe, m_Psi_R_ghosted, dof_handler_2, fe_2, constraints_2, m_Psi_C_ghosted );

    parallel::distributed::SolutionTransfer<dim,LA::MPI::Vector> solution_transfer(dof_handler_2);
    solution_transfer.prepare_serialization(m_Psi_C_ghosted);
    triangulation.save( filename.c_str() );
  }  

  template<int dim>
  void binR_to_binC<dim>::load( const string& filename )
  {
    make_grid();
    triangulation.load( filename.c_str() );
    
    setup_system();

    LA::MPI::Vector tmp;
    tmp.reinit (locally_owned_dofs, mpi_communicator);

    parallel::distributed::SolutionTransfer<dim,LA::MPI::Vector> solution_transfer(dof_handler);
    solution_transfer.deserialize(tmp);
    constraints.distribute(tmp);
    m_Psi_R_ghosted=tmp;
  }
} // end of namespace 

int main ( int argc, char *argv[] )
{
  cxxopts::Options options("binR_to_binC", "Converts real deal.ii binary format to complex deal.ii binary format.");
  
  options.add_options()
  ("p,params",  "parameter xml file" , cxxopts::value<std::string>()->default_value("params.xml") )
  ("positional", "Positional arguments: these are the arguments that are entered without an option", cxxopts::value<std::vector<std::string>>())
  ("help","Print help")
  ;
  
  options.parse_positional({"positional"});
  auto result = options.parse(argc, argv);

  std::string bin_filename, params_filename;
  try
  {
   if (result.count("") == 0)
    {
      std::cout << options.help({""}) << std::endl;
      return EXIT_FAILURE;
    }

    if( result.count("positional") > 0 )
    {
      bin_filename = result["positional"].as<std::vector<std::string>>()[0]; 
    }
    else
    {        
      std::cout << options.help({""}) << std::endl;
      return EXIT_FAILURE;
    }
    params_filename = result["p"].as<std::string>();
  }
  catch (const cxxopts::OptionException& e)
  {
    std::cout << "error parsing options: " << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  MyParameterHandler params(params_filename);
  int dim=0;

  try
  {
    dim = int(params.Get_Mesh("DIM",0));
  }
  catch (mu::Parser::exception_type &e)
  {
    std::cout << "Message:  " << e.GetMsg() << "\n";
    std::cout << "Formula:  " << e.GetExpr() << "\n";
    std::cout << "Token:    " << e.GetToken() << "\n";
    std::cout << "Position: " << e.GetPos() << "\n";
    std::cout << "Errc:     " << e.GetCode() << "\n";
  }

  dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv );
  {
    switch(dim)
    {
      case 2: { HelperPrograms::binR_to_binC<2> solver( params, bin_filename );
                solver.run();
                break; }
      case 3: { HelperPrograms::binR_to_binC<3> solver( params, bin_filename );
                solver.run();
                break; }
      default: std::cout << "You have found a new dimension!" << std::endl;
    }
  }
return EXIT_SUCCESS;
}
