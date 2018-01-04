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
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>

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
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/utilities.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
//#include <deal.II/tria.h>

#include <fstream>
#include <iostream>

#include "global.h"
#include "MyParameterHandler.h"
#include "MyComplexTools.h"
#include "cxxopts.hpp"
#include "muParser.h"

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

    MyParameterHandler &m_ph;
    Triangulation<dim> triangulation;
    FE_Q<dim> fe;
    DoFHandler<dim> dof_handler;
    FESystem<dim> fe_2;
    DoFHandler<dim> dof_handler_2;
    ConstraintMatrix constraints, constraints_2;

    string m_bin_filename;

    Vector<double> m_Psi_C;
    Vector<double> m_Psi_R;

    double m_xmin, m_xmax;
    //double m_ymin, m_ymax;
    unsigned m_global_refinement;
  };

  /*************************************************************************************************/

  /**
   * Constructor
   */
  template <int dim>
  binR_to_binC<dim>::binR_to_binC ( MyParameterHandler &ph, const string fn )
    :
    m_ph(ph),
    triangulation ( typename Triangulation<dim>::MeshSmoothing(Triangulation<dim>::limit_level_difference_at_vertices | Triangulation<dim>::eliminate_refined_inner_islands | Triangulation<dim>::smoothing_on_refinement | Triangulation<dim>::smoothing_on_coarsening)),
    fe (gl_degree_fe),
    dof_handler (triangulation),
    fe_2 (FE_Q<dim>(gl_degree_fe), 2),
    dof_handler_2 (triangulation),
    m_bin_filename(fn)
  {
    m_xmin = m_ph.Get_Mesh("xrange", 0);
    m_xmax = m_ph.Get_Mesh("xrange", 1);
    //m_ymin = m_ph.Get_Mesh("yrange",0);
    //m_ymax = m_ph.Get_Mesh("yrange",1);
    m_global_refinement = unsigned(m_ph.Get_Mesh("global_refinements", 0));
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
    Point<dim, double> pt1, pt2;

    double min[] = {m_xmin};
    double max[] = {m_xmax};

    for ( int i = 0; i < dim; i++ )
    {
      pt1(i) = min[i];
      pt2(i) = max[i];
    }

    GridGenerator::hyper_rectangle(triangulation, pt2, pt1);
    triangulation.refine_global(m_global_refinement);
  }

  template <int dim>
  void binR_to_binC<dim>::setup_system()
  {
    // stuff for the first dof handler
    dof_handler.distribute_dofs (fe);

    m_Psi_R.reinit (dof_handler.n_dofs());

    constraints.clear ();
    DoFTools::make_hanging_node_constraints (dof_handler, constraints);
    VectorTools::interpolate_boundary_values (dof_handler, 0, ZeroFunction<dim>(), constraints);
    constraints.close ();

    // stuff for the second dof handler
    dof_handler_2.distribute_dofs (fe_2);

    m_Psi_C.reinit (dof_handler_2.n_dofs());

    constraints_2.clear ();
    DoFTools::make_hanging_node_constraints (dof_handler_2, constraints_2);
    VectorTools::interpolate_boundary_values (dof_handler_2, 0, ZeroFunction<dim>(2), constraints_2);
    constraints_2.close ();
  }

  template <int dim>
  void binR_to_binC<dim>::run()
  {
    make_grid();
    setup_system();

    ifstream in(m_bin_filename);
    m_Psi_R.block_read(in);

    MyComplexTools::Interpolate_R_to_C( dof_handler, fe, m_Psi_R, dof_handler_2, fe_2, constraints_2, m_Psi_C );

    ofstream out("C" + m_bin_filename);
    m_Psi_C.block_write(out);
/*
    DataOut<1> data_out;
    data_out.attach_dof_handler (dof_handler_2);
    data_out.add_data_vector (m_Psi_C, "Psi");
    data_out.build_patches ();
    
    ofstream out2 ("C" + m_bin_filename + ".gnuplot");
    data_out.write_gnuplot ( out2 );
*/
  }
} // end of namespace

int main ( int argc, char *argv[] )
{
  cxxopts::Options options("binR_to_binC_1", "Converts real 1D deal.ii binary format to complex deal.ii binary format.");
  
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
  int dim = 0;

  try
  {
    dim = int(params.Get_Mesh("DIM", 0));
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
    switch (dim)
    {
      case 1:
      {
        HelperPrograms::binR_to_binC<1> solver( params, bin_filename );
        solver.run();
        break;
      }
      default:
        std::cout << "You have found a new dimension!" << std::endl;
    }
  }
  return EXIT_SUCCESS;
}
