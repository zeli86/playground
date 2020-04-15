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

/**
 * @file rtprop.cpp Real time propagation for the Gross--Pitaevskii equation (Cartesian coordinates)
 * @brief Real time propagation for the Gross--Pitaevskii equation (Cartesian coordinates)
 * @author Želimir Marojević
 */

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/timer.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparsity_tools.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_tools.h>    

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_q.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/function_parser.h>

#include <fstream>
#include <iostream>
#include <stdlib.h>

#include "global.h"
#include "MyParameterHandler.h"
#include "MyComplexTools.h"
#include "my_table.h"
#include "cxxopts.hpp"

///namespace realtime_propagation
namespace realtime_propagation
{
  using namespace std;
  using namespace dealii;

  class MySolver
  {
  public:
    MySolver( const std::string& );
    ~MySolver();

    void run ( const std::string& );

  protected:
    void make_grid();
    void setup_system();
    
    void solve();
    void output_results ( string );
    void load( string );
    void save( string );    

    MyParameterHandler m_ph;
    Triangulation<1> triangulation;
    FESystem<1> fe;
    DoFHandler<1> dof_handler;

    SparsityPattern sparsity_pattern;
    SparseMatrix<double> system_matrix;

    Vector<double> system_rhs;
    Vector<double> m_Psi; // Psi(t)
    Vector<double> m_workspace;
    Vector<double> m_error_per_cell;

    double m_gs;
    double m_t;
    double m_dt;
    double m_xmin;
    double m_xmax;

    /// total number of outputs
    unsigned m_NA;
    /// number of intermediate time steps between two outputs
    unsigned m_NK;
    
    unsigned m_global_refinement;
    
    MyTable m_table;
    string m_filename;
  };

  /** Default constructor
   */
  MySolver::MySolver ( const std::string& xmlfilename ) 
    : 
    m_ph(xmlfilename),
    triangulation (),
    fe (FE_Q<1>(gl_degree_fe), 2),
    dof_handler (triangulation)
  {
    try
    {
      m_gs = m_ph.Get_Physics("gs_1",0);
      
      m_xmin = m_ph.Get_Mesh("xrange",0);
      m_xmax = m_ph.Get_Mesh("xrange",1);
      m_global_refinement = unsigned(m_ph.Get_Mesh("global_refinements",0));

      m_NA = m_ph.Get_Algorithm("NA",0);
      m_NK = m_ph.Get_Algorithm("NK",0);
      m_dt = m_ph.Get_Algorithm("dt",0);
    }
    catch( const std::string info )
    {
      std::cerr << info << endl;
      exit(0);
    }    
    m_t = 0;
  }

  /** Default destructor
   */
  MySolver::~MySolver ()
  {
    dof_handler.clear ();
  }
  
  #include "utils_complex.h"
 
  void MySolver::make_grid ()
  {
    Point<1,double> pt1( m_xmin );
    Point<1,double> pt2( m_xmax );

    GridGenerator::hyper_rectangle(triangulation, pt2, pt1);
    triangulation.refine_global(m_global_refinement);
  }
  
  void MySolver::setup_system()
  {
    dof_handler.distribute_dofs (fe);
    //DoFRenumbering::component_wise (dof_handler);

    m_Psi.reinit (dof_handler.n_dofs());
    m_workspace.reinit (dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());
    m_error_per_cell.reinit(dof_handler.n_dofs());

    DynamicSparsityPattern dsp(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern (dof_handler, dsp);
    sparsity_pattern.copy_from(dsp);
    system_matrix.reinit (sparsity_pattern);
  }  

  /** Output of data to vtu files
   * @param[in] path specifies the output folder 
   */
  void MySolver::output_results ( string path ) 
  {
    ComputeIntensity<1> intensities;
    ComputePhase<1> phase;
    
    vector<std::string> solution_names;

    DataOut<1> data_out;
    data_out.attach_dof_handler (dof_handler);
    solution_names.push_back ("Re Psi");
    solution_names.push_back ("Im Psi");    
    data_out.add_data_vector (m_Psi, solution_names);
    data_out.add_data_vector (m_Psi , intensities);
    data_out.build_patches ();
    
    string filename = path + "solution-" + to_string(m_t) + ".gnuplot";
    ofstream output (filename);
    data_out.write_gnuplot ( output );
  }   

  void MySolver::save( string filename )
  {
  }

  void MySolver::load( string filename )
  {
    ifstream in(filename);
    m_Psi.block_read(in);
  }    
  
  void MySolver::solve ()
  {
    map<types::global_dof_index, double> boundary_values;
    VectorTools::interpolate_boundary_values (dof_handler, 0, ZeroFunction<1>(), boundary_values);
    MatrixTools::apply_boundary_values (boundary_values, system_matrix, m_Psi, system_rhs);    

    SparseDirectUMFPACK  A_direct;
    A_direct.initialize(system_matrix);
    A_direct.vmult (m_Psi, system_rhs);
  }

  void MySolver::run( const std::string& filename )
  {
    make_grid();
    setup_system();

    load( filename );

    double min_cell_diameter = GridTools::minimal_cell_diameter(triangulation);
    double max_cell_diameter = GridTools::maximal_cell_diameter(triangulation);
    cout << "min_cell_diameter = " << min_cell_diameter << "\n";
    cout << "max_cell_diameter = " << max_cell_diameter << "\n";
    cout << "dt/dx^2 == " << m_dt/(min_cell_diameter*min_cell_diameter) << endl;
    
    vector<double> p(1,0);
    vector<double> pos(1,0);

    output_results("");

    double N = MyComplexTools::Particle_Number( dof_handler, fe, m_Psi );
    //cout << "N == " << N << endl;
    MyComplexTools::Expectation_value_position( dof_handler, fe, m_Psi, pos );
    MyComplexTools::Expectation_value_momentum( dof_handler, fe, m_Psi, p );

    cout << "t == " << m_t << endl;
    cout << "N == " << N << endl;
    cout << "p == " << p[0]/N << endl;
    cout << "pos == " << pos[0]/N << endl;

    double dth = 0.5*m_dt;

    FunctionParser<1> Potential;
    Potential.initialize( FunctionParser<1>::default_variable_names(), m_ph.Get_String( "POTENTIAL", 0 ), m_ph.Get_Constants_Map());

    for( unsigned i=1; i<=m_NA; ++i )
    {
      MyComplexTools::AssembleSystem_LIN_Step( dof_handler, fe, m_Psi, dth, system_matrix, system_rhs );
      solve();
      for( unsigned j=1; j<m_NK; ++j )
      {
        MyComplexTools::AssembleSystem_NL_Step( dof_handler, fe, m_Psi, Potential, m_dt, m_gs,  system_matrix, system_rhs );
        solve();
        MyComplexTools::AssembleSystem_LIN_Step( dof_handler, fe, m_Psi, m_dt, system_matrix, system_rhs );
        solve();
      }
      MyComplexTools::AssembleSystem_NL_Step( dof_handler, fe, m_Psi, Potential, m_dt, m_gs,  system_matrix, system_rhs );
      solve();
      MyComplexTools::AssembleSystem_LIN_Step( dof_handler, fe, m_Psi, dth, system_matrix, system_rhs );
      solve();
      
      N = MyComplexTools::Particle_Number( dof_handler, fe, m_Psi );
      MyComplexTools::Expectation_value_position( dof_handler, fe, m_Psi, pos );
      MyComplexTools::Expectation_value_momentum( dof_handler, fe, m_Psi, p );

      cout << "N == " << N << endl;
      cout << "p == " << p[0]/N << endl;
      cout << "pos == " << pos[0]/N << endl;

      output_results("");
      save( "solution-" + to_string(m_t) + ".bin" );
    }
  }
} // end of namespace 

int main ( int argc, char *argv[] )
{
  using namespace dealii;
  deallog.depth_console (0);

  cxxopts::Options options("rt_prop_1", "real time propagation 1D");
  
  options.add_options()
  ("i,input", "input initial wave function" , cxxopts::value<std::string>()->default_value("Cfinal.bin") )
  ("p,params", "input parameter xml file" , cxxopts::value<std::string>()->default_value("params.xml") )
  ("help","Print help")
  ;
  
  auto result = options.parse(argc, argv);

  if (result.count("") == 0)
  {
    std::cout << options.help({""}) << std::endl;
    return EXIT_FAILURE;
  }

  std::string bin_filename, params_filename;
  try
  {
    bin_filename = result["i"].as<std::string>(); 
    params_filename = result["p"].as<std::string>();
  }
  catch (const cxxopts::OptionException& e)
  {
    std::cout << "error parsing options: " << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  realtime_propagation::MySolver solver(params_filename);
  solver.run(bin_filename);
return EXIT_SUCCESS;
}
