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

/*
    Želimir Marojević, Ertan Göklü, Claus Lämmerzahl
    (ZARM - Center of Applied Space Technology and Microgravity, Germany, http://www.zarm.uni-bremen.de/)

    Public use and modification of this code are allowed provided that the following papers are cited:

    Marojević, Želimir, Ertan Göklü, und Claus Lämmerzahl. "ATUS-PRO: A FEM-based solver for the time-dependent and stationary Gross–Pitaevskii equation",
    Computer Physics Communications, Vol 202, 2016, p. 216--232. doi:10.1016/j.cpc.2015.12.004.

    W. Bangerth and D. Davydov and T. Heister and L. Heltai and G. Kanschat and M. Kronbichler and M. Maier and B. Turcksin and D. Wells.
    "The \texttt{deal.II} Library, Version 8.4", Journal of Numerical Mathematics, vol 24, 2016.

    The authors would be grateful for all information and/or comments regarding the use of the code.
*/

#include <deal.II/lac/vector.h>

#include <deal.II/base/function.h>
#include <deal.II/base/function_parser.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>

#include <cstdlib>
#include <fstream>
#include <iostream>

#include "cxxopts.hpp"
#include "MyComplexTools.h"
#include "MyParameterHandler.h"
#include "functions.h"
#include "global.h"

namespace HelperPrograms
{
  using namespace std;
  using namespace dealii;

  template <int dim> 
  class WFGenerator
  {
  public:
    WFGenerator( const std::string& );
    virtual ~WFGenerator();

    void run( const string& );

  protected:
    double m_xmin, m_xmax;
    double m_ymin, m_ymax;

    unsigned m_global_refinement;

    void make_grid();
    void setup_system();
    void save(string);

    MyParameterHandler m_ph;
    Triangulation<dim> triangulation;
    FESystem<dim> fe;
    DoFHandler<dim> dof_handler;

    Vector<double> m_Psi;
    std::vector<std::string> m_wf_str;
  };

  /**
   * Constructor
   */
  template <int dim>
  WFGenerator<dim>::WFGenerator( const string& xmlfilename ) : 
    m_ph(xmlfilename),
    triangulation(),
    fe (FE_Q<dim>(gl_degree_fe), 2), 
    dof_handler(triangulation)
  {
    try
    {
      m_xmin = m_ph.Get_Mesh("xrange", 0);
      m_xmax = m_ph.Get_Mesh("xrange", 1);
      //m_ymin = m_ph.Get_Mesh("yrange", 0);
      //m_ymax = m_ph.Get_Mesh("yrange", 1);
      m_wf_str = m_ph.Get_AllStrings("WAVEFUNCTION_" + to_string(dim) + "D");
      m_global_refinement = unsigned(m_ph.Get_Mesh("global_refinements", 0));
    }
    catch (const std::string info)
    {
      std::cerr << info << endl;
    }
  }

  template <int dim>
  WFGenerator<dim>::~WFGenerator()
  {
    dof_handler.clear();
  }

  template <int dim> 
  void WFGenerator<dim>::make_grid()
  {
    Point<dim,double> pt1, pt2;

    double min[] = {m_xmin, m_ymin};
    double max[] = {m_xmax, m_ymax};

    for ( int i = 0; i < dim; i++ )
    {
      pt1(i) = min[i];
      pt2(i) = max[i];
    }

    GridGenerator::hyper_rectangle(triangulation, pt2, pt1);
    triangulation.refine_global(m_global_refinement);
  }

  template <int dim> 
  void WFGenerator<dim>::setup_system()
  {
    dof_handler.distribute_dofs(fe);

    m_Psi.reinit(dof_handler.n_dofs());

    cout << "no of dofs = " << dof_handler.n_dofs() << endl;
  }

  template <int dim> 
  void WFGenerator<dim>::run( const string& filename )
  {
    make_grid();
    setup_system();

    FunctionParser<dim> Psi(2);
    Psi.initialize( FunctionParser<dim>::default_variable_names(), m_wf_str, m_ph.Get_Constants_Map());

    VectorTools::interpolate(dof_handler, Psi, m_Psi);

    double N = MyComplexTools::Particle_Number( dof_handler, fe, m_Psi );
    cout << "N == " << N << endl;

    ofstream ofs(filename);
    m_Psi.block_write(ofs);
  }
} // end of namespace

int main(int argc, char *argv[])
{
  using namespace dealii;
  deallog.depth_console(0);

  cxxopts::Options options("binR_to_binC", "Converts real deal.ii binary format to complex deal.ii binary format.");
  
  options.add_options()
  ("p,params",  "input parameter xml file" , cxxopts::value<std::string>()->default_value("params.xml") )
  ("o,output",  "output 1D binary file" , cxxopts::value<std::string>()->default_value("Psi0.1d.bin") )
  ;
  
  options.parse_positional({"positional"});
  auto result = options.parse(argc, argv);

  std::string bin_filename, params_filename;
  try
  {
    bin_filename = result["o"].as<std::string>(); 
    params_filename = result["p"].as<std::string>();
  }
  catch (const cxxopts::OptionException& e)
  {
    std::cout << "error parsing options: " << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  HelperPrograms::WFGenerator<1> solver(params_filename);
  solver.run(bin_filename);
  return EXIT_SUCCESS;
}