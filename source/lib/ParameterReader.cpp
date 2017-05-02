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

#include "ParameterReader.h"
#include <deal.II/numerics/data_out.h>
#include <fstream>

using namespace dealii;

ParameterReader::ParameterReader(ParameterHandler &ph) : prm(ph)
{}


void ParameterReader::declare_parameters()
{  
  prm.enter_subsection ("Mesh & geometry parameters");
  {
      prm.declare_entry("global_refinements", "7", Patterns::Integer(), "level of global mesh refinements for the inital mesh");
      prm.declare_entry("max_refinements", "3", Patterns::Integer(), "maximal number of mesh refinements");
      prm.declare_entry("max_grid_level", "11", Patterns::Integer(), "");
      prm.declare_entry("min_grid_level", "7", Patterns::Integer(), "");
      prm.declare_entry("delta_refinement", "5", Patterns::Integer(), "frequency of mesh refinement");
      prm.declare_entry("top_fraction", "0.33", Patterns::Double(), "");
      prm.declare_entry("bottom_fraction", "0.11", Patterns::Double(), "");
#if POTENTIAL==1
      prm.declare_entry("xMin", "-15.0", Patterns::Double(), "");
      prm.declare_entry("xMax", "15.0", Patterns::Double(), "");
#endif
#if POTENTIAL==2
      prm.declare_entry("xMin", "0.0", Patterns::Double(), "");
      prm.declare_entry("xMax", "30.0", Patterns::Double(), "");
#endif
      prm.declare_entry("yMin", "-15.0", Patterns::Double(), "");
      prm.declare_entry("yMax", "15.0", Patterns::Double(), "");
      prm.declare_entry("zMin", "-15.0", Patterns::Double(), "");
      prm.declare_entry("zMax", "15.0", Patterns::Double(), "");
      prm.declare_entry("filename", "final.bin", Patterns::Anything(), "");
  }
  prm.leave_subsection ();

  prm.enter_subsection ("Physical constants");
  {
    prm.declare_entry("QN1_x", "1", Patterns::Integer(), "quantum number 1 x-direction");
    prm.declare_entry("QN1_y", "1", Patterns::Integer(), "quantum number 1 y-direction");
    prm.declare_entry("QN1_z", "1", Patterns::Integer(), "quantum number 1 z-direction");
    prm.declare_entry("QN2_x", "1", Patterns::Integer(), "quantum number 2 x-direction");
    prm.declare_entry("QN2_y", "1", Patterns::Integer(), "quantum number 2 y-direction");
    prm.declare_entry("QN2_z", "1", Patterns::Integer(), "quantum number 2 z-direction");
    prm.declare_entry("QN3_x", "1", Patterns::Integer(), "quantum number 3 x-direction");
    prm.declare_entry("QN3_y", "1", Patterns::Integer(), "quantum number 3 y-direction");
    prm.declare_entry("QN3_z", "1", Patterns::Integer(), "quantum number 3 z-direction");
    prm.declare_entry("QN4_x", "1", Patterns::Integer(), "quantum number 4 x-direction");
    prm.declare_entry("QN4_y", "1", Patterns::Integer(), "quantum number 4 y-direction");
    prm.declare_entry("QN4_z", "1", Patterns::Integer(), "quantum number 4 z-direction");
    prm.declare_entry("omega_x", "0.5", Patterns::Double(), "frequency of the harmonic trap for the x-direction");
    prm.declare_entry("omega_y", "0.5", Patterns::Double(), "frequency of the harmonic trap for the y-direction");
    prm.declare_entry("omega_z", "0.5", Patterns::Double(), "frequency of the harmonic trap for the z-direction");
    prm.declare_entry("s_x", "1", Patterns::Double(), "");
    prm.declare_entry("s_y", "1", Patterns::Double(), "");
    prm.declare_entry("s_z", "1", Patterns::Double(), "");
    prm.declare_entry("gf", "1.0", Patterns::Double(), "gravitational acceleration");
    prm.declare_entry("gs_11", "1.0", Patterns::Double(), "self interaction parameter");
    prm.declare_entry("gs_12", "1.0", Patterns::Double(), "self interaction parameter");
    prm.declare_entry("gs_22", "1.0", Patterns::Double(), "self interaction parameter");
    prm.declare_entry("gs_21", "1.0", Patterns::Double(), "self interaction parameter");
    prm.declare_entry("mu", "5.0", Patterns::Double(), "mu" );
    prm.declare_entry("mu2", "5.0", Patterns::Double(), "mu2" );
    prm.declare_entry("guess_fct", "exp(-0.5*x^2)*exp(-0.5*y^2)", Patterns::Anything(), "" );
    prm.declare_entry("sign", "-1.0", Patterns::Double(), "");
  }
  prm.leave_subsection ();

  prm.enter_subsection ("Algorithm control parameters");
  {
    prm.declare_entry("epsilon", "1e-5", Patterns::Double(), "convergency criterium" );
    prm.declare_entry("ti", "1.0", Patterns::Double(), "initial value for the point in function space" );
    prm.declare_entry("NA", "10", Patterns::Integer(), "frequency of data output" );
    prm.declare_entry("NK", "10", Patterns::Integer(), "number of intermediate steps" );
    prm.declare_entry("dmu", "0.5", Patterns::Double(), "delta mu" );
    prm.declare_entry("Ndmu", "10", Patterns::Integer(), "number of dmu steps" );
    prm.declare_entry("dt", "1e-3", Patterns::Double(), "dt for real time propagation" );
  }
  prm.leave_subsection ();
}

bool ParameterReader::read_parameters ( const std::string filename )
{
  declare_parameters();  
  return prm.read_input(filename);
}

void ParameterReader::print_parameters ( const std::string filename )
{
  std::ofstream out(filename);
  prm.print_parameters(out, ParameterHandler::Text);
  std::cout << "\nParameter file " << filename << " created.\n" << std::endl;
}



