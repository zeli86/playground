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


#include <cstdio>
#include <cstdlib>
#include <string>
#include <fstream>
#include <vector>
#include <iostream>

using namespace std;

int main ( int argc, char *argv[] )
{
  if( argc != 4 ) 
  {
    cout << argv[0] <<  " ti dt tf" << endl;
    return EXIT_FAILURE;
  }
  
  vector<string> args;
  for( int i=1; i<argc; ++i )
    args.push_back(argv[i]);
  
  const double ti = stod(args[0]);
  const double dt = stod(args[1]);
  const double tf = stod(args[2]);
  
  ofstream out("solution.pvd");
  out << "<?xml version=\"1.0\"?>\n<VTKFile type=\"Collection\" version=\"0.1\" ByteOrder=\"LittleEndian\">\n<Collection>\n";
  for( double t=ti; t<=tf; t+=dt )
    out << "<DataSet timestep=\"" << to_string(t) << "\" group=\"\" part=\"0\" file=\"step-" << to_string(t) << ".vtu\"/>\n";
  
  out << "</Collection>\n</VTKFile>\n";
return EXIT_SUCCESS;
}