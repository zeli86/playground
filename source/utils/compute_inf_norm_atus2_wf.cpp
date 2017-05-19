//
// ATUS2 - The ATUS2 package is atom interferometer Toolbox developed at ZARM
// (CENTER OF APPLIED SPACE TECHNOLOGY AND MICROGRAVITY), Germany. This project is
// founded by the DLR Agentur (Deutsche Luft und Raumfahrt Agentur). Grant numbers:
// 50WM0942, 50WM1042, 50WM1342.
// Copyright (C) 2017 Želimir Marojević, Ertan Göklü, Claus Lämmerzahl
//
// This file is part of ATUS2.
//
// ATUS2 is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// ATUS2 is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with ATUS2.  If not, see <http://www.gnu.org/licenses/>.
//

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <omp.h>
#include <fstream>
#include <iostream>
#include <memory>
#include "atus2.h"
#include "pugixml.hpp"

int main(int argc, char *argv[])
{
  std::string filename = argv[1];
  generic_header header = {};

  std::ifstream in;
  in.open(filename.c_str());
  if ( !in.is_open() ) return EXIT_FAILURE;
  in.read( reinterpret_cast<char*>(&header), sizeof(generic_header) );

  char *envstr = getenv( "MY_NO_OF_THREADS" );
  if ( envstr != nullptr )
  {
    int no_of_threads = atoi( envstr );
    omp_set_num_threads( no_of_threads );
  }   

  if ( header.bComplex == true )
  {
    const long long N = 2 * header.nDimX * header.nDimY * header.nDimZ;
    double * field = new double[N];
    in.read( reinterpret_cast<char*>(field), N*sizeof(double) );

    double retval=0;
    #pragma omp parallel for reduction(max:retval)
    for( long long i=0; i<N; i++ )
    {
      if( fabs(field[i]) > retval ) retval = fabs(field[i]); 
    }

    pugi::xml_document doc;
    pugi::xml_parse_result result = doc.load_file("info.xml");
    std::cerr << "Load result: " << result.description() << std::endl; 
        
    double mu = std::stod( doc.child( "INFO" ).child( "MU" ).child_value() );
    
    delete [] field;

    std::cout << mu << "\t" << retval << "\n";
  }
}
