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
Želimir Marojević
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <cstdlib>
#include <fstream>
#include <vector>

using namespace std;

int main( int argc, char *argv[] )
{
  vector<string> sub_folder_names;
  //sub_folder_names.push_back( "0_1_0" );
  sub_folder_names.push_back( "0_1_1" );
  sub_folder_names.push_back( "1_0_0" );
  sub_folder_names.push_back( "1_0_1" );
  sub_folder_names.push_back( "2_0_0" );
  sub_folder_names.push_back( "2_0_1" );

  ofstream bash_script( "run.sh");
  bash_script << "#!/bin/bash\n";
  
  for( auto sfn : sub_folder_names )
  {
    bash_script << "cd " << sfn << "\n";

    for( int i=0; i<91; i+=5 )
    {
      char buf[10] = {};
      sprintf( buf, "%.4d", i );
      string ssfn = buf;

      bash_script << "cd " << ssfn << "\n";
      bash_script << "pwd\n";
      bash_script << "mpirun -np 6 compute_lin_op\n";
      bash_script << "cd ..\n";      
    }

    bash_script << "cd ..\n";      
  }
return 0;
}


