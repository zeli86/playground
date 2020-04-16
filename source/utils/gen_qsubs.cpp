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

/*
Želimir Marojević
*/

#include <deal.II/base/parameter_handler.h>

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <cstdlib>
#include <fstream>

#include "ParameterReader.h"

using namespace std;

int main( int argc, char *argv[] )
{
  const string HomePath = getenv ("HOME");
  const string folder_name = "/GOST_2015-02-18_15:10";
  
  string shell_command, filename1, filename2;
  
  vector<string> sub_folder_names;
  sub_folder_names.push_back( "/0_0_0" );
  sub_folder_names.push_back( "/0_0_1" );
  sub_folder_names.push_back( "/0_1_0" );
  sub_folder_names.push_back( "/0_1_1" );
  sub_folder_names.push_back( "/1_0_0" );
  sub_folder_names.push_back( "/1_0_1" );
  sub_folder_names.push_back( "/2_0_0" );
  sub_folder_names.push_back( "/2_0_1" );

  vector<string> sub_sub_folder_names;
  //sub_sub_folder_names.push_back( "/0000" );
  //sub_sub_folder_names.push_back( "/0010" );
  sub_sub_folder_names.push_back( "/0020" );
  sub_sub_folder_names.push_back( "/0030" );
  //sub_sub_folder_names.push_back( "/0040" );
  sub_sub_folder_names.push_back( "/0050" );
  //sub_sub_folder_names.push_back( "/0060" );
  //sub_sub_folder_names.push_back( "/0070" );

  ParameterHandler prm;
  ParameterReader param(prm);
  param.declare_parameters();
  
  ofstream bash_script( HomePath + "/run.sh");
  bash_script << "#!/bin/bash" << endl;
  
  for( auto sfn : sub_folder_names )
  {
    for( auto ssfn : sub_sub_folder_names )
    {
      filename1 = HomePath + folder_name + sfn + "/params.prm";
      filename2 = HomePath + folder_name + sfn + ssfn + "/params.prm";
      param.read_parameters( filename1 );
      prm.enter_subsection("Algorithm control parameters");
      prm.set( "NA", long(2000) );
      prm.set( "NK", long(20) );
      prm.set( "dt", 0.01 );
      prm.leave_subsection();
      param.print_parameters( filename2 );
      
      filename1 = HomePath + folder_name + sfn + ssfn + "/job.pbs";
      ofstream pbs_file( filename1 );
      pbs_file << "#!/bin/bash\n";
      pbs_file << "#PBS -M zeli@zarm.uni-bremen.de\n";
      pbs_file << "#PBS -N " << sfn + ssfn << "\n";
      pbs_file << "#PBS -j oe\n";
      pbs_file << "#PBS -l walltime=600:00:00\n"; // our cluster has a max of 600:00:00
      pbs_file << "#PBS -l nodes=1:ppn=4\n";
      pbs_file << "#PBS -V\n";
      pbs_file << "cd " << HomePath + folder_name + sfn + ssfn << "\n";
      pbs_file << "mpirun rt_prop_mpi_cs\n";
      pbs_file.close();
      
      bash_script << "qsub -k oe " << filename1 << endl;
    }
  }
  bash_script.close();
  shell_command = "chmod +x " + HomePath + "/run.sh";
  system( shell_command.c_str() );  
return 0;
}


