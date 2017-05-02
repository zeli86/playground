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

#include <deal.II/base/parameter_handler.h>

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <cstdlib>
#include <fstream>

using namespace std;

int main( int argc, char *argv[] )
{
  const string HomePath = getenv ("HOME");
  
  string shell_command, filename1, filename2;
  
  vector<string> sub_folder_names;
  sub_folder_names.push_back( "/0_1_0" );
  sub_folder_names.push_back( "/0_1_1" );
  sub_folder_names.push_back( "/1_0_0" );
  sub_folder_names.push_back( "/1_0_1" );
  sub_folder_names.push_back( "/2_0_0" );
  sub_folder_names.push_back( "/2_0_1" );

  ofstream bash_script( HomePath + "/run.sh");
  bash_script << "#!/bin/bash" << endl;
  
  int counter = 0;

  for( auto sfn : sub_folder_names )
  {
    for( int i=0; i<90; i+=5 )
    {
      counter++;

      char buf[10] = {};
      sprintf( buf, "/%.4d", i );
      string ssfn = buf;

      string jobname = "matrixhell-" + to_string(counter);

      filename1 = jobname + "-job.sh";
      ofstream slurm_file( filename1 );
      slurm_file << "#!/bin/bash\n";
      slurm_file << "#SBATCH -J " << jobname << "\n";
      slurm_file << "#SBATCH -N 1 -n 8\n";
      slurm_file << "#SBATCH --mail-user=zeli@zarm.uni-bremen.de\n";
      slurm_file << "#SBATCH --mail-type=BEGIN --mail-type=END\n";
      slurm_file << "#SBATCH -o matrixhell-%j.out\n";

      slurm_file << "export MODULEPATH=$MODULEPATH:/home/zeli/local/modules/modulefiles\n";
      slurm_file << "module load GCC/5.4.0\n";
      slurm_file << "module load Cmake/3.7.0\n";
      slurm_file << "module load openmpi-2.0.2\n";
      slurm_file << "module load gsl-2.3\n";
      slurm_file << "module load p4est-1.1\n";
      slurm_file << "module load petsc-3.6.3\n";
      slurm_file << "module load deal.ii-8.4.1\n";
      slurm_file << "cd " << HomePath + sfn + ssfn << "\n";
      slurm_file << "mpirun -np 8 compute_lin_op final.bin\n";
      slurm_file.close();
      
      bash_script << "sbatch " << filename1 << endl;
    }
  }
  bash_script.close();
  shell_command = "chmod +x " + HomePath + "/run.sh";
  //system( shell_command.c_str() );  
return 0;
}


