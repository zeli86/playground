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

#include <cstdio>
#include <cstdlib>
#include <string.h>
#include <fstream>

using namespace std;

int main( int argc, char *argv[] )
{
  ofstream pbs_file( "compute_lin_op.pbs" );
  pbs_file << "#!/bin/bash\n";
  pbs_file << "#PBS -A hbp00035\n";
  pbs_file << "#PBS -j oe\n";
  pbs_file << "#PBS -l walltime=11:59:00\n"; 
  pbs_file << "#PBS -e my_job.$PBS_JOBID.err\n";
  pbs_file << "#PBS -o my_job.$PBS_JOBID.out\n";
  pbs_file << "#PBS -l feature=mpp\n";
  pbs_file << "#PBS -l nodes=15:ppn=24\n";
      
  for( int i=0; i<90; ++i )
  {
    char buf[10] = {};
    sprintf( buf, "%.4d", i );
    string subfolder = buf;    
    pbs_file << "cd " << subfolder << "\n";
    pbs_file << "aprun -n 4 compute_lin_op final.bin\n";
    pbs_file << "cd ..\n";
  }

  pbs_file << "wait\n" << endl;
return 0;
}


