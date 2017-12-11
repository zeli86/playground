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

#include "tinyxml2.h"

using namespace std;
using namespace tinyxml2;

const long N1[] = {4,4,4};
const long NA = 10000;
const long Ndmu = 40;
const double dmu = .1;
double omega[] = {0.5, 0.5, 0.5};

int main( int argc, char *argv[] )
{
  //const char * HomePath = getenv ("HOME");

  string tmpstr;

  char base_folder[255];
  char folder[255];
  char filename[255];
  char shellcmd[255];
  char datetime[255];

  time_t rawtime;
  struct tm * timeinfo;

  time(&rawtime);
  timeinfo = localtime(&rawtime);
  strftime(datetime,255,"%F_%R",timeinfo);

#if POTENTIAL==1
  sprintf( base_folder, "HTRAP_%s", datetime );
#endif
#if POTENTIAL==2
  sprintf( base_folder, "GOST_%s", datetime );
  omega[0] = 0.5039287608;
#endif
  
#if DIMENSION==2  
  for( long i1=0; i1<N1[0]; i1++ )
  {
    for( long j1=0; j1<N1[1]; j1++ )
    {
      sprintf( folder, "%ld_%ld", i1, j1 );
      sprintf( shellcmd, "mkdir -p %s/%s", base_folder, folder );
      system( shellcmd );

      sprintf( filename, "%s/%s/params.xml", base_folder, folder );
      FILE * fh = fopen( filename, "w" );

      XMLPrinter printer(fh);
      printer.OpenElement( "PARAMETER" );
        printer.OpenElement( "filename" );
          printer.PushText("final.bin");
        printer.CloseElement();

      
        // problem related parameter
        printer.OpenElement( "PHYSICS" );
          printer.OpenElement( "gs_1" );
            printer.PushText("1,1");
          printer.CloseElement();
        printer.CloseElement(); // close PHYSICS
        
        // mesh related parameter
        printer.OpenElement( "MESH" );
          printer.OpenElement( "xrange" );
            printer.PushText("0,10");
          printer.CloseElement();
          printer.OpenElement( "global_refinements" );
            printer.PushText( "11" );
          printer.CloseElement();          
        printer.CloseElement(); // close MESH

        // algorithm related parameter
        printer.OpenElement( "ALGORITHM" );
          printer.OpenElement( "ti" );
            printer.PushText( "1" );
          printer.CloseElement();
          printer.OpenElement( "NA" );
            tmpstr = to_string(NA);
            printer.PushText( tmpstr.c_str() );
          printer.CloseElement();
          printer.OpenElement( "NK" );
            printer.PushText( "10" );
          printer.CloseElement();
          printer.OpenElement( "dmu" );
            tmpstr = to_string(dmu);
            printer.PushText( tmpstr.c_str() );
          printer.CloseElement();
          printer.OpenElement( "Ndmu" );
            tmpstr = to_string(Ndmu);
            printer.PushText( tmpstr.c_str() );
          printer.CloseElement();
          printer.OpenElement( "dt" );
            printer.PushText( "0.001" );
          printer.CloseElement();
          printer.OpenElement( "epsilon" );
            printer.PushText( "1e-5,1e-10" );
          printer.CloseElement();
        printer.CloseElement(); // close ALGORITHM
      printer.CloseElement(); // close PARAMETER
      
/*
      sprintf( filename, "%s/%s.pbs", base_folder, folder );
      ofstream pbs_file( filename );
      pbs_file << "#!/bin/bash\n";
      pbs_file << "#PBS -M zeli@zarm.uni-bremen.de\n";
      pbs_file << "#PBS -N " << folder << "\n";
      pbs_file << "#PBS -j oe\n";
      pbs_file << "#PBS -l walltime=10:00:00\n"; // our cluster has a max of 600:00:00
      pbs_file << "#PBS -l nodes=1:ppn=8\n";
      pbs_file << "#PBS -V\n";
      pbs_file << "cd " << HomePath << "/" << base_folder << "/" << folder << "\n";
      pbs_file << "mpirun my_solver_mpi\n";
      //pbs_file << "mpirun --map-by core my_csolver_mpi\n";
      pbs_file.close();
*/
    }
  }
#endif
return 0;
}


