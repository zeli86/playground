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
#include <iostream>
#include <fstream>

#include "pugixml.hpp"

using namespace std;

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

      pugi::xml_document doc;
      pugi::xml_node parameter_node = doc.append_child("PARAMETER");
      pugi::xml_node physics_node = parameter_node.append_child("PHYSICS");
      pugi::xml_node mesh_node = parameter_node.append_child("MESH");
      pugi::xml_node algorithm_node = parameter_node.append_child("ALGORITHM");      

      // PARAMTER childs
      pugi::xml_node node = parameter_node.append_child("FILENAME");
      node.append_child(pugi::node_pcdata).set_value("final.bin");

      node = parameter_node.append_child("guess_fct");
      node.append_child(pugi::node_pcdata).set_value("exp(-x^2-y^2)");

      // PHYSICS childs
      node = physics_node.append_child("QN1");
      tmpstr = to_string(i1) + "," + to_string(j1);
      node.append_child(pugi::node_pcdata).set_value(tmpstr.c_str());

      node = physics_node.append_child("omega");
      tmpstr = to_string(omega[0]) + "," + to_string(omega[1]) + "," + to_string(omega[2]);
      node.append_child(pugi::node_pcdata).set_value(tmpstr.c_str());

      node = physics_node.append_child("gs_1");
      node.append_child(pugi::node_pcdata).set_value("1,1");

      node = physics_node.append_child("mu");
      node.append_child(pugi::node_pcdata).set_value("5,5");

      // MESH childs
      node = mesh_node.append_child("xrange");
      node.append_child(pugi::node_pcdata).set_value("0,30");

      node = mesh_node.append_child("yrange");
      node.append_child(pugi::node_pcdata).set_value("-15,15");

      node = mesh_node.append_child("yrange");
      node.append_child(pugi::node_pcdata).set_value("-15,15");

      node = mesh_node.append_child("global_refinements");
      node.append_child(pugi::node_pcdata).set_value("9");

      // ALGORTHM childs
      node = algorithm_node.append_child("ti");
      node.append_child(pugi::node_pcdata).set_value("1");

      node = algorithm_node.append_child("NA");
      node.append_child(pugi::node_pcdata).set_value("10");

      node = algorithm_node.append_child("NK");
      node.append_child(pugi::node_pcdata).set_value("10");

      node = algorithm_node.append_child("dt");
      node.append_child(pugi::node_pcdata).set_value("0.001");

      node = algorithm_node.append_child("epsilon");
      node.append_child(pugi::node_pcdata).set_value("1e-5,1e-10");

      node = algorithm_node.append_child("dmu");
      tmpstr = to_string(dmu);
      node.append_child(pugi::node_pcdata).set_value(tmpstr.c_str());

      node = algorithm_node.append_child("Ndmu");
      tmpstr = to_string(Ndmu);
      node.append_child(pugi::node_pcdata).set_value(tmpstr.c_str());

      // add param node before the description
      //pugi::xml_node param = node.insert_child_before("param", descr);

      // add attributes to param node
      //param.append_attribute("name") = "version";
      //param.append_attribute("value") = 1.1;
      //param.insert_attribute_after("type", param.attribute("name")) = "float";

      cout << "Saving result: " << doc.save_file(filename) << endl;
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


