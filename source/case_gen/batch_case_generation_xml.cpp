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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <map>
#include <functional>

#include "pugixml.hpp"
#include "cxxopts.hpp"

using namespace std;

const long N1[] = {4,4,4};
double omega[] = {0.5, 0.5, 0.5};

std::string qn_to_str_1D( const int l, const int /*maxN*/ )
{
  return std::to_string(l) + ",0,0"; 
}

std::string qn_to_str_2D( const int l, const int maxN )
{
  int I = l / maxN; 
  int J = l - I*maxN;
  return std::to_string(I) + "," + std::to_string(J) + ",0"; 
}

std::string qn_to_str_3D( const int l, const int maxN )
{
  int I = l / (maxN * maxN); 
  int J = l / maxN - I*maxN;
  int K = l - I * (maxN * maxN) - J * maxN;
  return std::to_string(I) + "," + std::to_string(J) + "," + std::to_string(K); 
}

std::string qn_to_path_1D( const int l, const int /*maxN*/ )
{
  return std::to_string(l); 
}

std::string qn_to_path_2D( const int l, const int maxN )
{
  int I = l / maxN; 
  int J = l - I*maxN;
  return std::to_string(I) + "_" + std::to_string(J); 
}

std::string qn_to_path_3D( const int l, const int maxN )
{
  int I = l / (maxN * maxN); 
  int J = l / maxN - I*maxN;
  int K = l - I * (maxN * maxN) - J * maxN;
  return std::to_string(I) + "_" + std::to_string(J) + "_" + std::to_string(K); 
}

int main( int argc, char *argv[] )
{
  //const string HomePath = getenv ("HOME");

  std::map<int,std::function<std::string(const int,const int)>> qn_map;
  qn_map[1] = qn_to_str_1D; 
  qn_map[2] = qn_to_str_2D; 
  qn_map[3] = qn_to_str_3D;
  std::map<int,std::function<std::string(const int,const int)>> path_map;
  path_map[1] = qn_to_path_1D; 
  path_map[2] = qn_to_path_2D; 
  path_map[3] = qn_to_path_3D;

  string tmpstr;

  char base_folder[255];
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

  cxxopts::Options options("batch_case_generation_xml", "Creates a set parameter files for the stationary solvers.");
  
  options.add_options()
  ("d,dim",  "spatial dimension (1,2,3)", cxxopts::value<int>()->default_value("2") )
  ("N,NQN",  "max quantum number" , cxxopts::value<int>()->default_value("4") )
  ("s,dmu", "delta mu", cxxopts::value<string>()->default_value("0.1"))
  ("k,Ndmu", "number of delta mu steps", cxxopts::value<string>()->default_value("10"))
  ("g,gr", "global refinement", cxxopts::value<string>()->default_value("8"))
  ("f,fn", "folder name", cxxopts::value<string>()->default_value(string(base_folder)))
  ("positional", "Positional arguments: these are the arguments that are entered without an option", cxxopts::value<std::vector<std::string>>())
  ;
  
  options.parse_positional({"positional"});
  auto result = options.parse(argc, argv);

  int N = result["N"].as<int>();
  int dim = result["dim"].as<int>();

  //cout << result["positional"].as<std::vector<std::string>>()[0] << endl;

  const int maxN = int(pow(N,dim));
  for( int s=0; s<maxN; s++ )
  {
    sprintf( shellcmd, "mkdir -p %s/%s", result["fn"].as<string>().c_str(), path_map[dim](s,N).c_str() );
    system( shellcmd );
    sprintf( filename, "%s/%s/params.xml", result["fn"].as<string>().c_str(), path_map[dim](s,N).c_str() );

    pugi::xml_document doc;
    pugi::xml_node parameter_node = doc.append_child("PARAMETER");
    pugi::xml_node physics_node = parameter_node.append_child("PHYSICS");
    pugi::xml_node mesh_node = parameter_node.append_child("MESH");
    pugi::xml_node algorithm_node = parameter_node.append_child("ALGORITHM");      

    // PARAMETER childs
    pugi::xml_node node = parameter_node.append_child("FILENAME");
    node.append_child(pugi::node_pcdata).set_value("final.bin");

    node = parameter_node.append_child("guess_fct");
    node.append_child(pugi::node_pcdata).set_value("exp(-x^2-y^2-z^2)");

    // PHYSICS childs
    node = physics_node.append_child("QN1");
    node.append_child(pugi::node_pcdata).set_value(qn_map[dim](s,N).c_str());

    node = physics_node.append_child("omega");
    tmpstr = to_string(omega[0]) + "," + to_string(omega[1]) + "," + to_string(omega[2]);
    node.append_child(pugi::node_pcdata).set_value(tmpstr.c_str());

    node = physics_node.append_child("gs_1");
    node.append_child(pugi::node_pcdata).set_value("1,1");

    node = physics_node.append_child("mu");
    node.append_child(pugi::node_pcdata).set_value("5,5");

    // MESH childs
    node = mesh_node.append_child("DIM");
    node.append_child(pugi::node_pcdata).set_value(to_string(dim).c_str());
#if POTENTIAL==1
    node = mesh_node.append_child("xrange");
    node.append_child(pugi::node_pcdata).set_value("-10,10");
#endif 
#if POTENTIAL==2
    node = mesh_node.append_child("xrange");
    node.append_child(pugi::node_pcdata).set_value("0,20");
#endif 
    node = mesh_node.append_child("yrange");
    node.append_child(pugi::node_pcdata).set_value("-10,10");

    node = mesh_node.append_child("zrange");
    node.append_child(pugi::node_pcdata).set_value("-10,10");

    node = mesh_node.append_child("global_refinements");
    node.append_child(pugi::node_pcdata).set_value(result["gr"].as<string>().c_str());

    // ALGORTHM childs
    node = algorithm_node.append_child("ti");
    node.append_child(pugi::node_pcdata).set_value("1");

    node = algorithm_node.append_child("NA");
    node.append_child(pugi::node_pcdata).set_value("1000");

    node = algorithm_node.append_child("NK");
    node.append_child(pugi::node_pcdata).set_value("10");

    node = algorithm_node.append_child("dt");
    node.append_child(pugi::node_pcdata).set_value("0.001");

    node = algorithm_node.append_child("epsilon");
    node.append_child(pugi::node_pcdata).set_value("1e-5,1e-10");

    node = algorithm_node.append_child("dmu");
    node.append_child(pugi::node_pcdata).set_value(result["dmu"].as<string>().c_str());

    node = algorithm_node.append_child("Ndmu");
    node.append_child(pugi::node_pcdata).set_value(result["Ndmu"].as<string>().c_str());

    // add param node before the description
    //pugi::xml_node param = node.insert_child_before("param", descr);

    // add attributes to param node
    //param.append_attribute("name") = "version";
    //param.append_attribute("value") = 1.1;
    //param.insert_attribute_after("type", param.attribute("name")) = "float";

      cout << "Saving result: " << doc.save_file(filename) << endl;
  }
return EXIT_SUCCESS;
}
