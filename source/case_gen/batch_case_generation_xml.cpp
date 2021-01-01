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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <functional>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

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

  static std::map<int,std::function<std::string(const int,const int)>> qn_map;
  qn_map[1] = qn_to_str_1D; 
  qn_map[2] = qn_to_str_2D; 
  qn_map[3] = qn_to_str_3D;
  static std::map<int,std::function<std::string(const int,const int)>> path_map;
  path_map[1] = qn_to_path_1D; 
  path_map[2] = qn_to_path_2D; 
  path_map[3] = qn_to_path_3D;

  std::string aux_sString;

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

  boost::program_options::options_description oOptionsDesc{"Options"};

  oOptionsDesc.add_options()
  ("dim", boost::program_options::value<int>()->default_value(2),  "spatial dimension (1,2,3)" )
  ("NQN", boost::program_options::value<int>()->default_value(4),  "max quantum number" )
  ("dmu", boost::program_options::value<std::string>()->default_value("0.1"), "delta mu")
  ("Ndmu", boost::program_options::value<std::string>()->default_value("10"), "number of delta mu steps")
  ("gr", boost::program_options::value<std::string>()->default_value("8"), "global refinement")
  ("fn", boost::program_options::value<std::string>()->default_value(std::string(base_folder)), "folder name")
  ("help","Print help");

  boost::program_options::variables_map oVarMap;
  boost::program_options::store(boost::program_options::parse_command_line(argc, argv, oOptionsDesc), oVarMap);
  boost::program_options::notify(oVarMap);
  
  int N = oVarMap["NQN"].as<int>();
  int dim = oVarMap["dim"].as<int>();

  const int maxN = int(pow(N,dim));
  
  for( int s=0; s<maxN; s++ )
  {
    boost::property_tree::ptree oPropertyTree;
        
    sprintf( shellcmd, "mkdir -p %s/%s", oVarMap["fn"].as<std::string>().c_str(), path_map[dim](s,N).c_str() );
    system( shellcmd );
    sprintf( filename, "%s/%s/params.xml", oVarMap["fn"].as<std::string>().c_str(), path_map[dim](s,N).c_str() );

    // header data
    oPropertyTree.put("atus.header.version","1.0.0");
    oPropertyTree.put("atus.header.creation_date", "");
    
    // general parameter
    oPropertyTree.put("atus.parameter.general.filename", "final.bin");
    
    // physics parameter
    oPropertyTree.put("atus.parameter.physics.qn1", qn_map[dim](s,N) );
    aux_sString = std::to_string(omega[0]) + "," + std::to_string(omega[1]) + "," + std::to_string(omega[2]);     
    oPropertyTree.put("atus.parameter.physics.omega", aux_sString );
    oPropertyTree.put("atus.parameter.physics.gs_1", "1,1");
    
    // mesh parameter
    oPropertyTree.put("atus.parameter.spatial_dimension", oVarMap["dim"].as<std::string>());
    oPropertyTree.put("atus.parameter.mesh.xrange","0,20");
    oPropertyTree.put("atus.parameter.mesh.yrange","-10,10");
    oPropertyTree.put("atus.parameter.mesh.zrange","-10,10");
    oPropertyTree.put("atus.parameter.mesh.global_refinements", oVarMap["gr"].as<std::string>());
     
    // algorithm parameter
    oPropertyTree.put("atus.parameter.algorithm.epsilon", "1e-5,1e-10");
    oPropertyTree.put("atus.parameter.algorithm.dt","0.001");
    oPropertyTree.put("atus.parameter.algorithm.dmu", oVarMap["dmu"].as<std::string>());
    oPropertyTree.put("atus.parameter.algorithm.Ndmu", oVarMap["Ndmu"].as<std::string>());
    oPropertyTree.put("atus.parameter.algorithm.NA","1000");
    oPropertyTree.put("atus.parameter.algorithm.NK","10");    
    oPropertyTree.put("atus.parameter.algorithm.ti","1");    
    
//     PARAMETER childsi
//     node = parameter_node.append_child("guess_fct");
//     node.append_child(pugi::node_pcdata).set_value("exp(-x^2-y^2-z^2)");
// 
//     PHYSICS childs
//     node = physics_node.append_child("QN1");
//     node.append_child(pugi::node_pcdata).set_value();
// 
//     node = physics_node.append_child("mu");
//     node.append_child(pugi::node_pcdata).set_value("5,5");
// 

    boost::property_tree::write_xml(filename,oPropertyTree);
   }
return EXIT_SUCCESS;
}
