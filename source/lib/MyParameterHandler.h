/* * atus-pro testing - atus-pro testing playgroung
 * Copyright (C) 2017 Želimir Marojević <zelimir.marojevic@gmail.com>
 *
 * This file is part of atus-pro testing.
 *
 * atus-pro testing is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * atus-pro testing is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with atus-pro testing.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef __class_MyParameterHandler__
#define __class_MyParameterHandler__

#include "pugixml.hpp"
#include <string>
#include <vector>
#include <map>
#include <list>

class MyParameterHandler
{
public:
  MyParameterHandler( const std::string );
  virtual ~MyParameterHandler(){};
  
  std::string Get_Parameter( const std::string );

  void Set_Physics( const std::string, const std::vector<double>& );
  double Get_Physics( const std::string, const int );
  std::vector<double> Get_Physics( const std::string );

  double Get_Mesh( const std::string, const int );
  std::vector<double> Get_Mesh( const std::string );

  double Get_Algorithm( const std::string, const int );
  std::vector<double> Get_Algorithm( const std::string );

  void SaveXMLFile( const std::string& );

  int Get_NA(); 
  int Get_NK(); 

  //void Setup_muParser( mu::Parser& );
protected:
  void populate_vconstants( const std::string, std::map<std::string,std::vector<double>>& );
  void populate_parameter();
  
  pugi::xml_document m_xml_doc;
  
  std::map<std::string,std::vector<double>> m_map_physics;
  std::map<std::string,std::vector<double>> m_map_mesh;
  std::map<std::string,std::vector<double>> m_map_algorithm;
  std::map<std::string,std::string> m_map_parameter;
};

#endif
