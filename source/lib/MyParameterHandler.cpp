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


#include "MyParameterHandler.h"
#include "strtk.hpp"
#include <iostream>
#include <string>
#include <regex>
#include <stdexcept> 

MyParameterHandler::MyParameterHandler( const std::string filename )
{
  if( !m_xml_doc.load_file(filename.c_str()) )
  {
    std::cerr << "Critical error occured during loading of file " << filename << std::endl;
    throw;
  }

  populate_vconstants( "PHYSICS", m_map_physics );
  populate_vconstants( "MESH", m_map_mesh );
  populate_vconstants( "ALGORITHM", m_map_algorithm );
  populate_constants();
  populate_stringlists();
  populate_parameter();
}

void MyParameterHandler::populate_constants()
{
  m_map_constants.clear();

  double val;  
  std::string tmp, str;

  std::string querystr = "/PARAMETER//CONSTANTS//*";
  pugi::xpath_node_set tools = m_xml_doc.select_nodes(querystr.c_str());

  for (pugi::xpath_node_set::const_iterator it = tools.begin(); it != tools.end(); ++it)
  {
    pugi::xpath_node node = *it;

    str = node.node().name();
    tmp = node.node().child_value();

    try
    {
      val = stod(std::regex_replace(tmp,std::regex("\\s+"), ""));
    }
    catch( const std::invalid_argument& ia )
    {
       std::cerr << "Error Parsing xml file: Unable to convert " << tmp << " to double for element <" << str << "> in section CONSTANTS\n";
    }
    m_map_constants[str] = val;
  }
}

void MyParameterHandler::populate_vconstants( const std::string section, std::map<std::string,std::vector<double>>& mymap )
{
  mymap.clear();

  double val;  
  std::string tmp, str;
  std::vector<std::string> vec;

  std::string querystr = "/PARAMETER//" + section + "//*";
  pugi::xpath_node_set tools = m_xml_doc.select_nodes(querystr.c_str());

  for (pugi::xpath_node_set::const_iterator it = tools.begin(); it != tools.end(); ++it)
  {
    pugi::xpath_node node = *it;

    vec.clear();
    str = node.node().name();
    tmp = node.node().child_value();
    strtk::parse(tmp,",",vec);
      
    for( auto i : vec )
    {
      try
      {
        val = stod(std::regex_replace(i,std::regex("\\s+"), ""));
      }
      catch( const std::invalid_argument& ia )
      {
         std::cerr << "Error Parsing xml file: Unable to convert " << i << " to double for element <" << str << "> in section " << section << '\n';
      }
      mymap[str].push_back(val);
    }
  }
}

void MyParameterHandler::populate_stringlists()
{
  m_map_strings.clear();

  std::string tmp, str;
  std::vector<std::string> vec;

  std::string querystr = "/PARAMETER//STRINGLISTS//*";
  pugi::xpath_node_set tools = m_xml_doc.select_nodes(querystr.c_str());

  for (pugi::xpath_node_set::const_iterator it = tools.begin(); it != tools.end(); ++it)
  {
    pugi::xpath_node node = *it;

    vec.clear();
    str = node.node().name();
    tmp = node.node().child_value();
    if(!strtk::parse(tmp,",",vec)) continue;
      
    m_map_strings[str] = vec;
  }
}

void MyParameterHandler::populate_parameter()
{
  m_map_parameter.clear();

  std::string tmp, str;
  pugi::xpath_node_set tools = m_xml_doc.select_nodes("/PARAMETER/*/text()");

  for (pugi::xpath_node_set::const_iterator it = tools.begin(); it != tools.end(); ++it)
  {
    pugi::xpath_node node = *it;
    m_map_parameter[node.parent().name()] = node.parent().child_value(); 
  }
}

/*
void MyParameterHandler::Setup_muParser( mu::Parser& mup )
{
  mup.ClearConst();
  //mup.ClearFun();
  
  for( auto i : m_map_constants )
  {  
    //std::cout << i.first << ", " << i.second << std::endl;
    mup.DefineConst( i.first.c_str(), i.second );
  }
  
  mup.DefineFun("Heaviside", Heaviside, false);
  mup.DefineFun("rect", rect, false);
  mup.DefineFun("sign", sign, false);
}
*/

std::string MyParameterHandler::Get_Parameter( const std::string k )
{
  auto it = m_map_parameter.find(k);
  if( it == m_map_parameter.end() ) throw std::string( "Error: Could not find the key: " + k + " in section PARAMETER." );  
return (*it).second;
}

std::string MyParameterHandler::Get_String( const std::string k, const int p )
{
  auto it = m_map_strings.find(k);
  if( it == m_map_strings.end() ) throw std::string( "Error: Could not find the key: " + k + " in section STRINGLISTS." );
  if( p >= (*it).second.size() ) throw std::string( "Error in section STRINGLISTS: Access of vector element for key " + k + " is out of bounds." );
return (*it).second[p];
}

std::vector<std::string> MyParameterHandler::Get_AllStrings( const std::string k )
{
  auto it = m_map_strings.find(k);
  if( it == m_map_strings.end() ) throw std::string( "Error: Could not find the key " + k + " in section STRINGLISTS." ); 
return (*it).second;
}

void MyParameterHandler::Set_Physics( const std::string k, const std::vector<double>& newdata )
{
  std::string newstr;
  auto it = m_map_physics.find(k);
  if( it == m_map_physics.end() ) throw std::string( "Error: Could not find the key " + k + " in section PHYSICS." ); 

  m_map_physics.at(k) = newdata;

  for( int i=0; i<newdata.size()-1; ++i )
    newstr += std::to_string(newdata[i]) + ",";
  newstr += std::to_string(newdata[newdata.size()-1]);

  std::string querystr = "/PARAMETER//PHYSICS/" + k;
  pugi::xpath_node_set tools = m_xml_doc.select_nodes(querystr.c_str());
  pugi::xpath_node node = *tools.begin();
  node.node().last_child().set_value(newstr.c_str());
}

std::vector<double> MyParameterHandler::Get_Physics( const std::string k )
{
  auto it = m_map_physics.find(k);
  if( it == m_map_physics.end() ) throw std::string( "Error: Could not find the key " + k + " in section PHYSICS." ); 
return (*it).second;
}

double MyParameterHandler::Get_Physics( const std::string k, const int p )
{
  auto it = m_map_physics.find(k);
  if( it == m_map_physics.end() ) throw std::string( "Error: Could not find the key " + k + " in section PHYSICS." ); 
  if( p >= (*it).second.size() ) throw std::string( "Error in section PHYSICS: Access of vector element for key " + k + " is out of bounds." );
return (*it).second[p];
}

std::vector<double> MyParameterHandler::Get_Mesh( const std::string k )
{
  auto it = m_map_mesh.find(k);
  if( it == m_map_mesh.end() ) throw std::string( "Error: Could not find the key " + k + " in section MESH." ); 
return (*it).second;
}

double MyParameterHandler::Get_Mesh( const std::string k, const int p )
{
  auto it = m_map_mesh.find(k);
  if( it == m_map_mesh.end() ) throw std::string( "Error: Could not find the key " + k + " in section MESH." ); 
  if( p >= (*it).second.size() ) throw std::string( "Error in section MESH: Access of vector element for key " + k + " is out of bounds." );
return (*it).second[p];
}

std::vector<double> MyParameterHandler::Get_Algorithm( const std::string k )
{
  auto it = m_map_algorithm.find(k);
  if( it == m_map_algorithm.end() ) throw std::string( "Error: Could not find the key " + k + " in section ALGORITHM." ); 
return (*it).second;
}

double MyParameterHandler::Get_Algorithm( const std::string k, const int p )
{
  auto it = m_map_algorithm.find(k);
  if( it == m_map_algorithm.end() ) throw std::string( "Error: Could not find the key " + k + " in section ALGORITHM." ); 
  if( p >= (*it).second.size() ) throw std::string( "Error in section ALGORITHM: Access of vector element for key " + k + " is out of bounds." );
return (*it).second[p];
}

double MyParameterHandler::Get_Constant( const std::string k )
{
  auto it = m_map_constants.find(k);
  if( it == m_map_constants.end() ) throw std::string( "Error: Could not find the key " + k + " in section ALGORITHM." ); 
return (*it).second;
}

std::map<std::string,double> MyParameterHandler::Get_Constants_Map()
{
return m_map_constants;
}


int MyParameterHandler::Get_NA() 
{ 
  int retval=10;
  auto it = m_map_algorithm.find("NA");
  if( it != m_map_algorithm.end() ) retval = int((*it).second[0]);
return retval;
}

int MyParameterHandler::Get_NK() 
{ 
  int retval=10;
  auto it = m_map_algorithm.find("NK");
  if( it != m_map_algorithm.end() ) retval = int((*it).second[0]);
return retval;
}

void MyParameterHandler::SaveXMLFile( const std::string& filename )
{
  m_xml_doc.save_file( filename.c_str() );
}
