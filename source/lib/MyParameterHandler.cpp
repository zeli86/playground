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
#include <boost/tokenizer.hpp>
#include <iostream>
#include <string>
#include <regex>
#include <stdexcept>

MyParameterHandler::MyParameterHandler(const std::string& sXMLFilename)
{
  PopulatePropertyTree(sXMLFilename);
}


void MyParameterHandler::PopulatePropertyTree(const std::string& sXMLFilename)
{
  boost::property_tree::read_xml(sXMLFilename, m_oPropertyTree);
}

void MyParameterHandler::GetParameter(const std::string sXMLNodeName, std::vector<double>& vRetval)
{
  vRetval = m_oPropertyTree.get<std::vector<double>>(std::string("atus.parameter.") + sXMLNodeName);
}

void MyParameterHandler::GetParameter(const std::string sXMLNodeName, std::vector<int>& vRetval)
{
  vRetval = m_oPropertyTree.get<std::vector<int>>(std::string("atus.parameter.") + sXMLNodeName);
}

void MyParameterHandler::GetParameter(const std::string sXMLNodeName, double& rRetval)
{
  rRetval = m_oPropertyTree.get<double>(std::string("atus.parameter.") + sXMLNodeName);
}

void MyParameterHandler::GetParameter(const std::string sXMLNodeName, int& iRetval)
{
  iRetval = m_oPropertyTree.get<int>(std::string("atus.parameter.") + sXMLNodeName);
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

