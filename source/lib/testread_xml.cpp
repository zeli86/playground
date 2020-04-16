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

#include "pugixml.hpp"

#include <iostream>

int main()
{

    pugi::xml_document doc;
    if (!doc.load_file("params.xml")) return -1;

    //pugi::xml_node tools = doc.child("Profile").child("Tools");
    pugi::xml_node root = doc.child("PARAMETER");

    pugi::xpath_node_set tools = doc.select_nodes("/PARAMETER//MESH//*");

    for (pugi::xpath_node_set::const_iterator it = tools.begin(); it != tools.end(); ++it)
    {
        std::cout << "* ";

        pugi::xpath_node node = *it;
        std::cout << node.node().name() << ": " <<  node.node().child_value() <<"\n";
    }

    std::cout << "****************************************************\n";
    tools = doc.select_nodes("/PARAMETER/*/text()");

    for (pugi::xpath_node_set::const_iterator it = tools.begin(); it != tools.end(); ++it)
    {
        std::cout << "* ";

        pugi::xpath_node node = *it;
        std::cout << node.parent().name() << ": " <<  node.parent().child_value() <<"\n";
    }

    std::cout << "****************************************************\n";

    pugi::xml_node node2 = doc.child( "PARAMETER" ).child( "PHYSICS" ).child( "omega" );

    node2.last_child().set_value("1,1,1");

    tools = doc.select_nodes("//PARAMETER/PHYSICS/omega");
    pugi::xpath_node node = *tools.begin();
    std::cout << node.node().name() << ": " <<  node.node().child_value() <<"\n";    
    node.node().last_child().set_value("2,2,2");

    tools = doc.select_nodes("//PARAMETER/PHYSICS/omega");
    node = *tools.begin();
    std::cout << node.node().name() << ": " <<  node.node().child_value() <<"\n";    
    
 }
