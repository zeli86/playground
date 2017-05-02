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

/** Želimir Marojević
 */

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/subscriptor.h>

#ifndef class_ParameterReader	
#define class_ParameterReader
using namespace dealii;

class ParameterReader : public Subscriptor
{
  public:
    ParameterReader(ParameterHandler &ph);
    bool read_parameters( const std::string = "params.prm" ); 
    void print_parameters( const std::string = "params.prm" );
    void declare_parameters();
   
    ParameterHandler &prm;
};
#endif