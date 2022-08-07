/* * atus-pro testing - atus-pro testing playgroung
 * Copyright (C) 2020 Želimir Marojević <zelimir.marojevic@gmail.com>
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

#pragma once

#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <list>
#include <map>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

namespace mylib
{
  namespace xml_translators
  {

    template<typename T> struct container
    {
      // types
      typedef T internal_type;
      typedef T external_type;

      boost::optional<T> get_value(const std::string& str) const
      {
        if (str.empty())
        {
          return boost::none;
        }

        T values;
        std::stringstream ss(str);

        typename T::value_type temp_value;
        while (ss >> temp_value)
        {
          values.insert(values.end(), temp_value);
        }

        return boost::make_optional(values);
      }

      boost::optional<std::string> put_value(const T& b)
      {
        std::stringstream ss;
        int32_t i = 0;
        for (auto v : b)
        {
          ss << (i++ < b.size() ? " " : "") << v;
        }
        return boost::make_optional(ss.str());
      }
    };

  }
}

namespace boost
{
  namespace property_tree
  {
    template<typename ch, typename traits, typename alloc, typename T>
    struct translator_between<std::basic_string<ch, traits, alloc>, std::vector<T>>
    {
      typedef mylib::xml_translators::container<std::vector<T>> type;
    };

    template<typename ch, typename traits, typename alloc, typename T>
    struct translator_between<std::basic_string<ch, traits, alloc>, std::list<T>>
    {
      typedef mylib::xml_translators::container<std::list<T>> type;
    };
  }
}

class MyParameterHandler
{
public:
  explicit MyParameterHandler(const std::string&);
  virtual ~MyParameterHandler() {};

  void get(const std::string, std::vector<double>&);
  void get(const std::string, std::vector<int>&);
  void get(const std::string, double&);
  void get(const std::string, int&);

  //void Setup_muParser( mu::Parser& );
protected:

  void PopulatePropertyTree(const std::string&);

  boost::property_tree::ptree m_oPropertyTree;
};


