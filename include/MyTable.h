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

#include <vector>
#include <string>
#include <cstdint>

#include <boost/config.hpp>
#include <boost/utility.hpp>
#include <boost/call_traits.hpp>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/composite_key.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/ordered_index.hpp>

class MyTableRow
{
  public:
    uint64_t m_iIndex = 0;
    double m_rResidum = 0;
    double m_rResidum2 = 0;
    double m_rMu = 0;
    double m_rMu2 = 0;
};

using tTable = boost::multi_index::multi_index_container<MyTableRow,
  boost::multi_index::indexed_by<
  boost::multi_index::ordered_unique<BOOST_MULTI_INDEX_MEMBER(MyTableRow,uint64_t,m_iIndex)>
>>;

class MyTable
{
  public:
    MyTable() = default;
    ~MyTable();

    //friend ostream& operator<<( ostream&, MyTable& );

    void clear();
    void load( const std::string& );

    void insert( const MyTableRow& );

    static const std::string COUNTER;
    static const std::string RES;
    static const std::string RESP;
    static const std::string RES_OVER_RESP;
    static const std::string RES_2;
    static const std::string RESP_2;
    static const std::string RES_OVER_RESP_2;
    static const std::string L2_NORM_PSI_REF;
    static const std::string INF_NORM_PSI_REF;
    static const std::string PARTICLE_NUMBER;
    static const std::string PARTICLE_NUMBER2;
    static const std::string MU;
    static const std::string MU2;
    static const std::string GS;
    static const std::string GS2;
    static const std::string l2norm_t;
    static const std::string Delta_l2norm_t;
    static const std::string t1;
    static const std::string t2;
    static const std::string t3;
    static const std::string t4;
    static const std::string ev1;
    static const std::string ev2;
    static const std::string ev3;
    static const std::string ev4;
    static const std::string total_no_cells;
    static const std::string total_no_active_cells;
    static const std::string STEPSIZE;
    static const std::string STATUS;
    static const std::string time;
    static const std::string ev_position_x;
    static const std::string ev_position_y;
    static const std::string ev_position_z;
    static const std::string ev_momentum_x;
    static const std::string ev_momentum_y;
    static const std::string ev_momentum_z;

    void save_txt( const std::string& );

    tTable m_miTable;
};

//ostream& operator<<( ostream& stream, MyTable& obj );
