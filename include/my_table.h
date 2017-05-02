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

#ifndef _MYTABLE_
#define _MYTABLE_

#include <vector>
#include <string>
#include <map>

using namespace std;

typedef map<string,double> columns;
typedef vector<columns> table;

class MyTable
{
  public:
    MyTable();
    ~MyTable();

    friend ostream& operator<<( ostream&, MyTable& );

    void clear();
    void load( string );
    columns& new_line();
    void insert( columns&, const string, const double );
    void dump_2_file( const string );
    
    static const string COUNTER;
    static const string RES;
    static const string RESP;
    static const string RES_OVER_RESP;
    static const string RES_2;
    static const string RESP_2;
    static const string RES_OVER_RESP_2;
    static const string L2_NORM_PSI_REF;
    static const string INF_NORM_PSI_REF;
    static const string PARTICLE_NUMBER;
    static const string PARTICLE_NUMBER2;
    static const string MU;
    static const string MU2;
    static const string GS;
    static const string GS2;
    static const string l2norm_t;
    static const string Delta_l2norm_t;
    static const string t1;
    static const string t2;
    static const string t3;
    static const string t4;
    static const string ev1;
    static const string ev2;
    static const string ev3;
    static const string ev4;
    static const string total_no_cells;
    static const string total_no_active_cells;
    static const string STEPSIZE;
    static const string STATUS;
    static const string time;
    static const string ev_position_x;
    static const string ev_position_y;
    static const string ev_position_z;
    static const string ev_momentum_x;
    static const string ev_momentum_y;
    static const string ev_momentum_z;

    void save_txt( string );

    table m_table;
  protected:
    vector<string> m_order;
};

ostream& operator<<( ostream& stream, MyTable& obj );
#endif