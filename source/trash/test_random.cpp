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

#include <random>
#include <iostream>
#include <vector>
#include <fstream>
 
int main()
{
    std::vector<double> t1, t2;
    std::random_device rd, rd2;
    std::mt19937 gen(rd()), gen2(rd2());
    std::uniform_real_distribution<> dis(-5, 5);
    std::uniform_real_distribution<> dis2(-5, 5);
    for( unsigned n=0; n<100; n++ ) 
    {
        t1.push_back(dis(gen));
        t2.push_back(dis2(gen2));
    }
    
    std::ofstream ofs( "test.txt" );
    for( unsigned n=0; n<t1.size(); n++ )
      ofs << t1[n] << "\t" << t2[n] << "\n";

return 0;
}