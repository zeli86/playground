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
//#include <boost/utility.hpp>
//#include <boost/call_traits.hpp>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/composite_key.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/ordered_index.hpp>

template <int iNoOfComponents>
class MyTableRow
{
public:
  uint64_t m_iIndex = 0;
  double m_rEigenvalue[iNoOfComponents] = {};
  double m_rNonlinearity[iNoOfComponents] = {};
  double m_rResidue[iNoOfComponents] = {};
  double m_rParticleNumber[iNoOfComponents] = {};
  double m_rExpectationValuePositionX[iNoOfComponents] = {};
  double m_rExpectationValuePositionY[iNoOfComponents] = {};
  double m_rExpectationValuePositionZ[iNoOfComponents] = {};
  double m_rExpectationValueMomentumX[iNoOfComponents] = {};
  double m_rExpectationValueMomentumY[iNoOfComponents] = {};
  double m_rExpectationValueMomentumZ[iNoOfComponents] = {};
  double m_rExpectationValueEnergyX[iNoOfComponents] = {};
  double m_rExpectationValueEnergyY[iNoOfComponents] = {};
  double m_rExpectationValueEnergyZ[iNoOfComponents] = {};
  double m_rL2Norm[iNoOfComponents] = {};
  double m_rInfNorm[iNoOfComponents] = {};
};

template <int iNoOfComponents>
class tTable : public boost::multi_index::multi_index_container<MyTableRow<iNoOfComponents>,
  boost::multi_index::indexed_by<boost::multi_index::ordered_unique<BOOST_MULTI_INDEX_MEMBER(MyTableRow<iNoOfComponents>, uint64_t, m_iIndex)>>>
  {};

template <int iNoOfComponents>
std::ostream& operator<<( std::ostream& oOutput, MyTableRow<iNoOfComponents>& oRow )
{
   oOutput << oRow.m_iIndex;
   for( int32_t i=0; i<iNoOfComponents; ++i )
   {
     oOutput << '\t' << oRow.m_rEigenvalue[i];
     oOutput << '\t' << oRow.m_rNonlinearity[i];
     oOutput << '\t' << oRow.m_rResidue[i];
     oOutput << '\t' << oRow.m_rParticleNumber[i];
     oOutput << '\t' << oRow.m_rExpectationValuePositionX[i];
     oOutput << '\t' << oRow.m_rExpectationValuePositionY[i];
     oOutput << '\t' << oRow.m_rExpectationValuePositionZ[i];
     oOutput << '\t' << oRow.m_rExpectationValueMomentumX[i];
     oOutput << '\t' << oRow.m_rExpectationValueMomentumY[i];
     oOutput << '\t' << oRow.m_rExpectationValueMomentumZ[i];
     oOutput << '\t' << oRow.m_rExpectationValueEnergyX[i];
     oOutput << '\t' << oRow.m_rExpectationValueEnergyY[i];
     oOutput << '\t' << oRow.m_rExpectationValueEnergyZ[i];
     oOutput << '\t' << oRow.m_rL2Norm[i];
     oOutput << '\t' << oRow.m_rInfNorm[i];    
   }
   oOutput << '\n';
}

template <int iNoOfComponents>
std::ostream& operator<<( std::ostream& oOutput, tTable<iNoOfComponents>& oTable )
{
  for( const auto& oRow : oTable )
  {
    oOutput << oRow;
  }
}
