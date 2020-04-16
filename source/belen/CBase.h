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


#include "MyParameterHandler.h"

#ifndef __CBASE_H__
#define __CBASE_H__

#define STR1(x) #x
#define STR2(x) STR1(x)

enum Status { SUCCESS, FAILED, ZERO_SOL, SLOW_CONV };

template <int dim> 
class CBase
{
  public:
    CBase( const std::string );
    virtual ~CBase() {};
 
    void dump_info_xml( const string="" );
    MPI_Comm mpi_communicator;
  protected:
    double m_xmin, m_xmax;
    double m_ymin, m_ymax;

    bool m_root;
    int m_rank;

    unsigned m_counter;
    unsigned m_global_refinement;
    unsigned m_total_no_cells;
    unsigned m_total_no_active_cells;    
    
    ofstream m_computing_timer_log;
    TimerOutput m_computing_timer;    
    MyParameterHandler m_ph;
    ConditionalOStream pcout;

    string m_mass_dist;
};

template <int dim>
CBase<dim>::CBase( const std::string xmlfilename  ) 
  : 
  mpi_communicator(MPI_COMM_WORLD), 
  m_computing_timer_log("benchmark.txt"),
  m_computing_timer(mpi_communicator, m_computing_timer_log, TimerOutput::summary, TimerOutput:: cpu_and_wall_times ), 
  m_ph(xmlfilename),
  pcout (cout, (Utilities::MPI::this_mpi_process(mpi_communicator) == 0))
{
  try 
  {

    m_xmin = m_ph.Get_Mesh("xrange",0);
    m_xmax = m_ph.Get_Mesh("xrange",1);
    m_ymin = m_ph.Get_Mesh("yrange",0);
    m_ymax = m_ph.Get_Mesh("yrange",1);

    m_global_refinement = m_ph.Get_Algorithm("global_refinement",0);
  }
  catch( const std::string info )
  {
    std::cerr << info << endl;
    MPI_Abort( mpi_communicator, 0 );
  }

  /*
  m_mu[0] = m_ph.get_double("mu");
  m_mu[1] = m_ph.get_double("mu2");
  m_guess_str = m_ph.get( "guess_fct" );
  */
  
  m_root = (Utilities::MPI::this_mpi_process(mpi_communicator) == 0);
  MPI_Comm_rank(mpi_communicator, &m_rank);
  
  m_counter=0;
  //m_final_error=0;
}

template <int dim>
void CBase<dim>::dump_info_xml ( const string path )
{
  string filename = path + "info.xml";

  wofstream fs2;
  fs2.open(filename);

  locale utf8_locale("en_US.UTF8");
  fs2.imbue(utf8_locale); 
  fs2 << L"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<INFO>\n";
  fs2 << L"<XMIN>" << m_xmin << L"</XMIN>\n";
  fs2 << L"<XMAX>" << m_xmax << L"</XMAX>\n";
  fs2 << L"<YMIN>" << m_ymin << L"</YMIN>\n";
  fs2 << L"<YMAX>" << m_ymax << L"</YMAX>\n";
  fs2 << L"<REVISION>" << STR2(GIT_SHA1) << L"</REVISION>\n";
  fs2 << L"</INFO>\n";
  fs2.close();
}
#endif