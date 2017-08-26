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

#ifndef __REF_PT_LIST__
#define __REF_PT_LIST__

#include <vector>
#include <cmath>
#include <algorithm>
#include <mpi.h>

namespace MyUtils
{
  template<int dim>
  class ref_pt_list_item
  {
    public:
      ref_pt_list_item()
      { 
        f = std::numeric_limits<double>::max();
        failed = true;        
      };
      Point<dim,double> ti;
      Point<dim,double> t;
      Point<dim,double> df;
      double f;
      bool failed;
      int status;

      void DumpItem( std::ostream& );
      double l2norm_t() const;
      double l2norm_df() const;
  };

  template <int dim>
  bool operator==( const ref_pt_list_item<dim>& rhs, const ref_pt_list_item<dim>& lhs )
  {
    return (fabs(rhs.l2norm_t()-lhs.l2norm_t()) < 1e-2);
  }

  template <int dim>
  bool operator<( const ref_pt_list_item<dim>& rhs, const ref_pt_list_item<dim>& lhs )
  {
    return (rhs.l2norm_t() < lhs.l2norm_t());
  }

  template <int dim>
  void ref_pt_list_item<dim>::DumpItem( std::ostream& out )
  {
    out << std::scientific;
    out << "|t|_l2 = " << l2norm_t() << ", |df|_l2 = " << l2norm_df() << ", f = " << f << ", " << status << ", " << ( failed == true ? "true" : "false" )  << std::endl;
  }

  template <int dim>
  double ref_pt_list_item<dim>::l2norm_t() const
  {
    double retval = 0;
    for( int i=0; i<dim; i++ )
      retval += t[i]*t[i];
  return sqrt(retval);
  }  

  template <int dim>
  double ref_pt_list_item<dim>::l2norm_df() const
  {
    double retval = 0;
    for( int i=0; i<dim; i++ )
      retval += df[i]*df[i];
  return sqrt(retval);
  }  
  
  template <int dim>
  class ref_pt_list
  {
    public:
      ref_pt_list( const unsigned int=10, const double=10.0 );
      
      std::vector<ref_pt_list_item<dim>> m_list;
      void reset( const unsigned int=10, const double=10.0 );
      void remove_f_non_zero( const double=1e-12 );
      void remove_zero_ref_points( const double=1.0 );
      void remove_duplicates();
      void broadcast( MPI_Comm, const int );
      void Dump( std::ostream& );
      void Dump( const string );
      const int Get_dim() const { return dim; };
    protected:
      void create_initial_list();
      void clear();

      std::vector<unsigned int> m_vertex;
      std::vector<double> m_coords;
      std::vector<double> m_shifted_coords;
      std::vector<std::vector<double>*> m_coords_array;
      unsigned int m_N;
      double m_ti_max;
      double m_dti;
  };

  template<int dim>
  void ref_pt_list<dim>::remove_f_non_zero( const double threshold )
  {
    typename std::vector<ref_pt_list_item<dim>>::iterator it = m_list.begin();
    for( ; it != m_list.end(); )
    {
      if( fabs((*it).f) > threshold )
        it = m_list.erase(it);
      else
        it++;
    }
  }

  template<int dim>
  void ref_pt_list<dim>::remove_zero_ref_points( const double threshold )
  {
    typename std::vector<ref_pt_list_item<dim>>::iterator it = m_list.begin();
    for( ; it != m_list.end(); )
    {
      if( (*it).l2norm_t() < threshold || (*it).failed )
        it = m_list.erase(it);
      else
        it++;
    }
  }
  
  template<int dim>
  void ref_pt_list<dim>::remove_duplicates()
  {
    std::sort( m_list.begin(), m_list.end() );
    m_list.erase( std::unique(m_list.begin(), m_list.end()), m_list.end());
  }

  template<int dim>
  ref_pt_list<dim>::ref_pt_list( unsigned N, double max )
  {
    m_N = N;
    m_ti_max = max;
    m_dti = m_ti_max / double(m_N);
    
    for( unsigned i=1; i<=m_N; i++ )
    {
      double x = double(i)*m_dti;
      m_coords.push_back(x);
    }

    for( unsigned i=0; i<m_N; i++ )
    {
      double x = double(2*i)*m_dti - m_ti_max + m_dti;
      m_shifted_coords.push_back(x);
    }

    m_coords_array.push_back(&m_coords);
    if( dim > 1 )
    {
      for( int i=1; i<dim; i++ )
        m_coords_array.push_back(&m_shifted_coords); 
    }
    create_initial_list();
  }
  
  template<int dim>
  void ref_pt_list<dim>::reset( unsigned N, double max )
  {
    clear();
    
    m_N = N;
    m_ti_max = max;
    m_dti = m_ti_max / double(m_N);
    
    for( unsigned i=1; i<=m_N; i++ )
    {
      double x = double(i)*m_dti;
      m_coords.push_back(x);
    }

    for( unsigned i=0; i<m_N; i++ )
    {
      double x = double(2*i)*m_dti - m_ti_max + m_dti;
      m_shifted_coords.push_back(x);
    }

    m_coords_array.push_back(&m_coords);
    if( dim > 1 )
    {
      for( int i=1; i<dim; i++ )
        m_coords_array.push_back(&m_shifted_coords); 
    }
    create_initial_list();
  }

  template<int dim>
  void ref_pt_list<dim>::clear()
  {
    m_coords.clear();
    m_shifted_coords.clear();
    m_coords_array.clear();
    m_list.clear();
  }
  
  template<int dim>
  void ref_pt_list<dim>::Dump( const string filename )
  { 
    ofstream out( filename );
    for( auto i : m_list )
    {
      i.DumpItem( out ); 
    }
  }

  template<int dim>
  void ref_pt_list<dim>::Dump( std::ostream& o  )
  { 
    ios init(NULL);
    init.copyfmt(cout);    
    for( auto i : m_list )
    {
      i.DumpItem( o ); 
    }
    std::cout.copyfmt(init);
  }
  
  template<int dim>
  void ref_pt_list<dim>::create_initial_list()
  {
    vector<int> perms(dim, 0);

    bool done;
    do
    {
      ref_pt_list_item<dim> item;
      for( int i=0; i<dim; i++ )
        item.ti[i] = (*m_coords_array[i])[perms[i]];
      m_list.push_back(item);

      done = true;
      for( int i=dim-1; i>=0; i-- ) 
      {
        if( ++perms[i] > m_N-1 ) 
        {
          perms[i] = 0;
          continue;
        } 
        else 
        {
          done = false; 
          break;
        }
      }
    } while (!done);
  }
  
  template<int dim> 
  void ref_pt_list<dim>::broadcast( MPI_Comm comm, const int rank )
  {
    double help1[dim];
    int len = m_list.size();
    MPI_Bcast( &len, 1, MPI_INT, 0, comm );

    if( rank == 0 )
    {
      typename std::vector<ref_pt_list_item<dim>>::iterator it = m_list.begin();
      for( ; it!=m_list.end(); it++ )
      {
        for( int i=0; i<dim; i++ ) help1[i] = (*it).ti[i];
        MPI_Bcast( help1, dim, MPI_DOUBLE, 0, comm );
        for( int i=0; i<dim; i++ ) help1[i] = (*it).t[i];
        MPI_Bcast( help1, dim, MPI_DOUBLE, 0, comm );
      }
    }
    else
    {
      m_list.clear();
      for( int l=0; l<len; l++ )
      {
        MyUtils::ref_pt_list_item<dim> newitem;
        MPI_Bcast( help1, dim, MPI_DOUBLE, 0, comm );
        for( int i=0; i<dim; i++ ) newitem.ti[i] = help1[i];
        MPI_Bcast( help1, dim, MPI_DOUBLE, 0, comm );
        for( int i=0; i<dim; i++ ) newitem.t[i] = help1[i];
        m_list.push_back( newitem );
      }
    }
  }
}
#endif