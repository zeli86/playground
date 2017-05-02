//
// atus-pro testing - atus-pro testing playgroung
// Copyright (C) 2017 Želimir Marojević <zelimir.marojevic@gmail.com>
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

/** Želimir Marojević
 */

#include <deal.II/lac/generic_linear_algebra.h>

namespace LA
{
  using namespace dealii::LinearAlgebraPETSc;
}

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/constraint_matrix.h>

#include <deal.II/lac/petsc_parallel_sparse_matrix.h>
#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_precondition.h>
#include <deal.II/lac/sparsity_tools.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_q.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/fe_field_function.h>

#include <deal.II/base/utilities.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/function_parser.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/timer.h>

#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/solution_transfer.h>

#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <cstdio>

#include "functions_cs.h"
#include "global.h"
#include "mpi.h"
#include "MyParameterHandler.h"
#include "my_table.h"
#include "pugixml.hpp"
#include "gsl/gsl_sf.h"

extern "C"
{
  #include "mmio.h"
}

namespace HelperPrograms
{
  using namespace std;
  using namespace dealii;

  bool file_exist(const std::string& fileName)
  {
      std::ifstream infile(fileName);
      return infile.good();
  }
/*
  typedef double (*FUNC)( double );

  class Box_EF : public Function<2>
  {
  public:
    Box_EF ( const double Lx, const double Ly, const unsigned QX, const unsigned QY, const unsigned QZ )  : Function<2>() 
    {
      m_fak = sqrt(2/Lx)*sqrt(2/Ly);       
      m_dkx = M_PI/Lx;       
      m_dky = M_PI/Ly;

      m_qx = double(QX);
      m_qy = double(QY);

      m_fkt_x = &sin;
      if( QZ == 0 )
      {
        m_fkt_y = &cos;
        m_s = 0.5;
      }
      else
      {
        m_fkt_y = &sin;
        m_s = 1;
      }
    }

    virtual double value (const Point<2> &p, const unsigned int component = 0) const
    {
      return m_fak * (*m_fkt_x)((m_qx+1)*m_dkx*p[0]) * (*m_fkt_y)((m_qy+m_s)*m_dky*p[1]);
    }
    
    virtual void value_list (const std::vector<Point<2> > &points, std::vector<double> &values, const unsigned int component = 0) const
    {

    }

    protected:
      FUNC m_fkt_x;
      FUNC m_fkt_y;
      double m_fak;
      double m_dkx;
      double m_dky;
      double m_qx;
      double m_qy;
      double m_s;
  };
*/

  std::string basename(const std::string& pathname)
  {
    return {std::string(std::find_if(pathname.rbegin(), pathname.rend(),[](char c) { return c == '/'; }).base(),pathname.end())};
  }  

  std::string dirname(std::string source)
  {
    source.erase(std::find(source.rbegin(), source.rend(), '/').base(), source.end());
    return source;
  }    

  template <int dim>
  class MySolver
  {
  public:
    MySolver( const std::string&, const std::string& );
    ~MySolver();

    void run ( string );

  protected:
    void make_grid();
    void setup_system();
    void setup_boundary_ids();
    void compute_matrix_element( const unsigned, const unsigned, const unsigned, const unsigned, double (&)[4] );
  
    void load( string );

    MyParameterHandler m_ph;
    MPI_Comm mpi_communicator;
    parallel::distributed::Triangulation<dim> triangulation;
    FE_Q<dim> fe;
    DoFHandler<dim> dof_handler;
    IndexSet locally_owned_dofs;
    IndexSet locally_relevant_dofs;
    ConstraintMatrix constraints;

    LA::MPI::Vector m_Psi;

    ConditionalOStream pcout;

    bool m_root;

    double m_xmin, m_xmax;
    double m_ymin, m_ymax;
    double m_zmin, m_zmax;
    double m_mu;
    double m_gs;
    vector<double> m_omega;

    unsigned m_QN1[3];
    
    int m_rank;
  };

  template <int dim>
  MySolver<dim>::MySolver ( const std::string& path, const std::string& xmlfilename ) 
    : 
    m_ph(path + xmlfilename),
    mpi_communicator (MPI_COMM_WORLD),
    triangulation (mpi_communicator, typename Triangulation<dim>::MeshSmoothing(Triangulation<dim>::smoothing_on_refinement|Triangulation<dim>::smoothing_on_coarsening)),
    fe (gl_degree_fe),
    dof_handler (triangulation),
    pcout (cout, (Utilities::MPI::this_mpi_process(mpi_communicator) == 0))
  {
    m_omega = m_ph.Get_Physics("omega");
    m_gs = m_ph.Get_Physics("gs_1",0);
    m_QN1[0] = int(m_ph.Get_Physics("QN1",0));
    m_QN1[1] = int(m_ph.Get_Physics("QN1",1));
    m_QN1[2] = int(m_ph.Get_Physics("QN1",2));

    m_xmin = m_ph.Get_Mesh("xrange",0);
    m_xmax = m_ph.Get_Mesh("xrange",1);
    m_ymin = m_ph.Get_Mesh("yrange",0);
    m_ymax = m_ph.Get_Mesh("yrange",1);
    m_zmin = m_ph.Get_Mesh("zrange",0);
    m_zmax = m_ph.Get_Mesh("zrange",1);

    pugi::xml_document doc;
    if (!doc.load_file("info.xml")) throw;
    
    m_mu = stod(doc.child( "INFO" ).child( "MU" ).child_value());

    m_root = (Utilities::MPI::this_mpi_process(mpi_communicator) == 0);
    MPI_Comm_rank(mpi_communicator, &m_rank);
  }

  template <int dim>
  MySolver<dim>::~MySolver ()
  {
    dof_handler.clear ();
  }

  template <int dim>
  void MySolver<dim>::make_grid ()
  {
#if DIMENSION==2
    Point<dim,double> pt1( m_xmin, m_ymin );
    Point<dim,double> pt2( m_xmax, m_ymax );
#endif
#if DIMENSION==3
    Point<dim,double> pt1( m_xmin, m_ymin, m_zmin );
    Point<dim,double> pt2( m_xmax, m_ymax, m_zmax );
#endif

    GridGenerator::hyper_rectangle(triangulation, pt2, pt1);
    triangulation.refine_global(1);
  }
  
  template <int dim>
  void MySolver<dim>::setup_system()
  {
    dof_handler.distribute_dofs (fe);

    locally_owned_dofs = dof_handler.locally_owned_dofs ();
    DoFTools::extract_locally_relevant_dofs (dof_handler, locally_relevant_dofs);

    m_Psi.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);

    vector<bool> mask (dof_handler.get_fe().n_components(), true );
    
    constraints.clear ();
    constraints.reinit (locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints (dof_handler, constraints);
    VectorTools::interpolate_boundary_values (dof_handler, 1, ZeroFunction<dim>(), constraints, ComponentMask(mask));
    constraints.close ();
  }

  template <int dim>
  void MySolver<dim>::compute_matrix_element( const unsigned q1, const unsigned q2, const unsigned q3, const unsigned q4, double (&retval)[4] )
  {
    const QGauss<dim> quadrature_formula(fe.degree);
    const QGauss<dim> quadrature_formula_2(fe.degree+1);
   
    FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_JxW_values|update_quadrature_points);
    FEValues<dim> fe_values_2 (fe, quadrature_formula_2, update_values|update_JxW_values|update_quadrature_points);
    
    unsigned QN1[] = {q1,q2,m_QN1[2]};
    unsigned QN2[] = {q3,q4,m_QN1[2]};
    CEigenfunctions Ef1( QN1, m_omega );
    CEigenfunctions Ef2( QN2, m_omega );    

    const unsigned dofs_per_cell = fe.dofs_per_cell;
    const unsigned n_q_points = quadrature_formula.size();
    const unsigned n_q_points_2 = quadrature_formula_2.size();

    vector<double> Psi(n_q_points);
    vector<double> Psi_2(n_q_points_2);
    Point<dim> pt;
    double tmp[4] = {};

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for ( ; cell!=endc; cell++ )
    {
      if( cell->is_locally_owned() )
      {
        fe_values.reinit (cell);
        fe_values.get_function_values( m_Psi, Psi );

        fe_values_2.reinit (cell);
        fe_values_2.get_function_values( m_Psi, Psi_2 );

        for( unsigned qp=0; qp<n_q_points; qp++ )
        {
          pt = fe_values.quadrature_point(qp);
          double tmp2 = fe_values.JxW(qp) * fabs(pt[1]) * Psi[qp] * Psi[qp] * Ef1.value(pt) * Ef2.value(pt);
          
          tmp[0] += tmp2;
          tmp[1] += 3 * tmp2;
        }

        for( unsigned qp=0; qp<n_q_points_2; qp++ )
        {
          pt = fe_values_2.quadrature_point(qp);
          double tmp2 = fe_values_2.JxW(qp) * fabs(pt[1]) * Psi_2[qp] * Psi_2[qp] * Ef1.value(pt) * Ef2.value(pt);
          
          tmp[2] += tmp2;
          tmp[3] += 3 * tmp2;
        }
      }
    }

    tmp[0] *= m_gs;
    tmp[1] *= m_gs;
    tmp[2] *= m_gs;
    tmp[3] *= m_gs;

    MPI_Allreduce( tmp, retval, 4, MPI_DOUBLE, MPI_SUM, mpi_communicator);
/*
    if( QN1[0] == QN2[0] && QN1[1] == QN2[1] && QN1[2] == QN2[2] )
    {
      double aizero = fabs(gsl_sf_airy_zero_Ai(QN1[0]+1));
      double E = pow( m_omega[0], 2.0/3.0 ) * aizero + m_omega[1]*(1+2*QN1[1]+QN1[2]);
      retval[0] += (E-m_mu);
      retval[1] += (E-m_mu);
      retval[2] += (E-m_mu);
      retval[3] += (E-m_mu);
    }
  */
  }

  template <int dim>
  void MySolver<dim>::run( string filename )
  {
    load(filename);

    const int N = 4;
    const int Ntot = N*N;

    double c[4];

    FullMatrix<double> block_1(Ntot, Ntot);    
    FullMatrix<double> block_2(Ntot, Ntot);    
    FullMatrix<double> block_1b(Ntot, Ntot);    
    FullMatrix<double> block_2b(Ntot, Ntot);    
    
    for( int i=0; i<Ntot; i++ )
    {
      for( int j=i; j<Ntot; j++ )
      {
        int I = i / N; 
        int J = i - I*N; 
        int K = j / N;
        int L = j - K*N; 
        compute_matrix_element ( I, J, K, L, c );

        block_1( i, j ) = c[0];
        block_1( j, i ) = c[0];
        block_2( i, j ) = c[1];
        block_2( j, i ) = c[1];

        block_1b( i, j ) = c[2];
        block_1b( j, i ) = c[2];
        block_2b( i, j ) = c[3];
        block_2b( j, i ) = c[3];
      }
      pcout << i << endl;
    }
    
    if( m_root )
    {
      MM_typecode matcode;
      mm_initialize_typecode(&matcode);
      mm_set_matrix(&matcode);
      mm_set_coordinate(&matcode);
      mm_set_real(&matcode);      
      mm_set_sparse(&matcode);   
      mm_set_general(&matcode);   

      //size_t pos = filename.find_last_of("/") + 1;
      //string filename2 = filename;
      //filename2.erase(pos,filename.length());

      //cout << filename2;

      string filename2 = "mat_" + to_string(Ntot) + "_1.mtx"; 

      FILE * fh = fopen( filename2.c_str(), "w+");

      mm_write_banner(fh, matcode); 
      mm_write_mtx_crd_size(fh, 2*Ntot, 2*Ntot, 2*Ntot*Ntot);

      for( int i=0; i<Ntot; i++ )
      {
        for( int j=0; j<Ntot; j++ )
        {
            fprintf(fh, "%d %d %.15e\n", i+Ntot+1, j+1, block_1(i,j));
        }    
      }
      for( int i=0; i<Ntot; i++ )
      {
        for( int j=0; j<Ntot; j++ )
        {
            fprintf(fh, "%d %d %.15e\n", i+1, j+Ntot+1, -block_2(i,j));
        }    
      }      
      fclose(fh);

      filename2 = "mat_" + to_string(Ntot) + "_1b.mtx"; 

      fh = fopen( filename2.c_str(), "w+");

      mm_write_banner(fh, matcode); 
      mm_write_mtx_crd_size(fh, 2*Ntot, 2*Ntot, 2*Ntot*Ntot);

      for( int i=0; i<Ntot; i++ )
      {
        for( int j=0; j<Ntot; j++ )
        {
            fprintf(fh, "%d %d %.15e\n", i+Ntot+1, j+1, block_1b(i,j));
        }    
      }
      for( int i=0; i<Ntot; i++ )
      {
        for( int j=0; j<Ntot; j++ )
        {
            fprintf(fh, "%d %d %.15e\n", i+1, j+Ntot+1, -block_2b(i,j));
        }    
      }      
      fclose(fh);
    }

    const unsigned lev = triangulation.n_global_levels ();
    const unsigned subs = unsigned(pow(2,double(lev-1)));

    const double dx = fabs(m_xmax-m_xmin)/double(subs);
    const double dy = fabs(m_ymax-m_ymin)/double(subs);
//    const double dkx = M_PI/fabs(m_xmax-m_xmin);
//    const double dky = M_PI/fabs(m_ymax-m_ymin);

    pcout << "(dx, dy) == (" << dx << ", " << dy << ")\n";
    pcout << "(dx, dy) == (" << dx << ", " << dy << ")\n";
  }

  template<int dim>
  void MySolver<dim>::load( string filename )
  {
    make_grid();
    triangulation.load( filename.c_str() );
    
    setup_system();

    parallel::distributed::SolutionTransfer<dim,LA::MPI::Vector> solution_transfer(dof_handler);
    solution_transfer.deserialize(m_Psi);
    
    LA::MPI::Vector tmp;
    tmp.reinit (locally_owned_dofs, mpi_communicator);
    
    tmp=m_Psi;
    constraints.distribute(tmp);
    m_Psi=tmp;
  }
} // end of namespace 

int main ( int argc, char *argv[] )
{
  using namespace dealii;
  deallog.depth_console (0);

  //std::string filename = argv[1];

  //size_t pos = filename.find_last_of("/") + 1;
  //string path = filename;
  //path.erase(pos,filename.length());

  //if( HelperPrograms::file_exist( path + "mat_100.mtx") ) 
  //  return EXIT_SUCCESS;

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv );
  {
    HelperPrograms::MySolver<DIMENSION> solver("", "params.xml");
    solver.run("final.bin");
  }  
return EXIT_SUCCESS;
}
