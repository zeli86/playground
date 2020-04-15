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

#include "global.h"
#include "mpi.h"
#include "MyParameterHandler.h"
#include "my_table.h"

namespace HelperPrograms
{
  std::string basename(const std::string& pathname)
  {
    return {std::string(std::find_if(pathname.rbegin(), pathname.rend(),[](char c) { return c == '/'; }).base(),pathname.end())};
  }  

  std::string dirname(std::string source)
  {
    source.erase(std::find(source.rbegin(), source.rend(), '/').base(), source.end());
    return source;
  }    
  
  enum Status { SUCCESS, FAILED };

  using namespace std;
  using namespace dealii;

  template <int dim>
  class MySolver
  {
  public:
    MySolver( const std::string& );
    ~MySolver();

    void run ( string );

  protected:
    void make_grid();
    void setup_system();
    void setup_boundary_ids();
    void compute_contributions( double (&retval)[5] );
  
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
    
    int m_rank;
  };

  template <int dim>
  MySolver<dim>::MySolver ( const std::string& xmlfilename ) 
    : 
    m_ph(xmlfilename),
    mpi_communicator (MPI_COMM_WORLD),
    triangulation (mpi_communicator, typename Triangulation<dim>::MeshSmoothing(Triangulation<dim>::smoothing_on_refinement|Triangulation<dim>::smoothing_on_coarsening)),
    fe (gl_degree_fe),
    dof_handler (triangulation),
    pcout (cout, (Utilities::MPI::this_mpi_process(mpi_communicator) == 0))
  {
    m_omega = m_ph.Get_Physics("omega");
    m_gs = m_ph.Get_Physics("gs_1",0);

    m_xmin = m_ph.Get_Mesh("xrange",0);
    m_xmax = m_ph.Get_Mesh("xrange",1);
    m_ymin = m_ph.Get_Mesh("yrange",0);
    m_ymax = m_ph.Get_Mesh("yrange",1);
    m_zmin = m_ph.Get_Mesh("zrange",0);
    m_zmax = m_ph.Get_Mesh("zrange",1);
  
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
  void MySolver<dim>::compute_contributions ( double (&retval)[5] )
  {
    const QGauss<dim> quadrature_formula(fe.degree+1);
   
    FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_gradients|update_JxW_values|update_quadrature_points);

    const unsigned dofs_per_cell = fe.dofs_per_cell;
    const unsigned n_q_points    = quadrature_formula.size();

    vector<double> Psi(n_q_points);
    vector<Tensor<1, dim> > Psi_grad(n_q_points);    

    double JxW, Psiq, x, y;
    double omqx = m_omega[0]*m_omega[0];
    double omqy = m_omega[1]*m_omega[1];
    
    double tmp[] = {0,0,0,0,0};
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for ( ; cell!=endc; ++cell )
    {
      if( cell->is_locally_owned() )
      {
        fe_values.reinit (cell);
        fe_values.get_function_values( m_Psi, Psi );
        fe_values.get_function_gradients( m_Psi, Psi_grad );

        for( unsigned qp=0; qp<n_q_points; ++qp )
        {
          x = fe_values.quadrature_point(qp)[0];
          y = fe_values.quadrature_point(qp)[1];
          JxW = fe_values.JxW(qp);

          Psiq = Psi[qp]*Psi[qp];
          tmp[0] += JxW*Psi_grad[qp]*Psi_grad[qp];
#if POTENTIAL==1
          tmp[1] += JxW*omqx*x*x*Psiq;
#endif
#if POTENTIAL==2
          tmp[1] += JxW*m_omega[0]*x*Psiq;
#endif  
          tmp[2] += JxW*omqy*y*y*Psiq;
          tmp[3] += JxW*Psiq*Psiq;
          tmp[4] += JxW*Psiq;
        }
      }
    }
    tmp[3] *= m_gs;
    MPI_Allreduce( tmp, retval, 5, MPI_DOUBLE, MPI_SUM, mpi_communicator);
  }

  template <int dim>
  void MySolver<dim>::run( string filename )
  {
    load(filename);

    double retval[] = {0,0,0,0,0};
    compute_contributions( retval );
    pcout << std::scientific;  
    pcout << "mu = " << (retval[0] + retval[1] + retval[2] + retval[3])/retval[4] << endl;
#if POTENTIAL==1
    pcout << "virial = " << 2*(retval[0]-retval[1]-retval[2])+0.5*double(dim)*retval[3] << endl;
#endif
#if POTENTIAL==2
    pcout << "virial = " << (2*retval[0]-retval[1]-2*retval[2]+0.5*double(dim)*retval[3]) << endl;
#endif
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

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv );
  {
    std::string filename = argv[1];
    HelperPrograms::MySolver<DIMENSION> solver("params.xml");
    solver.run(filename);
  }  
return EXIT_SUCCESS;
}