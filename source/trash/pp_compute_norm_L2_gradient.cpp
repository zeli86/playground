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
 * Purpose: convert deal.ii native file format files into vtu files
 */

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/parameter_handler.h>

#include <deal.II/lac/generic_linear_algebra.h>

namespace LA
{
  using namespace dealii::LinearAlgebraPETSc;
}

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/compressed_simple_sparsity_pattern.h>

#include <deal.II/lac/petsc_parallel_sparse_matrix.h>
#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_precondition.h>

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
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/solution_transfer.h>

#include <fstream>
#include <iostream>
#include <stdlib.h>

#include "mpi.h"
#include "functions_cs.h"
#include "ParameterReader.h"
#include "my_table.h"
#include "ref_pt_list.h"
#include "tinyxml2.h"

namespace HelperPrograms
{
  enum Status { SUCCESS, FAILED };

  using namespace std;
  using namespace dealii;

  template <int dim>
  class MySolver
  {
  public:
    MySolver( ParameterHandler & );
    ~MySolver();

    void run ( string );

  protected:
    void make_grid();
    void setup_system();
    void setup_boundary_ids();
    void set_im_to_zero();

    double compute_visibility();
    double compute_difference();
   
    void load( string );
    void load_2( string );
    void save( string );

    ParameterHandler &m_ph;
    MPI_Comm mpi_communicator;
    parallel::distributed::Triangulation<dim> triangulation;
    FESystem<dim> fe;
    DoFHandler<dim> dof_handler;
    IndexSet locally_owned_dofs;
    IndexSet locally_relevant_dofs;
    ConstraintMatrix constraints;

    LA::MPI::Vector m_Psi;
    LA::MPI::Vector m_Psi_0;
    LA::MPI::Vector m_workspace;
    LA::MPI::Vector system_rhs;
    Vector<double> m_error_per_cell;

    ConditionalOStream pcout;
    ConditionalOStream pcerr;

    bool m_root;

    double m_xmin;
    double m_xmax;
    double m_ymin;
    double m_ymax;
    double m_zmin;
    double m_zmax;
    double m_res;
    double m_omega[3];
    double m_mu;
    double m_gs;
    double m_t;
    
    int m_rank;
    unsigned int m_NA;
    unsigned int m_NK;
    unsigned int m_global_refinement;
    unsigned int m_QN1[3];
  };

  template <int dim>
  MySolver<dim>::MySolver ( ParameterHandler &ph ) 
    : 
    m_ph(ph),
    mpi_communicator (MPI_COMM_WORLD),
    triangulation (mpi_communicator, typename Triangulation<dim>::MeshSmoothing(Triangulation<dim>::smoothing_on_refinement|Triangulation<dim>::smoothing_on_coarsening)),
    fe (FE_Q<dim>(1), 2),
    dof_handler (triangulation),
    pcout (cout, (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)),
    pcerr (cerr, (Utilities::MPI::this_mpi_process(mpi_communicator) == 0))
  {
    m_ph.enter_subsection("Mesh & geometry parameters");    
    m_global_refinement = m_ph.get_integer("global_refinements");
    m_xmin = m_ph.get_double("xMin");
    m_xmax = m_ph.get_double("xMax");
    m_ymin = m_ph.get_double("yMin");
    m_ymax = m_ph.get_double("yMax");
    m_zmin = m_ph.get_double("zMin");
    m_zmax = m_ph.get_double("zMax");
    m_ph.leave_subsection();
    
    m_ph.enter_subsection("Physical constants");
    m_QN1[0] = m_ph.get_integer("QN1_x");
    m_QN1[1] = m_ph.get_integer("QN1_y");
    m_QN1[2] = m_ph.get_integer("QN1_z");
  #if POTENTIAL==1
    m_omega[0] = m_ph.get_double("omega_x");
  #endif
  #if POTENTIAL==2
    m_omega[0] = m_ph.get_double("gf");
  #endif
    m_omega[1] = m_ph.get_double("omega_y");
    m_omega[2] = m_ph.get_double("omega_z");
    m_gs = m_ph.get_double("gs_11");
    m_ph.leave_subsection();    

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
    m_Psi_0.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    m_workspace.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);    
    system_rhs.reinit(locally_owned_dofs, mpi_communicator);
    m_error_per_cell.reinit(triangulation.n_active_cells());
    
    system_rhs = 0;

    vector<bool> mask (dof_handler.get_fe().n_components(), true );
    
    constraints.clear ();
    constraints.reinit (locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints (dof_handler, constraints);
    VectorTools::interpolate_boundary_values (dof_handler, 1, ZeroFunction<dim>(), constraints, ComponentMask(mask));
    
    if( m_QN1[2] > 0 )
      VectorTools::interpolate_boundary_values (dof_handler, 0, ZeroFunction<dim>(), constraints, ComponentMask(mask));
    //if( m_QN1[1] & 1 )
    //  VectorTools::interpolate_boundary_values (dof_handler, 0, ZeroFunction<dim>(), constraints);

    constraints.close ();
  }
  
  template<int dim>
  double MySolver<dim>::compute_visibility()
  {
    double retval[]={0,0,0}, tmp1[]={0,0,0}, JxW;
    
    const QGauss<dim>  quadrature_formula(fe.degree+1);
    FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_quadrature_points|update_JxW_values);

    const unsigned int n_q_points = quadrature_formula.size();
    vector<Vector<double>> Psi_0(n_q_points,Vector<double>(2));
    vector<Vector<double>> Psi(n_q_points,Vector<double>(2));

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
      if( cell->is_locally_owned() )
      {
        fe_values.reinit (cell);
        fe_values.get_function_values( m_Psi_0, Psi_0 );
        fe_values.get_function_values( m_Psi, Psi );
        for( unsigned qp = 0; qp < n_q_points; ++qp )
        {
         JxW = fe_values.JxW(qp)*fabs(fe_values.quadrature_point(qp)[1]);
         tmp1[0] += JxW*(Psi_0[qp][0]*Psi[qp][0]+Psi_0[qp][1]*Psi[qp][1]);
         tmp1[1] += JxW*(Psi_0[qp][0]*Psi[qp][1]-Psi_0[qp][1]*Psi[qp][0]);
	 tmp1[2] += JxW*(Psi_0[qp][0]*Psi_0[qp][0]+Psi_0[qp][1]*Psi_0[qp][1]);
        }
      }
    }
    MPI_Allreduce( tmp1, retval, 3, MPI_DOUBLE, MPI_SUM, mpi_communicator);
    return sqrt(retval[0]*retval[0] + retval[1]*retval[1])/retval[2];
  }  
  
  template <int dim>
  double MySolver<dim>::compute_difference()
  {
    double retval=0, tmp1=0, JxW=0;

    const QGauss<dim>  quadrature_formula(fe.degree+1);
    FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_quadrature_points|update_JxW_values);

    const unsigned int n_q_points = quadrature_formula.size();
    vector<Vector<double>> Psi_0(n_q_points,Vector<double>(2));
    vector<Vector<double>> Psi(n_q_points,Vector<double>(2));

    
    const double re = cos(-m_mu*m_t);
    const double im = sin(-m_mu*m_t);
    
    double Psi_neu[] = {0,0};
    
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
      if( cell->is_locally_owned() )
      {
        fe_values.reinit (cell);
        fe_values.get_function_values( m_Psi_0, Psi_0 );
        fe_values.get_function_values( m_Psi, Psi );
        for( unsigned qp = 0; qp < n_q_points; ++qp )
        {
         JxW = fe_values.JxW(qp)*fabs(fe_values.quadrature_point(qp)[1]);
	 
	 Psi_neu[0] = (Psi_0[qp][0]*re - Psi_0[qp][1]*im) - Psi[qp][0];
	 Psi_neu[1] = (Psi_0[qp][0]*im + Psi_0[qp][1]*re) - Psi[qp][1];
	 
         tmp1 += JxW*(Psi_neu[0]*Psi_neu[0]+Psi_neu[1]*Psi_neu[1]);
        }
      }
    }
    
    MPI_Allreduce( &tmp1, &retval, 1, MPI_DOUBLE, MPI_SUM, mpi_communicator);
    return retval;
  }
  
  template <int dim>
  void MySolver<dim>::run( string filename )
  {
    tinyxml2::XMLDocument doc;
    doc.LoadFile( "info.xml" );
    tinyxml2::XMLElement* xml_element_mu = doc.FirstChildElement()->FirstChildElement( "MU" );
    xml_element_mu->QueryDoubleText( &m_mu );
    
    m_t=0; 
    double vis=1;
    double d=0;
    m_res=0;
    
    if( filename == "final.bin" )
    {
      load(filename);
    }
    else
    { 
      load_2( filename );
      string tmpstr=filename;
      std::size_t pos = tmpstr.find("-");
      tmpstr.erase(0,pos+1);
      tmpstr.erase(tmpstr.length()-4,4);
      m_t = stod(tmpstr);
    }
    
    /*
    double min_cell_diameter = GridTools::minimal_cell_diameter(triangulation);
    double max_cell_diameter = GridTools::maximal_cell_diameter(triangulation);
    pcerr << "min_cell_diameter = " << min_cell_diameter << endl;
    pcerr << "max_cell_diameter = " << max_cell_diameter << endl;
    */
    
    //compute_residual();
    if( filename != "final.bin" )
    { 
      vis = compute_visibility();
      d = compute_difference();
    }

    //pcout << "t = " << t << endl;
    //pcout << "m_mu = " << m_mu << endl;
    //pcout << "m_gs = " << m_gs << endl;
    //pcout << "m_res = " << m_res << endl;
    pcout << std::scientific;  
    pcout << m_t << "\t" << std::setprecision(5) << m_res << "\t" << vis << "\t" << d <<endl;
    //pcerr << m_t << "\t" << std::setprecision(5) << m_res << "\t" << vis << "\t" << d <<endl;
  }

  template<int dim>
  void MySolver<dim>::load( string filename )
  {
    make_grid();
    triangulation.load( filename.c_str() );
    setup_boundary_ids();
    
    setup_system();

    parallel::distributed::SolutionTransfer<dim,LA::MPI::Vector> solution_transfer(dof_handler);
    solution_transfer.deserialize(m_Psi);
    
    set_im_to_zero();
  }
  
  template<int dim>
  void MySolver<dim>::load_2( string filename )
  {
    make_grid();
    triangulation.load( filename.c_str() );
    setup_boundary_ids();
    
    setup_system();

    std::vector<LA::MPI::Vector*> x_system (2);
    x_system[0] = &m_Psi;
    x_system[1] = &m_Psi_0;

    parallel::distributed::SolutionTransfer<dim,LA::MPI::Vector> solution_transfer(dof_handler);
    solution_transfer.deserialize(x_system);
    
    LA::MPI::Vector tmp;
    tmp.reinit (locally_owned_dofs, mpi_communicator);
    
    tmp=m_Psi;
    constraints.distribute(tmp);
    m_Psi=tmp;
    
    tmp=m_Psi_0;
    constraints.distribute(tmp);
    m_Psi_0=tmp;  
  }

  template<int dim>
  void MySolver<dim>::setup_boundary_ids()
  {
    typename parallel::distributed::Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active(), endc = triangulation.end();
    for(; cell!=endc; ++cell)
    {
      for (unsigned int f=0; f < GeometryInfo<dim>::faces_per_cell; ++f)
      {
	const Point<dim> face_center = cell->face(f)->center();
	if (cell->face(f)->at_boundary() && !(face_center[1]==0) )    
	{
	  cell->face(f)->set_all_boundary_indicators(1);
	}
      }
    }
  }
  
  template<int dim>
  void MySolver<dim>::set_im_to_zero ()
  {
    //printf( "dof_handler.get_fe().n_components() = %d\n", dof_handler.get_fe().n_components() );
    
    vector<bool> selected_dofs ( dof_handler.n_locally_owned_dofs() );
    vector<bool> mask (dof_handler.get_fe().n_components(), false);
    
    mask[1] = true;
    DoFTools::extract_dofs (dof_handler, ComponentMask(mask), selected_dofs);

    vector<types::global_dof_index> indices(dof_handler.n_locally_owned_dofs());
    locally_owned_dofs.fill_index_vector ( indices );
    

    for( unsigned int i=0; i<indices.size(); i++ )
    {
      if( selected_dofs[i] == true )
      {
        m_Psi[indices[i]] = 0.0;
      }
    }
    m_Psi.compress(VectorOperation::insert);

    LA::MPI::Vector tmp;
    tmp.reinit (locally_owned_dofs, mpi_communicator);
    
    tmp=m_Psi;
    constraints.distribute(tmp);
    m_Psi=tmp;
  }
} // end of namespace 

int main ( int argc, char *argv[] )
{
  int myrank;

  using namespace dealii;
  deallog.depth_console (0);

  ParameterHandler  prm;
  ParameterReader   param(prm);

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv );
  {
    MPI_Comm_rank( MPI_COMM_WORLD, &myrank );
    if( !param.read_parameters("params.prm") )
    {
      MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
      return EXIT_FAILURE;
    }

    std::string filename = argv[1];
    HelperPrograms::MySolver<DIMENSION> solver(prm);
    solver.run(filename);
  }  
return EXIT_SUCCESS;
}
