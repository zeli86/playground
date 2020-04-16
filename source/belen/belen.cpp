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


#include <deal.II/lac/generic_linear_algebra.h>

namespace LA
{
  using namespace dealii::LinearAlgebraPETSc;
}

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>

#include <deal.II/lac/petsc_parallel_sparse_matrix.h>
#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_precondition.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/derivative_approximation.h>
#include <deal.II/numerics/fe_field_function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/function_parser.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/solution_transfer.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <locale>
#include <limits>
#include <cmath>
#include <array>
#include <unistd.h>

#include "global.h"
#include "mpi.h"
#include "functions_cs.h"
#include "my_table.h"
#include "MyParameterHandler.h"

namespace BreedSolver
{
  using namespace std;
  using namespace dealii;

  class MySolver 
  {
  public:
    MySolver( const std::string );
    virtual ~MySolver();

    void run ();

  protected:
    
    void make_grid_custom();
    void setup_system();
    void assemble_system_and_rhs();
    void solve();
    void save( const string& );
    void output_vec ( const string&, PETScWrappers::MPI::Vector& );

    double m_xmin, m_xmax;
    double m_ymin, m_ymax;

    MPI_Comm mpi_communicator;
    MyParameterHandler m_ph;
    parallel::distributed::Triangulation<2> triangulation;
    FE_Q<2> fe;
    DoFHandler<2> dof_handler;
    IndexSet locally_owned_dofs;
    IndexSet locally_relevant_dofs;
    ConstraintMatrix constraints;
    FunctionParser<2> m_mass_dist;
    FunctionParser<2> m_mass_dist_z;

    LA::MPI::SparseMatrix m_system_matrix;
    LA::MPI::Vector m_system_rhs;
    LA::MPI::Vector m_solution;
    LA::MPI::Vector m_workspace_1;
    LA::MPI::Vector m_workspace_2;

    ConditionalOStream pcout;

    bool m_root;
    int m_rank;

    vector<double> m_x;
    vector<double> m_y;
    gsl_interp_accel * m_acc;
    gsl_spline * m_spline;
    
    //vector<double> m_omega;

    MyTable m_table;
  };

/*************************************************************************************************/
/**
 * Constructor
 */

  MySolver::MySolver ( const std::string xmlfilename ) 
    : 
    mpi_communicator(MPI_COMM_WORLD), 
    m_ph(xmlfilename),
    pcout (cout, (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)),
    triangulation (mpi_communicator, Triangulation<2>::MeshSmoothing(Triangulation<2>::limit_level_difference_at_vertices|Triangulation<2>::eliminate_refined_inner_islands|Triangulation<2>::smoothing_on_refinement|Triangulation<2>::smoothing_on_coarsening)),
    fe (gl_degree_fe),
    dof_handler (triangulation),
    m_mass_dist(1),
    m_mass_dist_z(1)
  {
    std::map <std::string, double> constants;
    std::string variables = "x,y";
    std::string expression, expression2;

    try
    {
      m_xmin = m_ph.Get_Mesh("xrange",0);
      m_xmax = m_ph.Get_Mesh("xrange",1);
      m_ymin = m_ph.Get_Mesh("yrange",0);
      m_ymax = m_ph.Get_Mesh("yrange",1);
      //m_global_refinement = unsigned(m_ph.Get_Mesh("global_refinements",0));
      //m_NA = int(m_ph.Get_Algorithm("NA",0));

      constants["z_0"] = m_ph.Get_Physics("z_0",0);
      constants["M_L"] = m_ph.Get_Physics("M_L",0);
      constants["C"] = m_ph.Get_Physics("C",0);
      constants["mu_0"] = m_ph.Get_Physics("mu_0",0);
      constants["m"] = m_ph.Get_Physics("m",0);
      m_x = m_ph.Get_Physics("xdata");
      m_y = m_ph.Get_Physics("ydata");
      expression = m_ph.Get_Parameter("mass_dist");
      expression2 = m_ph.Get_Parameter("mass_dist_z");


    }
    catch( const std::string info )
    {
      std::cerr << info << endl;
      MPI_Abort( mpi_communicator, 0 );
    }
    m_root = (Utilities::MPI::this_mpi_process(mpi_communicator) == 0);
    MPI_Comm_rank(mpi_communicator, &m_rank);

    m_mass_dist.initialize(variables, expression, constants);
    m_mass_dist_z.initialize(variables, expression2, constants);

    // init spline stuff
    int N = m_x.size();
    m_acc = gsl_interp_accel_alloc ();
    m_spline = gsl_spline_alloc (gsl_interp_cspline, N);

    gsl_spline_init(m_spline, m_x.data(), m_y.data(), N);
  }

  MySolver::~MySolver ()
  {
    // free spline stuff
    gsl_spline_free (m_spline);
    gsl_interp_accel_free (m_acc);
    dof_handler.clear ();
  }

  void MySolver::assemble_system_and_rhs ()
  {
    const QGauss<2> quadrature_formula(fe.degree+1);
   
    double last_m_x = m_x[m_x.size()-1];

    double xi, yi;

    m_system_matrix=0;
    m_system_rhs=0;

    FEValues<2> fe_values (fe, quadrature_formula, update_values|update_gradients|update_JxW_values|update_quadrature_points);

    const unsigned dofs_per_cell = fe.dofs_per_cell;
    const unsigned n_q_points = quadrature_formula.size();

    Vector<double> cell_rhs (dofs_per_cell);
    FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);

    vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    Point<2> p;

    DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for ( ; cell!=endc; ++cell )
    {
      if( cell->is_locally_owned() )
      {
        cell_matrix=0;
        cell_rhs=0;
        fe_values.reinit (cell);

        for ( unsigned qp=0; qp<n_q_points; ++qp )
        {
          double JxW = fe_values.JxW(qp)*fe_values.quadrature_point(qp)[0];
          p = fe_values.quadrature_point(qp);

          xi = p[0];
          if( xi > last_m_x )
            yi = 0;
          else
            yi = gsl_spline_eval (m_spline, xi, m_acc);
          for ( unsigned i=0; i<dofs_per_cell; ++i )
          {
            cell_rhs(i) += JxW*((m_mass_dist.value(p)+m_mass_dist_z.value(p)*yi)*fe_values.shape_value(i,qp));

            for ( unsigned j=0; j<dofs_per_cell; ++j )
            {
              cell_matrix(i,j) += -JxW*(fe_values.shape_grad(i,qp)*fe_values.shape_grad(j,qp)); // grad h1 grad h2
            }
          }
        }
        cell->get_dof_indices (local_dof_indices);
        constraints.distribute_local_to_global(cell_matrix, cell_rhs, local_dof_indices, m_system_matrix, m_system_rhs);
      }
    }
    m_system_matrix.compress(VectorOperation::add);
    m_system_rhs.compress(VectorOperation::add);   
  }
  
  void MySolver::solve ()
  {
    pcout << "Solving..." << endl;
    SolverControl solver_control;
    
    PETScWrappers::SparseDirectMUMPS solver(solver_control, mpi_communicator);
    solver.solve(m_system_matrix, m_solution, m_system_rhs);
    constraints.distribute (m_solution);    
  }
  
  void MySolver::make_grid_custom ()
  {
    Point<2,double> pt1(0,m_ymin); 
    Point<2,double> pt2(m_xmax,m_ymax);
  
    GridGenerator::hyper_rectangle(triangulation, pt2, pt1);
    triangulation.refine_global(4);

    const double fakx = 0.8;
    const double faky = 0.8;
    double Lx = m_xmax-m_xmin;
    double Ly = m_ymax;

    for( unsigned step=0; step<4; step++ )
    {
      Lx *= fakx;
      Ly *= faky;
    
      parallel::distributed::Triangulation<2>::active_cell_iterator cell = triangulation.begin_active(), endc = triangulation.end();
      for( ; cell!=endc; ++cell )
        for ( unsigned v=0; v < GeometryInfo<2>::vertices_per_cell; v++ )
        {
          Point<2> p = cell->vertex(v);
          if( fabs(p[0]) < Lx && fabs(p[1]) < Ly )
          {
            cell->set_refine_flag ();
            break;
          }
        }
      triangulation.execute_coarsening_and_refinement ();
    }

    parallel::distributed::Triangulation<2>::active_cell_iterator cell = triangulation.begin_active(), endc = triangulation.end();
    for(; cell!=endc; ++cell)
    {
      for ( unsigned f=0; f < GeometryInfo<2>::faces_per_cell; ++f)
      {
        const Point<2> face_center = cell->face(f)->center();
        if (cell->face(f)->at_boundary() && face_center[0]==0 )  
        {  
          cell->face(f)->set_all_boundary_ids(1);
        }
      }
    }
/*
    GridOutFlags::Msh opt(true, true);
    string filename = "grid-" + Utilities::int_to_string(triangulation.locally_owned_subdomain(), 4);
    ofstream output ((filename + ".msh").c_str());
    GridOut grid_out;
    grid_out.set_flags(opt);
    grid_out.write_msh (triangulation,output);
*/
  }

  void MySolver::setup_system()
  {
    dof_handler.distribute_dofs (fe);
      
    locally_owned_dofs = dof_handler.locally_owned_dofs ();
    DoFTools::extract_locally_relevant_dofs (dof_handler, locally_relevant_dofs);
  
    m_solution.reinit (locally_owned_dofs, mpi_communicator); 
    m_system_rhs.reinit(locally_owned_dofs, mpi_communicator);
    m_workspace_1.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    m_workspace_2.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);

    //cout << "(" << myrank << ") locally_owned_dofs = " << m_Psi_1.local_size()  << endl;
    
    constraints.clear ();
    constraints.reinit (locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints (dof_handler, constraints);
      VectorTools::interpolate_boundary_values (dof_handler, 0, ZeroFunction<2>(), constraints);
    constraints.close ();

    DynamicSparsityPattern dsp (locally_relevant_dofs);
    DoFTools::make_sparsity_pattern (dof_handler, dsp, constraints, false);
    SparsityTools::distribute_sparsity_pattern (dsp, dof_handler.n_locally_owned_dofs_per_processor(), mpi_communicator, locally_relevant_dofs);
    m_system_matrix.reinit (locally_owned_dofs, locally_owned_dofs, dsp, mpi_communicator);    
  }

  void MySolver::output_vec ( const std::string& filename, PETScWrappers::MPI::Vector& vec )
  {
    constraints.distribute(vec);
    m_workspace_1=vec;

    VectorTools::interpolate (dof_handler, m_mass_dist, m_system_rhs );
    constraints.distribute(m_system_rhs);
    m_workspace_2 = m_system_rhs;

    Vector<float> subdomain (triangulation.n_active_cells());
    for ( unsigned i=0; i<subdomain.size(); ++i )
      subdomain(i) = triangulation.locally_owned_subdomain();

    DataOut<2> data_out;
    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (m_workspace_1, "vec");
    data_out.add_data_vector (m_workspace_2, "mass_dist");
    data_out.add_data_vector (subdomain, "subdomain");
    data_out.build_patches ();
    data_out.write_vtu_in_parallel (filename.c_str(), mpi_communicator);
  }

  void MySolver::save ( const std::string& filename )
  {
    constraints.distribute(m_solution);
    m_workspace_1=m_solution;

    parallel::distributed::SolutionTransfer<2,LA::MPI::Vector> solution_transfer(dof_handler);
    solution_transfer.prepare_serialization(m_workspace_1);

    triangulation.save( filename.c_str() );
  }

  void MySolver::run()
  {
    make_grid_custom();
    setup_system();
    
    assemble_system_and_rhs();
    solve();
  
    output_vec( "grav_pot.vtu", m_solution );
    save( "grav_pot.bin" );
  }
} // end of namespace 

int main ( int argc, char *argv[] )
{
  using namespace dealii;
  deallog.depth_console (0);

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv );
  {
    BreedSolver::MySolver solver("params.xml");
    solver.run ();
  }
return EXIT_SUCCESS;
}
