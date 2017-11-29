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

#include <deal.II/lac/generic_linear_algebra.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>

#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/petsc_parallel_sparse_matrix.h>
#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/petsc_precondition.h>
#include <deal.II/lac/petsc_solver.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_in.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
//#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_q.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/function_parser.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/grid_refinement.h>

#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <iomanip>

#include "global.h"
#include "mpi.h"
#include "muParser.h"
#include "MyParameterHandler.h"
#include "my_table.h"
#include "MyComplexTools.h"

namespace LA
{
  using namespace dealii::LinearAlgebraPETSc;
}

namespace edu
{
  using namespace std;
  using namespace dealii;

  template <int dim>
  class MySolver
  {
  public:
    MySolver( MyParameterHandler& );
    ~MySolver();

    void run ();

  protected:
    void refine_mesh();
    void setup_system();
    
    void Do_NL_Step();
    void Do_Lin_Step( const double );
    void output_results ( string );

    MyParameterHandler& m_ph;
    MPI_Comm mpi_communicator;
    parallel::distributed::Triangulation<dim> triangulation;
    FESystem<dim> fe;
    DoFHandler<dim> dof_handler;
    bool m_root;
    ConditionalOStream pcout;

    IndexSet locally_owned_dofs;
    IndexSet locally_relevant_dofs;
    ConstraintMatrix constraints;

    LA::MPI::SparseMatrix system_matrix;
    LA::MPI::Vector system_rhs; // no ghost cells
    LA::MPI::Vector sol; // no ghost cells
    LA::MPI::Vector m_Psi; 

    double m_gs;
    double m_t;
    double m_dt;
    double m_dth; // 0.5*m_dt
 
    int m_rank;
    unsigned m_NA;
    unsigned m_NK;

    std::string m_mesh_filename;
    std::vector<std::string> m_wf_str;
  };

  template <int dim>
  MySolver<dim>::MySolver ( MyParameterHandler& ph ) 
    : 
    m_ph(ph),
    mpi_communicator (MPI_COMM_WORLD),
    triangulation (mpi_communicator, typename Triangulation<dim>::MeshSmoothing(Triangulation<dim>::smoothing_on_refinement|Triangulation<dim>::smoothing_on_coarsening)),
    fe (FE_Q<dim>(gl_degree_fe), 2),
    dof_handler (triangulation),
    m_root(Utilities::MPI::this_mpi_process(mpi_communicator) == 0),
    pcout (cout, m_root )
  {
    try
    {
      m_gs = m_ph.Get_Physics("gs_1",0);

      m_NA = unsigned(m_ph.Get_Algorithm("NA",0));
      m_NK = unsigned(m_ph.Get_Algorithm("NK",0));
      m_dt = m_ph.Get_Algorithm("dt",0);    

      m_mesh_filename = m_ph.Get_String( "FILENAME_MESH", dim-2 );
      m_wf_str = m_ph.Get_AllStrings("WAVEFUNCTION_" + to_string(dim) + "D");
      
    }
    catch( const std::string info )
    {
      cerr << info << endl;
      MPI_Abort( mpi_communicator, 0 );
    }    
    m_dth = 0.5*m_dt;
    m_t = 0;
    MPI_Comm_rank(mpi_communicator, &m_rank);
  }

  template <int dim>
  MySolver<dim>::~MySolver ()
  {
    dof_handler.clear ();
  }
  
  template <int dim>
  void MySolver<dim>::setup_system()
  {
    dof_handler.distribute_dofs (fe);
   
    locally_owned_dofs = dof_handler.locally_owned_dofs ();
    DoFTools::extract_locally_relevant_dofs (dof_handler, locally_relevant_dofs);

    m_Psi.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    sol.reinit(locally_owned_dofs, mpi_communicator); // no ghost cells
    system_rhs.reinit(locally_owned_dofs, mpi_communicator); // no ghost cells

    constraints.clear ();
    constraints.reinit (locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints (dof_handler, constraints);
    VectorTools::interpolate_boundary_values (dof_handler, 0, ZeroFunction<dim>(2), constraints);
    constraints.close ();

    DynamicSparsityPattern dsp (locally_relevant_dofs);
    DoFTools::make_sparsity_pattern (dof_handler, dsp, constraints, false);
    SparsityTools::distribute_sparsity_pattern (dsp, dof_handler.n_locally_owned_dofs_per_processor(), mpi_communicator, locally_relevant_dofs);
    system_matrix.reinit (locally_owned_dofs, locally_owned_dofs, dsp, mpi_communicator);

    pcout << "total number of dofs = " << dof_handler.n_dofs () << endl;
  }  

  template <int dim>
  void MySolver<dim>::output_results ( string path ) 
  {
    string filename;
    
    vector<std::string> solution_names;
    solution_names.push_back ("Re_Psi");
    solution_names.push_back ("Im_Psi");    

    DataOut<dim> data_out;
    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector ( m_Psi, solution_names );
    data_out.build_patches (gl_degree_fe);
    
    filename = path + "solution-" + to_string(m_t) + ".vtu";
    data_out.write_vtu_in_parallel ( filename.c_str(), mpi_communicator );
  }   
  
  template<int dim> 
  void MySolver<dim>::Do_NL_Step()
  {
    pcout << "Computing nonlinear step, t = " << m_dt << endl;

    MyComplexTools::MPI::AssembleSystem_NL_Step<dim>( dof_handler, fe, constraints, m_Psi, m_gs*m_dt, system_matrix, system_rhs );

    SolverControl solver_control (sol.size(), 1e-15);
    PETScWrappers::SolverBicgstab solver (solver_control, mpi_communicator);

    PETScWrappers::PreconditionParaSails::AdditionalData adata;
    PETScWrappers::PreconditionParaSails preconditioner(system_matrix,adata);    

    solver.solve (system_matrix, sol, system_rhs, preconditioner);
    
    constraints.distribute (sol);
    m_Psi = sol;
  }
  
  template<int dim> 
  void MySolver<dim>::Do_Lin_Step( const double dt )
  {
    pcout << "Computing linear step (" << dt << "), t = " << m_t << endl;

    FunctionParser<dim> Potential;
    Potential.initialize( FunctionParser<dim>::default_variable_names(), m_ph.Get_String( "POTENTIAL", dim-2 ), m_ph.Get_Constants_Map());

    MyComplexTools::MPI::AssembleSystem_LIN_Step<dim>( dof_handler, fe, constraints, m_Psi, Potential, dt, system_matrix, system_rhs );

    SolverControl solver_control (sol.size(), 1e-15);
    PETScWrappers::SolverBicgstab solver (solver_control, mpi_communicator);
    
    PETScWrappers::PreconditionBlockJacobi preconditioner(system_matrix);    

    //PETScWrappers::PreconditionParaSails::AdditionalData adata;
    //PETScWrappers::PreconditionParaSails preconditioner(system_matrix,adata);    
    solver.solve (system_matrix, sol, system_rhs, preconditioner);
    
    constraints.distribute (sol);
    m_Psi = sol;    

    m_t += dt;
  }

  template <int dim>
  void MySolver<dim>::refine_mesh()
  {

    for( unsigned step=0; step<2; step++ )
    {
      typename parallel::distributed::Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active(), endc = triangulation.end();
      for( ; cell!=endc; ++cell )
      {
        if( cell->measure() > 0.1 ) 
        {
          cell->set_refine_flag();
        }
      }
      triangulation.execute_coarsening_and_refinement ();
      pcout << "global active cells = " << triangulation.n_global_active_cells() << endl;      
    }        

    for( unsigned step=0; step<1; step++ )
    {
      typename parallel::distributed::Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active(), endc = triangulation.end();
      for( ; cell!=endc; ++cell )
      {
        for (unsigned v=0; v < GeometryInfo<dim>::vertices_per_cell; ++v )
        {
          Point<dim> p = cell->vertex(v);
          
          if( p(dim-1) > 0 )
          {
            cell->set_refine_flag ((dim == 3) ? RefinementCase<dim>::cut_z : RefinementCase<dim>::cut_y );
            break;
          }
        }
      }
      triangulation.execute_coarsening_and_refinement ();
      pcout << "global active cells = " << triangulation.n_global_active_cells() << endl;      
    }    
  }
 
  template <int dim>
  void MySolver<dim>::run()
  {
    double T, N, W;

    GridIn<dim> grid_in;
    grid_in.attach_triangulation(triangulation);

    std::ifstream input_file(m_mesh_filename.c_str());
    Assert (input_file, ExcFileNotOpen(m_mesh_filename.c_str()));

    grid_in.read_msh(input_file);
    
    //refine_mesh();
    setup_system();
    
    FunctionParser<dim> wf0(2);
    wf0.initialize( FunctionParser<dim>::default_variable_names(), m_wf_str, m_ph.Get_Constants_Map());
    
    VectorTools::interpolate (dof_handler, wf0, system_rhs );
    constraints.distribute (system_rhs);
    m_Psi = system_rhs;

    output_results("");
   
    std::vector<double> pos(dim);
    std::vector<double> var(dim);

    N = MyComplexTools::MPI::Particle_Number( mpi_communicator, dof_handler, fe, m_Psi );
    pcout << "N == " << N << endl;

    MyComplexTools::MPI::Expectation_value_position( mpi_communicator, dof_handler, fe, m_Psi, pos );
    MyComplexTools::MPI::Expectation_value_width( mpi_communicator, dof_handler, fe, m_Psi, pos, var );

    pcout << "N == " << MyComplexTools::MPI::Particle_Number( mpi_communicator, dof_handler, fe, m_Psi ) << endl;
    pcout << "pos == " << pos[0]/N << ", " << pos[1]/N << ", " << pos[2]/N << endl;
    pcout << "var == " << var[0]/N << ", " << var[1]/N << ", " << var[2]/N << endl;

    for( unsigned i=1; i<=m_NA; i++ )
    {
      Do_Lin_Step( m_dth );
      for( unsigned j=2; j<=m_NK; j++ )
      {
	      Do_NL_Step();
        Do_Lin_Step( m_dt );
      }
      Do_NL_Step();
      Do_Lin_Step( m_dth );

      MyComplexTools::MPI::Expectation_value_position( mpi_communicator, dof_handler, fe, m_Psi, pos );
      MyComplexTools::MPI::Expectation_value_width( mpi_communicator, dof_handler, fe, m_Psi, pos, var );

      pcout << "N == " << MyComplexTools::MPI::Particle_Number( mpi_communicator, dof_handler, fe, m_Psi ) << endl;
      pcout << "pos == " << pos[0]/N << ", " << pos[1]/N << ", " << pos[2]/N << endl;
      pcout << "var == " << var[0]/N << ", " << var[1]/N << ", " << var[2]/N << endl;

      output_results("");
    }
  }
} // end of namespace 

int main ( int argc, char *argv[] )
{
  using namespace dealii;
  deallog.depth_console (0);

  if( argc != 2 )
  {
    cout << "missing xml file" << endl;
    return EXIT_FAILURE;
  }

  MyParameterHandler params(argv[1]);
  int dim=0;

  try
  {
    dim = int(params.Get_Mesh("DIM",0));
  }
  catch (mu::Parser::exception_type &e)
  {
    cout << "Message:  " << e.GetMsg() << "\n";
    cout << "Formula:  " << e.GetExpr() << "\n";
    cout << "Token:    " << e.GetToken() << "\n";
    cout << "Position: " << e.GetPos() << "\n";
    cout << "Errc:     " << e.GetCode() << "\n";
  }

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv );
  {
    switch(dim)
    {
      case 2: { edu::MySolver<2> solver(params);
                solver.run();
                break; }
      case 3: { edu::MySolver<3> solver(params);
                solver.run();
                break; }
      default: cout << "You have found a new dimension!" << endl;
    }
  }
return EXIT_SUCCESS;
}