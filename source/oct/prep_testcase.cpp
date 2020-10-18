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

/** 
 * Purpose: Real time propagation for the Gross-Pitaevskii equation (cartesian coordinates)
 * Method: Operator splitting - the linear step is solved with Crank-Nicolson 
 */

#include <deal.II/lac/generic_linear_algebra.h>

namespace LA
{
  using namespace dealii::LinearAlgebraPETSc;
}

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>

#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_precondition.h>
#include <deal.II/lac/sparsity_tools.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_out.h>

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
#include <iomanip>

#include "global.h"
#include "mpi.h"
#include "MyParameterHandler.h"
#include "my_table.h"
#include "MyComplexTools.h"

namespace realtime_propagation
{
  enum Status { SUCCESS, FAILED };

  using namespace std;
  using namespace dealii;

  template <int dim>
  class CPotential : public dealii::Function<dim> 
  {
    public:
      CPotential ( const std::vector<double>& omega, const double dom ) : Function<dim>(), m_omegaq(2)
      { 
        m_dom = dom;
        for( int i=0; i<dim; ++i )
          m_omegaq[i] = omega[i]*omega[i];
      }

      virtual double value ( const Point<dim> &p, const unsigned component = 0) const 
      {
        double t = this->get_time();
      return (1+2*sin(m_dom*t))*m_omegaq[0]*p[0]*p[0] + m_omegaq[1]*p[1]*p[1];
      }

      protected:
        std::vector<double> m_omegaq;
        double m_dom;
  };

  template <int dim>
  class MySolver
  {
  public:
    explicit MySolver( const std::string& );
    ~MySolver();

    void run ();

  protected:
    void make_grid();
    void setup_system();
    
    void Do_NL_Step( const double );
    void Do_Lin_Step();
    void output_results ( string );
    void load( string );
    void save( string );

    MyParameterHandler m_ph;
    MPI_Comm mpi_communicator;
    parallel::distributed::Triangulation<dim> triangulation;
    FESystem<dim> fe;
    DoFHandler<dim> dof_handler;
    IndexSet locally_owned_dofs;
    IndexSet locally_relevant_dofs;
    AffineConstraints<double> constraints;

    LA::MPI::SparseMatrix system_matrix;
    LA::MPI::Vector system_rhs;
    LA::MPI::Vector sol;
    LA::MPI::Vector m_Psi; 
    LA::MPI::Vector m_Psi_0; 

    ConditionalOStream pcout;
    ofstream m_computing_timer_log;
    TimerOutput m_computing_timer;    

    bool m_root;

    double m_gs;
    double m_t;
    double m_dt;
    double m_T;
    double m_domega;
    vector<double> m_omega;

    double m_xmin, m_xmax;
    double m_ymin, m_ymax;
    double m_zmin, m_zmax;

    int m_rank;

//    MyTable m_table;    
  };

  template <int dim>
  MySolver<dim>::MySolver ( const std::string& xmlfilename ) 
    : 
    m_ph(xmlfilename),
    mpi_communicator (MPI_COMM_WORLD),
    triangulation (mpi_communicator, typename Triangulation<dim>::MeshSmoothing(Triangulation<dim>::smoothing_on_refinement|Triangulation<dim>::smoothing_on_coarsening)),
    fe (FE_Q<dim>(gl_degree_fe), 2),
    dof_handler (triangulation),
    pcout (cout, (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)),
    m_computing_timer_log("benchmark.txt"),
    m_computing_timer(mpi_communicator, m_computing_timer_log, TimerOutput::summary, TimerOutput:: cpu_and_wall_times )
  {
    try
    {
      m_omega = m_ph.Get_Physics("omega");
      m_gs = m_ph.Get_Physics("gs_1",0);

      m_xmin = m_ph.Get_Mesh("xrange",0);
      m_xmax = m_ph.Get_Mesh("xrange",1);
      m_ymin = m_ph.Get_Mesh("yrange",0);
      m_ymax = m_ph.Get_Mesh("yrange",1);
      m_zmin = m_ph.Get_Mesh("zrange",0);
      m_zmax = m_ph.Get_Mesh("zrange",1);
    }
    catch( const std::string info )
    {
      std::cerr << info << endl;
      MPI_Abort( mpi_communicator, 0 );
    }    
    m_t = 0;
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
    TimerOutput::Scope timing_section(m_computing_timer, "");
    Point<dim,double> pt1;
    Point<dim,double> pt2;

    double min[] = {m_xmin, m_ymin, m_zmin};
    double max[] = {m_xmax, m_ymax, m_zmax};

    for( int i=0; i<dim; ++i )
    {
      pt1(i) = min[i];
      pt2(i) = max[i];
    }

    GridGenerator::hyper_rectangle(triangulation, pt2, pt1);
    
  }
  
  template <int dim>
  void MySolver<dim>::setup_system()
  {
    TimerOutput::Scope timing_section(m_computing_timer, "");

    dof_handler.distribute_dofs (fe);
   
    locally_owned_dofs = dof_handler.locally_owned_dofs ();
    DoFTools::extract_locally_relevant_dofs (dof_handler, locally_relevant_dofs);

    m_Psi.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    m_Psi_0.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    sol.reinit(locally_owned_dofs, mpi_communicator); // no ghosts
    system_rhs.reinit(locally_owned_dofs, mpi_communicator); // no ghosts

    constraints.clear ();
    constraints.reinit (locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints (dof_handler, constraints);
    VectorTools::interpolate_boundary_values (dof_handler, 0, ZeroFunction<dim>(2), constraints);
    constraints.close ();

    DynamicSparsityPattern dsp (locally_relevant_dofs);
    DoFTools::make_sparsity_pattern (dof_handler, dsp, constraints, false);
    SparsityTools::distribute_sparsity_pattern (dsp, dof_handler.n_locally_owned_dofs_per_processor(), mpi_communicator, locally_relevant_dofs);
    system_matrix.reinit (locally_owned_dofs, locally_owned_dofs, dsp, mpi_communicator);

    
  }  

  template <int dim>
  void MySolver<dim>::output_results ( string path ) 
  {
    TimerOutput::Scope timing_section(m_computing_timer, "");
    string filename;
    
    vector<std::string> solution_names;
    solution_names.push_back ("Re");
    solution_names.push_back ("Im");    

    DataOut<dim> data_out;
    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector ( m_Psi, solution_names );
    data_out.build_patches ();
    
    filename = path + "solution-" + to_string(m_t) + ".vtu";
    data_out.write_vtu_in_parallel ( filename.c_str(), mpi_communicator );

    
  }   

  template<int dim>
  void MySolver<dim>::save( string filename )
  {
    std::vector<const LA::MPI::Vector*> x_system (2);
    x_system[0] = &m_Psi_0;
    x_system[1] = &m_Psi;
    
    parallel::distributed::SolutionTransfer<dim,LA::MPI::Vector> solution_transfer(dof_handler);
    solution_transfer.prepare_for_serialization(x_system);
    triangulation.save( filename.c_str() );
  }

  template<int dim>
  void MySolver<dim>::load( string filename )
  {
    make_grid();
    triangulation.load( filename.c_str() );

    setup_system();

    parallel::distributed::SolutionTransfer<dim,LA::MPI::Vector> solution_transfer(dof_handler);
    solution_transfer.deserialize(system_rhs);

    constraints.distribute(system_rhs);
    m_Psi = system_rhs;
  }    

  template<int dim> 
  void MySolver<dim>::Do_NL_Step( double t )
  {
    TimerOutput::Scope timing_section(m_computing_timer, "");

    CPotential<dim> Potential (m_omega,m_domega);
    Potential.set_time(t);
    MyComplexTools::MPI::AssembleSystem_NL_Step<dim>( dof_handler, fe, constraints, m_Psi, Potential, m_dt, m_gs, system_matrix, system_rhs );

    sol=0;
    SolverControl solver_control (sol.size(), 1e-15);
    PETScWrappers::SolverCG solver (solver_control, mpi_communicator);
    PETScWrappers::PreconditionBlockJacobi::AdditionalData adata;
    PETScWrappers::PreconditionBlockJacobi preconditioner(system_matrix,adata);    
    solver.solve (system_matrix, sol, system_rhs, preconditioner);
    
    constraints.distribute (sol);
    m_Psi = sol;

       
  }
  
  template<int dim> 
  void MySolver<dim>::Do_Lin_Step()
  {
    TimerOutput::Scope timing_section(m_computing_timer, "");
    MyComplexTools::MPI::AssembleSystem_LIN_Step<dim>( dof_handler, fe, constraints, m_Psi, m_dt, system_matrix, system_rhs );

    sol=0;
    SolverControl solver_control (sol.size(), 1e-15);
    PETScWrappers::SolverGMRES solver (solver_control, mpi_communicator);
    PETScWrappers::PreconditionBlockJacobi::AdditionalData adata;
    PETScWrappers::PreconditionBlockJacobi preconditioner(system_matrix,adata);    
    solver.solve (system_matrix, sol, system_rhs, preconditioner);

/*
    SolverControl solver_control;
    PETScWrappers::SparseDirectMUMPS solver(solver_control, mpi_communicator);
    solver.set_symmetric_mode(false);
    solver.solve(system_matrix, sol, system_rhs);
*/
    constraints.distribute (sol);
    m_Psi = sol;    
       
  }
 
  template <int dim>
  void MySolver<dim>::run()
  {
    unsigned Nt = 101;
    m_T = 1;
    m_dt = m_T/double(Nt);
    m_domega = M_PI/m_T;

    pcout << "domega = " << m_domega << endl;

    load( "Cfinal.bin" );

    std::vector<double> p(dim);
    std::vector<double> pos(dim);
    std::vector<double> var(dim);

    double N = MyComplexTools::MPI::Particle_Number( mpi_communicator, dof_handler, fe, m_Psi );

    system_rhs *= 1.0/sqrt(N);
    constraints.distribute(system_rhs);
    m_Psi = system_rhs;
    m_Psi_0 = system_rhs;

    N = MyComplexTools::MPI::Particle_Number( mpi_communicator, dof_handler, fe, m_Psi );
    pcout << "N == " << N << endl;

    output_results("Cfinal");

    for( unsigned i=0; i<Nt; ++i )
    {
      Do_Lin_Step();
      Do_NL_Step(double(i)*m_dt);

      N = MyComplexTools::MPI::Particle_Number( mpi_communicator, dof_handler, fe, m_Psi );
      pcout << "N == " << N << endl;

      m_t += m_dt;

      output_results("Psi_");
    }
    save( "oct_0.bin" );
  }
} // end of namespace 

int main ( int argc, char *argv[] )
{
  using namespace dealii;
  deallog.depth_console (0);

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv );
  {
    realtime_propagation::MySolver<SPATIAL_DIM> solver("params_one.xml");
    solver.run();
  }  
return EXIT_SUCCESS;
}