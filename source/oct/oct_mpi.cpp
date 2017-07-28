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
 * DANGER Will Robinson!
 */

#include <deal.II/lac/generic_linear_algebra.h>

namespace LA
{
  using namespace dealii::LinearAlgebraPETSc;
}

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>

#include <deal.II/lac/petsc_parallel_sparse_matrix.h>
#include <deal.II/lac/petsc_parallel_vector.h>
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
#include <deal.II/base/function_cspline.h>
#include <deal.II/base/timer.h>

#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/solution_transfer.h>

#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <iomanip>
#include <array>

#include "global.h"
#include "mpi.h"
#include "MyParameterHandler.h"
#include "my_table.h"
#include "MyComplexTools.h"
#include "muParser.h"


namespace realtime_propagation
{
  enum Status { SUCCESS, FAILED };

  using namespace std;
  using namespace dealii;
  
  #include "potential.h"    

  template <int dim, int no_time_steps>
  class MySolver
  {
  public:
    MySolver( const std::string&, const double );
    ~MySolver();

    void run ();

  protected:
    CPotential<dim,no_time_steps> m_potential;

    void rt_propagtion_forward ( const int ); 
    void rt_propagtion_backward ( const int ); 
    void compute_correction ( const int ); 
    
    void make_grid();
    void setup_system();

    void assemble_system( const int );
    
    void solve();
    void solve_cg();
    void solve_eq1();
    void output_results ( string );
    void output_vec ( string, LA::MPI::Vector& );
    void load( string );
    void save( string );    

    MyParameterHandler m_ph;
    MPI_Comm mpi_communicator;
    parallel::distributed::Triangulation<dim> triangulation;
    FESystem<dim> fe;
    DoFHandler<dim> dof_handler;
    IndexSet locally_owned_dofs;
    IndexSet locally_relevant_dofs;
    ConstraintMatrix constraints;

    LA::MPI::SparseMatrix system_matrix;
    LA::MPI::Vector system_rhs;
    LA::MPI::Vector m_sol;
    LA::MPI::Vector m_Psi; // Psi(t)
    LA::MPI::Vector m_Psi_d; // desired state
    LA::MPI::Vector m_workspace;
    LA::MPI::Vector m_workspace_2;
    LA::MPI::Vector m_workspace_3;
    LA::MPI::Vector m_workspace_ng;
    LA::MPI::Vector m_p_D;
    LA::MPI::Vector m_p_L;
    Vector<double> m_error_per_cell;
    
    array<LA::MPI::Vector,no_time_steps> m_all_Psi;
    array<LA::MPI::Vector,no_time_steps> m_all_p;

    ConditionalOStream pcout;
    ofstream m_computing_timer_log;
    TimerOutput m_computing_timer;    

    bool m_root;

    double m_gs;
    double m_dt;
    double m_T;
    double m_dth; // 0.5*m_dt
    double m_res;
    double m_N;
    double m_N_pT;
    vector<double> m_norm_grad;
    vector<double> m_omega;

    double m_xmin, m_xmax;
    double m_ymin, m_ymax;
    double m_zmin, m_zmax;

    int m_rank;
    
    MyTable m_table;    
  };

  template <int dim, int no_time_steps>
  MySolver<dim,no_time_steps>::MySolver ( const std::string& xmlfilename, const double T ) 
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
    m_T = T;   
    m_dt = T/double(no_time_steps-1);
    m_dth = 0.5*m_dt;

    m_root = (Utilities::MPI::this_mpi_process(mpi_communicator) == 0);
    MPI_Comm_rank( mpi_communicator, &m_rank );    
  }

  template <int dim, int no_time_steps>
  MySolver<dim,no_time_steps>::~MySolver ()
  {
    dof_handler.clear ();
  }
  
  #include "eq1_mpi.h"
  #include "eq2_mpi.h"
  #include "eq3_mpi.h"

  template <int dim, int no_time_steps>
  void MySolver<dim,no_time_steps>::make_grid ()
  {
    m_computing_timer.enter_section(__func__);
    Point<dim,double> pt1;
    Point<dim,double> pt2;

    double min[] = {m_xmin, m_ymin, m_zmin};
    double max[] = {m_xmax, m_ymax, m_zmax};

    for( int i=0; i<dim; i++ )
    {
      pt1(i) = min[i];
      pt2(i) = max[i];
    }

    GridGenerator::hyper_rectangle(triangulation, pt2, pt1);
    m_computing_timer.exit_section();
  }
  
  template <int dim, int no_time_steps>
  void MySolver<dim,no_time_steps>::setup_system()
  {
    m_computing_timer.enter_section(__func__);
    dof_handler.distribute_dofs (fe);
    
    //DoFRenumbering::component_wise (dof_handler);

    locally_owned_dofs = dof_handler.locally_owned_dofs ();
    DoFTools::extract_locally_relevant_dofs (dof_handler, locally_relevant_dofs);
    
    for( int i=0; i<no_time_steps; i++ )
      m_all_Psi[i].reinit(locally_owned_dofs, mpi_communicator); // no ghosts
    for( int i=0; i<no_time_steps; i++ )
      m_all_p[i].reinit(locally_owned_dofs, mpi_communicator); // no ghosts

    m_sol.reinit(locally_owned_dofs, mpi_communicator); // no ghosts
    system_rhs.reinit(locally_owned_dofs, mpi_communicator); // no ghosts
    m_p_D.reinit(locally_owned_dofs, mpi_communicator); // no ghosts
    m_p_L.reinit(locally_owned_dofs, mpi_communicator); // no ghosts
    m_workspace.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    m_workspace_2.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    m_workspace_3.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    m_workspace_ng.reinit (locally_owned_dofs, mpi_communicator);
    m_Psi.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    m_Psi_d.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    m_error_per_cell.reinit(triangulation.n_active_cells());

    cout << "(" << m_rank << ") locally_owned_dofs = " << system_rhs.local_size()  << endl;

    constraints.clear ();
    constraints.reinit (locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints (dof_handler, constraints);
    VectorTools::interpolate_boundary_values (dof_handler, 0, ZeroFunction<dim>(2), constraints);
    constraints.close ();

    DynamicSparsityPattern dsp (locally_relevant_dofs);

    DoFTools::make_sparsity_pattern (dof_handler, dsp, constraints, false);
    SparsityTools::distribute_sparsity_pattern (dsp, dof_handler.n_locally_owned_dofs_per_processor(), mpi_communicator, locally_relevant_dofs);
    system_matrix.reinit (locally_owned_dofs, locally_owned_dofs, dsp, mpi_communicator);
    m_computing_timer.exit_section();
  }  

  template <int dim, int no_time_steps>
  void MySolver<dim,no_time_steps>::output_results ( std::string filename ) 
  {
    m_computing_timer.enter_section(__func__);

    vector<std::string> solution_names;
    vector<std::string> solution_names_2;

    constraints.distribute(m_Psi_d);
    constraints.distribute(m_all_Psi[no_time_steps-1]);
    m_workspace = m_Psi_d;
    m_workspace_2 = m_all_Psi[no_time_steps-1];

    DataOut<dim> data_out;
    data_out.attach_dof_handler (dof_handler);
    solution_names.push_back ("Re_Psi_d");
    solution_names.push_back ("Im_Psi_d");    
    data_out.add_data_vector (m_workspace, solution_names );
    solution_names_2.push_back ("Re_Psi_f");
    solution_names_2.push_back ("Im_Psi_f");    
    data_out.add_data_vector (m_workspace_2, solution_names_2 );
    data_out.build_patches ();
    
    //filename = path + "solution-" + to_string(m_timeindex) + ".vtu";
    data_out.write_vtu_in_parallel ( filename.c_str(), mpi_communicator );

    m_computing_timer.exit_section();
  }   

  template <int dim, int no_time_steps>
  void MySolver<dim,no_time_steps>::output_vec ( string filename, LA::MPI::Vector& vec ) 
  {
    m_computing_timer.enter_section(__func__);
    
    vector<string> solution_names;
    m_workspace_3 = vec;
    
    DataOut<dim> data_out;
    data_out.attach_dof_handler (dof_handler);
    solution_names.push_back ("Re");
    solution_names.push_back ("Im");    
    data_out.add_data_vector (m_workspace_3, solution_names);
    data_out.build_patches ();
    data_out.write_vtu_in_parallel ( filename.c_str(), mpi_communicator );

    m_computing_timer.exit_section();
  }     
  
  template <int dim, int no_time_steps>
  void MySolver<dim,no_time_steps>::save( string filename )
  {
    parallel::distributed::SolutionTransfer<dim,LA::MPI::Vector> solution_transfer(dof_handler);
    solution_transfer.prepare_serialization(m_Psi);

    triangulation.save( filename.c_str() );
  }

  template <int dim, int no_time_steps>
  void MySolver<dim,no_time_steps>::load( string filename )
  {
    make_grid();
    triangulation.load( filename.c_str() );
    setup_system();
    
    vector<LA::MPI::Vector*> x_system (2);
    x_system[0] = &m_sol;
    x_system[1] = &system_rhs;
    parallel::distributed::SolutionTransfer<dim,LA::MPI::Vector> solution_transfer(dof_handler);
    solution_transfer.deserialize(x_system);

    constraints.distribute(m_sol);
    constraints.distribute(system_rhs);

    m_all_Psi[0] = m_sol;
    m_Psi_d = system_rhs;
  }    

  template <int dim, int no_time_steps>
  void MySolver<dim,no_time_steps>::solve ()
  {
    m_computing_timer.enter_section(__func__);
    LA::MPI::Vector tmp_vec(locally_owned_dofs, mpi_communicator);
    SolverControl solver_control;
    
    PETScWrappers::SparseDirectMUMPS solver(solver_control, mpi_communicator);
    solver.set_symmetric_mode(false);
    solver.solve(system_matrix, tmp_vec, system_rhs);

    constraints.distribute (tmp_vec);
    m_Psi = tmp_vec;
    m_computing_timer.exit_section();
  }

  template <int dim, int no_time_steps>
  void MySolver<dim,no_time_steps>::run ()
  {
    map<string,double> con_map;
    con_map["omq_x"] = m_omega[0]*m_omega[0];
    con_map["omq_y"] = m_omega[1]*m_omega[1];
    con_map["omq_z"] = m_omega[2]*m_omega[2];

    // V(x,y;lambda,..) 
    vector<string> pot_str;
//    pot_str.push_back("(1+lam_0)*omq_x*x^2+omq_y*y^2"); // test case
//    pot_str.push_back("omq_x*x^2" ); // test case

    pot_str.push_back("omq_x*x^2 + omq_y*y^2 + lam_0*sin(lam_1+lam_2*(x+y))" );
    pot_str.push_back("sin(lam_1+lam_2*(x+y))");
    pot_str.push_back("lam_0*cos(lam_1+lam_2*(x+y))");
    pot_str.push_back("lam_0*cos(lam_1+lam_2*(x+y))*(x+y)");

    double domega = M_PI/m_T;

    pcout << "domega = " << domega << endl;

    // initial guess for lambda
    vector<string> lam_str;
    string str;
    str = "sin(" + to_string(domega) + "*t)"; // test case
    lam_str.push_back(str); // lam_0
    lam_str.push_back(str); // lam_1
    lam_str.push_back(str); // lam_2

    m_potential.init( lam_str, pot_str, con_map, m_T );
    m_potential.output( "lambda_guess.txt" );

    m_norm_grad.resize(m_potential.get_no_lambdas());

    load( "oct_0.bin" );

    double p[] = {0,0,0};
    double pos[] = {0,0,0};
    double var[] = {0,0,0};

    m_workspace = m_all_Psi[0];
    m_N = MyComplexTools::MPI::Particle_Number( mpi_communicator, dof_handler, fe, m_workspace );
    
    pcout << "N(Psi_0) == " << m_N << endl;
    pcout << "N(Psi_d) == " << MyComplexTools::MPI::Particle_Number( mpi_communicator, dof_handler, fe, m_Psi_d ) << endl;
    pcout << "dt == " << m_dt << endl;
//    Expectation_value_position( m_Psi, pos );
//    Expectation_value_momentum( m_Psi, p );

    //pcout << "p == " << p[0]/m_N << ", " << p[1]/m_N << ", " << p[2]/m_N << endl;
    //pcout << "pos == " << pos[0]/m_N << ", " << pos[1]/m_N << ", " << pos[2]/m_N << endl;
    
    for( int i=1; i<=400; i++ )
    {
//      pcout << "Step 1" << endl;
      rt_propagtion_forward(i);
//      m_N = MyComplexTools::MPI::Particle_Number( mpi_communicator, dof_handler, fe, m_workspace );      
//      pcout << "N == " << m_N << endl;
//      Expectation_value_momentum( m_Psi, p );
//      Expectation_value_position( m_Psi, pos );
//      pcout << "p == " << p[0]/m_N << ", " << p[1]/m_N << ", " << p[2]/m_N << endl;
//      pcout << "pos == " << pos[0]/m_N << ", " << pos[1]/m_N << ", " << pos[2]/m_N << endl;      
//      pcout << "Step 2" << endl;
      rt_propagtion_backward(i);
      if( m_N_pT < 1e-4 ) break;
//      pcout << "Step 3" << endl;
      compute_correction(i);
      //double cost = compute_costfunction();
      //if(m_root) printf( "cost = %g\n", cost );
      if(m_root) m_potential.output( "lambda_" + to_string(i) + ".txt" );
    }

    output_results( "oct_final.vtu");
  }
} // end of namespace 

int main ( int argc, char *argv[] )
{
  using namespace dealii;
  deallog.depth_console (0);

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);
  {
    realtime_propagation::MySolver<DIMENSION,101> solver("params_one.xml", 1);
    //realtime_propagation::MySolver<DIMENSION,6> solver("params_one.xml", 0.01);
    solver.run();
  }
return EXIT_SUCCESS;
}
