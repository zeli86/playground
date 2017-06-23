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
#include "muParser.h"
#include "MyComplexTools.h"


namespace realtime_propagation
{
  enum Status { SUCCESS, FAILED };

  using namespace std;
  using namespace dealii;
  
  template <int dim, int N>
  class CPotential : public Function<dim> 
  {
    public:
      // Constructor
      CPotential () 
        : 
        Function<dim>(),
        m_pos_val(3)
      { 
      }
      
      // Destructor
      virtual ~CPotential()
      {
        for( auto& i : m_lambdas )
        {
          delete i;
        }
        for( auto& i : m_pot )
        {
          delete i;
        }
      }

      void init(  const vector<string>& all_lambdas, const vector<string>& all_potential, const map<string,double>& constants, const double T )
      {
        m_dt = T/double(N-1);

        for( int i=0; i<N; i++ )
        {
          m_t[i] = double(i)*m_dt;
        }

        // Setup all initial lambda_i(t) with the guess from all_lambdas
        vector<double> tmpvec(N);
        FunctionParser<1> lam;
        for( auto str : all_lambdas )
        {
          lam.initialize( "lam", str, constants );
          
          for( int i=0; i<N; i++ )
          {
            Point<1> t(double(i)*m_dt);
            tmpvec[i] = lam.value(t);
          }

          m_lambdas.push_back( new Functions::CSpline<1>(m_t, tmpvec) );
        }

        m_lam_val.resize(all_lambdas.size());
        // Setup the potential and all derivatives with respect to lambda_i(t)
        for( auto str : all_potential )
        {
          m_pot.push_back( new mu::Parser() );
          mu::Parser * p = m_pot.back();
          
          p->SetExpr(str);

          p->DefineVar("x", &(m_pos_val.data()[0]));
          p->DefineVar("y", &(m_pos_val.data()[1]));
          p->DefineVar("z", &(m_pos_val.data()[2]));

          for( int i=0; i<all_lambdas.size(); i++ )
          {
            string tmp = "lam_" + to_string(i);
            p->DefineVar(tmp, &(m_lam_val.data()[i]));
          }
          for( auto i : constants )
          {
            p->DefineConst(i.first, i.second);
          }
        }
      }

      virtual double value ( const Point<dim> &p, const unsigned component = 0) 
      {
        for( int i=0; i<dim; i++ )
        {
          m_pos_val[i] = p[i]; // setting the spatial coordinate
        }
        int s=0;
        for( auto i : m_lambdas )
        {
          Point<1> pt(this->get_time());
          m_lam_val[s] = i->value(pt); 
          s++;
        }        
      return m_pot[component]->Eval();
      }      
      
      vector<Functions::CSpline<1>*> m_lambdas;

      protected:
        double m_dt;
        vector<mu::Parser*> m_pot; 
        vector<double> m_t;
        vector<double> m_pos_val;
        vector<double> m_lam_val;
  };
    
  /*************************************************************************************************/
  template <int dim, int no_time_steps>
  class MySolver
  {
  public:
    MySolver( const std::string& );
    ~MySolver();

    void run ();
    
    double Get_dt() { return m_dt; };
    void Set_dt( const double a ) { m_dt=a; m_dth=0.5*a; };
    void output_lambdas( const string& );
    void output_lambdas_grad( const string& );

  protected:
    void Scalar_product( LA::MPI::Vector&, LA::MPI::Vector&, double&, double& );

    CPotential<dim,no_time_steps> m_potential;

    void rt_propagtion_forward (); 
    void rt_propagtion_backward (); 
    void compute_correction (); 
    
    void make_grid();
    void setup_system();
    
    void assemble_system_2(); // required by rt_propagtion_backward
    
    void solve();
    void solve_eq1();
    void solve_eq2();
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
    LA::MPI::Vector m_Psi_t; // Psi trial
    LA::MPI::Vector m_workspace;
    LA::MPI::Vector m_workspace_2;
    LA::MPI::Vector m_workspace_ng;
    Vector<double> m_error_per_cell;
    
    array<LA::MPI::Vector,no_time_steps> m_all_Psi;
    array<LA::MPI::Vector,no_time_steps> m_all_p;

    ConditionalOStream pcout;
    ofstream m_computing_timer_log;
    TimerOutput m_computing_timer;    

    bool m_root;

    double m_gs;
    double m_dt;
    double m_dth; // 0.5*m_dt
    double m_res;
    double m_N;
    double m_overlap;
    vector<double> m_omega;

    double m_xmin, m_xmax;
    double m_ymin, m_ymax;
    double m_zmin, m_zmax;

    unsigned m_NA;
    unsigned m_NK;
    unsigned m_timeindex;
    
    MyTable m_table;    
  };

  template <int dim, int no_time_steps>
  MySolver<dim,no_time_steps>::MySolver ( const std::string& xmlfilename ) 
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

      m_NA = m_ph.Get_Algorithm("NA",0);
      m_NK = m_ph.Get_Algorithm("NK",0);
      m_dt = m_ph.Get_Algorithm("dt",0);    
    }
    catch( const std::string info )
    {
      std::cerr << info << endl;
      MPI_Abort( mpi_communicator, 0 );
    }    
    m_dth = 0.5*m_dt;

    m_root = (Utilities::MPI::this_mpi_process(mpi_communicator) == 0);
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
  void MySolver<dim,no_time_steps>::Scalar_product( LA::MPI::Vector& vec1, LA::MPI::Vector& vec2, double& re, double& im )
  {
    m_computing_timer.enter_section(__func__);
    double tmp[] = {0,0};
    
    constraints.distribute(vec1);
    m_workspace = vec1;

    constraints.distribute(vec2);
    m_workspace_2 = vec2;
    
    const QGauss<dim> quadrature_formula(fe.degree+1);
    FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_JxW_values);

    const unsigned n_q_points = quadrature_formula.size();
    vector<Vector<double>> vec_vals1(n_q_points,Vector<double>(2));
    vector<Vector<double>> vec_vals2(n_q_points,Vector<double>(2));

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for( ; cell!=endc; cell++ )
    {
      if( cell->is_locally_owned() )
      {
        fe_values.reinit (cell);
        fe_values.get_function_values( m_workspace, vec_vals1 );
        fe_values.get_function_values( m_workspace_2, vec_vals2 );
        for ( unsigned qp=0; qp<n_q_points; qp++ )
        {
          tmp[0] += fe_values.JxW(qp)*(vec_vals1[qp][0]*vec_vals2[qp][0]+vec_vals1[qp][1]*vec_vals2[qp][1]);
          tmp[1] += fe_values.JxW(qp)*(vec_vals1[qp][0]*vec_vals2[qp][1]-vec_vals1[qp][1]*vec_vals2[qp][0]);
        }
      }
    }

    double retval[2];
    MPI_Allreduce( tmp, retval, 2, MPI_DOUBLE, MPI_SUM, mpi_communicator);
    re = retval[0];
    im = retval[1];
    m_computing_timer.exit_section();
  }  
  
  template <int dim, int no_time_steps>
  void MySolver<dim,no_time_steps>::make_grid ()
  {
    m_computing_timer.enter_section(__func__);
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
    m_workspace.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    m_workspace_2.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    m_workspace_ng.reinit (locally_owned_dofs, mpi_communicator);
    m_Psi.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    m_Psi_d.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    m_Psi_t.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    m_error_per_cell.reinit(triangulation.n_active_cells());

    int myrank;
    MPI_Comm_rank( mpi_communicator, &myrank );    
    cout << "(" << myrank << ") locally_owned_dofs = " << system_rhs.local_size()  << endl;

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
    m_workspace = vec;
    
    DataOut<dim> data_out;
    data_out.attach_dof_handler (dof_handler);
    solution_names.push_back ("Re_vec");
    solution_names.push_back ("Im_vec");    
    data_out.add_data_vector (m_workspace, solution_names);
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

    m_Psi = m_sol;
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
    //TODO: setup m_Psi(0) and m_Psi_d and lambdas

    map<string,double> con_map;
    con_map["omq_x"] = 0.5;
    con_map["omq_y"] = 0.5;
    con_map["omq_z"] = 0.5;

    // V(x,y;lambda,..) 
    vector<string> pot_str;
    pot_str.push_back("omq_x*x^2 + omq_y*y^2 + lam_0*sin(lam_1*x) + lam_2*sin(lam_3*x)" );
    pot_str.push_back("sin(lam_1*x)" );
    pot_str.push_back("lam_0*cos(lam_1*x)" );
    pot_str.push_back("sin(lam_3*x)" );
    pot_str.push_back("lam_2*cos(lam_3*x)" );

    // initial guess for lambda
    vector<string> lam_str;
    lam_str.push_back("");

    m_potential.init( lam_str, pot_str, con_map, 5 );

    load( "oct_0.bin" );

    double p[] = {0,0,0};
    double pos[] = {0,0,0};
    double var[] = {0,0,0};

    output_results("bla.vtu");

    m_workspace = m_Psi;
    m_N = MPI::MyComplexTools::Particle_Number( mpi_communicator, dof_handler, fe, m_workspace );
    
    pcout << "N == " << m_N << endl;
//    Expectation_value_position( m_Psi, pos );
//    Expectation_value_momentum( m_Psi, p );

    pcout << "N == " << m_N << endl;
    pcout << "p == " << p[0]/m_N << ", " << p[1]/m_N << ", " << p[2]/m_N << endl;
    pcout << "pos == " << pos[0]/m_N << ", " << pos[1]/m_N << ", " << pos[2]/m_N << endl;
    
    if(m_root) output_lambdas( "lambda_0.txt" );
    
    m_all_Psi[0] = m_Psi;
    
    for( int i=1; i<=400; i++ )
    {
      m_Psi = m_all_Psi[0];
      pcout << "Step 1" << endl;
      rt_propagtion_forward();
//      Expectation_value_momentum( m_Psi, p );
//      Expectation_value_position( m_Psi, pos );
      pcout << "p == " << p[0]/m_N << ", " << p[1]/m_N << ", " << p[2]/m_N << endl;
      pcout << "pos == " << pos[0]/m_N << ", " << pos[1]/m_N << ", " << pos[2]/m_N << endl;      
      pcout << "Step 2" << endl;
      rt_propagtion_backward();
      pcout << "Step 3" << endl;
      compute_correction();
      //double cost = compute_costfunction();
      //if(m_root) printf( "cost = %g\n", cost );
      if(m_root) output_lambdas( "lambda_" + to_string(i) + ".txt" );
    }

    output_results( "oct_final.vtu");

  }
  
  template <int dim, int no_time_steps>
  void MySolver<dim,no_time_steps>::output_lambdas( const string& filename )
  {
    ofstream out( filename );
/*    
    for( int i=0; i<no_time_steps; i++ )
    {
      out << double(i)*m_dt << "\t";
      for( int j=0; j<no_lam; j++ ) 
      {
      }
    }    
*/
  }

  template <int dim, int no_time_steps>
  void MySolver<dim,no_time_steps>::output_lambdas_grad( const string& filename )
  {
    ofstream out( filename );
    /*
    for( int i=0; i<no_time_steps; i++ )
    {
      out << double(i)*m_dt << "\t";
      for( int j=0; j<no_lam; j++ ) 
      {
      }
    }
    */
  }  
} // end of namespace 

int main ( int argc, char *argv[] )
{
  using namespace dealii;
  deallog.depth_console (0);

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);
  {
    realtime_propagation::MySolver<DIMENSION,31> solver("params_one.xml");
    //solver.Set_dt(dt);
    //realtime_propagation::Setup_initial_lambdas<DIMENSION,NT,no_lam>(solver);    
    solver.run();
  }
return EXIT_SUCCESS;
}
