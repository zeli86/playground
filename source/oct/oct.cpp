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

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/function_cspline.h>

#include <deal.II/lac/generic_linear_algebra.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/constraint_matrix.h>

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
#include <deal.II/numerics/matrix_tools.h>
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
#include <iomanip>
#include <array>

#include "global.h"
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
     void assemble_system( const int ); // required by DoIter
    
    void solve();
    void solve_eq1();
    void output_results ( string );
    void output_vec ( string, Vector<double>& );

    MyParameterHandler m_ph;
    Triangulation<dim> triangulation;
    FESystem<dim> fe;
    DoFHandler<dim> dof_handler;
    ConstraintMatrix constraints;
    SparsityPattern sparsity_pattern;
    SparseMatrix<double> system_matrix;
 
    Vector<double> system_rhs;
    Vector<double> sol;
    Vector<double> m_Psi; // Psi(t)
    Vector<double> m_Psi_d; // desired state
    Vector<double> m_workspace;
    Vector<double> m_workspace_2;
    Vector<double> m_error_per_cell;
    
    array<Vector<double>,no_time_steps> m_all_Psi;
    array<Vector<double>,no_time_steps> m_all_p;

    double m_gs;
    double m_dt;
    double m_dth; // 0.5*m_dt
    vector<double> m_omega;
    vector<double> m_norm_grad;

    double m_N;
    double m_T;

    double m_xmin, m_xmax;
    double m_ymin, m_ymax;
   
    MyTable m_table;    
  };

  template <int dim, int no_time_steps>
  MySolver<dim,no_time_steps>::MySolver ( const std::string& xmlfilename, const double T ) 
    : 
    m_ph(xmlfilename),
    triangulation (),
    fe (FE_Q<dim>(gl_degree_fe), 2),
    dof_handler (triangulation)
  {
    try
    {
      m_omega = m_ph.Get_Physics("omega");
      m_gs = m_ph.Get_Physics("gs_1",0);
      
      m_xmin = m_ph.Get_Mesh("xrange",0);
      m_xmax = m_ph.Get_Mesh("xrange",1);
      m_ymin = m_ph.Get_Mesh("yrange",0);
      m_ymax = m_ph.Get_Mesh("yrange",1);
    }
    catch( const std::string info )
    {
      std::cerr << info << endl;
      exit(0);
    }    
    m_T = T;   
    m_dt = T/double(no_time_steps-1);
    m_dth = 0.5*m_dt;
  }

  template <int dim, int no_time_steps>
  MySolver<dim,no_time_steps>::~MySolver ()
  {
    dof_handler.clear ();
  }
  
  #include "eq1.h"
  #include "eq2.h"
  #include "eq3.h"
  
  template <int dim, int no_time_steps>
  void MySolver<dim,no_time_steps>::make_grid ()
  {
    Point<dim,double> pt1;
    Point<dim,double> pt2;

    double min[] = {m_xmin, m_ymin};
    double max[] = {m_xmax, m_ymax};

    for( int i=0; i<dim; i++ )
    {
      pt1(i) = min[i];
      pt2(i) = max[i];
    }

    GridGenerator::hyper_rectangle(triangulation, pt2, pt1);
    triangulation.refine_global(7);
    
    double min_cell_diameter = GridTools::minimal_cell_diameter(triangulation);
    double max_cell_diameter = GridTools::maximal_cell_diameter(triangulation);
    cout << "min_cell_diameter     = " << min_cell_diameter << "\n";
    cout << "max_cell_diameter     = " << max_cell_diameter << "\n";
    cout << "total_no_cells        = " << triangulation.n_cells() << "\n";
    cout << "total_no_active_cells = " << triangulation.n_active_cells() << "\n";
  }
  
  template <int dim, int no_time_steps>
  void MySolver<dim,no_time_steps>::setup_system()
  {
    dof_handler.distribute_dofs (fe);
    //DoFRenumbering::component_wise (dof_handler);

    for( int i=0; i<no_time_steps; i++ )
      m_all_Psi[i].reinit(dof_handler.n_dofs());
    for( int i=0; i<no_time_steps; i++ )
      m_all_p[i].reinit(dof_handler.n_dofs());

    sol.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());
    m_workspace.reinit (dof_handler.n_dofs());
    m_workspace_2.reinit (dof_handler.n_dofs());
    m_Psi.reinit (dof_handler.n_dofs());
    m_Psi_d.reinit (dof_handler.n_dofs());
    m_error_per_cell.reinit(triangulation.n_active_cells());

    constraints.clear ();
    DoFTools::make_hanging_node_constraints (dof_handler, constraints);
    VectorTools::interpolate_boundary_values (dof_handler, 0,ZeroFunction<dim>(2), constraints);
    constraints.close ();

    DynamicSparsityPattern dsp(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern (dof_handler, dsp);

    sparsity_pattern.copy_from(dsp);
    system_matrix.reinit (sparsity_pattern);
  }  

  template <int dim, int no_time_steps>
  void MySolver<dim,no_time_steps>::output_results ( string path ) 
  {
    string filename;
    
    vector<std::string> solution_names;

    DataOut<1> data_out;
    data_out.attach_dof_handler (dof_handler);
    solution_names.push_back ("Re Psi");
    solution_names.push_back ("Im Psi");    
    data_out.add_data_vector (m_Psi, solution_names);
    data_out.build_patches ();
    
    //filename = path + "solution-" + to_string(m_timeindex) + ".gnuplot";
    ofstream output (filename);
    data_out.write_gnuplot ( output );
  }   

  template <int dim, int no_time_steps>
  void MySolver<dim,no_time_steps>::output_vec ( string filename, Vector<double>& vec ) 
  {
    vector<string> solution_names;

    DataOut<1> data_out;
    data_out.attach_dof_handler (dof_handler);
    solution_names.push_back ("Re vec");
    solution_names.push_back ("Im vec");    
    data_out.add_data_vector (vec, solution_names);
    data_out.build_patches ();

    ofstream output (filename);
    data_out.write_gnuplot ( output );
  }     

  template <int dim, int no_time_steps>
  void MySolver<dim,no_time_steps>::solve ()
  {
    SparseDirectUMFPACK A_direct;
    A_direct.initialize(system_matrix);
    A_direct.vmult (m_Psi, system_rhs);
  }

  template <int dim, int no_time_steps>
  void MySolver<dim,no_time_steps>::run ()
  {
    make_grid();
    setup_system();
    
    ifstream in("Psi_0.bin");
    m_Psi.block_read(in);
    
    ifstream in2("Psi_d.bin");
    m_Psi_d.block_read(in2);
    
    output_vec( "Psi_0.gnuplot", m_Psi );
    output_vec( "Psi_d.gnuplot", m_Psi_d );

    map<string,double> con_map;
    con_map["omq_x"] = m_omega[0];
    con_map["omq_y"] = m_omega[1];
    con_map["omq_z"] = m_omega[2];

    // V(x,y;lambda,..) 
    vector<string> pot_str;
//    pot_str.push_back("omq_x*x^2 + omq_y*y^2 + lam_0*sin(lam_1*x) + lam_2*sin(lam_3*y)" );
//    pot_str.push_back("sin(lam_1*x)" );
//    pot_str.push_back("lam_0*cos(lam_1*x)" );
//    pot_str.push_back("sin(lam_3*y)" );
//    pot_str.push_back("lam_2*cos(lam_3*y)" );

    pot_str.push_back("(1+lam_0)*omq_x*x^2 + (1+lam_1)*y^2 + lam_2*x^3 + lam_3*y^3" );
    pot_str.push_back("x^2" );
    pot_str.push_back("y^2" );
    pot_str.push_back("x^3" );
    pot_str.push_back("y^3" );

    double domega = M_PI/m_T;

    // initial guess for lambda
    vector<string> lam_str;
    string str;
    str = "sin(lam*" + to_string(domega) + ")";
    lam_str.push_back(str);
    lam_str.push_back(str);
    str = "0";
    lam_str.push_back(str);
    lam_str.push_back(str);

    m_potential.init( lam_str, pot_str, con_map, m_T );
    m_potential.output( "lambda_guess.txt" );

    m_norm_grad.resize(m_potential.get_no_lambdas());

    double p[] = {0,0,0};
    double pos[] = {0,0,0};
    double var[] = {0,0,0};

    m_workspace = m_all_Psi[0];
    m_N = MyComplexTools::Particle_Number( dof_handler, fe, m_workspace );
    
    cout << "N == " << m_N << endl;
    cout << "dt == " << m_dt << endl;
//    Expectation_value_position( m_Psi, pos );
//    Expectation_value_momentum( m_Psi, p );

    //cout << "p == " << p[0]/m_N << ", " << p[1]/m_N << ", " << p[2]/m_N << endl;
    //cout << "pos == " << pos[0]/m_N << ", " << pos[1]/m_N << ", " << pos[2]/m_N << endl;
    
    for( int i=1; i<=10; i++ )
    {
      cout << "Step 1" << endl;
      rt_propagtion_forward(i);
      m_N = MyComplexTools::Particle_Number( dof_handler, fe, m_workspace );      
      cout << "N == " << m_N << endl;
//      Expectation_value_momentum( m_Psi, p );
//      Expectation_value_position( m_Psi, pos );
//      cout << "p == " << p[0]/m_N << ", " << p[1]/m_N << ", " << p[2]/m_N << endl;
//      cout << "pos == " << pos[0]/m_N << ", " << pos[1]/m_N << ", " << pos[2]/m_N << endl;      
      cout << "Step 2" << endl;
      rt_propagtion_backward(i);
      cout << "Step 3" << endl;
      compute_correction(i);
      m_potential.output( "lambda_" + to_string(i) + ".txt" );
    }

    //output_results( "oct_final.vtu");     
  } 
} // end of namespace 

int main ( int argc, char *argv[] )
{
  dealii::deallog.depth_console (0);

  realtime_propagation::MySolver<1,101> solver("params.xml",1);
  solver.run();

return EXIT_SUCCESS;
}
