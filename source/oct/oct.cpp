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
//#include <deal.II/lac/affine_constraints.h>

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

#include <deal.II/base/utilities.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/function_parser.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/grid_refinement.h>

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
#include "cxxopts.hpp"

namespace realtime_propagation
{
  enum Status { SUCCESS, FAILED };

  using namespace std;
  using namespace dealii;
  
  #include "potential.hpp"    

  template <int dim, int no_time_steps>
  class MySolver
  {
  public:
    MySolver( const std::string&, const double );
    ~MySolver();

    void run ( const std::string&, const std::string& );
    void run2 ( const std::string&, const std::string& );
    
  protected:
    CPotential<dim> m_potential;
    CPotential<dim> m_potential_backup;

    void rt_propagtion_forward ( const int ); 
    void rt_propagtion_backward ( const int ); 
    void assemble_system( const int );
    void compute_correction ( const int ); 
    void compute_cost( double& );
    void compute_difference( const LAPACKFullMatrix<double>&, const LAPACKFullMatrix<double>&, LAPACKFullMatrix<double>& );
    double compute_dot_product( const LAPACKFullMatrix<double>&, const LAPACKFullMatrix<double>& );

    void make_grid();
    void setup_system();
    void one_loop( const int );
    void solve();
    void output_results ( string );
    void output_vec ( string, Vector<double>& );
    void compute_beta();
    MyParameterHandler m_ph;
    Triangulation<dim> triangulation;
    FESystem<dim> fe;
    DoFHandler<dim> dof_handler;
    SparsityPattern sparsity_pattern;
    SparseMatrix<double> system_matrix;
 
    Vector<double> system_rhs;
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
    vector<double> m_beta;
    LAPACKFullMatrix<double> m_grad;
    LAPACKFullMatrix<double> m_old_grad;
    LAPACKFullMatrix<double> m_direction;
    LAPACKFullMatrix<double> m_old_direction;
    LAPACKFullMatrix<double> m_diff_grads;

    double m_N;
    double m_T;
    double m_cost;
    double m_old_cost;

    double m_xmin, m_xmax;
    double m_ymin, m_ymax;

    unsigned m_global_refinement;
   
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
      m_global_refinement = unsigned(m_ph.Get_Mesh("global_refinements",0));
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
  
  #include "eq1.hpp"
  #include "eq2.hpp"
  #include "eq3.hpp"
  
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
    triangulation.refine_global(m_global_refinement);
    
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

    system_rhs.reinit(dof_handler.n_dofs());
    m_workspace.reinit (dof_handler.n_dofs());
    m_workspace_2.reinit (dof_handler.n_dofs());
    m_Psi.reinit (dof_handler.n_dofs());
    m_Psi_d.reinit (dof_handler.n_dofs());
    m_error_per_cell.reinit(triangulation.n_active_cells());

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
    map<types::global_dof_index, double> boundary_values;
    VectorTools::interpolate_boundary_values (dof_handler, 0, ZeroFunction<dim>(), boundary_values);
    MatrixTools::apply_boundary_values (boundary_values, system_matrix, m_Psi, system_rhs);    

    SparseDirectUMFPACK A_direct;
    A_direct.initialize(system_matrix);
    A_direct.vmult (m_Psi, system_rhs);
  }

  template <int dim, int no_time_steps>
  void MySolver<dim,no_time_steps>::one_loop( const int ti )
  {
    rt_propagtion_forward(ti);
    rt_propagtion_backward(ti);
    compute_correction(ti);
  }

  template <int dim, int no_time_steps>
  void MySolver<dim,no_time_steps>::run ( const std::string& fn_i, const std::string& fn_d )
  {
    make_grid();
    setup_system();
    
    ifstream in(fn_i);
    m_Psi.block_read(in);
    
    ifstream in2(fn_d);
    m_Psi_d.block_read(in2);
    
    output_vec( "Psi_0.gnuplot", m_Psi );
    output_vec( "Psi_d.gnuplot", m_Psi_d );

    m_workspace = m_Psi_d;
    m_workspace -= m_Psi;

    double cost = MyComplexTools::Particle_Number( dof_handler, fe, m_workspace );
    printf( "initial cost = %g\n", cost );   
    printf( "N Psi_0 = %g\n", MyComplexTools::Particle_Number( dof_handler, fe, m_Psi ) );   
    printf( "N Psi_d = %g\n", MyComplexTools::Particle_Number( dof_handler, fe, m_Psi_d ) );   
    
    output_vec( "Psi_0.gnuplot", m_Psi );
    output_vec( "Psi_d.gnuplot", m_Psi_d );

    map<string,double> con_map;
    con_map["omq_x"] = m_omega[0];

    // V(x,y;lambda,..) 
    vector<string> pot_str;
    pot_str.push_back("(1+lam_0)*omq_x*(x-lam_1)^2 + lam_2*(x-lam_3)^3");
    pot_str.push_back("omq_x*lam_0*(x-lam_1)^2");
    pot_str.push_back("-2*(1+lam_0)*omq_x*(x-lam_1)");
    pot_str.push_back("(x-lam_3)^3");
    pot_str.push_back("-3*lam_2*(x-lam_3)^2");

    double domega = M_PI/m_T;

    // initial guess for lambda
    vector<string> lam_str;
    string str;
    str = "sin(t*" + to_string(domega) + ")";
    lam_str.push_back(str);
    str = "0";
    lam_str.push_back(str);
    lam_str.push_back(str);
    lam_str.push_back(str);
    
    m_potential.init( lam_str, pot_str, con_map, m_T, no_time_steps );
    m_potential.output( "lambda_guess.txt" );
    m_potential_backup = m_potential;
    m_potential_backup.output( "lambda_guess_backup.txt" );

    m_norm_grad.resize(m_potential.get_no_lambdas());
    m_grad.reinit(no_time_steps,m_potential.get_no_lambdas());
    m_old_grad.reinit(no_time_steps,m_potential.get_no_lambdas());
    m_direction.reinit(no_time_steps,m_potential.get_no_lambdas());
    m_old_direction.reinit(no_time_steps,m_potential.get_no_lambdas());
    m_diff_grads.reinit(no_time_steps,m_potential.get_no_lambdas());
    m_beta.resize(no_time_steps,0);

    m_all_Psi[0] = m_Psi;
    cout << "dt == " << m_dt << endl;
    
    for( int i=1; i<400; i++ )
    {
      one_loop(i);

      m_potential.add( 0.1, m_grad );

      printf("(%d) cost = %g\n", i, m_cost );

      m_potential.output( "lambda_" + to_string(i) + ".txt" );
      output_vec( "res_" + to_string(i) + ".gnuplot", m_all_Psi[no_time_steps-1] );
    }

    //output_results( "oct_final.vtu");     
  } 

  template <int dim, int no_time_steps>
  void MySolver<dim,no_time_steps>::run2 ( const std::string& fn_i, const std::string& fn_d )
  {
    make_grid();
    setup_system();
    
    ifstream in(fn_i);
    m_Psi.block_read(in);
    
    ifstream in2(fn_d);
    m_Psi_d.block_read(in2);
    
    output_vec( "Psi_0.gnuplot", m_Psi );
    output_vec( "Psi_d.gnuplot", m_Psi_d );

    map<string,double> con_map;
    con_map["omq_x"] = m_omega[0];

    // V(x,y;lambda,..) 
    vector<string> pot_str;
/*
    pot_str.push_back("(1+lam_0)*omq_x*x^2 + lam_1*x^3");
    pot_str.push_back("x^2");
    pot_str.push_back("x^3");
*/
    pot_str.push_back("omq_x*x + lam_0*sin(lam_1*x+lam_2)");
    pot_str.push_back("sin(lam_1*x+lam_2)");
    pot_str.push_back("lam_0*cos(lam_1*x+lam_2)*x");
    pot_str.push_back("lam_0*cos(lam_1*x+lam_2)");

    double domega = M_PI/m_T;

    // initial guess for lambda
    vector<string> lam_str;
    string str;
    str = "sin(t*" + to_string(domega) + ")";
    lam_str.push_back(str);
    str = "0";
    lam_str.push_back(str);

    m_potential.init( lam_str, pot_str, con_map, m_T, no_time_steps );
    m_potential.output( "lambda_guess.txt" );
    m_potential_backup = m_potential;
    m_potential_backup.output( "lambda_guess_backup.txt" );

    m_norm_grad.resize(m_potential.get_no_lambdas());
    m_grad.reinit(no_time_steps,m_potential.get_no_lambdas());
    m_old_grad.reinit(no_time_steps,m_potential.get_no_lambdas());
    m_direction.reinit(no_time_steps,m_potential.get_no_lambdas());
    m_old_direction.reinit(no_time_steps,m_potential.get_no_lambdas());
    m_beta.resize(no_time_steps,0);

    double p[3] = {};
    double pos[3] = {};
    double var[3] = {};

    m_all_Psi[0] = m_Psi;
    m_N = MyComplexTools::Particle_Number( dof_handler, fe, m_Psi );
    
    cout << "N == " << m_N << endl;
    cout << "dt == " << m_dt << endl;
//    Expectation_value_position( m_Psi, pos );
//    Expectation_value_momentum( m_Psi, p );

    //cout << "p == " << p[0]/m_N << ", " << p[1]/m_N << ", " << p[2]/m_N << endl;
    //cout << "pos == " << pos[0]/m_N << ", " << pos[1]/m_N << ", " << pos[2]/m_N << endl;
    
    one_loop(0); // compute the first gradient
    m_old_cost = m_cost;
    m_direction = m_grad;
    for( int i=1; i<40; i++ )
    {
      m_potential_backup = m_potential;

      double tau=1;
      for( int j=0; j<10; j++ )
      {
        m_potential.add( tau, m_direction );
        rt_propagtion_forward(0);
        if( m_cost - m_old_cost > tau*compute_dot_product(m_direction,m_grad) )
        {
          m_potential = m_potential_backup;
          tau *= 0.5;
        }
      }

      m_old_grad = m_grad;
      m_old_direction = m_direction;
      rt_propagtion_backward(0);
      compute_correction(0);
      compute_beta();

      for( int s=0; s<m_direction.n(); s++ )
        for( int ti=1; ti<m_direction.m()-1; ti++ )
          m_direction(ti,s) = m_grad(ti,s) - m_beta[s]*m_old_direction(ti,s);

      m_potential.output( "lambda_" + to_string(i) + ".txt" );
    }

    //output_results( "oct_final.vtu");     
  } 
} // end of namespace 

int main ( int argc, char *argv[] )
{
  dealii::deallog.depth_console (0);

  cxxopts::Options options("binR_to_atus2", "Converts real deal.ii binary format to atus2 binary format.");
  
  options.add_options()
  ("p,params", "input parameter xml file" , cxxopts::value<std::string>()->default_value("params_one.xml") )
  ("i,input", "input initial wave function" , cxxopts::value<std::string>() )
  ("d,desired", "input desired wave function" , cxxopts::value<std::string>() )
  ("help","Print help")
  ;
  
  auto result = options.parse(argc, argv);

  std::string bin_filename_i, bin_filename_d, params_filename;
  try
  {
    if (result.count("i") > 0 && result.count("d") > 0 && result.count("p") > 0 )
    {
      bin_filename_i = result["i"].as<std::string>(); 
      bin_filename_d = result["d"].as<std::string>(); 
      params_filename = result["p"].as<std::string>();
    }
    else
    {      std::cout << options.help({""}) << std::endl;
      return EXIT_FAILURE;
    }
  }
  catch (const cxxopts::OptionException& e)
  {
    std::cout << "error parsing options: " << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  realtime_propagation::MySolver<1,101> solver(params_filename,1);
  solver.run(bin_filename_i,bin_filename_d);

return EXIT_SUCCESS;
}
