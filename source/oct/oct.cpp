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

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/timer.h>

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

#include "MyParameterHandler.h"
#include "my_table.h"

namespace realtime_propagation
{
  enum Status { SUCCESS, FAILED };

  using namespace std;
  using namespace dealii;
  
  template <int no_time_steps, int no_lam>
  class CPotential : public Function<1> 
  {
    public:
      CPotential ( array<array<double,no_time_steps>,no_lam>& lam, const int timeindex, const int sel=0 ) : Function<1>(), m_lam(lam) 
      { 
        m_sel = sel;
        m_timeindex = timeindex;
      }
      
      ~CPotential()
      {
      }
      virtual double value ( const Point<1> &p, const unsigned component = 0) const;
      
      int m_sel;
      int m_timeindex;
    protected:
      array<array<double,no_time_steps>,no_lam>& m_lam;
  };
    
  /*************************************************************************************************/
  template <int no_time_steps, int no_lam>
  double CPotential<no_time_steps,no_lam>::value( const Point<1> &p, const unsigned component ) const
  {
    double retval;

    const double lam_1 = m_lam[0][m_timeindex];
    
    switch(m_sel) 
    {
      case 0: retval =  0.5*(lam_1-p(0))*(lam_1-p(0)); //0.5*(lam_1-p(0))*(lam_1-p(0)) + 0.5*p(1)*p(1);
      break;
      case 1: retval = lam_1-p(0);
        break;
//       case 2: retval = p(1)*p(1);
//         break;
    };
  return retval;
  }

  template <int no_time_steps, int no_lam>
  class MySolver
  {
  public:
    MySolver( const std::string& );
    ~MySolver();

    void run ();
    double Particle_Number( Vector<double>& );
    void Expectation_value_momentum( Vector<double>&, double* );
    void Expectation_value_position( Vector<double>&, double* );
    
    double Get_dt() { return m_dt; };
    void Set_dt( const double a ) { m_dt=a; m_dth=0.5*a; };
    void output_lambdas( const string& );
    void output_lambdas_grad( const string& );

    array<array<double,no_time_steps>,no_lam> m_all_lambdas;
    array<array<double,no_time_steps>,no_lam> m_all_lambdas_grad;
  protected:
    void Scalar_product( Vector<double>&, Vector<double>&, double&, double& );

    void rt_propagtion_forward (); 
    void rt_propagtion_backward (); 
    void compute_initial_p ();
    void compute_correction (); 
    
    void compute_all_lambdas_t();
    void compute_all_lambdas_tt();
    
    double compute_costfunction();

    void make_grid();
    void setup_system();
    void DoIter();
    
    void assemble_rhs(); // required by DoIter
    void assemble_system(); // required by DoIter
    void assemble_system_4_initial_p( const double, const double );
    void assemble_system_2(); // required by rt_propagtion_backward
    void assemble_system_3(); // required by compute_initial_p
    
    void solve();
    void solve_eq1();
    void solve_eq2();
    void output_results ( string );
    void output_vec ( string, Vector<double>& );

    MyParameterHandler m_ph;
    Triangulation<1> triangulation;
    FESystem<1> fe;
    DoFHandler<1> dof_handler;
    SparsityPattern sparsity_pattern;
    SparseMatrix<double> system_matrix;
 
    Vector<double> system_rhs;
    Vector<double> sol;
    Vector<double> m_Psi; // Psi(t)
    Vector<double> m_Psi_d; // desired state
    Vector<double> m_Psi_t; // Psi trial
    Vector<double> m_workspace;
    Vector<double> m_workspace_2;
    Vector<double> m_error_per_cell;
    
    array<Vector<double>,no_time_steps> m_all_Psi;
    array<Vector<double>,no_time_steps> m_all_p;
    array<array<double,no_time_steps>,no_lam> m_all_lambdas_t;

    double m_gs;
    double m_dt;
    double m_dth; // 0.5*m_dt
    vector<double> m_omega;
    double m_res;
    double m_res_old;
    double m_N;
    double m_overlap;

    double m_xmin;
    double m_xmax;

    unsigned m_NA;
    unsigned m_NK;
    
    int m_timeindex;
    
    MyTable m_table;    
  };

  template <int no_time_steps, int no_lam>
  MySolver<no_time_steps,no_lam>::MySolver ( const std::string& xmlfilename ) 
    : 
    m_ph(xmlfilename),
    triangulation (),
    fe (FE_Q<1>(1), 2),
    dof_handler (triangulation)
  {
    try
    {
      m_omega = m_ph.Get_Physics("omega");
      m_gs = m_ph.Get_Physics("gs_1",0);
      
      m_xmin = m_ph.Get_Mesh("xrange",0);
      m_xmax = m_ph.Get_Mesh("xrange",1);

      m_NA = m_ph.Get_Algorithm("NA",0);
      m_NK = m_ph.Get_Algorithm("NK",0);
      m_dt = m_ph.Get_Algorithm("dt",0);    
    }
    catch( const std::string info )
    {
      std::cerr << info << endl;
      exit(0);
    }    
    m_dth = 0.5*m_dt;
  }

  template <int no_time_steps, int no_lam>
  MySolver<no_time_steps,no_lam>::~MySolver ()
  {
    dof_handler.clear ();
  }
  
  #include "expectation_values.h"
  #include "eq1.h"
  #include "eq2.h"
  #include "eq3.h"
//   #include "cost.h"

  template <int no_time_steps, int no_lam>
  void MySolver<no_time_steps,no_lam>::Scalar_product( Vector<double>& vec1, Vector<double>& vec2, double& re, double& im )
  {
    re=0;
    im=0;
    
    const QGauss<1>  quadrature_formula(fe.degree+1);
    FEValues<1> fe_values (fe, quadrature_formula, update_values|update_JxW_values);

    const unsigned int n_q_points = quadrature_formula.size();
    vector<Vector<double>> vec_vals1(n_q_points,Vector<double>(2));
    vector<Vector<double>> vec_vals2(n_q_points,Vector<double>(2));

    DoFHandler<1>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
      fe_values.reinit (cell);
      fe_values.get_function_values( vec1, vec_vals1 );
      fe_values.get_function_values( vec2, vec_vals2 );
      for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
      {
        re += fe_values.JxW(q_point)*(vec_vals1[q_point][0]*vec_vals2[q_point][0]+vec_vals1[q_point][1]*vec_vals2[q_point][1]);
        im += fe_values.JxW(q_point)*(vec_vals1[q_point][0]*vec_vals2[q_point][1]-vec_vals1[q_point][1]*vec_vals2[q_point][0]);
      }
    }
  }  
  
  template <int no_time_steps, int no_lam>
  void MySolver<no_time_steps,no_lam>::make_grid ()
  {
    Point<1,double> pt1( m_xmin );
    Point<1,double> pt2( m_xmax );

    GridGenerator::hyper_rectangle(triangulation, pt2, pt1);
    triangulation.refine_global(1);
    
    double min_cell_diameter = GridTools::minimal_cell_diameter(triangulation);
    double max_cell_diameter = GridTools::maximal_cell_diameter(triangulation);
    cout << "min_cell_diameter     = " << min_cell_diameter << "\n";
    cout << "max_cell_diameter     = " << max_cell_diameter << "\n";
    cout << "total_no_cells        = " << triangulation.n_cells() << "\n";
    cout << "total_no_active_cells = " << triangulation.n_active_cells() << "\n";
  }
  
  template <int no_time_steps, int no_lam>
  void MySolver<no_time_steps,no_lam>::setup_system()
  {
    dof_handler.distribute_dofs (fe);
    DoFRenumbering::component_wise (dof_handler);

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
    m_Psi_t.reinit (dof_handler.n_dofs());
    m_error_per_cell.reinit(triangulation.n_active_cells());

    system_rhs=0;

    DynamicSparsityPattern dsp(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern (dof_handler, dsp);

    sparsity_pattern.copy_from(dsp);
    system_matrix.reinit (sparsity_pattern);
  }  

  template <int no_time_steps, int no_lam>
  void MySolver<no_time_steps,no_lam>::output_results ( string path ) 
  {
    string filename;
    
    vector<std::string> solution_names;

    DataOut<1> data_out;
    data_out.attach_dof_handler (dof_handler);
    solution_names.push_back ("Re Psi");
    solution_names.push_back ("Im Psi");    
    data_out.add_data_vector (m_Psi, solution_names);
    data_out.build_patches ();
    
    filename = path + "solution-" + to_string(m_timeindex) + ".gnuplot";
    ofstream output (filename);
    data_out.write_gnuplot ( output );
  }   

  template <int no_time_steps, int no_lam>
  void MySolver<no_time_steps,no_lam>::output_vec ( string filename, Vector<double>& vec ) 
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

  template <int no_time_steps, int no_lam>
  void MySolver<no_time_steps,no_lam>::solve ()
  {
    SparseDirectUMFPACK A_direct;
    A_direct.initialize(system_matrix);
    A_direct.vmult (m_Psi, system_rhs);
  }

  template <int no_time_steps, int no_lam>
  void MySolver<no_time_steps,no_lam>::run ()
  {
    make_grid();
    setup_system();
    
    ifstream in("Psi_0.bin");
    m_Psi.block_read(in);
    
    ifstream in2("Psi_d.bin");
    m_Psi_d.block_read(in2);
    
    output_vec( "Psi_0.gnuplot", m_Psi );
    output_vec( "Psi_d.gnuplot", m_Psi_d );
    
    double p, pos;

    output_results("");

    m_N = Particle_Number(m_Psi);
    cout << "N == " << m_N << endl;
    Expectation_value_position( m_Psi, &pos );
    Expectation_value_momentum( m_Psi, &p );

    cout << "N == " << m_N << endl;
    cout << "p == " << p/m_N << endl;
    cout << "pos == " << pos/m_N << endl;
    
//     CPotential<no_time_steps,no_lam> Potential ( m_all_lambdas, m_dt, 2.8 );
//     VectorTools::interpolate(dof_handler, Potential, m_workspace );;
//     output_vec( "pot.vtu", m_workspace );

    output_lambdas( "lambda_0.txt" );
    
    m_all_Psi[0] = m_Psi;
    
    for( int i=1; i<=1000; i++ )
    {
      printf( "----- counter == %d\n", i );
      m_Psi = m_all_Psi[0];
      rt_propagtion_forward();
      Expectation_value_momentum( m_Psi, &p );
      Expectation_value_position( m_Psi, &pos );
      cout << "p == " << p/m_N << endl;
      cout << "pos == " << pos/m_N << endl;      
      rt_propagtion_backward();
      compute_correction();
//       double cost = compute_costfunction();
//       printf( "cost = %g\n", cost );
      output_lambdas( "lambda_" + to_string(i) + ".txt" );
      if( fabs(m_overlap-m_N) < 0.01 ) break;
    }
  }
  
  template <int no_time_steps, int no_lam>
  void MySolver<no_time_steps,no_lam>::output_lambdas( const string& filename )
  {
    ofstream out( filename );
    
    for( int i=0; i<no_time_steps; i++ )
    {
      out << double(i)*m_dt << "\t";
      for( int j=0; j<no_lam; j++ ) 
      {
        out << m_all_lambdas[j][i] << ((j==no_lam-1)? "\n" : "\t");
      }
    }
  }

  template <int no_time_steps, int no_lam>
  void MySolver<no_time_steps,no_lam>::output_lambdas_grad( const string& filename )
  {
    ofstream out( filename );
    
    for( int i=0; i<no_time_steps; i++ )
    {
      out << double(i)*m_dt << "\t";
      for( int j=0; j<no_lam; j++ ) 
      {
        out << m_all_lambdas_grad[j][i] << ((j==no_lam-1)? "\n" : "\t");
      }
    }
  }
  
  template <int no_time_steps, int no_lam>
  void Setup_initial_lambdas( MySolver<no_time_steps,no_lam>& solver  )
  {
    const double dt=solver.Get_dt();
    const double T=(no_time_steps-1)*dt;
    const double alp=6.0/T;
    const double x0=-3.0;
    
    for( int s=0; s<no_lam; s++ )
      for( int i=0; i<no_time_steps; i++ )
      {
        double t = double(i)*dt;
        solver.m_all_lambdas[s][i] = t*alp+x0;
        //solver.m_all_lambdas[s][i] = 3*sin(0.523598775598*t);
      }
  }  
} // end of namespace 

int main ( int argc, char *argv[] )
{
  using namespace dealii;
  deallog.depth_console (0);

  const int no_lam=1;
  const int NT = 61;
  const double T=3.0;
  const double dt=T/(NT-1);
  
  realtime_propagation::MySolver<NT,no_lam> solver("params.xml");
  solver.Set_dt(dt);
  realtime_propagation::Setup_initial_lambdas<NT,no_lam>(solver);    
  solver.run();

return EXIT_SUCCESS;
}
