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
 * @file rtprop.cpp Real time propagation for the Gross--Pitaevskii equation (Cartesian coordinates)
 * @brief Real time propagation for the Gross--Pitaevskii equation (Cartesian coordinates)
 * @author Želimir Marojević
 */

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/timer.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparsity_tools.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
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
#include <deal.II/base/function_parser.h>

#include <fstream>
#include <iostream>
#include <stdlib.h>

#include "global.h"
#include "MyParameterHandler.h"
#include "MyComplexTools.h"
#include "my_table.h"
#include "cxxopts.hpp"

///namespace realtime_propagation
namespace realtime_propagation
{
  using namespace std;
  using namespace dealii;

  class MySolver
  {
  public:
    MySolver( const std::string& );
    ~MySolver();

    void run ( const std::string& );

  protected:
    void make_grid();
    void setup_system();
    void assemble_system();
    void assemble_rhs();
    
    void DoIter();
    void solve();
    void output_results ( string );
    void load( string );
    void save( string );    

    MyParameterHandler m_ph;
    Triangulation<1> triangulation;
    FESystem<1> fe;
    DoFHandler<1> dof_handler;

    SparsityPattern sparsity_pattern;
    SparseMatrix<double> system_matrix;

    Vector<double> system_rhs;
    Vector<double> newton_update;
    Vector<double> m_Psi; // Psi(t)
    Vector<double> m_Psi_t; // Psi(t)
    Vector<double> m_workspace;
    Vector<double> m_error_per_cell;

    double m_gs;
    double m_t;
    double m_dt;
    double m_res;
    double m_xmin;
    double m_xmax;

    /// total number of outputs
    unsigned m_NA;
    /// number of intermediate time steps between two outputs
    unsigned m_NK;
    
    unsigned m_global_refinement;
    
    MyTable m_table;
    string m_filename;
  };

  /** Default constructor
   */
  MySolver::MySolver ( const std::string& xmlfilename ) 
    : 
    m_ph(xmlfilename),
    triangulation (),
    fe (FE_Q<1>(gl_degree_fe), 2),
    dof_handler (triangulation)
  {
    try
    {
      m_gs = m_ph.Get_Physics("gs_1",0);
      
      m_xmin = m_ph.Get_Mesh("xrange",0);
      m_xmax = m_ph.Get_Mesh("xrange",1);
      m_global_refinement = unsigned(m_ph.Get_Mesh("global_refinements",0));

      m_NA = m_ph.Get_Algorithm("NA",0);
      m_NK = m_ph.Get_Algorithm("NK",0);
      m_dt = m_ph.Get_Algorithm("dt",0);
    }
    catch( const std::string info )
    {
      std::cerr << info << endl;
      exit(0);
    }    
    m_t = 0;
  }

  /** Default destructor
   */
  MySolver::~MySolver ()
  {
    dof_handler.clear ();
  }
  
  #include "utils_complex.h"
 
  void MySolver::make_grid ()
  {
    Point<1,double> pt1( m_xmin );
    Point<1,double> pt2( m_xmax );

    GridGenerator::hyper_rectangle(triangulation, pt2, pt1);
    triangulation.refine_global(m_global_refinement);
  }
  
  void MySolver::setup_system()
  {
    dof_handler.distribute_dofs (fe);
    //DoFRenumbering::component_wise (dof_handler);

    m_Psi.reinit (dof_handler.n_dofs());
    m_workspace.reinit (dof_handler.n_dofs());
    m_Psi_t.reinit (dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());
    newton_update.reinit(dof_handler.n_dofs());
    m_error_per_cell.reinit(dof_handler.n_dofs());

    DynamicSparsityPattern dsp(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern (dof_handler, dsp);
    sparsity_pattern.copy_from(dsp);
    system_matrix.reinit (sparsity_pattern);
  }  

  /** Output of data to vtu files
   * @param[in] path specifies the output folder 
   */
  void MySolver::output_results ( string path ) 
  {
    ComputeIntensity<1> intensities;
    ComputePhase<1> phase;
    
    vector<std::string> solution_names;

    DataOut<1> data_out;
    data_out.attach_dof_handler (dof_handler);
    solution_names.push_back ("Re Psi");
    solution_names.push_back ("Im Psi");    
    data_out.add_data_vector (m_Psi, solution_names);
    data_out.add_data_vector (m_Psi , intensities);
    data_out.build_patches ();
    
    string filename = path + "solution-" + to_string(m_t) + ".gnuplot";
    ofstream output (filename);
    data_out.write_gnuplot ( output );
  }   

  void MySolver::save( string filename )
  {
  }

  void MySolver::load( string filename )
  {
    ifstream in(filename);
    m_Psi.block_read(in);
  }    

  void MySolver::assemble_system ()
  {
    const QGauss<1> quadrature_formula(fe.degree+1);

    FunctionParser<1> Potential;
    Potential.initialize( FunctionParser<1>::default_variable_names(), m_ph.Get_String( "POTENTIAL", 0 ), m_ph.Get_Constants_Map());

    const FEValuesExtractors::Scalar rt (0);
    const FEValuesExtractors::Scalar it (1);

    system_matrix=0;

    FEValues<1> fe_values (fe, quadrature_formula, update_values|update_gradients|update_quadrature_points|update_JxW_values);

    const unsigned dofs_per_cell = fe.dofs_per_cell;
    const unsigned n_q_points    = quadrature_formula.size();

    FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);

    vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    vector<Vector<double>> Psi(n_q_points,Vector<double>(2));
    vector<vector<Tensor<1,1>>> Psi_grad(n_q_points, vector<Tensor<1,1>>(2));
    vector<Vector<double>> Psi_t(n_q_points,Vector<double>(2));
    vector<vector<Tensor<1,1>>> Psi_t_grad(n_q_points, vector<Tensor<1,1>>(2));
 
    double JxW, pot=0, tmp1a, tmp1b, tmp2, sum_re, sum_im, sum_req, sum_imq;

    const double fak2 = 0.5*m_dt;
    const double fak4 = 0.25*m_gs*m_dt;
    const double fak8 = 0.125*m_gs*m_dt;

    DoFHandler<1>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for( ; cell!=endc; ++cell )
    {
      cell_matrix=0;

      fe_values.reinit (cell);
      fe_values.get_function_values(m_Psi, Psi);
      fe_values.get_function_gradients(m_Psi, Psi_grad);
      fe_values.get_function_values(m_Psi_t, Psi_t);
      fe_values.get_function_gradients(m_Psi_t, Psi_t_grad);

      for( unsigned qp=0; qp<n_q_points; ++qp )
      {
        JxW = fe_values.JxW(qp);
        pot = Potential.value(fe_values.quadrature_point(qp)); 

        sum_re = Psi[qp][0]+Psi_t[qp][0];
        sum_im = Psi[qp][1]+Psi_t[qp][1];
        sum_req = sum_re*sum_re;
        sum_imq = sum_im*sum_im;
        tmp1a = fak8*(sum_req + 3*sum_imq);
        tmp1b = fak8*(sum_imq + 3*sum_req);
        tmp2 = fak4*sum_re*sum_im;
        
        for( unsigned i=0; i<dofs_per_cell; ++i )
        {
          for( unsigned j=0; j<dofs_per_cell; ++j )
          {
            cell_matrix(i,j) += JxW*((1.0-tmp2)*fe_values[rt].value(i,qp)*fe_values[rt].value(j,qp) + (1.0+tmp2)*fe_values[it].value(i,qp)*fe_values[it].value(j,qp)
                                     - tmp1a*fe_values[rt].value(i,qp)*fe_values[it].value(j,qp) - fak2*(fe_values[rt].gradient(i,qp)*fe_values[it].gradient(j,qp) + pot*fe_values[rt].value(i,qp)*fe_values[it].value(j,qp))
                                     + tmp1b*fe_values[it].value(i,qp)*fe_values[rt].value(j,qp) + fak2*(fe_values[it].gradient(i,qp)*fe_values[rt].gradient(j,qp) + pot*fe_values[it].value(i,qp)*fe_values[rt].value(j,qp)));
          }
        }
      }      
      cell->get_dof_indices (local_dof_indices);
      for ( unsigned i = 0; i < dofs_per_cell; ++i )
        for ( unsigned j = 0; j < dofs_per_cell; ++j )
          system_matrix.add (local_dof_indices[i], local_dof_indices[j], cell_matrix(i, j));
    }
  }
  
  void MySolver::assemble_rhs ()
  {
    const QGauss<1> quadrature_formula(fe.degree+1);

    FunctionParser<1> Potential;
    Potential.initialize( FunctionParser<1>::default_variable_names(), m_ph.Get_String( "POTENTIAL", 0 ), m_ph.Get_Constants_Map());

    const FEValuesExtractors::Scalar rt (0);
    const FEValuesExtractors::Scalar it (1);

    system_rhs=0;

    FEValues<1> fe_values (fe, quadrature_formula, update_values|update_gradients|update_quadrature_points|update_JxW_values);

    const unsigned dofs_per_cell = fe.dofs_per_cell;
    const unsigned n_q_points = quadrature_formula.size();

    Vector<double> cell_rhs (dofs_per_cell);
    vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    vector<Vector<double>> Psi(n_q_points,Vector<double>(2));
    vector<vector<Tensor<1,1>>> Psi_grad(n_q_points, vector<Tensor<1,1>>(2));
    vector<Vector<double>> Psi_t(n_q_points,Vector<double>(2));
    vector<vector<Tensor<1,1>>> Psi_t_grad(n_q_points, vector<Tensor<1,1>>(2));
 
    const double fak2 = 0.5*m_dt;
    const double fak8 = 0.125*m_gs*m_dt;

    DoFHandler<1>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
      cell_rhs = 0;
      fe_values.reinit (cell);
      fe_values.get_function_values(m_Psi, Psi);
      fe_values.get_function_gradients(m_Psi, Psi_grad);
      fe_values.get_function_values(m_Psi_t, Psi_t);
      fe_values.get_function_gradients(m_Psi_t, Psi_t_grad);

      cell->get_dof_indices (local_dof_indices);
      for( unsigned qp=0; qp<n_q_points; ++qp )
      {
        double JxW = fe_values.JxW(qp);
        double pot = Potential.value(fe_values.quadrature_point(qp)); 
        double sum_re = Psi[qp][0]+Psi_t[qp][0];
        double sum_im = Psi[qp][1]+Psi_t[qp][1];
        double tmp1 = fak8*(sum_re*sum_re + sum_im*sum_im);

        for( unsigned i=0; i<dofs_per_cell; ++i )
        {
          cell_rhs(i) += JxW*(-fak2*((Psi_grad[qp][1]+Psi_t_grad[qp][1])*fe_values[rt].gradient(i,qp) + pot*sum_im*fe_values[rt].value(i,qp))
                              +fak2*((Psi_grad[qp][0]+Psi_t_grad[qp][0])*fe_values[it].gradient(i,qp) + pot*sum_re*fe_values[it].value(i,qp))
                              +(Psi_t[qp][0]-Psi[qp][0])*fe_values[rt].value(i,qp) - tmp1*sum_im*fe_values[rt].value(i,qp)  
                              +(Psi_t[qp][1]-Psi[qp][1])*fe_values[it].value(i,qp) + tmp1*sum_re*fe_values[it].value(i,qp));
        }
      }
      for ( unsigned i = 0; i < dofs_per_cell; ++i )
        system_rhs(local_dof_indices[i]) += cell_rhs(i);
    }
    m_res = system_rhs.l2_norm(); 
  }
  
  void MySolver::solve ()
  {
    map<types::global_dof_index, double> boundary_values;
    VectorTools::interpolate_boundary_values (dof_handler, 0, ZeroFunction<1>(), boundary_values);
    MatrixTools::apply_boundary_values (boundary_values, system_matrix, newton_update, system_rhs);    

    SparseDirectUMFPACK  A_direct;
    A_direct.initialize(system_matrix);
    A_direct.vmult (newton_update, system_rhs);
  }

  void MySolver::DoIter()
  {
    m_Psi_t = 0; 
    m_res = 0;
    assemble_rhs();
    //cout << "m_res = " << m_res << endl;       
    do
    {
      assemble_system();
      solve();

      m_Psi_t.add( -1, newton_update );
      
      assemble_rhs();
      //cout << "m_res = " << m_res << endl;       
    }
    while( m_res > 1e-16 ); 
    m_t += m_dt;

    m_Psi = m_Psi_t;
  }

  void MySolver::run( const std::string& filename )
  {
    make_grid();
    setup_system();

    load( filename );

    double min_cell_diameter = GridTools::minimal_cell_diameter(triangulation);
    double max_cell_diameter = GridTools::maximal_cell_diameter(triangulation);
    cout << "min_cell_diameter = " << min_cell_diameter << "\n";
    cout << "max_cell_diameter = " << max_cell_diameter << "\n";
    cout << "dt/dx^2 == " << m_dt/(min_cell_diameter*min_cell_diameter) << endl;
    
    vector<double> p(1,0);
    vector<double> pos(1,0);

    output_results("");

    double N = MyComplexTools::Particle_Number( dof_handler, fe, m_Psi );
    //cout << "N == " << N << endl;
    MyComplexTools::Expectation_value_position( dof_handler, fe, m_Psi, pos );
    MyComplexTools::Expectation_value_momentum( dof_handler, fe, m_Psi, p );

    cout << "t == " << m_t << endl;
    cout << "N == " << N << endl;
    cout << "p == " << p[0]/N << endl;
    cout << "pos == " << pos[0]/N << endl;

    for( unsigned i=1; i<=m_NA; ++i )
    {
      for( unsigned j=1; j<=m_NK; ++j )
      {
        cout << "t == " << m_t << endl;
        DoIter();
      }
      
      N = MyComplexTools::Particle_Number( dof_handler, fe, m_Psi );
      MyComplexTools::Expectation_value_position( dof_handler, fe, m_Psi, pos );
      MyComplexTools::Expectation_value_momentum( dof_handler, fe, m_Psi, p );

      cout << "N == " << N << endl;
      cout << "p == " << p[0]/N << endl;
      cout << "pos == " << pos[0]/N << endl;

      output_results("");
      save( "solution-" + to_string(m_t) + ".bin" );
    }
  }
} // end of namespace 

int main ( int argc, char *argv[] )
{
  using namespace dealii;
  deallog.depth_console (0);

  cxxopts::Options options("rt_prop_1", "real time propagation 1D");
  
  options.add_options()
  ("i,input", "input initial wave function" , cxxopts::value<std::string>()->default_value("Cfinal.bin") )
  ("p,params", "input parameter xml file" , cxxopts::value<std::string>()->default_value("params.xml") )
  ("help","Print help")
  ;
  
  auto result = options.parse(argc, argv);

  if (result.count("") == 0)
  {
    std::cout << options.help({""}) << std::endl;
    return EXIT_FAILURE;
  }

  std::string bin_filename, params_filename;
  try
  {
    bin_filename = result["i"].as<std::string>(); 
    params_filename = result["p"].as<std::string>();
  }
  catch (const cxxopts::OptionException& e)
  {
    std::cout << "error parsing options: " << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  realtime_propagation::MySolver solver(params_filename);
  solver.run(bin_filename);
return EXIT_SUCCESS;
}
