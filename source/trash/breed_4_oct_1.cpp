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

/** Želimir Marojević
 */

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/function_parser.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/constraint_matrix.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_out.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_q.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>

#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <limits>

#include <deal.II/numerics/solution_transfer.h>

#include "ParameterReader.h"
#include "my_table.h"
#include "ef.h"
#include "functions_1.h"

namespace BreedSolver
{
  template <int dim>
  class MySolver;

  using namespace std;
  using namespace dealii;
  
  #include "my_solver_mpi_ortho_funcs.h"

#ifdef __variant_1__
  using namespace variant_1;
#endif
#ifdef __variant_2__
  using namespace variant_2;
#endif

  #include "CBase_no_mpi.h"

  template <int dim>
  class MySolver : public CBase<2>
  {
  public:
    MySolver( ParameterHandler & );
    ~MySolver();

    void run ();
    void run2 ();
    void run2b ();
    
    double m_T[2];
    double m_W[5];
    double m_I12; 
  private:
    void make_grid();
    void make_grid_custom();
    void setup_system();
    void assemble_system();
    void assemble_rhs();
    int DoIter(string="");
    
    void solve();
    void compute_contributions();
    void compute_E_lin( Vector<double>&, double&, double&, double& );
    void output_results ( string, string = "step" );
    void output_guess();
    void save( string );
    double Particle_Number( Vector<double>& );

    Triangulation<dim>   triangulation;
    DoFHandler<dim>      dof_handler;
    FE_Q<dim>            fe;
    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;

    Vector<double> newton_update;
    Vector<double> system_rhs;
    Vector<double> m_Psi_ref;
    Vector<double> m_Psi_0;
    Vector<double> m_Psi_1;
    Vector<double> m_workspace;
    Vector<double> m_error_per_cell;

    string m_Potential_str;
    
    MyTable m_table;
    MyTable m_results;
  };
  
/*************************************************************************************************/
  template <int dim>
  class BoundaryValues : public Function<dim>
  {
  public:
    BoundaryValues () : Function<dim>() {}
    virtual double value (const Point<dim> &p, const unsigned int  component = 0) const;
  };

  template <int dim>
  double BoundaryValues<dim>::value (const Point<dim> &p, const unsigned int /*component*/) const
  {
    return 0.0;
  }
/*************************************************************************************************/

  template <int dim>
  MySolver<dim>::MySolver ( ParameterHandler &ph ) 
    : 
    CBase<2>(ph),
    triangulation(typename Triangulation<dim>::MeshSmoothing(Triangulation<dim>::smoothing_on_refinement|Triangulation<dim>::smoothing_on_coarsening), true), 
    dof_handler (triangulation), 
    fe (1)
  {
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
    Point<dim,double> pt1( m_xmin );
    Point<dim,double> pt2( m_xmax );

    GridGenerator::hyper_rectangle(triangulation, pt2, pt1);
    triangulation.refine_global(m_global_refinement);

    double min_cell_diameter = GridTools::minimal_cell_diameter(triangulation);
    double max_cell_diameter = GridTools::maximal_cell_diameter(triangulation);
    cout << "min_cell_diameter = " << min_cell_diameter << "\n";
    cout << "max_cell_diameter = " << max_cell_diameter << "\n";
    
    m_total_no_cells = triangulation.n_cells();
    m_total_no_active_cells = triangulation.n_active_cells();
    
    
  }
  
  template <int dim>
  void MySolver<dim>::make_grid_custom ()
  {
    TimerOutput::Scope timing_section(m_computing_timer, "");

    Point<1,double> pt1( m_xmin );
    Point<1,double> pt2( m_xmax );

    GridGenerator::hyper_rectangle(triangulation, pt2, pt1);

    triangulation.refine_global(5);    
    //double isovalues[] = {28000,27000,26000,25000,24000,23000,22000,21000,20000,19000,18000,17000,16000,15000,14000};
    double isovalues[] = {28000,27000,26000,25000,24000,23000,22000,21000,20000,19000,18000,17000};
    
    for( unsigned int step=0; step<sizeof(isovalues)/sizeof(double); step++ )
    {
      typename Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active(), endc = triangulation.end();
      for( ; cell!=endc; ++cell )
        for (unsigned int v=0; v < GeometryInfo<dim>::vertices_per_cell; ++v )
        {
          Point<dim> p = cell->vertex(v);
          if( abs(p(0)) < isovalues[step] )
          {
            cell->set_refine_flag ();
            break;
          }
        }
      triangulation.execute_coarsening_and_refinement ();
    }

    m_total_no_cells = triangulation.n_cells();
    m_total_no_active_cells = triangulation.n_active_cells();
    
    double min_cell_diameter = GridTools::minimal_cell_diameter(triangulation);
    double max_cell_diameter = GridTools::maximal_cell_diameter(triangulation);
    cout << "min_cell_diameter = " << min_cell_diameter << "\n";
    cout << "max_cell_diameter = " << max_cell_diameter << "\n";
    
    
  }

  template<int dim>
  double MySolver<dim>::Particle_Number( Vector<double>& vec )
  {
    double retval=0;
    
    const QGauss<dim>  quadrature_formula(fe.degree+1);
    FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_JxW_values);

    const unsigned int n_q_points    = quadrature_formula.size();
    vector<double> vec_vals(n_q_points);

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
      fe_values.reinit (cell);
      fe_values.get_function_values( vec, vec_vals );      
      
      for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
        retval += fe_values.JxW(q_point)*vec_vals[q_point]*vec_vals[q_point];
    }
  return retval;
  }  
  
  template <int dim>
  void MySolver<dim>::compute_E_lin( Vector<double>& vec, double& T, double& N, double& W )
  {
    TimerOutput::Scope timing_section(m_computing_timer, "");
    
    map<string,double> constants;
    constants["pi"] = numbers::PI;
    FunctionParser<dim> Potential;
    Potential.initialize( FunctionParser<dim>::default_variable_names(), m_Potential_str, constants );
    
    const QGauss<dim>  quadrature_formula(fe.degree+1);
    FEValues<dim> fe_values (fe, quadrature_formula, update_gradients|update_values|update_JxW_values|update_quadrature_points);

    const unsigned int n_q_points = quadrature_formula.size();
    vector<double> vec_vals(n_q_points);
    vector<Tensor<1,dim> > vec_grad(n_q_points);

    double JxW, vec_val_q;
    double T1=0.0, N1=0.0, W1=0.0;
    
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for ( ; cell!=endc; ++cell)
    {
      fe_values.reinit (cell);
      fe_values.get_function_values(vec, vec_vals);
      fe_values.get_function_gradients(vec, vec_grad);

      for (unsigned int q_point = 0; q_point<n_q_points; ++q_point)
      {
        JxW = fe_values.JxW(q_point);
        vec_val_q = vec_vals[q_point]*vec_vals[q_point];
        T1 += JxW*( vec_grad[q_point]*vec_grad[q_point] + Potential.value(fe_values.quadrature_point(q_point))*vec_val_q);
        N1 += JxW*vec_val_q;
        W1 += JxW*vec_val_q*vec_val_q;
      }
    }
    T = T1;
    N = N1;
    W = W1;
    
  }
  
  template <int dim>
  void MySolver<dim>::compute_contributions()
  {
    TimerOutput::Scope timing_section(m_computing_timer, "");

    map<string,double> constants;
    constants["pi"] = numbers::PI;
    FunctionParser<dim> Potential;
    Potential.initialize( FunctionParser<dim>::default_variable_names(), m_Potential_str, constants );

    const QGauss<dim>  quadrature_formula(fe.degree+1);
    FEValues<dim> fe_values (fe, quadrature_formula, update_gradients|update_values|update_JxW_values|update_quadrature_points);

    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    const unsigned int n_q_points    = quadrature_formula.size();

    vector<double> Psi_1(n_q_points);
    vector<double> Psi_2(n_q_points);
    vector<Tensor<1, dim> > Psi_1_grad(n_q_points);
    vector<Tensor<1, dim> > Psi_2_grad(n_q_points);
    vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    
    double JxW, p12, p1q, p2q, Q;

    memset( m_T, 0, sizeof(m_T));
    memset( m_W, 0, sizeof(m_W));
    m_I12 = 0;
    
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
      fe_values.reinit (cell);
      fe_values.get_function_values( m_Psi_0, Psi_1 );
      fe_values.get_function_values( m_Psi_1, Psi_2 );
      fe_values.get_function_gradients(m_Psi_0, Psi_1_grad);
      fe_values.get_function_gradients(m_Psi_1, Psi_2_grad);
      
      for (unsigned int qpt = 0; qpt < n_q_points; ++qpt)
      {
        JxW = fe_values.JxW(qpt);
        p12 = Psi_1[qpt]*Psi_2[qpt];
        p1q = Psi_1[qpt]*Psi_1[qpt];
        p2q = Psi_2[qpt]*Psi_2[qpt];
	Q = Potential.value(fe_values.quadrature_point(qpt)) - m_mu;

        m_T[0] += JxW*(Psi_1_grad[qpt]*Psi_1_grad[qpt] + Q*p1q);
        m_T[1] += JxW*(Psi_2_grad[qpt]*Psi_2_grad[qpt] + Q*p2q);
        m_I12 += JxW*(Psi_1_grad[qpt]*Psi_2_grad[qpt] + Q*p12);
        m_W[0] += JxW*p1q*p1q;
        m_W[1] += JxW*p2q*p2q;
        m_W[2] += JxW*p1q*p2q;
        m_W[3] += JxW*p1q*p12;
        m_W[4] += JxW*p2q*p12;
      }  
    }
    
  }

  template <int dim>
  void MySolver<dim>::setup_system()
  {
    TimerOutput::Scope timing_section(m_computing_timer, "");

    dof_handler.distribute_dofs (fe);

    m_Psi_ref.reinit (dof_handler.n_dofs());
    m_Psi_0.reinit (dof_handler.n_dofs());
    m_Psi_1.reinit (dof_handler.n_dofs());
    newton_update.reinit (dof_handler.n_dofs());
    system_rhs.reinit (dof_handler.n_dofs());
    m_workspace.reinit(dof_handler.n_dofs());
    m_error_per_cell.reinit(triangulation.n_active_cells());

    CompressedSparsityPattern c_sparsity(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern (dof_handler, c_sparsity);

    sparsity_pattern.copy_from(c_sparsity);
    system_matrix.reinit (sparsity_pattern);
    
  }

  template <int dim>
  void MySolver<dim>::assemble_system ()
  {
    TimerOutput::Scope timing_section(m_computing_timer, "");

    vector<bool> boundary_dofs (dof_handler.n_dofs());
    DoFTools::extract_boundary_dofs (dof_handler, ComponentMask(), boundary_dofs);
    for (unsigned int i=0; i<dof_handler.n_dofs(); ++i)
      if (boundary_dofs[i] == true) m_Psi_ref(i) = 0.0;    
    
    map<string,double> constants;
    constants["pi"] = numbers::PI;
    FunctionParser<dim> Potential;
    Potential.initialize( FunctionParser<dim>::default_variable_names(), m_Potential_str, constants );

    const QGauss<dim> quadrature_formula(fe.degree+1);

    system_matrix=0;

    FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_gradients|update_JxW_values|update_quadrature_points);

    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    const unsigned int n_q_points    = quadrature_formula.size();

    FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);

    vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    vector<Tensor<1, dim> > Psi_ref_grad(n_q_points);
    vector<double> Psi_ref(n_q_points);
 
    double u2, JxW, pot, Q1, Q2;
    
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
      fe_values.reinit (cell);
      fe_values.get_function_values(m_Psi_ref, Psi_ref);
      fe_values.get_function_gradients(m_Psi_ref, Psi_ref_grad);

      cell_matrix=0;
 
      for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
      {
        JxW = fe_values.JxW(q_point);
        u2 = Psi_ref[q_point]*Psi_ref[q_point];
	Q2 = Potential.value(fe_values.quadrature_point(q_point)) - m_mu + 3.0*m_gs*u2;
        for (unsigned int i=0; i<dofs_per_cell; ++i)
        {
          for (unsigned int j=0; j<dofs_per_cell; ++j)
          {
            cell_matrix(i,j) += (fe_values.shape_grad(i,q_point)*fe_values.shape_grad(j,q_point) + Q2*fe_values.shape_value(i,q_point)*fe_values.shape_value(j,q_point))*JxW;
          }
        }
      }
      
      cell->get_dof_indices (local_dof_indices);
      for (unsigned int i=0; i<dofs_per_cell; ++i)
      {
        for (unsigned int j=0; j<dofs_per_cell; ++j) 
          system_matrix.add (local_dof_indices[i],local_dof_indices[j], cell_matrix(i,j));
      }
    }

    map<types::global_dof_index,double> boundary_values;
    VectorTools::interpolate_boundary_values (dof_handler, 0, ZeroFunction<dim>(), boundary_values);
    MatrixTools::apply_boundary_values (boundary_values, system_matrix, newton_update, system_rhs);

    
  }

  template <int dim>
  void MySolver<dim>::assemble_rhs ()
  {
    TimerOutput::Scope timing_section(m_computing_timer, "");

    vector<bool> boundary_dofs (dof_handler.n_dofs());
    DoFTools::extract_boundary_dofs (dof_handler, ComponentMask(), boundary_dofs);
    for (unsigned int i=0; i<dof_handler.n_dofs(); ++i)
      if (boundary_dofs[i] == true) m_Psi_ref(i) = 0.0;    
    
    map<string,double> constants;
    constants["pi"] = numbers::PI;
    FunctionParser<dim> Potential;
    Potential.initialize( FunctionParser<dim>::default_variable_names(), m_Potential_str, constants );
    
    const QGauss<dim> quadrature_formula(fe.degree+1);

    system_rhs=0;

    FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_gradients|update_JxW_values|update_quadrature_points);

    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    const unsigned int n_q_points    = quadrature_formula.size();

    vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    vector<Tensor<1, dim> > Psi_ref_grad(n_q_points);
    vector<double> Psi_ref(n_q_points);
 
    double u2, JxW, pot, Q1, Q2;
    
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
      fe_values.reinit (cell);
      fe_values.get_function_values(m_Psi_ref, Psi_ref);
      fe_values.get_function_gradients(m_Psi_ref, Psi_ref_grad);
      cell->get_dof_indices (local_dof_indices);

      for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
      {
        JxW = fe_values.JxW(q_point);
        u2 = Psi_ref[q_point]*Psi_ref[q_point];
        Q1 = Potential.value(fe_values.quadrature_point(q_point)) - m_mu + m_gs*u2;

	for (unsigned int i=0; i<dofs_per_cell; ++i)
        {
          system_rhs(local_dof_indices[i]) += (Psi_ref_grad[q_point]*fe_values.shape_grad(i,q_point) + Q1*Psi_ref[q_point]*fe_values.shape_value(i,q_point))*JxW;
        }
      }
    }

    map<types::global_dof_index,double> boundary_values;
    VectorTools::interpolate_boundary_values (dof_handler, 0, ZeroFunction<dim>(), boundary_values);
    MatrixTools::apply_boundary_values (boundary_values, system_matrix, newton_update, system_rhs);

    

    m_workspace = system_rhs;
    VectorTools::integrate_difference ( dof_handler,  m_workspace, ZeroFunction<dim>(), m_error_per_cell,  QGauss<dim>(fe.degree+2), VectorTools::L2_norm);
    m_res = m_error_per_cell.l2_norm();
  }
  
  template <int dim>
  void MySolver<dim>::solve ()
  {
    TimerOutput::Scope timing_section(m_computing_timer, "");

    SparseDirectUMFPACK  A_direct;
    A_direct.initialize(system_matrix);
    A_direct.vmult (newton_update, system_rhs);
    
    /*
    SolverControl solver_control (system_rhs.size(), 1e-16);
    SolverGMRES<Vector<double>> solver (solver_control);

    PreconditionSSOR<> preconditioner;
    preconditioner.initialize(system_matrix, 1.2);

    solver.solve(system_matrix, newton_update, system_rhs, preconditioner);
    */
    
  }

  template <int dim>
  int MySolver<dim>::DoIter ( string path )
  {
    int retval = Status::CONTINUE;
    
    m_table.clear();
    
    m_t[0] = m_ti;
    m_t[1] = m_ti;
    
    m_res=0;
    m_counter=0;
    m_Psi_ref=0;
    m_Psi_ref.add(m_t[0],m_Psi_0,m_t[1],m_Psi_1);
    assemble_rhs();
    m_res_old = m_res;
    do
    {
      cout << "--------------------------------------------------------------------------------" << endl;
      cout << "- " << m_counter << " - " << path << endl;

      assemble_system();
      solve();

      m_Psi_1.add( -m_t[1]/fabs(m_t[1]), newton_update); 

      if( m_counter % m_NA == 0 ) output_results(path);

#ifdef __variant_1__
      find_ortho_min();
#endif
#ifdef __variant_2__
      find_ortho_min(false);
#endif
      
      m_Psi_ref=0;
      m_Psi_ref.add(m_t[0],m_Psi_0,m_t[1],m_Psi_1);
      assemble_rhs();

      m_resp = m_res_old-m_res;
      m_res_over_resp = fabs( m_res/m_resp );
      m_res_old = m_res;

      columns& cols = m_table.new_line();
      m_table.insert( cols, MyTable::COUNTER, double(m_counter) );
      m_table.insert( cols, MyTable::RES, m_res );
      m_table.insert( cols, MyTable::RESP, m_resp );
      m_table.insert( cols, MyTable::RES_OVER_RESP, m_res_over_resp );
      m_table.insert( cols, MyTable::MU, m_mu );
      m_table.insert( cols, MyTable::GS, m_gs );
      m_table.insert( cols, MyTable::t1, m_t[0] );
      m_table.insert( cols, MyTable::t2, m_t[1] );
      m_table.insert( cols, MyTable::l2norm_t, l2norm_t() );
      m_table.insert( cols, MyTable::total_no_cells, double(m_total_no_cells) );
      m_table.insert( cols, MyTable::total_no_active_cells, double(m_total_no_active_cells) );

      cout << m_table;
      
      m_counter++;
      
      if( l2norm_t() < 1e-4 ) retval=Status::ZERO_SOL; 
      if( m_res < m_epsilon ) retval=Status::SUCCESS; 
    }while( retval==Status::CONTINUE );

    m_Psi_ref=0;
    m_Psi_ref.add(m_t[0],m_Psi_0,m_t[1],m_Psi_1);
    m_N = Particle_Number(m_Psi_ref);
    
    if( m_N < 1e-5 ) retval = Status::ZERO_SOL;
    
    string filename = path + "log.csv";
    m_table.dump_2_file(filename);
    return retval;
  }  

  /** @brief computes only one solution
   */
  template <int dim>
  void MySolver<dim>::run2 ()
  {
    int status;
    double T, N, W;

    make_grid();
    setup_system();

    m_Potential_str = "0.5*(x+3)^2";
    m_guess_str = "exp(-0.5*(x+3)^2)";    
    
    map<string,double> constants;
    constants["pi"] = numbers::PI;
    FunctionParser<dim> guess_fct;
    guess_fct.initialize( FunctionParser<dim>::default_variable_names(), m_guess_str , constants );

    VectorTools::interpolate (dof_handler, guess_fct, system_rhs );
    m_Psi_0 = system_rhs;
    m_Psi_0 *= 1.0/sqrt(Particle_Number(m_Psi_0));
    m_Psi_1 = 0;
    
    compute_E_lin( m_Psi_0, T, N, W );
    output_guess();

    m_mu = T/N+m_gs/fabs(m_gs)*m_dmu;
#ifdef NEHARI    
    m_ti = sqrt((m_mu*N-T)/(m_gs*W)); 
#endif
      
    cout << "T = " << T << endl;
    cout << "N = " << N << endl;
    cout << "W = " << W << endl;
    cout << "m_mu = " << m_mu << endl;
    cout << "m_gs = " << m_gs << endl;
    cout << "m_ti = " << m_ti << endl;

    status = DoIter();
    
    columns& cols = m_results.new_line();
    m_results.insert( cols, MyTable::MU, m_mu );
    m_results.insert( cols, MyTable::GS, m_gs );
    m_results.insert( cols, MyTable::PARTICLE_NUMBER, m_N );
    m_results.insert( cols, MyTable::COUNTER, double(m_counter) );
    m_results.insert( cols, MyTable::STATUS, double(status) );    

    cout << "Particle number: " << m_N << endl;

    if( status == Status::SUCCESS )
    {
      output_results("","Psi_0");
      save("Psi_0.bin");
    }

    ofstream ofs("log.txt");
    //ofs << m_table;
    
    m_Potential_str = "0.5*(x-3)^2";
    m_guess_str = "exp(-0.5*(x-3)^2)";    
    
    guess_fct.initialize( FunctionParser<dim>::default_variable_names(), m_guess_str , constants );

    VectorTools::interpolate (dof_handler, guess_fct, system_rhs );
    m_Psi_0 = system_rhs;
    m_Psi_0 *= 1.0/sqrt(Particle_Number(m_Psi_0));
    m_Psi_1 = 0;
    
    compute_E_lin( m_Psi_0, T, N, W );
    output_guess();

    m_mu = T/N+m_gs/fabs(m_gs)*m_dmu;
#ifdef NEHARI    
    m_ti = sqrt((m_mu*N-T)/(m_gs*W)); 
#endif
      
    cout << "T = " << T << endl;
    cout << "N = " << N << endl;
    cout << "W = " << W << endl;
    cout << "m_mu = " << m_mu << endl;
    cout << "m_gs = " << m_gs << endl;
    cout << "m_ti = " << m_ti << endl;

    status = DoIter();
    
    cols = m_results.new_line();
    m_results.insert( cols, MyTable::MU, m_mu );
    m_results.insert( cols, MyTable::GS, m_gs );
    m_results.insert( cols, MyTable::PARTICLE_NUMBER, m_N );
    m_results.insert( cols, MyTable::COUNTER, double(m_counter) );
    m_results.insert( cols, MyTable::STATUS, double(status) );    

    cout << "Particle number: " << m_N << endl;

    if( status == Status::SUCCESS )
    {
      output_results("","Psi_d");
      save("Psi_d.bin");
    }

    ofstream ofs2("log2.txt");
    ofs2 << m_table;

    m_results.dump_2_file( "results.csv" );
  }

  template <int dim>
  void MySolver<dim>::output_guess ()
  {
    TimerOutput::Scope timing_section(m_computing_timer, "");
    string filename = "guess.gnuplot";

    map<string,double> constants;
    constants["pi"] = numbers::PI;
    FunctionParser<dim> Potential;
    Potential.initialize( FunctionParser<dim>::default_variable_names(), m_Potential_str, constants );

    VectorTools::interpolate (dof_handler, Potential, m_workspace );
    
    DataOut<dim> data_out;
    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (m_Psi_0, "Psi_0");
    data_out.add_data_vector (m_Psi_1, "Psi_1");
    data_out.add_data_vector (m_workspace, "Potential");
    data_out.build_patches ();

    ofstream output (filename.c_str());
    data_out.write_gnuplot (output);
    
  }
  
  template <int dim>
  void MySolver<dim>::output_results ( string path, string prefix )
  {
    TimerOutput::Scope timing_section(m_computing_timer, "");
    string filename = path + prefix + "-" + Utilities::int_to_string (m_counter,5) + ".gnuplot";

    DataOut<dim> data_out;
    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (m_Psi_ref, "Psi_ref");
    data_out.add_data_vector (m_error_per_cell, "error_per_cell");
    data_out.build_patches ();

    ofstream output (filename.c_str());
    data_out.write_gnuplot (output);
    
  }
  
  template<int dim>
  void MySolver<dim>::save( string filename )
  {
    Vector<double> tmp( 2*dof_handler.n_dofs() );
    tmp=0;
    for( int i=0; i<dof_handler.n_dofs(); ++i )
    {
      tmp[i] = m_Psi_ref[i]; 
    }
    ofstream out(filename);
    tmp.block_write(out);
  }
} // end of namespace 

int main ( int argc, char *argv[] )
{
  using namespace dealii;
  using namespace BreedSolver;
  deallog.depth_console (0);
  
  ParameterHandler  prm;
  ParameterReader   param(prm);
  
  if( !param.read_parameters() )
    param.print_parameters();
    
  MySolver<1> solver(prm);
  solver.run2(); 
return EXIT_SUCCESS;
}
