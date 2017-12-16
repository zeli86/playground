/*
    This file is part of atus-pro package.

    atus-pro is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    atus-pro is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with atus-pro.  If not, see <http://www.gnu.org/licenses/>.
*/
/*
    Želimir Marojević, Ertan Göklü, Claus Lämmerzahl
    (ZARM - Center of Applied Space Technology and Microgravity, Germany, http://www.zarm.uni-bremen.de/)

    Public use and modification of this code are allowed provided that the following papers are cited:

    Marojević, Želimir, Ertan Göklü, und Claus Lämmerzahl. "ATUS-PRO: A FEM-based solver for the time-dependent and stationary Gross–Pitaevskii equation",
    Computer Physics Communications, Vol 202, 2016, p. 216--232. doi:10.1016/j.cpc.2015.12.004.

    W. Bangerth and D. Davydov and T. Heister and L. Heltai and G. Kanschat and M. Kronbichler and M. Maier and B. Turcksin and D. Wells.
    "The \texttt{deal.II} Library, Version 8.4", Journal of Numerical Mathematics, vol 24, 2016.

    The authors would be grateful for all information and/or comments regarding the use of the code.
*/

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>
//include <deal.II/base/function.h>
//#include <deal.II/base/function_parser.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_q.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>

#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>

#include "global.h"
#include "my_table.h"
#include "functions.h"
#include "MyParameterHandler.h"
#include "MyRealTools.h"
#include "muParser.h"

/// namespace BreedSolver_1
namespace BreedSolver_1
{
  using namespace std;
  using namespace dealii;

  #include "CBase_no_mpi.h"

  template <int dim>
  class MySolver : public CBase<dim,2>
  {
  public:
    MySolver( const std::string & );
    ~MySolver();

    void run ();

    double m_I[8];
  private:
    void make_grid();
    void setup_system();
    void assemble_system();
    void assemble_rhs();
    int DoIter(string = "");

    void solve();
    void compute_contributions();
    void compute_E_lin( Vector<double> &, double &, double &, double & );
    void output_results ( string, string = "step" );
    void output_guess();
    void save( string );
    void save_one( string );
    void estimate_error( double & );

    Triangulation<dim> triangulation;
    DoFHandler<dim> dof_handler;
    FE_Q<dim> fe;

    SparsityPattern m_sparsity_pattern;
    SparseMatrix<double> m_system_matrix;

    Vector<double> m_newton_update;
    Vector<double> m_system_rhs;
    Vector<double> m_Psi_ref;
    Vector<double> m_Psi_0;
    Vector<double> m_Psi_1;
    Vector<double> m_workspace;
    Vector<double> m_error_per_cell;

    MyTable m_table;
    MyTable m_results;

    using CBase<dim,2>::m_xmin;
    using CBase<dim,2>::m_xmax;
    using CBase<dim,2>::m_ymin;
    using CBase<dim,2>::m_ymax;
    using CBase<dim,2>::m_t;
    using CBase<dim,2>::m_ti;
    using CBase<dim,2>::m_N;
    using CBase<dim,2>::m_omega;
    using CBase<dim,2>::m_mu;
    using CBase<dim,2>::m_gs;
    using CBase<dim,2>::m_counter;
    using CBase<dim,2>::m_ph;
    using CBase<dim,2>::m_final_error;
    using CBase<dim,2>::m_NA;
    using CBase<dim,2>::m_Ndmu;
    using CBase<dim,2>::m_dmu;
    using CBase<dim,2>::m_QN1;
    using CBase<dim,2>::m_res;
    using CBase<dim,2>::m_resp;
    using CBase<dim,2>::m_res_old;
    using CBase<dim,2>::m_epsilon;
    using CBase<dim,2>::m_global_refinement;
    using CBase<dim,2>::m_total_no_cells;
    using CBase<dim,2>::m_total_no_active_cells;
  };

  template <int dim>
  MySolver<dim>::MySolver ( const std::string &xml_filename )
    :
    CBase<dim,2>(xml_filename),
    triangulation(),
    dof_handler (triangulation),
    fe (gl_degree_fe)
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
    Point<dim, double> pt1, pt2;

    double min[] = {m_xmin, m_ymin};
    double max[] = {m_xmax, m_ymax};

    for ( int i = 0; i < dim; i++ )
    {
      pt1(i) = min[i];
      pt2(i) = max[i];
    }

    GridGenerator::hyper_rectangle(triangulation, pt2, pt1);
    triangulation.refine_global(m_global_refinement);

    m_total_no_cells = triangulation.n_cells();
    m_total_no_active_cells = triangulation.n_active_cells();
  }

  template <int dim>
  void MySolver<dim>::compute_E_lin( Vector<double> &vec, double &T, double &N, double &W )
  {
    CPotential<dim> Potential(m_omega);
    const QGauss<dim>  quadrature_formula(fe.degree + 1);
    FEValues<dim> fe_values (fe, quadrature_formula, update_gradients | update_values | update_JxW_values | update_quadrature_points);

    const unsigned n_q_points = quadrature_formula.size();
    vector<double> vec_vals(n_q_points);
    vector<Tensor<1, dim>> vec_grad(n_q_points);

    T=0, N=0, W=0;

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for ( ; cell != endc; ++cell)
    {
      fe_values.reinit (cell);
      fe_values.get_function_values(vec, vec_vals);
      fe_values.get_function_gradients(vec, vec_grad);

      for ( unsigned qp = 0; qp < n_q_points; qp++ )
      {
        double JxW = fe_values.JxW(qp);
        double vec_val_q = vec_vals[qp] * vec_vals[qp];
        T += JxW * ( vec_grad[qp] * vec_grad[qp] + Potential.value(fe_values.quadrature_point(qp)) * vec_val_q);
        N += JxW * vec_val_q;
        W += JxW * vec_val_q * vec_val_q;
      }
    }
  }

  template <int dim>
  void MySolver<dim>::compute_contributions()
  {
    CPotential<dim> Potential(m_omega);
    const QGauss<dim>  quadrature_formula(fe.degree + 1);
    FEValues<dim> fe_values (fe, quadrature_formula, update_gradients | update_values | update_JxW_values | update_quadrature_points);

    const unsigned dofs_per_cell = fe.dofs_per_cell;
    const unsigned n_q_points    = quadrature_formula.size();

    vector<double> Psi_1(n_q_points);
    vector<double> Psi_2(n_q_points);
    vector<Tensor<1, dim>> Psi_1_grad(n_q_points);
    vector<Tensor<1, dim>> Psi_2_grad(n_q_points);
    vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    for( int i=0; i<8; i++ ) m_I[i] = 0;

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell != endc; cell++ )
    {
      fe_values.reinit (cell);
      fe_values.get_function_values( m_Psi_0, Psi_1 );
      fe_values.get_function_values( m_Psi_1, Psi_2 );
      fe_values.get_function_gradients(m_Psi_0, Psi_1_grad);
      fe_values.get_function_gradients(m_Psi_1, Psi_2_grad);

      for ( unsigned qp = 0; qp < n_q_points; qp++ )
      {
        double JxW = fe_values.JxW(qp);
        double p12 = Psi_1[qp] * Psi_2[qp];
        double p1q = Psi_1[qp] * Psi_1[qp];
        double p2q = Psi_2[qp] * Psi_2[qp];
        double Q = Potential.value(fe_values.quadrature_point(qp)) - m_mu;

        m_I[0] += JxW*p1q*p1q;
        m_I[1] += JxW*p2q*p2q;
        m_I[2] += JxW*p1q*p2q;
        m_I[3] += JxW*p1q*p12;
        m_I[4] += JxW*p2q*p12;
        m_I[5] += JxW*(Psi_1_grad[qp]*Psi_1_grad[qp] + Q*p1q);
        m_I[6] += JxW*(Psi_2_grad[qp]*Psi_2_grad[qp] + Q*p2q);
        m_I[7] += JxW*(Psi_1_grad[qp]*Psi_2_grad[qp] + Q*p12);
      }
    }
    for( int i=0; i<5; i++ ) m_I[i] *= m_gs;
  }

  template <int dim>
  void MySolver<dim>::setup_system()
  {
    // real
    dof_handler.distribute_dofs (fe);

    m_Psi_ref.reinit (dof_handler.n_dofs());
    m_Psi_0.reinit (dof_handler.n_dofs());
    m_Psi_1.reinit (dof_handler.n_dofs());
    m_newton_update.reinit (dof_handler.n_dofs());
    m_system_rhs.reinit (dof_handler.n_dofs());
    m_workspace.reinit(dof_handler.n_dofs());
    m_error_per_cell.reinit(triangulation.n_active_cells());

    DynamicSparsityPattern dsp(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern (dof_handler, dsp);
    m_sparsity_pattern.copy_from(dsp);
    m_system_matrix.reinit (m_sparsity_pattern);
  }

  template <int dim>
  void MySolver<dim>::assemble_system ()
  {
    CPotential<dim> Potential(m_omega);
    const QGauss<dim> quadrature_formula(fe.degree + 1);

    m_system_matrix = 0;

    FEValues<dim> fe_values (fe, quadrature_formula, update_values | update_gradients | update_JxW_values | update_quadrature_points);

    const unsigned dofs_per_cell = fe.dofs_per_cell;
    const unsigned n_q_points = quadrature_formula.size();

    FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);
    vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    vector<Tensor<1, dim>> Psi_ref_grad(n_q_points);
    vector<double> Psi_ref(n_q_points);

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell != endc; ++cell)
    {
      cell_matrix = 0;

      fe_values.reinit (cell);
      fe_values.get_function_values(m_Psi_ref, Psi_ref);
      fe_values.get_function_gradients(m_Psi_ref, Psi_ref_grad);

      for ( unsigned qp = 0; qp < n_q_points; qp++ )
      {
        double JxW = fe_values.JxW(qp);
        double Q2 = Potential.value(fe_values.quadrature_point(qp)) - m_mu + 3.0 * m_gs * Psi_ref[qp] * Psi_ref[qp];

        for ( unsigned i = 0; i < dofs_per_cell; i++ )
          for ( unsigned j = 0; j < dofs_per_cell; j++ )
            cell_matrix(i, j) += JxW * (fe_values.shape_grad(i, qp) * fe_values.shape_grad(j, qp) + Q2 * fe_values.shape_value(i, qp) * fe_values.shape_value(j, qp));
      }
      cell->get_dof_indices (local_dof_indices);
      for ( unsigned i = 0; i < dofs_per_cell; i++ )
        for ( unsigned j = 0; j < dofs_per_cell; j++ )
          m_system_matrix.add (local_dof_indices[i], local_dof_indices[j], cell_matrix(i, j));
    }
  }

  template <int dim>
  void MySolver<dim>::assemble_rhs ()
  {
    /*vector<bool> boundary_dofs (dof_handler.n_dofs());
    DoFTools::extract_boundary_dofs (dof_handler, ComponentMask(), boundary_dofs);
    for (unsigned int i = 0; i < dof_handler.n_dofs(); ++i)
      if (boundary_dofs[i] == true) m_Psi_ref(i) = 0.0;*/

    CPotential<dim> Potential(m_omega);
    const QGauss<dim> quadrature_formula(fe.degree + 1);

    m_system_rhs = 0;

    FEValues<dim> fe_values (fe, quadrature_formula, update_values | update_gradients | update_JxW_values | update_quadrature_points);

    const unsigned dofs_per_cell = fe.dofs_per_cell;
    const unsigned n_q_points = quadrature_formula.size();

    Vector<double> cell_rhs (dofs_per_cell);
    vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    vector<Tensor<1, dim>> Psi_ref_grad(n_q_points);
    vector<double> Psi_ref(n_q_points);

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell != endc; cell++)
    {
      cell_rhs=0;
      fe_values.reinit (cell);
      fe_values.get_function_values(m_Psi_ref, Psi_ref);
      fe_values.get_function_gradients(m_Psi_ref, Psi_ref_grad);

      for ( unsigned qp = 0; qp < n_q_points; qp++ )
      {
        double JxW = fe_values.JxW(qp);
        double Q1 = Potential.value(fe_values.quadrature_point(qp)) - m_mu + m_gs * Psi_ref[qp] * Psi_ref[qp];

        for ( unsigned i = 0; i < dofs_per_cell; i++ )
          cell_rhs(i) += JxW * (Psi_ref_grad[qp] * fe_values.shape_grad(i, qp) + Q1 * Psi_ref[qp] * fe_values.shape_value(i, qp));
      }
      cell->get_dof_indices (local_dof_indices);
      for ( unsigned i = 0; i < dofs_per_cell; i++ )
        m_system_rhs(local_dof_indices[i]) += cell_rhs(i);
    }
  }

  template <int dim>
  void MySolver<dim>::solve ()
  {
    map<types::global_dof_index, double> boundary_values;
    VectorTools::interpolate_boundary_values (dof_handler, 0, ZeroFunction<dim>(), boundary_values);
    MatrixTools::apply_boundary_values (boundary_values, m_system_matrix, m_newton_update, m_system_rhs);    

    SparseDirectUMFPACK A_direct;
    A_direct.initialize(m_system_matrix);
    A_direct.vmult (m_newton_update, m_system_rhs);
    m_res = m_newton_update.l2_norm();
  }

  template <int dim>
  int MySolver<dim>::DoIter ( string path )
  {
    CPotential<dim> Potential( m_omega );

    int retval = Status::CONTINUE;

    m_table.clear();

    m_t[0] = m_ti;
    m_t[1] = m_ti;

    m_res = 0;
    m_counter = 0;
    m_Psi_ref = 0;
    m_Psi_ref.add(m_t[0], m_Psi_0, m_t[1], m_Psi_1);
    assemble_rhs();
    m_res_old = m_res;
    do
    {
      cout << "--------------------------------------------------------------------------------" << endl;
      cout << "- " << m_counter << " - " << path << endl;

      assemble_system();
      solve();

      double tau;
      MyRealTools::compute_stepsize( dof_handler, fe, Potential, m_Psi_ref, m_newton_update, m_mu, m_gs, tau );
      m_Psi_1.add( tau*m_t[1] / fabs(m_t[1]), m_newton_update );

      if ( m_counter % m_NA == 0 ) output_results(path);

      this->find_ortho_min();

      m_Psi_ref = 0;
      m_Psi_ref.add(m_t[0], m_Psi_0, m_t[1], m_Psi_1);
      assemble_rhs();

      m_resp = m_res_old - m_res;
      m_res_old = m_res;

      columns &cols = m_table.new_line();
      m_table.insert( cols, MyTable::COUNTER, double(m_counter) );
      m_table.insert( cols, MyTable::RES, m_res );
      m_table.insert( cols, MyTable::RESP, m_resp );
      m_table.insert( cols, MyTable::MU, m_mu );
      m_table.insert( cols, MyTable::GS, m_gs );
      m_table.insert( cols, MyTable::t1, m_t[0] );
      m_table.insert( cols, MyTable::t2, m_t[1] );
      m_table.insert( cols, MyTable::l2norm_t, this->l2norm_t() );

      cout << m_table;

      m_counter++;

      if ( this->l2norm_t() < 1e-4 ) retval = Status::ZERO_SOL;
      if ( m_res < m_epsilon[0] ) retval = Status::SUCCESS;
    }
    while ( retval == Status::CONTINUE );

    m_Psi_ref = 0;
    m_Psi_ref.add(m_t[0], m_Psi_0, m_t[1], m_Psi_1);

    // Standard Newton
    do
    {
      cout << "--------------------------------------------------------------------------------" << endl;
      cout << "-- " << path << " - " << m_counter << endl;

      assemble_system();
      solve();

      m_Psi_ref.add( -1, m_newton_update );
      assemble_rhs();

      m_resp = m_res_old - m_res;
      m_res_old = m_res;

      if ( m_counter % m_NA == 0 ) output_results(path);

      columns &cols = m_table.new_line();
      m_table.insert( cols, MyTable::COUNTER, double(m_counter) );
      m_table.insert( cols, MyTable::RES, m_res );
      m_table.insert( cols, MyTable::RESP, m_resp );
      m_table.insert( cols, MyTable::MU, m_mu );
      m_table.insert( cols, MyTable::GS, m_gs );
      m_table.insert( cols, MyTable::t1, m_t[0] );
      m_table.insert( cols, MyTable::t2, m_t[1] );
      m_table.insert( cols, MyTable::l2norm_t, this->l2norm_t() );
      m_table.insert( cols, MyTable::PARTICLE_NUMBER, m_N );

      cout << m_table;
      if ( m_res < m_epsilon[1] )
      {
        retval = Status::SUCCESS;
        break;
      }

      m_counter++;
    }
    while ( true );

    m_N = MyRealTools::Particle_Number(dof_handler,fe,m_Psi_ref);

    if ( m_N < 1e-5 ) retval = Status::ZERO_SOL;

    string filename = path + "log.csv";
    m_table.dump_2_file(filename);
    return retval;
  }

  /**
  * estimate_error ( double& err ) interpolates the L_2 gradient onto the FE space and computes
  * the L_2 norm of the interpolated gradient.
  */
  template <int dim>
  void MySolver<dim>::estimate_error ( double &err )
  {
    CPotential<dim> Potential( m_omega );
    const QGauss<1> quadrature_formula(fe.degree + 1);

    m_system_rhs = 0;
    m_system_matrix = 0;

    FEValues<1> fe_values (fe, quadrature_formula, update_values | update_gradients | update_JxW_values | update_quadrature_points);

    const unsigned dofs_per_cell = fe.dofs_per_cell;
    const unsigned n_q_points = quadrature_formula.size();

    Vector<double> cell_rhs (dofs_per_cell);
    FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);

    vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    vector<double> vals(n_q_points);
    vector<Tensor<1, 1>> grads(n_q_points);

    DoFHandler<1>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for ( ; cell != endc; ++cell )
    {
      cell_rhs = 0;
      cell_matrix = 0;

      fe_values.reinit (cell);
      cell->get_dof_indices (local_dof_indices);

      fe_values.get_function_values(m_Psi_ref, vals);
      fe_values.get_function_gradients(m_Psi_ref, grads);

      for ( unsigned qp = 0; qp < n_q_points; qp++ )
      {
        double JxW = fe_values.JxW(qp);
        double Q1 = Potential.value(fe_values.quadrature_point(qp)) - m_mu + m_gs * (vals[qp] * vals[qp]);

        for ( unsigned i = 0; i < dofs_per_cell; i++ )
        {
          cell_rhs(i) += JxW * (grads[qp] * fe_values.shape_grad(i, qp) + Q1 * vals[qp] * fe_values.shape_value(i, qp));
          for ( unsigned j = 0; j < dofs_per_cell; j++ )
            cell_matrix(i, j) += JxW * (fe_values.shape_value(i, qp) * fe_values.shape_value(j, qp));
        }
      }
      for ( unsigned i = 0; i < dofs_per_cell; i++ )
      {
        m_system_rhs(local_dof_indices[i]) += cell_rhs(i);
        for ( unsigned j = 0; j < dofs_per_cell; j++ )
          m_system_matrix.add (local_dof_indices[i], local_dof_indices[j], cell_matrix(i, j));
      }    
    }
    solve();
    VectorTools::integrate_difference ( dof_handler, m_newton_update, ZeroFunction<dim>(), m_error_per_cell, QGauss<1>(fe.degree + 2), VectorTools::L2_norm);
    err = m_error_per_cell.l2_norm();
  }  

   /** @brief computes multiple solution
   */
  template <int dim>
  void MySolver<dim>::run ()
  {
    string path;
    char shellcmd[255];
    double T, N, W;
    int status;

    make_grid();
    setup_system();

    CEigenfunctions<dim> Ef1( m_QN1, m_omega );
    VectorTools::interpolate (dof_handler, Ef1, m_Psi_0 );

    m_Psi_0 *= 1.0 / sqrt(MyRealTools::Particle_Number<dim>(dof_handler,fe,m_Psi_0));
    m_Psi_1 = 0;

    compute_E_lin( m_Psi_0, T, N, W );
    double m_mu_0 = T / N;
    m_mu = ceil(10.0 * m_mu_0) / 10.0 + m_gs / fabs(m_gs) * m_dmu;

    cout << "T = " << T << endl;
    cout << "N = " << N << endl;
    cout << "W = " << W << endl;
    cout << "m_mu = " << m_mu << endl;
    cout << "m_gs = " << m_gs << endl;

    output_guess();
    m_results.clear();
    for ( unsigned i = 0; i < m_Ndmu; i++ )
    {
      sprintf( shellcmd, "mkdir %.4d/", i );
      system(shellcmd);
      sprintf( shellcmd, "%.4d/", i );
      path = shellcmd;

      // Nehari Reset
      // sqrt((m_mu*N-T)/(m_gs*W)); if m_Psi_1 == 0
      // sqrt((m_mu*N-T)/(4.0*m_gs*W)); if m_Psi_1 == m_Psi_0
#ifdef NEHARI
      m_ti = sqrt((m_mu * N - T) / (m_gs * W));
#endif

      status = DoIter(path);
      estimate_error(m_final_error);

      columns &cols = m_results.new_line();
      m_results.insert( cols, MyTable::MU, m_mu );
      m_results.insert( cols, MyTable::GS, m_gs );
      m_results.insert( cols, MyTable::PARTICLE_NUMBER, m_N );
      m_results.insert( cols, MyTable::COUNTER, double(m_counter) );
      m_results.insert( cols, MyTable::RES, m_final_error );
      m_results.insert( cols, MyTable::STATUS, double(status) );

      if ( status == Status::SUCCESS )
      {
        output_results(path, "final");
        save( path + "final.bin");
        save_one( path + "final_one.bin");
        
        vector<double> newgs = {m_gs*m_N}; 
        m_ph.Set_Physics( "gs_1", newgs );
        m_ph.SaveXMLFile( path + "params_one.xml" );
        newgs[0] = m_gs;
        m_ph.Set_Physics( "gs_1", newgs );
        m_ph.SaveXMLFile( path + "params.xml" );

        m_Psi_0 = m_Psi_ref;
        m_Psi_1 = 0;
      }
      else
      {
        break;
      }
      m_mu += m_gs / fabs(m_gs) * m_dmu;
      compute_E_lin( m_Psi_ref, T, N, W );
    }
    m_results.dump_2_file( "results.csv" );
  }

  template <int dim>
  void MySolver<dim>::output_guess ()
  {
    string filename = "guess.gnuplot";

    CPotential<dim> Pot( m_omega );
    VectorTools::interpolate (dof_handler, Pot, m_workspace );

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
    string filename = path + prefix + "-" + Utilities::int_to_string (m_counter, 5) + ".gnuplot";

    DataOut<dim> data_out;
    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (m_Psi_ref, "Psi_ref");
    data_out.add_data_vector (m_Psi_0, "Psi_1");
    data_out.add_data_vector (m_Psi_1, "Psi_2");
    data_out.add_data_vector (m_newton_update, "m_newton_update");
    data_out.add_data_vector (m_error_per_cell, "error_per_cell");
    data_out.build_patches ();

    ofstream output (filename.c_str());
    data_out.write_gnuplot (output);
  }

  template <int dim>
  void MySolver<dim>::save( string filename )
  {
    ofstream out(filename);
    m_Psi_ref.block_write(out);
  }

  template <int dim>
  void MySolver<dim>::save_one( string filename )
  {
    ofstream out(filename);
  
    m_workspace = m_Psi_ref;
    m_workspace *= 1.0/sqrt(MyRealTools::Particle_Number<dim>(dof_handler,fe,m_Psi_ref));
    m_workspace.block_write(out);
  }
} // end of namespace

int main ( int argc, char *argv[] )
{
  deallog.depth_console (0);

  BreedSolver_1::MySolver<1> solver("params.xml");
  solver.run();
return EXIT_SUCCESS;
}
