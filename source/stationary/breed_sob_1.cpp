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

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/vector.h>

#include <deal.II/base/function.h>
#include <deal.II/base/function_parser.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>

#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <cmath>
#include <limits>

#include "MyComplexTools.h"
#include "MyParameterHandler.h"
#include "MyRealTools.h"
#include "functions.h"
#include "global.h"
#include "my_table.h"

#define STR1(x) #x
#define STR2(x) STR1(x)

namespace BreedSolver_1
{
  using namespace std;
  using namespace dealii;

  enum Status { SUCCESS, FAILED, ZERO_SOL, SLOW_CONV };

  template <int dim> 
  class MySolver
  {
  public:
    explicit MySolver(const std::string&);
    virtual ~MySolver();

    void run();

  protected:
    double m_xmin=-5.0, m_xmax=5.0;
    double m_ymin=-5.0, m_ymax=5.0;
    double m_res, m_res_old, m_resp;
    double m_final_error = 0.0;
    double m_N = 0.0;

    double m_mu=0.0;
    double m_gs=1.0;
    vector<double> m_omega;
    vector<double> m_epsilon;

    unsigned m_counter=0;
    unsigned m_global_refinement;
    unsigned m_NA;

    int DoIter(string = "");

    void make_grid();
    void setup_system();
    void assemble_rhs();
    void assemble_system();
    void compute_Psi_sob();
    void compute_mu();
    void save(string);
    void Project_gradient();

    void solve();
    void compute_E_lin(Vector<double> &, double &, double &, double &);
    void estimate_error(double &);

    void output_results(string);

    MyParameterHandler m_ph;
    Triangulation<dim> triangulation;
    FE_Q<dim> fe;
    DoFHandler<dim> dof_handler;

    SparsityPattern m_sparsity_pattern;
    SparseMatrix<double> m_system_matrix;
    Vector<double> m_system_rhs;
    Vector<double> m_Psi;
    Vector<double> m_Psi_C_ghosted;
    Vector<double> m_sob_grad;
    Vector<double> m_Psi_sob;
    Vector<double> m_solution;
    Vector<double> m_error_per_cell;

    string m_guess_str;
    MyTable m_table;
    MyTable m_results;
  };

  /**
   * Constructor
   */
  template <int dim>
  MySolver<dim>::MySolver(const std::string& xmlfilename) : 
    m_ph(xmlfilename),
    triangulation(),
    fe(gl_degree_fe), 
    dof_handler(triangulation)
  {
    try
    {
      m_omega = m_ph.Get_Physics("omega");
      m_gs = m_ph.Get_Physics("gs_1", 0);

      m_xmin = m_ph.Get_Mesh("xrange", 0);
      m_xmax = m_ph.Get_Mesh("xrange", 1);
      m_ymin = m_ph.Get_Mesh("yrange", 0);
      m_ymax = m_ph.Get_Mesh("yrange", 1);
      m_global_refinement = unsigned(m_ph.Get_Mesh("global_refinements", 0));

      m_NA = int(m_ph.Get_Algorithm("NA", 0));
      m_epsilon = m_ph.Get_Algorithm("epsilon");
      m_guess_str = m_ph.Get_Parameter("guess_fct");
    }
    catch (const std::string& info)
    {
      std::cerr << info << endl;
    }

    m_counter = 0;
    m_final_error = 0;
  }

  template <int dim>
  MySolver<dim>::~MySolver()
  {
    dof_handler.clear();
  }

  template <int dim>
  void MySolver<dim>::compute_E_lin(Vector<double> &vec, double &T, double &N, double &W)
  {
    CPotential<dim> Potential(m_omega);
    const QGauss<dim> quadrature_formula(fe.degree + 1);
    FEValues<dim> fe_values(fe, quadrature_formula, update_gradients | update_values | update_JxW_values | update_quadrature_points);

    const unsigned n_q_points = quadrature_formula.size();
    vector<double> vals(n_q_points);
    vector<Tensor<1, dim>> grad(n_q_points);

    T = 0;
    N = 0;
    W = 0;

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell != endc; ++cell)
    {
      fe_values.reinit(cell);
      fe_values.get_function_values(vec, vals);
      fe_values.get_function_gradients(vec, grad);
      for (unsigned qp = 0; qp < n_q_points; ++qp)
      {
        const double vec_val_q = vals[qp] * vals[qp];
        const double JxW = fe_values.JxW(qp);
        T += JxW * (grad[qp] * grad[qp] + Potential.value(fe_values.quadrature_point(qp)) * vec_val_q);
        N += JxW * vec_val_q;
        W += JxW * vec_val_q * vec_val_q;
      }
    }
  }

  template <int dim> 
  void MySolver<dim>::Project_gradient()
  {
    const QGauss<dim> quadrature_formula(fe.degree + 1);
    FEValues<dim> fe_values(fe, quadrature_formula, update_values | update_JxW_values | update_quadrature_points);

    const unsigned n_q_points = quadrature_formula.size();
    vector<double> vals_1(n_q_points);
    vector<double> vals_2(n_q_points);
    vector<double> vals_3(n_q_points);

    double sum[] = {0, 0};

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell != endc; ++cell)
    {
      fe_values.reinit(cell);
      fe_values.get_function_values(m_Psi_sob, vals_1);
      fe_values.get_function_values(m_sob_grad, vals_2);
      fe_values.get_function_values(m_Psi, vals_3);

      for (unsigned qp = 0; qp < n_q_points; ++qp)
      {
        const double JxW = fe_values.JxW(qp);
        sum[0] += JxW * (vals_3[qp] * vals_2[qp]);
        sum[1] += JxW * (vals_3[qp] * vals_1[qp]);
      }
    }
    m_sob_grad.add(-sum[0] / sum[1], m_Psi_sob);
  }

  template <int dim> 
  void MySolver<dim>::estimate_error(double &err)
  {
    compute_mu();

    CPotential<dim> Potential( m_omega );
    const QGauss<dim> quadrature_formula(fe.degree+1);

    m_system_rhs=0;
    m_system_matrix=0;

    FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_gradients|update_JxW_values|update_quadrature_points);

    const unsigned dofs_per_cell = fe.dofs_per_cell;
    const unsigned n_q_points = quadrature_formula.size();

    Vector<double> cell_rhs (dofs_per_cell);
    FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);

    vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    vector<double> vals(n_q_points);
    vector<Tensor<1,dim>> grads(n_q_points);

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for ( ; cell!=endc; ++cell )
    {
      cell_rhs=0;
      cell_matrix=0;

      fe_values.reinit (cell);
      fe_values.get_function_values(m_Psi, vals);
      fe_values.get_function_gradients(m_Psi, grads);

      for ( unsigned qp=0; qp<n_q_points; ++qp )
      {
          const double JxW = fe_values.JxW(qp);
          const double Q1 = Potential.value(fe_values.quadrature_point(qp)) - m_mu + m_gs*(vals[qp]*vals[qp]);

          for ( unsigned i=0; i<dofs_per_cell; ++i )
          {
            cell_rhs(i) += JxW*(grads[qp]*fe_values.shape_grad(i,qp) + Q1*vals[qp]*fe_values.shape_value(i,qp));
            for ( unsigned j=0; j<dofs_per_cell; ++j )
                cell_matrix(i,j) += JxW*(fe_values.shape_value(i,qp)*fe_values.shape_value(j,qp));
          }
        }
        cell->get_dof_indices (local_dof_indices);

        for (unsigned i = 0; i < dofs_per_cell; ++i)
        {
          m_system_rhs(local_dof_indices[i]) += cell_rhs(i);
          for (unsigned j = 0; j < dofs_per_cell; ++j)
            m_system_matrix.add(local_dof_indices[i], local_dof_indices[j], cell_matrix(i, j));
        }      
      }

      solve();

      VectorTools::integrate_difference ( dof_handler, m_solution, ZeroFunction<dim>(), m_error_per_cell, QGauss<dim>(fe.degree+2), VectorTools::L2_norm);
      err = m_error_per_cell.l2_norm();
  }

  template <int dim> 
  void MySolver<dim>::assemble_rhs()
  {
    m_system_rhs = 0;

    CPotential<dim> Potential(m_omega);
    const QGauss<dim> quadrature_formula(fe.degree + 1);

    FEValues<dim> fe_values(fe, quadrature_formula, update_values | update_gradients | update_JxW_values | update_quadrature_points);

    const unsigned dofs_per_cell = fe.dofs_per_cell;
    const unsigned n_q_points = quadrature_formula.size();

    Vector<double> cell_rhs(dofs_per_cell);
    vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
    vector<double> vals(n_q_points);
    vector<Tensor<1, dim>> grads(n_q_points);

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                   endc = dof_handler.end();
    for (; cell != endc; ++cell)
    {
      cell_rhs = 0;
      fe_values.reinit(cell);
      fe_values.get_function_values(m_Psi, vals);
      fe_values.get_function_gradients(m_Psi, grads);

      for (unsigned qp = 0; qp < n_q_points; ++qp)
      {
        const double JxW = fe_values.JxW(qp);
        const double Q1 = Potential.value(fe_values.quadrature_point(qp)) + m_gs * (vals[qp] * vals[qp]);

        for (unsigned i = 0; i < dofs_per_cell; ++i)
          cell_rhs(i) +=  JxW * (grads[qp] * fe_values.shape_grad(i, qp) + Q1 * vals[qp] * fe_values.shape_value(i, qp));
      }
      cell->get_dof_indices(local_dof_indices);
      for ( unsigned i = 0; i < dofs_per_cell; ++i )
        m_system_rhs(local_dof_indices[i]) += cell_rhs(i);      
    }
  }

  template <int dim> 
  void MySolver<dim>::assemble_system()
  {
    const QGauss<dim> quadrature_formula(fe.degree + 1);

    m_system_matrix = 0;

    FEValues<dim> fe_values(fe, quadrature_formula, update_values | update_gradients | update_JxW_values | update_quadrature_points);

    const unsigned dofs_per_cell = fe.dofs_per_cell;
    const unsigned n_q_points = quadrature_formula.size();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell != endc; ++cell)
    {
      cell_matrix = 0;
      fe_values.reinit(cell);

      for (unsigned qp = 0; qp < n_q_points; ++qp)
      {
        const double JxW = fe_values.JxW(qp);

        for (unsigned i = 0; i < dofs_per_cell; ++i)
          for (unsigned j = 0; j < dofs_per_cell; ++j)
            cell_matrix(i, j) += JxW * (fe_values.shape_grad(i, qp) * fe_values.shape_grad(j, qp) + fe_values.shape_value(i, qp) * fe_values.shape_value(j, qp));
      }
      cell->get_dof_indices (local_dof_indices);
      for ( unsigned i = 0; i < dofs_per_cell; ++i )
        for ( unsigned j = 0; j < dofs_per_cell; ++j )
          m_system_matrix.add (local_dof_indices[i], local_dof_indices[j], cell_matrix(i, j));      
    }
  }

  template <int dim> 
  void MySolver<dim>::compute_Psi_sob()
  {
    m_system_rhs = 0;

    CPotential<dim> Potential(m_omega);
    const QGauss<dim> quadrature_formula(fe.degree + 1);

    FEValues<dim> fe_values(fe, quadrature_formula, update_values | update_gradients | update_JxW_values );

    const unsigned dofs_per_cell = fe.dofs_per_cell;
    const unsigned n_q_points = quadrature_formula.size();

    vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
    vector<double> vals(n_q_points);
    Vector<double> cell_rhs(dofs_per_cell);
    FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell != endc; ++cell)
    {
      cell_rhs = 0;
      cell_matrix = 0;
      fe_values.reinit(cell);
      fe_values.get_function_values(m_Psi, vals);
      cell->get_dof_indices(local_dof_indices);

      for (unsigned qp = 0; qp < n_q_points; ++qp)
      {
        const double JxW = fe_values.JxW(qp);

        for (unsigned i = 0; i < dofs_per_cell; ++i)
        {
          cell_rhs(i) += JxW * vals[qp] * fe_values.shape_value(i, qp);
          for (unsigned j = 0; j < dofs_per_cell; ++j)
            cell_matrix(i, j) += JxW * (fe_values.shape_grad(i, qp) * fe_values.shape_grad(j, qp) + fe_values.shape_value(i, qp) * fe_values.shape_value(j, qp));
        }
      }
      cell->get_dof_indices (local_dof_indices);
      for ( unsigned i = 0; i < dofs_per_cell; ++i )
      {
        m_system_rhs(local_dof_indices[i]) += cell_rhs(i);
        for ( unsigned j = 0; j < dofs_per_cell; ++j )
          m_system_matrix.add (local_dof_indices[i], local_dof_indices[j], cell_matrix(i, j));    
      } 
    }
    solve();
    m_Psi_sob = m_solution;
  }

  template <int dim> 
  void MySolver<dim>::compute_mu()
  {
    CPotential<dim> Potential(m_omega);
    const QGauss<dim> quadrature_formula(fe.degree + 1);
    FEValues<dim> fe_values(fe, quadrature_formula, update_gradients | update_values | update_JxW_values | update_quadrature_points);

    const unsigned n_q_points = quadrature_formula.size();
    vector<double> vals(n_q_points);
    vector<Tensor<1, dim>> grads(n_q_points);

    double psum = 0;
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell != endc; ++cell)
    {
      fe_values.reinit(cell);
      fe_values.get_function_values(m_Psi, vals);
      fe_values.get_function_gradients(m_Psi, grads);

      for (unsigned qp = 0; qp < n_q_points; ++qp)
      {
        const double uq = vals[qp] * vals[qp];
        psum += fe_values.JxW(qp) * (grads[qp] * grads[qp] + (Potential.value(fe_values.quadrature_point(qp)) + m_gs * uq) * uq);
      }
    }
  }

  template <int dim> 
  void MySolver<dim>::solve()
  {
    map<types::global_dof_index, double> boundary_values;
    VectorTools::interpolate_boundary_values (dof_handler, 0, ZeroFunction<dim>(), boundary_values);
    MatrixTools::apply_boundary_values (boundary_values, m_system_matrix, m_solution, m_system_rhs);    

    SparseDirectUMFPACK A_direct;
    A_direct.initialize(m_system_matrix);
    A_direct.vmult(m_solution, m_system_rhs);
  }

  template <int dim> 
  void MySolver<dim>::make_grid()
  {
    Point<dim, double> pt1, pt2;

    double min[] = {m_xmin, m_ymin};
    double max[] = {m_xmax, m_ymax};

    for ( int i = 0; i < dim; ++i )
    {
      pt1(i) = min[i];
      pt2(i) = max[i];
    }

    GridGenerator::hyper_rectangle(triangulation, pt2, pt1);
    triangulation.refine_global(m_global_refinement);
  }

  template <int dim> 
  void MySolver<dim>::setup_system()
  {
    dof_handler.distribute_dofs(fe);

    m_Psi.reinit(dof_handler.n_dofs());
    m_Psi_sob.reinit(dof_handler.n_dofs());
    m_system_rhs.reinit(dof_handler.n_dofs());
    m_solution.reinit(dof_handler.n_dofs());
    m_sob_grad.reinit(dof_handler.n_dofs());
    m_error_per_cell.reinit(triangulation.n_active_cells());

    cout << "no of dofs = " << dof_handler.n_dofs() << endl;

    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, dsp);
    m_sparsity_pattern.copy_from(dsp);
    m_system_matrix.reinit(m_sparsity_pattern);
  }

  template <int dim>
  void MySolver<dim>::output_results( string filename ) 
  {
    DataOut<dim> data_out;
    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (m_Psi, "Psi");
    data_out.add_data_vector (m_error_per_cell, "error_per_cell");
    data_out.build_patches ();
    
    ofstream output (filename.c_str());
    data_out.write_gnuplot (output);
  }

  template <int dim> 
  int MySolver<dim>::DoIter(string path)
  {
    int retval = Status::SUCCESS;

    m_table.clear();

    assemble_rhs();

    m_res = 0;
    m_counter = 0;
    do
    {
      cout << "--------------------------------------------------------------------------------" << endl;
      cout << "- " << m_counter << endl;

      assemble_system();
      solve();
      m_sob_grad = m_solution;
      
      compute_Psi_sob();
      Project_gradient();
      m_res = m_sob_grad.l2_norm();
      
      m_Psi.add( -0.1, m_sob_grad);
      
      // compute_mu(m_mu);
      m_N = MyRealTools::Particle_Number(dof_handler, fe, m_Psi);
      
      if (fabs(m_N - 1) > 1e-5)
        m_Psi *= 1 / sqrt(m_N);

      assemble_rhs();
        
      m_resp = m_res_old - m_res;
      m_res_old = m_res;

      columns &cols = m_table.new_line();
      m_table.insert(cols, MyTable::COUNTER, double(m_counter));
      m_table.insert(cols, MyTable::RES, m_res);
      m_table.insert(cols, MyTable::RESP, m_resp);
      m_table.insert(cols, MyTable::MU, m_mu);
      m_table.insert(cols, MyTable::GS, m_gs);
      m_table.insert(cols, MyTable::PARTICLE_NUMBER, m_N);

      cout << m_table;
      if (m_res < m_epsilon[1])
      {
        retval = Status::SUCCESS;
        break;
      }

      m_counter++;
    }
    while (true);

    m_N = MyRealTools::Particle_Number(dof_handler, fe, m_Psi);

    if (m_N < 1e-5)
      retval = Status::ZERO_SOL;

    string filename = path + "log.csv";
    m_table.dump_2_file(filename);

    return retval;
  }

  template <int dim> 
  void MySolver<dim>::run()
  {
    int status;
    double T, N, W;

    make_grid();
    setup_system();

    unsigned QN[] = {0, 0, 0};
    CEigenfunctions<dim> Ef1(QN, m_omega);
    VectorTools::interpolate(dof_handler, Ef1, m_Psi);

    compute_E_lin(m_Psi, T, N, W);
    m_Psi *= sqrt(1 / N);

    cout << setprecision(9);
    cout << "T = " << T << endl;
    cout << "N = " << N << endl;
    cout << "W = " << W << endl;

    status = DoIter("");

    if (status == Status::SUCCESS)
    {
      output_results("final_gs_one.gnuplot");
      save("final_gs_one.bin");
    }

    ofstream ofs("log.txt");
    ofs << m_table;
  }

  template <int dim> 
  void MySolver<dim>::save(string filename) 
  {
    ofstream ofs(filename);
    m_Psi.block_write(ofs);
  }
} // end of namespace

int main(int argc, char *argv[])
{
  using namespace dealii;
  deallog.depth_console(0);

  BreedSolver_1::MySolver<1> solver("params_one.xml");
  solver.run();
  return EXIT_SUCCESS;
}