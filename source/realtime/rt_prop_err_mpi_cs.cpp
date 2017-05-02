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
 * Purpose: Real time propagation for the Gross-Pitaevskii equation (cylinder symmetry)
 * Method: fully implicit Crank-Nicolson 
 */

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/parameter_handler.h>

#include <deal.II/lac/generic_linear_algebra.h>

namespace LA
{
  using namespace dealii::LinearAlgebraPETSc;
}

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>

#include <deal.II/lac/petsc_parallel_sparse_matrix.h>
#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_precondition.h>

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
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/solution_transfer.h>

#include <fstream>
#include <iostream>
#include <regex>
#include <stdlib.h>

#include "global.h"
#include "mpi.h"
#include "my_table.h"
#include "functions_cs.h"
#include "tinyxml2.h"
#include "MyParameterHandler.h"

#include "sys/types.h"
#include "sys/stat.h"
#include "dirent.h"

void Find_last_file( string path, double& t, string& filename )
{
  DIR *pDir;
  struct dirent *pDirEnt;
  struct stat fileinfo;
  string info(".info");
  string prefix("step-");
  string tmpstr;
  string filename2;
  
  pDir = opendir( path.c_str() );
  pDirEnt = readdir( pDir );
  while( pDirEnt != NULL ) 
  {
    if( stat( pDirEnt->d_name, &fileinfo ) == -1 ) { pDirEnt = readdir( pDir ); continue; }
    filename2 = pDirEnt->d_name;
    size_t pos = filename2.find("info");
    if( pos!=string::npos ) { pDirEnt = readdir( pDir ); continue; }
    pos = filename2.find(".bin");
    if( pos==string::npos ) { pDirEnt = readdir( pDir ); continue; }

    
    if( filename2.compare(0, prefix.length(), prefix) == 0 )
    {
      tmpstr=filename2;
      tmpstr.erase(0,prefix.length());
      tmpstr.erase(tmpstr.length()-4,4);
      double tmp = stod(tmpstr);
      if( tmp > t )
      {
        t=tmp;
        filename=filename2;
      }
    }
    pDirEnt = readdir( pDir );
  }
  closedir( pDir );
}

namespace realtime_propagation
{
  enum Status { SUCCESS, FAILED };

  using namespace std;
  using namespace dealii;
  
  template <int dim>
  class MySolver
  {
  public:
    MySolver( const std::string );
    ~MySolver();

    void run ();
    double Particle_Number( LA::MPI::Vector& );
    void Expectation_value_momentum( LA::MPI::Vector&, double* );
    void Expectation_value_position( LA::MPI::Vector&, double* );

  protected:
    void make_grid();
    void setup_system( const bool );
    void setup_boundary_ids();
    void assemble_system( LA::MPI::Vector& );
    void assemble_rhs( LA::MPI::Vector& );
    double compute_correlation();
    void orthogonalize( LA::MPI::Vector&, LA::MPI::Vector& );
    
    void compute_eta();
    
    void DoIter( LA::MPI::Vector& );
    void solve();
    void output_results ( string );
    
    void load( string );
    void load_2( string );
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
    LA::MPI::Vector newton_update;
    LA::MPI::Vector m_Psi; // Psi
    LA::MPI::Vector m_Psi_t; // Psi trial
    LA::MPI::Vector m_Psi_0;
    LA::MPI::Vector eta;
    LA::MPI::Vector m_workspace_1;
    LA::MPI::Vector m_workspace_2;
    Vector<double> m_error_per_cell;

    ConditionalOStream pcout;
    ofstream m_computing_timer_log;
    TimerOutput m_computing_timer;    

    bool m_root;

    double m_mu;
    double m_gs;
    double m_t;
    double m_dt;
    vector<double> m_omega;
    double m_res;

    double m_xmin;
    double m_xmax;
    double m_ymin;
    double m_ymax;
    double m_zmin;
    double m_zmax;

    int m_rank;
    unsigned int m_NA;
    unsigned int m_NK;
    unsigned int m_global_refinement;
    unsigned int m_QN1[3];

    MyTable m_table;    
  };

  template <int dim>
  MySolver<dim>::MySolver ( const std::string xmlfilename ) 
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
      m_QN1[0] = int(m_ph.Get_Physics("QN1",0));
      m_QN1[1] = int(m_ph.Get_Physics("QN1",1));
      m_QN1[2] = int(m_ph.Get_Physics("QN1",2));      

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
    m_t = 0;

    m_root = (Utilities::MPI::this_mpi_process(mpi_communicator) == 0);
    MPI_Comm_rank(mpi_communicator, &m_rank);
  }

  template <int dim>
  MySolver<dim>::~MySolver ()
  {
    dof_handler.clear ();
  }

  #include "shared_rt_prop_cs.h"

  template <int dim>
  void MySolver<dim>::compute_eta ()
  {
    m_computing_timer.enter_section(__func__);
    const QGauss<dim> quadrature_formula(fe.degree+1);

    CPotential Potential ( m_omega, m_QN1[2] );

    const FEValuesExtractors::Scalar rt (0);
    const FEValuesExtractors::Scalar it (1);
        
    system_matrix = 0;
    system_rhs = 0;

    FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_gradients|update_quadrature_points|update_JxW_values);

    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    const unsigned int n_q_points    = quadrature_formula.size();

    Vector<double> cell_rhs (dofs_per_cell);
    FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);

    vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    vector<Vector<double>> Psi0(n_q_points,Vector<double>(2));
    vector<vector<Tensor<1,dim>>> Psi0_grad(n_q_points, vector<Tensor<1,dim>>(2));
 
    double JxW, pot=0, tmp1a, tmp1b, tmp2, sum_re, sum_req, sum_im, sum_imq;

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
      if( cell->is_locally_owned() )
      {
        cell_matrix=0;
	      cell_rhs=0;

        fe_values.reinit (cell);
        fe_values.get_function_values(m_Psi_0, Psi0);
        fe_values.get_function_gradients(m_Psi_0, Psi0_grad);

        for( unsigned int qp=0; qp<n_q_points; qp++ )
        {
          JxW = fe_values.JxW(qp)*fabs(fe_values.quadrature_point(qp)[1]);
          pot = Potential.value(fe_values.quadrature_point(qp)) - m_mu + m_gs*Psi0[qp][0]*Psi0[qp][0];        

          for (unsigned int i=0; i<dofs_per_cell; i++ )
          {
            for (unsigned int j=0; j<dofs_per_cell; j++ )
            {
              cell_matrix(i,j) += JxW*( fe_values[rt].value(i,qp)*fe_values[rt].value(j,qp) +
                                        fe_values[rt].value(i,qp)*fe_values[it].value(j,qp) +
					fe_values[it].value(i,qp)*fe_values[rt].value(j,qp) + 
					fe_values[it].value(i,qp)*fe_values[it].value(j,qp) );
            }
            cell_rhs(i) += JxW*(Psi0_grad[qp][0]*fe_values[rt].gradient(i,qp) + pot*Psi0[qp][0]*fe_values[rt].value(i,qp) + Psi0_grad[qp][1]*fe_values[it].gradient(i,qp) + pot*Psi0[qp][1]*fe_values[it].value(i,qp) );            
	  }
        }
        cell->get_dof_indices (local_dof_indices);
        constraints.distribute_local_to_global( cell_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs );
      }
    }
    system_matrix.compress(VectorOperation::add);
    system_rhs.compress(VectorOperation::add);
    
    pcout << "HERE 1\n" << endl;
    SolverControl solver_control;
    
    PETScWrappers::SparseDirectMUMPS solver(solver_control, mpi_communicator);
    solver.set_symmetric_mode(false);
    solver.solve(system_matrix, newton_update, system_rhs);
    constraints.distribute(newton_update);
    pcout << "HERE 2\n" << endl;
    
    eta=newton_update;
    
    m_workspace_1 = system_rhs;
    constraints.distribute(m_workspace_1);
    
    VectorTools::integrate_difference ( dof_handler,  m_workspace_1, ZeroFunction<dim>(), m_error_per_cell,  QGauss<dim>(fe.degree+2), VectorTools::L2_norm);
    double tmp1 = m_error_per_cell.l2_norm();
    tmp1 = tmp1*tmp1;
    MPI_Allreduce( &tmp1, &m_res, 1, MPI_DOUBLE, MPI_SUM, mpi_communicator);
    printf( "err = %g\n", sqrt(m_res) );    
    
    m_computing_timer.exit_section();    
  }
  
  template <int dim>
  void MySolver<dim>::orthogonalize( LA::MPI::Vector& A, LA::MPI::Vector& B )
  {
    m_computing_timer.enter_section(__func__);
    double tmp[] = {0,0,0,0}, s[] = {0,0,0,0}, JxW;

    constraints.distribute(A);
    m_workspace_1 = A;
    constraints.distribute(B);
    m_workspace_2 = B;

    const QGauss<dim>  quadrature_formula(fe.degree+1);
    FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_quadrature_points|update_JxW_values);

    const unsigned int n_q_points = quadrature_formula.size();
    vector<Vector<double>> vec_vals_A(n_q_points,Vector<double>(2));
    vector<Vector<double>> vec_vals_B(n_q_points,Vector<double>(2));

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
      if( cell->is_locally_owned() )
      {
        fe_values.reinit (cell);
        fe_values.get_function_values( m_workspace_1, vec_vals_A );      
        fe_values.get_function_values( m_workspace_2, vec_vals_B );      

        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
        {
          JxW = fe_values.JxW(q_point)*fabs(fe_values.quadrature_point(q_point)[1]);
          tmp[0] += JxW*vec_vals_A[q_point][0]*vec_vals_B[q_point][0];
          tmp[1] += JxW*vec_vals_A[q_point][0]*vec_vals_A[q_point][0];
          tmp[2] += JxW*vec_vals_A[q_point][1]*vec_vals_B[q_point][1];
          tmp[3] += JxW*vec_vals_A[q_point][1]*vec_vals_A[q_point][1];
      	}
      }
    }
    MPI_Allreduce( tmp, s, 4, MPI_DOUBLE, MPI_SUM, mpi_communicator);

    double f1 = -s[0]/s[1];
    double f2 = -s[2]/s[3];
    
    if( std::isnan(f1) ) f1 = 0;
    if( std::isnan(f2) ) f2 = 0;

    vector<bool> selected_dofs (dof_handler.n_locally_owned_dofs());
    vector<bool> mask (dof_handler.get_fe().n_components(), false );
    mask[0] = true;

    DoFTools::extract_dofs (dof_handler, ComponentMask(mask), selected_dofs);

    vector<types::global_dof_index> indices(dof_handler.n_locally_owned_dofs());
    locally_owned_dofs.fill_index_vector ( indices );

    for( unsigned int i=0; i<indices.size(); i++ )
    {
      if( selected_dofs[i] == true )
        system_rhs[indices[i]] = B[indices[i]] + f1*A[indices[i]];
      else
        system_rhs[indices[i]] = B[indices[i]] + f2*A[indices[i]]; 
    }
    system_rhs.compress(VectorOperation::insert);
    constraints.distribute(system_rhs); 
    B = system_rhs;
    m_computing_timer.exit_section();
  }

  template <int dim>
  void MySolver<dim>::assemble_system ( LA::MPI::Vector& vec )
  {
    m_computing_timer.enter_section(__func__);
    const QGauss<dim> quadrature_formula(fe.degree+1);

    CPotential Potential ( m_omega, m_QN1[2] );

    const FEValuesExtractors::Scalar rt (0);
    const FEValuesExtractors::Scalar it (1);
        
    system_matrix = 0;

    constraints.distribute(vec);
    m_workspace_1=vec;
    constraints.distribute(m_Psi_t);
    m_workspace_2=m_Psi_t;    

    FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_gradients|update_quadrature_points|update_JxW_values);

    const unsigned dofs_per_cell = fe.dofs_per_cell;
    const unsigned n_q_points    = quadrature_formula.size();

    FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);

    vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    vector<Vector<double>> Psi(n_q_points,Vector<double>(2));
    vector<vector<Tensor<1,dim>>> Psi_grad(n_q_points, vector<Tensor<1,dim>>(2));
    vector<Vector<double>> Psi_t(n_q_points,Vector<double>(2));
    vector<vector<Tensor<1,dim>>> Psi_t_grad(n_q_points, vector<Tensor<1,dim>>(2));
    vector<Vector<double>> Psi0(n_q_points,Vector<double>(2));
 
    double JxW, pot=0, tmp1a, tmp1b, tmp2, sum_re, sum_req, sum_im, sum_imq;

    const double fak2 = 0.5*m_dt;
    const double fak4 = 0.25*m_gs*m_dt;
    const double fak8 = 0.125*m_gs*m_dt;

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
      if( cell->is_locally_owned() )
      {
        cell_matrix = 0;

        fe_values.reinit (cell);
        fe_values.get_function_values(m_workspace_1, Psi);
        fe_values.get_function_gradients(m_workspace_1, Psi_grad);
        fe_values.get_function_values(m_workspace_2, Psi_t);
        fe_values.get_function_gradients(m_workspace_2, Psi_t_grad);
        fe_values.get_function_values(m_Psi_0, Psi0);

        for( unsigned qp=0; qp<n_q_points; qp++ )
        {
          JxW = fe_values.JxW(qp)*fabs(fe_values.quadrature_point(qp)[1]);
          pot = Potential.value(fe_values.quadrature_point(qp)) - m_mu + m_gs*Psi0[qp][0]*Psi0[qp][0];        

          sum_re = Psi[qp][0]+Psi_t[qp][0];
          sum_im = Psi[qp][1]+Psi_t[qp][1];
          sum_req = sum_re*sum_re;
          sum_imq = sum_im*sum_im;
          tmp1a = fak8*(sum_req + 3*sum_imq);
          tmp1b = fak8*(sum_imq + 3*sum_req);
      	  tmp2 = fak4*sum_re*sum_im;
	  
          for (unsigned i=0; i<dofs_per_cell; i++ )
          {
            for (unsigned j=0; j<dofs_per_cell; j++ )
            {
	            // J00
              cell_matrix(i,j) += JxW*(1.0-tmp2+fak2*m_gs*Psi0[qp][0]*sum_im)*fe_values[rt].value(i,qp)*fe_values[rt].value(j,qp);
              // J01
              cell_matrix(i,j) -= JxW*tmp1a*fe_values[rt].value(i,qp)*fe_values[it].value(j,qp);
              cell_matrix(i,j) -= JxW*fak2*(fe_values[rt].gradient(i,qp)*fe_values[it].gradient(j,qp) + (pot-m_gs*Psi0[qp][0]*sum_re+m_gs*Psi0[qp][0]*Psi0[qp][0])*fe_values[rt].value(i,qp)*fe_values[it].value(j,qp));
              // J10
              cell_matrix(i,j) += JxW*tmp1b*fe_values[it].value(i,qp)*fe_values[rt].value(j,qp);
              cell_matrix(i,j) += JxW*fak2*(fe_values[it].gradient(i,qp)*fe_values[rt].gradient(j,qp) + (pot-3*m_gs*Psi0[qp][0]*sum_re-m_gs*Psi0[qp][0]*Psi0[qp][0])*fe_values[it].value(i,qp)*fe_values[rt].value(j,qp));
               // J11
	            cell_matrix(i,j) += JxW*(1.0+tmp2-fak2*m_gs*Psi0[qp][0]*sum_im)*fe_values[it].value(i,qp)*fe_values[it].value(j,qp);
           }
	        }
        }
        cell->get_dof_indices (local_dof_indices);
        constraints.distribute_local_to_global( cell_matrix, local_dof_indices, system_matrix );
      }
    }
    system_matrix.compress(VectorOperation::add);
    m_computing_timer.exit_section();
  }
  
  template <int dim>
  void MySolver<dim>::assemble_rhs ( LA::MPI::Vector& vec )
  {
    m_computing_timer.enter_section(__func__);
    const QGauss<dim> quadrature_formula(fe.degree+1);

    CPotential Potential ( m_omega, m_QN1[2] );

    const FEValuesExtractors::Scalar rt (0);
    const FEValuesExtractors::Scalar it (1);
        
    system_rhs = 0;
    constraints.distribute(vec);
    m_workspace_1=vec;
    constraints.distribute(m_Psi_t);
    m_workspace_2=m_Psi_t;

    FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_gradients|update_quadrature_points|update_JxW_values);

    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    const unsigned int n_q_points    = quadrature_formula.size();

    Vector<double> cell_rhs (dofs_per_cell);

    vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    vector<Vector<double>> Psi(n_q_points,Vector<double>(2));
    vector<vector<Tensor<1,dim>>> Psi_grad(n_q_points, vector<Tensor<1,dim>>(2));
    vector<Vector<double>> Psi_t(n_q_points,Vector<double>(2));
    vector<vector<Tensor<1,dim>>> Psi_t_grad(n_q_points, vector<Tensor<1,dim>>(2));
    vector<Vector<double>> Psi0(n_q_points,Vector<double>(2));
 
    double JxW, pot=0, tmp1, sum_re, sum_im;

    const double fak2 = 0.5*m_dt;
    const double fak8 = 0.125*m_gs*m_dt;

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
      if( cell->is_locally_owned() )
      {
        cell_rhs = 0;

        fe_values.reinit (cell);
        fe_values.get_function_values(m_workspace_1, Psi);
        fe_values.get_function_gradients(m_workspace_1, Psi_grad);
        fe_values.get_function_values(m_workspace_2, Psi_t);
        fe_values.get_function_gradients(m_workspace_2, Psi_t_grad);
        fe_values.get_function_values(m_Psi_0, Psi0);	

        for( unsigned int qp=0; qp<n_q_points; qp++ )
        {
          JxW = fe_values.JxW(qp)*fabs(fe_values.quadrature_point(qp)[1]);
          pot = Potential.value(fe_values.quadrature_point(qp)) - m_mu + m_gs*Psi0[qp][0]*Psi0[qp][0];        
	  
          sum_re = Psi[qp][0]+Psi_t[qp][0];
          sum_im = Psi[qp][1]+Psi_t[qp][1];	  
          tmp1 = fak8*(sum_re*sum_re + sum_im*sum_im);	  
	  
          for (unsigned int i=0; i<dofs_per_cell; i++ )
          {
            cell_rhs(i) -= JxW*fak2*((Psi_grad[qp][1]+Psi_t_grad[qp][1])*fe_values[rt].gradient(i,qp) + pot*sum_im*fe_values[rt].value(i,qp) - m_gs*Psi[qp][0]*sum_re*sum_im*fe_values[rt].value(i,qp) + m_gs*Psi[qp][0]*Psi[qp][0]*sum_im*fe_values[rt].value(i,qp) );
            cell_rhs(i) += JxW*fak2*((Psi_grad[qp][0]+Psi_t_grad[qp][0])*fe_values[it].gradient(i,qp) + pot*sum_re*fe_values[it].value(i,qp) - m_gs*Psi0[qp][0]*(1.5*sum_re*sum_re+0.5*sum_im*sum_im)*fe_values[it].value(i,qp) - m_gs*Psi0[qp][0]*Psi0[qp][0]*sum_re*fe_values[it].value(i,qp));
            cell_rhs(i) += JxW*((Psi_t[qp][0]-Psi[qp][0])*fe_values[rt].value(i,qp) - tmp1*sum_im*fe_values[rt].value(i,qp) + 
                                (Psi_t[qp][1]-Psi[qp][1])*fe_values[it].value(i,qp) + tmp1*sum_re*fe_values[it].value(i,qp));
          }
        }
        cell->get_dof_indices (local_dof_indices);
        constraints.distribute_local_to_global(cell_rhs, local_dof_indices, system_rhs);
      }
    }
    system_rhs.compress(VectorOperation::add);
    m_res = system_rhs.l2_norm();
    
    m_computing_timer.exit_section();
  }  
  

  template <int dim>
  void MySolver<dim>::solve ()
  {
    m_computing_timer.enter_section(__func__);
    SolverControl solver_control;
    
    PETScWrappers::SparseDirectMUMPS solver(solver_control, mpi_communicator);
    solver.set_symmetric_mode(false);
    solver.solve(system_matrix, newton_update, system_rhs);
    constraints.distribute(newton_update);
    
    m_computing_timer.exit_section();
  }
  

  /*  
  template <int dim>
  void MySolver<dim>::solve ()
  {
    m_computing_timer.enter_section(__func__);

    newton_update=0;
    
    PETScWrappers::PreconditionSOR preconditioner;
    PETScWrappers::PreconditionSOR::AdditionalData data;
    preconditioner.initialize(system_matrix, data);

    SolverControl solver_control (dof_handler.n_dofs(), 1e-15);
    PETScWrappers::SolverGMRES solver(solver_control, mpi_communicator);
    
    solver.solve(system_matrix, newton_update, system_rhs, preconditioner);
    constraints.distribute (newton_update);

    m_computing_timer.exit_section();
  }    
  */
  
  template<int dim> 
  void MySolver<dim>::DoIter( LA::MPI::Vector& vec )
  {
    m_Psi_t = m_Psi;
    m_res = 0;
    assemble_rhs(vec);
    pcout << "m_res = " << m_res << endl;
    do
    {
      assemble_system(vec);
      solve();

      m_Psi_t.add( -1, newton_update );
      constraints.distribute(m_Psi_t);

      assemble_rhs(vec);
      pcout << "m_res = " << m_res << endl;
    }
    while( m_res > 1e-6 );

    vec = m_Psi_t;
  }

  template <int dim>
  void MySolver<dim>::run()
  {
    double N;
    
    tinyxml2::XMLDocument doc;
    doc.LoadFile( "info.xml" );
    tinyxml2::XMLElement* xml_element_mu = doc.FirstChildElement()->FirstChildElement( "MU" );
    xml_element_mu->QueryDoubleText( &m_mu );
    pcout << "m_mu = " << m_mu << endl;

    load( "final.bin" );

    m_Psi_0 = m_Psi;

//     unsigned q[] = { 0, 0, m_QN1[2]};
//     double f[] = { 0.5, 0.5, 0.5};
//     CEigenfunctions<dim> Ef1( q, f );
//     
//     VectorTools::interpolate (dof_handler, Ef1, system_rhs );
//     constraints.distribute(system_rhs);
//     m_Psi = system_rhs;
//     m_Psi *= sqrt(1e-4)/sqrt(Particle_Number(m_Psi));

    // Define some constants that will be used by the function parser
    std::map<std::string,double> constants;
    constants["pi"] = numbers::PI;
    
    std::vector<std::string> expressions(2);
    expressions[0] = "exp(-0.1*(x^2+y^2))";
    expressions[1] = "exp(-0.1*(x^2+y^2))";
    FunctionParser<2> vector_function(2);
    vector_function.initialize( FunctionParser<2>::default_variable_names(), expressions, constants);    
    
    VectorTools::interpolate (dof_handler, vector_function, m_Psi );
    constraints.distribute(m_Psi);
    //m_Psi *= sqrt(1e-4)/sqrt(Particle_Number(m_Psi));

    //orthogonalize( m_Psi_0, m_Psi );
    
    //compute_eta();
    
    N = Particle_Number(m_Psi);
    
    output_results("");

    pcout << "t == " << m_t << endl;
    pcout << "N == " << N << endl;

    MyTable table;
    columns& cols = m_table.new_line();
    m_table.insert( cols, "t", m_t );
    m_table.insert( cols, "N", N );
    
    for( unsigned i=1; i<=m_NA; i++ )
    {
      for( unsigned j=1; j<=m_NK; j++ )
      {
        pcout << "t == " << m_t << endl;
        DoIter( m_Psi );
        m_t += m_dt;
      }

      N = Particle_Number(m_Psi);
      pcout << "N == " << N << endl;

      columns& cols = m_table.new_line();
      m_table.insert( cols, "t", m_t );
      m_table.insert( cols, "N", N );
      m_table.dump_2_file("err_log.csv");
      
      save( "step-" + to_string(m_t) + ".bin" );
    }
  }
  
  template<int dim>
  void MySolver<dim>::load_2( string filename )
  {
    make_grid();
    triangulation.load( filename.c_str() );
    setup_boundary_ids();

    setup_system(true);
    eta.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);;

    std::vector<LA::MPI::Vector*> x_system (2);
    x_system[0] = &m_Psi;
    x_system[1] = &m_Psi_0;

    parallel::distributed::SolutionTransfer<dim,LA::MPI::Vector> solution_transfer(dof_handler);
    solution_transfer.deserialize(x_system);
  }
} // end of namespace 

int main ( int argc, char *argv[] )
{
  int myrank;

  using namespace dealii;
  using namespace realtime_propagation;
  deallog.depth_console (0);

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv );
  {
    MySolver<DIMENSION> solver("params.xml");
    solver.run();
  }  
return EXIT_SUCCESS;
}
