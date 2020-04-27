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

#include <deal.II/lac/generic_linear_algebra.h>

namespace LA
{
  using namespace dealii::LinearAlgebraPETSc;
}

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>

#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_precondition.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/derivative_approximation.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/function_parser.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/solution_transfer.h>

#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <locale>
#include <limits>
#include <cmath>
#include <array>
#include <algorithm>

#include "global.h"
#include "mpi.h"
#include "functions.h"
#include "my_table.h"
//#include "ref_pt_list.h"
#include "MyParameterHandler.h"
#include "MyRealTools.h"
#include "muParser.h"

namespace BreedSolver
{
  using namespace std;
  using namespace dealii;

  #include "CBase.h"

  template <int dim>
  class MySolver : public CBase<dim,2>
  {
  public:
    explicit MySolver( const std::string );

    void run2b ();

    double m_I[8]; 
  protected:
    int DoIter( string="" );

    void make_grid();
    void make_grid_custom();

    void estimate_error( double& );

    bool solve();
    void compute_contributions();
    void compute_E_lin( LA::MPI::Vector&, double&, double&, double& );

    void output_vector ( LA::MPI::Vector&, string );

    MyTable m_table;
    MyTable m_results;

    double m_mu_punkt;

    using CBase<dim,2>::mpi_communicator;
    using CBase<dim,2>::m_root;
    using CBase<dim,2>::m_rank;
    using CBase<dim,2>::m_t;
    using CBase<dim,2>::m_ti;
    using CBase<dim,2>::m_N;
    using CBase<dim,2>::m_omega;
    using CBase<dim,2>::m_mu;
    using CBase<dim,2>::m_gs;
    using CBase<dim,2>::m_counter;
    using CBase<dim,2>::pcout;
    using CBase<dim,2>::m_computing_timer;
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
    using CBase<dim,2>::m_maxiter;
    using CBase<dim,2>::m_total_no_cells;
    using CBase<dim,2>::m_total_no_active_cells;
    using CBase<dim,2>::m_global_refinement;
  };

/*************************************************************************************************/
/**
 * Constructor
 */

  template <int dim>
  MySolver<dim>::MySolver ( const std::string xmlfilename ) : CBase<dim,2>(xmlfilename)
  {
  }
  
  #include "shared_1.h"
  #include "grid.h"  


  template <int dim>
  void MySolver<dim>::compute_contributions()
  {
    m_computing_timer.enter_section(__func__);

    this->update_workspace();
    
    CPotential<dim> Potential_fct ( m_omega );
    const QGauss<dim> quadrature_formula(this->m_FE.degree+1);
    FEValues<dim> fe_values (this->m_FE, quadrature_formula, update_gradients|update_values|update_JxW_values|update_quadrature_points);

    const unsigned dofs_per_cell = this->m_FE.dofs_per_cell;
    const unsigned n_q_points    = quadrature_formula.size();

    vector<double> Psi_1(n_q_points);
    vector<double> Psi_2(n_q_points);
    vector<Tensor<1,dim>> Psi_1_grad(n_q_points);
    vector<Tensor<1,dim>> Psi_2_grad(n_q_points);
    vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    double locint[8] = {};

    typename DoFHandler<dim>::active_cell_iterator cell = this->m_DOF_Handler.begin_active(), endc = this->m_DOF_Handler.end();
    for (; cell!=endc; ++cell)
    {
      if( cell->is_locally_owned() )
      {
        fe_values.reinit (cell);
        fe_values.get_function_values( this->m_Workspace[0], Psi_1 );
        fe_values.get_function_values( this->m_Workspace[1], Psi_2 );
        fe_values.get_function_gradients( this->m_Workspace[0], Psi_1_grad);
        fe_values.get_function_gradients( this->m_Workspace[1], Psi_2_grad);

        for ( unsigned qp=0; qp<n_q_points; ++qp )
        {
          double JxW = fe_values.JxW(qp);
          double p12 = Psi_1[qp]*Psi_2[qp];
          double p1q = Psi_1[qp]*Psi_1[qp];
          double p2q = Psi_2[qp]*Psi_2[qp];
          double Q = Potential_fct.value(fe_values.quadrature_point(qp)) - m_mu;

          locint[0] += JxW*p1q*p1q;
          locint[1] += JxW*p2q*p2q;
          locint[2] += JxW*p1q*p2q;
          locint[3] += JxW*p1q*p12;
          locint[4] += JxW*p2q*p12;
          locint[5] += JxW*(Psi_1_grad[qp]*Psi_1_grad[qp] + Q*p1q);
          locint[6] += JxW*(Psi_2_grad[qp]*Psi_2_grad[qp] + Q*p2q);
          locint[7] += JxW*(Psi_1_grad[qp]*Psi_2_grad[qp] + Q*p12);
        }  
      }
    }
    for( int i=0; i<5; ++i ) 
      locint[i] *= m_gs;

    MPI_Allreduce( locint, m_I, 8, MPI_DOUBLE, MPI_SUM, mpi_communicator);
    m_computing_timer.exit_section();
  }


  template <int dim>
  void MySolver<dim>::output_vector ( LA::MPI::Vector& vec, string filename )
  {
    m_computing_timer.enter_section(__func__);

    this->m_constraints.distribute(vec);
    this->m_Workspace[0]=vec;

    DataOut<dim> data_out;
    data_out.attach_this->m_DOF_Handler (this->m_DOF_Handler);
    data_out.add_data_vector (this->m_Workspace[0], "vec");
    data_out.build_patches (gl_subdivisions);
    data_out.write_vtu_in_parallel (filename.c_str(), mpi_communicator);

    m_computing_timer.exit_section();    
  }  

  template <int dim>
  int MySolver<dim>::DoIter ( string path )
  {
    int retval = Status::SUCCESS;
    
    m_table.clear();
    
    m_t[0] = m_ti;
    m_t[1] = m_ti;
    
    CPotential<dim> Potential( m_omega );

    this->do_linear_superposition();
    this->m_Workspace[0] = this->m_Psi_Ref;
    MyRealTools::MPI::AssembleRHS_L2gradient<dim>( this->m_DOF_Handler, this->m_FE, this->m_constraints, this->m_Workspace[0], Potential, m_mu, m_gs, m_res, this->m_System_RHS );
    m_res_old = m_res;
    m_counter=0;
    bool bempty;

    do
    {
      pcout << "--------------------------------------------------------------------------------" << endl;
      pcout << "- " << path << " - " << m_counter << endl;

      this->m_Workspace[0] = this->m_Psi_Ref;
      MyRealTools::MPI::AssembleSystem_Jacobian<dim>( this->m_DOF_Handler, this->m_FE, this->m_constraints, this->m_Workspace[0], Potential, m_mu, m_gs, this->m_System_Matrix);
      bool bsucc = solve();
      if( !bsucc ) { return Status::SINGULAR; }

      pcout << "m_Search_Direction.l2_norm() = " << this->m_Search_Direction.l2_norm() << endl;

/*
      this->m_Workspace[0] = m_Psi_1;
      this->m_Workspace[1] = this->m_Search_Direction;
      MyRealTools::MPI::orthonormalize(mpi_communicator, this->m_DOF_Handler, fe, this->m_constraints, this->m_Workspace[1], this->m_Workspace[0], this->m_Workspace_NG );
      this->m_Search_Direction = this->m_Workspace_NG;
*/
      double tau;
      this->m_Workspace[1] = this->m_Search_Direction;
      MyRealTools::MPI::compute_stepsize( mpi_communicator, this->m_DOF_Handler, this->m_FE, Potential, this->m_Workspace[0], this->m_Workspace[1], m_mu, m_gs, tau );

      this->m_Psi[1].add( tau*m_t[1]/fabs(m_t[1]), this->m_Search_Direction); 
      this->m_constraints.distribute(this->m_Psi[1]);

      bempty = this->find_ortho_min();

      if( bempty ) return Status::FAILED;
        
      this->do_linear_superposition();
      this->m_Workspace[0] = this->m_Psi_Ref;      
      MyRealTools::MPI::AssembleRHS_L2gradient<dim>( this->m_DOF_Handler, this->m_FE, this->m_constraints, this->m_Workspace[0], Potential, m_mu, m_gs, m_res, this->m_System_RHS );
      
      m_resp = m_res_old-m_res;
      m_res_old = m_res;

      if( m_counter % m_NA == 0 ) this->output_results(path);

      columns& cols = m_table.new_line();
      m_table.insert( cols, MyTable::COUNTER, double(m_counter) );
      m_table.insert( cols, MyTable::RES, m_res );
      m_table.insert( cols, MyTable::RESP, m_resp );
      m_table.insert( cols, MyTable::STEPSIZE, tau );
      m_table.insert( cols, MyTable::MU, m_mu );
      m_table.insert( cols, MyTable::GS, m_gs );
      m_table.insert( cols, MyTable::t1, m_t[0] );
      m_table.insert( cols, MyTable::t2, m_t[1] );
      m_table.insert( cols, MyTable::l2norm_t, this->l2norm_t() );

      m_counter++;
      if( m_root ) cout << m_table;
      if( m_res < m_epsilon[0] ) { retval=Status::SUCCESS; break; }
      if( this->l2norm_t() < 1e-4 ) { retval=Status::ZERO_SOL; break; }
      if( m_counter == m_maxiter ) { retval=Status::MAXITER; break; }
      if( isnan(m_res) ) { retval=Status::FAILED; break; }
   
    }while( true );
    
    // Standard Newton 
    do
    {
      pcout << "--------------------------------------------------------------------------------" << endl;
      pcout << "-- " << path << " - " << m_counter << endl;

      this->m_Workspace[0] = this->m_Psi_Ref;
      MyRealTools::MPI::AssembleSystem_Jacobian<dim>( this->m_DOF_Handler, this->m_FE, this->m_constraints, this->m_Workspace[0], Potential, m_mu, m_gs, this->m_System_Matrix);
      bool bsucc = solve();
      if( !bsucc ) { return Status::SINGULAR; }

      double tau;
      this->m_Workspace[1] = this->m_Search_Direction;
      MyRealTools::MPI::compute_stepsize( mpi_communicator, this->m_DOF_Handler, this->m_FE, Potential, this->m_Workspace[0], this->m_Workspace[1], m_mu, m_gs, tau );

      this->m_Psi_Ref.add( tau*m_t[1]/fabs(m_t[1]), this->m_Search_Direction); 
      this->m_constraints.distribute(this->m_Psi_Ref);

      this->m_Workspace[0] = this->m_Psi_Ref;
      MyRealTools::MPI::AssembleRHS_L2gradient<dim>( this->m_DOF_Handler, this->m_FE, this->m_constraints, this->m_Workspace[0], Potential, m_mu, m_gs, m_res, this->m_System_RHS );
      
      m_resp = m_res_old-m_res;
      m_res_old = m_res;

      if( m_counter % m_NA == 0 ) this->output_results(path);

      columns& cols = m_table.new_line();
      m_table.insert( cols, MyTable::COUNTER, double(m_counter) );
      m_table.insert( cols, MyTable::RES, m_res );
      m_table.insert( cols, MyTable::RESP, m_resp );
      m_table.insert( cols, MyTable::STEPSIZE, tau );
      m_table.insert( cols, MyTable::MU, m_mu );
      m_table.insert( cols, MyTable::GS, m_gs );
      m_table.insert( cols, MyTable::t1, m_t[0] );
      m_table.insert( cols, MyTable::t2, m_t[1] );
      m_table.insert( cols, MyTable::l2norm_t, this->l2norm_t() );
      m_table.insert( cols, MyTable::PARTICLE_NUMBER, m_N );

      m_counter++;

      if( m_root ) cout << m_table;
      if( m_res < m_epsilon[1] ) { retval=Status::SUCCESS; break; }
      if( m_resp < 0 ) { retval=Status::SUCCESS; break; }
      if( m_counter == m_maxiter ) { retval=Status::MAXITER; break; }
    }while( true );
    
    this->m_Workspace[0] = this->m_Psi_Ref;
    m_N = MyRealTools::MPI::Particle_Number( mpi_communicator, this->m_DOF_Handler, this->m_FE, this->m_Workspace[0] );
    
    if( m_N < 1e-5 ) retval = Status::ZERO_SOL;
    
    string filename = path + "log.csv";
    if( m_root ) m_table.dump_2_file(filename);
    
    return retval;
  }

  
  template <int dim>
  void MySolver<dim>::run2b ()
  {
    string path;
    char shellcmd[255];
    double T, N, W;
    int status;

    make_grid_custom();
    this->setup_system();

    CEigenfunctions<dim> Ef1( m_QN1, m_omega );
    //CEigenfunctions<dim> Ef2( m_QN2, m_omega );
    
    VectorTools::interpolate (this->m_DOF_Handler, Ef1, this->m_Psi[0] );

    this->m_Workspace[0] = this->m_Psi[0];
    this->m_Psi[0] *= 1.0/sqrt(MyRealTools::MPI::Particle_Number( mpi_communicator, this->m_DOF_Handler, this->m_FE, this->m_Workspace[0] ));
    this->m_Psi[1] = 0; 

    compute_E_lin( this->m_Psi[0], T, N, W );
    double m_mu_0 = T/N;
    m_mu = ceil(10.0*m_mu_0)/10.0 + m_gs/fabs(m_gs)*m_dmu;

    pcout << "T = " << T << endl;
    pcout << "N = " << N << endl;
    pcout << "W = " << W << endl;
    pcout << "m_mu = " << m_mu << endl;

    //output_guess();
    m_results.clear();
    for( int i=0; i<m_Ndmu; ++i )
    {
      sprintf( shellcmd, "mkdir %.4d/", i );
      if( m_root ) system(shellcmd);
      sprintf( shellcmd, "%.4d/", i );
      path = shellcmd;

      // nehari
      // sqrt((m_mu*N-T)/(m_gs*W)); if this->m_Psi[1] == 0
      // sqrt((m_mu*N-T)/(4.0*m_gs*W)); if this->m_Psi[1] == m_Psi_1
      m_ti = sqrt((m_mu*N-T)/(m_gs*W));
      
      status = DoIter(path);

      columns& cols = m_results.new_line();
      m_results.insert( cols, MyTable::MU, m_mu );
      m_results.insert( cols, MyTable::GS, m_gs );
      m_results.insert( cols, MyTable::PARTICLE_NUMBER, m_N );
      m_results.insert( cols, MyTable::COUNTER, double(m_counter) );
      m_results.insert( cols, MyTable::STATUS, double(status) );

      if( status == Status::SUCCESS )
      {
        estimate_error(m_final_error);

        this->save( path + "final.bin" );

        m_ph.SaveXMLFile( path + "params.xml" );       
        
        this->output_results(path,"final");
        this->dump_info_xml(path);
        this->m_Psi[0]  = this->m_Psi_Ref;
        this->m_Psi[1] = 0;
      }
      else if( status == Status::MAXITER || status == Status::SINGULAR )
      {
        this->m_Psi_Ref = this->m_Psi[0];
        this->m_Psi[1] = 0; 
      }
      else
      {
        //this->m_Psi[1] = 0; 
        break;
      }
      m_mu += m_gs/fabs(m_gs)*m_dmu;
      compute_E_lin( this->m_Psi_Ref, T, N, W ); // TODO: kommentier mich aus, falls ich kein nehari reset habe
      pcout << "T = " << T << endl;
      pcout << "N = " << N << endl;
      pcout << "W = " << W << endl;      
    }
    if( m_root ) m_results.dump_2_file( "results.csv" );
  }
} // end of namespace 

int main ( int argc, char *argv[] )
{
  using namespace dealii;
  //deallog.depth_console (2);

  MyParameterHandler params("params.xml");
  int dim=0;

  try
  {
    dim = int(params.Get_Mesh("DIM",0));
  }
  catch (mu::Parser::exception_type &e)
  {
    cout << "Message:  " << e.GetMsg() << "\n";
    cout << "Formula:  " << e.GetExpr() << "\n";
    cout << "Token:    " << e.GetToken() << "\n";
    cout << "Position: " << e.GetPos() << "\n";
    cout << "Errc:     " << e.GetCode() << "\n";
  }

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv );
  {
    switch(dim)
    {
      case 2: { BreedSolver::MySolver<2> solver("params.xml");
                solver.run2b();
                break; }
      case 3: { BreedSolver::MySolver<3> solver("params.xml");
                solver.run2b();
                break; }
      default: cout << "You have found a new dimension!" << endl;
    }
  }
return EXIT_SUCCESS;
}
