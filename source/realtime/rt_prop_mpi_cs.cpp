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
#include <deal.II/grid/grid_tools.h>

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
#include <deal.II/base/timer.h>

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
    void run_restart();
    double Particle_Number( LA::MPI::Vector& );
    void Expectation_value_momentum( LA::MPI::Vector&, double* );
    void Expectation_value_position( LA::MPI::Vector&, double* );

  protected:
    void make_grid();
    void setup_system( const bool );
    void setup_boundary_ids();
    void assemble_system( LA::MPI::Vector& );
    void assemble_rhs( LA::MPI::Vector& );
    void compute_visibility_and_difference( double&, double& );

    void DoIter( LA::MPI::Vector& );
    void solve();
    void output_vtu ();
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
    double m_res;
    vector<double> m_omega;

    double m_xmin;
    double m_xmax;
    double m_ymin;
    double m_ymax;
    double m_zmin;
    double m_zmax;

    int m_rank;
    unsigned m_NA;
    unsigned m_NK;
    unsigned m_QN1[3];

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

        for( unsigned qp=0; qp<n_q_points; qp++ )
        {
          double JxW = fe_values.JxW(qp)*fabs(fe_values.quadrature_point(qp)[1]);
          double pot = Potential.value(fe_values.quadrature_point(qp));

          double sum_re = Psi[qp][0]+Psi_t[qp][0];
          double sum_im = Psi[qp][1]+Psi_t[qp][1];
          double sum_req = sum_re*sum_re;
          double sum_imq = sum_im*sum_im;
          double tmp1a = fak8*(sum_req + 3*sum_imq);
          double tmp1b = fak8*(sum_imq + 3*sum_req);
          double tmp2 = fak4*sum_re*sum_im;

          for( unsigned i=0; i<dofs_per_cell; i++ )
          {
            for( unsigned j=0; j<dofs_per_cell; j++ )
            {
              cell_matrix(i,j) += JxW*(1.0-tmp2)*fe_values[rt].value(i,qp)*fe_values[rt].value(j,qp);
              cell_matrix(i,j) += JxW*(1.0+tmp2)*fe_values[it].value(i,qp)*fe_values[it].value(j,qp);
              cell_matrix(i,j) -= JxW*tmp1a*fe_values[rt].value(i,qp)*fe_values[it].value(j,qp);
              cell_matrix(i,j) -= JxW*fak2*(fe_values[rt].gradient(i,qp)*fe_values[it].gradient(j,qp) + pot*fe_values[rt].value(i,qp)*fe_values[it].value(j,qp));
              cell_matrix(i,j) += JxW*tmp1b*fe_values[it].value(i,qp)*fe_values[rt].value(j,qp);
              cell_matrix(i,j) += JxW*fak2*(fe_values[it].gradient(i,qp)*fe_values[rt].gradient(j,qp) + pot*fe_values[it].value(i,qp)*fe_values[rt].value(j,qp));
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

        for( unsigned int qp=0; qp<n_q_points; qp++ )
        {
          JxW = fe_values.JxW(qp)*fabs(fe_values.quadrature_point(qp)[1]);
          pot = Potential.value(fe_values.quadrature_point(qp));

          sum_re = Psi[qp][0]+Psi_t[qp][0];
          sum_im = Psi[qp][1]+Psi_t[qp][1];
          tmp1 = fak8*(sum_re*sum_re + sum_im*sum_im);

          for (unsigned int i=0; i<dofs_per_cell; i++ )
          {
            cell_rhs(i) -= JxW*fak2*((Psi_grad[qp][1]+Psi_t_grad[qp][1])*fe_values[rt].gradient(i,qp) + pot*sum_im*fe_values[rt].value(i,qp));
            cell_rhs(i) += JxW*fak2*((Psi_grad[qp][0]+Psi_t_grad[qp][0])*fe_values[it].gradient(i,qp) + pot*sum_re*fe_values[it].value(i,qp));
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

  /*
  template <int dim>
  void MySolver<dim>::solve ()
  {
    m_computing_timer.enter_section(__func__);
    SolverControl solver_control;

    PETScWrappers::SparseDirectMUMPS solver(solver_control, mpi_communicator);
    solver.set_symmetric_mode(false);
    solver.solve(system_matrix, newton_update, system_rhs);

    constraints.distribute (newton_update);
    m_computing_timer.exit_section();
  }
  */
  
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

  template<int dim>
  void MySolver<dim>::DoIter( LA::MPI::Vector& vec )
  {
    m_Psi_t = vec;
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
    while( m_res > 1e-15 );

    vec = m_Psi_t;
  }

  template <int dim>
  void MySolver<dim>::run()
  {
    double N, vis, d;

    pugi::xml_document doc;
    if (!doc.load_file("info.xml")) throw;
    
    m_mu = stod(doc.child( "INFO" ).child( "MU" ).child_value());
    pcout << "m_mu = " << m_mu << endl;

    load( "Cfinal.bin" );
    
    double min_cell_diameter = GridTools::minimal_cell_diameter(triangulation);
    double max_cell_diameter = GridTools::maximal_cell_diameter(triangulation);
    pcout << "min_cell_diameter = " << min_cell_diameter << "\n";
    pcout << "dt/dx^2 == " << m_dt/(min_cell_diameter*min_cell_diameter) << endl;     
    pcout << "max_cell_diameter = " << max_cell_diameter << "\n";
    pcout << "dt/dx^2 == " << m_dt/(max_cell_diameter*max_cell_diameter) << endl;        

    m_Psi_0=m_Psi;

//     output_results("");

    N = Particle_Number(m_Psi);

    pcout << "t == " << m_t << endl;
    pcout << "N == " << N << endl;

    for( unsigned i=1; i<=m_NA; i++ )
    {
      for( unsigned j=1; j<=m_NK; j++ )
      {
        pcout << "t == " << m_t << endl;
        DoIter( m_Psi );
        m_t += m_dt;

        compute_visibility_and_difference( vis, d );
        if( m_root )
        {
          ofstream err( "err_prop.txt", ofstream::out | ofstream::app ); 
          err << m_t << "\t" << vis << "\t" << d << "\n";
          err.close();
        }
      }

      pcout << "N == " << N << endl;
      pcout << "vis == " << vis << endl;
      pcout << "d == " << d << endl;

      output_vtu();
      save( "step-" + to_string(m_t) + ".bin" );
    }
  }

  template <int dim>
  void MySolver<dim>::run_restart()
  {
    double N, vis, d;

    pugi::xml_document doc;
    if (!doc.load_file("info.xml")) throw;
    
    m_mu = stod(doc.child( "INFO" ).child( "MU" ).child_value());
    pcout << "m_mu = " << m_mu << endl;

    string filename;
    Find_last_file( ".", m_t, filename );
    pcout << "filename == " << filename << endl;
    pcout << "m_t == " << m_t << endl;
    if(m_root) m_table.load( "rt_prop_log.csv" );

    load_2( filename );

    system_rhs = m_Psi;
    constraints.distribute(system_rhs);
    m_Psi = system_rhs;

    system_rhs = m_Psi_0;
    constraints.distribute(system_rhs);
    m_Psi_0 = system_rhs;

    N = Particle_Number(m_Psi);
    for( unsigned i=1; i<=m_NA; i++ )
    {
      for( unsigned j=1; j<=m_NK; j++ )
      {
        pcout << "t == " << m_t << endl;
        DoIter( m_Psi );
        m_t += m_dt;
      
        compute_visibility_and_difference( vis, d );
        if( m_root )
        {
          ofstream err( "err_prop.txt", ofstream::out | ofstream::app ); 
          err << m_t << "\t" << vis << "\t" << d << "\n";
          err.close();
        }
      }
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

    std::vector<LA::MPI::Vector*> x_system (2);
    x_system[0] = &m_Psi;
    x_system[1] = &m_Psi_0;

    parallel::distributed::SolutionTransfer<dim,LA::MPI::Vector> solution_transfer(dof_handler);
    solution_transfer.deserialize(x_system);
  }
  
  template <int dim>
  void MySolver<dim>::output_vtu ()
  {
    m_computing_timer.enter_section(__func__);

    string filename;

    m_workspace_1=m_Psi;
    
    DataOut<dim> data_out;
    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (m_workspace_1, "Psi");
    data_out.build_patches ();

    filename = "step-" + to_string(m_t) + ".vtu";
    data_out.write_vtu_in_parallel (filename.c_str(), mpi_communicator);

    m_computing_timer.exit_section();    
  }

  template <int dim>
  void MySolver<dim>::compute_visibility_and_difference( double& retval1, double& retval2 )
  {
    retval1=0;
    retval2=0;

    m_workspace_1=m_Psi_0;    
    constraints.distribute(m_Psi);
    m_workspace_2=m_Psi;
    
    double retval[4]={}, tmp1[4]={}, Psi_neu[2]={}, JxW;
    
    const double re = cos(-m_mu*m_t);
    const double im = sin(-m_mu*m_t);
    
    const QGauss<dim>  quadrature_formula(fe.degree+1);
    FEValues<dim> fe_values (fe, quadrature_formula, update_values|update_quadrature_points|update_JxW_values);

    const unsigned int n_q_points = quadrature_formula.size();
    vector<Vector<double>> Psi_0(n_q_points,Vector<double>(2));
    vector<Vector<double>> Psi(n_q_points,Vector<double>(2));

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
      if( cell->is_locally_owned() )
      {
        fe_values.reinit (cell);
        fe_values.get_function_values( m_workspace_1, Psi_0 );
        fe_values.get_function_values( m_workspace_2, Psi );
        for( unsigned qp=0; qp<n_q_points; qp++ )
        {
          JxW = fe_values.JxW(qp)*fabs(fe_values.quadrature_point(qp)[1]);

          Psi_neu[0] = (Psi_0[qp][0]*re - Psi_0[qp][1]*im) - Psi[qp][0];
          Psi_neu[1] = (Psi_0[qp][0]*im + Psi_0[qp][1]*re) - Psi[qp][1];
 
          tmp1[0] += JxW*(Psi_0[qp][0]*Psi[qp][0]+Psi_0[qp][1]*Psi[qp][1]);
          tmp1[1] += JxW*(Psi_0[qp][0]*Psi[qp][1]-Psi_0[qp][1]*Psi[qp][0]);
          tmp1[2] += JxW*(Psi_0[qp][0]*Psi_0[qp][0]+Psi_0[qp][1]*Psi_0[qp][1]);
          tmp1[3] += JxW*(Psi_neu[0]*Psi_neu[0]+Psi_neu[1]*Psi_neu[1]);	 
        }
      }
    }
    MPI_Allreduce( tmp1, retval, 4, MPI_DOUBLE, MPI_SUM, mpi_communicator);
    retval1 = sqrt(retval[0]*retval[0] + retval[1]*retval[1])/retval[2];
    retval2 = retval[3];    
  }  
} // end of namespace

int main ( int argc, char *argv[] )
{
  using namespace dealii;
  deallog.depth_console (0);

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv );
  {
    realtime_propagation::MySolver<DIMENSION> solver("params.xml");
    solver.run();
    //solver.run_restart();
  }
return EXIT_SUCCESS;
}
