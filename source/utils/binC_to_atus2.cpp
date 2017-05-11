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

#include <deal.II/lac/generic_linear_algebra.h>

namespace LA
{
  using namespace dealii::LinearAlgebraPETSc;
}

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/constraint_matrix.h>


#include <deal.II/lac/petsc_parallel_sparse_matrix.h>
#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_precondition.h>
#include <deal.II/lac/sparsity_tools.h>

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

#include <deal.II/base/utilities.h>
#include <deal.II/base/aligned_vector.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/function_parser.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/data_out_base.h>
#include <deal.II/base/table.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/solution_transfer.h>

#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <algorithm>

#include "global.h"
#include "mpi.h"
#include "MyParameterHandler.h"
#include "my_table.h"
#include "atus2.h"
#include "anyoption.h"

namespace HelperPrograms
{
  using namespace std;
  using namespace dealii;

  template <int dim>
  class MyDataOut : public DataOut<dim>
  {
    public:
      MyDataOut(){};
      virtual ~MyDataOut(){};

      void write_atus2( const std::string& filename, MPI_Comm& mpi_communicator, const std::vector<double>& data )
      {
        const std::vector<DataOutBase::Patch<dim,dim>> &patches = this->get_patches();
        const DataOutBase::Patch<dim,dim> &patch0 = patches[0];

        const unsigned subdiv = patch0.n_subdivisions;

        generic_header header = {};

        int rank, commsize;           
        MPI_Comm_rank(mpi_communicator, &rank);
        MPI_Comm_size(mpi_communicator, &commsize);
        
        const unsigned lev = this->triangulation->n_global_levels ();

        unsigned N = unsigned(pow(2,double(lev-1)));
        N = N / subdiv;

        header.nself = sizeof(generic_header);   
        header.nDatatyp = 2*sizeof(double); 
        //header.nself_and_data; 
        header.nDims = dim;
        header.nDimX = N;
        header.nDimY = ( dim == 2 || dim == 3 ) ? N : 1;
        header.nDimZ = ( dim == 3 ) ? N : 1;
        header.bAtom = true;
        header.bComplex = true;
        header.xMin = data[0];
        header.xMax = data[1];
        header.yMin = data[2];
        header.yMax = data[3];
        header.zMin = data[4];
        header.zMax = data[5];
        header.dx = fabs(data[1]-data[0])/double(N);
        header.dy = fabs(data[3]-data[2])/double(N);
        header.dz = fabs(data[5]-data[4])/double(N);
        header.dkx = 2*M_PI/fabs(data[1]-data[0]);
        header.dky = 2*M_PI/fabs(data[3]-data[2]);
        header.dkz = 2*M_PI/fabs(data[5]-data[4]);
        header.dt = 0.001;
 
        ptrdiff_t locN;

        switch(dim)
        {
          case 2: locN = N*N/commsize; break;
          case 3: locN = N*N*N/commsize; break;
        }

        double * wf = nullptr;
        MPI_Alloc_mem( 2*locN*sizeof(double), MPI_INFO_NULL, &wf );
        memset( wf, 0, 2*locN*sizeof(double));

        std::vector<int> coord_idx(3);
        std::vector<int> idx_shifts(3);
        idx_shifts[0] = ( header.xMin < 0 ) ? N/2 : 0;
        idx_shifts[1] = ( header.yMin < 0 ) ? N/2 : 0;
        idx_shifts[2] = ( header.zMin < 0 ) ? N/2 : 0; 
        
        std::vector<double> deltas(3);
        deltas[0] = header.dx;
        deltas[1] = header.dy;
        deltas[2] = header.dz;

        int target_rank, target_position, target_ijk;

        MPI_Win win;
        MPI_Win_create( wf, 2*locN, sizeof(double), MPI_INFO_NULL, mpi_communicator, &win );
        MPI_Win_fence( MPI_MODE_NOPRECEDE, win);

        for( unsigned i=0; i<patches.size(); i++ )
        {
          const DataOutBase::Patch<dim,dim>  &patch = patches[i];

          for( unsigned l=0; l<GeometryInfo<dim>::vertices_per_cell; l++ )
          {
            const Point<dim> &pt = patch.vertices[l];

            bool skip = false;
            for( unsigned k=0; k<dim; k++ )
            {
               coord_idx[k] = int(pt[k] / deltas[k]) + idx_shifts[k];
               skip = skip || (coord_idx[k]==N); 
            }
            if( skip ) continue;

            switch(dim)
            {
              case 2: target_ijk = coord_idx[1] + N * coord_idx[0]; break;
              case 3: target_ijk = coord_idx[2] + N * (coord_idx[1] + N * coord_idx[0]); break;
            }

            target_rank = target_ijk / locN;
            target_position = target_ijk-target_rank*locN;
            
            double data[] = {patch.data(0,l), patch.data(1,l)};
            MPI_Put( data, 2, MPI_DOUBLE, target_rank, 2*target_position, 2, MPI_DOUBLE, win );
          }
        }

        MPI_Win_fence( (MPI_MODE_NOSTORE|MPI_MODE_NOSUCCEED), win);
        MPI_Win_free(&win);

        MPI_Status status;
        MPI_File fh;
     
        MPI_Offset offset = sizeof(generic_header) + 2 * sizeof(double) * locN * ptrdiff_t(rank);
      
        MPI_File_open( mpi_communicator, const_cast<char*>(filename.c_str()), MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fh );

        if( rank == 0 )
	        MPI_File_write( fh, &header, sizeof(generic_header), MPI_BYTE, &status );
        MPI_File_write_at( fh, offset, wf, 2*locN, MPI_DOUBLE, MPI_STATUS_IGNORE );
        MPI_File_close( &fh );  

        MPI_Free_mem(wf);
      }       
  };

  template <int dim>
  class MySolver
  {
  public:
    MySolver( const std::string& );
    ~MySolver();

    void run ( string );

  protected:
    void make_grid();
    void setup_system();
    void setup_boundary_ids();
  
    void load( string );

    MyParameterHandler m_ph;
    MPI_Comm mpi_communicator;
    parallel::distributed::Triangulation<dim> triangulation;
    FESystem<dim> fe;
    DoFHandler<dim> dof_handler;
    IndexSet locally_owned_dofs;
    IndexSet locally_relevant_dofs;
    ConstraintMatrix constraints;

    LA::MPI::Vector m_Psi;

    ConditionalOStream pcout;

    bool m_root;

    double m_xmin, m_xmax;
    double m_ymin, m_ymax;
    double m_zmin, m_zmax;
    double m_mu;
    double m_gs;
    vector<double> m_omega;
    
    int m_rank;
  };

  template <int dim>
  MySolver<dim>::MySolver ( const std::string& xmlfilename ) 
    : 
    m_ph(xmlfilename),
    mpi_communicator (MPI_COMM_WORLD),
    triangulation (mpi_communicator, typename Triangulation<dim>::MeshSmoothing(Triangulation<dim>::smoothing_on_refinement|Triangulation<dim>::smoothing_on_coarsening)),
    fe (FE_Q<dim>(gl_degree_fe), 2),
    dof_handler (triangulation),
    pcout (cout, (Utilities::MPI::this_mpi_process(mpi_communicator) == 0))
  {
    m_omega = m_ph.Get_Physics("omega");
    m_gs = m_ph.Get_Physics("gs_1",0);

    m_xmin = m_ph.Get_Mesh("xrange",0);
    m_xmax = m_ph.Get_Mesh("xrange",1);
    m_ymin = m_ph.Get_Mesh("yrange",0);
    m_ymax = m_ph.Get_Mesh("yrange",1);
    m_zmin = m_ph.Get_Mesh("zrange",0);
    m_zmax = m_ph.Get_Mesh("zrange",1);
  
    m_root = (Utilities::MPI::this_mpi_process(mpi_communicator) == 0);
    MPI_Comm_rank(mpi_communicator, &m_rank);
  }

  template <int dim>
  MySolver<dim>::~MySolver ()
  {
    dof_handler.clear ();
  }

  template <int dim>
  void MySolver<dim>::make_grid ()
  {
#if DIMENSION==2
    Point<dim,double> pt1( m_xmin, m_ymin );
    Point<dim,double> pt2( m_xmax, m_ymax );
#endif
#if DIMENSION==3
    Point<dim,double> pt1( m_xmin, m_ymin, m_zmin );
    Point<dim,double> pt2( m_xmax, m_ymax, m_zmax );
#endif

    GridGenerator::hyper_rectangle(triangulation, pt2, pt1);
    triangulation.refine_global(1);
  }
  
  template <int dim>
  void MySolver<dim>::setup_system()
  {
    dof_handler.distribute_dofs (fe);

    locally_owned_dofs = dof_handler.locally_owned_dofs ();
    DoFTools::extract_locally_relevant_dofs (dof_handler, locally_relevant_dofs);

    m_Psi.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);

    constraints.clear ();
    constraints.reinit (locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints (dof_handler, constraints);
    constraints.close ();
  }

  template <int dim>
  void MySolver<dim>::run( string filename )
  {
    load(filename);

    std::string new_filename = "atus2_" + filename ;

    std::vector<double> add_data(6);
    add_data[0] = m_xmin;
    add_data[1] = m_xmax;
    add_data[2] = m_ymin;
    add_data[3] = m_ymax;
    add_data[4] = m_zmin;
    add_data[5] = m_zmax;

    MyDataOut<dim> data_out;
    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (m_Psi, "Psi");
    data_out.build_patches ();
    data_out.write_atus2 ( new_filename, mpi_communicator, add_data );
  }

  template<int dim>
  void MySolver<dim>::load( string filename )
  {
    make_grid();
    triangulation.load( filename.c_str() );
    
    setup_system();

    LA::MPI::Vector tmp;
    tmp.reinit (locally_owned_dofs, mpi_communicator);

    parallel::distributed::SolutionTransfer<dim,LA::MPI::Vector> solution_transfer(dof_handler);
    solution_transfer.deserialize(tmp);   
    constraints.distribute(tmp);
    m_Psi=tmp;
  }
} // end of namespace 

int main ( int argc, char *argv[] )
{
  using namespace dealii;
  deallog.depth_console (0);

  std::string filename;

  AnyOption * opt = new AnyOption();
  int dim;

  opt->noPOSIX(); 
  //opt->setVerbose();
  //opt->autoUsagePrint(true); 

  opt->addUsage( "" );
  opt->addUsage( "Usage: binR_to_atus2 [options] filename" );
  opt->addUsage( "" );
  opt->addUsage( " -h --help	Prints this help " );
  opt->addUsage( " -dim val (2 or 3)" );
  opt->addUsage( "" );
  opt->setFlag(  "help", 'h' );   

  opt->processCommandArgs( argc, argv );

  if( opt->getFlag( "help" ) || opt->getFlag( 'h' ) ) opt->printUsage();

  if( opt->getValue("-dim") != NULL ) 
  { 
    dim = atof(opt->getValue("-dim"));  
  };
  
  if( opt->getArgc() != 0 ) filename = opt->getArgv(0);
  else opt->printUsage();
  delete opt; 

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv );
  {
    if( dim == 2 )
    {
      std::string filename = argv[1];
      HelperPrograms::MySolver<2> solver("params.xml");
      solver.run(filename);
    }
    if( dim == 3 )
    {
      std::string filename = argv[1];
      HelperPrograms::MySolver<3> solver("params.xml");
      solver.run(filename);
    }
  }  
return EXIT_SUCCESS;
}