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

#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>
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
#include "cxxopts.hpp"
#include "muParser.h"

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
        idx_shifts[0] = int(-( header.xMin / header.dx ));
        idx_shifts[1] = int(-( header.yMin / header.dy ));
        idx_shifts[2] = int(-( header.zMin / header.dz ));
        
        std::vector<double> deltas(3);
        deltas[0] = header.dx;
        deltas[1] = header.dy;
        deltas[2] = header.dz;

        int target_rank, target_position, target_ijk;

        MPI_Win win;
        MPI_Win_create( wf, 2*locN, sizeof(double), MPI_INFO_NULL, mpi_communicator, &win );
        MPI_Win_fence( MPI_MODE_NOPRECEDE, win);

        for( unsigned i=0; i<patches.size(); ++i )
        {
          const DataOutBase::Patch<dim,dim>  &patch = patches[i];

          for( unsigned l=0; l<GeometryInfo<dim>::vertices_per_cell; l++ )
          {
            const Point<dim> &pt = patch.vertices[l];

            bool skip = false;
            for( unsigned k=0; k<dim; ++k )
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

  template <int dim, int comp>
  class MySolver
  {
  public:
    MySolver( const std::string& );
    ~MySolver();

    void run ( string, const bool, const int=1 );

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
    AffineConstraints<double> constraints;

    LA::MPI::Vector m_Psi;

    bool m_root;

    ConditionalOStream pcout;

    double m_xmin, m_xmax;
    double m_ymin, m_ymax;
    double m_zmin, m_zmax;
    double m_mu;
    
    int m_rank;
    int m_sel; // the selected wave function
  };

  template <int dim, int comp>
  MySolver<dim,comp>::MySolver ( const std::string& xmlfilename ) 
    : 
    m_ph(xmlfilename),
    mpi_communicator (MPI_COMM_WORLD),
    triangulation (mpi_communicator, typename Triangulation<dim>::MeshSmoothing(Triangulation<dim>::smoothing_on_refinement|Triangulation<dim>::smoothing_on_coarsening)),
    fe (FE_Q<dim>(gl_degree_fe), 2),
    dof_handler (triangulation),
    m_root(Utilities::MPI::this_mpi_process(mpi_communicator) == 0),
    pcout (cout, m_root),
    m_sel(0)
  {
    try
    {
      m_xmin = m_ph.Get_Mesh("xrange",0);
      m_xmax = m_ph.Get_Mesh("xrange",1);
      m_ymin = m_ph.Get_Mesh("yrange",0);
      m_ymax = m_ph.Get_Mesh("yrange",1);
      m_zmin = m_ph.Get_Mesh("zrange",0);
      m_zmax = m_ph.Get_Mesh("zrange",1);
    }
    catch( const std::string info )
    {
      std::cerr << info << endl;
      exit(0);
    }    

    MPI_Comm_rank(mpi_communicator, &m_rank);
  }

  template <int dim, int comp>
  MySolver<dim,comp>::~MySolver ()
  {
    dof_handler.clear ();
  }

  template <int dim, int comp>
  void MySolver<dim,comp>::make_grid ()
  {
    Point<dim,double> pt1;
    Point<dim,double> pt2;

    double min[] = {m_xmin, m_ymin, m_zmin};
    double max[] = {m_xmax, m_ymax, m_zmax};

    for( int i=0; i<dim; ++i )
    {
      pt1(i) = min[i];
      pt2(i) = max[i];
    }

    GridGenerator::hyper_rectangle(triangulation, pt2, pt1);
  }
  
  template <int dim, int comp>
  void MySolver<dim,comp>::setup_system()
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

  template <int dim, int comp>
  void MySolver<dim,comp>::run( string filename, const bool bvtu, const int wf )
  {
    assert( wf > 0 && wf <= comp );

    m_sel = wf-1;
    load(filename);

    std::string new_filename = "atus2_" + to_string(wf) + "_" + filename ;

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


    if( bvtu )
    {
      regex pat {"*.bin"};

    }
/*
    vector<std::string> solution_names;

    DataOut<dim> data_out;
    data_out.attach_dof_handler (dof_handler);
    solution_names.push_back ("Re Psi");
    solution_names.push_back ("Im Psi");    
    data_out.add_data_vector (m_Psi, solution_names);
    data_out.add_data_vector (m_Psi, intensities);
    
    solution_names.clear();
    solution_names.push_back ("Re Psi_0");
    solution_names.push_back ("Im Psi_0");    
    data_out.add_data_vector (m_Psi_0, solution_names);
    data_out.add_data_vector (m_Psi_0, intensities2);
    data_out.build_patches ();
    data_out.write_vtu_in_parallel ( filename.c_str(), mpi_communicator );
*/

  }

  template<int dim, int comp>
  void MySolver<dim,comp>::load( string filename )
  {
    make_grid();
    triangulation.load( filename.c_str() );

    setup_system();

    vector<LA::MPI::Vector*> all_vec;

    for( int i=0; i<comp; ++i )
    {
      all_vec.push_back(new LA::MPI::Vector());
    }

    for( auto& i : all_vec )
    {
      i->reinit (locally_owned_dofs, mpi_communicator);
    }

    parallel::distributed::SolutionTransfer<dim,LA::MPI::Vector> solution_transfer(dof_handler);
    solution_transfer.deserialize(all_vec);

    constraints.distribute(*all_vec[m_sel]);
    m_Psi=*all_vec[m_sel];

    for( auto& i : all_vec )
    {
      delete i;
    }
  }
} // end of namespace 

int main ( int argc, char *argv[] )
{

try
{
  string str { "Cfinal.bin" };
  regex pat {"(.*)(bin)\b"};
  string format {"$1\n"};
  cout << regex_replace(str,pat,format,regex_constants::format_no_copy);
}
catch( const std::regex_error& e )
{
  cout << e.what() << endl;
}

  return EXIT_SUCCESS;

  using namespace dealii;
  deallog.depth_console (0);

  cxxopts::Options options("binC_to_atus2", "Converts complex deal.ii binary format to atus2 binary format.");

  options.add_options()
  ("p,params",  "parameter xml file", cxxopts::value<std::string>()->default_value("params.xml") )
  ;

  options.add_options("wave function index")
  ("i,iwf",  "select the i-th wave function for export (1,2,..)", cxxopts::value<int>()->default_value("1")  )
  ("N,nowf",  "total number of wave functions stored in the binary file (1,..)", cxxopts::value<int>()->default_value("1")  )
  ("e,vtu", "output addionally a vtu file", cxxopts::value<bool>()->default_value("false") )
  ("positional", "Positional arguments: these are the arguments that are entered without an option", cxxopts::value<std::vector<std::string>>())
  ("help","Print help")
  ;
  
  options.parse_positional({"positional"});
  auto result = options.parse(argc, argv);

  std::string bin_filename, params_filename;
  try
  {
    if (result.count("") == 0)
    {
      std::cout << options.help({""}) << std::endl;
      return EXIT_FAILURE;
    }

    if( result.count("positional") > 0 )
    {
      bin_filename = result["positional"].as<std::vector<std::string>>()[0]; 
    }
    else
    {        
      std::cout << options.help({""}) << std::endl;
      return EXIT_FAILURE;
    }
    params_filename = result["p"].as<std::string>();
  }
  catch (const cxxopts::OptionException& e)
  {
    std::cout << "error parsing options: " << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  MyParameterHandler params(params_filename);
  int dim=0;

  try
  {
    dim = int(params.Get_Mesh("DIM",0));
  }
  catch (mu::Parser::exception_type &e)
  {
    std::cout << "Message:  " << e.GetMsg() << "\n";
    std::cout << "Formula:  " << e.GetExpr() << "\n";
    std::cout << "Token:    " << e.GetToken() << "\n";
    std::cout << "Position: " << e.GetPos() << "\n";
    std::cout << "Errc:     " << e.GetCode() << "\n";
  }

  int nowf = result["nowf"].as<int>();
  int wf =  result["iwf"].as<int>();
  bool bvtu = result["vtu"].as<bool>();
  
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv );
  {
    if( dim == 2 && nowf == 1 )
    {
      HelperPrograms::MySolver<2,1> solver(params_filename);
      solver.run(bin_filename, bvtu);
    }
    if( dim == 3 && nowf == 1 )
    {
      HelperPrograms::MySolver<3,1> solver(params_filename);
      solver.run(bin_filename,bvtu);
    }
    if( dim == 2 && nowf == 2 )
    {
      HelperPrograms::MySolver<2,2> solver(params_filename);
      solver.run(bin_filename,bvtu,wf);
    }
    if( dim == 3 && nowf == 2 )
    {
      HelperPrograms::MySolver<3,2> solver(params_filename);
      solver.run(bin_filename,bvtu,wf);
    }
  }  
return EXIT_SUCCESS;
}