/* * atus-pro testing - atus-pro testing playgroung
 * Copyright (C) 2020 Želimir Marojević <zelimir.marojevic@gmail.com>
 *
 * This file is part of atus-pro testing.
 *
 * atus-pro testing is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * atus-pro testing is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with atus-pro testing.  If not, see <http://www.gnu.org/licenses/>.
 */


  template <int dim>
  void MySolver<dim>::make_grid ()
  {
    // TimerOutput::Scope timing_section(m_computing_timer, "");

    // Point<dim,double> pt1;
    // Point<dim,double> pt2;

    // // double min[] = {m_ph.Get_Mesh("xrange",0), m_ph.Get_Mesh("yrange",0), m_ph.Get_Mesh("zrange",0)};
    // // double max[] = {m_ph.Get_Mesh("xrange",1), m_ph.Get_Mesh("yrange",1), m_ph.Get_Mesh("zrange",1)};

    // for( int i=0; i<dim; i++ )
    // {
    //   pt1(i) = min[i];
    //   pt2(i) = max[i];
    // }

    // GridGenerator::hyper_rectangle(this->m_Triangulation, pt2, pt1);
    // this->m_Triangulation.refine_global(m_global_refinement);

    // unsigned tmp1[2], tmp2[2];
    // tmp1[0] = this->m_Triangulation.n_cells();
    // tmp1[1] = this->m_Triangulation.n_active_cells();

    // MPI_Allreduce( tmp1, tmp2, 2, MPI_UNSIGNED, MPI_SUM, mpi_communicator);

    // m_total_no_cells = tmp2[0];
    // m_total_no_active_cells = tmp2[1];
    
  }
  
  template <int dim>
  void MySolver<dim>::make_grid_custom ()
  {
//     TimerOutput::Scope timing_section(m_computing_timer, "");

//     CPotential<dim> Potential_fct ( m_omega );
    
//     Point<dim,double> pt1;
//     Point<dim,double> pt2;

//     double min[] = {m_ph.Get_Mesh("xrange",0), m_ph.Get_Mesh("yrange",0), m_ph.Get_Mesh("zrange",0)};
//     double max[] = {m_ph.Get_Mesh("xrange",1), m_ph.Get_Mesh("yrange",1), m_ph.Get_Mesh("zrange",1)};

//     for( int i=0; i<dim; i++ )
//     {
//       pt1(i) = min[i];
//       pt2(i) = max[i];
//     }
    
//     GridGenerator::hyper_rectangle(this->m_Triangulation, pt2, pt1);
    
// #if SPATIAL_DIM==2
//     this->m_Triangulation.refine_global(5);    
//     //triangulation.refine_global(6);    
//     double isovalues[] = {34,32,30};
// #endif
// #if SPATIAL_DIM==3
//     this->m_Triangulation.refine_global(3);
//     double isovalues[] = {18,16,14};
// #endif
    
//     for( unsigned step=0; step<sizeof(isovalues)/sizeof(double); step++ )
//     {
//       typename parallel::distributed::Triangulation<dim>::active_cell_iterator cell = this->m_Triangulation.begin_active(), endc = this->m_Triangulation.end();
//       for( ; cell!=endc; ++cell )
//         for (unsigned v=0; v < GeometryInfo<dim>::vertices_per_cell; ++v )
//         {
//           Point<dim> p = cell->vertex(v);
//           if( Potential_fct.value(p)  < isovalues[step] )
//           {
//             cell->set_refine_flag ();
//             break;
//           }
//         }
//       this->m_Triangulation.execute_coarsening_and_refinement ();
//     }

//     unsigned tmp1[2], tmp2[2];
//     tmp1[0] = this->m_Triangulation.n_cells();
//     tmp1[1] = this->m_Triangulation.n_active_cells();

//     MPI_Allreduce( tmp1, tmp2, 2, MPI_UNSIGNED, MPI_SUM, mpi_communicator);

//     m_total_no_cells = tmp2[0];
//     m_total_no_active_cells = tmp2[1];
    
  }
