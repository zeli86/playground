/* * atus-pro testing - atus-pro testing playgroung
 * Copyright (C) 2017 Želimir Marojević <zelimir.marojevic@gmail.com>
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
    m_computing_timer.enter_section(__func__);

    Point<dim,double> pt1;
    Point<dim,double> pt2;

    double min[] = {m_xmin, m_ymin, m_zmin};
    double max[] = {m_xmax, m_ymax, m_zmax};

    for( int i=0; i<dim; i++ )
    {
      pt1(i) = min[i];
      pt2(i) = max[i];
    }

    GridGenerator::hyper_rectangle(triangulation, pt2, pt1);
    triangulation.refine_global(m_global_refinement);

    unsigned tmp1[2], tmp2[2];
    tmp1[0] = triangulation.n_cells();
    tmp1[1] = triangulation.n_active_cells();

    MPI_Allreduce( tmp1, tmp2, 2, MPI_UNSIGNED, MPI_SUM, mpi_communicator);

    m_total_no_cells = tmp2[0];
    m_total_no_active_cells = tmp2[1];
    m_computing_timer.exit_section();
  }
  
  template <int dim>
  void MySolver<dim>::make_grid_custom ()
  {
    m_computing_timer.enter_section(__func__);

    CPotential<dim> Potential_fct ( m_omega );
    
    Point<dim,double> pt1;
    Point<dim,double> pt2;

    double min[] = {m_xmin, m_ymin, m_zmin};
    double max[] = {m_xmax, m_ymax, m_zmax};

    for( int i=0; i<dim; i++ )
    {
      pt1(i) = min[i];
      pt2(i) = max[i];
    }
    
    GridGenerator::hyper_rectangle(triangulation, pt2, pt1);
    
#if DIMENSION==2
    triangulation.refine_global(5);    
    //triangulation.refine_global(6);    
    double isovalues[] = {34,32,30};
#endif
#if DIMENSION==3
    triangulation.refine_global(3);
    double isovalues[] = {18,16,14};
#endif
    
    for( unsigned step=0; step<sizeof(isovalues)/sizeof(double); step++ )
    {
      typename parallel::distributed::Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active(), endc = triangulation.end();
      for( ; cell!=endc; ++cell )
        for (unsigned v=0; v < GeometryInfo<dim>::vertices_per_cell; ++v )
        {
          Point<dim> p = cell->vertex(v);
          if( Potential_fct.value(p)  < isovalues[step] )
          {
            cell->set_refine_flag ();
            break;
          }
        }
      triangulation.execute_coarsening_and_refinement ();
    }

    unsigned tmp1[2], tmp2[2];
    tmp1[0] = triangulation.n_cells();
    tmp1[1] = triangulation.n_active_cells();

    MPI_Allreduce( tmp1, tmp2, 2, MPI_UNSIGNED, MPI_SUM, mpi_communicator);

    m_total_no_cells = tmp2[0];
    m_total_no_active_cells = tmp2[1];
    m_computing_timer.exit_section();
  }
