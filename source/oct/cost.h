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

  template <int dim, int no_time_steps, int no_lam>
  void MySolver<dim,no_time_steps,no_lam>::compute_all_lambdas_t()
  {
    m_computing_timer.enter_section(__func__);

    m_computing_timer.exit_section();
  }
  
  template <int dim, int no_time_steps, int no_lam>
  double MySolver<dim,no_time_steps,no_lam>::compute_costfunction()
  {
    double retval=0;
    return retval;
  }