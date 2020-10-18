/** atus-pro testing - atus-pro testing playgroung
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


#pragma once

#define STR1(x) #x
#define STR2(x) STR1(x)

#include <type_traits>

#include "ref_pt_list.h"
#include "nlopt.hpp"
#include "pugixml.hpp"

template <int dim>
class MySolver;


// template <int dim>
// double myfunc(const std::vector<double>& t, std::vector<double>& grad, void* data)
// {
//   MySolver<dim>* sol = reinterpret_cast<MySolver<dim>*>(data);

//   const double retval = 0.25 * sol->m_I[0] * pow(t[0], 4)
//                         + 0.5 * sol->m_I[5] * pow(t[0], 2)
//                         + 0.25 * sol->m_I[1] * pow(t[1], 4)
//                         + t[0] * sol->m_I[4] * pow(t[1], 3)
//                         + 1.5 * sol->m_I[2] * pow(t[0], 2) * pow(t[1], 2)
//                         + 0.5 * sol->m_I[6] * pow(t[1], 2)
//                         + pow(t[0], 3) * sol->m_I[3] * t[1]
//                         + t[0] * sol->m_I[7] * t[1];

//   if (!grad.empty())
//   {
//     grad[0] = sol->m_I[0] * pow(t[0], 3)
//               + sol->m_I[5] * t[0]
//               //+ 0.25 * sol->m_I[1] * pow(t[1], 4)
//               + sol->m_I[4] * pow(t[1], 3)
//               + 3 * sol->m_I[2] * t[0] * pow(t[1], 2)
//               //+ 0.5 * sol->m_I[6] * pow(t[1], 2)
//               + 3 * pow(t[0], 2) * sol->m_I[3] * t[1]
//               + sol->m_I[7] * t[1];

//     grad[1] = //0.25 * sol->m_I[0] * pow(t[0], 4)
//       //+ 0.5 * sol->m_I[5] * pow(t[0], 2)
//       + sol->m_I[1] * pow(t[1], 3)
//       + 3 * t[0] * sol->m_I[4] * pow(t[1], 2)
//       + 3 * sol->m_I[2] * pow(t[0], 2) * t[1]
//       + sol->m_I[6] * t[1]
//       + pow(t[0], 3) * sol->m_I[3]
//       + t[0] * sol->m_I[7];
//   }
//   return retval;
// }

template <int dim>
double myfunc(const std::vector<double>& t, std::vector<double>& grad, void* data)
{
  MySolver<dim>* sol = reinterpret_cast<MySolver<dim>*>(data);

  const double t0q = t[0] * t[0];
  const double t1q = t[1] * t[1];

  const auto& coeffs = sol->m_coeffs;

  const double retval = coeffs.at("t0^2") * t0q
                        + coeffs.at("t0_t1") * t[0] * t[1]
                        + coeffs.at("t1^2") * t1q
                        + coeffs.at("t0^4") * t0q * t0q
                        + coeffs.at("t0^1_t1^3") * t[0] * t[1] * t1q
                        + coeffs.at("t0^2_t1^2") * t0q * t1q
                        + coeffs.at("t0^3_t1^1") * t[0] * t[1] * t0q
                        + coeffs.at("t1^4") * t1q * t1q;

  if (!grad.empty())
  {
    grad[0] = 2 * coeffs.at("t0^2") * t[0]
              + coeffs.at("t0_t1") * t[1]
              + 4 * coeffs.at("t0^4") * t0q * t[0]
              + coeffs.at("t0^1_t1^3") * t[1] * t1q
              + 2 * coeffs.at("t0^2_t1^2") * t[0] * t1q
              + 3 * coeffs.at("t0^3_t1^1")  * t[1] * t0q;

    grad[1] =  coeffs.at("t0_t1") * t[0]
               + 2 * coeffs.at("t1^2") * t[1]
               + 3 * coeffs.at("t0^1_t1^3") * t[0] * t1q
               + 2 * coeffs.at("t0^2_t1^2") * t0q * t[1]
               + coeffs.at("t0^3_t1^1") * t[0] * t0q
               + 4 * coeffs.at("t1^4") * t1q * t[1];
  }

  return retval;
}

enum Status { SUCCESS, FAILED, ZERO_SOL, SLOW_CONV, MAXITER, SINGULAR };


template <int dim, int N, class T>
class cBaseMPI
{
public:
  explicit cBaseMPI(const std::string&, T& );
  virtual ~cBaseMPI()
  {
    m_DOF_Handler.clear();
  };

  int find_ortho_min();

  void dump_info_xml(const string = "");

  double l2norm_t();

  virtual void compute_contributions() = 0;

  void setup_system();

  void do_linear_superposition();

  void update_workspace();

  void save(const std::string&);

  void output_guess();

  void output_results(const std::string&, std::string = "step");

  MPI_Comm mpi_communicator;
protected:

  typedef typename dealii::LinearAlgebraPETSc::MPI::Vector aVector;
  typedef typename dealii::LinearAlgebraPETSc::MPI::SparseMatrix aSparseMatrix;
  typedef typename dealii::parallel::distributed::Triangulation<dim> aTriangulation;

  friend double myfunc<dim>(const std::vector<double>&, std::vector<double>&, void*);

  aVector m_Psi_Ref; // non ghosted
  aVector m_System_RHS; // non ghosted
  aVector m_Search_Direction; // non ghosted
  aVector m_Workspace_NG;
  dealii::Vector<double> m_error_per_cell;

  std::array<aVector, N> m_Psi; // non ghosted
  std::array<aVector, N> m_Workspace; // ghosted

  aSparseMatrix m_System_Matrix;

  AffineConstraints<double> m_constraints;
  IndexSet m_locally_owned_dofs;
  IndexSet m_locally_relevant_dofs;

  void screening();

  double m_t[N];
  double m_t_guess[N];

  double m_res;
  double m_res_old;
  double m_resp;
  double m_ti;
  double m_final_error;
  double m_N;
  double m_mu = 0;
  double m_dmu = 0.1;
  double m_gs = 1;
  vector<double> m_omega;
  vector<double> m_epsilon;

  int m_rank = -1;

  unsigned m_counter;
  unsigned m_maxiter = 500;
  unsigned m_global_refinement;
  unsigned m_total_no_cells;
  unsigned m_total_no_active_cells;
  unsigned m_NA;
  unsigned m_Ndmu;
  unsigned m_QN1[3];

  ofstream m_computing_timer_log;
  TimerOutput m_computing_timer;
  MyParameterHandler m_ph;
  bool m_root;
  ConditionalOStream pcout;

  aTriangulation m_Triangulation;
  T& m_FE;
  DoFHandler<dim> m_DOF_Handler;

  std::map<std::string, double> m_coeffs;

  MyUtils::ref_pt_list<N> m_ref_pt_list;
};


template <int dim, int N, class T>
cBaseMPI<dim, N, T>::cBaseMPI(const std::string& xmlfilename, T& FE)
  :
  mpi_communicator(MPI_COMM_WORLD),
  m_computing_timer_log("benchmark.txt"),
  m_computing_timer(mpi_communicator, m_computing_timer_log, TimerOutput::summary, TimerOutput:: cpu_and_wall_times),
  m_ph(xmlfilename),
  m_root(Utilities::MPI::this_mpi_process(mpi_communicator) == 0),
  pcout(cout, m_root),
  m_Triangulation(mpi_communicator, typename Triangulation<dim>::MeshSmoothing(Triangulation<dim>::limit_level_difference_at_vertices | Triangulation<dim>::eliminate_refined_inner_islands | Triangulation<dim>::smoothing_on_refinement | Triangulation<dim>::smoothing_on_coarsening)),
  m_FE(FE),
  m_DOF_Handler(m_Triangulation)
{
  // if constexpr ( std::is_same<T,dealii::FE_Q<dim>>::value )
  // {
  //   std::cout << "Hurra!\n";
  //   m_FE
  // }

  try
  {
    m_omega = m_ph.Get_Physics("omega");
    m_gs = m_ph.Get_Physics("gs_1", 0);
    m_QN1[0] = int(m_ph.Get_Physics("QN1", 0));
    m_QN1[1] = int(m_ph.Get_Physics("QN1", 1));
    m_QN1[2] = int(m_ph.Get_Physics("QN1", 2));

    m_global_refinement = unsigned(m_ph.Get_Mesh("global_refinements", 0));

    m_ti = m_ph.Get_Algorithm("ti", 0);
    m_epsilon = m_ph.Get_Algorithm("epsilon");
    m_t[0] = m_ti;
    m_t[1] = m_ti;
    m_t_guess[0] = m_ti;
    m_t_guess[1] = m_ti;

    m_NA = int(m_ph.Get_Algorithm("NA", 0));
    m_Ndmu = m_ph.Get_Algorithm("Ndmu", 0);
    m_dmu = m_ph.Get_Algorithm("dmu", 0);
  }
  catch (const std::string& info)
  {
    std::cerr << info << endl;
    MPI_Abort(mpi_communicator, 0);
  }

  MPI_Comm_rank(mpi_communicator, &m_rank);

  m_counter = 0;
  m_final_error = 0;
}

template <int dim, int N, class T>
double cBaseMPI<dim, N, T>::l2norm_t()
{
  double retval = 0;
  for (int i = 0; i < N; i++)
  {
    retval += m_t[i] * m_t[i];
  }
  return sqrt(retval);
}

template <int dim, int N, class T>
void cBaseMPI<dim, N, T>::screening()
{
  m_ref_pt_list.reset(5, 20);

  for (auto& it : m_ref_pt_list.m_list)
  {
    nlopt::opt opt(nlopt::LD_MMA, N);

    opt.set_xtol_rel(1e-10);
    opt.set_min_objective(myfunc<dim>, this);

    std::vector<double> x(2);
    double minf;

    for (int i = 0; i < N; i++)
    {
      x[i] = it.ti[i];
    }

    int status = opt.optimize(x, minf);

    it.status = status;
    it.failed = (status < 0);

    if (!isfinite(minf))
    {
      continue;
    }
    it.f = minf;
    for (int i = 0; i < N; i++)
    {
      it.t[i] = x[i];
    }
  }
}

template <int dim, int N, class T>
int cBaseMPI<dim, N, T>::find_ortho_min()
{
  compute_contributions();

  TimerOutput::Scope timing_section(m_computing_timer, "");

  if (m_root)
  {
    double l2_norm_t_old = 0;
    double min = std::numeric_limits<double>::max();
    for (int i = 0; i < N; i++)
    {
      l2_norm_t_old += m_t[i] * m_t[i];
    }

    l2_norm_t_old = sqrt(l2_norm_t_old);

    cBaseMPI<dim, N, T>::screening();

    m_ref_pt_list.condense();
    m_ref_pt_list.Dump(std::cout);

    for (auto it : m_ref_pt_list.m_list)
    {
      double tmp1 = fabs(it.l2norm_t() - l2_norm_t_old);
      if (tmp1 < min)
      {
        min = tmp1;
        for (int i = 0; i < N; i++)
        {
          m_t[i] = it.t[i];
        }
      }
    }
  }
  int retval = m_ref_pt_list.m_list.empty();
  MPI_Bcast(m_t, N, MPI_DOUBLE, 0, mpi_communicator);
  MPI_Bcast(&retval, 1, MPI_INT, 0, mpi_communicator);
  
  return retval;
}

template <int dim, int N, class T>
void cBaseMPI<dim, N, T>::dump_info_xml(const string path)
{
  string filename = path + "info.xml";

  pugi::xml_document doc;
  pugi::xml_node parameter_node = doc.append_child("INFO");

  pugi::xml_node node = parameter_node.append_child("MU");
  node.append_child(pugi::node_pcdata).set_value(to_string(m_mu).c_str());

  node = parameter_node.append_child("GS");
  node.append_child(pugi::node_pcdata).set_value(to_string(m_gs).c_str());

  node = parameter_node.append_child("N");
  node.append_child(pugi::node_pcdata).set_value(to_string(m_N).c_str());

  node = parameter_node.append_child("FINAL_ERROR");
  node.append_child(pugi::node_pcdata).set_value(to_string(m_final_error).c_str());

  node = parameter_node.append_child("REVISION");
  node.append_child(pugi::node_pcdata).set_value(STR2(GIT_SHA1));

  doc.save_file(filename.c_str());
}

template <int dim, int N, class T>
void cBaseMPI<dim, N, T>::setup_system()
{
  TimerOutput::Scope timing_section(m_computing_timer, "");

  m_DOF_Handler.distribute_dofs(m_FE);

  m_locally_owned_dofs = m_DOF_Handler.locally_owned_dofs();
  DoFTools::extract_locally_relevant_dofs(m_DOF_Handler, m_locally_relevant_dofs);

  m_Psi_Ref.reinit(m_locally_owned_dofs, m_locally_relevant_dofs, mpi_communicator);
  m_Search_Direction.reinit(m_locally_owned_dofs, mpi_communicator);
  m_System_RHS.reinit(m_locally_owned_dofs, mpi_communicator);
  m_Workspace_NG.reinit(m_locally_owned_dofs, mpi_communicator);
  m_error_per_cell.reinit(m_Triangulation.n_active_cells());

  for (int i = 0; i < N; ++i)
  {
    m_Psi[i].reinit(m_locally_owned_dofs, mpi_communicator);
    m_Workspace[i].reinit(m_locally_owned_dofs, m_locally_relevant_dofs, mpi_communicator);
  }

  m_constraints.clear();
  m_constraints.reinit(m_locally_relevant_dofs);
  DoFTools::make_hanging_node_constraints(m_DOF_Handler, m_constraints);
  VectorTools::interpolate_boundary_values(m_DOF_Handler, 0, ZeroFunction<dim>(), m_constraints);
  m_constraints.close();

  DynamicSparsityPattern csp(m_locally_relevant_dofs);
  DoFTools::make_sparsity_pattern(m_DOF_Handler, csp, m_constraints, false);
  SparsityTools::distribute_sparsity_pattern(csp, m_DOF_Handler.n_locally_owned_dofs_per_processor(), mpi_communicator, m_locally_relevant_dofs);
  m_System_Matrix.reinit(m_locally_owned_dofs, m_locally_owned_dofs, csp, mpi_communicator);

  
}

template <int dim, int N, class T>
void cBaseMPI<dim, N, T>::do_linear_superposition()
{
  m_Psi_Ref = 0;

  for (int i = 0; i < N; i++)
  {
    m_Psi_Ref.add(m_t[i], m_Psi[i]);
  }
  m_constraints.distribute(m_Psi_Ref);
}

template <int dim, int N, class T>
void cBaseMPI<dim, N, T>::update_workspace()
{
  for (int i = 0; i < N; i++)
  {
    m_constraints.distribute(m_Psi[i]);
    m_Workspace[i] = m_Psi[i];
  }
}

template <int dim, int N, class T>
void cBaseMPI<dim, N, T>::save(const std::string& filename)
{
  m_constraints.distribute(m_Psi_Ref);
  parallel::distributed::SolutionTransfer<dim, LA::MPI::Vector> solution_transfer(m_DOF_Handler);
  solution_transfer.prepare_for_serialization(m_Psi_Ref);

  m_Triangulation.save(filename.c_str());
}

template <int dim, int N, class T>
void cBaseMPI<dim, N, T>::output_results(const std::string& path, std::string prefix)
{
  TimerOutput::Scope timing_section(m_computing_timer, "");

  string filename;

  Vector<float> subdomain(m_Triangulation.n_active_cells());
  for (unsigned int i = 0; i < subdomain.size(); ++i)
  {
    subdomain(i) = m_Triangulation.locally_owned_subdomain();
  }

  DataOut<dim> data_out;
  data_out.attach_dof_handler(m_DOF_Handler);
  data_out.add_data_vector(m_Psi_Ref, "Psi_sol");
  data_out.add_data_vector(m_error_per_cell, "error per cell");
  data_out.add_data_vector(subdomain, "subdomain");
  data_out.build_patches();

  filename = path + prefix + "-" + Utilities::int_to_string(m_counter, 5) + ".vtu";
  data_out.write_vtu_in_parallel(filename.c_str(), mpi_communicator);

  
}

template <int dim, int N, class T>
void cBaseMPI<dim, N, T>::output_guess()
{
  TimerOutput::Scope timing_section(m_computing_timer, "");

  this->update_workspace();

  // CPotential Potential_fct ( m_omega, m_QN1[2] );
  // VectorTools::interpolate (this->m_DOF_Handler, Potential_fct, this->m_Workspace_NG );
  // this->m_constraints.distribute(this->m_Workspace_NG);
  // this->m_Psi_Ref=this->m_Workspace_NG;

  DataOut<dim> data_out;
  data_out.attach_dof_handler(this->m_DOF_Handler);
  data_out.add_data_vector(this->m_Workspace[0], "Psi_0");  // todo : loop
  data_out.add_data_vector(this->m_Workspace[1], "Psi_1");
  // data_out.add_data_vector (this->m_Psi_Ref, "m_Potential");
  data_out.build_patches();
  data_out.write_vtu_in_parallel("guess.vtu", mpi_communicator);

  
}
