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

template <int dim>
class MySolver;


template <int dim>
double myfunc2(const std::vector<double>& t, std::vector<double>& grad, void* data)
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
              //+ coeffs.at("t1^2") * t1q
              + 4 * coeffs.at("t0^4") * t0q * t[0]
              + coeffs.at("t0^1_t1^3") * t[1] * t1q
              + 2 * coeffs.at("t0^2_t1^2") * t[0] * t1q
              + 3 * coeffs.at("t0^3_t1^1")  * t[1] * t0q;
    //+ coeffs.at("t1^4") * t1q * t1q;

    grad[1] =  //coeffs.at("t0^2") * t0q
      + coeffs.at("t0_t1") * t[0]
      + 2 * coeffs.at("t1^2") * t[1]
      //+ coeffs.at("t0^4") * t0q * t0q
      + 3 * coeffs.at("t0^1_t1^3") * t[0] * t1q
      + 2 * coeffs.at("t0^2_t1^2") * t0q * t[1]
      + coeffs.at("t0^3_t1^1") * t[0] * t0q
      + 4 * coeffs.at("t1^4") * t1q * t[1];
  }

  return retval;
}

enum Status { SUCCESS, FAILED, ZERO_SOL, SLOW_CONV, MAXITER, SINGULAR, CONTINUE };

template <int dim, int N>
class cBase
{
public:
  explicit cBase(const std::string&);
  virtual ~cBase()
  {
    m_DOF_Handler.clear();
  };

  int find_ortho_min();

  void dump_info_xml(const string = "");

  double l2norm_t();

  virtual void compute_contributions() = 0;

  void setup_system();

  void do_linear_superposition();

  void save(const std::string&);

  void save_one(const std::string&);

  void output_guess();

  void output_results(const std::string&, std::string = "step");

protected:

  friend double myfunc2<dim>(const std::vector<double>&, std::vector<double>&, void*);

  typedef typename dealii::Vector<double> aVector;
  typedef typename dealii::SparseMatrix<double> aSparseMatrix;
  typedef typename dealii::Triangulation<dim> aTriangulation;

  aVector m_Psi_Ref;
  aVector m_System_RHS;
  aVector m_Search_Direction;
  dealii::Vector<double> m_error_per_cell;

  std::array<aVector, N> m_Psi;
  std::array<aVector, N> m_Workspace;

  aSparseMatrix m_System_Matrix;

  AffineConstraints<double> m_constraints;

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

  MyParameterHandler m_ph;

  aTriangulation m_Triangulation;
  FE_Q<dim> m_FE;
  DoFHandler<dim> m_DOF_Handler;

  std::map<std::string, double> m_coeffs;

  MyUtils::ref_pt_list<N> m_ref_pt_list;
};


template <int dim, int N>
cBase<dim, N>::cBase(const std::string& xmlfilename)
  :
  m_ph(xmlfilename),
  m_Triangulation(typename Triangulation<dim>::MeshSmoothing(Triangulation<dim>::limit_level_difference_at_vertices | Triangulation<dim>::eliminate_refined_inner_islands | Triangulation<dim>::smoothing_on_refinement | Triangulation<dim>::smoothing_on_coarsening)),
  m_FE(gl_degree_fe),
  m_DOF_Handler(m_Triangulation)
{
  try
  {
    // m_omega = m_ph.Get_Physics("omega");
    // m_gs = m_ph.Get_Physics("gs_1", 0);
    // m_QN1[0] = int(m_ph.Get_Physics("QN1", 0));
    // m_QN1[1] = int(m_ph.Get_Physics("QN1", 1));
    // m_QN1[2] = int(m_ph.Get_Physics("QN1", 2));

    // m_global_refinement = unsigned(m_ph.Get_Mesh("global_refinements", 0));

    // m_ti = m_ph.Get_Algorithm("ti", 0);
    // m_epsilon = m_ph.Get_Algorithm("epsilon");
    // m_t[0] = m_ti;
    // m_t[1] = m_ti;
    // m_t_guess[0] = m_ti;
    // m_t_guess[1] = m_ti;

    // m_NA = int(m_ph.Get_Algorithm("NA", 0));
    // m_Ndmu = m_ph.Get_Algorithm("Ndmu", 0);
    // m_dmu = m_ph.Get_Algorithm("dmu", 0);
  }
  catch (const std::string& info)
  {
    std::cerr << info << endl;
    //MPI_Abort(mpi_communicator, 0);
  }

  m_counter = 0;
  m_final_error = 0;
}

template <int dim, int N>
double cBase<dim, N>::l2norm_t()
{
  double retval = 0;
  for (int i = 0; i < N; i++)
  {
    retval += m_t[i] * m_t[i];
  }
  return sqrt(retval);
}

template <int dim, int N>
void cBase<dim, N>::screening()
{
  m_ref_pt_list.reset(5, 20);

  for (auto& it : m_ref_pt_list.m_list)
  {
    nlopt::opt opt(nlopt::LD_MMA, N);

    opt.set_xtol_rel(1e-10);
    opt.set_min_objective(myfunc2<dim>, this);

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

template <int dim, int N>
int cBase<dim, N>::find_ortho_min()
{
  compute_contributions();

  double l2_norm_t_old = 0;
  double min = std::numeric_limits<double>::max();
  for (int i = 0; i < N; i++)
  {
    l2_norm_t_old += m_t[i] * m_t[i];
  }

  l2_norm_t_old = sqrt(l2_norm_t_old);

  cBase<dim, N>::screening();

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

  int retval = m_ref_pt_list.m_list.empty();

  return retval;
}

template <int dim, int N>
void cBase<dim, N>::dump_info_xml(const string path)
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

template <int dim, int N>
void cBase<dim, N>::setup_system()
{
  m_DOF_Handler.distribute_dofs(m_FE);

  const auto ndofs = m_DOF_Handler.n_dofs();

  m_Psi_Ref.reinit(ndofs);
  m_Search_Direction.reinit(ndofs);
  m_System_RHS.reinit(ndofs);
  m_error_per_cell.reinit(m_Triangulation.n_active_cells());

  for (int i = 0; i < N; ++i)
  {
    m_Psi[i].reinit(ndofs);
    m_Workspace[i].reinit(ndofs);
  }

  // m_constraints.clear();
  // m_constraints.reinit(ndofs);
  // DoFTools::make_hanging_node_constraints(m_DOF_Handler, m_constraints);
  // VectorTools::interpolate_boundary_values(m_DOF_Handler, 0, ZeroFunction<dim>(), m_constraints);
  // m_constraints.close();

  DynamicSparsityPattern dsp(ndofs);
  DoFTools::make_sparsity_pattern(m_DOF_Handler, dsp);
  SparsityPattern sp;
  sp.copy_from(dsp);
  m_System_Matrix.reinit(sp);
}

template <int dim, int N>
void cBase<dim, N>::do_linear_superposition()
{
  m_Psi_Ref = 0;

  for (int i = 0; i < N; i++)
  {
    m_Psi_Ref.add(m_t[i], m_Psi[i]);
  }
  m_constraints.distribute(m_Psi_Ref);
}


template <int dim, int N>
void cBase<dim, N>::save(const std::string& filename)
{
  ofstream out(filename);
  this->m_Psi_Ref.block_write(out);
}


template <int dim, int N>
void cBase<dim, N>::save_one(const std::string& filename)
{
  ofstream out(filename);

  m_Workspace[0] = this->m_Psi_Ref;
  m_Workspace[0] *= 1.0 / sqrt(MyRealTools::Particle_Number<dim>(this->m_DOF_Handler, m_FE, this->m_Psi_Ref));
  m_Workspace[0].block_write(out);
}

template <int dim, int N>
void cBase<dim, N>::output_results(const std::string& path, std::string prefix)
{
}


template <int dim, int N>
void cBase<dim, N>::output_guess()
{
}
