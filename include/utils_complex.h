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
  class ComputeIntensity : public DataPostprocessorScalar<dim>
  {
  public:
    ComputeIntensity ( string="Intensity" );
    virtual void compute_derived_quantities_vector (const vector<Vector<double>> &uh,
                                                    const vector<vector<Tensor<1,dim>>> &duh,
                                                    const vector<vector<Tensor<2,dim>>> &dduh,
                                                    const vector<Point<dim>> &normals,
                                                    const vector<Point<dim> > &evaluation_points,
                                                    vector<Vector<double>> &computed_quantities) const;
  };

  template <int dim>
  ComputeIntensity<dim>::ComputeIntensity ( string name ) : DataPostprocessorScalar<dim> ( name, update_values)
  {}

  template <int dim> 
  void ComputeIntensity<dim>::compute_derived_quantities_vector ( const vector<Vector<double>> &uh,
                                                                  const vector<vector<Tensor<1,dim>>> &/*duh*/,
                                                                  const vector<vector<Tensor<2,dim>>> &/*dduh*/,
                                                                  const vector<Point<dim>> &/*normals*/,
                                                                  const vector<Point<dim>> &/*evaluation_points*/,
                                                                  vector<Vector<double>> &computed_quantities) const
  {
    Assert(computed_quantities.size() == uh.size(), ExcDimensionMismatch (computed_quantities.size(), uh.size()));

    for (unsigned int i=0; i<computed_quantities.size(); i++)
    {
      Assert(computed_quantities[i].size() == 1, ExcDimensionMismatch (computed_quantities[i].size(), 1));
      Assert(uh[i].size() == 2, ExcDimensionMismatch (uh[i].size(), 2));

      computed_quantities[i](0) = sqrt(uh[i](0)*uh[i](0) + uh[i](1)*uh[i](1));
    }
  }

  template <int dim>
  class ComputePhase : public DataPostprocessorScalar<dim>
  {
  public:
    ComputePhase ();
    virtual void compute_derived_quantities_vector (const vector<Vector<double>> &uh,
                                                    const vector<vector<Tensor<1,dim>>> &duh,
                                                    const vector<vector<Tensor<2,dim>>> &dduh,
                                                    const vector<Point<dim>> &normals,
                                                    const vector<Point<dim> > &evaluation_points,
                                                    vector<Vector<double>> &computed_quantities) const;
  };

  template <int dim>
  ComputePhase<dim>::ComputePhase () : DataPostprocessorScalar<dim> ("Phase", update_values)
  {}

  template <int dim> 
  void ComputePhase<dim>::compute_derived_quantities_vector ( const vector<Vector<double>> &uh,
                                                              const vector<vector<Tensor<1,dim>>> &/*duh*/,
                                                              const vector<vector<Tensor<2,dim>>> &/*dduh*/,
                                                              const vector<Point<dim>> &/*normals*/,
                                                              const vector<Point<dim>> &/*evaluation_points*/,
                                                              vector<Vector<double>> &computed_quantities) const
  {
    Assert(computed_quantities.size() == uh.size(), ExcDimensionMismatch (computed_quantities.size(), uh.size()));

    for (unsigned int i=0; i<computed_quantities.size(); i++)
    {
      Assert(computed_quantities[i].size() == 1, ExcDimensionMismatch (computed_quantities[i].size(), 1));
      Assert(uh[i].size() == 2, ExcDimensionMismatch (uh[i].size(), 2));

      computed_quantities[i](0) = atan2(uh[i](1),uh[i](0));
    }
  }

