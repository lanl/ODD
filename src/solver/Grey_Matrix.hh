//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   solver/Grey_Matrix.hh
 * \author Mathew Cleveland
 * \brief  Define class Grey_Matrix
 * \note   Copyright (C) 2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#ifndef odd_solver_Grey_Matrix_hh
#define odd_solver_Grey_Matrix_hh

#include "Interface_Data.hh"
#include "Orthogonal_Mesh.hh"

namespace odd_solver {

//================================================================================================//
/*!
 * \class Grey_Matrix
 * \brief
 *
 * This builder function constructs the base matrix data that will be used by the diffusion solver
 *
 */
//================================================================================================//

class Grey_Matrix {
public:
  //! Default constructors.
  Grey_Matrix(const Control_Data &control_data);

  //! Initialize solver data
  void initialize_solver_data(const Orthogonal_Mesh &mesh, const Mat_Data &mat_data,
                              const double dt);

  //! Build the solution matrix
  void build_matrix(const Orthogonal_Mesh &mesh, const double dt);

  //! solve using Gauss Siedel
  void gs_solver(const double eps, const size_t max_iter);

private:
  // boundary conditions
  const std::array<bool, 6> reflect_bnd;
  const std::array<double, 6> bnd_temp;

  // Local matrix data
  std::vector<double> fleck;
  std::vector<double> sigma_a;
  std::vector<std::vector<double>> face_D;

  //! Helper function for mass averaging
  double mass_average(const std::vector<double> &mass, const std::vector<double> &variable) const;

public:
  //! Fundamental matrix data and solution vectors
  Solver_Data solver_data;
};

} // end namespace odd_solver

#endif // odd_solver_Grey_Matrix_hh

//------------------------------------------------------------------------------------------------//
// end of solver/Grey_Matrix.hh
//------------------------------------------------------------------------------------------------//
