//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   solver/Mg_P1_Matrix.hh
 * \author Mathew Cleveland
 * \brief  Define class Mg_P1_Matrix
 * \note   Copyright (C) 2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#ifndef odd_solver_Mg_P1_Matrix_hh
#define odd_solver_Mg_P1_Matrix_hh

#include "Ghost_Comm.hh"
#include "Interface_Data.hh"
#include "Orthogonal_Mesh.hh"
#include <memory>

namespace odd_solver {

//================================================================================================//
/*!
 * \class Mg_P1_Matrix
 * \brief
 *
 * This builder function constructs the base matrix data that will be used by the diffusion solver
 *
 */
//================================================================================================//

class Mg_P1_Matrix {
public:
  //! Default constructors.
  Mg_P1_Matrix(const Control_Data &control_data);

  //! Initialize solver data
  void initialize_solver_data(const Orthogonal_Mesh &mesh, const Mat_Data &mat_data,
                              const double dt);

  //! Build the solution matrix
  void build_matrix(const Orthogonal_Mesh &mesh, const double dt);

  //! solve using Gauss Siedel
  void gs_solver(const double eps, const size_t max_iter, const bool print = false);

  //! Calculate the output data
  void calculate_output_data(const Orthogonal_Mesh &mesh, const Mat_Data &mat_data, const double dt,
                             Output_Data &output_data);

private:
  // boundary conditions
  const std::array<bool, 6> reflect_bnd;
  const std::array<double, 6> bnd_temp;
  const size_t ngroups;
  const std::vector<double> group_bounds;

  // Local matrix data
  double current_dt;
  std::vector<double> ext_imp_source;
  std::vector<double> ext_exp_source;
  std::vector<double> rad_source;
  std::vector<double> fleck;
  std::vector<double> cell_epsilon;
  std::vector<double> cell_correction_source;
  std::vector<std::vector<double>> cell_planck_weights;
  std::vector<double> volume;
  std::vector<std::vector<double>> sigma_a;
  std::vector<double> sigma_planck;
  std::vector<std::vector<std::vector<double>>> face_D;
  std::vector<std::vector<std::vector<double>>> face_sigma_tr;
  // Local Ghost data
  std::vector<double> ghost_face_D;

  // Ghost Data Communicator
  std::unique_ptr<Ghost_Comm> gcomm;

  //! Helper function for mass averaging
  double mass_average(const std::vector<double> &mass, const std::vector<double> &variable) const;

public:
  //! Fundamental matrix data and solution vectors
  MG_Solver_Data solver_data;
};

} // end namespace odd_solver

#endif // odd_solver_Mg_P1_Matrix_hh

//------------------------------------------------------------------------------------------------//
// end of solver/Mg_P1_Matrix.hh
//------------------------------------------------------------------------------------------------//
