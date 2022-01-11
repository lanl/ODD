//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   solver/Grey_Matrix.cc
 * \author Mathew Cleveland
 * \date   January 5th 2022
 * \brief  Build Matrix data to be used by the solver
 * \note   Copyright (C) 2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "Grey_Matrix.hh"
#include "Constants.hh"
#include "Opacity_Reader.hh"
#include "ds++/dbc.hh"
#include <cmath>
#include <numeric>
#include <sstream>

namespace odd_solver {

//================================================================================================//
/*!
 * \brief
 *
 * Grey Matrix constructor sets initial constructor data
 *
 *  \param[in] control_data interface control data
 *
 */
//================================================================================================//
Grey_Matrix::Grey_Matrix(const Control_Data &control_data)
    : reflect_bnd(control_data.reflect_bnd), bnd_temp(control_data.bnd_temp) {
  Insist(!control_data.multigroup, "Multigroup currently not supported");
}

//================================================================================================//
/*!
 * \brief
 *
 * Helper function to mass average variables
 *
 *  \param[in] mass of each material
 *  \param[in] variable to be averaged for each material
 *  \return mass averaged variable
 *
 */
//================================================================================================//
double Grey_Matrix::mass_average(const std::vector<double> &mass,
                                 const std::vector<double> &variable) const {
  Require(mass.size() == variable.size());
  const double total_mass = std::accumulate(mass.begin(), mass.end(), 0);
  const double accumulated_value =
      std::inner_product(mass.begin(), mass.end(), variable.begin(), 0.0);
  return total_mass > 0.0 ? accumulated_value / total_mass : 0.0;
}

//================================================================================================//
/*!
 * \brief
 *
 * Initialize solver data populates the fundimental physical data, from the interface, that will be
 * used to generate the solution matrix.
 *
 *  \param[in] mesh orthogonal mesh class
 *  \param[in] mat_data material interface data
 *  \param[in] dt time step size
 *
 */
//================================================================================================//
void Grey_Matrix::initialize_solver_data(const Orthogonal_Mesh &mesh, const Mat_Data &mat_data,
                                         const double dt) {
  Insist(!mesh.domain_decomposed(), "Domain decomposition currently not supported");
  Opacity_Reader opacity_reader(mat_data.ipcress_filename, mat_data.problem_matids);
  const auto ncells = mesh.number_of_local_cells();
  // Allocate solver cell data
  solver_data.diagonal.resize(ncells, 0.0);
  solver_data.source.resize(ncells, 0.0);
  solver_data.cell_density.resize(ncells, 0.0);
  solver_data.cell_cve.resize(ncells, 0.0);
  solver_data.cell_eden0.resize(ncells, 0.0);
  solver_data.cell_temperature0.resize(ncells, 0.0);
  solver_data.cell_eden.resize(ncells, 0.0);
  solver_data.cell_temperature.resize(ncells, 0.0);
  // Allocate solver face/neighbor data
  solver_data.off_diagonal.resize(ncells);
  solver_data.off_diagonal_id.resize(ncells);
  solver_data.face_flux0.resize(ncells);
  solver_data.face_flux.resize(ncells);
  // Resize local matrix data
  sigma_a.resize(ncells, 0.0);
  fleck.resize(ncells, 0.0);
  face_D.resize(ncells);
  // Loop over all cells
  for (size_t cell = 0; cell < ncells; cell++) {
    const auto cell_volume = mesh.cell_volume(cell);
    Check(cell_volume > 0.0);
    // Populate the homogenized cell material data
    const auto nmats = mat_data.number_of_cell_mats[cell];
    // Calculate the cell mat data
    std::vector<double> cell_mat_sigma_tr(nmats, 0.0);
    std::vector<double> cell_mat_sigma_a(nmats, 0.0);
    std::vector<double> cell_mat_mass(nmats, 0.0);
    std::vector<double> cell_mat_cv(nmats, 0.0);
    std::vector<double> cell_mat_T4(nmats, 0.0);
    std::vector<double> cell_mat_emission(nmats, 0.0);
    for (size_t mat = 0; mat < nmats; mat++) {
      size_t matid = mat_data.cell_mats[cell][mat];
      cell_mat_sigma_a[matid] =
          opacity_reader.mat_planck_abs_models[mat]->getOpacity(
              mat_data.cell_mat_temperature[cell][mat], mat_data.cell_mat_density[cell][mat]) *
          mat_data.cell_mat_density[cell][mat];
      cell_mat_sigma_tr[mat] =
          opacity_reader.mat_rosseland_total_models[matid]->getOpacity(
              mat_data.cell_mat_temperature[cell][mat], mat_data.cell_mat_density[cell][mat]) *
          mat_data.cell_mat_density[cell][mat];
      cell_mat_mass[mat] = mat_data.cell_mat_density[cell][mat] *
                           mat_data.cell_mat_vol_frac[cell][mat] * cell_volume;
      cell_mat_cv[mat] =
          mat_data.cell_mat_density[cell][mat] * mat_data.cell_mat_specific_heat[cell][mat];
      cell_mat_T4[mat] =
          mat_data.cell_mat_temperature[cell][mat] * mat_data.cell_mat_temperature[cell][mat] *
          mat_data.cell_mat_temperature[cell][mat] * mat_data.cell_mat_temperature[cell][mat];
      cell_mat_emission[mat] = cell_mat_T4[mat] * cell_mat_sigma_a[mat];
    }
    const double cell_mass = std::accumulate(cell_mat_mass.begin(), cell_mat_mass.end(), 0.0);
    Check(cell_mass > 0.0);
    const double cell_density = cell_mass / cell_volume;
    Check(cell_density > 0.0);
    const double cell_cve = std::inner_product(cell_mat_cv.begin(), cell_mat_cv.end(),
                                               mat_data.cell_mat_vol_frac[cell].begin(), 0.0);
    Check(cell_cve > 0.0);

    sigma_a[cell] = std::inner_product(cell_mat_sigma_a.begin(), cell_mat_sigma_a.end(),
                                       mat_data.cell_mat_vol_frac[cell].begin(), 0.0);

    // Calc mass averaged T4 (just in case the opacity is zero)
    const double mass_avg_T4 = mass_average(cell_mat_mass, cell_mat_T4);
    const double cell_T =
        sigma_a[cell] > 0.0
            ? std::pow(std::inner_product(cell_mat_emission.begin(), cell_mat_emission.end(),
                                          mat_data.cell_mat_vol_frac[cell].begin(), 0.0) /
                           sigma_a[cell],
                       0.25)
            : std::pow(mass_avg_T4, 0.25);

    solver_data.cell_density[cell] = cell_density;
    solver_data.cell_cve[cell] = cell_cve;
    Check(cell_cve > 0.0);
    fleck[cell] = 1.0 / (1.0 + constants::a * constants::c * sigma_a[cell] * 4.0 * cell_T * cell_T *
                                   cell_T * dt / cell_cve);
    // Set initial conditions
    solver_data.cell_temperature0[cell] = cell_T;
    solver_data.cell_temperature[cell] = cell_T;
    solver_data.cell_eden0[cell] = mat_data.cell_rad_eden[cell];
    solver_data.cell_eden[cell] = mat_data.cell_rad_eden[cell];
    // Allocate face data
    const auto nfaces = mesh.number_of_faces(cell);
    solver_data.off_diagonal[cell].resize(nfaces, 0.0);
    solver_data.face_flux0[cell].resize(nfaces, 0.0);
    solver_data.face_flux[cell].resize(nfaces, 0.0);
    solver_data.off_diagonal_id[cell].resize(nfaces, ncells);
    face_D[cell].resize(nfaces, 0.0);
    // Loop over faces
    for (size_t face = 0; face < nfaces; face++) {
      const auto ftype = mesh.face_type(cell, face);
      Check(ftype != FACE_TYPE::GHOST_FACE);
      if (ftype == FACE_TYPE::INTERNAL_FACE) {
        solver_data.off_diagonal_id[cell][face] = mesh.next_cell(cell, face);
      } else {
        Check(ftype == FACE_TYPE::BOUNDARY_FACE);
        Check(solver_data.off_diagonal_id[cell][face] == ncells);
      }
    }
  }
  // calculate the face diffusion coefficients
  for (size_t cell = 0; cell < ncells; cell++) {
    const auto cell_T4 = solver_data.cell_temperature[cell] * solver_data.cell_temperature[cell] *
                         solver_data.cell_temperature[cell] * solver_data.cell_temperature[cell];

    const auto nfaces = mesh.number_of_faces(cell);
    for (size_t face = 0; face < nfaces; face++) {
      const auto ftype = mesh.face_type(cell, face);
      Check(ftype != FACE_TYPE::GHOST_FACE);
      // calculate the average face temperature
      const auto nmats = mat_data.number_of_cell_mats[cell];
      if (ftype == FACE_TYPE::INTERNAL_FACE) {
        const auto next_cell = solver_data.off_diagonal_id[cell][face];
        const auto next_cell_T4 =
            solver_data.cell_temperature[next_cell] * solver_data.cell_temperature[next_cell] *
            solver_data.cell_temperature[next_cell] * solver_data.cell_temperature[next_cell];
        // Calculate the face averaged temperature
        const auto face_T = std::pow(0.5 * (next_cell_T4 + cell_T4), 0.25);
        // Calculate the face opacity
        std::vector<double> cell_mat_sigma_tr(nmats, 0.0);
        for (size_t mat = 0; mat < nmats; mat++) {
          size_t matid = mat_data.cell_mats[cell][mat];
          cell_mat_sigma_tr[mat] = opacity_reader.mat_rosseland_total_models[matid]->getOpacity(
                                       face_T, mat_data.cell_mat_density[cell][mat]) *
                                   mat_data.cell_mat_density[cell][mat];
        }
        const double sigma_tr =
            std::inner_product(cell_mat_sigma_tr.begin(), cell_mat_sigma_tr.end(),
                               mat_data.cell_mat_vol_frac[cell].begin(), 0.0);
        face_D[cell][face] = 1.0 / (3.0 * sigma_tr);
      } else {
        Check(ftype == FACE_TYPE::BOUNDARY_FACE);
        Check(solver_data.off_diagonal_id[cell][face] == ncells);
        Insist(reflect_bnd[face], "Vacuum/Source boundaries are not yet supported");
        // Calculate the face averaged temperature
        const auto face_T = solver_data.cell_temperature[cell];
        // Calculate the face opacity
        std::vector<double> cell_mat_sigma_tr(nmats, 0.0);
        for (size_t mat = 0; mat < nmats; mat++) {
          size_t matid = mat_data.cell_mats[cell][mat];
          cell_mat_sigma_tr[mat] = opacity_reader.mat_rosseland_total_models[matid]->getOpacity(
                                       face_T, mat_data.cell_mat_density[cell][mat]) *
                                   mat_data.cell_mat_density[cell][mat];
        }
        const double sigma_tr =
            std::inner_product(cell_mat_sigma_tr.begin(), cell_mat_sigma_tr.end(),
                               mat_data.cell_mat_vol_frac[cell].begin(), 0.0);
        face_D[cell][face] = 1.0 / (3.0 * sigma_tr);
      }
    }
  }
}

//================================================================================================//
/*!
 * \brief
 *
 * Build the grey solution matrix using the initialized matrix data and the mesh.
 *
 *  \param[in] mesh orthogonal mesh class
 *  \param[in] dt time step size
 *
 */
//================================================================================================//
void Grey_Matrix::build_matrix(const Orthogonal_Mesh &mesh, const double dt) {
  Insist(!mesh.domain_decomposed(), "Domain decomposition currently not supported");
  const auto ncells = mesh.number_of_local_cells();
  // Loop over all cells
  for (size_t cell = 0; cell < ncells; cell++) {
    // Populate the homogenized cell material data
    const auto cell_volume = mesh.cell_volume(cell);
    // Apply collision terms to the diagonal (census + absorption)
    solver_data.diagonal[cell] = 1.0 + constants::c * dt * fleck[cell] * sigma_a[cell];
    // Apply the emission source
    solver_data.source[cell] =
        solver_data.cell_eden0[cell] +
        constants::a * constants::c * 4.0 * fleck[cell] * solver_data.cell_temperature[cell] *
            solver_data.cell_temperature[cell] * solver_data.cell_temperature[cell] *
            solver_data.cell_temperature[cell] * dt;

    const auto nfaces = mesh.number_of_faces(cell);
    for (size_t face = 0; face < nfaces; face++) {
      const auto ftype = mesh.face_type(cell, face);
      Check(ftype != FACE_TYPE::GHOST_FACE);
      // calculate the fringe values
      if (ftype == FACE_TYPE::INTERNAL_FACE) {
        const auto face_area = mesh.face_area(cell, face);
        const auto next_cell = solver_data.off_diagonal_id[cell][face];
        const auto next_face = face % 2 == 0 ? face + 1 : face - 1;
        const auto D = face_D[cell][face];
        const auto next_D = face_D[next_cell][face];
        const auto half_width = mesh.distance_center_to_face(cell, face);
        const auto next_half_width = mesh.distance_center_to_face(next_cell, next_face);
        const double fring = -constants::c * dt * face_area * (D * next_D) /
                             (cell_volume * (half_width * next_D + next_half_width * D));
        solver_data.off_diagonal[cell][face] = fring;
        solver_data.diagonal[cell] -= fring;
      } else {
        Check(ftype == FACE_TYPE::BOUNDARY_FACE);
        Check(solver_data.off_diagonal_id[cell][face] == ncells);
        Insist(reflect_bnd[face], "Vacuum/Source boundaries are not yet supported");
        Insist(bnd_temp[0] < 1.0e6, "Junk check for bnd_temp to be used later");
      }
    }
  }
}

void Grey_Matrix::gs_solver(const double eps, const size_t max_iter) {
  Require(max_iter > 0);
  double max_error = 1.0;
  size_t count = 0;
  std::stringstream diagnostics;
  while (max_error > eps && count < max_iter) {
    max_error = 0.0;
    for (size_t i = 0; i < solver_data.diagonal.size(); i++) {
      const double diag = solver_data.diagonal[i];
      double b = solver_data.source[i];
      for (size_t f = 0; f < solver_data.off_diagonal_id[i].size(); f++) {
        const size_t next_cell = solver_data.off_diagonal_id[i][f];
        b -= solver_data.off_diagonal[i][f] * solver_data.cell_eden[next_cell];
      }
      Check(diag > 0.0);
      const double last_eden = solver_data.cell_eden[i];
      solver_data.cell_eden[i] = b / diag;
      const double eden_diff =
          solver_data.cell_eden[i] > 0.0
              ? std::abs(solver_data.cell_eden[i] - last_eden) / solver_data.cell_eden[i]
              : std::abs(solver_data.cell_eden[i] - last_eden);
      max_error = std::max(max_error, eden_diff);
    }
    diagnostics << "  Iteration = " << count << " error = " << max_error << std::endl;
    count++;
  }
  if (max_error > eps) {
    std::cout << "WARNING: Did not converge cell_eden max_error = " << max_error << std::endl;
    std::cout << diagnostics.str();
  } else {
    std::cout << "Converged eden -> Iteration = " << count << " error = " << max_error << std::endl;
  }
}

} // namespace odd_solver

//------------------------------------------------------------------------------------------------//
// end of Grey_Matrix.cc
//------------------------------------------------------------------------------------------------//
