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
#include "Correction.hh"
#include "Opacity_Reader.hh"
#include "c4/global.hh"
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
Grey_Matrix::Grey_Matrix(const Control_Data &control_data, const bool flux_limiter_)
    : flux_limiter(flux_limiter_), reflect_bnd(control_data.reflect_bnd),
      bnd_temp(control_data.bnd_temp) {
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
  Opacity_Reader opacity_reader(mat_data.ipcress_filename, mat_data.problem_matids);
  const auto ncells = mesh.number_of_local_cells();
  current_dt = dt;
  // declare some dd put data structures
  std::map<size_t, std::vector<double>> put_ghost_cell_temp;
  std::map<size_t, std::vector<double>> put_ghost_cell_eden;
  std::map<size_t, std::vector<double>> put_dist_center_to_face;
  if (mesh.domain_decomposed()) {
    gcomm = std::make_unique<Ghost_Comm>(Ghost_Comm(mesh));
    const auto nghost = gcomm->local_ghost_buffer_size;
    // Initialize dd receive (ghost) data structures
    solver_data.ghost_cell_temperature.resize(nghost, 0.0);
    solver_data.ghost_cell_eden.resize(nghost, 0.0);
    solver_data.ghost_dist_center_to_face.resize(nghost, 0.0);
    // Initialize dd put data structures
    for (auto &map : gcomm->put_buffer_size) {
      put_ghost_cell_temp[map.first] = std::vector<double>(map.second, 0.0);
      put_ghost_cell_eden[map.first] = std::vector<double>(map.second, 0.0);
      put_dist_center_to_face[map.first] = std::vector<double>(map.second, 0.0);
    }
  }
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
  solver_data.flux_source.resize(ncells);
  solver_data.face_type.resize(ncells);
  solver_data.off_diagonal_id.resize(ncells);
  solver_data.face_flux0.resize(ncells);
  // Resize local matrix data
  sigma_a.resize(ncells, 0.0);
  ext_imp_source.resize(ncells, 0.0);
  ext_exp_source.resize(ncells, 0.0);
  rad_source.resize(ncells, 0.0);
  fleck.resize(ncells, 0.0);
  cell_epsilon.resize(ncells, 0.0);
  cell_correction_source.resize(ncells, 0.0);
  volume.resize(ncells, 0.0);
  face_sigma_tr.resize(ncells);
  // Loop over all cells
  for (size_t cell = 0; cell < ncells; cell++) {
    const auto cell_volume = mesh.cell_volume(cell);
    volume[cell] = cell_volume;
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
      cell_mat_sigma_a[mat] =
          opacity_reader.mat_planck_abs_models[matid]->getOpacity(
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
    const double src_val = std::inner_product(mat_data.cell_mat_electron_source[cell].begin(),
                                              mat_data.cell_mat_electron_source[cell].end(),
                                              mat_data.cell_mat_vol_frac[cell].begin(), 0.0);

    ext_imp_source[cell] = src_val * (1.0 - fleck[cell]);
    ext_exp_source[cell] = src_val * fleck[cell];
    rad_source[cell] = mat_data.cell_rad_source[cell];

    // Set initial conditions
    solver_data.cell_temperature0[cell] = cell_T;
    solver_data.cell_temperature[cell] = cell_T;
    solver_data.cell_eden0[cell] = mat_data.cell_rad_eden[cell];
    solver_data.cell_eden[cell] = mat_data.cell_rad_eden[cell];
    // Allocate face data
    const auto nfaces = mesh.number_of_faces(cell);
    solver_data.off_diagonal[cell].resize(nfaces, 0.0);
    solver_data.flux_source[cell].resize(nfaces, 0.0);
    solver_data.face_type[cell].resize(nfaces, 0);
    solver_data.face_flux0[cell].resize(nfaces, 0.0);
    solver_data.off_diagonal_id[cell].resize(nfaces, ncells);
    face_sigma_tr[cell].resize(nfaces, 0.0);
    // Loop over faces
    for (size_t face = 0; face < nfaces; face++) {
      const auto ftype = mesh.face_type(cell, face);
      if (ftype == FACE_TYPE::INTERNAL_FACE || ftype == FACE_TYPE::GHOST_FACE) {
        solver_data.off_diagonal_id[cell][face] = mesh.next_cell(cell, face);
      } else {
        Check(ftype == FACE_TYPE::BOUNDARY_FACE);
        Check(solver_data.off_diagonal_id[cell][face] == ncells);
      }
    }
  }
  // collect ghost cell data
  if (gcomm) {
    for (auto &map : gcomm->put_map) {
      auto cell = map.first;
      for (auto &put_face : map.second) {
        auto rank = put_face.second.first;
        auto buffer_index = put_face.second.second;
        put_ghost_cell_temp[rank][buffer_index] = solver_data.cell_temperature[cell];
        put_ghost_cell_eden[rank][buffer_index] = solver_data.cell_eden[cell];
      }
    }
    gcomm->exchange_ghost_data(put_ghost_cell_temp, solver_data.ghost_cell_temperature);
    gcomm->exchange_ghost_data(put_ghost_cell_eden, solver_data.ghost_cell_eden);
  }
  // calculate the face diffusion coefficients
  for (size_t cell = 0; cell < ncells; cell++) {
    const auto cell_T4 = solver_data.cell_temperature[cell] * solver_data.cell_temperature[cell] *
                         solver_data.cell_temperature[cell] * solver_data.cell_temperature[cell];

    const auto nfaces = mesh.number_of_faces(cell);
    for (size_t face = 0; face < nfaces; face++) {
      const auto ftype = mesh.face_type(cell, face);
      // calculate the average face temperature
      const auto nmats = mat_data.number_of_cell_mats[cell];
      if (ftype == FACE_TYPE::GHOST_FACE) {
        const auto next_cell = solver_data.off_diagonal_id[cell][face];
        const auto next_face = face % 2 == 0 ? face + 1 : face - 1;
        // pull the ghost cell temperature
        const auto ghost_T =
            solver_data.ghost_cell_temperature[gcomm->ghost_map[next_cell][next_face]];
        const auto next_cell_T4 = ghost_T * ghost_T * ghost_T * ghost_T;
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
        face_sigma_tr[cell][face] = sigma_tr;
        // populate the put buffer
        const auto rank = gcomm->put_map[cell][face].first;
        const auto put_index = gcomm->put_map[cell][face].second;
        put_dist_center_to_face[rank][put_index] = mesh.distance_center_to_face(cell, face);

      } else if (ftype == FACE_TYPE::INTERNAL_FACE) {
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
        face_sigma_tr[cell][face] = sigma_tr;
      } else {
        Check(ftype == FACE_TYPE::BOUNDARY_FACE);
        Check(solver_data.off_diagonal_id[cell][face] == ncells);

        // Calculate the face averaged temperature
        auto face_T = solver_data.cell_temperature[cell];

        if (!reflect_bnd[face])
          face_T =
              std::pow(0.5 * (bnd_temp[face] * bnd_temp[face] * bnd_temp[face] * bnd_temp[face] +
                              face_T * face_T * face_T * face_T),
                       0.25);
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
        face_sigma_tr[cell][face] = sigma_tr;
      }
    }
  }
  if (gcomm) {
    // populate distance to center face data
    gcomm->exchange_ghost_data(put_dist_center_to_face, solver_data.ghost_dist_center_to_face);
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
  const auto ncells = mesh.number_of_local_cells();

  std::map<size_t, std::vector<double>> put_face_D;
  if (gcomm) {
    const auto nghost = gcomm->local_ghost_buffer_size;
    // Initialize dd receive (ghost) data structures
    ghost_face_D.resize(nghost, 0.0);
    // Initialize dd put data structures
    for (auto &map : gcomm->put_buffer_size) {
      put_face_D[map.first] = std::vector<double>(map.second, 0.0);
    }
  }
  // Build face diffusion coefficients
  face_D.resize(ncells);
  for (size_t cell = 0; cell < ncells; cell++) {
    const auto cell_eden = solver_data.cell_eden[cell];

    const auto nfaces = mesh.number_of_faces(cell);
    face_D[cell].resize(nfaces, 0.0);
    for (size_t face = 0; face < nfaces; face++) {
      const auto half_width = mesh.distance_center_to_face(cell, face);
      const auto ftype = mesh.face_type(cell, face);
      if (ftype == FACE_TYPE::GHOST_FACE) {
        const auto next_cell = solver_data.off_diagonal_id[cell][face];
        const auto next_face = face % 2 == 0 ? face + 1 : face - 1;
        const auto ghost_eden = solver_data.ghost_cell_eden[gcomm->ghost_map[next_cell][next_face]];
        const double deden = std::abs(cell_eden - ghost_eden) / half_width;
        // J.E. Morel,
        // Diffusion-limit asymptotics of the transport equation, the P1/3 equations, and two
        // flux-limited diffusion theories,
        // Journal of Quantitative Spectroscopy and Radiative Transfer,
        // Volume 65, Issue 5,
        // 2000,
        // Pages 769-778,
        const auto sigma_tr = face_sigma_tr[cell][face];
        face_D[cell][face] = flux_limiter
                                 ? 1.0 / std::sqrt(9.0 * sigma_tr * sigma_tr +
                                                   (deden * deden) / (cell_eden * cell_eden))
                                 : 1.0 / (3.0 * sigma_tr);
        // populate the put buffer
        const auto rank = gcomm->put_map[cell][face].first;
        const auto put_index = gcomm->put_map[cell][face].second;
        put_face_D[rank][put_index] = face_D[cell][face];

      } else if (ftype == FACE_TYPE::INTERNAL_FACE) {
        const auto next_cell = solver_data.off_diagonal_id[cell][face];
        const auto next_eden = solver_data.cell_eden[next_cell];
        const auto max_eden = std::max(cell_eden, next_eden);
        const double deden = std::abs(cell_eden - next_eden) / half_width;
        // J.E. Morel,
        // Diffusion-limit asymptotics of the transport equation, the P1/3 equations, and two
        // flux-limited diffusion theories,
        // Journal of Quantitative Spectroscopy and Radiative Transfer,
        // Volume 65, Issue 5,
        // 2000,
        // Pages 769-778,
        const double sigma_tr = face_sigma_tr[cell][face];
        face_D[cell][face] = flux_limiter ? 1.0 / std::sqrt(9.0 * sigma_tr * sigma_tr +
                                                            (deden * deden) / (max_eden * max_eden))
                                          : 1.0 / (3.0 * sigma_tr);
      } else {
        Check(ftype == FACE_TYPE::BOUNDARY_FACE);
        Check(solver_data.off_diagonal_id[cell][face] == ncells);

        const double bound_eden = constants::a * std::pow(bnd_temp[face], 4.0);
        const auto max_eden = std::max(cell_eden, bound_eden);
        const double deden = std::abs(cell_eden - bound_eden) / half_width;
        // J.E. Morel,
        // Diffusion-limit asymptotics of the transport equation, the P1/3 equations, and two
        // flux-limited diffusion theories,
        // Journal of Quantitative Spectroscopy and Radiative Transfer,
        // Volume 65, Issue 5,
        // 2000,
        // Pages 769-778,
        const double sigma_tr = face_sigma_tr[cell][face];
        face_D[cell][face] = flux_limiter ? 1.0 / std::sqrt(9.0 * sigma_tr * sigma_tr +
                                                            (deden * deden) / (max_eden * max_eden))
                                          : 1.0 / (3.0 * sigma_tr);
      }
      if (flux_limiter)
        face_D[cell][face] =
            std::min(face_D[cell][face], 1.0 / (20 * half_width + constants::c * dt));
    }
  }

  if (gcomm) {
    // populate remaining face data
    gcomm->exchange_ghost_data(put_face_D, ghost_face_D);
  }
  // Loop over all cells
  for (size_t cell = 0; cell < ncells; cell++) {
    // Populate the homogenized cell material data
    const auto cell_volume = mesh.cell_volume(cell);
    // Apply collision terms to the diagonal (census + absorption)
    solver_data.diagonal[cell] = 1.0 + constants::c * dt * fleck[cell] * sigma_a[cell];
    // Apply the cell source
    solver_data.source[cell] =
        solver_data.cell_eden0[cell] +
        constants::a * constants::c * fleck[cell] * sigma_a[cell] *
            solver_data.cell_temperature0[cell] * solver_data.cell_temperature0[cell] *
            solver_data.cell_temperature0[cell] * solver_data.cell_temperature0[cell] * dt +
        ext_imp_source[cell] + rad_source[cell];

    const auto nfaces = mesh.number_of_faces(cell);
    for (size_t face = 0; face < nfaces; face++) {
      const auto ftype = mesh.face_type(cell, face);
      solver_data.face_type[cell][face] = ftype;
      // calculate the fringe values
      if (ftype == FACE_TYPE::GHOST_FACE) {
        const auto face_area = mesh.face_area(cell, face);
        const auto next_cell = solver_data.off_diagonal_id[cell][face];
        const auto next_face = face % 2 == 0 ? face + 1 : face - 1;
        const auto D = face_D[cell][face];
        const auto next_D = ghost_face_D[gcomm->ghost_map[next_cell][next_face]];
        const auto half_width = mesh.distance_center_to_face(cell, face);
        const auto next_half_width =
            solver_data.ghost_dist_center_to_face[gcomm->ghost_map[next_cell][next_face]];
        const double fring = -constants::c * dt * face_area * (D * next_D) /
                             (cell_volume * (half_width * next_D + next_half_width * D));

        solver_data.off_diagonal[cell][face] = fring;
        solver_data.diagonal[cell] -= fring;
      } else if (ftype == FACE_TYPE::INTERNAL_FACE) {
        const auto face_area = mesh.face_area(cell, face);
        const auto next_cell = solver_data.off_diagonal_id[cell][face];
        const auto next_face = face % 2 == 0 ? face + 1 : face - 1;
        const auto D = face_D[cell][face];
        const auto next_D = face_D[next_cell][next_face];
        const auto half_width = mesh.distance_center_to_face(cell, face);
        const auto next_half_width = mesh.distance_center_to_face(next_cell, next_face);
        const double fring = -constants::c * dt * face_area * (D * next_D) /
                             (cell_volume * (half_width * next_D + next_half_width * D));

        solver_data.off_diagonal[cell][face] = fring;
        solver_data.diagonal[cell] -= fring;
      } else {
        Check(ftype == FACE_TYPE::BOUNDARY_FACE);
        Check(solver_data.off_diagonal_id[cell][face] == ncells);
        if (!reflect_bnd[face]) {
          const double E0 = constants::a * std::pow(bnd_temp[face], 4.0);
          const auto D = face_D[cell][face];
          const auto half_width = mesh.distance_center_to_face(cell, face);
          const auto face_area = mesh.face_area(cell, face);
          const double fring = constants::c * D * dt * face_area /
                               (cell_volume * half_width * (1.0 + 2.0 * D / half_width));
          const double flux_source = fring * E0;
          solver_data.source[cell] += flux_source;
          solver_data.diagonal[cell] += fring;
          solver_data.off_diagonal[cell][face] = fring;
          solver_data.flux_source[cell][face] = flux_source;
        }
      }
    }
  }
}

//================================================================================================//
/*!
 * \brief Solver the radiation energy vector using Gauss-Siedel
 *
 *  \param[in] eps max convergence error
 *  \param[in] max_iter maximum number of iterations
 *
 */
//================================================================================================//
void Grey_Matrix::gs_solver(const double eps, const size_t max_iter, const bool print) {
  Require(max_iter > 0);
  double max_error = 1.0;
  std::stringstream diagnostics;
  std::vector<double> last_ghost_eden;
  std::map<size_t, std::vector<double>> put_ghost_cell_eden;
  if (gcomm) {
    last_ghost_eden = solver_data.ghost_cell_eden;
    for (auto &map : gcomm->put_buffer_size) {
      put_ghost_cell_eden[map.first] = std::vector<double>(map.second, 0.0);
    }
  }
  size_t global_count = 0;
  size_t total_inners = 0;
  double local_max_error = 0.0;
  while (max_error > eps && global_count < max_iter) {
    size_t inner_count = 0;
    while (max_error > eps && inner_count < max_iter) {
      max_error = 0.0;
      for (size_t i = 0; i < solver_data.diagonal.size(); i++) {
        const double diag = solver_data.diagonal[i];
        double b = solver_data.source[i] + cell_correction_source[i];
        for (size_t f = 0; f < solver_data.off_diagonal_id[i].size(); f++) {
          const auto ftype = solver_data.face_type[i][f];
          const size_t next_cell = solver_data.off_diagonal_id[i][f];
          if (ftype == FACE_TYPE::INTERNAL_FACE) {
            b -= solver_data.off_diagonal[i][f] * solver_data.cell_eden[next_cell];
          } else if (ftype == FACE_TYPE::GHOST_FACE) {
            const auto next_face = f % 2 == 0 ? f + 1 : f - 1;
            b -= solver_data.off_diagonal[i][f] *
                 last_ghost_eden[gcomm->ghost_map[next_cell][next_face]];
          }
        }
        Check(diag > 0.0);
        const double last_eden = solver_data.cell_eden[i];
        solver_data.cell_eden[i] = b / diag;
        Correction::calc_correction(
            cell_epsilon[i], cell_correction_source[i], solver_data.cell_temperature[i],
            solver_data.cell_eden[i], sigma_a[i], sigma_a[i], fleck[i], solver_data.cell_cve[i],
            volume[i], current_dt, solver_data.cell_temperature0[i], ext_exp_source[i]);
        Check(!(solver_data.cell_eden[i] < 0.0));
        const double eden_diff =
            solver_data.cell_eden[i] > 0.0
                ? std::abs(solver_data.cell_eden[i] - last_eden) / solver_data.cell_eden[i]
                : std::abs(solver_data.cell_eden[i] - last_eden);
        max_error = std::max(max_error, eden_diff);
      }

      /* HANGS FOR DIFFERENT NUMBER OF INNER SOLVES
      // Inner loop communication increased mpi window calls but drastically reduces the total
      // number of iterations
      if (gcomm) {
        // update the put vectors
        for (auto &map : gcomm->put_map) {
          auto cell = map.first;
          for (auto &put_face : map.second) {
            auto rank = put_face.second.first;
            auto buffer_index = put_face.second.second;
            put_ghost_cell_eden[rank][buffer_index] = solver_data.cell_eden[cell];
          }
        }
        // write to the ghost window data
        gcomm->exchange_ghost_data(put_ghost_cell_eden, solver_data.ghost_cell_eden);
        for (size_t i = 0; i < solver_data.ghost_cell_eden.size(); i++) {
          const double eden_diff =
              solver_data.ghost_cell_eden[i] > 0.0
                  ? std::abs(solver_data.ghost_cell_eden[i] - last_ghost_eden[i]) /
                        solver_data.ghost_cell_eden[i]
                  : std::abs(solver_data.cell_eden[i] - last_ghost_eden[i]);
          max_error = std::max(max_error, eden_diff);
        }
        last_ghost_eden = solver_data.ghost_cell_eden;
      }
      */
      inner_count++;
      total_inners++;
    }

    if (gcomm) {
      rtt_c4::global_max(inner_count);
      // update the put vectors
      for (auto &map : gcomm->put_map) {
        auto cell = map.first;
        for (auto &put_face : map.second) {
          auto rank = put_face.second.first;
          auto buffer_index = put_face.second.second;
          put_ghost_cell_eden[rank][buffer_index] = solver_data.cell_eden[cell];
        }
      }
      // echange ghost data
      gcomm->exchange_ghost_data(put_ghost_cell_eden, solver_data.ghost_cell_eden);
      // calculate difference on exchange
      for (size_t i = 0; i < solver_data.ghost_cell_eden.size(); i++) {
        const double eden_diff =
            solver_data.ghost_cell_eden[i] > 0.0
                ? std::abs(solver_data.ghost_cell_eden[i] - last_ghost_eden[i]) /
                      solver_data.ghost_cell_eden[i]
                : std::abs(solver_data.ghost_cell_eden[i] - last_ghost_eden[i]);
        max_error = std::max(max_error, eden_diff);
      }
      last_ghost_eden = solver_data.ghost_cell_eden;
      local_max_error = max_error;
      rtt_c4::global_max(max_error);
    }
    global_count++;
  }
  diagnostics << "Node = " << rtt_c4::node() << " Glogal Iteration = " << global_count
              << " Total Inner Iteration = " << total_inners << " error = " << local_max_error
              << "\n";

  if (gcomm) {
    rtt_c4::global_max(total_inners);
  }
  if (max_error > eps) {
    std::cout << "WARNING: Did not converge cell_eden max_error = " << max_error << std::endl;
    std::cout << diagnostics.str();
  } else if (print && rtt_c4::node() == 0) {
    std::cout << "Converged eden -> Global_Iteration = " << global_count
              << " Total_Inner_Iterations = " << total_inners << " error = " << max_error
              << std::endl;
  }
}

//================================================================================================//
/*!
 * \brief Calculate the output data 
 *
 *  \param[in] mesh orthogonal mesh class
 *  \param[in] mat_data material interface data
 *  \param[in] dt time step size
 *  \param[inout] output_data the output interface data 
 *
 */
//================================================================================================//
void Grey_Matrix::calculate_output_data(const Orthogonal_Mesh &mesh, const Mat_Data &mat_data,
                                        const double dt, Output_Data &output_data) {
  const auto ncells = mesh.number_of_local_cells();
  Require(output_data.cell_mat_dedv.size() == ncells);
  Require(output_data.cell_rad_eden.size() == ncells);
  Require(output_data.face_flux.size() == ncells);
  Require(solver_data.off_diagonal.size() == ncells);
  Require(solver_data.flux_source.size() == ncells);
  Opacity_Reader opacity_reader(mat_data.ipcress_filename, mat_data.problem_matids);
  // Loop over all cells
  for (size_t cell = 0; cell < ncells; cell++) {
    // compute the flux
    const auto nfaces = mesh.number_of_faces(cell);
    const auto volume = mesh.cell_volume(cell);
    for (size_t face = 0; face < nfaces; face++) {
      const auto area = mesh.face_area(cell, face);
      const auto ftype = mesh.face_type(cell, face);
      if (ftype == FACE_TYPE::INTERNAL_FACE) {
        const auto next_cell = solver_data.off_diagonal_id[cell][face];
        // calculate the current flux
        output_data.face_flux[cell][face] =
            (solver_data.off_diagonal[cell][face] *
                 (solver_data.cell_eden[next_cell] - solver_data.cell_eden[cell]) +
             solver_data.flux_source[cell][face]) *
            volume / (area * dt);
      } else if (ftype == FACE_TYPE::GHOST_FACE) {
        const auto next_face = face % 2 == 0 ? face + 1 : face - 1;
        const auto next_cell = solver_data.off_diagonal_id[cell][face];
        // calculate the current flux
        output_data.face_flux[cell][face] =
            (solver_data.off_diagonal[cell][face] *
                 (solver_data.ghost_cell_eden[gcomm->ghost_map[next_cell][next_face]] -
                  solver_data.cell_eden[cell]) +
             solver_data.flux_source[cell][face]) *
            volume / (area * dt);
      } else if (ftype == FACE_TYPE::BOUNDARY_FACE) {
        output_data.face_flux[cell][face] =
            (solver_data.flux_source[cell][face] -
             solver_data.off_diagonal[cell][face] * solver_data.cell_eden[cell]) *
            volume / (area * dt);
      }
    }
    // compute the remaining output data
    output_data.cell_rad_eden[cell] = solver_data.cell_eden[cell];
    const double cell_vol_edep = (fleck[cell] + cell_epsilon[cell]) * sigma_a[cell] * constants::c *
                                 solver_data.cell_eden[cell] * dt;
    // Populate the homogenized cell material data
    const auto nmats = mat_data.number_of_cell_mats[cell];
    Check(output_data.cell_mat_dedv[cell].size() == nmats);
    // Split up the change in material energy between materials
    for (size_t mat = 0; mat < nmats; mat++) {
      const size_t matid = mat_data.cell_mats[cell][mat];
      const double mat_sigma_a =
          opacity_reader.mat_planck_abs_models[matid]->getOpacity(
              mat_data.cell_mat_temperature[cell][mat], mat_data.cell_mat_density[cell][mat]) *
          mat_data.cell_mat_density[cell][mat];
      const double mat_vol_edep = cell_vol_edep * mat_sigma_a / sigma_a[cell];
      const double mat_T4 =
          mat_data.cell_mat_temperature[cell][mat] * mat_data.cell_mat_temperature[cell][mat] *
          mat_data.cell_mat_temperature[cell][mat] * mat_data.cell_mat_temperature[cell][mat];
      const double mat_vol_emission = (fleck[cell] + cell_epsilon[cell]) * mat_sigma_a *
                                      constants::a * constants::c * mat_T4 * dt;
      output_data.cell_mat_dedv[cell][mat] =
          (mat_vol_edep - mat_vol_emission) +
          mat_data.cell_mat_electron_source[cell][mat] * (fleck[cell] + cell_epsilon[cell]) * dt;
    }
  }
}

} // namespace odd_solver

//------------------------------------------------------------------------------------------------//
// end of Grey_Matrix.cc
//------------------------------------------------------------------------------------------------//
