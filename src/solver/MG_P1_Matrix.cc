//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   solver/MG_P1_Matrix.cc
 * \author Mathew Cleveland
 * \date   January 5th 2022
 * \brief  Build Matrix data to be used by the solver
 * \note   Copyright (C) 2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "MG_P1_Matrix.hh"
#include "Constants.hh"
#include "Correction.hh"
#include "Opacity_Reader.hh"
#include "c4/global.hh"
#include "cdi/CDI.hh"
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
MG_P1_Matrix::MG_P1_Matrix(const Control_Data &control_data)
    : reflect_bnd(control_data.reflect_bnd), bnd_temp(control_data.bnd_temp),
      correction(control_data.correction) {
  Insist(control_data.multigroup, "Must be multigroup");
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
double MG_P1_Matrix::mass_average(const std::vector<double> &mass,
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
void MG_P1_Matrix::initialize_solver_data(const Orthogonal_Mesh &mesh, const Mat_Data &mat_data,
                                          const double dt) {
  Opacity_Reader opacity_reader(mat_data.ipcress_filename, mat_data.problem_matids);
  group_bounds = opacity_reader.group_bounds;
  ngroups = group_bounds.size() - 1;
  const auto ncells = mesh.number_of_local_cells();
  current_dt = dt;
  Require(mat_data.face_flux.size() == ncells);
  // declare some dd put data structures
  std::map<size_t, std::vector<double>> put_ghost_cell_temp;
  std::map<size_t, std::vector<double>> put_ghost_cell_eden;
  std::map<size_t, std::vector<double>> put_dist_center_to_face;
  if (mesh.domain_decomposed()) {
    gcomm = std::make_unique<Ghost_Comm>(Ghost_Comm(mesh));
    const auto nghost = gcomm->local_ghost_buffer_size;
    // Initialize dd receive (ghost) data structures
    solver_data.ghost_cell_temperature.resize(nghost, 0.0);
    solver_data.ghost_cell_eden.resize(nghost * ngroups, 0.0);
    solver_data.ghost_dist_center_to_face.resize(nghost, 0.0);
    // Initialize dd put data structures
    for (auto &map : gcomm->put_buffer_size) {
      put_ghost_cell_temp[map.first] = std::vector<double>(map.second, 0.0);
      put_ghost_cell_eden[map.first] = std::vector<double>(map.second * ngroups, 0.0);
      put_dist_center_to_face[map.first] = std::vector<double>(map.second, 0.0);
    }
  }
  // Allocate solver cell data
  solver_data.diagonal.resize(ncells, std::vector<double>(ngroups, 0.0));
  solver_data.source.resize(ncells, std::vector<double>(ngroups, 0.0));
  solver_data.cell_density.resize(ncells, 0.0);
  solver_data.cell_cve.resize(ncells, 0.0);
  solver_data.cell_eden0.resize(ncells, std::vector<double>(ngroups, 0.0));
  solver_data.cell_temperature0.resize(ncells, 0.0);
  solver_data.cell_eden.resize(ncells, std::vector<double>(ngroups, 0.0));
  solver_data.cell_temperature.resize(ncells, 0.0);
  // Allocate solver face/neighbor data
  solver_data.off_diagonal.resize(ncells);
  solver_data.flux_source.resize(ncells);
  solver_data.face_type.resize(ncells);
  solver_data.off_diagonal_id.resize(ncells);
  solver_data.face_flux0.resize(ncells);
  // Resize local matrix data
  sigma_a.resize(ncells, std::vector<double>(ngroups, 0.0));
  sigma_planck.resize(ncells, 0.0);
  ext_imp_source.resize(ncells, 0.0);
  ext_exp_source.resize(ncells, 0.0);
  rad_source.resize(ncells, 0.0);
  fleck.resize(ncells, 0.0);
  cell_epsilon.resize(ncells, 0.0);
  cell_correction_source.resize(ncells, 0.0);
  cell_planck_weights.resize(ncells, std::vector<double>(ngroups, 0.0));
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
    std::vector<std::vector<double>> cell_mat_sigma_tr(nmats, std::vector<double>(ngroups, 0.0));
    std::vector<std::vector<double>> cell_mat_sigma_a(nmats, std::vector<double>(ngroups, 0.0));
    std::vector<double> cell_mat_mass(nmats, 0.0);
    std::vector<double> cell_mat_cv(nmats, 0.0);
    std::vector<double> cell_mat_T4(nmats, 0.0);
    std::vector<double> cell_mat_emission(nmats, 0.0);
    std::vector<double> cell_mat_sigma_planck(nmats, 0.0);
    std::vector<double> planck_spec(ngroups, 0.0);
    for (size_t mat = 0; mat < nmats; mat++) {
      size_t matid = mat_data.cell_mats[cell][mat];
      cell_mat_sigma_a[mat] = opacity_reader.mat_mg_planck_abs_models[matid]->getOpacity(
          mat_data.cell_mat_temperature[cell][mat], mat_data.cell_mat_density[cell][mat]);
      std::transform(cell_mat_sigma_a[mat].begin(), cell_mat_sigma_a[mat].end(),
                     cell_mat_sigma_a[mat].begin(),
                     std::bind1st(std::multiplies<double>(), mat_data.cell_mat_density[cell][mat]));
      cell_mat_sigma_tr[mat] = opacity_reader.mat_mg_rosseland_total_models[matid]->getOpacity(
          mat_data.cell_mat_temperature[cell][mat], mat_data.cell_mat_density[cell][mat]);
      std::transform(cell_mat_sigma_tr[mat].begin(), cell_mat_sigma_tr[mat].end(),
                     cell_mat_sigma_tr[mat].begin(),
                     std::bind1st(std::multiplies<double>(), mat_data.cell_mat_density[cell][mat]));
      cell_mat_mass[mat] = mat_data.cell_mat_density[cell][mat] *
                           mat_data.cell_mat_vol_frac[cell][mat] * cell_volume;
      cell_mat_cv[mat] =
          mat_data.cell_mat_density[cell][mat] * mat_data.cell_mat_specific_heat[cell][mat];
      cell_mat_T4[mat] =
          mat_data.cell_mat_temperature[cell][mat] * mat_data.cell_mat_temperature[cell][mat] *
          mat_data.cell_mat_temperature[cell][mat] * mat_data.cell_mat_temperature[cell][mat];
      rtt_cdi::CDI::integrate_Planckian_Spectrum(
          group_bounds, mat_data.cell_mat_temperature[cell][mat], planck_spec);
      cell_mat_sigma_planck[mat] = std::inner_product(
          cell_mat_sigma_a[mat].begin(), cell_mat_sigma_a[mat].end(), planck_spec.begin(), 0.0);
      cell_mat_emission[mat] = cell_mat_T4[mat] * cell_mat_sigma_planck[mat];
    }
    const double cell_mass = std::accumulate(cell_mat_mass.begin(), cell_mat_mass.end(), 0.0);
    Check(cell_mass > 0.0);
    const double cell_density = cell_mass / cell_volume;
    Check(cell_density > 0.0);
    const double cell_cve = std::inner_product(cell_mat_cv.begin(), cell_mat_cv.end(),
                                               mat_data.cell_mat_vol_frac[cell].begin(), 0.0);
    Check(cell_cve > 0.0);
    sigma_a[cell] = std::vector<double>(ngroups, 0.0);
    const double mass_ave_sigma_p =
        std::inner_product(cell_mat_sigma_planck.begin(), cell_mat_sigma_planck.end(),
                           mat_data.cell_mat_vol_frac[cell].begin(), 0.0);

    // Calc mass averaged T4 (just in case the opacity is zero)
    const double mass_avg_T4 = mass_average(cell_mat_mass, cell_mat_T4);
    const double cell_T =
        mass_ave_sigma_p > 0.0
            ? std::pow(std::inner_product(cell_mat_emission.begin(), cell_mat_emission.end(),
                                          mat_data.cell_mat_vol_frac[cell].begin(), 0.0) /
                           mass_ave_sigma_p,
                       0.25)
            : std::pow(mass_avg_T4, 0.25);

    for (size_t mat = 0; mat < nmats; mat++) {
      for (size_t g = 0; g < ngroups; g++)
        sigma_a[cell][g] += cell_mat_sigma_a[mat][g] * mat_data.cell_mat_vol_frac[cell][mat];
    }

    // calculate the planck
    rtt_cdi::CDI::integrate_Planckian_Spectrum(group_bounds, cell_T, planck_spec);
    cell_planck_weights[cell] = planck_spec;
    sigma_planck[cell] =
        std::inner_product(sigma_a[cell].begin(), sigma_a[cell].end(), planck_spec.begin(), 0.0);

    solver_data.cell_density[cell] = cell_density;
    solver_data.cell_cve[cell] = cell_cve;
    Check(cell_cve > 0.0);
    fleck[cell] = 1.0 / (1.0 + constants::a * constants::c * sigma_planck[cell] * 4.0 * cell_T *
                                   cell_T * cell_T * dt / cell_cve);

    const double src_val = std::inner_product(mat_data.cell_mat_electron_source[cell].begin(),
                                              mat_data.cell_mat_electron_source[cell].end(),
                                              mat_data.cell_mat_vol_frac[cell].begin(), 0.0);

    ext_imp_source[cell] = src_val * (1.0 - fleck[cell]);
    ext_exp_source[cell] = src_val * fleck[cell];
    rad_source[cell] = mat_data.cell_rad_source[cell];

    // Set initial conditions
    solver_data.cell_temperature0[cell] = cell_T;
    solver_data.cell_temperature[cell] = cell_T;
    solver_data.cell_eden0[cell] = mat_data.cell_rad_mg_eden[cell];
    solver_data.cell_eden[cell] = mat_data.cell_rad_mg_eden[cell];
    // Allocate face data
    const auto nfaces = mesh.number_of_faces(cell);
    solver_data.off_diagonal[cell].resize(nfaces, std::vector<double>(ngroups, 0.0));
    solver_data.flux_source[cell].resize(nfaces, std::vector<double>(ngroups, 0.0));
    solver_data.face_type[cell].resize(nfaces, 0);
    solver_data.face_flux0[cell].resize(nfaces, std::vector<double>(ngroups, 0.0));
    solver_data.off_diagonal_id[cell].resize(nfaces, ncells);
    face_sigma_tr[cell].resize(nfaces, std::vector<double>(ngroups, 0.0));
    // Loop over faces
    for (size_t face = 0; face < nfaces; face++) {
      solver_data.face_flux0[cell][face] = mat_data.face_mg_flux[cell][face];
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

        for (size_t g = 0; g < ngroups; g++)
          put_ghost_cell_eden[rank][ngroups * buffer_index + g] = solver_data.cell_eden[cell][g];
      }
    }
    gcomm->exchange_ghost_data(put_ghost_cell_temp, solver_data.ghost_cell_temperature);
    gcomm->exchange_ghost_data(put_ghost_cell_eden, solver_data.ghost_cell_eden, ngroups);
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
        std::vector<double> cell_mat_sigma_tr(ngroups, 0.0);
        std::vector<double> sigma_tr(ngroups, 0.0);
        for (size_t mat = 0; mat < nmats; mat++) {
          size_t matid = mat_data.cell_mats[cell][mat];
          cell_mat_sigma_tr = opacity_reader.mat_mg_rosseland_total_models[matid]->getOpacity(
              face_T, mat_data.cell_mat_density[cell][mat]);
          for (size_t g = 0; g < ngroups; g++)
            sigma_tr[g] += cell_mat_sigma_tr[g] * mat_data.cell_mat_vol_frac[cell][mat] *
                           mat_data.cell_mat_density[cell][mat];
        }

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
        std::vector<double> cell_mat_sigma_tr(ngroups, 0.0);
        std::vector<double> sigma_tr(ngroups, 0.0);
        for (size_t mat = 0; mat < nmats; mat++) {
          size_t matid = mat_data.cell_mats[cell][mat];
          cell_mat_sigma_tr = opacity_reader.mat_mg_rosseland_total_models[matid]->getOpacity(
              face_T, mat_data.cell_mat_density[cell][mat]);
          for (size_t g = 0; g < ngroups; g++)
            sigma_tr[g] += cell_mat_sigma_tr[g] * mat_data.cell_mat_vol_frac[cell][mat] *
                           mat_data.cell_mat_density[cell][mat];
        }
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
        std::vector<double> cell_mat_sigma_tr(ngroups, 0.0);
        std::vector<double> sigma_tr(ngroups, 0.0);
        for (size_t mat = 0; mat < nmats; mat++) {
          size_t matid = mat_data.cell_mats[cell][mat];
          cell_mat_sigma_tr = opacity_reader.mat_mg_rosseland_total_models[matid]->getOpacity(
              face_T, mat_data.cell_mat_density[cell][mat]);
          for (size_t g = 0; g < ngroups; g++)
            sigma_tr[g] += cell_mat_sigma_tr[g] * mat_data.cell_mat_vol_frac[cell][mat] *
                           mat_data.cell_mat_density[cell][mat];
        }
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
void MG_P1_Matrix::build_matrix(const Orthogonal_Mesh &mesh, const double dt) {
  const auto ncells = mesh.number_of_local_cells();

  std::map<size_t, std::vector<double>> put_face_D;
  if (gcomm) {
    const auto nghost = gcomm->local_ghost_buffer_size;
    // Initialize dd receive (ghost) data structures
    ghost_face_D.resize(nghost * ngroups, 0.0);
    // Initialize dd put data structures
    for (auto &map : gcomm->put_buffer_size) {
      put_face_D[map.first] = std::vector<double>(map.second * ngroups, 0.0);
    }
  }
  // Build face diffusion coefficients
  face_D.resize(ncells);
  for (size_t cell = 0; cell < ncells; cell++) {
    const auto nfaces = mesh.number_of_faces(cell);
    face_D[cell].resize(nfaces, std::vector<double>(ngroups, 0.0));
    for (size_t face = 0; face < nfaces; face++) {
      const auto ftype = mesh.face_type(cell, face);
      if (ftype == FACE_TYPE::GHOST_FACE) {
        const auto sigma_tr = face_sigma_tr[cell][face];
        const auto rank = gcomm->put_map[cell][face].first;
        const auto put_index = gcomm->put_map[cell][face].second;
        for (size_t g = 0; g < ngroups; g++) {
          face_D[cell][face][g] =
              (constants::c * dt) / (3.0 * (1.0 + constants::c * dt * sigma_tr[g]));
          // populate the put buffer
          put_face_D[rank][put_index * ngroups + g] = face_D[cell][face][g];
        }

      } else if (ftype == FACE_TYPE::INTERNAL_FACE) {
        for (size_t g = 0; g < ngroups; g++)
          face_D[cell][face][g] = (constants::c * dt) /
                                  (3.0 * (1.0 + constants::c * dt * face_sigma_tr[cell][face][g]));

      } else {
        Check(ftype == FACE_TYPE::BOUNDARY_FACE);
        Check(solver_data.off_diagonal_id[cell][face] == ncells);
        for (size_t g = 0; g < ngroups; g++)
          face_D[cell][face][g] = (constants::c * dt) /
                                  (3.0 * (1.0 + constants::c * dt * face_sigma_tr[cell][face][g]));
      }
    }
  }

  if (gcomm) {
    // populate remaining face data
    gcomm->exchange_ghost_data(put_face_D, ghost_face_D, ngroups);
  }
  // Loop over all cells
  for (size_t cell = 0; cell < ncells; cell++) {
    // Populate the homogenized cell material data
    const auto cell_volume = mesh.cell_volume(cell);
    // Apply collision terms to the diagonal (census + absorption)
    std::vector<double> planck_spec(ngroups, 0.0);
    rtt_cdi::CDI::integrate_Planckian_Spectrum(group_bounds, solver_data.cell_temperature0[cell],
                                               planck_spec);
    const double total_emission =
        constants::a * constants::c * fleck[cell] * sigma_planck[cell] *
        solver_data.cell_temperature0[cell] * solver_data.cell_temperature0[cell] *
        solver_data.cell_temperature0[cell] * solver_data.cell_temperature0[cell] * dt;
    for (size_t g = 0; g < ngroups; g++) {
      solver_data.diagonal[cell][g] = 1.0 + constants::c * dt * fleck[cell] * sigma_a[cell][g];
      // Apply the cell source
      solver_data.source[cell][g] =
          solver_data.cell_eden0[cell][g] + total_emission * planck_spec[g] +
          ext_imp_source[cell] * planck_spec[g] + rad_source[cell] * planck_spec[g];
    }

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
        const auto next_D = std::vector<double>(
            ghost_face_D.begin() + gcomm->ghost_map[next_cell][next_face],
            ghost_face_D.begin() + gcomm->ghost_map[next_cell][next_face] + ngroups);
        const auto half_width = mesh.distance_center_to_face(cell, face);
        const auto next_half_width =
            solver_data.ghost_dist_center_to_face[gcomm->ghost_map[next_cell][next_face]];
        for (size_t g = 0; g < ngroups; g++) {
          const double fring = -constants::c * dt * face_area * (D[g] * next_D[g]) /
                               (cell_volume * (half_width * next_D[g] + next_half_width * D[g]));
          const double flux_source = -fring * 3.0 * (half_width + next_half_width) /
                                     (constants::c * dt) * solver_data.face_flux0[cell][face][g] /
                                     constants::c;

          solver_data.off_diagonal[cell][face][g] = fring;
          solver_data.flux_source[cell][face][g] = flux_source;
          solver_data.source[cell][g] -= flux_source;
          solver_data.diagonal[cell][g] -= fring;
        }
      } else if (ftype == FACE_TYPE::INTERNAL_FACE) {
        const auto face_area = mesh.face_area(cell, face);
        const auto next_cell = solver_data.off_diagonal_id[cell][face];
        const auto next_face = face % 2 == 0 ? face + 1 : face - 1;
        const auto D = face_D[cell][face];
        const auto next_D = face_D[next_cell][next_face];
        const auto half_width = mesh.distance_center_to_face(cell, face);
        const auto next_half_width = mesh.distance_center_to_face(next_cell, next_face);
        for (size_t g = 0; g < ngroups; g++) {
          const double fring = -constants::c * dt * face_area * (D[g] * next_D[g]) /
                               (cell_volume * (half_width * next_D[g] + next_half_width * D[g]));
          const double flux_source = -fring * 3.0 * (half_width + next_half_width) /
                                     (constants::c * dt) * solver_data.face_flux0[cell][face][g] /
                                     constants::c;

          solver_data.off_diagonal[cell][face][g] = fring;
          solver_data.flux_source[cell][face][g] = flux_source;
          solver_data.source[cell][g] -= flux_source;
          solver_data.diagonal[cell][g] -= fring;
        }
      } else {
        Check(ftype == FACE_TYPE::BOUNDARY_FACE);
        Check(solver_data.off_diagonal_id[cell][face] == ncells);
        if (!reflect_bnd[face]) {
          const double E0 = constants::a * std::pow(bnd_temp[face], 4.0);
          const auto D = face_D[cell][face];
          const auto half_width = mesh.distance_center_to_face(cell, face);
          const auto face_area = mesh.face_area(cell, face);
          for (size_t g = 0; g < ngroups; g++) {
            const double fring = constants::c * D[g] * dt * face_area /
                                 (cell_volume * half_width * (1.0 + 2.0 * D[g] / half_width));
            const double flux_source = fring * 3.0 * half_width / (constants::c * dt) *
                                       solver_data.face_flux0[cell][face][g] / constants::c;
            solver_data.source[cell][g] += fring * E0;
            solver_data.source[cell][g] += flux_source;
            solver_data.diagonal[cell][g] += fring;
            solver_data.off_diagonal[cell][face][g] = -fring;
            solver_data.flux_source[cell][face][g] = fring * E0 + flux_source;
          }
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
void MG_P1_Matrix::gs_solver(const double eps, const size_t max_iter, const bool print) {
  Require(max_iter > 0);
  double max_error = 1.0;
  std::stringstream diagnostics;
  std::vector<double> last_ghost_eden;
  std::map<size_t, std::vector<double>> put_ghost_cell_eden;
  if (gcomm) {
    last_ghost_eden = solver_data.ghost_cell_eden;
    for (auto &map : gcomm->put_buffer_size) {
      put_ghost_cell_eden[map.first] = std::vector<double>(map.second * ngroups, 0.0);
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
        const std::vector<double> last_mg_eden = solver_data.cell_eden[i];
        double new_cell_eden = 0.0;
        double sigma_dep = 0.0;
        for (size_t g = 0; g < ngroups; g++) {
          const double diag = solver_data.diagonal[i][g];
          double b =
              solver_data.source[i][g] + cell_correction_source[i] * cell_planck_weights[i][g];
          for (size_t f = 0; f < solver_data.off_diagonal_id[i].size(); f++) {
            const auto ftype = solver_data.face_type[i][f];
            const size_t next_cell = solver_data.off_diagonal_id[i][f];
            if (ftype == FACE_TYPE::INTERNAL_FACE) {
              b -= solver_data.off_diagonal[i][f][g] * solver_data.cell_eden[next_cell][g];
            } else if (ftype == FACE_TYPE::GHOST_FACE) {
              const auto next_face = f % 2 == 0 ? f + 1 : f - 1;
              b -= solver_data.off_diagonal[i][f][g] *
                   last_ghost_eden[(gcomm->ghost_map[next_cell][next_face]) * ngroups + g];
            }
          }
          Check(diag > 0.0);
          solver_data.cell_eden[i][g] = b / diag;
          new_cell_eden += solver_data.cell_eden[i][g];
          sigma_dep += solver_data.cell_eden[i][g] * sigma_a[i][g];
        }
        sigma_dep /= new_cell_eden;
        if (correction)
          Correction::calc_correction(
              cell_epsilon[i], cell_correction_source[i], solver_data.cell_temperature[i],
              new_cell_eden, sigma_dep, sigma_planck[i], fleck[i], solver_data.cell_cve[i],
              volume[i], current_dt, solver_data.cell_temperature0[i], ext_exp_source[i]);
        Check(!(new_cell_eden < 0.0));
        // Check convergence
        for (size_t g = 0; g < ngroups; g++) {
          const double eden_diff = solver_data.cell_eden[i][g] > 0.0
                                       ? std::abs(solver_data.cell_eden[i][g] - last_mg_eden[g]) /
                                             solver_data.cell_eden[i][g]
                                       : std::abs(solver_data.cell_eden[i][g] - last_mg_eden[g]);
          max_error = std::max(max_error, eden_diff);
        }
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
          for (size_t g = 0; g < ngroups; g++)
            put_ghost_cell_eden[rank][ngroups * buffer_index + g] = solver_data.cell_eden[cell][g];
        }
      }
      // echange ghost data
      gcomm->exchange_ghost_data(put_ghost_cell_eden, solver_data.ghost_cell_eden, ngroups);
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
void MG_P1_Matrix::calculate_output_data(const Orthogonal_Mesh &mesh, const Mat_Data &mat_data,
                                         const double dt, Output_Data &output_data) {
  const auto ncells = mesh.number_of_local_cells();
  Require(output_data.cell_mat_dedv.size() == ncells);
  Require(output_data.cell_rad_eden.size() == ncells);
  Require(output_data.cell_rad_mg_eden.size() == ncells);
  Require(output_data.face_flux.size() == ncells);
  Require(output_data.face_mg_flux.size() == ncells);
  Require(solver_data.off_diagonal.size() == ncells);
  Require(solver_data.flux_source.size() == ncells);
  Opacity_Reader opacity_reader(mat_data.ipcress_filename, mat_data.problem_matids);
  for (size_t cell = 0; cell < ncells; cell++) {
    // compute the flux
    const auto nfaces = mesh.number_of_faces(cell);
    const auto volume = mesh.cell_volume(cell);
    for (size_t face = 0; face < nfaces; face++) {
      const auto area = mesh.face_area(cell, face);
      const auto ftype = mesh.face_type(cell, face);
      output_data.face_flux[cell][face] = 0.0;
      if (ftype == FACE_TYPE::INTERNAL_FACE) {
        const auto next_cell = solver_data.off_diagonal_id[cell][face];
        // calculate the current flux
        for (size_t g = 0; g < ngroups; g++) {
          output_data.face_mg_flux[cell][face][g] =
              (solver_data.off_diagonal[cell][face][g] *
                   (solver_data.cell_eden[next_cell][g] - solver_data.cell_eden[cell][g]) +
               solver_data.flux_source[cell][face][g]) *
              volume / (area * dt);
          output_data.face_flux[cell][face] += output_data.face_mg_flux[cell][face][g];
        }
      } else if (ftype == FACE_TYPE::GHOST_FACE) {
        const auto next_face = face % 2 == 0 ? face + 1 : face - 1;
        const auto next_cell = solver_data.off_diagonal_id[cell][face];
        // calculate the current flux
        for (size_t g = 0; g < ngroups; g++) {
          output_data.face_mg_flux[cell][face][g] =
              (solver_data.off_diagonal[cell][face][g] *
                   (solver_data
                        .ghost_cell_eden[gcomm->ghost_map[next_cell][next_face] * ngroups + g] -
                    solver_data.cell_eden[cell][g]) +
               solver_data.flux_source[cell][face][g]) *
              volume / (area * dt);

          output_data.face_flux[cell][face] += output_data.face_mg_flux[cell][face][g];
        }
      } else if (ftype == FACE_TYPE::BOUNDARY_FACE) {
        for (size_t g = 0; g < ngroups; g++) {
          output_data.face_mg_flux[cell][face][g] =
              (solver_data.flux_source[cell][face][g] -
               solver_data.off_diagonal[cell][face][g] * solver_data.cell_eden[cell][g]) *
              volume / (area * dt);

          output_data.face_flux[cell][face] += output_data.face_mg_flux[cell][face][g];
        }
      }
    }
    // compute the remaining output data
    output_data.cell_rad_mg_eden[cell] = solver_data.cell_eden[cell];
    output_data.cell_rad_eden[cell] = std::accumulate(solver_data.cell_eden[cell].begin(),
                                                      solver_data.cell_eden[cell].end(), 0.0);
    const double cell_vol_edep =
        (fleck[cell] + cell_epsilon[cell]) * constants::c * dt *
        std::inner_product(solver_data.cell_eden[cell].begin(), solver_data.cell_eden[cell].end(),
                           sigma_a[cell].begin(), 0.0);
    // Populate the homogenized cell material data
    const auto nmats = mat_data.number_of_cell_mats[cell];
    std::vector<double> planck_spec(ngroups, 0.0);
    Check(output_data.cell_mat_dedv[cell].size() == nmats);
    // Split up the change in material energy between materials
    for (size_t mat = 0; mat < nmats; mat++) {
      const size_t matid = mat_data.cell_mats[cell][mat];
      auto mat_sigma_a = opacity_reader.mat_mg_planck_abs_models[matid]->getOpacity(
          mat_data.cell_mat_temperature[cell][mat], mat_data.cell_mat_density[cell][mat]);
      std::transform(mat_sigma_a.begin(), mat_sigma_a.end(), mat_sigma_a.begin(),
                     std::bind1st(std::multiplies<double>(), mat_data.cell_mat_density[cell][mat]));
      rtt_cdi::CDI::integrate_Planckian_Spectrum(
          group_bounds, mat_data.cell_mat_temperature[cell][mat], planck_spec);
      const double mat_planck =
          std::inner_product(mat_sigma_a.begin(), mat_sigma_a.end(), planck_spec.begin(), 0.0);

      const double mat_vol_edep = cell_vol_edep * mat_planck / sigma_planck[cell];
      const double mat_T4 =
          mat_data.cell_mat_temperature[cell][mat] * mat_data.cell_mat_temperature[cell][mat] *
          mat_data.cell_mat_temperature[cell][mat] * mat_data.cell_mat_temperature[cell][mat];
      const double mat_vol_emission = (fleck[cell] + cell_epsilon[cell]) * mat_planck *
                                      constants::a * constants::c * mat_T4 * dt;
      output_data.cell_mat_dedv[cell][mat] =
          (mat_vol_edep - mat_vol_emission) +
          mat_data.cell_mat_electron_source[cell][mat] * (fleck[cell] + cell_epsilon[cell]) * dt;
    }
  }
}

} // namespace odd_solver

//------------------------------------------------------------------------------------------------//
// end of MG_P1_Matrix.cc
//------------------------------------------------------------------------------------------------//
