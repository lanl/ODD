//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   api/test/tstInterface_Builder.cc
 * \author Mathew Cleveland
 * \date   October 21st 2021
 * \brief  Testing opacity reader
 * \note   Copyright (C) 2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "api/Arguments.hh"
#include "api/Interface_Builder.hh"
#include "solver/Interface_Data.hh"
#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/dbc.hh"

using namespace rtt_dsxx;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//
void test(rtt_dsxx::UnitTest &ut) {
  Arguments arg;
  // Should fail
  try {
    arg.control_data.check_arguments();
  } catch (...) {
    std::cout << "Caught the expected Control_Data::check_arguments assertions" << std::endl;
  }
  // fill control data with "valid junk" and check arguments
  std::string file = "not_empty";
  arg.control_data.opacity_file = &file[0];
  arg.control_data.dt = 0.1;
  arg.control_data.max_iter = 1;
  arg.control_data.min_tol = 1.0e-12;
  arg.control_data.check_arguments();

  try {
    arg.zonal_data.check_arguments();
  } catch (...) {
    std::cout << "Caught the expected Zonal_Data::check_arguments assertions" << std::endl;
  }
  // fill the zonal data with "valid junk" and check arguments
  arg.zonal_data.domain_decomposed = 0;
  arg.zonal_data.number_of_local_cells = 2;
  arg.zonal_data.number_of_global_cells = 2;
  arg.zonal_data.dimensions = 1;
  arg.zonal_data.coord_sys = 0;
  // 2 zones | 0 || 1 | with dx=0.5 dy=0 and dz=0
  std::vector<double> cell_position{0.25, 0.0, 0.0, 0.75, 0.0, 0.0};
  arg.zonal_data.cell_position = &cell_position[0];
  std::vector<double> cell_size{0.5, 0.0, 0.0, 0.5, 0.0, 0.0};
  arg.zonal_data.cell_size = &cell_size[0];
  std::vector<size_t> cell_global_id{0, 1};
  arg.zonal_data.cell_global_id = &cell_global_id[0];
  std::vector<size_t> face_type{
      odd_solver::FACE_TYPE::BOUNDARY_FACE, odd_solver::FACE_TYPE::INTERNAL_FACE,
      odd_solver::FACE_TYPE::INTERNAL_FACE, odd_solver::FACE_TYPE::BOUNDARY_FACE};
  arg.zonal_data.face_type = &face_type[0];
  std::vector<size_t> next_cell_id{2, 1, 0, 2};
  arg.zonal_data.next_cell_id = &next_cell_id[0];

  // global material data
  std::vector<size_t> matids = {19000, 190001};
  arg.zonal_data.number_of_mats = 2;
  arg.zonal_data.problem_matids = &matids[0];
  // cell wise material data
  // 1 material in cell 1 and 2 materials in cell 2
  std::vector<size_t> cell_number_of_mats{1, 2};
  arg.zonal_data.number_of_cell_mats = &cell_number_of_mats[0];
  // cell 1 (mat 1) cell 2 (mat 1 and mat 2)
  std::vector<size_t> cell_mats{0, 0, 1};
  std::vector<std::vector<size_t>> cell_mats_ref{{0}, {0, 1}};
  // cell 1 (1.0) cell 2 mat 1 (0.5) and mat 2 (0.5)
  std::vector<double> cell_mat_vol_frac{1.0, 0.5, 0.5};
  std::vector<std::vector<double>> cell_mat_vol_frac_ref{{1.0}, {0.5, 0.5}};
  // cell 1 (1.0) cell 2 mat 1 (10.0) and mat 2 (1.0)
  std::vector<double> cell_mat_temperature{1.0, 10.0, 1.0};
  std::vector<std::vector<double>> cell_mat_temperature_ref{{1.0}, {10.0, 1.0}};
  // cell 1 (10) cell 2 mat 1 (1.0) and mat 2 (10.0)
  std::vector<double> cell_mat_density{10.0, 1.0, 10.0};
  std::vector<std::vector<double>> cell_mat_density_ref{{10.0}, {1.0, 10.0}};
  // cell 1 (0.1) cell 2 mat 0 (0.1) and mat 1 (0.01)
  std::vector<double> cell_mat_specific_heat{0.1, 0.1, 0.01};
  std::vector<std::vector<double>> cell_mat_specific_heat_ref{{0.1}, {0.1, 0.01}};
  // cell velocity
  std::vector<double> cell_velocity{0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  std::vector<std::array<double, 3>> cell_velocity_ref{{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
  // cell radiation energy density
  std::vector<double> cell_rad_eden{32.0, 23.0};
  arg.zonal_data.cell_mats = &cell_mats[0];
  arg.zonal_data.cell_mat_vol_frac = &cell_mat_vol_frac[0];
  arg.zonal_data.cell_mat_temperature = &cell_mat_temperature[0];
  arg.zonal_data.cell_mat_density = &cell_mat_density[0];
  arg.zonal_data.cell_mat_specific_heat = &cell_mat_specific_heat[0];
  arg.zonal_data.cell_velocity = &cell_velocity[0];
  arg.zonal_data.cell_erad = &cell_rad_eden[0];

  // Check that the pointer matches the data
  FAIL_IF_NOT(arg.zonal_data.problem_matids[0] == matids[0]);
  FAIL_IF_NOT(arg.zonal_data.problem_matids[1] == matids[1]);
  arg.zonal_data.check_arguments();

  try {
    arg.output_data.check_arguments(arg.zonal_data.number_of_local_cells,
                                    arg.zonal_data.number_of_cell_mats);
  } catch (...) {
    std::cout << "Caught the expected Output_Data::check_arguments assertions" << std::endl;
  }
  std::vector<double> cell_erad = {17000.0, 17000.1};
  std::vector<double> cell_Trad = {18000.0, 18000.1};
  std::vector<double> cell_mat_delta_e = {-3.0, 4.0, -5.0};
  arg.output_data.cell_erad = &cell_erad[0];
  arg.output_data.cell_Trad = &cell_Trad[0];
  arg.output_data.cell_mat_delta_e = &cell_mat_delta_e[0];
  // Check that the pointer matches the data
  FAIL_IF_NOT(rtt_dsxx::soft_equiv(arg.output_data.cell_erad[0], cell_erad[0]));
  FAIL_IF_NOT(rtt_dsxx::soft_equiv(arg.output_data.cell_erad[1], cell_erad[1]));
  FAIL_IF_NOT(rtt_dsxx::soft_equiv(arg.output_data.cell_Trad[0], cell_Trad[0]));
  FAIL_IF_NOT(rtt_dsxx::soft_equiv(arg.output_data.cell_Trad[1], cell_Trad[1]));
  FAIL_IF_NOT(rtt_dsxx::soft_equiv(arg.output_data.cell_mat_delta_e[0], cell_mat_delta_e[0]));
  FAIL_IF_NOT(rtt_dsxx::soft_equiv(arg.output_data.cell_mat_delta_e[1], cell_mat_delta_e[1]));
  FAIL_IF_NOT(rtt_dsxx::soft_equiv(arg.output_data.cell_mat_delta_e[2], cell_mat_delta_e[2]));
  arg.output_data.check_arguments(arg.zonal_data.number_of_local_cells,
                                  arg.zonal_data.number_of_cell_mats);

  odd_solver::Interface_Data iface = odd_api::build_interface_data(arg);

  // FAIL_IF_NOT(iface.mat_data.ipcress_filename == arg.control_data.opacity_file);
  // Check Control data
  FAIL_IF(iface.control_data.multigroup);
  for (auto &rbc : iface.control_data.reflect_bnd)
    FAIL_IF_NOT(rbc);
  for (auto &bT : iface.control_data.bnd_temp)
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(bT, 0.0));

  // Check Mesh Data
  FAIL_IF(iface.mesh_data.domain_decomposed);
  FAIL_IF_NOT(iface.mesh_data.n_dims == 1);
  FAIL_IF_NOT(iface.mesh_data.coord_sys == odd_solver::COORDINATE_SYSTEM::CARTESIAN);
  FAIL_IF_NOT(iface.mesh_data.number_of_local_cells == 2);
  FAIL_IF_NOT(iface.mesh_data.number_of_ghost_cells == 0);
  FAIL_IF_NOT(iface.mesh_data.number_of_global_cells == 2);
  FAIL_IF_NOT(rtt_dsxx::soft_equiv(iface.mesh_data.cell_position.begin(),
                                   iface.mesh_data.cell_position.end(), cell_position.begin(),
                                   cell_position.end()));
  FAIL_IF_NOT(std::equal(iface.mesh_data.cell_global_id.begin(),
                         iface.mesh_data.cell_global_id.end(), cell_global_id.begin(),
                         cell_global_id.end()));
  FAIL_IF_NOT(std::equal(iface.mesh_data.face_types.begin(), iface.mesh_data.face_types.end(),
                         face_type.begin(), face_type.end()));
  FAIL_IF_NOT(std::equal(iface.mesh_data.next_cell_id.begin(), iface.mesh_data.next_cell_id.end(),
                         next_cell_id.begin(), next_cell_id.end()));

  // Check Material Data
  FAIL_IF_NOT(iface.mat_data.ipcress_filename == "not_empty");
  FAIL_IF_NOT(iface.mat_data.number_of_mats == 2);
  FAIL_IF_NOT(std::equal(iface.mat_data.problem_matids.begin(), iface.mat_data.problem_matids.end(),
                         matids.begin(), matids.end()));
  FAIL_IF_NOT(std::equal(iface.mat_data.number_of_cell_mats.begin(),
                         iface.mat_data.number_of_cell_mats.end(), cell_number_of_mats.begin(),
                         cell_number_of_mats.end()));
  FAIL_IF_NOT(std::equal(iface.mat_data.cell_mats.begin(), iface.mat_data.cell_mats.end(),
                         cell_mats_ref.begin(), cell_mats_ref.end()));
  FAIL_IF_NOT(std::equal(iface.mat_data.cell_mat_vol_frac.begin(),
                         iface.mat_data.cell_mat_vol_frac.end(), cell_mat_vol_frac_ref.begin(),
                         cell_mat_vol_frac_ref.end()));
  FAIL_IF_NOT(std::equal(iface.mat_data.cell_mat_temperature.begin(),
                         iface.mat_data.cell_mat_temperature.end(),
                         cell_mat_temperature_ref.begin(), cell_mat_temperature_ref.end()));
  FAIL_IF_NOT(std::equal(iface.mat_data.cell_mat_density.begin(),
                         iface.mat_data.cell_mat_density.end(), cell_mat_density_ref.begin(),
                         cell_mat_density_ref.end()));
  FAIL_IF_NOT(std::equal(iface.mat_data.cell_mat_specific_heat.begin(),
                         iface.mat_data.cell_mat_specific_heat.end(),
                         cell_mat_specific_heat_ref.begin(), cell_mat_specific_heat_ref.end()));
  FAIL_IF_NOT(std::equal(iface.mat_data.cell_velocity.begin(), iface.mat_data.cell_velocity.end(),
                         cell_velocity_ref.begin(), cell_velocity_ref.end()));
  FAIL_IF_NOT(std::equal(iface.mat_data.cell_rad_eden.begin(), iface.mat_data.cell_rad_eden.end(),
                         cell_rad_eden.begin(), cell_rad_eden.end()));

  if (ut.numFails == 0) {
    std::ostringstream m;
    m << "the ODD Interface Builder seem to be working";
    PASSMSG(m.str());
  }
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    // >>> UNIT TESTS
    test(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstArgumetns.cc
//------------------------------------------------------------------------------------------------//
