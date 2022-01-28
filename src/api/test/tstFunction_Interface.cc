//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   api/test/tstArguments.cc
 * \author Mathew Cleveland
 * \date   October 21st 2021
 * \brief  Testing opacity reader
 * \note   Copyright (C) 2021-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "api/Arguments.hh"
#include "api/Function_Interface.hh"
#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"

using namespace rtt_dsxx;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//
void test(rtt_dsxx::UnitTest &ut) {
  Arguments arg;

  // fill control data with "valid junk" and check arguments
  std::string file = "two-mats.ipcress";
  arg.control_data.opacity_file = &file[0];
  arg.control_data.dt = 0.1;
  arg.control_data.max_iter = 100;
  arg.control_data.min_tol = 1.0e-12;
  arg.control_data.check_arguments();

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
  std::vector<size_t> matids = {10001, 10002};
  arg.zonal_data.number_of_mats = 2;
  arg.zonal_data.problem_matids = &matids[0];
  // cell wise material data
  // 1 material in cell 1 and 2 materials in cell 2
  std::vector<size_t> cell_number_of_mats{1, 2};
  arg.zonal_data.number_of_cell_mats = &cell_number_of_mats[0];
  // cell 1 (mat 1) cell 2 (mat 0 and mat 1)
  std::vector<size_t> cell_mats{0, 0, 1};
  // cell 1 (1.0) cell 2 mat 0 (0.5) and mat 1 (0.5)
  std::vector<double> cell_mat_vol_frac{1.0, 0.5, 0.5};
  // cell 1 (1.0) cell 2 mat 0 (10.0) and mat 1 (1.0)
  std::vector<double> cell_mat_temperature{1.0, 10.0, 1.0};
  // cell 1 (10) cell 2 mat 0 (1.0) and mat 1 (10.0)
  std::vector<double> cell_mat_density{10.0, 1.0, 10.0};
  // cell 1 (0.1) cell 2 mat 0 (0.1) and mat 1 (0.01)
  std::vector<double> cell_mat_specific_heat{0.1, 0.1, 0.01};
  // cell velocity
  std::vector<double> cell_velocity{0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  // cell radiation energy density
  std::vector<double> cell_erad0{32.0, 23.0};
  arg.zonal_data.cell_mats = &cell_mats[0];
  arg.zonal_data.cell_mat_vol_frac = &cell_mat_vol_frac[0];
  arg.zonal_data.cell_mat_temperature = &cell_mat_temperature[0];
  arg.zonal_data.cell_mat_density = &cell_mat_density[0];
  arg.zonal_data.cell_mat_specific_heat = &cell_mat_specific_heat[0];
  arg.zonal_data.cell_velocity = &cell_velocity[0];
  arg.zonal_data.cell_erad = &cell_erad0[0];

  // setup the output data
  std::vector<double> cell_erad(arg.zonal_data.number_of_local_cells, 0.0);
  std::vector<double> cell_Trad(arg.zonal_data.number_of_local_cells, 0.0);
  std::vector<double> cell_mat_delta_e(cell_mats.size(), 0.0);
  arg.output_data.cell_erad = &cell_erad[0];
  arg.output_data.cell_Trad = &cell_Trad[0];
  arg.output_data.cell_mat_delta_e = &cell_mat_delta_e[0];

  // Call the solver on the fake arguments list
  std::cout << "Call solver" << std::endl;
  Odd_Diffusion_Solve(arg);

  // print out the opacity data
  size_t mat_index = 0;
  for (size_t i = 0; i < arg.zonal_data.number_of_local_cells; i++) {
    std::cout << "erad[" << i << "] = " << arg.output_data.cell_erad[i] << std::endl;
    std::cout << "Trad[" << i << "] = " << arg.output_data.cell_Trad[i] << std::endl;
    for (size_t m = 0; m < cell_number_of_mats[i]; m++, mat_index++)
      std::cout << "cell_mat_delta_e[" << i << "][" << m
                << "] = " << arg.output_data.cell_mat_delta_e[mat_index] << std::endl;
  }

  if (ut.numFails == 0) {
    std::ostringstream m;
    m << "the ODD API seems to be working";
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
// end of tstFunction_Interface.cc
//------------------------------------------------------------------------------------------------//
