//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   driver/Odd_Functions.hh
 * \author Mathew Cleveland
 * \brief  Define Odd Functions for problem setup
 * \note   Copyright (C) 2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#ifndef odd_driver_Odd_Functions_hh
#define odd_driver_Odd_Functions_hh

#include "api/Arguments.hh"
#include "cdi/EoS.hh"
#include <array>
#include <memory>
#include <string>
#include <vector>

namespace odd_driver {

struct Odd_Driver_Data {
  // simple mesh and single material data
  std::string opacity_file;
  size_t n_cycles;
  std::array<double, 3> mesh_size;
  std::array<size_t, 3> mesh_n_cells;
  std::unique_ptr<rtt_cdi::EoS> eos{nullptr};
  size_t matid;
  double temperature;
  double density;
  double specific_heat;
  double specific_heat_Tref{1.0};
  double specific_heat_Tpow{0.0};
  double rad_temperature;
  // Mesh data holder
  std::vector<double> cell_position;
  std::vector<double> cell_size;
  std::vector<size_t> cell_global_id;
  std::vector<size_t> face_type;
  std::vector<size_t> next_cell_id;
  // Material data holder
  std::vector<size_t> problem_matids;
  std::vector<size_t> number_of_cell_mats;
  std::vector<size_t> cell_mats;
  std::vector<double> cell_mat_vol_frac;
  std::vector<double> cell_mat_temperature;
  std::vector<double> cell_mat_density;
  std::vector<double> cell_mat_specific_heat;
  std::vector<double> cell_velocity;
  std::vector<double> cell_erad;
  // Output data
  std::vector<double> output_cell_erad;
  std::vector<double> output_cell_Trad;
  std::vector<double> output_cell_mat_delta_e;
  // ODD Conservation Data
  std::vector<double> mat_energy_density;
  double total_energy;
  double total_rad_energy;
  double total_mat_energy;
  // Boundary conditions
  std::vector<double> bnd_temp{0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  std::vector<size_t> reflect_bnd{1, 1, 1, 1, 1, 1};
};

void build_arguments_from_cmd(const std::vector<std::string> argv, Arguments &args,
                              Odd_Driver_Data &odd_data);

void energy_update(Arguments &args, Odd_Driver_Data &odd_data);

} // end namespace odd_driver

#endif // odd_driver_Odd_Functions_hh

//------------------------------------------------------------------------------------------------//
// end of driver/Odd_Functions.hh
//------------------------------------------------------------------------------------------------//
