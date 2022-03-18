#include "../Arguments.h"

int main(int argc, char *argv[]) {

  if (argc > 1)
    printf("%s", argv[1]);

  struct Arguments arg;

  arg.control_data.opacity_file = "two-mats.ipcress";
  arg.control_data.dt = 0.1;
  arg.control_data.max_iter = 20;
  arg.control_data.min_tol = 1.0e-12;
  arg.control_data.print = 1;
  double bnd_temp[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  arg.control_data.bnd_temp = &bnd_temp[0];
  size_t reflect_bnd[6] = {1, 1, 1, 1, 1, 1};
  arg.control_data.reflect_bnd = &reflect_bnd[0];

  // fill the zonal data with "valid junk" and check arguments
  arg.zonal_data.domain_decomposed = 0;
  arg.zonal_data.number_of_local_cells = 2;
  arg.zonal_data.number_of_global_cells = 2;
  arg.zonal_data.dimensions = 1;
  arg.zonal_data.coord_sys = 0;

  // 2 zones | 0 || 1 | with dx=0.5 dy=0 and dz=0
  double cell_position[6] = {0.25, 0.0, 0.0, 0.75, 0.0, 0.0};
  arg.zonal_data.cell_position = &cell_position[0];
  double cell_size[6] = {0.5, 0.0, 0.0, 0.5, 0.0, 0.0};
  arg.zonal_data.cell_size = &cell_size[0];
  size_t cell_global_id[2] = {0, 1};
  arg.zonal_data.cell_global_id = &cell_global_id[0];
  size_t face_type[4] = {1, 0, 0, 1};
  arg.zonal_data.face_type = &face_type[0];
  size_t next_cell_id[4] = {2, 1, 0, 2};
  arg.zonal_data.next_cell_id = &next_cell_id[0];

  // global material data
  size_t matids[2] = {10001, 10002};
  arg.zonal_data.number_of_mats = 2;
  arg.zonal_data.problem_matids = &matids[0];
  // cell wise material data
  // 1 material in cell 1 and 2 materials in cell 2
  size_t cell_number_of_mats[2] = {1, 2};
  arg.zonal_data.number_of_cell_mats = &cell_number_of_mats[0];
  // cell 1 (mat 1) cell 2 (mat 1 and mat 2)
  size_t cell_mats[3] = {0, 0, 1};
  // cell 1 (1.0) cell 2 mat 1 (0.5) and mat 2 (0.5)
  double cell_mat_vol_frac[3] = {1.0, 0.5, 0.5};
  // cell 1 (1.0) cell 2 mat 1 (10.0) and mat 2 (1.0)
  double cell_mat_temperature[3] = {1.0, 10.0, 1.0};
  // cell 1 (10) cell 2 mat 1 (1.0) and mat 2 (10.0)
  double cell_mat_density[3] = {10.0, 1.0, 10.0};
  // cell 1 (0.1) cell 2 mat 0 (0.1) and mat 1 (0.01)
  double cell_mat_specific_heat[3] = {0.1, 0.1, 0.01};
  // cell velocity
  double cell_velocity[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  // cell radiation energy density
  double cell_rad_eden[2] = {32.0, 23.0};
  arg.zonal_data.cell_mats = &cell_mats[0];
  arg.zonal_data.cell_mat_vol_frac = &cell_mat_vol_frac[0];
  arg.zonal_data.cell_mat_temperature = &cell_mat_temperature[0];
  arg.zonal_data.cell_mat_density = &cell_mat_density[0];
  arg.zonal_data.cell_mat_specific_heat = &cell_mat_specific_heat[0];
  arg.zonal_data.cell_velocity = &cell_velocity[0];
  arg.zonal_data.cell_erad = &cell_rad_eden[0];

  // Output data
  double cell_erad[2] = {17000.0, 17000.1};
  double cell_Trad[2] = {18000.0, 18000.1};
  double cell_mat_delta_e[3] = {-3.0, 4.0, -5.0};
  arg.output_data.cell_erad = &cell_erad[0];
  arg.output_data.cell_Trad = &cell_Trad[0];
  arg.output_data.cell_mat_delta_e = &cell_mat_delta_e[0];

  Odd_Diffusion_Solve(&arg);

  // print out the opacity data
  size_t mat_index = 0;
  for (size_t i = 0; i < arg.zonal_data.number_of_local_cells; i++) {
    printf("erad[%lu] = %e\n", i, arg.output_data.cell_erad[i]);
    printf("Trad[%lu] = %e\n", i, arg.output_data.cell_Trad[i]);
    for (size_t m = 0; m < cell_number_of_mats[i]; m++, mat_index++)
      printf("cell_mat_delta_e[%lu][%lu] = %e\n", i, m,
             arg.output_data.cell_mat_delta_e[mat_index]);
  }
}
