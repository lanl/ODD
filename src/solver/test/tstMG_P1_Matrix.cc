//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   solver/test/tstGrey_P1_Matrix.cc
 * \author Mathew Cleveland
 * \date   December 20th 2021
 * \brief  Testing Gray Matrix class
 * \note   Copyright (C) 2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "Test_Interface_Builder.hh"
#include "odd_release/Release.hh"
#include "solver/Constants.hh"
#include "solver/MG_P1_Matrix.hh"
#include "solver/Orthogonal_Mesh.hh"
#include "c4/ParallelUnitTest.hh"
#include "cdi/CDI.hh"
#include "ds++/Release.hh"
#include "ds++/dbc.hh"

using namespace rtt_dsxx;
using namespace odd_solver;
using namespace odd_solver_test;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//
void test_1d_matrix(rtt_dsxx::UnitTest &ut) {
  // setup test interface
  Interface_Data iface;
  Test_1D_Interface_Builder(iface);
  iface.control_data.multigroup = true;
  Orthogonal_Mesh mesh(iface.mesh_data);
  const double dt = 0.1;
  // Test Single Material Matrix
  {
    Test_Single_Mat_Builder(iface);
    Test_Output_Builder(iface);
    iface.mat_data.ipcress_filename = "const_mg_one_two.ipcress";
    iface.mat_data.ipcress_filename = ut.getTestSourcePath() + iface.mat_data.ipcress_filename;
    Opacity_Reader op_reader(iface.mat_data.ipcress_filename);
    size_t ngroups = op_reader.ngroups;
    Test_MG_Data_Builder(iface, ngroups);
    // store the planck PDF for reference T=3.0 [kev]
    std::vector<double> planck_spec(ngroups, 0.0);
    rtt_cdi::CDI::integrate_Planckian_Spectrum(op_reader.group_bounds, 3.0, planck_spec);

    MG_P1_Matrix matrix(iface.control_data);
    matrix.initialize_solver_data(mesh, iface.mat_data, dt);
    // Check initialized/homogenized data
    for (auto &t : matrix.solver_data.cell_temperature0)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(t, 3.0));
    for (auto &t : matrix.solver_data.cell_temperature)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(t, 3.0));
    for (auto &d : matrix.solver_data.cell_density)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(d, 3.0));
    for (auto &mg_e : matrix.solver_data.cell_eden0)
      for (auto &e : mg_e)
        FAIL_IF_NOT(rtt_dsxx::soft_equiv(e, 3.0));
    for (auto &mg_e : matrix.solver_data.cell_eden)
      for (auto &e : mg_e)
        FAIL_IF_NOT(rtt_dsxx::soft_equiv(e, 3.0));
    for (auto &cv : matrix.solver_data.cell_cve)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(cv, 9.0));

    matrix.build_matrix(mesh, dt);
    // Check the matrix values
    for (auto &mg_d : matrix.solver_data.diagonal)
      for (auto &d : mg_d)
        FAIL_IF_NOT(rtt_dsxx::soft_equiv(d, 19.8671, 1e-5));
    // source strength q_g = E0_g + emission*planck_spec[g]
    // emmission = fleck*sigma_p*a*c*T^4 = 3*a*c*3.0**4*dt
    // fleck =  1.0/(1.0 + 4.0*sigma_p*a*c*T^3*dt/cv) = 1.0/(1.0 + 4.0*3.0*a*c*3.0**3*dt/9.0)
    const double fleck = 1.0 / (1.0 + 36.0 * constants::a * constants::c * dt);
    const double emission = fleck * 243.0 * constants::a * constants::c * dt;
    for (auto &mg_b : matrix.solver_data.source)
      for (size_t g = 0; g < ngroups; g++)
        FAIL_IF_NOT(rtt_dsxx::soft_equiv(mg_b[g], 3.0 + emission * planck_spec[g]));

    // Check connectivity vector
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[0][0] == 2);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[0][1] == 1);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[1][0] == 0);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[1][1] == 2);
    // Reflecting boundaries should be zero
    for (auto &od : matrix.solver_data.off_diagonal[0][0])
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(od, 0.0));
    for (auto &od : matrix.solver_data.off_diagonal[1][1])
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(od, 0.0));
    // Internal leakage should match left==right
    for (size_t g = 0; g < matrix.solver_data.off_diagonal[0][1].size(); g++)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[0][1][g],
                                       matrix.solver_data.off_diagonal[1][0][g]));
    for (auto &od : matrix.solver_data.off_diagonal[0][1])
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(od, -13.1776, 1e-5));
    for (auto &od : matrix.solver_data.off_diagonal[1][0])
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(od, -13.1776, 1e-5));

    // Solve matrix using gauss siedel
    matrix.gs_solver(1.0e-12, 100);
    std::vector<double> ref_e{0.448461, 0.448461, 0.448461, 0.448461, 0.448461, 0.448461, 0.448461,
                              0.448461, 0.448461, 0.448461, 0.448461, 0.448461, 0.448461, 0.448461,
                              0.448461, 0.448463, 0.448506, 0.449344, 0.463573, 0.622761, 1.05409,
                              0.59761,  0.448545, 0.448461, 0.448461};
    for (auto &mg_e : matrix.solver_data.cell_eden)
      for (size_t g = 0; g < ngroups; g++)
        FAIL_IF_NOT(rtt_dsxx::soft_equiv(mg_e[g], ref_e[g], 1e-5));
    // Update the output data
    matrix.calculate_output_data(mesh, iface.mat_data, dt, iface.output_data);
    for (auto &mat_de : iface.output_data.cell_mat_dedv)
      for (auto &e : mat_de)
        FAIL_IF_NOT(rtt_dsxx::soft_equiv(e, 62.8433, 1e-5));
    for (auto &e : iface.output_data.cell_rad_eden)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(e, 12.1567, 1e-5));
    for (auto &mg_e : iface.output_data.cell_rad_mg_eden)
      for (size_t g = 0; g < ngroups; g++)
        FAIL_IF_NOT(rtt_dsxx::soft_equiv(mg_e[g], ref_e[g], 1e-5));
    // flux should be zero to within solver tolerance
    for (auto &f : iface.output_data.face_flux) {
      for (auto &v : f)
        FAIL_IF_NOT(rtt_dsxx::soft_equiv(v, 0.0, 1e-5));
    }
    for (auto &mg_f : iface.output_data.face_mg_flux) {
      for (auto &f : mg_f)
        for (auto &v : f)
          FAIL_IF_NOT(rtt_dsxx::soft_equiv(v, 0.0, 1e-5));
    }
  }

  /*
  // Test vacuum boundary condition
  {
    Test_Single_Mat_Builder(iface);
    Test_Output_Builder(iface);
    iface.mat_data.ipcress_filename = ut.getTestSourcePath() + iface.mat_data.ipcress_filename;
    iface.control_data.reflect_bnd[0] = false;
    MG_P1_Matrix matrix(iface.control_data);
    matrix.initialize_solver_data(mesh, iface.mat_data, dt);
    // Check initialized/homogenized data
    for (auto &t : matrix.solver_data.cell_temperature0)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(t, 3.0));
    for (auto &t : matrix.solver_data.cell_temperature)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(t, 3.0));
    for (auto &d : matrix.solver_data.cell_density)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(d, 3.0));
    for (auto &e : matrix.solver_data.cell_eden0)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(e, 3.0));
    for (auto &e : matrix.solver_data.cell_eden)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(e, 3.0));
    for (auto &cv : matrix.solver_data.cell_cve)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(cv, 9.0));

    matrix.build_matrix(mesh, dt);
    // Check the matrix values
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.diagonal[0], 30.7083, 1e-5));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.diagonal[1], 17.8768, 1e-5));
    for (auto &b : matrix.solver_data.source)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(b, 9.29023, 1e-5));
    // Check connectivity vector
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[0][0] == 2);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[0][1] == 1);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[1][0] == 0);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[1][1] == 2);
    // Reflecting boundaries should be zero
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[0][0], 12.8316, 1e-5));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[1][1], 0.0));
    // Internal leakage should match left==right
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[0][1],
                                     matrix.solver_data.off_diagonal[1][0]));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[0][1], -11.2167, 1e-5));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[1][0], -11.2167, 1e-5));

    // Solve matrix using gauss siedel
    matrix.gs_solver(1.0e-6, 100);
    FAIL_IF_NOT(matrix.solver_data.cell_eden[0] < matrix.solver_data.cell_eden[1]);

    // Update the output data
    matrix.calculate_output_data(mesh, iface.mat_data, dt, iface.output_data);
    // Just check for the expected slope
    FAIL_IF_NOT(iface.output_data.cell_rad_eden[0] < iface.output_data.cell_rad_eden[1]);
    FAIL_IF_NOT(iface.output_data.face_flux[0][0] < 0.0);
    FAIL_IF_NOT(iface.output_data.face_flux[0][1] < 0.0);
    FAIL_IF_NOT(iface.output_data.face_flux[1][0] > 0.0);
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(iface.output_data.face_flux[0][1],
                                     -iface.output_data.face_flux[1][0]));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(iface.output_data.face_flux[1][1], 0.0));
    // reset the boundary condition for the rest of the tests
    iface.control_data.reflect_bnd[0] = true;
  }

  // Test vacuum boundary condition on both sides
  {
    Test_Single_Mat_Builder(iface);
    Test_Output_Builder(iface);
    iface.mat_data.ipcress_filename = ut.getTestSourcePath() + iface.mat_data.ipcress_filename;
    iface.control_data.reflect_bnd[0] = false;
    iface.control_data.reflect_bnd[1] = false;
    MG_P1_Matrix matrix(iface.control_data);
    matrix.initialize_solver_data(mesh, iface.mat_data, dt);
    // Check initialized/homogenized data
    for (auto &t : matrix.solver_data.cell_temperature0)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(t, 3.0));
    for (auto &t : matrix.solver_data.cell_temperature)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(t, 3.0));
    for (auto &d : matrix.solver_data.cell_density)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(d, 3.0));
    for (auto &e : matrix.solver_data.cell_eden0)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(e, 3.0));
    for (auto &e : matrix.solver_data.cell_eden)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(e, 3.0));
    for (auto &cv : matrix.solver_data.cell_cve)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(cv, 9.0));

    matrix.build_matrix(mesh, dt);
    // Check the matrix values
    for (auto &d : matrix.solver_data.diagonal)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(d, 30.7083, 1e-5));
    for (auto &b : matrix.solver_data.source)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(b, 9.29023, 1e-5));
    // Check connectivity vector
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[0][0] == 2);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[0][1] == 1);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[1][0] == 0);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[1][1] == 2);
    // Reflecting boundaries should be zero
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[0][0], 12.8316, 1e-5));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[1][1], 12.8316, 1e-5));
    // Internal leakage should match left==right
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[0][1],
                                     matrix.solver_data.off_diagonal[1][0]));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[0][1], -11.2167, 1e-5));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[1][0], -11.2167, 1e-5));

    // Solve matrix using gauss siedel
    matrix.gs_solver(1.0e-6, 100);
    for (auto &e : matrix.solver_data.cell_eden)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(e, 0.476626, 1e-5));

    // Update the output data
    matrix.calculate_output_data(mesh, iface.mat_data, dt, iface.output_data);
    for (auto &mat_de : iface.output_data.cell_mat_dedv)
      for (auto &e : mat_de)
        FAIL_IF_NOT(rtt_dsxx::soft_equiv(e, -3.59249, 1e-5));
    for (auto &e : iface.output_data.cell_rad_eden)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(e, 0.476626, 1e-5));
    // reset the boundary condition for the rest of the tests
    iface.control_data.reflect_bnd[0] = true;
    iface.control_data.reflect_bnd[1] = true;
  }

  // Test boundary sources
  {
    Test_Single_Mat_Builder(iface);
    Test_Output_Builder(iface);
    iface.mat_data.ipcress_filename = ut.getTestSourcePath() + iface.mat_data.ipcress_filename;
    iface.control_data.reflect_bnd[0] = false;
    iface.control_data.reflect_bnd[1] = false;
    iface.control_data.bnd_temp[0] = 5.0;
    iface.control_data.bnd_temp[1] = 5.0;
    MG_P1_Matrix matrix(iface.control_data);
    matrix.initialize_solver_data(mesh, iface.mat_data, dt);
    // Check initialized/homogenized data
    for (auto &t : matrix.solver_data.cell_temperature0)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(t, 3.0));
    for (auto &t : matrix.solver_data.cell_temperature)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(t, 3.0));
    for (auto &d : matrix.solver_data.cell_density)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(d, 3.0));
    for (auto &e : matrix.solver_data.cell_eden0)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(e, 3.0));
    for (auto &e : matrix.solver_data.cell_eden)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(e, 3.0));
    for (auto &cv : matrix.solver_data.cell_cve)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(cv, 9.0));

    matrix.build_matrix(mesh, dt);
    // Check the matrix values
    for (auto &d : matrix.solver_data.diagonal)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(d, 30.7083, 1e-5));
    for (auto &b : matrix.solver_data.source)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(b, 119.322, 1e-5));
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[0][0] == 2);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[0][1] == 1);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[1][0] == 0);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[1][1] == 2);
    // Cheat and store the source energy density in the off diagonal
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[0][0], 12.8316, 1e-5));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[1][1], 12.8316, 1e-5));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.flux_source[0][0], 110.032, 1e-5));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.flux_source[1][1], 110.032, 1e-5));
    // Internal leakage should match left==right
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[0][1],
                                     matrix.solver_data.off_diagonal[1][0]));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[0][1], -11.2167, 1e-5));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[1][0], -11.2167, 1e-5));

    // Solve matrix using gauss siedel
    matrix.gs_solver(1.0e-6, 100);
    for (auto &e : matrix.solver_data.cell_eden)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(e, 6.12172, 1e-5));

    // Update the output data
    matrix.calculate_output_data(mesh, iface.mat_data, dt, iface.output_data);
    for (auto &mat_de : iface.output_data.cell_mat_dedv)
      for (auto &e : mat_de)
        FAIL_IF_NOT(rtt_dsxx::soft_equiv(e, 28.3591, 1e-5));
    for (auto &e : iface.output_data.cell_rad_eden)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(e, 6.12172, 1e-5));
    // reset the boundary condition for the rest of the tests
    iface.control_data.reflect_bnd[0] = true;
    iface.control_data.reflect_bnd[1] = true;
  }

  // Test Multi Material Matrix
  {
    Test_Multi_Mat_Builder(iface);
    Test_Output_Builder(iface);
    iface.mat_data.ipcress_filename = ut.getTestSourcePath() + iface.mat_data.ipcress_filename;
    MG_P1_Matrix matrix(iface.control_data);
    matrix.initialize_solver_data(mesh, iface.mat_data, dt);
    // Check initialized/homogenized data
    for (auto &t : matrix.solver_data.cell_temperature0)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(t, 4.27236, 1.0e-5));
    for (auto &t : matrix.solver_data.cell_temperature)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(t, 4.27236, 1.0e-5));
    for (auto &d : matrix.solver_data.cell_density)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(d, 4.0));
    for (auto &e : matrix.solver_data.cell_eden0)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(e, 4.0));
    for (auto &e : matrix.solver_data.cell_eden)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(e, 4.0));
    for (auto &cv : matrix.solver_data.cell_cve)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(cv, 17.0));

    matrix.build_matrix(mesh, dt);
    // Check the matrix values
    for (auto &d : matrix.solver_data.diagonal)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(d, 12.5616, 1e-5));
    for (auto &b : matrix.solver_data.source)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(b, 21.2717, 1e-5));
    // Check connectivity vector
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[0][0] == 2);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[0][1] == 1);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[1][0] == 0);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[1][1] == 2);
    // Reflecting boundaries should be zero
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[0][0], 0.0));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[1][1], 0.0));
    // Internal leakage should match left==right
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[0][1],
                                     matrix.solver_data.off_diagonal[1][0]));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[0][1], -7.78318, 1e-5));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[1][0], -7.78318, 1e-5));

    // Solve matrix using gauss siedel
    matrix.gs_solver(1.0e-6, 100);
    for (auto &e : matrix.solver_data.cell_eden)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(e, 4.45166, 1e-5));

    // Update the output data
    matrix.calculate_output_data(mesh, iface.mat_data, dt, iface.output_data);
    for (auto &mat_de : iface.output_data.cell_mat_dedv) {
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(mat_de[0], 13.541, 1e-5));
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(mat_de[1], -14.4443, 1e-5));
    }
    //    FAIL_IF_NOT(rtt_dsxx::soft_equiv(e, -4.52159, 1e-5));
    for (auto &e : iface.output_data.cell_rad_eden)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(e, 4.45166, 1e-5));
  }

  // Test Multi Material Matrix with internal and radiation source
  {
    Test_Multi_Mat_Builder(iface);
    Test_Output_Builder(iface);
    iface.mat_data.ipcress_filename = ut.getTestSourcePath() + iface.mat_data.ipcress_filename;
    // ADD sources
    iface.mat_data.cell_rad_source =
        std::vector<double>(iface.mesh_data.number_of_local_cells, 1.0);
    iface.mat_data.cell_mat_electron_source =
        std::vector<std::vector<double>>(iface.mesh_data.number_of_local_cells, {0.5, 1.0});
    MG_P1_Matrix matrix(iface.control_data);
    matrix.initialize_solver_data(mesh, iface.mat_data, dt);
    // Check initialized/homogenized data
    for (auto &t : matrix.solver_data.cell_temperature0)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(t, 4.27236, 1.0e-5));
    for (auto &t : matrix.solver_data.cell_temperature)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(t, 4.27236, 1.0e-5));
    for (auto &d : matrix.solver_data.cell_density)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(d, 4.0));
    for (auto &e : matrix.solver_data.cell_eden0)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(e, 4.0));
    for (auto &e : matrix.solver_data.cell_eden)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(e, 4.0));
    for (auto &cv : matrix.solver_data.cell_cve)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(cv, 17.0));

    matrix.build_matrix(mesh, dt);
    // Check the matrix values
    for (auto &d : matrix.solver_data.diagonal)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(d, 12.5616, 1e-5));
    for (auto &b : matrix.solver_data.source)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(b, 22.9851, 1e-5));
    // Check connectivity vector
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[0][0] == 2);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[0][1] == 1);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[1][0] == 0);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[1][1] == 2);
    // Reflecting boundaries should be zero
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[0][0], 0.0));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[1][1], 0.0));
    // Internal leakage should match left==right
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[0][1],
                                     matrix.solver_data.off_diagonal[1][0]));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[0][1], -7.78318, 1e-5));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[1][0], -7.78318, 1e-5));

    // Solve matrix using gauss siedel
    matrix.gs_solver(1.0e-6, 100);
    for (auto &e : matrix.solver_data.cell_eden)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(e, 4.81024, 1e-5));

    // Update the output data
    matrix.calculate_output_data(mesh, iface.mat_data, dt, iface.output_data);
    for (auto &mat_de : iface.output_data.cell_mat_dedv) {
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(mat_de[0], 14.997, 1e-5));
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(mat_de[1], -13.1834, 1e-5));
    }
    //    FAIL_IF_NOT(rtt_dsxx::soft_equiv(e, -4.52159, 1e-5));
    for (auto &e : iface.output_data.cell_rad_eden)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(e, 4.81024, 1e-5));
  }
*/
  if (ut.numFails == 0) {
    std::ostringstream m;
    m << "1D matrix test passed";
    PASSMSG(m.str());
  }
}

void test_1d_dd_matrix(rtt_dsxx::UnitTest &ut) {
  // setup test interface
  Interface_Data iface;
  bool dd = true;
  Test_1D_Interface_Builder(iface, dd);
  iface.control_data.multigroup = true;
  Orthogonal_Mesh mesh(iface.mesh_data);
  const double dt = 0.1;
  // Test Single Material Matrix
  {
    Test_Single_Mat_Builder(iface);
    Test_Output_Builder(iface);
    iface.mat_data.ipcress_filename = "const_mg_one_two.ipcress";
    iface.mat_data.ipcress_filename = ut.getTestSourcePath() + iface.mat_data.ipcress_filename;
    Opacity_Reader op_reader(iface.mat_data.ipcress_filename);
    size_t ngroups = op_reader.ngroups;
    Test_MG_Data_Builder(iface, ngroups);
    // store the planck PDF for reference T=3.0 [kev]
    std::vector<double> planck_spec(ngroups, 0.0);
    rtt_cdi::CDI::integrate_Planckian_Spectrum(op_reader.group_bounds, 3.0, planck_spec);

    MG_P1_Matrix matrix(iface.control_data);
    matrix.initialize_solver_data(mesh, iface.mat_data, dt);
    // Check initialized/homogenized data
    for (auto &t : matrix.solver_data.cell_temperature0)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(t, 3.0));
    for (auto &t : matrix.solver_data.cell_temperature)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(t, 3.0));
    for (auto &d : matrix.solver_data.cell_density)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(d, 3.0));
    for (auto &mg_e : matrix.solver_data.cell_eden0)
      for (auto &e : mg_e)
        FAIL_IF_NOT(rtt_dsxx::soft_equiv(e, 3.0));
    for (auto &mg_e : matrix.solver_data.cell_eden)
      for (auto &e : mg_e)
        FAIL_IF_NOT(rtt_dsxx::soft_equiv(e, 3.0));
    for (auto &cv : matrix.solver_data.cell_cve)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(cv, 9.0));

    matrix.build_matrix(mesh, dt);
    // Check the matrix values
    for (auto &mg_d : matrix.solver_data.diagonal)
      for (auto &d : mg_d)
        FAIL_IF_NOT(rtt_dsxx::soft_equiv(d, 19.8671, 1e-5));
    // source strength q_g = E0_g + emission*planck_spec[g]
    // emmission = fleck*sigma_p*a*c*T^4 = 3*a*c*3.0**4*dt
    // fleck =  1.0/(1.0 + 4.0*sigma_p*a*c*T^3*dt/cv) = 1.0/(1.0 + 4.0*3.0*a*c*3.0**3*dt/9.0)
    const double fleck = 1.0 / (1.0 + 36.0 * constants::a * constants::c * dt);
    const double emission = fleck * 243.0 * constants::a * constants::c * dt;
    for (auto &mg_b : matrix.solver_data.source)
      for (size_t g = 0; g < ngroups; g++)
        FAIL_IF_NOT(rtt_dsxx::soft_equiv(mg_b[g], 3.0 + emission * planck_spec[g]));

    // Check connectivity vector
    if (rtt_c4::node() == 0) {
      FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[0][0] == 1);
      FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[0][1] == 0);
    }
    if (rtt_c4::node() == 1) {
      FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[0][0] == 0);
      FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[0][1] == 1);
    }

    // Reflecting boundaries should be zero
    if (rtt_c4::node() == 0)
      for (auto &od : matrix.solver_data.off_diagonal[0][0])
        FAIL_IF_NOT(rtt_dsxx::soft_equiv(od, 0.0));
    if (rtt_c4::node() == 1)
      for (auto &od : matrix.solver_data.off_diagonal[0][1])
        FAIL_IF_NOT(rtt_dsxx::soft_equiv(od, 0.0));

    // Internal leakage should match left==right
    if (rtt_c4::node() == 0)
      for (auto &od : matrix.solver_data.off_diagonal[0][1])
        FAIL_IF_NOT(rtt_dsxx::soft_equiv(od, -13.1776, 1e-5));
    if (rtt_c4::node() == 1)
      for (auto &od : matrix.solver_data.off_diagonal[0][0])
        FAIL_IF_NOT(rtt_dsxx::soft_equiv(od, -13.1776, 1e-5));

    // Solve matrix using gauss siedel
    matrix.gs_solver(1.0e-12, 100);
    std::vector<double> ref_e{0.448461, 0.448461, 0.448461, 0.448461, 0.448461, 0.448461, 0.448461,
                              0.448461, 0.448461, 0.448461, 0.448461, 0.448461, 0.448461, 0.448461,
                              0.448461, 0.448463, 0.448506, 0.449344, 0.463573, 0.622761, 1.05409,
                              0.59761,  0.448545, 0.448461, 0.448461};
    for (auto &mg_e : matrix.solver_data.cell_eden)
      for (size_t g = 0; g < ngroups; g++)
        FAIL_IF_NOT(rtt_dsxx::soft_equiv(mg_e[g], ref_e[g], 1e-5));
    // Update the output data
    matrix.calculate_output_data(mesh, iface.mat_data, dt, iface.output_data);
    for (auto &mat_de : iface.output_data.cell_mat_dedv)
      for (auto &e : mat_de)
        FAIL_IF_NOT(rtt_dsxx::soft_equiv(e, 62.8433, 1e-5));
    for (auto &e : iface.output_data.cell_rad_eden)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(e, 12.1567, 1e-5));
    for (auto &mg_e : iface.output_data.cell_rad_mg_eden)
      for (size_t g = 0; g < ngroups; g++)
        FAIL_IF_NOT(rtt_dsxx::soft_equiv(mg_e[g], ref_e[g], 1e-5));
    // flux should be zero to within solver tolerance
    for (auto &f : iface.output_data.face_flux) {
      for (auto &v : f)
        FAIL_IF_NOT(rtt_dsxx::soft_equiv(v, 0.0, 1e-5));
    }
    for (auto &mg_f : iface.output_data.face_mg_flux) {
      for (auto &f : mg_f)
        for (auto &v : f)
          FAIL_IF_NOT(rtt_dsxx::soft_equiv(v, 0.0, 1e-5));
    }
  }
  if (ut.numFails == 0) {
    std::ostringstream m;
    m << "1D DD matrix test passed";
    PASSMSG(m.str());
  }
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_c4::ParallelUnitTest ut(argc, argv, rtt_odd::release);
  try {
    // >>> UNIT TESTS
    test_1d_matrix(ut);
    /*
    test_2d_matrix(ut);
    test_3d_matrix(ut);
    */
    if (rtt_c4::nodes() == 2) {
      test_1d_dd_matrix(ut);
    }
    /*else if (rtt_c4::nodes() == 3) {
      test_2d_dd_matrix(ut);
      test_3d_dd_matrix(ut);
    }
    */
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstOpacity_Reader.cc
//------------------------------------------------------------------------------------------------//
