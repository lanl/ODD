//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   solver/test/tstGrey_Matrix.cc
 * \author Mathew Cleveland
 * \date   December 20th 2021
 * \brief  Testing Gray Matrix class
 * \note   Copyright (C) 2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "Test_Interface_Builder.hh"
#include "solver/Grey_Matrix.hh"
#include "solver/Orthogonal_Mesh.hh"
#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
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
  Orthogonal_Mesh mesh(iface.mesh_data);
  const double dt = 0.1;
  // Test Single Material Matrix
  {
    Test_Single_Mat_Builder(iface);
    Test_Output_Builder(iface);
    iface.mat_data.ipcress_filename = ut.getTestSourcePath() + iface.mat_data.ipcress_filename;
    Grey_Matrix matrix(iface.control_data);
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
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(d, 17.9828, 1e-5));
    for (auto &b : matrix.solver_data.source)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(b, 9.29023, 1e-5));
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
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[0][1], -11.3227, 1e-5));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[1][0], -11.3227, 1e-5));

    // Solve matrix using gauss siedel
    matrix.gs_solver(1.0e-6, 100);
    for (auto &e : matrix.solver_data.cell_eden)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(e, 1.39491, 1e-5));

    // Update the output data
    matrix.calculate_output_data(iface.mat_data, dt, iface.output_data);
    for (auto &mat_de : iface.output_data.cell_mat_dedv)
      for (auto &e : mat_de)
        FAIL_IF_NOT(rtt_dsxx::soft_equiv(e, 1.60509, 1e-5));
    for (auto &e : iface.output_data.cell_rad_eden)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(e, 1.39491, 1e-5));
  }

  // Test Multi Material Matrix
  {
    Test_Multi_Mat_Builder(iface);
    Test_Output_Builder(iface);
    iface.mat_data.ipcress_filename = ut.getTestSourcePath() + iface.mat_data.ipcress_filename;
    Grey_Matrix matrix(iface.control_data);
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
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(d, 12.6124, 1e-5));
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
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[0][1], -7.83406, 1e-5));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[1][0], -7.83406, 1e-5));

    // Solve matrix using gauss siedel
    matrix.gs_solver(1.0e-6, 100);
    for (auto &e : matrix.solver_data.cell_eden)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(e, 4.45166, 1e-5));

    // Update the output data
    matrix.calculate_output_data(iface.mat_data, dt, iface.output_data);
    for (auto &mat_de : iface.output_data.cell_mat_dedv) {
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(mat_de[0], 13.541, 1e-5));
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(mat_de[1], -14.4443, 1e-5));
    }
    //    FAIL_IF_NOT(rtt_dsxx::soft_equiv(e, -4.52159, 1e-5));
    for (auto &e : iface.output_data.cell_rad_eden)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(e, 4.45166, 1e-5));
  }

  if (ut.numFails == 0) {
    std::ostringstream m;
    m << "1D matrix test passed";
    PASSMSG(m.str());
  }
}

void test_2d_matrix(rtt_dsxx::UnitTest &ut) {
  Interface_Data iface;
  Test_2D_Interface_Builder(iface);
  Orthogonal_Mesh mesh(iface.mesh_data);
  const double dt = 0.1;
  // Test Single Material Matrix
  {
    Test_Single_Mat_Builder(iface);
    Test_Output_Builder(iface);
    iface.mat_data.ipcress_filename = ut.getTestSourcePath() + iface.mat_data.ipcress_filename;
    Grey_Matrix matrix(iface.control_data);
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
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(d, 29.3054, 1e-5));
    for (auto &b : matrix.solver_data.source)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(b, 9.29023, 1e-5));
    // Check connectivity vector
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[0][0] == 4);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[0][1] == 1);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[0][2] == 4);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[0][3] == 2);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[1][0] == 0);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[1][1] == 4);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[1][2] == 4);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[1][3] == 3);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[2][0] == 4);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[2][1] == 3);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[2][2] == 0);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[2][3] == 4);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[3][0] == 2);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[3][1] == 4);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[3][2] == 1);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[3][3] == 4);

    // Reflecting boundaries should be zero
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[0][0], 0.0));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[0][2], 0.0));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[1][1], 0.0));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[1][2], 0.0));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[2][0], 0.0));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[2][3], 0.0));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[3][1], 0.0));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[3][3], 0.0));

    // Internal leakage should match left==right==up==down
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[0][1],
                                     matrix.solver_data.off_diagonal[1][0]));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[2][1],
                                     matrix.solver_data.off_diagonal[3][0]));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[0][3],
                                     matrix.solver_data.off_diagonal[2][2]));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[1][3],
                                     matrix.solver_data.off_diagonal[3][2]));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[0][1], -11.3227, 1e-5));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[0][3], -11.3227, 1e-5));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[1][0], -11.3227, 1e-5));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[1][3], -11.3227, 1e-5));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[2][1], -11.3227, 1e-5));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[2][2], -11.3227, 1e-5));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[3][0], -11.3227, 1e-5));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[3][2], -11.3227, 1e-5));

    // Solve matrix using gauss siedel
    matrix.gs_solver(1.0e-6, 100);
    for (auto &e : matrix.solver_data.cell_eden)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(e, 1.39491, 1e-5));

    // Update the output data
    matrix.calculate_output_data(iface.mat_data, dt, iface.output_data);
    for (auto &mat_de : iface.output_data.cell_mat_dedv)
      for (auto &e : mat_de)
        FAIL_IF_NOT(rtt_dsxx::soft_equiv(e, 1.60509, 1e-5));
    for (auto &e : iface.output_data.cell_rad_eden)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(e, 1.39491, 1e-5));
  }

  // Test Multi Material Matrix
  {
    Test_Multi_Mat_Builder(iface);
    Test_Output_Builder(iface);
    iface.mat_data.ipcress_filename = ut.getTestSourcePath() + iface.mat_data.ipcress_filename;
    Grey_Matrix matrix(iface.control_data);
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
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(d, 20.4465, 1e-5));
    for (auto &b : matrix.solver_data.source)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(b, 21.2717, 1e-5));
    // Check connectivity vector
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[0][0] == 4);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[0][1] == 1);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[0][2] == 4);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[0][3] == 2);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[1][0] == 0);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[1][1] == 4);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[1][2] == 4);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[1][3] == 3);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[2][0] == 4);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[2][1] == 3);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[2][2] == 0);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[2][3] == 4);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[3][0] == 2);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[3][1] == 4);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[3][2] == 1);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[3][3] == 4);

    // Reflecting boundaries should be zero
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[0][0], 0.0));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[0][2], 0.0));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[1][1], 0.0));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[1][2], 0.0));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[2][0], 0.0));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[2][3], 0.0));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[3][1], 0.0));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[3][3], 0.0));

    // Internal leakage should match left==right
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[0][1],
                                     matrix.solver_data.off_diagonal[1][0]));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[2][1],
                                     matrix.solver_data.off_diagonal[3][0]));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[0][3],
                                     matrix.solver_data.off_diagonal[2][2]));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[1][3],
                                     matrix.solver_data.off_diagonal[3][2]));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[0][1], -7.83406, 1e-5));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[0][3], -7.83406, 1e-5));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[1][0], -7.83406, 1e-5));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[1][3], -7.83406, 1e-5));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[2][1], -7.83406, 1e-5));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[2][2], -7.83406, 1e-5));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[3][0], -7.83406, 1e-5));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[3][2], -7.83406, 1e-5));

    // Solve matrix using gauss siedel
    matrix.gs_solver(1.0e-6, 100);
    for (auto &e : matrix.solver_data.cell_eden)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(e, 4.45166, 1e-5));

    // Update the output data
    matrix.calculate_output_data(iface.mat_data, dt, iface.output_data);
    for (auto &mat_de : iface.output_data.cell_mat_dedv) {
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(mat_de[0], 13.541, 1e-5));
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(mat_de[1], -14.4443, 1e-5));
    }
    //    FAIL_IF_NOT(rtt_dsxx::soft_equiv(e, -4.52159, 1e-5));
    for (auto &e : iface.output_data.cell_rad_eden)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(e, 4.45166, 1e-5));
  }

  if (ut.numFails == 0) {
    std::ostringstream m;
    m << "2D orthogonal mesh test passed";
    PASSMSG(m.str());
  }
}

void test_3d_matrix(rtt_dsxx::UnitTest &ut) {
  Interface_Data iface;
  Test_3D_Interface_Builder(iface);
  Orthogonal_Mesh mesh(iface.mesh_data);
  const double dt = 0.1;
  // Test Single Material Matrix
  {
    Test_Single_Mat_Builder(iface);
    Test_Output_Builder(iface);
    iface.mat_data.ipcress_filename = ut.getTestSourcePath() + iface.mat_data.ipcress_filename;
    Grey_Matrix matrix(iface.control_data);
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
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(d, 40.6281, 1e-5));
    for (auto &b : matrix.solver_data.source)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(b, 9.29023, 1e-5));
    // Check connectivity vector
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[0][0] == 8);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[0][1] == 1);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[0][2] == 8);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[0][3] == 2);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[0][4] == 8);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[0][5] == 4);

    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[1][0] == 0);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[1][1] == 8);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[1][2] == 8);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[1][3] == 3);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[1][4] == 8);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[1][5] == 5);

    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[2][0] == 8);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[2][1] == 3);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[2][2] == 0);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[2][3] == 8);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[2][4] == 8);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[2][5] == 6);

    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[3][0] == 2);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[3][1] == 8);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[3][2] == 1);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[3][3] == 8);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[3][4] == 8);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[3][5] == 7);

    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[4][0] == 8);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[4][1] == 5);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[4][2] == 8);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[4][3] == 6);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[4][4] == 0);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[4][5] == 8);

    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[5][0] == 4);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[5][1] == 8);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[5][2] == 8);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[5][3] == 7);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[5][4] == 1);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[5][5] == 8);

    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[6][0] == 8);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[6][1] == 7);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[6][2] == 4);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[6][3] == 8);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[6][4] == 2);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[6][5] == 8);

    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[7][0] == 6);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[7][1] == 8);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[7][2] == 5);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[7][3] == 8);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[7][4] == 3);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[7][5] == 8);

    // Reflecting boundaries should be zero
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[0][0], 0.0));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[0][2], 0.0));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[0][4], 0.0));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[1][1], 0.0));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[1][2], 0.0));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[1][4], 0.0));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[2][0], 0.0));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[2][3], 0.0));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[2][4], 0.0));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[3][1], 0.0));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[3][3], 0.0));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[3][4], 0.0));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[4][0], 0.0));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[4][2], 0.0));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[4][5], 0.0));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[5][1], 0.0));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[5][2], 0.0));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[5][5], 0.0));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[6][0], 0.0));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[6][3], 0.0));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[6][5], 0.0));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[7][1], 0.0));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[7][3], 0.0));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[7][5], 0.0));

    // Internal leakage should match left==right==up==down
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[0][1],
                                     matrix.solver_data.off_diagonal[1][0]));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[2][1],
                                     matrix.solver_data.off_diagonal[3][0]));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[0][3],
                                     matrix.solver_data.off_diagonal[2][2]));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[1][3],
                                     matrix.solver_data.off_diagonal[3][2]));

    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[4][1],
                                     matrix.solver_data.off_diagonal[5][0]));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[6][1],
                                     matrix.solver_data.off_diagonal[7][0]));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[4][3],
                                     matrix.solver_data.off_diagonal[6][2]));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[5][3],
                                     matrix.solver_data.off_diagonal[7][2]));

    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[0][5],
                                     matrix.solver_data.off_diagonal[4][4]));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[2][5],
                                     matrix.solver_data.off_diagonal[6][4]));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[0][5],
                                     matrix.solver_data.off_diagonal[6][4]));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[1][5],
                                     matrix.solver_data.off_diagonal[7][4]));

    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[0][1], -11.3227, 1e-5));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[0][3], -11.3227, 1e-5));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[0][5], -11.3227, 1e-5));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[1][0], -11.3227, 1e-5));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[1][3], -11.3227, 1e-5));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[1][5], -11.3227, 1e-5));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[2][1], -11.3227, 1e-5));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[2][2], -11.3227, 1e-5));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[2][5], -11.3227, 1e-5));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[3][0], -11.3227, 1e-5));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[3][2], -11.3227, 1e-5));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[3][5], -11.3227, 1e-5));

    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[4][1], -11.3227, 1e-5));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[4][3], -11.3227, 1e-5));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[4][4], -11.3227, 1e-5));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[5][0], -11.3227, 1e-5));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[5][3], -11.3227, 1e-5));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[5][4], -11.3227, 1e-5));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[6][1], -11.3227, 1e-5));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[6][2], -11.3227, 1e-5));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[6][4], -11.3227, 1e-5));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[7][0], -11.3227, 1e-5));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[7][2], -11.3227, 1e-5));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[7][4], -11.3227, 1e-5));

    // Solve matrix using gauss siedel
    matrix.gs_solver(1.0e-6, 100);
    for (auto &e : matrix.solver_data.cell_eden)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(e, 1.39491, 1e-5));

    // Update the output data
    matrix.calculate_output_data(iface.mat_data, dt, iface.output_data);
    for (auto &mat_de : iface.output_data.cell_mat_dedv)
      for (auto &e : mat_de)
        FAIL_IF_NOT(rtt_dsxx::soft_equiv(e, 1.60509, 1e-5));
    for (auto &e : iface.output_data.cell_rad_eden)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(e, 1.39491, 1e-5));
  }

  // Test Multi Material Matrix
  {
    Test_Multi_Mat_Builder(iface);
    Test_Output_Builder(iface);
    iface.mat_data.ipcress_filename = ut.getTestSourcePath() + iface.mat_data.ipcress_filename;
    Grey_Matrix matrix(iface.control_data);
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
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(d, 28.2806, 1e-5));
    for (auto &b : matrix.solver_data.source)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(b, 21.2717, 1e-5));
    // Check connectivity vector
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[0][0] == 8);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[0][1] == 1);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[0][2] == 8);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[0][3] == 2);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[0][4] == 8);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[0][5] == 4);

    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[1][0] == 0);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[1][1] == 8);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[1][2] == 8);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[1][3] == 3);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[1][4] == 8);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[1][5] == 5);

    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[2][0] == 8);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[2][1] == 3);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[2][2] == 0);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[2][3] == 8);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[2][4] == 8);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[2][5] == 6);

    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[3][0] == 2);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[3][1] == 8);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[3][2] == 1);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[3][3] == 8);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[3][4] == 8);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[3][5] == 7);

    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[4][0] == 8);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[4][1] == 5);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[4][2] == 8);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[4][3] == 6);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[4][4] == 0);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[4][5] == 8);

    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[5][0] == 4);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[5][1] == 8);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[5][2] == 8);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[5][3] == 7);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[5][4] == 1);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[5][5] == 8);

    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[6][0] == 8);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[6][1] == 7);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[6][2] == 4);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[6][3] == 8);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[6][4] == 2);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[6][5] == 8);

    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[7][0] == 6);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[7][1] == 8);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[7][2] == 5);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[7][3] == 8);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[7][4] == 3);
    FAIL_IF_NOT(matrix.solver_data.off_diagonal_id[7][5] == 8);

    // Reflecting boundaries should be zero
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[0][0], 0.0));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[0][2], 0.0));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[0][4], 0.0));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[1][1], 0.0));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[1][2], 0.0));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[1][4], 0.0));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[2][0], 0.0));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[2][3], 0.0));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[2][4], 0.0));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[3][1], 0.0));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[3][3], 0.0));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[3][4], 0.0));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[4][0], 0.0));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[4][2], 0.0));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[4][5], 0.0));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[5][1], 0.0));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[5][2], 0.0));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[5][5], 0.0));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[6][0], 0.0));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[6][3], 0.0));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[6][5], 0.0));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[7][1], 0.0));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[7][3], 0.0));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[7][5], 0.0));

    // Internal leakage should match left==right==up==down
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[0][1],
                                     matrix.solver_data.off_diagonal[1][0]));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[2][1],
                                     matrix.solver_data.off_diagonal[3][0]));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[0][3],
                                     matrix.solver_data.off_diagonal[2][2]));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[1][3],
                                     matrix.solver_data.off_diagonal[3][2]));

    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[4][1],
                                     matrix.solver_data.off_diagonal[5][0]));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[6][1],
                                     matrix.solver_data.off_diagonal[7][0]));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[4][3],
                                     matrix.solver_data.off_diagonal[6][2]));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[5][3],
                                     matrix.solver_data.off_diagonal[7][2]));

    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[0][5],
                                     matrix.solver_data.off_diagonal[4][4]));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[2][5],
                                     matrix.solver_data.off_diagonal[6][4]));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[0][5],
                                     matrix.solver_data.off_diagonal[6][4]));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[1][5],
                                     matrix.solver_data.off_diagonal[7][4]));

    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[0][1], -7.83406, 1e-5));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[0][3], -7.83406, 1e-5));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[0][5], -7.83406, 1e-5));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[1][0], -7.83406, 1e-5));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[1][3], -7.83406, 1e-5));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[1][5], -7.83406, 1e-5));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[2][1], -7.83406, 1e-5));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[2][2], -7.83406, 1e-5));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[2][5], -7.83406, 1e-5));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[3][0], -7.83406, 1e-5));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[3][2], -7.83406, 1e-5));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[3][5], -7.83406, 1e-5));

    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[4][1], -7.83406, 1e-5));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[4][3], -7.83406, 1e-5));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[4][4], -7.83406, 1e-5));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[5][0], -7.83406, 1e-5));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[5][3], -7.83406, 1e-5));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[5][4], -7.83406, 1e-5));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[6][1], -7.83406, 1e-5));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[6][2], -7.83406, 1e-5));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[6][4], -7.83406, 1e-5));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[7][0], -7.83406, 1e-5));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[7][2], -7.83406, 1e-5));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(matrix.solver_data.off_diagonal[7][4], -7.83406, 1e-5));

    // Solve matrix using gauss siedel
    matrix.gs_solver(1.0e-6, 100);
    for (auto &e : matrix.solver_data.cell_eden)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(e, 4.45166, 1e-5));

    // Update the output data
    matrix.calculate_output_data(iface.mat_data, dt, iface.output_data);
    for (auto &mat_de : iface.output_data.cell_mat_dedv) {
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(mat_de[0], 13.541, 1e-5));
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(mat_de[1], -14.4443, 1e-5));
    }
    //    FAIL_IF_NOT(rtt_dsxx::soft_equiv(e, -4.52159, 1e-5));
    for (auto &e : iface.output_data.cell_rad_eden)
      FAIL_IF_NOT(rtt_dsxx::soft_equiv(e, 4.45166, 1e-5));
  }

  if (ut.numFails == 0) {
    std::ostringstream m;
    m << "3D orthogonal mesh test passed";
    PASSMSG(m.str());
  }
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    // >>> UNIT TESTS
    test_1d_matrix(ut);
    test_2d_matrix(ut);
    test_3d_matrix(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstOpacity_Reader.cc
//------------------------------------------------------------------------------------------------//
