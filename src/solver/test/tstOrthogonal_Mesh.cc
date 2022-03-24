//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   solver/test/tstOrthogonal_Mesh.cc
 * \author Mathew Cleveland
 * \date   December 20th 2021
 * \brief  Testing Orthogonal Mesh class
 * \note   Copyright (C) 2021-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "Test_Interface_Builder.hh"
#include "odd_release/Release.hh"
#include "solver/Interface_Data.hh"
#include "solver/Orthogonal_Mesh.hh"
#include "c4/ParallelUnitTest.hh"
#include "ds++/Release.hh"
#include "ds++/dbc.hh"

using namespace rtt_dsxx;
using namespace odd_solver;
using namespace odd_solver_test;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//
void test_1d_mesh(rtt_dsxx::UnitTest &ut) {
  Interface_Data iface;
  Test_1D_Interface_Builder(iface);
  Orthogonal_Mesh mesh(iface.mesh_data);
  FAIL_IF_NOT(mesh.number_of_local_cells() == 2);
  FAIL_IF_NOT(mesh.number_of_global_cells() == 2);
  bool domain_replicated = !mesh.domain_decomposed();
  FAIL_IF_NOT(domain_replicated);
  FAIL_IF_NOT(mesh.number_of_ghost_cells() == 0);
  // Check cell data
  {
    // Cell 0
    size_t cell = 0;
    std::array<double, 3> cell_center = mesh.cell_center(cell);
    std::array<double, 3> cell_center_gold{0.25, 0.0, 0.0};
    FAIL_IF_NOT(soft_equiv(cell_center.begin(), cell_center.end(), cell_center_gold.begin(),
                           cell_center_gold.end()));
    FAIL_IF_NOT(soft_equiv(mesh.cell_volume(cell), 0.5));
    FAIL_IF_NOT(mesh.cell_global_id(cell) == 0);
    FAIL_IF_NOT(mesh.number_of_faces(cell) == 2);
    size_t face = 0;
    FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 1.0));
    FAIL_IF_NOT(mesh.face_type(cell, face) == odd_solver::FACE_TYPE::BOUNDARY_FACE);
    face = 1;
    FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 1.0));
    FAIL_IF_NOT(mesh.face_type(cell, face) == odd_solver::FACE_TYPE::INTERNAL_FACE);
    FAIL_IF_NOT(mesh.next_cell(cell, face) == 1);
  }

  // Check cell data
  {
    // Cell 1
    size_t cell = 1;
    std::array<double, 3> cell_center = mesh.cell_center(cell);
    std::array<double, 3> cell_center_gold{0.75, 0.0, 0.0};
    FAIL_IF_NOT(soft_equiv(cell_center.begin(), cell_center.end(), cell_center_gold.begin(),
                           cell_center_gold.end()));
    FAIL_IF_NOT(soft_equiv(mesh.cell_volume(cell), 0.5));
    FAIL_IF_NOT(mesh.cell_global_id(cell) == 1);
    FAIL_IF_NOT(mesh.number_of_faces(cell) == 2);
    size_t face = 0;
    FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 1.0));
    FAIL_IF_NOT(mesh.face_type(cell, face) == odd_solver::FACE_TYPE::INTERNAL_FACE);
    FAIL_IF_NOT(mesh.next_cell(cell, face) == 0);
    face = 1;
    FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 1.0));
    FAIL_IF_NOT(mesh.face_type(cell, face) == odd_solver::FACE_TYPE::BOUNDARY_FACE);
  }

  if (ut.numFails == 0) {
    std::ostringstream m;
    m << "1D orthogonal mesh test passed";
    PASSMSG(m.str());
  }
}

void test_1d_dd_mesh(rtt_dsxx::UnitTest &ut) {
  Interface_Data iface;
  const bool dd = true;
  Test_1D_Interface_Builder(iface, dd);
  Orthogonal_Mesh mesh(iface.mesh_data);
  FAIL_IF_NOT(mesh.number_of_local_cells() == 1);
  FAIL_IF_NOT(mesh.number_of_global_cells() == 2);
  bool domain_decomposed = mesh.domain_decomposed();
  FAIL_IF_NOT(domain_decomposed);
  FAIL_IF_NOT(mesh.number_of_ghost_cells() == 1);
  // Check cell data
  if (rtt_c4::node() == 0) {
    // Cell 0
    size_t cell = 0;
    std::array<double, 3> cell_center = mesh.cell_center(cell);
    std::array<double, 3> cell_center_gold{0.25, 0.0, 0.0};
    FAIL_IF_NOT(soft_equiv(cell_center.begin(), cell_center.end(), cell_center_gold.begin(),
                           cell_center_gold.end()));
    FAIL_IF_NOT(soft_equiv(mesh.cell_volume(cell), 0.5));
    FAIL_IF_NOT(mesh.cell_global_id(cell) == 0);
    FAIL_IF_NOT(mesh.number_of_faces(cell) == 2);
    size_t face = 0;
    FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 1.0));
    FAIL_IF_NOT(mesh.face_type(cell, face) == odd_solver::FACE_TYPE::BOUNDARY_FACE);
    face = 1;
    FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 1.0));
    FAIL_IF_NOT(mesh.face_type(cell, face) == odd_solver::FACE_TYPE::GHOST_FACE);
    FAIL_IF_NOT(mesh.next_cell(cell, face) == 0);
  }

  // Check cell data
  if (rtt_c4::node() == 1) {
    // Cell 1
    size_t cell = 0;
    std::array<double, 3> cell_center = mesh.cell_center(cell);
    std::array<double, 3> cell_center_gold{0.75, 0.0, 0.0};
    FAIL_IF_NOT(soft_equiv(cell_center.begin(), cell_center.end(), cell_center_gold.begin(),
                           cell_center_gold.end()));
    FAIL_IF_NOT(soft_equiv(mesh.cell_volume(cell), 0.5));
    FAIL_IF_NOT(mesh.cell_global_id(cell) == 1);
    FAIL_IF_NOT(mesh.number_of_faces(cell) == 2);
    size_t face = 0;
    FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 1.0));
    FAIL_IF_NOT(mesh.face_type(cell, face) == odd_solver::FACE_TYPE::GHOST_FACE);
    FAIL_IF_NOT(mesh.next_cell(cell, face) == 0);
    face = 1;
    FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 1.0));
    FAIL_IF_NOT(mesh.face_type(cell, face) == odd_solver::FACE_TYPE::BOUNDARY_FACE);
  }

  if (ut.numFails == 0) {
    std::ostringstream m;
    m << "1D DD orthogonal mesh test passed";
    PASSMSG(m.str());
  }
}

void test_2d_mesh(rtt_dsxx::UnitTest &ut) {
  Interface_Data iface;
  Test_2D_Interface_Builder(iface);
  Orthogonal_Mesh mesh(iface.mesh_data);
  FAIL_IF_NOT(mesh.number_of_local_cells() == 4);
  FAIL_IF_NOT(mesh.number_of_global_cells() == 4);
  bool domain_replicated = !mesh.domain_decomposed();
  FAIL_IF_NOT(domain_replicated);
  FAIL_IF_NOT(mesh.number_of_ghost_cells() == 0);
  // Check cell data
  {
    // Cell 0
    size_t cell = 0;
    std::array<double, 3> cell_center = mesh.cell_center(cell);
    std::array<double, 3> cell_center_gold{0.25, 0.25, 0.0};
    FAIL_IF_NOT(soft_equiv(cell_center.begin(), cell_center.end(), cell_center_gold.begin(),
                           cell_center_gold.end()));
    FAIL_IF_NOT(soft_equiv(mesh.cell_volume(cell), 0.25));
    FAIL_IF_NOT(mesh.cell_global_id(cell) == 0);
    FAIL_IF_NOT(mesh.number_of_faces(cell) == 4);
    size_t face = 0;
    FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 0.5));
    FAIL_IF_NOT(mesh.face_type(cell, face) == odd_solver::FACE_TYPE::BOUNDARY_FACE);
    face = 1;
    FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 0.5));
    FAIL_IF_NOT(mesh.face_type(cell, face) == odd_solver::FACE_TYPE::INTERNAL_FACE);
    FAIL_IF_NOT(mesh.next_cell(cell, face) == 1);
    face = 2;
    FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 0.5));
    FAIL_IF_NOT(mesh.face_type(cell, face) == odd_solver::FACE_TYPE::BOUNDARY_FACE);
    face = 3;
    FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 0.5));
    FAIL_IF_NOT(mesh.face_type(cell, face) == odd_solver::FACE_TYPE::INTERNAL_FACE);
    FAIL_IF_NOT(mesh.next_cell(cell, face) == 2);
  }

  {
    // Cell 1
    size_t cell = 1;
    std::array<double, 3> cell_center = mesh.cell_center(cell);
    std::array<double, 3> cell_center_gold{0.75, 0.25, 0.0};
    FAIL_IF_NOT(soft_equiv(cell_center.begin(), cell_center.end(), cell_center_gold.begin(),
                           cell_center_gold.end()));
    FAIL_IF_NOT(soft_equiv(mesh.cell_volume(cell), 0.25));
    FAIL_IF_NOT(mesh.cell_global_id(cell) == 1);
    FAIL_IF_NOT(mesh.number_of_faces(cell) == 4);
    size_t face = 0;
    FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 0.5));
    FAIL_IF_NOT(mesh.face_type(cell, face) == odd_solver::FACE_TYPE::INTERNAL_FACE);
    FAIL_IF_NOT(mesh.next_cell(cell, face) == 0);
    face = 1;
    FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 0.5));
    FAIL_IF_NOT(mesh.face_type(cell, face) == odd_solver::FACE_TYPE::BOUNDARY_FACE);
    face = 2;
    FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 0.5));
    FAIL_IF_NOT(mesh.face_type(cell, face) == odd_solver::FACE_TYPE::BOUNDARY_FACE);
    face = 3;
    FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 0.5));
    FAIL_IF_NOT(mesh.face_type(cell, face) == odd_solver::FACE_TYPE::INTERNAL_FACE);
    FAIL_IF_NOT(mesh.next_cell(cell, face) == 3);
  }

  {
    // Cell 2
    size_t cell = 2;
    std::array<double, 3> cell_center = mesh.cell_center(cell);
    std::array<double, 3> cell_center_gold{0.25, 0.75, 0.0};
    FAIL_IF_NOT(soft_equiv(cell_center.begin(), cell_center.end(), cell_center_gold.begin(),
                           cell_center_gold.end()));
    FAIL_IF_NOT(soft_equiv(mesh.cell_volume(cell), 0.25));
    FAIL_IF_NOT(mesh.cell_global_id(cell) == 2);
    FAIL_IF_NOT(mesh.number_of_faces(cell) == 4);
    size_t face = 0;
    FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 0.5));
    FAIL_IF_NOT(mesh.face_type(cell, face) == odd_solver::FACE_TYPE::BOUNDARY_FACE);
    face = 1;
    FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 0.5));
    FAIL_IF_NOT(mesh.face_type(cell, face) == odd_solver::FACE_TYPE::INTERNAL_FACE);
    FAIL_IF_NOT(mesh.next_cell(cell, face) == 3);
    face = 2;
    FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 0.5));
    FAIL_IF_NOT(mesh.face_type(cell, face) == odd_solver::FACE_TYPE::INTERNAL_FACE);
    FAIL_IF_NOT(mesh.next_cell(cell, face) == 0);
    face = 3;
    FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 0.5));
    FAIL_IF_NOT(mesh.face_type(cell, face) == odd_solver::FACE_TYPE::BOUNDARY_FACE);
  }

  {
    // Cell 3
    size_t cell = 3;
    std::array<double, 3> cell_center = mesh.cell_center(cell);
    std::array<double, 3> cell_center_gold{0.75, 0.75, 0.0};
    FAIL_IF_NOT(soft_equiv(cell_center.begin(), cell_center.end(), cell_center_gold.begin(),
                           cell_center_gold.end()));
    FAIL_IF_NOT(soft_equiv(mesh.cell_volume(cell), 0.25));
    FAIL_IF_NOT(mesh.cell_global_id(cell) == 3);
    FAIL_IF_NOT(mesh.number_of_faces(cell) == 4);
    size_t face = 0;
    FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 0.5));
    FAIL_IF_NOT(mesh.face_type(cell, face) == odd_solver::FACE_TYPE::INTERNAL_FACE);
    FAIL_IF_NOT(mesh.next_cell(cell, face) == 2);
    face = 1;
    FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 0.5));
    FAIL_IF_NOT(mesh.face_type(cell, face) == odd_solver::FACE_TYPE::BOUNDARY_FACE);
    face = 2;
    FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 0.5));
    FAIL_IF_NOT(mesh.face_type(cell, face) == odd_solver::FACE_TYPE::INTERNAL_FACE);
    FAIL_IF_NOT(mesh.next_cell(cell, face) == 1);
    face = 3;
    FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 0.5));
    FAIL_IF_NOT(mesh.face_type(cell, face) == odd_solver::FACE_TYPE::BOUNDARY_FACE);
  }

  if (ut.numFails == 0) {
    std::ostringstream m;
    m << "2D orthogonal mesh test passed";
    PASSMSG(m.str());
  }
}

void test_3d_mesh(rtt_dsxx::UnitTest &ut) {
  Interface_Data iface;
  Test_3D_Interface_Builder(iface);
  Orthogonal_Mesh mesh(iface.mesh_data);
  FAIL_IF_NOT(mesh.number_of_local_cells() == 8);
  FAIL_IF_NOT(mesh.number_of_global_cells() == 8);
  bool domain_replicated = !mesh.domain_decomposed();
  FAIL_IF_NOT(domain_replicated);
  FAIL_IF_NOT(mesh.number_of_ghost_cells() == 0);
  // Check cell data
  {
    // Cell 0
    size_t cell = 0;
    std::array<double, 3> cell_center = mesh.cell_center(cell);
    std::array<double, 3> cell_center_gold{0.25, 0.25, 0.25};
    FAIL_IF_NOT(soft_equiv(cell_center.begin(), cell_center.end(), cell_center_gold.begin(),
                           cell_center_gold.end()));
    FAIL_IF_NOT(soft_equiv(mesh.cell_volume(cell), 0.125));
    FAIL_IF_NOT(mesh.cell_global_id(cell) == 0);
    FAIL_IF_NOT(mesh.number_of_faces(cell) == 6);
    size_t face = 0;
    FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 0.25));
    FAIL_IF_NOT(mesh.face_type(cell, face) == odd_solver::FACE_TYPE::BOUNDARY_FACE);
    face = 1;
    FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 0.25));
    FAIL_IF_NOT(mesh.face_type(cell, face) == odd_solver::FACE_TYPE::INTERNAL_FACE);
    FAIL_IF_NOT(mesh.next_cell(cell, face) == 1);
    face = 2;
    FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 0.25));
    FAIL_IF_NOT(mesh.face_type(cell, face) == odd_solver::FACE_TYPE::BOUNDARY_FACE);
    face = 3;
    FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 0.25));
    FAIL_IF_NOT(mesh.face_type(cell, face) == odd_solver::FACE_TYPE::INTERNAL_FACE);
    FAIL_IF_NOT(mesh.next_cell(cell, face) == 2);
    face = 4;
    FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 0.25));
    FAIL_IF_NOT(mesh.face_type(cell, face) == odd_solver::FACE_TYPE::BOUNDARY_FACE);
    face = 5;
    FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 0.25));
    FAIL_IF_NOT(mesh.face_type(cell, face) == odd_solver::FACE_TYPE::INTERNAL_FACE);
    FAIL_IF_NOT(mesh.next_cell(cell, face) == 4);
  }

  {
    // Cell 1
    size_t cell = 1;
    std::array<double, 3> cell_center = mesh.cell_center(cell);
    std::array<double, 3> cell_center_gold{0.75, 0.25, 0.25};
    FAIL_IF_NOT(soft_equiv(cell_center.begin(), cell_center.end(), cell_center_gold.begin(),
                           cell_center_gold.end()));
    FAIL_IF_NOT(soft_equiv(mesh.cell_volume(cell), 0.125));
    FAIL_IF_NOT(mesh.cell_global_id(cell) == 1);
    FAIL_IF_NOT(mesh.number_of_faces(cell) == 6);
    size_t face = 0;
    FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 0.25));
    FAIL_IF_NOT(mesh.face_type(cell, face) == odd_solver::FACE_TYPE::INTERNAL_FACE);
    FAIL_IF_NOT(mesh.next_cell(cell, face) == 0);
    face = 1;
    FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 0.25));
    FAIL_IF_NOT(mesh.face_type(cell, face) == odd_solver::FACE_TYPE::BOUNDARY_FACE);
    face = 2;
    FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 0.25));
    FAIL_IF_NOT(mesh.face_type(cell, face) == odd_solver::FACE_TYPE::BOUNDARY_FACE);
    face = 3;
    FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 0.25));
    FAIL_IF_NOT(mesh.face_type(cell, face) == odd_solver::FACE_TYPE::INTERNAL_FACE);
    FAIL_IF_NOT(mesh.next_cell(cell, face) == 3);
    face = 4;
    FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 0.25));
    FAIL_IF_NOT(mesh.face_type(cell, face) == odd_solver::FACE_TYPE::BOUNDARY_FACE);
    face = 5;
    FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 0.25));
    FAIL_IF_NOT(mesh.face_type(cell, face) == odd_solver::FACE_TYPE::INTERNAL_FACE);
    FAIL_IF_NOT(mesh.next_cell(cell, face) == 5);
  }

  {
    // Cell 2
    size_t cell = 2;
    std::array<double, 3> cell_center = mesh.cell_center(cell);
    std::array<double, 3> cell_center_gold{0.25, 0.75, 0.25};
    FAIL_IF_NOT(soft_equiv(cell_center.begin(), cell_center.end(), cell_center_gold.begin(),
                           cell_center_gold.end()));
    FAIL_IF_NOT(soft_equiv(mesh.cell_volume(cell), 0.125));
    FAIL_IF_NOT(mesh.cell_global_id(cell) == 2);
    FAIL_IF_NOT(mesh.number_of_faces(cell) == 6);
    size_t face = 0;
    FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 0.25));
    FAIL_IF_NOT(mesh.face_type(cell, face) == odd_solver::FACE_TYPE::BOUNDARY_FACE);
    face = 1;
    FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 0.25));
    FAIL_IF_NOT(mesh.face_type(cell, face) == odd_solver::FACE_TYPE::INTERNAL_FACE);
    FAIL_IF_NOT(mesh.next_cell(cell, face) == 3);
    face = 2;
    FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 0.25));
    FAIL_IF_NOT(mesh.face_type(cell, face) == odd_solver::FACE_TYPE::INTERNAL_FACE);
    FAIL_IF_NOT(mesh.next_cell(cell, face) == 0);
    face = 3;
    FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 0.25));
    FAIL_IF_NOT(mesh.face_type(cell, face) == odd_solver::FACE_TYPE::BOUNDARY_FACE);
    face = 4;
    FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 0.25));
    FAIL_IF_NOT(mesh.face_type(cell, face) == odd_solver::FACE_TYPE::BOUNDARY_FACE);
    face = 5;
    FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 0.25));
    FAIL_IF_NOT(mesh.face_type(cell, face) == odd_solver::FACE_TYPE::INTERNAL_FACE);
    FAIL_IF_NOT(mesh.next_cell(cell, face) == 6);
  }

  {
    // Cell 3
    size_t cell = 3;
    std::array<double, 3> cell_center = mesh.cell_center(cell);
    std::array<double, 3> cell_center_gold{0.75, 0.75, 0.25};
    FAIL_IF_NOT(soft_equiv(cell_center.begin(), cell_center.end(), cell_center_gold.begin(),
                           cell_center_gold.end()));
    FAIL_IF_NOT(soft_equiv(mesh.cell_volume(cell), 0.125));
    FAIL_IF_NOT(mesh.cell_global_id(cell) == 3);
    FAIL_IF_NOT(mesh.number_of_faces(cell) == 6);
    size_t face = 0;
    FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 0.25));
    FAIL_IF_NOT(mesh.face_type(cell, face) == odd_solver::FACE_TYPE::INTERNAL_FACE);
    FAIL_IF_NOT(mesh.next_cell(cell, face) == 2);
    face = 1;
    FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 0.25));
    FAIL_IF_NOT(mesh.face_type(cell, face) == odd_solver::FACE_TYPE::BOUNDARY_FACE);
    face = 2;
    FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 0.25));
    FAIL_IF_NOT(mesh.face_type(cell, face) == odd_solver::FACE_TYPE::INTERNAL_FACE);
    FAIL_IF_NOT(mesh.next_cell(cell, face) == 1);
    face = 3;
    FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 0.25));
    FAIL_IF_NOT(mesh.face_type(cell, face) == odd_solver::FACE_TYPE::BOUNDARY_FACE);
    face = 4;
    FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 0.25));
    FAIL_IF_NOT(mesh.face_type(cell, face) == odd_solver::FACE_TYPE::BOUNDARY_FACE);
    face = 5;
    FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 0.25));
    FAIL_IF_NOT(mesh.face_type(cell, face) == odd_solver::FACE_TYPE::INTERNAL_FACE);
    FAIL_IF_NOT(mesh.next_cell(cell, face) == 7);
  }

  {
    // Cell 4
    size_t cell = 4;
    std::array<double, 3> cell_center = mesh.cell_center(cell);
    std::array<double, 3> cell_center_gold{0.25, 0.25, 0.75};
    FAIL_IF_NOT(soft_equiv(cell_center.begin(), cell_center.end(), cell_center_gold.begin(),
                           cell_center_gold.end()));
    FAIL_IF_NOT(soft_equiv(mesh.cell_volume(cell), 0.125));
    FAIL_IF_NOT(mesh.cell_global_id(cell) == 4);
    FAIL_IF_NOT(mesh.number_of_faces(cell) == 6);
    size_t face = 0;
    FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 0.25));
    FAIL_IF_NOT(mesh.face_type(cell, face) == odd_solver::FACE_TYPE::BOUNDARY_FACE);
    face = 1;
    FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 0.25));
    FAIL_IF_NOT(mesh.face_type(cell, face) == odd_solver::FACE_TYPE::INTERNAL_FACE);
    FAIL_IF_NOT(mesh.next_cell(cell, face) == 5);
    face = 2;
    FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 0.25));
    FAIL_IF_NOT(mesh.face_type(cell, face) == odd_solver::FACE_TYPE::BOUNDARY_FACE);
    face = 3;
    FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 0.25));
    FAIL_IF_NOT(mesh.face_type(cell, face) == odd_solver::FACE_TYPE::INTERNAL_FACE);
    FAIL_IF_NOT(mesh.next_cell(cell, face) == 6);
    face = 4;
    FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 0.25));
    FAIL_IF_NOT(mesh.face_type(cell, face) == odd_solver::FACE_TYPE::INTERNAL_FACE);
    FAIL_IF_NOT(mesh.next_cell(cell, face) == 0);
    face = 5;
    FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 0.25));
    FAIL_IF_NOT(mesh.face_type(cell, face) == odd_solver::FACE_TYPE::BOUNDARY_FACE);
  }

  {
    // Cell 5
    size_t cell = 5;
    std::array<double, 3> cell_center = mesh.cell_center(cell);
    std::array<double, 3> cell_center_gold{0.75, 0.25, 0.75};
    FAIL_IF_NOT(soft_equiv(cell_center.begin(), cell_center.end(), cell_center_gold.begin(),
                           cell_center_gold.end()));
    FAIL_IF_NOT(soft_equiv(mesh.cell_volume(cell), 0.125));
    FAIL_IF_NOT(mesh.cell_global_id(cell) == 5);
    FAIL_IF_NOT(mesh.number_of_faces(cell) == 6);
    size_t face = 0;
    FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 0.25));
    FAIL_IF_NOT(mesh.face_type(cell, face) == odd_solver::FACE_TYPE::INTERNAL_FACE);
    FAIL_IF_NOT(mesh.next_cell(cell, face) == 4);
    face = 1;
    FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 0.25));
    FAIL_IF_NOT(mesh.face_type(cell, face) == odd_solver::FACE_TYPE::BOUNDARY_FACE);
    face = 2;
    FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 0.25));
    FAIL_IF_NOT(mesh.face_type(cell, face) == odd_solver::FACE_TYPE::BOUNDARY_FACE);
    face = 3;
    FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 0.25));
    FAIL_IF_NOT(mesh.face_type(cell, face) == odd_solver::FACE_TYPE::INTERNAL_FACE);
    FAIL_IF_NOT(mesh.next_cell(cell, face) == 7);
    face = 4;
    FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 0.25));
    FAIL_IF_NOT(mesh.face_type(cell, face) == odd_solver::FACE_TYPE::INTERNAL_FACE);
    FAIL_IF_NOT(mesh.next_cell(cell, face) == 1);
    face = 5;
    FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 0.25));
    FAIL_IF_NOT(mesh.face_type(cell, face) == odd_solver::FACE_TYPE::BOUNDARY_FACE);
  }

  {
    // Cell 6
    size_t cell = 6;
    std::array<double, 3> cell_center = mesh.cell_center(cell);
    std::array<double, 3> cell_center_gold{0.25, 0.75, 0.75};
    FAIL_IF_NOT(soft_equiv(cell_center.begin(), cell_center.end(), cell_center_gold.begin(),
                           cell_center_gold.end()));
    FAIL_IF_NOT(soft_equiv(mesh.cell_volume(cell), 0.125));
    FAIL_IF_NOT(mesh.cell_global_id(cell) == 6);
    FAIL_IF_NOT(mesh.number_of_faces(cell) == 6);
    size_t face = 0;
    FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 0.25));
    FAIL_IF_NOT(mesh.face_type(cell, face) == odd_solver::FACE_TYPE::BOUNDARY_FACE);
    face = 1;
    FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 0.25));
    FAIL_IF_NOT(mesh.face_type(cell, face) == odd_solver::FACE_TYPE::INTERNAL_FACE);
    FAIL_IF_NOT(mesh.next_cell(cell, face) == 7);
    face = 2;
    FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 0.25));
    FAIL_IF_NOT(mesh.face_type(cell, face) == odd_solver::FACE_TYPE::INTERNAL_FACE);
    FAIL_IF_NOT(mesh.next_cell(cell, face) == 4);
    face = 3;
    FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 0.25));
    FAIL_IF_NOT(mesh.face_type(cell, face) == odd_solver::FACE_TYPE::BOUNDARY_FACE);
    face = 4;
    FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 0.25));
    FAIL_IF_NOT(mesh.face_type(cell, face) == odd_solver::FACE_TYPE::INTERNAL_FACE);
    FAIL_IF_NOT(mesh.next_cell(cell, face) == 2);
    face = 5;
    FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 0.25));
    FAIL_IF_NOT(mesh.face_type(cell, face) == odd_solver::FACE_TYPE::BOUNDARY_FACE);
  }

  {
    // Cell 7
    size_t cell = 7;
    std::array<double, 3> cell_center = mesh.cell_center(cell);
    std::array<double, 3> cell_center_gold{0.75, 0.75, 0.75};
    FAIL_IF_NOT(soft_equiv(cell_center.begin(), cell_center.end(), cell_center_gold.begin(),
                           cell_center_gold.end()));
    FAIL_IF_NOT(soft_equiv(mesh.cell_volume(cell), 0.125));
    FAIL_IF_NOT(mesh.cell_global_id(cell) == 7);
    FAIL_IF_NOT(mesh.number_of_faces(cell) == 6);
    size_t face = 0;
    FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 0.25));
    FAIL_IF_NOT(mesh.face_type(cell, face) == odd_solver::FACE_TYPE::INTERNAL_FACE);
    FAIL_IF_NOT(mesh.next_cell(cell, face) == 6);
    face = 1;
    FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 0.25));
    FAIL_IF_NOT(mesh.face_type(cell, face) == odd_solver::FACE_TYPE::BOUNDARY_FACE);
    face = 2;
    FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 0.25));
    FAIL_IF_NOT(mesh.face_type(cell, face) == odd_solver::FACE_TYPE::INTERNAL_FACE);
    FAIL_IF_NOT(mesh.next_cell(cell, face) == 5);
    face = 3;
    FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 0.25));
    FAIL_IF_NOT(mesh.face_type(cell, face) == odd_solver::FACE_TYPE::BOUNDARY_FACE);
    face = 4;
    FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 0.25));
    FAIL_IF_NOT(mesh.face_type(cell, face) == odd_solver::FACE_TYPE::INTERNAL_FACE);
    FAIL_IF_NOT(mesh.next_cell(cell, face) == 3);
    face = 5;
    FAIL_IF_NOT(soft_equiv(mesh.face_area(cell, face), 0.25));
    FAIL_IF_NOT(mesh.face_type(cell, face) == odd_solver::FACE_TYPE::BOUNDARY_FACE);
  }

  if (ut.numFails == 0) {
    std::ostringstream m;
    m << "3D orthogonal mesh test passed";
    PASSMSG(m.str());
  }
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_c4::ParallelUnitTest ut(argc, argv, rtt_odd::release);
  try {
    // >>> Replicated UNIT TESTS
    test_1d_mesh(ut);
    test_2d_mesh(ut);
    test_3d_mesh(ut);
    if (rtt_c4::nodes() == 2) {
      test_1d_dd_mesh(ut);
    }
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstOpacity_Reader.cc
//------------------------------------------------------------------------------------------------//
