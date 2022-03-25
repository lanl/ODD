//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   solver/test/tstGhost_Comm.cc
 * \author Mathew Cleveland
 * \date   December 20th 2021
 * \brief  Testing Orthogonal Mesh class
 * \note   Copyright (C) 2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "Test_Interface_Builder.hh"
#include "odd_release/Release.hh"
#include "solver/Ghost_Comm.hh"
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
void test_1d_dd_comm(rtt_dsxx::UnitTest &ut) {
  Interface_Data iface;
  const bool dd = true;
  Test_1D_Interface_Builder(iface, dd);
  Orthogonal_Mesh mesh(iface.mesh_data);
  FAIL_IF_NOT(mesh.number_of_local_cells() == 1);
  FAIL_IF_NOT(mesh.number_of_global_cells() == 2);
  bool domain_decomposed = mesh.domain_decomposed();
  FAIL_IF_NOT(domain_decomposed);
  FAIL_IF_NOT(mesh.number_of_ghost_cells() == 1);
  Ghost_Comm gcomm(mesh);
  for (auto &map : gcomm.put_buffer_size)
    std::cout << rtt_c4::node() << " BUFFER SIZE put_rank=" << map.first
              << " put_size=" << map.second << std::endl;
  for (auto &map : gcomm.put_map) {
    std::cout << rtt_c4::node() << " PUT MAP put_cell=" << map.first;
    for (auto &face_map : map.second)
      std::cout << " face_id=" << face_map.first << " rank_id=" << face_map.second.first
                << " put_buffer_id=" << face_map.second.second;
    std::cout << std::endl;
  }
  for (auto &map : gcomm.ghost_map) {
    std::cout << rtt_c4::node() << " GHOST MAP local_ghost_id=" << map.first;
    for (auto &face_map : map.second)
      std::cout << " face_id=" << face_map.first << " ghost_buffer_index=" << face_map.second;
    std::cout << std::endl;
  }
  std::vector<double> local_ghost_data(gcomm.local_ghost_buffer_size);
  std::map<size_t, std::vector<double>> local_put_data;
  for (auto &map : gcomm.put_buffer_size)
    local_put_data[map.first] = std::vector<double>(map.second, rtt_c4::node() + 3);

  gcomm.exchange_ghost_data(local_put_data, local_ghost_data);
  // Check cell data
  if (rtt_c4::node() == 0) {
    FAIL_IF_NOT(gcomm.local_ghost_buffer_size == 1);
    FAIL_IF_NOT(gcomm.put_buffer_size.size() == 1);
    FAIL_IF_NOT(gcomm.put_buffer_size[1] == 1);
    FAIL_IF_NOT(gcomm.put_map.size() == 1);
    FAIL_IF_NOT(gcomm.put_map[0].size() == 1);
    FAIL_IF_NOT(gcomm.put_map[0][1].first == 1);
    FAIL_IF_NOT(gcomm.put_map[0][1].second == 0);
    FAIL_IF_NOT(gcomm.ghost_map.size() == 1);
    FAIL_IF_NOT(gcomm.ghost_map[0].size() == 1);
    FAIL_IF_NOT(gcomm.ghost_map[0][1] == 0);
    FAIL_IF_NOT(local_ghost_data[0] == 4);
  }

  // Check cell data
  if (rtt_c4::node() == 1) {
    FAIL_IF_NOT(gcomm.local_ghost_buffer_size == 1);
    FAIL_IF_NOT(gcomm.put_buffer_size.size() == 1);
    FAIL_IF_NOT(gcomm.put_buffer_size[0] == 1);
    FAIL_IF_NOT(gcomm.put_map.size() == 1);
    FAIL_IF_NOT(gcomm.put_map[0][0].first == 0);
    FAIL_IF_NOT(gcomm.put_map[0][0].second == 0);
    FAIL_IF_NOT(gcomm.ghost_map.size() == 1);
    FAIL_IF_NOT(gcomm.ghost_map[0].size() == 1);
    FAIL_IF_NOT(gcomm.ghost_map[0][0] == 0);
    FAIL_IF_NOT(local_ghost_data[0] == 3);
  }

  if (ut.numFails == 0) {
    std::ostringstream m;
    m << "1D DD orthogonal mesh test passed";
    PASSMSG(m.str());
  }
}

//------------------------------------------------------------------------------------------------//
void test_2d_dd_comm(rtt_dsxx::UnitTest &ut) {
  Interface_Data iface;
  const bool dd = true;
  Test_2D_Interface_Builder(iface, dd);
  Orthogonal_Mesh mesh(iface.mesh_data);
  FAIL_IF_NOT(mesh.number_of_global_cells() == 4);
  bool domain_decomposed = mesh.domain_decomposed();
  FAIL_IF_NOT(domain_decomposed);
  Ghost_Comm gcomm(mesh);
  for (auto &map : gcomm.put_buffer_size)
    std::cout << rtt_c4::node() << " BUFFER SIZE put_rank=" << map.first
              << " put_size=" << map.second << std::endl;
  for (auto &map : gcomm.put_map) {
    std::cout << rtt_c4::node() << " PUT MAP put_cell=" << map.first;
    for (auto &face_map : map.second)
      std::cout << " face_id=" << face_map.first << " rank_id=" << face_map.second.first
                << " put_buffer_id=" << face_map.second.second;
    std::cout << std::endl;
  }
  for (auto &map : gcomm.ghost_map) {
    std::cout << rtt_c4::node() << " GHOST MAP local_ghost_id=" << map.first;
    for (auto &face_map : map.second)
      std::cout << " face_id=" << face_map.first << " ghost_buffer_index=" << face_map.second;
    std::cout << std::endl;
  }
  std::vector<double> local_ghost_data(gcomm.local_ghost_buffer_size);
  std::map<size_t, std::vector<double>> local_put_data;
  for (auto &map : gcomm.put_buffer_size)
    local_put_data[map.first] = std::vector<double>(map.second, rtt_c4::node() + 3);

  gcomm.exchange_ghost_data(local_put_data, local_ghost_data);
  // Check cell data
  if (rtt_c4::node() == 0) {
    FAIL_IF_NOT(mesh.number_of_local_cells() == 1);
    FAIL_IF_NOT(mesh.number_of_ghost_cells() == 2);
    FAIL_IF_NOT(gcomm.local_ghost_buffer_size == 2);
    FAIL_IF_NOT(gcomm.put_buffer_size.size() == 2);
    FAIL_IF_NOT(gcomm.put_buffer_size[1] == 1);
    FAIL_IF_NOT(gcomm.put_buffer_size[2] == 1);
    FAIL_IF_NOT(gcomm.put_map.size() == 1);
    // put cell 0 face 1
    FAIL_IF_NOT(gcomm.put_map[0].size() == 2);
    FAIL_IF_NOT(gcomm.put_map[0][1].first == 1);
    FAIL_IF_NOT(gcomm.put_map[0][1].second == 0);
    // put cell 0 face 3
    FAIL_IF_NOT(gcomm.put_map[0][3].first == 2);
    FAIL_IF_NOT(gcomm.put_map[0][3].second == 0);
    // ghost cell 0 (1 global_id) face 0
    FAIL_IF_NOT(gcomm.ghost_map.size() == 2);
    FAIL_IF_NOT(gcomm.ghost_map[0].size() == 1);
    FAIL_IF_NOT(gcomm.ghost_map[0][0] == 0);
    // ghost cell 1 (2 global_id) face 2
    FAIL_IF_NOT(gcomm.ghost_map[1].size() == 1);
    FAIL_IF_NOT(gcomm.ghost_map[1][2] == 1);

    FAIL_IF_NOT(local_ghost_data[0] == 4);
    FAIL_IF_NOT(local_ghost_data[1] == 5);
  }

  // Check cell data
  if (rtt_c4::node() == 1) {
    FAIL_IF_NOT(mesh.number_of_local_cells() == 1);
    FAIL_IF_NOT(mesh.number_of_ghost_cells() == 2);
    FAIL_IF_NOT(gcomm.local_ghost_buffer_size == 2);
    FAIL_IF_NOT(gcomm.put_buffer_size.size() == 2);
    FAIL_IF_NOT(gcomm.put_buffer_size[0] == 1);
    FAIL_IF_NOT(gcomm.put_buffer_size[2] == 1);
    // put cell 1 face 0
    FAIL_IF_NOT(gcomm.put_map.size() == 1);
    FAIL_IF_NOT(gcomm.put_map[0][0].first == 0);
    FAIL_IF_NOT(gcomm.put_map[0][0].second == 0);
    // put cell 1 face 3
    FAIL_IF_NOT(gcomm.put_map[0][3].first == 2);
    FAIL_IF_NOT(gcomm.put_map[0][3].second == 0);
    // ghost cell 0 (0 global_id) face 1
    FAIL_IF_NOT(gcomm.ghost_map.size() == 2);
    FAIL_IF_NOT(gcomm.ghost_map[0].size() == 1);
    FAIL_IF_NOT(gcomm.ghost_map[0][0] == 0);
    // ghost cell 1 (3 global_id) face 2
    FAIL_IF_NOT(gcomm.ghost_map[1].size() == 1);
    FAIL_IF_NOT(gcomm.ghost_map[1][2] == 1);

    FAIL_IF_NOT(local_ghost_data[0] == 3);
    FAIL_IF_NOT(local_ghost_data[1] == 5);
  }

  // Check cell data
  if (rtt_c4::node() == 2) {
    FAIL_IF_NOT(mesh.number_of_local_cells() == 2);
    FAIL_IF_NOT(mesh.number_of_ghost_cells() == 2);
    FAIL_IF_NOT(gcomm.local_ghost_buffer_size == 2);
    FAIL_IF_NOT(gcomm.put_buffer_size.size() == 2);
    FAIL_IF_NOT(gcomm.put_buffer_size[0] == 1);
    FAIL_IF_NOT(gcomm.put_buffer_size[1] == 1);
    // put cell 2 face 2
    FAIL_IF_NOT(gcomm.put_map.size() == 2);
    FAIL_IF_NOT(gcomm.put_map[0][2].first == 0);
    FAIL_IF_NOT(gcomm.put_map[0][2].second == 0);
    // put cell 3 face 2
    FAIL_IF_NOT(gcomm.put_map[1][2].first == 1);
    FAIL_IF_NOT(gcomm.put_map[1][2].second == 0);
    // ghost cell 0 (0 global_id) face 3
    FAIL_IF_NOT(gcomm.ghost_map.size() == 2);
    FAIL_IF_NOT(gcomm.ghost_map[0].size() == 1);
    FAIL_IF_NOT(gcomm.ghost_map[0][3] == 0);
    // ghost cell 0 (0 global_id) face 3
    FAIL_IF_NOT(gcomm.ghost_map[1].size() == 1);
    FAIL_IF_NOT(gcomm.ghost_map[1][3] == 1);

    FAIL_IF_NOT(local_ghost_data[0] == 3);
    FAIL_IF_NOT(local_ghost_data[1] == 4);
  }

  if (ut.numFails == 0) {
    std::ostringstream m;
    m << "2D DD orthogonal mesh test passed";
    PASSMSG(m.str());
  }
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_c4::ParallelUnitTest ut(argc, argv, rtt_odd::release);
  try {
    // >>> Replicated UNIT TESTS
    if (rtt_c4::nodes() == 2) {
      test_1d_dd_comm(ut);
    } else if (rtt_c4::nodes() == 3) {
      test_2d_dd_comm(ut);
    }
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstGhost_Comm.cc
//------------------------------------------------------------------------------------------------//
