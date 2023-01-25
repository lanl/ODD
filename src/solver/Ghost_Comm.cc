//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   solver/Ghost_Comm.cc
 * \author Mathew Cleveland
 * \date   March 21st 2022
 * \brief  Collect Ghost data implementation
 * \note   Copyright (C) 2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "Ghost_Comm.hh"
#include "c4/global.hh"
#include <cmath>

namespace odd_solver {

Ghost_Comm::Ghost_Comm(const Orthogonal_Mesh &mesh) {
  const auto nodes = static_cast<size_t>(rtt_c4::nodes());
  const auto node = static_cast<size_t>(rtt_c4::node());
  std::vector<int> global_ghost_per_proc(nodes * nodes, 0);
  std::vector<int> send_buffer_offset(nodes, 0);
  const auto ncells = mesh.number_of_local_cells();

  local_ghost_buffer_size = 0;
  // Loop over all cells and build up put maps
  for (size_t cell = 0; cell < ncells; cell++) {
    for (size_t face = 0; face < mesh.number_of_faces(cell); face++) {
      const auto ftype = mesh.face_type(cell, face);
      if (ftype == FACE_TYPE::GHOST_FACE) {
        size_t ghost_id = mesh.next_cell(cell, face);
        size_t ghost_proc = mesh.ghost_cell_proc(ghost_id);
        put_local_cell[ghost_proc].push_back(cell);
        put_global_cell[ghost_proc].push_back(mesh.cell_global_id(cell));
        const size_t size = put_global_cell[ghost_proc].size();
        put_face[ghost_proc].push_back(face);
        put_map[cell][face] = {ghost_proc, size - 1};
        put_rank_offset[ghost_proc] = size;
        put_buffer_size[ghost_proc] = size;
        global_ghost_per_proc[ghost_proc * nodes + node]++;
        local_ghost_buffer_size++;
      }
    }
  }
  rtt_c4::global_sum(&global_ghost_per_proc[0], nodes * nodes);
  for (auto &map : put_rank_offset) {
    auto recv_proc = map.first;
    size_t offset = 0;
    for (size_t send_proc = 0; send_proc < node; send_proc++) {
      offset += global_ghost_per_proc[recv_proc * nodes + send_proc];
    }
    map.second = offset;
  }
  build_ghost_map(mesh);
}

#ifdef C4_MPI
//------------------------------------------------------------------------------------------------//
// call MPI_put using a chunk style write to avoid error in MPI_put with large local buffers.
auto chunk_put_lambda = [](auto &put_rank, auto &put_offset, auto &put_buffer, auto &win,
                           MPI_Datatype mpi_data_type) {
  const size_t put_size = put_buffer.size();
  // This is dumb, but we need to write in chunks because MPI_Put writes
  // junk with large (>10,000) buffer sizes.
  int chunk_size = 1000;
  const auto nchunks =
      static_cast<int>(std::ceil(static_cast<double>(put_size) / static_cast<double>(chunk_size)));
  int nput = 0;
  for (int c = 0; c < nchunks; c++) {
    chunk_size = std::min(chunk_size, static_cast<int>(put_size) - nput);
    Check(chunk_size > 0);
    MPI_Put(&put_buffer[nput], chunk_size, mpi_data_type, static_cast<int>(put_rank),
            static_cast<int>(put_offset), chunk_size, mpi_data_type, win);
    nput += chunk_size;
  }
};
#endif

//------------------------------------------------------------------------------------------------//
/*!
 * \brief build_ghost_map data
 *
 *  Build up the ghost map data to make it easy to exchange ghost information 
 * 
 *
 *
 * \param[in] mesh
 */
void Ghost_Comm::build_ghost_map(const Orthogonal_Mesh &mesh) {
  std::vector<size_t> local_global_id_buffer(local_ghost_buffer_size, 32);
  exchange_ghost_data(put_global_cell, local_global_id_buffer);
  std::vector<size_t> local_face_buffer(local_ghost_buffer_size, 32);
  exchange_ghost_data(put_face, local_face_buffer);

  for (size_t i = 0; i < local_ghost_buffer_size; i++) {
    size_t buffer_global_id = local_global_id_buffer[i];
    size_t buffer_face_id = local_face_buffer[i];
    for (size_t g = 0; g < mesh.number_of_ghost_cells(); g++) {
      size_t gid = mesh.ghost_cell_global_id(g);
      if (gid == buffer_global_id)
        ghost_map[g][buffer_face_id] = i;
    }
  }
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief exchange ghost size_t data
 *
 * Collect ghost data for vector of 3 dimensional arrays. This function uses RMA and the local
 * put_window_map to allow each rank to independently fill in its data to ghost cells of other
 * ranks.
 *
 * \param[in] put_data Ghost data to be put on the other ranks
 * \param[inout] ghost_data vector to write the new ghost data to
 */
void Ghost_Comm::exchange_ghost_data(const std::map<size_t, std::vector<size_t>> &put_data,
                                     std::vector<size_t> &ghost_data) const {
#ifdef C4_MPI // temporary work around until RMA is available in c4
  Require(ghost_data.size() == local_ghost_buffer_size);
  Require(put_data.size() == put_buffer_size.size());

  // Use one sided MPI Put commands to fill up the ghost cell location data
  MPI_Win win;
  MPI_Win_create(ghost_data.data(), local_ghost_buffer_size * sizeof(size_t), sizeof(size_t),
                 MPI_INFO_NULL, MPI_COMM_WORLD, &win);
  // working from my local data put the ghost data on the other ranks
  Remember(int errorcode =) MPI_Win_fence(MPI_MODE_NOSTORE, win);
  Check(errorcode == MPI_SUCCESS);
  // temporary work around until RMA is available in c4
  // loop over all ranks we need to send this buffer too.
  for (auto &put : put_data) {
    const size_t put_rank = put.first;
    auto &put_buffer = put.second;
    // use map.at() to allow const access
    const size_t put_offset = put_rank_offset.at(put_rank);
    chunk_put_lambda(put_rank, put_offset, put_buffer, win, MPI_UNSIGNED_LONG);
  }
  Remember(errorcode =) MPI_Win_fence((MPI_MODE_NOSTORE | MPI_MODE_NOSUCCEED), win);
  Check(errorcode == MPI_SUCCESS);
  Remember(errorcode =) MPI_Win_fence((MPI_MODE_NOSTORE | MPI_MODE_NOSUCCEED), win);
  Check(errorcode == MPI_SUCCESS);
  MPI_Win_free(&win);
#endif
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief exchange ghost double data
 *
 * Collect ghost data for vector of 3 dimensional arrays. This function uses RMA and the local
 * put_window_map to allow each rank to independently fill in its data to ghost cells of other
 * ranks.
 *
 * \param[in] put_data Ghost data to be put on the other ranks
 * \param[inout] ghost_data vector to write the new ghost data to
 */
void Ghost_Comm::exchange_ghost_data(const std::map<size_t, std::vector<double>> &put_data,
                                     std::vector<double> &ghost_data) const {
#ifdef C4_MPI // temporary work around until RMA is available in c4
  Require(ghost_data.size() == local_ghost_buffer_size);
  Require(put_data.size() == put_buffer_size.size());

  // Use one sided MPI Put commands to fill up the ghost cell location data
  MPI_Win win;
  MPI_Win_create(ghost_data.data(), local_ghost_buffer_size * sizeof(double), sizeof(double),
                 MPI_INFO_NULL, MPI_COMM_WORLD, &win);
  // working from my local data put the ghost data on the other ranks
  Remember(int errorcode =) MPI_Win_fence(MPI_MODE_NOSTORE, win);
  Check(errorcode == MPI_SUCCESS);
  // temporary work around until RMA is available in c4
  // loop over all ranks we need to send this buffer too.
  for (auto &put : put_data) {
    const size_t put_rank = put.first;
    auto &put_buffer = put.second;
    // use map.at() to allow const access
    const size_t put_offset = put_rank_offset.at(put_rank);
    chunk_put_lambda(put_rank, put_offset, put_buffer, win, MPI_DOUBLE);
  }
  Remember(errorcode =) MPI_Win_fence((MPI_MODE_NOSTORE | MPI_MODE_NOSUCCEED), win);
  Check(errorcode == MPI_SUCCESS);
  Remember(errorcode =) MPI_Win_fence((MPI_MODE_NOSTORE | MPI_MODE_NOSUCCEED), win);
  Check(errorcode == MPI_SUCCESS);
  MPI_Win_free(&win);
#endif
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief exchange ghost double data packed with a fixed stride
 *
 * Collect ghost data for vector of 3 dimensional arrays. This function uses RMA and the local
 * put_window_map to allow each rank to independently fill in its data to ghost cells of other
 * ranks.
 *
 * \param[in] put_data Ghost data to be put on the other ranks
 * \param[inout] ghost_data vector to write the new ghost data to
 * \param[in] stride of packed data to be exchanged
 */
void Ghost_Comm::exchange_ghost_data(const std::map<size_t, std::vector<double>> &put_data,
                                     std::vector<double> &ghost_data, const size_t stride) const {
#ifdef C4_MPI // temporary work around until RMA is available in c4
  Require(ghost_data.size() == local_ghost_buffer_size * stride);
  Require(put_data.size() == put_buffer_size.size());

  // Use one sided MPI Put commands to fill up the ghost cell location data
  MPI_Win win;
  MPI_Win_create(ghost_data.data(), local_ghost_buffer_size * sizeof(double) * stride,
                 sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win);
  // working from my local data put the ghost data on the other ranks
  Remember(int errorcode =) MPI_Win_fence(MPI_MODE_NOSTORE, win);
  Check(errorcode == MPI_SUCCESS);
  // temporary work around until RMA is available in c4
  // loop over all ranks we need to send this buffer too.
  for (auto &put : put_data) {
    const size_t put_rank = put.first;
    auto &put_buffer = put.second;
    // use map.at() to allow const access
    const size_t put_offset = put_rank_offset.at(put_rank) * stride;
    chunk_put_lambda(put_rank, put_offset, put_buffer, win, MPI_DOUBLE);
  }
  Remember(errorcode =) MPI_Win_fence((MPI_MODE_NOSTORE | MPI_MODE_NOSUCCEED), win);
  Check(errorcode == MPI_SUCCESS);
  Remember(errorcode =) MPI_Win_fence((MPI_MODE_NOSTORE | MPI_MODE_NOSUCCEED), win);
  Check(errorcode == MPI_SUCCESS);
  MPI_Win_free(&win);
#endif
}

} // namespace odd_solver

//------------------------------------------------------------------------------------------------//
// end of <basename>.cc
//------------------------------------------------------------------------------------------------//
