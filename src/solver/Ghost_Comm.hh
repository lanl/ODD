//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   solver/Ghost_Comm.hh
 * \author Mathew Cleveland
 * \brief  Collect ghost data between processors
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#ifndef odd_solver_Ghost_Comm_hh
#define odd_solver_Ghost_Comm_hh

#include "Orthogonal_Mesh.hh"
#include <map>

namespace odd_solver {

//================================================================================================//
/*!
 * \class Ghost_Comm
 * \brief
 *
 * Communicate ghost information between ranks
 */
//================================================================================================//

class Ghost_Comm {
public:
  //! Default constructors.
  Ghost_Comm(const Orthogonal_Mesh &mesh);

  // DATA
  std::map<size_t, size_t> put_buffer_size;
  // put_mpa[local_cell_id][face_id] = {rank_id, put_buffer_id}
  std::map<size_t, std::map<size_t, std::pair<size_t, size_t>>> put_map;
  size_t local_ghost_buffer_size;
  // ghost_mpa[local_ghost_id][face_id][ghost_buffer_id]
  std::map<size_t, std::map<size_t, size_t>> ghost_map;

  // Interface functions
  void exchange_ghost_data(const std::map<size_t, std::vector<size_t>> &put_data,
                           std::vector<size_t> &ghost_data) const;
  void exchange_ghost_data(const std::map<size_t, std::vector<double>> &put_data,
                           std::vector<double> &ghost_data) const;

private:
  // Initialize ghost maps
  void build_ghost_map(const Orthogonal_Mesh &mesh);
  // Private comm data
  std::map<size_t, std::vector<size_t>> put_local_cell;
  std::map<size_t, std::vector<size_t>> put_face;
  std::map<size_t, std::vector<size_t>> put_global_cell;
  std::map<size_t, size_t> put_rank_offset;
};

} // namespace odd_solver

#endif // odd_solver_Ghost_Comm_hh

//------------------------------------------------------------------------------------------------//
// end of solver/Ghost_Comm.hh
//------------------------------------------------------------------------------------------------//
