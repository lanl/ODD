//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   solver/Orthogonal_Mesh.cc
 * \author Mathew Cleveland
 * \date   Dec. 15 2021
 * \brief  Functions to calculate mesh info from the mesh interface data
 * \note   Copyright (C) 2021 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "Orthogonal_Mesh.hh"
#include "ds++/dbc.hh"

namespace odd_solver {

//================================================================================================//
/*!
 * \brief Get face type
 *
 * \param[in] cell id (local)
 * \param[in] face id
 * 
 * \return face type (interior, boundary, or ghost
 */
//================================================================================================//
size_t Orthogonal_Mesh::face_type(const size_t cell, const size_t face) const {
  Require(cell >= 0 && cell < mesh_data.number_of_local_cells);
  Require(face >= 0 && face < number_of_faces(cell));
  const size_t index = cell * number_of_faces(cell) + face;
  Check(index >= 0 && index < mesh_data.face_types.size());
  const size_t face_t = mesh_data.face_types[index];
  Ensure(face_t >= 0 && face_t < FACE_TYPE::N_FACE_TYPES);
  return face_t;
}

//================================================================================================//
/*!
 * \brief Get the cell on the other side of a face
 *
 * \param[in] cell id (local)
 * \param[in] face id
 * 
 * \return cell_id (local or ghost)
 */
//================================================================================================//
size_t Orthogonal_Mesh::next_cell(const size_t cell, const size_t face) const {
  Require(cell >= 0 && cell < mesh_data.number_of_local_cells);
  Require(face >= 0 && face < number_of_faces(cell));
  Remember(size_t face_t = face_type(cell, face));
  // Next cell is only valid for internal or ghost faces
  Check(face_t == FACE_TYPE::INTERNAL_FACE || face_t == FACE_TYPE::GHOST_FACE);
  const size_t index = cell * number_of_faces(cell) + face;
  Check(index >= 0 && index < mesh_data.next_cell_id.size());
  const size_t cell_id = mesh_data.next_cell_id[index];
  Ensure(cell_id >= 0 && face_t == FACE_TYPE::INTERNAL_FACE
             ? cell_id < mesh_data.number_of_local_cells
             : cell_id < mesh_data.number_of_ghost_cells);
  return cell_id;
}

//================================================================================================//
/*!
 * \brief Calculate the cell volume
 *
 * \param[in] cell id (local)
 * 
 * \return volume of cell
 */
//================================================================================================//
double Orthogonal_Mesh::cell_volume(const size_t cell) const {
  Require(cell >= 0 && cell < mesh_data.number_of_local_cells);
  double volume = mesh_data.cell_size[cell * 3];
  for (size_t d = 1; d < mesh_data.n_dims; d++) {
    volume *= mesh_data.cell_size[cell * 3 + d];
  }
  Ensure(volume > 0.0);
  return volume;
}

//================================================================================================//
/*!
 * \brief Calculate the distance from the cell center to the face
 *
 * \param[in] cell id (local)
 * \param[in] face (local)
 * 
 * \return distance to cell face
 */
//================================================================================================//
double Orthogonal_Mesh::distance_center_to_face(const size_t cell, const size_t face) const {
  Require(cell >= 0 && cell < mesh_data.number_of_local_cells);
  Require(face >= 0 && face < number_of_faces(cell));
  size_t face_d = face < 2 ? 0 : (face < 4 ? 1 : 2);
  return 0.5 * mesh_data.cell_size[cell * 3 + face_d];
}

//================================================================================================//
/*!
 * \brief Get the face area of a cell
 *
 * \param[in] cell id (local)
 * \param[in] face id
 * 
 * \return area of a cell
 */
//================================================================================================//
double Orthogonal_Mesh::face_area(const size_t cell, const size_t face) const {
  Require(cell >= 0 && cell < mesh_data.number_of_local_cells);
  Require(face >= 0 && face < number_of_faces(cell));
  size_t face_d = face < 2 ? 0 : (face < 4 ? 1 : 2);
  Check(face_d < mesh_data.n_dims);
  double area = 1.0;
  for (size_t d = 0; d < mesh_data.n_dims; d++) {
    area *= face_d == d ? 1.0 : mesh_data.cell_size[cell * 3 + d];
  }
  return area;
}

//================================================================================================//
/*!
 * \brief Get the cell center
 *
 * \param[in] cell id (local)
 * 
 * \return center location of the cell
 */
//================================================================================================//
std::array<double, 3> Orthogonal_Mesh::cell_center(const size_t cell) const {
  Require(cell >= 0 && cell < mesh_data.number_of_local_cells);
  const std::array<double, 3> center{mesh_data.cell_position[cell * 3],
                                     mesh_data.cell_position[cell * 3 + 1],
                                     mesh_data.cell_position[cell * 3 + 2]};
  return center;
}

} // namespace odd_solver

//------------------------------------------------------------------------------------------------//
// end of Orthogonal_Mesh.cc
//------------------------------------------------------------------------------------------------//
