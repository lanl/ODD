//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   solver/Orthogonal_Mesh.hh
 * \author Mathew Cleveland
 * \brief  Define class Orthogonal_Mesh
 * \note   Copyright (C) 2021-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#ifndef solver_Orthogonal_Mesh_hh
#define solver_Orthogonal_Mesh_hh

#include "Interface_Data.hh"
#include <array>

namespace odd_solver {

//================================================================================================//
/*!
 * \class Orthogonal Mesh
 * \brief
 *
 * Simple class with minimum mesh data provided by the interface and a collection of data
 * construction functions (avoid storing temporary data where possible)
 *
 */
//================================================================================================//

class Orthogonal_Mesh {
public:
  //================================================================================================//
  /*!
 * \brief Orthogonal mesh constructor class
 *
 * \param[in] iface odd interface 
 *
 */
  //================================================================================================//
  Orthogonal_Mesh(const Mesh_Data mData) : mesh_data(std::move(mData)) {}

  //! get the number of local cells
  size_t number_of_local_cells() const { return mesh_data.number_of_local_cells; }

  //! get the number of local cells
  size_t number_of_global_cells() const { return mesh_data.number_of_global_cells; }

  //! get the number of local cells
  size_t number_of_ghost_cells() const { return mesh_data.number_of_ghost_cells; }

  //! Get the number of faces for each cell
  size_t number_of_faces(const size_t /*cell*/) const { return mesh_data.n_dims * 2; }

  //! Get the cells global id
  size_t cell_global_id(const size_t cell) const { return mesh_data.cell_global_id[cell]; }

  //! get face type
  size_t face_type(const size_t cell, const size_t face) const;

  //! get next cell based on current cell and face
  size_t next_cell(const size_t cell, const size_t face) const;

  //! get volume
  double cell_volume(const size_t cell) const;

  //! get the distance from the cell center to the face center
  double distance_center_to_face(const size_t cell, const size_t face) const;

  //! get face area
  double face_area(const size_t cell, const size_t face) const;

  //! get cell center
  std::array<double, 3> cell_center(const size_t cell) const;

  //! get the topology
  bool domain_decomposed() const { return mesh_data.domain_decomposed; }

private:
  // DATA
  const Mesh_Data mesh_data;
};

} // end namespace odd_solver

#endif // solver_Orthogonal_Mesh_hh

//------------------------------------------------------------------------------------------------//
// end of solver/Orthogonal_Mesh.hh
//------------------------------------------------------------------------------------------------//
