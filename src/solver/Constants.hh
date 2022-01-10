//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   solver/Constants.hh
 * \author Mathew Cleveland
 * \brief  Hold the physical constants that we need
 * \note   Copyright (C) 2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#ifndef odd_solver_constants_hh
#define odd_solver_constants_hh

namespace odd_solver {
namespace constants {

//! CGSH physical constants
constexpr rtt_units::PhysicalConstexprs<rtt_units::CGSH> phys_constants;

//! Speed of light
constexpr double c = phys_constants.speedOfLight(); // cm/sh

//! Radiation constant
constexpr double a = 4.0 * phys_constants.stefanBoltzmann() / phys_constants.speedOfLight();

} // namespace constants

} // end namespace odd_solver

#endif // odd_solver_constants_hh

//------------------------------------------------------------------------------------------------//
// end of solver/Constants.hh
//------------------------------------------------------------------------------------------------//
