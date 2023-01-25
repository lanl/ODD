//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   solver/Correction.hh
 * \author <user>
 * \brief  Define class Correction
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#ifndef odd_solver_Correction_hh
#define odd_solver_Correction_hh

#include <cmath>
#include <iostream>

namespace odd_solver {
namespace Correction {

void calc_correction(double &epsilon, double &correction_source, double &Tstar, double eden,
                     const double sigma_a, const double sigma_planck, const double fleck,
                     const double cv, const double volume, const double dt, const double T0,
                     const double ext_exp_source);
} // end namespace Correction

} // end namespace odd_solver

#endif // solver_Correction_hh

//------------------------------------------------------------------------------------------------//
// end of solver/Correction.hh
//------------------------------------------------------------------------------------------------//

