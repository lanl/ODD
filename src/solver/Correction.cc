//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   solver/Correction.hh
 * \author <user>
 * \brief  Define class Correction
 * \note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "Correction.hh"
#include "Constants.hh"
#include <cmath>
#include <iostream>

namespace odd_solver {
namespace Correction {

void calc_correction(double &epsilon, double &correction_source, double &Tstar, double eden,
                     const double sigma_a, const double sigma_planck, const double fleck,
                     const double cv, const double volume, const double dt, const double T0,
                     const double ext_exp_source) {
  double error = 1.0;
  const double convergence = 1e-6;
  const double dep = fleck * eden * sigma_a * constants::c * dt * volume;
  const double sampled =
      fleck * sigma_a * constants::a * constants::c * T0 * T0 * T0 * T0 * dt * volume;
  const double Sc = 1.0 / fleck * (dep + ext_exp_source - sampled);
  size_t count = 0;
  while (error > convergence && count < 1000) {
    double last_Tstar = Tstar;
    count++;
    const double Tc = std::pow(0.25 * (Tstar + T0) * (Tstar * Tstar + T0 * T0), 1.0 / 3.0);
    const double fleck_c =
        1.0 / (1.0 + constants::a * constants::c * sigma_planck * 4.0 * Tc * Tc * Tc * dt / cv);
    epsilon = fleck_c - fleck;
    const double delta_E = dep - sampled + ext_exp_source + epsilon * Sc;
    const double q = 1.0 / (fleck_c) * ((Tstar - T0) * cv - delta_E / volume);
    const double dqdT =
        4.0 * sigma_planck * constants::a * constants::c * dt * Tstar * Tstar * Tstar + cv;
    Tstar = (Tstar - q / dqdT);
    error = fabs(Tstar - last_Tstar) / Tstar;
    correction_source = -epsilon * Sc / volume;
  }
  std::cout << " T0 = " << T0;
  std::cout << " Tstar = " << Tstar;
  std::cout << " Tr0 = " << std::pow(eden / constants::a, 0.25);
  std::cout << " TrStar = " << std::pow((eden + correction_source) / constants::a, 0.25);
  std::cout << " correction_source " << correction_source;
  std::cout << " count = " << count;
  std::cout << " eps = " << epsilon;
  std::cout << std::endl;
  //Check(!(correction_source < 0.0));
}

} // end namespace Correction

} // end namespace odd_solver

//------------------------------------------------------------------------------------------------//
// end of solver/Correction.cc
//------------------------------------------------------------------------------------------------//

