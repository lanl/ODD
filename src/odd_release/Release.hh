//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   odd_release/Release.hh
 * \author Kelly Thompson <kgt@lanl.gov>
 * \date   June 14, 2011
 * \brief  Copyright function for the odd libraries.
 * \note   Copyright (C) 2011-2021 Triad National Security, LLCC., All rights
 * reserved */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_odd_release_Release_hh
#define rtt_odd_release_Release_hh

#include "odd_release/config.h"
#include <string>

///================================================================================================//
/*!
 * \namespace rtt_odd
 *
 * \brief Namespace that contains the odd package classes and variables.
 */
//================================================================================================//
namespace rtt_odd {

/*! Returns a version string of the form "ODD-0_1_0, build date 2022/07/25; build type: DEBUG;
 *      REPRODUCIBLE, GPU kernel: UNAVAILABLE, Roundoff Mode: ACCURATE, DBC: 7, and
 *      DRACO_DIAGNOSTICS: 0." */
std::string const release(std::string const &codename);

//! This signature used by UnitTest ctor().
std::string const release();

std::string const release(std::string const &package_name_version, std::string const &build_date,
                          std::string const &build_type);

//! Print a 1-line copyright notice.
std::string const copyright();

//! Print a list of authors (current and prior contributors)
std::string const authors(bool const use_doxygen_formatting = false);

/*! Returns a code header block of the form:
 *
 * \verbatim
 *
 * ================================================================================================
 * ODD : odd-0_1_0, build date 2022/07/25; build type: DEBUG; DBC: 7, and
 * DRACO_DIAGNOSTICS: 0.
 *
 * Odd API
 * ----------------------------------------
 * CCS-2 ODD Team: Mathew A. Cleveland, and Andrew Till.
 *
 * ================================================================================================
 *
 * \endverbatim
 */
void printBanner(std::ostream &os, std::string const &codename,
                 std::string const &desc = std::string(""));
} // namespace rtt_odd

#endif // rtt_odd_release_Release_hh

//------------------------------------------------------------------------------------------------//
// end of odd_release/Release.hh
//------------------------------------------------------------------------------------------------//
