//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   api/Interface_Builder.hh
 * \author Mathew Cleveland
 * \brief  Interface constructor build from the Arguments API
 * \note   Copyright (C) 2018-2020 Triad National Security, LLC.
 *         All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef odd_api_Interface_Builder_hh
#define odd_api_Interface_Builder_hh

#include "Arguments.hh"
#include "solver/Interface_Data.hh"

namespace odd_api {

odd_solver::Interface_Data build_interface_data(const Arguments &arg);

} // end namespace odd_api

#endif // odd_api_Interface_Builder_hh

//------------------------------------------------------------------------------------------------//
// end of api/Interface_Builder.hh
//------------------------------------------------------------------------------------------------//
