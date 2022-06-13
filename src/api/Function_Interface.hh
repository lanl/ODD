//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   api/Function_Interface.hh
 * \author Mathew Cleveland
 * \brief  Define the C interface functions
 * \note   Copyright (C) 2021-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#ifndef api_Function_Interface_hh
#define api_Function_Interface_hh

#include "api/Arguments.hh"

extern "C" {

//! Odd solver interface call
void Odd_Diffusion_Solve(Arguments &arg);

void MPI_Initialize(int argc, char *argv[]);

} // end extern "C"

#endif // api_Function_Interface_hh

//------------------------------------------------------------------------------------------------//
// end of api/Function_Interface.hh
//------------------------------------------------------------------------------------------------//
