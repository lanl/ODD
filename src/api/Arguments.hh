//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   api/Arguments.hh
 * \author Mathew Cleveland
 * \brief  API argument definitions 
 * \note   Copyright (C) 2018-2020 Triad National Security, LLC.
 *         All rights reserved. */
//------------------------------------------------------------------------------------------------//

// clang-format off

#ifndef api_arguments_hh
#define api_arguments_hh

#include <string>

namespace odd_api {

struct Control_Data{
    std::string opacity_file;
    // Constructor to initialize data
    Control_Data();
    void check_arguments();
};

struct Zonal_Data{
    size_t number_of_cells;
    size_t dimensions;
    double dx;
    double dy;
    double dz;
    size_t *matid;
    // constructor to initialize data
    Zonal_Data();
    void check_arguments();
};

struct Output_Data{
    double *opacity_data;
    // constructor to initialize data
    Output_Data();
    void check_arguments();
};

struct Arguments{
    Control_Data control_data;
    Zonal_Data zonal_data;
    Output_Data output_data;
    // constructor to initialize all data fields
    Arguments();
};


} // end namespace odd_api

#endif // api_arguments_hh

//------------------------------------------------------------------------------------------------//
// end of api/arguments.hh
//------------------------------------------------------------------------------------------------//

