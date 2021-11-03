//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   api/test/tstArguments.cc
 * \author Mathew Cleveland
 * \date   October 21st 2021
 * \brief  Testing opacity reader
 * \note   Copyright (C) 2021 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "api/Arguments.hh"
#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"

using namespace rtt_dsxx;
using namespace odd_api;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//
void test(rtt_dsxx::UnitTest &ut) {
  Arguments arg;
  // Should fail
  try {
    arg.control_data.check_arguments();
  } catch (...) {
    std::cout << "Caught the expected Control_Data::check_arguments assertions" << std::endl;
  }
  // fill control data with "valid junk" and check arguments
  arg.control_data.opacity_file = "not_empty";
  arg.control_data.check_arguments();

  try {
    arg.zonal_data.check_arguments();
  } catch (...) {
    std::cout << "Caught the expected Zonal_Data::check_arguments assertions" << std::endl;
  }
  // fill the zonal data with "valid junk" and check arguments
  arg.zonal_data.number_of_cells = 2;
  arg.zonal_data.dimensions = 1;
  arg.zonal_data.dx = 31.0;
  std::vector<size_t> matids = {19000, 190001};
  arg.zonal_data.matid = &matids[0];

  // Check that the pointer matches the data
  FAIL_IF_NOT((*arg.zonal_data.matid) == matids[0]);
  arg.zonal_data.matid++;
  FAIL_IF_NOT((*arg.zonal_data.matid) == matids[1]);
  arg.zonal_data.check_arguments();

  try {
    arg.output_data.check_arguments();
  } catch (...) {
    std::cout << "Caught the expected Output_Data::check_arguments assertions" << std::endl;
  }
  std::vector<double> opacity_data = {17000.0, 17000.1};
  arg.output_data.opacity_data = &opacity_data[0];
  // Check that the pointer matches the data
  FAIL_IF_NOT((*arg.output_data.opacity_data) == opacity_data[0]);
  arg.output_data.opacity_data++;
  FAIL_IF_NOT((*arg.output_data.opacity_data) == opacity_data[1]);
  arg.output_data.check_arguments();

  if (ut.numFails == 0) {
    std::ostringstream m;
    m << "the ODD API arguments seem to be working";
    PASSMSG(m.str());
  }
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    // >>> UNIT TESTS
    test(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstArgumetns.cc
//------------------------------------------------------------------------------------------------//
