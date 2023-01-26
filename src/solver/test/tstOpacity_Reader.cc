//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   solver/test/tstOpacity_Reader.cc
 * \author Mathew Cleveland
 * \date   October 21st 2021
 * \brief  Testing opacity reader
 * \note   Copyright (C) 2021 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "solver/Opacity_Reader.hh"
#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"

using namespace rtt_dsxx;
using namespace odd_solver;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//
void test_print(rtt_dsxx::UnitTest &ut) {
  Opacity_Reader op_reader(ut.getTestSourcePath() + "two-mats.ipcress");
  size_t ngroups = op_reader.ngroups;
  if (ngroups != 33)
    ITFAILS;
  op_reader.print_available_mats();
  op_reader.print_available_data(10001);
  op_reader.print_available_data(10002);
  if (ut.numFails == 0) {
    std::ostringstream m;
    m << "opacity reader print test passed";
    PASSMSG(m.str());
  }
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    // >>> UNIT TESTS
    test_print(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstOpacity_Reader.cc
//------------------------------------------------------------------------------------------------//
