//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   solver/test/tstopacity_reader.cc
 * \author Mathew Cleveland
 * \date   October 21st 2021
 * \brief  Testing opacity reader
 * \note   Copyright (C) 2018-2020 Triad National Security, LLC.
 *         All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "solver/opacity_reader.hh"

using namespace rtt_dsxx;
using namespace odd_solver;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//
void test_print(rtt_dsxx::UnitTest &ut){
  opacity_reader op_reader("two_mats.ipcress");
  op_reader.print_file_info();
  op_reader.print_mat_info(1);
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
// end of tstopacity_reader.cc
//------------------------------------------------------------------------------------------------//

