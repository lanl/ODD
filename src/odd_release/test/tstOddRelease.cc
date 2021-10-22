//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   odd_release/test/tstOddRelease.cc
 * \author Mathew Cleveland
 * \date   October 22nd 2021
 * \brief  Ensure that the Release.hh/cc files compile.
 * \note   Copyright (C) 2011-2021 Triad National Security, LLCC., All rights
 * reserved. */
//------------------------------------------------------------------------------------------------//

#include "ds++/DracoStrings.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"
#include "odd_release/Release.hh"
#include <iostream>
#include <sstream>

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_odd;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

void maintest(UnitTest &ut) {
  // Print the release information
  std::ostringstream const releaseString(release("Dummycode"));

  if (releaseString.str().length() > 0)
    ut.passes("releaseString len > 0");
  else
    ut.failure("releaseString len == 0");

  // Print the copyright statement and author list
  std::ostringstream const copyrightString(copyright());

  if (copyrightString.str().length() > 0)
    ut.passes("copyrightString len > 0");
  else
    ut.failure("copyrightString len == 0");

  // Check the output
  bool verbose(false);
  std::map<std::string, unsigned> wc = rtt_dsxx::get_word_count(releaseString, verbose);

  if (wc[string("DRACO_DIAGNOSTICS")] != 1)
    ITFAILS;
  if (wc[string("built")] != 1)
    ITFAILS;

  // Test the banner function
  printBanner(std::cout, "dummyCode", "This is a dummy code.");

  return;
}

//------------------------------------------------------------------------------------------------//

int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    maintest(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstOddRelease.cc
//------------------------------------------------------------------------------------------------//
