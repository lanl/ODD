//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   odd_release/Release.cc
 * \author Mathew Cleveland
 * \date   October 22nd 2021
 * \brief  Copyright function implementation for odd libraries.
 * \note   Copyright (C) 2011-2021 Triad National Security, LLC., All rights
 * reserved. */
//------------------------------------------------------------------------------------------------//

#include "Release.hh"
#include "ds++/DracoTerminal.hh"
#include "ds++/Release.hh"
#include <functional>
#include <iostream>
#include <map>
#include <sstream>

namespace rtt_odd {

//------------------------------------------------------------------------------------------------//
/*!
 * \return The release string, including information about rng and dbc.
 *
 * Function definition for Release, define the local version number for this
 * library in the form odd-\#_\#_\# in pkg_release variable
 *
 * No package specific CPP macros values in this function as it is also called
 * from api and mcgrid.
 *
 * \code
 * dummyCode: odd-7_11_20160729, built: 2016/07/29, Debug, DBC: 7,
 * and DRACO_DIAGNOSTICS: 0.
 * \endcode
 */
std::string const release(std::string const &package_name_version, std::string const &build_date,
                          std::string const &build_type) {
  using rtt_dsxx::fomdev;
  using rtt_dsxx::mmdevs;
  using rtt_dsxx::print_devs;
  std::ostringstream pkg_release;

  // build a list of information to print in the banner
  mmdevs info_items;

  info_items.insert(fomdev(100, std::string("built: ") + build_date));
  info_items.insert(fomdev(99, build_type));

#ifdef DBC
  std::ostringstream dbc_string;
  dbc_string << "DBC: " << DBC;
  info_items.insert(fomdev(80, dbc_string.str()));
#endif

#ifdef DRACO_DIAGNOSTICS
  std::ostringstream diag_string;
  diag_string << "DRACO_DIAGNOSTICS: " << DRACO_DIAGNOSTICS;
  info_items.insert(fomdev(70, diag_string.str()));
#endif

#ifdef DRACO_DIAGNOSTICS_LEVEL_3
#ifdef FPETRAP_SUPPORTED
  info_items.insert(fomdev(60, "FPE_TRAP: ON"));
#endif
#endif

  size_t const maxlinelen(90);
  pkg_release << print_devs(maxlinelen, package_name_version + ", ", info_items);

  return pkg_release.str();
}
//------------------------------------------------------------------------------------------------//
/*!
 * \return string of the release number
 *
 * Function definition for Release, define the local version number for this
 * library in the form odd-\#_\#_\# in pkg_release variable
 */
std::string const release(std::string const &codename) {
  std::ostringstream package_name_version;
  if (codename.length() > 0)
    package_name_version << codename << ": ";
  package_name_version << Term::ccolor(Term::style::bold)
                       << Term::ccolor(Term::fg::cyan) << "Odd-"
                       << Odd_VERSION_MAJOR << "_" << Odd_VERSION_MINOR << "_"
                       << Odd_VERSION_PATCH << Term::ccolor(Term::fg::reset)
                       << Term::ccolor(Term::style::reset);

  std::string const build_date(ODD_BUILD_DATE);
  std::string const build_type(ODD_BUILD_TYPE);

  return release(package_name_version.str(), build_date, build_type);
}

//------------------------------------------------------------------------------------------------//
std::string const release() { return release(std::string("")); }

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Returns an author list
 *
 * \param[in] use_doxygen_formatting If true, use extra decoration in the
 * output.
 *
 * Order of names is based on lines-of-code as a figure-of-merit (FOM) listed by
 * name by 'git log -numstat' for odd.  See draco/regression/alist.sh.
 */
std::string const authors(bool const use_doxygen_formatting) {

  using rtt_dsxx::fomdev;
  using rtt_dsxx::mmdevs;
  using rtt_dsxx::print_devs;
  std::ostringstream alist;

  // Order selected

  // 2. Active contributers by LOC
  mmdevs current_developers;

  current_developers.insert(fomdev(0, "Mathew A. Cleveland"));
  current_developers.insert(fomdev(0, "Andrew T. Till"));

  // 3. Inactive contributers by LOC

  //   mmdevs prior_developers;
  //   prior_developers.insert(fomdev(45217, "Allan B. Wollaber"));
  //   prior_developers.insert(fomdev(11742, "Robert B. Lowrie"));
  //   prior_developers.insert(fomdev(8486, "Timothy M. Kelley"));
  //   prior_developers.insert(fomdev(3830, "Gabriel M. Rockefeller"));
  //   prior_developers.insert(fomdev(983, "HyeongKae Park"));
  //   prior_developers.insert(fomdev(76, "Todd J. Urbatsch"));
  //   prior_developers.insert(fomdev(19, "Paul Talbot"));
  //   prior_developers.insert(fomdev(4, "Jeff D. Densmore"));
  //
  //   // Previous authors with no current LOC attribution:
  //   prior_developers.insert(fomdev(1, "Aimee K. Hungerford"));
  //   prior_developers.insert(fomdev(1, "Michael W. Buksas"));
  //   prior_developers.insert(fomdev(1, "Paul J. Henning"));
  //   prior_developers.insert(fomdev(1, "Ryan McClarren"));
  //   prior_developers.insert(fomdev(1, "Scott W. Mosher"));
  //   prior_developers.insert(fomdev(1, "Seth R. Johnson"));
  //   prior_developers.insert(fomdev(1, "Tom M. Evans"));

  size_t maxlinelen(100);
  std::string line_name("CCS-2 ODD Team: ");
  if (use_doxygen_formatting) {
    maxlinelen = 400;
    alist << "\\par " << line_name << "\n\n";
    line_name = "";
  } else {
    line_name = Term::ccolor(Term::style::bold) + Term::ccolor(Term::fg::cyan) + line_name +
                Term::ccolor(Term::fg::reset) + Term::ccolor(Term::style::reset);
  }
  alist << print_devs(maxlinelen, line_name, current_developers) << "\n";

  //   line_name = "Prior Contributors: ";
  //   if (use_doxygen_formatting) {
  //     alist << "\\par " << line_name << "\n\n";
  //     line_name = "";
  //   } else {
  //     line_name = Term::ccolor(Term::style::bold) +
  //     Term::ccolor(Term::fg::cyan) + line_name +
  //                 Term::ccolor(Term::fg::reset) +
  //                 Term::ccolor(Term::style::reset);
  //   }
  //   alist << print_devs(maxlinelen, line_name, prior_developers);

  return alist.str();
}

//------------------------------------------------------------------------------------------------//
//! Returns the odd copyright text block
std::string const copyright() {
  std::ostringstream msg;
  msg << "\n"
      << Term::ccolor(Term::fg::green) << "Copyright (C) 1998-2021 Triad National "
      << "Security, LLC. (C19029, LA-CC-16-001)\n     OD999 Export Controlled Software\n"
      << Term::ccolor(Term::fg::reset) << std::endl;
  return msg.str();
}

//------------------------------------------------------------------------------------------------//
//! Print a code banner
void printBanner(std::ostream &os, std::string const &codename,
                 std::string const &desc) {
  std::ostringstream out;
  std::string const divider1(100, '=');
  std::string const divider2(40, '-');

  out << divider1 << "\n" << release(codename) << "\n";
  if (desc.size() > 0)
    out << desc;
  out << "\n" << divider2 << "\n" << authors() << copyright() << "\n" << divider1;
  os << out.str() << std::endl;
  return;
}

} // namespace rtt_odd

//------------------------------------------------------------------------------------------------//
// end of Release.cc
//------------------------------------------------------------------------------------------------//
