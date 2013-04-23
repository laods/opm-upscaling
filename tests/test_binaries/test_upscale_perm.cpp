/*
  Copyright 2013 Statoil ASA

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "config.h"

#include <opm/core/utility/have_boost_redef.hpp>

#if defined(HAVE_DYNAMIC_BOOST_TEST)
#define BOOST_TEST_DYN_LINK
#endif

#define BOOST_TEST_MODULE UpscalePerm
#include <boost/test/unit_test.hpp>

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

#include "testBinary.hpp"

using namespace std;

BOOST_AUTO_TEST_CASE(upscale_perm_fixed_PeriodicTilted)	
{
    // --------------------------------------------------------------
    // INPUT DATA - edit this section to create new test
    // --------------------------------------------------------------

    // Specify accepted relative tolerance
    double relTol = 1e-6;

    // Specify number of rows and columns in result file
    int rows = 3;
    int cols = 3;
	
    // Specify relative paths from directory 'tests' in the build directory
    string gridPath = "input_data/grids/PeriodicTilted.grdecl";
    string refSolnPath = "input_data/reference_solutions/upscale_perm_fixed_PeriodicTilted.txt";
    string outputPath = "temp_results.txt";
	
    // Specify binary and input arguments that should be run
    // This stream should be equal to terminal command 
    stringstream runCommand;
    runCommand << "../bin/upscale_perm "
	       << "-output " << outputPath << " "
	       << gridPath;

    // --------------------------------------------------------------
    // RUN TEST 
    // --------------------------------------------------------------

    // Run executable
    system(runCommand.str().c_str());

    // Store solutions in Matrix object
    Opm::Matrix newSoln = Opm::readResultFile(outputPath, rows, cols);
    Opm::Matrix refSoln = Opm::readResultFile(refSolnPath, rows, cols);
	
    // Remove temporary result file
    string removeCommand = string("rm ") + outputPath;
    system(removeCommand.c_str());
	
    // Compare result with reference solution
    bool test = Opm::matrixAlmostEqual(refSoln, newSoln, relTol);
    BOOST_CHECK(test);
}
   
