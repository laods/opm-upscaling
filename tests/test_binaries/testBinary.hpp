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

#ifndef TEST_BINARY_HEADER
#define TEST_BINARY_HEADER

#include <string>
#include <iostream>
#include <fstream>

#include <opm/porsol/common/Matrix.hpp>

using namespace std;

namespace Opm {

    typedef FullMatrix<double, OwnData, COrdering> Matrix;

    /// @brief 
    ///    Read result file and store results in a matrix. 
    ///
    /// @param [in] filename
    ///    Path to result file to be read. Lines starting with '#' are ignored.
    ///
    /// @param [in] rows
    ///    Number of rows in result file
    ///
    /// @param [in] cols
    ///    Number of cols in result file
    Matrix readResultFile(string filename, int rows, int cols) {
	
	Matrix result(rows, cols, (const double*)0);
	
	// Check if filename is readable
	char* file = (char*)filename.c_str();
	ifstream resultfile(file, ios::in);
	if (resultfile.fail()) {
	    cerr << "Error: Filename " << filename << " not found or not readable." << endl;
	    exit(1);
	}
	
	// Read file and put results into matrix result
	double value;
	int row = 0;
	int col = 0;
	string line;
	while (getline(resultfile, line)) {
	    if (line[0] == '#') // Skip lines starting with #
		continue;
	    istringstream iss(line);
	    while (iss >> value) {
		// Check if input file contains more rows or columns than specified
		if (row >= rows) {
		    cerr << "Error: More rows in " << filename << " than specified." << endl;
		    exit(1);
		}
		else if (col >= cols) {
		    cerr << "Error: More columns in " << filename << " than specified." << endl;
		    exit(1);
		}
		else {
		    result(row, col) = value;
		    ++col;
		}
	    }
	    col = 0;
	    ++row;
	}
	
	resultfile.close();
	
	return result;	
    }

    /// @brief 
    ///    Test if two matrices are equal within a relative tolerance
    ///
    /// @param [in] refSoln
    ///    Reference solution stored as a matrix
    ///
    /// @param [in] newSoln
    ///    Solution to compare with
    ///
    /// @param [in] relTol
    ///    Relative tolerance
    bool matrixAlmostEqual(Matrix refSoln, Matrix newSoln, double relTol) {

	ASSERT(refSoln.numRows() == newSoln.numRows());
	ASSERT(refSoln.numCols() == newSoln.numCols());

	// Test element by element
	for (int row=0; row<refSoln.numRows(); ++row) {
	    for (int col=0; col<refSoln.numCols(); ++col) {
		bool equal = true;
		if (abs(refSoln(row,col)) < 1e-14) { // If refSoln close to zero do simple difference test
		    double absDiff = abs(refSoln(row,col) - newSoln(row,col));
		    if (absDiff > 1e-6) equal = false;
		}
		else { // If refSoln != 0 do relative difference test
		    double absRelDiff = abs( (refSoln(row,col) - newSoln(row,col)) / refSoln(row,col) );
		    if (absRelDiff > relTol) equal = false;
		}
		if (!equal) {
		    cout << endl << "Verification error: Calculated solution not equal to reference solution "
			 << "within a relative tolerance of " << relTol
			 << " (at position " << row << ", " << col << ")." << endl
			 << "Calculated solution:" << endl
			 << newSoln
			 << "Reference solution:" << endl
			 << refSoln << endl;
		    return false;
		}
	    }
	}
	
	return true;
    }

    
} // namespace Opm
#endif // TEST_BINARY_HEADER
