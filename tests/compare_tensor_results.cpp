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

#include <string>
#include <iostream>
#include <fstream>
#include <vector>

#include <dune/common/fmatrix.hh>

using namespace std;

typedef Dune::FieldMatrix<double, 3, 3> Matrix;


/// @brief
///    Read result file with tensors and store results in a vector of matrices.
///
/// @param [in] filename
///    Path to result file to be read. Lines starting with '#' are ignored.
vector<Matrix> readTensorFile(const char* filename) {
    
    // Check if filename is readable
    ifstream resultfile(filename, ios::in);
    if (resultfile.fail()) {
        cerr << "Error: Filename " << filename << " not found or not readable." << endl;
        exit(1);
    }
    
    vector<Matrix> result;
    Matrix nextMatrix(0.0);
    string line;
    double nextNumber;
    int row = 0;
    // Read file line by line
    while (getline(resultfile, line)) {
	if (line[0] == '#') // Skip lines starting with #
            continue;
        istringstream iss(line);
	for (int col=0; col<3; ++col) {
	    if (iss >> nextNumber) 
		nextMatrix[row][col] = nextNumber;
	    else { // This means we are missing an entry
		cerr << "Error: Number of columns is less than three in result file " 
		     << filename << "!" << endl;
		exit(1);
	    }
	}
	if (iss >> nextNumber) { // This means we have too many entries
	    cerr << "Error: Number of columns is bigger than three in result file " 
		 << filename << "!" << endl;
	    exit(1);
	}
	++row;
	if (row == 2) { // nextMatrix is filled up with results
	    row = 0;
	    result.push_back(nextMatrix);
	    nextMatrix = 0.0;
	}
    }

    return result;
}
    

/// @brief
///    Tests if two tensor result files are equal
///
/// Command input variables;
///    1) Path to reference solution file
///    2) Path to new solution file to compare with
///    3) Tolerance
///
/// The results in the input files must be represented as tensors (like the output from upscale_perm)
/// and tensors are compared with the Frobenius norm within the given tolerance. The result files
/// can contain multiple tensors (which is the case when running upscale_perm with different BCs).
int main(int varnum, char** vararg) {

    // Check if the correct number of variables are given
    if (varnum != 4) {
        cout << "Error: Wrong number of input variables, should be three!" << endl
             << "Usage: ./compare_tensor_result refSolnFile newSolnFile tol" << endl;
        exit(1);
    }

    // Process input
    const char* refSolnFile(vararg[1]);
    const char* newSolnFile(vararg[2]);
    double tol = atof(vararg[3]);

    // Read result files
    vector<Matrix> refSoln = readTensorFile(refSolnFile);
    vector<Matrix> newSoln = readTensorFile(newSolnFile);

    // Test if solutions contains equally many tensors
    if (refSoln.size() != newSoln.size()) {
	cerr << "Error: The result files contain different number of tensors!" << endl;
	exit(1);
    }  

    // Compare each tensor to each other with Frobenius norm
    for (int iTensor = 0; iTensor < refSoln.size(); ++iTensor) {
	Matrix diffMatrix = refSoln[iTensor];
	diffMatrix -= newSoln[iTensor];
	if (diffMatrix.frobenius_norm() > tol) {
	    cerr << "Tensor #" << (iTensor + 1) << " in " << refSolnFile 
		 << " is not equal to reference solution within a tolerance of " 
		 << tol << "!" << endl;
	    exit(1);
	}
    }
    
    return 0;
}
    

