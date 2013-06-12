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
///    Class to store upscaling results
class UpscalingResult {

public:
    UpscalingResult() : satPoints(0) {}
    UpscalingResult(const char* filename);
    vector<double> getCol1() {return col1;}
    vector<double> getCol2() {return col2;}
    Matrix getTensor(int i) {return tensors[i];}
    int getSatPoints() {return satPoints;}

private:
    vector<double> col1; // Typically capillary pressure or fractional flow
    vector<double> col2; // Typically saturation
    vector<Matrix> tensors; // E.g. rel.perm.
    int satPoints; // Number of saturation points in solution
    vector<pair<int,int> > vec2matIndex;
};

/// @brief
///    Construct object from result file
UpscalingResult::UpscalingResult(const char* filename) {

    // Check if filename is readable
    ifstream resultfile(filename, ios::in);
    if (resultfile.fail()) {
        cerr << "Error: Filename " << filename << " not found or not readable." << endl;
        exit(1);
    }

    // Create vector with correct matrix indexes
    vec2matIndex.push_back(pair<int,int>(0,0));
    vec2matIndex.push_back(pair<int,int>(1,1));
    vec2matIndex.push_back(pair<int,int>(2,2));
    vec2matIndex.push_back(pair<int,int>(1,2));
    vec2matIndex.push_back(pair<int,int>(0,2));
    vec2matIndex.push_back(pair<int,int>(0,1));
    vec2matIndex.push_back(pair<int,int>(2,1));
    vec2matIndex.push_back(pair<int,int>(2,0));
    vec2matIndex.push_back(pair<int,int>(1,0));

    Matrix nextMatrix(0.0);
    string line;
    double nextNumber;
    satPoints = 0;
    // Read file line by line
    while (getline(resultfile, line)) {
	if (line[0] == '#') // Skip lines starting with #
            continue;
        istringstream iss(line);
	// Get first row entry
	if (iss >> nextNumber)
	    col1.push_back(nextNumber);
	else {
	    cerr << "Error: Filename " << filename << " has to few columns!" << endl;
	    exit(1);
	}
	// Get second row entry
	if (iss >> nextNumber)
	    col2.push_back(nextNumber);
	else {
	    cerr << "Error: Filename " << filename << " has to few columns!" << endl;
	    exit(1);
	}
	// Read tensor
	int entry = 0;
	while (iss >> nextNumber) {
	    nextMatrix[vec2matIndex[entry].first, vec2matIndex[entry].second] = nextNumber;
	    ++entry;
	}
	if (entry != 3 || entry != 9) {
	    cerr << "Error: Filename " << filename << " has wrong number of tensor columns. " 
		 << "Should be either 3 or 9." << endl;
	    exit(1);
	}
	tensors.push_back(nextMatrix);
	nextMatrix=0.0;
	++satPoints;
    }
}


/// @brief
///    Calculate the 2-norm of a vector  
double two_norm(vector<double> vec) {
    double norm = 0;
    for (vector<double>::iterator it = vec.begin(); it != vec.end(); ++it)
	norm += *it * *it;
    return sqrt(norm);
}


/// @brief
///    Vector difference
vector<double> operator-(vector<double> vec1, vector<double> vec2) {
    vector<double> diffVec;
    int i=0;
    for (vector<double>::iterator it = vec1.begin(); it != vec1.end(); ++it) {
	diffVec.push_back(*it - vec2[i]);
	++i;
    }
    return diffVec;
}
  

/// @brief
///    Tests if two result files are equal
///
/// Command input variables:
///    1) Path to reference solution file
///    2) Path to new solution file to compare with
///    3) Tolerance
///
/// The results in the input files must be in list format, that is, one line for each saturation point.
/// The two first columns of the files are treated as vectors and compared to each other with the 2-norm.
/// The remaining entries on a line is treated as a tensor and compared to each other with the Frobenius
/// norm. Each line is supposed to contain either 3 or 9 tensor entries in voigt ordering. If only three
/// entries specified, these are assumed to be the diagonal entries and off-diagonal entries are set to 0.
int main(int varnum, char** vararg) {

    // Check if the correct number of variables are given
    if (varnum != 4) {
        cout << "Error: Wrong number of input variables, should be three!" << endl
             << "Usage: ./compare_list_result refSolnFile newSolnFile tol" << endl;
        exit(1);
    }

    // Process input
    const char* refSolnFile(vararg[1]);
    const char* newSolnFile(vararg[2]);
    double tol = atof(vararg[3]);

    // Read result files
    UpscalingResult refSoln(refSolnFile);
    UpscalingResult newSoln(newSolnFile);

    // Check number of points
    if (refSoln.getSatPoints() != newSoln.getSatPoints()) {
	cerr << "Error: The number of saturation points for the two result files are not equal!" 
	     << endl;
	exit(1);
    }
    // Check col1
    if (two_norm(refSoln.getCol1() - newSoln.getCol1()) > tol) {
	cerr << "Error: The first column of " << refSolnFile
	     << " is not equal to reference solution within a tolerance of " << tol << "!" << endl;

	exit(1);
    }
    // Check col2
    if (two_norm(refSoln.getCol2() - newSoln.getCol2()) > tol) {
	cerr << "Error: The second column of " << refSolnFile 
	     << " is not equal to reference solution within a tolerance of " << tol << "!" << endl;
	exit(1);
    }
    // Check tensors
    for (int satPoint = 0; satPoint < refSoln.getSatPoints(); ++satPoint) {
	Matrix diffMatrix = refSoln.getTensor(satPoint);
	diffMatrix -= newSoln.getTensor(satPoint);
	if (diffMatrix.frobenius_norm() > tol) {
	    cerr << "Tensor #" << (satPoint + 1) << " in " << refSolnFile 
		 << " is not equal to reference solution within a tolerance of " 
		 << tol << "!" << endl;
	    exit(1);
	}
    }

    return 0;
}
