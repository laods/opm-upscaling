#include <config.h>

#include <iostream>
#include <fstream>

#include <dune/grid/CpGrid.hpp>
//#include <opm/core/io/eclipse/EclipseGridParser.hpp>
#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/porsol/common/GridInterfaceEuler.hpp>
#include <opm/porsol/common/Rock.hpp>

using namespace std;


int main(int varnum, char** vararg) {

    // Read input
    const char* ECLIPSEFILENAME(vararg[1]);
    int dir = atoi(vararg[2]);
    double deltap = atoi(vararg[3]);

    // Test if filename exists and is readable
    ifstream eclipsefile(ECLIPSEFILENAME, ios::in);
    if (eclipsefile.fail()) {
        cerr << "Error: Filename " << ECLIPSEFILENAME << " not found or not readable." << endl;
        exit(1);
    }
    eclipsefile.close();
    
    // Parse grdecl file
    std::cout << "Parsing grid file '" << ECLIPSEFILENAME << "' ..." << std::endl;
    Opm::ParserPtr parser(new Opm::Parser);
    Opm::DeckConstPtr deck(parser->parseFile(ECLIPSEFILENAME));

    Dune::CpGrid grid;
    grid.processEclipseFormat(deck, 0, false);
    
    typedef Opm::GridInterfaceEuler<Dune::CpGrid> GridInterface;
    GridInterface gridinterf(grid);

    // Calculate average cell size
    typedef GridInterface::Vector Vector;
    typedef GridInterface::CellIterator CI;
    typedef CI::FaceIterator FI;
    int count = 0;
    double l  = 0.0;
    double min_coord = 1e10;
    double max_coord = -1e10;
    for (CI ci = gridinterf.cellbegin(); ci != gridinterf.cellend(); ++ci) {
        const Vector cellCenter = ci->centroid();
	for (FI fi = ci->facebegin(); fi != ci->faceend(); ++fi) {
	    if ( ! fi->boundary()) {
		Vector neighbourVector = fi->neighbourCell().centroid();
		neighbourVector -= cellCenter;
		const double distance = neighbourVector.two_norm();
		l += distance;
		++count;
	    }
	    else {
		min_coord = min(min_coord, fi->centroid()[dir]);
                max_coord = max(max_coord, fi->centroid()[dir]);
	    }
	}
    }
    l = l/count;


}
		
		
