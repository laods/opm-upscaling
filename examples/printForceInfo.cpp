/* Arguments:
 *   1) eclipse file
 *   2) surface tension (sigma) in N/m
 *   3) Pressure drop in Pa
 */

#include <config.h>

#include <iostream>
#include <fstream>

#include <dune/grid/CpGrid.hpp>
#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/porsol/common/GridInterfaceEuler.hpp>
#include <opm/porsol/common/Rock.hpp>

using namespace std;

int main(int varnum, char** vararg) {
    
    // Read input
    const char* ECLIPSEFILENAME(vararg[1]);
    int dir = atoi(vararg[2]);
    double sigma = atof(vararg[3]);
    double deltap = atof(vararg[4]);
    const string section_name = string(vararg[5]);
     
    // Cap pressure data rocks and map poro and satnum
    // OBS!!! HARD-CODED!
    double Pc_avg[20] = { };
    double Pb[20] = { };
    map<double,int> map_poro_satnum;
    if (section_name == "SN74" || section_name == "sn74") {
	Pc_avg = {23771,86713,1987,50694,28470,66404,26280,58528,4547,79796,61050,10623};
	Pb = {3038,12617,679,3712,5526,8420,4949,9626,1440,14447,8060,3386};
	map_poro_satnum[0.1933] = 1;
	map_poro_satnum[0.1033] = 2;
	map_poro_satnum[0.27] = 3;
	map_poro_satnum[0.1755] = 4;
	map_poro_satnum[0.2249] = 5;
	map_poro_satnum[0.1532] = 6;
	map_poro_satnum[0.203] = 7;
	map_poro_satnum[0.1701] = 8;
	map_poro_satnum[0.2734] = 9;
	map_poro_satnum[0.0801] = 10;
	map_poro_satnum[0.12] = 11;
	map_poro_satnum[0.24] = 12;
    }
    else if (section_name == "SN82" || section_name == "sn82") {
	Pc_avg = {23771,86713,50694,20412,2882,31040,4547,79796,20369,10623};
	Pb = {3038,12617,3712,1630,1153,10224,1440,14447,1573,3386};
	map_poro_satnum[0.1933] = 1;
	map_poro_satnum[0.1033] = 2;
	map_poro_satnum[0.1755] = 3;
	map_poro_satnum[0.215] = 4;
	map_poro_satnum[0.2348] = 5;
	map_poro_satnum[0.2098] = 6;
	map_poro_satnum[0.2734] = 7;
	map_poro_satnum[0.0801] = 8;
	map_poro_satnum[0.2149] = 9;
	map_poro_satnum[0.24] = 10;
    }
    else if (section_name == "SN91" || section_name == "sn91") {
	Pc_avg = {23771,86713,1987,50694,78336,3516,9222,66404,6282,5675,20412,26280,4547,79796,7158};
	Pb = {3038,12617,679,3712,17123,1327,1869,8420,1188,1401,1630,4949,1440,14447,1775};
	map_poro_satnum[0.1933] = 1;
	map_poro_satnum[0.1033] = 2;
	map_poro_satnum[0.27] = 3;
	map_poro_satnum[0.1755] = 4;
	map_poro_satnum[0.15] = 5;
	map_poro_satnum[0.2299] = 6;
	map_poro_satnum[0.2534] = 7;
	map_poro_satnum[0.1532] = 8;
	map_poro_satnum[0.2428] = 9;
	map_poro_satnum[0.2098] = 10;
	map_poro_satnum[0.215] = 11;
	map_poro_satnum[0.203] = 12;
	map_poro_satnum[0.2734] = 13;
	map_poro_satnum[0.0801] = 14;
	map_poro_satnum[0.2604] = 15;
    }
    else if (section_name == "Parallel" || section_name == "parallel") {
	Pc_avg = {1.2286e4,2.7473e5};
	Pb = {4554.9,1.0185e5};
	map_poro_satnum[0.2] = 1;
	map_poro_satnum[0.1] = 2;
    }
    else {
	cerr << "Model type (argument 5: '" << section_name << "') unknown. Either SN74, SN82, SN91 or Parallel.";
	exit(1);
    }

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
    
    Opm::Rock<3> rock;
    rock.init(deck, grid.globalCell());
   
    double accumulated_porevol = 0.0;
    double nabla_Pc_J = 0.0;
    double nabla_Pc_avg = 0.0;
    double nabla_Pb = 0.0;
    double min_coord = 1e16;
    double max_coord = -1e16;
    
    typedef GridInterface::CellIterator CI;
    typedef CI::FaceIterator FI;
    for (CI cell_iter = gridinterf.cellbegin(); cell_iter != gridinterf.cellend(); ++cell_iter) {
        int cell_index = cell_iter->index();
        const double cell_perm = rock.permeability(cell_index)(dir,dir);
        const double cell_poro = rock.porosity(cell_index);
        const double cell_porevol = cell_iter->volume()*cell_poro;
        const double cell_lambda = sqrt(cell_perm/cell_poro);
	if (map_poro_satnum.find(cell_poro) == map_poro_satnum.end()) {
	    cerr << "Can't map poro " << cell_poro << " to a satnum. Check map!" << endl;
	    exit(1);
	}
	const int cell_satnum = (map_poro_satnum.find(cell_poro))->second;
        const double cell_Pc_avg = Pc_avg[cell_satnum];
        const double cell_Pb = Pb[cell_satnum];
        for (FI face_iter = cell_iter->facebegin(); face_iter != cell_iter->faceend(); ++face_iter) {
            if ( ! face_iter->boundary()) {
                int neighbour_index = face_iter->neighbourCellIndex();
                const double neighbour_perm = rock.permeability(neighbour_index)(dir,dir);
                const double neighbour_poro = rock.porosity(neighbour_index);
                const double neighbour_porevol = face_iter->neighbourCellVolume()*neighbour_poro;
                const double neighbour_lambda = sqrt(neighbour_perm/neighbour_poro);
		if (map_poro_satnum.find(neighbour_poro) == map_poro_satnum.end()) {
		    cerr << "Can't map poro " << neighbour_poro << " to a satnum. Check map!" << endl;
		    exit(1);
		}
		const int neighbour_satnum = (map_poro_satnum.find(neighbour_poro))->second;
                const double neighbour_Pc_avg = Pc_avg[neighbour_satnum];
                const double neighbour_Pb = Pb[neighbour_satnum];
                const double length = (cell_iter->centroid() - face_iter->neighbourCell().centroid()).two_norm();
                const double porevol = cell_porevol + neighbour_porevol;
                accumulated_porevol += porevol;
                nabla_Pc_J += abs(1/cell_lambda - 1/neighbour_lambda)*porevol/length;
                nabla_Pc_avg += abs(cell_Pc_avg - neighbour_Pc_avg)*porevol/length;
                nabla_Pb += abs(cell_Pb - neighbour_Pb)*porevol/length;
            }
            else {
                min_coord = min(min_coord, face_iter->centroid()[dir]);
                max_coord = max(max_coord, face_iter->centroid()[dir]);
            }
        }
    }
    nabla_Pc_J = sigma*nabla_Pc_J/accumulated_porevol;
    nabla_Pc_avg = nabla_Pc_avg/accumulated_porevol;
    nabla_Pb = nabla_Pb/accumulated_porevol;
    
    double model_length = max_coord - min_coord;
    double nabla_visc_pressure = deltap/model_length;
    
    cout << scientific;
    cout << endl;
    cout << "Pressure drop:        " << deltap << endl;
    cout << "||nabla p||           " << nabla_visc_pressure << "\t(or " << 1/model_length << "*DeltaP)" << endl;
    cout << endl;
    cout << "With |J|=1:" << endl;
    cout << "||nabla Pc||          " << nabla_Pc_J << endl;
    cout << "Ca:                   " << nabla_visc_pressure/nabla_Pc_J << "\t(or " << 1/(model_length*nabla_Pc_J) << "*DeltaP)" << endl;
    cout << "Ca = 1 when Delta p = " << nabla_Pc_J * model_length << endl;
    cout << endl;
    cout << "With Pc_avg:" << endl;
    cout << "||nabla Pc||          " << nabla_Pc_avg << endl;
    cout << "Ca:                   " << nabla_visc_pressure/nabla_Pc_avg << "\t(or " << 1/(model_length*nabla_Pc_avg) << "*DeltaP)" << endl;
    cout << "Ca = 1 when Delta p = " << nabla_Pc_avg * model_length << endl;
    cout << endl;
    cout << "With Pb:" << endl;
    cout << "||nabla Pc||          " << nabla_Pb << endl;
    cout << "Ca:                   " << nabla_visc_pressure/nabla_Pb << "\t(or " << 1/(model_length*nabla_Pb) << "*DeltaP)" << endl;
    cout << "Ca = 1 when Delta p = " << nabla_Pb * model_length << endl;
    cout << endl;
}
