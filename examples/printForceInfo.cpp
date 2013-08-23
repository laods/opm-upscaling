/* Arguments:
 *   1) eclipse file
 *   2) surface tension (sigma) in N/m
 *   3) Pressure drop in Pa
 */

#include <config.h>

#include <iostream>
#include <fstream>

#include <dune/grid/CpGrid.hpp>
#include <opm/core/io/eclipse/EclipseGridParser.hpp>
#include <opm/porsol/common/GridInterfaceEuler.hpp>
#include <opm/porsol/common/Rock.hpp>

using namespace std;

int main(int varnum, char** vararg) {
    
    // Read input
    const char* ECLIPSEFILENAME(vararg[1]);
    int dir = atoi(vararg[2]);
    double sigma = atof(vararg[3]);
    double deltap = atof(vararg[4]);
     
    // Cap pressure data for SN91 rocks.
    // OBS!!! HARD-CODED!
    //const double Pc_avg[] = {23771,86712,1987,50694,78336,3516,9222,66404,6282,5675,20412,26280,4547,79796,7158};
    //const double Pb[] = {3038,12617,679,3712,17123,1327,1869,8420,1188,1401,1630,4949,1440,14447,1775};
    
    // Cap pressure data for parallel rocks
    // OBS!!! HARD-CODED!
    const double Pc_avg[] = {1.2286e4,2.7473e5};
    const double Pb[] = {4554.9,1.0185e5};
    
    // Test if filename exists and is readable
    ifstream eclipsefile(ECLIPSEFILENAME, ios::in);
    if (eclipsefile.fail()) {
        cerr << "Error: Filename " << ECLIPSEFILENAME << " not found or not readable." << endl;
        exit(1);
    }
    eclipsefile.close();
    
    Opm::EclipseGridParser * eclParser_p;
    try {
        eclParser_p = new Opm::EclipseGridParser(ECLIPSEFILENAME);
    }
    catch (...) {
        cout << "Error: Filename " << ECLIPSEFILENAME << " does not look like an eclipse grid file." << endl;
        exit(1);
    }
    Opm::EclipseGridParser& eclParser = *eclParser_p;

    Dune::CpGrid grid;
    grid.processEclipseFormat(eclParser, 0, false);
    
    typedef Opm::GridInterfaceEuler<Dune::CpGrid> GridInterface;
    GridInterface gridinterf(grid);
    
    Opm::Rock<3> rock;
    rock.init(eclParser, grid.globalCell());
    
    // OBS!!! HARD-CODED!
    vector<int> cell2satnum(grid.numCells(),0);
    for (GridInterface::CellIterator cell_iter = gridinterf.cellbegin(); cell_iter != gridinterf.cellend(); ++cell_iter) {
        int cell_index = cell_iter->index();
        const double cell_poro = rock.porosity(cell_index);
        /*
        if (cell_poro == 0.1933) {
            cell2satnum[cell_index] = 1;
        }
        else if (cell_poro == 0.1033) {
            cell2satnum[cell_index] = 2;
        }
        else if (cell_poro == 0.27) {
            cell2satnum[cell_index] = 3;
        }
        else if (cell_poro == 0.1755) {
            cell2satnum[cell_index] = 4;
        }
        else if (cell_poro == 0.15) {
            cell2satnum[cell_index] = 5;
        }
        else if (cell_poro == 0.2299) {
            cell2satnum[cell_index] = 6;
        }
        else if (cell_poro == 0.2534) {
            cell2satnum[cell_index] = 7;
        }
        else if (cell_poro == 0.1532) {
            cell2satnum[cell_index] = 8;
        }
        else if (cell_poro == 0.2428) {
            cell2satnum[cell_index] = 9;
        }
        else if (cell_poro == 0.2098) {
            cell2satnum[cell_index] = 10;
        }
        else if (cell_poro == 0.215) {
            cell2satnum[cell_index] = 11;
        }
        else if (cell_poro == 0.203) {
            cell2satnum[cell_index] = 12;
        }
        else if (cell_poro == 0.2734) {
            cell2satnum[cell_index] = 13;
        }
        else if (cell_poro == 0.0801) {
            cell2satnum[cell_index] = 14;
        }
        else if (cell_poro == 0.2604) {
            cell2satnum[cell_index] = 15;
        }
        else {
            cout << "Unknown porosity (" << cell_poro << ") at cell " << cell_index << endl;
            exit(1);
        }
        */
        if (cell_poro == 0.2) {
            cell2satnum[cell_index] = 1;
        }
        else if (cell_poro == 0.1) {
            cell2satnum[cell_index] = 2;
        }
        else {
            cout << "Unknown porosity (" << cell_poro << ") at cell " << cell_index << endl;
            exit(1);
        }
    }
    
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
        const int cell_satnum = cell2satnum[cell_index];
        const double cell_Pc_avg = Pc_avg[cell_satnum];
        const double cell_Pb = Pb[cell_satnum];
        for (FI face_iter = cell_iter->facebegin(); face_iter != cell_iter->faceend(); ++face_iter) {
            if ( ! face_iter->boundary()) {
                int neighbour_index = face_iter->neighbourCellIndex();
                const double neighbour_perm = rock.permeability(neighbour_index)(dir,dir);
                const double neighbour_poro = rock.porosity(neighbour_index);
                const double neighbour_porevol = face_iter->neighbourCellVolume()*neighbour_poro;
                const double neighbour_lambda = sqrt(neighbour_perm/neighbour_poro);
                const int neighbour_satnum = cell2satnum[neighbour_index];
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