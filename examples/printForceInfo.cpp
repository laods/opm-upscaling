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
    
    double accumulated_porevol = 0.0;
    double nabla_cap_pressure = 0.0;
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
        for (FI face_iter = cell_iter->facebegin(); face_iter != cell_iter->faceend(); ++face_iter) {
            if ( ! face_iter->boundary()) {
                int neighbour_index = face_iter->neighbourCellIndex();
                const double neighbour_perm = rock.permeability(neighbour_index)(dir,dir);
                const double neighbour_poro = rock.porosity(neighbour_index);
                const double neighbour_porevol = face_iter->neighbourCellVolume()*neighbour_poro;
                const double neighbour_lambda = sqrt(neighbour_perm/neighbour_poro);
                const double length = (cell_iter->centroid() - face_iter->neighbourCell().centroid()).two_norm();
                const double porevol = cell_porevol + neighbour_porevol;
                accumulated_porevol += porevol;
                nabla_cap_pressure += abs(1/cell_lambda - 1/neighbour_lambda)*porevol/length;
            }
            else {
                min_coord = min(min_coord, face_iter->centroid()[dir]);
                max_coord = max(max_coord, face_iter->centroid()[dir]);
            }
        }
    }
    nabla_cap_pressure = sigma*nabla_cap_pressure/accumulated_porevol;
    
    double model_length = max_coord - min_coord;
    
    double nabla_visc_pressure = deltap/model_length;
    
    cout << scientific;
    cout << "Pressure drop:                           " << deltap << endl;
    cout << "Characteristic size of capillary force:  " << nabla_cap_pressure << endl;
    cout << "Characteristic size of viscous force:    " << nabla_visc_pressure << "\t(or " << 1/model_length << "*DeltaP)" << endl;
    cout << "Capillary number (visc/cap):             " << nabla_visc_pressure/nabla_cap_pressure << "\t(or " << 1/(model_length*nabla_cap_pressure) << "*DeltaP)" << endl;
    cout << "Ca = 1 when pressure drop equals         " << nabla_cap_pressure * model_length << endl;
}