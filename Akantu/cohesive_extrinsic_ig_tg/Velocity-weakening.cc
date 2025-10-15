#include "solid_mechanics_model.hh"
#include "material_elastic.hh"
#include <iostream>
#include <omp.h>

using namespace akantu;

int main() {
    omp_set_num_threads(8);
    
    std::cout << "Using " << omp_get_max_threads() << " OpenMP threads.\n";



    const Int spatial_dimension = 2;
    const Int max_steps = 100;

    Mesh mesh(spatial_dimension);
    mesh.read("../../Models/50mm-PMMA.msh");

    if (mesh.getNbNodes() == 0) {
        std::cerr << "ERROR: Mesh has 0 nodes. Check if 'block.msh' exists in your directory.\n";
        return 1;
    }
    std::cout << "Mesh loaded: " << mesh.getNbNodes() << " nodes.\n";


    SolidMechanicsModel model(mesh);
    model.initFull(_analysis_method = _explicit_lumped_mass);

    return 0;
}