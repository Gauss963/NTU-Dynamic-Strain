#include "solid_mechanics_model.hh"
#include "material_elastic.hh"
#include <iostream>
#include <omp.h>

int main() {
    omp_set_num_threads(8);
    
    std::cout << "Using " << omp_get_max_threads() << " OpenMP threads.\n";
    
    const akantu::Int spatial_dimension = 3;
    const akantu::Int max_steps = 100;
    const std::string MESHFILE = "../../../Models/50mm-PMMA.msh";


    akantu::Mesh mesh(spatial_dimension);
    std::cout << "Reading mesh from: " << MESHFILE << std::endl;
    mesh.read(MESHFILE);

    if (mesh.getNbNodes() == 0) {
        std::cerr << "ERROR: Mesh has 0 nodes. Check if " << MESHFILE << " exists in your directory.\n";
        return 1;
    }
    std::cout << "Mesh loaded: " << mesh.getNbNodes() << " nodes.\n";


    akantu::SolidMechanicsModel model(mesh);
    model.initFull(akantu::_analysis_method = akantu::_explicit_lumped_mass);

    return 0;
}