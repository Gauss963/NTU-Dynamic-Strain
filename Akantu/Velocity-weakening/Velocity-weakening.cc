#include "solid_mechanics_model_cohesive.hh"
#include "mesh.hh"
#include "aka_common.hh"
#include <iostream>
#include <chrono>

using namespace akantu;

int main(int argc, char *argv[])
{
    constexpr Int sd = 3;
    const std::string mesh_file = "../../../Models/50mm-PMMA-CZM.msh";
    const std::string mat_file = "../../../Materials/material-mm-MPa.dat";

    initialize(mat_file, argc, argv);
    std::cout << "Initialized" << std::endl;

    Mesh mesh(sd);
    mesh.read(mesh_file);
    std::cout << "Load files successful." << std::endl;

    std::cout << "Cells (3D): " << mesh.getNbElement(mesh.getSpatialDimension()) << std::endl;
    std::cout << "Faces (2D): " << mesh.getNbElement(mesh.getSpatialDimension() - 1) << std::endl;
    std::cout << "Edges (1D): " << mesh.getNbElement(1) << std::endl;

    SolidMechanicsModelCohesive model(mesh);

    MaterialCohesiveRules rules{
        {{"friction_master", "friction_slave"}, "interface_mat"},
        {{"friction_slave", "friction_slave"}, "interface_mat"},
        {{"friction_master", "friction_master"}, "interface_mat"}};
    std::cout << "Got material" << std::endl;

    auto cohesive_selector = std::make_shared<MaterialCohesiveRulesSelector>(model, rules);
    auto bulk_selector = std::make_shared<MeshDataMaterialSelector<std::string>>("physical_names", model);
    std::cout << "Got physical names" << std::endl;

    cohesive_selector->setFallback(bulk_selector);
    bulk_selector->setFallback(model.getMaterialSelector());
    model.setMaterialSelector(cohesive_selector);
    std::cout << "Set material selector" << std::endl;

    
    // model.initFull(_analysis_method = _explicit_lumped_mass, _is_extrinsic = true);
    model.initFull(_analysis_method = _explicit_lumped_mass, _is_extrinsic = false);
    // This is where went WRONG if I use "_is_extrinsic = true" !!

    std::cout << "After model initialization" << std::endl;

    Real dt = model.getStableTimeStep() * 0.5;
    model.setTimeStep(dt);
    std::cout << "dt = " << dt << std::endl;

    model.assembleMassLumped();
    // model.setBaseName("czm_debug");
    model.addDumpFieldVector("displacement");
    model.addDumpFieldVector("velocity");
    model.addDumpFieldVector("internal_force");
    model.addDumpField("stress");
    // model.addDumpField("material_index");
    model.addDumpField("grad_u");
    // model.addDumpField("cohesive_opening");
    // model.addDumpField("cohesive_traction");

    std::cout << "Before getVelocity" << std::endl;
    auto &vel = model.getVelocity();
    std::cout << "Before getDisplacement" << std::endl;
    auto &disp = model.getDisplacement();
    std::cout << "Before set velocity to 0" << std::endl;
    vel.set(0.);
    disp.set(0.);
    std::cout << "After setting vel and disp" << std::endl;

    Vector<Real, 3> t_front{8.0, 0.0, 0.0}; // MPa traction (+X)
    Vector<Real, 3> t_left{ 0.0, 6.0, 0.0}; // MPa traction (+Y)

    model.applyBC(BC::Neumann::FromTraction(t_front), "moving-block-front");
    model.applyBC(BC::Neumann::FromTraction(t_left), "moving-block-left");

    model.applyBC(BC::Dirichlet::FixedValue(0., _y), "stationary-block-right");
    model.applyBC(BC::Dirichlet::FixedValue(0., _x), "stationary-block-back");

    std::cout << "set B.C. successful." << std::endl;

    const Int SIMULATION_TIME = 10;             // total simulation time in ms
    const Int max_steps = SIMULATION_TIME / dt; // total number of time steps
    std::cout << "Starting time integration for " << SIMULATION_TIME
              << " ms (" << max_steps << " steps)" << std::endl;
    auto start_time = std::chrono::high_resolution_clock::now();

    for (Int s = 0; s < max_steps; ++s)
    {
        model.solveStep();
        if (s % 100 == 0) {
            model.dump();
        }
        auto current_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = current_time - start_time;

        double time_per_iter = elapsed.count() / s;
        double estimated_total = time_per_iter * max_steps;
        double remaining = estimated_total - elapsed.count();

        std::cout << "Step " << s << "/" << max_steps
                  << " | Elapsed: " << elapsed.count() << " s"
                  << " | ETA: " << remaining << " s" << std::endl;
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> total_elapsed = end_time - start_time;
    std::cout << "Total elapsed time: " << total_elapsed.count() << " s" << std::endl;

    finalize();
    return 0;
}