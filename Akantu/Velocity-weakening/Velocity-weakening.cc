#include "solid_mechanics_model_cohesive.hh"
#include "mesh.hh"
#include "aka_common.hh"
#include <iostream>
#include <iomanip>
#include <chrono>


int main(int argc, char *argv[])
{
    constexpr akantu::Int sd = 3;
    const akantu::Real us = 1e-6;
    const akantu::Real ms = 1e-3;
    const std::string mesh_file = "../../../Models/100mm-PMMA-CZM.msh";
    const std::string mat_file = "../../../Materials/material-mm-MPa.dat";

    akantu::initialize(mat_file, argc, argv);
    std::cout << "Initialized" << std::endl;
    akantu::Mesh mesh(sd);
    mesh.read(mesh_file);
    std::cout << "Load files successful." << std::endl;
    std::cout << "Cells (3D): " << mesh.getNbElement(mesh.getSpatialDimension()) << std::endl;
    std::cout << "Faces (2D): " << mesh.getNbElement(mesh.getSpatialDimension() - 1) << std::endl;
    std::cout << "Edges (1D): " << mesh.getNbElement(1) << std::endl;

    akantu::SolidMechanicsModelCohesive model(mesh);

    akantu::MaterialCohesiveRules rules{
        {{"moving-block", "moving-block"}, "non_interface"},
        {{"stationary-block", "stationary-block"}, "non_interface"},
        {{"moving-block", "stationary-block"}, "interface"},
        {{"stationary-block", "moving-block"}, "interface"}
    };
    std::cout << "Got material" << std::endl;

    auto cohesive_selector = std::make_shared<akantu::MaterialCohesiveRulesSelector>(model, rules);
    auto bulk_selector = std::make_shared<akantu::MeshDataMaterialSelector<std::string>>("physical_names", model);
    std::cout << "Got physical names" << std::endl;

    // cohesive_selector->setFallback(bulk_selector);
    bulk_selector->setFallback(model.getMaterialSelector());
    model.setMaterialSelector(cohesive_selector);
    std::cout << "Set material selector" << std::endl;



    model.initFull(akantu::_analysis_method = akantu::_explicit_lumped_mass, akantu::_is_extrinsic = true);

    std::cout << "After model initialization" << std::endl;

    akantu::Real dt = model.getStableTimeStep() * 0.5;
    model.setTimeStep(dt);
    std::cout << "dt = " << dt << " s (" << dt / us << " us)\n";

    model.setBaseName("50mm-PMMA-CZM-velocity-weakening");
    model.assembleMassLumped();
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

    akantu::Vector<akantu::Real, 3> t_front{8.0, 0.0, 0.0}; // MPa traction (+X)
    akantu::Vector<akantu::Real, 3> t_left{ 0.0, 6.0, 0.0}; // MPa traction (+Y)

    model.applyBC(akantu::BC::Neumann::FromTraction(t_front), "moving-block-front");
    model.applyBC(akantu::BC::Neumann::FromTraction(t_left), "moving-block-left");

    model.applyBC(akantu::BC::Dirichlet::FixedValue(0., akantu::_y), "stationary-block-right");
    model.applyBC(akantu::BC::Dirichlet::FixedValue(0., akantu::_x), "stationary-block-back");

    std::cout << "set B.C. successful." << std::endl;

    const akantu::Real SIMULATION_TIME = 20.0 * ms;           // 20 ms = 0.02 s
    const akantu::Int max_steps = ceil(SIMULATION_TIME / dt); // total number of time steps
    std::cout << "Starting time integration for " << (SIMULATION_TIME / ms) << " ms (" << max_steps << " steps)\n";
    
    auto start_time = std::chrono::high_resolution_clock::now();

    for (akantu::Int s = 0; s < max_steps; ++s)
    {
        model.checkCohesiveStress();
        model.solveStep();
        if (s % 5 == 0) {
            model.dump();
        }
        auto current_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = current_time - start_time;

        double time_per_iter = elapsed.count() / s;
        double estimated_total = time_per_iter * max_steps;
        double remaining = estimated_total - elapsed.count();

        std::cout << "Step " << std::setw(8) << s << "/" << max_steps
                  << " | Elapsed: " << std::fixed << std::setprecision(1) << std::setw(8) << elapsed.count() << " s"
                  << " | ETA: " << std::fixed << std::setprecision(1) << std::setw(8) << remaining << " s\r";
        std::cout.flush();
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> total_elapsed = end_time - start_time;
    std::cout << "Total elapsed time: " << total_elapsed.count() << " s" << std::endl;

    akantu::finalize();
    return 0;
}