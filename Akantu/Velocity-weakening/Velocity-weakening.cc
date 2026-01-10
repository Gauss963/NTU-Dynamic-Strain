#include "ntn_friclaw_linear_slip_weakening.hh"
#include "solid_mechanics_model_cohesive.hh"
#include "aka_common.hh"
#include "mesh.hh"
#include <iostream>
#include <iomanip>
#include <chrono>


int main(int argc, char *argv[])
{
    constexpr akantu::Int sd = 3;
    constexpr akantu::Real us = 1e-6;
    constexpr akantu::Real ms = 1e-3;
    constexpr int PMMA_thickness = 50;
    constexpr akantu::Real time_factor = 0.5;
    const std::string mesh_file =  "../../../Models/" + std::to_string(PMMA_thickness) + "mm-PMMA-CZM.msh";
    const std::string mat_file = "../../../Materials/material-mm-MPa.dat";

    akantu::initialize(mat_file, argc, argv);
    std::cout << "Initialized" << std::endl;
    akantu::Mesh mesh(sd);

    // const auto &comm = akantu::Communicator::getStaticCommunicator();
    // akantu::Int prank = comm.whoAmI();
    // if (prank == 0) {
    //     mesh.read(mesh_file);
    // }
    // mesh.distribute();


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

    bulk_selector->setFallback(model.getMaterialSelector());
    model.setMaterialSelector(cohesive_selector);
    std::cout << "Set material selector" << std::endl;


    model.initFull(akantu::_analysis_method = akantu::_explicit_lumped_mass, akantu::_is_extrinsic = true);

    std::cout << "After model initialization" << std::endl;

    akantu::Real dt = model.getStableTimeStep() * time_factor;
    model.setTimeStep(dt);
    std::cout << "dt = " << dt << " s (" << dt / us << " us)\n";
    std::string SYS_PATH = "../../../../../../../../../../";
    // std::string DUMP_PATH = "Volumes/Gauss-T7/Dump/";
    std::string DUMP_PATH = "Volumes/Gauss-T7/Dump_5mm/";

    model.setBaseName(SYS_PATH + DUMP_PATH + std::to_string(PMMA_thickness) + "mm-PMMA-CZM-velocity-weakening");
    // model.setBaseName(std::to_string(PMMA_thickness) + "mm-PMMA-CZM-velocity-weakening");
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

    const akantu::Real normalStress  = 8.0;  // MPa
    const akantu::Real shearStress   = 3.0;  // MPa
    const akantu::Real displacement  = 3.0;  // mm
    akantu::Real pushValue           = 0.0;  // mm
    const akantu::Real riseEnd       = 0.70; // % of total time

    akantu::Vector<akantu::Real, 3> t_front{normalStress, 0.0, 0.0};  // MPa traction (+X)
    akantu::Vector<akantu::Real, 3> t_left{0.0, shearStress, 0.0};    // MPa traction (+Y)

    model.applyBC(akantu::BC::Neumann::FromTraction(t_front), "moving-block-front");
    // model.applyBC(akantu::BC::Neumann::FromTraction(t_left), "moving-block-left");

    model.applyBC(akantu::BC::Dirichlet::FixedValue(0., akantu::_y), "stationary-block-right");
    model.applyBC(akantu::BC::Dirichlet::FixedValue(0., akantu::_x), "stationary-block-back");

    std::cout << "set B.C. successful." << std::endl;

    // const akantu::Real SIMULATION_TIME = 20.0 * ms;
    const akantu::Real SIMULATION_TIME = 0.4 * ms;
    const akantu::Int MAX_STEPS = ceil(SIMULATION_TIME / dt);
    const akantu::Int TOTAL_FRAMES = 840;
    // const akantu::Int DUMP_INTERVAL = MAX_STEPS / TOTAL_FRAMES;
    const akantu::Int DUMP_INTERVAL = 3;

        std::cout << "Starting time integration for " << (SIMULATION_TIME / ms) << " ms (" << MAX_STEPS << " steps)\n";
    std::cout << "Dumping every " << DUMP_INTERVAL << " steps (" << TOTAL_FRAMES << " frames)" << std::endl;
    auto start_time = std::chrono::high_resolution_clock::now();

    for (akantu::Int s = 0; s < MAX_STEPS; ++s) {

        akantu::Real progress = static_cast<akantu::Real>(s) / MAX_STEPS;
        if (progress < riseEnd) {
            akantu::Real alpha = progress / riseEnd;
            pushValue = alpha * displacement;
        }
        else {
            pushValue = displacement;
        }
        model.applyBC(akantu::BC::Dirichlet::FixedValue(pushValue, akantu::_y), "moving-block-left");
        model.checkCohesiveStress();
        model.solveStep();
        if (s % DUMP_INTERVAL == 0) {
            model.dump();
        }


        auto current_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = current_time - start_time;

        double time_per_iter = elapsed.count() / s;
        double estimated_total = time_per_iter * MAX_STEPS;
        double remaining = estimated_total - elapsed.count();

        std::cout << "Step " << std::setw(8) << s << "/" << MAX_STEPS
                  << " | Elapsed: " << std::fixed << std::setprecision(1) << std::setw(8) << elapsed.count() << " s"
                  << " | ETA: " << std::fixed << std::setprecision(1) << std::setw(8) << remaining << " s = "
                  << std::fixed << std::setprecision(1) << std::setw(5) << remaining / 3600 << " h\r";
        std::cout.flush();
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> total_elapsed = end_time - start_time;
    std::cout << "Total elapsed time: " << total_elapsed.count() << " s = " << total_elapsed.count() / 3600 << " h" << std::endl;

    akantu::finalize();
    return 0;
}