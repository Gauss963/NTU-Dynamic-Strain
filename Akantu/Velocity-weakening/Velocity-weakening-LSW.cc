// #include "aka_common.hh"
#include "coupler_solid_contact.hh"
// #include "solid_mechanics_model.hh"
#include "mesh.hh"

// #include "ntn_contact.hh"
// #include "ntn_friclaw_linear_slip_weakening.hh"
// #include "ntn_fricreg_no_regularisation.hh"
// #include "ntn_friction.hh"

#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>

int main(int argc, char *argv[]) {
    constexpr akantu::Int sd = 3;
    constexpr akantu::Real us = 1e-6;
    constexpr akantu::Real ms = 1e-3;
    constexpr akantu::Real TIME_FACTOR = 0.0005;
    constexpr akantu::Real SIMULATION_TIME = 40 * us;
    constexpr int PMMA_thickness = 50;

    const std::string mesh_file = "../../../Models/" + std::to_string(PMMA_thickness) + "mm-BS-PMMA.msh";
    const std::string mat_file = "../../../Materials/NTN-LSW.dat";

    akantu::initialize(mat_file, argc, argv);

    akantu::Mesh mesh(sd);
    const auto &comm = akantu::Communicator::getStaticCommunicator();
    const akantu::Int prank = comm.whoAmI();

    if (prank == 0) {
        mesh.read(mesh_file);
    }

    mesh.distribute();

    akantu::CouplerSolidContact coupler(mesh);

    auto &solid = coupler.getSolidMechanicsModel();
    auto &contact = coupler.getContactMechanicsModel();
    auto &&selector = std::make_shared<akantu::MeshDataMaterialSelector<std::string>>("physical_names", solid);
    solid.setMaterialSelector(selector);

    coupler.initFull(akantu::_analysis_method = akantu::_explicit_lumped_mass);

    akantu::Real dt = coupler.getStableTimeStep() * TIME_FACTOR;
    if (prank == 0) {
        std::cout << "[SIM] dt = " << dt << " s (" << dt / us << " us)\n";
    }
    coupler.setTimeStep(dt);

    auto &&surface_selector = std::make_shared<akantu::PhysicalSurfaceSelector>(mesh);
    contact.getContactDetector().setSurfaceSelector(surface_selector);

    coupler.setBaseName(std::to_string(PMMA_thickness) + "mm-PMMA-NTN-LSW");
    coupler.addDumpFieldVector("displacement");
    coupler.addDumpFieldVector("velocity");
    coupler.addDumpFieldVector("internal_force");
    coupler.addDumpFieldVector("external_force");
    coupler.addDumpField("stress");
    coupler.addDumpField("grad_u");
    coupler.addDumpFieldVector("contact_force");
    coupler.addDumpFieldVector("normals");
    coupler.addDumpField("areas");

    coupler.getVelocity().set(0.0);
    coupler.getDisplacement().set(0.0);

    constexpr akantu::Real normalStress = 8.0; // MPa
    constexpr akantu::Real shearDisp = 3.0;    // mm
    constexpr akantu::Real riseEnd = 0.30;
    constexpr akantu::Int TOTAL_FRAMES = 2400;

    const akantu::Int SHEAR_STEPS = std::ceil(SIMULATION_TIME / dt);
    const akantu::Int PRESS_STEPS = 0.05 * SHEAR_STEPS;
    const akantu::Int DUMP_INTERVAL = std::max<akantu::Int>(1, SHEAR_STEPS / TOTAL_FRAMES);
    const akantu::Int rise_steps = std::ceil(riseEnd * SHEAR_STEPS);
    const akantu::Real dy = shearDisp / static_cast<akantu::Real>(rise_steps);

    // akantu::Vector<akantu::Real, 3> t_normal{normalStress, 0.0, 0.0};
    // solid.applyBC(akantu::BC::Neumann::FromTraction(t_normal), "moving-block-back");
    solid.applyBC(akantu::BC::Dirichlet::FixedValue(0., akantu::_x), "stationary-block-front");
    solid.applyBC(akantu::BC::Dirichlet::FixedValue(0., akantu::_y), "stationary-block-left");

    coupler.dump();

    auto start_time = std::chrono::high_resolution_clock::now();
    if (prank == 0) {
        std::cout << "[SIM] Starting time integration for compression phase for "
                  << PRESS_STEPS
                  << " steps\n";
    }

    // // Compression step
    // akantu::Vector<akantu::Real, 3> t_normal{normalStress, 0.0, 0.0};
    // solid.applyBC(akantu::BC::Neumann::FromTraction(t_normal), "moving-block-back");
    // for (akantu::Int step_now = 0; step_now < PRESS_STEPS; step_now++) {
    //     coupler.solveStep();
    //     if (step_now % DUMP_INTERVAL == 0) {
    //         coupler.dump();
    //     }
    // }

    // Compression step (DISPLACEMENT CONTROL)
    // Goal: close initial gap (~1 mm) and add a small overclosure to guarantee contact.
    constexpr akantu::Real initial_gap = 1.0; // mm (must match your gmsh initial_offdet)
    constexpr akantu::Real overclosure = 0.2; // mm, small extra to ensure contact is engaged
    const akantu::Real total_dx = initial_gap + overclosure;
    const akantu::Real dx_step = total_dx / static_cast<akantu::Real>(PRESS_STEPS);
    // IMPORTANT: do NOT apply the Neumann traction here
    // akantu::Vector<akantu::Real, 3> t_normal{normalStress, 0.0, 0.0};
    // solid.applyBC(akantu::BC::Neumann::FromTraction(t_normal), "moving-block-back");
    for (akantu::Int step_now = 0; step_now < PRESS_STEPS; ++step_now) {
        solid.applyBC(akantu::BC::Dirichlet::IncrementValue(dx_step, akantu::_x), "moving-block-back");
        coupler.solveStep();
        if (step_now % DUMP_INTERVAL == 0) {
            coupler.dump();
        }
    }

    if (prank == 0) {
        std::cout << "[SIM] Starting time integration for shear phase for "
                  << (SIMULATION_TIME / ms) << " ms (" << SHEAR_STEPS
                  << " steps)\n";
    }

    // Shear step
    for (akantu::Int step_now = 0; step_now < SHEAR_STEPS; ++step_now) {
        solid.applyBC(akantu::BC::Dirichlet::IncrementValue(dy, akantu::_y), "moving-block-right");

        coupler.solveStep();

        if (step_now % DUMP_INTERVAL == 0) {
            coupler.dump();
        }
        auto now = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = now - start_time;

        const double denom = (step_now > 0) ? static_cast<double>(step_now) : 1.0;
        const double time_per_iter = elapsed.count() / denom;
        const double remaining = time_per_iter * SHEAR_STEPS - elapsed.count();

        double percent = 100.0 * static_cast<double>(step_now) / static_cast<double>(SHEAR_STEPS);
        if (prank == 0 && step_now % 100 == 0) {
            std::cout << "[SIM] Step " << std::setw(8) << step_now << "/" << SHEAR_STEPS
                      << " | Elapsed: " << std::fixed << std::setprecision(1)
                      << std::setw(8) << elapsed.count() << " s"
                      << " | Progress: " << std::setw(6)
                      << std::fixed << std::setprecision(1) << percent << " %"
                      << " | ETA: " << std::fixed << std::setprecision(1)
                      << std::setw(8) << std::max(0.0, remaining) << " s\n";
        }
    }
}