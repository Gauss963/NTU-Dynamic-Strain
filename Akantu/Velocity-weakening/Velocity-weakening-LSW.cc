#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <unordered_set>
#include <vector>

#include "aka_common.hh"
#include "mesh.hh"
#include "solid_mechanics_model.hh"

#include "ntn_base_contact.hh"
#include "ntn_contact_solvercallback.hh"
#include "ntn_initiation_function.hh"

#include "ntn_friclaw_linear_slip_weakening.hh"
#include "ntn_fricreg_no_regularisation.hh"
#include "ntn_friction.hh"

// #include "faceRebuild.hh"

int main(int argc, char *argv[]) {
    constexpr akantu::Int sd = 3;
    // constexpr akantu::Real us = 1e-6;
    constexpr akantu::Real ms = 1e-3;

    constexpr akantu::Real TIME_FACTOR = 0.0020;
    constexpr akantu::Real SIMULATION_TIME = 0.2 * ms;
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

    akantu::SolidMechanicsModel model(mesh);

    const std::string slave_surface = "stationary-block-back";
    const std::string master_surface = "moving-block-front";
    constexpr akantu::Int normal_dir = 0;

    auto solver_ntn = std::make_unique<akantu::NTNContactSolverCallback>(model, slave_surface, master_surface, normal_dir, TIME_FACTOR);
    auto contact = solver_ntn->getContact();
    contact->initParallel();
    
    auto friction = solver_ntn->getFriction();
    auto &msh = model.getFEEngine().getMesh();

    contact->addSurfacePair(slave_surface, master_surface, normal_dir); // 20260203 Test NaN issue
    if (prank == 0) {
        std::cout << "[DBG] dt = " << model.getTimeStep() << std::endl;
        std::cout << "[DBG] nb_contact_nodes = " << contact->getNbContactNodes() << "\n";
        std::cout << "[DBG] masters size = " << contact->getMasters().size()
                  << " slaves size = " << contact->getSlaves().size() << "\n";
    }
    contact->updateNormals();
    contact->updateLumpedBoundary();
    contact->updateImpedance();

    const auto &s_nodes_aka = msh.getElementGroup(slave_surface).getNodeGroup().getNodes();
    const auto &m_nodes_aka = msh.getElementGroup(master_surface).getNodeGroup().getNodes();

    if (prank == 0) {
        std::cout << "[DBG] slave group nodes = " << s_nodes_aka.size()
                  << " | master group nodes = " << m_nodes_aka.size() << "\n";
    }

    std::vector<akantu::Idx> s_nodes(s_nodes_aka.begin(), s_nodes_aka.end());
    std::vector<akantu::Idx> m_nodes(m_nodes_aka.begin(), m_nodes_aka.end());
    std::sort(s_nodes.begin(), s_nodes.end());
    std::sort(m_nodes.begin(), m_nodes.end());

    std::vector<akantu::Idx> inter;
    inter.reserve(std::min(s_nodes.size(), m_nodes.size()));
    std::set_intersection(s_nodes.begin(), s_nodes.end(), m_nodes.begin(), m_nodes.end(), std::back_inserter(inter));

    if (prank == 0) {
        std::cout << "[DBG] node-id intersection = " << inter.size() << "\n";
        if (!inter.empty()) {
            std::cout << "[DBG] first shared node id = " << inter.front() << "\n";
        }
    }
    // ---------------------------------------------------------------

    if (prank == 0) {
        std::cout << "[DBG] nb_contact_nodes = " << contact->getNbContactNodes() << "\n";
        std::cout << "[DBG] masters size = " << contact->getMasters().size()
                  << " slaves size = " << contact->getSlaves().size() << "\n";
    }

    // contact->updateNormals();
    // contact->updateLumpedBoundary();
    // contact->updateImpedance();

    // akantu::Real dt = model.getStableTimeStep() * TIME_FACTOR;
    // if (prank == 0) {
    //     std::cout << "[SIM] dt = " << dt << " s (" << dt / us << " us)\n";
    // }
    // model.setTimeStep(dt);
    // constexpr akantu::Real dt = 1.0e3 * us;
    // akantu::Real dt = model.getStableTimeStep();
    akantu::Real dt = model.getTimeStep();

    model.setBaseName(std::to_string(PMMA_thickness) + "mm-PMMA-NTN-LSW");
    model.addDumpField("stress");
    model.addDumpField("mass");
    model.addDumpFieldVector("strain");
    model.addDumpFieldVector("displacement");
    model.addDumpFieldVector("internal_force");
    model.addDumpFieldVector("external_force");
    model.getVelocity().set(0.0);
    model.getDisplacement().set(0.0);

    constexpr akantu::Real shearDisp = 3.0;
    constexpr akantu::Real riseEnd = 0.30;
    constexpr akantu::Int TOTAL_FRAMES = 2400;
    constexpr akantu::Real initial_gap = 0.05;
    constexpr akantu::Real overclosure = 0.01;

    const akantu::Int SHEAR_STEPS = std::ceil(SIMULATION_TIME / dt);
    const akantu::Int PRESS_STEPS = SHEAR_STEPS;
    const akantu::Real total_dx = initial_gap + overclosure;
    const akantu::Real dx_step = total_dx / static_cast<akantu::Real>(PRESS_STEPS);

    const akantu::Int DUMP_INTERVAL = std::max<akantu::Int>(1, SHEAR_STEPS / TOTAL_FRAMES);
    const akantu::Int rise_steps = std::ceil(riseEnd * SHEAR_STEPS);
    const akantu::Real dy = shearDisp / static_cast<akantu::Real>(rise_steps);

    model.applyBC(akantu::BC::Dirichlet::FixedValue(0., akantu::_x), "stationary-block-front");
    model.applyBC(akantu::BC::Dirichlet::FixedValue(0., akantu::_y), "stationary-block-left");

    model.dump();

    //
    constexpr akantu::Real NORMAL_STRESS = 16.0;
    akantu::Vector<akantu::Real, 3> normal_traction{NORMAL_STRESS, 0.0, 0.0};
    model.applyBC(akantu::BC::Neumann::FromTraction(normal_traction), "moving-block-back");
    //

    auto start_time = std::chrono::high_resolution_clock::now();
    if (prank == 0) {
        std::cout << "[SIM] Starting time integration for compression phase for "
                  << PRESS_STEPS << " steps\n";
    }

    auto &displacement = model.getDisplacement();
    auto &velocity = model.getVelocity();
    auto &blocked = model.getBlockedDOFs();
    auto &mesh_ref = model.getFEEngine().getMesh();
    const auto &moving_nodes = mesh_ref.getElementGroup("moving-block-back").getNodeGroup().getNodes();

    auto &increment = model.getIncrement();

    //  /$$   /$$                                             /$$       /$$$$$$$  /$$
    // | $$$ | $$                                            | $$      | $$__  $$| $$
    // | $$$$| $$  /$$$$$$   /$$$$$$  /$$$$$$/$$$$   /$$$$$$ | $$      | $$  \ $$| $$$$$$$   /$$$$$$   /$$$$$$$  /$$$$$$
    // | $$ $$ $$ /$$__  $$ /$$__  $$| $$_  $$_  $$ |____  $$| $$      | $$$$$$$/| $$__  $$ |____  $$ /$$_____/ /$$__  $$
    // | $$  $$$$| $$  \ $$| $$  \__/| $$ \ $$ \ $$  /$$$$$$$| $$      | $$____/ | $$  \ $$  /$$$$$$$|  $$$$$$ | $$$$$$$$
    // | $$\  $$$| $$  | $$| $$      | $$ | $$ | $$ /$$__  $$| $$      | $$      | $$  | $$ /$$__  $$ \____  $$| $$_____/
    // | $$ \  $$|  $$$$$$/| $$      | $$ | $$ | $$|  $$$$$$$| $$      | $$      | $$  | $$|  $$$$$$$ /$$$$$$$/|  $$$$$$$
    // |__/  \__/ \______/ |__/      |__/ |__/ |__/ \_______/|__/      |__/      |__/  |__/ \_______/|_______/  \_______/

    for (akantu::Int step_now = 0; step_now < PRESS_STEPS; ++step_now) {
        // model.applyBC(akantu::BC::Dirichlet::IncrementValue(dx_step, akantu::_x), "moving-block-back");
        // model.solveStep(*solver_ntn, "explicit_lumped");

        // increment.set(0.0);
        // blocked.set(false);

        // for (auto n : moving_nodes) {
        //     displacement(n, akantu::_x) += dx_step;
        //     velocity(n, akantu::_x) = dx_step / dt;
        //     increment(n, akantu::_x) = dx_step;
        //     blocked(n, akantu::_x) = true;
        // }

        model.solveStep(*solver_ntn, "explicit_lumped");

        if (step_now % DUMP_INTERVAL == 0) {
            model.dump();
        }

        auto now = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = now - start_time;

        double percent_press = 100.0 * static_cast<double>(step_now) / static_cast<double>(PRESS_STEPS);
        const double denom = (step_now > 0) ? static_cast<double>(step_now) : 1.0;
        const double time_per_iter = elapsed.count() / denom;
        const double remaining = time_per_iter * SHEAR_STEPS - elapsed.count();

        if (prank == 0 && step_now % 100 == 0) {
            std::cout << "[SIM] Step " << std::setw(8) << step_now << "/"
                      << PRESS_STEPS
                      << " | Elapsed: " << std::fixed << std::setprecision(1)
                      << std::setw(8) << elapsed.count() << " s"
                      << " | Progress: " << std::setw(6) << std::fixed
                      << std::setprecision(1) << percent_press << " %"
                      << " | ETA: " << std::fixed << std::setprecision(1)
                      << std::setw(8) << std::max(0.0, remaining) << " s\n";
        }
    }

    if (prank == 0) {
        std::cout << "[SIM] Starting time integration for shear phase for "
                  << (SIMULATION_TIME / ms) << " ms (" << SHEAR_STEPS
                  << " steps)\n";
    }

    //   /$$$$$$  /$$                                           /$$$$$$$  /$$
    //  /$$__  $$| $$                                          | $$__  $$| $$
    // | $$  \__/| $$$$$$$   /$$$$$$   /$$$$$$   /$$$$$$       | $$  \ $$| $$$$$$$   /$$$$$$   /$$$$$$$  /$$$$$$
    // |  $$$$$$ | $$__  $$ /$$__  $$ |____  $$ /$$__  $$      | $$$$$$$/| $$__  $$ |____  $$ /$$_____/ /$$__  $$
    //  \____  $$| $$  \ $$| $$$$$$$$  /$$$$$$$| $$  \__/      | $$____/ | $$  \ $$  /$$$$$$$|  $$$$$$ | $$$$$$$$
    //  /$$  \ $$| $$  | $$| $$_____/ /$$__  $$| $$            | $$      | $$  | $$ /$$__  $$ \____  $$| $$_____/
    // |  $$$$$$/| $$  | $$|  $$$$$$$|  $$$$$$$| $$            | $$      | $$  | $$|  $$$$$$$ /$$$$$$$/|  $$$$$$$
    //  \______/ |__/  |__/ \_______/ \_______/|__/            |__/      |__/  |__/ \_______/|_______/  \_______/

    for (akantu::Int step_now = 0; step_now < SHEAR_STEPS; ++step_now) {
        model.applyBC(akantu::BC::Dirichlet::IncrementValue(dy, akantu::_y), "moving-block-right");

        model.solveStep(*solver_ntn, "explicit_lumped");

        if (step_now % DUMP_INTERVAL == 0) {
            model.dump();
        }

        auto now = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = now - start_time;

        const double denom = (step_now > 0) ? static_cast<double>(step_now) : 1.0;
        const double time_per_iter = elapsed.count() / denom;
        const double remaining = time_per_iter * SHEAR_STEPS - elapsed.count();

        double percent_shear =
            100.0 * static_cast<double>(step_now) / static_cast<double>(SHEAR_STEPS);

        if (prank == 0 && step_now % 100 == 0) {
            std::cout << "[SIM] Step " << std::setw(8) << step_now << "/"
                      << SHEAR_STEPS << " | Elapsed: " << std::fixed
                      << std::setprecision(1) << std::setw(8) << elapsed.count()
                      << " s"
                      << " | Progress: " << std::setw(6) << std::fixed
                      << std::setprecision(1) << percent_shear << " %"
                      << " | ETA: " << std::fixed << std::setprecision(1)
                      << std::setw(8) << std::max(0.0, remaining) << " s\n";
        }
    }

    akantu::finalize();
    return 0;
}