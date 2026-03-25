#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <set>
#include <string>
#include <unordered_set>
#include <vector>

#include "aka_common.hh"
#include "mesh.hh"
#include "mesh_partition_mesh_data.hh"
#include "mesh_utils.hh"

#include "solid_mechanics_model.hh"

#include "ntn_base_contact.hh"
#include "ntn_contact_solvercallback.hh"
#include "ntn_initiation_function.hh"

#include "ntn_friclaw_linear_slip_weakening.hh"
#include "ntn_fricreg_no_regularisation.hh"
#include "ntn_friction.hh"

int main(int argc, char *argv[]) {
    constexpr akantu::Int DIMENSION = 2;
    // constexpr akantu::Real us = 1e-6;
    constexpr akantu::Real ms = 1e-3;

    constexpr akantu::Real TIME_FACTOR = 0.02;
    constexpr akantu::Real SIMULATION_TIME = 20 * ms;
    constexpr int PMMA_thickness = 50;
    constexpr akantu::Int normal_dir = 0;
    const std::string slave_surface = "stationary-block-back";
    const std::string master_surface = "moving-block-front";

    const std::string mat_file = "../../../Materials/NTN-LSW.dat";
    std::string mesh_file;
    if (DIMENSION == 2) {
        mesh_file = "../../../Models/2D-PMMA.msh";
    } else {
        mesh_file = "../../../Models/" + std::to_string(PMMA_thickness) + "mm-BS-PMMA.msh";
    }

    akantu::initialize(mat_file, argc, argv);

    akantu::Mesh mesh(DIMENSION);
    const auto &comm = akantu::Communicator::getStaticCommunicator();
    const akantu::Int prank = comm.whoAmI();

#include "NTN_mpi.hh"

    // if (prank == 0) {
    //     mesh.read(mesh_file);
    // }
    // mesh.distribute();


    akantu::SolidMechanicsModel model(mesh);

    auto solver_ntn = std::make_unique<akantu::NTNContactSolverCallback>(model, slave_surface, master_surface, normal_dir, TIME_FACTOR);
    auto contact = solver_ntn->getContact();
    contact->initParallel();

    auto friction = solver_ntn->getFriction();
    auto &msh = model.getFEEngine().getMesh();

    const akantu::Real mu_s = friction->get("mu_s");
    const akantu::Real mu_k = friction->get("mu_k");
    const akantu::Real d_c = friction->get("d_c");
    std::cout << "[SIM] mu_s = " << mu_s << " | mu_k = " << mu_k << " | d_c = " << d_c << "mm\n";

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

    akantu::Real dt = model.getTimeStep();

    model.setBaseName(std::to_string(PMMA_thickness) + "mm-PMMA-NTN-LSW");

    if (DIMENSION == 3) {
        model.addDumpFieldVector("strain");
    } else {
        model.addDumpField("strain");
    }
    model.addDumpField("stress");
    model.addDumpField("mass");
    model.addDumpFieldVector("displacement");
    model.addDumpFieldVector("internal_force");
    model.addDumpFieldVector("external_force");

    model.getVelocity().set(0.0);
    model.getDisplacement().set(0.0);

    // constexpr akantu::Real shearDisp = 3.0;
    constexpr akantu::Real riseEnd = 0.30;
    constexpr akantu::Int TOTAL_FRAMES = 2400;

    const akantu::Int SHEAR_STEPS = std::ceil(SIMULATION_TIME / dt);
    const akantu::Int PRESS_STEPS = SHEAR_STEPS;
    const akantu::Int TAU_K_START_STEP = PRESS_STEPS * 0.99;

    const akantu::Int DUMP_INTERVAL = std::max<akantu::Int>(1, SHEAR_STEPS / TOTAL_FRAMES);

    model.applyBC(akantu::BC::Dirichlet::FixedValue(0., akantu::_x), "stationary-block-front");
    model.applyBC(akantu::BC::Dirichlet::FixedValue(0., akantu::_y), "stationary-block-left");

    model.dump();

    constexpr akantu::Real NORMAL_STRESS = 16.0;
    akantu::Vector<akantu::Real, DIMENSION> normal_traction;
    akantu::Vector<akantu::Real, DIMENSION> shear_traction_init;
    normal_traction.set(0.0);
    shear_traction_init.set(0.0);

    const akantu::Real tau_k = (500.0 / 145.0) * mu_k * NORMAL_STRESS;
    const akantu::Real tau_s = (500.0 / 145.0) * mu_s * NORMAL_STRESS;

    const akantu::Int rise_steps = std::max<akantu::Int>(1, std::ceil(riseEnd * SHEAR_STEPS));
    const akantu::Real dtau = (tau_s - tau_k) / static_cast<akantu::Real>(rise_steps);

    normal_traction(0) = NORMAL_STRESS;

    model.applyBC(akantu::BC::Neumann::FromTraction(normal_traction), "moving-block-back");

    auto start_time = std::chrono::high_resolution_clock::now();
    if (prank == 0) {
        std::cout << "[SIM] Starting time integration for compression phase for "
                  << PRESS_STEPS << " steps\n";
    }

    //  /$$   /$$                                             /$$       /$$$$$$$  /$$
    // | $$$ | $$                                            | $$      | $$__  $$| $$
    // | $$$$| $$  /$$$$$$   /$$$$$$  /$$$$$$/$$$$   /$$$$$$ | $$      | $$  \ $$| $$$$$$$   /$$$$$$   /$$$$$$$  /$$$$$$
    // | $$ $$ $$ /$$__  $$ /$$__  $$| $$_  $$_  $$ |____  $$| $$      | $$$$$$$/| $$__  $$ |____  $$ /$$_____/ /$$__  $$
    // | $$  $$$$| $$  \ $$| $$  \__/| $$ \ $$ \ $$  /$$$$$$$| $$      | $$____/ | $$  \ $$  /$$$$$$$|  $$$$$$ | $$$$$$$$
    // | $$\  $$$| $$  | $$| $$      | $$ | $$ | $$ /$$__  $$| $$      | $$      | $$  | $$ /$$__  $$ \____  $$| $$_____/
    // | $$ \  $$|  $$$$$$/| $$      | $$ | $$ | $$|  $$$$$$$| $$      | $$      | $$  | $$|  $$$$$$$ /$$$$$$$/|  $$$$$$$
    // |__/  \__/ \______/ |__/      |__/ |__/ |__/ \_______/|__/      |__/      |__/  |__/ \_______/|_______/  \_______/

    for (akantu::Int step_now = 0; step_now < PRESS_STEPS; ++step_now) {

        if (step_now == TAU_K_START_STEP) {
            akantu::Vector<akantu::Real, DIMENSION> shear_traction_tau_k;
            shear_traction_tau_k.set(0.0);
            shear_traction_tau_k(1) = tau_k;
            model.applyBC(akantu::BC::Neumann::FromTraction(shear_traction_tau_k), "moving-block-right");
        }

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
        // model.applyBC(akantu::BC::Dirichlet::IncrementValue(dy, akantu::_y), "moving-block-right");
        if (step_now < rise_steps) {
            akantu::Vector<akantu::Real, DIMENSION> shear_traction_increment;
            shear_traction_increment.set(0.0);
            shear_traction_increment(1) = dtau;
            model.applyBC(akantu::BC::Neumann::FromTraction(shear_traction_increment), "moving-block-right");
        }

        model.solveStep(*solver_ntn, "explicit_lumped");

        if (step_now % DUMP_INTERVAL == 0) {
            model.dump();
        }

        auto now = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = now - start_time;

        const double denom = (step_now > 0) ? static_cast<double>(step_now) : 1.0;
        const double time_per_iter = elapsed.count() / denom;
        const double remaining = time_per_iter * SHEAR_STEPS - elapsed.count();

        double percent_shear = 100.0 * static_cast<double>(step_now) / static_cast<double>(SHEAR_STEPS);

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