#include <chrono>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <ostream>
#include <sstream>
#include <string>

#include "aka_common.hh"
#include "dumpable_iohelper.hh"
#include "dumper_text.hh"
#include "dumper_variable.hh"
#include "mesh_utils.hh"
#include "ntn_base_contact.hh"
#include "ntn_contact.hh"
#include "ntn_contact_solvercallback.hh"
#include "ntn_friclaw_linear_slip_weakening.hh"
#include "ntn_fricreg_no_regularisation.hh"
#include "ntn_friction.hh"
#include "ntn_initiation_function.hh"
#include "solid_mechanics_model.hh"

int main(int argc, char *argv[]) {

    constexpr akantu::Int sd = 3;
    constexpr akantu::Real us = 1e-6;
    constexpr akantu::Real ms = 1e-3;
    constexpr akantu::Real TIME_FACTOR = 0.05;
    constexpr akantu::Real SIMULATION_TIME = 2 * ms;
    // constexpr int PMMA_thickness = 50;
    constexpr int PMMA_thickness = 2;

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

    return 0;
}