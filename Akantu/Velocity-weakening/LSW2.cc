#include "aka_common.hh"
#include "mesh.hh"
#include "coupler_solid_contact.hh"
#include "solid_mechanics_model.hh"

#include "ntn_contact.hh"
#include "ntn_friction.hh"
#include "ntn_friclaw_linear_slip_weakening.hh"
#include "ntn_fricreg_no_regularisation.hh"

#include "non_linear_solver.hh"

#include <iostream>
#include <iomanip>
#include <chrono>
#include <cmath>


int main(int argc, char *argv[])
{
    constexpr akantu::Int sd = 3;
    constexpr akantu::Real us = 1e-6;
    constexpr akantu::Real ms = 1e-3;
    constexpr akantu::Real time_factor = 0.05;
    constexpr akantu::Real SIMULATION_TIME = 20 * ms;
    constexpr int PMMA_thickness = 50;

    const std::string mesh_file = "../../../Models/" + std::to_string(PMMA_thickness) + "mm-BS-PMMA.msh";
    const std::string mat_file = "../../../Materials/NTN-LSW.dat";

    akantu::initialize(mat_file, argc, argv);

    akantu::Mesh mesh(sd);
    const auto &comm = akantu::Communicator::getStaticCommunicator();
    const akantu::Int prank = comm.whoAmI();

    if (prank == 0)
    {
        mesh.read(mesh_file);
    }

    mesh.distribute();

    akantu::CouplerSolidContact coupler(mesh);

    auto &solid = coupler.getSolidMechanicsModel();
    auto &contact = coupler.getContactMechanicsModel();
    auto &&selector = std::make_shared<akantu::MeshDataMaterialSelector<std::string>>("physical_names", solid);
    solid.setMaterialSelector(selector);

    coupler.initFull(akantu::_analysis_method = akantu::_explicit_lumped_mass);

    akantu::Real dt = solid.getStableTimeStep() * time_factor;
    coupler.setTimeStep(dt);
    if (prank == 0)
    {
        std::cout << "dt = " << dt << " s (" << dt / us << " us)\n";
    }

    auto &&surface_selector = std::make_shared<akantu::PhysicalSurfaceSelector>(mesh);
    contact.getContactDetector().setSurfaceSelector(surface_selector);

    coupler.setTimeStep(dt);
    coupler.setBaseName(std::to_string(PMMA_thickness) + "mm-PMMA-NTN-LSW");
    coupler.addDumpFieldVector("displacement");
    coupler.addDumpFieldVector("velocity");
    coupler.addDumpFieldVector("normals");
    coupler.addDumpFieldVector("tangents");
    coupler.addDumpFieldVector("contact_force");
    coupler.addDumpFieldVector("external_force");
    coupler.addDumpFieldVector("internal_force");
    coupler.addDumpField("areas");
    coupler.addDumpField("stress");
    coupler.addDumpField("blocked_dofs");

    coupler.getVelocity().set(0.);
    coupler.getDisplacement().set(0.);

}