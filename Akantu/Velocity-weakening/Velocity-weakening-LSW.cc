#include "aka_common.hh"
#include "mesh.hh"
#include "solid_mechanics_model.hh"

#include "ntn_contact.hh"
#include "ntn_friction.hh"
#include "ntn_friclaw_linear_slip_weakening.hh"
#include "ntn_fricreg_no_regularisation.hh"

#include <iostream>
#include <iomanip>
#include <chrono>
#include <cmath>

int main(int argc, char *argv[])
{
    constexpr akantu::Int sd = 3;
    constexpr akantu::Real us = 1e-6;
    constexpr akantu::Real ms = 1e-3;

    constexpr int PMMA_thickness = 50;
    constexpr akantu::Real time_factor = 0.5;

    const std::string mesh_file =
        "../../../Models/" + std::to_string(PMMA_thickness) + "mm-BS-PMMA.msh";
    const std::string mat_file = "../../../Materials/NTN-LSW.dat";

    akantu::initialize(mat_file, argc, argv);

    akantu::Mesh mesh(sd);
    mesh.read(mesh_file);
    mesh.createGroupsFromMeshData<std::string>("physical_names");

    std::cout << "Cells (3D): " << mesh.getNbElement(mesh.getSpatialDimension()) << "\n";
    std::cout << "Faces (2D): " << mesh.getNbElement(mesh.getSpatialDimension() - 1) << "\n";

    akantu::SolidMechanicsModel model(mesh);

    auto bulk_selector =
        std::make_shared<akantu::MeshDataMaterialSelector<std::string>>("physical_names", model);
    model.setMaterialSelector(bulk_selector);

    model.initFull(akantu::_analysis_method = akantu::_explicit_lumped_mass);
    model.assembleMassLumped();

    akantu::Real dt = model.getStableTimeStep() * time_factor;
    model.setTimeStep(dt);
    std::cout << "dt = " << dt << " s (" << dt / us << " us)\n";

    model.setBaseName(std::to_string(PMMA_thickness) + "mm-PMMA-NTN-LSW");

    model.addDumpFieldVector("displacement");
    model.addDumpFieldVector("velocity");
    model.addDumpFieldVector("internal_force");
    model.addDumpFieldVector("external_force");
    model.addDumpField("stress");
    model.addDumpField("grad_u");

    model.getVelocity().set(0.);
    model.getDisplacement().set(0.);

    akantu::NTNContact ntn_contact(model, "contact");
    ntn_contact.addSurfacePair("moving-block-front", "stationary-block-back", akantu::_x);
    ntn_contact.initParallel();
    ntn_contact.updateInternalData();

    using Reg = akantu::NTNFricRegNoRegularisation;
    using Fric = akantu::NTNFriction<akantu::NTNFricLawLinearSlipWeakening, Reg>;
    Fric ntn_friction(ntn_contact, "friction");

    const akantu::Real normalStress = 8.0; // MPa
    const akantu::Real shearDisp = 3.0;    // mm
    const akantu::Real riseEnd = 0.70;
    akantu::Real pushValue = 0.0;

    akantu::Vector<akantu::Real, 3> t_normal{normalStress, 0.0, 0.0};
    model.applyBC(akantu::BC::Neumann::FromTraction(t_normal), "moving-block-back");
    model.applyBC(akantu::BC::Dirichlet::FixedValue(0., akantu::_x), "stationary-block-front");
    model.applyBC(akantu::BC::Dirichlet::FixedValue(0., akantu::_y), "stationary-block-left");

    // ---------------------------------------------------------------------------
    // Time control
    // ---------------------------------------------------------------------------
    const akantu::Real SIMULATION_TIME = 0.4 * ms;
    const akantu::Int MAX_STEPS = static_cast<akantu::Int>(std::ceil(SIMULATION_TIME / dt));
    const akantu::Int DUMP_INTERVAL = 3;

    std::cout << "Starting time integration for " << (SIMULATION_TIME / ms)
              << " ms (" << MAX_STEPS << " steps)\n";

    auto start_time = std::chrono::high_resolution_clock::now();

    for (akantu::Int s = 0; s < MAX_STEPS; ++s)
    {
        const akantu::Real progress = static_cast<akantu::Real>(s) / MAX_STEPS;
        pushValue = (progress < riseEnd) ? (progress / riseEnd) * shearDisp : shearDisp;

        model.applyBC(akantu::BC::Dirichlet::FixedValue(pushValue, akantu::_y), "moving-block-right");

        // Contact + friction assembly BEFORE solveStep
        ntn_contact.updateInternalData();
        ntn_contact.computeContactPressure();
        ntn_contact.assembleGlobalContactPressure();

        ntn_friction.assembleGlobalFrictionTraction();

        model.solveStep();

        if (s % DUMP_INTERVAL == 0)
            model.dump();

        auto now = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = now - start_time;

        const double denom = (s > 0) ? static_cast<double>(s) : 1.0;
        const double time_per_iter = elapsed.count() / denom;
        const double remaining = time_per_iter * MAX_STEPS - elapsed.count();

        std::cout << "Step " << std::setw(8) << s << "/" << MAX_STEPS
                  << " | Elapsed: " << std::fixed << std::setprecision(1)
                  << std::setw(8) << elapsed.count() << " s"
                  << " | ETA: " << std::fixed << std::setprecision(1)
                  << std::setw(8) << std::max(0.0, remaining) << " s\r";
        std::cout.flush();
    }

    std::cout << "\n";
    akantu::finalize();
    return 0;
}