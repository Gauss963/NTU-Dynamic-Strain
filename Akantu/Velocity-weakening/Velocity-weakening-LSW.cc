#include "ntn_friclaw_linear_slip_weakening.hh"
#include "solid_mechanics_model.hh"
#include "ntn_contact.hh"
#include "aka_common.hh"
#include "mesh.hh"

#include <iostream>
#include <iomanip>
#include <chrono>
#include <cmath>
#include <algorithm>

int main(int argc, char *argv[])
{
    constexpr akantu::Int sd = 3;
    constexpr akantu::Real us = 1e-6;
    constexpr akantu::Real ms = 1e-3;

    constexpr int PMMA_thickness = 50;
    constexpr akantu::Real time_factor = 0.5;

    const std::string mesh_file = "../../../Models/" + std::to_string(PMMA_thickness) + "mm-BS-PMMA.msh";
    const std::string mat_file = "../../../Materials/NTN-LSW.dat";

    akantu::initialize(mat_file, argc, argv);

    akantu::Mesh mesh(sd);
    mesh.read(mesh_file);
    mesh.createGroupsFromMeshData<std::string>("physical_names");

    std::cout << "Cells (3D): " << mesh.getNbElement(mesh.getSpatialDimension()) << "\n";
    std::cout << "Faces (2D): " << mesh.getNbElement(mesh.getSpatialDimension() - 1) << "\n";

    akantu::SolidMechanicsModel model(mesh);

    auto bulk_selector = std::make_shared<akantu::MeshDataMaterialSelector<std::string>>("physical_names", model);
    model.setMaterialSelector(bulk_selector);

    model.initFull(akantu::_analysis_method = akantu::_explicit_lumped_mass);
    model.assembleMassLumped();

    akantu::Real dt = model.getStableTimeStep() * time_factor;
    model.setTimeStep(dt);
    std::cout << "dt = " << dt << " s (" << dt / us << " us)\n";

    std::string SYS_PATH = "../../../../../../../../../../";
    std::string DUMP_PATH = "Volumes/Gauss-T7/Dump_5mm/";
    model.setBaseName(std::to_string(PMMA_thickness) + "mm-PMMA-NTN-LSW");

    model.addDumpFieldVector("displacement");
    model.addDumpFieldVector("velocity");
    model.addDumpFieldVector("internal_force");
    model.addDumpField("stress");
    model.addDumpField("grad_u");

    model.getVelocity().set(0.);
    model.getDisplacement().set(0.);

    akantu::NTNContact ntn_contact(model, "contact");
    ntn_contact.addSurfacePair("moving-block-front", "stationary-block-back", akantu::_x);
    ntn_contact.updateInternalData();

    
    
    const akantu::Real normalStress = 8.0; // MPa
    const akantu::Real shearDisp = 3.0;    // mm
    const akantu::Real riseEnd = 0.70;     // Fullfill Deip. at 70%
    akantu::Real pushValue = 0.0;          // mm

    akantu::Vector<akantu::Real, 3> t_front{normalStress, 0.0, 0.0};
    model.applyBC(akantu::BC::Neumann::FromTraction(t_front), "moving-block-back");
    model.applyBC(akantu::BC::Dirichlet::FixedValue(0., akantu::_y), "stationary-block-left");
    model.applyBC(akantu::BC::Dirichlet::FixedValue(0., akantu::_x), "stationary-block-front");


    const akantu::Real SIMULATION_TIME = 0.4 * ms;
    const akantu::Int MAX_STEPS = ceil(SIMULATION_TIME / dt);
    const akantu::Int DUMP_INTERVAL = 3;

    std::cout << "Starting time integration for " << (SIMULATION_TIME / ms)
              << " ms (" << MAX_STEPS << " steps)\n";

    auto start_time = std::chrono::high_resolution_clock::now();

    for (akantu::Int s = 0; s < MAX_STEPS; ++s)
    {
        akantu::Real progress = static_cast<akantu::Real>(s) / MAX_STEPS;
        if (progress < riseEnd)
        {
            pushValue = (progress / riseEnd) * shearDisp;
        }
        else
        {
            pushValue = shearDisp;
        }

        model.applyBC(akantu::BC::Dirichlet::FixedValue(pushValue, akantu::_y), "moving-block-right");

        // --------------------------
        // Contact + friction contribution BEFORE solveStep
        // --------------------------
        // 常見流程：
        // 1) ntn_contact.computeContact();   // 算 gap / normal traction / active set
        // 2) ntn_contact.applyFriction();    // 依 friction law 更新剪力 traction
        // 3) model.solveStep();              // explicit step uses updated internal/external forces
        //
        // 函式名稱依版本不同，保留骨架即可：
        // ntn_contact.computeContact();
        // ntn_contact.computeFriction();  // 或 updateFriction(), applyContactForces() 等

        ntn_contact.updateInternalData();
        ntn_contact.computeContactPressure();
        model.solveStep();

        if (s % DUMP_INTERVAL == 0)
            model.dump();

        auto now = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = now - start_time;
        double time_per_iter = (s > 0) ? (elapsed.count() / s) : 0.0;
        double remaining = (s > 0) ? (time_per_iter * MAX_STEPS - elapsed.count()) : 0.0;

        std::cout << "Step " << std::setw(8) << s << "/" << MAX_STEPS
                  << " | Elapsed: " << std::fixed << std::setprecision(1)
                  << std::setw(8) << elapsed.count() << " s"
                  << " | ETA: " << std::fixed << std::setprecision(1)
                  << std::setw(8) << remaining << " s\r";
        std::cout.flush();
    }

    std::cout << "\n";
    akantu::finalize();
    return 0;
}