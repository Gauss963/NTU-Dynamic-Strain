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

    akantu::SolidMechanicsModel model(mesh);

    auto bulk_selector = std::make_shared<akantu::MeshDataMaterialSelector<std::string>>("physical_names", model);
    model.setMaterialSelector(bulk_selector);
    model.initFull(akantu::_analysis_method = akantu::_explicit_lumped_mass);
    model.assembleMassLumped();

    akantu::Real dt = model.getStableTimeStep() * time_factor;
    if (prank == 0)
    {
        std::cout << "dt = " << dt << " s (" << dt / us << " us)\n";
    }
    model.setTimeStep(dt);
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
    // ntn_contact.addSurfacePair("stationary-block-back", "moving-block-front", akantu::_x);
    ntn_contact.addSurfacePair("moving-block-front", "stationary-block-back", akantu::_x);


    

    int local_pairs = static_cast<int>(ntn_contact.getNbContactNodes());
    if (prank == 0)
    {
        std::cout << "[CHK] local_pairs(rank0)=" << local_pairs  << "\n";
    }




    ntn_contact.initParallel();
    ntn_contact.updateInternalData();

    using Regular = akantu::NTNFricRegNoRegularisation;
    using Friction = akantu::NTNFriction<akantu::NTNFricLawLinearSlipWeakening, Regular>;
    Friction ntn_friction(ntn_contact, "friction");

    const akantu::Real normalStress = 8.0; // MPa
    const akantu::Real shearDisp = 8.0;    // mm
    const akantu::Real riseEnd = 0.30;
    // akantu::Real pushValue = 0.0;

    akantu::Vector<akantu::Real, 3> t_normal{normalStress, 0.0, 0.0};
    model.applyBC(akantu::BC::Neumann::FromTraction(t_normal), "moving-block-back");
    model.applyBC(akantu::BC::Dirichlet::FixedValue(0., akantu::_x), "stationary-block-front");
    model.applyBC(akantu::BC::Dirichlet::FixedValue(0., akantu::_y), "stationary-block-left");

    const akantu::Int TOTAL_FRAMES = 2400;
    const akantu::Int MAX_STEPS = static_cast<akantu::Int>(std::ceil(SIMULATION_TIME / dt));
    const akantu::Int DUMP_INTERVAL = MAX_STEPS / TOTAL_FRAMES;

    const akantu::Int rise_steps = static_cast<akantu::Int>(std::ceil(riseEnd * MAX_STEPS));
    const akantu::Real dy = shearDisp / rise_steps;

    if (prank == 0)
    {
        std::cout << "Starting time integration for " << (SIMULATION_TIME / ms)
                  << " ms (" << MAX_STEPS << " steps)\n";
    }

    auto start_time = std::chrono::high_resolution_clock::now();

    for (akantu::Int s = 0; s < MAX_STEPS; ++s)
    {
        // const akantu::Real progress = static_cast<akantu::Real>(s) / MAX_STEPS;
        // pushValue = (progress < riseEnd) ? (progress / riseEnd) * shearDisp : shearDisp;

        // model.applyBC(akantu::BC::Dirichlet::FixedValue(pushValue, akantu::_y), "moving-block-right");
        model.applyBC(akantu::BC::Dirichlet::IncrementValue(dy, akantu::_y), "moving-block-right");

        ntn_contact.updateInternalData();
        ntn_contact.computeContactPressure();
        ntn_contact.assembleGlobalContactPressure();

        ntn_friction.assembleGlobalFrictionTraction();

        //
        auto &f_ext = model.getExternalForce();
        auto &f_c = ntn_contact.getGlobalContactPressure();
        const auto dim = model.getSpatialDimension();
        for (akantu::Idx n = 0; n < mesh.getNbNodes(); ++n)
        {
            for (akantu::Int d = 0; d < dim; ++d)
            {
                f_ext(n, d) += static_cast<akantu::Real>(f_c(n, d));
            }
        }
        model.solveStep();
        for (akantu::Idx n = 0; n < mesh.getNbNodes(); ++n)
        {
            for (akantu::Int d = 0; d < dim; ++d)
            {
                f_ext(n, d) -= static_cast<akantu::Real>(f_c(n, d));
            }
        }
        //

        // model.solveStep();

        if (s % DUMP_INTERVAL == 0)
        {
            model.dump();
        }
        auto now = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = now - start_time;

        const double denom = (s > 0) ? static_cast<double>(s) : 1.0;
        const double time_per_iter = elapsed.count() / denom;
        const double remaining = time_per_iter * MAX_STEPS - elapsed.count();
        if (prank == 0 && s % 200 == 0)
        {
            std::cout << "Step " << std::setw(8) << s << "/" << MAX_STEPS
                      << " | Elapsed: " << std::fixed << std::setprecision(1)
                      << std::setw(8) << elapsed.count() << " s"
                      << " | ETA: " << std::fixed << std::setprecision(1)
                      << std::setw(8) << std::max(0.0, remaining) << " s\n";
        }

        if (prank == 0 && s % 200 == 0)
        {
            auto &normals_arr = ntn_contact.getNormals().getArray(); // Array<Real>
            akantu::Real nx_avg = 0.0;

            for (akantu::Idx i = 0; i < normals_arr.size(); ++i)
            {
                nx_avg += normals_arr(i, 0);
            }
            nx_avg /= std::max<akantu::Idx>(1, normals_arr.size());
            akantu::Real nx_min = 1e100;
            akantu::Real nx_max = -1e100;

            for (akantu::Idx i = 0; i < normals_arr.size(); ++i)
            {
                const auto nx = normals_arr(i, 0);
                nx_min = std::min(nx_min, nx);
                nx_max = std::max(nx_max, nx);
            }
            std::cout << "[DBG] nx[min,max]=[" << nx_min << "," << nx_max << "]\n";
        }
    }

    if (prank == 0)
    {
        std::cout << "\n";
    }
    akantu::finalize();
    return 0;
}