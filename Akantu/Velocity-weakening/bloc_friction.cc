#include "coupler_solid_contact.hh"
#include "non_linear_solver.hh"
#include <fstream>
#include <iomanip>
#include <vector>

using namespace akantu;

int main(int argc, char *argv[])
{
    initialize("material.dat", argc, argv);

    Mesh mesh(2);
    mesh.read("bloc.msh");

    // Create the model
    CouplerSolidContact coupler(mesh);

    // Initialize each model
    auto &solid = coupler.getSolidMechanicsModel();
    auto &contact = coupler.getContactMechanicsModel();
    auto &&selector = std::make_shared<MeshDataMaterialSelector<std::string>>("physical_names", solid);
    solid.setMaterialSelector(selector);

    // Initialize the coupler
    coupler.initFull(_analysis_method = _explicit_lumped_mass);
    Real time_step = solid.getStableTimeStep() * 0.1;
    coupler.setTimeStep(time_step);
    std::cout << "Time step: " << time_step << std::endl;

    // Setup the contact
    auto &&surface_selector = std::make_shared<PhysicalSurfaceSelector>(mesh);
    contact.getContactDetector().setSurfaceSelector(surface_selector);

    // Configuration of the dumper
    coupler.setBaseName("bloc_friction");
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

    // Add the boundary conditions
    solid.applyBC(BC::Dirichlet::FixedValue(0.0, _x), "XFixed");
    solid.applyBC(BC::Dirichlet::FixedValue(0.0, _y), "YFixed");
    solid.applyBC(BC::Dirichlet::FixedValue(0.0, _x), "loading");
    solid.applyBC(BC::Dirichlet::FixedValue(0.0, _y), "loading");

    // Register velocity and gaps for future damping
    auto &velocity = solid.getVelocity();
    auto &gaps = contact.getGaps();

    // Dump the initial state
    coupler.dump();

    // First loop : compression of the bloc against the wall
    for (int s = 0; s < 10000; s++)
    {
        solid.applyBC(BC::Dirichlet::IncrementValue(-1.0 / 10000, _y), "loading");

        coupler.solveStep();

        // damping velocities only along the contacting zone
        for (auto &&tuple : zip(gaps, make_view(velocity, 2)))
        {
            auto &gap = std::get<0>(tuple);
            auto &vel = std::get<1>(tuple);
            if (gap > 0)
            {
                vel *= 0.99;
            }
        }

        if (s % 100 == 0)
        {
            std::cout << "Step " << s << "\t\r" << std::flush;
            coupler.dump();
        }
    }

    std::cout << "Compression done !" << std::endl;

    // Second loop : sliding of the bloc against the wall
    for (int s = 0; s < 10000; s++)
    {
        solid.applyBC(BC::Dirichlet::IncrementValue(2.0 / 10000, _x), "loading");

        coupler.solveStep();

        // damping velocities only along the contacting zone
        for (auto &&tuple : zip(gaps, make_view(velocity, 2)))
        {
            auto &gap = std::get<0>(tuple);
            auto &vel = std::get<1>(tuple);
            if (gap > 0)
            {
                vel *= 0.99;
            }
        }

        if (s % 100 == 0)
        {
            std::cout << "Step " << s << "\t\r" << std::flush;
            coupler.dump();
        }
    }

    std::cout << "Sliding done !" << std::endl;

    // Third loop : stabilization of the system (no external force)
    for (int s = 0; s < 10000; s++)
    {
        coupler.solveStep();

        // damping velocities only along the contacting zone
        for (auto &&tuple : zip(gaps, make_view(velocity, 2)))
        {
            auto &gap = std::get<0>(tuple);
            auto &vel = std::get<1>(tuple);
            if (gap > 0)
            {
                vel *= 0.99;
            }
        }

        if (s % 100 == 0)
        {
            std::cout << "Step " << s << "\t\r" << std::flush;
            coupler.dump();
        }
    }

    return 0;
}