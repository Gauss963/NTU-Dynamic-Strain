#include "solid_mechanics_model.hh"
#include "mesh.hh"
#include "aka_common.hh"

#include <omp.h>
#include <iostream>
#include <chrono>
#include <cmath>
#include <algorithm>

int main(int argc, char *argv[])
{
    omp_set_num_threads(12);
    constexpr akantu::Int dim = 3;

    const std::string mesh_file = "../../../Models/Ball-Drop.msh";
    const std::string mat_file = "../../../Materials/material-Ball-Drop.dat";

    akantu::initialize(mat_file, argc, argv);
    std::cout << "[✓] Akantu initialized\n";

    akantu::Mesh mesh(dim);
    mesh.read(mesh_file);
    std::cout << "[✓] Mesh loaded\n";

    const auto &nodes = mesh.getNodes();

    akantu::SolidMechanicsModel model(mesh);
    model.initFull(akantu::_analysis_method = akantu::_explicit_lumped_mass);

    akantu::Real dt = model.getStableTimeStep() * 0.5;
    model.setTimeStep(dt);
    std::cout << "[✓] Time step set: " << dt << " s\n";

    model.assembleMassLumped();
    model.addDumpFieldVector("displacement");
    model.addDumpFieldVector("velocity");
    model.addDumpFieldVector("external_force");
    model.addDumpField("stress");

    auto &vel = model.getVelocity();
    auto &disp = model.getDisplacement();
    auto &f_ext = model.getExternalForce();

    vel.set(0.);
    disp.set(0.);
    f_ext.set(0.);

    // === Apply BC: fix Y direction of 4 corner bottom nodes ===
    std::vector<akantu::Int> bottom_corner_nodes;
    const akantu::Real eps_z = 1e-3;
    const akantu::Real eps_corner = 1.0;

    for (akantu::Int node = 0; node < mesh.getNbNodes(); ++node)
    {
        akantu::Real x = nodes(node, 0);
        akantu::Real y = nodes(node, 1);
        akantu::Real z = nodes(node, 2);

        if (std::abs(z) > eps_z)
            continue;

        bool is_corner =
            (std::abs(x - 0.) < eps_corner && std::abs(y - 0.) < eps_corner) ||
            (std::abs(x - 300.) < eps_corner && std::abs(y - 0.) < eps_corner) ||
            (std::abs(x - 0.) < eps_corner && std::abs(y - 300.) < eps_corner) ||
            (std::abs(x - 300.) < eps_corner && std::abs(y - 300.) < eps_corner);

        if (is_corner)
            bottom_corner_nodes.push_back(node);
    }

    // === 手動設定底部四角節點的 Y 方向固定 ===
    auto &blocked_dofs = model.getBlockedDOFs();
    blocked_dofs.set(false); // 先全部解除

    for (const auto &node : bottom_corner_nodes)
    {
        blocked_dofs(node, 1) = true; // 1 = y-direction
    }

    std::cout << "[✓] Fixed Y on " << bottom_corner_nodes.size() << " corner nodes (manually via getBlockedDOFs)\n";

    std::cout << "[✓] Fixed Y on " << bottom_corner_nodes.size() << " corner nodes\n";

    // === Find top-center node
    akantu::Real x_center = 150.;
    akantu::Real y_center = 150.;
    akantu::Real z_top = 20.; // mm

    akantu::Int node_closest = -1;
    akantu::Real min_dist = 1e10;

    for (akantu::Int node = 0; node < mesh.getNbNodes(); ++node)
    {
        akantu::Real x = nodes(node, 0);
        akantu::Real y = nodes(node, 1);
        akantu::Real z = nodes(node, 2);

        if (std::abs(z - z_top) > 1e-3)
            continue;

        akantu::Real dx = x - x_center;
        akantu::Real dy = y - y_center;
        akantu::Real dist2 = dx * dx + dy * dy;

        if (dist2 < min_dist)
        {
            min_dist = dist2;
            node_closest = node;
        }
    }

    std::cout << "[✓] Central top node ID: " << node_closest << "\n";

    // === Apply force pulse
    const akantu::Real SIM_TIME = 0.01; // sec
    const akantu::Real t_total = 0.002; // sec
    const akantu::Real F_peak = 1e6;    // N
    const akantu::Int max_steps = static_cast<akantu::Int>(SIM_TIME / dt);

    std::cout << "[✓] Starting time loop: " << max_steps << " steps\n";
    auto start_time = std::chrono::high_resolution_clock::now();

    for (akantu::Int step = 0; step < max_steps; ++step)
    {
        akantu::Real t = step * dt;

        f_ext.set(0.);
        if (t < t_total)
        {
            akantu::Real fz = -F_peak * std::sin(3.14159 * t / t_total);
            f_ext(node_closest, 2) = fz; // z-direction
        }

        model.solveStep();

        if (step % 50 == 0)
            model.dump();

        auto now = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = now - start_time;
        double eta = elapsed.count() / (step + 1) * (max_steps - step - 1);
        std::cout << "Step " << step << "/" << max_steps
                    << " | Time: " << t * 1e3 << " ms"
                    << " | ETA: " << eta << " s\n";
    }

    akantu::finalize();
    return 0;
}