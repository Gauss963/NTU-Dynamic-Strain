#include "solid_mechanics_model.hh"
#include "mesh.hh"
#include "aka_common.hh"

#include <omp.h>
#include <iostream>
#include <chrono>
#include <cmath>
#include <algorithm>
#include <vector>
#include <limits>

namespace hertz
{
    // (McLaskey 2009) Hertzian STF：f(t) = fmax * sin(pi t / tc)^(3/2), 0<t<tc
    inline akantu::Real delta(akantu::Real E, akantu::Real nu)
    {
        return (1.0 - nu * nu) / (M_PI * E);
    }

    struct Params
    {
        akantu::Real rho;     // kg/m^3  (block / target)
        akantu::Real R;       // m       (sphere radius)
        akantu::Real v;       // m/s     (impact velocity)
        akantu::Real E1, nu1; // sphere (steel)
        akantu::Real E2, nu2; // block  (PMMA)
    };

    struct Precomp
    {
        akantu::Real tc;   // s
        akantu::Real fmax; // N
    };

    inline Precomp precompute(const Params &p)
    {
        const akantu::Real del1 = delta(p.E1, p.nu1);
        const akantu::Real del2 = delta(p.E2, p.nu2);
        const akantu::Real A = (4.0 * p.rho * M_PI * (del1 + del2) / 3.0);
        const akantu::Real tc = 4.53 * std::pow(A, 0.4) * p.R * std::pow(p.v, -0.2);
        const akantu::Real fmax = 1.917 * std::pow(p.rho, 0.6) * std::pow(del1 + del2, -0.4) * (p.R * p.R) * std::pow(p.v, 1.2);
        return {tc, fmax};
    }

    inline akantu::Real stf(akantu::Real t, const Precomp &c)
    {
        if (t <= 0.0 || t >= c.tc)
            return 0.0;
        akantu::Real s = std::sin(M_PI * t / c.tc);
        return c.fmax * std::pow(s, 1.5);
    }
} // namespace hertz

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

    akantu::Real dt_stable = model.getStableTimeStep();
    akantu::Real dt = std::min(dt_stable * 0.5, 1e-7);
    model.setTimeStep(dt);
    std::cout << "[✓] Time step set: " << dt << " s (stable*0.5=" << dt_stable * 0.5 << ")\n";

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

    // === Apply BC: fix Z direction of 4 corner bottom nodes ===
    std::vector<akantu::Int> bottom_corner_nodes;
    const akantu::Real eps_z = 1e-3;     // mm
    const akantu::Real eps_corner = 1.0; // mm

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

    auto &blocked_dofs = model.getBlockedDOFs();
    blocked_dofs.set(false);
    for (const auto &node : bottom_corner_nodes)
        blocked_dofs(node, 2) = true; // z-dir
    std::cout << "[✓] Fixed Z on " << bottom_corner_nodes.size() << " bottom-corner nodes\n";

    // === Find top-center node (z ≈ 20 mm; 若找不到嚴格 z=20，退而求其次以最高 z 為準) ===
    akantu::Real x_center = 150.;
    akantu::Real y_center = 150.;
    akantu::Real z_top_nominal = 20.; // mm

    akantu::Int node_closest = -1;
    akantu::Real min_dist = std::numeric_limits<akantu::Real>::max();

    for (akantu::Int node = 0; node < mesh.getNbNodes(); ++node)
    {
        akantu::Real x = nodes(node, 0);
        akantu::Real y = nodes(node, 1);
        akantu::Real z = nodes(node, 2);
        if (std::abs(z - z_top_nominal) > 1e-3)
            continue;
        akantu::Real dx = x - x_center, dy = y - y_center;
        akantu::Real d2 = dx * dx + dy * dy;
        if (d2 < min_dist)
        {
            min_dist = d2;
            node_closest = node;
        }
    }
    if (node_closest < 0)
    {
        akantu::Real zmax = -1e9;
        for (akantu::Int node = 0; node < mesh.getNbNodes(); ++node)
        {
            zmax = std::max(zmax, nodes(node, 2));
        }
        min_dist = std::numeric_limits<akantu::Real>::max();
        for (akantu::Int node = 0; node < mesh.getNbNodes(); ++node)
        {
            if (std::abs(nodes(node, 2) - zmax) > 1e-6)
                continue;
            akantu::Real dx = nodes(node, 0) - x_center;
            akantu::Real dy = nodes(node, 1) - y_center;
            akantu::Real d2 = dx * dx + dy * dy;
            if (d2 < min_dist)
            {
                min_dist = d2;
                node_closest = node;
            }
        }
    }
    std::cout << "[✓] Central top node ID: " << node_closest << "\n";

    // === Hertzian STF 參數（請依你的材料檔更新 PMMA 參數） ===
    // 目標材料（block, PMMA）
    const akantu::Real rho_block = 1.19e-9 * 1e12; // = 1190.0 kg/m^3
    const akantu::Real E_block = 3.2e3 * 1e6; // = 3.2e9 Pa
    const akantu::Real nu_block = 0.35; // 無單位

    // 球（鋼）
    const akantu::Real E_ball = 208.197e9; // Pa
    const akantu::Real nu_ball = 0.286;

    // 球半徑與落高 → 撞擊速度
    const akantu::Real R_sphere = 1e-3; // m
    const akantu::Real h_drop = 0.300;  // m
    const akantu::Real g = 9.80665;     // m/s^2
    const akantu::Real v_impact = std::sqrt(2.0 * g * h_drop);

    hertz::Params hp{rho_block, R_sphere, v_impact,
                     E_ball, nu_ball, E_block, nu_block};
    auto hc = hertz::precompute(hp);

    // 可選：如需以 M0 一致性作幅度調整，可設 scale_factor（預設 1，不額外縮放）
    // akantu::Real M0_scaling = 1.748 * hc.fmax * hc.tc / M_PI; // 依你 Python 的計算
    const akantu::Real scale_factor = 1.0;

    std::cout << "[i] Impact v = " << v_impact << " m/s\n";
    std::cout << "[i] Hertz tc = " << hc.tc * 1e6 << " μs, fmax = " << hc.fmax << " N\n";

    // === 總模擬時間：200 μs ===
    const akantu::Real SIM_TIME = 200e-6; // 200 μs
    const akantu::Int max_steps = static_cast<akantu::Int>(std::ceil(SIM_TIME / dt));
    std::cout << "[✓] Starting time loop: " << max_steps << " steps (~" << SIM_TIME * 1e6 << " μs)\n";

    auto t0 = std::chrono::high_resolution_clock::now();

    for (akantu::Int step = 0; step < max_steps; ++step)
    {
        akantu::Real t = step * dt;

        f_ext.set(0.);

        akantu::Real fz = -scale_factor * hertz::stf(t, hc); // 作用在 -Z（向下）
        if (node_closest >= 0)
            f_ext(node_closest, 2) = fz;

        model.solveStep();

        if (step % 1 == 0)
            model.dump();

        auto now = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = now - t0;
        double eta = elapsed.count() / (step + 1) * (max_steps - step - 1);
        std::cout << "Step " << step << "/" << max_steps
                    << " | t = " << t * 1e6 << " μs"
                    << " | dt = " << dt * 1e6 << " μs"
                    << " | fz = " << fz << " N"
                    << " | ETA: " << eta << " s\n";
    }

    akantu::finalize();
    return 0;
}