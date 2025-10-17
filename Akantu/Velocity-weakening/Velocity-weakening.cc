#include "solid_mechanics_model.hh"
#include "contact_mechanics_model.hh"
#include "coupler_solid_contact.hh"
#include "mesh.hh"
#include "aka_common.hh"
#include <omp.h>

#include <iostream>
#include <algorithm>
#include <cmath>

using namespace akantu;

int main(int argc, char *argv[])
{
    // omp_set_num_threads(8);
    // std::cout << "Using " << omp_get_max_threads() << " OpenMP threads.\n";

    constexpr Int sd = 3;

    const std::string MESHFILE = "../../../Models/50mm-PMMA.msh";
    const std::string MATERIALFILE = "../../../Materials/material.dat";

    initialize(MATERIALFILE, argc, argv);
    Mesh mesh(sd);
    mesh.read(MESHFILE);

    CouplerSolidContact coupler(mesh);
    auto &solid = coupler.getSolidMechanicsModel();
    auto &contact = coupler.getContactMechanicsModel();

    // 用 Gmsh 的 physical surfaces 當接觸面（需 mesh 內有 friction_master / friction_slave）
    auto surf_sel = std::make_shared<PhysicalSurfaceSelector>(mesh);
    contact.getContactDetector().setSurfaceSelector(surf_sel);

    // 顯式初始化
    coupler.initFull(_analysis_method = _explicit_lumped_mass);

    // 時間步
    Real dt = solid.getStableTimeStep() * 0.5;
    coupler.setTimeStep(dt);
    std::cout << "dt = " << dt << std::endl;

    // 場變數
    Array<Real> &vel = solid.getVelocity();            // [nb_nodes x 3]
    Array<Real> &force = solid.getExternalForce();     // [nb_nodes x 3]
    Array<Real> &disp = solid.getDisplacement();       // [nb_nodes x 3]
    const Array<Real> &gaps = contact.getGaps();       // [nb_nodes x 1]（非 slave 結點通常為 0）
    const Array<Real> &normals = contact.getNormals(); // [nb_nodes x 3]（非 slave 結點通常為 0）

    // slip-weakening 參數（請依需求調整）
    const Real mu_s = 0.60;   // 初始/靜摩擦
    const Real mu_d = 0.20;   // 殘餘/動摩擦
    const Real Dc = 0.20e-3;  // 臨界滑移 (m)
    const Real epsN = 1.0e10; // 要與 material.dat 的 epsilon_n 一致
    const Real Aeff = 1.0e-6; // 節點等效面積（可日後精算成每節點不同）

    // 為每個節點建立 slip 累積量（與 gaps 尺寸一致）
    Array<Real> slip(gaps.size(), 1);
    slip.set(0.);

    // 初始條件
    vel.set(0.);
    disp.set(0.);
    solid.getBlockedDOFs().set(false);


    // Boundary Conditions Here:
    // I want: stationary-block-right   to be fixed on Y
    //         stationary-block-back    to be fixed on X
    //         moving-block-bottom      to be fixed on Z
    //         stationary-block-bottom  to be fixed on Z
    //
    //         A 16MPa acting on moving-block-front on +X
    //         A 16MPa acting on moving-block-left  on +Y


    // Add normal and shearing stress
    Vector<Real, 3> t_front{16e6, 0., 0.}; // +X
    Vector<Real, 3> t_left{0., 16e16, 0.};  // +Y
    solid.applyBC(BC::Neumann::FromTraction(t_front), "moving-block-front");
    solid.applyBC(BC::Neumann::FromTraction(t_left), "moving-block-left");

    // Fix some surfaces
    solid.applyBC(BC::Dirichlet::FixedValue(0., _y), "stationary-block-right");
    solid.applyBC(BC::Dirichlet::FixedValue(0., _x), "stationary-block-back");
    solid.applyBC(BC::Dirichlet::FixedValue(0., _z), "moving-block-bottom");
    solid.applyBC(BC::Dirichlet::FixedValue(0., _z), "stationary-block-bottom");

    const Int max_steps = 200000;

    for (Int s = 0; s < max_steps; ++s)
    {
        force.set(0.);

        Real sum_vt = 0.;
        Int cnt_vt = 0;

        const Int N = gaps.size();

        for (Int i = 0; i < N; ++i)
        {
            const Real g = gaps(i);
            if (g <= 0.)
                continue;

            const Real vx = vel(i, 0), vy = vel(i, 1), vz = vel(i, 2);
            const Real nx = normals(i, 0), ny = normals(i, 1), nz = normals(i, 2);

            // vn = v · n
            const Real vn = vx * nx + vy * ny + vz * nz;
            // vt = v - vn * n
            const Real vtx = vx - vn * nx;
            const Real vty = vy - vn * ny;
            const Real vtz = vz - vn * nz;
            const Real vtn = std::sqrt(vtx * vtx + vty * vty + vtz * vtz);

            slip(i) += vtn * dt;

            const Real mu_now = mu_d + (mu_s - mu_d) * std::max(Real(0.), Real(1.) - slip(i) / Dc);

            const Real Fn = epsN * g * Aeff;
            const Real Ft_max = mu_now * Fn;

            if (vtn > 0.)
            {
                force(i, 0) += (-Ft_max / vtn) * vtx;
                force(i, 1) += (-Ft_max / vtn) * vty;
                force(i, 2) += (-Ft_max / vtn) * vtz;
            }

            sum_vt += vtn;
            ++cnt_vt;
        }

        coupler.solveStep();

        if (s % 100 == 0)
        {
            Real slip_avg = 0.;
            for (Int i = 0; i < N; ++i)
                slip_avg += slip(i);
            if (N)
                slip_avg /= N;

            const Real vbar = (cnt_vt ? sum_vt / cnt_vt : 0.);
            std::cout << "step " << s << "  vbar=" << vbar << "  <slip>=" << slip_avg << std::endl;
            // Maybe: coupler.dump();
        }
    }

    finalize();
    return 0;
}