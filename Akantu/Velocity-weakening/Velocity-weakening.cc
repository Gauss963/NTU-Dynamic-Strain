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

    constexpr Int sd = 3;

    const std::string MESHFILE = "../../../Models/50mm-PMMA.msh";
    const std::string MATERIALFILE = "../../../Materials/material.dat";

    initialize(MATERIALFILE, argc, argv);
    Mesh mesh(sd);
    mesh.read(MESHFILE);

    CouplerSolidContact coupler(mesh);
    auto &solid = coupler.getSolidMechanicsModel();
    auto &contact = coupler.getContactMechanicsModel();

    // 用 Gmsh 的 physical surfaces 當接觸面（mesh 需有 friction_master / friction_slave）
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
    const Array<Real> &gaps = contact.getGaps();       // [nb_nodes x 1]（非 slave 多為 0）
    const Array<Real> &normals = contact.getNormals(); // [nb_nodes x 3]（非 slave 多為 0）

    // slip-weakening 參數
    const Real mu_s = 0.60;   // 靜摩擦
    const Real mu_d = 0.20;   // 動摩擦
    const Real Dc = 0.20e-3;  // 臨界滑移 (m)
    const Real epsN = 1.0e10; // 要與 material.dat 的 epsilon_n 一致（單位自洽）
    const Real Aeff = 1.0e-6; // 接觸節點等效面積（先用常數，之後可精算）

    // 累積滑移量（與 gaps 尺寸一致）
    Array<Real> slip(gaps.size(), 1);
    slip.set(0.);

    // 初始條件
    vel.set(0.);
    disp.set(0.);
    solid.getBlockedDOFs().set(false);

    // 邊界條件（Neumann traction + Dirichlet 夾持）
    // 16 MPa traction
    Vector<Real, 3> t_front{16, 0., 0.}; // +X on moving-block-front
    Vector<Real, 3> t_left{0., 16, 0.};  // +Y on moving-block-left   ← 修正：16e6 (不是 16e16)
    solid.applyBC(BC::Neumann::FromTraction(t_front), "moving-block-front");
    solid.applyBC(BC::Neumann::FromTraction(t_left), "moving-block-left");

    // 固定面
    solid.applyBC(BC::Dirichlet::FixedValue(0., _y), "stationary-block-right");
    solid.applyBC(BC::Dirichlet::FixedValue(0., _x), "stationary-block-back");
    solid.applyBC(BC::Dirichlet::FixedValue(0., _z), "moving-block-bottom");
    solid.applyBC(BC::Dirichlet::FixedValue(0., _z), "stationary-block-bottom");

    const Int max_steps = 200000;

    // 自訂摩擦力的暫存（避免動到 external force 原本由 Neumann 組裝的部分）
    Array<Real> fric(force.size(), 3); // [nb_nodes x 3]
    const Real v_eps = 1e-12;

    for (Int s = 0; s < max_steps; ++s)
    {
        // 關鍵：不要清 external force！只清自己的暫存
        fric.set(0.);

        // 以「壓入才有接觸」的慣例：g < 0 視為接觸；Fn = max(0, -epsN * g) * Aeff
        const Int N = gaps.size();

        Real sum_vt = 0.;
        Int cnt_vt = 0;

        for (Int i = 0; i < N; ++i)
        {
            const Real g = gaps(i);
            if (g >= 0.)
                continue; // 分離(或非接觸點)就略過

            const Real nx = normals(i, 0), ny = normals(i, 1), nz = normals(i, 2);
            const Real vx = vel(i, 0), vy = vel(i, 1), vz = vel(i, 2);

            // 法向/切向分解
            const Real vn = vx * nx + vy * ny + vz * nz;
            const Real vtx = vx - vn * nx;
            const Real vty = vy - vn * ny;
            const Real vtz = vz - vn * nz;
            const Real vtn = std::sqrt(vtx * vtx + vty * vty + vtz * vtz);

            // 累積滑移
            slip(i) += vtn * dt;

            // 線性弱化
            const Real weakening = std::max(Real(0.), Real(1.) - slip(i) / Dc);
            const Real mu_now = mu_d + (mu_s - mu_d) * weakening;

            // 正向力（僅壓入）
            const Real Fn = std::max(0., -epsN * g) * Aeff;
            const Real Ft_max = mu_now * Fn;

            // 切向摩擦力方向 = - v_t / |v_t|
            if (vtn > v_eps && Ft_max > 0.)
            {
                fric(i, 0) += (-Ft_max / vtn) * vtx;
                fric(i, 1) += (-Ft_max / vtn) * vty;
                fric(i, 2) += (-Ft_max / vtn) * vtz;
            }

            sum_vt += vtn;
            ++cnt_vt;
        }

        // 疊加到 external force（不清 force，直接加）
        for (Int i = 0; i < force.size(); ++i)
        {
            force(i, 0) += fric(i, 0);
            force(i, 1) += fric(i, 1);
            force(i, 2) += fric(i, 2);
        }

        // 時間推進
        coupler.solveStep();

        // 監控輸出
        if (s % 2 == 0)
        {
            Real slip_avg = 0.;
            for (Int i = 0; i < N; ++i)
                slip_avg += slip(i);
            if (N)
                slip_avg /= N;

            // 也印個最大速度確認有在動
            Real vmax = 0.;
            for (Int i = 0; i < vel.size(); ++i)
            {
                const Real v = std::sqrt(vel(i, 0) * vel(i, 0) + vel(i, 1) * vel(i, 1) + vel(i, 2) * vel(i, 2));
                if (v > vmax)
                    vmax = v;
            }

            const Real vbar = (cnt_vt ? sum_vt / cnt_vt : 0.);
            std::cout << "step " << s
                      << "  vmax=" << vmax
                      << "  vbar=" << vbar
                      << "  <slip>=" << slip_avg
                      << std::endl;
        }
    }

    finalize();
    return 0;
}