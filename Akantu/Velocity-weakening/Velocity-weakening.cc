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
    // auto surf_sel = std::make_shared<PhysicalSurfaceSelector>(mesh);
    // contact.getContactDetector().setSurfaceSelector(surf_sel);

    // 顯式初始化
    coupler.initFull(_analysis_method = _explicit_lumped_mass);
    auto surf_sel = std::make_shared<PhysicalSurfaceSelector>(mesh);
    contact.getContactDetector().setSurfaceSelector(surf_sel);

    coupler.setBaseName("contact-debug");
    coupler.addDumpField("gaps");
    coupler.addDumpFieldVector("normal_force");
    coupler.addDumpFieldVector("external_force");
    coupler.addDumpFieldVector("displacement");

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

    // ==== 單位：mm–MPa ====
    // slip-weakening 參數（線性弱化）
    const Real mu_s = 0.60; // 靜摩擦
    const Real mu_d = 0.20; // 動摩擦
    const Real Dc = 0.20;   // 臨界滑移 (mm)  ← mm–MPa 制，請用 mm
    const Real Aeff = 1.0;  // 節點等效面積 (mm^2) 先給 1 mm^2 方便測試
    // 懲罰法正向剛度：epsN ≈ E / h （MPa/mm）
    // 若 PMMA E≈3000 MPa、元素尺寸 h≈1 mm，則 epsN~3000 MPa/mm 是合理級距
    const Real epsN = 3.0e3; // MPa/mm 先給一個保守值，之後可依網格調整

    // 累積滑移量
    Array<Real> slip(gaps.size(), 1);
    slip.set(0.);

    // 初始條件
    vel.set(0.);
    disp.set(0.);
    solid.getBlockedDOFs().set(false);

    // 邊界條件（Neumann traction + Dirichlet 夾持）
    Vector<Real, 3> t_front {2.0, 0.0, 0.0}; // +X on moving-block-front,  16 MPa
    Vector<Real, 3> t_left  {0.0, 5.0, 0.0}; // +Y on moving-block-left,   16 MPa
    solid.applyBC(BC::Neumann::FromTraction(t_front), "moving-block-front");
    solid.applyBC(BC::Neumann::FromTraction(t_left), "moving-block-left");

    // 固定面Ｆ
    solid.applyBC(BC::Dirichlet::FixedValue(0., _y), "stationary-block-right");
    solid.applyBC(BC::Dirichlet::FixedValue(0., _x), "stationary-block-back");
    // solid.applyBC(BC::Dirichlet::FixedValue(0., _z), "moving-block-bottom");
    // solid.applyBC(BC::Dirichlet::FixedValue(0., _z), "stationary-block-bottom");

    const Int max_steps = 200000;

    // 自訂摩擦力的暫存（避免動到 external force 原本由 Neumann 組裝的部分）
    Array<Real> fric(force.size(), 3); // [nb_nodes x 3]
    const Real v_eps = 1e-12;

    // 自動偵測 gap 號誌：第一個步驟統計 g<0 與 g>0 哪個多，就用那個當「壓入」
    bool penetration_is_negative = true; // default
    bool decided_sign = false;

    for (Int s = 0; s < max_steps; ++s)
    {
        fric.set(0.); // 只清自己的暫存，不要清 external force

        // 第 0 步做 gap 號誌偵測 + 基本診斷
        if (!decided_sign)
        {
            Int cnt_neg = 0, cnt_pos = 0, cnt_nrm = 0;
            Real gmin = +1e30, gmax = -1e30;
            const Int N0 = gaps.size();
            for (Int i = 0; i < N0; ++i)
            {
                const Real g = gaps(i);
                if (g < gmin)
                    gmin = g;
                if (g > gmax)
                    gmax = g;
                if (g < 0)
                    ++cnt_neg;
                else if (g > 0)
                    ++cnt_pos;

                const Real nx = normals(i, 0), ny = normals(i, 1), nz = normals(i, 2);
                const Real n2 = nx * nx + ny * ny + nz * nz;
                if (n2 > 1e-20)
                    ++cnt_nrm;
            }
            penetration_is_negative = (cnt_neg >= cnt_pos);
            decided_sign = true;

            std::cout << "[contact diag] gmin=" << gmin
                      << " gmax=" << gmax
                      << "  (#g<0,#g>0,#|n|>0)=("
                      << cnt_neg << "," << cnt_pos << "," << cnt_nrm << ")  "
                      << "penetration_is_negative=" << (penetration_is_negative ? "true" : "false")
                      << std::endl;
        }

        const auto in_contact = [&](Real g) -> bool
        {
            return penetration_is_negative ? (g < 0.) : (g > 0.);
        };

        // 接觸摩擦（線性弱化）
        const Int N = gaps.size();
        Real sum_vt = 0.;
        Int cnt_vt = 0;

        for (Int i = 0; i < N; ++i)
        {
            const Real g = gaps(i);
            if (!in_contact(g))
                continue; // 只有「壓入」視為接觸

            const Real nx = normals(i, 0), ny = normals(i, 1), nz = normals(i, 2);
            const Real vx = vel(i, 0), vy = vel(i, 1), vz = vel(i, 2);

            const Real vn = vx * nx + vy * ny + vz * nz;
            const Real vtx = vx - vn * nx;
            const Real vty = vy - vn * ny;
            const Real vtz = vz - vn * nz;
            const Real vtn = std::sqrt(vtx * vtx + vty * vty + vtz * vtz);

            // 累積滑移（mm）
            slip(i) += vtn * dt;

            // 線性弱化 μ(s) = μ_d + (μ_s − μ_d) * max(0, 1 − s/Dc)
            const Real weakening = std::max(Real(0.), Real(1.) - slip(i) / Dc);
            const Real mu_now = mu_d + (mu_s - mu_d) * weakening;

            // 正向力（只在壓入）：Fn = epsN * |g| * Aeff
            const Real Fn = epsN * std::abs(g) * Aeff; // (MPa/mm * mm * mm^2) = MPa*mm^2 = N
            const Real Ft_max = mu_now * Fn;

            if (vtn > v_eps && Ft_max > 0.)
            {
                fric(i, 0) += (-Ft_max / vtn) * vtx;
                fric(i, 1) += (-Ft_max / vtn) * vty;
                fric(i, 2) += (-Ft_max / vtn) * vtz;
            }

            sum_vt += vtn;
            ++cnt_vt;
        }

        // 疊加到 external force（保留 Neumann）
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
            Real slip_avg = 0., vmax = 0., umax = 0.;
            for (Int i = 0; i < N; ++i)
                slip_avg += slip(i);
            if (N)
                slip_avg /= N;

            for (Int i = 0; i < vel.size(); ++i)
            {
                const Real v = std::sqrt(vel(i, 0) * vel(i, 0) + vel(i, 1) * vel(i, 1) + vel(i, 2) * vel(i, 2));
                if (v > vmax)
                    vmax = v;
                const Real u = std::sqrt(disp(i, 0) * disp(i, 0) + disp(i, 1) * disp(i, 1) + disp(i, 2) * disp(i, 2));
                if (u > umax)
                    umax = u;
            }

            const Real vbar = (cnt_vt ? sum_vt / cnt_vt : 0.);
            std::cout << "step " << s
                      << "  vmax=" << vmax
                      << "  umax=" << umax
                      << "  vbar(contact)=" << vbar
                      << "  <slip>=" << slip_avg
                      << "  (cnt_contact=" << cnt_vt << ")"
                      << std::endl;
        }
        solid.dump();
        contact.dump();
    }

    finalize();
    return 0;
}