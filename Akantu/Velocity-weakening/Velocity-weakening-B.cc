#include "aka_common.hh"
#include "coupler_solid_contact.hh"
#include <iostream>

// using namespace akantu;

int main(int argc, char **argv)
{
    const akantu::Int sd = 3;

    const std::string meshfile = "../../../Models/50mm-PMMA.msh";
    const std::string matfile = "../../../Materials/material.dat";

    akantu::initialize(matfile, argc, argv);

    akantu::Mesh mesh(sd);
    mesh.read(meshfile);

    // 建立耦合器：裡面同時有 SolidMechanicsModel 與 ContactMechanicsModel
    akantu::CouplerSolidContact coupler(mesh);
    auto &solid = coupler.getSolidMechanicsModel();
    auto &contact = coupler.getContactMechanicsModel();

    // 若用 Gmsh 的 physical_names 指派材料，建議掛 MaterialSelector
    auto selector = std::make_shared<akantu::MeshDataMaterialSelector<std::string>>("physical_names", solid);
    solid.setMaterialSelector(selector);

    // 告訴 Contact 要用「物理面」當作接觸面來源（對應你在 Gmsh 取的名字）
    auto surface_selector = std::make_shared<akantu::PhysicalSurfaceSelector>(mesh);
    contact.getContactDetector().setSurfaceSelector(surface_selector); // v5 官方作法
    // ↑ master/slave 的實際名稱由 material.dat 的 contact_detector 區塊決定

    // 初始化（顯式）
    coupler.initFull(akantu::_analysis_method = akantu::_explicit_lumped_mass);

    // time step（保守地拿固體穩定步長的 10%）
    akantu::Real dt = solid.getStableTimeStep() * 0.1;
    coupler.setTimeStep(dt);
    std::cout << "dt = " << dt << "\n";

    // 速度弱化參數（你可以改）
    const akantu::Real mu_s = 0.60;
    const akantu::Real mu_d = 0.20;
    const akantu::Real v1 = 0.50;
    const std::string iface = "friction_slave";

    // 取出需要的陣列
    auto &vel = solid.getVelocity();
    auto &gaps = contact.getGaps();       // > 0 表示在 active set（有「穿入」量） [oai_citation:4‡Akantu](https://akantu.readthedocs.io/en/latest/manual/contactmechanicsmodel.html?utm_source=chatgpt.com)
    auto &normals = contact.getNormals(); // slave 節點法向量（與 gaps 對齊） [oai_citation:5‡Akantu](https://akantu.readthedocs.io/en/latest/manual/contactmechanicsmodel.html?utm_source=chatgpt.com)

    // （可選）初始條件
    vel.set(0.);
    solid.getDisplacement().set(0.);
    solid.getBlockedDOFs().set(false);

    // 這裡省略外力/邊界條件：請在每步迴圈中自行施加位移或力到 solid
    // 例如：solid.applyBC(BC::Dirichlet::IncrementValue(...), "loading");

    const akantu::Int max_steps = 1000;
    for (akantu::Int s = 0; s < max_steps; ++s)
    {

        akantu::Real sum_vt = 0.;
        akantu::Int cnt = 0;

        for (auto &&tpl : zip(gaps, make_view(vel, sd), make_view(normals, sd)))
        {
            const akantu::Real gap = std::get<0>(tpl);
            const auto &v = std::get<1>(tpl);
            const auto &n = std::get<2>(tpl);
            if (gap > 0)
            {
                akantu::Real vn = v.dot(n);
                akantu::Vector<akantu::Real> vt = v - vn * n;
                sum_vt += vt.norm();
                ++cnt;
            }
        }
        const akantu::Real vbar = (cnt > 0) ? (sum_vt / cnt) : 0.;

        // === 速度弱化：更新目前的 μ ===
        akantu::Real mu_now = mu_d + (mu_s - mu_d) * std::max(0., 1. - vbar / v1);

        // v5 用 ParameterRegistry 路徑設定 contact_resolution 的參數
        // 兩種路徑字串擇一會生效（都呼叫不影響正確性）
        // contact.set("contact_resolution:friction_slave:mu", mu_now);
        contact.set("contact_resolution.friction_slave.mu", mu_now);

        // === 你的載入條件（位移/力）請放在 solveStep() 前後合適位置 ===
        // solid.applyBC(...); // ← 依你的案例加

        // 前進一步
        coupler.solveStep();

        if (s % 50 == 0)
        {
            std::cout << "step " << s << "  vbar=" << vbar << "  mu=" << mu_now << "\n";
            // 可選：coupler.dump();
        }
    }

    akantu::finalize();
    return 0;
}