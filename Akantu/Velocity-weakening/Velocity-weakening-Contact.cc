#include "solid_mechanics_model_cohesive.hh"
#include "coupler_solid_cohesive_contact.hh"
#include "mesh.hh"
#include "aka_common.hh"

#include <iostream>
#include <iomanip>
#include <chrono>
#include <memory>

int main(int argc, char *argv[])
{

    constexpr akantu::Int sd = 3;
    constexpr akantu::Real us = 1e-6;
    constexpr akantu::Real ms = 1e-3;
    constexpr akantu::Int PMMA_thickness = 500;
    constexpr akantu::Real time_factor = 0.5;

    const std::string mesh_file = "../../../Models/" + std::to_string(PMMA_thickness) + "mm-PMMA-CZM.msh";
    const std::string mat_file = "../../../Materials/material-mm-MPa.dat";
    const std::string SYS_PATH = "../../../../../../../../../../";
    const std::string DUMP_PATH = "Volumes/Gauss-T7/Dump_5mm/";

    akantu::initialize(mat_file, argc, argv);
    std::cout << "Initialized\n";

    akantu::Mesh mesh(sd);
    mesh.read(mesh_file);

    std::cout << "Load files successful.\n";
    std::cout << "Cells (3D): " << mesh.getNbElement(mesh.getSpatialDimension()) << "\n";
    std::cout << "Faces (2D): " << mesh.getNbElement(mesh.getSpatialDimension() - 1) << "\n";
    std::cout << "Edges (1D): " << mesh.getNbElement(1) << "\n";

    // 以「CZM×接觸」耦合器建立兩個 model：SolidMechanicsModelCohesive 與 ContactMechanicsModel
    akantu::CouplerSolidCohesiveContact coupler(mesh);
    auto &solid = coupler.getSolidMechanicsModelCohesive();
    auto &contact = coupler.getContactMechanicsModel();

    // 材料選擇（bulk 與 cohesive 規則仍讀自 material.dat）
    auto bulk_selector = std::make_shared<akantu::MeshDataMaterialSelector<std::string>>("physical_names", solid);
    solid.setMaterialSelector(bulk_selector);

    // 接觸面選擇：本題要在兩個 Vol 交界（你的 friction_master/slave）與可能的斷裂面上都偵測
    // 若僅限物理面：PhysicalSurfaceSelector(mesh)
    // 若要同時支援物理面與斷裂後的 cohesive 面：AllSurfaceSelector(mesh.getMeshFacets())
    auto surface_selector = std::make_shared<akantu::AllSurfaceSelector>(mesh.getMeshFacets());
    contact.getContactDetector().setSurfaceSelector(surface_selector);

    // 初始化（顯式 + extrinsic）
    coupler.initFull(akantu::_analysis_method = akantu::_explicit_lumped_mass, akantu::_is_extrinsic = true);
    std::cout << "After model initialization\n";


    akantu::Real dt = solid.getStableTimeStep() * time_factor;
    coupler.setTimeStep(dt);
    std::cout << "dt = " << dt << " s (" << dt / us << " us)\n";


    coupler.setBaseName(SYS_PATH + DUMP_PATH + std::to_string(PMMA_thickness) + "mm-PMMA-CZM-contact");
    coupler.addDumpFieldVector("displacement");
    coupler.addDumpFieldVector("velocity");
    coupler.addDumpFieldVector("internal_force");
    coupler.addDumpField("stress");
    coupler.addDumpField("grad_u");
    coupler.addDumpField("gaps");
    coupler.addDumpField("areas");
    coupler.addDumpFieldVector("normal_force");
    coupler.addDumpFieldVector("external_force");

    auto &vel = solid.getVelocity();
    auto &disp = solid.getDisplacement();
    vel.set(0.0);
    disp.set(0.0);

    
    const akantu::Real normalStress = 8.0; // MPa
    const akantu::Real targetDisp = 3.0;   // mm 在 y 方向
    const akantu::Real riseEnd = 0.70;     // 步入式位移比


    akantu::Vector<akantu::Real, 3> t_front{normalStress, 0.0, 0.0}; // (+X) traction
    solid.applyBC(akantu::BC::Neumann::FromTraction(t_front), "moving-block-front");
    solid.applyBC(akantu::BC::Dirichlet::FixedValue(0., akantu::_y), "stationary-block-right");
    solid.applyBC(akantu::BC::Dirichlet::FixedValue(0., akantu::_x), "stationary-block-back");

    std::cout << "set B.C. successful.\n";

    const akantu::Real SIMULATION_TIME = 0.4 * ms;
    const akantu::Int MAX_STEPS = static_cast<akantu::Int>(std::ceil(SIMULATION_TIME / dt));
    const akantu::Int DUMP_INTERVAL = 3;

    std::cout << "Starting time integration for " << (SIMULATION_TIME / ms)
              << " ms (" << MAX_STEPS << " steps)\n";

    auto t0 = std::chrono::high_resolution_clock::now();

    for (akantu::Int s = 0; s < MAX_STEPS; ++s)
    {
        const akantu::Real progress = static_cast<akantu::Real>(s) / MAX_STEPS;
        const akantu::Real alpha = (progress < riseEnd) ? (progress / riseEnd) : 1.0;
        const akantu::Real pushY = alpha * targetDisp;

        solid.applyBC(akantu::BC::Dirichlet::FixedValue(pushY, akantu::_y), "moving-block-left");
        coupler.solveStep();

        solid.checkCohesiveStress();

        if (s % DUMP_INTERVAL == 0)
            coupler.dump();



        auto t1 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> el = t1 - t0;
        double time_per_iter = (s > 0) ? el.count() / s : 0.0;
        double estimated_total = time_per_iter * MAX_STEPS;
        double remaining = estimated_total - el.count();

        std::cout << "Step " << std::setw(8) << s << "/" << MAX_STEPS
                  << " | Elapsed: " << std::fixed << std::setprecision(1) << std::setw(8) << el.count() << " s"
                  << " | ETA: " << std::fixed << std::setprecision(1) << std::setw(8) << remaining << " s = "
                  << std::fixed << std::setprecision(1) << std::setw(5) << remaining / 3600.0 << " h\r";
        std::cout.flush();
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> total = t1 - t0;
    std::cout << "\nTotal elapsed time: " << total.count() << " s = " << total.count() / 3600.0 << " h\n";

    akantu::finalize();
    return 0;
}