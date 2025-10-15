#include "cohesive_extrinsic_ig_tg.hh"
#include "material_elastic.hh"
#include "material_cohesive_linear.hh"
#include <iostream>
#include <omp.h>

using namespace akantu;

// === 自定義 Velocity functor ===
Velocity::Velocity(SolidMechanicsModel & model, Real vel, BC::Axis ax)
    : DirichletFunctor(ax), model(model), vel(vel) {
    disp = vel * model.getTimeStep();
}

inline void Velocity::operator()(Idx node, VectorProxy<bool> &, VectorProxy<Real> & disp, const VectorProxy<const Real> & coord) {
    Real sign = std::signbit(coord(axis)) ? -1. : 1.;
    disp(axis) += sign * this->disp;
    model.getVelocity()(node, axis) = sign * vel;
}

int main(int argc, char * argv[]) {
    omp_set_num_threads(16);
    std::cout << "Using " << omp_get_max_threads() << " OpenMP threads.\n";

    const Int spatial_dimension = 2;
    const Int max_steps = 10000;

    // === 1. 讀取網格 ===
    Mesh mesh(spatial_dimension);
    mesh.read("square.msh");

    if (mesh.getNbNodes() == 0) {
        std::cerr << "ERROR: Mesh has 0 nodes. Did you forget to place 'square.msh' in the build directory?\n";
        return 1;
    }
    std::cout << "Mesh read succeeded.\n";

    // === 2. 建立模型 ===
    SolidMechanicsModelCohesive model(mesh);

    // === 3. 設定材料選擇器（含 fallback）===
    MaterialCohesiveRules rules{
        {{"btop", "bbottom"}, "tg_cohesive"},
        {{"btop", "btop"},     "ig_cohesive"},
        {{"bbottom", "bbottom"}, "ig_cohesive"},
        // fallback 給所有沒有明確 match 的 interface
        { {".*", ".*"}, "default" }
    };

    auto cohesive_selector = std::make_shared<MaterialCohesiveRulesSelector>(model, rules);
    auto bulk_selector = std::make_shared<MeshDataMaterialSelector<std::string>>("material.dat", model);

    cohesive_selector->setFallback(bulk_selector);
    bulk_selector->setFallback(model.getMaterialSelector());
    model.setMaterialSelector(cohesive_selector);
    // === 4. 初始化 ===
    model.initFull(_analysis_method = _explicit_lumped_mass, _is_extrinsic = true);

    Real time_step = model.getStableTimeStep() * 0.05;
    model.setTimeStep(time_step);
    std::cout << "Time step: " << time_step << "\n";

    model.assembleMassLumped();

    // === 5. Apply boundary conditions ===
    auto & position = mesh.getNodes();
    auto & velocity = model.getVelocity();

    model.applyBC(BC::Dirichlet::FlagOnly(_y), "top");
    model.applyBC(BC::Dirichlet::FlagOnly(_y), "bottom");
    model.applyBC(BC::Dirichlet::FlagOnly(_x), "left");
    model.applyBC(BC::Dirichlet::FlagOnly(_x), "right");

    // === 6. 輸出欄位 ===
    model.setBaseName("extrinsic");
    model.addDumpFieldVector("displacement");
    model.addDumpField("velocity");
    model.addDumpField("acceleration");
    model.addDumpField("internal_force");
    model.addDumpField("stress");
    model.addDumpField("grad_u");
    model.addDumpField("material_index");
    model.dump();

    // === 7. 初始速度設定 ===
    Real loading_rate = 0.1;
    Real VI = loading_rate * 2 * 0.5;
    for (auto && [pos, vel] : zip(make_view(position, spatial_dimension),
                                  make_view(velocity, spatial_dimension))) {
        vel = loading_rate * pos;
    }

    model.dump();

    // === 8. 動態施加 BC（velocity）===
    Velocity vely(model, VI, _y);
    Velocity velx(model, VI, _x);

    for (Int s = 1; s <= max_steps; ++s) {
        model.applyBC(vely, "top");
        model.applyBC(vely, "bottom");
        model.applyBC(velx, "left");
        model.applyBC(velx, "right");

        model.checkCohesiveStress();
        model.solveStep();

        if (s % 10 == 0) {
            model.dump();
            std::cout << "passing step " << s << "/" << max_steps << "\n";
        }
    }

    return 0;
}