#include "solid_mechanics_model_cohesive.hh"
#include "mesh.hh"
#include "aka_common.hh"
#include <iostream>

using namespace akantu;

int main(int argc, char *argv[])
{
    constexpr Int sd = 3;

    const std::string mesh_file = "../../../Models/50mm-PMMA.msh";
    const std::string mat_file = "../../../Materials/material.dat";

    std::cout << "Load files successful." << std::endl;

    std::cout << "Before initialize" << std::endl;
    initialize(mat_file, argc, argv);
    std::cout << "After initialize" << std::endl;

    Mesh mesh(sd);
    mesh.read(mesh_file);
    std::cout << "Mesh spatial dimension: " << mesh.getSpatialDimension() << std::endl;
    std::cout << "Number of elements (bulk): " << mesh.getNbElement(0) << std::endl;

    SolidMechanicsModelCohesive model(mesh);

    MaterialCohesiveRules rules{
        {{"friction_master", "friction_slave"}, "interface_mat"}};
    std::cout << "After getting material" << std::endl;

    auto cohesive_selector = std::make_shared<MaterialCohesiveRulesSelector>(model, rules);
    auto bulk_selector = std::make_shared<MeshDataMaterialSelector<std::string>>("physical_names", model);

    cohesive_selector->setFallback(bulk_selector);
    bulk_selector->setFallback(model.getMaterialSelector());
    model.setMaterialSelector(cohesive_selector);

    std::map<std::string, Int> material_counts;

    for (const auto &[name, count] : material_counts)
    {
        std::cout << "Material \"" << name << "\" assigned to " << count << " elements." << std::endl;
    }

    model.initFull(_analysis_method = _explicit_lumped_mass, _is_extrinsic = true);

    Real dt = model.getStableTimeStep() * 0.5;
    model.setTimeStep(dt);
    std::cout << "dt = " << dt << std::endl;

    model.setBaseName("czm_debug");
    model.addDumpFieldVector("displacement");
    model.addDumpFieldVector("velocity");
    model.addDumpFieldVector("internal_force");
    model.addDumpField("stress");
    model.addDumpField("material_index");
    model.addDumpField("grad_u");
    model.addDumpField("cohesive_opening");
    model.addDumpField("cohesive_traction");

    // auto &vel = model.getVelocity();
    // auto &disp = model.getDisplacement();
    // vel.set(0.);
    // disp.set(0.);

    std::cout << "Before getVelocity" << std::endl;
    auto &vel = model.getVelocity();
    std::cout << "Before getDisplacement" << std::endl;
    auto &disp = model.getDisplacement();
    std::cout << "Before set velocity to 0" << std::endl;
    vel.set(0.);
    disp.set(0.);
    std::cout << "After setting vel and disp" << std::endl;

    Vector<Real, 3> t_front{2.0, 0.0, 0.0}; // MPa traction (+X)
    Vector<Real, 3> t_left{0.0, 5.0, 0.0};  // MPa traction (+Y)

    model.applyBC(BC::Neumann::FromTraction(t_front), "moving-block-front");
    model.applyBC(BC::Neumann::FromTraction(t_left), "moving-block-left");

    model.applyBC(BC::Dirichlet::FixedValue(0., _y), "stationary-block-right");
    model.applyBC(BC::Dirichlet::FixedValue(0., _x), "stationary-block-back");

    std::cout << "set B.C. successful." << std::endl;

    const Int max_steps = 200000;
    for (Int s = 0; s < max_steps; ++s)
    {
        model.solveStep();

        if (s % 100 == 0)
        {
            std::cout << "Step: " << s << " / " << max_steps << std::endl;
            model.dump();
        }
    }

    finalize();
    return 0;
}