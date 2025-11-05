#include "solid_mechanics_model_cohesive.hh"
#include "mesh.hh"
#include "aka_common.hh"
#include <iostream>

using namespace akantu;

int main(int argc, char *argv[])
{
    constexpr Int sd = 3;
    const std::string mesh_file = "../../../Models/50mm-PMMA-CZM.msh";
    const std::string mat_file = "../../../Materials/material.dat";

    initialize(mat_file, argc, argv);
    std::cout << "Initialized" << std::endl;

    Mesh mesh(sd);
    mesh.read(mesh_file);
    std::cout << "Load files successful." << std::endl;

    std::cout << "Cells (3D): " << mesh.getNbElement(mesh.getSpatialDimension()) << std::endl;
    std::cout << "Faces (2D): " << mesh.getNbElement(mesh.getSpatialDimension() - 1) << std::endl;
    std::cout << "Edges (1D): " << mesh.getNbElement(1) << std::endl;

    SolidMechanicsModelCohesive model(mesh);

    MaterialCohesiveRules rules{
        {{"friction_master", "friction_slave"}, "interface_mat"},
        {{"friction_slave", "friction_master"}, "interface_mat"}};
    std::cout << "Got material" << std::endl;

    auto cohesive_selector = std::make_shared<MaterialCohesiveRulesSelector>(model, rules);
    auto bulk_selector = std::make_shared<MeshDataMaterialSelector<std::string>>("physical_names", model);
    std::cout << "Got physical names" << std::endl;

    cohesive_selector->setFallback(bulk_selector);
    bulk_selector->setFallback(model.getMaterialSelector());
    model.setMaterialSelector(cohesive_selector);
    std::cout << "Set material selector" << std::endl;

    
    model.initFull(_analysis_method = _explicit_lumped_mass, _is_extrinsic = true); 
    // This is where went WRONG!!


    std::cout << "After model initialization" << std::endl;

    Real dt = model.getStableTimeStep() * 0.5;
    model.setTimeStep(dt);
    std::cout << "dt = " << dt << std::endl;


    

    std::map<std::string, akantu::Int> material_counts;
    const int dim = mesh.getSpatialDimension();
    for (auto type : mesh.elementTypes(dim))
    {
        const auto &mat_by_el = model.getMaterialByElement(type);
        for (akantu::Idx e = 0; e < mat_by_el.size(); ++e)
        {
            akantu::Idx mid = mat_by_el(e);
            const auto &mat = model.getMaterial(mid);
            material_counts[mat.getName()]++;
        }
    }
    for (auto &&kv : material_counts)
    {
        std::cout << "Material \"" << kv.first << "\" assigned to "
                  << kv.second << " bulk elements.\n";
    }
    akantu::Int sum = 0;
    for (auto &&kv : material_counts)
        sum += kv.second;
    if (sum != mesh.getNbElement(dim))
    {
        std::cerr << "[WARN] material count != total 3D elements (" << sum
                  << " vs " << mesh.getNbElement(dim) << ")\n";
    }




    model.assembleMassLumped();
    // model.setBaseName("czm_debug");
    model.addDumpFieldVector("displacement");
    model.addDumpFieldVector("velocity");
    model.addDumpFieldVector("internal_force");
    model.addDumpField("stress");
    // model.addDumpField("material_index");
    model.addDumpField("grad_u");
    // model.addDumpField("cohesive_opening");
    // model.addDumpField("cohesive_traction");

    std::cout << "Before getVelocity" << std::endl;
    auto &vel = model.getVelocity();
    std::cout << "Before getDisplacement" << std::endl;
    auto &disp = model.getDisplacement();
    std::cout << "Before set velocity to 0" << std::endl;
    vel.set(0.);
    disp.set(0.);
    std::cout << "After setting vel and disp" << std::endl;

    Vector<Real, 3> t_front{32.0, 0.0, 0.0}; // MPa traction (+X)
    Vector<Real, 3> t_left{  0.0, 5.0, 0.0}; // MPa traction (+Y)

    model.applyBC(BC::Neumann::FromTraction(t_front), "moving-block-front");
    model.applyBC(BC::Neumann::FromTraction(t_left), "moving-block-left");

    model.applyBC(BC::Dirichlet::FixedValue(0., _y), "stationary-block-right");
    model.applyBC(BC::Dirichlet::FixedValue(0., _x), "stationary-block-back");

    std::cout << "set B.C. successful." << std::endl;

    const Int max_steps = 200000;
    // for (Int s = 0; s < max_steps; ++s)
    // {
    //     model.solveStep();

    //     if (s % 100 == 0)
    //     {
    //         std::cout << "Step: " << s << " / " << max_steps << std::endl;
    //         model.dump();
    //     }
    // }
    bool czm_fields_added = false;
    auto has_elem_type = [&](const akantu::Mesh &m, akantu::ElementType t) -> bool
    {
        const int sd = m.getSpatialDimension();
        for (int d = 0; d <= sd; ++d)
            for (auto et : m.elementTypes(d))
                if (et == t)
                    return true;
        return false;
    };
    for (Int s = 0; s < max_steps; ++s)
    {
        model.checkCohesiveStress();
        model.solveStep();

        if (!czm_fields_added)
        {
            int n6 = has_elem_type(mesh, _cohesive_3d_6) ? mesh.getNbElement(_cohesive_3d_6, _not_ghost) : 0;
            int n8 = has_elem_type(mesh, _cohesive_3d_8) ? mesh.getNbElement(_cohesive_3d_8, _not_ghost) : 0;
            if (n6 + n8 > 0)
            {
                model.addDumpField("cohesive_opening");
                model.addDumpField("cohesive_traction");
                model.addDumpField("material_index");
                czm_fields_added = true;
                std::cout << "[INFO] cohesive inserted: tri6=" << n6 << ", quad8=" << n8 << "\n";
            }
        }

        if (s % 100 == 0)
        {
            model.dump();
            std::cout << "Step: " << s << "\n";
        }
    }

    finalize();
    return 0;
}