#include "solid_mechanics_model_cohesive.hh"
#include <iostream>

using namespace akantu;

/* -------------------------------------------------------------------------- */
class Velocity : public BC::Dirichlet::DirichletFunctor
{
public:
    explicit Velocity(SolidMechanicsModel &model, Real vel, BC::Axis ax = _x)
        : DirichletFunctor(ax), model(model), vel(vel)
    {
        disp = vel * model.getTimeStep();
    }

public:
    inline void operator()(Idx node, VectorProxy<bool> & /*flags*/,
                           VectorProxy<Real> &disp,
                           const VectorProxy<const Real> &coord) override
    {
        Real sign = std::signbit(coord(axis)) ? -1. : 1.;
        disp(axis) += sign * this->disp;
        model.getVelocity()(node, axis) = sign * vel;
    }

private:
    SolidMechanicsModel &model;
    Real vel, disp;
};

/* -------------------------------------------------------------------------- */
int main(int argc, char *argv[])
{
    initialize("material.dat", argc, argv);

    const Int spatial_dimension = 2;
    const Int max_steps = 1000;

    Mesh mesh(spatial_dimension);
    mesh.read("square.msh");

    SolidMechanicsModelCohesive model(mesh);
    MaterialCohesiveRules rules{{{"btop", "bbottom"}, "tg_cohesive"},
                                {{"btop", "btop"}, "ig_cohesive"},
                                {{"bbottom", "bbottom"}, "ig_cohesive"}};

    /// model initialization
    auto cohesive_material_selector =
        std::make_shared<MaterialCohesiveRulesSelector>(model, rules);
    auto bulk_material_selector =
        std::make_shared<MeshDataMaterialSelector<std::string>>("physical_names",
                                                                model);
    auto &&current_selector = model.getMaterialSelector();

    cohesive_material_selector->setFallback(bulk_material_selector);
    bulk_material_selector->setFallback(current_selector);

    model.setMaterialSelector(cohesive_material_selector);

    model.initFull(_analysis_method = _explicit_lumped_mass,
                   _is_extrinsic = true);

    Real time_step = model.getStableTimeStep() * 0.05;
    model.setTimeStep(time_step);
    std::cout << "Time step: " << time_step << "\n";

    model.assembleMassLumped();

    auto &position = mesh.getNodes();
    auto &velocity = model.getVelocity();

    model.applyBC(BC::Dirichlet::FlagOnly(_y), "top");
    model.applyBC(BC::Dirichlet::FlagOnly(_y), "bottom");

    model.applyBC(BC::Dirichlet::FlagOnly(_x), "left");
    model.applyBC(BC::Dirichlet::FlagOnly(_x), "right");

    model.setBaseName("extrinsic");
    model.addDumpFieldVector("displacement");
    model.addDumpField("velocity");
    model.addDumpField("acceleration");
    model.addDumpField("internal_force");
    model.addDumpField("stress");
    model.addDumpField("grad_u");
    model.addDumpField("material_index");
    model.dump();

    /// initial conditions
    Real loading_rate = 0.1;
    // bar_height  = 2
    Real VI = loading_rate * 2 * 0.5;
    for (auto &&[pos, vel] : zip(make_view(position, spatial_dimension),
                                 make_view(velocity, spatial_dimension)))
    {
        vel = loading_rate * pos;
    }

    model.dump();

    Velocity vely(model, VI, _y);
    Velocity velx(model, VI, _x);

    /// Main loop
    for (Int s = 1; s <= max_steps; ++s)
    {

        model.applyBC(vely, "top");
        model.applyBC(vely, "bottom");

        model.applyBC(velx, "left");
        model.applyBC(velx, "right");

        model.checkCohesiveStress();

        model.solveStep();

        if (s % 10 == 0)
        {
            model.dump();
            std::cout << "passing step " << s << "/" << max_steps << "\n";
        }
    }

    return 0;
}
