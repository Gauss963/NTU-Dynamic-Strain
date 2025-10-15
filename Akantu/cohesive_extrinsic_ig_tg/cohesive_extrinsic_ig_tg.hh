#ifndef COHESIVE_EXTRINSIC_IG_TG_HH
#define COHESIVE_EXTRINSIC_IG_TG_HH

#include "solid_mechanics_model_cohesive.hh"

using namespace akantu;

class Velocity : public BC::Dirichlet::DirichletFunctor {
public:
    explicit Velocity(SolidMechanicsModel & model, Real vel, BC::Axis ax = _x);

    inline void operator()(Idx node, VectorProxy<bool> & flags,
                            VectorProxy<Real> & disp,
                            const VectorProxy<const Real> & coord) override;

private:
    SolidMechanicsModel & model;
    Real vel, disp;
};

#endif