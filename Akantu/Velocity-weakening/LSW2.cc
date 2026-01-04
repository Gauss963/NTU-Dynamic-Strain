#include "aka_common.hh"
#include "mesh.hh"
#include "solid_mechanics_model.hh"

#include "ntn_contact.hh"
#include "ntn_friction.hh"
#include "ntn_friclaw_linear_slip_weakening.hh"
#include "ntn_fricreg_no_regularisation.hh"

#include <iostream>
#include <iomanip>
#include <chrono>
#include <cmath>
#include <type_traits>
#include <utility>

#ifdef AKANTU_USE_MPI
#include <mpi.h>
#endif

using namespace akantu;

// -----------------------------
// Utilities: C++17 detection
// -----------------------------
namespace detail
{

    template <class, class = void>
    struct has_getGlobalFrictionTraction : std::false_type
    {
    };

    template <class T>
    struct has_getGlobalFrictionTraction<T, std::void_t<decltype(std::declval<T &>().getGlobalFrictionTraction())>>
        : std::true_type
    {
    };

    template <class, class = void>
    struct has_getGlobalFrictionTractions : std::false_type
    {
    };

    template <class T>
    struct has_getGlobalFrictionTractions<T, std::void_t<decltype(std::declval<T &>().getGlobalFrictionTractions())>>
        : std::true_type
    {
    };

    // add src nodal vector into dst external force
    template <class NodalArrayLike>
    inline void addNodal(Array<Real> &dst_ext, const NodalArrayLike &src, Int dim)
    {
        const Idx nmax = std::min<Idx>(dst_ext.size(), src.size());
        for (Idx n = 0; n < nmax; ++n)
        {
            for (Int d = 0; d < dim; ++d)
            {
                dst_ext(n, d) += static_cast<Real>(src(n, d));
            }
        }
    }

    template <class NodalArrayLike>
    inline void subNodal(Array<Real> &dst_ext, const NodalArrayLike &src, Int dim)
    {
        const Idx nmax = std::min<Idx>(dst_ext.size(), src.size());
        for (Idx n = 0; n < nmax; ++n)
        {
            for (Int d = 0; d < dim; ++d)
            {
                dst_ext(n, d) -= static_cast<Real>(src(n, d));
            }
        }
    }

} // namespace detail

int main(int argc, char *argv[])
{
    constexpr Int sd = 3;
    constexpr Real us = 1e-6;
    constexpr Real ms = 1e-3;

    constexpr Real time_factor = 0.05;
    constexpr Real SIMULATION_TIME = 20 * ms;

    constexpr int PMMA_thickness = 50;

    const std::string mesh_file =
        "../../../Models/" + std::to_string(PMMA_thickness) + "mm-BS-PMMA.msh";
    const std::string mat_file = "../../../Materials/NTN-LSW.dat";

    akantu::initialize(mat_file, argc, argv);

    const auto &comm = akantu::Communicator::getStaticCommunicator();
    const Int prank = comm.whoAmI();

    // -----------------------------
    // Mesh (read on rank0 + distribute)
    // -----------------------------
    akantu::Mesh mesh(sd);
    if (prank == 0)
    {
        mesh.read(mesh_file);
    }
    mesh.distribute();

    // -----------------------------
    // Solid model
    // -----------------------------
    akantu::SolidMechanicsModel model(mesh);

    auto bulk_selector =
        std::make_shared<akantu::MeshDataMaterialSelector<std::string>>("physical_names", model);
    model.setMaterialSelector(bulk_selector);

    model.initFull(akantu::_analysis_method = akantu::_explicit_lumped_mass);
    model.assembleMassLumped();

    const Real dt = model.getStableTimeStep() * time_factor;
    model.setTimeStep(dt);

    if (prank == 0)
    {
        std::cout << "dt = " << dt << " s (" << dt / us << " us)\n";
    }

    model.setBaseName(std::to_string(PMMA_thickness) + "mm-PMMA-NTN-LSW");
    model.addDumpFieldVector("displacement");
    model.addDumpFieldVector("velocity");
    model.addDumpFieldVector("internal_force");
    model.addDumpFieldVector("external_force");
    model.addDumpField("stress");
    model.addDumpField("grad_u");

    model.getVelocity().set(0.);
    model.getDisplacement().set(0.);

    // -----------------------------
    // NTN contact + friction (acts like "contact model" inside a coupler)
    // -----------------------------
    akantu::NTNContact ntn_contact(model, "contact");

    // 你原本的 pair（請確保方向/法向定義正確）
    // ntn_contact.addSurfacePair("stationary-block-back", "moving-block-front", akantu::_x);
    ntn_contact.addSurfacePair("moving-block-front", "stationary-block-back", akantu::_x);

    // Parallel init MUST come after mesh.distribute() and model.initFull()
    ntn_contact.initParallel();
    ntn_contact.updateInternalData();

    using Regular = akantu::NTNFricRegNoRegularisation;
    using Friction = akantu::NTNFriction<akantu::NTNFricLawLinearSlipWeakening, Regular>;
    Friction ntn_friction(ntn_contact, "friction");

    // -----------------------------
    // Quick sanity check of contact nodes (after initParallel + updateInternalData)
    // -----------------------------
    int local_pairs = static_cast<int>(ntn_contact.getNbContactNodes());
    int global_pairs = local_pairs;

#ifdef AKANTU_USE_MPI
    // NOTE: usually OK; if your Akantu communicator is not MPI_COMM_WORLD,
    // you can replace this with the correct communicator in your environment.
    MPI_Allreduce(&local_pairs, &global_pairs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif

    if (prank == 0)
    {
        std::cout << "[CHK] contact_pairs local(rank0)=" << local_pairs
                  << " global(sum)=" << global_pairs << "\n";
    }

    // -----------------------------
    // Boundary conditions (yours)
    // -----------------------------
    const Real normalStress = 8.0; // MPa
    const Real shearDisp = 8.0;    // mm
    const Real riseEnd = 0.30;

    Vector<Real, 3> t_normal{normalStress, 0.0, 0.0};
    model.applyBC(BC::Neumann::FromTraction(t_normal), "moving-block-back");

    model.applyBC(BC::Dirichlet::FixedValue(0., _x), "stationary-block-front");
    model.applyBC(BC::Dirichlet::FixedValue(0., _y), "stationary-block-left");

    // -----------------------------
    // Time loop controls
    // -----------------------------
    const Int TOTAL_FRAMES = 2400;
    const Int MAX_STEPS = static_cast<Int>(std::ceil(SIMULATION_TIME / dt));
    const Int DUMP_INTERVAL = std::max<Int>(1, MAX_STEPS / TOTAL_FRAMES);

    const Int rise_steps = std::max<Int>(1, static_cast<Int>(std::ceil(riseEnd * MAX_STEPS)));
    const Real dy = shearDisp / rise_steps;

    if (prank == 0)
    {
        std::cout << "Starting time integration for " << (SIMULATION_TIME / ms)
                  << " ms (" << MAX_STEPS << " steps)\n";
    }

    auto start_time = std::chrono::high_resolution_clock::now();

    // -----------------------------
    // "Coupler-like" loop:
    //  - update contact data
    //  - assemble contact pressure (+ friction traction if available)
    //  - temporarily add to model external force
    //  - solveStep
    //  - remove contact contributions
    // -----------------------------
    auto &f_ext = model.getExternalForce();
    const Int dim = model.getSpatialDimension();

    for (Int s = 0; s < MAX_STEPS; ++s)
    {
        // loading: ramp then keep constant (你原本用 IncrementValue)
        if (s < rise_steps)
        {
            model.applyBC(BC::Dirichlet::IncrementValue(dy, _y), "moving-block-right");
        }
        else
        {
            // no-op (keep)
            // 若你想要 lock 在 shearDisp，可改成 FixedValue(shearDisp, _y) 但要注意每步重複 applyBC 的行為
        }

        // (1) contact kinematics/state update
        ntn_contact.updateInternalData();

        // (2) normal contact
        ntn_contact.computeContactPressure();
        ntn_contact.assembleGlobalContactPressure();
        const auto &f_cn = ntn_contact.getGlobalContactPressure(); // nodal vector-like

        // (3) friction
        // Many implementations assemble internally; still call this every step
        ntn_friction.assembleGlobalFrictionTraction();

        // try to fetch global friction traction if this API exists in your build
        bool has_friction_array = false;

        // (4) inject contact/friction into solid external force (coupler duty)
        detail::addNodal(f_ext, f_cn, dim);

        if constexpr (detail::has_getGlobalFrictionTraction<Friction>::value)
        {
            const auto &f_ct = ntn_friction.getGlobalFrictionTraction();
            detail::addNodal(f_ext, f_ct, dim);
            has_friction_array = true;
        }
        else if constexpr (detail::has_getGlobalFrictionTractions<Friction>::value)
        {
            const auto &f_ct = ntn_friction.getGlobalFrictionTractions();
            detail::addNodal(f_ext, f_ct, dim);
            has_friction_array = true;
        }

        // (5) time integration step
        model.solveStep();

        // (6) remove injected contributions (so f_ext does not accumulate)
        if (has_friction_array)
        {
            if constexpr (detail::has_getGlobalFrictionTraction<Friction>::value)
            {
                const auto &f_ct = ntn_friction.getGlobalFrictionTraction();
                detail::subNodal(f_ext, f_ct, dim);
            }
            else if constexpr (detail::has_getGlobalFrictionTractions<Friction>::value)
            {
                const auto &f_ct = ntn_friction.getGlobalFrictionTractions();
                detail::subNodal(f_ext, f_ct, dim);
            }
        }
        detail::subNodal(f_ext, f_cn, dim);

        // dump
        if (s % DUMP_INTERVAL == 0)
        {
            model.dump();
        }

        // progress / debug
        auto now = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = now - start_time;

        const double denom = (s > 0) ? static_cast<double>(s) : 1.0;
        const double time_per_iter = elapsed.count() / denom;
        const double remaining = time_per_iter * MAX_STEPS - elapsed.count();

        if (prank == 0 && s % 200 == 0)
        {
            std::cout << "Step " << std::setw(8) << s << "/" << MAX_STEPS
                      << " | Elapsed: " << std::fixed << std::setprecision(1)
                      << std::setw(8) << elapsed.count() << " s"
                      << " | ETA: " << std::fixed << std::setprecision(1)
                      << std::setw(8) << std::max(0.0, remaining) << " s\n";
        }

        if (prank == 0 && s % 200 == 0)
        {
            auto &normals_arr = ntn_contact.getNormals().getArray();
            Real nx_min = 1e100, nx_max = -1e100, nx_avg = 0.0;
            for (Idx i = 0; i < normals_arr.size(); ++i)
            {
                const auto nx = normals_arr(i, 0);
                nx_avg += nx;
                nx_min = std::min(nx_min, nx);
                nx_max = std::max(nx_max, nx);
            }
            nx_avg /= std::max<Idx>(1, normals_arr.size());
            std::cout << "[DBG] nx[min,avg,max]=[" << nx_min << "," << nx_avg << "," << nx_max << "]\n";
        }
    }

    if (prank == 0)
    {
        std::cout << "\n";
    }

    akantu::finalize();
    return 0;
}