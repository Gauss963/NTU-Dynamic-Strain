#include "aka_common.hh"
#include "mesh.hh"

static void rebuild_contact_groups_xplane_by_normal(akantu::Mesh &mesh,
                                                    akantu::Real x_plane,
                                                    akantu::Real tol,
                                                    const std::string &slave_name,
                                                    const std::string &master_name,
                                                    akantu::Int normal_dir = 0) {
    const akantu::Int sd = mesh.getSpatialDimension();
    const akantu::Int face_dim = sd - 1;

    auto &slave = mesh.createElementGroup(slave_name, face_dim);
    auto &master = mesh.createElementGroup(master_name, face_dim);

    const auto &nodes = mesh.getNodes();

    auto get_point = [&](akantu::Idx n) -> akantu::Vector<akantu::Real> {
        akantu::Vector<akantu::Real> p(sd);
        for (akantu::Int d = 0; d < sd; ++d)
            p(d) = nodes(n, d);
        return p;
    };

    for (auto type : mesh.elementTypes(face_dim, akantu::_not_ghost, akantu::_ek_regular)) {
        const auto nb = mesh.getNbElement(type, akantu::_not_ghost);
        for (akantu::Idx e = 0; e < nb; ++e) {
            akantu::Element face{type, e, akantu::_not_ghost};

            const auto bary = mesh.getBarycenter(face);
            if (std::abs(bary(normal_dir) - x_plane) > tol)
                continue;

            // 取前三個節點計算法向（tri/quad 都可用前三點近似）
            const auto conn = mesh.getConnectivity(face); // should be iterable
            const akantu::Idx n0 = conn(0);
            const akantu::Idx n1 = conn(1);
            const akantu::Idx n2 = conn(2);

            auto p0 = get_point(n0);
            auto p1 = get_point(n1);
            auto p2 = get_point(n2);

            // normal = (p1-p0) x (p2-p0)
            akantu::Vector<akantu::Real> v1(sd), v2(sd), n(sd);
            for (akantu::Int d = 0; d < sd; ++d) {
                v1(d) = p1(d) - p0(d);
                v2(d) = p2(d) - p0(d);
            }

            // cross in 3D only; 你是 3D (sd=3)
            n(0) = v1(1) * v2(2) - v1(2) * v2(1);
            n(1) = v1(2) * v2(0) - v1(0) * v2(2);
            n(2) = v1(0) * v2(1) - v1(1) * v2(0);

            // 用 normal 的 x 分量分 side
            if (n(normal_dir) >= 0.)
                master.add(face);
            else
                slave.add(face);
        }
    }
}