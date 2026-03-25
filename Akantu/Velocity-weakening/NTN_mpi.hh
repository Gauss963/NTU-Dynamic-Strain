std::shared_ptr<akantu::MeshPartition> partitioner;
if (prank == 0) {
    mesh.read(mesh_file);

    const auto psize = comm.getNbProc();

    akantu::Int partition_dir = akantu::_x;
    if constexpr (DIMENSION == 2) {
        partition_dir = (normal_dir == akantu::_x) ? akantu::_y : akantu::_x;
    } else {
        // For your current setup with normal_dir == _x, keep the split along y.
        partition_dir = (normal_dir == akantu::_x) ? akantu::_y : akantu::_x;
    }

    const auto lower = mesh.getLowerBounds()(partition_dir);
    const auto upper = mesh.getUpperBounds()(partition_dir);
    const auto slice =
        (upper - lower) / static_cast<akantu::Real>(psize);

    auto partition =
        std::make_shared<akantu::ElementTypeMapArray<akantu::Idx>>("partition");
    partition->initialize(mesh, akantu::_spatial_dimension = DIMENSION,
                          akantu::_with_nb_element = true);

    std::set<akantu::Idx> sanity_proc_check;

    akantu::for_each_element(
        mesh,
        [&](auto && element) {
            auto barycenter = mesh.getBarycenter(element);

            auto raw = static_cast<akantu::Idx>(
                std::floor((barycenter(partition_dir) - lower) / slice));
            auto proc = std::max<akantu::Idx>(
                0, std::min<akantu::Idx>(psize - 1, raw));

            sanity_proc_check.insert(proc);
            (*partition)(element) = proc;
        },
        akantu::_spatial_dimension = DIMENSION,
        akantu::_ghost_type = akantu::_not_ghost);

    AKANTU_DEBUG_ASSERT(sanity_proc_check.size() == std::size_t(psize),
                        "There are too many MPI ranks for this mesh partition");

    partitioner =
        std::make_shared<akantu::MeshPartitionMeshData>(mesh, partition, DIMENSION);
    partitioner->partitionate(psize);
}
mesh.distribute(partitioner);
