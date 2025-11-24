import gmsh
import math

msh_file = "../Models/50mm-PMMA-CZM.msh"
tol = 1e-6

targets = [
    (173.057, 346.667, 41.6667),
    (301.607, 254.944, 33.3333),
]

gmsh.initialize()
gmsh.open(msh_file)

# 先拿到所有 2D element 的節點與 ID
types, tags, nodeTags = gmsh.model.mesh.getElements(dim=2)

print("2D element types:", types)

# 通常只會有一種 type，例如 _triangle_3 -> Gmsh type 2
for etype, etags, enodes in zip(types, tags, nodeTags):
    # 每個三角形 3 個節點
    num_per_elem = int(len(enodes) / len(etags))
    for i, elem_tag in enumerate(etags):
        node_ids = enodes[i * num_per_elem:(i + 1) * num_per_elem]
        # 取出節點座標
        coords = []
        for nid in node_ids:
            x, y, z = gmsh.model.mesh.getNode(nid)[0]
            coords.append((x, y, z))
        # 算 barycenter
        bx = sum(c[0] for c in coords) / 3.0
        by = sum(c[1] for c in coords) / 3.0
        bz = sum(c[2] for c in coords) / 3.0

        for tx, ty, tz in targets:
            if (abs(bx - tx) < tol and
                abs(by - ty) < tol and
                abs(bz - tz) < tol):
                print("\n=== Found suspect facet ===")
                print(f"Element tag: {elem_tag}, type: {etype}")
                print(f"Barycenter: ({bx:.6f}, {by:.6f}, {bz:.6f})")
                print(f"Node IDs: {node_ids}")

                # 查它屬於哪些 entity（surface）、有哪些 physical group
                # 先從 node 去查 entity 不方便，改用 bounding box + getEntitiesInBoundingBox
                bb_tol = 1e-3
                ents = gmsh.model.getEntitiesInBoundingBox(
                    bx - bb_tol, by - bb_tol, bz - bb_tol,
                    bx + bb_tol, by + bb_tol, bz + bb_tol,
                    dim=2
                )
                print("Geometrical entities (dim=2) around it:", ents)

                for dim, tag in ents:
                    phys = gmsh.model.getPhysicalGroupsForEntity(dim, tag)
                    print(f"  Surface (dim={dim}, tag={tag}) physical groups:", phys)

                # 再看它連到哪些 3D elements（tetra）
                conn = gmsh.model.mesh.getElementAdjacent(dim=2, tag=elem_tag, targetDim=3)
                print("  Adjacent 3D elements:", conn)

gmsh.finalize()