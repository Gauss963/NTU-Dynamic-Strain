import gmsh
import SFunctions as Functions

def main():
    PMMA_THICKNESSES = [50, 100]
    mesh_size = 5

    for t in PMMA_THICKNESSES:
        gmsh.initialize()
        gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
        gmsh.option.setNumber("Mesh.ElementOrder", 1)
        gmsh.model.add("ContactModel")

        Lx1, Ly1, Lz = 200.0, 500.0, float(t)
        Lx2, Ly2      = 145.0, 550.0
        x_split       = 200.0

        blk1 = Functions.create_block((0, 0, 0),       (Lx1, Ly1, Lz), mesh_size, "moving-block")
        blk2 = Functions.create_block((x_split, 0, 0), (Lx2, Ly2, Lz), mesh_size, "stationary-block")

        gmsh.model.mesh.setSize(gmsh.model.getEntities(0), mesh_size)

        gmsh.model.occ.fragment([(3, blk1["volume"])], [(3, blk2["volume"])])
        gmsh.model.occ.synchronize()

        vols = [tag for (dim, tag) in gmsh.model.getEntities(3)]
        left_vols, right_vols = [], []
        for v in vols:
            xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.occ.getBoundingBox(3, v)
            if abs(xmin - 0.0) < 1e-3 and abs(xmax - x_split) < 1e-3:
                left_vols.append(v)
            elif abs(xmin - x_split) < 1e-3 and abs(xmax - (x_split + Lx2)) < 1e-3:
                right_vols.append(v)
        assert left_vols and right_vols, "找不到 fragment 後的左右體積"

        # fragment 後再掛 Physical Volume（交給 Gmsh 分配 tag）
        pgL = gmsh.model.addPhysicalGroup(3, left_vols)
        gmsh.model.setPhysicalName(3, pgL, "moving-block")
        pgR = gmsh.model.addPhysicalGroup(3, right_vols)
        gmsh.model.setPhysicalName(3, pgR, "stationary-block")

        # 共享面 = 左右邊界交集
        left_faces, right_faces = set(), set()
        for v in left_vols:
            left_faces |= {t for (d, t) in gmsh.model.getBoundary([(3, v)], oriented=False) if d == 2}
        for v in right_vols:
            right_faces |= {t for (d, t) in gmsh.model.getBoundary([(3, v)], oriented=False) if d == 2}
        iface = sorted(left_faces & right_faces)
        assert iface, "沒有任何共享面（介面不共形）"

        # 單一 interface 名稱（extrinsic 最穩）
        pgI = gmsh.model.addPhysicalGroup(2, iface)
        gmsh.model.setPhysicalName(2, pgI, "interface_zone")

        gmsh.model.mesh.generate(3)
        gmsh.write(f"../Models/{t}mm-PMMA.msh")
        gmsh.finalize()

if __name__ == "__main__":
    main()