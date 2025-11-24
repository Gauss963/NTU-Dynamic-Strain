import gmsh, sys

PMMA_THK = [50, 100, 500]
LC       = 20          # 全域網格尺寸

gmsh.initialize()
for thk in PMMA_THK:
    gmsh.model.reset()
    gmsh.model.add(f"CZM_{thk}mm")

    # ---- 幾何 ----
    blk1 = gmsh.model.occ.addBox(  0, 0, 0, 200, 500, thk)
    blk2 = gmsh.model.occ.addBox(200, 0, 0, 145, 550, thk)
    gmsh.model.occ.fragment([(3, blk1), (3, blk2)], [])
    gmsh.model.occ.removeAllDuplicates()
    gmsh.model.occ.synchronize()

    # ---- Physical groups: volume ----
    vols = gmsh.model.getEntities(dim=3)
    v1, v2 = vols[0][1], vols[1][1]
    gmsh.model.addPhysicalGroup(3, [v1], 11)
    gmsh.model.setPhysicalName(3, 11, "moving-block")
    gmsh.model.addPhysicalGroup(3, [v2], 21)
    gmsh.model.setPhysicalName(3, 21, "stationary-block")

    # ---- Physical group: interface surface ----
    s_int = list(set(gmsh.model.getAdjacencies(3, v1)[1]) &
                 set(gmsh.model.getAdjacencies(3, v2)[1]))
    if not s_int:
        gmsh.fatal("❌ 介面面交集為空，請檢查幾何")
    gmsh.model.addPhysicalGroup(2, s_int, 55)
    gmsh.model.setPhysicalName(2, 55, "czm_interface")

    # ---- 網格 ----
    gmsh.model.mesh.setSize(gmsh.model.getEntities(0), LC)
    gmsh.model.mesh.generate(3)

    gmsh.write(f"../Models/{thk}mm-PMMA-CZM.msh")
    print(f"✅ 產生完成：{thk} mm")

gmsh.finalize()