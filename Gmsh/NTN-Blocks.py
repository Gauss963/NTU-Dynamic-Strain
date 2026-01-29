import gmsh
import math


def _classify_box_faces(box_tag, origin, dims, tol=1e-8):
    """
    Return dict of face tags for a box volume:
      back/front  : x = x0 / x0+dx
      right/left  : y = y0 / y0+dy   (注意：沿用你原本 Functions.py 的命名：y=y0 -> right, y=y0+dy -> left)
      bottom/top  : z = z0 / z0+dz
    """
    x0, y0, z0 = origin
    dx, dy, dz = dims

    faces = gmsh.model.getBoundary([(3, box_tag)], oriented=False)
    out = {}

    for dim, s in faces:
        com = gmsh.model.occ.getCenterOfMass(dim, s)
        x, y, z = com

        if abs(x - x0) < tol:
            out["back"] = s
        elif abs(x - (x0 + dx)) < tol:
            out["front"] = s
        elif abs(y - y0) < tol:
            out["right"] = s     # 沿用你舊 mapping
        elif abs(y - (y0 + dy)) < tol:
            out["left"] = s      # 沿用你舊 mapping
        elif abs(z - z0) < tol:
            out["bottom"] = s
        elif abs(z - (z0 + dz)) < tol:
            out["top"] = s

    return out


def _set_transfinite_all(mesh_size):
    gmsh.model.occ.synchronize()

    curves = gmsh.model.getEntities(1)
    for _, ctag in curves:
        bnd = gmsh.model.getBoundary([(1, ctag)], oriented=False)
        if len(bnd) != 2:
            continue

        p0 = bnd[0][1]
        p1 = bnd[1][1]

        x0, y0, z0 = gmsh.model.getValue(0, p0, [])
        x1, y1, z1 = gmsh.model.getValue(0, p1, [])

        L = ((x1 - x0)**2 + (y1 - y0)**2 + (z1 - z0)**2) ** 0.5
        n = max(2, int(round(L / mesh_size)) + 1)

        gmsh.model.mesh.setTransfiniteCurve(ctag, n)

    for _, stag in gmsh.model.getEntities(2):
        gmsh.model.mesh.setTransfiniteSurface(stag)
        gmsh.model.mesh.setRecombine(2, stag)

    for _, vtag in gmsh.model.getEntities(3):
        gmsh.model.mesh.setTransfiniteVolume(vtag)
        gmsh.model.mesh.setRecombine(3, vtag)


def build_mesh(PMMA_thickness_mm, mesh_size, scaling_factor=0.1, initial_offset_mm=0.0):
    """
    Generate split stationary block (500 + 50) with contact patch 500 high.
    Writes:
      ../Models/{thickness}mm-BS-PMMA.msh
      ../Models/{thickness}mm-BS-PMMA.brep
    """
    sf = scaling_factor
    th = PMMA_thickness_mm * sf

    # Dimensions (scaled)
    mov_origin = (0.0, 0.0, 0.0)
    mov_dims   = (200.0 * sf, 500.0 * sf, th)

    # stationary: split into lower (500) + upper (50)
    off = initial_offset_mm * sf
    sta_origin_low = (200.0 * sf + off, 0.0, 0.0)
    sta_dims_low   = (145.0 * sf, 500.0 * sf, th)

    sta_origin_up  = (200.0 * sf + off, 500.0 * sf, 0.0)
    sta_dims_up    = (145.0 * sf, 50.0 * sf, th)

    gmsh.initialize()
    gmsh.option.setNumber("General.NumThreads", 10)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.ElementOrder", 1)
    gmsh.option.setNumber("Mesh.Algorithm3D", 1)

    # IMPORTANT: 不要用 Netgen optimize（可能把 structured mesh 搞成 tet）
    gmsh.option.setNumber("Mesh.Optimize", 0)
    gmsh.option.setNumber("Mesh.OptimizeNetgen", 0)

    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.Binary", 0)
    gmsh.model.add("ContactModel")

    # geometry
    mov = gmsh.model.occ.addBox(*mov_origin, *mov_dims)
    sta_low = gmsh.model.occ.addBox(*sta_origin_low, *sta_dims_low)
    sta_up  = gmsh.model.occ.addBox(*sta_origin_up,  *sta_dims_up)

    gmsh.model.occ.synchronize()

    # global mesh size
    gmsh.model.mesh.setSize(gmsh.model.getEntities(0), mesh_size)
    gmsh.option.setNumber("Mesh.MeshSizeMin", mesh_size)
    gmsh.option.setNumber("Mesh.MeshSizeMax", mesh_size)

    # --- classify faces ---
    tol = 1e-7
    mov_faces = _classify_box_faces(mov, mov_origin, mov_dims, tol=tol)
    sta_low_faces = _classify_box_faces(sta_low, sta_origin_low, sta_dims_low, tol=tol)
    sta_up_faces  = _classify_box_faces(sta_up,  sta_origin_up,  sta_dims_up,  tol=tol)

    # --- physical groups: volumes ---
    # keep your exact naming/tag pattern
    # moving-block volume tag = 11
    gmsh.model.addPhysicalGroup(3, [mov], tag=11)
    gmsh.model.setPhysicalName(3, 11, "moving-block")

    # stationary-block volume tag = 21 (contains both low+up volumes)
    gmsh.model.addPhysicalGroup(3, [sta_low, sta_up], tag=21)
    gmsh.model.setPhysicalName(3, 21, "stationary-block")

    # --- physical groups: surfaces ---
    # moving surfaces: bottom12 top13 back14 front15 right16 left17
    gmsh.model.addPhysicalGroup(2, [mov_faces["bottom"]], tag=12); gmsh.model.setPhysicalName(2, 12, "moving-block-bottom")
    gmsh.model.addPhysicalGroup(2, [mov_faces["top"]],    tag=13); gmsh.model.setPhysicalName(2, 13, "moving-block-top")
    gmsh.model.addPhysicalGroup(2, [mov_faces["back"]],   tag=14); gmsh.model.setPhysicalName(2, 14, "moving-block-back")
    gmsh.model.addPhysicalGroup(2, [mov_faces["front"]],  tag=15); gmsh.model.setPhysicalName(2, 15, "moving-block-front")  # contact master
    gmsh.model.addPhysicalGroup(2, [mov_faces["right"]],  tag=16); gmsh.model.setPhysicalName(2, 16, "moving-block-right")
    gmsh.model.addPhysicalGroup(2, [mov_faces["left"]],   tag=17); gmsh.model.setPhysicalName(2, 17, "moving-block-left")

    # stationary surfaces: bottom22 top23 back24 front25 right26 left27
    # IMPORTANT:
    #   - stationary-block-back 只放 "下半塊" 的 back 面 (500 高) => NTN 才能 node-to-node 配對
    #   - 其他面可合併 low+up 兩塊的外表面
    gmsh.model.addPhysicalGroup(2, [sta_low_faces["bottom"]], tag=22); gmsh.model.setPhysicalName(2, 22, "stationary-block-bottom")

    # top: both low+up have a top face (z=max), include both
    gmsh.model.addPhysicalGroup(2, [sta_low_faces["top"], sta_up_faces["top"]], tag=23); gmsh.model.setPhysicalName(2, 23, "stationary-block-top")

    # back: ONLY lower back is contact slave
    gmsh.model.addPhysicalGroup(2, [sta_low_faces["back"]], tag=24); gmsh.model.setPhysicalName(2, 24, "stationary-block-back")

    # front/left/right: union of low+up external faces
    gmsh.model.addPhysicalGroup(2, [sta_low_faces["front"], sta_up_faces["front"]], tag=25); gmsh.model.setPhysicalName(2, 25, "stationary-block-front")
    gmsh.model.addPhysicalGroup(2, [sta_low_faces["right"], sta_up_faces["right"]], tag=26); gmsh.model.setPhysicalName(2, 26, "stationary-block-right")
    gmsh.model.addPhysicalGroup(2, [sta_low_faces["left"],  sta_up_faces["left"]],  tag=27); gmsh.model.setPhysicalName(2, 27, "stationary-block-left")

    # --- structured mesh (node-to-node on contact patch) ---
    _set_transfinite_all(mesh_size)

    gmsh.model.occ.synchronize()
    gmsh.model.mesh.generate(3)

    # sanity print
    types3d, _, _ = gmsh.model.mesh.getElements(3)
    print(f"[mesh] thickness={PMMA_thickness_mm}mm  3D element types = {types3d}")

    gmsh.write(f"../Models/{PMMA_thickness_mm}mm-BS-PMMA.msh")
    gmsh.write(f"../Models/{PMMA_thickness_mm}mm-BS-PMMA.brep")

    # gmsh.fltk.run()
    gmsh.finalize()


def main():
    PMMA_THICKNESSES = [50, 100, 500]
    mesh_size = 2.0
    # scaling_factor = 1 / 10
    scaling_factor = 1

    initial_offset_mm = 0.0

    for th in PMMA_THICKNESSES:
        build_mesh(th, mesh_size, scaling_factor=scaling_factor, initial_offset_mm=initial_offset_mm)


if __name__ == "__main__":
    main()