import gmsh

def create_block(origin, dimensions, mesh_size, tag_prefix=1):
    """
    建立一個六面體的 block 並使用 Transfinite 結構網格（保證節點對齊），適用於接觸分析。

    Parameters:
    - origin: (x, y, z) 起點座標
    - dimensions: (Lx, Ly, Lz) 長寬高
    - mesh_size: 用於節點密度控制的基準長度
    - tag_prefix: 用來區分各個 block 的 base ID
    """

    x, y, z = origin
    Lx, Ly, Lz = dimensions

    # 計算各方向節點數量（確保 >= 2）
    nx = max(2, round(Lx / mesh_size))
    ny = max(2, round(Ly / mesh_size))
    nz = max(2, round(Lz / mesh_size))

    # Tag generator
    pt_tag = lambda i: tag_prefix * 100 + i
    line_tag = lambda i: tag_prefix * 1000 + i
    loop_tag = lambda i: tag_prefix * 2000 + i
    surface_tag = lambda i: tag_prefix * 3000 + i
    volume_tag = tag_prefix * 4000

    # === 1. Points ===
    # Bottom face
    gmsh.model.geo.addPoint(x     , y     , z     , mesh_size, pt_tag(1))
    gmsh.model.geo.addPoint(x+Lx  , y     , z     , mesh_size, pt_tag(2))
    gmsh.model.geo.addPoint(x+Lx  , y+Ly , z     , mesh_size, pt_tag(3))
    gmsh.model.geo.addPoint(x     , y+Ly , z     , mesh_size, pt_tag(4))

    # Top face
    gmsh.model.geo.addPoint(x     , y     , z+Lz  , mesh_size, pt_tag(5))
    gmsh.model.geo.addPoint(x+Lx  , y     , z+Lz  , mesh_size, pt_tag(6))
    gmsh.model.geo.addPoint(x+Lx  , y+Ly , z+Lz  , mesh_size, pt_tag(7))
    gmsh.model.geo.addPoint(x     , y+Ly , z+Lz  , mesh_size, pt_tag(8))

    # === 2. Lines ===
    # Bottom
    gmsh.model.geo.addLine(pt_tag(1), pt_tag(2), line_tag(1))
    gmsh.model.geo.addLine(pt_tag(2), pt_tag(3), line_tag(2))
    gmsh.model.geo.addLine(pt_tag(3), pt_tag(4), line_tag(3))
    gmsh.model.geo.addLine(pt_tag(4), pt_tag(1), line_tag(4))

    # Top
    gmsh.model.geo.addLine(pt_tag(5), pt_tag(6), line_tag(5))
    gmsh.model.geo.addLine(pt_tag(6), pt_tag(7), line_tag(6))
    gmsh.model.geo.addLine(pt_tag(7), pt_tag(8), line_tag(7))
    gmsh.model.geo.addLine(pt_tag(8), pt_tag(5), line_tag(8))

    # Vertical
    gmsh.model.geo.addLine(pt_tag(1), pt_tag(5), line_tag(9))
    gmsh.model.geo.addLine(pt_tag(2), pt_tag(6), line_tag(10))
    gmsh.model.geo.addLine(pt_tag(3), pt_tag(7), line_tag(11))
    gmsh.model.geo.addLine(pt_tag(4), pt_tag(8), line_tag(12))

    # === 3. Transfinite Lines ===
    # x direction: 1-2, 5-6, 4-3, 8-7
    gmsh.model.geo.mesh.setTransfiniteCurve(line_tag(1), nx)
    gmsh.model.geo.mesh.setTransfiniteCurve(line_tag(2), ny)
    gmsh.model.geo.mesh.setTransfiniteCurve(line_tag(3), nx)
    gmsh.model.geo.mesh.setTransfiniteCurve(line_tag(4), ny)

    gmsh.model.geo.mesh.setTransfiniteCurve(line_tag(5), nx)
    gmsh.model.geo.mesh.setTransfiniteCurve(line_tag(6), ny)
    gmsh.model.geo.mesh.setTransfiniteCurve(line_tag(7), nx)
    gmsh.model.geo.mesh.setTransfiniteCurve(line_tag(8), ny)

    gmsh.model.geo.mesh.setTransfiniteCurve(line_tag(9), nz)
    gmsh.model.geo.mesh.setTransfiniteCurve(line_tag(10), nz)
    gmsh.model.geo.mesh.setTransfiniteCurve(line_tag(11), nz)
    gmsh.model.geo.mesh.setTransfiniteCurve(line_tag(12), nz)

    # === 4. Curve Loops + Surfaces ===
    # Bottom
    gmsh.model.geo.addCurveLoop([line_tag(1), line_tag(2), line_tag(3), line_tag(4)], loop_tag(1))
    gmsh.model.geo.addPlaneSurface([loop_tag(1)], surface_tag(1))

    # Top
    gmsh.model.geo.addCurveLoop([line_tag(5), line_tag(6), line_tag(7), line_tag(8)], loop_tag(2))
    gmsh.model.geo.addPlaneSurface([loop_tag(2)], surface_tag(2))

    # Front (x-z)
    gmsh.model.geo.addCurveLoop([line_tag(1), line_tag(10), -line_tag(5), -line_tag(9)], loop_tag(3))
    gmsh.model.geo.addPlaneSurface([loop_tag(3)], surface_tag(3))

    # Right (y-z)
    gmsh.model.geo.addCurveLoop([line_tag(2), line_tag(11), -line_tag(6), -line_tag(10)], loop_tag(4))
    gmsh.model.geo.addPlaneSurface([loop_tag(4)], surface_tag(4))

    # Back (x-z)
    gmsh.model.geo.addCurveLoop([line_tag(3), line_tag(12), -line_tag(7), -line_tag(11)], loop_tag(5))
    gmsh.model.geo.addPlaneSurface([loop_tag(5)], surface_tag(5))

    # Left (y-z)
    gmsh.model.geo.addCurveLoop([line_tag(4), line_tag(9), -line_tag(8), -line_tag(12)], loop_tag(6))
    gmsh.model.geo.addPlaneSurface([loop_tag(6)], surface_tag(6))

    # === 5. Transfinite Surfaces + Recombine ===
    for i in range(1, 7):
        gmsh.model.geo.mesh.setTransfiniteSurface(surface_tag(i))
        gmsh.model.geo.mesh.setRecombine(2, surface_tag(i))

    # === 6. Volume ===
    gmsh.model.geo.addSurfaceLoop([surface_tag(i) for i in range(1, 7)], loop_tag(7))
    gmsh.model.geo.addVolume([loop_tag(7)], volume_tag)
    gmsh.model.geo.mesh.setTransfiniteVolume(volume_tag)