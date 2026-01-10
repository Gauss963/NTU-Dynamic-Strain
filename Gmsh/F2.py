import gmsh
import math

def create_block(origin, dimensions, mesh_size, block_name, tag_prefix=1):
    x, y, z = origin
    dx, dy, dz = dimensions
    
    # 建立幾何
    box = gmsh.model.occ.addBox(x, y, z, dx, dy, dz)
    gmsh.model.occ.synchronize()

    # --- [關鍵修正] 設定 Transfinite (結構化網格) ---
    # 這確保節點位置是數學上可預測的，從而讓兩個方塊的接觸面節點完美重疊
    
    # 計算每個方向需要的節點數 (Ceiling 確保至少有一格)
    nx = int(math.ceil(dx / mesh_size)) + 1
    ny = int(math.ceil(dy / mesh_size)) + 1
    nz = int(math.ceil(dz / mesh_size)) + 1

    # 取得 Box 的所有 Lines (Curves)
    # Box 的線條順序不一定，所以我們透過 Bounding Box 來判斷線條方向
    lines = gmsh.model.getBoundary([(3, box)], combined=False, oriented=False, recursive=True)
    
    for _, line_tag in lines:
        # 取得線條的邊界框來判斷它是長、寬還是高
        bbox = gmsh.model.getBoundingBox(1, line_tag)
        lx = abs(bbox[3] - bbox[0])
        ly = abs(bbox[4] - bbox[1])
        lz = abs(bbox[5] - bbox[2])

        # 設定線條上的節點數
        if lx > ly and lx > lz: # X方向
            gmsh.model.mesh.setTransfiniteCurve(line_tag, nx)
        elif ly > lx and ly > lz: # Y方向
            gmsh.model.mesh.setTransfiniteCurve(line_tag, ny)
        else: # Z方向 (或厚度方向)
            gmsh.model.mesh.setTransfiniteCurve(line_tag, nz)

    # 取得 Box 的所有 Surfaces
    surfaces = gmsh.model.getBoundary([(3, box)], combined=True, oriented=False)
    for _, surf_tag in surfaces:
        # 設定表面為 Transfinite，並指定排列方式 (Right)
        gmsh.model.mesh.setTransfiniteSurface(surf_tag)
        # 為了獲得更好的四邊形網格，通常會 recombine (變成 Hex/Quad)
        # 但 Akantu 的 Tetrahedron 也可以運作，只要節點對齊即可。
        # 如果你想要純六面體 (Hex) 網格，請取消下面這行的註解：
        # gmsh.model.mesh.setRecombine(2, surf_tag)

    # 設定體積為 Transfinite (這會生成非常整齊的網格)
    gmsh.model.mesh.setTransfiniteVolume(box)
    # 如果上面 recombine 了表面，這裡也要 recombine 體積
    # gmsh.model.mesh.setRecombine(3, box)

    # ----------------------------------------------------

    # 設定 Physical Groups (這部分維持你的邏輯)
    faces = gmsh.model.getBoundary([(3, box)], oriented=False)
    face_tags = {}
    face_phys = {}

    gmsh.model.addPhysicalGroup(3, [box], tag=tag_prefix * 10 + 1)
    gmsh.model.setPhysicalName(3, tag_prefix * 10 + 1, block_name)

    tolerance = 1e-2
    for dim, tag in faces:
        com = gmsh.model.occ.getCenterOfMass(dim, tag)
        name = ""
        tag_val = 0

        if abs(com[0] - x) < tolerance:
            name = f"{block_name}-back"
            tag_val = tag_prefix * 10 + 4
        elif abs(com[0] - (x + dx)) < tolerance:
            name = f"{block_name}-front"
            tag_val = tag_prefix * 10 + 5
        elif abs(com[1] - y) < tolerance:
            name = f"{block_name}-right"
            tag_val = tag_prefix * 10 + 6
        elif abs(com[1] - (y + dy)) < tolerance:
            name = f"{block_name}-left"
            tag_val = tag_prefix * 10 + 7
        elif abs(com[2] - z) < tolerance:
            name = f"{block_name}-bottom"
            tag_val = tag_prefix * 10 + 2
        elif abs(com[2] - (z + dz)) < tolerance:
            name = f"{block_name}-top"
            tag_val = tag_prefix * 10 + 3
        else:
            continue

        gmsh.model.addPhysicalGroup(2, [tag], tag=tag_val)
        gmsh.model.setPhysicalName(2, tag_val, name)
        face_tags[name] = tag
        face_phys[name] = tag_val

    return {"volume": box, "faces_geo": face_tags, "faces_phys": face_phys}