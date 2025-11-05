# Functions.py
import gmsh

def create_block(origin, dimensions, mesh_size, block_name):
    x, y, z = origin
    dx, dy, dz = dimensions

    box = gmsh.model.occ.addBox(x, y, z, dx, dy, dz)
    gmsh.model.occ.synchronize()
    gmsh.model.mesh.setSize(gmsh.model.getEntities(0), mesh_size)

    # 這裡只回傳 Volume tag，不要在這裡掛 Physical Group（避免 fragment 後失效）
    # 網格尺寸全域設一次即可（在主程式）
    return {"volume": box, "name": block_name}