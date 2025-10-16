import gmsh

def create_block(origin, dimensions, mesh_size, tag_prefix=1):
    x, y, z = origin
    dx, dy, dz = dimensions
    box = gmsh.model.occ.addBox(x, y, z, dx, dy, dz)

    gmsh.model.occ.synchronize()

    gmsh.model.mesh.setTransfiniteVolume(box)
    gmsh.model.mesh.setRecombine(3, box)

    faces = gmsh.model.getBoundary([(3, box)], oriented=False)
    for dim, tag in faces:
        gmsh.model.mesh.setTransfiniteSurface(tag)
        gmsh.model.mesh.setRecombine(2, tag)

    gmsh.model.addPhysicalGroup(3, [box], tag=tag_prefix * 10 + 1)
    gmsh.model.setPhysicalName(3, tag_prefix * 10 + 1, "PMMA")

    for dim, tag in faces:
        com = gmsh.model.occ.getCenterOfMass(dim, tag)            
        if abs(com[0] - x) < 1e-6:
            gmsh.model.addPhysicalGroup(2, [tag], tag=tag_prefix * 10 + 4)
            gmsh.model.setPhysicalName(2, tag_prefix * 10 + 4, "front")

        elif abs(com[0] - (x + dx)) < 1e-6:
            gmsh.model.addPhysicalGroup(2, [tag], tag=tag_prefix * 10 + 5)
            gmsh.model.setPhysicalName(2, tag_prefix * 10 + 5, "back")
            
        elif abs(com[1] - y) < 1e-6:
            gmsh.model.addPhysicalGroup(2, [tag], tag=tag_prefix * 10 + 6)
            gmsh.model.setPhysicalName(2, tag_prefix * 10 + 6, "left")
        elif abs(com[1] - (y + dy)) < 1e-6:
            gmsh.model.addPhysicalGroup(2, [tag], tag=tag_prefix * 10 + 7)
            gmsh.model.setPhysicalName(2, tag_prefix * 10 + 7, "right")
            
        elif abs(com[2] - z) < 1e-6:
            gmsh.model.addPhysicalGroup(2, [tag], tag=tag_prefix * 10 + 2)
            gmsh.model.setPhysicalName(2, tag_prefix * 10 + 2, "bottom")
        elif abs(com[2] - (z + dz)) < 1e-6:
            gmsh.model.addPhysicalGroup(2, [tag], tag=tag_prefix * 10 + 3)
            gmsh.model.setPhysicalName(2, tag_prefix * 10 + 3, "top")

    return box