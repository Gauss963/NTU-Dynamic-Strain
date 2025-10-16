import gmsh

def create_block(origin, dimensions, mesh_size, tag_prefix=1):
    x, y, z = origin
    dx, dy, dz = dimensions
    box = gmsh.model.occ.addBox(x, y, z, dx, dy, dz)

    gmsh.model.occ.synchronize()

    gmsh.model.mesh.setTransfiniteVolume(box)
    gmsh.model.mesh.setRecombine(3, box)

    faces = gmsh.model.getBoundary([(3, box)], oriented=False)
    face_tags = {}

    for dim, tag in faces:
        gmsh.model.mesh.setTransfiniteSurface(tag)
        gmsh.model.mesh.setRecombine(2, tag)

    gmsh.model.addPhysicalGroup(3, [box], tag=tag_prefix * 10 + 1)
    gmsh.model.setPhysicalName(3, tag_prefix * 10 + 1, "PMMA")

    for dim, tag in faces:
        com = gmsh.model.occ.getCenterOfMass(dim, tag)

        if abs(com[0] - x) < 1e-6:
            name = "front"
            tag_val = tag_prefix * 10 + 4
        elif abs(com[0] - (x + dx)) < 1e-6:
            name = "back"
            tag_val = tag_prefix * 10 + 5
        elif abs(com[1] - y) < 1e-6:
            name = "left"
            tag_val = tag_prefix * 10 + 6
        elif abs(com[1] - (y + dy)) < 1e-6:
            name = "right"
            tag_val = tag_prefix * 10 + 7
        elif abs(com[2] - z) < 1e-6:
            name = "bottom"
            tag_val = tag_prefix * 10 + 2
        elif abs(com[2] - (z + dz)) < 1e-6:
            name = "top"
            tag_val = tag_prefix * 10 + 3
        else:
            continue

        gmsh.model.addPhysicalGroup(2, [tag], tag=tag_val)
        gmsh.model.setPhysicalName(2, tag_val, name)
        face_tags[name] = tag

    return {
        "volume": box,
        "faces": face_tags
    }