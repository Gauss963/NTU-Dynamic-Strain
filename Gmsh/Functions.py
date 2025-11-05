import gmsh

def create_block(origin, dimensions, mesh_size, block_name, tag_prefix=1):
    x, y, z = origin
    dx, dy, dz = dimensions
    box = gmsh.model.occ.addBox(x, y, z, dx, dy, dz)

    gmsh.model.occ.synchronize()

    gmsh.model.mesh.setTransfiniteVolume(box)
    # gmsh.model.mesh.setRecombine(3, box)

    faces = gmsh.model.getBoundary([(3, box)], oriented=False)
    face_tags = {}
    face_phys = {}

    for dim, tag in faces:
        gmsh.model.mesh.setTransfiniteSurface(tag)
        # gmsh.model.mesh.setRecombine(2, tag)

    gmsh.model.addPhysicalGroup(3, [box], tag=tag_prefix * 10 + 1)
    gmsh.model.setPhysicalName(3, tag_prefix * 10 + 1, block_name)

    gmsh.model.occ.synchronize()
    gmsh.model.mesh.setSize(gmsh.model.getEntities(0), mesh_size)

    tolerance = 1e-2

    for dim, tag in faces:
        com = gmsh.model.occ.getCenterOfMass(dim, tag)

        if abs(com[0] - x) < tolerance:
            name = f"{block_name}-front"
            tag_val = tag_prefix * 10 + 4
        elif abs(com[0] - (x + dx)) < tolerance:
            name = f"{block_name}-back"
            tag_val = tag_prefix * 10 + 5
        elif abs(com[1] - y) < tolerance:
            name = f"{block_name}-left"
            tag_val = tag_prefix * 10 + 6
        elif abs(com[1] - (y + dy)) < tolerance:
            name = f"{block_name}-right"
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

    gmsh.model.occ.synchronize()
    return {"volume": box, "faces_geo": face_tags, "faces_phys": face_phys}