import gmsh

def create_block(origin, dimensions, mesh_size, block_name, tag_prefix=1):
    x, y, z = origin
    dx, dy, dz = dimensions
    box = gmsh.model.occ.addBox(x, y, z, dx, dy, dz)

    gmsh.model.occ.synchronize()

    gmsh.model.mesh.setTransfiniteVolume(box)
    gmsh.model.mesh.setRecombine(3, box)

    faces = gmsh.model.getBoundary([(3, box)], oriented=False)
    face_tags = {}
    face_phys = {}

    for dim, tag in faces:
        pass
        gmsh.model.mesh.setTransfiniteSurface(tag)
        gmsh.model.mesh.setRecombine(2, tag)


    gmsh.model.addPhysicalGroup(3, [box], tag=tag_prefix * 10 + 1)
    gmsh.model.setPhysicalName(3, tag_prefix * 10 + 1, block_name)

    gmsh.model.occ.synchronize()
    gmsh.model.mesh.setSize(gmsh.model.getEntities(0), mesh_size)
    gmsh.option.setNumber("Mesh.MeshSizeMin", mesh_size)
    gmsh.option.setNumber("Mesh.MeshSizeMax", mesh_size)

    tolerance = 1e-2

    for dim, tag in faces:
        com = gmsh.model.occ.getCenterOfMass(dim, tag)

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

    gmsh.model.occ.synchronize()
    return {"volume": box, "faces_geo": face_tags, "faces_phys": face_phys}

def create_block_2d_quad(origin, dimensions, mesh_size, block_name, tag_prefix=1):
    x, y, z = origin
    dx, dy, _ = dimensions

    # Rectangle
    surf = gmsh.model.occ.addRectangle(x, y, z, dx, dy)
    gmsh.model.occ.synchronize()

    # --- Get boundary curves ---
    curves = gmsh.model.getBoundary([(2, surf)], oriented=False)

    # Number of nodes along each direction
    nx = int(round(dx / mesh_size)) + 1
    ny = int(round(dy / mesh_size)) + 1

    # Tolerance for edge identification
    tol = 1e-6

    edge_tags = {}
    edge_phys = {}

    for dim, tag in curves:
        com = gmsh.model.occ.getCenterOfMass(dim, tag)

        # Vertical edges → ny nodes
        if abs(com[0] - x) < tol or abs(com[0] - (x + dx)) < tol:
            gmsh.model.mesh.setTransfiniteCurve(tag, ny)

            if abs(com[0] - x) < tol:
                name = f"{block_name}-back"
                tag_val = tag_prefix * 10 + 4
            else:
                name = f"{block_name}-front"
                tag_val = tag_prefix * 10 + 5

        # Horizontal edges → nx nodes
        elif abs(com[1] - y) < tol or abs(com[1] - (y + dy)) < tol:
            gmsh.model.mesh.setTransfiniteCurve(tag, nx)

            if abs(com[1] - y) < tol:
                name = f"{block_name}-right"
                tag_val = tag_prefix * 10 + 6
            else:
                name = f"{block_name}-left"
                tag_val = tag_prefix * 10 + 7
        else:
            continue

        gmsh.model.addPhysicalGroup(1, [tag], tag=tag_val)
        gmsh.model.setPhysicalName(1, tag_val, name)

        edge_tags[name] = tag
        edge_phys[name] = tag_val

    # --- Transfinite surface & recombine ---
    gmsh.model.mesh.setTransfiniteSurface(surf)
    gmsh.model.mesh.setRecombine(2, surf)

    # --- Physical surface ---
    gmsh.model.addPhysicalGroup(2, [surf], tag=tag_prefix * 10 + 1)
    gmsh.model.setPhysicalName(2, tag_prefix * 10 + 1, block_name)

    gmsh.model.occ.synchronize()

    return {
        "surface": surf,
        "edges_geo": edge_tags,
        "edges_phys": edge_phys,
        "nx": nx,
        "ny": ny
    }