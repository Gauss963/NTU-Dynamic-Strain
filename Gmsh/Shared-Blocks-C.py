import gmsh

def main():

    PMMA_THICKNESSES = [50, 100, 500]
    mesh_size = 5
    mesh_size = 1

    for PMMA_thickness in PMMA_THICKNESSES:

        gmsh.initialize()
        gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
        gmsh.option.setNumber("Mesh.Algorithm3D", 1)
        gmsh.option.setNumber("Mesh.ElementOrder", 1)
        gmsh.option.setNumber("Mesh.SubdivisionAlgorithm", 1) 
        gmsh.option.setNumber("Mesh.Binary", 0)
        gmsh.option.setNumber("Mesh.MeshSizeMin", mesh_size)
        gmsh.option.setNumber("Mesh.MeshSizeMax", mesh_size)
        gmsh.option.setNumber("Mesh.Smoothing", 20)

        gmsh.model.add("ContactModel")

        blk1_origin = (0, 0, 0)
        blk1_dims = (200, 500, PMMA_thickness)

        blk2_origin = (200, 0, 0)
        blk2_dims   = (145, 550, PMMA_thickness)

        blk1_tag = gmsh.model.occ.addBox(
            blk1_origin[0], blk1_origin[1], blk1_origin[2],
            blk1_dims[0], blk1_dims[1], blk1_dims[2]
        )

        blk2_tag = gmsh.model.occ.addBox(
            blk2_origin[0], blk2_origin[1], blk2_origin[2],
            blk2_dims[0], blk2_dims[1], blk2_dims[2]
        )

        all_vols = [(3, blk1_tag), (3, blk2_tag)]
        ov, ovm = gmsh.model.occ.fragment(all_vols, [])
        gmsh.model.occ.synchronize()

        new_blk1_tags = [entry[1] for entry in ovm[0]]
        new_blk2_tags = [entry[1] for entry in ovm[1]]

        gmsh.model.mesh.setSize(gmsh.model.getEntities(0), mesh_size)
        gmsh.model.occ.synchronize()
        gmsh.model.addPhysicalGroup(3, new_blk1_tags, 11)
        gmsh.model.setPhysicalName(3, 11, "moving-block")
        
        gmsh.model.addPhysicalGroup(3, new_blk2_tags, 21)
        gmsh.model.setPhysicalName(3, 21, "stationary-block")

        # ===== Contact Interface Handling =====
        tol = 1e-3
        # search contact candidate surfaces near x=200 over full y-range of stationary block
        y_min = 0.0 - tol
        y_max = 550.0 + tol
        contact_candidates = gmsh.model.getEntitiesInBoundingBox(
            200.0 - tol, y_min, 0.0 - tol,
            200.0 + tol, y_max, PMMA_thickness + tol,
            dim=2
        )

        if not contact_candidates:
            print("Error, cannot find CZM interface face candidates!")
        else:
            contact_tags = [tag for dim, tag in contact_candidates]
            print(f"[DEBUG] Found {len(contact_tags)} contact candidate surfaces: {contact_tags}")

            # collect boundary surfaces of each new volume (may be multiple after fragment)
            blk1_bound = {t for (d, t) in gmsh.model.getBoundary([(3, tag) for tag in new_blk1_tags], oriented=False) if d == 2}
            blk2_bound = {t for (d, t) in gmsh.model.getBoundary([(3, tag) for tag in new_blk2_tags], oriented=False) if d == 2}

            print(f"[DEBUG] blk1_bound: {blk1_bound}")
            print(f"[DEBUG] blk2_bound: {blk2_bound}")

            moving_front = [t for t in contact_tags if t in blk1_bound]
            stationary_back = [t for t in contact_tags if t in blk2_bound]
            unassigned = [t for t in contact_tags if t not in moving_front and t not in stationary_back]

            print(f"[DEBUG] moving_front tags: {moving_front}")
            print(f"[DEBUG] stationary_back tags: {stationary_back}")
            print(f"[DEBUG] unassigned tags: {unassigned}")

            if moving_front:
                gmsh.model.addPhysicalGroup(2, moving_front, 15)
                gmsh.model.setPhysicalName(2, 15, "moving-block-front")
                print(f"✓ Created moving-block-front with {len(moving_front)} surface(s)")
            else:
                print("⚠️ no moving-block-front surfaces found")

            if stationary_back:
                gmsh.model.addPhysicalGroup(2, stationary_back, 24)
                gmsh.model.setPhysicalName(2, 24, "stationary-block-back")
                print(f"✓ Created stationary-block-back with {len(stationary_back)} surface(s)")
            else:
                print("⚠️ no stationary-block-back surfaces found")

            # in case some faces are neither (safety)
            if unassigned:
                gmsh.model.addPhysicalGroup(2, unassigned, 99)
                gmsh.model.setPhysicalName(2, 99, "contact-interface-unassigned")
                print(f"⚠️ {len(unassigned)} unassigned contact surface(s)")

        # ===== Boundary Conditions =====
        front_face = gmsh.model.getEntitiesInBoundingBox(
            0 - tol, 0 - tol, 0 - tol,
            0 + tol, 500 + tol, PMMA_thickness + tol,
            dim=2
        )
        if front_face:
            gmsh.model.addPhysicalGroup(2, [front_face[0][1]], 14)
            gmsh.model.setPhysicalName(2, 14, "moving-block-back")
            print("✓ Created moving-block-back (x=0 face)")

        right_face = gmsh.model.getEntitiesInBoundingBox(
            0 - tol, 0 - tol, 0 - tol,
            200 + tol, 0 + tol, PMMA_thickness + tol,
            dim=2
        )
        if right_face:
            gmsh.model.addPhysicalGroup(2, [right_face[0][1]], 13)
            gmsh.model.setPhysicalName(2, 13, "moving-block-right")
            print("✓ Created moving-block-right (y=0 face)")

        left_face = gmsh.model.getEntitiesInBoundingBox(
            0 - tol, 500 - tol, 0 - tol,
            200 + tol, 500 + tol, PMMA_thickness + tol,
            dim=2
        )
        if left_face:
            gmsh.model.addPhysicalGroup(2, [left_face[0][1]], 16)
            gmsh.model.setPhysicalName(2, 16, "moving-block-left")
            print("✓ Created moving-block-left (y=500 face)")

        back_faces = gmsh.model.getEntitiesInBoundingBox(
            345 - tol, 0 - tol, 0 - tol,
            345 + tol, 550 + tol, PMMA_thickness + tol,
            dim=2
        )
        if back_faces:
            back_face_tags = [tag for dim, tag in back_faces]
            gmsh.model.addPhysicalGroup(2, back_face_tags, 25)
            gmsh.model.setPhysicalName(2, 25, "stationary-block-front")
            print(f"✓ Created stationary-block-front (x=345 face) with {len(back_face_tags)} surface(s)")
        else:
            print("⚠️ cannot find stationary-block-front")
        
        right_faces = gmsh.model.getEntitiesInBoundingBox(
            200 - tol, 0 - tol, 0 - tol,
            345 + tol, 0 + tol, PMMA_thickness + tol,
            dim=2
        )
        if right_faces:
            right_face_tags = [tag for dim, tag in right_faces]
            gmsh.model.addPhysicalGroup(2, right_face_tags, 26)
            gmsh.model.setPhysicalName(2, 26, "stationary-block-right")
            print(f"✓ Created stationary-block-right (y=0 face) with {len(right_face_tags)} surface(s)")
        else:
            print("⚠️ cannot find stationary-block-right")

        left_upper_faces = gmsh.model.getEntitiesInBoundingBox(
            200 - tol, 550 - tol, 0 - tol,
            345 + tol, 550 + tol, PMMA_thickness + tol,
            dim=2
        )
        if left_upper_faces:
            left_upper_face_tags = [tag for dim, tag in left_upper_faces]
            gmsh.model.addPhysicalGroup(2, left_upper_face_tags, 27)
            gmsh.model.setPhysicalName(2, 27, "stationary-block-left")
            print(f"✓ Created stationary-block-left (y=550 face) with {len(left_upper_face_tags)} surface(s)")
        else:
            print("⚠️ cannot find stationary-block-left")

        gmsh.model.occ.removeAllDuplicates()
        gmsh.model.occ.synchronize()
        gmsh.model.mesh.generate(3)
        gmsh.model.mesh.removeDuplicateNodes()
        gmsh.model.mesh.removeDuplicateElements()
        
        gmsh.write(f"../Models/{PMMA_thickness}mm-PMMA-CZM.msh")
        gmsh.write(f"../Models/{PMMA_thickness}mm-PMMA.brep")
        print(f"\n✓ Generated {PMMA_thickness}mm-PMMA-CZM.msh")
        print("✓ Contact interface (moving-block-front / stationary-block-back) has been created on conformal mesh.\n")
        
        # gmsh.fltk.run()
        gmsh.finalize()

if __name__ == "__main__":
    main()
