import gmsh

def main():

    PMMA_THICKNESSES = [50, 100, 500]
    mesh_size = 5

    for PMMA_thickness in PMMA_THICKNESSES:

        gmsh.initialize()
        gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
        gmsh.option.setNumber("Mesh.Algorithm3D", 1)
        gmsh.option.setNumber("Mesh.ElementOrder", 1)
        gmsh.option.setNumber("Mesh.SubdivisionAlgorithm", 1) 
        gmsh.option.setNumber("Mesh.Binary", 0)

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

        new_blk1_tag = ovm[0][0][1]
        new_blk2_tag = ovm[1][0][1]

        gmsh.model.mesh.setSize(gmsh.model.getEntities(0), mesh_size)
        gmsh.model.occ.synchronize()
        gmsh.model.addPhysicalGroup(3, [new_blk1_tag], 11)
        gmsh.model.setPhysicalName(3, 11, "moving-block")
        
        gmsh.model.addPhysicalGroup(3, [new_blk2_tag], 21)
        gmsh.model.setPhysicalName(3, 21, "stationary-block")

        tol = 1e-3
        czm_face = gmsh.model.getEntitiesInBoundingBox(
            200 - tol, 0 - tol, 0 - tol,
            200 + tol, 500 + tol, PMMA_thickness + tol,
            dim=2
        )
        if not czm_face:
            print("Error, cannot find CZM interface face!")
        else:
            czm_face_tag = czm_face[0][1]
            
            slave_tag = 15 # tag_prefix * 10 + 5 (moving-block-back)
            master_tag = 24 # tag_prefix * 10 + 4 (stationary-block-front)
            shear_tag = 55
            
            gmsh.model.addPhysicalGroup(2, [czm_face_tag], slave_tag)
            gmsh.model.setPhysicalName(2, slave_tag, "friction-surface")
            
        front_face = gmsh.model.getEntitiesInBoundingBox(
            0 - tol, 0 - tol, 0 - tol,
            0 + tol, 500 + tol, PMMA_thickness + tol,
            dim=2
        )
        if front_face:
            gmsh.model.addPhysicalGroup(2, [front_face[0][1]], 14) # tag 14
            gmsh.model.setPhysicalName(2, 14, "moving-block-back")

        left_face = gmsh.model.getEntitiesInBoundingBox(
            0 - tol, 0 - tol, 0 - tol,
            200 + tol, 0 + tol, PMMA_thickness + tol,
            dim=2
        )
        if left_face:
            gmsh.model.addPhysicalGroup(2, [left_face[0][1]], 13) # tag 13
            gmsh.model.setPhysicalName(2, 13, "moving-block-right")
            
        back_faces = gmsh.model.getEntitiesInBoundingBox(
            345 - tol, 0 - tol, 0 - tol,
            345 + tol, 550 + tol, PMMA_thickness + tol,
            dim=2
        )
        if back_faces:
            back_face_tags = [tag for dim, tag in back_faces]
            gmsh.model.addPhysicalGroup(2, back_face_tags, 25)  # tag 25
            gmsh.model.setPhysicalName(2, 25, "stationary-block-front")
        else:
            print("⚠️ cannot find stationary-block-front")
        
        right_faces = gmsh.model.getEntitiesInBoundingBox(
            200 - tol, 550 - tol, 0 - tol,
            345 + tol, 550 + tol, PMMA_thickness + tol,
            dim=2
        )
        if right_faces:
            right_face_tags = [tag for dim, tag in right_faces]
            gmsh.model.addPhysicalGroup(2, right_face_tags, 26)  # tag 26
            gmsh.model.setPhysicalName(2, 26, "stationary-block-left")
        else:
            print("⚠️ cannot find stationary-block-left")

        gmsh.model.occ.removeAllDuplicates()
        gmsh.model.occ.synchronize()
        gmsh.model.mesh.generate(3)
        gmsh.model.mesh.removeDuplicateNodes()
        gmsh.model.mesh.removeDuplicateElements()
        
        gmsh.write(f"../Models/{PMMA_thickness}mm-PMMA-CZM.msh")
        gmsh.write(f"../Models/{PMMA_thickness}mm-PMMA.brep")
        print(f"Generated {PMMA_thickness}mm-PMMA-CZM.msh")
        print("CZM interface (Slave/Master) has been created on the conformal mesh.")
        
        gmsh.fltk.run()
        gmsh.finalize()

if __name__ == "__main__":
    main()