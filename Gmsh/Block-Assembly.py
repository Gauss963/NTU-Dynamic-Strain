import gmsh
import Functions

def main():

    PMMA_THICKNESSES = [50, 100, 500]
    mesh_size = 2

    initial_offdet = 1 # initial gap, in mm

    for PMMA_thickness in PMMA_THICKNESSES:

        gmsh.initialize()
        gmsh.option.setNumber("General.NumThreads", 10)
        gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
        gmsh.option.setNumber("Mesh.ElementOrder", 1)
        gmsh.option.setNumber("Mesh.Algorithm3D", 1)

        # gmsh.option.setNumber("Mesh.RecombineAll", 1)
        gmsh.option.setNumber("Mesh.SubdivisionAlgorithm", 1)
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)

        gmsh.option.setNumber("Mesh.Binary", 0)
        gmsh.option.setNumber("Mesh.RecombineAll", 0)
        gmsh.model.add("ContactModel")

        blk1 = Functions.create_block(origin=(0, 0, 0), dimensions=(200, 500, PMMA_thickness), mesh_size=mesh_size, block_name="moving-block", tag_prefix=1)
        blk2 = Functions.create_block(origin=(200 + initial_offdet, 0, 0), dimensions=(145, 550, PMMA_thickness), mesh_size=mesh_size, block_name="stationary-block", tag_prefix=2)

        mov_face = blk1["faces_geo"]["moving-block-front"]
        sta_face = blk2["faces_geo"]["stationary-block-back"]

        gmsh.model.addPhysicalGroup(2, [mov_face], tag=1001)
        gmsh.model.setPhysicalName(2, 1001, "contact-moving")

        gmsh.model.addPhysicalGroup(2, [sta_face], tag=1002)
        gmsh.model.setPhysicalName(2, 1002, "contact-stationary")
        

        # blk1 = Functions.create_block(
        #     origin=(0, 0, 0),
        #     dimensions=(200, 500, PMMA_thickness),
        #     mesh_size=mesh_size,
        #     block_name="moving-block",
        #     tag_prefix=1
        # )
        # blk2a = Functions.create_block(
        #     origin=(200, 0, 0),
        #     dimensions=(145, 500, PMMA_thickness),
        #     mesh_size=mesh_size,
        #     block_name="stationary-block",
        #     tag_prefix=2
        # )
        # blk2b = Functions.create_block(
        #     origin=(200, 500, 0),
        #     dimensions=(145, 50, PMMA_thickness),
        #     mesh_size=mesh_size,
        #     block_name="stationary-block-upper",
        #     tag_prefix=3
        # )

        # gmsh.model.occ.synchronize()
        # stationary_pg = gmsh.model.addPhysicalGroup(3, [blk2a["volume"], blk2b["volume"]], tag=221)
        # gmsh.model.setPhysicalName(3, stationary_pg, "stationary-block")


        gmsh.option.setNumber("Mesh.Optimize", 1)
        gmsh.option.setNumber("Mesh.OptimizeNetgen", 1)
        gmsh.option.setNumber("Mesh.Smoothing", 120)
        gmsh.model.occ.synchronize()
        gmsh.model.mesh.generate(3)
        gmsh.model.mesh.optimize("Netgen")

        types3d, _, _ = gmsh.model.mesh.getElements(3)
        print("3D element types =", types3d)
        print("3D element types =", types3d)
        print("3D element types =", types3d)
        print("3D element types =", types3d)
        print("3D element types =", types3d)
        print("3D element types =", types3d)
        print("3D element types =", types3d)
        print("3D element types =", types3d)
        print("3D element types =", types3d)
        print("3D element types =", types3d)
        print("3D element types =", types3d)
        print("3D element types =", types3d)
        print("3D element types =", types3d)

        gmsh.write(f"../Models/{PMMA_thickness}mm-BS-PMMA.msh")
        gmsh.write(f"../Models/{PMMA_thickness}mm-BS-PMMA.brep")
        # gmsh.fltk.run()
        gmsh.finalize()

if __name__ == "__main__":
    main()