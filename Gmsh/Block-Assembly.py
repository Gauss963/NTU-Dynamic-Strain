import gmsh
import Functions

def main():

    PMMA_THICKNESSES = [50, 100]
    mesh_size = 5

    for PMMA_thickness in PMMA_THICKNESSES:

        gmsh.initialize()
        # gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
        gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
        gmsh.option.setNumber("Mesh.Algorithm3D", 1)
        gmsh.option.setNumber("Mesh.RecombineAll", 1)
        gmsh.option.setNumber("Mesh.ElementOrder", 1)
        gmsh.option.setNumber("Mesh.SubdivisionAlgorithm", 1)
        gmsh.option.setNumber("Mesh.Binary", 0)

        gmsh.model.add("ContactModel")

        blk1 = Functions.create_block(origin=(0, 0, 0), dimensions=(200, 500, PMMA_thickness), mesh_size=mesh_size, block_name="moving-block", tag_prefix=1)
        blk2 = Functions.create_block(origin=(200, 0, 0), dimensions=(145, 550, PMMA_thickness), mesh_size=mesh_size, block_name="stationary-block", tag_prefix=2)

        gmsh.model.setPhysicalName(2, blk1["faces_phys"]["moving-block-back"],  "friction_master")
        gmsh.model.setPhysicalName(2, blk2["faces_phys"]["stationary-block-front"], "friction_slave")


        slave_pg  = gmsh.model.addPhysicalGroup(2, [blk1["faces_geo"]["moving-block-back"]])
        master_pg = gmsh.model.addPhysicalGroup(2, [blk2["faces_geo"]["stationary-block-front"]])
        gmsh.model.setPhysicalName(2, slave_pg,  "friction_slave")
        gmsh.model.setPhysicalName(2, master_pg, "friction_master")


        gmsh.model.occ.synchronize()

        gmsh.model.mesh.generate(3)
        gmsh.write(f"../Models/{PMMA_thickness}mm-PMMA.msh")
        # gmsh.fltk.run()
        gmsh.finalize()

if __name__ == "__main__":
    main()