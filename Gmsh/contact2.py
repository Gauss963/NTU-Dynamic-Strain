import gmsh
import Functions
import Test
import Test2

def main():

    # PMMA_THICKNESSES = [50, 100]
    PMMA_THICKNESSES = [100]

    for PMMA_thickness in PMMA_THICKNESSES:

        gmsh.initialize()
        gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
        gmsh.option.setNumber("Mesh.Algorithm3D", 1)
        gmsh.option.setNumber("Mesh.RecombineAll", 1)
        gmsh.option.setNumber("Mesh.ElementOrder", 1)
        gmsh.option.setNumber("Mesh.SubdivisionAlgorithm", 1)
        gmsh.option.setNumber("Mesh.Binary", 0)  # optional, ASCII for debug

        gmsh.model.add("ContactModel")

        mesh_size = 20

        blk1 = Test2.create_block(origin=(0, 0, 0), dimensions=(200, 500, PMMA_thickness), mesh_size=mesh_size, tag_prefix=1)
        blk2 = Test2.create_block(origin=(200, 0, 0), dimensions=(145, 550, PMMA_thickness), mesh_size=mesh_size, tag_prefix=2)

        gmsh.model.addPhysicalGroup(2, [blk1["faces"]["back"]], tag=100)
        gmsh.model.setPhysicalName(2, 100, "friction_master")

        gmsh.model.addPhysicalGroup(2, [blk2["faces"]["front"]], tag=101)
        gmsh.model.setPhysicalName(2, 101, "friction_slave")

        # gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
        # gmsh.option.setNumber("Mesh.Algorithm3D", 1)
        # gmsh.option.setNumber("Mesh.RecombineAll", 1)
        # gmsh.option.setNumber("Mesh.ElementOrder", 1)

        gmsh.model.mesh.generate(3)
        gmsh.write(f"../Models/{PMMA_thickness}mm-PMMA.msh")
        gmsh.finalize()

if __name__ == "__main__":
    main()