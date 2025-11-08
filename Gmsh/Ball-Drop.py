import gmsh
import Functions

def main():

    PMMA_THICKNESSES = [20]
    mesh_size = 5

    for PMMA_thickness in PMMA_THICKNESSES:

        gmsh.initialize()
        gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
        # gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
        gmsh.option.setNumber("Mesh.Algorithm3D", 1)
        gmsh.option.setNumber("Mesh.RecombineAll", 1)
        gmsh.option.setNumber("Mesh.ElementOrder", 1)
        gmsh.option.setNumber("Mesh.SubdivisionAlgorithm", 1)
        gmsh.option.setNumber("Mesh.Binary", 0)

        gmsh.model.add("ContactModel")

        Sample = Functions.create_block(origin=(0, 0, 0), dimensions=(300, 300, PMMA_thickness), mesh_size=mesh_size, block_name="Sample", tag_prefix=1)
        


        gmsh.model.occ.synchronize()

        gmsh.model.mesh.generate(3)
        gmsh.write(f"../Models/Ball-Drop.msh")
        gmsh.write(f"../Models/Ball-Drop.brep")
        gmsh.fltk.run()
        gmsh.finalize()

if __name__ == "__main__":
    main()