import gmsh
import Functions

def main():

    mesh_size = 5

    gmsh.initialize()
    gmsh.option.setNumber("General.NumThreads", 10)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.ElementOrder", 1)
    gmsh.option.setNumber("Mesh.RecombineAll", 1)
    gmsh.option.setNumber("Mesh.Algorithm", 8)
    gmsh.option.setNumber("Mesh.Binary", 0)

    gmsh.model.add("ContactModel-2D-Quad")

    Functions.create_block_2d_quad(
        origin=(0, 0, 0),
        dimensions=(200, 500, 0),
        mesh_size=mesh_size,
        block_name="moving-block",
        tag_prefix=1
    )

    Functions.create_block_2d_quad(
        origin=(200, 0, 0),
        dimensions=(145, 550, 0),
        mesh_size=mesh_size,
        block_name="stationary-block",
        tag_prefix=2
    )

    gmsh.model.occ.synchronize()
    gmsh.model.mesh.generate(2)

    gmsh.write("../Models/2D-PMMA.msh")
    gmsh.write("../Models/2D-PMMA.brep")
    gmsh.fltk.run()
    gmsh.finalize()

if __name__ == "__main__":
    main()