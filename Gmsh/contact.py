import gmsh
import Functions

def main():
    gmsh.initialize()
    gmsh.model.add("ContactModel")

    mesh_size = 20

    Functions.create_block(origin=(0, 0, 0), dimensions=(200, 500, 50), mesh_size=mesh_size, tag_prefix=1)
    Functions.create_block(origin=(200, 0, 0), dimensions=(145, 550, 50), mesh_size=mesh_size, tag_prefix=2)

    gmsh.model.geo.synchronize()

    gmsh.option.setNumber("Mesh.Algorithm3D", 1)
    gmsh.option.setNumber("Mesh.RecombineAll", 1)
    gmsh.option.setNumber("Mesh.ElementOrder", 1)

    gmsh.model.mesh.generate(3)
    gmsh.write("../Models/5cm-PMMA.msh")
    gmsh.fltk.run()
    gmsh.finalize()

if __name__ == "__main__":
    main()