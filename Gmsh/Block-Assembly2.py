import gmsh
import F2 as Functions

def main():
    PMMA_THICKNESSES = [50] # 測試用
    mesh_size = 5

    for PMMA_thickness in PMMA_THICKNESSES:
        gmsh.initialize()
        gmsh.option.setNumber("Mesh.MshFileVersion", 2.2) 
        gmsh.option.setNumber("Mesh.ElementOrder", 1)
        gmsh.option.setNumber("Mesh.Algorithm3D", 1) 
        gmsh.model.add("ContactModel")

        blk1 = Functions.create_block(
            origin=(0, 0, 0),
            dimensions=(200, 500, PMMA_thickness),
            mesh_size=mesh_size,
            block_name="moving-block",
            tag_prefix=1
        )
        blk2a = Functions.create_block(
            origin=(200, 0, 0),
            dimensions=(145, 500, PMMA_thickness),
            mesh_size=mesh_size,
            block_name="stationary-block",
            tag_prefix=2
        )
        blk2b = Functions.create_block(
            origin=(200, 500, 0),
            dimensions=(145, 50, PMMA_thickness),
            mesh_size=mesh_size,
            block_name="stationary-block-upper",
            tag_prefix=3
        )

        gmsh.model.occ.synchronize()
        stationary_pg = gmsh.model.addPhysicalGroup(3, [blk2a["volume"], blk2b["volume"]], tag=221)
        gmsh.model.setPhysicalName(3, stationary_pg, "stationary-block")
        gmsh.model.mesh.generate(3)
        
        gmsh.write(f"../Models/{PMMA_thickness}mm-BS-PMMA.msh")
        gmsh.finalize()

if __name__ == "__main__":
    main()