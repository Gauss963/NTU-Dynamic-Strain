import gmsh
gmsh.initialize()
gmsh.open("../Models/50mm-BS-PMMA.msh")

types, elemTags, _ = gmsh.model.mesh.getElements(3)
for t, tags in zip(types, elemTags):
    name = gmsh.model.mesh.getElementProperties(t)[0]
    print(f"type {t:2d}: {name:20s}  count={len(tags)}")

gmsh.finalize()