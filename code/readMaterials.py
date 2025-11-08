import FolderActions

mat_file = "../Materials/material-mm-MPa.dat"
materials = FolderActions.read_materials(mat_file)

print(materials)

Gamma = materials["interface_mat"]["parameters"]["G_c"]  # Fracture energy (J/m^2)
E = materials["moving-block"]["parameters"]["E"]         # Young's modulus (MPa)
nu = materials["moving-block"]["parameters"]["nu"]       # Poisson's ratio


print(f"Fracture energy (Gamma): {Gamma} J/m^2")
print(f"Young's modulus (E): {E} MPa")
print(f"Poisson's ratio (nu): {nu}")


FolderActions.delete_pycache()