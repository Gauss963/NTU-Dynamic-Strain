import CohesiveCrackPY
import FolderActions

mat_file = "../Materials/material-mm-MPa.dat"
materials = FolderActions.read_materials(mat_file)

print(materials)

Gamma = materials["interface"]["parameters"]["G_c"]      # Fracture energy (J/m^2)
E = materials["moving-block"]["parameters"]["E"]         # Young's modulus (MPa)
nu = materials["moving-block"]["parameters"]["nu"]       # Poisson's ratio
rho = materials["moving-block"]["parameters"]["rho"]     # Density (tonne/mm^3)


Vp = 3500  # P-wave velocity (m/s) - Example value
Vs = 1700  # S-wave velocity (m/s) - Example value

E, nu = CohesiveCrackPY.compute_E_nu_from_VpVsRho(Vp, Vs, rho * 1e12)  # Convert density to kg/m^3

print(f"Fracture energy (Gamma): {Gamma} J/m^2")
print(f"Young's modulus (E): {E/1e9} GPa")
print(f"Poisson's ratio (nu): {nu}")
print(f"Density (rho): {rho} tonne/mm^3")


FolderActions.delete_pycache()