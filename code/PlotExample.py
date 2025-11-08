import numpy as np
import matplotlib.pyplot as plt

import FolderActions
import CohesiveCrack

def main():

    materials = FolderActions.read_materials("../Materials/material-mm-MPa.dat")

    Gamma = materials["interface_mat"]["parameters"]["G_c"]  # Fracture energy (J/m^2)
    E = materials["moving-block"]["parameters"]["E"]         # Young's modulus (MPa)
    nu = materials["moving-block"]["parameters"]["nu"]       # Poisson's ratio


    Gamma = 0.21  # Fracture energy (J/m^2)
    E = 51e9      # Young's modulus (Pa)
    nu = 0.25     # Poisson's ratio
    C_f = 2404    # Rupture speed (m/s)
    C_s = 2760    # Shear wave speed (m/s)
    C_d = 4790    # Longitudinal wave speed (m/s)
    X_c = 13.8e-3 # Cohesive zone size (m)
    y_values = [1e-8, 0.1e-3, 0.5e-3, 1.0e-3, 2.0e-3, 5e-3, 10e-3, 15e-3]

    x = np.linspace(-50e-3, 50e-3, 8192)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True)

    for i, y in enumerate(y_values):
        delta_sigma_xx = CohesiveCrack.delta_sigma_xx(x, y, X_c, C_f, C_s, C_d, nu, Gamma, E)
        axes[0].plot(x * 1000, delta_sigma_xx / 1e5 + i * 5, '-.', label=f'y = {y * 1e3:.1f} mm')

    axes[0].set_xlabel('Rupture tip position x (mm)')
    axes[0].set_ylabel('Shear stress fluctuation $\\Delta \\sigma_{xx}$ (MPa)')
    axes[0].set_title('Normal stress fluctuation (xx)')
    axes[0].axvline(0, color='k', linestyle='--', linewidth=1)
    axes[0].legend(loc = "lower right")
    axes[0].grid()

    for i, y in enumerate(y_values):
        delta_sigma_xy = CohesiveCrack.delta_sigma_xy(x, y, X_c, C_f, C_s, C_d, nu, Gamma, E)
        axes[1].plot(x * 1000, delta_sigma_xy / 1e5 + i * 5, '-.', label=f'y = {y * 1e3:.1f} mm')

    axes[1].set_xlabel('Rupture tip position x (mm)')
    axes[1].set_ylabel('Shear stress fluctuation $\\Delta \\sigma_{xy}$ (MPa)')
    axes[1].set_title('Shear stress fluctuation (xy)')
    axes[1].axvline(0, color='k', linestyle='--', linewidth=1)
    axes[1].legend(loc = "lower right")
    axes[1].grid()

    plt.suptitle(f'Stress fluctuations along the fault | E = {E/1e9}GPa | ν = {nu} | $C_f$ = {C_f}m/s, $C_s$ = {C_s}m/s, $C_d$ = {C_d}m/s', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig('./Plot/example_xx_xy.pdf', dpi=900)
    # plt.show()
    
    return 0

if __name__ == "__main__":
    main()
    FolderActions.delete_pycache()