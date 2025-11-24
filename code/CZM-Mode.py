def convert_modeII_to_akantu(tau_c, G_IIc, beta):
    """
    Convert Mode II cohesive properties (tau_c, G_IIc)
    into Akantu parameters (sigma_c, fracture energy, beta).

    Inputs:
        tau_c : Mode II critical shear strength (MPa)
        G_IIc : Mode II fracture energy (N/mm)
        beta  : Akantu mode-mixity parameter (given by user)

    Output:
        sigma_c : Akantu's critical cohesive strength (MPa)
        G_c     : fracture energy = G_IIc (N/mm)
        gamma_c : computed shear displacement at failure (mm) for checking
    """

    # Akantu defines effective stress: sigma_eff = sqrt( σ_n^2 + τ^2/β^2 )
    # For pure Mode II: τ_c = β * σ_c
    sigma_c = tau_c / beta

    # Akantu's fracture energy is the Mode II fracture energy
    G_c = G_IIc

    # For verification: G_IIc = 0.5 * τ_c * γ_c
    gamma_c = 2 * G_IIc / tau_c

    return sigma_c, G_c, beta, gamma_c


# Example usage:
tau_c = 3.0     # MPa
G_IIc = 0.12     # N/mm
beta  = 0.01     # chosen for Akantu

sigma_c, G_c, beta_out, gamma_c = convert_modeII_to_akantu(tau_c, G_IIc, beta)

print("Converted parameters for Akantu:")
print(f"  sigma_c  = {sigma_c:.4f} MPa")
print(f"  G_c      = {G_c:.4f} N/mm  (Akantu fracture energy)")
print(f"  beta     = {beta_out}")
print(f"  gamma_c  = {gamma_c:.6f} mm (for verification only)")