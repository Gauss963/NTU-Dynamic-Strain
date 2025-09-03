import numpy as np
import matplotlib.pyplot as plt

def alpha_s(C_f, C_s):
    return np.sqrt(1 - (C_f / C_s) ** 2)

def alpha_d(C_f, C_d):
    return np.sqrt(1 - (C_f / C_d) ** 2)

def D(alpha_s, alpha_d):
    return 4 * alpha_s * alpha_d - (1 + alpha_s ** 2) ** 2

def M_of_z(tau_p, X_c, z):
    return (2 / np.pi) * tau_p * ((1 + z / X_c) * np.arctan(1 / np.sqrt(z / X_c)) - np.sqrt(z / X_c))

def compute_A2(C_f, C_s, nu, D_value):
    alpha_s_value = alpha_s(C_f, C_s)
    psfactor = 1 / (1 - nu)
    return (C_f ** 2 * alpha_s_value * psfactor) / (C_s ** 2 * D_value)

def compute_K2(Gamma, E, nu, A2):
    return np.sqrt((Gamma * E) / ((1 - nu ** 2) * A2))

def compute_tau_p(K2, X_c):
    return K2 * np.sqrt(9 * np.pi / (32 * X_c))

def compute_stress_components(M_z_d, M_z_s, alpha_s_value, alpha_d_value):
    Sxx_tmp = (1 + 2 * alpha_d_value ** 2 - alpha_s_value ** 2) * M_z_d - (1 + alpha_s_value ** 2) * M_z_s
    Syy_tmp = M_z_d - M_z_s
    Sxy_tmp = 4 * alpha_s_value * alpha_d_value * M_z_d - (1 + alpha_s_value ** 2) ** 2 * M_z_s
    return Sxx_tmp, Syy_tmp, Sxy_tmp

def compute_stresses(Sxx_tmp, Syy_tmp, Sxy_tmp, alpha_s_value, D_value):
    Sxx = 2 * alpha_s_value / D_value * Sxx_tmp.imag
    Syy = -2 * alpha_s_value * (1 + alpha_s_value ** 2) / D_value * Syy_tmp.imag
    Sxy = Sxy_tmp.real / D_value
    return Sxx, Syy, Sxy

def delta_sigma_xy(x, y, X_c, C_f, C_s, C_d, nu, Gamma, E):
    alpha_s_value = alpha_s(C_f, C_s)
    alpha_d_value = alpha_d(C_f, C_d)
    D_value = D(alpha_s_value, alpha_d_value)
    A2 = compute_A2(C_f, C_s, nu, D_value)
    K2 = compute_K2(Gamma, E, nu, A2)
    tau_p = compute_tau_p(K2, X_c)
    
    z_d_value = x + 1j * alpha_d_value * y
    z_s_value = x + 1j * alpha_s_value * y
    
    M_z_d = M_of_z(tau_p, X_c, z_d_value)
    M_z_s = M_of_z(tau_p, X_c, z_s_value)
    
    Sxx_tmp, Syy_tmp, Sxy_tmp = compute_stress_components(M_z_d, M_z_s, alpha_s_value, alpha_d_value)
    
    Sxx, Syy, Sxy = compute_stresses(Sxx_tmp, Syy_tmp, Sxy_tmp, alpha_s_value, D_value)
    
    delta_sigma = Sxy
    
    return delta_sigma

def delta_sigma_xx(x, y, X_c, C_f, C_s, C_d, nu, Gamma, E):
    alpha_s_value = alpha_s(C_f, C_s)
    alpha_d_value = alpha_d(C_f, C_d)
    D_value = D(alpha_s_value, alpha_d_value)
    A2 = compute_A2(C_f, C_s, nu, D_value)
    K2 = compute_K2(Gamma, E, nu, A2)
    tau_p = compute_tau_p(K2, X_c)
    
    z_d_value = x + 1j * alpha_d_value * y
    z_s_value = x + 1j * alpha_s_value * y
    
    M_z_d = M_of_z(tau_p, X_c, z_d_value)
    M_z_s = M_of_z(tau_p, X_c, z_s_value)
    
    Sxx_tmp, Syy_tmp, Sxy_tmp = compute_stress_components(M_z_d, M_z_s, alpha_s_value, alpha_d_value)
    
    Sxx, Syy, Sxy = compute_stresses(Sxx_tmp, Syy_tmp, Sxy_tmp, alpha_s_value, D_value)
    
    delta_sigma = Sxx
    
    return delta_sigma