import numpy as np
import matplotlib.pyplot as plt

import FolderActions
import CohesiveCrack  # pybind11 模組
import CohesiveCrackPY

def main():
    # === 包裝 C++ 函數為 NumPy vectorized ===
    cpp_delta_sigma_xx = np.vectorize(CohesiveCrack.delta_sigma_xx)
    cpp_delta_sigma_xy = np.vectorize(CohesiveCrack.delta_sigma_xy)

    # === 參數設定 ===
    Gamma = 0.21
    E = 51e9
    nu = 0.25
    C_f = 2404
    C_s = 2760
    C_d = 4790
    X_c = 13.8e-3
    y_values = [1e-8, 0.1e-3, 0.5e-3, 1.0e-3, 2.0e-3, 5e-3, 10e-3, 15e-3]
    x = np.linspace(-50e-3, 50e-3, 8192)

    fig, axes = plt.subplots(2, 2, figsize=(14, 10), sharex=True, sharey=True)

    # ======== 第一列：Python 版 ========
    for i, y in enumerate(y_values):
        delta_sigma_xx_py = CohesiveCrackPY.delta_sigma_xx(x, y, X_c, C_f, C_s, C_d, nu, Gamma, E)
        axes[0][0].plot(x * 1000, delta_sigma_xx_py / 1e5 + i * 5, '--', label=f'y = {y * 1e3:.1f} mm')

    axes[0][0].set_title('Python: $\\Delta \\sigma_{xx}$')
    axes[0][0].set_ylabel('Stress (MPa)')
    axes[0][0].axvline(0, color='k', linestyle='--')
    axes[0][0].grid(True)
    axes[0][0].legend(loc='lower right')

    for i, y in enumerate(y_values):
        delta_sigma_xy_py = CohesiveCrackPY.delta_sigma_xy(x, y, X_c, C_f, C_s, C_d, nu, Gamma, E)
        axes[0][1].plot(x * 1000, delta_sigma_xy_py / 1e5 + i * 5, '--', label=f'y = {y * 1e3:.1f} mm')

    axes[0][1].set_title('Python: $\\Delta \\sigma_{xy}$')
    axes[0][1].axvline(0, color='k', linestyle='--')
    axes[0][1].grid(True)
    axes[0][1].legend(loc='lower right')

    # ======== 第二列：C++ pybind11 版 ========
    for i, y in enumerate(y_values):
        delta_sigma_xx_cpp = cpp_delta_sigma_xx(x, y, X_c, C_f, C_s, C_d, nu, Gamma, E)
        axes[1][0].plot(x * 1000, delta_sigma_xx_cpp / 1e5 + i * 5, '-', label=f'y = {y * 1e3:.1f} mm')

    axes[1][0].set_title('C++ (pybind11): $\\Delta \\sigma_{xx}$')
    axes[1][0].set_xlabel('x (mm)')
    axes[1][0].set_ylabel('Stress (MPa)')
    axes[1][0].axvline(0, color='k', linestyle='--')
    axes[1][0].grid(True)

    for i, y in enumerate(y_values):
        delta_sigma_xy_cpp = cpp_delta_sigma_xy(x, y, X_c, C_f, C_s, C_d, nu, Gamma, E)
        axes[1][1].plot(x * 1000, delta_sigma_xy_cpp / 1e5 + i * 5, '-', label=f'y = {y * 1e3:.1f} mm')

    axes[1][1].set_title('C++ (pybind11): $\\Delta \\sigma_{xy}$')
    axes[1][1].set_xlabel('x (mm)')
    axes[1][1].axvline(0, color='k', linestyle='--')
    axes[1][1].grid(True)

    # === 全圖設定 ===
    fig.suptitle(f'Stress fluctuation comparison: Python vs C++', fontsize=16, fontweight='bold')
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig('../Plot/compare_py_cpp_xx_xy.pdf', dpi=900)
    # plt.show()

    FolderActions.delete_pycache()

if __name__ == "__main__":
    main()