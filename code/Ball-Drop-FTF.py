import numpy as np
import matplotlib.pyplot as plt
import CohesiveCrackPY
import FolderActions

def stf_herzian_mclaskey2009(t, rho, R, v, E1, nu1, E2, nu2):
    def get_delta(E, nu):
        return (1 - nu ** 2) / (np.pi * E)

    del1 = get_delta(E1, nu1)
    del2 = get_delta(E2, nu2)

    tc = 4.53 * ((4 * rho * np.pi * (del1 + del2) / 3) ** (2 / 5)) * R * (v ** (-1 / 5))
    fmax = 1.917 * (rho ** (3 / 5)) * ((del1 + del2) ** (-2 / 5)) * (R ** 2) * (v ** (6 / 5))

    stf = np.zeros_like(t)
    for i, tt in enumerate(t):
        if 0 < tt < tc:
            stf[i] = fmax * (np.sin(np.pi * tt / tc)) ** (3 / 2)
    return stf, fmax, tc

# ========================================
# 1. 時間設定
# ========================================
dt = 1e-7  # 0.1 μs
tvec = np.arange(np.ceil(20e-6 / dt)) * dt  # 20 μs 時間軸

# ========================================
# 2. 掉落參數：h 決定撞擊速度
# ========================================
h = 0.300  # [m]
g = 9.80665
v = np.sqrt(2 * g * h)  # 撞擊速度

# ========================================
# 3. 載入 PMMA 材料參數（你的 moving-block）
# ========================================
materials = FolderActions.read_materials("../Materials/material-mm-MPa.dat")
rho_pmma = materials["moving-block"]["parameters"]["rho"] * 1e9      # g/mm³ → kg/m³
E_pmma   = materials["moving-block"]["parameters"]["E"]   * 1e6      # MPa → Pa
nu_pmma  = materials["moving-block"]["parameters"]["nu"]             # Poisson's ratio

# 可選：計算剪力波速度（僅用於確認材料性質）
C_s = CohesiveCrackPY.get_Cs(E_pmma, nu_pmma, rho_pmma)
C_d = CohesiveCrackPY.get_Cd(E_pmma, nu_pmma, rho_pmma)

# ========================================
# 4. 鋼球參數（McLaskey 用的鋼球）
# ========================================
E_ball  = 208.197e9  # Pa
nu_ball = 0.286

# ========================================
# 5. 計算 Hertzian 接觸的 STF
# ========================================
ft_herz, fmax, tc = stf_herzian_mclaskey2009(
    tvec,
    rho=rho_pmma, R=1e-3, v=v,
    E1=E_ball, nu1=nu_ball,     # 球
    E2=E_pmma, nu2=nu_pmma      # PMMA 材料（你的 block）
)

# ========================================
# 6. Scaling → 得到 M₀ 一致的幅值
# ========================================
M0_scaling = 1.748 * fmax * tc / np.pi

# ========================================
# 7. 顯示資訊
# ========================================
print(f"v     = {v:.3f} m/s")
print(f"C_s   = {C_s:.1f} m/s, C_d = {C_d:.1f} m/s")
print(f"fmax  = {fmax:.2f} N")
print(f"tc    = {tc*1e6:.2f} μs")
print(f"Number of non-zero STF = {np.count_nonzero(ft_herz)}")

# ========================================
# 8. 畫圖
# ========================================
fig, ax = plt.subplots(figsize=(5, 4))
ax.plot(tvec * 1e6, ft_herz * M0_scaling, "b--", label="Hertzian STF (PMMA)")
ax.set_xlim([0, 20])
ax.set_ylim([0, None])
ax.set_xlabel("Time [μs]")
ax.set_ylabel("Force [N]")
ax.legend()
plt.tight_layout()
plt.show()