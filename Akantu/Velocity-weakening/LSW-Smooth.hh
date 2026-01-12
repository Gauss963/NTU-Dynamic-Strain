#ifndef LSW_SMOOTH_HH
#define LSW_SMOOTH_HH

// ============================================================
//  Linear Slip Weakening (LSW) friction law utilities
// ============================================================

// --- Probability / activation utilities ---

/// Standard normal cumulative distribution function Φ(x)
double normal_cdf(double x);

/// Gaussian Error Linear Unit
/// GELU(x) = x * Φ(x)
double GELU(double x);

// --- Linear Slip Weakening friction laws ---

/// Classic LSW friction coefficient (piecewise linear, non-smooth)
/// μ(δ) = μ_s + (μ_k - μ_s) * min(δ, D_c) / D_c
double LSW_mu(double delta);

/// Smoothed LSW friction coefficient
/// Uses GELU-based smoothing near δ = D_c
double LSW_mu_smooth(double delta);

#endif // LSW_SMOOTH_HH