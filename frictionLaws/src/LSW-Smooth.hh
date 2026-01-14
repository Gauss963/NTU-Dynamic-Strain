#pragma once
#include <algorithm>
#include <cmath>

double normal_cdf(double x) {
    return 0.5 * (1.0 + std::erf(x / std::sqrt(2.0)));
}

double GELU(double x) {
    return x * normal_cdf(x);
}

double LSW_mu(double delta, double mu_s, double mu_k, double D_c) {
    double delta_eff = std::min(delta, D_c);
    return delta_eff * (mu_k - mu_s) / D_c + mu_s;
}

double LSW_mu_smooth(double delta, double mu_s, double mu_k, double D_c) {
    
    const double eps = 0.01 * D_c;
    double x = (delta - D_c) / eps;
    double delta_eff = delta - eps * GELU(x);

    return mu_s + (mu_k - mu_s) * delta_eff / D_c;
}