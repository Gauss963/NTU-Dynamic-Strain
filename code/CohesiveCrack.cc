#include <iostream>
#include <complex>
#include <cmath>
#include <vector>
#include <chrono>
#include <iomanip>

class StressAnalysis {
private:
    static constexpr double PI = M_PI;
    
public:
    static double alpha_s(double C_f, double C_s) {
        return std::sqrt(1.0 - (C_f / C_s) * (C_f / C_s));
    }
    
    static double alpha_d(double C_f, double C_d) {
        return std::sqrt(1.0 - (C_f / C_d) * (C_f / C_d));
    }
    
    static double D(double alpha_s_val, double alpha_d_val) {
        double term = 1.0 + alpha_s_val * alpha_s_val;
        return 4.0 * alpha_s_val * alpha_d_val - term * term;
    }
    
    static std::complex<double> M_of_z(double tau_p, double X_c, const std::complex<double>& z) {
        std::complex<double> ratio = z / X_c;
        std::complex<double> sqrt_ratio = std::sqrt(ratio);
        std::complex<double> one_plus_ratio = 1.0 + ratio;
        
        std::complex<double> arctan_arg = 1.0 / sqrt_ratio;
        std::complex<double> arctan_result = std::atan(arctan_arg);
        
        return (2.0 / PI) * tau_p * (one_plus_ratio * arctan_result - sqrt_ratio);
    }
    
    static double compute_A2(double C_f, double C_s, double nu, double D_value) {
        double alpha_s_value = alpha_s(C_f, C_s);
        double psfactor = 1.0 / (1.0 - nu);
        return (C_f * C_f * alpha_s_value * psfactor) / (C_s * C_s * D_value);
    }
    
    static double compute_K2(double Gamma, double E, double nu, double A2) {
        return std::sqrt((Gamma * E) / ((1.0 - nu * nu) * A2));
    }
    
    static double compute_tau_p(double K2, double X_c) {
        return K2 * std::sqrt(9.0 * PI / (32.0 * X_c));
    }
    
    static void compute_stress_components(
        const std::complex<double>& M_z_d,
        const std::complex<double>& M_z_s,
        double alpha_s_value,
        double alpha_d_value,
        std::complex<double>& Sxx_tmp,
        std::complex<double>& Syy_tmp,
        std::complex<double>& Sxy_tmp
    ) {
        double alpha_s_sq = alpha_s_value * alpha_s_value;
        double alpha_d_sq = alpha_d_value * alpha_d_value;
        double term1 = 1.0 + alpha_s_sq;
        
        Sxx_tmp = (1.0 + 2.0 * alpha_d_sq - alpha_s_sq) * M_z_d - term1 * M_z_s;
        Syy_tmp = M_z_d - M_z_s;
        Sxy_tmp = 4.0 * alpha_s_value * alpha_d_value * M_z_d - term1 * term1 * M_z_s;
    }
    
    static void compute_stresses(
        const std::complex<double>& Sxx_tmp,
        const std::complex<double>& Syy_tmp,
        const std::complex<double>& Sxy_tmp,
        double alpha_s_value,
        double D_value,
        double& Sxx,
        double& Syy,
        double& Sxy
    ) {
        double alpha_s_sq = alpha_s_value * alpha_s_value;
        
        Sxx = 2.0 * alpha_s_value / D_value * Sxx_tmp.imag();
        Syy = -2.0 * alpha_s_value * (1.0 + alpha_s_sq) / D_value * Syy_tmp.imag();
        Sxy = Sxy_tmp.real() / D_value;
    }
    
    static double delta_sigma_xy(
        double x, double y, double X_c, double C_f, double C_s, 
        double C_d, double nu, double Gamma, double E
    ) {
        double alpha_s_value = alpha_s(C_f, C_s);
        double alpha_d_value = alpha_d(C_f, C_d);
        double D_value = D(alpha_s_value, alpha_d_value);
        double A2 = compute_A2(C_f, C_s, nu, D_value);
        double K2 = compute_K2(Gamma, E, nu, A2);
        double tau_p = compute_tau_p(K2, X_c);
        
        std::complex<double> z_d_value(x, alpha_d_value * y);
        std::complex<double> z_s_value(x, alpha_s_value * y);
        
        std::complex<double> M_z_d = M_of_z(tau_p, X_c, z_d_value);
        std::complex<double> M_z_s = M_of_z(tau_p, X_c, z_s_value);
        
        std::complex<double> Sxx_tmp, Syy_tmp, Sxy_tmp;
        compute_stress_components(M_z_d, M_z_s, alpha_s_value, alpha_d_value, 
                                Sxx_tmp, Syy_tmp, Sxy_tmp);
        
        double Sxx, Syy, Sxy;
        compute_stresses(Sxx_tmp, Syy_tmp, Sxy_tmp, alpha_s_value, D_value, 
                        Sxx, Syy, Sxy);
        
        return Sxy;
    }
    
    static double delta_sigma_xx(
        double x, double y, double X_c, double C_f, double C_s, 
        double C_d, double nu, double Gamma, double E
    ) {
        double alpha_s_value = alpha_s(C_f, C_s);
        double alpha_d_value = alpha_d(C_f, C_d);
        double D_value = D(alpha_s_value, alpha_d_value);
        double A2 = compute_A2(C_f, C_s, nu, D_value);
        double K2 = compute_K2(Gamma, E, nu, A2);
        double tau_p = compute_tau_p(K2, X_c);
        
        std::complex<double> z_d_value(x, alpha_d_value * y);
        std::complex<double> z_s_value(x, alpha_s_value * y);
        
        std::complex<double> M_z_d = M_of_z(tau_p, X_c, z_d_value);
        std::complex<double> M_z_s = M_of_z(tau_p, X_c, z_s_value);
        
        std::complex<double> Sxx_tmp, Syy_tmp, Sxy_tmp;
        compute_stress_components(M_z_d, M_z_s, alpha_s_value, alpha_d_value, 
                                Sxx_tmp, Syy_tmp, Sxy_tmp);
        
        double Sxx, Syy, Sxy;
        compute_stresses(Sxx_tmp, Syy_tmp, Sxy_tmp, alpha_s_value, D_value, 
                        Sxx, Syy, Sxy);
        
        return Sxx;
    }

    static void benchmark_test() {
        std::cout << "Running benchmark test..." << std::endl;
        
        double x = 1.0, y = 2.0, X_c = 10.0;
        double C_f = 1000.0, C_s = 3000.0, C_d = 5000.0;
        double nu = 0.3, Gamma = 100.0, E = 200000.0;
        
        auto start = std::chrono::high_resolution_clock::now();
        
        const int iterations = 100000;
        double sum_xy = 0.0, sum_xx = 0.0;
        
        for (int i = 0; i < iterations; ++i) {
            double xi = x + i * 0.001;
            double yi = y + i * 0.001;
            sum_xy += delta_sigma_xy(xi, yi, X_c, C_f, C_s, C_d, nu, Gamma, E);
            sum_xx += delta_sigma_xx(xi, yi, X_c, C_f, C_s, C_d, nu, Gamma, E);
        }
        
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        
        std::cout << std::fixed << std::setprecision(6);
        std::cout << "Benchmark results:" << std::endl;
        std::cout << "Iterations: " << iterations << std::endl;
        std::cout << "Total time: " << duration.count() / 1000.0 << " ms" << std::endl;
        std::cout << "Time per iteration: " << duration.count() / (double)iterations << " μs" << std::endl;
        std::cout << "Sum xy: " << sum_xy << std::endl;
        std::cout << "Sum xx: " << sum_xx << std::endl;
    }
};

int main() {
    std::cout << "=== Stress Analysis C++ Implementation ===" << std::endl;
    
    double x = 1.0, y = 2.0, X_c = 10.0;
    double C_f = 1000.0, C_s = 3000.0, C_d = 5000.0;
    double nu = 0.3, Gamma = 100.0, E = 200000.0;
    
    std::cout << "\nSingle calculation test:" << std::endl;
    std::cout << "Parameters: x=" << x << ", y=" << y << ", X_c=" << X_c << std::endl;
    std::cout << "C_f=" << C_f << ", C_s=" << C_s << ", C_d=" << C_d << std::endl;
    std::cout << "nu=" << nu << ", Gamma=" << Gamma << ", E=" << E << std::endl;
    
    double result_xy = StressAnalysis::delta_sigma_xy(x, y, X_c, C_f, C_s, C_d, nu, Gamma, E);
    double result_xx = StressAnalysis::delta_sigma_xx(x, y, X_c, C_f, C_s, C_d, nu, Gamma, E);
    
    std::cout << "\nResults:" << std::endl;
    std::cout << "delta_sigma_xy: " << result_xy << std::endl;
    std::cout << "delta_sigma_xx: " << result_xx << std::endl;
    
    std::cout << "\n" << std::endl;
    StressAnalysis::benchmark_test();
    
    return 0;
}