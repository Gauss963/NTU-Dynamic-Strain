#include <TAxis.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TROOT.h>
#include <iostream>
#include <vector>

#include "LSW-Smooth.hh"

int main() {
    gROOT->SetBatch(true);

    constexpr double mu_s = 0.8;
    constexpr double mu_k = 0.6;
    constexpr double D_c = 5.5;

    const int N = 8192;
    const double TOTAL_SLIP = 10.0;
    std::vector<double> x(N), mu(N), mu_smooth(N);

    for (int i = 0; i < N; ++i) {
        x[i] = i * TOTAL_SLIP / (N - 1);
        mu_smooth[i] = LSW_mu_smooth(x[i], mu_s, mu_k, D_c);
        mu[i] = LSW_mu(x[i], mu_s, mu_k, D_c);
    }

    TCanvas canvas("c", "LSW", 800, 600);
    TGraph graph(N, x.data(), mu_smooth.data());

    graph.SetLineWidth(2);
    graph.SetTitle("Linear Slip Weakening;Slip #delta;Friction #mu");

    graph.Draw("AL");
    canvas.SaveAs("../../Plot/lsw_mu.pdf");

    std::cout << "Done" << std::endl;

    return 0;
}