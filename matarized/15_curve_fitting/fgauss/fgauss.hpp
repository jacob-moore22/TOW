// Gaussian basis functions with derivatives (Numerical Recipes FGAUSS).
//
// Evaluates y = sum of Gaussians and partial derivatives dyda for
// Levenberg-Marquardt fitting. Parameters are grouped in triplets:
//   a(i)   = amplitude
//   a(i+1) = center
//   a(i+2) = width
// Each Gaussian: a(i) * exp(-((x - a(i+1)) / a(i+2))^2)

#pragma once
#include <cmath>
#include <matar.h>

using namespace mtr;

inline void fgauss(double x, DFMatrixKokkos<double>& a, double& y,
                   DFMatrixKokkos<double>& dyda, int na)
{
    y = 0.0;
    for (int i = 1; i <= na - 1; i += 3) {
        double arg = (x - a.host(i + 1)) / a.host(i + 2);
        double ex  = std::exp(-arg * arg);
        double fac = a.host(i) * ex * 2.0 * arg;
        y += a.host(i) * ex;
        dyda.host(i)     = ex;
        dyda.host(i + 1) = fac / a.host(i + 2);
        dyda.host(i + 2) = fac * arg / a.host(i + 2);
    }
}
