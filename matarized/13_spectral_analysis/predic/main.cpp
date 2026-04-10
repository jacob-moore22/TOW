// Driver for MATAR predic -- linear prediction.
// Uses known AR(1) coefficients to predict a simple sequence and
// compares against direct computation.

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <matar.h>
#include "predic.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        // AR(1): x(n) = 0.8 * x(n-1), starting at x(1) = 1.0
        const int NDATA = 10;
        const int NPOLES = 1;
        const int NFUT = 8;

        DFMatrixKokkos<double> data(NDATA);
        double val = 1.0;
        for (int j = 1; j <= NDATA; j++) {
            data.host(j) = val;
            val *= 0.8;
        }
        data.update_device();

        DFMatrixKokkos<double> d(NPOLES);
        d.host(1) = 0.8;
        d.update_device();

        DFMatrixKokkos<double> future(NFUT);
        for (int j = 1; j <= NFUT; j++) future.host(j) = 0.0;
        future.update_device();

        predic(data, NDATA, d, NPOLES, future, NFUT);

        future.update_host();
        data.update_host();

        std::printf("Linear Prediction (AR(1), a=0.8)\n");
        std::printf("  Step   Predicted    Expected\n");

        double max_err = 0.0;
        double last = data.host(NDATA);
        for (int j = 1; j <= NFUT; j++) {
            double expected = last * std::pow(0.8, j);
            double err = std::fabs(future.host(j) - expected);
            if (err > max_err) max_err = err;
            std::printf("  %3d   %12.6f %12.6f\n", j, future.host(j), expected);
        }
        std::printf("\n  Max |predicted - expected| = %.2e\n\n", max_err);
    }
    MATAR_FINALIZE();
    return 0;
}
