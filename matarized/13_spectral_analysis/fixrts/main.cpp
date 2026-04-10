// Driver for MATAR fixrts -- stabilise LP roots.
// Creates LP coefficients with a known unstable root, fixes them, and
// verifies all roots lie inside the unit circle.

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <matar.h>
#include "fixrts.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        // Build LP coefficients from roots: r1 = 0.5, r2 = 1.5 (unstable)
        // Polynomial: (z - 0.5)(z - 1.5) = z^2 - 2z + 0.75
        // LP coefficients d(1) = 2, d(2) = -0.75  (negated polynomial, reversed)
        const int NPOLES = 2;
        DFMatrixKokkos<double> d(NPOLES);
        d.host(1) =  2.0;
        d.host(2) = -0.75;
        d.update_device();

        std::printf("Before fixrts:\n");
        d.update_host();
        for (int j = 1; j <= NPOLES; j++)
            std::printf("  d(%d) = %12.6f\n", j, d.host(j));

        fixrts(d, NPOLES);

        d.update_host();
        std::printf("\nAfter fixrts:\n");
        for (int j = 1; j <= NPOLES; j++)
            std::printf("  d(%d) = %12.6f\n", j, d.host(j));

        // Verify roots via quadratic formula
        // Polynomial: z^2 - d(1)*z - d(2) = 0
        double disc = d.host(1) * d.host(1) + 4.0 * d.host(2);
        if (disc >= 0) {
            double r1 = (d.host(1) + std::sqrt(disc)) / 2.0;
            double r2 = (d.host(1) - std::sqrt(disc)) / 2.0;
            std::printf("\nRoots: r1 = %.6f (|r1| = %.6f)\n", r1, std::fabs(r1));
            std::printf("       r2 = %.6f (|r2| = %.6f)\n", r2, std::fabs(r2));
            bool stable = (std::fabs(r1) <= 1.0 + 1e-6) && (std::fabs(r2) <= 1.0 + 1e-6);
            std::printf("Stable: %s\n\n", stable ? "YES" : "NO");
        } else {
            double re = d.host(1) / 2.0;
            double im = std::sqrt(-disc) / 2.0;
            double mag = std::sqrt(re * re + im * im);
            std::printf("\nRoots: %.6f +/- %.6fi  (|r| = %.6f)\n", re, im, mag);
            bool stable = (mag <= 1.0 + 1e-6);
            std::printf("Stable: %s\n\n", stable ? "YES" : "NO");
        }
    }
    MATAR_FINALIZE();
    return 0;
}
