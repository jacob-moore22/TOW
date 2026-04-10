#include <cstdio>
#include <cmath>
#include <matar.h>
#include "anneal.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        std::printf("Simulated annealing for TSP (anneal)\n");
        std::printf("=====================================\n\n");

        // 10 cities arranged in a circle — optimal tour follows the circle
        constexpr int ncity = 10;
        double x[ncity], y[ncity];
        int iorder[ncity];

        for (int i = 0; i < ncity; i++) {
            double theta = 2.0 * M_PI * i / ncity;
            x[i] = std::cos(theta);
            y[i] = std::sin(theta);
            iorder[i] = i;
        }

        // Shuffle the initial tour
        anneal_detail::RNG rng(123);
        for (int i = ncity - 1; i > 0; i--) {
            int j = rng.irand(i + 1);
            int tmp = iorder[i];
            iorder[i] = iorder[j];
            iorder[j] = tmp;
        }

        // Compute initial path length
        double init_path = 0.0;
        for (int i = 0; i < ncity - 1; i++)
            init_path += anneal_detail::alen(x[iorder[i]], x[iorder[i + 1]],
                                             y[iorder[i]], y[iorder[i + 1]]);
        init_path += anneal_detail::alen(x[iorder[ncity - 1]], x[iorder[0]],
                                         y[iorder[ncity - 1]], y[iorder[0]]);

        std::printf("Initial path length: %.6f\n\n", init_path);
        std::printf("Annealing:\n");

        anneal(x, y, iorder, ncity);

        // Compute final path length
        double final_path = 0.0;
        for (int i = 0; i < ncity - 1; i++)
            final_path += anneal_detail::alen(x[iorder[i]], x[iorder[i + 1]],
                                              y[iorder[i]], y[iorder[i + 1]]);
        final_path += anneal_detail::alen(x[iorder[ncity - 1]], x[iorder[0]],
                                          y[iorder[ncity - 1]], y[iorder[0]]);

        // Optimal circle tour length: 10 * 2*sin(pi/10) ≈ 6.18034
        double optimal = ncity * 2.0 * std::sin(M_PI / ncity);

        std::printf("\nFinal path length:   %.6f\n", final_path);
        std::printf("Optimal (circle):    %.6f\n", optimal);
        std::printf("Ratio:               %.4f\n\n", final_path / optimal);

        std::printf("Tour: ");
        for (int i = 0; i < ncity; i++) std::printf("%d ", iorder[i]);
        std::printf("\n\n");

        bool pass = final_path < 1.05 * optimal;
        std::printf("Test %s\n", pass ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
