#pragma once
#include <cmath>
#include <cstdio>
#include <matar.h>

using namespace mtr;

// Simulated annealing for the Travelling Salesman Problem (TSP).
//
// Fortran COMMON blocks and external routines (RAN3, IRBIT1, TRNCST, METROP,
// TRNSPT, REVCST, REVERS) are replaced by self-contained inline helpers and
// a simple LCG random number generator.
//
// x, y:    city coordinates (length ncity, 0-based).
// iorder:  permutation of 0..ncity-1 giving the tour. Modified in place.
// ncity:   number of cities.

namespace anneal_detail {

constexpr int NMAX = 2048;

struct RNG {
    unsigned long state;
    explicit RNG(unsigned long seed = 12345UL) : state(seed) {}
    double operator()() {
        state = state * 6364136223846793005ULL + 1442695040888963407ULL;
        return static_cast<double>(state >> 33) / static_cast<double>(1ULL << 31);
    }
    int irand(int n) { return static_cast<int>((*this)() * n); }
    bool bit() { return irand(2) != 0; }
};

inline double alen(double x1, double x2, double y1, double y2)
{
    double dx = x2 - x1, dy = y2 - y1;
    return std::sqrt(dx * dx + dy * dy);
}

// Cost change for reversing the segment [n0..n1] (cyclic) in the tour
inline double revcst(const double* x, const double* y, const int* iorder,
                     int ncity, int n0, int n1)
{
    int prev0 = (n0 + ncity - 1) % ncity;
    int next1 = (n1 + 1) % ncity;
    return -alen(x[iorder[n0]], x[iorder[prev0]],
                 y[iorder[n0]], y[iorder[prev0]])
           -alen(x[iorder[n1]], x[iorder[next1]],
                 y[iorder[n1]], y[iorder[next1]])
           +alen(x[iorder[prev0]], x[iorder[n1]],
                 y[iorder[prev0]], y[iorder[n1]])
           +alen(x[iorder[next1]], x[iorder[n0]],
                 y[iorder[next1]], y[iorder[n0]]);
}

// Reverse the segment [n0..n1] (cyclic) in the tour
inline void revers(int* iorder, int ncity, int n0, int n1)
{
    int seg_len = ((n1 - n0 + ncity) % ncity) + 1;
    int swaps = seg_len / 2;
    int i = n0, j = n1;
    for (int k = 0; k < swaps; k++) {
        int tmp = iorder[i];
        iorder[i] = iorder[j];
        iorder[j] = tmp;
        i = (i + 1) % ncity;
        j = (j + ncity - 1) % ncity;
    }
}

// Cost change for transporting segment [n0..n1] to follow position n2
inline double trncst(const double* x, const double* y, const int* iorder,
                     int ncity, int n0, int n1, int n2)
{
    int prev0 = (n0 + ncity - 1) % ncity;
    int next1 = (n1 + 1) % ncity;
    int next2 = (n2 + 1) % ncity;

    return -alen(x[iorder[prev0]], x[iorder[n0]],
                 y[iorder[prev0]], y[iorder[n0]])
           -alen(x[iorder[n1]],   x[iorder[next1]],
                 y[iorder[n1]],   y[iorder[next1]])
           -alen(x[iorder[n2]],   x[iorder[next2]],
                 y[iorder[n2]],   y[iorder[next2]])
           +alen(x[iorder[prev0]], x[iorder[next1]],
                 y[iorder[prev0]], y[iorder[next1]])
           +alen(x[iorder[n2]],   x[iorder[n0]],
                 y[iorder[n2]],   y[iorder[n0]])
           +alen(x[iorder[n1]],   x[iorder[next2]],
                 y[iorder[n1]],   y[iorder[next2]]);
}

// Transport segment [n0..n1] (cyclic) to follow position n2 in the tour.
// n2 must lie outside the segment, in the arc from n1+1 to n0-1 (cyclic).
inline void trnspt(int* iorder, int ncity, int n0, int n1, int n2)
{
    int jorder[NMAX];

    int seg_len = ((n1 - n0 + ncity) % ncity) + 1;
    int seg[NMAX];
    for (int i = 0; i < seg_len; i++)
        seg[i] = iorder[(n0 + i) % ncity];

    // Rebuild: rest-before-n2 | segment | rest-after-n2
    int pos = 0;
    int j = (n1 + 1) % ncity;
    int stop = (n2 + 1) % ncity;

    // Copy the part of the rest from n1+1 through n2
    while (j != stop) {
        jorder[pos++] = iorder[j];
        j = (j + 1) % ncity;
    }

    // Insert the transported segment
    for (int i = 0; i < seg_len; i++)
        jorder[pos++] = seg[i];

    // Copy the remaining rest from n2+1 through n0-1
    while (pos < ncity) {
        jorder[pos++] = iorder[j];
        j = (j + 1) % ncity;
    }

    for (int i = 0; i < ncity; i++)
        iorder[i] = jorder[i];
}

inline bool metrop(double de, double t, RNG& rng)
{
    return de < 0.0 || rng() < std::exp(-de / t);
}

} // namespace anneal_detail

inline void anneal(const double* x, const double* y, int* iorder, int ncity)
{
    using namespace anneal_detail;

    int nover  = 100 * ncity;
    int nlimit = 10 * ncity;
    double tfactr = 0.9;

    double path = 0.0;
    for (int i = 0; i < ncity - 1; i++)
        path += alen(x[iorder[i]], x[iorder[i + 1]],
                     y[iorder[i]], y[iorder[i + 1]]);
    path += alen(x[iorder[ncity - 1]], x[iorder[0]],
                 y[iorder[ncity - 1]], y[iorder[0]]);

    double t = 0.5;
    RNG rng(42);

    for (int j = 0; j < 100; j++) {
        int nsucc = 0;
        for (int k = 0; k < nover; k++) {
            int n0 = rng.irand(ncity);
            int n1 = rng.irand(ncity - 1);
            if (n1 >= n0) n1++;

            int nn = 1 + ((n0 - n1 + ncity - 1) % ncity);
            if (nn < 3) continue;

            if (rng.bit()) {
                int n2 = (n1 + 1 + rng.irand(nn - 2)) % ncity;
                double de = trncst(x, y, iorder, ncity, n0, n1, n2);
                if (metrop(de, t, rng)) {
                    nsucc++;
                    path += de;
                    trnspt(iorder, ncity, n0, n1, n2);
                }
            } else {
                double de = revcst(x, y, iorder, ncity, n0, n1);
                if (metrop(de, t, rng)) {
                    nsucc++;
                    path += de;
                    revers(iorder, ncity, n0, n1);
                }
            }
            if (nsucc >= nlimit) break;
        }

        std::printf("  T = %10.6f   Path Length = %10.4f   Moves: %d\n",
                    t, path, nsucc);
        t *= tfactr;
        if (nsucc == 0) return;
    }
}
