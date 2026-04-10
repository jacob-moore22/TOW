#pragma once
#include <cmath>
#include <stdexcept>
#include <matar.h>

using namespace mtr;

// Quicksort. Sorts arr(1..n) into ascending order using the quicksort
// algorithm. Uses a randomized pivot selection and falls back to insertion sort
// for small subarrays.
inline void qcksrt(int n, DFMatrixKokkos<double>& arr)
{
    constexpr int    M      = 7;
    constexpr int    NSTACK = 50;
    constexpr double FM     = 7875.0;
    constexpr double FA     = 211.0;
    constexpr double FC     = 1663.0;
    constexpr double FMI    = 1.2698413e-4;

    int istack[NSTACK + 1];
    int jstack = 0;
    int l  = 1;
    int ir = n;
    double fx = 0.0;

    while (true) {
        if (ir - l < M) {
            for (int j = l + 1; j <= ir; j++) {
                double a = arr.host(j);
                int i;
                for (i = j - 1; i >= 1; i--) {
                    if (arr.host(i) <= a) break;
                    arr.host(i + 1) = arr.host(i);
                }
                arr.host(i + 1) = a;
            }
            if (jstack == 0) return;
            ir = istack[jstack];
            l  = istack[jstack - 1];
            jstack -= 2;
        } else {
            int i = l;
            int j = ir;
            fx = std::fmod(fx * FA + FC, FM);
            int iq = l + static_cast<int>((ir - l + 1) * (fx * FMI));
            double a = arr.host(iq);
            arr.host(iq) = arr.host(l);

            while (true) {
                while (j > 0 && a < arr.host(j)) {
                    j--;
                }
                if (j <= i) {
                    arr.host(i) = a;
                    break;
                }
                arr.host(i) = arr.host(j);
                i++;
                while (i <= n && a > arr.host(i)) {
                    i++;
                }
                if (j <= i) {
                    arr.host(j) = a;
                    i = j;
                    break;
                }
                arr.host(j) = arr.host(i);
                j--;
            }

            jstack += 2;
            if (jstack > NSTACK) {
                throw std::runtime_error("NSTACK must be made larger in qcksrt.");
            }
            if (ir - i >= i - l) {
                istack[jstack]     = ir;
                istack[jstack - 1] = i + 1;
                ir = i - 1;
            } else {
                istack[jstack]     = i - 1;
                istack[jstack - 1] = l;
                l = i + 1;
            }
        }
    }
}
