#pragma once
#include <cmath>
#include <cstdio>
#include <matar.h>

using namespace mtr;

// Find all root-bracketing sub-intervals of func in [x1,x2].
// Divides [x1,x2] into n equal sub-intervals and returns pairs where sign
// changes occur. On entry nb_max is the max number of brackets to store;
// on exit nb is the number found. xb1[i],xb2[i] give the i-th bracket.
template<typename Func>
inline void zbrak(Func func, double x1, double x2, int n,
                  double* xb1, double* xb2, int nb_max, int& nb)
{
    nb = 0;
    double dx = (x2 - x1) / n;
    double x  = x1;
    double fp = func(x);

    for (int i = 0; i < n; i++) {
        x += dx;
        double fc = func(x);
        if (fc * fp < 0.0) {
            xb1[nb] = x - dx;
            xb2[nb] = x;
            nb++;
            if (nb == nb_max) return;
        }
        fp = fc;
    }
}
