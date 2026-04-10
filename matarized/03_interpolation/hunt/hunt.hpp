// Search with initial guess (Numerical Recipes HUNT).
// Given array xx[1..n] and a value x, and an initial guess jlo,
// returns updated jlo such that x is between xx(jlo) and xx(jlo+1).
// jlo=0 or jlo=n indicates out of range.

#pragma once
#include <matar.h>

using namespace mtr;

inline void hunt(DFMatrixKokkos<double>& xx, int n, double x, int& jlo)
{
    bool ascnd = (xx.host(n) > xx.host(1));

    if (jlo <= 0 || jlo > n) {
        jlo = 0;
        int jhi = n + 1;
        while (jhi - jlo != 1) {
            int jm = (jhi + jlo) / 2;
            if ((x > xx.host(jm)) == ascnd) {
                jlo = jm;
            } else {
                jhi = jm;
            }
        }
        return;
    }

    int inc = 1;
    int jhi;
    if ((x >= xx.host(jlo)) == ascnd) {
        jhi = jlo + inc;
        while (true) {
            if (jhi > n) {
                jhi = n + 1;
                break;
            } else if ((x >= xx.host(jhi)) == ascnd) {
                jlo = jhi;
                inc += inc;
                jhi = jlo + inc;
            } else {
                break;
            }
        }
    } else {
        jhi = jlo;
        jlo = jhi - inc;
        while (true) {
            if (jlo < 1) {
                jlo = 0;
                break;
            } else if ((x < xx.host(jlo)) == ascnd) {
                jhi = jlo;
                inc += inc;
                jlo = jhi - inc;
            } else {
                break;
            }
        }
    }

    while (jhi - jlo != 1) {
        int jm = (jhi + jlo) / 2;
        if ((x > xx.host(jm)) == ascnd) {
            jlo = jm;
        } else {
            jhi = jm;
        }
    }
}
