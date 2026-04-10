// Jacobi eigenvalue method for real symmetric matrices (Numerical Recipes JACOBI).
// Computes all eigenvalues and eigenvectors of a real symmetric matrix a[1..n][1..n].
// On output, d[1..n] contains eigenvalues, v[1..n][1..n] contains eigenvectors
// (columns), and nrot is the number of Jacobi rotations performed.
// The upper triangle of a is destroyed.

#pragma once
#include <cmath>
#include <cstdio>
#include <matar.h>

using namespace mtr;

inline void jacobi(DFMatrixKokkos<double>& a, int n,
                   DFMatrixKokkos<double>& d, DFMatrixKokkos<double>& v,
                   int& nrot)
{
    constexpr int NMAX = 100;
    double b[NMAX + 1], z[NMAX + 1];

    for (int ip = 1; ip <= n; ip++) {
        for (int iq = 1; iq <= n; iq++) {
            v.host(ip, iq) = 0.0;
        }
        v.host(ip, ip) = 1.0;
    }

    for (int ip = 1; ip <= n; ip++) {
        b[ip] = a.host(ip, ip);
        d.host(ip) = b[ip];
        z[ip] = 0.0;
    }

    nrot = 0;

    for (int i = 1; i <= 50; i++) {
        double sm = 0.0;
        for (int ip = 1; ip <= n - 1; ip++) {
            for (int iq = ip + 1; iq <= n; iq++) {
                sm += std::fabs(a.host(ip, iq));
            }
        }
        if (sm == 0.0) return;

        double tresh;
        if (i < 4) {
            tresh = 0.2 * sm / (n * n);
        } else {
            tresh = 0.0;
        }

        for (int ip = 1; ip <= n - 1; ip++) {
            for (int iq = ip + 1; iq <= n; iq++) {
                double g = 100.0 * std::fabs(a.host(ip, iq));
                if (i > 4 &&
                    std::fabs(d.host(ip)) + g == std::fabs(d.host(ip)) &&
                    std::fabs(d.host(iq)) + g == std::fabs(d.host(iq))) {
                    a.host(ip, iq) = 0.0;
                } else if (std::fabs(a.host(ip, iq)) > tresh) {
                    double h = d.host(iq) - d.host(ip);
                    double t;
                    if (std::fabs(h) + g == std::fabs(h)) {
                        t = a.host(ip, iq) / h;
                    } else {
                        double theta = 0.5 * h / a.host(ip, iq);
                        t = 1.0 / (std::fabs(theta) + std::sqrt(1.0 + theta * theta));
                        if (theta < 0.0) t = -t;
                    }
                    double c = 1.0 / std::sqrt(1.0 + t * t);
                    double s = t * c;
                    double tau = s / (1.0 + c);
                    h = t * a.host(ip, iq);
                    z[ip] -= h;
                    z[iq] += h;
                    d.host(ip) -= h;
                    d.host(iq) += h;
                    a.host(ip, iq) = 0.0;

                    for (int j = 1; j <= ip - 1; j++) {
                        g = a.host(j, ip);
                        h = a.host(j, iq);
                        a.host(j, ip) = g - s * (h + g * tau);
                        a.host(j, iq) = h + s * (g - h * tau);
                    }
                    for (int j = ip + 1; j <= iq - 1; j++) {
                        g = a.host(ip, j);
                        h = a.host(j, iq);
                        a.host(ip, j) = g - s * (h + g * tau);
                        a.host(j, iq) = h + s * (g - h * tau);
                    }
                    for (int j = iq + 1; j <= n; j++) {
                        g = a.host(ip, j);
                        h = a.host(iq, j);
                        a.host(ip, j) = g - s * (h + g * tau);
                        a.host(iq, j) = h + s * (g - h * tau);
                    }
                    for (int j = 1; j <= n; j++) {
                        g = v.host(j, ip);
                        h = v.host(j, iq);
                        v.host(j, ip) = g - s * (h + g * tau);
                        v.host(j, iq) = h + s * (g - h * tau);
                    }
                    nrot++;
                }
            }
        }

        for (int ip = 1; ip <= n; ip++) {
            b[ip] += z[ip];
            d.host(ip) = b[ip];
            z[ip] = 0.0;
        }
    }

    std::fprintf(stderr, "jacobi: 50 iterations should never happen\n");
}
