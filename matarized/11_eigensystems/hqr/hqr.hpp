// QR algorithm for eigenvalues of an upper Hessenberg matrix
// (Numerical Recipes HQR).
// On input, a[1..n][1..n] is an upper Hessenberg matrix (output of elmhes).
// On output, wr[1..n] and wi[1..n] contain the real and imaginary parts
// of the eigenvalues. The matrix a is destroyed.

#pragma once
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <matar.h>

using namespace mtr;

inline void hqr(DFMatrixKokkos<double>& a, int n,
                DFMatrixKokkos<double>& wr, DFMatrixKokkos<double>& wi)
{
    double p, q, r, s, t, w, x, y, z;
    double anorm, u, v;
    int i, j, k, l, m, nn, its, mmin;

    anorm = std::fabs(a.host(1, 1));
    for (i = 2; i <= n; i++) {
        for (j = i - 1; j <= n; j++) {
            anorm += std::fabs(a.host(i, j));
        }
    }

    nn = n;
    t = 0.0;

    while (nn >= 1) {
        its = 0;

        while (true) {
            for (l = nn; l >= 2; l--) {
                s = std::fabs(a.host(l - 1, l - 1)) + std::fabs(a.host(l, l));
                if (s == 0.0) s = anorm;
                if (std::fabs(a.host(l, l - 1)) + s == s) break;
            }

            x = a.host(nn, nn);

            if (l == nn) {
                wr.host(nn) = x + t;
                wi.host(nn) = 0.0;
                nn--;
                break;
            }

            y = a.host(nn - 1, nn - 1);
            w = a.host(nn, nn - 1) * a.host(nn - 1, nn);

            if (l == nn - 1) {
                p = 0.5 * (y - x);
                q = p * p + w;
                z = std::sqrt(std::fabs(q));
                x += t;
                if (q >= 0.0) {
                    z = p + std::copysign(z, p);
                    wr.host(nn) = x + z;
                    wr.host(nn - 1) = wr.host(nn);
                    if (z != 0.0) wr.host(nn) = x - w / z;
                    wi.host(nn) = 0.0;
                    wi.host(nn - 1) = 0.0;
                } else {
                    wr.host(nn) = x + p;
                    wr.host(nn - 1) = wr.host(nn);
                    wi.host(nn) = z;
                    wi.host(nn - 1) = -z;
                }
                nn -= 2;
                break;
            }

            if (its == 30) {
                std::fprintf(stderr, "hqr: too many iterations\n");
                return;
            }

            if (its == 10 || its == 20) {
                t += x;
                for (i = 1; i <= nn; i++) {
                    a.host(i, i) -= x;
                }
                s = std::fabs(a.host(nn, nn - 1)) + std::fabs(a.host(nn - 1, nn - 2));
                x = 0.75 * s;
                y = x;
                w = -0.4375 * s * s;
            }

            its++;

            for (m = nn - 2; m >= l; m--) {
                z = a.host(m, m);
                r = x - z;
                s = y - z;
                p = (r * s - w) / a.host(m + 1, m) + a.host(m, m + 1);
                q = a.host(m + 1, m + 1) - z - r - s;
                r = a.host(m + 2, m + 1);
                s = std::fabs(p) + std::fabs(q) + std::fabs(r);
                p /= s;
                q /= s;
                r /= s;
                if (m == l) break;
                u = std::fabs(a.host(m, m - 1)) * (std::fabs(q) + std::fabs(r));
                v = std::fabs(p) * (std::fabs(a.host(m - 1, m - 1)) +
                    std::fabs(z) + std::fabs(a.host(m + 1, m + 1)));
                if (u + v == v) break;
            }

            for (i = m + 2; i <= nn; i++) {
                a.host(i, i - 2) = 0.0;
                if (i != m + 2) a.host(i, i - 3) = 0.0;
            }

            for (k = m; k <= nn - 1; k++) {
                if (k != m) {
                    p = a.host(k, k - 1);
                    q = a.host(k + 1, k - 1);
                    r = 0.0;
                    if (k != nn - 1) r = a.host(k + 2, k - 1);
                    x = std::fabs(p) + std::fabs(q) + std::fabs(r);
                    if (x != 0.0) {
                        p /= x;
                        q /= x;
                        r /= x;
                    }
                }
                s = std::copysign(std::sqrt(p * p + q * q + r * r), p);
                if (s != 0.0) {
                    if (k == m) {
                        if (l != m) a.host(k, k - 1) = -a.host(k, k - 1);
                    } else {
                        a.host(k, k - 1) = -s * x;
                    }
                    p += s;
                    x = p / s;
                    y = q / s;
                    z = r / s;
                    q /= p;
                    r /= p;
                    for (j = k; j <= nn; j++) {
                        p = a.host(k, j) + q * a.host(k + 1, j);
                        if (k != nn - 1) {
                            p += r * a.host(k + 2, j);
                            a.host(k + 2, j) -= p * z;
                        }
                        a.host(k + 1, j) -= p * y;
                        a.host(k, j) -= p * x;
                    }
                    mmin = std::min(nn, k + 3);
                    for (i = l; i <= mmin; i++) {
                        p = x * a.host(i, k) + y * a.host(i, k + 1);
                        if (k != nn - 1) {
                            p += z * a.host(i, k + 2);
                            a.host(i, k + 2) -= p * r;
                        }
                        a.host(i, k + 1) -= p * q;
                        a.host(i, k) -= p;
                    }
                }
            }
        }
    }
}
