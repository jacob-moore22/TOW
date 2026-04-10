#include <cstdio>
#include <cmath>
#include <matar.h>
#include "quad3d.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        std::printf("3D integration (quad3d) tests\n\n");
        std::printf("%-40s %18s %18s %14s\n",
                    "Integral", "Computed", "Exact", "Error");

        // Constant y/z limits for the unit cube
        auto y1 = [](double)           { return 0.0; };
        auto y2 = [](double)           { return 1.0; };
        auto z1 = [](double, double)   { return 0.0; };
        auto z2 = [](double, double)   { return 1.0; };

        // Test 1: integral of 1 over [0,1]^3 = 1
        auto f_one = [](double, double, double) { return 1.0; };
        double r1 = quad3d(f_one, 0.0, 1.0, y1, y2, z1, z2);
        double e1 = 1.0;
        std::printf("%-40s %18.12f %18.12f %14.2e\n",
                    "1 over [0,1]^3", r1, e1, std::fabs(r1 - e1));

        // Test 2: integral of x*y*z over [0,1]^3 = 1/8
        auto f_xyz = [](double x, double y, double z) { return x * y * z; };
        double r2 = quad3d(f_xyz, 0.0, 1.0, y1, y2, z1, z2);
        double e2 = 1.0 / 8.0;
        std::printf("%-40s %18.12f %18.12f %14.2e\n",
                    "x*y*z over [0,1]^3", r2, e2, std::fabs(r2 - e2));

        // Test 3: integral of x^2+y^2+z^2 over [0,1]^3 = 1
        auto f_sum2 = [](double x, double y, double z) {
            return x * x + y * y + z * z;
        };
        double r3 = quad3d(f_sum2, 0.0, 1.0, y1, y2, z1, z2);
        double e3 = 1.0;
        std::printf("%-40s %18.12f %18.12f %14.2e\n",
                    "x^2+y^2+z^2 over [0,1]^3", r3, e3, std::fabs(r3 - e3));

        // Test 4: non-rectangular region — tetrahedron 0<=z<=y<=x<=1
        // integral of 1 = 1/6
        auto yt1 = [](double)          { return 0.0; };
        auto yt2 = [](double x)        { return x; };
        auto zt1 = [](double, double)  { return 0.0; };
        auto zt2 = [](double, double y) { return y; };

        double r4 = quad3d(f_one, 0.0, 1.0, yt1, yt2, zt1, zt2);
        double e4 = 1.0 / 6.0;
        std::printf("%-40s %18.12f %18.12f %14.2e\n",
                    "1 over tetrahedron 0<=z<=y<=x<=1", r4, e4, std::fabs(r4 - e4));

        bool pass = (std::fabs(r1 - e1) < 1e-6) &&
                    (std::fabs(r2 - e2) < 1e-6) &&
                    (std::fabs(r3 - e3) < 1e-6) &&
                    (std::fabs(r4 - e4) < 1e-6);
        std::printf("\nTest %s\n", pass ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
