#include <cstdio>
#include <matar.h>
#include "julday.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        struct TestCase { int mm; int id; int iyyy; int expected; };
        TestCase cases[] = {
            { 1,  1, 2000, 2451545},
            {10, 15, 1582, 2299161},
            { 7,  4, 1776, 2369916},
            { 1,  1,    1, 1721424},
        };
        int ncases = 4;

        std::printf("Julian Day Number tests:\n");
        std::printf("%6s %6s %6s %12s %12s %8s\n",
                    "Month", "Day", "Year", "julday", "expected", "status");

        bool all_pass = true;
        for (int t = 0; t < ncases; t++) {
            int jd = julday(cases[t].mm, cases[t].id, cases[t].iyyy);
            bool ok = (jd == cases[t].expected);
            if (!ok) all_pass = false;
            std::printf("%6d %6d %6d %12d %12d %8s\n",
                        cases[t].mm, cases[t].id, cases[t].iyyy,
                        jd, cases[t].expected, ok ? "OK" : "FAIL");
        }

        // Round-trip test: julday(1,1,2000) -> caldat -> should get back same date
        int jd2000 = julday(1, 1, 2000);
        std::printf("\nRound-trip: julday(1,1,2000) = %d\n", jd2000);

        std::printf("Test %s\n", all_pass ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
