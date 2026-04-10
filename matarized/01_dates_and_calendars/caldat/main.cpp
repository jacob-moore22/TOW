#include <cstdio>
#include <cmath>
#include <matar.h>
#include "caldat.hpp"
#include "../julday/julday.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        std::printf("Calendar date from Julian Day tests:\n");

        // Round-trip test: julday(1,1,2000) -> caldat -> should give 1/1/2000
        int jd = julday(1, 1, 2000);
        int mm, id, iyyy;
        caldat(jd, mm, id, iyyy);

        std::printf("julday(1,1,2000) = %d\n", jd);
        std::printf("caldat(%d) = %d/%d/%d\n", jd, mm, id, iyyy);

        bool pass1 = (mm == 1 && id == 1 && iyyy == 2000);

        // More round-trip tests
        struct TestCase { int m; int d; int y; };
        TestCase cases[] = {
            { 7,  4, 1776},
            {12, 25, 1999},
            { 3, 15,   44},
            { 6, 15, 2026},
        };
        int ncases = 4;
        bool all_pass = pass1;

        std::printf("\nRound-trip tests (julday -> caldat):\n");
        std::printf("%6s %6s %6s %6s %6s %6s %8s\n",
                    "M_in", "D_in", "Y_in", "M_out", "D_out", "Y_out", "match");

        for (int t = 0; t < ncases; t++) {
            int jd_t = julday(cases[t].m, cases[t].d, cases[t].y);
            int m_out, d_out, y_out;
            caldat(jd_t, m_out, d_out, y_out);
            bool ok = (m_out == cases[t].m && d_out == cases[t].d);
            // Year comparison needs care for BC dates
            if (cases[t].y > 0)
                ok = ok && (y_out == cases[t].y);
            if (!ok) all_pass = false;
            std::printf("%6d %6d %6d %6d %6d %6d %8s\n",
                        cases[t].m, cases[t].d, cases[t].y,
                        m_out, d_out, y_out, ok ? "OK" : "FAIL");
        }

        std::printf("\nTest %s\n", all_pass ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
