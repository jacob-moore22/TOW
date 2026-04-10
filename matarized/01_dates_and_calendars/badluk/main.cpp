#include <cstdio>
#include <matar.h>

using namespace mtr;

#include "../julday/julday.hpp"
#include "../caldat/caldat.hpp"
#include "../flmoon/flmoon.hpp"
#include "badluk.hpp"

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        double timzon = -5.0 / 24.0;

        // Known result: Sep 13, 2019 is a Friday the 13th with a full moon
        std::printf("=== Friday the 13th full moons (1900-2100) ===\n\n");
        badluk(1900, 2100, timzon);

        // Verify a known Friday the 13th
        int jd = julday(9, 13, 2019);
        int idwk = (jd + 1) % 7;
        std::printf("\nVerification: Sep 13, 2019 is day-of-week %d (5=Friday)\n", idwk);
        std::printf("Test %s\n", idwk == 5 ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
