#include <cstdio>
#include <matar.h>
#include "ran4.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        std::printf("ran4 — STUB (depends on des, not yet implemented)\n");
        int idum = -1;
        Ran4State state;
        double x = ran4(idum, state);
        std::printf("  ran4(-1) = %.6f  (stub returns 0)\n", x);
        std::printf("  Test : SKIPPED\n");
    }
    MATAR_FINALIZE();
    return 0;
}
