#include <cstdio>
#include <matar.h>
#include "desks.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        std::printf("desks — STUB (DES key schedule not yet implemented)\n");
        unsigned long key[2] = {0, 0};
        desks(key);
        std::printf("  Test : SKIPPED\n");
    }
    MATAR_FINALIZE();
    return 0;
}
