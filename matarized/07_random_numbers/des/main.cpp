#include <cstdio>
#include <matar.h>
#include "des.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        std::printf("des — STUB (DES encryption not yet implemented)\n");
        unsigned long lword = 0, irword = 0;
        des(lword, irword, 1);
        std::printf("  Test : SKIPPED\n");
    }
    MATAR_FINALIZE();
    return 0;
}
