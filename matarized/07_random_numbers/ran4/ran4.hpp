#pragma once
#include <cmath>
#include <matar.h>
#include "des.hpp"

// STUB: ran4 depends on the DES cipher for its hash.
// Returns 0.0 until des is implemented.
struct Ran4State {
    int iff;
    unsigned long jflone, jflmsk;
    Ran4State() : iff(0), jflone(0), jflmsk(0) {}
};

inline double ran4(int& idum, Ran4State& state)
{
    (void)state;
    (void)idum;
    return 0.0;
}
