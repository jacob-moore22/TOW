#pragma once
#include <cmath>
#include <matar.h>

// State for ran2: single LCG with Bays-Durham shuffle table.
struct Ran2State {
    int iff;
    int ir[97];
    int iy;
    Ran2State() : iff(0), iy(0) {
        for (int i = 0; i < 97; i++) ir[i] = 0;
    }
};

// Long-period random number generator using a single LCG with
// Bays-Durham shuffle.
// Call with idum < 0 to initialize; subsequent calls leave idum alone.
inline double ran2(int& idum, Ran2State& state)
{
    constexpr int M = 714025, IA = 1366, IC = 150889;
    constexpr double RM = 1.4005112e-6;

    if (idum < 0 || state.iff == 0) {
        state.iff = 1;
        idum = (IC - idum) % M;
        for (int j = 0; j < 97; j++) {
            idum = (IA * idum + IC) % M;
            state.ir[j] = idum;
        }
        idum = (IA * idum + IC) % M;
        state.iy = idum;
    }
    int j = (97 * state.iy) / M;
    if (j > 96) j = 96;
    if (j < 0)  j = 0;
    state.iy = state.ir[j];
    double result = state.iy * RM;
    idum = (IA * idum + IC) % M;
    state.ir[j] = idum;
    return result;
}
