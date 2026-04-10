#pragma once
#include <cmath>
#include <matar.h>

// State for ran1: three-LCG combination with Bays-Durham shuffle table.
struct Ran1State {
    int ix1, ix2, ix3;
    double r[97];
    int iff;
    Ran1State() : ix1(0), ix2(0), ix3(0), iff(0) {
        for (int i = 0; i < 97; i++) r[i] = 0.0;
    }
};

// Portable random number generator combining three linear congruential
// generators with a Bays-Durham shuffle.
// Call with idum < 0 to initialize; subsequent calls use idum = 1.
inline double ran1(int& idum, Ran1State& state)
{
    constexpr int M1 = 259200, IA1 = 7141, IC1 = 54773;
    constexpr double RM1 = 3.8580247e-6;
    constexpr int M2 = 134456, IA2 = 8121, IC2 = 28411;
    constexpr double RM2 = 7.4373773e-6;
    constexpr int M3 = 243000, IA3 = 4561, IC3 = 51349;

    if (idum < 0 || state.iff == 0) {
        state.iff = 1;
        state.ix1 = (IC1 - idum) % M1;
        state.ix1 = (IA1 * state.ix1 + IC1) % M1;
        state.ix2 = state.ix1 % M2;
        state.ix1 = (IA1 * state.ix1 + IC1) % M1;
        state.ix3 = state.ix1 % M3;
        for (int j = 0; j < 97; j++) {
            state.ix1 = (IA1 * state.ix1 + IC1) % M1;
            state.ix2 = (IA2 * state.ix2 + IC2) % M2;
            state.r[j] = (static_cast<double>(state.ix1) +
                          static_cast<double>(state.ix2) * RM2) * RM1;
        }
        idum = 1;
    }
    state.ix1 = (IA1 * state.ix1 + IC1) % M1;
    state.ix2 = (IA2 * state.ix2 + IC2) % M2;
    state.ix3 = (IA3 * state.ix3 + IC3) % M3;
    int j = (97 * state.ix3) / M3;
    if (j > 96) j = 96;
    if (j < 0)  j = 0;
    double result = state.r[j];
    state.r[j] = (static_cast<double>(state.ix1) +
                  static_cast<double>(state.ix2) * RM2) * RM1;
    return result;
}
