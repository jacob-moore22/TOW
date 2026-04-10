#pragma once
#include <cmath>
#include <matar.h>

// State for ran0: Bays-Durham shuffle of a Park-Miller LCG.
struct Ran0State {
    int iff;
    int iseed;
    double v[97];
    double y;
    Ran0State() : iff(0), iseed(0), y(0.0) {
        for (int i = 0; i < 97; i++) v[i] = 0.0;
    }
};

// Park-Miller minimal standard LCG, replacing Fortran system RAN()
inline double park_miller_ran(int& iseed)
{
    constexpr int IA = 16807;
    constexpr int IM = 2147483647;
    constexpr double AM = 1.0 / static_cast<double>(IM);
    int k = iseed / 127773;
    iseed = IA * (iseed - k * 127773) - 2836 * k;
    if (iseed < 0) iseed += IM;
    return AM * iseed;
}

// Minimal random number generator with Bays-Durham shuffle.
// Call with idum < 0 to initialize; subsequent calls use idum = 1.
inline double ran0(int& idum, Ran0State& state)
{
    if (idum < 0 || state.iff == 0) {
        state.iff = 1;
        state.iseed = (idum < 0) ? -idum : idum;
        if (state.iseed == 0) state.iseed = 1;
        idum = 1;
        for (int j = 0; j < 97; j++)
            park_miller_ran(state.iseed);
        for (int j = 0; j < 97; j++)
            state.v[j] = park_miller_ran(state.iseed);
        state.y = park_miller_ran(state.iseed);
    }
    int j = static_cast<int>(97.0 * state.y);
    if (j > 96) j = 96;
    if (j < 0)  j = 0;
    state.y = state.v[j];
    state.v[j] = park_miller_ran(state.iseed);
    return state.y;
}
