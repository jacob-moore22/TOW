#pragma once
#include <cmath>
#include <matar.h>

// State for ran3: Knuth's subtractive method.
struct Ran3State {
    int iff;
    int inext, inextp;
    int ma[55];
    Ran3State() : iff(0), inext(0), inextp(0) {
        for (int i = 0; i < 55; i++) ma[i] = 0;
    }
};

// Knuth's subtractive random number generator.
// Call with idum < 0 to initialize; subsequent calls use idum = 1.
inline double ran3(int& idum, Ran3State& state)
{
    constexpr int MBIG  = 1000000000;
    constexpr int MSEED = 161803398;
    constexpr int MZ    = 0;
    constexpr double FAC = 1.0e-9;

    if (idum < 0 || state.iff == 0) {
        state.iff = 1;
        int mj = MSEED - (idum < 0 ? -idum : idum);
        mj = mj % MBIG;
        if (mj < 0) mj += MBIG;
        state.ma[54] = mj;
        int mk = 1;
        for (int i = 1; i <= 54; i++) {
            int ii = (21 * i) % 55;      // Fortran 1-based index (1..54)
            state.ma[ii - 1] = mk;       // convert to 0-based
            mk = mj - mk;
            if (mk < MZ) mk += MBIG;
            mj = state.ma[ii - 1];
        }
        for (int k = 0; k < 4; k++) {
            for (int i = 0; i < 55; i++) {
                state.ma[i] -= state.ma[(i + 31) % 55];
                if (state.ma[i] < MZ) state.ma[i] += MBIG;
            }
        }
        state.inext  = -1;
        state.inextp = 30;
        idum = 1;
    }
    state.inext++;
    if (state.inext == 55) state.inext = 0;
    state.inextp++;
    if (state.inextp == 55) state.inextp = 0;
    int mj = state.ma[state.inext] - state.ma[state.inextp];
    if (mj < MZ) mj += MBIG;
    state.ma[state.inext] = mj;
    return mj * FAC;
}
