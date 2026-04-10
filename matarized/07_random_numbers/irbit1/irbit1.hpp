#pragma once
#include <matar.h>

// Returns a random bit (0 or 1) using a linear-feedback shift register
// with primitive polynomial taps at bits 1, 2, 5, 18.
// iseed is modified in place (serves as the LFSR state).
KOKKOS_INLINE_FUNCTION
int irbit1(int& iseed)
{
    constexpr int IB1 = 1, IB2 = 2, IB5 = 16, IB18 = 131072;
    bool newbit = (iseed & IB18) != 0;
    if ((iseed & IB5) != 0) newbit = !newbit;
    if ((iseed & IB2) != 0) newbit = !newbit;
    if ((iseed & IB1) != 0) newbit = !newbit;
    int result = 0;
    iseed = (iseed << 1) & (~IB1);
    if (newbit) {
        result = 1;
        iseed = iseed | IB1;
    }
    return result;
}
