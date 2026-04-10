#pragma once
#include <matar.h>

// Returns a random bit (0 or 1) using a linear-feedback shift register
// with XOR-based feedback at bits 1, 2, 5, 18.
// iseed is modified in place (serves as the LFSR state).
KOKKOS_INLINE_FUNCTION
int irbit2(int& iseed)
{
    constexpr int IB1 = 1, IB2 = 2, IB5 = 16, IB18 = 131072;
    constexpr int MASK = IB1 + IB2 + IB5;
    if ((iseed & IB18) != 0) {
        iseed = ((iseed ^ MASK) << 1) | IB1;
        return 1;
    } else {
        iseed = (iseed << 1) & (~IB1);
        return 0;
    }
}
