#pragma once
#include <cmath>
#include <matar.h>
#include "ran1.hpp"

// Exponential deviate: returns a random deviate drawn from an
// exponential distribution with unit mean, using ran1 as the
// uniform source.
inline double expdev(int& idum, Ran1State& state)
{
    return -log(ran1(idum, state));
}
