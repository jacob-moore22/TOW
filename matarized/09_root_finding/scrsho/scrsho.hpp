#pragma once
#include <cmath>
#include <cstdio>
#include <matar.h>

using namespace mtr;

// Simple text-mode screen plot of a function over [x1,x2].
// Produces a character grid showing the function shape and zero crossings.
template<typename Func>
inline void scrsho(Func func, double x1, double x2)
{
    constexpr int ISCR = 60;
    constexpr int JSCR = 21;

    char scr[ISCR][JSCR];
    double y[ISCR];

    // Initialize border
    for (int j = 0; j < JSCR; j++) {
        scr[0][j]        = '|';
        scr[ISCR - 1][j] = '|';
    }
    for (int i = 1; i < ISCR - 1; i++) {
        scr[i][0]        = '-';
        scr[i][JSCR - 1] = '-';
        for (int j = 1; j < JSCR - 1; j++)
            scr[i][j] = ' ';
    }

    // Evaluate function and find range
    double dx = (x2 - x1) / (ISCR - 1);
    double xv = x1;
    double ybig = 0.0, ysml = 0.0;

    for (int i = 0; i < ISCR; i++) {
        y[i] = func(xv);
        if (y[i] < ysml) ysml = y[i];
        if (y[i] > ybig) ybig = y[i];
        xv += dx;
    }

    if (ybig == ysml) ybig = ysml + 1.0;

    double dyj = (JSCR - 1) / (ybig - ysml);
    int jz = static_cast<int>(1 - ysml * dyj);
    if (jz < 0)        jz = 0;
    if (jz >= JSCR)    jz = JSCR - 1;

    for (int i = 0; i < ISCR; i++) {
        scr[i][jz] = '-';
        int j = static_cast<int>((y[i] - ysml) * dyj);
        if (j < 0)     j = 0;
        if (j >= JSCR) j = JSCR - 1;
        scr[i][j] = 'x';
    }

    // Print (top row = ybig, bottom row = ysml)
    std::printf(" %10.3e ", ybig);
    for (int i = 0; i < ISCR; i++) std::printf("%c", scr[i][JSCR - 1]);
    std::printf("\n");

    for (int j = JSCR - 2; j >= 1; j--) {
        std::printf("            ");
        for (int i = 0; i < ISCR; i++) std::printf("%c", scr[i][j]);
        std::printf("\n");
    }

    std::printf(" %10.3e ", ysml);
    for (int i = 0; i < ISCR; i++) std::printf("%c", scr[i][0]);
    std::printf("\n");

    std::printf("            %-10.3e%*s%10.3e\n", x1, ISCR - 20, "", x2);
}
