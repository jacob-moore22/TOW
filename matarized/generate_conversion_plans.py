#!/usr/bin/env python3
"""
Generate CONVERSION_PLAN.md files for all examples under matarized/.

Parses Fortran 77 sources and auto-generated Makefiles from the companion
fortran/ directory to produce structured conversion plans for porting each
Numerical Recipes example to performance-portable C++ using MATAR.

Usage:
    python3 generate_conversion_plans.py            # generate all plans
    python3 generate_conversion_plans.py --dry-run   # print stats, don't write
    python3 generate_conversion_plans.py <example>   # generate for one example only
"""

import os
import re
import sys
from pathlib import Path
from collections import defaultdict

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR  = Path(__file__).resolve().parent
ROOT_DIR    = SCRIPT_DIR.parent
FORTRAN_DIR = ROOT_DIR / "fortran"

# ---------------------------------------------------------------------------
# Chapter metadata
# ---------------------------------------------------------------------------
CHAPTER_TITLES = {
    "01": "Dates and Calendars",
    "02": "Linear Algebra",
    "03": "Interpolation",
    "04": "Integration",
    "05": "Evaluation of Functions",
    "06": "Special Functions",
    "07": "Random Numbers",
    "08": "Sorting",
    "09": "Root Finding",
    "10": "Optimization",
    "11": "Eigensystems",
    "12": "Fourier Transform",
    "13": "Spectral Analysis",
    "14": "Statistics",
    "15": "Curve Fitting",
    "16": "ODE Integration",
    "17": "Boundary Value Problems",
    "19": "Partial Differential Equations",
}

CHAPTER_TAGS = {
    "01": ["utility", "calendar"],
    "02": ["linear-algebra", "matrix"],
    "03": ["interpolation", "approximation"],
    "04": ["integration", "quadrature"],
    "05": ["function-evaluation", "polynomial", "chebyshev"],
    "06": ["special-function", "mathematical-function"],
    "07": ["random-number", "stochastic"],
    "08": ["sorting", "ordering"],
    "09": ["root-finding", "nonlinear-equation"],
    "10": ["optimization", "minimization"],
    "11": ["eigensystem", "eigenvalue", "matrix"],
    "12": ["fourier-transform", "fft", "signal-processing"],
    "13": ["spectral-analysis", "signal-processing"],
    "14": ["statistics", "hypothesis-testing"],
    "15": ["curve-fitting", "regression", "least-squares"],
    "16": ["ode", "differential-equation", "time-integration"],
    "17": ["boundary-value", "differential-equation"],
    "19": ["pde", "differential-equation", "stencil"],
}

DESCRIPTIONS = {
    "julday": "Compute Julian Day Number from a calendar date.",
    "caldat": "Compute calendar date from a Julian Day Number.",
    "flmoon": "Calculate date and time of phases of the moon.",
    "badluk": "Find all Friday-the-13th dates in a given range.",
    "gaussj": "Gauss-Jordan elimination for matrix inversion and solving linear systems Ax=b.",
    "ludcmp": "LU decomposition of a square matrix using Crout's algorithm with partial pivoting.",
    "lubksb": "Solve Ax=b by LU back-substitution given the LU decomposition from ludcmp.",
    "mprove": "Iterative improvement of a solution vector obtained by LU decomposition.",
    "svdcmp": "Singular Value Decomposition of a rectangular matrix.",
    "svbksb": "Solve Ax=b using SVD back-substitution given the decomposition from svdcmp.",
    "tridag": "Solve a tridiagonal system of linear equations.",
    "toeplz": "Solve a Toeplitz system of linear equations using Levinson's method.",
    "sparse": "Solve sparse linear systems using the biconjugate gradient method.",
    "vander": "Solve a Vandermonde system of linear equations.",
    "polint": "Polynomial interpolation using Neville's algorithm.",
    "ratint": "Rational function interpolation using diagonal rational approximation.",
    "spline": "Compute cubic spline interpolation coefficients.",
    "splint": "Evaluate a cubic spline interpolation at a given point.",
    "splie2": "Construct 2D cubic spline coefficients from tabulated function values.",
    "splin2": "Evaluate a 2D bicubic spline interpolation at a given point.",
    "locate": "Binary search to locate a value in an ordered table.",
    "hunt": "Search an ordered table for a value with a cached initial guess.",
    "polcoe": "Compute polynomial coefficients from tabulated values.",
    "polcof": "Compute polynomial coefficients using successive synthetic division.",
    "polin2": "2D polynomial interpolation on a grid using repeated 1D interpolation.",
    "bcucof": "Compute coefficients for bicubic interpolation on a grid cell.",
    "bcuint": "Bicubic interpolation within a grid cell using bcucof.",
    "trapzd": "Trapezoidal rule integration (one refinement step).",
    "qtrap": "Integration using iterated trapezoidal rule to convergence.",
    "qsimp": "Integration using Simpson's rule to convergence.",
    "qromb": "Romberg integration using Richardson extrapolation on the trapezoidal rule.",
    "midpnt": "Midpoint rule integration (one refinement step, open formula).",
    "midinf": "Midpoint rule integration mapped for semi-infinite intervals.",
    "qromo": "Romberg integration using open midpoint-type rules.",
    "qgaus": "Gauss-Legendre 10-point quadrature on a finite interval.",
    "gauleg": "Compute Gauss-Legendre quadrature weights and abscissas.",
    "quad3d": "Three-dimensional integration by repeated Gaussian quadrature.",
    "chebft": "Chebyshev fit: compute Chebyshev coefficients of a function.",
    "chebev": "Evaluate a Chebyshev series approximation at a given point.",
    "chder": "Compute the derivative of a Chebyshev series.",
    "chint": "Compute the integral of a Chebyshev series.",
    "chebpc": "Convert Chebyshev coefficients to polynomial coefficients.",
    "pcshft": "Shift a Chebyshev polynomial from [-1,1] to an arbitrary interval.",
    "ddpoly": "Evaluate a polynomial and all its derivatives at a given point.",
    "poldiv": "Polynomial division producing quotient and remainder.",
    "eulsum": "Euler's transformation for accelerating alternating series.",
    "gammln": "Log of the gamma function using the Lanczos approximation.",
    "factrl": "Factorial function n! using table lookup and gammln.",
    "factln": "Log-factorial ln(n!) using gammln.",
    "bico": "Binomial coefficient C(n,k) using factln.",
    "beta": "Beta function B(z,w) using gammln.",
    "gammp": "Incomplete gamma function P(a,x).",
    "gammq": "Complementary incomplete gamma function Q(a,x) = 1 - P(a,x).",
    "gser": "Incomplete gamma function by series expansion (used by gammp).",
    "gcf": "Incomplete gamma function by continued fraction (used by gammq).",
    "erf": "Error function erf(x) using the incomplete gamma function.",
    "erfc": "Complementary error function erfc(x) = 1 - erf(x).",
    "erfcc": "Complementary error function using rational Chebyshev approximation.",
    "betai": "Incomplete beta function I_x(a,b).",
    "betacf": "Continued fraction evaluation for the incomplete beta function.",
    "bessj0": "Bessel function J0(x) using rational approximation.",
    "bessj1": "Bessel function J1(x) using rational approximation.",
    "bessj": "Bessel function Jn(x) for integer order n by Miller's downward recurrence.",
    "bessy0": "Bessel function Y0(x) using rational approximation.",
    "bessy1": "Bessel function Y1(x) using rational approximation.",
    "bessy": "Bessel function Yn(x) for integer order n by upward recurrence from Y0, Y1.",
    "bessi0": "Modified Bessel function I0(x) using polynomial approximation.",
    "bessi1": "Modified Bessel function I1(x) using polynomial approximation.",
    "bessi": "Modified Bessel function In(x) for integer order n by downward recurrence.",
    "bessk0": "Modified Bessel function K0(x) using polynomial approximation.",
    "bessk1": "Modified Bessel function K1(x) using polynomial approximation.",
    "bessk": "Modified Bessel function Kn(x) for integer order n by upward recurrence.",
    "plgndr": "Associated Legendre polynomial P_l^m(x).",
    "sncndn": "Jacobian elliptic functions sn(u|m), cn(u|m), dn(u|m).",
    "cel": "General complete elliptic integral.",
    "el2": "General elliptic integral of the second kind.",
    "ran0": "Minimal Park-Miller random number generator.",
    "ran1": "Minimal random number generator with Bays-Durham shuffle.",
    "ran2": "Long-period random number generator combining two linear congruential generators.",
    "ran3": "Knuth's subtractive random number generator.",
    "ran4": "Random number generator using the Data Encryption Standard.",
    "expdev": "Generate exponentially distributed random deviates.",
    "gasdev": "Generate Gaussian (normal) random deviates via Box-Muller transform.",
    "gamdev": "Generate gamma-distributed random deviates.",
    "poidev": "Generate Poisson-distributed random deviates.",
    "bnldev": "Generate binomially distributed random deviates.",
    "irbit1": "Generate random bit sequence using a linear shift register.",
    "irbit2": "Generate random bit sequence using a primitive polynomial modulo 2.",
    "des": "Data Encryption Standard (DES) block cipher implementation.",
    "desks": "DES key schedule setup.",
    "piksrt": "Straight insertion sort for an array.",
    "piksr2": "Insertion sort of one array while rearranging a second in parallel.",
    "shell": "Shell sort for an array.",
    "sort": "Heapsort for an array.",
    "sort2": "Heapsort of one array while rearranging a second in parallel.",
    "sort3": "Sort by constructing an index table, rearranging multiple arrays.",
    "indexx": "Construct an index table (indirect sort) via heapsort.",
    "rank": "Construct a rank table from an index table.",
    "qcksrt": "Quicksort for an array.",
    "eclass": "Determine equivalence classes from pair labels.",
    "eclazz": "Determine equivalence classes using an equivalence-testing function.",
    "mdian1": "Find the median of an array via full sort.",
    "mdian2": "Find the median via partial quicksort (selection algorithm).",
    "zbrac": "Bracket a root of a function by expanding outward.",
    "zbrak": "Search for root-bracketing intervals on a uniform grid.",
    "rtbis": "Find a root by bisection.",
    "rtflsp": "Find a root by false position (regula falsi).",
    "rtsec": "Find a root by the secant method.",
    "rtnewt": "Find a root by Newton-Raphson iteration.",
    "rtsafe": "Find a root by Newton-Raphson with bisection safeguard.",
    "zbrent": "Find a root by Brent's method (inverse quadratic interpolation).",
    "laguer": "Find a root of a polynomial using Laguerre's method.",
    "zroots": "Find all roots of a polynomial by successive Laguerre deflation.",
    "qroot": "Find a quadratic factor of a polynomial by Bairstow's method.",
    "mnewt": "Multi-dimensional Newton-Raphson root finding for nonlinear systems.",
    "scrsho": "Display a function graph for visual root identification.",
    "golden": "Minimize a function of one variable by golden section search.",
    "brent": "Minimize a function of one variable by Brent's method (parabolic interpolation).",
    "dbrent": "Minimize a function using Brent's method with first-derivative information.",
    "mnbrak": "Bracket a minimum of a function of one variable.",
    "linmin": "One-dimensional line minimization along a direction in N-D space.",
    "f1dim": "Utility for linmin: project an N-D function onto a line.",
    "df1dim": "Utility for linmin with derivatives: project the gradient onto a line.",
    "frprmn": "Minimize in N-D using Fletcher-Reeves-Polak-Ribiere conjugate gradient.",
    "dfpmin": "Minimize in N-D using Broyden-Fletcher-Goldfarb-Shanno (BFGS) quasi-Newton.",
    "powell": "Minimize in N-D using Powell's direction-set method (no derivatives).",
    "amoeba": "Minimize in N-D using the downhill simplex (Nelder-Mead) method.",
    "simplx": "Linear programming using the simplex method.",
    "simp1": "Simplex utility: locate the pivot column.",
    "simp2": "Simplex utility: locate the pivot row.",
    "simp3": "Simplex utility: perform a pivot exchange.",
    "anneal": "Combinatorial optimization by simulated annealing (TSP example).",
    "link": "Linked-list utility for the simulated-annealing TSP solver.",
    "jacobi": "Eigenvalues and eigenvectors of a real symmetric matrix by Jacobi rotations.",
    "eigsrt": "Sort eigenvalues and eigenvectors into descending order.",
    "tred2": "Householder reduction of a real symmetric matrix to tridiagonal form.",
    "tqli": "Eigenvalues and eigenvectors of a symmetric tridiagonal matrix by QL iteration.",
    "balanc": "Balance a nonsymmetric matrix to improve eigenvalue accuracy.",
    "elmhes": "Reduce a general matrix to upper Hessenberg form by elimination.",
    "hqr": "Eigenvalues of an upper Hessenberg matrix by the QR algorithm.",
    "four1": "Cooley-Tukey radix-2 FFT for complex data (in-place, interleaved real/imag).",
    "realft": "FFT of real-valued data using the complex FFT (four1).",
    "twofft": "Simultaneous FFT of two real functions packed into one complex FFT.",
    "sinft": "Discrete sine transform using realft.",
    "cosft": "Discrete cosine transform using realft.",
    "fourn": "Multi-dimensional FFT (N-dimensional generalization of four1).",
    "convlv": "FFT-based convolution and deconvolution of real data.",
    "correl": "FFT-based cross-correlation of two real data sets.",
    "spctrm": "Power spectrum estimation using windowed and segmented FFT.",
    "memcof": "Maximum entropy (Burg) method: compute spectral coefficients.",
    "evlmem": "Evaluate the power spectrum from maximum entropy coefficients.",
    "fixrts": "Stabilize a polynomial by moving roots inside the unit circle.",
    "predic": "Linear prediction of future data using maximum entropy coefficients.",
    "smooft": "Smooth data using an FFT-based low-pass filter.",
    "moment": "Compute mean, variance, skewness, and kurtosis of a data set.",
    "avevar": "Compute mean and variance of a data set in one pass.",
    "ttest": "Student's t-test for difference of means (equal variances assumed).",
    "tutest": "Student's t-test for difference of means (unequal variances).",
    "tptest": "Paired-sample Student's t-test.",
    "ftest": "F-test for comparing two variances.",
    "chsone": "Chi-square test: compare binned data to an expected distribution.",
    "chstwo": "Chi-square test: compare two binned data sets.",
    "ksone": "Kolmogorov-Smirnov test: compare data to a known distribution function.",
    "kstwo": "Kolmogorov-Smirnov test: compare two unbinned data sets.",
    "probks": "Kolmogorov-Smirnov probability function Q_KS.",
    "cntab1": "Chi-square contingency table analysis with Cramér's V.",
    "cntab2": "Contingency table analysis using entropy-based measures.",
    "pearsn": "Pearson's linear correlation coefficient and significance.",
    "spear": "Spearman rank-order correlation coefficient.",
    "crank": "Replace data values with their ranks, handling ties.",
    "kendl1": "Kendall's tau rank correlation for two data sets.",
    "kendl2": "Kendall's tau rank correlation for a contingency table.",
    "fit": "Weighted least-squares fit to a straight line y = a + bx.",
    "lfit": "General linear least-squares fit to a set of basis functions.",
    "svdfit": "Least-squares fit using singular value decomposition.",
    "svdvar": "Compute covariance matrix of SVD fit parameters.",
    "covsrt": "Re-sort covariance matrix to include frozen (held) parameters.",
    "fpoly": "Evaluate polynomial basis functions for lfit or svdfit.",
    "fleg": "Evaluate Legendre polynomial basis functions for lfit or svdfit.",
    "fgauss": "Evaluate a sum-of-Gaussians model and its derivatives for mrqmin.",
    "mrqmin": "Levenberg-Marquardt nonlinear least-squares fitting driver.",
    "mrqcof": "Compute the Hessian matrix and gradient vector for Levenberg-Marquardt.",
    "medfit": "Robust least-absolute-deviation (L1-norm) straight line fit.",
    "rofunc": "Evaluate the L1 objective function (used by medfit).",
    "rk4": "Fourth-order Runge-Kutta single-step ODE integrator.",
    "rkdumb": "Integrate ODEs using fixed-step fourth-order Runge-Kutta.",
    "rkqc": "Adaptive step-size fifth-order Runge-Kutta-Cash-Karp ODE step.",
    "odeint": "Driver for ODE integration with adaptive step-size control.",
    "mmid": "Modified midpoint method for ODE integration (used by bsstep).",
    "bsstep": "Bulirsch-Stoer ODE step with Richardson extrapolation.",
    "pzextr": "Polynomial extrapolation for the Bulirsch-Stoer method.",
    "rzextr": "Rational extrapolation for the Bulirsch-Stoer method.",
    "shoot": "Solve a two-point BVP by shooting from one boundary.",
    "shootf": "Solve a two-point BVP by shooting from both boundaries to a fitting point.",
    "sfroid": "Solve the spheroidal wave equation BVP by relaxation.",
    "solvde": "General relaxation solver for two-point boundary value problems.",
    "bksub": "Back-substitution step for the banded system in the relaxation solver.",
    "difeq": "Finite difference equations for the sfroid example.",
    "pinvs": "Diagonalize and back-substitute a section of the banded matrix.",
    "red": "Reduce columns of the banded matrix using results from pinvs.",
    "sor": "Successive Over-Relaxation solver for 2D elliptic PDEs with red-black ordering.",
    "adi": "Alternating Direction Implicit method for 2D parabolic/elliptic PDEs.",
}

CONVERTED_EXAMPLES = {
    ("12_fourier_transform", "four1"),
    ("12_fourier_transform", "realft"),
    ("12_fourier_transform", "twofft"),
    ("13_spectral_analysis", "convlv"),
    ("19_partial_differential_equations", "sor"),
}

LIBRARY_ONLY = {
    ("04_integration", "midinf"),
    ("06_special_functions", "betacf"),
    ("10_optimization", "link"),
    ("10_optimization", "simp1"),
    ("10_optimization", "simp2"),
    ("10_optimization", "simp3"),
    ("17_boundary_value_problems", "bksub"),
    ("17_boundary_value_problems", "difeq"),
    ("17_boundary_value_problems", "pinvs"),
    ("17_boundary_value_problems", "red"),
    ("17_boundary_value_problems", "solvde"),
}

INTRINSICS = {
    'ABS', 'ACOS', 'AIMAG', 'AINT', 'ALOG', 'ALOG10', 'AMAX0', 'AMAX1',
    'AMIN0', 'AMIN1', 'AMOD', 'ANINT', 'ASIN', 'ATAN', 'ATAN2',
    'CABS', 'CCOS', 'CEXP', 'CHAR', 'CLOG', 'CMPLX', 'CONJG', 'COS',
    'COSH', 'CSIN', 'CSQRT', 'DABS', 'DACOS', 'DASIN', 'DATAN', 'DATAN2',
    'DBLE', 'DCOS', 'DCOSH', 'DDIM', 'DEXP', 'DIM', 'DINT', 'DLOG',
    'DLOG10', 'DMAX1', 'DMIN1', 'DMOD', 'DNINT', 'DPROD', 'DSIGN',
    'DSIN', 'DSINH', 'DSQRT', 'DTAN', 'DTANH', 'EXP', 'FLOAT',
    'IABS', 'ICHAR', 'IDIM', 'IDINT', 'IDNINT', 'IFIX', 'INDEX',
    'INT', 'ISIGN', 'LEN', 'LGE', 'LGT', 'LLE', 'LLT', 'LOG',
    'LOG10', 'MAX', 'MAX0', 'MAX1', 'MIN', 'MIN0', 'MIN1', 'MOD',
    'NINT', 'REAL', 'SIGN', 'SIN', 'SINH', 'SNGL', 'SQRT', 'TAN',
    'TANH',
}

# ===================================================================
# Fortran preprocessor
# ===================================================================

def preprocess_fortran(text):
    """Strip comments, join continuation lines.  Returns list of code lines."""
    raw = text.split('\n')
    lines = []
    for line in raw:
        if not line.rstrip():
            continue
        fc = line[0:1]
        if fc in ('C', 'c', '*', '!'):
            continue
        if fc == '\t':
            after = line[1:]
            if after and after[0] in '123456789':
                if lines:
                    lines[-1] += ' ' + after[1:].strip()
                continue
            stmt = line.lstrip('\t').strip()
            if stmt:
                lines.append(stmt)
            continue
        if len(line) > 5 and line[5] not in (' ', '0'):
            if lines:
                lines[-1] += ' ' + line[6:].strip()
            continue
        stmt = line[6:].strip() if len(line) > 6 else line.strip()
        if stmt:
            lines.append(stmt)
    return lines


# ===================================================================
# Fortran source parser
# ===================================================================

_RE_SUB  = re.compile(r'SUBROUTINE\s+(\w+)\s*\(([^)]*)\)', re.I)
_RE_FUNC = re.compile(r'(?:(\w[\w*]*)\s+)?FUNCTION\s+(\w+)\s*\(([^)]*)\)', re.I)
_RE_DIM  = re.compile(r'DIMENSION\s+(.*)', re.I)
_RE_IMPL = re.compile(r'IMPLICIT\s+(.*)', re.I)
_RE_PARAM = re.compile(r'PARAMETER\s*\((.*)\)', re.I)
_RE_COMMON = re.compile(r'COMMON\s*/\s*(\w*)\s*/\s*(.*)', re.I)
_RE_CALL = re.compile(r'\bCALL\s+(\w+)', re.I)
_RE_DO   = re.compile(
    r'DO\s+(\d+)\s+(\w+)\s*=\s*([^,]+),\s*([^,]+)(?:,\s*(.+))?', re.I)
_RE_ARRAY_DECL = re.compile(r'(\w+)\s*\(([^)]+)\)')
_RE_TYPE = re.compile(
    r'(REAL\*8|REAL\*4|REAL|DOUBLE\s+PRECISION|INTEGER\*?\d*|'
    r'COMPLEX\*?\d*|CHARACTER[^A-Z]*|LOGICAL)\s+(.*)', re.I)
_RE_ACCUM = re.compile(r'(\w+)\s*=\s*\1\s*[+*]', re.I)
_RE_WRITE = re.compile(r'(\w+)\s*\(([^)]*)\)\s*=', re.I)


def _parse_decl_list(text):
    """Parse 'A(NP,NP),B(NP,MP),C' into [(name, [dims])]."""
    items, depth, cur = [], 0, ''
    for ch in text:
        if ch == '(':
            depth += 1; cur += ch
        elif ch == ')':
            depth -= 1; cur += ch
        elif ch == ',' and depth == 0:
            items.append(cur.strip()); cur = ''
        else:
            cur += ch
    if cur.strip():
        items.append(cur.strip())
    result = []
    for it in items:
        m = _RE_ARRAY_DECL.match(it)
        if m:
            result.append((m.group(1).upper(), [d.strip() for d in m.group(2).split(',')]))
        else:
            result.append((it.strip().upper(), []))
    return result


def _implicit_type(name, implicit_rules):
    """Determine type of a variable from IMPLICIT rules and Fortran defaults."""
    first = name[0].upper()
    for rule_type, letters in implicit_rules:
        if first in letters:
            return rule_type
    if first in 'IJKLMN':
        return 'INTEGER'
    return 'REAL'


def parse_fortran_source(filepath):
    """Parse a .f file.  Returns dict with keys:
       name, kind, params, variables, loops, calls, accumulators, lines, implicit_rules
    """
    text = filepath.read_text(errors='replace')
    raw_lines = text.split('\n')
    code = preprocess_fortran(text)

    info = {
        'name': '', 'kind': '', 'params': [],
        'variables': {},      # name -> {type, dims, role, is_param, value}
        'loops': [],
        'calls': set(),
        'accumulators': [],
        'array_writes': [],
        'source_lines': len(raw_lines),
        'implicit_rules': [],
        'common_blocks': [],
        'parameters': {},     # name -> value
    }

    for ln in code:
        m = _RE_SUB.search(ln)
        if m:
            info['name'] = m.group(1).upper()
            info['kind'] = 'subroutine'
            info['params'] = [p.strip().upper() for p in m.group(2).split(',') if p.strip()]
            continue
        m = _RE_FUNC.search(ln)
        if m:
            info['name'] = m.group(2).upper()
            info['kind'] = 'function'
            ret_type = m.group(1) or ''
            info['params'] = [p.strip().upper() for p in m.group(3).split(',') if p.strip()]
            if ret_type:
                info['variables'][info['name']] = {
                    'type': ret_type.upper(), 'dims': [], 'role': 'return', 'is_param': False, 'value': ''}
            continue

    for ln in code:
        m = _RE_IMPL.match(ln)
        if m:
            body = m.group(1).strip().upper()
            tm = re.match(r'(REAL\*8|REAL|INTEGER|DOUBLE\s+PRECISION|COMPLEX)\s*\(([^)]+)\)', body, re.I)
            if tm:
                ftype = tm.group(1).upper()
                ranges = tm.group(2)
                letters = set()
                for part in ranges.split(','):
                    part = part.strip()
                    if '-' in part:
                        a, b = part.split('-')
                        for c in range(ord(a.strip().upper()), ord(b.strip().upper()) + 1):
                            letters.add(chr(c))
                    else:
                        letters.add(part.upper())
                info['implicit_rules'].append((ftype, letters))

    for ln in code:
        m = _RE_TYPE.match(ln)
        if m:
            ftype = m.group(1).upper().strip()
            rest = m.group(2)
            for vname, dims in _parse_decl_list(rest):
                if vname in info['variables']:
                    info['variables'][vname]['type'] = ftype
                    if dims:
                        info['variables'][vname]['dims'] = dims
                else:
                    role = 'parameter' if vname in info['params'] else 'local'
                    info['variables'][vname] = {
                        'type': ftype, 'dims': dims, 'role': role,
                        'is_param': False, 'value': ''}

    for ln in code:
        m = _RE_DIM.match(ln)
        if m:
            for vname, dims in _parse_decl_list(m.group(1)):
                if vname in info['variables']:
                    info['variables'][vname]['dims'] = dims
                else:
                    role = 'parameter' if vname in info['params'] else 'local'
                    info['variables'][vname] = {
                        'type': '', 'dims': dims, 'role': role,
                        'is_param': False, 'value': ''}

    for ln in code:
        m = _RE_PARAM.match(ln)
        if m:
            for part in m.group(1).split(','):
                part = part.strip()
                if '=' in part:
                    pn, pv = part.split('=', 1)
                    pn = pn.strip().upper()
                    pv = pv.strip()
                    info['parameters'][pn] = pv
                    if pn in info['variables']:
                        info['variables'][pn]['is_param'] = True
                        info['variables'][pn]['value'] = pv
                    else:
                        ptype = _implicit_type(pn, info['implicit_rules'])
                        info['variables'][pn] = {
                            'type': ptype, 'dims': [], 'role': 'constant',
                            'is_param': True, 'value': pv}

    for ln in code:
        m = _RE_COMMON.match(ln)
        if m:
            bname = m.group(1).upper() or '(blank)'
            vnames = [v.strip().upper() for v in m.group(2).split(',') if v.strip()]
            info['common_blocks'].append({'name': bname, 'vars': vnames})

    for ln in code:
        for m in _RE_CALL.finditer(ln):
            name = m.group(1).upper()
            if name not in INTRINSICS:
                info['calls'].add(name)

    for ln in code:
        m = _RE_DO.match(ln)
        if m:
            info['loops'].append({
                'label': m.group(1),
                'var': m.group(2).upper(),
                'start': m.group(3).strip(),
                'end': m.group(4).strip(),
                'step': (m.group(5) or '1').strip(),
            })

    for ln in code:
        m = _RE_ACCUM.search(ln)
        if m:
            info['accumulators'].append(m.group(1).upper())

    for ln in code:
        for m in _RE_WRITE.finditer(ln):
            arr = m.group(1).upper()
            if arr not in ('WRITE', 'READ', 'OPEN', 'CLOSE', 'FORMAT', 'IF', 'MOD', 'MAX', 'MIN'):
                info['array_writes'].append({'array': arr, 'indices': m.group(2)})

    for p in info['params']:
        if p not in info['variables']:
            ptype = _implicit_type(p, info['implicit_rules'])
            info['variables'][p] = {
                'type': ptype, 'dims': [], 'role': 'parameter',
                'is_param': False, 'value': ''}

    for vname, vinfo in info['variables'].items():
        if not vinfo['type']:
            vinfo['type'] = _implicit_type(vname, info['implicit_rules'])

    return info


# ===================================================================
# Dependency graph from Makefiles
# ===================================================================

def _chapter_num(chapter_dir_name):
    m = re.match(r'(\d+)', chapter_dir_name)
    return m.group(1) if m else ''


def parse_makefile_deps(makefile_path):
    """Return list of {'name', 'chapter', 'path', 'is_self'}."""
    text = makefile_path.read_text(errors='replace')
    text = re.sub(r'\\\s*\n', ' ', text)
    m = re.search(r'^DEPS\s*=\s*(.*)$', text, re.M)
    if not m:
        return []
    tokens = m.group(1).split()
    results = []
    example_dir = makefile_path.parent
    for tok in tokens:
        tok = tok.strip()
        if not tok:
            continue
        abs_p = (example_dir / tok).resolve()
        try:
            rel = abs_p.relative_to(FORTRAN_DIR)
        except ValueError:
            continue
        parts = rel.parts
        if len(parts) >= 2:
            results.append({
                'name': parts[1],
                'chapter': parts[0],
                'path': str(rel),
                'is_self': (parts[1] == example_dir.name),
            })
    return results


def discover_examples():
    """Scan matarized/ for all example directories.  Returns list of dicts."""
    examples = []
    for chapter_dir in sorted(SCRIPT_DIR.iterdir()):
        if not chapter_dir.is_dir():
            continue
        cnum = _chapter_num(chapter_dir.name)
        if cnum not in CHAPTER_TITLES:
            continue
        for ex_dir in sorted(chapter_dir.iterdir()):
            if not ex_dir.is_dir():
                continue
            if ex_dir.name.startswith('.') or ex_dir.name.startswith('build'):
                continue
            examples.append({
                'name': ex_dir.name,
                'chapter': chapter_dir.name,
                'chapter_num': cnum,
                'matarized_dir': ex_dir,
                'fortran_dir': FORTRAN_DIR / chapter_dir.name / ex_dir.name,
            })
    return examples


def build_dependency_graphs(examples):
    """Build forward and reverse dependency maps.

    Returns (forward, reverse) where each is:
        dict of (chapter, name) -> list of (dep_chapter, dep_name)
    """
    forward = {}
    for ex in examples:
        key = (ex['chapter'], ex['name'])
        mf = ex['fortran_dir'] / 'Makefile'
        if mf.exists():
            raw = parse_makefile_deps(mf)
            deps = [(d['chapter'], d['name']) for d in raw if not d['is_self']]
            deps = list(dict.fromkeys(deps))
            forward[key] = deps
        else:
            forward[key] = []

    reverse = defaultdict(list)
    for key, deps in forward.items():
        for dep in deps:
            reverse[dep].append(key)
    return forward, dict(reverse)


def compute_conversion_order(forward, reverse):
    """Return dict (chapter, name) -> int (1 = highest priority)."""
    all_keys = set(forward.keys())
    scored = []
    for k in all_keys:
        rc = len(reverse.get(k, []))
        fc = len(forward.get(k, []))
        scored.append((k, rc, fc))
    scored.sort(key=lambda t: (-t[1], t[2], t[0][1]))
    return {k: i + 1 for i, (k, _, _) in enumerate(scored)}


# ===================================================================
# MATAR type mapper
# ===================================================================

def map_to_matar(vname, vinfo, is_subroutine_param=False):
    """Return recommended MATAR type string for a Fortran variable."""
    ft = vinfo['type'].upper()
    dims = vinfo['dims']

    if vinfo.get('is_param'):
        if 'INTEGER' in ft:
            return f"constexpr int {vname} = {vinfo['value']};"
        return f"constexpr double {vname} = {vinfo['value']};"

    if 'COMPLEX' in ft:
        cpp = 'double'
    elif 'INTEGER' in ft:
        cpp = 'int'
    else:
        cpp = 'double'

    if not dims:
        return cpp

    dims_str = ', '.join(dims)
    return f"DFMatrixKokkos<{cpp}>({dims_str})"


# ===================================================================
# Thread-safety classifier
# ===================================================================

def classify_loops(parsed):
    """Classify each loop for thread safety.  Returns list of dicts."""
    loop_vars = {lp['var'] for lp in parsed['loops']}
    accum_vars = set(parsed['accumulators'])
    results = []
    for lp in parsed['loops']:
        safety = 'safe'
        macro = 'DO_ALL'
        note_parts = []
        if accum_vars:
            safety = 'reduction'
            macro = 'DO_REDUCE_SUM'
            note_parts.append(f"Accumulates: {', '.join(sorted(accum_vars))}")
        has_nonloop_write = False
        has_stencil = False
        nonloop_arrays = set()
        for aw in parsed['array_writes']:
            idx = aw['indices'].upper()
            if not any(v in idx for v in loop_vars):
                has_nonloop_write = True
                nonloop_arrays.add(aw['array'])
            if f'{lp["var"]}+' in idx or f'{lp["var"]}-' in idx:
                has_stencil = True
        if has_nonloop_write:
            safety = 'unsafe_review'
            note_parts.append(
                f"Array write(s) not indexed by loop variable: "
                f"{', '.join(sorted(nonloop_arrays))}. Verify thread safety.")
        if has_stencil and safety == 'safe':
            note_parts.append("Stencil access pattern detected -- verify neighbor independence.")
        results.append({
            'loop': lp,
            'safety': safety,
            'macro': macro,
            'notes': ' '.join(note_parts) if note_parts else 'None',
        })
    return results


# ===================================================================
# CMake template
# ===================================================================

def generate_cmake(example_name, chapter_dir_name, dep_includes):
    """Generate a CMakeLists.txt string.

    dep_includes: list of {'chapter_dir': '12_fourier_transform', 'example': 'four1'}
    """
    project = f"{example_name}_matar_parallel"

    inc_vars = []
    inc_lines = []
    seen = set()
    for d in dep_includes:
        cdir = d['chapter_dir']
        if cdir not in seen:
            seen.add(cdir)
            vname = cdir.split('_', 1)[1].upper().replace('_', '') + '_DIR'
            inc_vars.append((vname, cdir))
        vname = cdir.split('_', 1)[1].upper().replace('_', '') + '_DIR'
        inc_lines.append(f"    ${{{vname}}}/{d['example']}")

    set_lines = ''
    if inc_vars:
        set_lines = '\n'.join(
            f"set({v:<20s} ${{MATARIZED_ROOT}}/{c})" for v, c in inc_vars)

    inc_block = ''
    if inc_lines:
        joined = '\n'.join(inc_lines)
        inc_block = f"""target_include_directories({example_name} PRIVATE
{joined}
)"""

    return f"""cmake_minimum_required(VERSION 3.18)
project({project} CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(FetchContent)

# --- Kokkos backend selection (Serial is always on) ---
set(Kokkos_ENABLE_SERIAL ON CACHE BOOL "Enable Kokkos serial backend")

option(ENABLE_OPENMP "Enable OpenMP backend" OFF)
option(ENABLE_CUDA   "Enable CUDA backend"   OFF)
option(ENABLE_HIP    "Enable HIP backend"    OFF)

if(ENABLE_OPENMP)
    set(Kokkos_ENABLE_OPENMP ON CACHE BOOL "")
endif()
if(ENABLE_CUDA)
    set(Kokkos_ENABLE_CUDA        ON CACHE BOOL "")
    set(Kokkos_ENABLE_CUDA_LAMBDA ON CACHE BOOL "")
endif()
if(ENABLE_HIP)
    set(Kokkos_ENABLE_HIP ON CACHE BOOL "")
endif()

# --- Fetch Kokkos ---
FetchContent_Declare(
    kokkos
    GIT_REPOSITORY https://github.com/kokkos/kokkos.git
    GIT_TAG        4.5.01
    GIT_SHALLOW    TRUE
)
FetchContent_MakeAvailable(kokkos)

# --- Fetch MATAR (header-only -- bypass its CMakeLists.txt) ---
FetchContent_Declare(
    matar
    GIT_REPOSITORY https://github.com/lanl/MATAR.git
    GIT_TAG        main
    GIT_SHALLOW    TRUE
)
FetchContent_GetProperties(matar)
if(NOT matar_POPULATED)
    FetchContent_Populate(matar)
endif()

add_library(matar_lib INTERFACE)
target_include_directories(matar_lib INTERFACE ${{matar_SOURCE_DIR}}/src/include)
target_link_libraries(matar_lib INTERFACE Kokkos::kokkos)
target_compile_definitions(matar_lib INTERFACE HAVE_KOKKOS=1)

# --- Cross-chapter dependency headers ---
set(MATARIZED_ROOT ${{CMAKE_CURRENT_SOURCE_DIR}}/../..)
{set_lines}

# --- Build the {example_name.upper()} example ---
add_executable({example_name} main.cpp)
target_link_libraries({example_name} matar_lib)
{inc_block}
"""


# ===================================================================
# Markdown generator
# ===================================================================

def _mermaid_dep_graph(example_name, forward_deps):
    """Generate a mermaid graph TD string for forward deps."""
    if not forward_deps:
        return f"```mermaid\ngraph TD\n    {example_name}\n```"
    lines = [f"graph TD"]
    for dep_ch, dep_nm in forward_deps:
        lines.append(f"    {example_name} --> {dep_nm}")
    return "```mermaid\n" + "\n".join(lines) + "\n```"


def _variable_table(parsed, info):
    """Generate a markdown table for the variable catalog."""
    rows = []
    for vname in sorted(parsed['variables']):
        v = parsed['variables'][vname]
        dims_str = ', '.join(v['dims']) if v['dims'] else '(scalar)'
        role = v['role']
        if vname in [p.upper() for p in info.get('params', [])]:
            role = 'parameter (input)'
        matar = map_to_matar(vname, v)
        notes = ''
        if v.get('is_param'):
            notes = f"constant = {v['value']}"
        rows.append(f"| `{vname}` | `{v['type']}` | {dims_str} | {role} | `{matar}` | {notes} |")
    if not rows:
        return "_No variables extracted._"
    header = "| Name | Fortran Type | Shape | Role | MATAR Type | Notes |\n"
    header += "|------|-------------|-------|------|-----------|-------|\n"
    return header + "\n".join(rows)


def _loop_analysis_section(loop_info):
    """Generate the compute kernel analysis section."""
    if not loop_info:
        return "_No DO loops detected in source._"
    parts = []
    for i, li in enumerate(loop_info, 1):
        lp = li['loop']
        step_str = f", step {lp['step']}" if lp['step'] != '1' else ''
        parts.append(
            f"### K{i}: DO {lp['label']}  {lp['var']}={lp['start']},{lp['end']}{step_str}\n\n"
            f"- **Thread safety:** `{li['safety']}`\n"
            f"- **Recommended macro:** `{li['macro']}`\n"
            f"- **Notes:** {li['notes'] or 'None'}\n"
        )
    return "\n".join(parts)


def _complexity(parsed, fwd_count):
    loc = parsed.get('source_lines', 0)
    if loc > 60 or fwd_count > 3:
        return 'high'
    if loc > 30 or fwd_count > 1:
        return 'medium'
    return 'low'


def generate_markdown(ex, parsed, forward, reverse, ordering):
    """Produce the full CONVERSION_PLAN.md content."""
    key = (ex['chapter'], ex['name'])
    fwd = forward.get(key, [])
    rev = reverse.get(key, [])
    order_pos = ordering.get(key, 999)
    cnum = ex['chapter_num']
    ctitle = CHAPTER_TITLES.get(cnum, ex['chapter'])
    desc = DESCRIPTIONS.get(ex['name'], "<!-- TODO: add description -->")
    is_converted = key in CONVERTED_EXAMPLES
    is_lib = key in LIBRARY_ONLY
    status = 'converted' if is_converted else ('library_only' if is_lib else 'not_started')
    cpx = _complexity(parsed, len(fwd))
    tags = list(CHAPTER_TAGS.get(cnum, []))
    if is_lib:
        tags.append('library')
    if not fwd:
        tags.append('leaf')
    if any(fc != ex['chapter'] for fc, _ in fwd):
        tags.append('cross-chapter')

    dep_names = [n for _, n in fwd]
    rev_names = [f"{n} ({c})" for c, n in rev]
    prereqs = [n for _, n in fwd]

    fwd_deps_list = "\n".join(f"  - `{n}` ({c})" for c, n in fwd) if fwd else "  (none)"
    rev_deps_list = "\n".join(f"  - `{n}` ({c})" for c, n in rev) if rev else "  (none)"

    # Source files
    f_src = ex['fortran_dir'] / f"{ex['name']}.f"
    f_dem = ex['fortran_dir'] / f"{ex['name']}.dem"
    data_files = []
    if ex['fortran_dir'].exists():
        for p in sorted(ex['fortran_dir'].iterdir()):
            if p.suffix.upper() == '.DAT' or p.name == 'fort.5':
                data_files.append(p.name)

    src_section = f"- **Fortran source:** `fortran/{ex['chapter']}/{ex['name']}/{ex['name']}.f`"
    if f_src.exists():
        src_section += f" ({parsed['source_lines']} lines)"
    src_section += "\n"
    if f_dem.exists():
        src_section += f"- **Driver/demo:** `fortran/{ex['chapter']}/{ex['name']}/{ex['name']}.dem`\n"
    else:
        src_section += "- **Driver/demo:** _(none -- library-only or program-in-.f)_\n"
    if data_files:
        src_section += f"- **Data files:** {', '.join(f'`{d}`' for d in data_files)}\n"
    src_section += f"- **Target:** `matarized/{ex['chapter']}/{ex['name']}/`\n"

    # Variable table
    vtable = _variable_table(parsed, {'params': parsed['params']})

    # Loop analysis
    loop_class = classify_loops(parsed)
    loop_section = _loop_analysis_section(loop_class)

    # Mermaid
    mermaid = _mermaid_dep_graph(ex['name'], fwd)

    # Function signature
    if parsed['kind'] == 'subroutine' and parsed['params']:
        cpp_params = []
        for p in parsed['params']:
            v = parsed['variables'].get(p, {})
            if v.get('dims'):
                cpp_params.append(f"DFMatrixKokkos<double>& {p.lower()}")
            elif 'INTEGER' in v.get('type', ''):
                cpp_params.append(f"int {p.lower()}")
            else:
                cpp_params.append(f"double {p.lower()}")
        sig = f"inline void {ex['name']}({', '.join(cpp_params)})"
    elif parsed['kind'] == 'function':
        ret = 'double'
        if 'INTEGER' in parsed['variables'].get(parsed['name'], {}).get('type', ''):
            ret = 'int'
        params_str = ', '.join(
            f"double {p.lower()}" for p in parsed['params'])
        sig = f"inline {ret} {ex['name']}({params_str})"
    else:
        sig = f"inline void {ex['name']}(/* parameters */)"

    is_header = len(rev) > 0 or is_lib
    file_type = ".hpp header" if is_header else ".cpp with main()"

    # CMake dep includes
    cmake_deps = []
    for dc, dn in fwd:
        cmake_deps.append({'chapter_dir': dc, 'example': dn})
    cmake_content = generate_cmake(ex['name'], ex['chapter'], cmake_deps)

    # Performance
    perf_notes = []
    perf_notes.append("- **FMatrix to CArray migration:** The initial translation uses `DFMatrixKokkos` "
                      "(column-major, 1-based) for Fortran compatibility.  For GPU targets, converting "
                      "to `DCArrayKokkos` (row-major, 0-based) with reordered loops will improve "
                      "coalesced memory access.")
    if parsed['loops']:
        perf_notes.append("- **Loop ordering:** Verify innermost parallel index matches the fastest-varying "
                          "array dimension for the chosen layout.")
    if parsed['accumulators']:
        perf_notes.append("- **Reduction fusion:** If multiple reductions share the same loop bounds, "
                          "consider fusing them into a single pass to reduce kernel launch overhead.")
    perf_notes.append("- **Fence elimination:** After conversion, audit `MATAR_FENCE()` placement.  "
                      "Remove fences between independent kernels that do not share data.")
    if len(parsed['loops']) > 2:
        perf_notes.append("- **Hierarchical parallelism:** For deeply nested loops, consider "
                          "`FOR_FIRST`/`FOR_SECOND` team-thread decomposition for better occupancy.")

    # Validation
    has_dem = f_dem.exists()
    validation = []
    validation.append("### Reference Output\n")
    if has_dem:
        validation.append(f"Build and run the Fortran version to capture reference output:\n")
        validation.append(f"```bash\n"
                          f"cd fortran/{ex['chapter']}/{ex['name']}\n"
                          f"make run > reference_output.txt 2>&1\n"
                          f"```\n")
    else:
        validation.append("_No standalone Fortran driver (.dem) exists.  "
                          "Validate by calling this routine from a dependent example._\n")
    validation.append("\n### Serial Validation\n")
    validation.append(f"```bash\n"
                      f"cd matarized/{ex['chapter']}/{ex['name']}\n"
                      f"mkdir -p build && cd build\n"
                      f"cmake .. && make\n"
                      f"./{ex['name']} > serial_output.txt 2>&1\n"
                      f"diff <(head -50 serial_output.txt) <(head -50 ../../../../fortran/{ex['chapter']}/{ex['name']}/reference_output.txt)\n"
                      f"```\n")
    validation.append("\n### Parallel Validation (OpenMP)\n")
    validation.append(f"```bash\n"
                      f"cd matarized/{ex['chapter']}/{ex['name']}\n"
                      f"mkdir -p build-omp && cd build-omp\n"
                      f"cmake .. -DENABLE_OPENMP=ON && make\n"
                      f"OMP_NUM_THREADS=1 ./{ex['name']} > omp1_output.txt 2>&1\n"
                      f"OMP_NUM_THREADS=4 ./{ex['name']} > omp4_output.txt 2>&1\n"
                      f"# Verify: omp1 output must exactly match serial output\n"
                      f"diff serial_output.txt omp1_output.txt\n"
                      f"# Verify: omp4 output must match within floating-point tolerance\n"
                      f"```\n")
    validation.append("\n### Pass Criteria\n")
    validation.append("- Max absolute difference vs. Fortran reference: **< 1e-10** (double precision)\n")
    validation.append("- OpenMP results must be deterministic across repeated runs\n")
    validation.append("- No runtime errors, memory leaks, or Kokkos warnings\n")

    # Conversion steps
    steps = []
    steps.append("1. **Translate data structures** -- replace Fortran arrays with `DFMatrixKokkos` "
                 "(see variable catalog below)")
    steps.append(f"2. **Translate routine** -- convert `{parsed['name'] or ex['name'].upper()}` "
                 f"to a C++ function as a `{file_type}`")
    steps.append("3. **Replace loops** -- convert DO loops to `DO_ALL` / `DO_REDUCE_*` macros "
                 "(see kernel analysis below)")
    steps.append("4. **Add synchronization** -- insert `MATAR_FENCE()` between dependent kernels; "
                 "add `update_host()`/`update_device()` for Dual types")
    if has_dem:
        steps.append("5. **Create driver** -- translate the `.dem` test program to `main.cpp` with "
                     "`MATAR_INITIALIZE` / `MATAR_FINALIZE` boilerplate")
    steps.append(f"6. **Generate CMakeLists.txt** -- use the template below (based on convlv reference)")
    steps.append("7. **Validate** -- follow the validation plan below")

    # --- Assemble markdown ---
    md = f"""---
example: {ex['name']}
chapter: {ex['chapter']}
chapter_title: "{ctitle}"
status: {status}
complexity: {cpx}
conversion_order: {order_pos}
priority: {len(rev)}
tags: [{', '.join(tags)}]
dependencies: [{', '.join(dep_names)}]
reverse_dependencies: [{', '.join(n for _, n in rev)}]
---

# {ex['name'].upper()} -- {ctitle}

## 1. Overview

| Field | Value |
|-------|-------|
| **Example** | `{ex['name']}` |
| **Chapter** | {cnum} -- {ctitle} |
| **Purpose** | {desc} |
| **Status** | `{status}` |
| **Complexity** | `{cpx}` |
| **Fortran LOC** | {parsed['source_lines']} |
| **Subroutine** | `{parsed['name'] or '(none)'}` ({parsed['kind'] or 'N/A'}) |

## 2. Source Files

{src_section}

## 3. Dependency Graph

### Forward Dependencies (this example depends on)

{fwd_deps_list}

### Diagram

{mermaid}

### Cross-Chapter Dependencies

{chr(10).join(f'- `{n}` from chapter {_chapter_num(c)}' for c, n in fwd if c != ex['chapter']) or '(none)'}

## 4. Reverse Dependencies (examples that depend on this)

{rev_deps_list}

> **Conversion note:** {"This routine is depended on by " + str(len(rev)) + " other examples and should be converted early." if rev else "No other examples depend on this routine."}

## 5. Fortran Variable Catalog

{vtable}

### MATAR Type Mapping Rationale

- **Layout:** `FMatrix` (column-major) preserves Fortran memory layout for correctness.
- **Index base:** `Matrix` (1-based) matches Fortran indexing with `DO_ALL` inclusive ranges.
- **Residence:** `Dual` (`DFMatrixKokkos`) enables both host I/O and device computation.
- **Ownership:** Owning types at call site; consider `ViewFMatrix` for sub-array slices.

## 6. Compute Kernel Analysis

{loop_section}

### Thread-Safety Legend

| Classification | Meaning | Action |
|---------------|---------|--------|
| `safe` | No write conflicts | Parallelize directly with `DO_ALL` |
| `reduction` | Accumulation to scalar | Use `DO_REDUCE_SUM` / `DO_REDUCE_MAX` |
| `unsafe_review` | Potential race condition | Restructure: inner serial loop or phased approach |
| `inherently_serial` | Sequential data dependency | Keep as serial `for` inside parallel region |

## 7. Conversion Strategy

### Proposed C++ Signature

```cpp
{sig}
```

### Output Format

- **{file_type}** {'(included by other examples via `#include`)' if is_header else '(standalone executable)'}

### Steps

{chr(10).join(steps)}

## 8. CMake Configuration

Based on the [convlv CMakeLists.txt](../../13_spectral_analysis/convlv/CMakeLists.txt) reference template.

```cmake
{cmake_content}```

## 9. Performance Improvements

{chr(10).join(perf_notes)}

## 10. Validation Plan

{chr(10).join(validation)}

## 11. Agent Metadata

| Field | Value |
|-------|-------|
| **Conversion order** | {order_pos} of {len(ordering)} |
| **Priority score** | {len(rev)} (reverse dependency count) |
| **Estimated effort** | {cpx} ({parsed['source_lines']} Fortran LOC, {len(fwd)} dependencies) |
| **Prerequisite conversions** | {', '.join(f'`{n}`' for n in prereqs) or '(none -- leaf node)'} |
| **Tags** | {', '.join(f'`{t}`' for t in tags)} |
| **MATAR reference sections** | {'Sec 5 (parallel loops), Sec 6 (reductions)' if parsed['accumulators'] else 'Sec 5 (parallel loops)'}{', Sec 15 (Fortran interop)' if parsed['kind'] else ''} |
"""
    return md


# ===================================================================
# Main
# ===================================================================

def main():
    dry_run = '--dry-run' in sys.argv
    target = None
    for a in sys.argv[1:]:
        if not a.startswith('-'):
            target = a

    print("=== CONVERSION_PLAN.md Generator ===\n")

    print("Discovering examples ...")
    examples = discover_examples()
    print(f"  Found {len(examples)} examples across "
          f"{len({e['chapter'] for e in examples})} chapters\n")

    print("Building dependency graphs ...")
    forward, reverse = build_dependency_graphs(examples)
    total_deps = sum(len(v) for v in forward.values())
    print(f"  {total_deps} forward dependency edges")
    print(f"  {sum(len(v) for v in reverse.values())} reverse dependency edges\n")

    print("Computing conversion ordering ...")
    ordering = compute_conversion_order(forward, reverse)

    generated = 0
    skipped = 0

    for ex in examples:
        if target and ex['name'] != target:
            skipped += 1
            continue

        f_src = ex['fortran_dir'] / f"{ex['name']}.f"
        if f_src.exists():
            parsed = parse_fortran_source(f_src)
        else:
            parsed = {
                'name': ex['name'].upper(), 'kind': '', 'params': [],
                'variables': {}, 'loops': [], 'calls': set(),
                'accumulators': [], 'array_writes': [],
                'source_lines': 0, 'implicit_rules': [],
                'common_blocks': [], 'parameters': {},
            }

        md = generate_markdown(ex, parsed, forward, reverse, ordering)

        out_path = ex['matarized_dir'] / 'CONVERSION_PLAN.md'
        if dry_run:
            print(f"  [DRY-RUN] {ex['chapter']}/{ex['name']}/CONVERSION_PLAN.md  "
                  f"({len(md)} chars)")
        else:
            out_path.write_text(md)
            print(f"  Wrote {out_path.relative_to(ROOT_DIR)}")
        generated += 1

    print(f"\n{'Would generate' if dry_run else 'Generated'} {generated} plans"
          f"{f' (skipped {skipped})' if skipped else ''}.")

    if not dry_run and not target:
        top = sorted(ordering.items(), key=lambda t: t[1])[:15]
        print("\nTop-15 conversion priority (most depended-on first):")
        for (ch, nm), pos in top:
            rc = len(reverse.get((ch, nm), []))
            print(f"  {pos:3d}. {nm:12s} ({_chapter_num(ch)}) -- {rc} consumers")


if __name__ == '__main__':
    main()
