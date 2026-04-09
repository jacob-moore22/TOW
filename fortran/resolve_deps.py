#!/usr/bin/env python3
"""
Resolve per-example Fortran subroutine/function dependencies and generate
self-contained Makefiles that reference exact .f files via relative paths.

Usage:
    python3 resolve_deps.py          # from the fortran/ directory
    python3 resolve_deps.py --dry-run  # print what would be generated
"""

import os
import re
import sys
from collections import defaultdict
from pathlib import Path

FORTRAN_DIR = Path(__file__).resolve().parent

# Regex patterns for parsing Fortran 77 sources
RE_SUB_DEF = re.compile(
    r'^\s+SUBROUTINE\s+(\w+)', re.IGNORECASE | re.MULTILINE
)
RE_FUNC_DEF = re.compile(
    r'^\s+(?:\w+\s+)*FUNCTION\s+(\w+)', re.IGNORECASE | re.MULTILINE
)
RE_PROGRAM_DEF = re.compile(
    r'^\s+PROGRAM\s+(\w+)', re.IGNORECASE | re.MULTILINE
)
RE_CALL = re.compile(
    r'\bCALL\s+(\w+)', re.IGNORECASE
)
RE_EXTERNAL = re.compile(
    r'^\s+EXTERNAL\s+(.+)', re.IGNORECASE | re.MULTILINE
)

# Intrinsics and builtins that should never be resolved as library deps
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
    'TANH', 'FABS', 'FMAX', 'FMIN', 'WRITE', 'READ', 'OPEN', 'CLOSE',
    'PRINT', 'FORMAT', 'STOP', 'PAUSE', 'RETURN', 'CONTINUE', 'GOTO',
    'IF', 'THEN', 'ELSE', 'ENDIF', 'DO', 'ENDDO', 'END',
    'MPI_INIT', 'MPI_FINALIZE', 'MPI_COMM_RANK', 'MPI_COMM_SIZE',
}


def strip_comments(text):
    """Remove Fortran 77 comment lines (C, c, or * in column 1)."""
    lines = []
    for line in text.split('\n'):
        if line and line[0] in ('C', 'c', '*', '!'):
            continue
        lines.append(line)
    return '\n'.join(lines)


def find_definitions(text):
    """Return set of SUBROUTINE/FUNCTION names defined in the source."""
    code = strip_comments(text)
    names = set()
    for m in RE_SUB_DEF.finditer(code):
        names.add(m.group(1).upper())
    for m in RE_FUNC_DEF.finditer(code):
        names.add(m.group(1).upper())
    return names


def find_program(text):
    """Return PROGRAM name if present, else None."""
    code = strip_comments(text)
    m = RE_PROGRAM_DEF.search(code)
    return m.group(1).upper() if m else None


def find_calls(text):
    """Return set of names that appear in CALL statements."""
    code = strip_comments(text)
    return {m.group(1).upper() for m in RE_CALL.finditer(code)}


def find_externals(text):
    """Return set of names declared EXTERNAL."""
    code = strip_comments(text)
    names = set()
    for m in RE_EXTERNAL.finditer(code):
        for name in re.split(r'[,\s]+', m.group(1).strip()):
            name = name.strip().upper()
            if name:
                names.add(name)
    return names


def find_function_refs(text, known_functions):
    """Find references to known library functions used as `NAME(` without CALL.

    Excludes single-character names (F, G, H, etc.) which produce false
    positives against array accesses and format specifiers.
    """
    code = strip_comments(text)
    refs = set()
    for fname in known_functions:
        if len(fname) <= 1:
            continue
        pattern = re.compile(r'\b' + re.escape(fname) + r'\s*\(', re.IGNORECASE)
        if pattern.search(code):
            refs.add(fname)
    return refs


def scan_all_sources(fortran_dir):
    """
    Scan all .f files to build:
      - name_to_file: maps uppercase routine name -> relative .f path
      - file_deps:    maps relative .f path -> set of routine names it calls
      - known_functions: set of names defined as FUNCTION (not SUBROUTINE)
    """
    name_to_file = {}
    file_deps = {}
    known_functions = set()

    for fpath in sorted(fortran_dir.rglob('*.f')):
        rel = fpath.relative_to(fortran_dir)
        text = fpath.read_text(errors='replace')
        code = strip_comments(text)

        defs = find_definitions(text)
        for name in defs:
            name_to_file[name] = str(rel)

        for m in RE_FUNC_DEF.finditer(code):
            known_functions.add(m.group(1).upper())

        calls = find_calls(text)
        externals = find_externals(text)
        deps = (calls | externals) - defs - INTRINSICS
        file_deps[str(rel)] = deps

    return name_to_file, file_deps, known_functions


def resolve_transitive(start_deps, name_to_file, file_deps, known_functions,
                       provided_names=None):
    """
    Given a set of direct dependency names, resolve the transitive closure
    of .f files needed. Returns a list of relative .f paths (sorted).

    provided_names: names already supplied by the driver (inline defs) that
    should not be resolved to external .f files.
    """
    provided = provided_names or set()
    needed_files = []
    seen_names = set()
    queue = list(start_deps)

    while queue:
        name = queue.pop(0)
        if name in seen_names or name in INTRINSICS or name in provided:
            continue
        seen_names.add(name)

        fpath = name_to_file.get(name)
        if fpath is None:
            continue

        if fpath not in needed_files:
            needed_files.append(fpath)

        for dep_name in file_deps.get(fpath, set()):
            if dep_name not in seen_names:
                queue.append(dep_name)

        # Also resolve function references within this .f file
        try:
            text = (FORTRAN_DIR / fpath).read_text(errors='replace')
            frefs = find_function_refs(text, known_functions)
            defs_in_file = find_definitions(text)
            for fname in frefs - defs_in_file - seen_names:
                queue.append(fname)
        except FileNotFoundError:
            pass

    return sorted(needed_files)


def analyze_driver(driver_path, name_to_file, file_deps, known_functions):
    """
    Analyze a .dem or special .f PROGRAM file to determine its dependencies.
    Returns (inline_defs, needed_f_files).
    """
    text = driver_path.read_text(errors='replace')

    inline_defs = find_definitions(text)
    calls = find_calls(text)
    externals = find_externals(text)
    func_refs = find_function_refs(text, known_functions)

    direct_deps = (calls | externals | func_refs) - inline_defs - INTRINSICS
    needed = resolve_transitive(direct_deps, name_to_file, file_deps,
                                known_functions, provided_names=inline_defs)
    return inline_defs, needed


def make_relative(from_dir, to_path):
    """Compute relative path from from_dir to to_path."""
    return os.path.relpath(to_path, from_dir)


def generate_makefile(prog_dir, prog_name, driver_name, dep_files, is_f_program=False):
    """Generate Makefile content for a single example."""
    rel_deps = []
    for dep in dep_files:
        abs_dep = FORTRAN_DIR / dep
        rel = make_relative(prog_dir, abs_dep)
        rel_deps.append(rel)

    lines = []
    lines.append('FC      = gfortran')
    lines.append('FFLAGS  = -O2 -std=legacy -w -fno-automatic')
    lines.append('')
    lines.append(f'PROG    = {prog_name}')
    lines.append(f'DRIVER  = {driver_name}')
    lines.append('')

    if rel_deps:
        lines.append('# Subroutine sources (transitive dependencies)')
        lines.append(f'DEPS    = {rel_deps[0]}')
        for dep in rel_deps[1:]:
            lines[-1] += ' \\'
            lines.append(f'          {dep}')
    else:
        lines.append('DEPS    =')

    lines.append('')
    if is_f_program:
        lines.append('$(PROG): $(DRIVER) $(DEPS)')
        lines.append('\t$(FC) $(FFLAGS) $(DRIVER) $(DEPS) -o $@')
    else:
        lines.append('$(PROG): $(DRIVER) $(DEPS)')
        lines.append('\t$(FC) $(FFLAGS) -x f77 $(DRIVER) -x none $(DEPS) -o $@')

    lines.append('')
    lines.append('.PHONY: clean run')
    lines.append('clean:')
    lines.append('\trm -f $(PROG)')
    lines.append('run: $(PROG)')
    lines.append('\t./$(PROG)')
    lines.append('')

    return '\n'.join(lines)


def main():
    dry_run = '--dry-run' in sys.argv

    print("Scanning all .f files for definitions and dependencies...")
    name_to_file, file_deps, known_functions = scan_all_sources(FORTRAN_DIR)
    print(f"  Found {len(name_to_file)} subroutine/function definitions")
    print(f"  Found {len(known_functions)} FUNCTION definitions")
    print(f"  Scanned {len(file_deps)} .f files for call dependencies")

    # Find all driver programs
    drivers = []
    for dem_path in sorted(FORTRAN_DIR.rglob('*.dem')):
        drivers.append(('dem', dem_path))

    # Special PROGRAM .f files (badluk, sfroid)
    for fpath in sorted(FORTRAN_DIR.rglob('*.f')):
        text = fpath.read_text(errors='replace')
        prog = find_program(text)
        if prog is not None:
            # Check it's not also a .dem-driven example
            dem_sibling = fpath.with_suffix('.dem')
            if not dem_sibling.exists():
                drivers.append(('f_program', fpath))

    print(f"\nProcessing {len(drivers)} driver programs...")

    generated = 0
    errors = []

    for kind, driver_path in drivers:
        prog_dir = driver_path.parent
        prog_name = driver_path.stem
        driver_name = driver_path.name

        try:
            inline_defs, needed_files = analyze_driver(
                driver_path, name_to_file, file_deps, known_functions
            )
        except Exception as e:
            errors.append((str(driver_path), str(e)))
            continue

        is_f_program = (kind == 'f_program')
        content = generate_makefile(
            prog_dir, prog_name, driver_name, needed_files, is_f_program
        )

        makefile_path = prog_dir / 'Makefile'
        if dry_run:
            print(f"\n--- {makefile_path.relative_to(FORTRAN_DIR)} ---")
            dep_summary = ', '.join(Path(f).name for f in needed_files) if needed_files else '(none)'
            print(f"  deps: {dep_summary}")
        else:
            makefile_path.write_text(content)

        generated += 1

    print(f"\n{'Would generate' if dry_run else 'Generated'} {generated} Makefiles")

    if errors:
        print(f"\n{len(errors)} errors:")
        for path, err in errors:
            print(f"  {path}: {err}")

    # Print a summary of unresolved externals for debugging
    print("\nChecking for potentially unresolved dependencies...")
    unresolved_count = 0
    for kind, driver_path in drivers:
        text = driver_path.read_text(errors='replace')
        inline_defs = find_definitions(text)
        calls = find_calls(text)
        externals = find_externals(text)
        func_refs = find_function_refs(text, known_functions)
        all_refs = (calls | externals | func_refs) - inline_defs - INTRINSICS
        missing = {n for n in all_refs if n not in name_to_file}
        if missing:
            rel = driver_path.relative_to(FORTRAN_DIR)
            print(f"  {rel}: unresolved -> {missing}")
            unresolved_count += len(missing)

    if unresolved_count == 0:
        print("  All dependencies resolved!")
    else:
        print(f"  {unresolved_count} total unresolved references (may be inline or user-supplied)")


if __name__ == '__main__':
    main()
