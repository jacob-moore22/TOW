FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic
AR      = ar
ARFLAGS = rcs

BUILDDIR = build

# Every .f that defines a SUBROUTINE or FUNCTION (not PROGRAM) goes into the library.
# badluk.f and sfroid.f contain PROGRAM and are excluded.
LIB_SRCS := $(filter-out %/badluk.f %/sfroid.f, \
               $(shell find . -path ./$(BUILDDIR) -prune -o -name '*.f' -print))
LIB_OBJS := $(addprefix $(BUILDDIR)/, $(notdir $(LIB_SRCS:.f=.o)))

VPATH := $(sort $(dir $(LIB_SRCS)))

LIB = libnrf77.a

# Discover all per-program Makefiles (depth 3: ./<chapter>/<prog>/Makefile)
PROG_DIRS := $(sort $(shell find . -mindepth 3 -maxdepth 3 -name Makefile -exec dirname {} \;))

# ── Targets ──────────────────────────────────────────────────────────────────

.PHONY: all lib programs clean

all: lib programs

lib: $(LIB)

$(LIB): $(LIB_OBJS)
	$(AR) $(ARFLAGS) $@ $^

$(BUILDDIR)/%.o: %.f | $(BUILDDIR)
	$(FC) $(FFLAGS) -c -o $@ $<

$(BUILDDIR):
	mkdir -p $@

programs: $(LIB)
	@fail=0; \
	for d in $(PROG_DIRS); do \
		$(MAKE) --no-print-directory -C $$d || fail=$$((fail+1)); \
	done; \
	echo ""; \
	if [ $$fail -gt 0 ]; then \
		echo "*** $$fail program(s) failed to build ***"; \
	else \
		echo "All programs built successfully"; \
	fi

clean:
	rm -rf $(BUILDDIR) $(LIB)
	@for d in $(PROG_DIRS); do $(MAKE) --no-print-directory -C $$d clean 2>/dev/null; done
