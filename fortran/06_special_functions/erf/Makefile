FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic

PROG    = erf
DRIVER  = erf.dem

# Subroutine sources (transitive dependencies)
DEPS    = erf.f \
          ../gammln/gammln.f \
          ../gammp/gammp.f \
          ../gcf/gcf.f \
          ../gser/gser.f

$(PROG): $(DRIVER) $(DEPS)
	$(FC) $(FFLAGS) -x f77 $(DRIVER) -x none $(DEPS) -o $@

.PHONY: clean run
clean:
	rm -f $(PROG)
run: $(PROG)
	./$(PROG)
