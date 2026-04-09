FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic

PROG    = linmin
DRIVER  = linmin.dem

# Subroutine sources (transitive dependencies)
DEPS    = ../brent/brent.f \
          ../f1dim/f1dim.f \
          linmin.f \
          ../mnbrak/mnbrak.f

$(PROG): $(DRIVER) $(DEPS)
	$(FC) $(FFLAGS) -x f77 $(DRIVER) -x none $(DEPS) -o $@

.PHONY: clean run
clean:
	rm -f $(PROG)
run: $(PROG)
	./$(PROG)
