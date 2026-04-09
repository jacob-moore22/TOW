FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic

PROG    = dfpmin
DRIVER  = dfpmin.dem

# Subroutine sources (transitive dependencies)
DEPS    = ../../06_special_functions/bessj0/bessj0.f \
          ../../06_special_functions/bessj1/bessj1.f \
          ../brent/brent.f \
          dfpmin.f \
          ../f1dim/f1dim.f \
          ../linmin/linmin.f \
          ../mnbrak/mnbrak.f

$(PROG): $(DRIVER) $(DEPS)
	$(FC) $(FFLAGS) -x f77 $(DRIVER) -x none $(DEPS) -o $@

.PHONY: clean run
clean:
	rm -f $(PROG)
run: $(PROG)
	./$(PROG)
