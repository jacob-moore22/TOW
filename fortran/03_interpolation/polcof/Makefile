FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic

PROG    = polcof
DRIVER  = polcof.dem

# Subroutine sources (transitive dependencies)
DEPS    = polcof.f \
          ../polint/polint.f

$(PROG): $(DRIVER) $(DEPS)
	$(FC) $(FFLAGS) -x f77 $(DRIVER) -x none $(DEPS) -o $@

.PHONY: clean run
clean:
	rm -f $(PROG)
run: $(PROG)
	./$(PROG)
