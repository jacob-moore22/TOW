FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic

PROG    = factln
DRIVER  = factln.dem

# Subroutine sources (transitive dependencies)
DEPS    = factln.f \
          ../gammln/gammln.f

$(PROG): $(DRIVER) $(DEPS)
	$(FC) $(FFLAGS) -x f77 $(DRIVER) -x none $(DEPS) -o $@

.PHONY: clean run
clean:
	rm -f $(PROG)
run: $(PROG)
	./$(PROG)
