FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic

PROG    = factrl
DRIVER  = factrl.dem

# Subroutine sources (transitive dependencies)
DEPS    = factrl.f \
          ../gammln/gammln.f

$(PROG): $(DRIVER) $(DEPS)
	$(FC) $(FFLAGS) -x f77 $(DRIVER) -x none $(DEPS) -o $@

.PHONY: clean run
clean:
	rm -f $(PROG)
run: $(PROG)
	./$(PROG)
