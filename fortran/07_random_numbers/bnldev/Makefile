FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic

PROG    = bnldev
DRIVER  = bnldev.dem

# Subroutine sources (transitive dependencies)
DEPS    = ../../06_special_functions/gammln/gammln.f \
          bnldev.f \
          ../ran1/ran1.f

$(PROG): $(DRIVER) $(DEPS)
	$(FC) $(FFLAGS) -x f77 $(DRIVER) -x none $(DEPS) -o $@

.PHONY: clean run
clean:
	rm -f $(PROG)
run: $(PROG)
	./$(PROG)
