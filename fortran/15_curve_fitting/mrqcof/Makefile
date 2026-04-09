FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic

PROG    = mrqcof
DRIVER  = mrqcof.dem

# Subroutine sources (transitive dependencies)
DEPS    = ../../06_special_functions/beta/beta.f \
          ../../06_special_functions/gammln/gammln.f \
          ../../07_random_numbers/gasdev/gasdev.f \
          ../../07_random_numbers/ran1/ran1.f \
          ../fgauss/fgauss.f \
          mrqcof.f

$(PROG): $(DRIVER) $(DEPS)
	$(FC) $(FFLAGS) -x f77 $(DRIVER) -x none $(DEPS) -o $@

.PHONY: clean run
clean:
	rm -f $(PROG)
run: $(PROG)
	./$(PROG)
