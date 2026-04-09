FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic

PROG    = mrqmin
DRIVER  = mrqmin.dem

# Subroutine sources (transitive dependencies)
DEPS    = ../../02_linear_algebra/gaussj/gaussj.f \
          ../../06_special_functions/beta/beta.f \
          ../../06_special_functions/gammln/gammln.f \
          ../../07_random_numbers/gasdev/gasdev.f \
          ../../07_random_numbers/ran1/ran1.f \
          ../covsrt/covsrt.f \
          ../fgauss/fgauss.f \
          ../mrqcof/mrqcof.f \
          mrqmin.f

$(PROG): $(DRIVER) $(DEPS)
	$(FC) $(FFLAGS) -x f77 $(DRIVER) -x none $(DEPS) -o $@

.PHONY: clean run
clean:
	rm -f $(PROG)
run: $(PROG)
	./$(PROG)
