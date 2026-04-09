FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic

PROG    = ksone
DRIVER  = ksone.dem

# Subroutine sources (transitive dependencies)
DEPS    = ../../06_special_functions/erf/erf.f \
          ../../06_special_functions/gammln/gammln.f \
          ../../06_special_functions/gammp/gammp.f \
          ../../06_special_functions/gcf/gcf.f \
          ../../06_special_functions/gser/gser.f \
          ../../07_random_numbers/gasdev/gasdev.f \
          ../../07_random_numbers/ran1/ran1.f \
          ../../08_sorting/sort/sort.f \
          ksone.f \
          ../probks/probks.f

$(PROG): $(DRIVER) $(DEPS)
	$(FC) $(FFLAGS) -x f77 $(DRIVER) -x none $(DEPS) -o $@

.PHONY: clean run
clean:
	rm -f $(PROG)
run: $(PROG)
	./$(PROG)
