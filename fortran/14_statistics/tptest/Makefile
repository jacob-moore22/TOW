FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic

PROG    = tptest
DRIVER  = tptest.dem

# Subroutine sources (transitive dependencies)
DEPS    = ../../06_special_functions/betacf/betacf.f \
          ../../06_special_functions/betai/betai.f \
          ../../06_special_functions/gammln/gammln.f \
          ../../07_random_numbers/gasdev/gasdev.f \
          ../../07_random_numbers/ran1/ran1.f \
          ../avevar/avevar.f \
          tptest.f

$(PROG): $(DRIVER) $(DEPS)
	$(FC) $(FFLAGS) -x f77 $(DRIVER) -x none $(DEPS) -o $@

.PHONY: clean run
clean:
	rm -f $(PROG)
run: $(PROG)
	./$(PROG)
