FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic

PROG    = pearsn
DRIVER  = pearsn.dem

# Subroutine sources (transitive dependencies)
DEPS    = ../../06_special_functions/betacf/betacf.f \
          ../../06_special_functions/betai/betai.f \
          ../../06_special_functions/gammln/gammln.f \
          pearsn.f

$(PROG): $(DRIVER) $(DEPS)
	$(FC) $(FFLAGS) -x f77 $(DRIVER) -x none $(DEPS) -o $@

.PHONY: clean run
clean:
	rm -f $(PROG)
run: $(PROG)
	./$(PROG)
