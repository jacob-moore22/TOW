FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic

PROG    = spear
DRIVER  = spear.dem

# Subroutine sources (transitive dependencies)
DEPS    = ../../06_special_functions/betacf/betacf.f \
          ../../06_special_functions/betai/betai.f \
          ../../06_special_functions/erfcc/erfcc.f \
          ../../06_special_functions/gammln/gammln.f \
          ../../08_sorting/sort2/sort2.f \
          ../crank/crank.f \
          spear.f

$(PROG): $(DRIVER) $(DEPS)
	$(FC) $(FFLAGS) -x f77 $(DRIVER) -x none $(DEPS) -o $@

.PHONY: clean run
clean:
	rm -f $(PROG)
run: $(PROG)
	./$(PROG)
