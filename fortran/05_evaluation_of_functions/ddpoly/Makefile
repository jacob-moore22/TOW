FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic

PROG    = ddpoly
DRIVER  = ddpoly.dem

# Subroutine sources (transitive dependencies)
DEPS    = ddpoly.f \
          ../../06_special_functions/factrl/factrl.f \
          ../../06_special_functions/gammln/gammln.f

$(PROG): $(DRIVER) $(DEPS)
	$(FC) $(FFLAGS) -x f77 $(DRIVER) -x none $(DEPS) -o $@

.PHONY: clean run
clean:
	rm -f $(PROG)
run: $(PROG)
	./$(PROG)
