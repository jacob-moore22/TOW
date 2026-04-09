FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic

PROG    = odeint
DRIVER  = odeint.dem

# Subroutine sources (transitive dependencies)
DEPS    = ../../06_special_functions/bessj/bessj.f \
          ../../06_special_functions/bessj0/bessj0.f \
          ../../06_special_functions/bessj1/bessj1.f \
          odeint.f \
          ../rk4/rk4.f \
          ../rkqc/rkqc.f

$(PROG): $(DRIVER) $(DEPS)
	$(FC) $(FFLAGS) -x f77 $(DRIVER) -x none $(DEPS) -o $@

.PHONY: clean run
clean:
	rm -f $(PROG)
run: $(PROG)
	./$(PROG)
