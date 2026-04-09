FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic

PROG    = shootf
DRIVER  = shootf.dem

# Subroutine sources (transitive dependencies)
DEPS    = ../../02_linear_algebra/lubksb/lubksb.f \
          ../../02_linear_algebra/ludcmp/ludcmp.f \
          ../../16_ode_integration/odeint/odeint.f \
          ../../16_ode_integration/rk4/rk4.f \
          ../../16_ode_integration/rkqc/rkqc.f \
          shootf.f

$(PROG): $(DRIVER) $(DEPS)
	$(FC) $(FFLAGS) -x f77 $(DRIVER) -x none $(DEPS) -o $@

.PHONY: clean run
clean:
	rm -f $(PROG)
run: $(PROG)
	./$(PROG)
