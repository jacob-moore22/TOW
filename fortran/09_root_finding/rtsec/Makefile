FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic

PROG    = rtsec
DRIVER  = rtsec.dem

# Subroutine sources (transitive dependencies)
DEPS    = ../../06_special_functions/bessj0/bessj0.f \
          rtsec.f \
          ../zbrak/zbrak.f

$(PROG): $(DRIVER) $(DEPS)
	$(FC) $(FFLAGS) -x f77 $(DRIVER) -x none $(DEPS) -o $@

.PHONY: clean run
clean:
	rm -f $(PROG)
run: $(PROG)
	./$(PROG)
