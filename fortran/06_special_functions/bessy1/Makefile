FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic

PROG    = bessy1
DRIVER  = bessy1.dem

# Subroutine sources (transitive dependencies)
DEPS    = ../bessj1/bessj1.f \
          bessy1.f

$(PROG): $(DRIVER) $(DEPS)
	$(FC) $(FFLAGS) -x f77 $(DRIVER) -x none $(DEPS) -o $@

.PHONY: clean run
clean:
	rm -f $(PROG)
run: $(PROG)
	./$(PROG)
