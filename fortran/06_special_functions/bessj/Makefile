FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic

PROG    = bessj
DRIVER  = bessj.dem

# Subroutine sources (transitive dependencies)
DEPS    = bessj.f \
          ../bessj0/bessj0.f \
          ../bessj1/bessj1.f

$(PROG): $(DRIVER) $(DEPS)
	$(FC) $(FFLAGS) -x f77 $(DRIVER) -x none $(DEPS) -o $@

.PHONY: clean run
clean:
	rm -f $(PROG)
run: $(PROG)
	./$(PROG)
