FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic

PROG    = quad3d
DRIVER  = quad3d.dem

# Subroutine sources (transitive dependencies)
DEPS    = quad3d.f

$(PROG): $(DRIVER) $(DEPS)
	$(FC) $(FFLAGS) -x f77 $(DRIVER) -x none $(DEPS) -o $@

.PHONY: clean run
clean:
	rm -f $(PROG)
run: $(PROG)
	./$(PROG)
