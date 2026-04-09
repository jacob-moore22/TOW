FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic

PROG    = tred2
DRIVER  = tred2.dem

# Subroutine sources (transitive dependencies)
DEPS    = tred2.f

$(PROG): $(DRIVER) $(DEPS)
	$(FC) $(FFLAGS) -x f77 $(DRIVER) -x none $(DEPS) -o $@

.PHONY: clean run
clean:
	rm -f $(PROG)
run: $(PROG)
	./$(PROG)
