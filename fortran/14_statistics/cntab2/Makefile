FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic

PROG    = cntab2
DRIVER  = cntab2.dem

# Subroutine sources (transitive dependencies)
DEPS    = cntab2.f

$(PROG): $(DRIVER) $(DEPS)
	$(FC) $(FFLAGS) -x f77 $(DRIVER) -x none $(DEPS) -o $@

.PHONY: clean run
clean:
	rm -f $(PROG)
run: $(PROG)
	./$(PROG)
