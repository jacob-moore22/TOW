FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic

PROG    = covsrt
DRIVER  = covsrt.dem

# Subroutine sources (transitive dependencies)
DEPS    = covsrt.f

$(PROG): $(DRIVER) $(DEPS)
	$(FC) $(FFLAGS) -x f77 $(DRIVER) -x none $(DEPS) -o $@

.PHONY: clean run
clean:
	rm -f $(PROG)
run: $(PROG)
	./$(PROG)
