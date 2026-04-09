FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic

PROG    = eulsum
DRIVER  = eulsum.dem

# Subroutine sources (transitive dependencies)
DEPS    = eulsum.f

$(PROG): $(DRIVER) $(DEPS)
	$(FC) $(FFLAGS) -x f77 $(DRIVER) -x none $(DEPS) -o $@

.PHONY: clean run
clean:
	rm -f $(PROG)
run: $(PROG)
	./$(PROG)
