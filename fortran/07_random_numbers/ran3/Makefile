FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic

PROG    = ran3
DRIVER  = ran3.dem

# Subroutine sources (transitive dependencies)
DEPS    = ran3.f

$(PROG): $(DRIVER) $(DEPS)
	$(FC) $(FFLAGS) -x f77 $(DRIVER) -x none $(DEPS) -o $@

.PHONY: clean run
clean:
	rm -f $(PROG)
run: $(PROG)
	./$(PROG)
