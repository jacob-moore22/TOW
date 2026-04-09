FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic

PROG    = ran2
DRIVER  = ran2.dem

# Subroutine sources (transitive dependencies)
DEPS    = ran2.f

$(PROG): $(DRIVER) $(DEPS)
	$(FC) $(FFLAGS) -x f77 $(DRIVER) -x none $(DEPS) -o $@

.PHONY: clean run
clean:
	rm -f $(PROG)
run: $(PROG)
	./$(PROG)
