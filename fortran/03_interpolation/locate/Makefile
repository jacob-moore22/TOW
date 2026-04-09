FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic

PROG    = locate
DRIVER  = locate.dem

# Subroutine sources (transitive dependencies)
DEPS    = locate.f

$(PROG): $(DRIVER) $(DEPS)
	$(FC) $(FFLAGS) -x f77 $(DRIVER) -x none $(DEPS) -o $@

.PHONY: clean run
clean:
	rm -f $(PROG)
run: $(PROG)
	./$(PROG)
