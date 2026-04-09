FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic

PROG    = vander
DRIVER  = vander.dem

# Subroutine sources (transitive dependencies)
DEPS    = vander.f

$(PROG): $(DRIVER) $(DEPS)
	$(FC) $(FFLAGS) -x f77 $(DRIVER) -x none $(DEPS) -o $@

.PHONY: clean run
clean:
	rm -f $(PROG)
run: $(PROG)
	./$(PROG)
