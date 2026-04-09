FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic

PROG    = balanc
DRIVER  = balanc.dem

# Subroutine sources (transitive dependencies)
DEPS    = balanc.f

$(PROG): $(DRIVER) $(DEPS)
	$(FC) $(FFLAGS) -x f77 $(DRIVER) -x none $(DEPS) -o $@

.PHONY: clean run
clean:
	rm -f $(PROG)
run: $(PROG)
	./$(PROG)
