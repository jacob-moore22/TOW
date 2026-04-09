FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic

PROG    = des
DRIVER  = des.dem

# Subroutine sources (transitive dependencies)
DEPS    = des.f \
          ../desks/desks.f

$(PROG): $(DRIVER) $(DEPS)
	$(FC) $(FFLAGS) -x f77 $(DRIVER) -x none $(DEPS) -o $@

.PHONY: clean run
clean:
	rm -f $(PROG)
run: $(PROG)
	./$(PROG)
