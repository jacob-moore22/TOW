FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic

PROG    = gasdev
DRIVER  = gasdev.dem

# Subroutine sources (transitive dependencies)
DEPS    = gasdev.f \
          ../ran1/ran1.f

$(PROG): $(DRIVER) $(DEPS)
	$(FC) $(FFLAGS) -x f77 $(DRIVER) -x none $(DEPS) -o $@

.PHONY: clean run
clean:
	rm -f $(PROG)
run: $(PROG)
	./$(PROG)
