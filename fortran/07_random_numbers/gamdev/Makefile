FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic

PROG    = gamdev
DRIVER  = gamdev.dem

# Subroutine sources (transitive dependencies)
DEPS    = gamdev.f \
          ../ran1/ran1.f

$(PROG): $(DRIVER) $(DEPS)
	$(FC) $(FFLAGS) -x f77 $(DRIVER) -x none $(DEPS) -o $@

.PHONY: clean run
clean:
	rm -f $(PROG)
run: $(PROG)
	./$(PROG)
