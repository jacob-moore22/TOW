FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic

PROG    = twofft
DRIVER  = twofft.dem

# Subroutine sources (transitive dependencies)
DEPS    = ../four1/four1.f \
          twofft.f

$(PROG): $(DRIVER) $(DEPS)
	$(FC) $(FFLAGS) -x f77 $(DRIVER) -x none $(DEPS) -o $@

.PHONY: clean run
clean:
	rm -f $(PROG)
run: $(PROG)
	./$(PROG)
