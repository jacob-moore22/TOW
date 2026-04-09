FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic

PROG    = cosft
DRIVER  = cosft.dem

# Subroutine sources (transitive dependencies)
DEPS    = cosft.f \
          ../four1/four1.f \
          ../realft/realft.f

$(PROG): $(DRIVER) $(DEPS)
	$(FC) $(FFLAGS) -x f77 $(DRIVER) -x none $(DEPS) -o $@

.PHONY: clean run
clean:
	rm -f $(PROG)
run: $(PROG)
	./$(PROG)
