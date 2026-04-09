FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic

PROG    = chebev
DRIVER  = chebev.dem

# Subroutine sources (transitive dependencies)
DEPS    = chebev.f \
          ../chebft/chebft.f

$(PROG): $(DRIVER) $(DEPS)
	$(FC) $(FFLAGS) -x f77 $(DRIVER) -x none $(DEPS) -o $@

.PHONY: clean run
clean:
	rm -f $(PROG)
run: $(PROG)
	./$(PROG)
