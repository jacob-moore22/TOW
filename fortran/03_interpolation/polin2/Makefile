FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic

PROG    = polin2
DRIVER  = polin2.dem

# Subroutine sources (transitive dependencies)
DEPS    = polin2.f \
          ../polint/polint.f

$(PROG): $(DRIVER) $(DEPS)
	$(FC) $(FFLAGS) -x f77 $(DRIVER) -x none $(DEPS) -o $@

.PHONY: clean run
clean:
	rm -f $(PROG)
run: $(PROG)
	./$(PROG)
