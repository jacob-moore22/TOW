FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic

PROG    = crank
DRIVER  = crank.dem

# Subroutine sources (transitive dependencies)
DEPS    = ../../08_sorting/sort2/sort2.f \
          crank.f

$(PROG): $(DRIVER) $(DEPS)
	$(FC) $(FFLAGS) -x f77 $(DRIVER) -x none $(DEPS) -o $@

.PHONY: clean run
clean:
	rm -f $(PROG)
run: $(PROG)
	./$(PROG)
