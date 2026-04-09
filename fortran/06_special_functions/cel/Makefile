FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic

PROG    = cel
DRIVER  = cel.dem

# Subroutine sources (transitive dependencies)
DEPS    = ../../04_integration/qsimp/qsimp.f \
          ../../04_integration/trapzd/trapzd.f \
          cel.f \
          ../../07_random_numbers/ran3/ran3.f

$(PROG): $(DRIVER) $(DEPS)
	$(FC) $(FFLAGS) -x f77 $(DRIVER) -x none $(DEPS) -o $@

.PHONY: clean run
clean:
	rm -f $(PROG)
run: $(PROG)
	./$(PROG)
