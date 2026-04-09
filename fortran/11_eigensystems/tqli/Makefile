FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic

PROG    = tqli
DRIVER  = tqli.dem

# Subroutine sources (transitive dependencies)
DEPS    = tqli.f \
          ../tred2/tred2.f

$(PROG): $(DRIVER) $(DEPS)
	$(FC) $(FFLAGS) -x f77 $(DRIVER) -x none $(DEPS) -o $@

.PHONY: clean run
clean:
	rm -f $(PROG)
run: $(PROG)
	./$(PROG)
