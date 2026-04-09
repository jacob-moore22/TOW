FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic

PROG    = hqr
DRIVER  = hqr.dem

# Subroutine sources (transitive dependencies)
DEPS    = ../balanc/balanc.f \
          ../elmhes/elmhes.f \
          hqr.f

$(PROG): $(DRIVER) $(DEPS)
	$(FC) $(FFLAGS) -x f77 $(DRIVER) -x none $(DEPS) -o $@

.PHONY: clean run
clean:
	rm -f $(PROG)
run: $(PROG)
	./$(PROG)
