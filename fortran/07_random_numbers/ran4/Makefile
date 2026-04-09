FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic

PROG    = ran4
DRIVER  = ran4.dem

# Subroutine sources (transitive dependencies)
DEPS    = ../des/des.f \
          ../desks/desks.f \
          ran4.f \
          ../../14_statistics/avevar/avevar.f

$(PROG): $(DRIVER) $(DEPS)
	$(FC) $(FFLAGS) -x f77 $(DRIVER) -x none $(DEPS) -o $@

.PHONY: clean run
clean:
	rm -f $(PROG)
run: $(PROG)
	./$(PROG)
