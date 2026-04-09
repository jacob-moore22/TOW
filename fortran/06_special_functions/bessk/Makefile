FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic

PROG    = bessk
DRIVER  = bessk.dem

# Subroutine sources (transitive dependencies)
DEPS    = ../bessi0/bessi0.f \
          ../bessi1/bessi1.f \
          bessk.f \
          ../bessk0/bessk0.f \
          ../bessk1/bessk1.f

$(PROG): $(DRIVER) $(DEPS)
	$(FC) $(FFLAGS) -x f77 $(DRIVER) -x none $(DEPS) -o $@

.PHONY: clean run
clean:
	rm -f $(PROG)
run: $(PROG)
	./$(PROG)
