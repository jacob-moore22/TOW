FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic

PROG    = bessk1
DRIVER  = bessk1.dem

# Subroutine sources (transitive dependencies)
DEPS    = ../bessi1/bessi1.f \
          bessk1.f

$(PROG): $(DRIVER) $(DEPS)
	$(FC) $(FFLAGS) -x f77 $(DRIVER) -x none $(DEPS) -o $@

.PHONY: clean run
clean:
	rm -f $(PROG)
run: $(PROG)
	./$(PROG)
