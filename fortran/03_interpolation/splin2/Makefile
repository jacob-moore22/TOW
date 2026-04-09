FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic

PROG    = splin2
DRIVER  = splin2.dem

# Subroutine sources (transitive dependencies)
DEPS    = ../splie2/splie2.f \
          splin2.f \
          ../spline/spline.f \
          ../splint/splint.f

$(PROG): $(DRIVER) $(DEPS)
	$(FC) $(FFLAGS) -x f77 $(DRIVER) -x none $(DEPS) -o $@

.PHONY: clean run
clean:
	rm -f $(PROG)
run: $(PROG)
	./$(PROG)
