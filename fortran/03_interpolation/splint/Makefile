FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic

PROG    = splint
DRIVER  = splint.dem

# Subroutine sources (transitive dependencies)
DEPS    = ../spline/spline.f \
          splint.f

$(PROG): $(DRIVER) $(DEPS)
	$(FC) $(FFLAGS) -x f77 $(DRIVER) -x none $(DEPS) -o $@

.PHONY: clean run
clean:
	rm -f $(PROG)
run: $(PROG)
	./$(PROG)
