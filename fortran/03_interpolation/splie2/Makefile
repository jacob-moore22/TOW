FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic

PROG    = splie2
DRIVER  = splie2.dem

# Subroutine sources (transitive dependencies)
DEPS    = splie2.f \
          ../spline/spline.f

$(PROG): $(DRIVER) $(DEPS)
	$(FC) $(FFLAGS) -x f77 $(DRIVER) -x none $(DEPS) -o $@

.PHONY: clean run
clean:
	rm -f $(PROG)
run: $(PROG)
	./$(PROG)
