FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic

PROG    = qromb
DRIVER  = qromb.dem

# Subroutine sources (transitive dependencies)
DEPS    = ../../03_interpolation/polint/polint.f \
          qromb.f \
          ../trapzd/trapzd.f

$(PROG): $(DRIVER) $(DEPS)
	$(FC) $(FFLAGS) -x f77 $(DRIVER) -x none $(DEPS) -o $@

.PHONY: clean run
clean:
	rm -f $(PROG)
run: $(PROG)
	./$(PROG)
