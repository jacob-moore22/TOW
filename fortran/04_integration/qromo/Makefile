FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic

PROG    = qromo
DRIVER  = qromo.dem

# Subroutine sources (transitive dependencies)
DEPS    = ../../03_interpolation/polint/polint.f \
          ../midinf/midinf.f \
          qromo.f

$(PROG): $(DRIVER) $(DEPS)
	$(FC) $(FFLAGS) -x f77 $(DRIVER) -x none $(DEPS) -o $@

.PHONY: clean run
clean:
	rm -f $(PROG)
run: $(PROG)
	./$(PROG)
