FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic

PROG    = erfc
DRIVER  = erfc.dem

# Subroutine sources (transitive dependencies)
DEPS    = erfc.f \
          ../gammln/gammln.f \
          ../gammp/gammp.f \
          ../gammq/gammq.f \
          ../gcf/gcf.f \
          ../gser/gser.f

$(PROG): $(DRIVER) $(DEPS)
	$(FC) $(FFLAGS) -x f77 $(DRIVER) -x none $(DEPS) -o $@

.PHONY: clean run
clean:
	rm -f $(PROG)
run: $(PROG)
	./$(PROG)
