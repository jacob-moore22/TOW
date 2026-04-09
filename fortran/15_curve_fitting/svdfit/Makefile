FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic

PROG    = svdfit
DRIVER  = svdfit.dem

# Subroutine sources (transitive dependencies)
DEPS    = ../../02_linear_algebra/svbksb/svbksb.f \
          ../../02_linear_algebra/svdcmp/svdcmp.f \
          ../../07_random_numbers/gasdev/gasdev.f \
          ../../07_random_numbers/ran1/ran1.f \
          ../fleg/fleg.f \
          ../fpoly/fpoly.f \
          svdfit.f \
          ../svdvar/svdvar.f

$(PROG): $(DRIVER) $(DEPS)
	$(FC) $(FFLAGS) -x f77 $(DRIVER) -x none $(DEPS) -o $@

.PHONY: clean run
clean:
	rm -f $(PROG)
run: $(PROG)
	./$(PROG)
