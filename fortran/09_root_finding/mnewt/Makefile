FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic

PROG    = mnewt
DRIVER  = mnewt.dem

# Subroutine sources (transitive dependencies)
DEPS    = ../../02_linear_algebra/lubksb/lubksb.f \
          ../../02_linear_algebra/ludcmp/ludcmp.f \
          ../../06_special_functions/beta/beta.f \
          ../../06_special_functions/gammln/gammln.f \
          mnewt.f

$(PROG): $(DRIVER) $(DEPS)
	$(FC) $(FFLAGS) -x f77 $(DRIVER) -x none $(DEPS) -o $@

.PHONY: clean run
clean:
	rm -f $(PROG)
run: $(PROG)
	./$(PROG)
