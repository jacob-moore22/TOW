FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic

PROG    = mprove
DRIVER  = mprove.dem

# Subroutine sources (transitive dependencies)
DEPS    = ../lubksb/lubksb.f \
          ../ludcmp/ludcmp.f \
          mprove.f \
          ../../07_random_numbers/ran3/ran3.f

$(PROG): $(DRIVER) $(DEPS)
	$(FC) $(FFLAGS) -x f77 $(DRIVER) -x none $(DEPS) -o $@

.PHONY: clean run
clean:
	rm -f $(PROG)
run: $(PROG)
	./$(PROG)
