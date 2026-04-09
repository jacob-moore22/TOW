FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic

PROG    = lubksb
DRIVER  = lubksb.dem

# Subroutine sources (transitive dependencies)
DEPS    = lubksb.f \
          ../ludcmp/ludcmp.f

$(PROG): $(DRIVER) $(DEPS)
	$(FC) $(FFLAGS) -x f77 $(DRIVER) -x none $(DEPS) -o $@

.PHONY: clean run
clean:
	rm -f $(PROG)
run: $(PROG)
	./$(PROG)
