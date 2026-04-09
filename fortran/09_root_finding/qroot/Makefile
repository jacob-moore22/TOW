FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic

PROG    = qroot
DRIVER  = qroot.dem

# Subroutine sources (transitive dependencies)
DEPS    = ../../05_evaluation_of_functions/poldiv/poldiv.f \
          qroot.f

$(PROG): $(DRIVER) $(DEPS)
	$(FC) $(FFLAGS) -x f77 $(DRIVER) -x none $(DEPS) -o $@

.PHONY: clean run
clean:
	rm -f $(PROG)
run: $(PROG)
	./$(PROG)
