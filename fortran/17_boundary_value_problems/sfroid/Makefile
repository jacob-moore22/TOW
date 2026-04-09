FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic

PROG    = sfroid
DRIVER  = sfroid.f

# Subroutine sources (transitive dependencies)
DEPS    = ../../06_special_functions/plgndr/plgndr.f \
          ../bksub/bksub.f \
          ../difeq/difeq.f \
          ../pinvs/pinvs.f \
          ../red/red.f \
          ../solvde/solvde.f

$(PROG): $(DRIVER) $(DEPS)
	$(FC) $(FFLAGS) $(DRIVER) $(DEPS) -o $@

.PHONY: clean run
clean:
	rm -f $(PROG)
run: $(PROG)
	./$(PROG)
