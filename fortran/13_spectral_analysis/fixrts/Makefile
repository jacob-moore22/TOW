FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic

PROG    = fixrts
DRIVER  = fixrts.dem

# Subroutine sources (transitive dependencies)
DEPS    = ../../09_root_finding/laguer/laguer.f \
          ../../09_root_finding/zroots/zroots.f \
          fixrts.f

$(PROG): $(DRIVER) $(DEPS)
	$(FC) $(FFLAGS) -x f77 $(DRIVER) -x none $(DEPS) -o $@

.PHONY: clean run
clean:
	rm -f $(PROG)
run: $(PROG)
	./$(PROG)
