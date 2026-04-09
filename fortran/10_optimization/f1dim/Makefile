FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic

PROG    = f1dim
DRIVER  = f1dim.dem

# Subroutine sources (transitive dependencies)
DEPS    = ../../09_root_finding/scrsho/scrsho.f \
          f1dim.f

$(PROG): $(DRIVER) $(DEPS)
	$(FC) $(FFLAGS) -x f77 $(DRIVER) -x none $(DEPS) -o $@

.PHONY: clean run
clean:
	rm -f $(PROG)
run: $(PROG)
	./$(PROG)
