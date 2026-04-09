FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic

PROG    = qtrap
DRIVER  = qtrap.dem

# Subroutine sources (transitive dependencies)
DEPS    = qtrap.f \
          ../trapzd/trapzd.f

$(PROG): $(DRIVER) $(DEPS)
	$(FC) $(FFLAGS) -x f77 $(DRIVER) -x none $(DEPS) -o $@

.PHONY: clean run
clean:
	rm -f $(PROG)
run: $(PROG)
	./$(PROG)
