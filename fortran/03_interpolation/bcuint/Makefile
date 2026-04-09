FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic

PROG    = bcuint
DRIVER  = bcuint.dem

# Subroutine sources (transitive dependencies)
DEPS    = ../bcucof/bcucof.f \
          bcuint.f

$(PROG): $(DRIVER) $(DEPS)
	$(FC) $(FFLAGS) -x f77 $(DRIVER) -x none $(DEPS) -o $@

.PHONY: clean run
clean:
	rm -f $(PROG)
run: $(PROG)
	./$(PROG)
