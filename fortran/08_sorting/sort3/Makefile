FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic

PROG    = sort3
DRIVER  = sort3.dem

# Subroutine sources (transitive dependencies)
DEPS    = ../indexx/indexx.f \
          sort3.f

$(PROG): $(DRIVER) $(DEPS)
	$(FC) $(FFLAGS) -x f77 $(DRIVER) -x none $(DEPS) -o $@

.PHONY: clean run
clean:
	rm -f $(PROG)
run: $(PROG)
	./$(PROG)
