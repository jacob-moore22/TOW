FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic

PROG    = chebpc
DRIVER  = chebpc.dem

# Subroutine sources (transitive dependencies)
DEPS    = ../chebft/chebft.f \
          chebpc.f

$(PROG): $(DRIVER) $(DEPS)
	$(FC) $(FFLAGS) -x f77 $(DRIVER) -x none $(DEPS) -o $@

.PHONY: clean run
clean:
	rm -f $(PROG)
run: $(PROG)
	./$(PROG)
