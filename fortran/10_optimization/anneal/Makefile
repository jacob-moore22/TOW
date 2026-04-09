FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic

PROG    = anneal
DRIVER  = anneal.dem

# Subroutine sources (transitive dependencies)
DEPS    = ../../07_random_numbers/irbit1/irbit1.f \
          ../../07_random_numbers/ran3/ran3.f \
          anneal.f \
          ../link/link.f

$(PROG): $(DRIVER) $(DEPS)
	$(FC) $(FFLAGS) -x f77 $(DRIVER) -x none $(DEPS) -o $@

.PHONY: clean run
clean:
	rm -f $(PROG)
run: $(PROG)
	./$(PROG)
