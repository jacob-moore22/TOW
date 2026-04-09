FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic

PROG    = irbit2
DRIVER  = irbit2.dem

# Subroutine sources (transitive dependencies)
DEPS    = irbit2.f

$(PROG): $(DRIVER) $(DEPS)
	$(FC) $(FFLAGS) -x f77 $(DRIVER) -x none $(DEPS) -o $@

.PHONY: clean run
clean:
	rm -f $(PROG)
run: $(PROG)
	./$(PROG)
