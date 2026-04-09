FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic

PROG    = eigsrt
DRIVER  = eigsrt.dem

# Subroutine sources (transitive dependencies)
DEPS    = eigsrt.f \
          ../jacobi/jacobi.f

$(PROG): $(DRIVER) $(DEPS)
	$(FC) $(FFLAGS) -x f77 $(DRIVER) -x none $(DEPS) -o $@

.PHONY: clean run
clean:
	rm -f $(PROG)
run: $(PROG)
	./$(PROG)
