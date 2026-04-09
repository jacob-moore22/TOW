FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic

PROG    = spctrm
DRIVER  = spctrm.dem

# Subroutine sources (transitive dependencies)
DEPS    = ../../12_fourier_transform/four1/four1.f \
          spctrm.f

$(PROG): $(DRIVER) $(DEPS)
	$(FC) $(FFLAGS) -x f77 $(DRIVER) -x none $(DEPS) -o $@

.PHONY: clean run
clean:
	rm -f $(PROG)
run: $(PROG)
	./$(PROG)
