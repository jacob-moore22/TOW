FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic

PROG    = smooft
DRIVER  = smooft.dem

# Subroutine sources (transitive dependencies)
DEPS    = ../../07_random_numbers/gasdev/gasdev.f \
          ../../07_random_numbers/ran1/ran1.f \
          ../../12_fourier_transform/four1/four1.f \
          ../../12_fourier_transform/realft/realft.f \
          smooft.f

$(PROG): $(DRIVER) $(DEPS)
	$(FC) $(FFLAGS) -x f77 $(DRIVER) -x none $(DEPS) -o $@

.PHONY: clean run
clean:
	rm -f $(PROG)
run: $(PROG)
	./$(PROG)
