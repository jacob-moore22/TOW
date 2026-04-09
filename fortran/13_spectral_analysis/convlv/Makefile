FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic

PROG    = convlv
DRIVER  = convlv.dem

# Subroutine sources (transitive dependencies)
DEPS    = ../../12_fourier_transform/four1/four1.f \
          ../../12_fourier_transform/realft/realft.f \
          ../../12_fourier_transform/twofft/twofft.f \
          convlv.f

$(PROG): $(DRIVER) $(DEPS)
	$(FC) $(FFLAGS) -x f77 $(DRIVER) -x none $(DEPS) -o $@

.PHONY: clean run
clean:
	rm -f $(PROG)
run: $(PROG)
	./$(PROG)
