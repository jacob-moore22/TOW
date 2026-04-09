FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic

PROG    = chstwo
DRIVER  = chstwo.dem

# Subroutine sources (transitive dependencies)
DEPS    = ../../06_special_functions/gammln/gammln.f \
          ../../06_special_functions/gammq/gammq.f \
          ../../06_special_functions/gcf/gcf.f \
          ../../06_special_functions/gser/gser.f \
          ../../07_random_numbers/expdev/expdev.f \
          ../../07_random_numbers/ran1/ran1.f \
          chstwo.f

$(PROG): $(DRIVER) $(DEPS)
	$(FC) $(FFLAGS) -x f77 $(DRIVER) -x none $(DEPS) -o $@

.PHONY: clean run
clean:
	rm -f $(PROG)
run: $(PROG)
	./$(PROG)
