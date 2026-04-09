FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic

PROG    = kendl2
DRIVER  = kendl2.dem

# Subroutine sources (transitive dependencies)
DEPS    = ../../06_special_functions/erfcc/erfcc.f \
          ../../07_random_numbers/irbit1/irbit1.f \
          ../../07_random_numbers/irbit2/irbit2.f \
          kendl2.f

$(PROG): $(DRIVER) $(DEPS)
	$(FC) $(FFLAGS) -x f77 $(DRIVER) -x none $(DEPS) -o $@

.PHONY: clean run
clean:
	rm -f $(PROG)
run: $(PROG)
	./$(PROG)
