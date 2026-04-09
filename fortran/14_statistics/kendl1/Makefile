FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic

PROG    = kendl1
DRIVER  = kendl1.dem

# Subroutine sources (transitive dependencies)
DEPS    = ../../06_special_functions/erfcc/erfcc.f \
          ../../07_random_numbers/des/des.f \
          ../../07_random_numbers/desks/desks.f \
          ../../07_random_numbers/ran0/ran0.f \
          ../../07_random_numbers/ran1/ran1.f \
          ../../07_random_numbers/ran2/ran2.f \
          ../../07_random_numbers/ran3/ran3.f \
          ../../07_random_numbers/ran4/ran4.f \
          kendl1.f

$(PROG): $(DRIVER) $(DEPS)
	$(FC) $(FFLAGS) -x f77 $(DRIVER) -x none $(DEPS) -o $@

.PHONY: clean run
clean:
	rm -f $(PROG)
run: $(PROG)
	./$(PROG)
