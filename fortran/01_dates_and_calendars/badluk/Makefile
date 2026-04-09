FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic

PROG    = badluk
DRIVER  = badluk.f

# Subroutine sources (transitive dependencies)
DEPS    = ../flmoon/flmoon.f \
          ../julday/julday.f

$(PROG): $(DRIVER) $(DEPS)
	$(FC) $(FFLAGS) $(DRIVER) $(DEPS) -o $@

.PHONY: clean run
clean:
	rm -f $(PROG)
run: $(PROG)
	./$(PROG)
