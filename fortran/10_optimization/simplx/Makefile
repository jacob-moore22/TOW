FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic

PROG    = simplx
DRIVER  = simplx.dem

# Subroutine sources (transitive dependencies)
DEPS    = ../simp1/simp1.f \
          ../simp2/simp2.f \
          ../simp3/simp3.f \
          simplx.f

$(PROG): $(DRIVER) $(DEPS)
	$(FC) $(FFLAGS) -x f77 $(DRIVER) -x none $(DEPS) -o $@

.PHONY: clean run
clean:
	rm -f $(PROG)
run: $(PROG)
	./$(PROG)
