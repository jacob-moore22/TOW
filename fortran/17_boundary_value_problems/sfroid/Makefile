FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic
LIB     = ../../libnrf77.a

PROG    = sfroid
SRC     = sfroid.f

$(PROG): $(SRC) $(LIB)
	$(FC) $(FFLAGS) -o $@ $(SRC) $(LIB)

.PHONY: clean run
clean:
	rm -f $(PROG)
run: $(PROG)
	./$(PROG)
