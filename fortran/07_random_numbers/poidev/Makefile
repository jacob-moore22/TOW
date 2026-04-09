FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic
LIB     = ../../libnrf77.a

PROG    = poidev
SRC     = poidev.dem

$(PROG): $(SRC) $(LIB)
	$(FC) $(FFLAGS) -x f77 $(SRC) -x none $(LIB) -o $@

.PHONY: clean run
clean:
	rm -f $(PROG)
run: $(PROG)
	./$(PROG)
