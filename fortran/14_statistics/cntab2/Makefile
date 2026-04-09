FC      = gfortran
FFLAGS  = -O2 -std=legacy -w -fno-automatic
LIB     = ../../libnrf77.a

PROG    = cntab2
SRC     = cntab2.dem

$(PROG): $(SRC) $(LIB)
	$(FC) $(FFLAGS) -x f77 $(SRC) -x none $(LIB) -o $@

.PHONY: clean run
clean:
	rm -f $(PROG)
run: $(PROG)
	./$(PROG)
