.PHONY: all fortran fortran-lib fortran-programs clean clean-fortran

all: fortran

fortran:
	$(MAKE) -C fortran all

fortran-lib:
	$(MAKE) -C fortran lib

fortran-programs:
	$(MAKE) -C fortran programs

clean: clean-fortran

clean-fortran:
	$(MAKE) -C fortran clean
