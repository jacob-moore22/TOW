PROG_DIRS := $(sort $(shell find . -mindepth 3 -maxdepth 3 -name Makefile -exec dirname {} \;))

.PHONY: all clean run

all:
	@pass=0; fail=0; \
	for d in $(PROG_DIRS); do \
		if $(MAKE) --no-print-directory -C $$d 2>/dev/null; then \
			pass=$$((pass+1)); \
		else \
			echo "FAIL: $$d"; \
			fail=$$((fail+1)); \
		fi; \
	done; \
	echo ""; \
	echo "$$pass passed, $$fail failed out of $$((pass+fail)) programs"; \
	if [ $$fail -gt 0 ]; then exit 1; fi

run:
	@pass=0; fail=0; \
	for d in $(PROG_DIRS); do \
		if $(MAKE) --no-print-directory -C $$d run 2>/dev/null; then \
			pass=$$((pass+1)); \
		else \
			echo "FAIL: $$d"; \
			fail=$$((fail+1)); \
		fi; \
	done; \
	echo ""; \
	echo "$$pass passed, $$fail failed out of $$((pass+fail)) programs"

clean:
	@for d in $(PROG_DIRS); do $(MAKE) --no-print-directory -C $$d clean 2>/dev/null; done
