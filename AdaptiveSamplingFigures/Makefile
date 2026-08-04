PYTHON ?= python3
MPLBACKEND ?= Agg
export MPLBACKEND

.PHONY: figures verify verify-strict

figures:
	$(PYTHON) src/generate_all.py

verify:
	$(PYTHON) src/verify_outputs.py

verify-strict:
	$(PYTHON) src/verify_outputs.py --strict
