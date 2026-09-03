PYTHON ?= python3

.PHONY: verify test clean

# run every verification script and compare against the expected values
verify:
	$(PYTHON) run_all.py

# unit tests of the exact-arithmetic checker only
test:
	$(PYTHON) -m pytest verify/ -q

clean:
	rm -rf results .pytest_cache
	find . -name __pycache__ -type d -prune -exec rm -rf {} +
