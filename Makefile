# From https://github.com/biocore/qurro/blob/master/Makefile

.PHONY: test stylecheck style

# The -B in the invocation of python prevents this from creating pycache stuff.
test:
	python3 -B -m pytest strainflye/tests --cov strainflye

stylecheck:
	flake8 --ignore=W503 strainflye/ misc/ setup.py
	black --check -l 79 strainflye/ misc/ setup.py

style:
	black -l 79 strainflye/ misc/ setup.py
