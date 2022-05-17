# From https://github.com/biocore/qurro/blob/master/Makefile
# (and https://github.com/fedarko/pleuk/blob/master/Makefile, but that one
# isn't public on GitHub yet)

.PHONY: test stylecheck style

# The -B in the invocation of python prevents this from creating pycache stuff.
# --cov-report xml is needed to make this visible to Codecov, and
# --cov-report term is needed in order to print a fancy table on the terminal
# (useful when testing locally).
test:
	python3 -B -m pytest strainflye/tests --cov-report xml --cov-report term --cov strainflye

stylecheck:
	flake8 --ignore=W503 strainflye/ misc/ setup.py
	black --check -l 79 strainflye/ misc/ setup.py

style:
	black -l 79 strainflye/ misc/ setup.py
