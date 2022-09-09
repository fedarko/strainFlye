# From https://github.com/biocore/qurro/blob/master/Makefile
# (and https://github.com/fedarko/pleuk/blob/master/Makefile, but that one
# isn't public on GitHub yet)

.PHONY: test stylecheck style

# The -B in the invocation of python prevents this from creating pycache stuff.
# --cov-report xml is needed to make this visible to Codecov, and
# --cov-report term is needed in order to print a fancy table on the terminal
# (useful when testing locally).
#
# Also, for reference -- if you wanna see all print output during testing, add
# -s to the end of this command. Useful when debugging stuff.
test:
	python3 -B -m pytest --doctest-modules strainflye/ --cov-report xml --cov-report term --cov-report html --cov strainflye

stylecheck:
	@# Ignoring E203 per https://github.com/psf/black/issues/315
	flake8 --ignore=W503,E203 strainflye/ setup.py
	black --check -l 79 strainflye/ setup.py

style:
	black -l 79 strainflye/ setup.py
