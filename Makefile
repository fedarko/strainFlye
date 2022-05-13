# From https://github.com/biocore/qurro/blob/master/Makefile
# (and https://github.com/fedarko/pleuk/blob/master/Makefile, but that one
# isn't public on GitHub yet)

.PHONY: test stylecheck style

# The -B in the invocation of python prevents this from creating pycache stuff.
# The --cov-report xml seems to be needed to make this visible to Codecov.
test:
	python3 -B -m pytest strainflye/tests --cov strainflye -s

stylecheck:
	flake8 --ignore=W503 strainflye/ misc/ setup.py
	black --check -l 79 strainflye/ misc/ setup.py

style:
	black -l 79 strainflye/ misc/ setup.py
