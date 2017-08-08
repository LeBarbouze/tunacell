pipinstall:
	pip install -e .

virtualenv:
	virtualenv --python=python2.7 venv
	venv/bin/pip install flake8 pytest

flake8:
	venv/bin/flake8 tuna

test:
	venv/bin/py.test -v tuna/tests/