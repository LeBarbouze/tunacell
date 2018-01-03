pipinstall:
	pip install -e .

virtualenv:
	virtualenv --system-site-packages venv
	venv/bin/pip install flake8 pytest

flake8:
	venv/bin/flake8 tunacell

test:
	venv/bin/py.test -v tunacell/tests/

simu:
	python bin/tunasimu -s 42

tuto: ~/tmptunacell/simutest
	python scripts/tutorial.py
