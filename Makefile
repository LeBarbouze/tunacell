pipinstall:
	pip install -e .

virtualenv:
	virtualenv venv
	venv/bin/pip install flake8 pytest

flake8:
	venv/bin/flake8 tunacell

test:
	venv/bin/py.test -v tunacell/tests/

simu:
	python bin/tunasimu -s 42

tuto: ~/tmptunacell/simutest
	python scripts/tutorial.py -i

plotting-demo:
	python scripts/plotting-samples.py -i

analysis-demo:
	python scripts/univariate-analysis.py -i
	python scripts/univariate-analysis-2.py -i
	python scripts/bivariate-analysis.py -i

full-demo:
	python bin/tunasimu -s 42 -f
	python scripts/tutorial.py --time .5  # waiting time minimal for efficacy
	python scripts/plotting-samples.py --time .5
	python scripts/univariate-analysis.py --time .5
	python scripts/univariate-analysis-2.py --time .5
	python scripts/bivariate-analysis.py --time .5
