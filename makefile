test:
	python -m pytest --cov=dalmisc tests -vx 2>&1 | tee errors.err
