tests:
	python -m pytest -vx 2>&1 | tee errors.err
