.PHONY: help install test lint format clean

help:
	@echo "HPC-ScaleTest Development Tasks"
	@echo "================================"
	@echo "install    - Install package and dependencies"
	@echo "test       - Run unit tests"
	@echo "lint       - Run linting checks"
	@echo "format     - Format code with black"
	@echo "clean      - Remove build artifacts"

install:
	pip install -e .
	pip install -r requirements.txt

test:
	python -m pytest tests/ -v --cov=. --cov-report=html

lint:
	flake8 core/ backends/ engine/ utils/ tests/
	mypy core/ backends/ engine/ utils/

format:
	black core/ backends/ engine/ utils/ tests/ scaletest.py

clean:
	rm -rf build/ dist/ *.egg-info
	find . -type d -name __pycache__ -exec rm -rf {} +
	find . -type f -name '*.pyc' -delete
	find . -type f -name '*.pyo' -delete
	rm -rf .pytest_cache .coverage htmlcov/
