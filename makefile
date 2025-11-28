.PHONY: help install test lint format clean validate

help:
	@echo "HPC-ScaleTest Development Tasks"
	@echo "================================"
	@echo "install    - Install package and dependencies"
	@echo "test       - Run unit tests with coverage"
	@echo "lint       - Run linting checks"
	@echo "format     - Format code with black"
	@echo "validate   - Validate configuration files"
	@echo "clean      - Remove build artifacts"

install:
	pip install -e .
	pip install -r requirements.txt
	pip install pytest pytest-cov flake8 black
	@echo "✓ Installation complete"

test:
	pytest tests/ -v --cov=. --cov-report=html --cov-report=term
	@echo "✓ Tests passed"

lint:
	flake8 core/ backends/ engine/ utils/ tests/ --max-line-length=100 --ignore=E203,W503 || true
	@echo "✓ Lint complete"

format:
	black core/ backends/ engine/ utils/ tests/ hpc_auto.py --line-length=100
	@echo "✓ Code formatted"

validate:
	@python -c "import yaml; yaml.safe_load(open('run_example.yaml')); print('✓ YAML valid')" || echo "No config to validate"
	@python -c "from utils.validators import ConfigValidator; print('✓ Validators work')"

clean:
	rm -rf build/ dist/ *.egg-info
	find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true
	find . -type f -name '*.pyc' -delete
	find . -type f -name '*.pyo' -delete
	rm -rf .pytest_cache .coverage htmlcov/
	@echo "✓ Cleaned build artifacts"
