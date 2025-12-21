.PHONY: help install test lint format clean validate check examples

help:
	@echo "HPC-ScaleTest Development Tasks"
	@echo "================================"
	@echo "install    - Install package and dependencies"
	@echo "test       - Run unit tests with coverage"
	@echo "lint       - Run linting checks"
	@echo "format     - Format code with black"
	@echo "validate   - Validate configuration files"
	@echo "check      - Run system check utility"
	@echo "examples   - Validate example configurations"
	@echo "clean      - Remove build artifacts"
	@echo "clean-all  - Remove all generated files including outputs"

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
	black core/ backends/ engine/ utils/ tests/ hpc_auto.py submit_jobs.py check_system.py --line-length=100
	@echo "✓ Code formatted"

validate:
	@python -c "import yaml; yaml.safe_load(open('run.yaml')); print('✓ run.yaml valid')" || echo "○ No run.yaml to validate"
	@python -c "import yaml; yaml.safe_load(open('run.generic.yaml')); print('✓ run.generic.yaml valid')"
	@python -c "from utils.validators import ConfigValidator; print('✓ Validators work')"
	@echo "✓ Validation complete"

check:
	@python3 check_system.py
	@echo ""

examples:
	@echo "Validating example configurations..."
	@python -c "import yaml; yaml.safe_load(open('examples/example_2d_weak_scaling.yaml')); print('✓ 2D weak scaling example')"
	@python -c "import yaml; yaml.safe_load(open('examples/example_3d_weak_scaling.yaml')); print('✓ 3D weak scaling example')"
	@python -c "import yaml; yaml.safe_load(open('examples/example_strong_scaling.yaml')); print('✓ Strong scaling example')"
	@python -c "import yaml; yaml.safe_load(open('examples/run.template.yaml')); print('✓ Run template')"
	@echo "✓ All example configurations valid"

clean:
	rm -rf build/ dist/ *.egg-info
	find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true
	find . -type f -name '*.pyc' -delete
	find . -type f -name '*.pyo' -delete
	rm -rf .pytest_cache .coverage htmlcov/
	@echo "✓ Cleaned build artifacts"

clean-all: clean
	rm -rf output/ results/ logs/
	rm -f *.out *.err slurm-*.out job_ids.txt
	rm -rf cloned_repos/ source/
	@echo "✓ Cleaned all generated files"
