SHELL := /bin/sh
.DEFAULT_GOAL := help
.DELETE_ON_ERROR:

PYTHON ?= python3
IMAGE ?= chimera:local

.PHONY: help install install-dev format format-check lint typecheck test coverage \
	build check-dist smoke security pre-commit container ci clean

.PHONY: check

help: ## Show the available development and release targets.
	@awk 'BEGIN {FS = ":.*## "} /^[a-zA-Z0-9_-]+:.*## / {printf "  %-16s %s\n", $$1, $$2}' $(MAKEFILE_LIST)

install: ## Install CHIMERA into the active Python environment.
	$(PYTHON) -m pip install .

install-dev: ## Install CHIMERA and all development/release dependencies.
	$(PYTHON) -m pip install --editable ".[dev]"

format: ## Format source and tests with Ruff.
	$(PYTHON) -m ruff check --fix .
	$(PYTHON) -m ruff format .

format-check: ## Verify formatting without changing files.
	$(PYTHON) -m ruff format --check .

lint: ## Run the complete configured Ruff rule set.
	$(PYTHON) -m ruff check .

typecheck: ## Run strict static type checking.
	$(PYTHON) -m mypy

test: ## Run the complete test suite.
	$(PYTHON) -m pytest

coverage: ## Run tests with branch coverage and enforce the configured threshold.
	$(PYTHON) -m pytest --cov=chimera --cov-report=term-missing --cov-report=xml

build: ## Build the wheel and source distribution.
	$(PYTHON) -m build

check-dist: build ## Validate built distribution metadata and long descriptions.
	$(PYTHON) -m twine check dist/*

smoke: check-dist ## Install the wheel in a temporary environment and test the CLI.
	@set -eu; \
	smoke_dir="$$(mktemp -d)"; \
	trap 'rm -rf "$$smoke_dir"' EXIT HUP INT TERM; \
	$(PYTHON) -m venv "$$smoke_dir"; \
	set -- dist/*.whl; \
	[ "$$#" -eq 1 ] || { echo "Expected exactly one wheel in dist/, found $$#" >&2; exit 1; }; \
	"$$smoke_dir/bin/python" -m pip install --disable-pip-version-check "$$1"; \
	"$$smoke_dir/bin/python" -m pip check; \
	"$$smoke_dir/bin/chimera" --version; \
	"$$smoke_dir/bin/chimera" --help; \
	"$$smoke_dir/bin/python" -c 'from chimera.schema_resources import validate_packaged_schemas; validate_packaged_schemas()'; \
	"$$smoke_dir/bin/chimera" schema bundle >/dev/null; \
	"$$smoke_dir/bin/chimera" suite --config examples/tiny/chimera.toml --outdir "$$smoke_dir/bundle"; \
	"$$smoke_dir/bin/chimera" validate "$$smoke_dir/bundle"

security: ## Audit installed runtime and development dependencies.
	$(PYTHON) -m pip_audit

pre-commit: ## Run every pre-commit hook against the repository.
	$(PYTHON) -m pre_commit run --all-files

container: ## Build the local non-root OCI image.
	docker build --file containers/Dockerfile --tag $(IMAGE) .

ci: format-check lint typecheck coverage check-dist ## Reproduce stable CI gates locally.

check: ci ## Alias for the complete local CI gate.

clean: ## Remove only generated local build, test, and analysis artifacts.
	rm -rf build dist htmlcov .coverage coverage.xml .pytest_cache .mypy_cache .ruff_cache
