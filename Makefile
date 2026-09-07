# Compatibility entry point. The research programme is intentionally nested
# under research/; preserve the established top-level make interface.

REPOSITORY_ROOT := $(abspath .)
RESEARCH_PYTHON := $(REPOSITORY_ROOT)/research/.venv/bin/python
LEGACY_PYTHON := $(REPOSITORY_ROOT)/.venv/bin/python
PYTHON ?= $(if $(wildcard $(RESEARCH_PYTHON)),$(RESEARCH_PYTHON),$(LEGACY_PYTHON))
SYSTEM_PYTHON ?= python3

.DEFAULT_GOAL := check
.PHONY: check

check:
	+$(MAKE) -C research PYTHON="$(PYTHON)" SYSTEM_PYTHON="$(SYSTEM_PYTHON)" check

.DEFAULT:
	+$(MAKE) -C research PYTHON="$(PYTHON)" SYSTEM_PYTHON="$(SYSTEM_PYTHON)" "$@"
