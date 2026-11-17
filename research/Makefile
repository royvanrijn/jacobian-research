SHELL := /bin/bash
SYSTEM_PYTHON ?= python3
.DEFAULT_GOAL := check

.PHONY: check check-navigation render-status render-navigation
check: check-navigation
	$(SYSTEM_PYTHON) -m compileall -q scripts

check-navigation:
	$(SYSTEM_PYTHON) scripts/research.py check
	$(SYSTEM_PYTHON) -m unittest discover -s tests -q

render-status render-navigation:
	$(SYSTEM_PYTHON) scripts/research.py render
