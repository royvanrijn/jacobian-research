SYSTEM_PYTHON ?= python3
.DEFAULT_GOAL := check

.PHONY: check check-navigation render-status render-navigation
check check-navigation render-status render-navigation:
	+$(MAKE) -C research SYSTEM_PYTHON="$(SYSTEM_PYTHON)" "$@"
