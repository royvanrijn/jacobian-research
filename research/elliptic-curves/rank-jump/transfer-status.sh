#!/bin/sh
# Status is read-only. Pass "stop" to stop the owned detached controller.
exec python3 "$(dirname "$0")/run_fresh_constructor_transfer.py" "${1:-status}"
