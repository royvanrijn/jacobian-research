SHELL := /bin/bash
PYTHON ?= .venv/bin/python
SYSTEM_PYTHON ?= python3
.DEFAULT_GOAL := check

.PHONY: check render-status verify-elliptic-curves check-singular-runner
check: check-navigation
	$(SYSTEM_PYTHON) -m compileall -q scripts elliptic-curves/cas elliptic-curves/ecsearch elliptic-curves/scripts elliptic-curves/tests

render-status:
	$(SYSTEM_PYTHON) scripts/render_status.py

# Metadata and documentation only: no CAS, certificate replay or inventory rebuild.
.PHONY: render-navigation check-navigation
render-navigation:
	$(SYSTEM_PYTHON) scripts/render_status.py
	$(SYSTEM_PYTHON) elliptic-curves/cas/render_main_readme_curves.py --from-inventory
	$(SYSTEM_PYTHON) scripts/research.py render

check-navigation:
	$(SYSTEM_PYTHON) scripts/audit_status.py
	$(SYSTEM_PYTHON) scripts/research.py check --require-partial-review
	$(SYSTEM_PYTHON) scripts/check_markdown_links.py
	$(SYSTEM_PYTHON) scripts/audit_repository_hygiene.py
	$(SYSTEM_PYTHON) scripts/audit_cleanup_snapshot.py
	$(SYSTEM_PYTHON) scripts/audit_programme_archive.py
	$(SYSTEM_PYTHON) elkies-k3/scripts/analyze_process_ledger.py --check-document
	$(SYSTEM_PYTHON) elliptic-curves/cas/render_main_readme_curves.py --check --from-inventory
	$(SYSTEM_PYTHON) -m unittest discover -s tests -p 'test_research*.py' -q

# Explicit finite arithmetic witness, independent of the original Sage producer.
.PHONY: verify-ns0031-arithmetic
verify-ns0031-arithmetic:
	$(SYSTEM_PYTHON) elkies-k3/scripts/verify_ns0031_marking_arithmetic.py
	$(SYSTEM_PYTHON) -m unittest discover -s tests -p 'test_ns0031_marking_replay.py' -q

# Explicit broad EC suite. Inspect its cost and prerequisites before invoking.
verify-elliptic-curves:
	$(PYTHON) elliptic-curves/scripts/audit_artifact_catalog.py
	PYTHONPATH=elliptic-curves:elliptic-curves/cas $(PYTHON) -m unittest discover -s elliptic-curves/tests -v
	$(PYTHON) elliptic-curves/scripts/verify_family_data.py
	$(PYTHON) elliptic-curves/scripts/verify_benchmarks.py
	$(PYTHON) elliptic-curves/scripts/verify_fermigier_rank_certificates.py
	$(PYTHON) elliptic-curves/scripts/verify_fermigier_rank20_near_miss.py
	$(PYTHON) elliptic-curves/scripts/verify_crt_lattice_calibration.py
	$(PYTHON) elliptic-curves/scripts/verify_fermigier_crt_seed.py
	$(PYTHON) elliptic-curves/scripts/verify_kihara_rank14.py
	$(PYTHON) elliptic-curves/scripts/verify_e29_independence.py
	$(SYSTEM_PYTHON) elliptic-curves/scripts/verify_k3_chain_ledger.py
	PYTHONPATH=elliptic-curves/cas $(PYTHON) elliptic-curves/cas/build_structural_search_groundwork.py --check
	$(PYTHON) elliptic-curves/cas/verify_icarm_curve273_rank30.py --check
	PYTHONPATH=elliptic-curves/cas $(PYTHON) elliptic-curves/cas/check_icarm_curve302_rank31_pinned.py
	$(PYTHON) elliptic-curves/cas/analyze_icarm_7fff_zip_sequence.py --check

check-singular-runner:
	$(SYSTEM_PYTHON) -m unittest discover -s tests -p 'test_singular_runner.py' -v
