"""Explicit runtime amendment of frozen historical half-lattice protocols.

Candidate lists, centre-selection rules, and budgets retain their own checks.
Only listed implementation replacements may use the pinned historical source
hash. New records bind the current sources and coordinates separately, and
must not resume or merge historical results under the amended search policy.
"""
from functools import lru_cache
from hashlib import sha256
from pathlib import Path
import subprocess

import pointed_quartic_search as pointed

ROOT = pointed.ROOT
REGRESSION_REVISION = "d30a742133f0658185c3bd4c99f0b0f815f2f74b"
MIGRATED = (
    "elliptic-curves/cas/run_curve385_iterated_half_lattice_search.sage",
    "elliptic-curves/cas/run_r17_refresh_jump_ladder_blind.sage",
    "elliptic-curves/cas/run_r17_refresh_jump_ladder_blind_v2.sage",
    "elliptic-curves/cas/run_mw17_jump_v2.sage",
    "elliptic-curves/cas/run_mw17_jump_v2_zero_gain_rescue.sage",
    "elliptic-curves/cas/run_curve398_mw16_adaptive_half_lattice_search.sage",
    "elliptic-curves/cas/run_icarm_mw16_parent_ladder_blind.sage",
    "elliptic-curves/cas/run_icarm_mw16_curve400_adaptive_calibration.sage",
    "elliptic-curves/cas/run_curve385_sparse_quotient_rank32_search.sage",
    "elliptic-curves/cas/run_curve385_sparse_quotient_rank32_search_v2.sage",
    "elliptic-curves/cas/run_a1_mw16_blind_parameter_experiment.sage",
    "elliptic-curves/cas/run_a1_mw16_target_free_parameter_search.sage",
    "elliptic-curves/cas/run_icarm_mw16_nagao_finalist_half_lattice.sage",
    "elliptic-curves/cas/run_mw16_sensitivity.sage",
    "elkies-k3/scripts/run_r17_prospective_crt_half_lattice_search.sage",
)


def validate_frozen_sources(expected):
    for name, digest in expected.items():
        if sha256((ROOT/name).read_bytes()).hexdigest() == digest:
            continue
        if name not in MIGRATED:
            raise ArithmeticError("frozen source changed outside the backend amendment: "+name)
        original = subprocess.check_output(["git", "show", REGRESSION_REVISION+":"+name], cwd=ROOT)
        if sha256(original).hexdigest() != digest:
            raise ArithmeticError("historical protocol source does not match the pinned regression revision: "+name)


@lru_cache(maxsize=1)
def runtime_search():
    paths = (*MIGRATED, str(Path(__file__).resolve().relative_to(ROOT)))
    return {"backend": pointed.BACKEND_NAME, "coordinate_policy": pointed.CoordinatePolicy().record(),
            "regression_revision": REGRESSION_REVISION,
            "source_hashes": {**pointed.sources(), **{p: sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}},
            "amendment": "Use PointedQuarticSearch instead of the historical search implementation; finite coordinate boxes change. Candidate and centre selection and rank certification remain governed by the source protocol."}


def require_runtime(record):
    if record.get("runtime_search") != runtime_search():
        raise ArithmeticError("checkpoint belongs to another backend/coordinate policy or source revision; use a fresh output path")
