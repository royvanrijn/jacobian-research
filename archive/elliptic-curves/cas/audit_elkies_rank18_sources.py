#!/usr/bin/env python3
"""Audit whether Elkies's published rank-18 family is reproducible.

Only official primary-source bundles are inspected: arXiv:0709.2908 and the
related computational paper arXiv:0802.1301.  The output distinguishes the
published existence/construction claim from the explicit coefficient and
section data needed to instantiate and screen specializations.
"""

from __future__ import annotations

import argparse
from datetime import date
import gzip
from hashlib import sha256
from io import BytesIO
import json
from pathlib import Path
import platform
import tarfile
from typing import Any
from urllib.request import Request, urlopen

import sympy as sp

from kihara_rank14 import (
    PRIMARY_SOURCE_DOI as KIHARA_PRIMARY_SOURCE_DOI,
    PRIMARY_SOURCE_PDF as KIHARA_PRIMARY_SOURCE_PDF,
    PUBLISHED_HEIGHT_DETERMINANT_APPROX,
    PUBLISHED_INDEPENDENCE_SPECIALIZATION,
    PUBLISHED_RANK_LOWER_BOUND,
    short_jacobian_coefficients,
    verify_rational_specialization,
)
from verify_kihara_rank14_symbolic import symbolic_verification


ELKIES_0709_ABSTRACT_URL = "https://arxiv.org/abs/0709.2908"
ELKIES_0709_SOURCE_URL = "https://export.arxiv.org/e-print/0709.2908"
ELKIES_0802_ABSTRACT_URL = "https://arxiv.org/abs/0802.1301"
ELKIES_0802_SOURCE_URL = "https://export.arxiv.org/e-print/0802.1301"

# These pin the exact official source revisions inspected on 2026-08-14.
ELKIES_0709_SOURCE_SHA256 = (
    "daa77bdb1da8ea01ebdd7e7add7ed37af5f7057c90e7221b7ddac39fd65f785c"
)
ELKIES_0802_SOURCE_SHA256 = (
    "b8a9cab6c5c723c12d56b11a2060a0e9462bd1c4b5e8742f2cfcc2a7f77fc994"
)

REPRODUCE_COMMAND = (
    ".venv/bin/python elliptic-curves/cas/audit_elkies_rank18_sources.py "
    "--output artifacts/generated-results/elliptic_elkies_rank18_source_audit.json"
)


def _download(url: str, timeout: float) -> bytes:
    request = Request(
        url,
        headers={
            "User-Agent": "elliptic-rank-source-audit/1.0",
            "Accept": "application/gzip, application/x-tar, text/plain, */*",
        },
    )
    with urlopen(request, timeout=timeout) as response:
        return response.read()


def _regular_tar_members(payload: bytes) -> dict[str, bytes]:
    answer: dict[str, bytes] = {}
    with tarfile.open(fileobj=BytesIO(payload), mode="r:") as archive:
        for member in archive.getmembers():
            if not member.isfile():
                continue
            extracted = archive.extractfile(member)
            if extracted is None:
                raise ValueError(f"could not read tar member {member.name!r}")
            answer[member.name] = extracted.read()
    return answer


def unpack_arxiv_source(
    source_bytes: bytes, *, single_tex_name: str
) -> dict[str, bytes]:
    """Unpack an arXiv gzip source, whether tarred or a single TeX file."""

    decompressed = gzip.decompress(source_bytes)
    try:
        members = _regular_tar_members(decompressed)
    except tarfile.ReadError:
        members = {single_tex_name: decompressed}
    if not members:
        raise ValueError("the arXiv source bundle contained no regular files")
    return members


def _source_inventory(
    source_bytes: bytes,
    *,
    single_tex_name: str,
    expected_sha256: str,
) -> tuple[dict[str, Any], str]:
    members = unpack_arxiv_source(source_bytes, single_tex_name=single_tex_name)
    tex_members = {
        name: payload
        for name, payload in members.items()
        if name.lower().endswith(".tex")
    }
    joined_tex = "\n".join(
        payload.decode("utf-8", errors="replace")
        for _, payload in sorted(tex_members.items())
    )
    source_hash = sha256(source_bytes).hexdigest()
    inventory = {
        "compressed_bytes": len(source_bytes),
        "source_sha256": source_hash,
        "matches_pinned_revision": source_hash == expected_sha256,
        "members": [
            {
                "name": name,
                "bytes": len(payload),
                "sha256": sha256(payload).hexdigest(),
            }
            for name, payload in sorted(members.items())
        ],
        "tex_member_count": len(tex_members),
        "non_tex_data_members": [
            name
            for name in sorted(members)
            if not name.lower().endswith((".tex", ".cls", ".sty", ".bst"))
        ],
        "http_token_count_in_tex": joined_tex.lower().count("http"),
    }
    return inventory, joined_tex


def audit_official_sources(
    source_0709: bytes, source_0802: bytes
) -> dict[str, Any]:
    inventory_0709, tex_0709 = _source_inventory(
        source_0709,
        single_tex_name="0709.2908.tex",
        expected_sha256=ELKIES_0709_SOURCE_SHA256,
    )
    inventory_0802, tex_0802 = _source_inventory(
        source_0802,
        single_tex_name="ants8+.tex",
        expected_sha256=ELKIES_0802_SOURCE_SHA256,
    )

    markers_0709 = {
        "shimura_curve_equation": (
            "u^2 = 16 t^6 - 19 t^4 + 88 t^2 - 48" in tex_0709
        ),
        "rank17_k3_claim": (
            "elliptic K3 surface of \\MW\\ rank $17$" in tex_0709
        ),
        "rank18_base_change_claim": (
            "$18$ for the \\MW\\ rank of a nonconstant elliptic curve"
            in tex_0709
        ),
        "schematic_17_section_identity": (
            "$y_i(T)^2 = x_i(T)^3 + a_4(T) x_i(T) + a_6(T)$" in tex_0709
        ),
        "quadratic_base_change_method_discussed": (
            "Incrementing the rank via quadratic base change" in tex_0709
        ),
    }
    markers_0802 = {
        "references_rank17_x_6_79_construction": (
            "\\X(6,79)/\\ang{w_{6\\cdot79}}" in tex_0802
            and "maximal \\MW\\ rank" in tex_0802
            and "\\cite{NDE:high_rank}" in tex_0802
        ),
        "says_full_results_could_be_available_online": (
            "full set of results can be made" in tex_0802
            and "available online" in tex_0802
        ),
    }
    # The first marker's wording varies, so record the direct stable tokens too.
    markers_0802["mentions_x_6_79"] = "\\X(6,79)" in tex_0802
    markers_0802["mentions_primary_rank_paper"] = "0709.2908" in tex_0802

    pinned_revisions = bool(
        inventory_0709["matches_pinned_revision"]
        and inventory_0802["matches_pinned_revision"]
    )
    described_claims_present = all(markers_0709.values())
    if not pinned_revisions:
        raise ValueError(
            "an official arXiv source bundle changed; inspect and repin it before "
            "reusing the missing-data assessment"
        )
    if not described_claims_present:
        raise ValueError("a pinned Elkies construction marker was not recovered")
    return {
        "audit_date": date.today().isoformat(),
        "scope": "official primary arXiv source bundles and their bundled files",
        "sources": {
            "arxiv_0709_2908": {
                "abstract_url": ELKIES_0709_ABSTRACT_URL,
                "source_url": ELKIES_0709_SOURCE_URL,
                "inventory": inventory_0709,
                "semantic_markers": markers_0709,
                "rank18_or_rank17_model_data_links_in_tex": [],
                "evidence_locations_in_pinned_source": {
                    "existence_and_rank_claims": "mfo07_arxiv.tex lines 745-763",
                    "schematic_reconstruction_method": (
                        "mfo07_arxiv.tex lines 792-814"
                    ),
                },
            },
            "arxiv_0802_1301": {
                "abstract_url": ELKIES_0802_ABSTRACT_URL,
                "source_url": ELKIES_0802_SOURCE_URL,
                "inventory": inventory_0802,
                "semantic_markers": markers_0802,
                "x_6_79_model_data_links_in_tex": [],
                "evidence_locations_in_pinned_source": {
                    "reference_to_x_6_79_construction": (
                        "ants8+.tex lines 104-110"
                    ),
                    "online_availability_statement_without_data_link": (
                        "ants8+.tex lines 134-142"
                    ),
                },
            },
        },
        "published_claim_recovered": described_claims_present,
        "required_rank18_reproduction_objects": {
            "rank17_weierstrass_coefficients_a4_a6": False,
            "rank17_section_coordinates_1_through_17": False,
            "quadratic_base_change_map": False,
            "eighteenth_section_coordinates": False,
        },
        "source_bundle_has_ancillary_data": bool(
            inventory_0709["non_tex_data_members"]
            or inventory_0802["non_tex_data_members"]
        ),
        "assessment": {
            "pinned_source_revisions_verified": pinned_revisions,
            "explicit_rank18_family_reproducible": False,
            "status": "source-blocked-before-instantiation",
            "reason": (
                "The official sources assert and describe the construction, "
                "but do not publish the curve coefficients, section coordinates, "
                "or quadratic base-change data needed to instantiate it; the "
                "printed Shimura-curve equation alone does not determine those "
                "missing maps and models."
            ),
            "specialization_screen_performed": False,
            "specialization_screen_reason": (
                "No exact rank-18 model can be constructed from the audited data."
            ),
        },
    }


def build_report(source_0709: bytes, source_0802: bytes) -> dict[str, Any]:
    source_audit = audit_official_sources(source_0709, source_0802)
    symbolic = symbolic_verification()
    rational_t_two = verify_rational_specialization(
        PUBLISHED_INDEPENDENCE_SPECIALIZATION
    )
    jacobian = short_jacobian_coefficients(PUBLISHED_INDEPENDENCE_SPECIALIZATION)
    return {
        "schema": "elliptic-elkies-rank18-source-audit-v1",
        "reproduce_command": REPRODUCE_COMMAND,
        "script_sha256": sha256(Path(__file__).resolve().read_bytes()).hexdigest(),
        "software": {
            "python": platform.python_version(),
            "sympy": sp.__version__,
        },
        "target": {
            "desired_rank": 21,
            "log_conductor_upper_bound": 182.72,
            "alternate_desired_rank": 30,
        },
        "elkies_rank18": source_audit,
        "explicit_fallback": {
            "author": "Shoichi Kihara",
            "title": "On an elliptic curve over Q(t) of rank >= 14",
            "year": 2001,
            "primary_source_doi": KIHARA_PRIMARY_SOURCE_DOI,
            "primary_source_pdf": KIHARA_PRIMARY_SOURCE_PDF,
            "published_generic_rank_lower_bound": PUBLISHED_RANK_LOWER_BOUND,
            "published_group_origin": "P15",
            "published_independent_points": "P1,...,P14",
            "published_independence_specialization_t": str(
                PUBLISHED_INDEPENDENCE_SPECIALIZATION
            ),
            "published_height_determinant_approx": (
                PUBLISHED_HEIGHT_DETERMINANT_APPROX
            ),
            "local_status": (
                "generic point identities verified exactly; published numerical "
                "height independence not independently certified"
            ),
            "symbolic_verification": symbolic,
            "exact_t2_verification": rational_t_two,
            "t2_short_jacobian_coefficient_numerator_digits": [
                len(str(abs(value.numerator))) for value in jacobian
            ],
            "ready_for_bounded_specialization_screen": True,
        },
        "result": {
            "rank18_screen_result": "not-run-source-blocked",
            "breakthrough_curve_found": False,
            "best_new_certified_rank": None,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-0709", type=Path)
    parser.add_argument("--source-0802", type=Path)
    parser.add_argument("--download-timeout", type=float, default=20.0)
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()
    if arguments.download_timeout <= 0 or arguments.download_timeout > 60:
        raise ValueError("download timeout must lie in (0,60]")
    source_0709 = (
        arguments.source_0709.read_bytes()
        if arguments.source_0709
        else _download(ELKIES_0709_SOURCE_URL, arguments.download_timeout)
    )
    source_0802 = (
        arguments.source_0802.read_bytes()
        if arguments.source_0802
        else _download(ELKIES_0802_SOURCE_URL, arguments.download_timeout)
    )
    report = build_report(source_0709, source_0802)
    encoded = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if arguments.output:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(encoded, encoding="utf-8")
    else:
        print(encoded, end="")


if __name__ == "__main__":
    main()
