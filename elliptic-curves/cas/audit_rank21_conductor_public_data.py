#!/usr/bin/env python3
"""Audit public exact data for the rank-21 / log-conductor frontier.

This is a bounded source audit, not a proof that no better curve exists.  It
pins the authoritative snapshots inspected on 2026-08-14, replays every
printed point on Dujella's rank-21 through rank-29 pages, and locally computes
exact conductors whenever the public/local factorization data make that replay
fast and deterministic.  The remaining three record curves are far outside
the target by Dujella's machine table and are retained as source-reported
values rather than silently promoted to local computations.

The artifact deliberately distinguishes:

* a published claim that displayed points are independent;
* exact point membership, which this script checks itself;
* exact finite-reduction rank certificates already present in this repository;
* a source-reported conductor logarithm; and
* an exact conductor recomputed here with PARI/GP.
"""

from __future__ import annotations

import argparse
import ast
from dataclasses import dataclass
from decimal import Decimal
from fractions import Fraction
from hashlib import sha256
import html
from io import BytesIO
import json
from pathlib import Path
import platform
import re
import shutil
import subprocess
from typing import Any, Iterable
from urllib.request import Request, urlopen

from pypdf import PdfReader, __version__ as pypdf_version


AUDIT_DATE = "2026-08-14"
STRICT_LOG_CONDUCTOR_TARGET = Decimal("182.72")
DEFAULT_OUTPUT = Path(
    "artifacts/generated-results/"
    "elliptic_rank21_conductor_public_data_audit.json"
)
REPRODUCE_COMMAND = (
    ".venv/bin/python "
    "elliptic-curves/cas/audit_rank21_conductor_public_data.py "
    "--source-dir /private/tmp/rank21-public-source-audit "
    "--artifact-dir artifacts/generated-results "
    "--output artifacts/generated-results/"
    "elliptic_rank21_conductor_public_data_audit.json"
)


@dataclass(frozen=True)
class SourceSpec:
    filename: str
    url: str
    sha256: str
    authority: str
    role: str


DUJELLA_BASE = "https://web.math.pmf.unizg.hr/~duje/tors"
SOURCE_SPECS: dict[str, SourceSpec] = {
    "dujella_torsion_frontier": SourceSpec(
        filename="tors.html",
        url=f"{DUJELLA_BASE}/tors.html",
        sha256=(
            "186d212466787027a130cde6f259bb6bf850cedf9561d9759f38f8304ae354d8"
        ),
        authority="Andrej Dujella's author-maintained high-rank catalogue",
        role="current rank lower bounds for every Mazur torsion group",
    ),
    "dujella_rank_history": SourceSpec(
        filename="rankhist.html",
        url=f"{DUJELLA_BASE}/rankhist.html",
        sha256=(
            "7cc15035e97d6c9d7c46f2031b1fa338017b3895351002fbe174b3d7715e4835"
        ),
        authority="Andrej Dujella's author-maintained rank-record catalogue",
        role="historical E/Q rank-record index and primary references",
    ),
    "dujella_coefficients": SourceSpec(
        filename="coeff.txt",
        url=f"{DUJELLA_BASE}/coeff.txt",
        sha256=(
            "170a2168984f7cd5cf38598d78b7dda100e8ffbb6115d786b5b98bc29e04ced3"
        ),
        authority="Andrej Dujella's downloadable machine table",
        role="machine-readable Weierstrass models by torsion group",
    ),
    "dujella_conductors": SourceSpec(
        filename="conductor.txt",
        url=f"{DUJELLA_BASE}/conductor.txt",
        sha256=(
            "0d3c0f501a014943b3e6c958cbde9250833879ace4448dde31cc6478d5c97e08"
        ),
        authority="Andrej Dujella's downloadable machine table",
        role="source-reported natural logarithms of conductors",
    ),
    "lmfdb_stats": SourceSpec(
        filename="lmfdb_stats.html",
        url="https://www.lmfdb.org/EllipticCurve/Q/stats",
        sha256=(
            "2eb81a54a5539c0723bb81e34d8aaae5e6067eaf83ed78b32bea4e8075976880"
        ),
        authority="official LMFDB statistics page",
        role="current E/Q database size, conductor ceiling, and rank distribution",
    ),
    "elkies_small_conductor_table": SourceSpec(
        filename="elkies_womack2001.html",
        url="https://people.math.harvard.edu/~elkies/womack2001.html",
        sha256=(
            "cbc48c2d7bc10edbe8a7c604806c022a5d49fadd589d90a72d750098caa9f919"
        ),
        authority="Noam Elkies's Harvard author site",
        role="author-hosted best-conductor table, whose displayed scope ends at rank 11",
    ),
    "elkies_watkins_paper": SourceSpec(
        filename="elkies_watkins_antsvi.pdf",
        url="https://magma.maths.usyd.edu.au/users/watkins/papers/antsVI.pdf",
        sha256=(
            "e55dec449c92ac1b952d4c7c6ea8c585ecec2c34249d2c19be93e3b1149374dc"
        ),
        authority="Mark Watkins's University of Sydney author directory",
        role="primary small-conductor paper and model/conductor tables for ranks 5 through 11",
    ),
    "fermigier_1997": SourceSpec(
        filename="fermigier_1997.pdf",
        url="https://matwbn.icm.edu.pl/ksiazki/aa/aa82/aa8243.pdf",
        sha256=(
            "9e0455228382c74b0e558b80b28346d5440531eb797b2539d3d379f1c86d77e4"
        ),
        authority="Acta Arithmetica official archive",
        role="primary rank-22 construction, specialization, model, and points",
    ),
}

for _rank, _digest in {
    21: "ec81b94ad1b2a36295c9c2eec312b746477549cb4972100645880d877b3ddd56",
    22: "711624e3c11a7e6a9413cd25872697fb51d3ca6060380beab045380cb5bb4c4c",
    23: "c9aaa9d28391f8c15aeebbf44f86206bfe560bf09a380a4bfb5b9fefd3c5ed46",
    24: "2407b9d77ccc3cdfe7c964fbc255aaff8b135ca953aec2c6d3afd72201e434cc",
    25: "3c24b0da3dbb0a96aee924932d8ff00c19e0fa09424009df7e2e8b88d5a4105c",
    26: "5f2d8c4c58b8db74142fbd29391c59a5370ccc9bc05768fd5402babef8404c77",
    27: "a6045f181da0489fa1d02f08b1f000971fb4eb85d74a52f94912ca715a4560bb",
    28: "9b8e9aa08b93f8da5054d371ddb7b230f60016e75635a618baacc99b6e8e6137",
    29: "94b673dd0655dfa2d8eae9247488c5475ae904c3a230f382d39e4116b76956aa",
}.items():
    SOURCE_SPECS[f"dujella_rank{_rank}"] = SourceSpec(
        filename=f"rk{_rank}.html",
        url=f"{DUJELLA_BASE}/rk{_rank}.html",
        sha256=_digest,
        authority="Andrej Dujella's author-maintained rank-record catalogue",
        role=f"printed rank-{_rank} model and {_rank} claimed independent points",
    )


EXPECTED_AUTHORS = {
    21: "Nagao - Kouya (1994)",
    22: "Fermigier (1997)",
    23: "Martin - McMillen (1998)",
    24: "Martin - McMillen (2000)",
    25: "Elkies (2006)",
    26: "Elkies (2006)",
    27: "Elkies (2016)",
    28: "Elkies (2006)",
    29: "Elkies - Klagsbrun (2024)",
}

# Complete minimal-discriminant factorizations for the exact local replay
# subset.  Every listed base is independently checked by GP's ``isprime``;
# Python checks the product before GP is allowed to compute a conductor.
DISCRIMINANT_FACTORIZATIONS: dict[int, dict[int, int]] = {
    21: {
        2: 13,
        3: 2,
        5: 4,
        7: 4,
        13: 4,
        17: 3,
        23: 3,
        47: 4,
        4507: 1,
        115482611374267602141168398241396608699381902319617225736297616061235976719: 1,
    },
    22: {
        2: 2,
        3: 9,
        5: 2,
        7: 6,
        13: 6,
        17: 4,
        37: 3,
        47293: 1,
        270704849145149791: 1,
        60794657878864337775664712674231370427122734380997: 1,
    },
    23: {
        2: 4,
        3: 8,
        5: 6,
        11: 2,
        13: 2,
        17: 2,
        19: 5,
        23: 4,
        199: 2,
        17858193374841750901974789649: 1,
        1006218106634655545344494448610726356220703995276273: 1,
    },
    24: {
        2: 2,
        3: 9,
        5: 2,
        11: 6,
        13: 2,
        17: 2,
        29: 2,
        31: 3,
        41: 2,
        458619970494582607679296750333015081: 1,
        264240973182971699094661154229360236070105974082503: 1,
    },
    25: {
        2: 16,
        3: 13,
        5: 6,
        7: 7,
        11: 3,
        13: 4,
        157: 1,
        283: 2,
        72551: 1,
        895613: 1,
        39466373175449369221626360280383658386480681928740709773906484427935361743867171: 1,
    },
    29: {
        2: 19,
        3: 7,
        5: 7,
        7: 4,
        11: 5,
        13: 3,
        17: 4,
        31: 3,
        41: 2,
        43: 2,
        61: 2,
        233: 1,
        241: 2,
        4139: 1,
        678146849364709860535420504397393: 1,
        159788990966780131363155786084695062643236502969: 1,
        4402149008473369392540402625019227412319473055901: 1,
    },
}

EXPECTED_EXACT_CONDUCTORS = {
    21: 26112076915897777815571388664430310998157918697219343275140810790098571234096793308930,
    22: 22720638514787473197194583889675055980109503436060704437972911338086049759883790,
    23: 113964706497998622564257410134352232583764612126645513242486361984405434455972952111271430,
    24: 325763846420202806803685538805031963849862314963649942610049440603050996133862755155678267754410,
    25: 3421625153333343825198002959770328244842456471613462044687992344664286334910845403447307652180335890,
    29: 188691278153039213628440510802703336227022485808859170041743892297765880547867276243424847942033952274950770795963060482541371967126966027206485174370,
}

EXPECTED_ROOT_NUMBERS = {21: -1, 22: 1, 23: -1, 24: 1, 25: -1, 29: -1}


def _download(url: str, timeout: float) -> bytes:
    request = Request(
        url,
        headers={
            "User-Agent": "elliptic-rank21-conductor-public-audit/1.0",
            "Accept": "text/html, application/pdf, text/plain, */*",
        },
    )
    with urlopen(request, timeout=timeout) as response:
        return response.read()


def verify_payload(name: str, payload: bytes) -> dict[str, Any]:
    """Check a source hash before any generated artifact can be written."""

    spec = SOURCE_SPECS[name]
    actual = sha256(payload).hexdigest()
    if actual != spec.sha256:
        raise ValueError(
            f"source {name!r} changed: expected {spec.sha256}, got {actual}; "
            "inspect and repin before reusing this audit"
        )
    return {
        "url": spec.url,
        "authority": spec.authority,
        "role": spec.role,
        "bytes": len(payload),
        "sha256": actual,
        "matches_pinned_revision": True,
    }


def load_sources(
    *, source_dir: Path | None, timeout: float
) -> tuple[dict[str, bytes], dict[str, dict[str, Any]]]:
    payloads: dict[str, bytes] = {}
    inventories: dict[str, dict[str, Any]] = {}
    for name, spec in SOURCE_SPECS.items():
        payload = (
            (source_dir / spec.filename).read_bytes()
            if source_dir is not None
            else _download(spec.url, timeout)
        )
        inventories[name] = verify_payload(name, payload)
        payloads[name] = payload
    return payloads, inventories


def _plain_html(payload: bytes, *, encoding: str = "utf-8") -> str:
    text = payload.decode(encoding, errors="replace")
    text = re.sub(r"<sup>(.*?)</sup>", r"^\1", text, flags=re.I | re.S)
    text = re.sub(r"<sub>(.*?)</sub>", r"_\1", text, flags=re.I | re.S)
    text = re.sub(r"<br\s*/?>", "\n", text, flags=re.I)
    text = re.sub(r"<[^>]+>", " ", text)
    return re.sub(r"\s+", " ", html.unescape(text)).strip()


def _signed(sign: str, value: str) -> int:
    return int(value) if sign == "+" else -int(value)


def parse_record_page(payload: bytes, rank: int) -> dict[str, Any]:
    """Recover one exact generalized model and all displayed points."""

    text = _plain_html(payload)
    equation_match = re.search(
        r"(y\^2.*?)Independent points", text, flags=re.I | re.S
    )
    if equation_match is None:
        raise ValueError(f"rank-{rank} equation was not recovered")
    equation = equation_match.group(1).strip()
    lhs, rhs = (part.strip() for part in equation.split("=", 1))
    a1 = int("+ xy" in lhs)
    a3 = int(re.search(r"\+ y$", lhs) is not None)
    rhs_match = re.fullmatch(
        r"x\^3(?: ([+-]) x\^2)? ([+-]) (\d+)x ([+-]) (\d+)", rhs
    )
    if rhs_match is None:
        raise ValueError(f"rank-{rank} right side was not parsed: {rhs!r}")
    a2 = (
        0
        if rhs_match.group(1) is None
        else (1 if rhs_match.group(1) == "+" else -1)
    )
    a4 = _signed(rhs_match.group(2), rhs_match.group(3))
    a6 = _signed(rhs_match.group(4), rhs_match.group(5))
    coefficients = (a1, a2, a3, a4, a6)

    points = []
    for match in re.finditer(
        r"P_(\d+)\s*=\s*\[\s*([-+]?\d+(?:/\d+)?)\s*,\s*"
        r"([-+]?\d+(?:/\d+)?)\s*\]",
        text,
    ):
        points.append(
            (int(match.group(1)), Fraction(match.group(2)), Fraction(match.group(3)))
        )
    indices = [index for index, _, _ in points]
    if indices != list(range(1, rank + 1)):
        raise AssertionError(
            f"rank-{rank} page has point indices {indices}, expected 1..{rank}"
        )
    bad_points = []
    for index, x_value, y_value in points:
        lhs_value = y_value**2 + a1 * x_value * y_value + a3 * y_value
        rhs_value = (
            x_value**3 + a2 * x_value**2 + a4 * x_value + a6
        )
        if lhs_value != rhs_value:
            bad_points.append(index)
    if bad_points:
        raise AssertionError(f"rank-{rank} off-curve points: {bad_points}")

    author = EXPECTED_AUTHORS[rank]
    if author not in text:
        raise AssertionError(f"rank-{rank} author/year marker changed")
    return {
        "rank_lower_bound_printed_by_source": rank,
        "source_claim": f"{rank} independent points of infinite order",
        "author_and_year": author,
        "equation": equation,
        "weierstrass_coefficients": [str(value) for value in coefficients],
        "printed_point_count": len(points),
        "printed_point_indices": indices,
        "exact_membership_checks_passed": len(points),
        "independence_reproved_from_this_page_alone": False,
    }


def _parse_sectioned_lines(
    payload: bytes, *, value_parser
) -> dict[str, list[Any]]:
    sections: dict[str, list[Any]] = {}
    current: str | None = None
    for raw_line in payload.decode("utf-8", errors="replace").splitlines():
        line = raw_line.strip()
        if re.fullmatch(r"Z(?:\d+|\d+\*Z\d+)", line):
            current = line
            sections[current] = []
            continue
        if current is None or not line:
            continue
        parsed = value_parser(line)
        if parsed is not None:
            sections[current].append(parsed)
    return sections


def parse_machine_tables(
    coefficients_payload: bytes,
    conductor_payload: bytes,
    records: dict[int, dict[str, Any]],
) -> dict[str, Any]:
    def coefficient_parser(line: str) -> tuple[int, ...] | None:
        if not line.startswith("["):
            return None
        # The pinned table contains one stray legacy-encoding byte after a
        # closing bracket in the Z6 block.  Parse the bracketed vector itself
        # and record the resulting row-count discrepancy below.
        closing = line.find("]")
        if closing < 0:
            raise ValueError(f"unterminated coefficient vector {line!r}")
        value = ast.literal_eval(line[: closing + 1])
        if not isinstance(value, list) or len(value) != 5:
            raise ValueError(f"malformed coefficient vector {line!r}")
        return tuple(int(entry) for entry in value)

    def log_parser(line: str) -> str | None:
        return line if re.fullmatch(r"\d+\.\d+", line) else None

    coefficient_sections = _parse_sectioned_lines(
        coefficients_payload, value_parser=coefficient_parser
    )
    conductor_sections = _parse_sectioned_lines(
        conductor_payload, value_parser=log_parser
    )
    if len(coefficient_sections["Z1"]) != 6:
        raise AssertionError("the pinned Z1 coefficient block no longer has six models")
    if len(conductor_sections["Z1"]) != 6:
        raise AssertionError("the pinned Z1 conductor block no longer has six values")

    z1_records = []
    for index, rank in enumerate(range(29, 23, -1)):
        coefficients = coefficient_sections["Z1"][index]
        page_coefficients = tuple(
            int(value) for value in records[rank]["weierstrass_coefficients"]
        )
        if coefficients != page_coefficients:
            raise AssertionError(f"machine/page model mismatch at rank {rank}")
        log_conductor = conductor_sections["Z1"][index]
        z1_records.append(
            {
                "rank": rank,
                "weierstrass_coefficients": [str(value) for value in coefficients],
                "source_reported_log_conductor": log_conductor,
                "margin_above_strict_target": str(
                    Decimal(log_conductor) - STRICT_LOG_CONDUCTOR_TARGET
                ),
            }
        )
    return {
        "coefficient_sections": {
            name: len(values) for name, values in coefficient_sections.items()
        },
        "conductor_sections": {
            name: len(values) for name, values in conductor_sections.items()
        },
        "coefficient_vector_count": sum(map(len, coefficient_sections.values())),
        "conductor_log_count": sum(map(len, conductor_sections.values())),
        "known_z6_row_count_mismatch_in_pinned_download": (
            len(coefficient_sections["Z6"]) == 165
            and len(conductor_sections["Z6"]) == 166
        ),
        "z1_current_and_predecessor_records": z1_records,
    }


def parse_torsion_frontier(payload: bytes) -> dict[str, Any]:
    source = payload.decode("utf-8", errors="replace")
    ranks = {
        filename.removesuffix(".html"): int(rank)
        for filename, rank in re.findall(
            r'<a href="(z[^"/]+\.html)">(\d+)</a>', source, flags=re.I
        )
    }
    if ranks.get("z1") != 29 or ranks.get("z2") != 20:
        raise AssertionError("Dujella torsion frontier markers changed")
    nontrivial = {name: rank for name, rank in ranks.items() if name != "z1"}
    return {
        "rank_lower_bounds_by_page": ranks,
        "trivial_torsion_frontier": ranks["z1"],
        "largest_nontrivial_torsion_frontier": max(nontrivial.values()),
        "nontrivial_torsion_page_attaining_that_frontier": sorted(
            name for name, rank in nontrivial.items() if rank == max(nontrivial.values())
        ),
        "any_nontrivial_torsion_frontier_at_least_21": any(
            rank >= 21 for rank in nontrivial.values()
        ),
    }


def parse_rank_history(payload: bytes) -> dict[str, Any]:
    source = payload.decode("utf-8", errors="replace")
    linked = [
        int(rank)
        for rank, label in re.findall(
            r'href="rk(\d+)\.html"[^>]*>(\d+)</a>', source, flags=re.I
        )
        if rank == label
    ]
    unique = sorted(set(linked))
    text = _plain_html(payload)
    return {
        "linked_record_ranks": unique,
        "highest_linked_record_rank": max(unique),
        "rank30_page_linked": 30 in unique,
        "current_record_rank29_marker": (
            "current record is an example of elliptic curve with rank >= 29"
            in text.replace("≥", ">=")
        ),
    }


def parse_lmfdb_stats(payload: bytes) -> dict[str, Any]:
    text = payload.decode("utf-8", errors="replace")
    summary = re.search(
        r"includes \$([0-9{},]+)\$.*?in \$([0-9{},]+)\$.*?"
        r"conductor.*?at most \$([0-9{},]+)\$",
        text,
        flags=re.I | re.S,
    )
    if summary is None:
        raise ValueError("LMFDB summary counts were not parsed")
    ranks = {
        int(rank): int(count)
        for rank, count in re.findall(
            r"href='/EllipticCurve/Q/\?rank=(\d+)'>(\d+)</a>", text
        )
    }
    if not ranks:
        raise ValueError("LMFDB rank distribution was not parsed")
    number = lambda value: int(value.replace("{", "").replace("}", "").replace(",", ""))
    return {
        "curve_count": number(summary.group(1)),
        "isogeny_class_count": number(summary.group(2)),
        "largest_conductor_in_database": number(summary.group(3)),
        "rank_counts": {str(rank): count for rank, count in sorted(ranks.items())},
        "largest_rank_displayed": max(ranks),
        "rank_at_least_21_entry_present": any(rank >= 21 for rank in ranks),
        "scope_note": "absence here is database-bounded, not a global nonexistence theorem",
    }


def parse_elkies_small_conductor_table(payload: bytes) -> dict[str, Any]:
    source = payload.decode("iso-8859-1", errors="replace")
    rows = []
    for row in re.findall(r"<tr[^>]*>(.*?)</tr>", source, flags=re.I | re.S):
        cells = []
        for cell in re.findall(r"<td[^>]*>(.*?)</td>", row, flags=re.I | re.S):
            cell_text = re.sub(r"<[^>]+>", " ", cell)
            cells.append(re.sub(r"\s+", " ", html.unescape(cell_text)).strip())
        if len(cells) == 5 and cells[0].isdigit():
            rows.append(
                {
                    "rank": int(cells[0]),
                    "model": cells[1],
                    "conductor": cells[2],
                    "log_conductor": cells[3],
                    "source": cells[4],
                }
            )
    if [row["rank"] for row in rows] != list(range(12)):
        raise AssertionError("Elkies author table no longer has exactly ranks 0..11")
    return {
        "row_count": len(rows),
        "displayed_rank_range": [rows[0]["rank"], rows[-1]["rank"]],
        "highest_displayed_rank": rows[-1]["rank"],
        "rank21_in_scope": False,
        "rank11_row": rows[-1],
    }


def _pdf_text(payload: bytes) -> tuple[str, int]:
    reader = PdfReader(BytesIO(payload))
    return "\n".join(page.extract_text() or "" for page in reader.pages), len(reader.pages)


def parse_pdf_sources(
    fermigier_payload: bytes,
    watkins_payload: bytes,
) -> dict[str, Any]:
    fermigier_text, fermigier_pages = _pdf_text(fermigier_payload)
    watkins_text, watkins_pages = _pdf_text(watkins_payload)
    assignments = re.findall(
        r"t\s*=\s*([0-9]+)\s*/\s*([0-9]+)", fermigier_text
    )
    specializations = [f"{numerator}/{denominator}" for numerator, denominator in assignments]
    score_rows = {}
    for rank in range(19, 23):
        match = re.search(
            rf"E\s*{rank}\s+((?:\d+\.\d+\s*){{8}})", fermigier_text
        )
        if match is None:
            raise ValueError(f"Fermigier E{rank} score row was not parsed")
        score_rows[f"E{rank}"] = re.findall(r"\d+\.\d+", match.group(1))
    if specializations != ["19754/39"]:
        raise AssertionError("Fermigier specialization inventory changed")
    printed_a4 = "940299517776391362903023121165864"
    printed_a6 = "10707363070719743033425295515449274534651125011362"
    if printed_a4 not in fermigier_text or printed_a6 not in fermigier_text:
        raise AssertionError("Fermigier printed model was not recovered")
    if not re.search(r"r\s*=\s*6\s*,\s*7\s*,\s*\.\.\.\s*,\s*11", watkins_text):
        raise AssertionError("Elkies-Watkins abstract rank range changed")
    if not re.search(
        r"Table\s*2\.\s*Low conductor records for ranks\s*5\s*[–-]\s*11",
        watkins_text,
    ):
        raise AssertionError("Elkies-Watkins Table 2 scope marker changed")
    return {
        "fermigier": {
            "page_count": fermigier_pages,
            "rational_t_specializations_printed": specializations,
            "only_printed_specialization_is_record_fiber": True,
            "record_specialization_normalized_repository_parameter": "39508/39",
            "printed_rank22_model_coefficients_a4_a6": [
                f"-{printed_a4}",
                printed_a6,
            ],
            "comparison_score_rows": score_rows,
            "adjacent_E19_E20_E21_exact_models_or_parameters_printed": False,
            "adjacent_rows_are_score_only": True,
        },
        "elkies_watkins": {
            "page_count": watkins_pages,
            "abstract_search_rank_range": [5, 11],
            "table2_title": "Low conductor records for ranks 5-11",
            "rank21_in_scope": False,
        },
    }


def curve_discriminant(coefficients: Iterable[int | Fraction]) -> Fraction:
    a1, a2, a3, a4, a6 = map(Fraction, coefficients)
    b2 = a1**2 + 4 * a2
    b4 = a1 * a3 + 2 * a4
    b6 = a3**2 + 4 * a6
    b8 = a1**2 * a6 + 4 * a2 * a6 - a1 * a3 * a4 + a2 * a3**2 - a4**2
    return -b2**2 * b8 - 8 * b4**3 - 27 * b6**2 + 9 * b2 * b4 * b6


def j_invariant(coefficients: Iterable[int | Fraction]) -> Fraction | None:
    a1, a2, a3, a4, a6 = map(Fraction, coefficients)
    b2 = a1**2 + 4 * a2
    b4 = a1 * a3 + 2 * a4
    c4 = b2**2 - 24 * b4
    discriminant = curve_discriminant((a1, a2, a3, a4, a6))
    return None if discriminant == 0 else c4**3 / discriminant


def _gp_version() -> str:
    executable = shutil.which("gp")
    if executable is None:
        raise FileNotFoundError("PARI/GP executable 'gp' was not found")
    result = subprocess.run(
        [executable, "-q"],
        input="print(version());\nquit\n",
        text=True,
        capture_output=True,
        check=True,
        timeout=10,
    )
    return result.stdout.strip()


def exact_conductor_replays(
    records: dict[int, dict[str, Any]], *, timeout: float
) -> tuple[dict[int, dict[str, Any]], str]:
    """Run one pre-factored exact GP batch for ranks 21--25 and 29."""

    executable = shutil.which("gp")
    if executable is None:
        raise FileNotFoundError("PARI/GP executable 'gp' was not found")
    all_primes = sorted(
        {prime for factors in DISCRIMINANT_FACTORIZATIONS.values() for prime in factors}
    )
    commands = [
        "default(realprecision,60);",
        f"P=[{','.join(map(str, all_primes))}];",
        'print("PRIMES_BEGIN");print(vecsum(vector(#P,i,isprime(P[i]))));print(#P);print("PRIMES_END");',
        "addprimes(P);",
    ]
    expected_discriminants: dict[int, int] = {}
    for rank in sorted(DISCRIMINANT_FACTORIZATIONS):
        coefficients = tuple(
            int(value) for value in records[rank]["weierstrass_coefficients"]
        )
        discriminant = curve_discriminant(coefficients)
        if discriminant.denominator != 1:
            raise AssertionError("an integral source model had nonintegral discriminant")
        expected_discriminants[rank] = int(discriminant)
        factor_product = 1
        for prime, exponent in DISCRIMINANT_FACTORIZATIONS[rank].items():
            factor_product *= prime**exponent
        if factor_product != abs(int(discriminant)):
            raise AssertionError(f"rank-{rank} factorization product mismatch")
        vector = ",".join(map(str, coefficients))
        commands.extend(
            [
                f"E=ellminimalmodel(ellinit([{vector}]));G=ellglobalred(E);",
                f'print("RANK_{rank}_BEGIN");',
                "print(E.a1);print(E.a2);print(E.a3);print(E.a4);print(E.a6);",
                "print(E.disc);print(G[1]);print(log(G[1]));print(ellrootno(E));",
                f'print("RANK_{rank}_END");',
            ]
        )
    commands.append("quit")
    result = subprocess.run(
        [executable, "-q", "-f", "-s", "1000000000"],
        input="\n".join(commands) + "\n",
        text=True,
        capture_output=True,
        timeout=timeout,
    )
    if result.returncode != 0 or "***" in result.stderr:
        raise RuntimeError(f"PARI/GP conductor batch failed: {result.stderr.strip()}")
    lines = [line.strip() for line in result.stdout.splitlines() if line.strip()]

    def block(name: str) -> list[str]:
        start = lines.index(f"{name}_BEGIN") + 1
        end = lines.index(f"{name}_END")
        return lines[start:end]

    prime_check = block("PRIMES")
    if int(prime_check[0]) != int(prime_check[1]):
        raise AssertionError("at least one pinned discriminant factor is not prime")
    output: dict[int, dict[str, Any]] = {}
    for rank in sorted(DISCRIMINANT_FACTORIZATIONS):
        values = block(f"RANK_{rank}")
        model = tuple(int(value) for value in values[:5])
        source_model = tuple(
            int(value) for value in records[rank]["weierstrass_coefficients"]
        )
        if model != source_model:
            raise AssertionError(f"rank-{rank} page model was not already minimal")
        discriminant = int(values[5])
        conductor = int(values[6])
        log_conductor = values[7]
        root_number = int(values[8])
        if discriminant != expected_discriminants[rank]:
            raise AssertionError(f"rank-{rank} GP discriminant mismatch")
        if conductor != EXPECTED_EXACT_CONDUCTORS[rank]:
            raise AssertionError(f"rank-{rank} exact conductor changed")
        if root_number != EXPECTED_ROOT_NUMBERS[rank]:
            raise AssertionError(f"rank-{rank} root number changed")
        output[rank] = {
            "status": "exact_local_pari_gp_replay",
            "minimal_model": [str(value) for value in model],
            "minimal_discriminant": str(discriminant),
            "complete_discriminant_factorization": [
                {"prime": str(prime), "exponent": exponent}
                for prime, exponent in DISCRIMINANT_FACTORIZATIONS[rank].items()
            ],
            "factor_product_verified_before_gp": True,
            "all_factor_bases_proved_prime_by_gp": True,
            "conductor": str(conductor),
            "log_conductor": log_conductor,
            "root_number": root_number,
            "below_strict_log_conductor_target": (
                Decimal(log_conductor) < STRICT_LOG_CONDUCTOR_TARGET
            ),
            "margin_above_strict_target": str(
                Decimal(log_conductor) - STRICT_LOG_CONDUCTOR_TARGET
            ),
        }
    return output, _gp_version()


def _fraction_or_none(value: Any) -> Fraction | None:
    if isinstance(value, bool) or value is None:
        return None
    try:
        return Fraction(value)
    except (TypeError, ValueError, ZeroDivisionError):
        return None


def scan_local_artifact_duplicates(
    records: dict[int, dict[str, Any]],
    *,
    artifact_directory: Path,
    output_path: Path,
) -> dict[int, list[dict[str, Any]]]:
    """Find isomorphic source models via exact j-invariant fingerprints."""

    rank_by_j = {
        j_invariant(tuple(int(x) for x in record["weierstrass_coefficients"])): rank
        for rank, record in records.items()
    }
    hits: dict[int, dict[str, list[str]]] = {rank: {} for rank in records}

    def walk(value: Any, path: str, filename: str) -> None:
        if isinstance(value, dict):
            for key, child in value.items():
                walk(child, f"{path}.{key}" if path else key, filename)
            return
        if not isinstance(value, list):
            return
        if len(value) == 5:
            coefficients = [_fraction_or_none(entry) for entry in value]
            if all(entry is not None for entry in coefficients):
                try:
                    fingerprint = j_invariant(coefficients)  # type: ignore[arg-type]
                except (ArithmeticError, ValueError):
                    fingerprint = None
                rank = rank_by_j.get(fingerprint)
                if rank is not None:
                    hits[rank].setdefault(filename, []).append(path)
        for index, child in enumerate(value):
            walk(child, f"{path}[{index}]", filename)

    output_resolved = output_path.resolve()
    files = sorted(artifact_directory.glob("*.json"))
    for path in files:
        if path.resolve() == output_resolved:
            continue
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        walk(payload, "", path.name)
    answer: dict[int, list[dict[str, Any]]] = {}
    for rank, by_file in hits.items():
        answer[rank] = [
            {
                "path": str((artifact_directory / filename).as_posix()),
                "sha256": sha256((artifact_directory / filename).read_bytes()).hexdigest(),
                "coefficient_vector_paths_with_same_j_invariant": sorted(paths),
            }
            for filename, paths in sorted(by_file.items())
        ]
    return answer


def build_report(
    payloads: dict[str, bytes],
    inventories: dict[str, dict[str, Any]],
    *,
    artifact_directory: Path,
    output_path: Path,
    gp_timeout: float,
) -> dict[str, Any]:
    records = {
        rank: parse_record_page(payloads[f"dujella_rank{rank}"], rank)
        for rank in range(21, 30)
    }
    machine = parse_machine_tables(
        payloads["dujella_coefficients"],
        payloads["dujella_conductors"],
        records,
    )
    source_logs = {
        item["rank"]: item["source_reported_log_conductor"]
        for item in machine["z1_current_and_predecessor_records"]
    }
    exact_replays, gp_version = exact_conductor_replays(
        records, timeout=gp_timeout
    )
    for rank in (24, 25, 29):
        if abs(
            Decimal(exact_replays[rank]["log_conductor"])
            - Decimal(source_logs[rank])
        ) >= Decimal("1e-16"):
            raise AssertionError(f"rank-{rank} local/source log mismatch")

    duplicates = scan_local_artifact_duplicates(
        records,
        artifact_directory=artifact_directory,
        output_path=output_path,
    )
    if {rank for rank, hits in duplicates.items() if hits} != {21, 22, 29}:
        raise AssertionError("the local high-rank duplicate inventory changed")

    candidate_inventory = []
    for rank in range(21, 30):
        local = exact_replays.get(rank)
        source_log = source_logs.get(rank)
        effective_log = local["log_conductor"] if local is not None else source_log
        if effective_log is None:
            raise AssertionError(f"no conductor evidence for rank {rank}")
        candidate_inventory.append(
            {
                "rank_lower_bound_printed_by_source": rank,
                "author_and_year": records[rank]["author_and_year"],
                "weierstrass_coefficients": records[rank]["weierstrass_coefficients"],
                "printed_points_exactly_checked_on_curve": records[rank][
                    "exact_membership_checks_passed"
                ],
                "independence_status": (
                    "exact finite-reduction certificate already local"
                    if rank in (22, 29)
                    else "printed independent-point claim; membership replay only in this audit"
                ),
                "conductor_evidence": (
                    local
                    if local is not None
                    else {
                        "status": "dujella_machine_table_log_only",
                        "source_reported_log_conductor": source_log,
                        "margin_above_strict_target": str(
                            Decimal(source_log) - STRICT_LOG_CONDUCTOR_TARGET
                        ),
                        "exact_local_conductor_replayed_here": False,
                        "reason_not_promoted": (
                            "the pinned page/table supplies no complete factorization; "
                            "the source value is already more than 65 natural-log units "
                            "above the target"
                        ),
                    }
                ),
                "below_strict_log_conductor_target": (
                    Decimal(effective_log) < STRICT_LOG_CONDUCTOR_TARGET
                ),
                "local_artifact_isomorphism_matches": duplicates[rank],
                "excluded_as_already_local": bool(duplicates[rank]),
                "new_to_local_artifacts_at_audit_time": not duplicates[rank],
                "plausible_target_candidate_after_conductor_gate": False,
            }
        )

    pdf_results = parse_pdf_sources(
        payloads["fermigier_1997"],
        payloads["elkies_watkins_paper"],
    )
    script_path = Path(__file__).resolve()
    return {
        "schema_version": 1,
        "artifact_kind": "bounded_elliptic_rank21_conductor_public_data_audit",
        "status": "no_new_public_rank21_subtarget_conductor_curve_recovered",
        "audit_date": AUDIT_DATE,
        "strict_target": {
            "rank_at_least": 21,
            "natural_log_conductor_strictly_less_than": "182.72",
        },
        "reproduction": {
            "command": REPRODUCE_COMMAND,
            "script_sha256": sha256(script_path.read_bytes()).hexdigest(),
            "all_source_hashes_verified_before_write": True,
            "source_mode": "offline pinned files supported; network refresh hash-guarded",
        },
        "software": {
            "python": platform.python_version(),
            "pari_gp": gp_version,
            "pypdf": pypdf_version,
        },
        "source_inventories": inventories,
        "source_results": {
            "dujella_torsion_frontier": parse_torsion_frontier(
                payloads["dujella_torsion_frontier"]
            ),
            "dujella_rank_history": parse_rank_history(
                payloads["dujella_rank_history"]
            ),
            "dujella_machine_tables": machine,
            "lmfdb": parse_lmfdb_stats(payloads["lmfdb_stats"]),
            "elkies_author_small_conductor_table": parse_elkies_small_conductor_table(
                payloads["elkies_small_conductor_table"]
            ),
            "primary_pdfs": pdf_results,
        },
        "record_page_replays": {str(rank): records[rank] for rank in records},
        "candidate_inventory": candidate_inventory,
        "local_duplicate_policy": {
            "method": (
                "recursive scan of every length-5 numeric vector in local JSON "
                "artifacts, compared by exact j-invariant"
            ),
            "ranks_already_local": [21, 22, 29],
            "ranks_new_to_local_artifacts_but_conductor_disqualified": [23, 24, 25, 26, 27, 28],
        },
        "fermigier_adjacent_specialization_result": {
            "exact_new_adjacent_specializations_recovered": [],
            "only_exact_rational_parameter_printed_in_primary_paper": "19754/39",
            "repository_normalization_of_that_parameter": "39508/39",
            "E19_E20_E21_comparison_data_type": "eight Nagao-score values each; no model or parameter",
            "archive_blocker": (
                "the paper does not publish the screened parameter list or exact "
                "neighboring fibers, so there is no adjacent specialization dataset to replay"
            ),
        },
        "assessment": {
            "new_public_model_meeting_both_rank_and_conductor_gates": False,
            "finite_reduction_certificate_triggered_for_new_curve": False,
            "reason_no_new_certificate_was_built": (
                "every newly recovered rank-23 through rank-28 record model is "
                "already above the conductor cutoff; the rank-21, rank-22, and "
                "rank-29 models were already present locally"
            ),
            "closest_public_certified_near_miss": {
                "rank_lower_bound": 22,
                "curve": "Fermigier E22",
                "log_conductor": exact_replays[22]["log_conductor"],
                "excess_over_strict_target": exact_replays[22][
                    "margin_above_strict_target"
                ],
                "exact_local_rank_certificate_artifact": (
                    "artifacts/generated-results/elliptic_fermigier_rank22_points.json"
                ),
            },
            "public_table_frontier": (
                "Dujella rank pages give exact models and printed point sets "
                "through rank 29; no audited table contains a rank>=21 curve "
                "with log conductor below 182.72"
            ),
            "epistemic_status": (
                "bounded negative public-source audit; not a theorem excluding "
                "unpublished curves or models outside the pinned catalogues"
            ),
            "breakthrough_result_found": False,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path)
    parser.add_argument("--artifact-dir", type=Path, default=Path("artifacts/generated-results"))
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--download-timeout", type=float, default=60.0)
    parser.add_argument("--gp-timeout", type=float, default=120.0)
    args = parser.parse_args()
    payloads, inventories = load_sources(
        source_dir=args.source_dir, timeout=args.download_timeout
    )
    report = build_report(
        payloads,
        inventories,
        artifact_directory=args.artifact_dir,
        output_path=args.output,
        gp_timeout=args.gp_timeout,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        f"status={report['status']} "
        f"sources={len(report['source_inventories'])} "
        "new_target_candidates=0"
    )


if __name__ == "__main__":
    main()
