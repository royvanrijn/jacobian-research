#!/usr/bin/env python3
"""Audit public first-party sources for a reproducible rank-30 curve.

This is a bounded source audit, not a theorem that no rank-30 curve exists.
It pins the official/author-controlled snapshots inspected on 2026-08-14,
recovers the public rank-29 frontier, and records the exact missing objects
that prevent replaying either a rank-30 example or the rank-17 K3 fibration
used in the Elkies--Klagsbrun search.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from datetime import datetime
from email.utils import parsedate_to_datetime
from hashlib import sha256
import html
import json
from pathlib import Path
import platform
import re
from typing import Any
from urllib.request import Request, urlopen
import xml.etree.ElementTree as ET


AUDIT_DATE = "2026-08-14"
DEFAULT_OUTPUT = Path(
    "artifacts/generated-results/elliptic_rank30_public_source_audit.json"
)
REPRODUCE_COMMAND = (
    ".venv/bin/python elliptic-curves/cas/audit_rank30_public_sources.py "
    "--source-dir /private/tmp/rank30-public-source-audit "
    "--output artifacts/generated-results/"
    "elliptic_rank30_public_source_audit.json"
)


@dataclass(frozen=True)
class SourceSpec:
    filename: str
    url: str
    sha256: str
    authority: str
    role: str


SOURCE_SPECS: dict[str, SourceSpec] = {
    "nmbrthry_index": SourceSpec(
        filename="nmbrthry_index.html",
        url=(
            "https://listserv.nodak.edu/cgi-bin/wa.exe?"
            "A0=NMBRTHRY&TOC=&S=b"
        ),
        sha256=(
            "d08ac1964de817faddbc08e6e3269d7315aa67a9bd67d4cb8ef8e0473fc7b4c0"
        ),
        authority="official NMBRTHRY LISTSERV archive",
        role="complete archive index snapshot through the audit date",
    ),
    "nmbrthry_rank29": SourceSpec(
        filename="nmbrthry_rank29.html",
        url=(
            "https://listserv.nodak.edu/cgi-bin/wa.exe?"
            "A2=NMBRTHRY;b9d018b1.2409&S=b"
        ),
        sha256=(
            "93151103437b6e1c190882963e1dcb2785258bf42c7a9eb60e2bb5a1fa42266a"
        ),
        authority="official NMBRTHRY LISTSERV archive; author announcement",
        role="primary announcement of the rank-29 curve and search method",
    ),
    "dujella_rank_history": SourceSpec(
        filename="dujella_rank_history.html",
        url="https://web.math.pmf.unizg.hr/~duje/tors/rankhist.html",
        sha256=(
            "7cc15035e97d6c9d7c46f2031b1fa338017b3895351002fbe174b3d7715e4835"
        ),
        authority="Andrej Dujella's author-maintained rank-record catalogue",
        role="current public record index",
    ),
    "dujella_rank29": SourceSpec(
        filename="dujella_rank29.html",
        url="https://web.math.pmf.unizg.hr/~duje/tors/rk29.html",
        sha256=(
            "94b673dd0655dfa2d8eae9247488c5475ae904c3a230f382d39e4116b76956aa"
        ),
        authority="Andrej Dujella's author-maintained rank-record catalogue",
        role="public rank-29 model and 29 point coordinates",
    ),
    "arxiv_rank30": SourceSpec(
        filename="arxiv_rank30.atom",
        url=(
            "https://export.arxiv.org/api/query?"
            "search_query=all:%22rank%2030%22%20AND%20"
            "all:%22elliptic%20curve%22&start=0&max_results=100&"
            "sortBy=submittedDate&sortOrder=descending"
        ),
        sha256=(
            "f1ff392e51a9e9229b04ae7c163797884311d0b401c89bea02aa4bfa0cae4d5a"
        ),
        authority="official arXiv API",
        role="exact-phrase rank-30 metadata query",
    ),
    "arxiv_elkies": SourceSpec(
        filename="arxiv_elkies.atom",
        url=(
            "https://export.arxiv.org/api/query?search_query=au:Elkies&"
            "start=0&max_results=100&sortBy=submittedDate&"
            "sortOrder=descending"
        ),
        sha256=(
            "49e6be35fa090d8e4009d8a6475bcab757672b44cf5c34f56de2a1d676738337"
        ),
        authority="official arXiv API",
        role="Noam Elkies author feed",
    ),
    "arxiv_klagsbrun": SourceSpec(
        filename="arxiv_klagsbrun.atom",
        url=(
            "https://export.arxiv.org/api/query?search_query=au:Klagsbrun&"
            "start=0&max_results=100&sortBy=submittedDate&"
            "sortOrder=descending"
        ),
        sha256=(
            "873ecbcc30109756790f9b69519575efa887c33dcdd099662a6f3a13a3c3fb83"
        ),
        authority="official arXiv API",
        role="Zev Klagsbrun author feed",
    ),
    "github_owner_repos": SourceSpec(
        filename="github_owner_repos.json",
        url=(
            "https://api.github.com/users/zevklagsbrun/repos?"
            "per_page=100&type=owner&sort=updated"
        ),
        sha256=(
            "463661ab1e803026a6ef0e256153fd6ce92179258cfa800b5afdb4f323c0ccc3"
        ),
        authority="official GitHub API; author's public account",
        role="public owner-repository inventory",
    ),
    "github_elkiescurve_tree": SourceSpec(
        filename="github_elkiescurve_tree.json",
        url=(
            "https://api.github.com/repos/zevklagsbrun/"
            "ElkiesCurve/git/trees/HEAD?recursive=1"
        ),
        sha256=(
            "a429c3dc5518325889f8f6b7abcdadb767738c9162cda51af29b2f3d9ec0ba59"
        ),
        authority="official GitHub API; author's public repository",
        role="recursive default-branch file inventory",
    ),
    "github_elkiescurve_branches": SourceSpec(
        filename="github_elkiescurve_branches.json",
        url=(
            "https://api.github.com/repos/zevklagsbrun/"
            "ElkiesCurve/branches?per_page=100"
        ),
        sha256=(
            "be6fff488590aec408fed7ca74d3619deb4ea089e84ddb27fbec589a5e966849"
        ),
        authority="official GitHub API; author's public repository",
        role="public branch inventory",
    ),
    "github_elkiescurve_tags": SourceSpec(
        filename="github_elkiescurve_tags.json",
        url=(
            "https://api.github.com/repos/zevklagsbrun/"
            "ElkiesCurve/tags?per_page=100"
        ),
        sha256=(
            "2ba33ca0557f1bb5b7ba88d67f9d0093c7185a36ec51fe2b7bd9372d3e001d6d"
        ),
        authority="official GitHub API; author's public repository",
        role="public tag inventory",
    ),
    "github_owner_gists": SourceSpec(
        filename="github_owner_gists.json",
        url="https://api.github.com/users/zevklagsbrun/gists?per_page=100",
        sha256=(
            "2ba33ca0557f1bb5b7ba88d67f9d0093c7185a36ec51fe2b7bd9372d3e001d6d"
        ),
        authority="official GitHub API; author's public account",
        role="public gist inventory",
    ),
    "zenodo_rank29": SourceSpec(
        filename="zenodo_rank29.json",
        url=(
            "https://zenodo.org/api/records?page=1&"
            "q=%22elliptic+curve%22+AND+%22rank+29%22&"
            "size=25&sort=bestmatch"
        ),
        sha256=(
            "a22dceaf5456db7a2327ad0d064261a1a161bb043538c9b9dd5c7c687d416c32"
        ),
        authority="official Zenodo API",
        role="rank-29 data-deposit query",
    ),
    "zenodo_rank30": SourceSpec(
        filename="zenodo_rank30.json",
        url=(
            "https://zenodo.org/api/records?page=1&"
            "q=%22elliptic+curve%22+AND+%22rank+30%22&"
            "size=25&sort=bestmatch"
        ),
        sha256=(
            "e23329b42f281ef222f8a497e07ac5bf3c5c3ff162589f8ac0ba301ec621f878"
        ),
        authority="official Zenodo API",
        role="rank-30 data-deposit query",
    ),
    "zenodo_elkies_klagsbrun": SourceSpec(
        filename="zenodo_elkies_klagsbrun.json",
        url=(
            "https://zenodo.org/api/records?page=1&"
            "q=Elkies+AND+Klagsbrun&size=25&sort=bestmatch"
        ),
        sha256=(
            "26fab1997e59f86d17722dd62f2bd6c91268be8b2267fa5bd962bfbca2426b0f"
        ),
        authority="official Zenodo API",
        role="author-pair data-deposit query",
    ),
    "elkies_k3_20si": SourceSpec(
        filename="elkies_k3_20si.html",
        url="https://people.math.harvard.edu/~elkies/K3_20SI.html",
        sha256=(
            "36ff78e6c487d4a8ab1f85d55c82fdad2667df063595302dd33566f73057d412"
        ),
        authority="Noam Elkies's Harvard author site",
        role="candidate explicit K3 data page found in the author archive",
    ),
    "elkies_k3_i18": SourceSpec(
        filename="elkies_k3_i18.html",
        url="https://people.math.harvard.edu/~elkies/M263.21/K3_I18.html",
        sha256=(
            "d340e3f919b59425728e2dd351ca6cc12cda2f736bd78830ff88e943c7c1193f"
        ),
        authority="Noam Elkies's Harvard author site",
        role="candidate explicit K3 page found in the author archive",
    ),
}


def _download(url: str, timeout: float) -> bytes:
    request = Request(
        url,
        headers={
            "User-Agent": "elliptic-rank30-public-source-audit/1.0",
            "Accept": "application/json, application/atom+xml, text/html, */*",
        },
    )
    with urlopen(request, timeout=timeout) as response:
        return response.read()


def verify_payload(name: str, payload: bytes) -> dict[str, Any]:
    """Verify one source snapshot before any artifact can be written."""

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
        inventory = verify_payload(name, payload)
        payloads[name] = payload
        inventories[name] = inventory
    return payloads, inventories


def _plain_html(payload: bytes) -> str:
    text = payload.decode("utf-8", errors="replace")
    text = re.sub(r"<br\s*/?>", "\n", text, flags=re.IGNORECASE)
    text = re.sub(r"<[^>]+>", " ", text)
    return re.sub(r"\s+", " ", html.unescape(text)).strip()


def _archive_subjects(payload: bytes) -> list[str]:
    source = payload.decode("utf-8", errors="replace")
    matches = re.findall(
        r'<div class="archive forcewrap"><a [^>]*>(.*?)</a></div>',
        source,
        flags=re.IGNORECASE | re.DOTALL,
    )
    return [_plain_html(match.encode()) for match in matches]


RANK_29_PATTERN = re.compile(
    r"(?:\bz\s*\^\s*\{?\s*29\s*\}?(?!\d)|"
    r"\brank(?:\s+(?:is|of))?\s*(?:at\s+least|>=|≥)?\s*29\b)",
    re.I,
)
RANK_30_PATTERN = re.compile(
    r"(?:\bz\s*\^\s*\{?\s*30\s*\}?(?!\d)|"
    r"\brank(?:\s+(?:is|of))?\s*(?:at\s+least|>=|≥)?\s*30\b)",
    re.I,
)


def audit_nmbrthry(index: bytes, announcement: bytes) -> dict[str, Any]:
    index_text = _plain_html(index)
    announcement_text = _plain_html(announcement)
    subjects = _archive_subjects(index)
    date_strings = re.findall(
        r"(?:Mon|Tue|Wed|Thu|Fri|Sat|Sun), \d{1,2} [A-Z][a-z]{2} "
        r"20\d{2} \d{2}:\d{2}:\d{2} [+-]\d{4}",
        index.decode("utf-8", errors="replace"),
    )
    parsed_dates: list[tuple[datetime, str]] = []
    for value in date_strings:
        try:
            parsed_dates.append((parsedate_to_datetime(value), value))
        except (TypeError, ValueError):
            continue
    latest = max(parsed_dates, key=lambda item: item[0])[1]
    rank29_subjects = [subject for subject in subjects if RANK_29_PATTERN.search(subject)]
    rank30_subjects = [subject for subject in subjects if RANK_30_PATTERN.search(subject)]
    announcement_markers = {
        "unconditional_rank_at_least_29": "has rank at least 29" in announcement_text,
        "rank17_fibration_named": "rank-17 fibration" in announcement_text,
        "twelve_extra_points_named": (
            "found 12 more independent points" in announcement_text
        ),
        "writeup_and_k3_computation_deferred": (
            "We intend to write up these results, including my computation "
            "of the K3 surface" in announcement_text
        ),
        "twenty_nine_x_coordinates_supplied": (
            "there are 29 independent points with X-coordinates"
            in announcement_text
        ),
        "gp_height_determinant_code_supplied": (
            "matdet(ellheightmatrix(E29,PT))" in announcement_text
        ),
    }
    if not all(announcement_markers.values()):
        missing = [key for key, value in announcement_markers.items() if not value]
        raise ValueError(f"rank-29 announcement markers missing: {missing}")
    return {
        "archive_subject_count": len(subjects),
        "latest_message_date_in_snapshot": latest,
        "rank29_subject_matches": rank29_subjects,
        "rank30_subject_matches": rank30_subjects,
        "rank30_phrase_occurrences_in_archive_text": len(
            RANK_30_PATTERN.findall(index_text)
        ),
        "announcement_semantic_markers": announcement_markers,
        "announcement_rank30_phrase_present": bool(
            RANK_30_PATTERN.search(announcement_text)
        ),
        "bounded_result": (
            "The pinned full-index snapshot contains the rank-29 announcement "
            "and its Cremona reply, but no subject or indexed excerpt matching "
            "the declared rank-30 patterns."
        ),
    }


def audit_dujella(rank_history: bytes, rank29: bytes) -> dict[str, Any]:
    history_source = rank_history.decode("utf-8", errors="replace")
    history_text = _plain_html(rank_history)
    rank29_source = rank29.decode("utf-8", errors="replace")
    point_indices = sorted(
        int(value)
        for value in re.findall(r"P<sub>(\d+)</sub>\s*=", rank29_source)
    )
    markers = {
        "current_record_rank_at_least_29": (
            "The current record is an example of elliptic curve with rank ≥ 29"
            in history_text
        ),
        "attributes_record_to_elkies_klagsbrun_2024": (
            "found by Elkies and Klagsbrun in 2024" in history_text
        ),
        "rank29_page_linked": 'href="rk29.html"' in history_source,
        "rank30_page_linked": bool(
            re.search(r'href=["\']rk30\.html["\']', history_source, re.I)
        ),
    }
    if not all(value for key, value in markers.items() if key != "rank30_page_linked"):
        raise ValueError("current rank-history markers were not recovered")
    if point_indices != list(range(1, 30)):
        raise ValueError("the pinned rank-29 page does not list P1 through P29")
    return {
        "rank_history_markers": markers,
        "rank29_page": {
            "weierstrass_model_present": all(
                token in rank29_source
                for token in (
                    "y<sup>2</sup> + xy = x<sup>3</sup>",
                    "27006183241630922218434652145297453784768054621836357954737385",
                    "55258058551342376475736699591118191821521067032535079608372404779149413277716173425636721497",
                )
            ),
            "point_indices": point_indices,
            "point_count": len(point_indices),
        },
        "bounded_result": (
            "The pinned current catalogue identifies rank >=29 as the record "
            "and publishes exactly P1,...,P29; it has no rk30.html link."
        ),
    }


ATOM = "http://www.w3.org/2005/Atom"
OPENSEARCH = "http://a9.com/-/spec/opensearch/1.1/"


def parse_arxiv_atom(payload: bytes) -> dict[str, Any]:
    root = ET.fromstring(payload)
    entries = []
    for entry in root.findall(f"{{{ATOM}}}entry"):
        entries.append(
            {
                "id": entry.findtext(f"{{{ATOM}}}id"),
                "title": " ".join(
                    (entry.findtext(f"{{{ATOM}}}title") or "").split()
                ),
                "published": entry.findtext(f"{{{ATOM}}}published"),
                "updated": entry.findtext(f"{{{ATOM}}}updated"),
            }
        )
    total_text = root.findtext(f"{{{OPENSEARCH}}}totalResults")
    return {
        "feed_title": root.findtext(f"{{{ATOM}}}title"),
        "feed_updated": root.findtext(f"{{{ATOM}}}updated"),
        "total_results": int(total_text or 0),
        "entry_count_in_snapshot": len(entries),
        "entries": entries,
    }


def _post_announcement_entries(feed: dict[str, Any]) -> list[dict[str, Any]]:
    cutoff = "2024-08-29T00:00:00Z"
    return [
        entry
        for entry in feed["entries"]
        if max(entry.get("published") or "", entry.get("updated") or "") >= cutoff
    ]


def _record_related(entries: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [
        entry
        for entry in entries
        if "rank" in entry["title"].lower()
        and "elliptic" in entry["title"].lower()
    ]


def audit_arxiv(
    rank30: bytes, elkies: bytes, klagsbrun: bytes
) -> dict[str, Any]:
    rank30_feed = parse_arxiv_atom(rank30)
    elkies_feed = parse_arxiv_atom(elkies)
    klagsbrun_feed = parse_arxiv_atom(klagsbrun)
    if rank30_feed["total_results"] != 0:
        raise ValueError("the pinned exact-phrase arXiv query is no longer empty")
    elkies_post = _post_announcement_entries(elkies_feed)
    klagsbrun_post = _post_announcement_entries(klagsbrun_feed)
    return {
        "exact_phrase_query": {
            key: value for key, value in rank30_feed.items() if key != "entries"
        },
        "author_feeds": {
            "elkies": {
                "total_results": elkies_feed["total_results"],
                "published_or_updated_after_rank29_announcement": elkies_post,
                "post_announcement_record_related_titles": _record_related(
                    elkies_post
                ),
            },
            "klagsbrun": {
                "total_results": klagsbrun_feed["total_results"],
                "published_or_updated_after_rank29_announcement": klagsbrun_post,
                "post_announcement_record_related_titles": _record_related(
                    klagsbrun_post
                ),
            },
        },
        "scope_limit": (
            "The exact-phrase query and two author feeds are metadata checks, "
            "not an exhaustive search of every arXiv paper or every wording."
        ),
    }


def audit_github(
    repos_payload: bytes,
    tree_payload: bytes,
    branches_payload: bytes,
    tags_payload: bytes,
    gists_payload: bytes,
) -> dict[str, Any]:
    repos = json.loads(repos_payload)
    tree_document = json.loads(tree_payload)
    branches = json.loads(branches_payload)
    tags = json.loads(tags_payload)
    gists = json.loads(gists_payload)
    paths = sorted(item["path"] for item in tree_document["tree"])
    rank_numbers = []
    for path in paths:
        rank_numbers.extend(
            int(value)
            for value in re.findall(
                r"(?:rank|rels|points|verify)[^0-9]{0,4}(\d+)", path, re.I
            )
        )
    terms = ("rank29", "rank30", "e29", "e30", "fibration", "k3")
    term_hits = {
        term: [path for path in paths if term in path.lower()] for term in terms
    }
    return {
        "owner_repository_count": len(repos),
        "owner_repositories": [
            {
                "full_name": repo["full_name"],
                "html_url": repo["html_url"],
                "updated_at": repo["updated_at"],
                "default_branch": repo["default_branch"],
            }
            for repo in repos
        ],
        "elkiescurve_repository": {
            "tree_truncated": tree_document["truncated"],
            "path_count": len(paths),
            "paths": paths,
            "branches": [branch["name"] for branch in branches],
            "tags": [tag["name"] for tag in tags],
            "largest_rank_number_in_path_names": max(rank_numbers),
            "rank_numbers_in_path_names": sorted(set(rank_numbers)),
            "candidate_term_hits": term_hits,
        },
        "public_gist_count": len(gists),
        "bounded_result": (
            "The author's public GitHub account has one relevant repository; "
            "its sole public branch contains rank-labelled files up to 28, "
            "but no rank-29/rank-30 or K3/fibration-named file."
        ),
    }


def audit_zenodo(
    rank29: bytes, rank30: bytes, authors: bytes
) -> dict[str, Any]:
    reports = {}
    for name, payload in (
        ("rank29_query", rank29),
        ("rank30_query", rank30),
        ("author_pair_query", authors),
    ):
        document = json.loads(payload)
        reports[name] = {
            "self_url": document["links"]["self"],
            "total_results": document["hits"]["total"],
        }
    return {
        "queries": reports,
        "data_deposit_recovered": any(
            item["total_results"] for item in reports.values()
        ),
        "scope_limit": "Three declared metadata queries, not all possible keywords.",
    }


def audit_author_k3_pages(k3_20si: bytes, k3_i18: bytes) -> dict[str, Any]:
    shioda_text = _plain_html(k3_20si)
    i18_text = _plain_html(k3_i18)

    def record_markers(text: str) -> dict[str, bool]:
        return {
            "mordell_weil_rank17": bool(
                re.search(r"(?:Mordell-Weil|MW)\s+rank[- ]?17\b", text, re.I)
            ),
            "record_lattice_discriminant_948": bool(
                re.search(
                    r"(?:\(\s*19\s*,\s*948\s*\)|"
                    r"disc(?:riminant)?[^0-9]{0,20}948\b)",
                    text,
                    re.I,
                )
            ),
        }

    return {
        "K3_20SI": {
            "page_title": "Mordell-Weil generators for singular Shioda-Inose surfaces over Q",
            "explicitly_infinite_cyclic_mordell_weil": (
                "infinite cyclic Mordell-Weil group" in shioda_text
            ),
            "record_rank17_markers_present": record_markers(shioda_text),
            "classification": (
                "Explicit singular Shioda-Inose rank-one fibrations; not the "
                "record rank-17 Mordell-Weil fibration."
            ),
        },
        "K3_I18": {
            "i18_fibre_marker": "I_{18}" in i18_text or "I18" in i18_text,
            "root_lattice_parameter_count_marker": (
                "20 - (2+17) = 1" in i18_text
            ),
            "record_rank17_markers_present": record_markers(i18_text),
            "classification": (
                "An explicit I18 reducible-fibre construction, where 17 is "
                "the fibre root-lattice contribution; not the record rank-17 "
                "Mordell-Weil fibration."
            ),
        },
    }


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def local_frontier() -> dict[str, Any]:
    repository = Path(__file__).resolve().parents[2]
    rank29_path = (
        repository
        / "artifacts/generated-results/"
        "elliptic_elkies_klagsbrun_rank29_certificate.json"
    )
    rank18_path = (
        repository
        / "artifacts/generated-results/elliptic_elkies_rank18_source_audit.json"
    )
    rank29 = _load_json(rank29_path)
    rank18 = _load_json(rank18_path)
    search_paths = sorted(
        (repository / "artifacts/generated-results").glob(
            "elliptic_elkies_klagsbrun_rank30*.json"
        )
    )
    searches = []
    for path in search_paths:
        document = _load_json(path)
        searches.append(
            {
                "path": str(path.relative_to(repository)),
                "sha256": sha256(path.read_bytes()).hexdigest(),
                "status": document["status"],
            }
        )
    if rank29["status"] != "exact_unconditional_rank_at_least_29":
        raise ValueError("local rank-29 certificate status changed")
    rank18_assessment = rank18["elkies_rank18"]["assessment"]
    if rank18_assessment["explicit_rank18_family_reproducible"]:
        raise ValueError("local rank-18 source blocker unexpectedly disappeared")
    return {
        "rank29_exact_certificate": {
            "path": str(rank29_path.relative_to(repository)),
            "sha256": sha256(rank29_path.read_bytes()).hexdigest(),
            "status": rank29["status"],
            "certified_algebraic_rank_lower_bound": rank29["claim"][
                "certified_algebraic_rank_lower_bound"
            ],
        },
        "rank18_source_audit": {
            "path": str(rank18_path.relative_to(repository)),
            "sha256": sha256(rank18_path.read_bytes()).hexdigest(),
            "status": rank18_assessment["status"],
            "explicit_rank18_family_reproducible": rank18_assessment[
                "explicit_rank18_family_reproducible"
            ],
        },
        "existing_rank30_bounded_searches": searches,
    }


def build_report(
    payloads: dict[str, bytes], inventories: dict[str, dict[str, Any]]
) -> dict[str, Any]:
    nmbrthry = audit_nmbrthry(
        payloads["nmbrthry_index"], payloads["nmbrthry_rank29"]
    )
    dujella = audit_dujella(
        payloads["dujella_rank_history"], payloads["dujella_rank29"]
    )
    arxiv = audit_arxiv(
        payloads["arxiv_rank30"],
        payloads["arxiv_elkies"],
        payloads["arxiv_klagsbrun"],
    )
    github = audit_github(
        payloads["github_owner_repos"],
        payloads["github_elkiescurve_tree"],
        payloads["github_elkiescurve_branches"],
        payloads["github_elkiescurve_tags"],
        payloads["github_owner_gists"],
    )
    zenodo = audit_zenodo(
        payloads["zenodo_rank29"],
        payloads["zenodo_rank30"],
        payloads["zenodo_elkies_klagsbrun"],
    )
    author_k3 = audit_author_k3_pages(
        payloads["elkies_k3_20si"], payloads["elkies_k3_i18"]
    )
    local = local_frontier()
    required = {
        "rank30_weierstrass_model": False,
        "thirty_rational_points_on_one_curve": False,
        "exact_point_membership_checks_for_thirty_points": False,
        "independence_certificate_for_thirty_points": False,
        "record_rank17_k3_weierstrass_coefficients": False,
        "record_rank17_k3_section_coordinates_1_through_17": False,
        "record_rank17_specialization_parameter_or_map": False,
        "public_record_search_data_beyond_the_rank29_announcement": False,
    }
    return {
        "schema": "elliptic-rank30-public-source-audit-v1",
        "artifact_kind": "bounded_primary_source_recovery_audit",
        "status": "source_audit_no_public_rank30_data_recovered",
        "audit_date": AUDIT_DATE,
        "script_sha256": sha256(Path(__file__).resolve().read_bytes()).hexdigest(),
        "software": {"python": platform.python_version()},
        "reproduce_command": REPRODUCE_COMMAND,
        "declared_scope": {
            "included": [
                "official NMBRTHRY full-index snapshot and primary announcement",
                "Dujella's author-maintained current rank catalogue",
                "official arXiv exact-phrase query and both authors' feeds",
                "Zev Klagsbrun's public GitHub owner repos, branches, tags, and gists",
                "official Zenodo metadata queries",
                "two explicit K3 pages located on Noam Elkies's author site",
                "the existing exact local rank-29 and rank-18 source audits",
            ],
            "excluded_or_not_proved": [
                "private, unindexed, access-controlled, or unpublished data",
                "all conceivable keyword variants and every repository on the web",
                "a mathematical nonexistence theorem for rank >=30",
            ],
        },
        "source_inventories": inventories,
        "source_results": {
            "nmbrthry": nmbrthry,
            "dujella": dujella,
            "arxiv": arxiv,
            "github": github,
            "zenodo": zenodo,
            "author_k3_pages": author_k3,
        },
        "local_existing_work": local,
        "recovered_frontier": {
            "highest_publicly_reproducible_algebraic_rank_lower_bound": 29,
            "rank29_exact_local_certificate_available": True,
            "rank30_curve_recovered": False,
            "explicit_record_rank17_k3_family_recovered": False,
            "published_search_mechanism_recovered": (
                "rank-17 fibration of the same K3 used for rank 28, then sieve "
                "specializations and search for points outside generic Z^17"
            ),
            "mechanism_replayable_from_public_data": False,
        },
        "required_rank30_reproduction_objects": required,
        "assessment": {
            "all_pinned_source_hashes_verified_before_write": all(
                source["matches_pinned_revision"]
                for source in inventories.values()
            ),
            "public_rank30_or_higher_curve_reproducible": False,
            "public_record_k3_search_family_reproducible": False,
            "breakthrough_result_found": False,
            "reason": (
                "The declared official/author-controlled source scope still "
                "ends at a fully reproducible rank-29 lower bound.  The primary "
                "announcement identifies a rank-17 K3 fibration and explicitly "
                "defers its writeup/computation; no audited source supplies that "
                "model, its 17 sections, a specialization map, a 30th point, or "
                "a rank-30 independence certificate."
            ),
            "epistemic_status": (
                "bounded negative source-recovery result; not evidence that a "
                "rank-30 curve does not exist"
            ),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source-dir",
        type=Path,
        help="read all pinned snapshots by SourceSpec.filename instead of downloading",
    )
    parser.add_argument("--download-timeout", type=float, default=20.0)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    arguments = parser.parse_args()
    if arguments.download_timeout <= 0 or arguments.download_timeout > 60:
        raise ValueError("download timeout must lie in (0,60]")
    payloads, inventories = load_sources(
        source_dir=arguments.source_dir, timeout=arguments.download_timeout
    )
    report = build_report(payloads, inventories)
    encoded = json.dumps(report, indent=2, sort_keys=True) + "\n"
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(encoded, encoding="utf-8")
    print(f"wrote {arguments.output}")
    print("public reproducible frontier: rank >= 29")
    print("rank >= 30 recovered: false")
    print("record rank-17 K3 model recovered: false")


if __name__ == "__main__":
    main()
