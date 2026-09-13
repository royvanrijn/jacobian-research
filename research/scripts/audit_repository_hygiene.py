#!/usr/bin/env python3
"""Audit active-note identities and repository-only generated debris.

Already-tracked files in ``artifacts/generated-results`` are deliberate pinned
certificates.  The directory remains ignored so a local replay cannot add a
large output set accidentally; tracked pinned files are therefore the one
allowed tracked/ignored class.
"""

from __future__ import annotations

import re
import subprocess
from collections import defaultdict
from pathlib import Path
from research_documents import document_paths


ROOT = Path(__file__).resolve().parents[1]
REPOSITORY_ROOT = next(
    parent for parent in (ROOT, *ROOT.parents) if (parent / ".git").exists()
)
PROJECT_PREFIX = ROOT.relative_to(REPOSITORY_ROOT).as_posix()
EXCLUDED_PARTS = {
    ".cache",
    ".git",
    ".idea",
    ".lake",
    ".venv",
    "archive",
    "output",
    "tmp",
}
DEFINITION_PATTERNS = (
    re.compile(
        r"\*\*(?:Theorem|Corollary|Proposition|Lemma)\s+`?"
        r"([A-Z][A-Z0-9-]*\d[A-Z0-9-]*)"
    ),
    re.compile(r"\*\*`?([A-Z][A-Z0-9-]*\d[A-Z0-9-]*)`?\.\*\*"),
)


definitions: dict[str, set[Path]] = defaultdict(set)
for document in document_paths():
    if ROOT not in document.parents:
        continue
    relative = document.relative_to(ROOT)
    if EXCLUDED_PARTS.intersection(relative.parts):
        continue
    for line in document.read_text().splitlines():
        for pattern in DEFINITION_PATTERNS:
            match = pattern.search(line)
            if match:
                definitions[match.group(1)].add(relative)
                break

duplicates = {
    identifier: paths
    for identifier, paths in definitions.items()
    if len(paths) > 1
}
if duplicates:
    details = []
    for identifier, paths in sorted(duplicates.items()):
        details.append(f"{identifier}: {', '.join(map(str, sorted(paths)))}")
    raise SystemExit(
        "headline theorem identifiers are defined in multiple active notes:\n"
        + "\n".join(details)
    )

tracked_ignored: list[str] = []
tracked_files: list[str] = []
result = subprocess.run(
    ["git", "ls-files", "-ci", "--exclude-standard"],
    cwd=REPOSITORY_ROOT,
    check=True,
    capture_output=True,
    text=True,
)
tracked_ignored = result.stdout.splitlines()
tracked_files = subprocess.run(
    ["git", "ls-files"],
    cwd=REPOSITORY_ROOT,
    check=True,
    capture_output=True,
    text=True,
).stdout.splitlines()
generated_results_prefix = f"{PROJECT_PREFIX}/artifacts/generated-results/"
unexpected_tracked_ignored = [
    path for path in tracked_ignored
    if not path.startswith(generated_results_prefix)
]
if unexpected_tracked_ignored:
    raise SystemExit(
        "ignored files outside the pinned generated-results tree are still tracked:\n"
        + "\n".join(unexpected_tracked_ignored)
    )

# Archive-only certificate citations in active Mestre notes must be explicit.
# Otherwise an agent can mistake historical evidence for a missing active output
# and repeat a completed experiment.
archive_output_citation_failures: list[str] = []
archive_output_citation_checks = 0
active_ec_outputs = ROOT / "artifacts/generated-results/elliptic-curves"
archived_ec_outputs = ROOT / "archive/elliptic-curves/artifacts/generated-results"
for document in ROOT.glob("elliptic-curves/notes/*.md"):
    for number, line in enumerate(document.read_text().splitlines(), 1):
        for name in set(re.findall(r"`([^`/]+\.json(?:\.gz)?)`", line)):
            if (active_ec_outputs / name).is_file() or not (archived_ec_outputs / name).is_file():
                continue
            archive_output_citation_checks += 1
            target = f"../../archive/elliptic-curves/artifacts/generated-results/{name}"
            if target not in line:
                archive_output_citation_failures.append(
                    f"{document.relative_to(ROOT)}:{number}: {name} is archive-only; link its retained archive path"
                )
if archive_output_citation_failures:
    raise SystemExit(
        "archive-only EC output citations must link directly to historical evidence:\n"
        + "\n".join(archive_output_citation_failures)
    )

# Completed proof notes may retain their experiment history, but the exact
# headings below used to make a finished or unscheduled route look like the
# programme's current execution instruction. Keep that wording out of all
# active EC/K3 notes, including exploratory references that are not canonical
# status sources; a source may instead state an explicit historical or
# separately scoped boundary. This deliberately does not reject a genuinely
# specified open theorem gate with different wording.
stale_execution_heading = re.compile(
    r"^#{1,6}\s+(?:"
    r"next experiment|next decisive computation|"
    r"authorized(?:\s+[\w-]+){0,2}\s+continuation|"
    r"(?:interpretation and|actionable|boundary and|ranked interpretation and|"
    r"revised mechanism ranking and)\s+next experiment"
    r")\s*$",
    re.IGNORECASE,
)
active_ec_k3_notes = []
for document in document_paths():
    if ROOT not in document.parents:
        continue
    relative = document.relative_to(ROOT)
    if not relative.parts or relative.parts[0] not in {"elliptic-curves", "elkies-k3"}:
        continue
    if EXCLUDED_PARTS.intersection(relative.parts):
        continue
    active_ec_k3_notes.append(document)
stale_execution_headings: list[str] = []
for document in sorted(active_ec_k3_notes):
    for number, line in enumerate(document.read_text().splitlines(), 1):
        if stale_execution_heading.fullmatch(line.strip()):
            stale_execution_headings.append(f"{document.relative_to(ROOT)}:{number}: {line.strip()}")
if stale_execution_headings:
    raise SystemExit(
        "active EC/K3 notes must not present stale execution headings:\n"
        + "\n".join(stale_execution_headings)
    )

# These dated handoffs retain reproducible commands and intermediate research
# history.  Their current-boundary labels are therefore part of the navigation
# contract: losing one could make an historical command look like a current
# compute request.  Keep this intentionally small and explicit rather than
# applying a generic word ban to legitimate open theorem gates.
historical_handoff_markers = {
    "elkies-k3/RESEARCH_UPDATE_2026-08-20.md": (
        "research update (2026-08-20, historical)",
        "direct rootless endpoint gate are now complete. Use the current programme",
    ),
    "elkies-k3/RESEARCH_UPDATE_2026-08-22.md": (
        "Research update — 2026-08-22 (historical)",
        "not the current K3 queue",
    ),
    "elkies-k3/Q80_LOW_Q_ALTERNATE_2026-08-22.md": (
        "Q80 low-q alternate corridor — 2026-08-22 (historical)",
        "not the current\nconstruction queue",
    ),
    "elkies-k3/LOWER_Q_FIBRATION_PATH_SEARCH_2026-08-20.md": (
        "historical",
        "not a current work request",
    ),
    "elkies-k3/E6_MW3_PROGRESS_2026-08-20.md": (
        "Historical snapshot",
        "Historical search/reconstruction strategy (not scheduled)",
    ),
    "elkies-k3/E6_MW3_ATTACK.md": (
        "historical",
        "Historical reproducible pipeline (not a current work request)",
    ),
    "elkies-k3/R17_RATIONAL_QUADRATIC_MW20_SEARCH_2026-09-04.md": (
        "Operational continuation is paused",
        "Remaining theorem gate (unscheduled)",
    ),
    "elkies-k3/RANK7_AUXILIARY_CATALOGUE_2026-09-01.md": (
        "Remaining catalogue frontier (unscheduled)",
        "No whole-backend or factory enumeration is current work.",
    ),
    "elkies-k3/NS0024_MW4_MODULAR_RECOVERY_2026-09-01.md": (
        "historical geometric frontier",
        "Historical extension-field question (not scheduled)",
    ),
    "elkies-k3/NS0024_NEW_ROOTLESS_SOURCE_ROUTE_2026-09-03.md": (
        "historical geometric route",
        "not a QQ arithmetic continuation",
    ),
    "elkies-k3/LATTICE_FOUNDRY_EMPIRICAL_SOURCE_RANKING_2026-09-02.md": (
        "historical scoring snapshot",
        "not a current construction queue",
    ),
    "elkies-k3/NS0024_EDGE1_COMPILER_PREPARATION_2026-09-01.md": (
        "historical geometric regression",
        "Historical reproduction (not a current work request)",
    ),
    "elliptic-curves/notes/FULL11952_LATE_BAND_TRIAL_2026-09-06.md": (
        "Historical second-prime-band trial",
        "This completed dated trial retains",
    ),
    "elliptic-curves/notes/MW16_OUTER_PARAMETER_BANDS_2026-09-06.md": (
        "Historical MW16 outer-parameter campaign",
        "This completed campaign retains",
    ),
    "elkies-k3/OTHER_RANK17_GATE_B_AUDIT_2026-08-31.md": (
        "historical Gate-B audit",
        "does not select a construction route",
    ),
    "elkies-k3/Q80_TO_ROOTLESS_PATH_2026-08-21.md": (
        "Historical generic q80-to-rootless MW17 path",
        "not the current equation",
    ),
    "elliptic-curves/notes/CURVE385_ITERATED_HALF_LATTICE_RECOVERY_2026-09-04.md": (
        "historical completed campaign",
        "sparse rank-32 no-growth campaign are sealed evidence.",
    ),
    "elkies-k3/E6_P2_REDUCTION_2026-08-20.md": (
        "historical GF(31) search state",
        "not a current construction queue.",
    ),
    "elkies-k3/GOLAY_OCTAD_LATTICE_DESIGN_2026-09-01.md": (
        "historical geometric design",
        "outside the arithmetic MW17 queue",
    ),
}
historical_handoff_failures: list[str] = []
for relative, markers in historical_handoff_markers.items():
    handoff = ROOT / relative
    content = handoff.read_text()
    for marker in markers:
        if marker not in content:
            historical_handoff_failures.append(
                f"{relative}: missing historical-boundary marker {marker!r}"
            )
if historical_handoff_failures:
    raise SystemExit(
        "historical EC/K3 handoffs must preserve their current-boundary labels:\n"
        + "\n".join(historical_handoff_failures)
    )

backup_pattern = re.compile(r"(?:~|\.orig|\.rej|\.bak(?:[-.]|$))")
active_backups = [
    path
    for path in tracked_files
    if "archive" not in Path(path).parts and backup_pattern.search(Path(path).name)
]
if active_backups:
    raise SystemExit(
        "backup snapshots must be moved under an explicit archive directory:\n"
        + "\n".join(active_backups)
    )

print(
    f"PASS repository hygiene: {len(definitions)} active headline identifiers are "
    f"file-unique; {len(tracked_ignored)} pinned generated artifacts are tracked "
    f"under the ignore guard; {archive_output_citation_checks} archive-only EC output citations "
    f"link directly to retained evidence; {len(historical_handoff_markers)} historical handoffs preserve "
    "their current-boundary labels; no stale active execution headings or other ignored files or "
    "active-tree backup snapshots are tracked"
)
