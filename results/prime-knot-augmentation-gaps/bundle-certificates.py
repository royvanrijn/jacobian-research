"""Package the exact inputs and replay scripts cited by the paper."""

from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile, ZipInfo

PAPER = Path(__file__).resolve().parent
PROOFS = PAPER.parents[1] / "research" / "meridional-rank" / "representations"
ARCHIVE = PAPER / "prime_knot_certificates.zip"

CERTIFICATES = (
    "protected_five_bridge_uniform_family",
    "protected_six_bridge_uniform_family",
    "protected_universal_parent_candidate",
    "threaded_curl_uniform_transfer",
)
IMPORTS = (
    "protected_five_bridge_prime_inputs.json",
    "verify_protected_five_bridge_prime.py",
    "verify_new_seventeen_crossing_knot.py",
    "verify_unbounded_prime_augmentation_gap_audit.py",
    "sources/uniform5-snappy332-init.py",
    "sources/uniform5-spherogram241-invariants.py",
    "sources/uniform5-spherogram241-links-base.py",
)


def main():
    files = {"README.md": PAPER / "certificate-replay.md"}
    for stem in CERTIFICATES:
        files[f"checks/{stem}_inputs.json"] = PROOFS / f"{stem}_inputs.json"
        files[f"checks/verify_{stem}.py"] = PROOFS / f"verify_{stem}.py"
    for name in IMPORTS:
        files[f"checks/{name}"] = PROOFS / name

    with ZipFile(ARCHIVE, "w") as archive:
        for name, source in sorted(files.items()):
            member = ZipInfo(f"prime-knot-certificates/{name}")
            member.external_attr = 0o644 << 16
            archive.writestr(member, source.read_bytes(), compress_type=ZIP_DEFLATED)


if __name__ == "__main__":
    main()
