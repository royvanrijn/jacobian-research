from __future__ import annotations

from io import BytesIO
import gzip
from hashlib import sha256
import json
from pathlib import Path
import sys
import tarfile
import unittest


CAS_DIRECTORY = Path(__file__).resolve().parents[1] / "cas"
if str(CAS_DIRECTORY) not in sys.path:
    sys.path.insert(0, str(CAS_DIRECTORY))

from audit_elkies_rank18_sources import (  # noqa: E402
    ELKIES_0709_SOURCE_SHA256,
    ELKIES_0802_SOURCE_SHA256,
    unpack_arxiv_source,
)


ARTIFACT = (
    Path(__file__).resolve().parents[2]
    / "artifacts"
    / "generated-results"
    / "elliptic_elkies_rank18_source_audit.json"
)
SCRIPT = CAS_DIRECTORY / "audit_elkies_rank18_sources.py"


class ElkiesRank18SourceAuditTests(unittest.TestCase):
    def test_unpack_handles_single_tex_and_tarred_arxiv_sources(self) -> None:
        single = gzip.compress(b"single tex source")
        self.assertEqual(
            unpack_arxiv_source(single, single_tex_name="paper.tex"),
            {"paper.tex": b"single tex source"},
        )

        buffer = BytesIO()
        with tarfile.open(fileobj=buffer, mode="w") as archive:
            for name, payload in (
                ("paper.tex", b"tex"),
                ("style.cls", b"class"),
            ):
                info = tarfile.TarInfo(name=name)
                info.size = len(payload)
                archive.addfile(info, BytesIO(payload))
        packed = gzip.compress(buffer.getvalue())
        self.assertEqual(
            unpack_arxiv_source(packed, single_tex_name="unused.tex"),
            {"paper.tex": b"tex", "style.cls": b"class"},
        )

    def test_pinned_artifact_records_the_source_level_blocker(self) -> None:
        report = json.loads(ARTIFACT.read_text(encoding="utf-8"))
        audit = report["elkies_rank18"]
        source_0709 = audit["sources"]["arxiv_0709_2908"]
        source_0802 = audit["sources"]["arxiv_0802_1301"]
        self.assertEqual(
            source_0709["inventory"]["source_sha256"],
            ELKIES_0709_SOURCE_SHA256,
        )
        self.assertEqual(
            source_0802["inventory"]["source_sha256"],
            ELKIES_0802_SOURCE_SHA256,
        )
        self.assertEqual(
            [member["name"] for member in source_0709["inventory"]["members"]],
            ["mfo07_arxiv.tex", "owrart.cls"],
        )
        self.assertEqual(source_0709["inventory"]["non_tex_data_members"], [])
        self.assertEqual(source_0802["inventory"]["non_tex_data_members"], [])
        self.assertTrue(all(source_0709["semantic_markers"].values()))
        self.assertTrue(all(source_0802["semantic_markers"].values()))
        self.assertFalse(
            any(audit["required_rank18_reproduction_objects"].values())
        )
        self.assertFalse(audit["assessment"]["explicit_rank18_family_reproducible"])
        self.assertEqual(
            audit["assessment"]["status"], "source-blocked-before-instantiation"
        )
        self.assertFalse(audit["assessment"]["specialization_screen_performed"])

    def test_reproduction_metadata_is_pinned(self) -> None:
        report = json.loads(ARTIFACT.read_text(encoding="utf-8"))
        self.assertEqual(
            report["script_sha256"], sha256(SCRIPT.read_bytes()).hexdigest()
        )
        self.assertEqual(
            report["software"], {"python": "3.14.6", "sympy": "1.14.0"}
        )

    def test_fallback_artifact_records_exact_identity_scope(self) -> None:
        report = json.loads(ARTIFACT.read_text(encoding="utf-8"))
        fallback = report["explicit_fallback"]
        self.assertEqual(fallback["published_generic_rank_lower_bound"], 14)
        self.assertEqual(fallback["published_group_origin"], "P15")
        self.assertTrue(
            fallback["symbolic_verification"]["all_fifteen_sections_exact"]
        )
        self.assertTrue(fallback["exact_t2_verification"]["all_points_exact"])
        self.assertIn("not independently certified", fallback["local_status"])
        self.assertFalse(report["result"]["breakthrough_curve_found"])


if __name__ == "__main__":
    unittest.main()
