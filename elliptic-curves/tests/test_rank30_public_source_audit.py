from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
import platform
import sys
import unittest


CAS_DIRECTORY = Path(__file__).resolve().parents[1] / "cas"
if str(CAS_DIRECTORY) not in sys.path:
    sys.path.insert(0, str(CAS_DIRECTORY))

from audit_rank30_public_sources import (  # noqa: E402
    RANK_30_PATTERN,
    SOURCE_SPECS,
    parse_arxiv_atom,
    verify_payload,
)


REPOSITORY = Path(__file__).resolve().parents[2]
ARTIFACT = (
    REPOSITORY
    / "artifacts/generated-results/elliptic_rank30_public_source_audit.json"
)
SCRIPT = CAS_DIRECTORY / "audit_rank30_public_sources.py"


class Rank30PublicSourceAuditTests(unittest.TestCase):
    def test_pinned_source_hashes_and_reproduction_metadata(self) -> None:
        report = json.loads(ARTIFACT.read_text(encoding="utf-8"))
        inventories = report["source_inventories"]
        self.assertEqual(set(inventories), set(SOURCE_SPECS))
        for name, spec in SOURCE_SPECS.items():
            self.assertEqual(inventories[name]["sha256"], spec.sha256)
            self.assertTrue(inventories[name]["matches_pinned_revision"])
            self.assertEqual(inventories[name]["url"], spec.url)
        self.assertEqual(
            report["script_sha256"], sha256(SCRIPT.read_bytes()).hexdigest()
        )
        self.assertEqual(report["software"], {"python": platform.python_version()})
        self.assertTrue(
            report["assessment"]["all_pinned_source_hashes_verified_before_write"]
        )

    def test_artifact_pins_the_public_frontier_and_exact_blockers(self) -> None:
        report = json.loads(ARTIFACT.read_text(encoding="utf-8"))
        self.assertEqual(
            report["status"], "source_audit_no_public_rank30_data_recovered"
        )
        frontier = report["recovered_frontier"]
        self.assertEqual(
            frontier["highest_publicly_reproducible_algebraic_rank_lower_bound"],
            29,
        )
        self.assertTrue(frontier["rank29_exact_local_certificate_available"])
        self.assertFalse(frontier["rank30_curve_recovered"])
        self.assertFalse(frontier["explicit_record_rank17_k3_family_recovered"])
        self.assertFalse(frontier["mechanism_replayable_from_public_data"])
        self.assertFalse(
            any(report["required_rank30_reproduction_objects"].values())
        )
        self.assertFalse(report["assessment"]["breakthrough_result_found"])
        self.assertIn(
            "bounded negative source-recovery result",
            report["assessment"]["epistemic_status"],
        )

    def test_official_archive_and_catalogue_have_rank29_but_no_rank30_hit(self) -> None:
        report = json.loads(ARTIFACT.read_text(encoding="utf-8"))
        source_results = report["source_results"]
        nmbrthry = source_results["nmbrthry"]
        self.assertEqual(nmbrthry["archive_subject_count"], 252)
        self.assertEqual(
            nmbrthry["latest_message_date_in_snapshot"],
            "Thu, 16 Jul 2026 10:09:54 +0000",
        )
        self.assertEqual(
            nmbrthry["rank29_subject_matches"],
            ["Z^29 in E(Q)", "Re: Z^29 in E(Q)"],
        )
        self.assertEqual(nmbrthry["rank30_subject_matches"], [])
        self.assertEqual(nmbrthry["rank30_phrase_occurrences_in_archive_text"], 0)
        self.assertTrue(all(nmbrthry["announcement_semantic_markers"].values()))

        dujella = source_results["dujella"]
        self.assertEqual(dujella["rank29_page"]["point_count"], 29)
        self.assertEqual(
            dujella["rank29_page"]["point_indices"], list(range(1, 30))
        )
        self.assertTrue(dujella["rank29_page"]["weierstrass_model_present"])
        self.assertFalse(
            dujella["rank_history_markers"]["rank30_page_linked"]
        )

    def test_archive_deposit_and_author_code_queries_record_no_new_data(self) -> None:
        report = json.loads(ARTIFACT.read_text(encoding="utf-8"))
        source_results = report["source_results"]
        arxiv = source_results["arxiv"]
        self.assertEqual(arxiv["exact_phrase_query"]["total_results"], 0)
        for author in ("elkies", "klagsbrun"):
            self.assertEqual(
                arxiv["author_feeds"][author][
                    "post_announcement_record_related_titles"
                ],
                [],
            )

        github = source_results["github"]
        repository = github["elkiescurve_repository"]
        self.assertEqual(repository["branches"], ["master"])
        self.assertEqual(repository["tags"], [])
        self.assertEqual(repository["rank_numbers_in_path_names"], [20, 21, 22, 23, 24, 27, 28])
        self.assertEqual(repository["largest_rank_number_in_path_names"], 28)
        self.assertTrue(
            all(not paths for paths in repository["candidate_term_hits"].values())
        )
        self.assertEqual(github["public_gist_count"], 0)

        zenodo = source_results["zenodo"]
        self.assertFalse(zenodo["data_deposit_recovered"])
        self.assertTrue(
            all(
                query["total_results"] == 0
                for query in zenodo["queries"].values()
            )
        )

    def test_author_k3_pages_are_not_the_record_fibration(self) -> None:
        report = json.loads(ARTIFACT.read_text(encoding="utf-8"))
        pages = report["source_results"]["author_k3_pages"]
        self.assertTrue(pages["K3_20SI"]["explicitly_infinite_cyclic_mordell_weil"])
        self.assertFalse(
            any(pages["K3_20SI"]["record_rank17_markers_present"].values())
        )
        self.assertTrue(pages["K3_I18"]["i18_fibre_marker"])
        self.assertTrue(pages["K3_I18"]["root_lattice_parameter_count_marker"])
        self.assertFalse(
            any(pages["K3_I18"]["record_rank17_markers_present"].values())
        )

    def test_local_existing_work_is_not_duplicated_or_overstated(self) -> None:
        report = json.loads(ARTIFACT.read_text(encoding="utf-8"))
        local = report["local_existing_work"]
        self.assertEqual(
            local["rank29_exact_certificate"]["status"],
            "exact_unconditional_rank_at_least_29",
        )
        self.assertEqual(
            local["rank29_exact_certificate"][
                "certified_algebraic_rank_lower_bound"
            ],
            29,
        )
        self.assertFalse(
            local["rank18_source_audit"]["explicit_rank18_family_reproducible"]
        )
        for key in ("rank29_exact_certificate", "rank18_source_audit"):
            item = local[key]
            self.assertEqual(
                item["sha256"],
                sha256((REPOSITORY / item["path"]).read_bytes()).hexdigest(),
            )
        self.assertTrue(local["existing_rank30_bounded_searches"])
        self.assertTrue(
            all(
                item["status"] == "bounded_search_no_certified_30th_point"
                for item in local["existing_rank30_bounded_searches"]
            )
        )
        for item in local["existing_rank30_bounded_searches"]:
            self.assertEqual(
                item["sha256"],
                sha256((REPOSITORY / item["path"]).read_bytes()).hexdigest(),
            )

    def test_atom_parser_and_hash_guard_are_independently_exercised(self) -> None:
        fixture = b"""<?xml version='1.0' encoding='UTF-8'?>
<feed xmlns='http://www.w3.org/2005/Atom'
      xmlns:opensearch='http://a9.com/-/spec/opensearch/1.1/'>
  <title>fixture</title><updated>2026-08-14T00:00:00Z</updated>
  <opensearch:totalResults>1</opensearch:totalResults>
  <entry><id>urn:one</id><title> A   title </title>
    <published>2026-01-01T00:00:00Z</published>
    <updated>2026-01-02T00:00:00Z</updated></entry>
</feed>"""
        parsed = parse_arxiv_atom(fixture)
        self.assertEqual(parsed["total_results"], 1)
        self.assertEqual(parsed["entries"][0]["title"], "A title")
        for spelling in (
            "Z^30 in E(Q)",
            "Z^{30} in E(Q)",
            "rank 30",
            "rank at least 30",
            "rank is at least 30",
            "rank >= 30",
            "rank ≥ 30",
        ):
            self.assertIsNotNone(RANK_30_PATTERN.search(spelling), spelling)
        with self.assertRaisesRegex(ValueError, "changed"):
            verify_payload("arxiv_rank30", b"not the pinned response")


if __name__ == "__main__":
    unittest.main()
