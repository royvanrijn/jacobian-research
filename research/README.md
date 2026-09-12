# Elliptic-curve research

The active programme is elliptic curves, including the K3 constructions used to
produce high-rank families. Start with the [current status](STATUS.md) and the
canonical source for the claim you need. [MATH_STATUS.json](MATH_STATUS.json)
remains the sole mathematical-status authority.

| Need | Entry point |
|---|---|
| Elliptic curves and rank jumps | [Programme map](elliptic-curves/README.md) · [Curve inventory](elliptic-curves/INVENTORY.md) |
| Supporting K3 constructions | [K3 map](elkies-k3/README.md) · [Process atlas](elkies-k3/ELKIES_K3_PROCESS_ATLAS.md) |
| Prior methods and failed approaches | [Research memory](KNOWLEDGE_BASE.md) |
| Remaining obligations and cleanup | [Work ledger](knowledge/WORK_LEDGER.md) · [Source reviews](knowledge/PARTIAL_REVIEW.md) |
| What this cleanup completed | [Dated report and remaining audit work](CLEANUP_REPORT_2026-09-12.md) |
| Exact results and provenance | [Claim catalogue](index/README.md) · [Structured records](index/resources.md) |
| Reproduction and discoveries | [Replay guide](REPRODUCE.md) · [Timeline](RESEARCH_TIMELINE.md) |

Before computing, read the full scope, replacements, retained inputs and failure
reason. Reuse the completed calculation when it answers the same question.
Missing local artifacts do not authorize reconstruction.

```sh
python3 research/scripts/research.py search "class group"
python3 research/scripts/research.py routes --area elliptic-curves
python3 research/scripts/research.py show METHOD-EC-CACHED-CONTINUATION
python3 research/scripts/research.py work --area elkies-k3
make check-navigation
```

Search and work lists default to this programme. `--history` includes the
[archived projects](archive/non-elliptic/README.md); `show ID` retrieves any
retained claim. Archiving preserves mathematical states and unfinished
obligations. It does not mark them proved or completed.

<!-- status-consumer: EC-K3-ELKIES-2026-R17 9208e67f51fc8c97 -->
<!-- status-consumer: OP-EC-NEXT b86e37cc3775f627 -->
<!-- status-consumer: EC-K3-R17-ALTERNATE-Q80-ARITHMETIC-RANK17 a304934727bb3f87 -->
<!-- status-consumer: EC-K3-R17-NONCYCLIC-4A1-DIRECT-EQUATION f657620e07f8f3f0 -->
<!-- status-consumer: EC-K3-R17-NORM12-11952-DIRECT-Q80-EQUATION 077c6409d76cbe63 -->
<!-- status-consumer: EC-K3-R17-NORM12-11952-CONTROL-J-PREIMAGES 1ef38474a0d7f629 -->
<!-- status-consumer: EC-K3-NS0024-QQ-MARKING-OBSTRUCTION b7f0cf002c0411fe -->
<!-- status-consumer: EC-K3-NS0031-MARKED-FORMAL-BRANCH b31e99bce4edac0a -->
