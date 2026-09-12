# Research

Start with the [generated status](STATUS.md), then the canonical source linked
from the claim you need. [MATH_STATUS.json](MATH_STATUS.json) is the sole
mathematical-status authority. The [catalogue](index/README.md) indexes every
registered result, including partial, falsified and replaced claims.
The [work ledger](knowledge/WORK_LEDGER.md) gives remaining gates and suggested
next steps; [structured records](index/resources.md) expose process lessons,
support conditions and curve evidence outside the prose notes.

| Programme | Entry point |
|---|---|
| Keller constructions, cancellation and arithmetic | [Verified core](verified/README.md) · [complete catalogue](index/core.md) |
| Gaussian moments, GVC, SIC and extended geometry | [Programme map](extended-geometry/README.md) |
| HC4, Hessian and Schur reductions | [Current claims and canonical sources](index/hessian.md) |
| Plane Jacobian programme | [Programme map](plane-jc/README.md) |
| Elliptic curves and rank jumps | [Programme map](elliptic-curves/README.md) · [current inventory](elliptic-curves/INVENTORY.md) |
| K3 constructions and lattice methods | [Programme map](elkies-k3/README.md) |
| Formal verification and papers | [Formal projects](formal/README.md) · [papers](papers/README.md) |

## Before another calculation

1. Search the [method memory](KNOWLEDGE_BASE.md) and current claims.
2. Read the full scope, dependencies, replacements and canonical proof.
3. Check retained inputs, completed exposure and the exact failure reason.
4. Choose the smallest replay or new calculation that answers the remaining question.

From the repository root:

```sh
python3 research/scripts/research.py search "class group"
python3 research/scripts/research.py routes --area elliptic-curves
python3 research/scripts/research.py show METHOD-EC-CACHED-CONTINUATION
python3 research/scripts/research.py work --area elliptic-curves
make check-navigation
```

Search includes archived notes and labels their authority. It reads live files;
it never runs a checker, point search, factorization or descent.

## Reproduction and history

The [short replay guide](REPRODUCE.md) explains environments, narrow checks and
the full command catalogue. The [timeline](RESEARCH_TIMELINE.md) records major
discoveries and corrections; the [archive](archive/README.md) preserves dated
investigations and superseded routes. [UNIFYING_THESIS.md](UNIFYING_THESIS.md)
gives the conceptual marked-root and boundary-construction framework.

<!-- status-consumer: EC-K3-ELKIES-2026-R17 9208e67f51fc8c97 -->
<!-- status-consumer: OP-EC-NEXT 50b9aeeb557b4df9 -->
<!-- status-consumer: EC-K3-R17-ALTERNATE-Q80-ARITHMETIC-RANK17 a304934727bb3f87 -->
<!-- status-consumer: EC-K3-R17-NONCYCLIC-4A1-DIRECT-EQUATION f657620e07f8f3f0 -->
<!-- status-consumer: EC-K3-R17-NORM12-11952-DIRECT-Q80-EQUATION 077c6409d76cbe63 -->
<!-- status-consumer: EC-K3-R17-NORM12-11952-CONTROL-J-PREIMAGES 1ef38474a0d7f629 -->
<!-- status-consumer: EC-K3-NS0024-QQ-MARKING-OBSTRUCTION b7f0cf002c0411fe -->
<!-- status-consumer: EC-K3-NS0031-MARKED-FORMAL-BRANCH b31e99bce4edac0a -->
