# Normalized planar shape exclusions through degree 150

The table records the exact scope of nine David Niedbala Giraudin Hexagon
preprints. Each linked PDF and source archive is the deposited version named in
the table. The archive contains the named ancillary script and its author's
expected output. A local script replay checks those exact calculations; it does
not audit the Guccione–Guccione–Valqui normalization, prove completeness of the
census, or independently review the paper proofs. In particular these entries
do not establish JC2 through degree 150. The already retained
[F2](sources/niedbala-giraudin-f2-2611.00685v1.pdf) and
[other degree-147](sources/niedbala-giraudin-dz147-2611.00897v1.pdf)
exclusions have different standing data.

Here a starting triple lists its corner, middle point and direction. The
successor condition in each Standing Data 1.1 is part of the normalization;
the shorthand in this table does not remove it. The census labels follow the
papers' readings of GGHV Section 6 and still require a complete independent
row-by-row comparison with the actual census.

| Census row or chain | Normalization and root pattern | Applicable result | Exact certificate |
|---|---|---|---|
| (126,84), `(12,30) → (16/3,10) → (11/6,3)`; also the stated degree-144 chain through `(16/3,10)` | Standard (3,2), intermediate regular corner in `L(3)`; no additional edge-root pattern in the theorem | [2611.00709v1](sources/niedbala-giraudin-2611.00709v1.pdf), Theorem 1.1 and Corollary 1.2 | [source archive](sources/niedbala-giraudin-2611.00709v1.zip), `anc/verify_126_84.py` |
| (84,126), `(7,35) → (19/7,5)` | Standard (2,3), starting triple `((7,35),(2,0),(7,-1))`; `h(t)=(t-c)^5`, `c≠0` | [2611.00710v1](sources/niedbala-giraudin-2611.00710v1.pdf), Standing Data 1.1 and Theorem 1.2 | [source archive](sources/niedbala-giraudin-2611.00710v1.zip), `anc/verify_84_126.py` |
| Same (84,126) starting data, root partition `4+1` | Regular corner `(18/7,4)` in `L(7)` is excluded; with the arithmetic for multiplicity three, Corollary 1.2 reduces to the preceding single-root row | [2611.00714v1](sources/niedbala-giraudin-2611.00714v1.pdf), Theorem 1.1 and Corollary 1.2 | [source archive](sources/niedbala-giraudin-2611.00714v1.zip), `anc/verify_7_35.py` |
| (88,132), `(11,33) → (19/4,8)` | Standard (2,3), starting triple `((11,33),(3,1),(4,-1))`; degree-eight `h`, whose single root of multiplicity eight is derived in Lemma 2.1 | [2611.00724v2](sources/niedbala-giraudin-2611.00724v2.pdf), Standing Data 1.1, Lemma 2.1, Theorem 1.2 | [source archive](sources/niedbala-giraudin-2611.00724v2.zip), `anc/verify_88_132.py` |
| (135,90) and (90,135), `(9,36) → (17/9,4)` | Standard (3,2) or (2,3), starting triple `((9,36),(1,0),(9,-2))`; `h(t)=(t-c)^4`, `c≠0` | [2611.00733v2](sources/niedbala-giraudin-2611.00733v2.pdf), Standing Data 1.1 and Theorem 1.2 | [source archive](sources/niedbala-giraudin-2611.00733v2.zip), `anc/verify_135.py` |
| (90,135), `(9,36) → (9,24) → (11/3,8)` | Minimal standard (2,3) complete chain; vertical-edge quartic `p(y)` has at most three distinct roots; multiplicities 2, 3 and 4 are treated separately | [2611.00944v1](sources/niedbala-giraudin-2611.00944v1.pdf), Standing Data 1.1 and Theorem 1.2 | [source archive](sources/niedbala-giraudin-2611.00944v1.zip), `anc/verify_936.py`; `anc/ggv4_chains_check.py` is an additional census cross-check |
| (56,140), `(7,21) → (13/7,3)` (F11, `j=0`) | Standard (2,5), starting triple `((7,21),(1,0),(7,-2))`; Lemma 2.1 derives `h(t)=(t-c)^3` | [2611.00717v1](sources/niedbala-giraudin-2611.00717v1.pdf), Standing Data 1.1, Lemma 2.1, Theorem 1.2 | [source archive](sources/niedbala-giraudin-2611.00717v1.zip), `anc/verify_56_140.py` |
| (42,147)/(147,42), F7, `(6,15) → (7/3,4)`; (63,147)/(147,63), F8, `(6,15) → (8/3,5)` | Standard `(m,n)∈{(2,7),(7,2),(3,7),(7,3)}`, starting triple `((6,15),(1,0),(3,-1))`; F7 root partition `4+1`, F8 single root of multiplicity five | [2611.00896v1](sources/niedbala-giraudin-2611.00896v1.pdf), Standing Data 1.1 and Theorem 1.2 | [source archive](sources/niedbala-giraudin-2611.00896v1.zip), `anc/verify_F7F8.py` |
| (150,100)/(100,150), `(10,40) → (18/5,8) → (8/5,3)` | Standard (3,2) or exchanged (2,3); intermediate regular corner in `L(5)`; no additional edge-root partition in Theorem 1.1 | [2611.00713v1](sources/niedbala-giraudin-2611.00713v1.pdf), Theorem 1.1 and Corollary 1.2 | [source archive](sources/niedbala-giraudin-2611.00713v1.zip), `anc/verify_150.py` |

## Replay and remaining coverage

From the repository root set `repo_python="$(pwd)/.venv/bin/python"`, extract
each exact source archive, then run `"$repo_python" verify_*.py` from that
archive's `anc` directory, using the precise script name in the table.
The available
environment is Python 3.12, SymPy 1.14, python-flint 0.9 and Singular 4.3.2.
The [replay receipt](SHAPE_EXCLUSION_REPLAY.json) gives each exact archive hash,
command, exit result and last status line. The scripts' reported passes are
internal arithmetic or ideal checks. They do not verify the reductions from
an arbitrary Keller pair to each standing datum. The degree-135 complete-chain
argument relies on completeness of GGHV's enumerated list at its `k=4` branch;
the ancillary chain script is a cross-check and does not replace that premise.

The next census task is to extract every GGHV Section 6 row with maximum
degree at most 150 from the original source, record its `(m,n)`, complete
corner chain and allowed root partitions, and match it to precisely one
applicable theorem or an open row. Check the normalization and successor
condition for each match, including polynomial versus Laurent source
admission. Only an exhaustive table with those bridges can support a new
global degree bound.

The F7/F8 paper, Lemmas 2.2–4.1, recognizes its polynomial source subring
using constraints from both Newton sides before layered elimination. For the
[exceptional polar fold](../../hessian/hc4-triangular-graphs/EXCEPTIONAL_FOLD_POLYNOMIAL_ADMISSION_AND_MARKED_TRACE.md),
the concrete transfer test is to find two independent Newton-side conditions
that characterize simultaneous membership of both proposed target
coordinates in the *original* polynomial source. The F7/F8 lattice argument
is a method comparison, not a theorem about fold realization.
