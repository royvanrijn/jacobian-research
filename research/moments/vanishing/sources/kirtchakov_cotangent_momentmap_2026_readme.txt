# Certify: replayable certificates for machine-checked mathematics

**Archival records (Zenodo):** all versions [doi:10.5281/zenodo.21799111](https://doi.org/10.5281/zenodo.21799111) · Part A (v0.1.x) [doi:10.5281/zenodo.21799112](https://doi.org/10.5281/zenodo.21799112), erratum (v0.1.2) [doi:10.5281/zenodo.21831894](https://doi.org/10.5281/zenodo.21831894) · Part B (v0.2.0) [doi:10.5281/zenodo.21799780](https://doi.org/10.5281/zenodo.21799780), update (v0.2.1) [doi:10.5281/zenodo.21831995](https://doi.org/10.5281/zenodo.21831995) · Part C (v0.3.0) [doi:10.5281/zenodo.21816010](https://doi.org/10.5281/zenodo.21816010) · Part D (v0.4.0) [doi:10.5281/zenodo.21816018](https://doi.org/10.5281/zenodo.21816018) · Part E (v0.5.0) [doi:10.5281/zenodo.21831896](https://doi.org/10.5281/zenodo.21831896) · Part F (v0.6.0) [doi:10.5281/zenodo.21832028](https://doi.org/10.5281/zenodo.21832028) · Part G (v0.7.0) [doi:10.5281/zenodo.21890619](https://doi.org/10.5281/zenodo.21890619) · Part H (v0.8.0) [doi:10.5281/zenodo.21895825](https://doi.org/10.5281/zenodo.21895825) · Part I (v0.9.0) [doi:10.5281/zenodo.21897011](https://doi.org/10.5281/zenodo.21897011) · Part J (v0.10.0) [doi:10.5281/zenodo.21898266](https://doi.org/10.5281/zenodo.21898266) · Part K (v0.11.0) [doi:10.5281/zenodo.21898343](https://doi.org/10.5281/zenodo.21898343) · Part L (v0.12.0) [doi:10.5281/zenodo.21898722](https://doi.org/10.5281/zenodo.21898722) · [erratum v0.12.1: doi:10.5281/zenodo.21910159](https://doi.org/10.5281/zenodo.21910159) · Part M (v0.13.0) [doi:10.5281/zenodo.21898996](https://doi.org/10.5281/zenodo.21898996) · Part N (v0.16.0) [doi:10.5281/zenodo.23792064](https://doi.org/10.5281/zenodo.23792064) · Part L correction and addendum (v0.17.0) [doi:10.5281/zenodo.24123940](https://doi.org/10.5281/zenodo.24123940)

Author: Daniel Kirtchakov (Independent researcher, Half Ounce Research) — daniel@halfounce.io. Repository: https://github.com/05oz/certify. Date of this snapshot: 2026-10-29.

> **Part A erratum — v0.1.2 (2026-08-06).** The Theorem D fiber count of the note over the pullback of `{Δ₂ = 0}` was wrong (stated as 2). The correct set-theoretic fiber sizes are **{3, 1, 0}**; over `{Δ₂ = 0}` the fiber has **three distinct** points (an *apparent* branch locus). No other claim depended on it. v0.1.2 supersedes v0.1.0/v0.1.1 on this point only. Details and independent re-verification: [`paper/ERRATUM-v0.1.2.md`](paper/ERRATUM-v0.1.2.md), scripts in [`scripts/erratum-check/`](scripts/erratum-check/).

Two independent bodies of work live here, sharing a method rather than a subject: **produce an artifact a skeptic can re-check without trusting the tool that made it, then say exactly what still has to be believed.**

| | Part A — Alpöge Keller map (v0.1.2) | Part B — quantum code distances (v0.2.1) |
|---|---|---|
| Subject | Degree minimality in the equivariant class of the Alpöge Keller map; moment-map structure of its cotangent lift | Certified minimum distances of eleven stabilizer codes, including the exact **d = 18** of IBM's [[288,12,18]] and a first lower bound for [[360,12,≤24]] |
| Artifact | 8 msolve Gröbner unit-ideal certificates + SymPy verification scripts | witness pairs + LRAT unsatisfiability proofs + ZX-duality permutations |
| Checker | `scripts/` (SymPy) | `qec-scripts/` — four checkers, 1,144 lines of Python, standard library only |
| Paper | `paper/preprint-dixmier-poisson.*` | `paper/preprint-qec-distances.*` |
| Certificates | `certificates/` | `qec-certificates/` |
| Provenance | [PROVENANCE.md](PROVENANCE.md) §§1–4 | [PROVENANCE.md](PROVENANCE.md) §5, [SWEEP-RECORD-QEC-2026-08-04.md](SWEEP-RECORD-QEC-2026-08-04.md) |

Every mathematical claim in either preprint maps to a script whose `assert` statements pass, or to a stored certificate a checker accepts. Nothing is conjectural unless labeled so.

**Since v0.3.0 two further parts live here under the same method:** **Part C** (v0.3.0) — the first certified determination of the tournament packing numbers ν₃(9) = 9 and ν₃(10) = 12, extending the verified range of Yuster's 2004 formula from n ≤ 8 to n ≤ 10 (`tt3-paper/`, `tt3-certificates/`, `tt3-scripts/`); **Part D** (v0.4.0) — a certificate-backed automorphism exclusion for the [[14,3,5]] quantum code existence question, which remains open (the [[14,3]] entry records lower bound 4 and upper bound 5; codetables.de, retrieved 2026-08-05): any such code has monomial automorphism group of order 2^a·3^b·5^c (`qec1435-paper/`, `qec1435-certificates/`, `qec1435-scripts/`). Both notes passed a three-lens adversarial review (claims-vs-artifacts, priority against primary sources, and a replay audit with negative controls) before release; the decision logs ship as `FIXLOG.md` in each paper directory.

**Parts E and F (v0.5.0, v0.6.0, released 2026-08-06), same method and same review pipeline:** **Part E** — an explicit 5×25 circular Florentine rectangle establishing **F_c(25) ≥ 5**, one more than the lower bound of 4 recorded in Table 62.27 (p. 677) of the *Handbook of Combinatorial Designs*, 2nd ed. (2006); verified exhaustively and exactly by a self-contained standard-library checker; no priority is claimed over H.-Y. Song's 2000 paper, which we could not access (`cfr-paper/`, `cfr-certificates/`, [SWEEP-RECORD-CFR-2026-08-06.md](SWEEP-RECORD-CFR-2026-08-06.md)). **Part F** — a bond-dimension-2 matrix-product state with integer transfer matrices that is an exact zero-energy eigenstate, at every length L, of the periodic spin-½ chain H = −Σᵢ(I+Xᵢ)(Xᵢ₊₁+Zᵢ₊₁), proved by a sixteen-equation integer telescoping certificate (technique due to Derrida–Evans–Hakim–Pasquier 1993, Gehrmann–Essler, and Garre Rubio et al., credited at point of use); no claim is made about the remainder of the spectrum (`mps-paper/`, `mps-certificates/`, [SWEEP-RECORD-MPS-2026-08-06.md](SWEEP-RECORD-MPS-2026-08-06.md)).

---

# Part A — the Alpöge Keller map

**Degree minimality in the equivariant class of the Alpöge Keller map, and the moment-map structure of its cotangent lift — certificates, verification scripts, and preprint.** All computations carried out with Claude (Fable 5), SymPy 1.14, msolve 0.10.1. Everything here is in characteristic zero.

## Credit — read this first

- **The map.** The degree-7 Keller counterexample F to the Jacobian conjecture in dimension 3 is **Levent Alpöge's**, announced **2026-07-19**. The problem was posed by **Akhil Mathew**; per the announcement, the search that produced the map was run with Claude Fable 5. Everything in this repository is downstream of that example. See T. Tao's expository post (terrytao.wordpress.com, 2026-07-21) and the Secret Blogging Seminar thread (sbseminar.wordpress.com, 2026-07-20).
- **The Dixmier corollary AND the symplectic cotangent lift.** Both are **W. G. P. Mayner's**, from the same repository and the same day. The statement that DC_n is false for n ≥ 3, with the explicit Weyl-algebra endomorphism Ψ_F as witness, was first stated publicly in a Secret Blogging Seminar comment (2026-07-20, 22:08), backed by a GitHub note prepared with Claude Fable 5 (github.com/wmayner/dixmier-counterexample, two commits on 2026-07-21). The accompanying `REPORT.md` §8 in that same repository contains the cotangent lift Φ(q,p) = (F(q), G(q)p), machine-verified there to preserve the standard symplectic form *exactly*, to have det DΦ = 1, and to be non-injective; his priority section §12 lists "the symplectic lift to C⁶" explicitly. **We credit both fully and claim neither.** An earlier version of this README credited Mayner for the Weyl endomorphism only — that was an attribution error, now corrected.
- **Most of the structural material here was anticipated.** The master equation, the anchor lemma, the parked square, the S₃/discriminant computation, the trace identity, and the exact image theorem were all obtained independently here and then found to have earlier public sources: T. Shaska, *Graded Keller maps and the Jacobian Conjecture*, arXiv:2607.20210 (v1 2026-07-22, v2 2026-07-25); the anonymous note at ulam.ai/research/jacobian.pdf (2026-07-20); MathOverflow 513387 (2026-07-20); A. Lou (2026-07-20); and Mayner. [PROVENANCE.md](PROVENANCE.md) §2b gives the claim-by-claim table with dates and theorem numbers. We retain those results for self-containedness and credit priority to those sources.
- **What is actually new here** (as of 2026-08-04): (1) the **degree-minimality theorem** in the weight-(1,−1,−2) equivariant class — no Keller lift of degree ≤ 6 exists at all in the sector containing Alpöge's map, and every degree-≤6 lift elsewhere in the class is an automorphism, by eight msolve unit-ideal certificates over Q reproduced mod 32003; (2) the **moment-map identity** ν∘Φ = μ and the **no-go lemma** (s-free C with s-affine A, B ⟹ automorphism, all weights k, all degrees); and (3) a **reframing**, not a new object: Mayner's Φ identified as an explicit PC_n witness via Adjamagbo–van den Essen, with degrees, momenta, and the quantization identity gr Ψ_F = Φ*. The k = 1, 3 certificates are in progress and are **not** claimed.
- **Open, and not claimed settled here:** JC_2, DC_1 (Zheglov's claimed proof is under review), DC_2, and unconditional minimality of degree 7 among all counterexamples in C^3. Note also that PC_n false for n > 2 is *not* our result — it follows formally from Mayner's ¬DC_3, and is already asserted in secondary sources.

## Errata

- **v0.1.2 (2026-08-06) — Theorem D fiber count over `{Δ₂ = 0}`.** Releases v0.1.0 and v0.1.1 stated that the fiber of `F` drops to **two** points over the pullback of `{Δ₂ = 0}` ("two sheets merge"). This is wrong: that fiber has **three distinct** points. `{Δ₂ = 0}` is an *apparent* branch locus — the invariant `u = 1+xy` fails to separate two of the three unramified sheets, which is exactly why `Δ₂` appears **squared** in `disc = −4·Δ₁·Δ₂²` (claims-table row 6, unaffected). The achievable set-theoretic fiber sizes are exactly **{3, 1, 0}**; `2` never occurs. Independently re-verified in exact arithmetic from the raw map: [`scripts/erratum-check/`](scripts/erratum-check/) (`fibre_check.py`, `exhibit2.py`, `structural.py`). No other claim in the note used the erroneous value — generic degree 3, the image `im F = ℂ³ ∖ Γ` (row 8), and the `S₃` monodromy (row 6) are all unaffected. Full write-up: [`paper/ERRATUM-v0.1.2.md`](paper/ERRATUM-v0.1.2.md). **This supersedes v0.1.0/v0.1.1 on this point only.**

## The claims table

Each row: claim → certificate/artifact → how to re-run. All scripts are pure SymPy and terminate with all asserts passing; the msolve certificates are re-run from their stored input files. **The "Origin" column records who first stated the claim publicly** — "ours" means we have found no earlier public source; everything else is credited. Verification status and priority are independent: a row marked *Mayner* or *Shaska* is still machine-checked here.

| # | Claim | Origin | Certificate / artifact | Re-run |
|---|-------|--------|------------------------|--------|
| 1 | F is Keller: det JF = −2 identically; matches the published form; degrees (7,6,4) | Alpöge (2026-07-19) | `scripts/core_verify.py` (PASS lines) | `python scripts/core_verify.py` |
| 2 | Rational triple collision F(1,−3/2,13/2) = F(−1,3/2,13/2) = F(0,0,−1/4) = (−1/4,0,0) | Alpöge (2026-07-19) | `scripts/core_verify.py` | same |
| 3 | Master equation (all weights k): det JG = C^k · [C·J(B,A) + k·A·J(B,C) + B·J(C,A)]; example bracket = 2 | **Shaska**, arXiv:2607.20210v2 Thm. 8.3 (2026-07-25); divisorial form v1 Thm. 6.1 (07-22) | `scripts/core_verify.py` | same |
| 4 | Critical-line restriction, Wronskian = 2; no-go lemma layer identities | no-go lemma **ours**; critical-line contraction in Shaska Prop. 5.1(1) | `scripts/core_verify.py` | same |
| 5 | Minimal cubic of the cover, trace identity (no u² term: fiber u-values sum to 0) | phenomenon in **Mayner** §4.2 (coord. x) and **Shaska** Rem. 5.4 (coord. s); u-coordinate version ours | `scripts/cover_verify.py` | `python scripts/cover_verify.py` |
| 6 | disc = −4·Δ₁·Δ₂², Δ₁ irreducible ⇒ monodromy full S₃ (cover not Galois) | **MathOverflow 513387** and **Lou** (both 2026-07-20); also Mayner §4.4, Shaska Thm. 4.4(2) | `scripts/cover_verify.py` | same |
| 7 | s ∈ Q(u,p,q) ⇒ F is generically exactly 3:1 | **ulam.ai** Thm. 4.2 / **MO 513387** (2026-07-20) | `scripts/cover_verify.py` | same |
| 8 | Exact image: F misses precisely the punctured curve (4/(27t²), 4/(3t), t), t ≠ 0 | **ulam.ai** Thm. 4.2 (2026-07-20); also Mayner §4.3, Shaska Prop. 5.1 | `scripts/cover_verify.py` (+ preprint §5) | same |
| 9 | Ψ_F is a well-defined unital endomorphism of W₃: (A) [D_j,F_i] = δ_ij, (B) [D_i,D_j] = 0 (all 9 + 9 identities) | **Mayner** (2026-07-20/21) | `scripts/weyl_verify.py`; expanded operators in `scripts/weyl_endomorphism.txt` | `python scripts/weyl_verify.py` |
| 10 | Cotangent lift Φ(q,p) = (F(q), (JF)^{-T} p) is polynomial, MᵀΩM = Ω exactly, det M = 1 | **Mayner**, `REPORT.md` §8 (2026-07-21) | `scripts/dixmier_symplectic_verify.py` | `python scripts/dixmier_symplectic_verify.py` |
| 10b | Component degrees of Φ are (7,6,4,9,10,12) | ours | `scripts/dixmier_symplectic_verify.py` | same |
| 11 | Φ has a rational triple collision in C⁶ (generically 3:1); Φ* is therefore an injective non-surjective Poisson endomorphism — a PC_n witness | non-injectivity of Φ is **Mayner's**; PC_n false for n > 2 follows formally from Mayner's ¬DC_3 via Adjamagbo–van den Essen and is already asserted in secondary sources; the **explicit identification, the momenta, and gr Ψ_F = Φ*** are ours | `scripts/dixmier_symplectic_verify.py` | same |
| 12 | Moment-map preservation: ν∘Φ = μ with μ = xp₁ − yp₂ − 2zp₃ (Hamiltonian C*-equivariance) | **ours** — no occurrence of "moment map" found in the 2026-07/08 literature on this example | `scripts/dixmier_symplectic_verify.py` | same |
| 13 | Degree minimality in the equivariant class: no Keller lift of degree ≤ 6 in the s-in-C sector; eight unit-ideal certificates over Q | **ours** — cf. Shaska Thm. 10.10 (weaker: A, B degree one in the invariants; general case explicitly open) and Mayner §7 (different ansatz, covering degree) | `certificates/ms_*_c0.ms` → `certificates/out_*_q.txt` (each basis = `[1]`) | `msolve -g 2 -f certificates/ms_I-A2c-g1_c0.ms -o out.txt` etc. (8 files) |
| 14 | Same eight certificates reproduced mod 32003 | ours | `certificates/ms_*_c32003.ms` → `certificates/out_*_p.txt` (each basis = `[1]`) | same, on the `_c32003` inputs |
| 15 | Positive control: the degree-7 example satisfies the identical normalized degree-7 system (pipeline provably contains the counterexample) | ours | `scripts/min_verify.py` part `d7control` (exact scaled substitution, asserted); computational record in `certificates/D7CONTROL-NEGATIVE-RESULT.md` (the msolve reduced-basis run on `certificates/ms_D7control_c32003.ms` did not terminate in 600 s — expected for a nonempty variety) | `python scripts/min_verify.py d7control` |

The eight emptiness certificates are the six Branch-I systems `I-A2c-g1`, `I-A2c-g2`, `I-A2L-g1`, `I-A2L-g2`, `I-B2-g1`, `I-B2-g2` and the two full-generality Branch-II systems `II-f0`, `II-f1` (18 and 19 unknowns respectively). Naming: `*_c0` = characteristic 0 input, `*_c32003` = mod 32003 input; `out_*_q.txt` = char-0 output, `out_*_p.txt` = mod-32003 output. A stored output whose reduced Gröbner basis is `[1]` is a Nullstellensatz certificate that the corresponding system is empty.

## Quickstart

```sh
python3 -m venv venv
venv/bin/pip install sympy

venv/bin/python scripts/core_verify.py                 # map, det, collisions, master equation, no-go identities
venv/bin/python scripts/cover_verify.py                # cubic, trace, discriminant, S3, 3:1, image locus
venv/bin/python scripts/weyl_verify.py                 # Weyl endomorphism: CCR (A) and (B), writes weyl_endomorphism.txt
venv/bin/python scripts/dixmier_symplectic_verify.py   # symplectomorphism, C^6 collision, moment map
venv/bin/python scripts/min_verify.py ident            # anchor lemma, layer identities, image certificate at (7/3, 4/27)
venv/bin/python scripts/min_verify.py d7control        # degree-7 positive control (exact scaled substitution)
venv/bin/python scripts/min_verify.py kdet             # det JF = -bracket_k upstairs, k = 1,2,3
venv/bin/python scripts/min_verify.py axis             # axis-target uniqueness certificate
venv/bin/python scripts/min_verify.py I                # SymPy GB reproduction, six Branch-I leaves, mod 32003
venv/bin/python scripts/min_verify.py II               # SymPy GB reproduction, two Branch-II leaves, mod 32003 (may TIMEOUT; see below)

# msolve (https://github.com/algebraic-solving/msolve), v0.10.1 used here:
brew install msolve          # macOS; or build from source
for f in certificates/ms_I-*_c0.ms certificates/ms_II-*_c0.ms;         do msolve -g 2 -f "$f" -o "${f%.ms}.out"; done   # 8 unit ideals over Q
for f in certificates/ms_I-*_c32003.ms certificates/ms_II-*_c32003.ms; do msolve -g 2 -f "$f" -o "${f%.ms}.out"; done   # mod-32003 reproductions
```

Every Python script must end with its `PASS` lines and no assertion failures — with one documented exception: `min_verify.py II` re-derives the two largest systems in SymPy under a 150 s per-leaf alarm, and on many machines prints `[II-f0] TIMEOUT after 150s` / `[II-f1] TIMEOUT after 150s` instead of a result. A timeout is **inconclusive, not a failure**: the proof of Branch-II emptiness is the stored msolve certificates (`ms_II-*` → reduced basis `[1]`), which the loop below replays; the SymPy run is a convenience cross-check. To give SymPy longer over ℚ, `venv/bin/python scripts/min_verify.py QQ II-f0 II-f1` uses a 240 s alarm. Every `ms_I-*`/`ms_II-*` msolve run must output the reduced basis `[1]`, matching the stored `out_*_q.txt` (char 0) and `out_*_p.txt` (mod 32003). Expected runtimes: the six Branch-I systems finish in well under a second each; the two `ms_II-*_c0.ms` char-0 runs take on the order of **5 minutes each** — do not interrupt them. Do not loop `ms_D7control_c32003.ms` in: that reduced-basis run does not terminate in reasonable time (see `certificates/D7CONTROL-NEGATIVE-RESULT.md`); the asserted positive control is `scripts/min_verify.py d7control`.

---

# Part B — certified quantum code distances

**Replayable minimum-distance certificates for stabilizer codes, with no solver and no proof assistant in the trusted base: the bivariate-bicycle family, and the exact distance of [[288,12,18]].** Paper: [`paper/preprint-qec-distances.pdf`](paper/preprint-qec-distances.pdf).

This is a **verification contribution.** The exact distance values are, with one exception, already known and are credited below; what did not exist is a standalone artifact anyone can replay without a SAT solver and without a proof assistant.

> **Update (2026-08-06).** [[288,12,18]] is upgraded from the lower bound `d_X ≥ 14` to the **exact value `d = 18`**, certified end to end by a profile-normalisation (`prof`) encoding plus two on-paper lemmas — Lemma P (all X-logicals even; Okada–Kasai) and Lemma S (the encoding's completeness; a standard lex-leader symmetry break) — and the shipped duality certificate. This **confirms** IBM's uncertified ILP value; it does **not** correct it. A **first lower bound of any kind for [[360,12,≤24]]** is added (`16 ≤ d ≤ 24`: the lower end certified here, the upper end cited from Bravyi et al.). The `prof` certificates are checked by a fourth standard-library checker, [`qec-scripts/check_prof.py`](qec-scripts/check_prof.py). The encoding is characterized only relative to a lex-leader baseline over the same group; it is **not** benchmarked against automated symmetry-breaking tools (BreakID, satsuma), so no novelty-versus-tools claim is made. The one honest dependency: the weight-16 exclusion (the rung that reaches 18 rather than 16) is proved only in the `prof` encoding; an independent encoding corroborates the ladder to `d_X ≥ 12`.

## Lead with the audit

The corpus was re-checked by a separate agent instance with no access to the generating pipeline. Its report is [INDEPENDENT-VERIFICATION.md](INDEPENDENT-VERIFICATION.md).

| | |
|---|---|
| Certificate checks run | **47** |
| Passed / failed | **47 / 0** |
| LRAT bytes replayed | **5,459,315,046 (5.08 GiB)**, in **pure Python** — no compiled code in the trusted base |
| Peak resident set, worst case | **79 MB** (the 2.94 GB proof) |
| IBM BB parity-check matrices vs. arXiv:2308.07915 | **byte-identical, all five codes** |
| `manifest.json` SHA-256 entries re-hashed | **182 / 182 match** |
| Negative controls | **11 / 11 correctly rejected** |
| Brute-force cross-checks (codes small enough) | **6 / 6 agree** |

## Results

| code | n, k | certified | lower-bound proof | in this repo? |
|---|---|---|---|---|
| Steane [[7,1,3]] | 7, 1 | **d = 3** | 699 B each sector | yes |
| five-qubit [[5,1,3]] (non-CSS) | 5, 1 | **d = 3** | 1,514 B | yes |
| rotated surface d=3 / 5 / 7 | 9 / 25 / 49, 1 | **d = 3 / 5 / 7** | 518 B – 75 kB | yes |
| Golay [[23,1,7]] | 23, 1 | **d = 7** | 226 kB each | yes |
| IBM BB [[72,12,6]] | 72, 12 | **d = 6** | 1.5 / 1.9 MB | yes |
| IBM BB [[90,8,10]] | 90, 8 | **d = 10** | 128 / 151 MB | yes (gzipped) |
| IBM BB [[108,8,10]] | 108, 8 | **d = 10** | 72 / 82 MB | yes (gzipped) |
| **IBM gross [[144,12,12]]** | 144, 12 | **d = 12** | 124 MB symmetry-broken; **868 / 672 MB symmetry-FREE** | symmetry-broken yes; symmetry-free regenerable |
| **IBM BB [[288,12,18]]** | 288, 12 | **d = 18** (exact) | 358 MB (`prof` ladder: 48 MB + 310 MB) | K14 rung ships; K16-exact regenerable |
| IBM BB [[360,12,≤24]] | 360, 12 | **d ≥ 16** (`d ≤ 24` cited, not certified here) | 43 MB (`prof` K14) + 14.7 MB (K12) | yes (gzipped); replays `d_X ≥ 16` — no duality certificate ships for n = 360 |

Five proofs (79 MB–646 MB compressed) are too large for git — the four `check_lower` proofs and the 310 MB exact-weight-16 `prof` proof that supplies the top rung at n = 288. [`qec-certificates/REGENERATE.md`](qec-certificates/REGENERATE.md) gives the exact CaDiCaL command, expected byte count, and expected SHA-256 for each. Everything else ships, so a reader with nothing but CPython can replay a certified **d = 12 for the gross code**, **d_X ≥ 16 at n = 288** (`prof` K14, `check_prof.py`), and **d_X ≥ 16 at n = 360** (`prof` K14, `check_prof.py`) — the final step to the exact **d = 18** is one regenerated proof away. At n = 360 the passage from `d_X ≥ 16` to `d ≥ 16` rests on the ZX-duality lemma, verified for that code in exact F₂ arithmetic but — unlike n = 288 — not shipped as a standalone certificate.

## Credit — read this first

- **The codes and the distances are IBM's.** The bivariate-bicycle family, including the gross code, is Bravyi, Cross, Gambetta, Maslov, Rall and Yoder, *Nature* **627** (2024) 778 / arXiv:2308.07915. Their distances were computed there by the MIP method of Landahl, Anderson and Rice (arXiv:1108.5738). The gross-code value `d = 12` was confirmed exactly, at MIP gap 0, by Cruz-Benito, Cross, Kremer and Faro (IBM, arXiv:2606.02418, 1 Jun 2026). **We claim no distance value.**
- **`d_X = d_Z` for BB codes is Bravyi et al.'s lemma**, from their supplemental material. Only the explicit permutation, packaged as a ~15 ms checkable certificate, is ours.
- **Machine-checked quantum distance proofs are LEAN-QEC's** (arXiv:2605.16523). Their *paper* dispatches the gross code to `cvc5` outside the Lean kernel and calls kernel replay "the next concrete engineering target" — but **their repository has moved past their paper**: commit `c73827d` (2026-07-10) records a full [[144,12,12]] verification via `bv_decide` in about 30 minutes, including kernel replay. **We claim no priority for a machine-checked gross-code distance.** What differs, at that commit: their `BB144.lean` carries two `sorry`s (`BB144_X_ker_rank` L69, `BB144_Z_ker_rank` L72) that `BB144_dist_12` routes through; three lemmas use `native_decide`, which their own paper notes extends the trusted base with Lean's compiler; no LRAT artifact is committed for BB144; and their encoding is symmetry-broken only. Ours has no admitted lemmas, ships the artifacts, includes symmetry-**free** proofs, and needs no proof assistant. Their kernel-checked ladder should **not** be described as reaching n = 108 either: `BB108.lean` carries `sorry` at L120 and L132 with `--bv_decide` commented out.
- **[[288,12,18]], stated correctly.** Bravyi et al. assert `d = 18` **exactly**, by ILP, without shipping a checkable artifact — their Table 3 lists it with no "≤", unlike [[360,12,≤24]]. The value `d = 18` is theirs; **we claim no distance value, only its first replayable determination.** We now certify `d = 18` end to end (upper: weight-18 witness; lower: the `prof` ladder + Lemma P + duality), **confirming** their value rather than correcting it. Chen, Jafari and Lai (arXiv:2606.12445) report `d ≥ 11` solver-asserted with no proof artifact and with every configuration timing out; ours is the first machine-checkable determination of this distance on record. Honest caveat: the weight-16 exclusion is single-encoding (`prof` only), corroborated by an independent encoding to `d_X ≥ 12`.
- **The even-weight lemma and symmetry breaking are prior art.** Lemma P (odd `H_Z` column weights ⟹ all-even `Z`-kernel ⟹ `d_X` even) is **Okada–Kasai's** (arXiv:2607.14091, Sec. V-A); we use it, claim it nowhere. The `prof` symmetry break (Lemma S) is an instance of **lexicographic-leader** symmetry breaking — **Crawford–Ginsberg–Luks–Roy** (KR 1996), automated in **BreakID** (Devriendt et al., SAT 2016) and **satsuma** (Anders et al., SAT 2024). We measure `prof` only against a lex-leader over the same translation group, not against those tools (not installable here), so novelty-versus-tools stays open and unclaimed.
- **[[360,12,≤24]], stated correctly.** Bravyi et al. give only the upper bound `d ≤ 24` (Table 3, marked "≤"); no lower bound of any kind had been reported. We certify `d ≥ 16` by the same `prof` + Lemma P + duality method — the **only new distance information** in this update — and cite their `d ≤ 24`.
- **Also prior art, cited at point of use:** QDistRnd (JOSS 2022, upper bounds only, "no performance guarantee"); Stim's `search_for_undetectable_logical_errors` (documented verbatim as "THIS IS A HEURISTIC METHOD"); the Webster–Jacob–Higgott survey (arXiv:2603.22532); PBLean (arXiv:2602.08692); Heule's `drat-trim`/LRAT; Biere's CaDiCaL; Sinz's cardinality encoding; Tseitin's gate encoding.

The dated adversarial sweep behind these statements is [SWEEP-RECORD-QEC-2026-08-04.md](SWEEP-RECORD-QEC-2026-08-04.md). It broke two claims of an earlier draft; both were withdrawn.

## Quickstart (Part B)

No virtual environment, no packages, no compiled binary — the checkers import only the Python standard library.

```sh
# small, instant
python3 qec-scripts/check_witness.py  qec-certificates/steane/witness_X.json
python3 qec-scripts/check_lower.py    qec-certificates/steane/lower_X_K2.json
python3 qec-scripts/check_lower.py    qec-certificates/golay/lower_X_K6.json

# the gross code, d = 12, from artifacts in this repository
# check_lower.py descriptors name the UNCOMPRESSED .lrat, so decompress first;
# check_prof.py descriptors name the .lrat.gz, so never gunzip a *_prof_*.lrat.gz
gunzip qec-certificates/bb144/lower_*.lrat.gz
python3 qec-scripts/check_witness.py  qec-certificates/bb144/witness_X.json     # d_X <= 12
python3 qec-scripts/check_lower.py    qec-certificates/bb144/lower_X_K11_sym.json  # d_X >= 12
python3 qec-scripts/check_duality.py  qec-certificates/bb144/duality.json       # d_X = d_Z

# n = 288: the exact d = 18, via the profile-normalisation ladder
# (no gunzip step here: check_prof.py reads its LRAT gzipped)
python3 qec-scripts/check_prof.py     qec-certificates/bb288/bb288_prof_K14.json   # d_X >= 16 (with Lemma P)
python3 qec-scripts/check_witness.py  qec-certificates/bb288/witness_X.json        # d_X <= 18
python3 qec-scripts/check_duality.py  qec-certificates/bb288/duality.json          # d_X = d_Z
# the final rung d_X >= 18 (weight-16 exclusion) is the one proof too large for git;
# regenerate its 310 MB LRAT from the shipped CNF per qec-certificates/REGENERATE.md, then:
#   python3 qec-scripts/check_prof.py qec-certificates/bb288/bb288_prof_K16_exact.json

# n = 360: first lower bound of any kind; what this artifact replays is d_X >= 16
python3 qec-scripts/check_prof.py     qec-certificates/bb360/bb360_prof_K14.json   # d_X >= 16 (with Lemma P)

# integrity: re-hash every audited artifact against the manifest
python3 qec-scripts/verify_manifest.py
```

Expected runtimes: the small checks are seconds; the bb144 `lower_X_K11_sym` replay is
~10–15 s in pure Python, `bb288_prof_K14` ~62 s, and `bb360_prof_K14` ~40 s. The
`prof` certificates are checked by `check_prof.py`, which rebuilds `H_X, H_Z` from each
code's polynomial spec, regenerates the CNF clause-for-clause, and replays the LRAT — so
it trusts neither the shipped matrices nor the shipped `.cnf`.

`verify_manifest.py` reports `172 match, 0 mismatch, 10 absent` on a fresh clone: the
ten "absent" are the six proofs that ship gzipped (decompress them and they match, since
the manifest hashes the *uncompressed* bytes) and the four that are too large for git.
After `gunzip qec-certificates/*/lower_*.lrat.gz` the count is 178 match, 4 absent. That
glob is deliberately restricted to the `check_lower` proofs: the three shipped
`*_prof_*.lrat.gz` are named in compressed form by their own descriptors, are not manifest
entries, and must stay gzipped.

Note that `qec-scripts/manifest.py` is the *pipeline's* manifest generator, not a verifier:
it expects a `certificates/` directory beside itself and it **overwrites** `manifest.json`.
It is included because it is part of the pipeline, and the pipeline is explicitly not trusted.

`check_lower.py` never reads the shipped `.cnf`: it regenerates the CNF from the raw parity-check matrices and replays the LRAT against its own clause list. Corrupting a shipped `.cnf` changes nothing, and the audit confirmed that by doing it.

## Known gaps in Part B

1. `bb288/duality.json` was generated **after** the audit closed — it passes `check_duality.py` and its permutation was independently re-verified, but it is outside the 47 audited checks and absent from `manifest.json`.
2. `manifest.json` covers the audited corpus including the large proofs this repository does not carry. The `prof` certificates (the `d = 18` ladder at n = 288 and the n = 360 bound) **postdate the audit and the manifest** — they are checked directly by `check_prof.py`, which regenerates every CNF from the polynomial spec, and each was replayed independently before release, but they are not among the 47 audited checks.
3. **The `d ≥ 18` rung at n = 288 is single-encoding.** The weight-16 exclusion (the step that reaches 18 rather than 16) is proved only in the `prof` encoding; an independent, differently structured encoding corroborates the ladder to `d_X ≥ 12`, and the retained symmetry-broken Sinz ladder audited earlier reaches `d_X ≥ 14`. The completeness of `prof` is the on-paper Lemma S (machine-checked hypotheses; conclusion tested against brute force on 31 small codes and against the known weight-18 logical at n = 288, but not by a second encoding at the decisive rung). The `prof` encoding is **not** benchmarked against automated symmetry-breaking tools, so no novelty-versus-tools claim is made.
4. The shipped `check_lower.py` is 488 lines; the auditor read 419. The difference is an optional totalizer cardinality encoding that no `check_lower` certificate in this release selects. The whole audited corpus was re-run against the shipped 488-line checker before release — 20 witnesses, 5 duality certificates, 23 lower-bound certificates, **48 in all, every one accepted**, including a second pure-Python replay of the 2.94 GB proof.
5. `run_all.sh` expects a `tools-drat-trim/lrat-check` binary that is not vendored. The pure-Python path needs nothing but CPython and is what every number in the paper reports.

Three further defects the audit found are reported verbatim in §8 of the paper, including a latent soundness hole in a checker branch that no shipped certificate exercises. A paper about trusted bases that suppresses its own audit findings is not one.

---

# Part C — tournament packing numbers ν₃(9), ν₃(10)

**The first certified determination of ν₃(9) = 9 and ν₃(10) = 12**, where ν₃(n) is the minimum over all n-vertex tournaments of the maximum number of arc-disjoint transitive triples, confirming Yuster's conjectured formula ⌈n(n−1)/6 − n/3⌉ at n = 9, 10. Paper: [`tt3-paper/note.pdf`](tt3-paper/note.pdf).

- **Credit.** The quantity and the formula are R. Yuster's (2004), who verified n ≤ 8 (n ≤ 7 by direct argument, n = 8 by computer, uncertified); the matching upper-bound constructions are also **Yuster's (2004)** — the note's new content is the certified integral **lower bounds** over all isomorphism classes (191,536 at n = 9; 9,733,056 at n = 10), plus SAT optimality certificates on the minimizing tournaments. Kabiya–Yuster 2008 supplies the fractional strengthening and is credited at point of use.
- **Declared trust assumption:** `gentourng` (nauty) enumerates completely — cross-validated exhaustively for n ≤ 7 and count-validated against OEIS A000568 and exact Burnside numbers at n = 9, 10.
- **Replay:** exact commands, file inventory, and all 32 MD5 pins are in [`tt3-paper/note.md`](tt3-paper/note.md) §6 ("Artifacts"); certificates in `tt3-certificates/` (LRAT optimality proofs `min9_ge10.lrat`, `min10_ge13.lrat` replay against `tt3-scripts/` with stock Python), full sweep logs included. Expected runtimes: the n = 9 sweep verification is seconds; the full n = 10 sweep verification over all sixteen gzipped slices takes ~2–3 minutes.
- The adversarial review that preceded release re-verified the sweeps with independently written code (fresh enumeration, an independent exact solver on both minimizers, negative controls); the dated novelty sweep is [SWEEP-RECORD-TT3-2026-08-05.md](SWEEP-RECORD-TT3-2026-08-05.md).

# Part D — [[14,3,5]]: a certificate-backed automorphism exclusion

**If a [[14,3,5]] qubit stabilizer code exists, its monomial automorphism group has order 2^a·3^b·5^c** — no automorphism of order divisible by 7, 11, or 13 is possible. The existence question itself **remains open and is not claimed**: the [[14,3]] entry records lower bound 4, from a construction dated June 30, 2005, and upper bound 5 (codetables.de, retrieved 2026-08-05). Paper: [`qec1435-paper/note.pdf`](qec1435-paper/note.pdf).

- **Credit.** The open entry is recorded in M. Grassl's codetables.de; the automorphism question descends from Ball–Centelles–Huber 2020 (Research Problem 1). The **CSS case is settled by Koh et al., arXiv:2601.20927** (exhaustive CSS enumeration at n ≤ 14; their Table VI gives max CSS [[14,3]] distance 4) — that result is theirs, cited and not claimed; an earlier in-house CSS derivation is subsumed and appears only as a remark. Cross–Vandeth arXiv:2501.17447 covers general stabilizer enumeration at n ≤ 9.
- **What ships:** 43 certificate files (`qec1435-certificates/`, SHA-256-pinned in the paper), the generators and the independent checker (`qec1435-scripts/`, incl. `check1435.c` and `verify_1435.py`), and the classical code tables used (`qec1435-scripts/data/`). Every candidate of every nonzero symmetry class was distance-checked; an exact-rational Krawtchouk LP lemma closes the fixed-qubit branches (a second LP lemma belongs to the order-5 work in progress, outside the theorem). Scope, gaps, and unarchived intermediate runs are disclosed in the paper itself (§5, §7).
- **Replay:** commands in [`qec1435-paper/note.md`](qec1435-paper/note.md) §7; quick control: `python3 qec1435-scripts/verify_1435.py qec1435-scripts/data/ct_14_3_stab.txt` reproduces the d = 4 control verdict. Dated novelty sweep: [SWEEP-RECORD-1435-2026-08-05.md](SWEEP-RECORD-1435-2026-08-05.md).

---

## Part G — the oriented Ramsey value k(3,4) = 21 (v0.7.0)

**A previously unknown value in Erdős Problem #112, determined and certified.** Let k(n,m) be the
least N such that every oriented graph on N vertices contains an independent set of size n or a
transitive tournament on m vertices. The published bounds for k(3,4) were 9 ≤ k(3,4) ≤ 25
(Ihringer–Rajendraprasad–Weinert 2021). We determine **k(3,4) = r(I₃,L₄) = 21**: an explicit
20-vertex oriented graph containing neither pattern (`k34-certificates/witness_sat_3_4_20.json`,
checkable in milliseconds), together with an exhaustion of the 21-vertex case — a neighbourhood
decomposition into 346 SAT instances, every one refuted with an LRAT certificate replayed by an
independently written standard-library checker. From the same campaign, `29 ≤ k(6,3) ≤ 33` (the
lower bound new, by a vertex-transitive 28-vertex witness; upper bound IRW's). Paper:
[`k34-paper/note.pdf`](k34-paper/note.pdf).

- **Reproducible without trusting us.** `k34-scripts/gen_cnf.py` builds the base propositional
  encoding directly from the definitions (no symmetry breaking on the UNSAT path), so the encoding
  the proofs are checked against is re-derivable from the problem statement alone.
  `k34-scripts/lrat_check.py` is standard-library Python and replays any LRAT certificate against
  its CNF. The LRAT proofs themselves (~245 GB) are a regenerable cache, deleted
  after verification; [`k34-certificates/CERTLOG.txt`](k34-certificates/CERTLOG.txt) records every
  certificate's verdict, checked-step count, and SHA-256.
- **What regenerating the 346 case CNFs takes.** The case CNFs come from `make_structured.py`
  applied to the 70 block-class representatives, and the solver command lines live in queue files;
  neither the representatives nor the queue files are deposited — they are part of the private
  regeneration kit, tested there end to end (it regenerates `s21_587_1.cnf` to its recorded
  `cnf_sha256`). From this repository alone a reader replays the two witnesses and any certificate
  whose CNF and proof they hold, but cannot reproduce a recorded `cnf_sha256`.
  [`k34-certificates/REGENERATE.md`](k34-certificates/REGENERATE.md) states the split artifact by
  artifact and gives the kit procedure.
- **Replay the witness now:** `python3 k34-scripts/verify_witness.py k34-certificates/witness_sat_3_4_20.json 3 4`.
- Every component was confirmed by an adversarial referee writing fresh code throughout; the
  decision log ships as `k34-paper/FIXLOG.md`.

## Part H — certified sub-threshold logical error brackets (v0.8.0)

**Superseded by Part K (v0.11.0):** the two-sided bracket below is replaced by the exact rational value of the logical error probability for these configurations; see the Part K section.

**Exact uncorrectable-set counts and a two-sided rational bracket on the logical error
probability of the rotated surface code, re-verifiable in the standard library.** Mullan,
Weippert, and Brown (arXiv:2607.27153) name the deep sub-threshold regime as inaccessible to
direct Monte Carlo simulation and answer it with a sampler. We answer it with a *certificate*.
For the distance-3 and distance-5 rotated surface codes under one round of circuit-level
depolarizing noise and a fixed lookup-table (coset-leader) decoder, we compute in exact
arithmetic the integer counts `A_w` of uncorrectable weight-`w` fault sets up to a truncation
weight `WMAX`, and convert them into a two-sided exact-rational bracket **`L ≤ P_L ≤ U`** on the
decoder's logical error probability for Stim's independent-mechanism detector error model, with
width equal to an exactly-computed Poisson-binomial tail `T = P(W ≥ WMAX+1)`. Paper:
[`wedge-paper/note.pdf`](wedge-paper/note.pdf).

- **Both verdicts, stated honestly.** Deep sub-threshold, at `p = 10^-3`, the bracket is far
  tighter than a `10^7`-shot Monte Carlo interval — about **18,500×** at `d = 3` and **626×** at
  `d = 5` (matching it by sampling would take ~`3.9×10^12` shots). The advantage is a
  deep-sub-threshold phenomenon that degrades as the expected number of firing mechanisms grows:
  at `p = 10^-2` the bracket still wins at `d = 3` (~2.3×) but **loses** at `d = 5` (~20×), because
  the weight-truncation tail is fat there — **not** because a threshold has been crossed (the
  logical rate still falls with distance, `P_L(d=5) < P_L(d=3)` at `p = 10^-2`).
- **Reproducible without trusting us.** The public unit is *certificate JSON + standard-library
  checker only.* `wedge-certificates/check_wedge.py` (d=3) and `check_wedge_d5.py` (d=5) rebuild
  the decoder, re-enumerate the fault sets, reproduce every count `A_w`, and recompute `L`, `T`,
  `U` in exact rational arithmetic from the certificate alone, failing loudly (`CHECK FAIL`, exit
  nonzero) on any mismatch. Both import only the Python standard library; no signals, subprocesses,
  network, or wall-clock. Verified on system `python3` 3.9.6.
- **Trust root, stated plainly.** Each checker certifies `L ≤ P_L ≤ U` *given* the `(det, obs, p)`
  mechanism list embedded in the certificate; it does not (and from the public artifact cannot)
  re-verify that this list equals Stim's DEM at the stated `p` — that binding lives in the private
  generator, which is **not** in this repository. The bracket bounds the independent-mechanism DEM,
  not the physical depolarizing circuit; "circuit-level" names the DEM's origin. Stim is a
  generator only and is not trusted.
- **Replay a bracket now:** `python3 wedge-certificates/check_wedge_d5.py
  wedge-certificates/certificate_d5_r1_p1over1000.json` (≈34 s → `CHECK PASS`, prints
  `2.6135024167e-05 ≤ P_L ≤ 2.6144956889e-05`). The d=3 checks are sub-second; the optional
  tighter `WMAX = 6` certificate (`certificate_d5_r1_p1over1000_w6.json`) also verifies but is
  heavy (order ten-plus minutes) — it is the edge of what a portable pure-Python checker reaches.
- **The wall (open frontier).** Exact re-verification at weight 7 (`C(77,7) = 2,404,808,340`
  subsets plus tens of millions of big-integer terms) exceeds a pure-Python checker — the next
  engine's target. The three-lens review (claims-vs-artifacts, replay with tamper controls,
  novelty/priority) found no MUST-level defect; the decision log ships as `wedge-paper/FIXLOG.md`
  and the dated sweep as [SWEEP-RECORD-WEDGE-2026-08-11.md](SWEEP-RECORD-WEDGE-2026-08-11.md).

## Part I — Kelmans' 1984 problem, verified through 22 vertices (v0.9.0)

**The first recorded computational verification of Kelmans' 1984 problem at any order.**
Let λ(G) be the maximum number of vertex-disjoint 3-vertex paths in G; counting vertices
gives λ(G) ≤ ⌊v(G)/3⌋. Kelmans asked in 1984 (Problem 1.10 of arXiv:0910.2766) whether
equality holds for every cubic 3-connected graph; at orders divisible by 3 this is the
Akiyama–Kano P₃-factor conjecture, and a positive answer would give Reed's domination
conjecture for cubic 3-connected graphs. The problem is open, and no computational
verification of it at any order appears in the record. We supply one: for **every**
3-connected cubic graph on at most 22 vertices — all **6,339,157** of them —
**λ(G) = ⌊v(G)/3⌋**, together with the applicable strong forms of Kelmans' equivalence
theorem ((z2),(z3),(z7),(z8) at orders 6, 12, 18; (t2) at 8, 14, 20; (f1),(f2) at 4, 10, 16
and (f1) at 22). Those equivalences are constructive, so a single failure at any of these
orders would have yielded an explicit cubic 3-connected graph of order divisible by 6 with
no Λ-factor. None was found. Paper: [`kelmans-paper/note.pdf`](kelmans-paper/note.pdf).

- **Reproducible without trusting us.** The public unit is *certificates + standard-library
  checkers only*. **53,356** Λ-factor certificates ship — all 43,580 through order 22 and
  9,776 at order 24 — each carrying only a `graph6` string, the vertex triples forming the
  paths, and the avoided vertices. Two checkers written independently of each other and of
  the searchers, [`kelmans-scripts/verify_cert.py`](kelmans-scripts/verify_cert.py) and
  [`kelmans-scripts/refcert.py`](kelmans-scripts/refcert.py), re-derive everything from the
  `graph6` string: own decoder, cubicity, 3-connectivity re-proved by exhaustive vertex-pair
  deletion, each triple a path, triples plus avoided vertices a partition with
  `|avoided| = n mod 3`. Both import only the Python standard library; verified on system
  CPython 3.9.6 and 3.14.2. **Every one of the 53,356 certificates passes both, 0 rejected.**
- **Replay the whole corpus now:** `python3 kelmans-scripts/verify_cert.py
  kelmans-certificates/certs_n*.txt` (about 35 s → `VERIFIED 53356 certificates`).
- **Two pipelines, no shared code.** Different connectivity tests (bitmask BFS vs.
  union–find), different search orders, different failure caches; they agree on every count
  at every order and report zero failures. Generated counts match OEIS A002851 and
  Brinkmann–Goedgebeur–McKay; filtered counts match OEIS A204198 at every order and
  McKay–Royle at orders 10–20. Per-slice summaries for both pipelines ship in
  [`kelmans-certificates/summaries/`](kelmans-certificates/summaries/).
- **Both failure paths are exercised.** With the 3-connectivity filter disabled the sweep
  finds the unique sub-3-connected base-claim failure at orders 10–16 (`O???E?oBEAWOKGK_@o?W_`,
  λ = 4 < 5), and the strong-form paths fire at orders 10, 12, 14, 16 with a recorded
  per-type breakdown. Eight distinct classes of doctored certificate are each rejected by
  the gate that targets them: `refcert.py`, run with `--g6set` and `--check3c`, rejects all
  eight; `verify_cert.py` carries no membership gate, so it rejects seven and accepts the
  certificate written for a non-canonical relabelling of the graph. Two
  controls-on-the-controls confirm the diagnosis. See
  [`kelmans-certificates/controls/`](kelmans-certificates/controls/).
- **Stated at exactly the strength the recount supports.** Orders ≤ 20 carry the referee's
  signed verdict of 2026-08-06; order 22's independent recount completed 2026-08-11 (read
  7,319,447, kept 5,909,292, zero failures, all 5,904 certificates cross-checked with
  membership and 3-connectivity). **Order 24 is search-side complete only** — 98,101,019
  graphs, zero failures, counts matching the published enumeration — with **no independent
  recount**, and is reported at that strength and no higher, in the paper and in
  [`kelmans-certificates/verdict-n22-24.md`](kelmans-certificates/verdict-n22-24.md).
- **Not in this deposit:** the two searchers. A reader can re-check every positive answer
  they gave and reproduce the enumeration counts, but cannot re-run a sweep from this
  deposit alone; the paper says so in §4. Replay instructions, exact commands and SHA-256
  hashes: [`kelmans-certificates/REGENERATE.md`](kelmans-certificates/REGENERATE.md).
- The three-lens review (claims-vs-artifacts, replay with tamper controls,
  novelty/priority) produced five must-fix entries in the fix pass's log: two about the
  cited source (the smallest member of Kelmans' R_s family is cyclically 5-connected, not
  6-connected; the page range of the Kostochka–Stodolsky reference), two about this
  deposit's own records (the per-type breakdown of the order-12 negative control, and an
  n = 22 referee verdict that had been left an unsigned stub), and one about the tense of
  the deposit sentence. None of the five changed a theorem or a proposition; M-1 corrected the per-type
  breakdown of a negative control and added an order-16 line the draft had omitted. The decision log ships as
  [`kelmans-paper/FIXLOG.md`](kelmans-paper/FIXLOG.md) and the dated sweep as
  [SWEEP-RECORD-KELMANS-2026-08-11.md](SWEEP-RECORD-KELMANS-2026-08-11.md).

## Part J — the k(3,4) extremal graph is not unique (v0.10.0)

Part G left open (its Question 8.1) whether the 20-vertex extremal graph for k(3,4) = 21 is
unique up to isomorphism. **It is not.** At least **thirteen** pairwise non-isomorphic
{I₃,TT₄}-free oriented graphs on 20 vertices exist, and every one of them is *rigid* — trivial
automorphism group. Each is verified free of both patterns by an exhaustive test
(`k34add-scripts/verify_witnesses.py`) over all C(20,3) = 1140 triples, none independent, and all
C(20,4) = 4845 quadruples, none transitive. That tester is not code-independent of the private
structure-mining dig; the rigidity and non-isomorphism verdicts around it rest on a
Weisfeiler–Leman refinement and canonical form that are.

The note also proves that the Paley tournament QR₇ inside these graphs is **forced, not designed**:
a vertex's non-neighbourhood induces a tournament with no transitive quadruple, such tournaments
have at most 7 vertices, and on 7 vertices there is exactly one — so any vertex with seven
non-neighbours necessarily carries a QR₇ block (`k34add-scripts/verify_qr7_lemma.py`). A third
script shows that the largest {I₃,TT₄}-free tournament blow-up is QR₇[I₂] on 14 vertices, six
short of the extremal order 20, so the algebraic route yields only k(3,4) ≥ 15 against the truth
21 and the extremal family is substantially non-algebraic (`k34add-scripts/blowup_bound.py`).

Paper: [`k34add-paper/note.pdf`](k34add-paper/note.pdf). Witnesses and verifiers:
`k34add-certificates/`, `k34add-scripts/`.

## Part K — the exact logical error probability (v0.11.0)

**The Part H bracket, collapsed to a single exact rational — and the weight-7 wall broken at
d≤5.** Part H (v0.8.0, DOI [10.5281/zenodo.21895825](https://doi.org/10.5281/zenodo.21895825))
certified two-sided rational brackets `L ≤ P_L ≤ U` on the logical error probability of the
d=3 and d=5 rotated surface codes (one round, circuit-level depolarizing noise, lookup-table
coset-leader decoder) and named its own frontier: exact re-verification at weight 7 exceeds a
pure-Python checker, because it enumerates `C(77,7) = 2,404,808,340` fault sets. This part
removes both the truncation and the wall: a syndrome-space character sum (a MacWilliams-type
Walsh–Hadamard argument, Theorems 2.1–2.2 of the note) computes **every** uncorrectable count
`A_w` in one `O(2^n n)` pass — delivering the previously unreachable **`A_7 = 832,441,445`**
and the full spectrum through `A_77 = 1` — and its probability-weighted variant computes the
**exact** dyadic-rational `P_L`, eliminating the bracket. At all four published operating
points the exact value lies strictly inside the Part H bracket (positions 0.36–0.52 of the
width), inside the tighter optional WMAX=6 bracket at d=5, p=1/1000, and inside both 10^7-shot
MC confidence intervals. Paper: [`wedge2-paper/note.pdf`](wedge2-paper/note.pdf).

- **Supersedes our own Part H core result, and nothing external.** The brackets survive inside
  the new certificates in one role: as published independent enclosures the exact values are
  verified to satisfy. The DEMs, decoder, hashes, and trust architecture are Part H's,
  byte-identical.
- **Reproducible without trusting us.** The public unit is *certificates + standard-library
  checker only*: [`wedge2-certificates/check_wedge2.py`](wedge2-certificates/check_wedge2.py)
  rebuilds the full BFS decoder, verifies syndrome-space spanning, runs the Walsh–Hadamard
  transforms on `array('q')` buffers, and re-derives the full `A_w` and `N_w` spectra, the
  circuit-level distance, two closed-form invariants (`sum A_w = 2^{m-1}`,
  `sum N_w = 2^{m-n-1}`), the exact `P_L` numerator/denominator, and the Part H bracket
  containment. Imports only `hashlib, json, sys, array, fractions`; no signals, subprocesses,
  network, or wall-clock. **All four certificates: `CHECK PASS`.**
- **Replay now:** `python3 wedge2-certificates/check_wedge2.py
  wedge2-certificates/certificate_d3_r1_p1over1000_exact.json` (0.15 s → `CHECK PASS`); the
  d=5 certificates verify in ≈129 s each within 0.85 GB. Identity self-test:
  `python3 wedge2-scripts/identity_selftest.py` (both theorems vs. exact brute force on random
  small DEMs). Tamper battery: `python3 wedge2-scripts/tamper_demo_w2.py
  wedge2-certificates/check_wedge2.py wedge2-certificates/certificate_d3_r1_p1over1000_exact.json`
  (8/8 corruption classes rejected).
- **Verification battery behind the release:** at d=3, three independent exact routes
  (character sum, CRT-reconstructed syndrome convolution, Gray-code full 2^23 enumeration)
  agree digit-for-digit at both operating points; at d=5 the character-sum values were checked
  against the independent convolution route modulo 24 distinct 25-bit primes (~600 bits of
  agreement) and the spectrum against every Part H certified count plus both invariants.
  Log: [`wedge2-paper/FIXLOG.md`](wedge2-paper/FIXLOG.md).
- **The honest trade, stated in §6 of the note:** the `C(m,w)` enumeration wall is exchanged
  for a `2^n` syndrome-space wall. One-round d=7 has n≈48 detectors (`2^48` syndromes), so d=7
  remains out of reach on a laptop by this method; the obstruction (the global
  minimum-cardinality decoder blocks localisation) is named precisely, and the note's
  questions ask for the structure that would break it.
- **Trust root unchanged:** the mechanism list (hash-pinned, identical to Part H) and its
  private binding to Stim's DEM; the certified `P_L` is the DEM's, not the physical
  circuit's. Dated sweep:
  [SWEEP-RECORD-WEDGE2-2026-08-12.md](SWEEP-RECORD-WEDGE2-2026-08-12.md).

## Part L — certified Newell demagnetization-tensor reference table (v0.12.0; erratum v0.12.1; correction and addendum v0.17.0)

**Two-sided rational enclosures of the Newell demagnetization-tensor entries, and a rigorous
map of where double-precision micromagnetics loses its digits.** OOMMF and other finite-difference
micromagnetic simulators build their demagnetizing field from Newell's demagnetization tensor of a
pair of uniformly magnetized rectangular cells (Newell–Williams–Dunlop, *J. Geophys. Res.* 98
(1993) 9551); some evaluate Newell's closed form, while MuMax3 computes the tensor by numerical
integration. It is
documented but uncertified that the closed-form evaluation loses all significant digits to
catastrophic cancellation once the cells are more than a few hundred cell widths apart
(Chernyshenko–Fangohr, arXiv:1403.1978: relative error `~10^-15 r^6`, no significant digits past
`~300` cells). We pin each entry with a certified interval `[N_lo, N_hi]` of dyadic rationals,
`N_lo ≤ N_true ≤ N_hi`, obtained by evaluating Newell's formulas in outward-rounded interval
arithmetic (77-digit working precision; `sqrt`, `atan`, `log` enclosed by rigorously-truncated
series), and measure the floating-point failure against it. Paper:
[`demag-paper/note.pdf`](demag-paper/note.pdf) (revised in v0.17.0, [doi:10.5281/zenodo.24123940](https://doi.org/10.5281/zenodo.24123940)).

- **The pathology, certified.** For the canonical cube on-axis `Nxx`, the naive double-precision
  analytic value falls by ~6 correct decimal digits per decade of separation — from **15.2** digits
  at one cell to **0.4** near `n = 300` (no correct significant figure, the documented breakdown)
  to **−8.9** at `n = 10^4`, where it returns `+1.2×10^-4` for a true value of `−1.6×10^-13`
  (wrong sign, nine orders of magnitude too large). The breakdown radius is **not** a universal
  300: across common geometries it ranges from `~100` cells (elongated cells, off-diagonal) to
  `~2000` (a thin film's out-of-plane `Nzz`); the certificate maps it for each of 50
  geometry/component pairs. A truncated asymptotic expansion is the mirror image, poor at short
  range and good at long, and the crossover is bracketed rigorously for each truncation. *(Through
  v0.16.0 this entry called the three-term form "OOMMF's asymptotic expansion". It is a three-term
  truncation taken from the demagderiv notebook; OOMMF's `demagcoef.cc` carries a fourth term and,
  since 2.0a2, chooses its switching radii from assumed error and speed models. See note §4 and CORRECTIONS.md.)*
- **Addendum (v0.17.0): a proven bound on the asymptotic truncation error.** For two identical
  axis-parallel cells with diagonal `D` and centre offset `r > D`, every tensor component and every
  odd truncation order `n ≥ 3`, the series truncated at order `n` has error at most
  `V D^(n−3) / (2π r^(n−2) (r² − D²))`, with a sharper shape-dependent form. This is the form of
  OOMMF's assumed error model with the cell diagonal in place of the longest edge and `1/(2π)` in
  place of `1/5`; the constant is attained in the limit of a thin needle, and the bound is loose for compact
  cells (about 105 times the true leading error for a cube at `n = 11`). Proof: note §8. Replay:
  `python3 demag-scripts/check_remainder_bound.py demag-certificates/demag_certificate_six_cells.json`
  tests both forms on all 1,162 certified entries with `r > D`. A second certificate,
  `demag_certificate_six_cells.json` (1,310 entries), adds the thin-film cells `(1,1,1/10)` and
  `(1,1.3,1/10)` and passes the unchanged checker; its other 862 entries are identical to the first.
- **Tighter than double precision everywhere.** At the worst point (`n = 10^4`) the enclosure is
  ~54 digits tight where double precision has none. Pre-registered kill condition (DEAD if
  enclosures cannot beat double precision anywhere in the regime): **not triggered — LIVE.**
- **Reproducible without trusting us.** The public unit is *certificate JSON + standard-library
  checker only.* [`demag-certificates/check_demag.py`](demag-certificates/check_demag.py)
  recomputes every enclosure by interval arithmetic and Newell evaluation that is not
  independent of the private engine (see note §5),
  verifies containment (`N_lo ≤ N_true ≤ N_hi`) with a width-sanity bound, recomputes each naive
  double (to a relative 1e-4; bit for bit with `--strict-naive`), recomputes every rigorous
  digit-loss bracket, and re-tests the tensor's own
  identities (trace encloses 0 at all 136 mutual points; each self-term's three diagonals sum to
  enclose 1). It imports only `sys, json, math, hashlib, fractions, collections`; no signals,
  subprocesses, network, or wall-clock. **All 862 entries: `CHECK PASS`.**
- **Replay now:** `python3 demag-certificates/check_demag.py
  demag-certificates/demag_certificate.json --sample 40` (seconds → `CHECK PASS` on a
  deterministic sample); the full 862-entry verification is `CHECK PASS` in ≈7.5 min. Independent
  anchor: `python3 demag-scripts/anchor_check.py` confirms that the 16 anchor rows taken from the
  31-row OOMMF/Maple 50-digit check table in `demagcoef.cc` agree with the recomputed enclosures
  to ≥49.6 digits; `python3 demag-scripts/tamper_demo.py`
  runs six corruption controls, all rejected.
- **Trust root, stated plainly.** The certificate bounds the value of the *analytic Newell tensor
  entry* the simulators compute; whether that entry is the right physical kernel for a given
  discretization is Newell's modelling choice, cited not claimed. The generator and the
  interval/Newell engine that built the certificate are **not** in this repository; the checker
  recomputes every enclosure, naive value and digit-loss bracket from the certificate alone
  (`demag-scripts/four_term_digits.py` recomputes the three-term truncation behind the stored
  asymptotic value) but shares 212 of its 558
  executable lines with them (counted on the v0.15.0 file; 568 since v0.16.0); the external anchor is 16 rows of OOMMF's Maple table, and the consistency
  identities cover all 862 entries. Build/verification
  log: [`demag-paper/FIXLOG.md`](demag-paper/FIXLOG.md); dated sweep:
  [SWEEP-RECORD-DEMAG-2026-08-12.md](SWEEP-RECORD-DEMAG-2026-08-12.md).

## Part M — certified ZEFOZ brackets for ¹⁶⁷Er³⁺:Y₂SiO₅ (v0.13.0)

**Existence, location, and curvature of the published ZEFOZ points, certified in exact rational
arithmetic — and a pre-registered kill condition that fired on the completeness question.** The
zero-first-order-Zeeman (ZEFOZ) points of ¹⁶⁷Er³⁺:Y₂SiO₅ — the magnetic fields at which a hyperfine
transition frequency is stationary — set the predicted coherence times of a leading solid-state
quantum-memory platform. The published atlas (Matsuura et al., *Phys. Rev. B* **113**, 085421 (2026);
arXiv:2412.10126) is produced by Newton iteration from finite grids, carries no error bounds on its
locations, frequencies, or curvatures, and its authors state that the number of points found depends
on the initial grid. We certify the points; we do not certify the list. Paper:
[`zefoz-paper/note.pdf`](zefoz-paper/note.pdf).

- **The twenty published nonzero-field points, certified.** For each of the twenty entries of v3
  Table 5 (ten per crystallographic site) the certificate proves, at an exactly specified dyadic
  rational field: eigenvalue brackets of width 2e-10 MHz (hence 4e-10 MHz on the transition), a
  gradient bound `|∇f| ≤ 3.2e-37` MHz/mT, two-sided brackets of width ≤ 2.1e-15 MHz/mT² on all three
  eigenvalues of the 3×3 frequency Hessian with certified signs, and — by a Krawczyk contraction with
  a rigorous third-derivative bound — **existence and local uniqueness of an exact stationary point
  within 2.9e-14 mT** (worst contraction ratio 0.096). The certified signatures settle the
  stationary-point types: **none of the twenty is a local minimum** of its transition frequency —
  thirteen are saddles, seven are local maxima.
- **Zero field, exactly.** A signed-permutation time-reversal identity (`M conj(H₀) Mᵀ = H₀`,
  `M conj(Z_k) Mᵀ = −Z_k`) holds exactly over Q[√7,√12,√15] for both sites, so the spectrum is even
  in **B** and all **120 transitions of both sites are stationary at B = 0** — the folklore statement
  upgraded to a machine-checked identity — with certified curvature brackets for the ten published
  zero-field pairs.
- **The completeness question is reported dead, by a pre-registered rule.** Kill condition K2 was
  registered in advance: if measured statistics project more than 200 laptop-hours to close the box
  `‖B‖∞ ≤ 100` mT, kill the claim. Six stratified branch-and-bound sample chunks give a **strict**
  lower bound above 7.9e5 laptop-hours per site — more than three orders of magnitude over budget.
  **K2 fired.** The obstruction is named and quantified rather than left as a gap: quasi-degenerate
  hyperfine doublets force the per-level spectral-gap machinery below its validity radius. The
  measured statistics ship as run logs in `zefoz-scripts/kill-logs/` and are **not** part of the
  certified surface.
- **Three errata in the reference, recorded with version history.** A load-bearing sign error in the
  printed site-1 quadrupole matrix (`Q₂₃ = +15.5` MHz) stood through arXiv v1 and v2 for eleven
  months and was corrected upstream in v3 — credited, not claimed; two v3 Table 5 frequencies are
  inconsistent with their own stationary points by 2.74 and 4.34 MHz (certified restatements given);
  and one printed field vector carries an inconsistent sign pairing. Any reader who rationalized v1
  or v2 as printed inherited a wrong Hamiltonian.
- **Reproducible without trusting us.** The public unit is *certificate JSON + standard-library
  checker only.* [`zefoz-certificates/zefoz_checker2.py`](zefoz-certificates/zefoz_checker2.py)
  re-derives every claim from the certificate alone as inequalities among exact rationals. It imports
  only `sys, json, math, fractions`; no floating point, no eigensolver, no signals, subprocesses,
  network, or wall-clock. 247 of its 541 executable lines are shared with the private
  certification engine; the interval arithmetic and the inertia count are not. Containment of
  all 352 eigenvalue brackets is separately confirmed by `zefoz-scripts/anchor_check.py`, which
  shares 14 of 100 lines with that engine. **23-object certificate: exit 0.**
- **Replay now:** `python3 zefoz-certificates/zefoz_checker2.py zefoz-certificates/certificate2.json`
  (1.9 MB certificate, ≈40 s → `CERTIFICATE VERIFIED`, exit 0; also passes on CPython 3.9.6); the re-verified pilot is
  `python3 zefoz-certificates/zefoz_checker_pilot.py zefoz-certificates/certificate_pilot.json`
  (≈47 s, exit 0). Independent anchor: `python3 zefoz-scripts/anchor_check.py` confirms containment
  against an independent 60-digit diagonalization; `python3 zefoz-scripts/tamper_demo.py` runs six
  corruption controls plus an untampered control, all correctly decided.
- **Trust root, stated plainly.** The certificate is unconditional about the mathematics of the
  exact effective spin Hamiltonian and **silent about the spectroscopy** — whether that Hamiltonian
  describes the physical crystal is a separate, experimental question, and the printed matrix
  entries are adopted exactly as published. The identification of certified points with rows of
  Table 5 is by proximity of the published rounded coordinates. The search engine and the
  interval-arithmetic generator are **not** in this repository. Build/verification log:
  [`zefoz-paper/FIXLOG.md`](zefoz-paper/FIXLOG.md); dated sweep:
  [SWEEP-RECORD-ZEFOZ-2026-08-12.md](SWEEP-RECORD-ZEFOZ-2026-08-12.md).

## Part N — circuit-level distances of the Kishony–Fowler colour-code circuit, d = 3 to 15 (v0.16.0)

**d_circ = 2, 4, 6, 7, 9, 11, 12 at d = 3, 5, 7, 9, 11, 13, 15, for d rounds and for one round,
which equals d − ⌊(d+3)/6⌋ at every point.** The circuit is the single-auxiliary colour-code
syndrome-extraction circuit of Kishony and Fowler (arXiv:2603.28852), built by their own companion
code (Classiq/classiq-library at commit a43d3bd4) with their default schedule and two-qubit
depolarising noise after every CNOT, as a Z-basis memory experiment. Every value is two-sided: an
explicit fault set, and an LRAT proof that no smaller set exists. A standard-library checker
re-derives the fault mechanisms from the hash-pinned circuit text without Stim, regenerates each
formula and replays each proof; at d = 15 drat-trim's `lrat-check` also accepted the refutation.
Paper: [`colorcode-dcirc-paper/note.pdf`](colorcode-dcirc-paper/note.pdf).

- **Credit.** The circuit, its schedule and the formula are Kishony and Fowler's; they state that
  the formula holds numerically up to d = 13. The extended code through which the lower bounds at
  d ≥ 11 go is Strikis, Browne and Beverland's (arXiv:2603.05481), and the projection check is the
  one Kornstein makes for bivariate-bicycle cycles (doi:10.5281/zenodo.23394575). The values at
  rounds = d for d = 3–9 and at one round for d = 5–11 reproduce qsim-lab's
  (github.com/dylanneve1/qsim-lab). Later results on the same CNOT schedule are Crocker's
  (doi:10.5281/zenodo.23633600) and Issa and Martonosi's (arXiv:2610.33347). What was and was not
  found elsewhere, with its as-of date and the sources read, is stated in
  [`colorcode-dcirc-certificates/README.md`](colorcode-dcirc-certificates/README.md) and
  [SWEEP-RECORD-COLORCODE-DCIRC-2026-10-20.md](SWEEP-RECORD-COLORCODE-DCIRC-2026-10-20.md).
- **Replay now:** from `colorcode-dcirc-certificates/`, `python3 run_all.py --fast` (about 10–15 s)
  or `python3 run_all.py` (everything shipped plus the tamper battery, about 1.5–2.5 minutes).
  Standard-library Python 3.9 or later.
- **Not in git:** five proofs larger than 50 MB, among them the refutations at d = 13 and d = 15.
  They are attached to the GitHub release v0.16.0, the two larger than 2 GiB in parts.
  [`colorcode-dcirc-certificates/REGENERATE.md`](colorcode-dcirc-certificates/REGENERATE.md) gives
  the joining commands, the SHA-256 of each file and the route to regenerate them. Without them,
  the shipped proofs cover every value up to d = 11.
- Build and verification log: [`colorcode-dcirc-paper/FIXLOG.md`](colorcode-dcirc-paper/FIXLOG.md)
  and `colorcode-dcirc-certificates/VERIFY.log`. Dated record:
  [`colorcode-dcirc-paper/PROVENANCE.md`](colorcode-dcirc-paper/PROVENANCE.md).

## Layout

```
README.md            this file
PROVENANCE.md        timeline, what is new, how to verify, what is not claimed (§§1-4 Part A, §5 Part B)
CITATION.cff         citation metadata
.zenodo.json         Zenodo deposit metadata
LICENSE-CODE         Apache-2.0 (code and machine-readable certificates)
LICENSE-DOCS         CC-BY-4.0 (prose and paper)
paper/               both preprints (LaTeX + PDF + readable Markdown mirror)

  -- Part A: Alpoge Keller --
certificates/        17 msolve input files (ms_*.ms) + stored outputs (out_*.txt) + D7 control record
scripts/             the verification scripts (incl. min_verify.py) + expanded Weyl operators

  -- Part B: quantum code distances --
qec-certificates/    the certificate corpus, by code; manifest.json; REGENERATE.md
qec-scripts/         the four checkers (check_witness / check_lower / check_duality /
                     check_prof), verify_manifest.py, and the (untrusted) generating pipeline
INDEPENDENT-VERIFICATION.md   the audit: 47/47 checks, 5.08 GiB replayed in pure Python
SWEEP-RECORD-QEC-2026-08-04.md   the dated adversarial prior-art sweep

  -- Part C: tournament packing (v0.3.0) --
tt3-paper/           note (LaTeX + PDF + Markdown mirror) + FIXLOG.md review log
tt3-certificates/    minimizers, SAT optimality certificates (LRAT), full sweep logs, MD5SUMS.txt
tt3-scripts/         tt3pack.c searcher, CNF encoder, stdlib-only independent verifiers
SWEEP-RECORD-TT3-2026-08-05.md   dated novelty sweep

  -- Part D: [[14,3,5]] automorphism exclusion (v0.4.0) --
qec1435-paper/       note (LaTeX + PDF + Markdown mirror) + FIXLOG.md review log
qec1435-certificates/  43 SHA-256-pinned certificates by symmetry class
qec1435-scripts/     generators, check1435.c, verify_1435.py, data/ classical code tables
SWEEP-RECORD-1435-2026-08-05.md  dated novelty sweep

  -- Part E: circular Florentine rectangles, F_c(25) >= 5 (v0.5.0) --
cfr-paper/           note (LaTeX + PDF + Markdown mirror) + FIXLOG.md review log
cfr-certificates/    CFR_5_25.json (the explicit 5x25 rectangle) + stdlib checker
                     (verify_cfr525.py, exhaustive and exact)
SWEEP-RECORD-CFR-2026-08-06.md   dated novelty sweep

  -- Part F: exact all-length zero-energy MPS eigenstate (v0.6.0) --
mps-paper/           note (LaTeX + PDF + Markdown mirror) + FIXLOG.md review log
mps-certificates/    object.json (the bond-dimension-2 MPS + telescoping certificate)
                     + stdlib checker (reverify.py, exact integers/Fractions);
                     xcheck.py is an auxiliary numpy cross-check, NOT trusted base
SWEEP-RECORD-MPS-2026-08-06.md   dated novelty sweep

  -- Part G: the oriented Ramsey value k(3,4) = 21 (v0.7.0) --
k34-paper/           note (LaTeX + PDF + Markdown mirror) + FIXLOG.md review log
                     + REFEREE-k63.md
k34-certificates/    the two lower-bound witnesses + CERTLOG.txt (per-certificate
                     verdict, checked-step count, SHA-256) + REGENERATE.md;
                     the ~245 GB LRAT corpus is a regenerable cache, not shipped
k34-scripts/         gen_cnf.py encoder, make_structured.py (structured-instance
                     builder), lrat_check.py and verify_witness.py (stdlib
                     checkers), audit_cnf.py; audit_multiset.py requires the
                     private regeneration kit and does not run from this
                     repository alone — see k34-certificates/REGENERATE.md
SWEEP-RECORD-K34-2026-08-11.md   dated novelty sweep

  -- Part H: certified sub-threshold logical error brackets (v0.8.0) --
wedge-paper/         note (LaTeX + PDF + Markdown mirror) + FIXLOG.md review log
wedge-certificates/  5 certificate JSONs (d=3, d=5; WMAX 4/5, plus optional WMAX 6)
                     + two stdlib checkers (check_wedge.py, check_wedge_d5.py); NO generator
SWEEP-RECORD-WEDGE-2026-08-11.md  dated novelty + verification sweep

  -- Part I: Kelmans' 1984 problem, verified through 22 vertices (v0.9.0) --
kelmans-paper/       note (LaTeX + PDF + Markdown mirror) + FIXLOG.md review log
kelmans-certificates/  53,356 Lambda-factor certificates by order + per-slice summaries
                     + negative controls + both referee verdict records + REGENERATE.md
kelmans-scripts/     the two independent stdlib checkers (verify_cert.py, refcert.py),
                     the 3-connectivity recount driver, the control builder; NO searchers
SWEEP-RECORD-KELMANS-2026-08-11.md  dated novelty + verification sweep

  -- Part J: the k(3,4) extremal graph is not unique (v0.10.0) --
k34add-paper/        note (LaTeX + PDF + Markdown mirror)
k34add-certificates/ the 13 pairwise non-isomorphic rigid witnesses (w01_W.json,
                     w02-w13.json) + README.md inventory
k34add-scripts/      the three stdlib verifiers: verify_witnesses.py (both patterns,
                     all 1140 triples and all 4845 quadruples), verify_qr7_lemma.py,
                     blowup_bound.py; NO searchers
SWEEP-RECORD-K34ADD-2026-08-11.md  dated novelty + verification sweep

  -- Part K: the exact logical error probability (v0.11.0) --
wedge2-paper/        note (LaTeX + PDF + Markdown mirror) + FIXLOG.md build/verification log
wedge2-certificates/ 4 exact-P_L certificate JSONs (d=3, d=5; p=1/1000, 1/100; full A_w and
                     N_w spectra + exact rational P_L) + stdlib checker (check_wedge2.py);
                     NO generator, NO transform engine
wedge2-scripts/      identity_selftest.py (both theorems vs. brute force, stdlib) +
                     tamper_demo_w2.py (8 controls)
SWEEP-RECORD-WEDGE2-2026-08-12.md  dated novelty + verification sweep

  -- Part L: certified Newell demagnetization-tensor table (v0.12.0) --
demag-paper/         note (LaTeX + PDF + Markdown mirror) + FIXLOG.md build/verification log
demag-certificates/  demag_certificate.json (862 certified enclosures) +
                     demag_certificate_six_cells.json (1,310; v0.17.0) + stdlib checker
                     (check_demag.py); NO generator, NO interval/Newell engine
demag-scripts/       anchor_check.py (16-row anchor subset of the OOMMF/Maple table)
                     + tamper_demo.py (6 controls); v0.17.0: series_terms.py (the asymptotic
                     series to any order), check_remainder_bound.py (tests the note's §8
                     bound on every certified entry with r > D), four_term_digits.py (Table 1)
SWEEP-RECORD-DEMAG-2026-08-12.md  dated novelty + verification sweep
SWEEP-RECORD-DEMAG-BOUND-2026-11-10.md  dated sweep for the v0.17.0 addendum

  -- Part M: certified ZEFOZ brackets, 167Er3+:Y2SiO5 (v0.13.0) --
zefoz-paper/         note (LaTeX + PDF + Markdown mirror) + FIXLOG.md build/verification log
zefoz-certificates/  certificate2.json (time-reversal identity; 2 zero-field points; 20
                     Krawczyk existence points with gradient/Hessian/signature brackets)
                     + stdlib checker (zefoz_checker2.py); the re-verified pilot
                     certificate_pilot.json + zefoz_checker_pilot.py; NO generator
zefoz-scripts/       anchor_check.py (independent 60-digit containment, mpmath) +
                     tamper_demo.py (6 tampers + control) + kill-logs/ (the measured
                     branch-and-bound statistics that fired kill condition K2)
SWEEP-RECORD-ZEFOZ-2026-08-12.md  dated novelty + verification sweep

  -- Part N: Kishony–Fowler colour-code circuit distances, d = 3 to 15 (v0.16.0) --
colorcode-dcirc-paper/        note (Markdown + PDF) + FIXLOG.md + PROVENANCE.md
colorcode-dcirc-certificates/ 24 objects, 25 witnesses, 16 lower-bound records, 11 shipped LRAT
                              proofs (5 more attached to the v0.16.0 release and regenerable,
                              REGENERATE.md) + stdlib checkers (run_all.py, certs/check_dcirc.py,
                              certs/lemma_b_check.py) + certs/tamper_demo.py (15 controls) +
                              tools/ (CNF emitter, binary-LRAT converter, untrusted witness
                              lifter); the build's own CNF encoder and the solver drivers
                              are not included
SWEEP-RECORD-COLORCODE-DCIRC-2026-10-20.md  dated novelty sweep
```

For Part A, the reduction library and the system generators are intentionally not part of this repository; the published claims are the certificates themselves plus the verification scripts, which are self-contained. For Part B, the generating pipeline **is** included (`qec-scripts/certify.py`, `qec_lib.py`, `run_all.sh`) precisely because it is *not* trusted: it can be deleted and every certificate still verifies.

## Licensing

Dual license by content type:

- **Code and machine-readable certificate files** — everything under `scripts/`, `certificates/`, `qec-scripts/`, `qec-certificates/`, and the corresponding `*-scripts/` and `*-certificates/` directories of the later parts (`tt3-`, `qec1435-`, `cfr-`, `mps-`, `k34-`, `wedge-`, `kelmans-`, `k34add-`, `wedge2-`, `demag-`, `zefoz-`, `colorcode-dcirc-`) — are licensed under the **Apache License 2.0** ([LICENSE-CODE](LICENSE-CODE)).
- **Documentation and the paper** — `paper/`, `README.md`, `PROVENANCE.md`, and all other prose — are licensed under **CC BY 4.0** ([LICENSE-DOCS](LICENSE-DOCS)).

## Platforms

The **trusted base is operating-system-independent**: every checker in this repository is Python,
with no OS-specific calls on any default verification path — the one opt-in exception is noted
below — and verification requires nothing else. **Part A** (`scripts/`) is the one part whose
checkers need a third-party package — standard library plus SymPy. **The checkers of every other
part, B through N, are standard library only.** Checked import-by-import: Part B
`qec-scripts/check_witness.py`, `check_lower.py`, `check_duality.py`, `check_prof.py`,
`verify_manifest.py`; Part C `tt3-scripts/verify_minimizer.py`, `verify_sweep.py`; Part D
`qec1435-scripts/verify_1435.py`; Part E `cfr-certificates/verify_cfr525.py`; Part F
`mps-certificates/reverify.py`; Part G `k34-scripts/verify_witness.py`, `lrat_check.py`; Part H
`wedge-certificates/check_wedge.py`, `check_wedge_d5.py`; Part I `kelmans-scripts/verify_cert.py`,
`refcert.py`; Part J `k34add-scripts/verify_witnesses.py`, `verify_qr7_lemma.py`, `blowup_bound.py`;
Part K `wedge2-certificates/check_wedge2.py`; Part L `demag-certificates/check_demag.py`, `demag-scripts/check_remainder_bound.py`, `series_terms.py`; Part M
`zefoz-certificates/zefoz_checker2.py`, `zefoz_checker_pilot.py`; Part N
`colorcode-dcirc-certificates/run_all.py`, `certs/check_dcirc.py`, `certs/dem_from_circuit.py`,
`certs/extcode_from_circuit.py`, `certs/purex_from_object.py`, `certs/lemma_b_check.py`. Between
them they import nothing outside the Python standard library (re-confirmed by AST scan over all
23 files, 2026-08-12; Part N's six on 2026-10-29; Part L's `check_remainder_bound.py` and
`series_terms.py` on 2026-11-10).

Outside Part A, five scripts use third-party packages, and none is in any trusted base. Two are
auxiliary cross-checks: `mps-certificates/xcheck.py` (numpy), an independent cross-check of Part F's
stdlib checker, and `zefoz-scripts/anchor_check.py` (mpmath), Part M's independent 60-digit anchor.
Three are generators this README already labels untrusted: `qec-scripts/certify.py` and
`qec-scripts/qec_lib.py` (numpy), and `qec1435-scripts/gen_generic.py` (sympy). Deleting all five
leaves every certified claim verifiable. One further caveat, documented since the Part B audit: the
only non-portable call anywhere in a checker is the `os.mkfifo` in `check_lower.py`'s **opt-in**
`--external BINARY` path, which hands LRAT replay to a compiled third-party checker that this
repository deliberately does not vendor. That path is outside the trusted base by construction (see
[INDEPENDENT-VERIFICATION.md](INDEPENDENT-VERIFICATION.md) §5.1), was never invoked for any published
number, and is not needed to verify anything shipped here. On Windows, substitute the platform's usual
forms: `venv\Scripts\python` for `venv/bin/python`, backslash paths, and
`python -m gzip -d <file>.lrat.gz` where the quickstart uses `gunzip`. The `for` loops in the
quickstarts are POSIX-shell; on Windows run the listed commands individually or use WSL or
Git Bash. C sources (`tt3-scripts/tt3pack.c`, `qec1435-scripts/check1435.c`) are portable C
and are **not** part of any trusted base — the Python verifiers stand alone.

**Regeneration** (as opposed to verification) uses solver toolchains — msolve, CaDiCaL,
nauty/gentourng — that are routinely built on Linux and macOS; on Windows we recommend WSL
for regeneration. Verifying the shipped certificates never requires them.

Development and the pre-release replay audit were performed on macOS. **Windows: tested.** On
2026-08-07 an independent replay on native Windows CPython 3.14.0 (not WSL) ran nine checks
across Parts A–D — the Part B witness/lower-bound/duality checkers, the Part C minimizer
verifiers, the Part D control (`verify_1435.py`, expected distance-4 verdict with exit 1 by
design), and the Part A SymPy suite including `min_verify.py I` after the SIGALRM-portability
fix — all nine at their documented exit codes. Not yet covered on Windows: the large bb288/bb360
LRAT replays and the full n = 10 tournament sweep. That replay predates Parts G–M and did not include Parts E and F,
which were released the day before it; those checkers are standard-library-only and
OS-independent by inspection but have not been independently replayed on Windows. Part N was run on
2026-10-28 on GitHub-hosted Linux and Windows runners under Python 3.9 and 3.14; see
`colorcode-dcirc-certificates/README.md`, Platforms. Reports from other platforms are welcome and will be credited.

## Citing

See [CITATION.cff](CITATION.cff). Archival DOIs are minted per release on Zenodo — concept DOI for all versions [10.5281/zenodo.21799111](https://doi.org/10.5281/zenodo.21799111); Part A (v0.1.x) [10.5281/zenodo.21799112](https://doi.org/10.5281/zenodo.21799112); Part B (v0.2.0) [10.5281/zenodo.21799780](https://doi.org/10.5281/zenodo.21799780); Part B update (v0.2.1) [10.5281/zenodo.21831995](https://doi.org/10.5281/zenodo.21831995); Part C (v0.3.0) [10.5281/zenodo.21816010](https://doi.org/10.5281/zenodo.21816010); Part D (v0.4.0) [10.5281/zenodo.21816018](https://doi.org/10.5281/zenodo.21816018); Part E (v0.5.0) [10.5281/zenodo.21831896](https://doi.org/10.5281/zenodo.21831896); Part F (v0.6.0) [10.5281/zenodo.21832028](https://doi.org/10.5281/zenodo.21832028); Part G (v0.7.0) [10.5281/zenodo.21890619](https://doi.org/10.5281/zenodo.21890619); Part H (v0.8.0) [10.5281/zenodo.21895825](https://doi.org/10.5281/zenodo.21895825); Part I (v0.9.0) [10.5281/zenodo.21897011](https://doi.org/10.5281/zenodo.21897011); Part J (v0.10.0) [10.5281/zenodo.21898266](https://doi.org/10.5281/zenodo.21898266); Part K (v0.11.0) [10.5281/zenodo.21898343](https://doi.org/10.5281/zenodo.21898343); Part L (v0.12.0) [10.5281/zenodo.21898722](https://doi.org/10.5281/zenodo.21898722); Part L erratum (v0.12.1) [10.5281/zenodo.21910159](https://doi.org/10.5281/zenodo.21910159); Part M (v0.13.0) [10.5281/zenodo.21898996](https://doi.org/10.5281/zenodo.21898996); Part N (v0.16.0) [10.5281/zenodo.23792064](https://doi.org/10.5281/zenodo.23792064); Part L correction and addendum (v0.17.0) [10.5281/zenodo.24123940](https://doi.org/10.5281/zenodo.24123940); the Part A erratum (v0.1.2) is [10.5281/zenodo.21831894](https://doi.org/10.5281/zenodo.21831894). To cite an individual result, cite its note and the matching version DOI. External timestamps for this repository's claims begin at the first public push and the Zenodo deposits — not at local file dates (see PROVENANCE.md §3).

## Contact

Daniel Kirtchakov — daniel@halfounce.io
