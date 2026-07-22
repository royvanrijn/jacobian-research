# Candidate frontier from 125 through 150

> **Bounded computational regression.**  The full complete-chain algorithm is
> published as pseudocode in Guccione--Guccione--Horruitiner--Valqui,
> arXiv:1708.07936.  No author-supplied implementation or source repository was
> located.  The local script therefore does not claim an independent
> reimplementation: it encodes every relevant row of the paper's accepted
> Section 7 tables, recomputes \((m(a+b),n(a+b))\), checks the printed maximum
> and coprimality, quotients by coordinate reversal, and asserts the exact
> distinct output.  This is a deterministic table regression.

Run:

```bash
python3 plane-jc/cas/frontier_125_150.py
```

It returns 13 normalized unordered pairs:

\[
\begin{gathered}
(75,125),(84,126),(96,128),(88,132),(90,135),\\
(56,140),(84,140),(96,144),(108,144),\\
(42,147),(63,147),(98,147),(100,150).
\end{gathered}
\]

In particular, 125 is an admissible **maximum** in the enumeration, through
\((75,125)\); it is not a claim that a Keller counterexample exists.

## Compiler handoff and loose ends — 2026-07-22

The reusable middle end is now implemented and tested:

```text
published exhaustive Laurent normal form
    -> Laurent lattice supports and bracket layers
    -> weighted-Wronskian block
    -> superelliptic character eigenspace
    -> canonical de Rham obstruction coordinates
    -> optional exact local-system reuse certificate
```

For `(72,108)`, one typed normal-form certificate expands the published chain
to both Laurent polygon cases.  Both compile to the same genus-three first
block, and the six compatibility equations are certified coordinates in the
six-dimensional compact de Rham eigenspace.  The implementation and current
boundary are documented in
[`LAURENT_BAND_FRONTEND.md`](LAURENT_BAND_FRONTEND.md) and
[`NEWTON_DERHAM_COMPILER.md`](NEWTON_DERHAM_COMPILER.md).

Open work, in dependency order:

1. **Coordinate before starting the F2 derivation.**  Another active agent may
   already be working on the `(75,125)` F2 `j=1` normal form; inspect current
   worktrees/branches and avoid duplicating that derivation.
2. **Derive or obtain an exhaustive F2 `j=1` normal-form certificate.**  The
   required output is every Laurent polygon branch, the normalizing
   transformations, vertex nonvanishing, coefficient normalizations, full
   supports, and transformed bracket monomial.  The compiler deliberately
   rejects the current corner-only record.
3. **Reconcile the `(96,144)` sources.**  Settle the relation between the 2016
   lower-side exclusion of `(8,28),(8,40)` and the 2017 raw complete-chain row
   before treating repeated-tail reuse as a live exclusion problem.
4. **Test reuse only after a new band derivation exists.**  Derive the new
   curve polynomial `A`, give an explicit base/coordinate map, pass the exact
   cyclic-curve identity check, and then compare the forcing differential,
   primitive support, residual scaling, and gauge/exact class.  Equal coarse
   fingerprints alone are not evidence of reuse.
5. **Optional engineering simplification.**  Emit one machine-readable report
   per normal-form branch containing provenance, polygons, bands, bracket
   layers, Wronskian data, de Rham certificate, and any reuse certificate.

Scope guard: this handoff records a compiler programme, not a new family
exclusion or bound.  The externally published bound remains 108; the local
125 statement retains its explicit dependence on the published normal-form
reduction and local `(72,108)` certificate.

## Ranked worklist

“System” below means the coefficient system still to be derived; none is
silently treated as solved.

| Pair | gcd / ratio | first corner and chain data | Remaining obstruction | Expected difficulty | Does the \((72,108)\) proof apply unchanged? |
| --- | --- | --- | --- | --- | --- |
| \((75,125)\) | 25; \(3:5\) | family F2, \(A_0=(5,20)\), final \((7/5,2)\) | derive its Laurent polygons and coefficient systems | high, but shortest next family case | no: different ratio, corner, bands, and approximate-root exponents |
| \((84,126)\) | 42; \(2:3\) | \(A_0=(7,35)\to(19/7,5)\), and a second realization \((12,30)\to(16/3,10)\to(11/6,3)\) | two chain realizations must both be excluded | high | no: the ratio agrees but neither chain is \((8,28)\to(11/4,7)\) |
| \((96,128)\) | 32; \(3:4\) | F24, \((8,24)\to(14/4,6)\to(5/4,0)\to(19/8,3)\) | longest family-derived corner chain at this maximum | very high | no |
| \((88,132)\) | 44; \(2:3\) | \((11,33)\to(19/4,8)\) | new \(2:3\) coefficient system | high | no |
| \((90,135)\) | 45; \(2:3\) | \((9,36)\to(17/9,4)\) with \((m,n)=(3,2),(2,3)\); \((9,36)\to(9,24)\to(11/3,8)\); and \((12,33)\to(11/3,8)\), both with \((2,3)\) | four published realizations; their Laurent polygons and coefficient systems remain to be derived | very high | no |
| \((56,140)\) | 28; \(2:5\) | F11, \((7,21)\to(13/7,3)\) | different approximate-root multiplicities | high | no |
| \((84,140)\) | 28; \(3:5\) | F9, \((7,21)\to(11/7,2)\) | family companion to the previously excluded \((56,84)\) | medium-to-high; prior \((7,21)\) work may help | no, but older \((7,21)\) reductions may be reusable |
| \((96,144)\) | 48; \(2:3\) | \((8,40)\to(8,28)\to(11/4,7)\), \((m,n)=(3,2)\); four \((2,3)\) chains from \((12,36)\), through respectively \((12,33)\to(11/3,8)\), \((9,24)\to(11/3,8)\), \((21/4,9)\to(19/4,8)\), and \((21/4,9)\to(12/4,5)\); and \((12,36)\to(12,30)\to(16/3,10)\to(11/6,3)\), \((m,n)=(3,2)\) | six rows in the 2017 table; no Laurent systems derived.  The repeated-tail row must first be reconciled with the 2016 lower-side remark that discards \(B_0=(8,28),B_1=(8,40)\) via the impossible corner \((8,4)\) | source reconciliation first | not unchanged; reuse is presently an architecture experiment, not a live exclusion |
| \((108,144)\) | 36; \(3:4\) | \((8,28)\to(7/4,3)\) | same first corner as the audited case but a different final corner and exponent ratio | high | no |
| \((42,147)\) | 21; \(2:7\) | F7, \((6,15)\to(7/3,4)\) | large degree ratio and approximate-root exponent | high | no |
| \((63,147)\) | 21; \(3:7\) | F8, \((6,15)\to(8/3,5)\) | companion chain to F7 | high | no |
| \((98,147)\) | 49; \(2:3\) | \((7,42)\to(13/7,6)\), both orientations | new \(2:3\) chain | high | no |
| \((100,150)\) | 50; \(2:3\) | \((10,40)\to(16/5,6)\to(23/10,3)\) and \((10,40)\to(18/5,8)\to(8/5,3)\), both with \((m,n)=(3,2)\) | two distinct Laurent/residual systems expected but not yet constructed | very high | no |

The table records every one of the 24 published chain realizations in this
window.  A “remaining system” is deliberately recorded as *not yet constructed*
when the cited paper supplies corner data rather than explicit Laurent
coefficient equations; the number of coefficient systems need not equal the
number of chains until the polygonal branching is derived.

## Does the method eliminate a family?

No family theorem follows from the reproduction.  The pair
\((72,108)=(2d,3d)\) with \(d=36\) is controlled not merely by its gcd or
ratio, but by the specific admissible tail

\[
 (8,28)\longrightarrow(11/4,7),
\]

the resulting two Laurent polygons, their band widths, and the degree-35
first-block field.  The frontier contains many other \(2:3\) pairs, and their
corner chains differ.  Even \((108,144)\), which shares \(A_0=(8,28)\), uses
ratio \(3:4\) and final corner \((7/4,3)\), so the coefficient system changes.

A defensible parameterized compiler conjecture is:

> A repeated admissible tail should permit reuse of its valuation/band descent
> after the earlier chain has been separately normalized.

The raw \((96,144)\) row through
\((8,40)\to(8,28)\to(11/4,7)\) initially looked like the first test.  A
source-level rereading changes that status.  In the 2016 lower-side paper,
the remark following Proposition 3.29 says that
\(B_0=(8,28),B_1=(8,40)\) leads to the impossible last lower corner
\((8,4)\) and can be discarded.  The 2017 Section 7 table nevertheless lists
the length-two complete chain.  The formal proposition proves the
impossibility of the relevant class of last lower corners; the application to
this particular `B_0,B_1` pair is stated in the remark as straightforward.

The two source statements must therefore be reconciled before the row is
treated as a live frontier case.  It remains a useful synthetic test of a
source-aware compiler: [`cas/newton_derham_compiler.py`](cas/newton_derham_compiler.py)
records the five missing band-level inputs and refuses to copy the old
genus-three block into the new row.  See
[`NEWTON_DERHAM_COMPILER.md`](NEWTON_DERHAM_COMPILER.md).
