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

1. **Exploit the completed `(72,108)` log replay.**  The primary proof yields the
   smooth-branch rays `(2,1),(3,1),(4,1)`, while the Laurent transformations
   have longer adapted base ideals `(t,x^4),(t,x^6),(t,x^8)`.  Their residual
   source chains compile in `4,6,8` blowups, and additive composition proves
   that all three cases share the order-four eight-blowup graph.  Its
   `F_4` transition and affine-plane fill give a unimodular `10 x 10` source
   boundary passing adjunction.  Alternative residue labels are encoded
   symbolically, and the unselected order-three factor is proved to separate
   from the common order-four center and the filled divisor.  Target-infinity
   pole orders are exact on all eight exceptionals.  The full common-graph
   pole vector has no dicritical component, so at least one additional
   global cluster is forced.  The first weighted-Wronskian equation now
   selects a base-order-ten crossing blowup at `E3 intersect E4` followed by
   ten simple children, and excludes the formerly numerical smooth-`E3`
   candidate.  In terminal Case 2, two further target-infinity blowups have
   poles `12,0` and produce a degree-twelve dicritical; the combined intrinsic
   package passes with remaining self-intersection `29`.  On the final
   plane-return edge, the top bracket forces `A=a*r^2,C=c*r^3` for a quartic
   `r`.  The edge alone permits five homogeneous root partitions, and every
   corresponding intrinsic model passes.  The primary split-factor formula
   is sharper: it forces `r=(s-beta)^4`.  Selecting the other order-three
   factor gives the exact adapted transverse parameter and resolves the
   second blowup with poles `12,0` and degree-twelve residue.  Thus Case 1
   selects the same 23-component numerical package as Case 2.  The
   chain-to-boundary translation is complete, but the resulting package is
   allowed by every current intrinsic gate;
   see [`FRONTIER_LOG_SCALE_AUDIT.md`](FRONTIER_LOG_SCALE_AUDIT.md).
   The transformed pair has `[P,Q]=X^2`, so its actual ramification is
   `K+3H+div(X^2)`, not the boundary-supported `K+3H` representative.
   The corrected dicritical normal indices are `3` in Case 1 and `5` in
   Case 2; the total ramification intersection remains `35`.  The exact
   degree-twelve residue has a priori normalization-cover degree `1,2,4`,
   with image degrees `12,6,3`.  Case 1 retains all three rows, with
   degree-29 valuation contributions `3,6,12`; its composition sieve first
   requires the omitted bands `P:z^-6,...,z^-8` and
   `Q:z^-5,...,z^-12`, or a proof that they cannot affect the remainder
   equations.  In Case 2, general
   polynomial-right-component remainder ideals exclude cover degrees `2`
   and `4` already from `(J3),(J2)` and the determined part of `(J1)`.
   The remaining degree-twelve row is also empty: after localization at its
   forced endpoint `G_12 != 0`, the seven residual `(J1)` compatibility
   cubics generate the unit ideal, without `(J0)`.  Thus the former
   contribution `5` and remainder `24` are no longer live rows.  The earlier
   bottom-equation reduction, writing
   `H=gcd(C',G')`, `C'=H*c`, `G'=H*g` reduces `(J0),(J1)` to
   `B=K*c,F=K*g` and
   `2H(Ag-cE)+K^2(cg'-c'g)=0`, with `deg(K)<=deg(H)+1` and `K=0` or `t|K`.
   Its initial coefficients force `t|H`, eliminating gcd degree zero.  The
   maximal gcd degree `7` is also excluded: three coefficients of
   `remainder(G',C')` plus the terminal `(J0)` coefficient generate the unit
   ideal over the exact degree-35 field, with no residual `(J1)`
   compatibility equation.  The same extreme-stratum method excludes gcd
   degree `6` using a reconstructed linear cofactor, two terminal
   `G' mod H` coefficients, and `(J0)_{19}`.  Its remaining exact gcd
   degrees `1,...,5` are now retired as pre-compatibility strata.  The
   degree-one row is already reduced to the single origin pattern
   `ord(B),ord(E),ord(G'),ord(F)=(1,2,3,3)`; its only other valuation pattern
   is excluded by `(J1)`.  After those coefficient checks comes an accounting
   of the remaining divisorial valuation degrees.
2. **Finish the F2 lower-boundary theorem.**  The chain arithmetic, Puiseux
   translation, transformed bracket, forced vertices, common-power band, and
   terminal type-I normalization are now exact and machine-readable; see
   [`F2_75_125_DERIVATION.md`](F2_75_125_DERIVATION.md).  What remains is to
   prove the `gamma` branches exhaustive and determine every lower Laurent
   boundary and band through bracket layer 4.  The compiler deliberately
   rejects the partial record.
3. **Retire the excluded `(96,144)` repeated-tail branch.**  The source
   statements are now reconciled in
   [`FRONTIER_CLOSING_ATTACKS.md`](FRONTIER_CLOSING_ATTACKS.md).  The 2017
   table is a necessary over-approximation whose source leaves the
   lower-corner filter commented out; the 2016 full transition is stated only
   in a remark.  Reusing the proved common-tail part of the 2022 `(8,32)`
   calculation gives `q1=d0=4` and reduces the `(8,40)` vertical factor to a
   cubic.  Every cubic root partition containing a simple root yields the
   forbidden corner `(8,4)`.  The remaining factor
   `R=kappa*x^2*y^7*(y-lambda)^3` translates the initial edge to
   `(8,40)->(8,12)`.  Its maximum complete-chain length is three; even with a
   permissive superset of possible lower corners, the open-chain counts are
   `1,6,3,0` and no final corner occurs.  This repeated-tail row is excluded
   before Laurent compilation.  The other five `(96,144)` rows remain live.
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
| \((96,144)\) | 48; \(2:3\) | \((8,40)\to(8,28)\to(11/4,7)\), \((m,n)=(3,2)\); four \((2,3)\) chains from \((12,36)\), through respectively \((12,33)\to(11/3,8)\), \((9,24)\to(11/3,8)\), \((21/4,9)\to(19/4,8)\), and \((21/4,9)\to(12/4,5)\); and \((12,36)\to(12,30)\to(16/3,10)\to(11/6,3)\), \((m,n)=(3,2)\) | The repeated-tail row is excluded: simple-root partitions give `(8,4)`, and the triple-root translation `(8,40)->(8,12)` has no complete chain.  The other five rows still have no Laurent systems. | five high-difficulty rows remain | no |
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

The reconciliation is now complete.  The common tail forces a cubic residual
factor.  Simple-root partitions give `(8,4)`; the sole triple-root partition
gives the translated edge `(8,40)->(8,12)`, whose permissive complete-chain
enumeration has counts `1,6,3,0` and no final corner.  The source-aware
compiler records the row as pre-excluded and refuses to manufacture a Laurent
block.  See [`NEWTON_DERHAM_COMPILER.md`](NEWTON_DERHAM_COMPILER.md) and
[`cas/complete_chain_no_escape.py`](cas/complete_chain_no_escape.py).
