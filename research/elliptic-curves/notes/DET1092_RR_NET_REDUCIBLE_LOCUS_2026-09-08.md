# The fixed first-centre RR net: complete reducible locus

## Result

**Verified application and new deduction.** The reducible members of the
historical first-centre RR net are exactly:

1. **22 unordered pairs of known generic sections:** one pair of heights0,10
   and21 pairs of heights4,6;
2. the projective pencil `C_min+F_tau`, consisting of the prior rational
   bisection and one fibre.

The22 section pairs supply no new independent point on specialization.
At302 the old bisection is nonsplit. When `tau=0`, the vertical component is
the entire302 fibre: it contains all its points but does not construct one.
No claim is made that this old bisection is useless at other parameters.

**Unknown.** Irreducible singular members, including members with rational
or genus-one normalization, remain possible. This result completes the
reducible part of one fixed net; it does not complete the arithmetic
explanation of302's rank gain.

The [v2 certificate](../../artifacts/generated-results/elliptic-curves/det1092_rr_net_reducible_locus_v2.json)
lists all44 minimum vectors and all22 explicit pairs as integer words in the
generic17 basis. The [independent replay](../../artifacts/generated-results/elliptic-curves/det1092_rr_net_reducible_locus_replay_v2.json)
uses a separately implemented exact closed-ball enumeration in a different
column order. The two enumerations visit1719 and1721 nodes respectively.

## The finite classification

**Established input, verified application.** The full geometric MW lattice
`M` and every section are rational; the surface has only irreducible fibres.
The selected word `w` has norm10. The net's divisor class is

\[
D=2O+5F+\phi(w)=C_{\min}+F=O+P_w.
\]

For a section `P_x`, the known divisor formula is

\[
P_x=O+\frac{\langle x,x\rangle}{2}F+\phi(x),
\]

also valid at `x=0`, where the section is `O`.

**New deduction.** Suppose a member splits horizontally as two sections
plus `k` fibres. Its two words must be `x,w-x`, and

\[
\langle x,x\rangle+\langle w-x,w-x\rangle=10-2k.
\]

Completing the square gives

\[
\boxed{\langle2x-w,2x-w\rangle=10-4k.}
\]

The exact minimum of the single coset `w+2M` is10, so `k=0`. Each minimum
vector `v` gives `x=(w+v)/2`; the pair `v,-v` gives the same unordered
section pair. The closed radius10 ball contains exactly44 vectors of this
parity, hence exactly22 pairs. There are no even-word double-section members
because `w` is not even.

For an irreducible horizontal bisection with `k` vertical fibres, its class
is `D-kF` and its arithmetic genus is `2-2k`. Nonnegative arithmetic genus
forces `k<=1`; at `k=1` the unique member is the already constructed
`C_min`. This exhausts all possibilities for a reducible degree-two divisor.
No larger parity census is required.

**Constructive check.** The member `O+P_w` has an explicit vertical-chord
equation through `T=P_{-w}`:

\[
-\operatorname{num}x_T+\operatorname{den}x_T\,x=0.
\]

Its exact coefficients in the RR basis `A,tA,B` are in the certificate and
are independently checked. The22 other word descriptions are constructive:
the generic section formulas and elliptic group law determine both points
of each pair. Completeness is a lattice/divisor theorem, not a bounded
rational-point-search claim.

## Consequence for the first successful coordinate

**New deduction from the completed independent-point replay.** The
[RR-net bridge](DET1092_FIRST_UNLOCK_RR_NET_2026-09-08.md) sends the first
saved successful coordinate to a finite value `u0` such that

\[
B+(u_0+vt)A
\]

contains the first exceptional point above zero for every rational `v`.
None of these affine members can be one of the22 section pairs: that would
place the first point on a specialized generic section, contradicting its
independent rank18 certificate. Nor can an affine member belong to the
projective line `span(A,tA)`, since its coefficient of `B` is1.
Consequently these affine members are geometrically irreducible.

This is a **retrospective structural deduction**, not a source-only rule
for choosing `u0`. No curve was selected or factored using that coordinate
in this experiment. It makes the remaining issue precise: useful lower-genus
members through the first point, if present, must be irreducible degenerations
of this pencil, rather than the known reducible members just classified.

The [genus-nine singular-member gate](DET1092_RR_NET_SINGULAR_MEMBER_GATE_2026-09-08.md)
controls the net's nondegenerate singular-incidence component. A subsequent
[whole-pencil obstruction](DET1092_FIRST_WITNESS_PENCIL_GENUS_GATE_2026-09-08.md)
now proves that every rational affine member through this first witness has
genus2. Thus no lower-genus member of this fixed pencil supplies the point;
other points, representatives and divisor systems remain unresolved.

## Evidence boundary correction and replay

**Rejected claim, retained evidence.** The immutable v1 certificate's final
sentence incorrectly excluded useful points on *every* smooth fibre.
That overbroad sentence is rejected. The v2 wrapper retains the identical
enumeration and equations and corrects the scope: only the22 section pairs
are uniformly old-subgroup points; the repeated bisection is known nonsplit
at302, not at every parameter. The independent replay explicitly records
this rejection. Only v2 supports the current theorem.

Limits were one parity coset,200000 exact nodes per enumeration, and one
canonical RR member check. All jobs finished in seconds. No point search,
parameter sweep, whole-orbit census, class-group calculation, or modification
to V3 was performed.

```sh
sage -python research/elliptic-curves/cas/verify_det1092_rr_net_reducible_locus.sage
```
