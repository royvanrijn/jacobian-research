# F7 quadruple-root pole chart:24 targets, with fibres seven and eight

14 November 2026. Pole-chart lemma for **GAOF7FOUR1**; see the [full theorem](PROOF.md).

Use exactly Gao F7 over C on c!=0, with normalized targets
Y=(cA,c^2 B,c^3 D,c^4 Z). Here s denotes a marked eliminant root, not
the target coordinate D. Exactly24 distinct normalized targets have a
root of multiplicity at least four satisfying

    D0(s)=s*(3s^3-2)=0.

One is the known fivefold origin Y=0 and has a seven-point fibre. The
other23 have a single quadruple root and eight other simple roots,
hence eight-point fibres. There are six marked incidence points at each
of the four complex roots of D0. In the original c!=0 chart these give
24 disjoint C* curves through the normalization formula.

The [exact parameterization](exceptional-parameterization.json) uses a
squarefree degree 24 polynomial in z=s+a and rational recovery formulas.
The [current verifier](verify.py) passed
[1275-check pole-only snapshot](history/exceptional-receipt.json) in 9.22 seconds, including
the complementary slope-coordinate identities. The
[certificate](exceptional-certificate.json) is 74296 bytes.

This excludes(4,3,2) only when its quadruple root is on D0=0. The
generic D0!=0 branch and the other remaining three-point list(3,2,2,2)
are **OPEN**. No full quadruple-root classification or global three-point
exclusion is claimed. The statements count geometric points, without
asserting a scheme structure, independent end-to-end or formal review,
external review or novelty priority.

## Why this exceptional chart is necessary

For the critical-root parameterization X(s,a,b) below, direct
coefficient comparison gives

    v=a'(s;X)
     =2*(3s^3-1)*(a^2+a)+2*s*(3s^3-2)*b
        +6s^6+4s^5-2s^3-3s^2.

The derivative holds the target X fixed. Thus on D0!=0 the change
(s,a,b)->(s,a,v) is birational, with inverse

    b=(v-2*(3s^3-1)*(a^2+a)-6s^6-4s^5+2s^3+3s^2)/(2D0).

After this substitution every critical target coordinate is at most
quadratic in a; the earlier cubic terms cancel. The verifier checks
the exact derivative identity and this degree assertion for all four
target coordinates. It also checks that all five
[generic slope equations](slope-triple-input.json) come from the
[raw residual-triple equations](residual-triple-input.json), with only
nonzero constants and powers of s and3s^3-2 removed. This supplies a
verified working system on D0!=0, rather than a global exclusion.

The two source equations in slope coordinates have a-degree three,
instead of the earlier degrees four and five. These coordinate gains
do not justify discarding any pole or claiming a faster exact lift.
The present theorem completes the entire D0=0 piece required by that
change of variables.

## Full root-incidence coverage and the residual octic

From [Gao 2608.00222v2](https://arxiv.org/html/2608.00222v2) use

    b(t)=Y4-3Y2*t^2+2Y1*t^3-t^4,
    a(t)=Y3-2Y2*t+Y1*t^2-t^3-t^2*b(t),
    R=a(t)^2+a(t)+t*b(t)+t^3+Y1*t-Y2.

Every simple root lifts uniquely to a source point on c!=0.
At a marked multiple root put a=a(s), b=b(s), and define

    G2=a^2+a+s*b+s^3,
    G3=2s*G2+a+s^3+s^2*b,
    G4=3s^2*G2+b+s^4,
    X1=-det(d(G2,G3,G4)/d(s,a,b)),
    X=(X1,G2+s*X1,G3+s^2*X1,G4+s^3*X1).

The root equation and definitions force
Y=(xi,G2+s*xi,G3+s^2*xi,G4+s^3*xi), with xi=Y1. Direct substitution
gives R'(s;Y)=xi-X1, so every multiple root is covered by X. There is
no variable denominator in this reconstruction. The checker replays
these identities by integer polynomial arithmetic.

The two primitive polynomials f2,f3 in
[quadruple-input.json](quadruple-input.json) are nonzero rational
multiples of the Taylor coefficients of orders two and three in
R(s+x;X). Coefficients zero and one vanish identically. The verifier
also reconstructs the complete monic quotient identity

    R(t;X)=(t-s)^4*Q8(t;s,a,b)
                   +c2*f2*(t-s)^2+c3*f3*(t-s)^3,

with c2,c3 nonzero rational constants and Q8 monic of degree eight.
Consequently f2=f3=0 is equivalent to a root of multiplicity at least
four, for all source values. The quotient has 892 coefficient terms.

For the remaining global(4,3,2) gate, the raw system is f2=f3=0 plus
Q8(s+d)=Q8'(s+d)=Q8''(s+d)=0. All d, including zero, are retained.
The slope transformation is checked against precisely these equations.

## Reduce the whole pole carrier before eliminating

Reduce f2,f3 modulo D0(s), then clear only nonzero rational constants.
The two retained polynomials F2,F3 in
[exceptional-source.json](exceptional-source.json) are LINEAR in b:

    F2=B2(s,a)*b+C2(s,a),
    F3=B3(s,a)*b+C3(s,a).

Every point on D0=0 with f2=f3=0 satisfies these equations, and the
converse holds there. Form B2*C3-B3*C2 and reduce it modulo D0, again
clearing only rational constants. Call the resulting integer polynomial
H(s,a); its a-degree is six. This linear-resultant identity is replayed
from the source by the checker, without a symbolic resultant backend.

Set z=s+a and h(s,z)=H(s,z-s). The determinant

    U(z)=Res_s(D0(s),h(s,z))

is a nonzero rational constant times the primitive integer polynomial
P24(z) in the certificate. The s-degrees are4,9 and the z-degrees are0,6,
so the13x13 Sylvester determinant has z-degree at most24. The checker
verifies the full integer polynomial at 25 distinct integers with both
FLINT determinants and a plain-Python fraction-free Bareiss algorithm.
It preserves nominal degrees at all specializations.

The coefficient rows of this matrix annihilate(s^12,...,s,1) at any
common finite s-root. Use the signed cofactors deleting row 9 and columns
11,12, denoted N(z),T(z). They have degrees17,18 and satisfy

    N(z)-s*T(z)=0

at every common root. At rank 12 the adjugate column is proportional to
the kernel vector; at rank at most11 both cofactors vanish, which proves
the same relation. Each cofactor has degree at most18, since deleting
the first h-row leaves only three z-dependent rows of degree at most six.
Both signed polynomial determinants are checked at 19 distinct integers,
in both engines.

## Finite-free units give characteristic-zero conclusions

Apply the finite-free good-reduction lemma proved in the
[complete fivefold-root note](../gao-f7-quintuple/PROOF.md). For an integer
polynomial P of actual degree n whose leading coefficient is nonzero
modulo ell=1000000007, the algebra

    Z_(ell)[z]/(P/leading_coefficient(P))

is free of rank n. A nonzero modular determinant of multiplication by an
element makes that element a unit over Z_(ell), hence nonzero at every
complex root of P. Previously certified coordinate inverses stay in the
same algebra, so their reduction is legitimate. No assumption of reducedness
is needed for this lemma.

The checker retains the actual integer degrees and proves each leading
coefficient an ell-unit, proves ell prime by trial division, and verifies
these multiplication determinants in FLINT and independent Python:

| Algebra | Tested element | Dimension | Residue |
|---|---|---:|---:|
| Q[z]/P24 | T |24|167006305|
| Q[z]/P24 | P24' |24|932964196|
| Q[z]/P24 | B2(s,a) |24|400935167|
| Q[z]/P24 | Res_t(Q8,Q8') |24|121642958|
| Q[z]/P23 | Q8(s) |23|358538936|

Here P24=z*P23, P23(0)!=0, s=N/T, a=z-s and b=-C2/B2.
Thus both coordinate denominators are units, and no zero-denominator
point is discarded. A genuine marked source point projects to a root
of P24 and the cofactor relation recovers s; the linear F2 equation then
recovers b.

Every modular source evaluation and multiplication matrix is separately
reconstructed with Python coefficient-list polynomial arithmetic.
The finite resultant of Q8,Q8' is a 15x15 determinant with coefficient
z-degrees at most23, hence degree at most345. It is verified at 346
distinct field values in two determinant engines. This checks the entire
finite-field polynomial identity, not a sampled numerical test.
The finite-free lemma then makes its norm a conclusive characteristic-zero
unit certificate. In particular Q8 is squarefree at every retained point.

A general modular unit-ideal report would not justify this passage to C.
The fixed finite FREE presentation and its checked leading coefficient
are indispensable.

## Positive coverage and the exact count24

The checker also proves the homogeneous substitutions

    T^4*D0(N/T)=0 mod P24,
    T^9*h(N/T,z)=0 mod P24

as exact INTEGER polynomial identities, by division and full coefficient
multiplication; the multiplication is also replayed in plain Python.
Since T is a unit, every root z of P24 really gives D0(s)=H(s,a)=0.
Since B2 is a unit, setting b=-C2/B2 gives F2=0 and

    F3=(B2*C3-B3*C2)/B2=0.

Thus every retained root is an actual f2=f3=0 point, not merely a
necessary projection value. The squarefreeness norm of P24 proves its
24 complex roots distinct, and the recovery formulas give one point
for each root.

If two distinct marked source points gave the same Y, either their s
values differ, producing another root of multiplicity at least four
in the squarefree residual Q8, or their s values agree, in which case
root-incidence reconstruction forces the same a,b and therefore the
same z. Both possibilities contradict distinctness. Hence there are
exactly 24 normalized TARGETS as well as24 marked incidence points.

D0 has four distinct complex roots. For any fixed one of them,
H(s,a) is a nonzero polynomial of degree at most six: if it vanished
identically, the nonzero resultant U(z) would vanish for every z.
There are therefore at most six points over each s. The exact total24
forces equality at every one of the four values. This yields six
points at s=0 and eighteen over the three roots of 3s^3=2.

At z=0 the exact recovery is s=a=b=0 and Y=0. Here

    R=t^5*H7, H7=t^7-2t^4+2t-1.

H7(0)=-1 and disc(H7)=713, giving the known seven-point fibre.
For the other23 roots, the last norm gives Q8(s)!=0. Each has one
quadruple root and eight distinct other simple roots, so its fibre has
eight source points. No scheme-structure assertion follows from this
geometric point count.

## Assurance, retained attempts and the remaining gate

The origin axis fibre belongs to Gao's prior work; the earlier exact
degree-five eight-point witness in GAOF7CARD1 also lies on s=0. The
present statement identifies the COMPLETE24-point marked pole chart,
without a novelty-priority claim. The local source and arXiv record
remain [Gao v2](https://arxiv.org/abs/2608.00222v2), dated 4 November 2026;
the original user paper remains
[de Bondt v1](https://arxiv.org/abs/2610.24889v1), dated 17 October 2026,
at the 14 November version check. No new external theorem is imported.

With [pinned dependencies](requirements.txt), from the repository root:

```sh
timeout 120s .venv/bin/python research/jacobian-followup/gao-f7-quadruple/verify.py
```

Python3.12.3, SymPy1.14.0 and python-flint0.9.0 were used. The current full-theorem
receipt binds the checker and its source/certificate inputs by
SHA256. The replay takes about 9.22 seconds, including the verified
complementary source transformations.

The [campaign record](campaign.json) retains the inconclusive global
attempts. The direct four-variable system produced a modular unit ideal
in 9.54 seconds. Two function-field lifts and both modular/rational
D0^2 lift attempts reached20 second caps. No identity was exported.
A local modular DAG hit its30 MiB raw-output cap after7.87 seconds;
streaming compression addressed that specific I/O failure, but the
bounded retry reached90 seconds with 1343 nodes and no terminal unit.
It has about 2.39million polynomial terms. No completed DAG or
characteristic-zero inference is claimed. Singular runs had1536 MiB
address-space and10 MiB output caps; DAG pilots had4 GiB and30 MiB caps.
All traces, inputs, logs and partial evidence are retained. No automatic
prime-generation campaign is scheduled.

Modulo32003, reduction puts D0^2, though not D0, in the slope-system
ideal. This is triage only. A literal rational identity expressing
D0^2 in that ideal, or a different exact argument or example on D0!=0,
is still needed. Together with this completed pole chart it would
settle the remaining(4,3,2) gate. No unchanged capped retry is proposed.
The other three-point list(3,2,2,2), global three-point existence, full
strata, F5 and JC2 remain open. There is no independent end-to-end,
formal or external review. MATH_STATUS.json is the status authority.

The pole-only snapshot used GAOF7QEX1 as an internal identifier; it was never
registered. The current theorem is GAOF7FOUR1 with [receipt.json](receipt.json).
