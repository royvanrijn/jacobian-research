# Gao F7 quadruple-root incidence: generic degree eight and finite degenerations

14 November 2026. Registered claim **GAOF7FOUR1**.

Current result: [GAOF7FOUR2](COMPLETE_PROOF.md) now classifies the entire
quadruple-root locus:912 six-point and118 seven-point exceptional targets,
and excludes(4,3,2). This document retains the predecessor lemma used by
that proof. References below to an uncomputed projector or undecided
(4,3,2) describe this lemma's original scope; current status is in the
linked complete theorem.

For exactly Gao F7 over C, on c!=0, mark an eliminant root s of
multiplicity at least four. Its incidence has finite fibres over EVERY
s-value. It has generic degree eight over the s-line. The generic
normalized targets in this locus have exactly one quadruple root and eight other
simple roots, hence eight F7 source points.

The nonordinary normalized-target locus is finite. In particular, the
still-undecided three-point type(4,3,2) can occur at only finitely many
normalized targets. This theorem does not compute that complete finite
set or exclude it. The other remaining type(3,2,2,2), global three-point
existence, full F7 strata, F5 and JC2 remain **OPEN**.

The complementary slope coordinate has pole carrier D0=s*(3s^3-2).
On D0=0 the complete marked locus consists of 24 distinct normalized
targets: the known fivefold origin with seven-point fibre, and23
quadruple targets with eight-point fibres. There are six marked points
at each of the four carrier values. This pole chart has an
[exact parameterization](exceptional-parameterization.json).

The [separate verifier](verify.py) passed [1535 checks](receipt.json) in
10.30 seconds. It replays source coverage, the entire pole chart, a
coefficient coprimality determinant and a good finite-free fibre of
rank eight, with separate native and Python arithmetic. The
[generic model](generic-model.json) specifies the plane curve and norm
functions. No independent end-to-end or formal review, external review
or novelty priority is claimed. Counts are geometric, without a
scheme-structure assertion.

## Source equations and full root coverage

Use the monic degree-twelve eliminant in
[Gao 2608.00222v2](https://arxiv.org/html/2608.00222v2):

    b(t)=Y4-3Y2*t^2+2Y1*t^3-t^4,
    a(t)=Y3-2Y2*t+Y1*t^2-t^3-t^2*b(t),
    R(t;Y)=a(t)^2+a(t)+t*b(t)+t^3+Y1*t-Y2.

A simple root lifts uniquely to a source point. The full critical-root
parameterization X(s,a,b) and its converse are reconstructed in the
[pole-chart lemma](EXCEPTIONAL_PROOF.md): no nonconstant denominator
is used to cover all multiple roots.

The verifier proves that the two primitive polynomials f2,f3 in
[quadruple-input.json](quadruple-input.json) are nonzero constant
multiples of the Taylor coefficients of orders two and three in
R(s+x;X). Coefficients zero and one vanish identically, and exactly

    R(t;X)=(t-s)^4*Q8(t;s,a,b)
                  +c2*f2*(t-s)^2+c3*f3*(t-s)^3,

where c2,c3 are nonzero rational constants and Q8 is monic of degree
eight. Thus f2=f3=0 is precisely the marked quadruple-or-higher incidence.
The checker replays the full integer polynomial identity, not merely
sampled derivatives.

## The generic plane model

Write f2=A*b^2+B*b+C, with A=4s^2*(3s^3-2)^2. On D0!=0 set

    L1=P1*b+Q1=(A*f3-[b^2]f3*f2)/(4D0),
    E1=A*Q1^2-B*P1*Q1+C*P1^2.

The literal identity

    E1=P1^2*f2+L1*(A*Q1-A*P1*b-B*P1)

is checked from the fresh source. Clearing only rational constants and
powers of s and3s^3-2 gives the primitive integer polynomial G(s,a)
stored as equation3 in the
[retained bivariate input](../gao-f7-linear-cube/high-root-bivariate-input.json).
The checker verifies every normalization. Its a-degree is eight, and its
leading coefficient is exactly

    ell(s)=2560*s^4*(3s^3-2)^2*(3s^3-1)^2*(3s^3+1).

Whenever D0 and P1 are nonzero, the source equations are equivalent to
G=0 and b=-Q1/P1. Necessity follows by cancellation. Conversely, G=0
implies E1=0, so the displayed identity and L1=0 give f2=0. The
cancellation equation then gives f3=0 because A=4D0^2 is nonzero.
Zeros of P1 are not discarded globally: they remain in the source
incidence and are handled by the fibre bound below. The good-fibre test
will prove P1 invertible only over Q(s) and on its specified sample.

## No vertical component and finite fibres at every s

Write G=sum_(j=0..8) g_j(s)*a^j. Let u=g_0 and v=sum_(j=1..8)g_j.
The checker constructs these exact integer polynomials from G. Both
have degree 27. Their full54x54 INTEGER Sylvester matrix has determinant

    779574954 mod1000000007,

nonzero in both FLINT and plain-Python Gaussian elimination. Therefore
its integer determinant is nonzero and u,v have no common complex root.
In particular all the g_j can never vanish simultaneously. This checks
every s-value, including a-leading-coefficient zeros; it is not a
generic coefficient-gcd report.

For fixed s with D0!=0, every marked source satisfies the nonzero
polynomial G(s,a)=0, so at most eight a-values occur. For each a,
f2 is quadratic in b with nonzero leading coefficient A(s), so at most
two b-values occur. Hence there are at most16 marked points in such a
fibre. No sharpness claim is made for that bound.

On D0=0 the pole-chart lemma classifies the full24-point set and proves
exactly six points at each of the four values. Thus the marked-root
projection has finite fibres EVERYWHERE. There can be no positive-
dimensional component confined to one s-value. No properness assumption
about this projection is used.

## A good rank-eight fibre at s=1

Let K=Q(s) and A_K=K[a]/(G/ell). This is a degree-eight finite FREE
algebra. It need not initially be assumed irreducible or reduced.
At s=1, D0=1 and ell(1) is nonzero modulo p=1000000007. The
[good-fibre certificate](generic-certificate.json) uses the exact integer
specialization G(1,a), still of degree eight. Its four multiplication
norm determinants are:

| Element at s=1 | Rank | Residue modulo p |
|---|---:|---:|
| partial_a G |8|338024181|
| P1 |8|777603188|
| Res_t(Q8,Q8') |8|731948363|
| Q8(1) |8|133730184|

Each actual integer leading coefficient is checked a p-unit. The
finite-free lemma over Z_(p), used in the pole-chart proof, therefore
makes these conclusive characteristic-zero unit tests in the sample
algebra. Matrix entries, coordinate inverses and source evaluations
are checked with FLINT and separate Python coefficient-list arithmetic;
every determinant is checked in native and Python engines. The
prime is proved by complete trial division.

The sample Q8 coefficients have a-degree at most seven. The 15x15
Sylvester determinant with Q8' consequently has degree at most105 in a.
It is verified at 106 distinct finite-field values in both determinant
engines. This proves the entire finite polynomial computation. The
finite-free unit lemma supplies its characteristic-zero meaning.

The source identities above prove that every root of G(1,a) really
gives f2=f3=0 after b=-Q1/P1. The four units prove these eight roots
distinct, their marked roots exactly quadruple, and all eight residual
roots simple. Thus this is an exact good fibre, not a numerical or
finite-box example.

## Why one good fibre proves a generic statement

This passage is deliberately weaker than a global unit assertion.
Work first over

    B=Z_(p)[s,1/(D0(s)*ell(s))].

Since D0(1) and ell(1) are p-units, specialization s->1 followed by
reduction modulo p is defined. The algebra B[a]/(G/ell) is free of
rank eight. Multiplication determinants are rational functions of s.
The first two sample units show that the determinants for partial_a G
and P1 are not identically zero: a zero rational function would have
zero specialization wherever the presentation is defined.

Hence over K=Q(s), G is separable of degree eight and P1 is a unit.
After inverting its nonzero norm, b=-Q1/P1 is defined in this algebra.
The source Q8 coefficients then lie in the same localized finite FREE
presentation, still specializing regularly at s=1. The last two sample
units show that the rational functions

    Norm_(A_K/K)(Res_t(Q8,Q8')),
    Norm_(A_K/K)(Q8(s))

are also not identically zero.

Remove the finitely many zeros and poles of these norm functions,
D0 and ell. On the resulting nonempty Zariski-open subset U of the
s-line, G has eight distinct roots, P1 is nonzero at every one, and
both source directions of the model hold. Every one of these eight
marked points has one quadruple root and eight other simple roots.
Thus the incidence has generic degree eight, and the associated
normalized F7 targets all have eight-point fibres.

This does NOT say that the norm functions are units at every s.
Their exceptional zeros and poles have not been computed in full.
A nonzero sample can prove generic nonvanishing, but cannot remove
those complex exceptional values.

## Finiteness of the remaining quadruple-root degenerations

Every point whose eliminant is not of type4+1^8 lies over the finite
complement of U. Every fibre over that complement is finite by the
coefficient-resultant argument and the pole-chart classification.
Therefore all such marked points form a finite set. Their image in
normalized target space is finite as well. In particular, targets of
three-point type(4,3,2) are a finite, still-undecided subset.

The statement concerns normalized Y on c!=0. In original coordinates
one normalized point gives the C* curve
(c,Y1/c,Y2/c^2,Y3/c^3,Y4/c^4); finiteness of normalized exceptions does
not mean finiteness of their original five-dimensional target points.
The marked quadruple-root locus itself is generically one-dimensional,
with infinitely many ordinary targets. Only its nonordinary locus is
proved finite here.

The complete pole piece is 24 normalized targets: one seven-point
origin and23 eight-point targets, with six points at each D0-root.
The [parameterization](exceptional-parameterization.json) and
[integer certificate](exceptional-certificate.json) retain all of them.
They are not asserted to be the complete global degeneration list;
only the origin is nonordinary within this pole piece.

## Remaining computation and assurance boundary

The [verified generic slope system](slope-triple-input.json) is a
concrete necessary-and-sufficient working system on D0!=0 for a
quadruple marked root with a triple root in its residual octic. It
includes every root separation d, including zero. Modulo32003,
D0^2, though not D0, reduces to zero in its ideal. This remains triage.
The bounded rational and modular lifts exported no identity, and
neither local modular DAG pilot reached a terminal unit. The full
[campaign record](campaign.json) retains every cap and failure reason.
No characteristic-zero generic exclusion is claimed from those runs.

A next exact route is to expand or otherwise certify the nonzero
relative norm projector, construct its finite exceptional presentation
with proved degrees, and test the residual triple-root equations there.
Alternatively, a literal rational D0^2 membership identity for the
verified slope system, together with the pole lemma, would exclude
(4,3,2). Neither route is completed. An exact example remains equally
valid evidence. The separate(3,2,2,2) gate, full F7 strata, F5 and JC2
are open.

With [pinned dependencies](requirements.txt), replay from the repository root:

```sh
timeout 120s .venv/bin/python research/jacobian/gao-f7-quadruple/verify.py
```

Python3.12.3, SymPy1.14.0 and python-flint0.9.0 were used. The final
receipt binds every proof input and the generic curve source by SHA256.
The 1,535-check replay takes about 10.30 seconds, including all source
transformations. The proof uses no Groebner search, symbolic resultant
backend or factorizer. Arithmetic independence is not independent
end-to-end or formal review; no external review, author endorsement or
novelty priority is claimed. The source boundary and earlier witnesses
retain the credit recorded in the pole-chart lemma. MATH_STATUS.json
is the status authority.
