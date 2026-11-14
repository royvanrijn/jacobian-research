# F7: a partial three-point obstruction and an elliptic incidence surface

14 November 2026. Registered claims **GAOF7THREE1** and **GAOF7CUSP1**.

For exactly Gao's F7 over C, three-point fibres remain **OPEN**. This note
exhausts their possible multiplicity patterns, excludes one coefficient
subfamily by an explicit rational polynomial identity, and determines the
elliptic geometry of the auxiliary triple-root incidence surface. It does
not exclude either full factorization family or settle JC2. No independent
end-to-end review or novelty priority is claimed.

The subsequent [GAOF7CUBE1](../gao-f7-three-point-certificate/PROOF.md)
excludes the full cubic-cube family, including repeated/shared roots.
This note retains the original partial identity and elliptic calculation;
the later [GAOF7MULT1 proof](../gao-f7-linear-cube/HIGH_ROOT_PROOF.md)
excludes(7,2) by the sharp maximum-five bound. The subsequent
[GAOF7FIVE1 proof](../gao-f7-quintuple/PROOF.md) classifies all118 fivefold
targets and leaves only(4,3,2),(3,2,2,2) and global three-point existence open.

## Inputs and notation

The map and simple-root lifting criterion are from
[Gao, arXiv:2608.00222v2](https://arxiv.org/html/2608.00222v2).
[GAOF7BND1](../gao-f7-fibres/PROOF.md) counts every c=0 fibre: it has one
or two points, so three points can occur only on c!=0. In this chart put
Y=(cA,c^2 B,c^3 D,c^4 Z) and use the monic degree-twelve eliminant

    b(t)=Y4-3Y2*t^2+2Y1*t^3-t^4,
    a(t)=Y3-2Y2*t+Y1*t^2-t^3-t^2*b(t),
    R(t;Y)=a(t)^2+a(t)+t*b(t)+t^3+Y1*t-Y2.

Exactly its simple roots lift. GAOF7CARD1 supplies ten occurring global
cardinalities, and [GAOF7SURJ1](../gao-f7-surjectivity/PROOF.md) excludes
zero. Eleven is impossible, leaving only three undecided.

## Exhausting three-point multiplicity types

Let S3 be the monic product of the three simple roots. It is squarefree
and disjoint from every nonsimple factor. The remaining multiplicities
sum to nine and are at least two. Their odd entries number either one
or three. If there is one, subtract three from it and halve every
remaining multiplicity. The resulting exponents sum to three, giving

    R=(t-r)^3*T3^2*S3,                              (I)

where T3 is any monic cubic; its roots may repeat or coincide with r.
If there are three odd entries, they are all three, giving

    R=Q3^3*S3.                                     (II)

Here Q3 is monic and squarefree. The eight possibilities for the
nonsimple multiplicities are

    (9), (7,2), (6,3), (5,4), (5,2,2), (4,3,2),
    (3,2,2,2), (3,3,3).

The last belongs to (II), the other seven to (I). Conversely, a
factorization of either type with S3 squarefree and disjoint from the
other factors has exactly three simple roots. This reduction introduces
no genericity restriction.

## A coefficient gate for the cubic-cube family

Write

    Q3=t^3+u*t^2+v*t+w,
    S3=t^3+z*t^2+h*t+k,
    Delta=(u-z)^2+2(v-h).

We prove that **Delta=0 is impossible whenever R=Q3^3*S3**. This statement
holds even without the squarefreeness and disjointness conditions.

For P=Q3^3*S3, comparison of degrees 11,10,9 forces

    Y1=-[t^11]P/4,
    Y2=([t^10]P-4Y1^2)/6,
    k=-12Y1*Y2-2-([t^9]P)|k=0.

After this k substitution, degrees eight and six force

    Y4=(6Y1+9Y2^2-[t^8]P)/2,
    Y3=([t^6]P-14Y1*Y2+6Y2*Y4-2)/2.

Only nonzero constant denominators occur. The seven remaining coefficient
equations, in degrees 0,1,2,3,4,5,7, are the primitive integer polynomials
in [cubic_cube-input.json](cubic_cube-input.json). The checker reconstructs
them from R, matches each up to a nonzero rational constant, and checks
that every eliminated coefficient vanishes identically.

In the saved primitive degree-seven equation the coefficient of w is
576 Delta. On Delta=0 make the invertible substitution

    m=u-z, z=u-m, h=v+m^2/2.

The resulting seven necessary equations g_i in (u,v,m,w), again made
primitive over Z, are in
[cubic-delta-zero-input.json](cubic-delta-zero-input.json). The
[certificate](cubic-delta-zero-macaulay-certificate.json) gives seven
rational polynomials lambda_i, with 1397 terms in total, such that

    sum_i lambda_i*g_i = 1.                         (*)

This is a literal identity over Q. It is replayed both in FLINT and in
the retained standard-library Fraction sparse-polynomial implementation.
Thus the seven equations have no simultaneous complex zero, proving the
claimed exclusion. The fixed t coordinate is part of the statement;
no assertion of invariance under arbitrary coordinate changes is needed.

This1397-term identity alone excludes only the Delta=0 slice of (II).
The full family (II), including Delta!=0, is subsequently excluded by
[GAOF7CUBE1](../gao-f7-three-point-certificate/PROOF.md). The five
lists in (I) were not excluded by GAOF7CUBE1 alone; GAOF7MULT1 subsequently
excludes(7,2); GAOF7FIVE1 subsequently excludes both fivefold lists,
leaving only(4,3,2) and(3,2,2,2) unresolved.

## The auxiliary triple-root incidence surface

The next calculation describes the discriminant's incidence geometry;
it concerns an auxiliary critical sweep, not a locus of ramification
points in the etale source of F7. Introduce independent coordinates
(s,a,b) and set

    G2=a^2+a+s*b+s^3,
    G3=2s*G2+a+s^3+s^2*b,
    G4=3s^2*G2+b+s^4,
    X1=-det(d(G2,G3,G4)/d(s,a,b)),
    X=(X1,G2+s*X1,G3+s^2*X1,G4+s^3*X1).

Substitution gives R(s;X)=R'(s;X)=0 identically. Every derivative here
is taken with Y held fixed, before setting t=s and Y=X. Consequently

    f(s,a,b) = R''(s;X)/2 = 0                       (T)

is the triple-or-higher-root incidence condition on this parameterization.
The condition includes higher multiplicities. Define the surface by (T).

The polynomial f is quadratic in b, say f=Ab*b^2+Bb*b+Cb, with

    Ab=4s^2*(3s^3-2)^2.

Let H(s,a)=disc_b(f)/4 and Z=Ab*b+Bb/2. The exact identity

    Z^2-H = Ab*f

gives an isomorphism over K=C(s) between (T) and Z^2=H(s,a). The polynomial
H is quartic in a, with leading coefficient

    16s^2*(3s^3-2)^2*(3s^3-1)^2
      = (4s*(3s^3-2)*(3s^3-1))^2.

Its discriminant is nonzero. Thus the smooth projective generic fibre
has genus one and two K-rational points at infinity. Choosing either
as origin makes it an elliptic curve, isomorphic to its Jacobian.
The gcd in C[s] of all coefficients of f as a polynomial in a,b is a
nonzero constant. Hence there are no extra vertical components. Since
H is squarefree over C(s), Z^2-H is geometrically irreducible over that
field; Gauss's lemma then makes the whole incidence surface irreducible.

The exact quartic coefficients are retained in
[elliptic-cusp-input.json](elliptic-cusp-input.json), not reconstructed
from numerical samples. At s=1 its discriminant is 5799391526912.

## Invariants and the relatively minimal completion

For H=A*a^4+B*a^3+C*a^2+D*a+E use

    I=12AE-3BD+C^2,
    J=72ACE+9BCD-27AD^2-27B^2E-2C^3.

The Jacobian model is

    y^2=x^3-27I(s)*x-27J(s).

These formulas and the genus-one interpretation are the standard
binary-quartic construction, as in Section 2 of
[Fisher, On binary quartics and the Cassels-Tate pairing](https://www.dpmms.cam.ac.uk/~taf1000/papers/bq-ctp.pdf).
Our D4 below denotes the ordinary polynomial discriminant of H;
Fisher's genus-one discriminant is 16D4. With this convention,

    D4=(4I^3-J^2)/27
      = -4096s^4*(3s^3-2)^4*Q64(s),
    disc_Weierstrass = 16*3^12*D4.

The saved Q64 has degree 64, is squarefree, and is coprime to
s*(3s^3-2). The saved I and J have degrees 23 and 40. Exact gcd checks
show gcd(I,D4) and gcd(I,J) are nonzero constants. These facts determine
every finite bad fibre of the relatively minimal elliptic surface:

| Geometric place | Minimal discriminant order | Kodaira fibre |
| --- | ---: | --- |
| s=0 and the three roots of 3s^3-2 | 4 at each | four I4 |
| the 64 roots of Q64 | 1 at each | 64 I1 |

At these finite places c4 is a unit, so the displayed model is minimal.
At infinity put T=1/s and scale x=T^-14*X, y=T^-21*Y. The resulting
integral coefficients and discriminant have orders

    ord_T(a4)=28-23=5,
    ord_T(a6)=42-40=2,
    ord_T(disc)=84-80=4.

The a6 order is below six, so this model is minimal. In characteristic
zero these orders give type IV. The table and numerical invariants
follow the Kodaira criteria in Section 4.4 and the surface formulas in
Sections 6.7 and 6.10 of
[Schuett and Shioda, Elliptic surfaces](https://arxiv.org/pdf/0907.0298).
In particular the full geometric fibre configuration is

    4 I4 + 64 I1 + IV.

It has Euler number 4*4+64+4=84 and chi(O)=84/12=7. The j-function is
nonconstant: I^3 has degree 69 and D4 has degree 80 with no common
finite root. The base is P1 and there is a section. Hence q=0 and
p_g=chi(O)-1=6. Geometric genus is a birational invariant of smooth
projective surfaces, so this completion, and therefore the affine
incidence surface (T), is **not rational**.

The Kodaira descriptions refer to the relatively minimal projective
elliptic completion. They do not describe the uncompleted affine
fibres of f=0 as already consisting of those complete configurations.
No Mordell-Weil rank, K3 property, three-point exclusion or plane
Jacobian consequence follows from this calculation alone.

## Reproducibility, failed attempts and remaining gate

From the repository root:

```sh
timeout 60s .venv/bin/python research/jacobian-followup/gao-f7-three-point/verify.py
```

This performs 78 exact checks, including two arithmetic replays of (*),
independent reconstruction of all necessary equations, and the invariant,
discriminant and local-valuation calculations. The geometric implications
use the cited theorems and written arguments. Neither those arguments nor
the entire proof have received independent end-to-end or formal review.
Sources and the retained Fraction implementation are SHA256-pinned in
[sources/manifest.json](sources/manifest.json).

The certificate finder uses modular linear algebra only to select support,
then solves the 1397-by-1397 supported system over Q and accepts only after
an exact polynomial identity replay. Its weighted bound was 18. The
finder is not part of default verification.

Raw and weighted Singular attempts used 30-second time, 1536MB address-space
and 1MB output caps. Their inputs, logs and receipts are retained, including
an initial linear-family script rejected because a ring name collided with
a variable name, followed by its corrected timed-out run. msolve reports
and finite-field unit ideals are triage only; no theorem here is inferred
from them. Other bounded Macaulay attempts, including Delta-power targets
for the then-unresolved branch, found no identity at their retained bounds.
Failure to find an identity is inconclusive. They are not scheduled for
unchanged repetition.

After the separate full cubic-cube exclusion, a complete answer still
requires an exact exclusion or a certified example for the two surviving
lists(4,3,2),(3,2,2,2), after GAOF7MULT1 and the complete GAOF7FIVE1
fivefold-root classification. The elliptic model provides a structural constraint on
multiple-root incidence and rules out rational parametrization of the
entire surface, but no map from that fact to the remaining coefficient
ideal is asserted. Full F7 strata and F5 remain
open. MATH_STATUS.json is the mathematical-status authority.
