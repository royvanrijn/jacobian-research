# F7 triple-root pole curves: complete7/8/9-point strata

14 November 2026. Claim **GAOF7TRP1**.

For precisely Gao F7 over C on c!=0, normalize the target and mark an
eliminant root s of multiplicity at least three. On the entire carrier

    D0=s*(3s^3-2)=0,

the marked locus is four rational curves parametrized by a!= -1/2.
Its complete fibre classification is:

| Actual F7 fibre | Normalized targets on these curves | Eliminant type |
|---|---:|---|
| 7 points | **288**,72 at each carrier root |287 of type3+2+1^7; the fivefold origin5+1^7 |
| 8 points | **23** |4+1^8 |
| 9 points | All others; infinite ordinary open on each curve |3+1^9 |

Thus there are exactly **311 nonordinary normalized targets** in this
pole locus. Every target has a unique triple-or-higher root, so none
is counted twice. In particular the remaining three-point type
**(3,2,2,2) cannot lie on D0=0**. The generic D0!=0 case, global
three-point existence, full F7 strata,F5 and JC2 remain open.

[verify.py](verify.py) passed12884 checks in45.28seconds; its
[receipt](receipt.json) pins every mathematical input. It checks all
norm coefficients at1510 integer points using FLINT and separate
Python Bareiss determinants, and four fixed integer-resultant
nonvanishing tests modulo a trial-proved prime in two engines.
Independent end-to-end, formal, external and novelty review are not claimed.

## Exact source and complete rational coverage

Use the monic degree-twelve eliminant from
[Gao 2608.00222v2](https://arxiv.org/html/2608.00222v2):

    b(t)=Y4-3Y2*t^2+2Y1*t^3-t^4,
    a(t)=Y3-2Y2*t+Y1*t^2-t^3-t^2*b(t),
    R(t;Y)=a(t)^2+a(t)+t*b(t)+t^3+Y1*t-Y2.

Exactly its simple roots lift, once each, to actual F7 source points.
The full multiple-root critical parameterization X(s,a,b), source
coverage and its inverse a=a(s;Y),b=b(s;Y) are proved in
[GAOF7FOUR1](../gao-f7-quadruple/PROOF.md). In particular this is an
auxiliary marked incidence, not ramification in the etale F7 source.

The new checker freshly reconstructs X and R(s+d;X), proves its constant
and linear coefficients zero and its a,b recovery, and checks

    R(s+d;X)=c2*f2(s,a,b)*d^2+d^3*Q9(d;s,a,b),

where c2 is a nonzero rational constant and Q9 is monic of degree nine.
The exact source and all1213 Q9 terms are in
[pole-source-candidate.json](pole-source-candidate.json). Its legacy
candidate label grants no authority; the receipt and replay do.
Thus f2=0 is exactly the triple-or-higher incidence.

Write f2=A*b^2+B*b+C. The source identity gives A=4D0^2. On D0=0,

    B=2*(2a+1)*(6s^3-1).

The last factor is -1 at s=0 and3 at3s^3=2, hence is never zero.
If a=-1/2, the exact remainder is

    f2(s,-1/2,b)=W(s)/108 mod D0,
    W=27-181s-60s^2-144s^3.

The complete7x7 integer Sylvester matrix of D0,W has determinant
-2098783467, nonzero in native and Python exact arithmetic. Hence W
never vanishes at a carrier root. There are no source points with
B=0 and no omitted free-b component. Consequently every actual pole
point, and every a!= -1/2 at each of the four carrier values, is given by

    b=-C(s,a)/B(s,a).

This proves both source directions and the entire affine curve domain.

For each fixed carrier root this affine source curve is isomorphic to
Gm: write C(a)=c+(2a+1)*T(a), where c=C(-1/2)!=0. Then f2=0 expresses
the inverse of2a+1 as a polynomial in a,b. Its coordinate ring is
C[a,1/(2a+1)]. The recovered a(s;Y),b(s;Y) are polynomial in Y for
fixed s. Hence the coordinate-ring map of X is surjective and its
restriction is a closed embedding of this Gm curve into normalized
target space. The gcd bound below makes the four images disjoint.
These are the affine fibres; the I4 configurations in
[GAOF7CUSP1](../gao-f7-three-point/PROOF.md) concern its minimal
projective elliptic completion.
No finite source target is represented by an omitted point at infinity:
a(s;Y),b(s;Y) recover finite a,b from any finite Y and fixed s.

## Polynomial residuals without coordinate denominators

Since Q9 is quadratic in b, form

    U=B^2*[b^0]Q9-B*C*[b^1]Q9+C^2*[b^2]Q9.

At the recovered source point, U=B^2*Q9. On each constant base
s=0 or3s^3=2, coefficient reduction and constant clearing give a
polynomial F9(d,a,s) of d-degree9,a-degree8 and s-degree below the
base degree. The checker verifies the entire identity

    clearing*U-content*F9 = base_polynomial*quotient

by exact coefficient division, not samples. At s=0, F9=U/4; on the
cubic base F9=729*U/4. Their d-leading coefficients are respectively

    (2a+1)^2,       6561*(2a+1)^2.

They are nonzero at every retained point, so F9 has exactly the root
multiplicities of the actual residual Q9. Only rational constants and
the already excluded coordinate pole were removed.

## A rank bound rules out two residual repeated roots

For F9 degree9 and F9' degree8, take the polynomials

    d^k*F9  (k=0..6),       d^k*F9' (k=0..7).

Their15 coefficient rows in the16-dimensional space of polynomials
of degree at most15 form a rectangular Sylvester map. If their gcd
has degree at least two, every row is divisible by that gcd, so the
row space has dimension at most14. Every15x15 minor must vanish.
In particular the two minors using the first14 descending coefficient
columns and one of the last two must both vanish.

On s=0 these are integer polynomials in a. On3s^3=2 take their norms
over K=Q[s]/(3s^3-2). For u=x+y*s+z*s^2, multiplication by u in the
basis1,s,s^2, scaled by3, is

    [ 3x  2z  2y ]
    [ 3y  3x  2z ]
    [ 3z  3y  3x ].

Replace every entry of a15x15 minor by this3x3 block. The resulting
45x45 integer polynomial matrix has determinant
3^45 times the field norm of the original minor. Thus if a minor
vanishes at any one complex cubic-base value, this integer norm
vanishes. Taking separate norms may introduce extra common zeros,
but cannot discard a necessary zero.

The two norm polynomials on each base are in
[lower-order-norms-candidate.json](lower-order-norms-candidate.json).
After removing only a nonzero integer content and a power of2a+1,
their degrees and coprimality certificates are:

| Base | Two stripped degrees | Full integer Sylvester size | Determinant mod1000000007 |
|---|---:|---:|---:|
| s=0 |60,58|118|861345170|
| 3s^3=2 |180,180|360|932400796|

These are the Sylvester matrices of the ACTUAL stored integer
polynomials. A nonzero residue of their fixed integer determinant
proves that determinant nonzero over Z, hence coprimality over C.
The checker constructs every entry from the verified coefficients,
proves the prime by trial division, and checks each determinant in
FLINT and separate Python Gaussian arithmetic. This is not a
characteristic-zero inference from a modular unit ideal or a presumed
module rank.

Since2a+1 was already excluded, the two minors cannot vanish together
at any retained point. Therefore

    deg gcd(Q9,Q9') <= 1

EVERYWHERE on the four curves. In particular Q9 has no triple root
and cannot have two double roots. Its only possible repeated-root
pattern is one double root, with all remaining roots simple.
This stronger result includes the desired exclusion of three residual
double roots at every carrier value.

## Exact discriminant projection and the288 targets

The full17x17 Sylvester matrix of F9,F9' detects a residual repeated
root. Its determinant, and its51x51 scaled regular-representation norm
on the cubic base, factor exactly as

    nonzero_constant*(2a+1)^16*P72(a)      (s=0),
    nonzero_constant*(2a+1)^48*P216(a)     (3s^3=2).

P72,P216 are primitive integer polynomials of their indicated ACTUAL
degrees. The checker verifies these complete coefficient products.
It also proves neither has a root at -1/2 by exact homogeneous evaluation.
Their squarefreeness is proved by the full fixed integer Sylvester
matrices with their derivatives:

| Projection | Sylvester size | Determinant mod1000000007 |
|---|---:|---:|
| P72 |143|460917130|
| P216 |431|592069571|

Both determinants are nonzero in both arithmetic engines, so these
are respectively72 and216 distinct complex a-values, without any
heuristic root finding or factorization assumption.

For s=0 each P72 root gives the actual source b=-C/B and a residual
double root. For P216, the base is the three distinct constant roots
of3s^3-2, and the norm is a product of the three specialized
resultants. Every root therefore gives at least one actual base point.
Its norm multiplicity is one, so exactly one specialization vanishes
there and does so simply. Thus the216 projection roots count exactly
216 actual marked points, not merely an upper bound or a necessary
finite cover. Their distribution is72 per cubic base root: the
irreducible rational cubic has conjugate defining polynomials, so
all three specialized residual degrees and root counts agree.

The rank bound shows these288 points each have exactly one residual
double root. On the D0=0 quadruple-or-higher subset,
GAOF7FOUR1 completely classifies24 targets: the fivefold origin and23
ordinary quadruple targets. Therefore a residual double root cannot
coexist with a simply marked quadruple root: that would give type4+2+1^6,
which does not occur in that pole lemma. The only case where the double
is at the marked root is the known fivefold origin. All other287
points have exactly type3+2+1^7 and seven-point F7 fibres. The origin
also has seven points. It is the a=0 root of P72.

The23 ordinary quadruple targets have squarefree residual Q9 with a
simple root at d=0, giving eight-point fibres. Every other pole point
has an exact triple marked root and squarefree Q9 avoiding d=0, hence
nine simple eliminant roots and a nine-point fibre. Each curve loses
only finitely many a-values, so its ordinary nine-point locus is infinite.

No target can have a second triple-or-higher root: it would be a
triple root of Q9 and violate the proved gcd-degree bound. Combined
with unique source-coordinate recovery, this makes the marked
parameter-to-target map injective, even across different carrier values.
The counts288+23=311 are therefore distinct normalized targets.
The [exact presentation](parameterization.json) supplies both projection
polynomials, recovery equations and the old quadruple pole list.
Original unnormalized targets include their C* scaling curves. These
are complex affine geometric counts, without a scheme-length assertion.

## Why the determinant identities are exact

The generator uses integer interpolation but is not trusted. The checker
constructs every rectangular coefficient map independently from the
verified F9 coefficients and every scaled regular-representation block.
All entries have a-degree at most eight. Thus an N-by-N determinant has
a-degree at most8N, independently of any observed cancellation.
For each proposed polynomial it checks8N+1 distinct integer values in
both FLINT and separate Python Bareiss arithmetic. The bounds are:

| Calculation | Matrix size | Degree bound | Integer points per polynomial |
|---|---:|---:|---:|
| s=0 order1 minors, two |15|120|121|
| cubic order1 norms, two |45|360|361|
| s=0 full resultant |17|136|137|
| cubic full resultant norm |51|408|409|

There are1510 points altogether. Complete interpolation proves every
integer polynomial identity. Specializations where the d-leading
coefficient vanishes remain part of the nominal determinant identity;
the source point a=-1/2 is excluded by a separate exact source argument.
No multivariate resultant, symbolic discriminant backend or factorizer
is used by the checker. The four large modular determinants concern
fixed actual-degree integer matrices, with the characteristic-zero
argument stated above.

## Global remaining case and retained attempts

The [campaign](campaign.json) retains the successful bounded generators,
the earlier order2-minor exploration, the complete replay and all
subsequent global working inputs. The old six-variable factorization
runs were capped and are not restarted unchanged.

A fresh eliminant coefficient calculation in
[derive_global_fivevar.py](derive_global_fivevar.py) rechecks the full
factorization R=(t-r)^3*T3^2*S3 and all old equations before introducing

    delta=-2h-3r*u+3r*z-u*z+2v+z^2,
    h=v+((z-u)*(3r+z)-delta)/2.

Its degree-seven coefficient is linear in w, with coefficient a
nonzero rational constant times delta. On delta!=0 solve for w and
retain every clearing factor; on delta=0 retain w. This gives two
five-variable working systems, with1677 and13708 total terms.
Both branches are retained. The delta-zero branch returned a modular
unit in7.94seconds, while the saturated delta-nonzero pilot hit its
20second cap. These are triage only; no rational identity or global
exclusion follows from them. The30second characteristic-zero lift of the delta-zero branch timed out
without exporting an identity. A second support target dropped the653-term
constant-coefficient equation; its remaining six equations returned a
modular unit in8.99seconds, but did not improve the measured unit-search
time. No unchanged rational retry followed. Neither branch is a registered
theorem.

GAOF7FOUR2 already excludes(4,3,2). This proof removes every D0=0
case of the sole remaining list(3,2,2,2). The D0!=0 global problem
still needs a literal rational identity, a proved finite exceptional
presentation or a certified example. Nonrational elliptic geometry
and finite-field unit outputs do not decide it. Full F7 strata,F5
and JC2 remain open.

Reproduce this theorem from the repository root:

```sh
timeout 120s .venv/bin/python research/jacobian/gao-f7-triple-doubles/verify.py
```

Dependencies are pinned in [requirements.txt](requirements.txt). The
[primary source review](source-review.json) checked Gao v2 and de Bondt
v1 again on14 November; no theorem from a newer paper was imported.
The [status authority](../../MATH_STATUS.json) records the scope and
assurance separately from the retained generation candidates.
