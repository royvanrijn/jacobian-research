# Complete F7 quadruple-root strata and the exclusion of (4,3,2)

14 November 2026. Claim **GAOF7FOUR2**. Mathematical status is maintained
in [MATH_STATUS.json](../../MATH_STATUS.json).

For precisely Gao's F7 over C, on c!=0 and after weighted target
normalization, every target whose eliminant has a root of multiplicity
at least four has exactly one of these types:

| Complete degree-twelve root list | Normalized targets | Actual F7 fibre |
|---|---:|---:|
| 4+1^8 | Infinite, ordinary open | 8 points |
| 4+2+1^6 | **912** | 6 points |
| 5+1^7 | **118** | 7 points |

Thus there are exactly **1030 nonordinary normalized quadruple-root
targets**, and the three-point type **(4,3,2) does not occur**.
Together with GAOF7THREE1, GAOF7CUBE1, GAOF7MULT1 and GAOF7FIVE1,
only **(3,2,2,2)** remains for a three-point fibre. Its existence, the
full F7 strata, F5 and JC2 remain open.

The [complete replay](verify_complete.py) runs the source and high-root
dependencies, a certified integer norm reconstruction, the residual
discriminant identity, and every presentation-pole chart. Its
[receipt](complete-receipt.json) binds the component checkers and receipts.
The proof combines exact arithmetic with the local norm argument below.
Two arithmetic engines are used; independent end-to-end review, formal
verification, external review and novelty priority are not claimed.

## Source and the finite plane chart

Use the degree-twelve monic eliminant and source reconstruction in
[Gao 2608.00222v2](https://arxiv.org/html/2608.00222v2) and the fully
replayed [GAOF7FOUR1 source lemma](PROOF.md). A simple eliminant root
lifts to exactly one F7 source point; multiple roots do not lift.
On the marked quadruple incidence,

    f2(s,a,b)=f3(s,a,b)=0,
    R(t;X(s,a,b))=(t-s)^4*Q8(t;s,a,b),

where Q8 is monic of degree eight. Every multiple root is covered by
the critical parameterization X; no root-coordinate division is used.
The [source input](quadruple-input.json) stores f2,f3 and all Q8 coefficients.

Write D=s*(3s^3-2), f2=A*b^2+B*b+C, A=4D^2, and

    L1=P1*b+Q1=(A*f3-[b^2]f3*f2)/(4D),
    E1=A*Q1^2-B*P1*Q1+C*P1^2.

Removing only nonzero constants and D-supported factors from E1 gives
G(s,a), of a-degree eight, with leading coefficient

    ell=2560*s^4*(3s^3-2)^2*(3s^3-1)^2*(3s^3+1).

On D*P1!=0 the source equations are equivalent to G=0 and b=-Q1/P1.
All source identities and this converse are proved in GAOF7FOUR1.
Its good fibre also proves G separable over Q(s). The marked incidence
has finite fibres over every s; none of its components is vertical.

## A residual octic with a polynomial discriminant

Set d=t-s and Q=Q8(s+d;s,a,b). Remove its quadratic b term:

    LQ=A*Q-[b^2]Q*f2=PQ*b+QQ,
    C_raw=P1*QQ-Q1*PQ.

The fresh integer polynomial replay proves exactly

    C_raw=8D*C0(d,a,s),
    LC_d(C0)=D*P1/2.

On f2=f3=0 with D*P1!=0 this gives

    C0(d,a,s)=(D*P1/2)*Q8(s+d;s,a,b).

C0 has degrees (8,7,27) in (d,a,s), with1175 terms; it is stored as
equation1 in [generic-triple-input.json](generic-triple-input.json).
No denominator-zero chart is discarded by this identity.

Let delta(a,s)=Disc_d(C0). The [integer coefficient list](residual-discriminant-candidate.json)
has a-degree94 and s-degree305. Despite its legacy candidate filename,
its entire identity is checked in [verify_discriminant.py](verify_discriminant.py):

    Res_d(C0,partial_d C0)=LC_d(C0)*delta.

The left side is the nominal15x15 Sylvester determinant and has degrees
at most105 in a and405 in s. The right side obeys the same bounds.
At every point of the106-by406 integer grid, the determinant equals
the stored polynomial in both FLINT and separate Python Bareiss
arithmetic. These43036 points prove the complete polynomial identity,
including specializations whose leading coefficients vanish. A nonzero
polynomial LC_d(C0) can then be cancelled in the polynomial ring.
No symbolic discriminant backend is used by the checker.

## Exact norm and its factorization

Reduce delta in a modulo G by87 fraction-free steps. Each step is an
explicit identity

    ell*R_i-leading_a(R_i)*a^(deg_a(R_i)-8)*G
                =content_i(s)*R_(i+1).

The initial content is also retained. Every content has only s,
3s^3-2,3s^3-1,3s^3+1 as nonconstant factors. The checker divides each
literal polynomial identity exactly; it performs no gcd search.
The [final remainder](primitive-remainder-candidate.json) has a-degree7,
s-degree1435 and11465 terms. On the complement of the stated pole
carrier, R_final is a unit multiple of delta in the algebra defined by G.

Define the integer polynomial

    N(s)=Res_a(G(s,a),R_final(s,a)).

Its15x15 Sylvester matrix gives the degree bound11669. If || ||_1
means the sum of absolute values of every integer coefficient, each
coefficient of N has absolute value at most

    B=||G||_1^7*||R_final||_1^8.

Indeed expansion of the determinant is bounded coefficientwise by the
product of the row coefficient l1 norms; summing over all row choices
also bounds the subset of choices that are permutations. B has10308
bits. The product M of345 distinct trial-proved primes has10315 bits
and satisfies M>2B.

The [saved integer norm](crt-norm-candidate.json) is verified independently
of its CRT generator by [replay_norm.py](replay_norm.py) and
[verify_norm_modular.c](verify_norm_modular.c). For EACH prime, the full
nominal Sylvester determinant equals the proposed N at11670 distinct
field points, in FLINT and a separate elementary C Gaussian elimination.
Thus each modular polynomial identity holds by the degree bound. The
candidate coefficients also satisfy the same height bound, so congruence
modulo M>2B proves the integer identity coefficient by coefficient.
The replay performs8052300 determinant comparisons. It uses univariate
product-tree evaluation and no multivariate resultant routine.

Exact coefficient multiplication verifies

    N = kappa*s^1645*(3s^3-2)^805*(3s^3-1)^805
              *(3s^3+1)^448*B87(s)^28*H912(s),

where kappa is a nonzero integer and B87,H912 are the explicit primitive
integer polynomials in the
[factor certificate](crt-norm-factors-candidate.json). The degree is11167.
The [discriminant replay](discriminant-receipt.json) proves H912 has
actual degree912, is squarefree over C, and is coprime to B87 and all
four displayed pole factors. Each test uses its actual degree912,
a p-unit leading coefficient at p=1000000007, and literal inverses in
Fp[s]/H912, replayed with native and separate Python coefficient lists.
The finite-free unit argument gives the characteristic-zero meaning.
This is a proved integer norm, not an inference from a modular unit ideal.

## Why a simple norm zero excludes residual triple roots

Work off the D,ell and P1 pole loci. G/ell gives a finite free algebra
of rank eight over the base line. Since G is generically separable,
its curve is reduced on this chart, although special fibres may be
ramified or nonreduced. R_final and delta differ by a base unit in
this algebra. The resultant N and the algebra norm of delta therefore
differ by a base unit as well.

For a degree-eight polynomial with nonzero leading coefficient, a
triple root or two distinct double roots give a gcd of degree at least
two with its derivative. Its15x15 Sylvester matrix has corank at least
two. All its first cofactors vanish, so the resultant and all its first
coefficient derivatives vanish. After dividing by the nonzero leading
coefficient, the same is true of its discriminant. Consequently, at
any such point (s0,a0), the composed delta belongs to

    (s-s0,a-a0)^2.

This forces its norm to vanish to order at least two at s0. To see this
without assuming an unramified plane projection, pass to the
normalization of each local curve branch. Both coordinate differences
have positive valuation there, so delta has valuation at least two.
For a finite separable function-field extension over C((s-s0)), the
valuation of a norm is the sum of branch valuations, weighted by
residue degrees; the residue fields here are C. The determinant norm
has the same value. Ramification does not divide this order. Other
branches contribute nonnegative valuations because delta is regular.
The same argument shows two distinct geometric points with delta=0
would make the norm vanish to order at least two.

Off the pole carriers, every zero of N is a simple zero of H912.
Thus each such s-value has exactly one geometric point with delta=0,
and its Q8 has exactly one double root and six other simple roots.
No residual triple root or two residual double roots can occur there.
This is the step that turns a discriminant projection into an exclusion;
nonzero generic norm alone would not suffice.

## Every presentation-pole chart is covered

The D=0 chart is the complete24-target classification in
[EXCEPTIONAL_PROOF.md](EXCEPTIONAL_PROOF.md): the known fivefold origin
has a seven-point fibre;23 targets have one quadruple root and eight
other simple roots. Six marked points occur at each of the four D roots.

For P1=0 on D!=0, L1=0 forces Q1=0. The full8x8 integer Sylvester
matrix of P1,Q1 in a gives

    Res_a(P1,Q1)=constant*s^2*(3s^3-2)^2*(3s^3-1)^3*B87(s).

B87 is exactly the degree87 factor in N. After the leading-pole values
are handled separately, every coordinate-pole point lies over B87=0.
Signed cofactors of row5, columns6,7 give a=N(s)/T(s).
Their full polynomial identities are checked at119 integer points
for the resultant and103 for each cofactor, with two determinant engines.
The raw denominator is a unit on B87, so the adjugate has rank seven
and every common a-root is this recovered value. Exact homogeneous
identities prove P1=Q1=0 for every retained root, not just necessity.
The base is squarefree, avoids D and ell, and all recovery denominators
are units. The actual monic f2 equation in b gives a finite-free algebra
of rank87*2=174. The original cancellation identity then gives f3=0.

At3s^3-1=0 and3s^3+1=0 the actual a-degree of G drops to seven,
with leading coefficients -38440/9 and -10728 respectively.
The degree-three base and monic degree-seven a equation give two
finite-free presentations of rank21. P1 is a unit in each, so b=-Q1/P1
and the source equivalence still apply. No a-leading zero is omitted.

In EACH of these three presentations, the
[chart certificate](remaining-charts-candidate.json) gives literal
inverses for the defining-polynomial derivative and Q8(s). It also
gives a polynomial inverse to partial_t Q8 modulo the monic Q8.
[verify_remaining_charts.py](verify_remaining_charts.py) reconstructs
every source coefficient and checks all identities with native nested
quotients and separate Python coefficient lists. The actual monic
presentations over Z_(p) are free of ranks174,21,21; adjoining monic Q8
remains free. Invertibility modulo p therefore makes the corresponding
multiplication determinants p-units, hence units in characteristic zero.
All216 targets have one quadruple root and eight other simple roots.
Squarefree base and extension derivatives give the asserted geometric
counts. Their unique quadruple root makes the source-to-target map
injective, and the three charts are disjoint. With D=0, all240
presentation-pole targets are covered: one seven-point and239 eight-point.

## Exact count and the remaining three-point gate

Every one of the912 distinct roots of H912 avoids all poles. The
specialized determinant is zero, so there is a common root of G and
R_final, and the plane model gives an actual marked source point.
The simple norm zero proves there is exactly one such geometric point.
Its residual octic has one double root and otherwise simple roots.
GAOF7MULT1 bounds root multiplicities by five, and GAOF7FIVE1 proves
every target with a fivefold root has type5+1^7. The residual double
root therefore cannot be the marked root, and no fifth-order marked
root occurs here. These targets have exactly type4+2+1^6.

Each has a unique quadruple root, so different s-values cannot give
the same normalized target. This proves the count912, not merely a
projection-degree upper bound. The explicit polynomial and unique
common-root recovery are recorded in
[finite-degenerations.json](finite-degenerations.json).

At every other pole-free point the residual octic is squarefree.
It either avoids the marked root (type4+1^8) or contains it simply
(type5+1^7). The pole charts have already been classified. The complete
fivefold theorem counts exactly118 of the latter targets, all with
seven-point fibres. The two exceptional lists are disjoint, giving
912+118=1030. Only simple roots lift, which gives the table's actual
F7 fibre counts.

In the original unnormalized target space, each normalized point
corresponds to a C* curve under weighted scaling. The counts above
are normalized complex geometric targets, without a scheme-length
assertion. They do not settle all F7 fibre strata. The separate
triple-root type(3,2,2,2), global cardinality three and JC2 remain open.

## Reproduction, provenance and limits

From the repository root, with Python3.12.3, SymPy1.14.0,
python-flint0.9.0, a C compiler, and matched FLINT3.3.1 headers/library:

```sh
timeout 360s .venv/bin/python research/jacobian-followup/gao-f7-quadruple/verify_complete.py
```

The modular norm helper uses the already installed Julia FLINT3.3.1
and OpenBLAS artifacts; set FLINT_ROOT and FLINT_BLAS_LIB for other
installations. Its compiler invocation and versions are recorded in
[norm-grid-receipt.json](norm-grid-receipt.json). The C helper is
GPL-2.0-or-later and uses FLINT under its upstream licence. The
Python witness generator is not trusted by the replay.

The [campaign](campaign.json) retains the direct integer norm timeout,
CRT reconstruction, serialization and representation failures, and all
earlier capped identity attempts. The failure of the direct integer
resultant motivated a finite, coefficient-bounded CRT reconstruction;
no cap was enlarged to continue the same failed search. Files retaining
candidate names are accepted only through their hash-bound complete
replays, never by their filenames or generation logs.

The immediate [source review](source-review.json) confirms Gao v2
(4 November) and de Bondt v1 (17 October) are still current as checked
on14 November. The restricted search found no newer imported theorem
for this computation; it is not an exhaustive novelty or priority review.
Credit for the eliminant and simple-root mechanism remains with Gao;
the original user paper is
[de Bondt 2610.24889](https://arxiv.org/abs/2610.24889).
This proof addresses the explicit F7 family and does not prove JC2.
