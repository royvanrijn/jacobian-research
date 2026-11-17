# Ambient height bounds for locally soluble fibres of the fixed map

This note applies quantitative results on thin sets to the same fixed quintic
Keller map `Phi`. The [ambient arithmetic criterion](FIXED_MAP_AMBIENT_LOCAL_CLASSIFICATION.md)
already decides local solubility at every rational target. Here we count those
targets in their actual projective height. The argument uses the exact
hypotheses of the external theorems identified below; finite group and
polynomial controls do not certify those imported theorems.

Write a target as `y=(u,b,c)`, put `pi=-u`, and retain its inverse polynomial

\[
 E(S)=\pi^5S^5-\frac32\pi^4S^4+\frac32\pi S^3
       +\left(-\frac54\pi-\frac9{32}b\right)S^2
       +\frac9{16}S-\frac9{32}c.                         \tag{1}
\]

The height throughout is `H(y)=H([1:u:b:c])`. Let `A(B),L(B),G(B)` count,
respectively, Hasse-failing, everywhere locally soluble, and rationally
soluble fibres at rational targets of height at most `B`. These are ambient
target counts, without quotienting targets by isomorphism of their fibres.
Exactly `L(B)=G(B)+A(B)`.

## Proposition ABC1: generic monodromy and the square-discriminant restriction

The monic polynomial `f=E/pi^5` has arithmetic and geometric generic group
`S_5` over `Q(pi,b,c)`, and its normal closure is regular over `Q`.
Let `U` be the open target set `pi Disc(E)!=0`. Every Hasse-failing target
lies in `U`; at such a target its root action is the intransitive action

\[
 H=(S_2\times S_3)\cap A_5\simeq S_3
 \quad\text{on a block of two and a block of three}.       \tag{2}
\]

In particular `Disc(E)` is a nonzero rational square there. For an
irreducible factorization `f=q h` of degrees two and three, this square
condition is equivalent to the common quadratic-resolvent condition in FMA3.
The union of the conjugates of (2) has exactly 36 elements of `S_5`, namely
the identity, the twenty three-cycles, and the fifteen double transpositions.

**Proof.** Over `K=Q(pi)`, the polynomial has the form `h_0(S)+B S^2-C`,
where `B,C` are independent variables, `deg h_0=5`, and
`h_0'(0)=9/(16pi^5)!=0`. The [quadratic-tilt lemma UAQTILT1](../jacobian/universal-keller-structure/PROOF.md)
gives both generic groups `S_5`, even after algebraically closing `K`.
It follows in particular that the geometric group over `Q` is already
`S_5`, equal to the arithmetic group; their constant-field quotient is
trivial. This is the regularity assertion.

FMA1 places every Hasse failure in `U`, and FMA2 gives a cubic field together
with its quadratic resolvent. A transposition on the cubic roots also swaps
the two quadratic roots, whereas a cubic three-cycle fixes them. This is
precisely the subgroup (2). Its elements are therefore all even. Alternatively,
the discriminant identity for monic coprime factors is

\[
 \operatorname{Disc}(q h)=\operatorname{Disc}(q)
       \operatorname{Disc}(h)\operatorname{Res}(q,h)^2.    \tag{3}
\]

Thus the discriminant is a square exactly when the factor discriminants have
the same squareclass. Replacing `f` by `E=pi^5 f` multiplies its discriminant
by `pi^40`, another square. The usual Vandermonde description of the
discriminant gives the same parity conclusion. The stated class sizes are
`1,20,15`, and all occur in conjugates of (2). \(\square\)

The square-discriminant condition alone is not sufficient for Hasse failure:
one still needs the factorization and the finite local tests in FMA3.

## Proposition ABC2: factor-selection covers of degrees ten and twenty

Over `U`, the scheme `X_2` selecting a monic quadratic factor of `f` is a
finite etale, geometrically integral cover of degree ten. Every rational
Hasse-failing target belongs to the image of `X_2(Q)`.

There is also a geometrically integral finite etale cover `X_H` of degree
twenty, with stabilizer (2), whose rational image contains every Hasse
failure. It selects a quadratic factor and a square root of the full
discriminant. Consequently the ambient Hasse-failing set is thin of type II.

**Proof.** The ordered-root cover over `U` has group `S_5`. Quotient by the
stabilizer `S_2 x S_3` of a two-element root subset. The geometric action on
the ten subsets is transitive, proving geometric integrality; its stabilizer
has order twelve, so the degree is ten. Equivalently its fibre consists of
the ten monic quadratic factors over an algebraic closure. A rational
quadratic factor is a rational point of that fibre. FMA3 supplies one at
every Hasse failure.

Intersect that stabilizer with the sign kernel. The new stabilizer has
order six, so the degree is twenty. On the quadratic-factor cover, adjoining
a square root of `Disc(f)` gives this connected double cover. At a Hasse
target the discriminant is a nonzero rational square by ABC1, so this cover
also has a rational point. The construction on the whole etale open, rather
than only a coefficient chart, includes the centered factors on the special
lines. No Hasse target has been lost at the boundary. \(\square\)

## Theorem ABC3: ambient rational bounds and sharp local/global orders

There are fixed positive constants such that, as `B` tends to infinity,

\[
 \boxed{\frac{B^2}{\log B}\ll A(B)\ll B^3},\qquad
 \boxed{G(B)\asymp B^3,\quad L(B)\asymp B^3}.           \tag{4}
\]

More precisely the lower limits satisfy

\[
 \liminf_{B\to\infty}\frac{A(B)\log B}{B^2}\ge C_\Phi,
 \qquad
 G(B)\ge\frac4{\zeta(3)}B^3+O(B^2).                   \tag{5}
\]

The total number of rational targets in `A^3` at this height is
`8 B^4/zeta(4)+O(B^3)`. Thus the proportion with locally soluble fibres
has sharp order `1/B`. The proportion of Hasse-failing targets is at most
`O(1/B)`; no corresponding sharp ambient Hasse order or leading constant
is asserted.

**Proof.** We use the following exact external interface: Buggenhout,
Cluckers, Salberger, Santens and Vermeulen, *Serre's question on thin sets in
projective space*, [arXiv:2506.13471v2](https://arxiv.org/pdf/2506.13471v2),
20 January 2026, Theorem 1.1, PDF p. 4. For any integral `X` and a fixed
quasi-finite dominant `Q`-morphism `g:X -> P^n` of degree at least five,
it bounds the number of image points of standard height at most `B` by
`O_g(B^n)`. Its source is not required to be projective or to map onto
every point of projective space.

Apply it to `X_2 -> U -> P^3`. The first morphism is finite etale of degree
ten, the second is an open immersion, and ABC2 verifies integrality and
containment of every Hasse target. The resulting morphism is quasi-finite,
dominant, and has degree ten, so `A(B)<<B^3`.

For rationally soluble targets, use the fixed map itself:
`A^3 --Phi--> A^3 -> P^3`. It is quasi-finite because it is etale, and
has degree five by the fixed-map construction. The same theorem gives
`G(B)<<B^3`, and hence `L(B)<<B^3` by the exact partition.

Every rational target in the plane `u=0` has a rational source point, by
FMA1. Its primitive projective representatives are `[q:0:m:n]` with
`1<=q<=B`, `|m|,|n|<=B`, and `gcd(q,m,n)=1`. Mobius inversion gives

\[
 \sum_{d\le B}\mu(d)\lfloor B/d\rfloor
                    (2\lfloor B/d\rfloor+1)^2
 =\frac4{\zeta(3)}B^3+O(B^2).                         \tag{6}
\]

Indeed the main term uses the absolutely convergent sum `sum mu(d)/d^3`;
the floor errors sum to `O(B^2 sum d^-2)`, and the completed main-term
tail costs `O(B)`. This proves the lower bound for `G`, and hence `L`.
The original-line Hasse count supplies (5) for `A` in exactly the same
height. Repeating the primitive count with three signed numerator
coordinates gives `8B^4/zeta(4)+O(B^3)`, proving the density statements.
\(\square\)

The external theorem is used with the degrees of finite field extensions
(five, ten, twenty), not the coordinate degree of `Phi` or the polynomial
degree of a defining equation for a cover. The two original arithmetic
lines have no integral targets, so rational and integral counting domains
must also remain distinct.

## Theorem ABC4: stronger integral-box bounds from the larger sieve

Let `A_Z(B),L_Z(B)` count the corresponding fibres at integer targets
`(u,b,c) in Z^3` with `max(|u|,|b|,|c|)<=B`. Then

\[
 A_{\mathbb Z}(B)\ll B^{23/10}\log B,\qquad
 L_{\mathbb Z}(B)\ll B^{5/2}\log B.                    \tag{7}
\]

These are upper bounds on integral *targets*, not on integral source
points or on the integer parameters of the original target line.

**Proof.** Apply Zywina, *Hilbert's irreducibility theorem and the larger
sieve*, [author PDF](https://pi.math.cornell.edu/~zywina/papers/Quantitative-HIT.pdf),
30 November 2010, Theorems 1.3--1.4, PDF p. 2. Theorem 1.4 says that for
a monic irreducible polynomial over `k(T_1,...,T_n)` whose normal closure
is regular over `k`, and a conjugacy-stable subset `C` of its generic
group `G`, the number of integral specializations with group contained
in `C` is `O(B^([k:Q](n-1+|C|/|G|)) log B)`, away from poles and
nonseparable fibres. Theorem 1.3 instead gives the usual large-sieve
bound `O(B^([k:Q](n-1/2)) log B)` for any proper such subset.

ABC1 verifies these hypotheses for `f=E/pi^5`, `k=Q,n=3`, and the
conjugacy-stable set of size 36 in `S_5`. Every Hasse specialization lies
in it and is away from the exceptional set. Since `36/120=3/10`, the
exponent is `3-1+3/10=23/10`.

For locally soluble targets in `U`, FMA3 says there is either a rational
root or a Hasse failure. Rational roots force the specialized group to be
contained in the conjugacy-stable proper subset of permutations with a
fixed point; Hasse groups are contained there too. Theorem 1.3 gives
`O(B^(5/2) log B)`. Outside `U`, a fixed nonzero polynomial
`pi Disc(E)` vanishes, so there are `O(B^2)` integer targets there by the
elementary polynomial zero bound. This is absorbed by (7). \(\square\)

For comparison, the rational-height part of Zywina's Theorem 1.4 alone
would give `A(B)<<B^(33/10) log B`. ABC3 improves it to `O(B^3)` using
the degree-ten cover and the more recent projective theorem. No affine
leading-form hypotheses from that theorem are presumed in (7).

## Verification and remaining questions

The exact controls in [verify_fixed_map_ambient_bounds.py](verify_fixed_map_ambient_bounds.py)
enumerate all 120 permutations, check both stabilizers, their conjugacy
union, the factor discriminant identity and the primitive plane-count
identity. They do not prove generic monodromy, the imported sieve or
projective theorems, or an asymptotic by finite extrapolation. The
[receipt](fixed_map_ambient_bounds_checks.json) pins these controls.

Run from the repository root with Python 3 and SymPy:

```sh
.venv/bin/python research/arithmetic/verify_fixed_map_ambient_bounds.py --write-receipt
```

The two primary PDFs are retained under `sources/`. The exact applied
statements were inspected. Their deeper inputs remain imported theorem
dependencies; no formal verification, independent whole-proof replay or
external review is claimed for this application. The ambient Hasse order
between (5) and (4), and leading constants for the full ambient local and
global counts, remain open.
