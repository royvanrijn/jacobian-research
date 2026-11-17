# A square-root logarithmic error for every fixed finite local sieve

On the existing elliptic slice, each fixed finite-place target count
has error `O(sqrt(log B))`. The additional input is a lower bound for
two p-adic elliptic logarithms. It makes the shortest vector of a deep
congruence kernel grow exponentially with its depth, so the lattice
errors of the local annuli can be summed without an iterated-logarithm
factor.

## 1. Scope and the primary logarithm bound

Retain `E`, `Gamma=E(Q)`, its rational torsion subgroup, the four-class
map `alpha`, the target height and `C_E` from
[EFA1--EFA2](ELLIPTIC_EXCEPTION_FINITE_LOCAL_ASYMPTOTICS.md).
Write `N_(S,a)(B)` for the full field targets of class `a` passing
the fixed finite set of places `S`, and `mu_(S,a)` for their normalized
Haar coefficient on the actual joint closure. The
[quantitative finite-place proof](ELLIPTIC_EXCEPTION_QUANTITATIVE_FINITE_LOCAL_ASYMPTOTICS.md),
EQA1--EQA3, supplies the open joint closure, rationality of this
coefficient and the explicit local annuli.

Fix an integral basis `P_1,P_2` of the free quotient of `Gamma`, only
for the proof. Thus every point is uniquely

\[
 P=mP_1+nP_2+\tau,
 \qquad m,n\in\mathbb Z,\quad \tau\in\Gamma_{\rm tors}.
                                                               \tag{1.1}
\]

The basis is not claimed to be one of the previously displayed pairs.
All constants below may depend on this fixed basis, hence on the
fixed elliptic data. Put `H(m,n)=max(3,|m|,|n|)`.

The new primary input is N. Hirata-Kohno and R. Takada,
*Linear forms in two elliptic logarithms in the p-adic case*,
Kyushu Journal of Mathematics **64** (2010), 239--260,
DOI `10.2206/kyushujm.64.239`, Corollary 2.1, printed p.243
(PDF p.5), with Assumptions 2.1--2.2 on pp.242--243:
[published primary PDF](https://www.jstage.jst.go.jp/article/kyushujm/64/2/64_2_239/_pdf),
[retained copy](sources/hirata-kohno-takada-2010-padic-elliptic-logarithms.pdf).
For fixed nonzero elliptic logarithms `u_1,u_2` with rational
exponentials, their bound implies

\[
 |\beta u_1-u_2|_p>
       \exp\{-C\max(1,h(\beta))\}
       \quad(\beta\in\mathbb Q^*,\ \beta u_1-u_2\ne0),           \tag{1.2}
\]

where `C` depends on the fixed curve, points and prime. Here `h(beta)`
is the absolute logarithmic height. The logarithms must belong to
the domain `|u|_p<p^(-lambda_p)`, with `lambda_2=3` and
`lambda_p=1/(p-1)` for odd `p`. These domain conditions will be
ensured by fixed multiples. No complex multiplication hypothesis
is present in the cited result.

**Lemma ESQ1 (formal depth is logarithmic in the coefficients).**
For each fixed prime `p`, a nontorsion point (1.1) of exact formal
depth `k>=1` satisfies

\[
 k\le C_{E,p}\log H(m,n).                           \tag{1.3}
\]

**Proof.** Choose a fixed integer `M` killing `Gamma_tors` and
putting both `MP_i` sufficiently deep in the formal group. Choose
it also with large enough `p`-adic valuation that `MP` is in the
domain above for every `P in E_1(p)`. Then

\[
 u_i=\log_E(MP_i),\qquad
 L=mu_1+nu_2=\log_E(MP).
\]

Both `u_i` are nonzero. The domain and rationality assumptions of
(1.2) hold because `exp_E(u_i)=MP_i in E(Q)`. In the notation of
the primary text, its two algebraic functions are `1/X(MP_i)`
and `-X(MP_i)/y(MP_i)`, so they are rational as required.
Its parameters other than `h(beta)` are fixed.

The form `L` is nonzero: its exponential is `MP`, and vanishing
would give the rational relation `mMP_1+nMP_2=O`, contradicting
the free-basis property. If both coefficients are nonzero, apply
(1.2) with `beta=-m/n`. Since `h(beta)<=log H(m,n)` and
`|n|_p>=H(m,n)^(-1)`, this gives

\[
 |L|_p=|n|_p|\beta u_1-u_2|_p
                    >H(m,n)^{-(C+1)}.              \tag{1.4}
\]

If one coefficient is zero, the same type of bound follows from
`|m u_1|_p>=|u_1|_p/H(m,n)`, or its analogue for `u_2`.

On this small domain the logarithm preserves the absolute value
of `t=-X/y`: this follows from the exponential description and
the identities `|phi(u)|_p=|u|_p`, `|phi'(u)|_p=1` on printed
p.242 of the same text. Formal multiplication in the integral
short model gives

\[
 v_p(L)=v_p(t(MP))=k+v_p(M).
\]

The depth identity is also the short-model filtration calculation
used in EQA1 from Milne, Chapter II, Theorem 4.1 and Lemma 5.9.
Taking valuations in (1.4), and enlarging a fixed constant for the
one-coefficient cases, proves (1.3). \(\square\)

In particular, write `X(P)=A(P)/C(P)^2` in reduced form with
`C(P)>0`. For every fixed prime `p`, admissible points with actual
target height at most `B` satisfy

\[
 v_p(C(P))\ll_{E,p}1+\log\log B\qquad(B\ge e^e).    \tag{1.5}
\]

Indeed a positive denominator valuation is exactly the formal depth,
and the canonical-height comparison bounds the free coefficient norm
by `O_E(sqrt(1+log B))`.

## 2. A lattice bound that improves with congruence depth

Use the canonical norm `||v||_q=sqrt(q_x(v))` on the fixed free
rank-two lattice. For a finite-index subgroup `Lambda subset Gamma`,
let `L_Lambda` be its projected free lattice, and let
`lambda_1(L_Lambda)` be the shortest nonzero length in that norm.

**Lemma ESQ2 (short vectors and uniform coset counting).** Fix a
nonempty finite set `T` of primes. If

\[
 \Lambda\subseteq\bigcap_{p\in T}E_{n_p}(p),\qquad n_p\ge1,
 \qquad K=\max_{p\in T}n_p,
\]

then there are constants `c,eta>0`, depending only on the fixed data
and `T`, such that

\[
 \lambda_1(L_\Lambda)\ge c\exp(\eta K).             \tag{2.1}
\]

For every coset `C` of any finite-index subgroup `Lambda` of `Gamma`,
its actual target-height count satisfies, uniformly in that subgroup
and coset,

\[
 \#\{P\in C:H(y(P))\le B\}
 =\frac{2C_E}{[\Gamma:\Lambda]}\log B
    +O_E\left(\frac{\sqrt{\log B}}{
                       \lambda_1(L_\Lambda)}+1\right).         \tag{2.2}
\]

**Proof of (2.1).** For each nonzero vector `v` of `L_Lambda`,
choose its actual torsion lift `P in Lambda`. Lemma ESQ1 applies
to this nontorsion point at every prime of `T`. Its coefficient
height is at most `3+C_E'||v||_q`, by equivalence of the fixed
norms. Taking the largest depth yields

\[
 K\le C_{E,T}\log(3+C_E'\|v\|_q).
\]

Inverting this inequality gives (2.1) for all large `K`. For bounded
`K`, the fixed ambient free lattice has a positive shortest length;
decreasing `c` proves it for every `K`. The torsion lift does not
change the canonical norm.

**Proof of (2.2).** Apply the fixed linear transformation taking the
canonical ellipse to a Euclidean disk. In the transformed lattice
choose a shortest primitive vector `u` of length `lambda_1` and a
complementary basis vector `v`. Subtract a multiple of `u` so that
the parallel component of `v` has absolute value at most
`lambda_1/2`. Minimality gives `|v|>=lambda_1`, so its perpendicular
row spacing `h` satisfies `h>=sqrt(3)lambda_1/2`.

For a disk of radius `R`, rows have spacing `h` and points on each
row have spacing `lambda_1`. Rounding the count on each row costs
`O(R/h+1)`. The compactly supported chord-length function has total
variation `4R`, so its shifted row sum differs from its integral
divided by `h` by `O(R)`. After dividing by `lambda_1`, the total
discrepancy is

\[
 O\left(R/\lambda_1+R/h+1\right)
                         =O(R/\lambda_1+1),         \tag{2.3}
\]

uniformly in the lattice, its skew and the translated coset. This
also covers `lambda_1>R`.

Each projected free point has one or two torsion lifts in a given
coset, a constant equal to `|Lambda intersect Gamma_tors|`.
The index formula gives the main fraction in (2.2), with that
factor at most two in its error. The existing bounded comparison

\[
 \log H(y(P))=3q_x(P)+O_E(1)
\]

brackets the actual height set by ellipses whose squared radii differ
by `O_E(1)`. Their area main terms differ by `O_E(1)` uniformly in
the index; both have error (2.3) with `R` of order `sqrt(log B)`.
This proves (2.2). Finitely many global excluded points can be
removed with a further uniform `O_E(1)`. \(\square\)

The branch centers in the local sieves need not be rational. Their
rectangles pull back to cosets, but the corresponding subgroup is
always a homogeneous formal-congruence kernel. Thus (2.1) uses only
the two logarithms of the fixed basis points. No lower bound for
an inhomogeneous form involving a third branch-center logarithm
is needed.

## 3. Summing the annuli

**Theorem ESQ3 (square-root logarithmic finite-place error).** For
every fixed finite set of places `S` and every `a in A`,

\[
 \boxed{N_{S,a}(B)=\mu_{S,a}C_E\log B
                              +O_{E,S}(\sqrt{\log B}).}        \tag{3.1}
\]

The coefficient is the rational number of EQA2. It is positive for
the certified class at every fixed finite `S`; the theorem permits
zero for other classes. The error constant is not asserted to be
uniform as `S` grows.

**Proof.** Omit the places with automatic local points, leaving
the active set `T` of EQA3, and put `r=|T|`, `L=log B`, `R=sqrt(L)`.
If `r=0`, the fixed-class coset count already gives (3.1).
Suppose `r>=1`.

Use the inner and outer local approximations from EQA3, with
`N_p=ceil(log_p R)+O_{E,p}(1)`. Their local measure gap is
`O_{E,p}(p^(-N_p))`. Their balls can be chosen with at most
`O_{E,p}(1)` pieces at each subgroup depth: the local descriptions
have finitely many centers, periodic depth conditions modulo three
and a fixed angular refinement at three. Finitely many initial balls
are kept at fixed depths, without refinement to the deepest level.
Each approximation therefore has `O_{E,p}(N_p+1)` pieces as before.

Fix the class `a` and intersect a product rectangle with the actual
joint closure. A nonempty intersection is a coset of its intersection
with the homogeneous product of the defining formal subgroups.
Density of the rational image makes its pullback a single nonempty
coset of a finite-index subgroup `Lambda` of `Gamma`. If the
rectangle depths are `(n_p)`, then `Lambda subset E_(n_p)(p)`
for every `p in T`. Lemma ESQ2 bounds its lattice error by

\[
 O_{E,S}\left(R\exp(-\eta\max_p n_p)+1\right).      \tag{3.2}
\]

There are only `O_{E,S}(1)` rectangles with any fixed depth vector.
The number of vectors with maximum `k` is `O_r((k+1)^(r-1))`.
Hence the sum of the first terms in (3.2) is at most

\[
 O_{E,S}\left(R\sum_{k\ge1}(k+1)^{r-1}
                                    e^{-\eta k}\right)=O_{E,S}(R).
                                                               \tag{3.3}
\]

The sum of the second terms is `O_{E,S}((1+log R)^r)`, which is
`O_{E,S}(R)` for fixed `r`. Thus both inner and outer point counts
have total lattice error `O_{E,S}(R)`.

Their difference in measure is at most the sum of the local gaps,
using the marginal Haar measures and the union bound. Its main-term
contribution is

\[
 O_{E,S}\left(L\sum_{p\in T}p^{-N_p}\right)=O_{E,S}(R).
\]

The two approximations therefore squeeze the point count to its
Haar main term with error `O_{E,S}(R)`. Remove the finite global
inadmissible and reducible-cubic exceptions. The local conditions
and the global class are invariant under negation. Each target has
exactly the two preimages `P,-P` away from finitely many exceptions,
so division by two gives (3.1). \(\square\)

## 4. The explicit three-place constants with this error

**Corollary ESQ4 (three-place counts and Hasse upper bounds).** Let
`S={2,3,17}`. The coefficients calculated in
[J2173--J2174](ELLIPTIC_EXCEPTION_JOINT_2_3_17_SIEVE.md) satisfy

\[
 N_{S,a}(B)=c_aC_E\log B+O_E(\sqrt{\log B}),
\]

with the following exact values, normalized on the whole joint group:

| Class `a`, with `q_0=2709173` | Quadratic field | `c_a` |
| --- | --- | --- |
| `1` | `Q(sqrt(-3))` | `205/111748` |
| `2` | `Q(sqrt(-6))` | `205/2456` |
| `5q_0` | `Q(sqrt(-40637595))` | `1/84` |
| `10q_0` | `K_0=Q(sqrt(-81275190))` | `1/8` |

In particular the total full field-target count passing these three
places is

\[
 M_{2,3,17}(B)=\frac{5321}{23946}C_E\log B
                                      +O_E(\sqrt{\log B}).      \tag{4.1}
\]

The all-places Hasse-target counts therefore satisfy

\[
 \boxed{A_E(B)\le\frac{5321}{23946}C_E\log B
                                     +O_E(\sqrt{\log B}),}
                                                               \tag{4.2}
\]
\[
 \boxed{A_{K_0}(B)\le\frac18 C_E\log B
                                     +O_E(\sqrt{\log B}).}      \tag{4.3}
\]

**Proof.** Apply ESQ3 to the exact four coefficients of J2174 and
sum them; their sum is the coefficient of J2173. Every all-places
Hasse target passes these three places, giving (4.2)--(4.3).
\(\square\)

The same deduction gives the square-root logarithmic error for every
previously computed fixed finite-place coefficient. Equations
(4.2)--(4.3) are upper bounds for the all-places counts; the leading
asymptotics here count only the specified finite local tests.

## 5. Scope of the conclusion

The primary logarithm bound, not a finite experiment, supplies the
exponential estimate (2.1). The reduced-basis argument and weighted
annulus sum are the written deductions above. The other imports are
the previously retained local tests, canonical-height comparison,
finite global exceptions and target quotient. No new numerical
checker is needed for this deduction, and no formal verification
or external review is claimed.

An explicit finite quotient is still required to evaluate a particular
Haar coefficient, as in EQA2. Existence of a basis in the proof does
not certify a displayed basis or provide numerical values for the
error constants. Most importantly, (3.1) is for each fixed finite
`S`. It gives no uniform error for a growing prime set and does not
resolve the all-primes measure or tail condition of
[EFA3](ELLIPTIC_EXCEPTION_FINITE_LOCAL_ASYMPTOTICS.md).
