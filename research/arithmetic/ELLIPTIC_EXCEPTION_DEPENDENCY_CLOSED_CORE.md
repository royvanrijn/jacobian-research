# A finite dependency core and safe-prime local families

For the two quadratic classes whose positivity at every finite set
is unresolved, the prime `q_0=2709173` is already automatic. The
remaining bad-prime conditions use only the primary factors
`2,3,5,7,1889`. Exposing those five factors first removes the large
unspecified cutoff from the triangular extension argument.

Every finite set consisting of the four core tests and primes with
no supported fatal branch is feasible in all four classes. This is
a restricted finite-set theorem; it does not prove feasibility for
arbitrary finite sets or produce an everywhere locally soluble point.

## 1. The core and its exact dependencies

Retain
\[
 E:y^2=x^3+ax^2+bx,\quad a=-10656,\quad b=27091730,
 \qquad F(x)=128x-25865,
\]
\[
 q_0=2709173,\qquad d=51730,\qquad
 A=\{1,2,5q_0,10q_0\}.
\]
Write `Gamma=E(Q)`. Its rank is two and its rational torsion is
`C_2`; put `T=(0,0)`. The four-class map is `alpha`, with quadratic
field `Q(sqrt(-3 alpha))`. Define the test core and its primary
dependency set by
\[
 C_0=\{2,3,5,1889\},\qquad
 D=\{2,3,5,7,1889\}.                              \tag{1.1}
\]
A test depends only on `D` if its pullback to
\[
 \widehat\Gamma=C_2\times\prod_\ell\mathbb Z_\ell^2
\]
is measurable in the torsion and factors indexed by `D`. This
formulation retains every finite correlation among the tests.

**Theorem EDC1 (explicit core and all-four feasibility).**
The core tests and the class map depend only on `D`. For every
`alpha in A`, their common full-field target coefficient satisfies
\[
                         \mu_{C_0,\alpha}>0.       \tag{1.2}
\]
For the two classes `alpha=2,5q_0`, every local test outside the
core is either automatic or a test at a good prime `p>=11`,
`p!=1889,q_0`, with `p=2 mod3`. In particular the bad prime
`q_0` imposes no test in those two classes.

**Proof of the dependency assertion.** The curve discriminant is
\[
 \Delta_E=16b^2(a^2-4b)=2^9 5^2 7^3 1889 q_0^2.     \tag{1.3}
\]
The actual rational two-adic closure is `<T>+E_4(2)`, by ETF1 in
[the two/five-adic classification](ELLIPTIC_EXCEPTION_TWO_FIVE_ADIC_CLASSIFICATION.md).
It is pro-two, so its map uses only the two-primary data. The
actual three-adic closure is `E(Q_3)[2]+E_2(3)`, by ETC3a in
[the complete three-adic classification](ELLIPTIC_EXCEPTION_COMPLETE_THREE_ADIC_CLASSIFICATION.md).
Its finite part is a two-group and its formal part is pro-three.
The global class map also factors through the two-primary data.

At five the reduced cubic has its unique node at `(0,0)`.
Here `v_5(b)=1` and `a` is a unit. Let `E_0(5)` be the subgroup
of points reducing nonsingularly and `E_1(5)` its formal kernel.
Every point outside `E_0(5)`, apart from `T` itself, has
`k=v_5(x)>=1`. Addition by `T` gives the exact formula
\[
          (x,y)+T=(b/x,-by/x^2).                  \tag{1.4}
\]
If `k=1`, the translated x-coordinate is a unit, so the translated
point reduces nonsingularly. If `k>=2`, the term `bx` has strictly
smaller valuation than `ax^2` and `x^3`. Thus `v_5(y^2)=k+1`,
which makes `k` odd and at least three; the translated coordinate
has valuation `1-k<0`, so it is formal. Hence all the singularly
reducing points lie in `T+E_0(5)`. Since `T` itself reduces to the
node, the quotient `E(Q_5)/E_0(5)` has exactly two elements.
The tangent coefficient at the node is `a=4 mod5`, a square.
The nonsingular nodal reduction group therefore has order four:
\[
                 |E(\mathbb Q_5)/E_1(5)|=8.       \tag{1.5}
\]
Together with its formal pro-five subgroup, the local map uses
only the two- and five-primary factors.

At `p=1889`, write
\[
 x^3+ax^2+bx=x\bigl((x-5328)^2-1295854\bigr),
                 \qquad1295854=1889\cdot686.      \tag{1.6}
\]
The unique reduced node is `(5328,0)`, and `5328` is a unit.
A point reducing to this node would have `x` a unit and
`x-5328` divisible by `1889`. The right side of (1.6) would then
have valuation one, impossible for a square. Thus every local
point reduces nonsingularly: `E(Q_1889)=E_0(1889)`.
After translating the node, its tangent coefficient is `5328`,
a nonsquare modulo `1889`. Its nonsingular group has order
`1889+1`, giving
\[
 |E(\mathbb Q_{1889})/E_1(1889)|
                   =1890=2\cdot3^3\cdot5\cdot7.  \tag{1.7}
\]
Its formal part is pro-1889. This proves the dependency set.
These arguments use nodal reduction and elementary valuations,
not a Tate-algorithm output.

The local interface is Milne, *Elliptic Curves* (2006), Chapter II,
section 3, printed pp.57--59 (split and nonsplit nodal groups), and
Theorem 4.1, pp.62--63 (nonsingular reduction quotient and formal
filtration): [primary book](https://www.jmilne.org/math/Books/ectext6.pdf),
[retained copy](sources/milne-elliptic-curves-2006.pdf).
The integral models at five and 1889 are minimal: their
discriminant valuations are two and one, both less than twelve.
The translation to the short model is integral at these primes.

**Proof of feasibility and erasure.** The concrete points
\[
 P=(373248/49,59459616/343),\quad
 Q=(3442688/529,304670016/12167),\qquad U=2P+3Q+T
\]
are the exact full-field witnesses of
[EFO1](ELLIPTIC_EXCEPTION_FINITE_OBSTRUCTION_CERTIFICATES.md).
`Q` has class two and `U` has class `5q_0`. They pass `2,3,5`.
In addition,
\[
                    v_{1889}(F(Q))=v_{1889}(F(U))=0. \tag{1.8}
\]
ELV1 applies at 1889 without a good-reduction assumption, so both
pass that test. Local solubility is open at these admissible
full-field points, so EQA2a gives positive coefficients. The
existing all-finite-set results ECO1 for class one and the certified
seed for class `10q_0` prove (1.2) for the other two classes.

By EPE2 in
[the eligibility classification](ELLIPTIC_EXCEPTION_PRIME_ELIGIBILITY.md),
`q_0` is automatic for class two because `q_0=5 mod8`, and for
class `5q_0` because that class is ramified there. Primes one modulo
three and the real place are automatic. The other bad prime seven
is one modulo three. Equation (1.3) gives the last assertion.
\(\square\)

## 2. Fresh formal factors outside the core

For `p notin D`, let `F^D_<p` be the sigma-algebra generated by
`C_2` and the primary factors indexed by `D` or less than `p`.
For classes `2,5q_0`, omit the automatic tests before forming a
smaller-prime prefix.

**Theorem EDC2 (core-adjusted triangular extension).** Fix
`alpha=2` or `5q_0`. At every remaining active prime `p>=11`,
`p!=1889,q_0`, the core event, class condition, remaining tests at
smaller primes, and reduction modulo `p` are `F^D_<p`-measurable.
The formal `p`-coordinate is independent of this sigma-algebra and
Haar-uniform in `E_(e_p)(p)`. The ETP4 residue-hazard formula holds
with this filtration and no numerical cutoff `p_*`.

**Proof.** A good prime `q>=7` has even `N_q=#E(F_q)<2q`,
by rational two-torsion and Hasse's bound. Hence every prime
divisor of `N_q` is smaller than `q`, as in ETP1. A good local
map at `q` therefore uses only its own primary factor and lower
ones. EDC1 exposes all core dependencies initially; the remaining
bad-prime test `q_0` is erased, not treated as an independent map.
At the new prime, reduction uses only divisors of `N_p`, all less
than `p`. Its formal map uses only the independent factor `Z_p^2`,
whose image is `E_(e_p)(p)` by EHD3.

In particular, when `p<1889`, the core condition at 1889 is not
lower-primary measurable in the original ETP filtration. It is
measurable here, since the 1889 factor was exposed initially and
the finite quotient (1.7) has no prime factor `p>=11`. Equivalently
one can condition on all non-`p` primary factors; the `p` factor
remains independent. This proves the asserted independence without
assuming independence of any reduction conditions.

At these active primes the branch and cubic constants are units,
and the zeros of `F` have unit y-coordinates when rational over
`Q_p`. ETP4's proof uses exactly this local geometry, the explicit
ELV1 valuation criterion and the fresh formal factor. Repeating it
gives its residue table here. In particular the only wholly
insoluble supported residues are the finite zero residues whose
offset depth satisfies
\[
                         k_p<e_p,\qquad3\nmid k_p. \tag{2.1}
\]
There are at most two. The pole residue always has a positive
soluble fraction. \(\square\)

## 3. Finite safe-prime families and the first possible failure

For either class `alpha=2,5q_0`, call a prime outside `C_0`
**safe for that class** if it is automatic for the class, or if it
is a remaining good active prime for which none of the zero
residues supported by that class satisfies (2.1). Here supported
means in the image under reduction of the actual class coset in
`widehat Gamma`. Let `G_alpha` be this set and put
`G_safe=G_2 intersection G_(5q_0)`. These sets use the actual
rational image and formal depth, not a displayed subgroup presumed
saturated.
Sufficient conditions at an active prime include `e_p=1`, absence
of rational zero residues, absence of those residues from the
rational reduction image, `k_p>=e_p`, or `3|k_p`.

**Theorem EDC3 (positive safe families and fatal-prime reduction).**
For every class `alpha in A` and every finite set of places `S`
whose finite primes lie in `C_0 union G_safe`,
\[
                            \mu_{S,\alpha}>0.       \tag{3.1}
\]
Its coefficient is rational, and the full-field target count is
\[
 N_{S,\alpha}(B)=\mu_{S,\alpha}C_E\log B
                              +O_{E,S}(\sqrt{\log B}). \tag{3.2}
\]
For either class `2,5q_0`, if any finite collection is infeasible,
add `C_0`, erase automatic primes and process the others increasingly.
Its first killing prime must have supported finite zero residues
satisfying (2.1). In particular
\[
 e_p>=2,\quad k_p<e_p,\quad3\nmid k_p,\quad
                    p\in\mathcal P_{\alpha,d}.     \tag{3.3}
\]
The entire positive prefix support must lie on those fatal residues.
This also holds for the largest new prime of a minimal obstruction
relative to the fixed feasible core.

**Proof.** Classes one and `10q_0` are already feasible at every
finite set. For the other two classes, enlarge `S` by `C_0` and
start with EDC1's positive event. Process the nonautomatic primes
outside the core increasingly. At each safe prime EDC2 supplies
an independent formal coordinate. Every supported reduction has
strictly positive conditional survival. There are finitely many
reductions, so integrating over any positive prefix leaves positive
mass. Finite induction proves (3.1). Rationality and (3.2) are
EQA2 and ESQ3, for each fixed finite set.

At a first killing prime, the same argument forces the prefix
onto reductions with zero conditional survival, exactly (2.1).
Since `k_p>=1`, this implies `e_p>=2`. EPE2 also makes a numerator
obstruction impossible outside `P_(alpha,d)`. For a minimal
obstruction relative to the core, deleting its largest new prime
leaves a positive core-prefix, proving the last assertion.
\(\square\)

For one fixed class, the same proof requires only that the added
primes belong to its own `G_alpha`.

The positive-density prime set of
[EKX1--EKX3](ELLIPTIC_EXCEPTION_KUMMER_BRANCH_EXCLUSION.md)
has no supported finite zero residues. Its primes outside the core
are therefore safe by the same local argument. It gives an
unconditional positive-density collection that may be adjoined in
arbitrary finite subsets, simultaneously in all four classes.
No distribution hypothesis about elliptic Wieferich primes is used.

Safe conditions cannot be erased from arbitrary sets containing
later fatal primes: their reductions can influence later support.
The result does not show that fatal primes exist, that any one kills
a prefix, or that all finite sets in the two unresolved classes are
feasible. It isolates their possible first failure to (3.3), without
an unspecified initial numerical cutoff.

## 4. Fixed safe pools of exact positive density

Let `P_(alpha,d)` be the numerator-eligible prime set of EPE2,
and fix the Chebotarev prime set `P_Kum` of EKX1, of density
`delta_Kum>0`. Its primes belong to both eligible sets. Define
\[
 T_\alpha=C_0\ \cup\ (\{\text{primes}\}\setminus
                  \mathcal P_{\alpha,d})\ \cup\mathcal P_{\rm Kum}
       \quad(\alpha=2,5q_0),                              \tag{4.1}
\]
\[
 T=C_0\ \cup\bigl(\{\text{primes}\}\setminus
        (\mathcal P_{2,d}\cup\mathcal P_{5q_0,d})\bigr)
                                      \ \cup\mathcal P_{\rm Kum}. \tag{4.2}
\]

**Corollary EDC4 (universal finite safe pools).** These fixed
sets have exact natural prime densities
\[
 \boxed{\operatorname{dens}(T_\alpha)=\frac78+\delta_{\rm Kum},
 \qquad\operatorname{dens}(T)=\frac{13}{16}+\delta_{\rm Kum}.} \tag{4.3}
\]
Every finite subset of `T_alpha` has positive coefficient in class
`alpha`. Every finite subset of `T` has positive coefficient in
all four classes, with the rational positive coefficient and sharp
fixed-set asymptotic (3.2).

**Proof.** Outside `P_(alpha,d)`, EPE2 excludes a numerator
obstruction in that class. An automatic prime needs no test; at a
remaining good prime the only possible obstruction is the pole,
whose soluble fraction is always positive. EDC2 therefore extends
every positive earlier prefix at such a prime. On `P_Kum`, EKX1
excludes both finite zero residues from the whole rational image,
so the same assertion holds. Thus `T_alpha` lies in the core
together with `G_alpha`, and `T` lies in the core together with
`G_safe`. EDC3 gives positivity.

The four quadratic characters of `-3,2,5q_0,d` are independent,
as proved in EPE3. Each numerator-eligible set has density `1/8`;
their intersection prescribes all four signs and has density
`1/16`. Their union has density `3/16`. Since `P_Kum` is a
known-density subset of their intersection, its union with either
complement in (4.1)--(4.2) is disjoint up to finitely many primes.
This proves (4.3). \(\square\)

The pools are fixed once the auxiliary EKX division field is fixed.
Their exact densities involve that field's positive constant, not a
computed numerical value. No coefficient lower bound uniform in the
finite subset is asserted; positivity of every finite stage does not
imply positive all-primes limiting measure or Hasse infinitude.

[The bounded depth certificate](ELLIPTIC_EXCEPTION_BOUNDED_FATAL_PRIME_CERTIFICATE.md)
shows that every prime through `100000` is safe for both unresolved
classes: every potentially fatal eligible prime there has depth one.
Consequently every finite local set supported below that inclusive
bound is feasible in all four classes, and may be enlarged by any
finite collection of the appropriate safe primes. A first killing
prime in the core-adjusted increasing order must exceed the bound.
No assertion about larger fatal depths is inferred from this check.

## 5. Exact replay and scope

Run from the repository root:

```sh
.venv/bin/python research/arithmetic/verify_elliptic_exception_dependency_closed_core.py
```

The checker passes 45 named exact controls and replays the existing
37-control Q/U baseline certificate. It checks
the nodal identities, exact discriminant and quotient orders,
nonsplit tangent at 1889, dependency sets and `F`-unit controls.
Its complete reduced-curve scans are at the two primes five and
1889; the density arithmetic enumerates all sixteen sign patterns.
Its [receipt](elliptic_exception_dependency_closed_core_checks.json)
pins the source, checker and baseline inputs. These exact finite
controls support the written local proofs; they do not establish
the unrestricted induction by experimentation, compute a full
Mordell--Weil basis or exhibit a finite obstruction. No formal
verification, external review or independent replay is claimed.
