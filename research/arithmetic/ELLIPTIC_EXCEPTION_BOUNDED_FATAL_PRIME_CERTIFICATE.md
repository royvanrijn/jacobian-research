# No fatal-prime depth through 100000 on the elliptic slice

Every prime that could create a first finite local obstruction in
either unresolved quadratic class has formal depth one through the
inclusive bound `100000`. One fixed rational point supplies every
depth-one witness. Consequently every finite set of local conditions
supported on primes at most this bound is feasible in all four
classes. This is a finite-range result; it makes no assertion about
larger simultaneous elliptic Wieferich primes or everywhere local
solubility.

## 1. Complete prime domain and the depth-one certificate

Use the short model and known rational points
\[
 E:y^2=X^3+AX+B,\quad A=-10758382,\quad B=6600759744,
 \qquad X=x-3552,
\]
\[
 P=(199200/49,59459616/343),\qquad
 Q=(1563680/529,304670016/12167).
\]
The original model is
`y^2=x^3-10656x^2+27091730x`, with
`q_0=2709173`, branch class `d=51730` and four global classes
`1,2,5q_0,10q_0`. For a good prime let `N_p=#E(F_p)`.
The actual rational formal depth is
\[
 e_p=\min_{W\in E(\mathbb Q),\ W\text{ nontorsion}}
                  v_p(t([N_p]W)),\qquad t=-X/y,
                                                               \tag{1.1}
\]
as proved in EHD3. All good primes `p>=7` have `p not dividing N_p`.

Let `L=100000` and define the complete selected prime set
\[
 \begin{split}
 \mathcal R_L=\{p\le L:\;&p\text{ prime},\ p\ge11,
       \ p\equiv2\pmod3,\ p\notin\{1889,q_0\},\\
       &(d/p)=1,\quad (2/p)=1\text{ or }(5q_0/p)=1\}.
 \end{split}                                                   \tag{1.2}
\]
Thus these are exactly the noncore good active primes in the union
of the numerator-eligible sets of the two unresolved classes.

**Theorem EBF1 (complete depth-one certificate).** The set
`R_L` has exactly `1722` primes. At every one of them,
\[
                         v_p(t([N_p]Q))=1,
 \qquad\boxed{e_p=1\quad(p\in\mathcal R_L).}       \tag{1.3}
\]
No generator-basis or saturation hypothesis is required.

**Certificate and proof.** The
[complete exact certificate](elliptic_exception_bounded_fatal_prime_certificate.json)
contains every selected prime, its exact curve order and the two
coefficients
\[
 u_P(p)=\frac{t([N_p]P)}p\bmod p,\qquad
 u_Q(p)=\frac{t([N_p]Q)}p\bmod p.                  \tag{1.4}
\]
The [standard-library verifier](verify_elliptic_exception_bounded_fatal_primes.py)
regenerates every prime through `L` by an Eratosthenes sieve and
applies the exact Euler-character tests in (1.2). It obtains the
following complete domain counts:

| Domain | Number of primes |
| --- | ---: |
| All primes at most `100000` | 9592 |
| Good active primes at least 11 outside the core | 4804 |
| Those with `(d/p)=1` | 2336 |
| Those satisfying the additional class-eligibility union | 1722 |

It recomputes each curve order by scanning every affine
x-coordinate. Explicitly, for the original integral cubic `f`,
\[
 N_p=p+1+\sum_{x\in\mathbb F_p}\chi_p(f(x)).        \tag{1.5}
\]
The character table is constructed from all nonzero square classes,
so this step makes no point-counting oracle assumption. It then
computes both multiples in `E(Q_p)/E_2(p)` using the exact chart
arithmetic proved in section 2. Every result is formal, and the
retained coefficient `u_Q(p)` is nonzero in every row. The equality
in (1.3) follows. Since a single rational witness bounds the
minimum (1.1) above by one, and all such minima are at least one,
this proves `e_p=1` without using the subgroup index of `P,Q`.
\(\square\)

The only selected prime for which `u_P(p)=0` is `401`:
\[
 (p,N_p,u_P(p),u_Q(p))=(401,408,0,288).             \tag{1.6}
\]
This is a higher-depth event for the point `P`, not for the actual
rational group. The point `Q` still proves `e_401=1`.

## 2. Exact group arithmetic modulo the second formal subgroup

For odd good `p>=11`, represent a class in `E(Q_p)/E_2(p)` by
one of two charts:

- `a(X,y)`: finite affine reduction, with both coordinates modulo
  `p^2`;
- `f(t)`: reduction at infinity, with `t in pZ/p^2Z`.

In the formal chart `t=0` represents `E_2(p)`. These are quotient
classes, not assertions that the underlying local point is the
identity. Rational affine input is reduced directly when its
denominators are units; input with a pole is represented by
`t=-X/y`. In particular the rational denominator of `Q` at `23`
is handled by the formal chart.

Write `f(X)=X^3+AX+B` and use all formulas below modulo `p^2`.
The formal and mixed laws are
\[
 f(t_1)+f(t_2)=f(t_1+t_2),
\]
\[
 a(X,y)+f(t)=a(X+2yt,\ y+f'(X)t).                \tag{2.1}
\]
For two affine points `a(X,y),a(U,V)`, use the usual chord
formula when `U-X` is a unit, with
\[
 \lambda=(V-y)/(U-X),\quad
 X'=\lambda^2-X-U,\quad y'=-y+\lambda(X-X').     \tag{2.2}
\]
If `U=X mod p` and `y+V` is a unit, use the same formula with
the unit-denominator slope
\[
                 \lambda=(X^2+XU+U^2+A)/(y+V).   \tag{2.3}
\]
If `U=X mod p`, `y+V=0 mod p` and `V-y` is a unit, the sum is
formal and
\[
                         t=(U-X)/(V-y).          \tag{2.4}
\]
The last case is reduction to the same two-torsion point, where
both `y,V` are divisible by `p`. Good reduction makes `f'(X)` a
unit, and the sum is formal with
\[
                         t=(y+V)/f'(X).          \tag{2.5}
\]
Thus every division used by the algorithm is by a unit.

For completeness, (2.1) is the linear part of the integral formal
group law and translation. Its tangent vector is
`(dX,dy)=(2y,f'(X)) dt`; quadratic terms vanish modulo `p^2`
since `t in pZ_p`. It also follows by inserting the short-model
expansions `X(t)=t^(-2)+O(t^2)`, `y(t)=-t^(-3)+O(t)` into
the chord formula. Formula (2.3) is the rationalized secant slope,
using `V^2-y^2=(U-X)(X^2+XU+U^2+A)`; it also covers doubling.
In (2.4), the inverse slope is in `pZ_p` and the chord formula
gives `t=lambda^(-1)+O(lambda^(-3))`, proving that formula.

For (2.5), let `r` be the unique local lift of the reduced
two-torsion root. Hensel applies because `f'(r)` is a unit.
Near `(r,0)`,
\[
 X-r=y^2/f'(r)+O(y^4),\qquad
 t(a(X,y)+(r,0))=y/f'(r)+O(y^3).
\]
Translate both inputs by `(r,0)` and add in the formal chart.
Their `X` coordinates agree with `r` modulo `p^2`, giving (2.5).
These derivations also cover exact opposite points and exact
two-torsion points by the regularity of the group law.

The local filtration and formal parameter are those of Milne,
*Elliptic Curves* (2006), Chapter II, Theorem 4.1, printed
pp.62--63, and the short-model addition expansion of Lemma 5.9,
p.67: [primary book](https://www.jmilne.org/math/Books/ectext6.pdf),
[retained copy](sources/milne-elliptic-curves-2006.pdf).
The explicit formulas above prove the arithmetic used by the
checker; no projective group-law formula with a nonunit inverse
or uncontrolled cancellation is used.

As an implementation control separate from the selected-prime
calculation, all `7623` additions on the ordered pairs of
`O, +/-P,...,+/-8P, +/-Q,...,+/-8Q` at the seven primes
`11,17,23,29,41,53,101` are compared with exact rational chord
arithmetic before reduction. All six chart branches are exercised.
This finite cross-check supplements the written formula proof.

## 3. Feasibility of every finite set below the bound

**Corollary EBF2 (complete bounded finite-local feasibility).**
For every class `alpha in {1,2,5q_0,10q_0}` and every finite set
of places `S` all of whose finite primes are at most `100000`,
\[
                           \mu_{S,\alpha}>0.       \tag{3.1}
\]
Its coefficient is rational and
\[
 N_{S,\alpha}(B)=\mu_{S,\alpha}C_E\log B
                                  +O_{E,S}(\sqrt{\log B}). \tag{3.2}
\]
The same is true after adjoining any finite collection of the
class-specific safe primes of EDC3. When noncore tests are processed
in increasing prime order after the fixed core, any first killing
prime in either unresolved class must exceed `100000`.

**Proof.** The classes one and `10q_0` are already feasible for
every finite set. For either other class start with the positive
test core `C_0={2,3,5,1889}` of EDC1, which lies below the bound.
Every other bad or automatic test is erased by EPE2. At a
remaining good active prime, a finite zero obstruction requires
membership in `P_(alpha,d)`. Such a prime at most `100000` lies
in (1.2), so EBF1 proves `e_p=1`. The fatal configuration
`k_p<e_p`, `3 not dividing k_p` is impossible. Outside that
eligible set, no finite zero obstruction occurs in the class.
Thus every noncore prime through the bound is safe for each
unresolved class. EDC2--EDC3 give (3.1), including any finitely
many further safe primes. EQA2 and ESQ3 supply rationality and
(3.2). The same argument excludes a killing step through the
bound in any core-adjusted increasing prefix. \(\square\)

Here the permitted safe primes in the extension are understood for
the class being counted; for a simultaneous all-four assertion use
their common intersection. No coefficient lower bound uniform in
the prime set is claimed. This does not exhibit a single rational
target passing every prime, or prove positivity for arbitrary
finite sets involving larger potential fatal primes.

## 4. Replay and assurance

From the repository root run

```sh
.venv/bin/python research/arithmetic/verify_elliptic_exception_bounded_fatal_primes.py
```

The [receipt](elliptic_exception_bounded_fatal_prime_checks.json)
records complete domain counts, every point-count and modular
multiple row covered, exact-rational control counts and source/input
hashes. The retained replay uses only Python standard-library
integers and fractions. The finite certificate and the written
induction establish exactly the bounded statements above; no claim
about unbounded formal depths, a saturated displayed basis, new
all-places Hasse points, external review or formal verification is
made.
