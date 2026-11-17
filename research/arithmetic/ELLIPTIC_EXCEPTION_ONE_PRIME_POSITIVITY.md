# Every quadratic class passes each individual prime on the elliptic slice

All four quadratic fields on the existing elliptic slice have a positive
target coefficient at every single prescribed place. A three-residue
argument proves this for all but finitely many primes, and small exact
witnesses handle the exceptions. More generally every feasible finite
set of local conditions remains feasible after adding all but finitely
many new primes. These statements do not prove simultaneous feasibility
for every finite set in the two remaining quadratic classes.

## 1. Single-prime positivity in all four classes

Use the notation of
[EFA2](ELLIPTIC_EXCEPTION_FINITE_LOCAL_ASYMPTOTICS.md):
`Gamma=E_0(Q)`, `A={1,2,5q_0,10q_0}`, `q_0=2709173`, the global
squareclass map `alpha`, matching-target coefficient `C_E`, and
finite-place coefficients `mu_{S,a}`. The translated model is

\[
 E_0:y^2=x^3-10656x^2+27091730x,\qquad X=x-3552,
 \quad F(x)=128x-25865.
\]

Retain the exact points

\[
 P=\left(\frac{373248}{49},\frac{59459616}{343}\right),\quad
 Q=\left(\frac{3442688}{529},\frac{304670016}{12167}\right),\quad
 T=(0,0),
\]
\[
 S=\left(\frac{92179201}{6400},
                   \frac{553221611649}{512000}\right),
 \qquad R=P+T.                                      \tag{1.1}
\]

Their squareclasses are `alpha(P)=alpha(Q)=2`, `alpha(T)=10q_0`,
`alpha(S)=1`, and `alpha(R)=5q_0`.

**Theorem EOP1 (all four classes at every single place).** For every
finite prime `p` and every `a in A`,

\[
 0<\mu_{\{p\},a}\in\mathbb Q,\qquad
 N_{\{p\},a}(B)=\mu_{\{p\},a}C_E\log B
                              +O_{E,p}(\sqrt{\log B}).         \tag{1.2}
\]

The same positivity holds at the real place, where the coefficient is
`1/4`. In particular no individual prime obstructs either of the two
classes whose feasibility for every finite set remains undecided.

**The three-residue argument.** Let `p>=11` be good and `p=2 mod3`.
Rational two-torsion and Hasse's bound give `p` prime to
`N_p=#E_0(F_p)`: indeed `N_p` is even and `N_p<2p`. By EHD3 in
[the ambient Haar diagnostic](ELLIPTIC_EXCEPTION_AMBIENT_HAAR_DIAGNOSTIC.md),
the rational local image is its finite reduction image times a formal
subgroup `E_{e_p}(p)`. In the profinite completion of `Gamma`, the
reduction depends only on prime-primary factors different from `p`,
while the formal image is the Haar image of the fresh `p`-primary
factor. The global squareclass also depends only on the two-primary
factor.

Only the two finite zero residues `R_+,R_-` of `F` can be entirely
insoluble inside a projected residue coset. Off the three residues
`O,R_+,R_-`, `F` is a unit and the product is soluble. The pole residue
`O` always contains an open soluble part: arbitrarily deep formal
depths divisible by three exist in `E_{e_p}(p)`. These are the exact
local valuation tests of ELV1, with no assumption that `e_p=1`.

Fix a rational representative `P_a` of class `a`. The positive
profinite coset `P_a+2 widehat(Gamma)` preserves that class, and its
reduction contains the cyclic coset `bar P_a+2 Z bar P`. If
`[4]bar P!=O`, the latter has at least three elements: the order of
`[2]bar P` is greater than two. At least one of these residues is a
pole or a nonbranch residue. Conditional on that residue and the
coset, the fresh formal `p`-factor is unrestricted, since `p` is odd.
It has a positive open soluble subset there. Thus the class-specific
single-prime soluble set has positive intrinsic Haar measure.

**The complete finite exceptional set.** Exact duplication gives the
following reduced denominator of `x([4]P)`:

\[
 2^{12}3^6 5^2 7^2 139^2 163^2 122011^2 68819^2
                  11069^2 437900039^2 227812349^2.              \tag{1.3}
\]

Every factor displayed is prime. At a good prime reduction of `[4]P`
is `O` exactly when that denominator is divisible by `p`. The only
inert primes at least eleven in (1.3) are

\[
 11069,\quad68819,\quad227812349,\quad437900039.     \tag{1.4}
\]

The remaining bad-reduction primes of the elliptic curve that can
matter here are `1889,q_0`; the discriminant is
`2^9 5^2 7^3*1889*q_0^2`. At each of these six primes, exact rational
arithmetic gives `v_p(F(P))=v_p(F(R))=0`. ELV1 therefore supplies
an admissible soluble witness in class 2 or `5q_0`, respectively.

For the small primes use these exact witnesses:

| Class | Prime | Witness and local certificate |
| --- | --- | --- |
| `2` | `2` | `P in T+E_4(2)`, since `R=P+T` has depth four |
| `2` | `3` | `Q` has `X=2 mod3` and `v_3(y)=2` |
| `2` | `5` | `P` has integral `x` |
| `5q_0` | `2` | `R-S` has `v_2(x)=-12`, hence formal depth six |
| `5q_0` | `3` | `R` has depth three and `t(R)/27=-1 mod9` |
| `5q_0` | `5` | `R` has integral `x` |

The local assertions are ETF2--ETF3 and ETC1. All points used as
local witnesses are admissible; their curve equations, classes and nonvanishing
discriminants and resultant are checked in the retained certificate.
Each actual soluble point lies in an open local soluble set, so its
class has positive Haar coefficient. A globally irreducible cubic is
not needed as a witness: removing the finite global reducible-cubic
exceptions still leaves the full-field count.

Class 1 is positive for every finite set by
[ECO1](ELLIPTIC_EXCEPTION_FINITE_LOCAL_CLASS_ONE.md), and class
`10q_0` is positive by the certified Hasse seed and EFA2. This handles
all four classes. Primes one modulo three and the real place are
automatic. Rationality follows from EQA2, and the sharp error in
(1.2) is ESQ3 in
[the sharp finite-place note](ELLIPTIC_EXCEPTION_SHARP_FINITE_LOCAL_ERROR.md).
This proves EOP1. \(\square\)

## 2. Cofinitely many one-prime extensions of any feasible finite set

**Corollary EOP2 (cofinite extension from a certified coset).** Fix
`S,a` with `mu_{S,a}>0`. There are a rational point `P_0` and a
positive even integer `M` giving a profinite coset

\[
 \mathscr C=P_0+[M]\widehat\Gamma
\]

entirely in the class-`a` soluble set for the places in `S`.
For every prime outside the finite computable set

\[
 \mathcal E_{S,a}=\{2,3,5,7\}
   \cup\{p:p\mid\Delta_E M\}
   \cup\{p:p\mid\operatorname{den}(x([2M]P))\},     \tag{2.1}
\]

one has

\[
 \boxed{\mu_{S\cup\{p\},a}>0.}                    \tag{2.2}
\]

Given the coset certificate `(P_0,M)`, the exceptional set in (2.1)
is computable by finite integer arithmetic. Some listed primes can
also be harmless; no minimality is claimed.

**Proof.** EFA1 makes `U_{S,a}` open. Choose a basic clopen coset
contained in it and meeting the dense rational image. Its pullback
in `Gamma` contains `P_0+M Gamma` for some even `M`, and its full
profinite pullback contains `mathscr C`. Thus the image of the entire
profinite coset lies in the chosen clopen neighbourhood. This choice
is stronger than taking the closure of an arbitrary rational subset
of an open set.

For `p` outside (2.1), `p` is good, `p` does not divide `M`, and
`[2M]bar P!=O`. The reduction of `mathscr C` consequently has at
least three residues, since it contains
`bar P_0+M Z bar P`. All finite group restrictions defining
`mathscr C` involve only prime-primary factors dividing `2M`.
Thus its `p`-primary factor is unrestricted. At this good prime,
reduction is prime to `p` as in EOP1, so reduction is independent
of this fresh formal factor, even if `S` contains primes larger
than `p`.

Choose a supported pole or nonbranch residue, avoiding the at most
two wholly insoluble finite branch residues. Its intersection with
`mathscr C` has positive measure, and its formal soluble subset has
positive conditional measure for every `e_p`. Their intersection is
already soluble at every place of `S` and at `p`. Pushforward to the
actual joint closure gives (2.2). Primes one modulo three are
automatic in any event. \(\square\)

The exceptional set changes when the starting finite conditions or
coset change. EOP2 therefore does not inductively prove positivity
for every finite set in class 2 or `5q_0`. Neither theorem supplies
an all-places Hasse target in either class.

## 3. Uniform absolute loss at a new large prime

Write `R_p=im(Gamma -> E_0(F_p))` at a good prime. All Haar measures
below use the intrinsic rational joint closure of EFA2, normalized to
have total mass one; they are not measures in the full ambient product.

**Theorem EOP3 (large reduction image and vanishing absolute loss).**
There are constants `c_E,K_E>0` and `p_E` such that every good prime
`p>p_E` satisfies

\[
 |R_p|\ge c_E\log p.                              \tag{3.1}
\]

For every finite prime set `S`, every `a in A`, and every such prime,

\[
 0\le\mu_{S,a}-\mu_{S\cup\{p\},a}
                   \le\frac{K_E}{\log p}.          \tag{3.2}
\]

The constants in (3.2) are independent of `S` and `a`. In particular,
for each fixed feasible prefix `S,a`,

\[
 \frac{\mu_{S\cup\{p\},a}}{\mu_{S,a}}\longrightarrow1
 \quad(p\longrightarrow\infty),\qquad
 \mu_{\{p\},a}=\frac14+O_E\!\left(\frac1{\log p}\right).
                                                        \tag{3.3}
\]

The second limit concerns the finite-place leading coefficient: for
one fixed `p`, first let target height `B` tend to infinity in EOP1.
No uniformity of its counting error as `p` grows is asserted.

**Proof of the reduction-image bound.** The two displayed points `P,Q`
in (1.1) are independent modulo torsion; this is proved by the two
finite-prime certificates in ETC3a of
[the complete three-adic classification](ELLIPTIC_EXCEPTION_COMPLETE_THREE_ADIC_CLASSIFICATION.md).
The `P_1,P_2` there have short coordinates `X=x-3552`, so they are
exactly these points on the translated model.

Use the canonical height `q_x` and its normalization from
[the elliptic target count](ELLIPTIC_EXCEPTION_COUNT.md). There are
constants `C_0,C_1>0`, depending only on this fixed curve and these
fixed points, such that for every rational nonidentity point `D` and
every pair of integers `m,n`,

\[
 \left|\log H_X(D)-q_x(D)\right|\le C_0,
 \qquad q_x(mP+nQ)\le C_1(m^2+n^2).                \tag{3.4}
\]

At a good prime, a nonidentity point reducing to `O` has
`v_p(X(D))<=-2`. Writing `X(D)=A/B^2` in lowest terms therefore gives

\[
 \log H_X(D)\ge2\log p.                           \tag{3.5}
\]

Choose a fixed `epsilon>0` with `2C_1 epsilon^2<1`, and put
`L=floor(epsilon sqrt(log p))`. Consider the `(L+1)^2` points

\[
 mP+nQ,\qquad 0\le m,n\le L.
\]

If two had the same reduction, their nonzero difference
`D=uP+vQ` would satisfy `|u|,|v|<=L` and reduce to `O`.
Independence ensures that this difference is a nonidentity point.
Equations (3.4)--(3.5) would give

\[
 2\log p\le\log H_X(D)
       \le 2C_1L^2+C_0<\log p+C_0,
\]

which is impossible for sufficiently large `p`. The reductions of the
box are distinct, and `(L+1)^2` is bounded below by a positive constant
times `log p`. This proves (3.1); it does not require a full
Mordell--Weil basis or a saturation calculation.

**Proof of the loss bound.** At a sufficiently large good prime
`p=2 mod3`, ELV1 confines local failure, outside the finite excluded
base points, to the at most three reduction residues
`O,bar R_+,bar R_-` used in EOP1. The pushforward of intrinsic Haar
measure under reduction is uniform on the actual finite image `R_p`.
Thus the *unconditioned* local failure mass is at most

\[
 \frac3{|R_p|}\le\frac{K_E}{\log p}.              \tag{3.6}
\]

The finitely many excluded base points have intrinsic Haar measure
zero by EFA1. At primes `p=1 mod3`, local solubility is automatic, so
the failure mass is zero. The event removed when passing from
`mu_{S,a}` to `mu_{S union {p},a}` is an intersection with the same
unconditioned local failure event. Its measure is bounded by (3.6),
regardless of correlations with `S` or the class `a`. This proves
(3.2). Divide by the fixed positive `mu_{S,a}` for the first limit in
(3.3). Taking `S` empty, whose class mass is `1/4`, proves the second.
\(\square\)

EOP3 strengthens the asymptotic conclusion of EOP2, while EOP2 retains
its separate computable exceptional set from a specified coset. The
absolute error in (3.2) does not control the conditional loss uniformly
when `mu_{S,a}` becomes small. It supplies neither a growing-prime
counting estimate nor a summable all-primes tail, and does not decide
positivity of the all-places Haar measure.

## 4. Exact replay and proof scope

The [certificate](elliptic_exception_one_prime_certificate.json)
retains the rational generators, denominator factorization and finite
witness claims. The
[standard-library checker](verify_elliptic_exception_one_prime_positivity.py)
verifies exact group arithmetic, primality and complete multiplication
of (1.3), all exceptional-prime unit checks, the local depths and cube
unit, and admissibility and classes of the witnesses. Its receipt is
[elliptic_exception_one_prime_positivity_checks.json](elliptic_exception_one_prime_positivity_checks.json).

The all-prime assertions come from the three-residue proof and the
explicit complete exceptional set; the finite computations do not
establish an unrestricted statement by extrapolation. The cofinite
extension argument is a written group and openness proof. EOP3 follows
from the canonical-height box argument and a uniform bound on the
unconditioned failure mass; its height proof is not a finite test. No full
Mordell--Weil basis, formal verification or external review is claimed.
