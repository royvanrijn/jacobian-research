# A three-adic upper sieve on the elliptic exceptional slice

The existing elliptic slice has a strict three-adic loss among its
matching-resolvent targets. This gives a smaller upper coefficient for
Hasse failures. The argument uses finite quotients of the local elliptic
group and does not assume that displayed rational points generate the
Mordell--Weil group.

## 1. Scope and the two upper bounds

Retain the curve, parameter, target height and matching-target coefficient
from [the elliptic height count](ELLIPTIC_EXCEPTION_COUNT.md):

\[
 E:T^2=X^3-10758382X+6600759744,\qquad
 x=X+3552,\qquad v=\frac{3(x+432)}{6400},
 \qquad C_E=\frac{\pi}{3\sqrt{R_x}}.                 \tag{1.1}
\]

Here `R_x` is the determinant of the canonical quadratic-height matrix in
the normalization of EEC1. Write `A_E(B)` for the number of admissible
Hasse-failing targets on this slice of actual projective height at most
`B`. The known unrestricted matching count is
`C_E log B + O_E(sqrt(log B))`.

**Theorem ETS1 (finite-quotient upper bounds).** For every fixed integer
`J>=1`,

\[
 A_E(B)\le
 \left(\frac{29}{52}+\frac9{13\,27^J}\right)C_E\log B
                 +O_{E,J}(\sqrt{\log B}).           \tag{1.2}
\]

In particular the first finite quotient gives

\[
 \boxed{A_E(B)\le\frac7{12}C_E\log B
                    +O_E(\sqrt{\log B}).}          \tag{1.3}
\]

**Corollary ETS2 (limiting upper proportion).** One has

\[
 \boxed{\limsup_{B\to\infty}
       \frac{A_E(B)}{C_E\log B}\le\frac{29}{52}.}    \tag{1.4}
\]

No uniform error in `J` is asserted. In particular (1.4) does not assert
an `O(sqrt(log B))` remainder with coefficient `29/52`, a density-zero
result, an asymptotic for `A_E`, or infinitude of Hasse failures.

## 2. Two three-adic obstruction regions

The translated equation is

\[
 T^2=x(x^2-10656x+27091730).                         \tag{2.1}
\]

It has good reduction at three, and its four reduced points are the
identity and the three points with `T=0`, `X=0,1,2`. The reduction group
is `V_4`. Denote the reduction of the rational two-torsion point by
`bar T_0=(X,T)=(0,0)`, equivalently `(x,T)=(0,0)`.

**Lemma ETS3 (the two-torsion residue disc is insoluble).** Every
admissible matching target whose elliptic point reduces to `bar T_0`
has no source point over `Q_3`.

**Proof.** In this residue disc the second factor of (2.1) is a unit.
Thus `v_3(x)=2v_3(T)>=2`, apart from the exact two-torsion point, which
is inadmissible. Formula (1.1) gives `v in 27 Z_3`.

The exact shifted cubic factor is

\[
 h(W+1)=W^3+\frac35W^2+
       \left(\frac{341787}{819200}-v\right)W
       +\frac{10v^2}{9}-\frac{851209v}{245760}
                         +\frac{769467}{819200}.    \tag{2.2}
\]

When `v in 27 Z_3`, each nonleading coefficient has valuation exactly
one: both displayed constant rational terms in the last two
coefficients have valuation one, while their changes have valuation
at least two. This cubic is Eisenstein. Its discriminant has valuation
three, since for a monic cubic with `v_3(a)=v_3(b)=v_3(c)=1`, the term
`-4b^3` in the discriminant has valuation three and all other terms have
strictly larger valuation. The cubic therefore has no `Q_3`-root, and
the matching quadratic discriminant has odd valuation and is not a
square. Neither field component supplies a local source point.
\(\square\)

For `n>=1`, put

\[
 E_n=\{P\in E_1(\mathbb Q_3):v_3(-X(P)/T(P))\ge n\},
\]

where `E_1` is the kernel of reduction. The identity belongs to every
`E_n`. A nonidentity point in `E_n minus E_{n+1}` has formal depth `n`.

**Lemma ETS4 (obstructed formal annuli).** An admissible matching
target from a point of formal depth `n` has no `Q_3` source point when
`3` does not divide `n`.

**Proof.** The standard local parameter gives
`v_3(X)=-2n`, `v_3(T)=-3n`. Hence (1.1) gives
`v_3(v)=1-2n`. In the unshifted cubic

\[
 h(W)=W^3-\frac{12}{5}W^2+(s_0-v)W+
       \frac{10}{9}v^2-\frac{10s_0}{9}v+\frac{125}{1024},
 \qquad s_0=\frac{1816347}{819200},                 \tag{2.3}
\]

the coefficient valuations in degrees zero through three are
`-4n, 1-2n, 1, 0`. The intermediate points are strictly above the
endpoint Newton edge `(0,-4n)--(3,0)`. Every root consequently has
valuation `-4n/3`. If `3` does not divide `n`, none is in `Q_3`.
The quadratic discriminant is `D=-3x/1600`, of valuation `1-2n`,
which is odd. Thus the quadratic has no local point either.
\(\square\)

Only necessity is used: depths divisible by three are retained by this
upper sieve, without being declared locally soluble.

## 3. The rational image in the finite local quotients

The exact local group input is Milne, *Elliptic Curves* (2006),
Chapter II, Theorem 4.1, printed pp. 62--63, for
`E/E_1=E(F_3)`, the parameter valuations and
`E_n/E_{n+1}=F_3`; and Lemma 5.9, printed p. 67, for the parameter
addition estimate on a short integral Weierstrass model:
[primary book](https://www.jmilne.org/math/Books/ectext6.pdf),
[retained copy](sources/milne-elliptic-curves-2006.pdf).

In particular, if `t=-X/T` and `P,Q in E_n`, that estimate gives

\[
 t(P+Q)=t(P)+t(Q)\pmod {3^{5n}}.                  \tag{3.1}
\]

For a point of exact depth `n>=1`, doubling preserves depth and (3.1)
then gives `t([3]P)=3t(P)+O(3^{5n})`, of exact depth `n+1`.

Use the rational points

\[
 T_0=(-3552,0),\qquad
 P_2=\left(\frac{1563680}{529},
                  \frac{304670016}{12167}\right).  \tag{3.2}
\]

Both satisfy the curve equation exactly. Their reductions are distinct
nonidentity points `(0,0)` and `(2,0)` in `E(F_3)`, so the rational
reduction image is the whole group of order four. Exact duplication
gives

\[
 X(2P_2)=\frac{89497017424304041335649}{1002118368684773376},
\]
\[
 T(2P_2)=
 -\frac{26756035205526466434703545465981295}
        {1003179235240703925530394624}.             \tag{3.3}
\]

Their valuations are `-4,-6`, so `2P_2` has formal depth two.
It follows from (3.1) and the quotient sizes that, for each `N>=3`,
`2P_2` generates `E_2/E_N`, of order `3^{N-2}`.

Let `G_N=E(Q_3)/E_N`, of order `4*3^{N-1}`. The subgroup `W_N`
generated by the images of `T_0,P_2` has order `4*3^{N-2}`: the
point `P_2` has order `2*3^{N-2}` in this quotient, and `T_0` is not
in its cyclic subgroup because their nonzero reductions are distinct.
Thus `W_N` has index three in `G_N`.

The image `H_N` of the entire group `E(Q)` contains `W_N`. Since this
index is prime, there are just two possibilities,

\[
 H_N=W_N\quad\hbox{or}\quad H_N=G_N.               \tag{3.4}
\]

This is a statement about finite quotient images. No global generator,
regulator, or saturation computation enters it.

## 4. Counting the retained residue classes

First take `N=3`. If `H_3=W_3`, it has twelve classes. The forbidden
`bar T_0` disc has three of them, and the identity disc has three;
its two nonidentity classes have depth two and are forbidden by ETS4.
At most seven of twelve classes remain. If `H_3=G_3`, there are
thirty-six classes. Nine are in the forbidden `bar T_0` disc, and eight
of the nine identity-disc classes have depth one or two. At most
nineteen of thirty-six remain, a smaller fraction than `7/12`.

For the general finite bound take `N=3J`. In the smaller-image case,
the identity disc begins at `E_2`. Conditional on that disc, its depth
`n`, for `2<=n<N`, has fraction

\[
 \frac23\,3^{-(n-2)},
\]

and the as yet untested tail `E_N` has fraction `3^{-(N-2)}`.
Two of the four residue discs are retained in full, `bar T_0` is
discarded, and the identity disc retains only depths divisible by
three and the tail. The resulting retained fraction is

\[
 \frac12+\frac14\left(
  \sum_{k=1}^{J-1}\frac23\,3^{-(3k-2)}
                         +3^{-(3J-2)}\right)
       =\frac{29}{52}+\frac9{13\,27^J}.            \tag{4.1}
\]

In the full-image case the identity disc begins at `E_1`; its analogous
fraction is

\[
 \frac{27}{52}+\frac3{13\,27^J},                  \tag{4.2}
\]

which is smaller than (4.1).

Each class in `H_N` pulls back to one coset of the same finite-index
subgroup of `E(Q)`. The canonical-height lattice argument of EEC1,
including its rational torsion subgroup, counts each such coset with
exactly `1/|H_N|` of the full point main term and error
`O_{E,N}(sqrt(log B))`. The retained classes are symmetric under
negation. Except at finitely many parameters, targets have the two
preimages `P,-P`, so dividing by two preserves the retained fraction
of the matching-target main term. The finite inadmissible and reducible
cubic exceptions do not change the estimate. This proves ETS1.

For each fixed `J`, take the upper limit in `B` in (1.2). Letting `J`
tend to infinity only after that upper limit gives ETS2. This order of
limits requires no uniform error term.

The inequalities are upper bounds for the existing elliptic slice.
The other two residue discs and the retained formal depths may contain
further local failures. No claim of exact three-adic density or of
control at all the varying bad primes is made.

## 5. Exact replay and assurance

The [standard-library checker](verify_elliptic_exception_three_adic_sieve.py)
verifies the rational points, exact duplication, the shifted cubic,
bounded valuation and Newton controls, the complete finite quotient
counts, the rational coefficients (4.1)--(4.2), and corruption rejection.
Its receipt is
[elliptic_exception_three_adic_sieve_checks.json](elliptic_exception_three_adic_sieve_checks.json).
The unbounded statements follow from the written local-group and lattice
proofs above; bounded controls do not establish their quantifiers.
This is internal exact replay and contextual mathematical review, with
no formal verification or external review claimed.
