# Exact inverse reachability for canonical polynomial densities

Internal research note, 22 November 2026. All assertions below concern the
canonical polynomial densities themselves. Equality of their finite moment
polynomials is a weaker condition; Section 5 gives an explicit positive-density
example separating the two notions. The algebraic identities hold over every
characteristic-zero field, using polynomial integration as coefficient division
and endpoint evaluation. The algorithm is effective over the rationals or any
coefficient field with exact arithmetic and equality tests. No SFC closure is claimed.

## 1. Definitions and exact density degrees

Fix nonnegative integers `r,u`, terminal order `n>=1`, and `s=r+u`.
Write `f_0=Beta(nr+1,nu+1)` and
`f_-1=Beta((n+1)r+1,(n+1)u+1)` as normalized polynomial densities on
`[0,1]`. At depth `ell` the represented moment degree is `k=n-ell`.
For a normalized previous/current pair `g,f`, define

\[
 a=(k+1)\int tg-k\int tf,
\]

\[
 K(t)=\int_0^t(t-x)(g(x)-f(x))\,dx
          -\frac1{k+1}\int_0^t(x-a)f(x)\,dx,
 \qquad \beta_k=k(k+1)\int_0^1 K.
\]

If `beta_k!=0`, the next canonical density is
`k(k+1)K/beta_k`. Signed densities and negative nonzero multipliers are
allowed. A canonical normalized chain stops if a required multiplier is zero.

**Degree theorem.** Wherever the normalized chain exists,

\[
 \deg f_0=ns,\qquad
 \deg f_\ell=(n+1)s+2\ell\quad(\ell\ge1).
 \tag{1}
\]

Moreover

\[
 f_\ell=t^{nr+\ell}(1-t)^{nu+\ell}H_\ell,
 \qquad\deg H_\ell=s\quad(\ell\ge1).
 \tag{2}
\]

No claim that `H_ell(0)` or `H_ell(1)` is nonzero is needed.

For `s>0`, the upper initial density has degree `(n+1)s`, exceeding
the degree `ns` of all lower-density terms in the first kernel. Its double
primitive therefore gives exact first-kernel degree `(n+1)s+2`.
For `s=0`, both initial densities equal one and
`K=t(1-t)/(2(n+1))`, so `beta_n=n/12` and `f_1=6t(1-t)`.
At this first step the kernel has the factors
`t^(nr+1)(1-t)^(nu+1)`: the defining integrals show these orders at zero;
the analogous right-endpoint formula gives the orders at one. An exponent-zero
axis can give only a simple endpoint zero, so one must not impose `K'=0`
there at the initial step.

For every later step, put `N=deg f_ell` and let `b` be its leading
coefficient. The previous density has degree at most `N-2` (exactly
`N-2` from depth two onward). Since the current
density vanishes at both endpoints, ordinary polynomial differentiation gives

\[
 K''=f_{\ell-1}-\frac{k+2}{k+1}f_\ell
                    -\frac{t-a}{k+1}f_\ell'.
\]

Its leading coefficient is `-(N+k+2)b/(k+1)`, which is nonzero.
Consequently

\[
 [t^{N+2}]f_{\ell+1}
   =-\frac{k(N+k+2)}{\beta_k(N+1)(N+2)}[t^N]f_\ell\ne0.
 \tag{3}
\]

The defining kernel has `K=K'=0` at both endpoints at these later steps.
The right side of its second-derivative equation has factors
`t^(nr+ell-1)(1-t)^(nu+ell-1)`. Integration from either endpoint gives
the factors in (2) for `ell+1`. Its exact total degree proves exact modifier
degree `s`. This proves (1)--(3) inductively, with no positivity assumption.

## 2. Finite candidate enumeration from a single density

Suppose a nonzero normalized polynomial `f` is proposed as a descended
canonical density of represented degree `k>=1`. Let

\[
 D=\deg f,\quad a_0=\operatorname{ord}_0 f,\quad
 a_1=\operatorname{ord}_1 f,\quad\delta=D-a_0-a_1.
\]

For any possible depth `ell>=1`, (1)--(2) force

\[
 D=(k+\ell+1)s+2\ell=(k+1)s+\ell(s+2),
 \tag{4}
\]

\[
 \boxed{\delta\le s\le
       \left\lfloor\frac{D-2}{k+2}\right\rfloor.}
 \tag{5}
\]

The left inequality follows because the required endpoint factors have total
degree `D-s`, while the actual endpoint orders have sum `D-delta`.
The right inequality uses `ell>=1` in (4). For every integer `s` in
(5), the only possible depth and terminal order are

\[
 \ell=\frac{D-(k+1)s}{s+2},\qquad n=k+\ell.
 \tag{6}
\]

Retain only positive integer depths. For each retained candidate, try the
finitely many integer splits `r+u=s`, and reject those violating
`nr+ell<=a_0`, `nu+ell<=a_1`, or the corresponding endpoint divisibilities
of the proposed previous density. The proposed previous degree must be
`ns` if `ell=1`, and `D-2` if `ell>=2`.

This produces a complete finite list, independent of any numerical cutoff.
There are at most `(S+1)(S+2)/2` triples before endpoint filtering, where
`S=floor((D-2)/(k+2))`; if `D<2`, no descended candidate exists.
The initial case `ell=0` is tested separately: `n=k`, `s=D/k` must be a
nonnegative integer, and each of the finitely many splits `r+u=s` gives
two prescribed Beta densities for comparison.

## 3. A uniquely determined backward step

Let `f` be the current density of represented degree `k`, and `g` the
previous density of represented degree `m=k+1`. At a candidate depth at
least two, write

\[
 D=\deg f,\quad d=D-2=\deg g,\quad
 b=[t^D]f,\ b_1=[t^{D-1}]f,\quad
 c=[t^d]g,\ c_1=[t^{d-1}]g.
\]

Here `D>=4` and `d>=2`. Solving the differentiated forward relation for
the earlier density `h` gives

\[
 h=\frac{(m+2)g+(t-\alpha)g'+(\beta/m)f''}{m+1}.
 \tag{7}
\]

An actual earlier density has degree at most `d-2`, including the last
backward step to the lower initial Beta density. The top two coefficients
in (7) must therefore vanish. These conditions determine uniquely

\[
 \boxed{\beta=-\frac{m(m+D)c}{D(D-1)b}},
 \tag{8}
\]

\[
 \boxed{\alpha=
 \frac{(m+d+1)c_1+(\beta/m)(D-1)d\,b_1}{dc}.}
 \tag{9}
\]

All denominators are nonzero in characteristic zero. In particular the
backward multiplier (8) is automatically nonzero; its sign is unrestricted.

The degree cancellation alone does not certify reachability. At every
backward step one must check all required endpoint factors and the exact
earlier degree: `D-4` while the earlier depth is positive, and `ns` when
the earlier depth is zero. The latter can demand cancellation of additional
coefficients beyond the first two.

These backward relations recover the actual canonical operation, not merely
a differential identity. The candidate endpoint factors imply `g(0)=g(1)=0`
and `f(0)=f(1)=f'(0)=f'(1)=0`. If `f,g` have mass one, integration by parts
in (7) gives

\[
 \int h=1,\qquad
 (m+1)\int th-m\int tg=\alpha.
\]

The forward kernel formed from `h,g` consequently satisfies
`K''=beta f''/(m(m+1))`. Its value and first derivative vanish at zero,
as do those of `f`; hence `K=beta f/(m(m+1))`. Its mass gives exactly
the multiplier (8). Thus the reconstructed forward normalization is `f`.

## 4. A complete canonical-density membership test

Enumerate the candidates in Section 2. For each one, apply the unique
backward step until depth one, checking the degrees, endpoint factors and
normalizations at every stage. When the backward recursion reaches a
depth-one pair, its previous density must equal `Beta(nr+1,nu+1)`
exactly. Form the first kernel from that Beta density and
`Beta((n+1)r+1,(n+1)u+1)`. Reject the candidate if this kernel has zero
mass; otherwise compare its normalized density with the proposed depth-one
density. This explicit rejection does not assume a separate uniform theorem
about the first multiplier.

**A candidate is accepted if and only if the proposed pair occurs in its
nondefective canonical density chain.** Necessity follows from (1)--(9).
For sufficiency, the exact initial comparison supplies the first step,
and the integration-by-parts argument in Section 3 rebuilds each subsequent
canonical forward step. The algorithm is finite by (5)--(6).

This is an exact membership theorem for a supplied polynomial pair. It does
not prove that all actual chains remain nondefective, since the construction
starts from an already supplied normalized pair.

## 5. Why this does not classify moment-polynomial reachability

The map `f -> integral (z-t)^k f(t)dt` sees only the first `k` moments.
The canonical density carries additional data. This distinction persists
even among strictly positive densities.

For `k=2`, set

\[
 g(t)=1,\qquad f(t)=6t(1-t),
\]

\[
 \Delta(t)=\frac{d^3}{dt^3}[t^4(1-t)^4]
 =24t(1-t)(1-9t+21t^2-14t^3),
 \quad\widetilde f=f+\frac{\Delta}{1000}.
\]

Three integrations by parts show
`integral t^j Delta=0` for `j=0,1,2`. Thus `f` and `f_tilde`
have exactly the same represented quadratic. Both are normalized and
strictly positive on `(0,1)`: the crude bound
`|24(1-9t+21t^2-14t^3)|<=1080` makes the bracket in
`f_tilde=t(1-t)[6+24(1-9t+21t^2-14t^3)/1000]` at least `123/25`.

The pair `g,f` is an actual canonical pair at terminal order `n=3`,
profile `r=u=0`, depth one. The pair `g,f_tilde` represents the same
moment-polynomial pair, yet `f_tilde` cannot be any descended canonical
density of represented degree two: it has `D=5`, endpoint orders `1,1`,
and `delta=3`, contradicting `delta<=floor((D-2)/(k+2))=0`.
It is not an initial density either, since initial degree `2s` is even.

Consequently excluding a supplied density pair does not by itself exclude
its finite moment-polynomial pair. Actual common-zero or resultant claims
must be checked at the polynomial level. The present theorem provides
finite exact canonical-density reachability, not a proof of SFC.

## 6. Verification

The companion `verify_inverse_reachability.py` independently implements the
forward kernel and backward coefficient cancellation. Its finite actual
controls include unequal profiles, an axis, the zero profile, first depth
and deeper backward steps. It also verifies the exact positive-density
noninjectivity example. These controls check the implementation and written
identities; the uniform membership theorem is the proof in Sections 1--4.
