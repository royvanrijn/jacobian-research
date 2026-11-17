# All-rank integrated monodromy and coordinate rigidity

This completes the group calculation for every integrated-power
field of KMFI109, rather than only rank4. It also bounds the
coefficient field for actual affine completions and identifies
exactly which coordinate equivalences preserve normalized sweep
lines. These are construction-class theorems. Exceptional11 and
the full faithful-transitive Keller action spectrum remain open.

## 1. Complete geometric and arithmetic groups

Let k have characteristic zero, let h in k[t] have degree d>=1,
and put r>=1, e=r+1, N=rd+1. Up to the polynomial target sign
and the spectators in KMFI109, the field is
\[
 k(A,B)\subset k(A,t),\qquad
                  B=I_r(t,A)=\int_0^t(A-h(s))^r\,ds.       \tag{1}
\]

**KMFI112.** If d1, the geometric action is regular cyclic C_e.
The arithmetic group is
C_e semidirect Gal(k(zeta_e)/k), with exact constants k(zeta_e).
It has no Keller model, even after common independent stabilization
and arbitrary birational source/target recharts.

For every d>=2 the complete action is
\[
 G_{\rm geom}=
 \begin{cases}
 S_N,&r\text{ odd},\\
 A_N,&r\text{ even}.
 \end{cases}                                               \tag{2}
\]
For odd r, arithmetic is S_N and the closure has no new constants.
For even r, define
\[
                 \Delta=(-1)^{rd/2}(rd+1).                 \tag{3}
\]
Arithmetic is A_N iff Delta is a square in k, otherwise S_N.
The exact normal-closure constants are k(sqrt(Delta)).

Here is an unbounded proof. At generic A the d roots h(t)=A
are simple, each with ramification index e. The critical-value
parametrization a=h(t), b=I_r(t,h(t)) satisfies
\[
                    \frac{d^r b}{da^r}=r!t.               \tag{4}
\]
Indeed the first r successive endpoint contributions vanish;
before the r-th step the remaining integrand still vanishes
at s=t. Differentiation in the finite separable extension
k(t)/k(a) preserves the subfield k(a,b), since the derivation
of k(a) has a unique extension to each finite separable field.
Thus t belongs to k(a,b). The parametrization is birational,
and distinct generic critical points have distinct values.
This says nothing about every specialized A.

We need primitivity even when e is composite. Over an algebraic
closure of k(A), suppose I_r=U composed with V, with both degrees
greater than1. Choose a critical point c of U. Every point of
V^{-1}(c) is critical for the composite and has the same value.
Distinct critical values force V^{-1}(c) to be one point, so
V=a(t-t0)^b+c. If U had another critical point c', its inverse
under V would contain b distinct points, another contradiction.
Therefore U has just the critical point c and is
a'(z-c)^a+c0. The composite is a monomial up to affine changes
and has just one critical point. Our polynomial has d>=2
critical points, a contradiction.

For completeness, this indecomposability is equivalent to
primitive monodromy here. An intermediate curve dominated by
P1 has genus0 by Riemann--Hurwitz. The unique place over infinity
in a polynomial cover gives a unique intermediate place there.
Choosing it as infinity makes both intermediate functions
polynomials. Thus a nontrivial block system would give precisely
the polynomial decomposition just excluded.

Each finite inertia is a single e-cycle fixing r(d-1) sheets.
If r(d-1)>=3, the primitive-cycle theorem gives G>=A_N
([Jones, Corollary1.3](https://arxiv.org/pdf/1209.5169)).
This is the composite-cycle result, whose proof uses the
classification of finite simple groups; the prime-cycle theorem
alone does not suffice for all r.

The small cases are elementary. For r1, the branch transpositions
generate a transitive group, hence S_N by the connected-edge graph.
For r2,d2 there are two3-cycles on five sheets. Transitivity forces
their supports to overlap in one sheet, and they generate A5;
the checker exhausts the representative (123),(345).
These cases cover every r(d-1)<3. The sign of an e-cycle is
(-1)^r, proving (2). A relative geometric S_N already is full.
For even r, the following square class bounds the full geometric
group over algebraically closed k by A_N, so no enlargement
from motion in A is hidden.

In fact, P=I_r(t,A)-B has leading coefficient h_d^r/N for even r,
and P'=(A-h)^r. The exact discriminant identity is
\[
 \operatorname{Disc}(P)=
 (-1)^{rd/2}N
 \left(\frac{\operatorname{Res}(P,A-h)^{r/2}}
                   {h_d^{r/2}}\right)^2.                 \tag{5}
\]
Here N is odd and N(N-1)/2 has parity rd/2.
The relative group A_N and (5) prove (3), the arithmetic statement
and the exact constants. For odd r the geometric group is
already S_N, so arithmetic cannot enlarge it.

When d1, writing w=A-h(t) integrates (1) to an affine expression
in w^e plus a polynomial in A. A triangular target change gives
exactly w^e=Z with independent A,Z. Its splitting field is
k(A,Z)(zeta_e,w), giving the claimed arithmetic semidirect
product and constants. The nontrivial regular geometric action
is absolutely excluded by KMFI2, not only by this constructor.

## 2. Actual completions with a uniform coefficient-degree bound

**KMFI113.** For d>=2 the KMFI109 input has an actual Keller
completion on ordinary A^(r+2), with exact action (2).
All coefficients lie in an explicitly chosen extension k'/k
of degree at most r(d-1). Over k', the arithmetic group and
constants are still given by KMFI112 with its square test in k'.

Here is the exact marked-point construction, not an existence
claim about an unspecified formal solution. Choose t0 in k with
h'(t0)!=0 and with distinct critical values over A0=h(t0).
There are only finitely many bad A and k is infinite, so such
a t0 in k exists by (4). Use p*=(t0,0,...,0), and compute its
target Y*=X(p*). Form
\[
 Q(t)=\frac{I_r(t,A0)-I_r(t0,A0)}{(t-t0)^{r+1}}.            \tag{6}
\]
Its degree is r(d-1), and it is coprime to A0-h(t).
The multiplicity at t0 is exactly r+1; any other root of h=A0
has a different critical value. Hence every root alpha of Q
has gamma0=A0-h(alpha)!=0.

Choose any irreducible factor of Q and any root alpha in its
degree-at-most-r(d-1) field. Recover the other p coordinates
successively from the triangular first r target equations.
Equation(6) verifies the last target equation. The entire source
point is over k'=k(alpha), is regular, and maps to Y*.

At p*, the first r coordinates of DX have determinant h'(t0),
so DX is injective. Since P=[-p|I]DX is nilpotent, kerP has
dimension1: DX maps it injectively into the one-dimensional
kernel of [-p|I]. Thus P is a single regular nilpotent Jordan
block. A cyclic vector and its P powers give a Jordan frame
over k, without adjoining constants. KMFI103--104 then use
only finite field operations over k': the polynomial first stage,
the inverse linear jet, the truncated curve and the lower
unitriangular correction. They specify every polynomial coefficient,
give the explicit source inverse, and prove the whole nonzero
constant Jacobian. The field pair after one common variable
is reconstructed in both directions.

For h=t^d one may always choose t0=1, with no exceptional-value
search. Write
\[
 c_{r,d}=I_r(1,1)
 =\sum_{j=0}^r\frac{(-1)^j\binom rj}{dj+1}
 =\frac{r!d^r}{\prod_{j=0}^r(dj+1)}\ne0.                 \tag{7}
\]
For a d-th root zeta of1, I_r(zeta,1)=zeta*c_(r,d).
Thus the d critical values are distinct and (6) is explicit over Q.
This generalizes the quartic-field rank4,d2 control, not its
already-known A9 abstract action.

There is also an exact geometric normalization consequence.
For r>=2, apply the specified higher-weight completion with
generic centering as in KMFI106--108, which may be chosen above.
The original finite normalization is A^(r+1), with class group0.
Its inertia is one e-cycle and N-e fixed sheets; its e-th power
is the identity. Over algebraically closed constants the new C0
has N split primes, precisely two retained and N-2 omitted
unramified. There is one other omitted prime of indexe.
Consequently
\[
 \operatorname{Cl}(Z_F)=\mathbb Z^{N-1},\qquad
                 [K_{Z_F}]=(r,0,\ldots,0).                \tag{8}
\]
This refers to that completed representative after algebraic
closure, not its class group over k' or every model of the action.
For r1 the nonlinear exceptional stage needs a separate calculation;
(8) is not asserted for r1.

## 3. Coordinate changes that really preserve normalized sweep lines

**KMFI114.** Let S_X and S_X' be normalized sweeps of the same
rank, both with determinant gamma^r. Suppose polynomial
automorphisms T of the target and Psi of the source satisfy
T composed with S_X = S_X' composed with Psi, and Psi and its
inverse carry gamma-lines onto gamma-lines. Then T is affine,
Psi has no gamma-translation, and for nonzero constant lambda
and invertible constant B,
\[
 T(Y)=
 \begin{pmatrix}\lambda&0\\v&B\end{pmatrix}Y+c,\qquad
 \Psi(\gamma,p)=
             (\lambda\gamma,(v+Bp)/\lambda).              \tag{9}
\]
In particular the induced change of p is affine.

The fiber-preserving polynomial automorphism has
p'=phi(p), gamma'=a(p)gamma+b(p), where phi is a polynomial
automorphism and a is a polynomial unit, hence constant lambda.
If the highest homogeneous degree of T were m>1, the coefficient
of gamma^m in T(X(p)+gamma(1,p)) would be T_m(1,p).
It vanishes identically because the other side is affine in gamma.
A homogeneous polynomial vanishing on the full affine chart
(1,p) is zero, a contradiction. Thus T=M Y+c.
Comparison of directions gives M(1,p)=lambda(1,phi(p)),
forcing the matrix form and the affine phi in (9).
Finally the critical divisor of either sweep is exactly gamma0.
The source automorphism carries it onto itself, so b(p)=0.
This proves (9); it does not classify arbitrary equivalences
that change the sweep-line foliation.

A nonfinite rank4 input cannot become finite triangular by any
polynomial source and target automorphisms, even without the
line-preserving hypothesis: finiteness itself is invariant.
More strongly, KMFI112 excludes exceptional11 through every
rank of this integrated field class, even after common independent
stabilization and arbitrary birational changes, because degree
and the exact permutation action are invariants of a field pair.
In degree11 the possible (r,d) are (1,10),(2,5),(5,2),(10,1);
their actions are respectively S11,A11,S11,C11.

One tempting coordinate theorem is inapplicable. Polynomial
nilpotency of P alone does not make P a Nijenhuis tensor.
With the convention
\[
 (N_P)^i{}_{jk}
 =\sum_l\left(P^l_j\partial_lP^i_k-P^l_k\partial_lP^i_j
       -P^i_l(\partial_jP^l_k-\partial_kP^l_j)\right),
\]
the full-rank rank4,h=t^2 example already has
(N_P)^1_(12)=4t-2, which is nonzero. A local nilpotent-Nijenhuis
normal-form theorem cannot supply a general triangularization.
A smooth local form would also require a separate global
polynomial argument.

## 4. Evidence and the unresolved objective

The checker tests exact higher critical-value identities, even-rank
discriminant signs, the marked cofactor and coefficient bounds,
higher-rank constructor identities, the small A5 generation,
and the nonzero torsion countercontrol.
The arguments for all ranks and all degrees are written above;
finite controls are not an unbounded search certificate.
The proof uses Jones's composite-cycle Corollary1.3, including its
classification input, and KMFI103--104's written whole-map proof.
No new expanded final determinant or independent whole-proof
review is claimed.

This decides this entire integrated/triangular architecture, with
actual affine completions for its nonregular groups and an
absolute obstruction for its regular cyclic members.
It does not prove that every normalized full-rank sweep has such
a presentation. It supplies neither the nonfinite rank4 PSL11
field nor an absolute exceptional11 action obstruction.
The full faithful-transitive spectrum remains open.

[Checker](verify_all_rank_integrated.py) ·
[exact controls](all_rank_integrated_model.json) ·
[receipt](all_rank_integrated_receipt.json) ·
[primary input](sources/all-rank-integrated-manifest.json).
