# Algebraic infinity limits and complete linear growth exclusions

14 November 2026. **QVL1.** The proper-power boundary in
[QVM1](CYCLIC_INFINITY_MATE_OBSTRUCTION.md) can be closed when the
remaining scalar coefficients are bounded at infinity. In the
quadratic chart [a,v]=h/v, every finite operator
S=v*H(a)+S0 with deg H>=2 and bounded coefficients in S0 has no
finite Weyl mate, over any finite algebraic coefficient cover.
No root-multiplicity condition on H is needed.

The new ingredient is a transported infinity obstruction. After an
algebraic position change and momentum shift, a nonlinear
constant-coefficient limit excludes every finite-order mate.
Repeated momentum roots supply exactly such a change. A related
Euler limit excludes even a nonconstant linear limit because an
algebraic coefficient field cannot supply a logarithmic conjugate.

These results constrain scalar finite operators. They do not exclude
every Long representative, matrices, arbitrary projective carriers or
formal series of unbounded operator order.

## A nonlinear ordinary infinity limit has no finite mate

Let C be algebraically closed of characteristic zero, let h be a
nonzero differential constant, and let F/C(u) be any finite
algebraic function field with derivation partial_u. At a point
over u=infinity, use the leading Puiseux exponent as deg_u.
Consider a finite scalar operator
\[
 P=\sum_j p_j(u)E^j,\qquad [E,f]=h f_u,\qquad
 \deg_u p_j\le0,\qquad
 P_\infty(E)=\sum_j \lim_{u\to\infty}p_j(u)E^j .
                                                               \tag{1}
\]
If deg P_infinity>=2, there is no finite Q over any finite
algebraic coefficient cover such that [P,Q]=kappa is a nonzero
differential constant. Higher-order coefficients of P are allowed
to tend to zero; P_infinity need not have P's generic order.

Indeed, write the leading coefficient-valued term of Q as
u^M*K(E), with M rational and K a nonzero polynomial over C.
Normal ordering gives the coefficient at exponent M-1:
\[
          [P,Q]_{M-1}=hM\,P_\infty'(E)K(E),\qquad
                  \deg_u[P,Q]\le M-1.                        \tag{2}
\]
All decaying terms in P and all lower coefficient exponents in Q
contribute strictly below M-1. If M>1, the displayed coefficient
is nonzero and has positive u degree. If M<1, the upper bound
is negative. Thus M=1, and a constant commutator would require
P_infinity'(E)*K(E) to be a nonzero constant. Its degree is
positive, a contradiction. This proves the assertion for all
partner orders at once.

The restriction to finite scalar differential operators matters:
the inverse of P_infinity'(E) in a division ring is not such
an operator. The positive pair E,u shows why a linear ordinary
limit is outside the exclusion.

## The Euler version detects the logarithmic boundary

Use instead [E,f]=h*v*f_v and the Puiseux valuation at v=infinity.
Suppose P has bounded coefficients and nonconstant limit F(E).
For a mate with leading term v^M*K(E), the leading coefficient is
\[
       [P,Q]_M=
              \bigl(F(E+hM)-F(E)\bigr)K(E).                   \tag{3}
\]
For M nonzero this polynomial is nonzero. A constant commutator
would force M=0. But then the leading polynomials commute and
every remaining coefficient has negative v exponent, again
excluding a nonzero constant. Thus an Euler limit of any positive
degree has no finite algebraic-coefficient mate.

This is not an ordinary-coordinate assertion at a finite point.
The identity [v*partial_v,log v]=1 suggests the missing conjugate,
but log v belongs to no finite algebraic extension of C(v).
The proof (3) uses only finite differences and Puiseux exponents.

## Transport through algebraic canonical changes

Return to a coefficient field with derivation delta and
[a,f]=h*delta(f). After a finite algebraic extension choose
f nonzero, alpha and u in that field such that
f*delta(u)=1. Then
\[
          E=f(a-\alpha),\qquad
                 a=\alpha+f^{-1}E,\qquad [E,u]=h .            \tag{4}
\]
Both directions preserve finite scalar differential order.
Any hypothetical finite mate remains a finite mate after this
change. Consequently (1)--(2) apply whenever u tends to infinity
and the transformed coefficients have the stated limits.
A change with f*delta=v*partial_v instead uses (3).

These changes need not be regular at h=0. In particular a root
alpha may have an h pole. This is legitimate for a generic-fiber
nonexistence proof: a finite regular mate would still become a
finite generic mate. It is not a bounded formal normalization
theorem for arbitrary quantum lifts.

## Repeated roots produce a nonlinear limit

Let q=v^d/d with d>=2, so delta=v^(1-d)*partial_v.
Suppose
\[
       S=vH(a)+S_0,\qquad H\in C[a],\qquad
              \deg_v\text{ every coefficient of }S_0\le0 .
                                                               \tag{5}
\]
Let alpha be a root of H of multiplicity k>=2, and put
\[
 E=v^{1/k}(a-\alpha),\qquad
 u=\frac{v^{\,d-1/k}}{d-1/k},\qquad [E,u]=h .                 \tag{6}
\]
This takes only a finite algebraic coefficient extension, and
u tends to infinity with v. Write b=a-alpha.
If H(alpha+b)=sum_(j>=k)c_j*b^j, then in the new normal order
the leading coefficient term of b^j is v^(-j/k)*E^j.
Each internal contraction decreases its coefficient v exponent
by the additional positive amount d-1/k.

Thus all coefficients of v*H(a) become bounded; its limit is
c_k*E^k. Every term from S0 remains bounded and contributes
only a scalar to the limit. If sigma0 is the value at alpha
of the polynomial made from S0's coefficient limits, then
\[
                          S_\infty(E)=c_kE^k+\sigma_0 .      \tag{7}
\]
The degree is k>=2, so (1)--(2) exclude every finite mate.
This covers repeated roots whether or not their multiplicities
have gcd one.

Combining with QVM1 gives the complete conclusion for (5):
\[
             \deg H\ge2\quad\Longrightarrow\quad
       \text{no finite algebraic-coefficient Weyl mate for }S.
                                                               \tag{8}
\]
If H is squarefree, QVM1 applies. If it is not squarefree,
choose a repeated root and use (6)--(7). These two cases
exhaust every polynomial H of degree at least two.

The boundedness requirement on S0 is precise. QVM1's primitive
case permits coefficient exponents strictly below one. The new
proper-power conclusion assumes exponents at most zero.
It does not silently cover arbitrary fractional positive growth.

For the natural Laurent-polynomial scalar chart, the two statements
coincide at this boundary: below leading v degree one, all other
integer exponents are at most zero. Hence proper powers are no
longer an escape for that entire v-linear Laurent class.

## A finite Newton support test beyond bounded lower terms

There is a useful sufficient criterion for more general finite
Laurent-polynomial operators. After a constant shift b=a-alpha,
write S=sum c_(j,l)*v^j*b^l. Choose rational rho with 0<rho<d.
If every nonzero support point satisfies
\[
 j-\rho l\le0,\qquad
           F(E)=\sum_{j=\rho l}c_{j,l}E^l,\qquad \deg F\ge2,
                                                               \tag{9}
\]
then E=v^rho*b and u=v^(d-rho)/(d-rho) transport S to (1).
The coefficient exponent of each contraction drops further by
d-rho>0, so its limiting polynomial is exactly F(E).
No finite mate is possible.

This test can include positive intermediate v powers absent from
(5), when their momentum vanishing order is sufficient. It needs
a finite support certificate, not a partner-order search.

At rho=d the change has Euler derivation rather than an
algebraic ordinary coordinate. Under the same support inequality,
the limiting polynomial is
\[
 F_{\rm Euler}(E)=
    \sum_{j=dl}c_{j,l}
                    \prod_{t=0}^{l-1}(E-thd).                \tag{10}
\]
If this is nonconstant, (3) excludes every finite mate.
The product records all ordering corrections; replacing it by
E^l would omit terms at the same infinity degree.

For the quadratic chart d=2, (10) and the repeated-root
normalization also give a second complete bounded-lower class:
S=v^2*H(a)+S0 has no finite mate whenever H is nonconstant.
A simple root uses E=v^2*(a-alpha) and the Euler limit.
A root of multiplicity k>=2 uses E=v^(2/k)*(a-alpha),
an ordinary coordinate tending to infinity, and a degree-k limit.

This quadratic-growth corollary retains the bounded-coefficient
hypothesis on S0. It does not cover every higher-v-degree
deformation of Long: its existing 2v term would be an unbounded
lower term unless additional root-vanishing conditions such as
(9) were actually verified.

## The proper-power Long repair is excluded in every degree

A natural replacement for QVM1's rejected cubic coefficient is
\[
             S_k=2v(1+ha)^k+4a-3ra^2,\qquad k\ge2 .         \tag{11}
\]
It has the correct classical Long S and a proper-power leading
polynomial. The preceding primitive-root test alone did not
exclude it.

Its root alpha=-1/h has multiplicity k. Formula (6) gives
\[
             (S_k)_\infty=
                2h^k E^k-\frac4h-\frac{3r}{h^2}.             \tag{12}
\]
Thus no S_k has a finite Weyl mate, on any finite algebraic
coefficient cover, in any partner order. Arbitrary finite
corrections with bounded v coefficients are still excluded.

For k=3 the spectral connection (S_3-s) also passes the earlier
minimum-rank and branch-splitting tests. In the basis 1,a,a^2 its
ordinary v matrix is polynomial of degree one, because the
leading a coefficient is 2h^3*v and multiplication by v/h clears
its v denominator. QKR1 therefore gives pushforward rank six
and inverse-Fourier rank three, with the usual algebraic branch
projectors. The missing full Weyl action is still impossible by
(12). There is no need to compute more quantum corrections
inside this proper-power family.

## Positive controls and the remaining boundary

For d=2, S=v*(lambda*a+mu) and T=v/lambda satisfy [S,T]=h.
Also S=c*v^2+G(a), T=-a/(2c) is an exact pair for c nonzero.
They preserve the degree-one and constant-H boundaries of
(8) and its quadratic-growth corollary.

A different infinity face can support a genuine proper-power pair.
For d=2 put
\[
 E=v^3a,\qquad u=-v^{-1},\qquad
 P=E^2+u=v^6a^2+3hv^4a-v^{-1},\qquad Q=-E .
                                                               \tag{13}
\]
Then [E,u]=h and [P,Q]=h exactly. Here u tends to a finite point
when v tends to infinity, and the lower v^4*a term is unbounded.
Neither the transported ordinary infinity hypothesis nor (5)
applies. This prevents an unjustified claim about all proper-power
leading faces.

Formal mates also remain possible. With [E,u]=h,
\[
 P_h=E+hE^2,\qquad
       Q_h=u(1+2hE)^{-1}
                    =\sum_{n\ge0}(-2h)^n uE^n
                                                               \tag{14}
\]
satisfy [P_h,Q_h]=h formally. The mate has unbounded momentum
order and lies outside the finite Ore ring. Thus this infinity
obstruction is compatible with QFL1's formal existence theorem.

## A surviving linear growth Long model costs at least rank four

For a rational Laurent-coefficient operator in (5), let m be its
generic a order and p=deg H. If m>p, the coefficient of a^m is
bounded at infinity, while some coefficient of order at most p
has v degree one. The commutative Newton polynomial therefore has
a root with positive v exponent: otherwise its elementary symmetric
functions, after division by the leading coefficient, would all
be bounded.

For such a growing a branch, a=(h/v)*partial_v gives a v
exponential slope strictly greater than two. On the quadratic
q-line its slope is strictly greater than one. This is a formal
Newton-polygon assertion: derivative contractions are lower at
every positive a exponent, so they do not remove that upper edge.
The strictly positive infinity contribution rho_infinity is an
integer and is at least one.

The quadratic pushforward of (S-s) has rank2m and a Kummer
self-twist. QKR1 gives its q=0 contribution at least m,
whether or not the v=0 connection is nonsingular. Hence
\[
                       r(E)\ge m+1 .                        \tag{15}
\]
Any other finite singularities add nonnegative cost.

For a bounded full scalar Long lift with exact undeformed R,
[QSO2](QUADRATIC_OPERATOR_MATE_CLASSIFICATION.md) gives m>=3.
Equation (8) forces p<=1 in the v-linear bounded-lower class.
Thus m>p and (15) forces r(E)>=4. This closes the rank-three
proposal throughout that class; it does not exclude the remaining
higher-rank models or change the scalar-reduction hypothesis.

For example S=v*(2+h*a)+4a-3ra^2+h^2*a^3 has m=3 and p=1.
It has one bounded momentum branch and two growing like
v^(1/2). Its v connection is nonsingular at every finite point;
the q infinity contribution is one and its inverse rank is
exactly four. No full mate for this example is supplied.

## A commuting Darboux algebra still cannot supply a Weyl generator

The bounded-coefficient setting also occurs in bispectral research;
compare [Horozov, Section2](https://arxiv.org/html/math-ph/0206029).
Bispectrality and a finite Weyl mate impose different requirements,
and no classification theorem from that paper is used here.

For an explicit rational Darboux example, with [E,u]=h, set
\[
 \begin{split}
 A&=E-h/u,\\
 L_2&=E^2-2h^2/u^2,\\
 L_3&=E^3-3h^2u^{-2}E+3h^3u^{-3}.
 \end{split}                                                  \tag{16}
\]
Direct multiplication gives L2*A=A*E^2 and L3*A=A*E^3.
Cancellation in the Ore domain gives [L2,L3]=0 and L3^2=L2^3.
Every nonconstant element of C[L2,L3] has a nonlinear ordinary
infinity limit in C[E^2,E^3], so (1)--(2) exclude every finite
Weyl mate, of every order and on every finite coefficient cover.
The limit map is injective: in the normal form p(L2)+L3*q(L2),
the even and odd highest E orders cannot cancel.

This rules out choosing a scalar Weyl generator inside this whole
commuting Darboux algebra, not only its quadratic member. It does
not exclude a full action using other noncommuting operators in a
larger idealizer, or identify this algebra with every projective
Long carrier.

The remaining scalar constructions lie outside the excluded
v-linear Laurent class or have leading H of degree at most one
and satisfy this stronger rank cost, together with the full
right action and global admission.
Arbitrary nonpositional first images and projective or matrix
carriers remain separate possibilities. No all-representative
Long nonlifting is inferred.

For the leading-bracket framework compare
[Han--Tan, Theorem2.8](https://arxiv.org/html/2210.00257) and
[Guccione--Guccione--Valqui, Remark1.13 and Section2](https://arxiv.org/html/1111.6100).
The infinity limits, transport and whole-curve assertions above
are proved directly; no classification of bispectral operators
or centralizers is imported.

[Exact controls](check_infinity_limits.py) verify the finite
canonical changes, complete transformed operators and limits,
a proper-power family, a repeated-root polynomial with two roots,
ordinary and Euler Newton examples, the cubic spectral matrix,
positive finite pairs, a commuting Darboux algebra and a formal-mate control. The universal
all-degree statements use the written Puiseux and Ore proofs
and the registered QVM1 theorem. No independent theorem replay,
formal verification, external review, novelty priority, actual
finite-left Long carrier, bounded full Long lift, scalar or
matrix counterexample, all-representative nonlifting or DC2
solution is claimed.
