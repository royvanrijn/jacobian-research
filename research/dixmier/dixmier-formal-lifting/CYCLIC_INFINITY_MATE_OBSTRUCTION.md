# A leading term obstruction to finite Weyl mates on cyclic charts

14 November 2026. **QVM1.** A simple infinity valuation excludes
every finite-order Weyl mate for an entire higher-momentum class,
including an h-small cubic repair of Long's quadratic S.
The exclusion allows arbitrary finite algebraic coefficient covers,
arbitrary partner orders and poles, and arbitrary lower-growth
corrections to S.

Let q=v^d/d with d>=2, and let [a,f]=h*v^(1-d)*f_v.
If the leading v term of S at infinity is v*H(a), where H has
degree at least two and is not a proper polynomial power, then
S has no finite algebraic-coefficient T with nonzero constant
commutator. This covers S=v*(2+h*a^3)+4a-3ra^2, whose generic
order is three and whose classical reduction is Long's S.

This is a statement about the specified infinity face. A bounded
Long lift can still change that face, make H a proper power, or
leave the exact-R scalar category. The result is not an
all-representative nonlifting theorem.

## Valuation and the exact hypothesis

Let C be algebraically closed of characteristic zero, let h be
a nonzero differential constant, and let F/C(v) be any finite
algebraic function-field extension. Extend
\[
             \delta=v^{1-d}\partial_v,\qquad d\ge2
\]
uniquely to F, and work in the finite scalar Ore ring F[a;h*delta].
At a point over v=infinity, use Puiseux series in v^(-1/e).
Define deg_v(f) as its leading v exponent, and for an operator
sum f_j*a^j use the maximum of those coefficient exponents.
All such exponents are rational.

Assume S has
\[
 S=vH(a)+S_{<1},\qquad
 H\in C[a],\quad \deg H\ge2,\quad \deg_v S_{<1}<1 .            \tag{1}
\]
There is no order bound on S_<1. Its coefficients may have
arbitrary algebraic poles away from infinity. Write
H=c*product_i(a-alpha_i)^(n_i), with c nonzero, and require
\[
                            \gcd_i n_i=1.                   \tag{2}
\]
Thus H is not a proper power after a constant rescaling over C.
In particular every squarefree H of degree at least two qualifies.

**Theorem.** Under (1)--(2), no finite T in any such algebraic
coefficient field satisfies
\[
                             [S,T]=\kappa,\qquad
                           0\ne\kappa\in C .                \tag{3}
\]
The constant fields may first be extended algebraically; the
result descends to any characteristic-zero constant field.

## The leading commutator

Suppose T has leading term v^M*K(a), with M rational and K nonzero
in C[a]. Since delta decreases v exponent by d, the Ore product gives
\[
 \deg_v[S,T]\le M+1-d,\qquad
 [S,T]_{\rm leading}
   =h\,v^{M+1-d}\bigl(MH'K-HK'\bigr).                        \tag{4}
\]
The displayed coefficient is interpreted as the coefficient at
that exponent; it may vanish. All repeated contractions and all
terms from S_<1 or lower terms of T have strictly smaller exponent.
The bound and formula are valid on every finite Puiseux cover.

If M>d-1, equation (3) forces
\[
                         MH'K-HK'=0 .
\]
Comparing residues of the logarithmic derivatives K'/K and H'/H
at the roots of H gives
\[
             \operatorname{ord}_{\alpha_i}K=M n_i.
\]
There are no other roots of K, because K'/K=M H'/H.
Condition (2) forces M to be an integer, and hence
\[
                              K=c_0 H^M.                    \tag{5}
\]
Here M>d-1>=1, so S^M is a finite operator. Replace T by
T-c0*S^M. This preserves (3) and strictly lowers deg_v(T).
Every exponent still greater than d-1 must again be an integer
by the same argument. Thus this reduction terminates after finitely
many subtractions; no infinite subtraction or completion is used.

We may consequently assume M<=d-1. If M<d-1, (4) makes the
commutator have negative v degree, contradicting (3). Thus
M=d-1, and its leading coefficient must obey
\[
                     (d-1)H'K-HK'=\kappa/h .                \tag{6}
\]

## The polynomial leading equation has no solution

Put n=deg H>=2 and t=deg K. Unless t=(d-1)n, the leading term
on the left of (6) is nonzero and has degree n+t-1>=1.
It cannot be a nonzero constant.

If t=(d-1)n, subtract a constant multiple of H^(d-1) from K to
reduce its degree. This does not alter the left side of (6).
If the result is zero the left side is zero. Otherwise its degree
is now below (d-1)n, so the preceding nonconstant leading-term
argument applies. Both outcomes contradict (6), proving the theorem.

The proof gives an all-order criterion from one infinity face.
Neither coefficient-degree bounds nor an enumeration of finite
covers or partner orders enters it.

## Why the restrictions matter

For degree-one H=alpha*a+beta, with alpha nonzero, the exact pair
\[
       S=v(\alpha a+\beta),\qquad
       T=\frac{v^{d-1}}{\alpha(d-1)}
\]
satisfies [S,T]=h. Thus the degree threshold cannot be discarded.

For H=(a^2+1)^2, the homogeneous leading equation in (4) has the
fractional-exponent solution M=3/2, K=(a^2+1)^3.
Such a leading term cannot be removed by subtracting an integer
power of S. This demonstrates the boundary of the reduction proof;
it does not give a finite Weyl mate for this proper-power H.

The assumption deg_v S_<1<1 is also essential to the method.
Terms on a different top face can change both (4) and the finite
reduction. Ordinary polynomial coordinate equivalence or formal
Hamiltonian equivalence does not preserve this specified face.

For the classical leading-bracket and proportional-power framework, see
[Han--Tan, Theorem2.8](https://arxiv.org/html/2210.00257) and
[Guccione--Guccione--Valqui, Remark1.13 and Section2](https://arxiv.org/html/1111.6100).
Those polynomial and fractional-Weyl results are context; the proof above
establishes the stated whole-curve, arbitrary-Puiseux-cover version directly.

The [quadratic finite-mate theorem](QUADRATIC_OPERATOR_MATE_CLASSIFICATION.md)
gives the complementary low-order input.
This proof uses only the finite Ore calculus, polynomial roots and
characteristic-zero Puiseux expansions. No novelty priority or
independent theorem review is asserted.

## An h-small cubic candidate passes the Fourier rank test

The [QKR1 rank budget](FOURIER_CYCLIC_RANK_BUDGET.md) makes a
generic inverse-Fourier rank-two proposal insufficient for a
bounded full exact-R scalar Long lift. The next natural local
test is
\[
              S_*=v(2+h a^3)+4a-3ra^2,\qquad q=v^2/2 .
                                                               \tag{7}
\]
Its special operator symbol is the retained Long S, while its
generic differential order is three. It therefore escapes the
order-two hypothesis of QSO2.

This candidate also passes the algebraic branch and minimum-rank
tests. For a source eigenvalue s, the connection S_*-s has rank
three over C(v). In the basis 1,a,a^2 its ordinary v derivative
matrix is
\[
 A(v)=
 \begin{pmatrix}
 0&0&(s-2v)/h^2\\
 v/h&0&-4/h^2\\
 0&v/h&3r/h^2
 \end{pmatrix}.                                               \tag{8}
\]
It is polynomial of degree at most one, so the connection is
nonsingular at every finite v and has slopes at infinity at most
two. Its quadratic pushforward has rank six, local characters
1^3 direct_sum K_(1/2)^3 at q=0, no other finite singularities,
and slopes at infinity at most one. The Fourier rank formula
therefore gives exactly
\[
                       N=6,\qquad r(E_*)=3.                 \tag{9}
\]
It saturates the quadratic rank budget.

For completeness, write A(v)=A0+v*A1 and put
C0=h*A1, Cminus=h*A0. Its a-connection is
C(v)=C0+Cminus/v. The pushed-forward q matrix is
\[
 H_*=
 \begin{pmatrix}
 C0&Cminus\\
 Cminus/(2q)&C0+hI/(2q)
 \end{pmatrix}.
\]
With V=[[0,2qI],[I,0]], the two algebraic projectors
(I+V/v)/2 and (I-V/v)/2 are horizontal after q=v^2/2.
Their image connections are C0+Cminus/v and C0-Cminus/v.
Thus the rank-three inverse model genuinely permits extraction;
a missing algebraic projector is not its failure.

However H(a)=2+h*a^3 in (7) is squarefree over the generic constant
field: H'=3h*a^2 and H(0)=2. It has degree three.
Equations (1)--(6) exclude every finite T satisfying [S_*,T]=h,
over all finite algebraic coefficient covers. The exclusion
continues to hold after arbitrary finite lower-v-degree corrections
to S_*. It is independent of the inverse-Fourier rank computation.

This is an exact failed candidate and an all-degree failed family.
It is not a full finite-left Long graph carrier. No fourth-image
action or global finite-left lattice is inferred from (8)--(9).

## Current scope

The leading-face obstruction above is independent of the inverse-Fourier
rank calculation in the explicit control. The
[proper-power infinity limits](ALGEBRAIC_INFINITY_LIMITS.md) and
[zero-mode transport](QUANTUM_ZERO_MODE_TRANSPORT.md) supply the remaining
cases used in the completed scalar classification. Other coefficient faces,
deformed first images and nonscalar carriers remain outside this theorem.
No formal verification or independent external review is claimed.
