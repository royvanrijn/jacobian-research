# A uniform quadratic chart for nonreciprocal paired critical algebras

Every good paired-critical odd Laurent pencil has a canonical quadratic
presentation over its squared-critical-value algebra. Its paired product
is a square in that algebra, even when it is not a constant. This yields
exact rational coefficient equations for every degree, without choosing
critical roots or imposing reciprocal critical pairing.

The chart below is an implicit algebraic locus in rational affine
coordinates. It is not a rational parametrization of its solutions and
does not establish a rational degree-12 or degree-26 point. The existing
degree-26 p-adic construction is used only as a scoped consistency control.

## Two numerator parities in one notation

Let k be a field of characteristic zero, let h>=1 and epsilon in {0,1},
and put

\[
 b=2h+\epsilon,\quad n=2b,\quad L=2b-1,\qquad
                         f(z)=N(z^2)/z.
\]

The numerator N is monic of degree b. The good paired locus has

\[
             2tN'(t)-N(t)=L(t-1)^\epsilon S(t),               \tag{1}
\]

where S is monic squarefree of degree 2h, S(0)!=0, N is nonzero at
every S-root, and N(t)^2/t takes h distinct nonzero values there,
each twice.

For epsilon=1 additionally require N(1)=0 and S(1)!=0.
Equation (1) then gives N'(1)=0 and 2N''(1)=L S(1)!=0, so the zero
at one has exactly order two. Its two critical points provide the
zero branch level. For epsilon=0, all branch levels are nonzero.
Thus both cases give n/2 finite branch levels with two simple
critical points at each.

This concerns the even-numerator locus H(z)+a=N(z^2). It does not
assert that arbitrary Laurent pencils admit this odd symmetry.

## Theorem 1. The canonical quadratic algebra

Given such N, put

\[
 A=k[t]/S,\quad u=t\bmod S,\quad \nu=N(u),\quad
                   M=\nu^2/u,\qquad B=k[M]\subset A.
\]

Then B is etale of rank h and A is free of rank two over B.
There are unique units a,rho in B such that

\[
 \boxed{
 A=B[u]/\left(u^2-\left(2\rho+\frac{M}{a^2}\right)u+\rho^2\right),
 \qquad \nu=a(u-\rho).}
\tag{2}
\]

In particular

\[
 \boxed{\operatorname{Norm}_{A/B}(u)=\rho^2,\qquad
 \rho=-\operatorname{Norm}_{A/B}(\nu)/M.}
\tag{3}
\]

The discriminant of the quadratic equation has square class
M(M+4rho a^2) in B; this element is a unit.

**Proof.** Over a separable closure, B consists of functions on the h
distinct squared critical values, and A consists of functions on their
2h distinct critical roots. Each B point has precisely two A points.
This description descends the pair exchange to an involution tau of A,
with B=A^tau.

The element u-tau(u) is a unit, because different critical roots are
distinct. Also nu-tau(nu) is a unit. Indeed, equality of its two values
on any pair, together with the nonzero relation nu^2/u=M, would force
u=tau(u) on that pair. Therefore

\[
                   a=\frac{\nu-\tau(\nu)}{u-\tau(u)}
\]

is an invariant unit. The elements 1,u are a B basis: on each pair
their evaluation determinant is the nonzero difference of the two
u values. Write nu=a u+c with c in B. Multiplying nu^2=M u by its
conjugate shows that u tau(u) is a square in B.

More precisely, put rho=-nu tau(nu)/M, an invariant unit. From
tau(nu)^2=M tau(u) one obtains

\[
 a(u-\rho)
 =\frac{\nu-\tau(\nu)}{u-\tau(u)}
       \left(u+\frac{\nu\tau(\nu)}M\right)=\nu.
\]

Consequently c=-a rho. Reduction of nu^2=M u in the basis 1,u gives

\[
       u\tau(u)=\rho^2,\qquad u+\tau(u)=2\rho+M/a^2.
\]

These identities prove (2)-(3) and uniqueness. The quadratic discriminant
is

\[
                  (u-\tau(u))^2
                    =\frac{M(M+4\rho a^2)}{a^4};
\]

it is a unit because the critical roots are distinct. ∎

The units a and rho are elements of B, not necessarily scalars in k.
Thus (2) permits a different paired product at each squared critical
value. Requiring rho^2 to lie in k recovers the reciprocal ansatz.

## Theorem 2. Complete rational coefficient equations

Choose a monic polynomial q(Y) of degree h with q(0)Disc(q)!=0.
Let B=k[Y]/q, write w=Y mod q, and choose elements a,rho in B.
Impose the unit conditions

\[
       a,\ \rho,\ \Delta=w(w+4\rho a^2)\quad\hbox{are units}.    \tag{4}
\]

Define

\[
 A=B[u]/(u^2-\theta u+\rho^2),\quad
             \theta=2\rho+w/a^2,\qquad \nu=a(u-\rho).
\tag{5}
\]

Let S(T) be the characteristic polynomial of multiplication by u on
the 2h-dimensional k algebra A. Impose

\[
            \operatorname{Disc}(S)\ne0,\qquad
                 S(1)\ne0\quad\hbox{if epsilon=1}.             \tag{6}
\]

Define N, with no unknown integration constants, by

\[
 R(T)=L(T-1)^\epsilon S(T)=\sum_j R_jT^j,\qquad
                         N(T)=\sum_j\frac{R_j}{2j-1}T^j.
\tag{7}
\]

The exact remaining equations are

\[
                       N(u)=\nu\quad\hbox{in }A,
       \qquad N(1)=0\quad\hbox{if epsilon=1}.                  \tag{8}
\]

These are 2h+epsilon scalar equations in 3h rational affine
coordinates: the h coefficients of q and the h coordinates of each
of a,rho in the power basis of B. Together with (4)-(6), they are
equivalent to every good paired model in the stated normalization.
All constructions and their inverse are rational over k.

**Proof.** Conditions (4) make A etale of rank two over B. Also u is a
unit since its norm is rho^2. Direct calculation in (5) gives

\[
                              \nu^2=u w.                     \tag{9}
\]

In particular nu is a unit. The discriminant condition in (6) says
that all 2h geometric u values are distinct, so u generates A and
A=k[u]. Formula (7) is monic of degree b and has
2T N'-N=R. Equation (8) identifies its critical value element N(u)^2/u
with w. The squarefree polynomial q supplies precisely h different
nonzero values, with two critical roots above each. For epsilon=1
the extra equation and S(1)!=0 give exactly the double zero and
simple distinguished critical root already stated. This proves the
good paired passport and all its required nondegeneracies.

Conversely apply Theorem 1 to a good N. Let q be the minimal polynomial
of M on B, equivalently the unique monic square root of its
characteristic polynomial on A. Then q has degree h and
q(0)Disc(q)!=0. Theorem 1 provides a,rho and all unit conditions.
The reconstructed characteristic polynomial S is the original S,
and the invertibility of 2T d/dT-1 on polynomials proves that (7)
recovers the original N. Hence (8) holds.

There is an explicit rational inverse without labeling roots. In A
form the ordered basis matrix with columns

\[
 1,M,\ldots,M^{h-1},\ u,uM,\ldots,uM^{h-1}.                     \tag{10}
\]

It is invertible by the pairwise determinant argument above.
Expressing M^h in this basis recovers q; expressing nu recovers
c(M)+a(M)u. Then rho=-c/a is computed by inversion in B.
Characteristic polynomials, this basis solve and multiplication
inversions are rational operations with nonzero determinants on
the indicated open set. All divisions in (7) are by fixed nonzero
integers. Thus no root choices or additional descent conditions
are hidden. ∎

For degree twelve, h=3, epsilon=0: this is a nine-variable system
with six equations. For normalized degree twenty-six, h=6,
epsilon=1: it is an eighteen-variable system with thirteen equations.
These equation counts do not prove finite solution sets or rationality.
They describe the full paired even-H locus, including nonreciprocal
critical pairings.

## Corollary 3. A uniform arithmetic norm condition

Every good model above satisfies

\[
 S(0)=\operatorname{Norm}_{B/k}(\rho)^2,\qquad
 \boxed{N(0)=(-1)^{b+1}(n-1)\operatorname{Norm}_{B/k}(\rho)^2.} \tag{11}
\]

Indeed S has even degree, so its constant is Norm_A/k(u).
Apply the norm tower to (3), then evaluate (1) at zero.

Over Q this gives a necessary squareclass condition on the pole
coefficient in both numerator parities. In degree twelve,
-N(0)/11 must be a rational square. In the admissible degree-26
case N(0) must be a rational square. This is not a sufficiency theorem.

Moreover the residual critical-point polynomial S(z^2) has literal
square discriminant:

\[
 \operatorname{Disc}_z(S(z^2))
 =4^{2h}S(0)\operatorname{Disc}(S)^2
 =\left(2^{2h}\operatorname{Norm}_{B/k}(\rho)
                           \operatorname{Disc}(S)\right)^2.   \tag{12}
\]

The first equality follows by multiplying the derivatives
2zS'(z^2) at its roots. This concerns the critical-point polynomial,
not the discriminant of the pencil N(z^2)-sz.

If the paired product is a nonsquare scalar lambda in k, (3)
places a square root of lambda in B. In particular, if h is odd
and k=Q, such a scalar cannot be nonsquare: every component field
of B would have even degree. The earlier reciprocal parity
obstruction is therefore visible directly in the critical-value algebra.

## Corollary 4. A rational matrix test for the entire real locus

For real chart data satisfying Theorem 2, define the symmetric h-by-h
trace matrices

\[
 (H_g)_{ij}=\operatorname{Tr}_{B/\mathbb R}(g w^{i+j}),
                           \qquad 0\le i,j<h.
\tag{13}
\]

Every finite critical point of f is real if and only if both

\[
                         H_w>0,\qquad H_{w+4\rho a^2}>0
\tag{14}
\]

are positive definite. For rational chart data, this is equivalent to
positivity of their 2h rational leading principal minors. Thus the
full all-real locus in this chart has an exact test requiring no
algebraic root choices.

**Proof.** Write the real etale algebra as a product of real and complex
fields. On a complex component, for g=alpha+i beta != 0 the quadratic
trace form at z=x+i y is

\[
 \operatorname{Tr}_{\mathbb C/\mathbb R}(gz^2)
       =2\alpha(x^2-y^2)-4\beta xy.
\]

Its determinant is -4(alpha^2+beta^2)<0. Therefore H_w positive
definite forces B to split as R^h; its h weights w_i must all be
positive. In this split algebra, the invertible Vandermonde evaluation
matrix makes H_g congruent to diag(g_1,...,g_h). Thus (14) is exactly

\[
                    w_i>0,\qquad w_i+4\rho_i a_i^2>0
                    \quad\hbox{for every }i.
\]

The quadratic u discriminant is w_i(w_i+4rho_i a_i^2)/a_i^4, so it is
positive. Also rho_i>-w_i/(4a_i^2), and hence its trace
2rho_i+w_i/a_i^2>w_i/(2a_i^2)>0. Its norm rho_i^2 is positive.
Both u roots are consequently positive, so the corresponding critical
points z=+-sqrt(u) are all real. The distinguished roots +-1 when
epsilon=1 are already real.

Conversely, positive critical u roots make their squared values
w_i=N(u_i)^2/u_i positive, so B splits over R. Each quadratic pair
has positive discriminant, giving the second strict inequality.
This proves (14). The leading-minor test is Sylvester's criterion,
obtained by successive completion of squares for a real symmetric
matrix. ∎

## Exact controls and remaining arithmetic problem

The [checker](verify_paired_critical_trace_chart.py) recovers the chart
and reconstructs the entire model for the existing rational degree-eight
and degree-ten seeds, including positive and negative exact controls for
the Hermite real-locus criterion. The degree-eight paired product is nonconstant;
the degree-ten data recover rho=2/15 exactly.

For the retained unrestricted degree-26 lift modulo nine, the same
rational algebra identities are valid because all matrix and element
inverses used are units. The reconstructed paired product is

\[
                         \rho^2=1+3w^2+3w^3\pmod9.
\]

It is nonconstant and so records the departure from constant reciprocity
in the first lift. This controls a finite shadow of the
[existing p-adic family](DEGREE_TWENTY_SIX_THREE_ADIC_PENCILS.md).
The checker does not divide nonunit Euler constants modulo nine and
does not infer rational descent from local lifting.

The [receipt](paired_critical_trace_chart_verification.json) records
all q,a,rho data and exact roundtrips; the
[source ledger](paired_critical_trace_chart_source_ledger.json) pins
the proof inputs. The algebraic equivalence has received internal
critical review. No uniform rational point theorem, rational
parametrization of the solution locus, new degree-26 rational model,
formal verification, independent external replay or external human
review is claimed.
