# Two traces close arbitrary plane boundary curves

14 November 2026. The single-trace Hodge test has a useful global
consequence: one pure trace forces every parallel trace of a hypothetical
nonautomorphic unit-Hessian polynomial to have a nonpositive pure
signature. This closes a new all-degree y-graph family with arbitrary
plane boundary curves, including positive genus, singularities and
multiplicities. A second corollary closes cuspidal boundaries from one
trace alone.

The only external mathematical inputs are the retained Tang alternative,
the retained affine-component theorem and the retained basic
Hodge--Deligne curve formula. Neither plane Keller endpoint is resolved.

## 1. A pure trace fixes the whole parallel profile (HC4MYHE0)

Work first over C, with actual affine source coordinates (s,z1,z2,z3).
Put t=uv. Suppose f has nonzero constant four-variable Hessian and its
gradient is not an automorphism. If at any scalar s=a the generic value
fibre has

\[
                    E(f(a,-)=\lambda)=t^2+\kappa t,          \tag{1}
\]

then, at **every** scalar s=b, its generic value fibre has

\[
                    E(f(b,-)=\lambda)=t^2-m_b t,
                       \qquad m_b\in\mathbb Z_{\ge0}.       \tag{2}
\]

The integers m_b vanish for all but finitely many b. The full relative
critical scheme is a finite disjoint union of reduced affine lines,
each contained in one scalar s hyperplane. In particular, no trace has
an isolated critical point.

Indeed, adapt the constant jointly affine two-plane supplied by
[HCTANG2](../tang-hc4-audit/STRUCTURAL_REVIEW.md) to this same base s.
In rank one, [HC4MYHC0](HODGE_TRACE_CLOSURE.md#2-a-single-trace-criterion-with-no-critical-point-hypothesis-hc4myhc0)
gives, with one fixed smooth plane Keller curve C=B=0,

\[
 E(f(b,-)=\lambda)=t^2+(d_b-c)t+G(u+v)+(R-c).              \tag{3}
\]

Here c,G,R count its components, total completion genus and total
punctures; d_b is the generic degree of D+bA on its nonconstant
components. Constant components cause no problem: avoid their finitely
many values when choosing lambda. Thus G and R-c are independent of b.
Equation (1) forces G=0 and R=c. Every component is an actual closed
A1, and HC4MYGR0 makes the plane Keller pair, hence the full gradient,
an automorphism. This contradicts the hypothesis.

The remaining rank-zero form is

\[
 f=D(s,U)+P F(s,U)+Q G(s,U),\qquad J(F,G)\in C^\times.      \tag{4}
\]

At b, the two specialized rows are not both zero. Complete their
primitive row over C[U], retaining the gcd g_b. The trace is polynomially
equivalent to D_b(U)+g_b(U)V, with W free, so its generic fibre is
(A1 minus m_b distinct roots) times A1. This proves (2).

The common zeros of the plane Keller pair (F,G) are finite and reduced.
At each such (s,U), F_U,G_U are not both zero and the last relative
critical equation is one nonconstant affine-linear equation in P,Q.
It gives one reduced A1. There are no other relative critical points.
This proves the asserted whole-scheme description and finite support
of m_b. It neither makes the plane pair proper nor identifies its
coefficient ideal from one trace.

Consequently a single pure trace together with **any** parallel trace
whose generic Hodge polynomial is outside the set (2) forces gradient
invertibility under unit Hessian. Two traces whose genus or constant
Hodge coefficients differ also suffice, by (3) and the rank-zero formula.
All statements concern whole value fibres in the same actual affine
fibration. No nonlinear straightening of source divisors is used.

## 2. A linear trace cannot have a nonpositive pure signature (HC4MYHE1)

Let

\[
              h(U,V,W)=H(U)+V a(U,W),\qquad
              \deg H=D\ge2,\quad \deg_W a\ge1.             \tag{5}
\]

For generic lambda, its whole value fibre **never** has Hodge polynomial
t^2-mt with m>=0. Neither smoothness, squarefreeness, monicity,
properness of U on the zero curve, nor a genus bound is assumed.

Let C=V(a)_red in the (U,W) plane, and let d be the summed generic
degree of H|C over its nonconstant components. Off C solve V. Over the
d distinct points H=lambda on C, V is free. Constant components are
avoided. These are actual product strata, so

\[
                         E(h=\lambda)=t^2-E(C)+d t.         \tag{6}
\]

Normalize the c irreducible components of C. Write G for their total
completion genus and R for their total punctures. For every singular
point x of the reduced curve let b_x be its number of normalization
preimages, and set e=sum_x(b_x-1)>=0. Additivity, applied outside the
finite singular set, gives

\[
 E(C)=ct-G(u+v)+c-R-e.                                    \tag{7}
\]

This formula counts branches rather than conductor lengths. In
particular a cusp contributes zero to e; a node or an intersection of
different components contributes positively. Repeated factors of a do
not change C or (6).

Suppose (6) were t^2-mt. Equations (6),(7) give

\[
                G=0,\qquad R+e=c,\qquad d+m=c.              \tag{8}
\]

Each normalization is affine and has at least one puncture. Hence all
normalizations are A1, e=0, and distinct components are disjoint. The
positive W degree of a supplies at least one nonvertical component
C_i: its equation depends on W. On its A1 normalization the function U
is a nonconstant polynomial and is therefore surjective onto A1.

There cannot be a vertical component U=rho: surjectivity would make
C_i intersect that line, contradicting e=0. Thus every component is
nonvertical. If e_i>=1 is the degree of U on its normalization, then
the generic affine degree of H|C_i is D e_i. Normalization changes
only finitely many points, whose H values are avoided. Therefore

\[
                      d=D\sum_i e_i\ge2c,                  \tag{9}
\]

contradicting d+m=c. This proves the claim in all degrees.

Combining it with HC4MYHE0: if one trace of a four-variable polynomial
is polynomially equivalent to (5), and one parallel trace has a pure
signature (1), **unit Hessian forces a polynomial gradient inverse**.
Each trace equivalence must have an actual polynomial inverse, but the
two equivalences need not extend to one global coordinate change.

## 3. Arbitrary q-independent r boundaries with one constant slice (HC4MYHE2)

Use the scaled Meng--Yang potential Psi=L A^2+M A+N B, L*N nonzero.
Take any actual graph y=R(x,p,q,r) satisfying

\[
 R(0,p,q,r)=Z(p,r),\quad \deg_r Z\ge1,\qquad
 R(a,p,q,r)=b\in K,\quad a\in K^\times.                  \tag{10}
\]

After **every** actual polynomial retained-three-dual automorphism over
K[x], unit full Hessian forces a polynomial gradient inverse.

At x=0 the trace is

\[
                  Lp^2+Mp+4NpZ^2+NqZ.
\]

The actual polynomial shear V=N(q+4pZ), U=p, W=r carries it to
H(U)+V Z(U,W), H=LU^2+MU. Thus HC4MYHE1 applies to this entire
trace, irrespective of the zero curve's genera, singularities or
multiplicities.

At x=a the physical base point is (a,b), with a nonzero. The coefficient
rows of the ambient linear forms A,B have cross product

\[
 \big(2a^2(2ab+3),\;-a(1+ab)(2+ab),\;b(1+ab)^3\big).       \tag{11}
\]

It is nonzero: the first entry can vanish only if ab=-3/2, when the
second entry is a/4. Complete A,B to a linear coordinate triple (A,B,C).
The further triangular change (A,B,C) to (A,H(A)+NB,C) is an actual
polynomial automorphism because N is a field unit. The trace is a
coordinate, and its generic value fibre is A2, of Hodge polynomial t^2.
HC4MYHE0 and HC4MYHE1 now contradict hypothetical nonautomorphy.

All off-slice terms are unrestricted. Equivalently, every polynomial

\[
 R=(1-x/a)Z(p,r)+(x/a)b+x(x-a)T(x,p,q,r)                   \tag{12}
\]

is covered for arbitrary T and arbitrary positive r degree Z. Conversely
(10) implies (12) by division by the two coprime scalar factors x,x-a.
This permits arbitrary old-q dependence away from the two slices, and
arbitrary graph, boundary and automorphism word degrees. Retained-dual
automorphisms preserve each whole trace fibre; the Hessian condition
is always imposed on the actual four-variable pullback. No nonlinear
Hessian invariance is assumed.

## 4. One-place unibranch boundaries need only one trace (HC4MYHE3)

In (5), suppose every irreducible component of C=V(a)_red has
normalization A1, the components are disjoint, and every singularity
is unibranch. Then E(C)=ct. The surjectivity argument above excludes
vertical components, and (6),(9) give

\[
                   E(h=\lambda)=t^2+\kappa t,
                      \qquad \kappa=d-c\ge c>0.             \tag{13}
\]

The single-trace HC4MYHC0 criterion applies. In particular every y graph
with boundary Z(p,r) satisfying these reduced zero-curve hypotheses
is closed after all retained-dual automorphisms, with **no second-slice
condition** and arbitrary off-x-zero terms. The zero curve need not be
smooth or a plane coordinate fibre. For example

\[
                    Z=(r^m-p^n)^k,
                \quad m,n\ge1,\;\gcd(m,n)=1,\;k\ge1        \tag{14}
\]

has normalization p=t^m,r=t^n, one puncture and one possible unibranch
singularity at the origin. Its generic trace polynomial is
t^2+(2m-1)t. The exponent k does not enter. When m,n>=2 this explicitly
includes singular cusps beyond the earlier smooth coordinate-boundary
class. Equation (13) asserts gradient invertibility under unit Hessian;
it does not assert existence or nonexistence of such members.

## Boundaries, descent and controls

One pure trace alone does not prove invertibility: the unresolved
relative-coordinate/unit-ideal endpoint has pure traces. Nor does a
general positive-genus boundary without the second trace satisfy the
positive single-trace criterion. A component with several punctures
is not relabeled A1. No assertion excluding C* components of arbitrary
plane Keller fibres is used.

Local mapping degrees at distinct points can be summed only within
one fixed value fibre; summing them across different critical values
would give a false bound. The present proof uses whole Hodge profiles
and normalization strata instead of that shortcut.

For any characteristic-zero field, the finitely many coefficients of
a candidate and its stipulated inverses embed into C. Polynomial
gradient invertibility descends by uniqueness. In HC4MYHE3 the geometric
curve hypotheses are taken over the algebraic closure.

[verify_two_trace.py](verify_two_trace.py) and its
[receipt](receipt-two-trace.json) check normalization/Hodge identities,
singular and multiple-component controls, trace shears, base-row
independence, interpolation, retained-dual inverses and source pins.
The arbitrary-degree normalization, surjectivity and Tang arguments
are written proofs, not finite-search conclusions. Full HC4, the
unit-ideal surjectivity endpoint and the distinct singleton-zero
properness endpoint remain open. No formal verification, independent
replay, external review or novelty priority is claimed.
