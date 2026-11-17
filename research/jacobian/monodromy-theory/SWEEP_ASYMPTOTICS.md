# Asymptotic directions and rational multiplier constraints

Work over algebraically closed constants of characteristic zero.
The hypotheses and conclusions concern normalized sweep constructions;
they do not classify arbitrary Keller maps.

**KMFI32.** Let X(w,z) be a polynomial three-vector and

\[
 S(\gamma,w,z)=X(w,z)+\gamma(1,w,z),\qquad \det JS=\gamma^2.
\]

If X is generically of rank two, its generic degree delta onto its
image is at most two. In particular, if S is finite, its affine
inertia has type 3^delta 1^(N-3delta), with delta=1 or2.
The exceptional degree-eleven PSL2(11) type 3^3 1^2 is impossible.
The natural degree-twelve action is also impossible for an actual
Keller target twist, since its order-three elements have no fixed sheet.

This closes the finite tautological four-dimensional Gao route
without any polynomial or rational multiplier-chart assumption.
The same result holds after a polynomial automorphism in the two
direction coordinates. A constant-Jacobian direction map which has
not been proved invertible does not supply that coordinate change.

**KMFI33.** Suppose a nonzero polynomial lambda and a rational
multiplier rho in k(w)[z] give a polynomial nilpotent matrix

\[
 P=\lambda\begin{pmatrix}-\rho&1\\-\rho^2&\rho\end{pmatrix}
\]

satisfying the integrability equation

\[
 \lambda_{ww}+2(\rho\lambda)_{wz}+(\rho^2\lambda)_{zz}=0. \tag{1}
\]

Then either rho belongs to k[w], or

\[
 \rho={\kappa\over w-a}\,z+q(w),\qquad
 \lambda_{\rm top}=C(w-a)^r z^s,\quad r>=2,\ s>=0,
 \tag{2}
\]
\[
 r(r-1)+2\kappa(r-1)(s+1)+\kappa^2(s+1)(s+2)=0. \tag{3}
\]

Here q has no pole away from a, and its pole at a has order at most
floor(r/2); its polynomial part is not classified by this statement.
In particular every polynomial multiplier rho(w,z) is actually
independent of z. A rational multiplier depending only on w
must be polynomial. Together with KMFI28, every primitive finite
sweep with a nonconstant polynomial multiplier is alternating.

**KMFI34.** For the homogeneous rational class rho=kappa z/w with kappa!=0,
equation (1) has either no nonzero polynomial solution with
polynomial matrix P, or its entire solution space is
lambda=Cw^r z^s for one pair r>=2,s>=0 satisfying (3).
The resulting sweep field pair is an independently stabilized
polynomial curve cover of degree N=r+s+1, with geometric group
A_N for odd N and S_N for even N.
For N=3 this means the regular C3 action and no absolute Keller map.
This complete homogeneous class cannot realize either PSL2(11) action.
The centered, zero-integration-constant sweep has no origin preimage
with gamma!=0, so its declared positive-weight Gao stage cannot exist.
This last origin statement is about that centered architecture;
translations or different affine presentations are separate.

The [checker](verify_sweep_asymptotics.py),
[model](sweep_asymptotics_model.json) and
[receipt](sweep_asymptotics_receipt.json) give exact controls.
Universal arguments are written below.
The full action classification and the abstract PSL2(11) cases remain open.

## 1. Nilpotency and an asymptotic tangent direction

Put Delta=(1,w,z), N=[[-w,1,0],[-z,0,1]], and P=N JX.
Since N JDelta=I, the Jacobian condition is

\[
 \det(P+\gamma I)=\gamma^2.
\]

Thus P has trace and determinant zero.
At a general point where JX has rank two, P has rank one:
if P=0 both columns of JX would be proportional to Delta.
For this rank-one nilpotent matrix, im P=ker P.

Tangency gives Delta=JX v for a unique nonzero source vector v.
Applying N shows Pv=0. Hence v belongs to im P.
Let ell be a left kernel vector of P and n=ell N a normal to X.
Then ell v=0 and n JDelta v=ell v=0.
Differentiate n dot Delta=0 in direction v:

\[
 (D_vn)\cdot\Delta+n\cdot J\Delta\,v=0,
 \quad\text{so}\quad (D_vn)\cdot JXv=0. \tag{4}
\]

This says the tangent direction Delta is a zero of the second
fundamental form of the image surface.

## 2. At most two preimages on the branch surface

The second fundamental form is a homogeneous quadratic on the
two-dimensional tangent space of the image surface D.
It is not identically zero at the generic point.
Otherwise the differential of the projective Gauss normal would
vanish. In characteristic zero that normal would be constant,
making D a plane.
A constant nonzero normal cannot annihilate (1,w,z) for all
independent source parameters w,z, so this is impossible.

A nonzero quadratic on a projective line has at most two zeros.
Over a general smooth point of D, distinct preimages under X
give distinct projective Delta-directions: their first coordinate
is one, so equality of directions would give the same w,z.
Equation (4) therefore proves delta<=2.

If S is finite, X maps gamma0 finitely onto a surface.
Thus rank two is automatic, and the full normalization is present.
KMFI27 makes the ramification index three and counts exactly delta
ramified sheets. Its inertia is 3^delta times fixed sheets.
This yields the PSL2(11) conclusions in KMFI32.

The geometric degree bound itself needs only rank two of X.
The full inertia conclusion needs the finite-normalization hypothesis:
missing branches of a nonfinite sweep are not controlled by (4).
The argument does not turn a bound for this sweep into a dimension
bound for an arbitrary Keller self-map.

The signed family gives a sharp control.
For X(q,sigma)=(sigma,b(q),d(q)-sigma^2/2), b'=q a',d'=a b',
a normal is (sigma,-a,1), and the second fundamental form is

\[
 d\sigma^2-q(a')^2dq^2.
\]

Its two directions correspond to w=+sqrt(q) and -sqrt(q).
The univariate-polynomial class studied in KMFI28 is a special case;
no multiplier chart was used in the proof above.

## 3. Polynomial multipliers cannot depend on the second parameter

Let s=deg_z lambda and f(w) its leading coefficient.
First suppose rho has z-degree d>=2, allowing rational coefficients
in w. The term (rho^2 lambda)_(zz) has z-degree s+2d-2,
strictly larger than the other two terms in (1).
Its coefficient is nonzero, a contradiction.

If rho=p(w)z+q(w), its highest coefficient satisfies

\[
 f''+2(s+1)(pf)' +(s+1)(s+2)p^2 f=0. \tag{5}
\]

For a nonzero polynomial p, the p^2 f term has strictly larger
w-degree than the derivative terms. Thus p cannot be polynomial.
This proves that a polynomial rho must be independent of z.

If rho=h(w) is rational and independent of z, the top coefficient
in (1) gives f''=0.
Write h in reduced form. Polynomiality of h^2 lambda forces
the squared denominator of h to divide f.
A nonconstant denominator is incompatible with f being linear
or constant. Hence h is polynomial.

## 4. The single-pole theorem for a rational z coefficient

Now take rho=p(w)z+q(w), p!=0 rational.
Polynomiality of P forces the squared reduced denominator of rho
to divide lambda, hence each coefficient, including f.

Equation (5) first shows that p has neither a nonzero polynomial
part nor order at least two of decay at infinity.
For a polynomial part, p^2 f has dominant degree.
If p decays faster than 1/w, f'' has dominant degree; f has degree
at least two because a denominator square divides it.
Therefore p behaves as kappa_infinity/w at infinity.

Every finite pole of p is simple.
At a pole of order t>1, write f of order r>=2t there.
The p^2 f term has order r-2t, strictly less than the derivative
terms r-t-1 and r-2, and cannot cancel.
At each simple pole a_i of residue kappa_i, let r_i>=2 be
the zero order of f. Equation (5) gives (3), with r=r_i.
At infinity the same equation gives (3) with r=R=deg f
and kappa=kappa_infinity.

Each residue is quadratic algebraic over Q.
The discriminant of its equation is

\[
 -4(r-1)(s+1)(r+s+1)<0,
\]

and its normalized rational trace is -(r-1)/(s+2).
Since kappa_infinity is the sum of the finite residues,
taking normalized traces in a common number field gives

\[
 R-1=\sum_i(r_i-1).
\]

If there are ell poles, R=sum r_i-ell+1.
But f is polynomial and R>=sum r_i, so ell<=1.
There is at least one pole, so ell=1 and f has no other zeros.
Thus p=kappa/(w-a) and f=C(w-a)^r.

Any pole of q elsewhere would force another zero of f through
the squared-denominator divisibility. At a its pole order cannot
exceed floor(r/2). This proves KMFI33.

This is a necessary normal form, not a classification of all
nonhomogeneous q or all their lower lambda coefficients.
It leaves a specific rational, nonhomogeneous search space.

## 5. Complete homogeneous normalization

For rho=kappa z/w, polynomiality of P requires w^2|lambda.
On a monomial w^r z^s, the left side of (1) is

\[
 \{r(r-1)+2\kappa(r-1)(s+1)
       +\kappa^2(s+1)(s+2)\}\,w^{r-2}z^s.
\]

The monomials do not mix.
For fixed kappa, at most one pair r,s can have zero coefficient.
Indeed the quadratic in kappa is irreducible over Q.
Two pairs sharing a root must have proportional quadratics.
Their ratios determine
t=(s+1)/r and u=(s+1)(s+2)/(r(r-1)), and then

\[
 r=(u+t)/(u-t^2),\qquad s+1=tr.
\]

The denominator is nonzero.
Thus the solution space is one dimensional when nonzero.
A nonzero scalar lambda and the integration constants of X
are removed by scaling and target translation for field purposes.
Take lambda=w^r z^s with zero constants.

Put t0=w^(r-1)z^(s+1) and

\[
 a=-\kappa-r/(s+1),\quad
 \beta=1/(s+1),\quad
 \chi={2\kappa(s+1)+r\over(s+1)(s+2)}.
\]

Direct integration gives

\[
 X=(a t_0,\ (a+\beta)t_0w,\ (a+\chi)t_0z).
 \tag{6}
\]

The quadratic residue relation (3) is exactly its integrability
condition. Both beta and chi are nonzero and distinct.

For Y=S(gamma,w,z), set xi=t0/Y1. Then

\[
 Y_2=wY_1(1+\beta\xi),\quad
 Y_3=zY_1(1+\chi\xi),\quad
 T={Y_2^{r-1}Y_3^{s+1}\over Y_1^{r+s+1}}
   =\xi(1+\beta\xi)^{r-1}(1+\chi\xi)^{s+1}.
 \tag{7}
\]

This gives the fixed polynomial f(xi) of degree N=r+s+1.

## 6. Exact independent stabilization and the full curve group

There is no hidden torus or constant extension in (7).
Ignoring rational xi factors, the exponents of (w,z) in
Y1,Y2,Y3 are the columns

\[
 (r-1,s+1),\quad(r,s+1),\quad(r-1,s+2).
\]

Their minors have gcd gcd(r-1,s+1,N)=1.
The primitive kernel vector is (-N,r-1,s+1).
A unimodular target monomial basis (T,U,V) therefore maps to
a source basis in w,z after xi is adjoined.
Consequently the exact field pair is

\[
 k(T,U,V)\subset k(\xi,U,V),\quad T=f(\xi),
\]

with U,V independent. This is the independent curve class, not
just a compatible curve slice.

The residue relation gives

\[
 f'=(1+\beta\xi)^{r-2}(1+\chi\xi)^s(1-a\xi)^2.
\]

The last critical point is distinct from the three zeros of f
and has nonzero branch value. Its inertia is one three-cycle.
The zero branch has cycles (r-1),(s+1),1; infinity is an N-cycle.
The finite index sum is (r-2)+s+2=N-1.

This group is primitive.
In a nontrivial block system a single three-cycle cannot permute
blocks of size greater than one, so its block action is trivial.
The block quotient would be cyclic, generated by the zero-branch
permutation. That permutation has a fixed sheet and therefore a
fixed block, which is incompatible with a nontrivial transitive
cyclic action. Hence no proper block system exists.

The normal closure of the three-cycle is transitive by primitivity.
The three-cycle support lemma gives A_N inside the group.
The zero-branch permutation is even exactly when N is odd.
Thus the full group is A_N for odd N and S_N for even N.
This proves the field-monodromy part of KMFI34, including N=3.

## 7. Centered affine-origin obstruction and controls

In (6), the first sweep equation at Y=0 gives gamma=-a t0.
The second then gives beta*w*t0=0.
Since r>=2,s>=0, either w=0 or t0=0 forces t0=0,
and then gamma=0.
Thus no preimage of the origin has gamma!=0.
The positive-weight Gao stage, which needs such a marked point,
cannot lift this centered homogeneous sweep.

The statement does not prohibit a different presentation of
the same A_N or S_N field; those actions have actual realizations
elsewhere in the theory. It also does not extend the origin
obstruction to arbitrary translated/deformed rational multipliers.

A genuinely nonhomogeneous rational control remains:
for any polynomial R(w), put kappa^2+kappa+1=0,
rho=(kappa*z+R(w))/w and lambda=w^2.
With zero-constant primitives I'=R, J'=wR, K'=R^2, take

\[
 X_1=-(\kappa+2)wz+wR+2\kappa I,\qquad
 X_2=\kappa^2w^2z+w^2R+2(\kappa-1)J,
\]
\[
 X_3=-wz^2+wRz-K.
\]

Differentiation gives the required polynomial nilpotent matrix
and det JS=gamma^2. But X(0,z)=0 for the entire line,
so S is nonfinite. This is a solution of the integrability equation,
not an affine Keller realization or an instance of the finite
inertia conclusion. The checker evaluates R=w^2+w.

Exact controls include the sharp signed delta-two surface,
polynomial degree balances, residue/trace identities, the complete
homogeneous kernels through selected degrees including eleven,
and a genuine nonhomogeneous rational-multiplier solution.
No formal verification or independent whole-proof review is claimed.

The [rational sweep classification](RATIONAL_SWEEPS.md), KMFI35--37,
now settles the entire z-dependent rational-multiplier slice with
lambda independent of z, including nonfinite sweeps and allowed
higher pole orders. Its nonhomogeneous groups are alternating;
the homogeneous cases have the groups above. An explicit
nonfinite rational slice also has complete polynomial affine lifts.
Positive z degree in lambda remains outside that classification.

[KMFI38--40](TAUTOLOGICAL_BOUNDARY.md) now show that every genuine
w-pole collapses a source line and that every z-dependent multiplier
in this rational chart has a birational ramification surface, even
with positive z degree in lambda. The full z-dependent rational class is now settled by
[KMFI41--43](RATIONAL_CHART_CLASSIFICATION.md), including every
positive z degree in lambda. The natural degree12
tautological Gao route is excluded without assuming finiteness;
an exceptional degree11 proposal must supply two additional ramified
normalization branches over its branch surface outside the source.

The latter classification proves finiteness on Y2!=0 and puts all
omitted affine normalization branches over a separate plane. Hence
the formerly required two-branch coincidence over the ramification
surface cannot occur in this z-dependent rational chart.
