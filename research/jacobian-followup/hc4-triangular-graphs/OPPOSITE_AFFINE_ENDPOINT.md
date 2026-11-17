# All-degree opposite affine-p graphs: only the base-only endpoint remains

For the original potential and notation in [CRITICAL_COMPONENTS.md](CRITICAL_COMPONENTS.md),
let `L*N != 0`, work in characteristic zero, and put

\[
 f(s,p,q,r)=\Psi(a(s)p+b(s),s,p,q,r),\qquad a,b\in K[s].
\]

Allow every actual polynomial automorphism of `(p,q,r)` over `K[s]`, with
polynomial inverse. There is no bound on coefficient degrees or automorphism
word length. The Hessian in the conclusion is always the actual Hessian after
that change.

**HC4MYGS2. If `a` is nonzero, no member of this class has nonzero constant
Hessian determinant.** The earlier [HC4MYGR3](CRITICAL_COMPONENTS.md) handles
`a(0) != 0`. The proof below handles the entire remaining boundary, without
successively increasing a degree cutoff.

**HC4MYGS3. Allowing arbitrary `a,b`, existence of a unit-Hessian member is
exactly existence of a unit-ideal plane Keller pair.** All such members must
have `a=0` and `b(0) != 0`; [HC4MYGU1/2/3](RELATIVE_COORDINATE_ENDPOINT.md)
proves the equivalence for that base-only class, and the fixed graph `x=1`
already realizes the whole endpoint. Equivalently, nonexistence for the whole
opposite affine-p family is the assertion that every complex plane Keller map
is surjective. That assertion remains open here. This is distinct from the
singleton-zero properness endpoint of the nonlinear p chart.

## 1. Preserve the whole gradient on the critical target line (HC4MYGS0)

For any polynomial `f(s,z)` and a relative polynomial automorphism
`Phi(s,z)=(s,phi(s,z))`, set `g=f composed with Phi`. The chain rule gives

\[
 \nabla_z g=(D_z\phi)^t(\nabla_z f)\circ\Phi,\qquad
 g_s=f_s\circ\Phi+\sum_i(f_{z_i}\circ\Phi)\phi_{i,s}.
\]

Since `D_z phi` is invertible over the polynomial ring, Phi identifies the
entire relative critical schemes `C_g` and `C_f`, and on those schemes
`g_s=f_s composed with Phi`. In particular it preserves every fibre, including
empty fibres, of the gradient over the target line

\[
 \ell=\{(\tau,0,0,0):\tau\in\mathbb C\}.
\]

This adds the base-gradient function to [HC4MYGK0](RELATIVE_CRITICAL_LOCUS.md).
The existing Hessian congruence also holds on `C_g`: the correction terms are
multiples of the relative partial derivatives, because the first coordinate
of Phi is exactly s. Thus if `delta=det D_z phi` is its nonzero constant
Jacobian, `det Hess(g)=delta^2 (det Hess(f)) composed with Phi` on `C_g`.
There is no assertion of nonlinear Hessian invariance away from this scheme.

It suffices to work over C. All uses of roots and poles below can be made
after passage to an algebraic closure; characteristic-zero finitely generated
coefficient fields embed in C, and nonexistence descends.

## 2. The entire critical curve in the remaining boundary

Assume henceforth `a != 0`, `a(0)=0`, `b(0) != 0`. At a zero of a, the dual
derivative rows of A,B are independent unless the old base is `(0,0)`:
[HC4MYGU1](RELATIVE_COORDINATE_ENDPOINT.md) computes their minors ideal
exactly `(x,y)`. At s=0 the old x is b(0), so no such point is critical.
At any other zero of a, s is nonzero. Also x=0 cannot be critical when a
is nonzero there: the q derivative is `N*s`. These arguments cover points
outside the rational chart; none are discarded by localization.

On `a(s)*x != 0`, use the relative computation coordinates `(x,T=A,q)`.
They have relative Jacobian `a(s)*x^3`. Put

\[
 p=(x-b)/a,\quad u=1+xs,\quad
 r=(u^3p+3xu^2q-T)/x^3.
\]

The potential is literally

\[
 F=LT^2+cT+e+dq,\quad
 c=M-N(2-3xs)/x^2,\quad
 e=N(1+xs)(2+xs)(x-b)/(ax^2),\quad
 d=2N(3+2xs)/x.
\]

The q equation forces `x_c=-3/(2s)` and hence s nonzero. The T equation and
then x equation give

\[
 T_c=\frac{13N/(2x_c^2)-M}{2L},\qquad
 p_c=\frac{x_c-b}{a},\qquad
 q_c=\frac{17T_c}{12x_c}+\frac{x_c-2b}{24ax_c},
\]

with r recovered above. The x equation solves uniquely because
`d_x=-6N/x^2` is nonzero. Consequently the **whole** relative critical scheme
is the reduced curve

\[
 C_f\simeq\operatorname{Spec}\mathbb C[s,1/(s a(s))].                 \tag{1}
\]

Projection to s is an isomorphism. In particular there are no extra vertical
components at roots of a, or points over s=0. These exclusions were proved in
the original polynomial coordinates, before using the rational chart.

Its critical value and base gradient are

\[
 D(s)=\frac{Ns(3+2sb(s))}{18a(s)}
       -\frac{(M-26Ns^2/9)^2}{4L},\qquad f_s|_{C_f}=D'(s).          \tag{2}
\]

The relative Hessian of F has determinant `-72*L*N^2/x^4` on its critical
curve. Schur complementation along that curve gives its full determinant as
this expression times D''. The native relative-coordinate Jacobian squared
therefore gives

\[
 (\det\operatorname{Hess} f)|_{C_f}
       =-162LN^2\frac{a(s)^2}{s^2}D''(s).                         \tag{3}
\]

All nonlinear coordinate corrections vanish here. If any permitted relative
automorphism produced a nonzero constant full Hessian, (3) and Section 1
would force, for a nonzero constant K,

\[
                     D''=K s^2/a^2.                              \tag{4}
\]

## 3. An all-degree pole classification (HC4MYGS1)

Write `k=ord_0(a)`, `d=deg(a)`. At a nonzero root rho of a, of multiplicity
m, the right side of (4) has pole order 2m. D has pole order at most m. A
pole of order n of a rational function has second derivative of order n+2
in characteristic zero; a regular function contributes no pole. Thus
`n=2m-2`, `1 <= n <= m`, forcing **m=2**.

At zero, (2) has leading term `N*s^(1-k)/(6*a_k)`, since b is regular and
`3+2sb` has nonzero constant term. If k>=2, (4) forces `k+1=2k-2`, hence
**k=3**. Otherwise k=1. In particular d is odd.

If d>=3, the right side of (4) decays at infinity with order `2d-2`.
The polynomial part of D must then have degree at most one. Its nonzero
proper part starts with `s^(-n)`, where `n=2d-4`. A rational function with
denominator dividing a has `n<=d`. Hence d<=4, and oddness forces **d=3**.
If d=1 it remains available separately. There are only three possibilities:

\[
 a=cs,\qquad a=cs^3,\qquad a=cs(s-\rho)^2\quad(c\rho\ne0).         \tag{5}
\]

For `a=cs^3`, the principal part of D at zero is
`N*s^-2/(6c)+N*b(0)*s^-1/(9c)`. Equation (4) has only an order-four pole;
the nonzero order-three term in D'' contradicts it. This case is impossible.

The remaining necessary packets are explicit; neither is asserted to
satisfy the full Hessian identity away from the critical curve.

**Packet I.** For `a=cs`, equation (4) is equivalent to

\[
 b=b_0+b_1s+\frac{169cN}{9L}s^3,\qquad
 b_0\ne0,\quad b_1+13cM/L\ne0,
\]

and

\[
 D=q_2s^2+q_1s+q_0,\quad
 q_2=\frac{N}{9c}(b_1+13cM/L),\quad q_1=\frac{Nb_0}{9c},\quad
 q_0=\frac{N}{6c}-\frac{M^2}{4L}.
\]

Here `C_f=C*`. On the target line, the gradient has exactly one preimage
for every `tau != q1`, and none for `tau=q1`.

**Packet II.** For `a=cs(s-rho)^2`, (4) gives

\[
 D=\frac{e_0}{(s-\rho)^2}+A_0s+B_0,\qquad e_0\ne0,\qquad
 B_0=-\frac{M^2}{4L}+\frac{N/(6c)-e_0}{\rho^2}.
\]

Conversely this is a polynomial b packet, since

\[
 b=\frac{(18c/N)\{e_0+(A_0s+B_0)(s-\rho)^2-E(s)(s-\rho)^2\}-3}{2s},
 \quad E(s)=-\frac{(M-26Ns^2/9)^2}{4L}
\]

has numerator zero at s=0. Retain the original condition
`b(0)=9c*rho^2*(A0+2e0/rho^3)/N-3/rho != 0`.
Here `C_f=C minus {0,rho}` and

\[
                         D'=A_0-2e_0/(s-\rho)^3.
\]

The gradient on the target line has three preimages for every
`tau` outside `{A0,A0+2e0/rho^3}`, two at the second listed value, and none
at A0. This is an exact fibre count, not a formal local degree or a claim
of properness. The full native determinant along C_f is the constant
`-972*L*N^2*c^2*e0`.

## 4. A general critical target-line degree obstruction (HC4MYGS0)

The argument in this section applies to any four-variable unit-Hessian
potential: if a gradient target line has both attained and omitted points,
its generic fibre cardinality equals the geometric degree of a plane Keller
map. In particular it cannot be one or three. Section 1 makes this criterion
usable after every retained-dual automorphism without computing its Hessian.

We use two classical inputs, specified and byte-pinned in the
[import manifest](sources/opposite-affine-import.json):

1. Nguyen Van Chau, [arXiv:0711.3894v1](https://arxiv.org/html/0711.3894v1),
   **Theorem 4**: a nonsingular polynomial plane map cannot have an
   irreducible nonproper-value component isomorphic to A1. In particular a
   plane Keller map's nonproper locus contains no affine line.
2. S. Yu. Orevkov, [On three-sheeted polynomial mappings of C2](https://www.math.univ-toulouse.fr/~orevkov/jc86.pdf),
   **Theorem 1.1**: no complex plane Keller map has geometric degree three.
   The first paragraph of Section 1 also records the classical birational
   Keller criterion: geometric degree one forces a polynomial inverse.
   The degree-three theorem is already registered as [PCC3](../../archive/non-elliptic/plane-jc/JC2_FINITE_NORMALIZATION_FRONTIER.md).

Suppose a permitted pullback g had unit Hessian. Both packets give an
omitted point of its gradient target, invariant by Section 1. Thus its
gradient is not an automorphism. Apply the retained Tang alternative
[HCTANG2](../tang-hc4-audit/STRUCTURAL_REVIEW.md). A constant affine change
of source coordinates puts the potential in cotangent form

\[
               h(U,V,S,T)=D_0(U,V)+S F(U,V)+T G(U,V),
\]

with `(F,G)` a plane Keller pair. Indeed its Hessian determinant is
`J(F,G)^2`, so the residual Jacobian is a nonzero constant. For each output,
the last two gradient coordinates specify `(F,G)`, and the first two solve
uniquely for `(S,T)` for every preimage `(U,V)`. Therefore its gradient
fibre over any target point has **exactly the same cardinality** as the
fibre of `(F,G)` over the projected target point.

The image of the original target line under the constant gradient-coordinate
change projects affinely to the `(F,G)` target plane. This projection cannot
be a point: a fixed residual pair value would make every point of the line
attained, or every point omitted. Each packet has both kinds. It is therefore
an affine line lambda, with the original parameter mapping bijectively to it.

By Chau's theorem, lambda is not contained in the plane nonproper locus.
It meets that algebraic curve in finitely many points (or the locus is
empty). Off those points the plane map is finite etale and its fibres have
its geometric degree. The packet's generic line fibre count is consequently
the actual residual plane geometric degree. Packet I forces degree one,
so the residual pair and full gradient have polynomial inverses, contradicting
the omitted point. Packet II forces degree three, contradicting Orevkov.
This excludes both packets and proves the remaining boundary case.

For completeness: `a(0) != 0` was excluded by HC4MYGR3. If `a(0)=b(0)=0`,
the entire s=0 trace is `H(p)` and [HC4MYGU0](RELATIVE_COORDINATE_ENDPOINT.md)
excludes unit Hessian, including after the allowed dual substitutions.
Together these prove HC4MYGS2 for every nonzero a, and the introductory base-only
endpoint proves HC4MYGS3 for the entire opposite affine-p family.

## 5. Verification and exact boundary

[verify_opposite_affine.py](verify_opposite_affine.py) checks symbolic
critical equations, coordinate and Hessian factors, both rational packets,
their fibre polynomials, polynomial divisibility and source hashes.
[The receipt](receipt-opposite-affine.json) reports its actual controls.
The pole-order classification and geometric exclusion are written proofs;
finite examples do not certify the all-degree quantifiers or the imports.
Tang, Chau and Orevkov are explicit mathematical inputs. No formal
verification, independent external review or novelty priority is claimed.

Still open are the unit-ideal plane-surjectivity endpoint, the separate
singleton-zero properness endpoint, nonlinear y=R(x,p) with a nonzero
no-simple-root boundary avoiding the critical p value, q/r-dependent omitted
base graphs, and genuine mixing of the free base with the duals. This theorem
does not close unrestricted HC4 or JC2. Earlier receipts retain their
then-current boundary statements and are superseded here only in scope.
