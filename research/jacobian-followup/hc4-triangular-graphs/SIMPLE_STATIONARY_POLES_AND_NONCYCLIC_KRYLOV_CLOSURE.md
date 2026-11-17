# Stationary pole groups and noncyclic orthogonal Krylov closure

15 November 2026. Claims **HC4MYJH0/1**. Compatible simple
determinant roots in an affine quadratic pencil can be counted
exactly: one pole represents a whole group of finite critical points.
This extends the JF0 reciprocal-Hessian obstruction beyond cyclic
offsets. In particular every invertible simple-spectrum pencil,
with every affine offset and every nonlinear outer polynomial,
has no unit-Hessian polynomial source orbit.

The original HC4 problem, both plane Keller endpoints and the
remaining noncyclic families are not resolved.

## 1. Setup and the new boundary-complete residue test (JH0)

Over C in n>=1 variables let A=A^T be invertible, B=B^T, and use
the actual source potential and rational stationarity data of JF:

    h=q_A+H(g),     g=q_B+2v^T A x+c,     T=A^(-1)B,
    p=H',          m=deg p>=1,
    P(s)=det(I+sT), x(s)=-s(I+sT)^(-1)v,
    R(s)=g(x(s)),  F(s)=s-p(R(s)),  W(s)=1/(P(s)F(s)).     (1)

The exact regular-parameter critical determinant is

    J(x(s)):=det Hess h(x(s))
             =2^n det(A)P(s)F'(s).                       (2)

This is the retained, explicitly derived [JF0 identity](AFFINE_QUADRATIC_PENCIL_RESIDUES_AND_ENDPOINTS.md).
No source inverse or common center is assumed.

Suppose F is not identically zero, W=O(s^(-2)) at infinity,
and each root alpha of P is of one of these two types:

- Blocked: (I+alpha T)x=-alpha v is inconsistent and W is
  regular at alpha.
- Simple stationary: P'(alpha)!=0 and that affine equation
  is consistent.

Then, if every native critical point is Morse, the complete
critical set is finite and nonempty and

    sum_{x in Crit(h)} 1/J(x)=0.                         (3)

Consequently no polynomial source automorphism Phi makes
h composed with Phi have nonzero constant Hessian determinant.

The second type is the new one. The proof below includes every
critical point above it; it does not discard singular pencils
or identify arbitrary determinant residues with finite masses.

## 2. One compatible simple pole is a whole critical group

Let alpha be of the second type, and put lambda=-1/alpha.
A simple root of P has a one-dimensional nondegenerate
A-eigenspace for lambda, orthogonal to the other primary spaces.
In constant coordinates normalize this direction to y with
q_A=y^2+q_A'(z). Consistency forces its component of v to vanish.
Thus the regular stationary solution has y=0, R extends regularly
to alpha, and the entire affine stationary set at alpha is

    z=z_alpha,     y free,
    g=R0+lambda y^2,     R0=R(alpha).                    (4)

Every root r of p(r)=alpha gives y^2=(r-R0)/lambda.
If F(alpha)=0, the root r=R0 gives a critical point with
a zero Hessian row in the y direction. If p'(r)=0 at a root,
the critical determinant also vanishes. Hence the Morse assumption
makes F(alpha)!=0 and all m roots simple and different from R0.
There are exactly 2m critical points in this group.

The singular diagonal entry of 2(A+alpha B) is zero.
Expanding its rank-one update by p'(r) grad(g) grad(g)^T gives

    J(x_{r,+})=J(x_{r,-})
       =2^n det(A) * 2(r-R0)p'(r)P'(alpha).              (5)

The formula is in the original source coordinates. Under the
constant normalization both J and det(A) acquire the same
squared determinant factor, which cancels in (5).

The proper one-variable partial-fraction identity

    1/(p(z)-alpha)
       =sum_{p(r)=alpha} 1/[p'(r)(z-r)]

at z=R0 gives

    sum_{p(r)=alpha} 1/[(r-R0)p'(r)]=1/F(alpha).

Keeping both signs of y in (5), the complete group therefore has

    sum_group 1/J
       =1/[2^n det(A)P'(alpha)F(alpha)]
       =Res_alpha(W ds)/[2^n det(A)].                   (6)

In particular one determinant pole is not one critical point.

### Completing the global sum

Away from P=0, zeros of F are in bijection with critical points.
Under the Morse assumption they are simple and their W coefficients
are 2^n det(A)/J by (2). Blocked roots have neither source critical
points nor W poles. At each compatible simple root, (6) accounts
for its whole group and its simple W pole. These are all the
finite poles and all the finite source critical points.

Partial fractions and W=O(s^(-2)) make the sum of the coefficients
of all finite simple poles zero. There is at least one pole:
a nonzero rational function with no finite pole would be a
polynomial and could not have that decay. Each pole corresponds
to a nonempty group or a regular critical point. This proves (3).

If h composed with Phi had unit Hessian, its critical points
would all be Morse and their native J values would equal one
common C!=0. The retained critical chain rule gives one constant
factor (det D Phi)^2 over the whole critical set. For N>=1 points,
(3) would say N/C=0, impossible in characteristic zero.
This proves JH0 without Tang, properness or an imported
multivariable Euler--Jacobi theorem.

A compatible multiple determinant root is not covered by (6).
It may carry an escaping contribution with no finite critical
point. Section 5 gives an explicit control of that boundary.

## 3. An orthogonal Krylov certificate (JH1)

Suppose T is invertible and let

    K=span(v,Tv,...,T^(n-1)v).

Assume the A-pairing restricted to K is nondegenerate, including
the zero-dimensional case. This does not require A(v,v)!=0:
the vector v itself may be isotropic across distinct eigenspaces.
Cayley-Hamilton makes K T-stable;
self-adjointness makes K-perp T-stable. Thus

    C^n=K direct-sum K-perp

is an actual constant orthogonal decomposition, with v in K.
This is linear algebra, not a polynomial coordinate or action
rectification claim.

Call a nonzero eigenvalue active if it occurs on K. Write l_lambda
for its algebraic multiplicity on K and k_lambda for its
multiplicity on the whole source. Suppose

    m(l_lambda+1)>=k_lambda           for every active lambda, (7)

and every inactive primary space is semisimple.
Then no source reparameterization of h has unit Hessian.

On K the vector v is cyclic by definition. Its nondegenerate
self-adjoint primary chains have nonzero top pairings, as
proved in JF1. Thus at an active root alpha=-1/lambda,
R has pole order l_lambda+1, its singular stationary equation
is inconsistent already after projection to K, and (7) makes
W regular at alpha. These roots are blocked.

At an inactive eigenvalue, v has no component in its primary
space. If that semisimple primary has dimension at least two,
the entire stationary set at alpha contains the quadric

    q_A(y)=(r-R0)/lambda,        p(r)=alpha.

It is nonempty and positive-dimensional for every such r over C.
This alone forbids unit Hessian in every source orbit.
This pre-exclusion does not require the other pole bounds.

Otherwise all inactive primary spaces are simple and their
determinant roots are simple stationary poles. JH0 applies:
T invertible makes R bounded at infinity and
F=s+O(1), W=O(s^(-n-1)). F is nonzero; either an active pole
makes that immediate or v=0 gives F=s-p(c).
This proves JH1.

### Complete simple-spectrum and semisimple consequences

If T is invertible with simple spectrum, the hypotheses hold
for every vector v and every nonlinear H. Each eigenspace is
a nondegenerate A-line. The Krylov space is the orthogonal sum
of exactly the lines where the projection of v is nonzero.
Its pairing is automatically nondegenerate; active l=k=1
satisfies (7) for all m>=1; inactive roots are simple.

Therefore:

    For every invertible simple-spectrum T, every affine offset v,c
    and every deg H>=2, no polynomial source orbit has unit Hessian. (8)

No generic-offset qualification is needed. If all native critical
points are Morse in this case, their number is exactly

    N=1+2mn.                                            (9)

Indeed each active root gives R a double pole and F a pole of
order 2m. F has one pole at infinity, so there are 1+2m*a
regular zeros, where a counts the active eigenvalues.
Every inactive eigenvalue supplies another 2m critical points.
No regular zero coincides with an inactive determinant root
under the Morse assumption.

More generally if T is semisimple, nondegeneracy of K says exactly
that every nonzero projected v_lambda has
A(v_lambda,v_lambda)!=0. Then l_lambda=1 and (7) reads

    2(deg H-1)>=k_lambda          for active eigenvalues. (10)

In four variables every deg H>=3 meets (10), because
k_lambda<=4. Thus every invertible semisimple pencil with
nondegenerate offset Krylov space is excluded in those outer
degrees. Every nonlinear H is excluded when active
multiplicities are at most two. Inactive repeated eigenvalues
are already excluded by their critical quadrics.

For actual graphs X=P0(g) in
Q=q_A+a0 X^2+b0 X+c0, a0!=0, P0 nonconstant,
H=a0 P0^2+b0 P0+c0 has degree 2 deg P0>=2.
Every closed parameterization factors by a source automorphism.
Thus (8) covers every such graph/source degree and every affine
offset for an invertible simple-spectrum pencil. The other
JH1 cases keep the exact bounds (7)/(10). Certified quadratic
similitudes and the literal HM0 transfer preserve these
nonexistence statements for genuine Meng--Yang images.

## 4. Why this is an all-degree extension

The critical-group formula (6) handles all m roots at once by
one partial fraction, rather than solving them or increasing
the graph degree. The Krylov certificate uses a constant
subspace, its actual pairing and eigenvalue multiplicities.
No source automorphism degree is bounded.

The prior JF0 theorem remains its all-blocked special case.
JG complete cyclic inverse closure and JE common-center
classification remain unchanged. The new result fills a
different gap: noncyclic offsets with actual finite critical
groups at simple determinant roots.

## 5. Exact limits and a genuine boundary residue

Take in two variables

    q_A=2xy,       T=[[1,0],[1,1]],       v=0, H(r)=r^2.

Then q_B=x^2+2xy, P=(1+s)^2, R=0, F=s, and

    W=1/[s(1+s)^2]
     =1/s-1/(1+s)-1/(1+s)^2.

At the determinant root s=-1, the stationary equation forces
x=0, so q_B=0 and H'(q_B)=0!= -1. There is no finite critical
point there, although W has a nonzero residue -1.
Thus it would be false to extend (6) automatically to a
compatible multiple root. It is an escaping boundary contribution,
not an omitted finite critical group.

Degenerate Krylov spaces are also outside JH1. For example
an isotropic nonzero offset in a repeated semisimple eigenspace
can make the inner R pole disappear while stationarity is
inconsistent; the determinant residue is then not covered by
the active regularity estimate.

The hypotheses are sufficient certificates, not a classification
of every failed certificate. Compatible multiple stationary
poles, degenerate Krylov spaces, uncovered high-multiplicity
active poles, singular-pencil infinity resonances, and the
two actual endpoints remain separate. The general JF0/JH0
rational tests and JG0 geometric criterion may still close
individual cases outside the stated corollaries.

Within HC4 the remaining affine-image candidates must escape
these results as well as the previous common-center, cyclic
and nilpotent classifications. Arbitrary image geometry,
full w-dependent recovery, higher helper traces, intrinsic
ambient alignment and unrestricted HC4/JC2 remain open.

For classical context see the already retained primary
[Soprunov record](https://arxiv.org/abs/math/0511684v2).
No toric or multivariable residue theorem is imported; the
complete finite group identity and all required pole accounting
are proved here. The checker supplies exact finite controls,
not a proof by testing outer degrees or source maps.
