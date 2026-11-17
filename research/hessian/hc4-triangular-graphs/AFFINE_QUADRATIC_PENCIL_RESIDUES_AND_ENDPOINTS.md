# Affine quadratic pencils: reciprocal Hessian sums and the two old endpoints

15 November 2026. Claims **HC4MYJF0/1/2**. The inner quadratic may
have an independent linear offset; no common center is assumed.
A rational pencil test excludes unit Hessian in an entire polynomial
source orbit, in all outer degrees. A cyclic-vector certificate makes
that test uniform on a nonempty open family. Every nilpotent affine
pencil reduces, for the universal inverse question, to the two
existing plane Keller endpoints. Their resolution and unrestricted
HC4/JC2 remain open.

## 1. The actual affine stationarity equations

Work over C in n variables. Let A=A^T be invertible, B=B^T, and put

    T=A^(-1)B,
    q_A(x)=x^T A x,
    g(x)=x^T B x+2 v^T A x+c,
    h(x)=q_A(x)+H(g(x)).                                 (1)

This permits every affine quadratic inner polynomial at the unique
center of q_A. Define rational functions of s:

    M(s)=I+sT,             P(s)=det M(s),
    x(s)=-s M(s)^(-1)v,
    R(s)=g(x(s)),
    F(s)=s-H'(R(s)).                                     (2)

Away from P=0, the critical points of h correspond bijectively to
zeros of F, through x(s). Indeed

    grad h=2A[(I+H'(g)T)x+H'(g)v].

At a critical point set s=H'(g(x)); the equation is
M(s)x=-s v. Conversely F(s)=0 makes x(s) critical.
Two parameter values cannot give the same critical point, since
s=H'(g(x)).

One useful identity, including its constants and sign, is

    R(s)=d/ds [s c-s^2 v^T A M(s)^(-1)v].                (3)

The bracket is q_A(x(s))+s g(x(s)); differentiation uses its
stationarity in x, so its derivative is g(x(s)).

At a regular-pencil critical parameter s, direct differentiation
and the rank-one determinant identity give

    det Hess h(x(s))=2^n det(A) P(s) F'(s).               (4)

For detail, write N=A+sB and z=grad g(x(s)). Differentiating the
stationarity equation gives x'=-N^(-1)z/2, hence
R'=-z^T N^(-1)z/2. The critical Hessian is
2N+H''(R) z z^T. Its determinant is
2^n det(N)(1-H''(R)R'), exactly (4).
This is a literal native critical Hessian identity.

## 2. JF0: a complete rational boundary test

Suppose F is not identically zero and both conditions hold:

1. For every root alpha of P, the equation
   (I+alpha T)x=-alpha v has no solution.
2. The rational function W=1/(P F) is regular at every root of P,
   and W(s)=O(s^(-2)) at infinity.

Then no polynomial source automorphism Phi makes h composed with
Phi have nonzero constant Hessian determinant.

Here is a proof by partial fractions, without importing a
multivariable residue theorem. Condition 1 excludes every critical
point at a singular pencil. By (2), all remaining critical points
are zeros of F. At such a zero, (4) says degeneracy is equivalent
to F'=0. Thus a hypothetical unit-Hessian source orbit would make
all these zeros simple.

It would also make their native Hessian determinants equal to
one common C!=0. If det D Phi=delta, the critical chain rule is
Hess(h composed with Phi)=D Phi^T Hess(h) D Phi; delta is one
constant, not a separate rescaling at different points.
This reproduces the retained critical-scheme argument of
[HG0](GLOBAL_VALUE_ORBITS.md#1-global-invariants-and-the-whole-value-profile-hc4myhg0).

By condition 2 all finite poles of W are therefore simple zeros
s_i of F outside P=0. Their partial-fraction coefficients are

    1/(P(s_i)F'(s_i))=2^n det(A)/C.                      (5)

There is at least one such pole: a nonzero rational function with
no finite pole would be a polynomial and could not decay as
s^(-2). Subtract from W the sum of (5)/(s-s_i). The difference
has no finite pole and tends to zero, so it is the zero polynomial.
Comparing its coefficient of 1/s at infinity makes the sum of
all coefficients (5) zero.
If there are N>=1 poles, the sum is N*2^n det(A)/C!=0.
This contradiction proves JF0. The argument uses neither Tang,
properness of a gradient map nor an unresolved Keller theorem.

Both boundary conditions are genuine. In particular a critical
point at P=0 may not be omitted from the count. The estimate at
infinity also cannot be discarded: a nilpotent centered example
below has W=1/s and actual unit Hessian.

A convenient sufficient finite-pole check is this. If R has a pole
of order r_alpha at a root alpha of P of multiplicity k_alpha,
and m=deg H'>=1, then F has pole order m r_alpha there. Hence
m r_alpha>=k_alpha makes W regular there. Nonsolvability in
condition 1 must still be checked; pole orders alone do not
silently supply it.

## 3. JF1: a cyclic vector excludes generic independent centers

Assume deg H>=2 and v is cyclic for T:

    det[v,Tv,...,T^(n-1)v]!=0.                            (6)

The following each imply that no source parameterization of h
has unit Hessian:

- T is invertible;
- deg H>=3;
- zero is an eigenvalue whose algebraic multiplicity is at least two.

This holds in every n>=1 and every admitted degree. In particular
invertible cyclic pencils cover all nonlinear H, including degree
two. In four variables the conditions det T!=0 and (6) define
a nonempty open set in the affine-quadratic coefficients. For
A=I, T=diag(1,2,3,4), v=(1,1,1,1), the determinant in (6) is 12.
The constant c and every coefficient of nonlinear H are unrestricted.

To prove the claim, use orthogonal generalized eigenspaces of the
A-self-adjoint T. Cyclicity means there is one Jordan block on
each primary space, and the projection v_lambda is cyclic there.
Put N_lambda=T-lambda I and let its primary dimension be l.
The Gram matrix on the chain v_lambda,N_lambda v_lambda,... is
anti-triangular with top coefficient

    a_lambda=A(v_lambda,N_lambda^(l-1)v_lambda)!=0.       (7)

Nondegeneracy of that primary space proves the inequality, as in
the chain argument of [JD0](NILPOTENT_QUADRATIC_PENCILS_AND_ENDPOINT_ORBITS.md).
This is finite-dimensional linear algebra, not an ambient
polynomial rectification assumption.

At alpha=-1/lambda for lambda!=0, the inverse expansion is

    M(s)^(-1)=sum_{j=0}^{l-1}
             (-s)^j N_lambda^j/(1+lambda s)^(j+1).

By (7), v^T A M^(-1)v has an exact pole of order l, and (3)
makes R have pole order l+1. Also v_lambda is not in the image
of N_lambda: its top iterate is nonzero. The singular-pencil
stationarity equation is therefore inconsistent. Since P has
order l there, W is regular for every m=deg H'>=1.

Only infinity remains. If T is invertible, x(s) has a finite
limit and R is bounded; F=s+O(1) and W=O(s^(-n-1)).
If the zero primary space has dimension l0>=1, the same chain
calculation gives a polynomial part of R of exact degree l0.
For m l0>=2, F has exact degree m l0 at infinity, while
deg P=n-l0. Thus

    W=O(s^(-n-(m-1)l0)).

The exponent is at least two in the latter two bullets above.
JF0 applies in every case.

### The exact remaining cyclic boundary in four variables

A cyclic unit-Hessian source orbit not excluded above must have
H quadratic and zero an eigenvalue of algebraic multiplicity one.
Write H'(r)=a r+b, a!=0. In this case

    R(s)=alpha s+beta+O(s^(-2)),     deg P=3.             (8)

There is no 1/s term. On each nonzero primary space the limit
of x(s) is a critical point of g; its first-order contribution
to R vanishes. Equivalently differentiate (3) at infinity.
On the one-dimensional zero primary space R is linear.

Unless both

    1-a alpha=0,              a beta+b=0,                (9)

hold, F has a nonzero linear or constant leading term, so W
still decays at least as s^(-3). JF0 excludes it. Consequently
(9) is a necessary double infinity resonance for any remaining
four-variable cyclic unit-Hessian orbit. It is not an existence
certificate or an inverse theorem for the resonant case.

### Actual graph images

For
Q(u,X)=q_A(u)+a0 X^2+b0 X+c0, a0!=0, and
X=P0(g(u)) with P0 nonconstant, set
H=a0 P0^2+b0 P0+c0. Every closed parameterization of that image
is its graph map composed with a polynomial source automorphism.
The preceding nonexistence conclusions therefore hold in every
source degree. They cover all nonconstant P0 for invertible
cyclic T, and every deg P0>=2 for cyclic T even when B is singular.
No common-centering certificate is needed.

Certified polynomial quadratic similitudes and the literal
Meng--Yang transfer of HP0/HM0 preserve the actual source function
up to a scalar and a constant, so the same conclusions apply to
their image saturations. No automatic symmetry rectification is
asserted.

## 4. JF2: every nilpotent affine pencil has only the two old endpoints

Now let n=4, T nilpotent, and allow every polynomial H.
Then P=1 and x(s), R(s), F(s) are polynomials, by the finite
geometric inverse of I+sT. The stationarity description is
global, including all source points.

- If deg F>=2, JF0 excludes unit Hessian in every source orbit,
  since W=1/F=O(s^(-2)).
- If F is identically zero, every x(s) is critical.
  Since H'(R(s))=s, x(s) is nonconstant. Its image has a
  positive-dimensional critical closure, again excluding
  unit Hessian.
- If F is a nonzero constant, there are no critical points.
- If F has degree one, there is exactly one critical point,
  reduced by (4).

Thus any hypothetical nonautomorphic unit-Hessian member in
this entire nilpotent affine-quadratic source family yields
respectively the existing unit-ideal or singleton-zero plane
Keller endpoint. This uses retained Tang/HG0 and
[HG2](GLOBAL_VALUE_ORBITS.md#3-the-two-known-endpoints-are-full-source-orbits-hc4myhg2):
the full critical scheme of a nonautomorphic cotangent form
is the reduced zero fibre of its plane Keller pair.
No simultaneous centering or unproved normal form for each
individual degree-one F potential is used.

Both endpoints are actually represented by fixed graph images
in this family, over all polynomial source parameterizations:

    Q=x^2+y^2+z^2+w^2-X^2,       X=x+1.

The potential is y^2+z^2+w^2-2x-1, a polynomial coordinate
with polynomial recovery
x=(y^2+z^2+w^2-h-1)/2. Its entire source orbit is the coordinate
orbit of HG2. Here B=0, v=(1/2,0,0,0), c=1, H=-t^2, and F=2.

For the other endpoint take

    Q=2xy+z^2+w^2+X^2,          X=x^2.

Here T is nonzero nilpotent, v=0, R=0, F=s.
The shear y->y+x^3/2 gives the nondegenerate quadratic source
orbit. Its native Hessian determinant is -16. In particular
W=1/s shows exactly why the infinity condition in JF0 is needed.

The two ambient forms are constant-linearly equivalent complex
nondegenerate quadratics, so these are genuine instances of the
same quadratic-ambient setup and transfer to Meng--Yang by HM0.
Their graph maps and every source automorphism are closed
embeddings. Therefore the universal gradient-inverse question
for all nilpotent affine-quadratic-composite graph images is
equivalent in both directions to resolving both existing plane
Keller endpoints. Neither endpoint is solved here.

For deg H>=3, any nonconstant R already makes deg F>=2.
A remaining nilpotent orbit in those outer degrees must have
R constant and hence F of degree one. For graph degree
deg P0>=2 this removes the empty-critical-point alternative;
it is not a search through progressively higher degrees.

## 5. Literature and remaining scope

The selected primary record
[Soprunov, arXiv:math/0511684v2](https://arxiv.org/abs/math/0511684v2)
provides context for global residues and interpolation.
No toric or multivariable Euler--Jacobi theorem is imported,
and no missing-zero-at-infinity hypothesis is assumed.
The required one-variable partial-fraction identity and all
pencil boundary conditions are proved in Sections 1--3.

The all-centered family of JE0/1 is retained. The results here
extend to independent linear offsets on the stated cyclic and
nilpotent families. Noncyclic affine pencils, the cyclic double
infinity resonance, arbitrary image geometry, full w-dependent
recovery, intrinsic ambient alignment, both actual endpoints
and unrestricted HC4/JC2 remain open.
