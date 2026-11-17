# A genus-one cone trace closes every cyclic affine quadratic pencil

15 November 2026. Claims **HC4MYJG0/1**. The cyclic quadratic
double-infinity resonance left by JF1 is now inverse-closed.
Together with JF1, every four-variable cyclic affine-quadratic
composite graph with nonconstant outer graph polynomial is
inverse-closed in every admitted graph and source degree.
This does not classify noncyclic affine pencils or solve either
old plane Keller endpoint or unrestricted HC4/JC2.

## 1. JG0: a whole-source-orbit cone-trace criterion

Let b(z) be a nondegenerate ternary homogeneous quadratic over C,
a(z) any homogeneous quadratic, and ell(z) a linear form. Put

    f(z,t)=a(z)+ell(z)+t b(z).                            (1)

Let D={b=0} in P2, a smooth conic. Assume:

- The restriction a_D is not a scalar multiple of ell_D^2.
  This includes a_D!=0 when ell=0.
- At every multiple zero of a_D on D, ell is nonzero.

Then for every polynomial source automorphism Phi of A4,
unit Hessian of f composed with Phi forces its polynomial gradient
inverse. No bound on the source parameterization degree is used.

Here is a direct proof of the curve input, including the boundary.
The conic D is P1: after a constant normalization of b its
parameterization is [p:q] -> [p^2:pq:q^2].
The restrictions of a and ell are a binary quartic A and a binary
quadratic L. The hypotheses say A is not in C L^2 and any common
zero of A and L is a simple zero of A.

For generic kappa!=0 the binary quartic

    Delta_kappa=L^2+4kappa A                             (2)

has four distinct zeros on P1. If L=0, the second hypothesis
already makes A squarefree. Otherwise the fixed common-zero
divisor of A and L^2 is reduced by the second hypothesis.
Away from it, the nonconstant rational function A/L^2 has only
finitely many critical values in characteristic zero: its
derivative is a nonzero rational function. Excluding those values
and the finitely many values which make a moving zero meet a
fixed zero makes every zero of (2) simple. There are four in all,
including zeros at infinity of either chosen P1 chart.

Consider the actual projective curve

    Cbar_kappa={b(z)=0,
               a(z)+ell(z)Z0-kappa Z0^2=0} in P3.        (3)

It misses the cone vertex [0:0:0:1], because kappa!=0.
Projection to D is a degree-two finite cover. Locally choose
a representative v of a conic direction; its equation is the
monic-up-to-unit quadratic
a(v)+ell(v)tau-kappa tau^2=0. Thus its branch divisor is (2).
Simple discriminant zeros make the actual projective curve smooth,
not merely its normalization. The cover is connected, since four
simple branch zeros prevent its discriminant from being a square.
The retained characteristic-zero
[Riemann-Hurwitz formula](https://stacks.math.columbia.edu/tag/0C1B)
gives

    2G-2=2*(-2)+4=0,         G=1.                        (4)

The affine curve C_kappa={b=0,a+ell=kappa} is exactly (3) with
Z0=0 removed. Write r for the number of distinct zeros of a_D.
It removes exactly r points, 1<=r<=4, not their intersection
multiplicities. In particular a tangential infinity intersection
does not by itself make the curve singular: ell at such a point
is nonzero. The retained affine-curve Hodge formula gives, with
L0=uv and U=u+v,

    E(C_kappa)=L0-U+1-r.                                (5)

On b!=0, (1)=kappa solves for t and gives A3 minus the cone.
The cone has E=L0^2: in the normalization b=XZ-Y^2, the X!=0
part is Gm times A1 and X=0 is one affine line.
On b=0 the whole value fibre is C_kappa times A1_t.
Therefore the entire generic value fibre has

    E(f=kappa)=L0^3-L0 U+(1-r)L0.                        (6)

The nonzero genus term contradicts the
[HG0 whole-value profile](GLOBAL_VALUE_ORBITS.md#1-global-invariants-and-the-whole-value-profile-hc4myhg0)
L0^3-N L0 for a nonautomorphic unit-Hessian source orbit,
for every critical count N. This proves JG0.

All product factors above are actual solved strata. The degree-two
cover is handled by (4), not by multiplying the E-polynomial of
its base by two. The Hodge and Riemann-Hurwitz foundations are the
explicit retained inputs of HC0 and IY0. No new generic-smoothness,
Bertini or complete-intersection genus theorem is imported.

## 2. A cyclic sufficient certificate for the cone criterion

Suppose a=q_A=z^T A z and b=q_B=z^T B z are both nondegenerate
in three variables. Set S=A^(-1)B and let u be cyclic for S:

    det[u,Su,S^2u]!=0.
    ell(z)=-2 u^T A z.                                  (7)

Then both hypotheses of JG0 hold.

First suppose q_A on D were gamma ell^2. The quadratic
q_A-gamma ell^2 vanishes on the irreducible conic q_B=0, so
it equals mu q_B. If mu=0, nondegeneracy of A would fail.
Otherwise A-mu B has rank at most one, giving an eigenspace
of S of dimension at least two. This contradicts cyclicity
in dimension three. The argument includes gamma=0.

Next a multiple zero p of q_A on D satisfies
A p=tau B p, with tau!=0, hence p is an isotropic eigenvector
of S. Distinct generalized eigenspaces are A-orthogonal.
On its primary space, cyclicity gives one Jordan chain with
cyclic vector u_lambda. The Gram matrix on
u_lambda,(S-lambda)u_lambda,... is anti-triangular; its top
pairing is nonzero by nondegeneracy. Consequently

    u^T A p !=0,

since p spans the last eigenline of that chain. A block of
length one has no isotropic eigenline at all. Thus ell(p)!=0
at every multiple zero, proving the second hypothesis.
This is the same elementary self-adjoint cyclic pairing used
in JF1, with no polynomial rectification premise.

The geometric criterion can also apply without cyclicity.
No converse to (7) is asserted, and failure of either geometric
hypothesis is not declared an HC4 counterexample.

## 3. JG1: normalize the entire double-infinity resonance at once

Use the affine-quadratic setup of JF:

    h(x)=x^T A x+H(g(x)),
    g=x^T B x+2v^T A x+c,
    T=A^(-1)B,        det[v,Tv,T^2v,T^3v]!=0.             (8)

Let deg H>=2. We claim every unit-Hessian polynomial source
reparameterization of h has a polynomial gradient inverse.

JF1 already excludes unit Hessian unless H is quadratic,
zero is a simple eigenvalue of T, and both infinity resonances
hold. Split the A-orthogonal primary spaces and normalize the
one-dimensional zero block. In constant coordinates write

    q_A=q_{A3}(z)+w^2,       T=S direct-sum 0,
    v=(v3,d),               d!=0,
    S invertible,           v3 cyclic.                  (9)

Here A3 and B3=A3 S are nondegenerate. Put

    u=S^(-1)v3,       beta=c-v3^T A3 u.

The unique zero primary contributes -2d^2 s to the stationary
inner value R(s). Its constant term at infinity is beta.
Thus the double resonance, with H'(r)=a r+b, is exactly

    a=-1/(2d^2),      b=beta/(2d^2),
    H(r)=H(beta)-(r-beta)^2/(4d^2).                     (10)

There is an explicit polynomial source automorphism from new
coordinates (z,t) to the old coordinates:

    z_old=z-u,
    w_old=-d t-q_{B3}(z)/(4d).                          (11)

Its inverse is
z=z_old+u,
t=-w_old/d-q_{B3}(z_old+u)/(4d^2).
Direct substitution, keeping every constant, gives

    h composed with (11)
       =q_{A3}(z)-2u^T A3 z+t q_{B3}(z)+C,
    C=q_{A3}(u)+H(beta).                                (12)

Since S is invertible, u=S^(-1)v3 is cyclic as well.
Section 2 verifies the cone hypotheses for the potential
on the right of (12), so JG0 applies.

Equation (12) is an equality of full polynomial source orbits.
It is not Hessian invariance under the shear. Composing a
given unit-Hessian reparameterization with the inverse of
(11) merely represents its very same actual source function
as a source reparameterization of the cone-trace potential.
The harmless constant C shifts generic values. HG0 then
excludes nonautomorphy of that actual gradient.

This proves JG1 for every cyclic four-variable affine quadratic
pencil and every nonlinear H. The former small stationary-
numerator/congruence proposal is not needed: the genus-one
argument handles all resonance coefficients simultaneously.

## 4. Actual graph images and current limits

For the nondegenerate ambient quadratic

    Q(u,X)=q_A(u)+a0 X^2+b0 X+c0,  a0!=0,

and actual graph X=P0(g(u)) with P0 nonconstant, the native
potential is q_A+H(g), where
H=a0 P0^2+b0 P0+c0 has degree 2 deg P0>=2.
Every closed parameterization of that same graph is its graph
map composed with a polynomial source automorphism.
If the affine offset v is cyclic for A^(-1)B, JG1 therefore
forces gradient invertibility under actual unit Hessian in
every graph and source degree. Most cases retain the stronger
JF1 nonexistence statement; the former quadratic resonance is
inverse-closed here, not asserted nonexistent.

Certified polynomial quadratic-similitude image saturation
and the literal HM0 transfer to original Meng--Yang embeddings
preserve the result because the actual source function is
only rescaled and translated. No ambient symmetry rectification
or transfer of normality between different images is assumed.

The whole common-center family remains classified by JE.
All nilpotent affine offsets remain reduced by JF2 to the
two old endpoint alternatives for universal inversion.
A new affine-pencil image mechanism beyond those results must
therefore be nonnilpotent and noncyclic. This restriction is
within the admitted quadratic-composite image family, not the
entire unrestricted HC4 frontier.

Both actual plane endpoints, noncyclic nonnilpotent affine
pencils, different image geometry, full w-dependent recovery,
higher helper traces, intrinsic ambient alignment and
unrestricted HC4/JC2 remain open. No endpoint existence or
new nonautomorphic unit-Hessian example is constructed.

For related classical context, the inspected
[UT Austin lecture notes](https://web.ma.utexas.edu/users/ikmartin/notes/algebraic-curves-08/)
give the genus-one smooth two-quadric example. We do not import
their adjunction argument or any converse embedding assertion.
Here the cone, tangencies and missing points are accounted for
by the explicit double cover and retained Riemann-Hurwitz.
The checker supplies finite exact controls of this written
geometric proof, rather than a proof by testing Jordan types
or bounded source degrees.
