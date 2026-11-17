# Affine polar curves close every invertible semisimple pencil

15 November 2026. Claims **HC4MYJI0/1**. The pencil is invertible
and semisimple; the affine offset, nonlinear outer polynomial and
polynomial source reparameterization have unrestricted degree.
No nondegeneracy assumption on the offset Krylov space remains.

The conclusion is gradient invertibility whenever a source
reparameterization has unit Hessian. It is not a claim that all these
source orbits have no unit-Hessian member. Singular pencils,
nonsemisimple affine pencils outside earlier certificates, both
actual plane Keller endpoints and unrestricted HC4/JC2 remain open.

## 1. Statements and retained foundations

Over C let A=A^T be invertible in four variables, B=B^T, and

    h=q_A+H(g),     g=q_B+2v^T A x+c,     T=A^(-1)B,
    d=deg H>=2,    L=uv,                U=u+v.

**JI1.** If T is invertible and diagonalizable, every polynomial
source automorphism Phi for which h composed with Phi has nonzero
constant Hessian determinant gives a polynomial gradient inverse.
This holds for every v,c,H, including isotropic projected offsets
and repeated eigenvalues. The actual graph family

    Q(x,X)=q_A(x)+a0 X^2+b0 X+c0,   a0!=0,   X=P0(g),

with P0 nonconstant has H=a0 P0^2+b0 P0+c0 of degree at least two.
Every closed parameterization of this fixed graph factors by a
source automorphism. Thus JI1 covers all its admitted graph/source
degrees, certified quadratic-similitude images and literal original
Meng--Yang transfers.

**JI0** is the all-offset plane-curve calculation in Section 2,
including its precise ramification and puncture counts. We use its
total genus, not an unproved assertion of connectedness.

The retained inputs are:

- [HG0](GLOBAL_VALUE_ORBITS.md): a hypothetical nonautomorphic
  unit-Hessian reparameterization in four variables has whole generic
  value-fibre polynomial E=L^3-NL, where N is its full critical count.
  Critical schemes and full value fibres survive actual source
  automorphisms; native Hessians at critical points change by one
  constant congruence factor.
- [HC0](HODGE_TRACE_CLOSURE.md): additivity, literal product strata,
  and E(C)=c_C L-G_C U+c_C-r_C for a smooth affine curve, where
  c_C is its component count, G_C its total compact genus, and r_C
  its number of normalization points at infinity.
- The already retained characteristic-zero
  [Riemann--Hurwitz formula](https://stacks.math.columbia.edu/tag/0C1B),
  explicitly used in [IY0](CONTRACTIBLE_RADICAL_CURVES_AND_PENCIL_DEFECT.md).
- [JH0/1](SIMPLE_STATIONARY_POLES_AND_NONCYCLIC_KRYLOV_CLOSURE.md):
  the complete simple stationary groups and the nondegenerate
  Krylov certificate with active budget 2(d-1)>=k_lambda for
  semisimple active primaries.
- The exact source/ambient scalar transfers
  [HP0](QUADRATIC_ISOMETRY_SATURATION.md) and
  [HM0](AMBIENT_QUADRIC_GRAPHS.md).

There is no finite-cover E multiplicativity, arbitrary source Hessian
invariance, proper-gradient theorem, or new unproved external input.
All values called generic avoid a finite set of values for finitely
many fixed polynomials. That set may depend on the actual offsets.

## 2. Plane curves in every affine offset (JI0)

### Two distinct nonzero eigenvalues

After a constant change write the plane potential as

    f(z,w)=a z^2+b w^2+e z+j w+H(z^2+w^2),
    a*b*(a-b)!=0.

This includes every two-variable affine-quadratic composite whose
nondegenerate quadratic pencil has distinct nonzero eigenvalues:
translate the inner quadratic to its center, simultaneously
diagonalize its two forms, and normalize the inner form.
Let C_kappa be the smooth affine curve f=kappa for generic kappa.
Put epsilon=d mod 2. If C_kappa has c_C components and total
compact genus G_C, then

    deg(t:C_kappa -> A1)=4,   t=z^2+w^2,
    finite ramification=4d,
    infinity ramification=2 epsilon,
    r_C=4-2 epsilon,
    G_C=2d-4+epsilon+c_C >= 2d-3+epsilon > 0.       (1)

These statements do not require C_kappa to be connected.

For finiteness, combine the equations t=z^2+w^2 and f=kappa.
Since a!=b, they solve both z^2 and w^2 as affine-linear
expressions in z,w with coefficients in C[t]. The two monic
relations have relatively prime leading monomials z^2,w^2
in total (z,w) degree. Their standard basis is
1,z,w,zw, so the map is finite free of degree four.

Its affine polar equation is

    J(t,f)=4(b-a)zw+2jz-2ew=0.                     (2)

It is a translated hyperbola or the union of two distinct
transverse lines. Its two points at infinity are [1:0] and [0:1],
where t is nonzero as a homogeneous quadratic. On a hyperbola
parameter r in C*, t has a pole of order two at each of r=0,infinity.
Since d>=2, f has pole order 2d at each end: its degree-two part
and affine terms cannot cancel H(t). Thus its map to kappa has
degree 4d. In the two-line case each line restriction has degree
2d. A generic kappa avoids the vertex and the finite critical
values of these restrictions. Hence C_kappa meets the smooth polar
transversely in exactly 4d distinct points. On the smooth curve
this intersection multiplicity is exactly the zero order of dt.
All finite ramification is simple, with total 4d.

Here is the full infinity calculation. For d=2e0, set
t=tau^(-1), z=tau^(-e0)Z, w=tau^(-e0)W.
The two scaled equations at tau=0 are

    Z^2+W^2=0,       a Z^2+b W^2=-h_d,

where h_d is the leading coefficient of H. They have four
solutions with ZW!=0 and Jacobian 4(b-a)ZW!=0. The implicit
function theorem gives four distinct branches, each with t a
simple pole. Affine terms have positive tau order and change
neither this count nor smoothness of the normalized branches.

For d odd use t=tau^(-2), z=tau^(-d)Z,
w=tau^(-d)W. The same four limiting solutions are paired by
tau -> -tau, (Z,W) -> (-Z,-W). There are two actual infinity
points, each of t-ramification index two. The finite free degree
four shows these are all the points above infinity; there are no
additional branches. Applying Riemann--Hurwitz to every compact
component and summing gives

    2G_C-2c_C=-8+4d+2 epsilon.

This proves (1), with actual punctures rather than their
intersection multiplicities.

### Scalar nonzero pencil in two variables

Write f=q+H(mu q+2<eta,z>+c), where q is nondegenerate in
two variables and mu!=0. Its generic smooth curve has these
three exact forms:

1. eta=0: d disjoint nonzero levels of q, so E=d(L-1).
2. eta nonzero and isotropic: in Witt coordinates
   q=2zw and g=2mu zw+2z+c. For kappa!=H(c), put t=g.
   Then 2z=t+mu H(t)-mu kappa-c. This is a degree-d
   polynomial with d distinct roots for generic kappa.
   Necessarily z!=0, and w is uniquely recovered.
   The whole curve is A1 minus d points, with E=L-d.
3. eta nonisotropic: normalize q=z^2+w^2 and the affine term
   to 2a z, a!=0. The exact polynomial presentation is

       z=(t+mu H(t)-mu kappa-c)/(2a),
       w^2=kappa-H(t)-z(t)^2.                     (3)

   The polynomial on the right has degree 2d. It has simple
   roots for generic kappa. Indeed every critical point of f
   has w=0: the alternative 1+mu H'(g)=0 makes its z
   derivative 2a H'(g) nonzero. The remaining z derivative
   is a nonzero polynomial of degree 2d-1, so f has only
   finitely many critical points. Repeated roots in (3)
   would give a singular point of its actual fibre.
   Thus the compact double cover is connected of genus d-1
   with two punctures, and

       E=L-(d-1)U-1.                              (4)

Equations (1)--(4) hold for arbitrary affine offsets.

## 3. A nonzero isotropic projected offset

Diagonalizable self-adjoint T gives an orthogonal sum of
nondegenerate eigenspaces. Suppose the projection v_lambda
to one eigenvalue lambda is nonzero and isotropic. Since T
is invertible, lambda!=0. Choose a constant Witt pair in
this eigenspace so that, for two remaining variables z,

    q_A=2xy+q(z),   g=2lambda xy+2x+g0(z),         (5)

where q is nondegenerate and the remaining projection to
lambda is zero. Write f0(z)=q(z)+H(g0(z)).
With t=g and

    M(t,z)=t+lambda H(t)+lambda q(z)-g0(z)
                                      -lambda kappa,

the fibre h=kappa has the following disjoint actual strata:

    x!=0:   2x=M(t,z),   y=(t-2x-g0(z))/(2lambda x),
            isomorphic to A3 minus V(M);
    x=0:    f0(z)=kappa, y arbitrary.

These maps and their inverses are regular on the stated strata.
Consequently the complete fibre polynomial is

    E(h=kappa)=L^3-E(M=0)+L E(f0=kappa).           (6)

If lambda has multiplicity k, where k is 2,3 or 4, the
quadratic part of lambda q-g0 has rank r=4-k. Its linear
part on its kernel is zero. Complete the square on its
nonzero eigenspaces. The exact generic base profiles are

    r=2: E(M=0)=L^2+(d-1)L;
    r=1: E(M=0)=L E(D),  D: y^2=t+lambda H(t)+c1-lambda kappa;
    r=0: E(M=0)=d L^2.                            (7)

For r=2 use a constant hyperbolic frame: M=0 is uv=K(t)
with d distinct roots. Its u!=0 stratum is C* times A1;
at each of the d roots the u=0 stratum is an affine line.
For r=1, D is the connected smooth double cover with
genus a_D=(d+epsilon-2)/2 and 2-epsilon punctures.
For r=0 the d roots give d separate affine planes.

We now compare (6) with HG0 under hypothetical gradient
nonautomorphy.

- k=2. The remaining pencil is invertible on two variables.
  If its eigenvalues are distinct, (1) gives a strictly
  positive coefficient G_C of -U in its curve profile.
  If it is scalar with a nonisotropic remaining offset,
  (4) again gives positive genus. Neither can cancel against
  the first profile in (7).
  If the remaining scalar offset is zero, substitution gives

      E=L^3+(d-1)L^2-(2d-1)L,

  which has a nonzero L^2 coefficient.
  If the remaining scalar offset is nonzero isotropic, it gives

      E=L^3-(2d-1)L.

  In this last case both marked offsets are isotropic eigenvectors
  in different orthogonal eigenspaces. The rational stationarity
  data have R(s)=c and F=s-H'(c). Away from determinant roots
  there is exactly one critical point. At a determinant root
  its nonzero isotropic projected offset makes stationarity
  inconsistent. Thus there are at most one critical point
  in total; a native degeneracy already forbids unit Hessian.
  Since 2d-1>=3, the displayed profile cannot be L^3-NL.
- k=3. The remaining two eigenvalues are lambda and mu!=lambda,
  both nonzero. By (1) their curve's total genus is at least
  2d-3+epsilon. This strictly exceeds a_D for every d>=2.
  The U terms in (6) therefore cannot cancel.
- k=4. The remaining scalar offset is zero. Equations (6),(7)
  give E=L^3-dL. Here R=c and F=s-H'(c) again give at most
  one critical point, including the possibility of none.
  Since d>=2 the HG0 profile is impossible.

Thus every case of a nonzero isotropic eigenspace projection
is inverse-closed under arbitrary unit-Hessian source
reparameterization. No orthogonal split of a degenerate Krylov
space has been assumed; (5) uses a genuine nondegenerate
Witt pair in the full eigenspace.

## 4. The large nonisotropic active primary

Suppose an eigenvalue lambda has multiplicity at least three
and its projected offset is nonisotropic. In a constant frame
put

    q_A=2uw+x^2+z^2,
    g=2lambda uw+lambda x^2+mu z^2+2a x+2b z+c,
    a!=0.

For multiplicity three, mu!=lambda and mu!=0. For multiplicity
four, mu=lambda and b=0. On the whole fibre h=kappa, t=g gives

    x(t,z)=[t+lambda H(t)-lambda kappa
                         -(mu-lambda)z^2-2b z-c]/(2a),
    2uw=S(t,z):=kappa-H(t)-x(t,z)^2-z^2.           (8)

This is a global isomorphism of the whole fibre with the
hypersurface 2uw=S(t,z), with polynomial inverse using (8).
It is a fibre identity, not a claim of polynomial coordinate
rectification of h or Hessian invariance.
Solving for w when u!=0, and retaining u=0 in full, gives

    E(h=kappa)=L^3-L^2+L E(S=0).                  (9)

The curve S=0 is globally isomorphic to

    x^2+z^2+H(lambda x^2+mu z^2+2a x+2b z+c)=kappa.

For mu!=lambda, (1) gives positive total genus. For mu=lambda,
(3),(4) give genus d-1>0. Hence (9) cannot have HG0's
generic profile. This closes all these primaries in every d,
including the quadratic outer degree where the JH budget
could fail.

## 5. Completion of the semisimple branch

If any projected offset is nonzero isotropic, Section 3 applies.
Otherwise the offset Krylov space is the orthogonal sum of
nondegenerate marked lines, hence is nondegenerate. If every
active eigenvalue satisfies 2(d-1)>=k_lambda, JH1 excludes
unit Hessian. If an inactive eigenvalue is repeated, its
stationary critical quadric already excludes unit Hessian,
independently of other active bounds.

Any remaining failed active bound in dimension four has d=2
and an active multiplicity three or four. Its offset is
nonisotropic, so Section 4 applies. These alternatives exhaust
all invertible semisimple pencils, with every affine offset.
In each admitted case a unit-Hessian source member must
therefore have automorphic gradient, proving JI1.

This is not the stronger statement that every unit-Hessian
source orbit is empty. The simple-spectrum nonexistence
theorem JH1 remains intact; no converse existence statement
is inferred for the larger family.

## 6. Boundary and next structural questions

Invertibility is used in the plane lemmas and in excluding
zero-eigenvalue infinity configurations. Diagonalizability
is used for the actual orthogonal eigenspace frames.
Neither condition is silently removed.

There is a concrete reason to keep singular pencils separate.
An exceptional rank-one plane potential can reduce to

    f(z,w)=w(z^2-1)+alpha z.

This occurs inside the original affine-quadratic family, not
merely in an unrelated example. Take

    h0=2xy+z^2+w^2+H(z^2+2a z+2w+a^2),
    H(t)=-t^2/4+t/2.

The actual source automorphism Z=z+a,
W=-w-Z^2/4+5/4 has polynomial inverse
z=Z-a, w=-W-Z^2/4+5/4. Substitution gives

    h0=2xy+W(Z^2-1)-2a Z+a^2+5/4.

Here T is semisimple but singular, of rank one.
Then 2xy+f has two reduced critical points with the same
nonzero critical Hessian determinant. Its generic plane curve
is A1 minus two points, so its whole generic four-variable
profile is exactly L^3-2L. It passes the critical-count and
HG0 tests; no unit-Hessian source member or Keller failure
is asserted. This is a residual source-orbit gate, not a
counterexample or a new theorem of nonexistence.

The zero-radical shear suggested previously remains available
as an exploratory reduction, but does not resolve this gate.
Nonsemisimple independently centered pencils must likewise
retain their actual multiple-pole boundary and infinity terms.
Both actual endpoint problems, full w-dependent global
recovery, different image geometry and unrestricted HC4/JC2
remain open.

For related suspension context we inspected the primary
[Asok--Dubouloz--Ostvaer record](https://arxiv.org/abs/2112.08241v2).
No A1-homotopy, cancellation, connectedness or suspension
classification theorem is imported. Equations (6)--(9) are
proved by the explicit product strata above.
