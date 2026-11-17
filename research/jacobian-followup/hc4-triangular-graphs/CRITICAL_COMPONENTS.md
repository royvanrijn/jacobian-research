# Affine critical components close nonlinear boundary-root graphs

14 November 2026. **HC4MYGR2** closes every y=R(x,p) whose boundary
R(0,p) has a simple zero: after every polynomial retained-dual
automorphism over K[x], a nonzero constant Hessian determinant forces
an invertible gradient. The polynomial R has arbitrary degree in
both variables. In many of these cases unit Hessian is impossible.

**HC4MYGR3** also forbids unit Hessian after all retained-dual
automorphisms for every opposite graph x=a(y)p+b(y) with a(0) nonzero.
This closes the earlier example with matching critical determinants;
its critical components carry an additional geometric obstruction.

The new global argument uses the reproduced Tang alternative and two
explicit classical inputs: closed affine-line rectification in the
plane, and Gwozdziewicz's line-injectivity theorem. The
[new import record](sources/critical-components-import.json) pins the
rectification statement as reproduced in
[Jelonek, Theorem 2.4](https://arxiv.org/html/2011.03472v3#S2.Thmdefi4).
The line theorem is already retained as
[HC4MYGL0](CRITICAL_COORDINATE_JET.md). No rectification in three
dimensions and no rational-curve substitute for an actual affine
line are assumed.

## 1. A plane Keller fibre with an affine-line component

**HC4MYGR0.** If a plane Keller pair (A,C) has a closed component
Gamma of C=0 isomorphic to the affine line, the pair is a polynomial
automorphism.

Take a polynomial isomorphism gamma:A1->Gamma. Its tangent never
vanishes. Since D(A,C) is invertible and C composed with gamma is
constant, the one-variable polynomial A composed with gamma has
derivative nowhere zero. Over C that derivative is a nonzero
constant, so the pair is injective on Gamma. Rectify this actual
closed embedded affine line by Abhyankar--Moh--Suzuki, and apply
HC4MYGL0 to the resulting plane Keller pair. This proves the claim.
In particular the entire C=0 fibre is then one affine line.

## 2. A four-variable geometric criterion

Let f(s,z1,z2,z3) have nonzero constant Hessian determinant. Let
Cf=V(f_z1,f_z2,f_z3), the relative critical scheme. Its three
Jacobian rows are rows of the invertible full Hessian. Thus Cf is
smooth of pure dimension one whenever nonempty.

**HC4MYGR1.** Suppose Cf has a closed irreducible component
isomorphic to A1. If either

1. some component of Cf has nonconstant projection to s; or
2. a generic three-dual fibre f(s,-)=lambda has Euler characteristic
   different from one,

then grad f is polynomially invertible. If Cf is additionally
reducible, these hypotheses forbid unit Hessian altogether.

For a hypothetical nonautomorphic unit-Hessian gradient, HCTANG2
supplies a jointly affine constant two-plane W. Its projection to
the one-dimensional free base has rank zero or one. Only constant
affine changes preserving the free-base fibration are used below.

### Rank zero

The potential has form

\[
 f=D(s,U)+P F(s,U)+Q G(s,U),\qquad
 \operatorname{Jac}(F,G)\in C^\times.             \tag{1}
\]

The common zeros of F,G are finite: a curve of common zeros would
make their differentials dependent at a smooth point. Every
relative critical point lies over such a zero (s,U). Thus Cf has
only finitely many base values, contradicting condition 1.

For every other base value, F(s,U),G(s,U) are coprime in C[U].
A Bezout completion of that row makes f(s,-) a polynomial coordinate
of the three-dual ring. All its value fibres are A2 and have Euler
characteristic one, contradicting condition 2. The completion is
used to identify a fibre variety, not to preserve a Hessian under
a nonlinear source change.

### Rank one

The normal form is

\[
 f=D(U,Z)+sA(U,Z)+VC(U,Z),\qquad
 \operatorname{Jac}(A,C)\in C^\times.             \tag{2}
\]

Its relative critical equations are C=0 and

\[
 \nabla D+s\nabla A+V\nabla C=0.
\]

The last two equations solve for s,V as polynomial functions of
U,Z, since their coefficient determinant is a field unit. Projection
therefore identifies the **whole critical scheme** with C=0 in the
two-dimensional residual plane. The affine-line component of Cf
becomes a closed embedded affine line in that plane. HC4MYGR0 makes
(A,C) an automorphism, and cotangent reconstruction makes grad f
an automorphism, contradicting the assumed nonautomorphy. This
proves HC4MYGR1's invertibility assertion.

Finally, if grad f is an automorphism, Cf is its inverse image of the
target coordinate line where the three dual outputs vanish. Hence
Cf is itself an irreducible affine line. Reducibility contradicts
this and proves the stated no-unit refinement. Invertibility alone
does not assert unit-Hessian nonexistence without that refinement.

All component hypotheses transport through actual retained-dual
automorphisms by [HC4MYGK0](RELATIVE_CRITICAL_LOCUS.md). A punctured
line, a normalization isomorphic to A1, or one branch at infinity
does not meet the closed affine-line hypothesis.

## 3. Every y=R(x,p) has a noncoordinate generic dual fibre

Assume deg_p R>0. For generic fixed x=s nonzero, the change
(p,q,r)->(p,T=A,q) is a polynomial automorphism of the three-dual
fibre: the factor 1/s^3 is now a field constant. The exact formula
from HC4MYGK1 becomes

\[
 f_s=L T^2+c(p)T+e(p)+d(p)q,\qquad
 d(p)=2N(3+2sR(s,p))/s.                          \tag{3}
\]

Complete the square in T. Let r>=1 be the number of distinct roots
of d(p); multiplicities are unrestricted. For a generic value lambda,
projection of f_s=lambda to p has the following two disjoint strata.

- Off those r roots, q is uniquely determined and T is free. This
  stratum is (A1 minus r points) times A1 and has Euler characteristic
  1-r.
- At each root, the nonzero quadratic coefficient L gives two distinct
  T values for generic lambda, and q is free. These are two affine
  lines, contributing two per root.

Use compactly supported complex Euler characteristic. Additivity and
multiplicativity give the exact formula

\[
 \chi(f_s^{-1}(\lambda))=1+r\geq2.               \tag{4}
\]

No rationality, irreducibility or degree bound is inferred from this
Euler value. It is enough to exclude the coordinate fibres in rank
zero of HC4MYGR1. All retained-dual automorphisms specialize to
polynomial automorphisms of these fibres and preserve (4).

## 4. A simple boundary zero supplies an actual affine critical line

Put R0(p)=R(0,p) and H(p)=Lp^2+Mp. At x=0,

\[
 f_0=H(p)+4NpR_0^2+NqR_0,
\]
\[
 (f_p,f_q,f_r)|_{x=0}
  =\bigl(H'(p)+4NR_0^2+N(8pR_0+q)R_0',\ NR_0,\ 0\bigr). \tag{5}
\]

For a simple root rho of R0, set v=R0'(rho) nonzero. The literal
closed line

\[
 x=0,\quad p=\rho,\quad q=-H'(\rho)/(Nv),\quad r\text{ free} \tag{6}
\]

is relatively critical. On it the native full Hessian determinant is
4*N^4*v^2, so the critical scheme is smooth of dimension one there;
the line is an irreducible component, not just a subset of a larger
critical surface.

Equations (4) and (6) meet condition 2 of HC4MYGR1. This proves
**HC4MYGR2**: every y=R(x,p) with a simple zero of R(0,p) has an
invertible gradient whenever its pullback under any retained-three-dual
polynomial automorphism over K[x] has unit Hessian. Degrees in both
x and p, coefficient sizes and parameter word lengths are unrestricted.

There are useful no-unit refinements. Two simple boundary roots give
two disjoint critical-line components. Also, any nonempty component
of 2xR+3=0 on x*R_p nonzero gives the additional critical curve
classified in HC4MYGK1. Its projection to x is nonconstant and it
cannot meet (6), since 2xR+3 equals 3 at x=0. Either circumstance
makes the critical scheme reducible, so HC4MYGR1 forbids unit Hessian.

A multiple root rho with H'(rho)=0 gives a different direct obstruction:
(5) vanishes on the whole two-plane x=0,p=rho, with q,r free.
This cannot occur in the smooth one-dimensional relative critical
scheme of a unit-Hessian potential, even after a retained-dual
automorphism. Thus this case also has no unit-Hessian pullback.

The case of an entirely repeated generic critical curve is retained.
For example R=3p+(3/2)xp^2 has
2xR+3=3(1+xp)^2 and R_p=3(1+xp), so the open critical curve
x*R_p nonzero is empty. Nevertheless R0=3p gives the boundary line,
and the generic dual fibre in (4) has Euler characteristic two. The
invertibility conclusion still applies. We do not assert no unit
Hessian from a second component that is absent in this example.

If a nonautomorphic unit-Hessian y=R(x,p) remains outside these
exclusions, R0 must be nonzero, have no simple root, and satisfy
R0(-M/(2L)) nonzero. Nonzero constant boundary traces and polynomials
whose roots are all multiple and avoid that critical p value are
not closed by this theorem. The identically zero trace was already
excluded by HC4MYGU0.

## 5. The opposite affine-p chart: a(0) nonzero is now excluded

Take x=a(y)p+b(y), with a(0)=a0 nonzero. Put b0=b(0),
rho=-b0/a0 and h0=H'(rho). At y=0 there is the closed critical line

\[
 p=\rho,\qquad r=-h_0(1+3a_0q)/(2Na_0),\qquad q\text{ free}. \tag{7}
\]

The native full Hessian determinant along it is 4*N^4*a0^2,
independently of all higher coefficients of a,b. Thus it is a
smooth irreducible critical component.

There is also a critical curve over the nonempty open y*a(y) nonzero.
Use x=-3/(2y), T=(13N/(2x^2)-M)/(2L), and
p=(x-b(y))/a(y). In the rational computation coordinates (x,T,q),
the q coefficient is 6N/x+4Ny. Its x derivative is -6N/x^2,
nonzero there, so the remaining x derivative uniquely determines q;
r is recovered from A. This curve projects nonconstantly to y and
cannot meet the boundary line at y=0 in the affine source.

Consequently the relative critical scheme is reducible, contains a
closed affine-line component, and has a component dominating the
free base. HC4MYGR1, now using condition 1 rather than Euler, forbids
unit Hessian after every retained-three-dual automorphism over K[y].
This proves **HC4MYGR3** in every degree.

In particular the earlier graph x=p+(1013/54)y^2, with L=N=1,M=0,
is now excluded even though its computed critical determinants are
all 4. The line-plus-dominating-curve geometry supplies the obstruction
that the determinant-value test did not provide.

If a(0)=0 and b(0)=0, HC4MYGU0 already excludes the whole zero slice.
The affine-p opposite-chart case with a(0)=0,b(0) nonzero remains
outside this new result; it includes the unresolved base-only endpoint.

## Reproduction, imports and remaining scope

The proofs are first over C. A characteristic-zero unit-Hessian
counterexample over another coefficient field would survive a
finite-coefficient-field embedding into C; polynomial inverses descend
uniquely in the invertible cases. Boundary-root conditions are taken
over the algebraic closure of that coefficient field.

The [checker](verify_critical_components.py) verifies both boundary
critical-line equations and full determinants, both projection normal
forms, the whole rank-one critical-scheme reconstruction, Euler
stratum counts and repeated-root controls, and the opposite-chart
dominating curve and prior matching-determinant example. It pins the
actual Tang, rectification and line-theorem sources. The all-degree
exhaustion and topology argument are written proofs, not inferred
from finite checks. [Actual receipt](receipt-critical-components.json).

```sh
timeout 60s .venv/bin/python research/jacobian-followup/hc4-triangular-graphs/verify_critical_components.py
```

The no-simple-root boundary cases above, q/r-dependent omitted-base
graphs, opposite-chart cases with a(0)=0,b(0) nonzero, nonlinear
free-base/dual mixing, the unit-ideal plane-surjectivity endpoint and
the separate singleton-zero properness endpoint remain unresolved.
No unrestricted HC4 or JC2 solution, formal verification, independent
external review or novelty priority is claimed.
