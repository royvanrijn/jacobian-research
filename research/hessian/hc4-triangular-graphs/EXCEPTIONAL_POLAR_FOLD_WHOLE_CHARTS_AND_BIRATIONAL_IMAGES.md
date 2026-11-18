# Whole charts and birational images of exceptional polar folds

An omitted ramified prime with source values (-2,1) and normal
ramification index two has an entire plane chart in the actual finite
normalization. Its target image has normalization A1 and residue degree
one. These conclusions use the actual finite target map; the linear
Jacobian of an abstract chart alone is insufficient.

The chart yields a concrete realization constraint for the exceptional
type (p,r,n)=(2,1,5) in
[the ramified-fold theorem](RAMIFIED_FOLD_VALUATIONS_AND_COMMON_COORDINATE_POLES.md).
It does not exclude that type. A smooth affine geometric control below
realizes the entire chart and the regular primitive, while a different
family realizes finite polynomial maps with linear Jacobian in every
generic degree at least two. Neither supplies the simultaneous source
polynomial admission required of an actual Keller pair.

## 1. Exact hypotheses and conclusion

Let (P,Q):X=A2_(a,b)->Y=A2 be polynomial with constant
J=Jac_(a,b)(P,Q) in C*. Put B=C[P,Q], L=C(a,b), and let S be the
FULL finite normalization of B in L, with Z=Spec S. The marked source
X embeds openly in Z. Suppose an omitted prime E of Z has
```
e_E=2,       v_E(a)=-2,       v_E(b)=1.                  (1)
```
Here a,b are polynomial coordinates of the given source, after shifting
b by its constant residue. The case of interest has a equal, up to
affine change, to the common coordinate of the unramified completion.
Once(1) is given, the present theorem does NOT require that completion
to have a nonempty boundary.

After an affine change of a and an admissible residue parameter t,
the valuation has the exact field chart
```
s=b,       a=s^-2+s*t,       L=C(s,t).                   (2)
```
For this normalized source pair the target Jacobian is still a nonzero
constant, denoted J. Substitution in(2) gives polynomials
```
P_E(t,s), Q_E(t,s) in C[t,s],
Jac_(t,s)(P_E,Q_E)=J*s.                                 (3)
```
Moreover every element of S lies in C[t,s], and the resulting morphism
```
phi:A2_(t,s)->Z
```
is an OPEN IMMERSION. Its image is exactly
```
(X minus {b=0}) union E,
```
and its line s=0 maps isomorphically onto the ENTIRE E. In particular
E is smooth A1 and meets no other omitted prime of Z.

Writing Gamma for E's target image, its normalization is A1 and
```
[C(E):C(Gamma)]=1.                                     (4)
```
The target curve need not be smooth or embedded as an affine line.
At any normalization parameter t0 where its polynomial parameterization
is nonimmersive, suitable affine linear target coordinates on that
individual branch have leading orders
```
m, m+1,       m>=2.                                    (5)
```
Other normalization parameters may map to the same target point.

## 2. The five point blowups give an exact chart

The canonical calculation in the
[full-normalization proof](CANONICAL_BOUNDARY_AND_STABLE_COMPLETE_INTERSECTION_OBSTRUCTION.md)
gives ord_E(da wedge db)=e_E-1=1. Work after the first source
infinity blowup, at the finite-b point with local coordinates
u=1/a,b. The local point multiplicities satisfy
```
sum m_j=2*2+(2-1)=5.
```
Their first value is min(2,1)=1 and they never increase, so the complete
sequence consists of five ones. This conclusion does not require the
unramified anchor inequality.

Put s=b. The five blowup coordinates near the generic final divisor
can be written exactly as
```
u=s^2*(kappa+A*s+B0*s^2+s^3*T),       kappa!=0,          (6)
```
where A,B0,kappa are constants and T is a separating affine parameter
on the final exceptional divisor. Indeed the first two substitutions
are u=s*u1, u1=s*u2; one then successively writes
u2=kappa+s*v, v=A+s*w, w=B0+s*T.
Thus
```
a=kappa^-1*s^-2-A*kappa^-2*s^-1
  +(A^2*kappa^-3-B0*kappa^-2)+O(s).
```

The ACTUAL target form theta=J^-1*P*dQ is regular at the generic point
of E and has derivative da wedge db. On the source polynomial plane
the polynomial Poincare lemma gives
```
theta=a*db+dF,       F in C[a,b].
```
On a formal transverse arc the derivative dF has no ds/s coefficient,
whereas a*db has coefficient -A*kappa^-2. Consequently A=0.

Put c=kappa^-1 and alpha=-B0*kappa^-2, and define the rational
function
```
t0=(a-alpha-c*s^-2)/s.
```
Its residue is -kappa^-2*T, so it is a residue-field generator.
The valuation of C(s,t0) is exactly the s-adic valuation, with t0
transcendental in the residue field. The identity
a=alpha+c*s^-2+s*t0 is exact, not an infinite formal approximation.
Replacing a by (a-alpha)/c and t0 by t=t0/c gives(2).

In particular, for every f in C[a,b], the substitution is a Laurent
polynomial in s with polynomial t-coefficients, and
```
v_E(f)>=0  iff  f(s^-2+s*t,s) belongs to C[t,s].          (7)
```
The residue generator ensures no unexplained cancellation between
different t-polynomials. No whole ramified component chart has been
assumed in deriving(7).

## 3. Finite target compatibility supplies the WHOLE plane

Actual P,Q belong to S, so they satisfy(7). This proves their
polynomiality in(3); the chain rule gives its displayed Jacobian.

More generally, if f belongs to S, it is integral over B. Under(2)
it is integral over C[t,s], since B maps into that normal ring.
It belongs to its fraction field L, so f belongs to C[t,s].
This defines phi over Y.

On s!=0, formula(2) is an isomorphism onto the original-source
open b!=0. At the generic point of s=0, phi has centre E, by(7)
and the original valuation centre. The restrictions
```
P0(t)=P_E(t,0),       Q0(t)=Q_E(t,0)
```
are not both constant: E is a divisor of the finite surface map
Z->Y and therefore has a curve as image.

A nonconstant polynomial map from A1_t to its affine curve image
is finite. One nonconstant coordinate polynomial already makes
C[t] integral and finite over a polynomial subring, and hence finite
over C[P0,Q0]. Thus the map s=0 -> Y is finite.
Its lift s=0 -> Z is also finite: its graph in A1_t times_Y Z is
closed because Z->Y is separated, and the projection to Z is a base
change of the finite map A1_t->Y. Its image is therefore closed
and contains the generic point of E, so it is all of E.

Every fibre of phi has finitely many points: s!=0 is an open immersion
into X, and s=0 has just been shown finite over Z. Thus phi is
quasi-finite. It is birational because its function field is L.
Zariski--Main for a separated quasi-finite birational map to a normal
target now makes phi an open immersion. Since its zero line is finite
onto the whole E, it is an isomorphism there. This proves the image
description and smoothness of E.

If a different omitted prime met E, its generic point would belong
to the same open plane, since an open neighbourhood of an intersection
contains that generization. But the only points of this plane outside
X are its zero line. Thus E is disjoint from every other omitted prime.

The target-finiteness step cannot be dropped. The polynomial map
```
(P_E,Q_E)=(s*t,s),       Jac_(t,s)=s
```
contracts s=0. It has the same displayed linear Jacobian but does NOT
satisfy the nonconstant zero-line hypothesis furnished by the actual
finite-normalization divisor.

## 4. A linear Jacobian forces birationality on its critical line

Here is a separate polynomial lemma.

Suppose U(t,s),V(t,s) are polynomials,
Jac_(t,s)(U,V)=C*s with C!=0, and their restrictions U0,V0
to s=0 are not both constant. Then
```
C(t)=C(U0,V0).                                         (8)
```

First, the normalization of the affine curve image is A1. To see this
without assuming a plane embedding theorem, C[t] is finite over
C[U0,V0]. Normalize the latter in its fraction field. The resulting
normal affine curve is finitely covered by A1_t. Its smooth
projective completion is dominated by P1_t, so Riemann--Hurwitz gives
genus zero. Every point missing from the affine normalization has
inverse image contained in the single point infinity of P1_t.
There is therefore exactly one missing point; the affine normalization
is P1 minus that point, hence A1. A generator h of its coordinate ring
pulls back to a polynomial h(t). We can write
```
U0=A(h), V0=B(h),       deg(h)=[C(t):C(U0,V0)].
```

Expand U=sum s^j*Uj(t), V=sum s^j*Vj(t).
Let G=gcd(A',B'), chosen nonzero, and put
```
A'=G*I,       B'=G*K,       gcd(I,K)=1.
```
This includes a constant one of A,B, using gcd(0,f)=f.
The coefficient of s^0 in the Jacobian is
```
h'*G(h)*(I(h)*V1-U1*K(h))=0.
```
Bezout for I,K, composed with h, shows that
```
U1=k(t)*I(h),       V1=k(t)*K(h)
```
for a polynomial k(t). The coefficient of s^1 is then
```
C=k^2*h'*(I'(h)*K(h)-I(h)*K'(h))
  +2*h'*G(h)*(I(h)*V2-U2*K(h)).                          (9)
```
It is divisible by h'. Thus h' is a nonzero constant and deg(h)=1,
proving(8). This calculation is valid for arbitrarily high degrees.

Apply the lemma to(3). Since the whole E has coordinate t, (8) is
exactly residue degree one and proves(4).

For the branch refinement, now let D=gcd(P0',Q0'), and write
P0'=D*I,Q0'=D*K with gcd(I,K)=1. The same coefficient identity reads
```
J=k^2*(I'*K-I*K')+2*D*(I*Q2-P2*K).                      (10)
```
At a root t0 of D, both k and the Wronskian I'*K-I*K' are nonzero.
Choose a linear target coordinate whose primitive derivative is
nonzero at t0. Subtract its suitable scalar multiple from the other
coordinate. The remaining primitive derivative has a SIMPLE zero
because the Wronskian is nonzero. If D has order m-1 at t0, integration
gives leading orders m and m+1, proving(5).
This concerns the individual branch, not all singularities of Gamma.

## 5. Exact polynomial admission remains an additional constraint

An arbitrary polynomial pair satisfying Jac=C*s need not come from
the marked original source. Both target coordinates must ALSO satisfy
```
P_E((a*b^2-1)/b^3,b),
Q_E((a*b^2-1)/b^3,b) in C[a,b].                         (11)
```
Equations(3),(8),(11), with the whole chart in the full S, are a
concrete exceptional-fold realization problem. No argument above
discharges(11) for an arbitrary linear-Jacobian pair.

The retained finite family in
[regular source-coordinate completions, Section5.1](REGULAR_SOURCE_COORDINATE_COMPLETIONS_AND_PRIMITIVE_RESIDUES.md)
already supplies linear-Jacobian maps of unbounded generic degree.
After the equivalent sign changes, its equations can be written
```
U=s+f(t),       V=t*(s+f(t))-F(t),       F'=f.
```
The cited proof establishes finiteness and generic degree deg(f)+1,
with Jacobian -s. We use that family as a control here; we do not
repeat or claim its unrestricted finite-degree theorem.
Its critical-line image is (f,t*f-F); for f=t^m its nonimmersive
branch has leading orders (m,m+1).

These examples FAIL(11). If f has degree d>=1, the smallest b-order
in U((a*b^2-1)/b^3,b) is -3d, with nonzero leading coefficient.
Thus U is not a polynomial in the original source frame. The family
is an exact finite-target control for the linear-Jacobian lemma,
not a Keller counterexample or a realization of the original fold.

## 6. A smooth affine whole-chart control for the exceptional type

The entire marked plane geometry and exact primitive can also occur
without a supplied Keller target. Define N by the following SEVEN
equations in b,u,w,k,l:
```
b^2*w=u*(u-1),      b*k=u*w,           k^2=u*l,
w*k=b*l,           w^3=(u-1)*l,       b^2*l=u*w^2,
(u-1)*k=b*w^2.                                         (12)
```
Its source map is
```
u=a*b^2,
w=a*(a*b^2-1),
k=b*a^2*(a*b^2-1),
l=a^3*(a*b^2-1)^2.
```
On D(u-1) put a=w/(u-1);(12) gives exactly these source
formulas. This chart is C[a,b,(a*b^2-1)^-1].
On D(u) put t=k/u^2. The equations give
```
u=1+b^3*t, w=b*u*t, k=u^2*t, l=u^3*t^2,
```
so this chart is C[b,t,(1+b^3*t)^-1].
The two charts cover N, are smooth integral domains, and have a
nonempty common open. Therefore N is smooth affine and integral.

The source map identifies the WHOLE A2_(a,b) with
D(b) union D(u-1), since b and a*b^2-1 have no common zero.
Its complement is E={b=0,u=1,w=0,l=k^2}, an affine line.
Similarly D(b) union D(u) is the WHOLE plane A2_(t,s), s=b,
with transition(2); its complement is the source line b=0.
Thus the two entire planes glue along b!=0 as claimed.
In this geometry,
```
da wedge db=s*dt wedge ds,
a*db+d(a*b)=3*s*t*ds+s^2*dt.
```
Both the volume and the same primitive are globally regular on N.
They have respectively order one and no pole along E. This realizes
(p,r,n)=(2,1,5) in whole noncontracting affine plane geometry.

The last relation in(12) cannot be omitted from the displayed
presentation: its omission permits b=u=w=l=0,k=epsilon with epsilon^2=0.
The spurious tangent direction is not part of the stated integral
surface.

No finite Keller target pair is supplied for N. In particular this
control does not produce an unramified completion boundary or satisfy
the additional simultaneous-D constraints used in the common-pencil
anchor. It shows that the exceptional curve's own affine geometry,
whole chart and exact primitive alone do not exclude it.

## 7. Evidence and remaining obligation

The checker controls the five exact blowups, the primitive residue,
the exact field chart and target Jacobian chain factor; the coefficient
identity(9), its nonbirational and contracting negative controls; the
unbounded-degree finite-target family and its failure of source
admission; and all seven affine equations, both localizations, both
whole planes and the global primitive.

The written arguments supply the unrestricted geometry and curve
normalization steps. The classical inputs are normality in the fraction
field, finiteness and localization of affine normalization, the
quasi-finite birational form of Zariski--Main, elementary finite
polynomial maps of A1, Riemann--Hurwitz for smooth projective curves,
and the polynomial Poincare lemma. The existing canonical/source-open
and fold proofs retain their stated dependency limits.

The remaining realization must combine actual polynomial admission,
the finite ordinary target trace pairing and every other original
boundary prime. No all-fold exclusion, automatic common-coordinate
integrality, finite target for the geometric control, or HC4/JC2
inversion is proved. No formal verification or external review is
asserted.
