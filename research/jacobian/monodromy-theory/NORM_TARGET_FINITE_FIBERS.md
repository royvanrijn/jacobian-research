# The cofactor norm target is not unirational

This completes the finite-fiber gate left by KMFI91. The norm double
cover has a nonzero regular three-form on a smooth projective model.
Consequently no dominant rational-target base change of the rational
m12 cofactor-marked Klein family can select an eleven-sheet marked
component with full PSL2(11) root action.

This closes that specified descent architecture. It does **not**
exclude the exceptional action for all polynomial Keller maps.

Use the notation and normalization of
[KMFI90--91](NORM_TARGET_REDUCTION.md): f is monic of degree11,
f(0)=0, f'(0)=1, the finite branch values are0,beta with beta!=0,
and their ramification profiles are \(2^4 1^3\) and \(3^3 1^2\).
The second value need not be1. Put
\[
 Q_h(T)=h+rT^2+q(T^3-cT^2),\qquad
 R_h=\operatorname{Res}_T(f(T)-b,Q_h(T)).
\]
The norm field E is k(b,q,r,Y), with \(Y^2=R_1\).
Its proper normal double-cover model has branch bidegree(4,12).
The candidate is
\[
                    \eta=q^3\,db\wedge dq\wedge dr/Y.       \tag{1}
\]

**KMFI92.** The form (1) extends to a nonzero regular three-form on
a smooth projective model of E. Its geometric genus is exactly3,
with basis given by numerators q^3,q^2*r,q*r^2 in (1).
Hence E is not unirational and is not stably rational.

**KMFI93.** In the rational relative model
\(a^{12}=f'(t)Q_1(t),\,f(t)=b\), no field extension of
k(b,q,r) that is rational over k can select an eleven-sheet marked
component with the full exceptional PSL2(11) root action.
This includes arbitrary additional target variables and birational
coordinate models of dominant base changes. It covers both A5 actions.
Non-dominant specializations, forgetting the mark, and unrelated
Keller architectures remain outside this statement.

KMFI91 already settles infinity and every ordinary fiber. Below we
prove extension at the two finite critical fibers, including their
coalescences, the h=0 contact, and the q=0 boundary.

We use the canonicality of compound Du Val threefold singularities
from [Corti--Haskins--Nordstrom--Pacini, Definition3.25 and
Proposition3.26](https://msp.org/gt/2013/17-4/gt-v17-n4-p02-p.pdf).
The finite covers and weighted-chart calculations are given below.

[Exact checker](verify_norm_finite_fibers.py) ·
[model](norm_finite_fibers_model.json) ·
[receipt](norm_finite_fibers_receipt.json) ·
[retained primary inputs](sources/norm-finite-fibers-manifest.json).

## 1. A finite cover that removes the cofactor-root collisions

On the projective coefficient chart q=1 write
\[
 Q_h(T)=T^3+(r-c)T^2+h.
\]
Order its three roots u,v,w. They satisfy uv+uw+vw=0.
The quadratic cone of ordered roots has the finite double cover
\[
 u=s z,\qquad v=s(s-z),\qquad w=-z(s-z).                    \tag{2}
\]
Indeed \(s^2=u+v,\ z^2=u+w,\ sz=u\), and the invariant ring for
(s,z)->(-s,-z) is exactly this cone. The ordered roots are integral
over the coefficient ring, so (2) is finite over the q=1 plane.
Set
\[
 D=sz(s-z),\quad B=(2z-s)(2s-z)(s+z).
\]
Then
\[
 h=D^2,\qquad r=c-(s^2-sz+z^2),\qquad
 \det\frac{\partial(h,r)}{\partial(s,z)}=-2DB.              \tag{3}
\]

Now take the normal integral complete intersection V
\[
 y_1^2=f(u)-b,\qquad y_2^2=f(v)-b,\qquad y_3^2=f(w)-b.
                                                               \tag{4}
\]
It is finite over A3_(b,s,z). The three branch primes f(u_i)-b
are distinct, and each is smooth at its generic point. Thus their
square classes are independent, giving an integral cover; (4) is
Cohen--Macaulay and regular in codimension1, hence normal.

In the q=1 chart of E the equation is \(Z^2=hR_h\).
Since deg f=11 and the cofactor is monic cubic,
\(R_h=-\prod_i(f(u_i)-b)\). Consequently
\[
                         Z=\sqrt{-1}\,D\,y_1y_2y_3
\]
defines a finite dominant map V -> E. The q-chart version of (1) is
\(-db\wedge dh\wedge dr/Z\), so (3) gives
\[
          \eta|_V=\text{nonzero constant}\,
               B\,\Omega_V,\qquad
          \Omega_V=\frac{db\wedge ds\wedge dz}{y_1y_2y_3}. \tag{5}
\]
Up to a nonzero constant, Omega_V is the global residue generator
of the dualizing sheaf of the complete intersection (4).
The multiplier B is a polynomial. In particular the root-ordering
cover has introduced no poles in (5).

A pole of a top form remains a pole after a finite dominant cover:
at an extension of a divisorial valuation its order becomes
e*a+(e-1), negative when the original integer a<=-1.
Thus it suffices to extend (5) on resolutions at the points below.
This test includes exceptional divisors, not only generic fiber points.

## 2. The origin and the two-zero-root axes

At s=z=0 all three roots in (2) are0.
For b!=0 all y_i are units and V is smooth.
For b=0 eliminate b from (4). The tangent cone consists of
two quadrics in (s,z,y_1,y_2,y_3):
\[
 y_1^2-y_3^2-2sz+z^2=0,\qquad
 y_2^2-y_3^2-s^2+z^2=0.                                  \tag{6}
\]
The determinant of their pencil is, up to a nonzero scalar,
\[
             \lambda\mu(\lambda+\mu)
                         (\lambda^2+\lambda\mu+\mu^2).
\]
All five roots are distinct. The projective intersection is therefore
a smooth degree4 del Pezzo surface: diagonalizing a nonsingular member
of the pencil gives five distinct eigenvalues, which rules out any
singular base point.

The ordinary blowup of the complete intersection vertex has
discrepancy \(5-1-2-2=0\), computed directly by substituting the
exceptional parameter in its residue form. Its exceptional surface
is smooth, so the strict transform is smooth along it.
Hence Omega_V extends at this point. Higher terms have order at
least4, because f'(0)=1 and u,v,w are quadratic.

At a nonzero point of the cone where one root is zero, exactly two
are zero. After relabeling they are u,v, and w is a unit.
The cone is smooth here; use u,w as local coordinates and
\[
                 v=-uw/(u+w)=-u+O(u^2).
\]
For b=0, the first two equations in (4) give
\[
 y_1^2-y_2^2=f(u)-f(v)=2u+O(u^2).
\]
Eliminate u and then b. To quadratic order
\(b=-(y_1^2+y_2^2)/2\).
The remaining equation is a hypersurface
\[
                  y_3^2=f(w)-b(y_1,y_2,w).
\]
It is smooth when f'(w)!=0, or has quadratic rank3 when w is
critical. In the latter case it is compound A with order2 or3
in the remaining w direction, hence canonical.
For b!=0 the same equations have unit y_1,y_2 locally and cause
no new pole; if the third root is critical, the same implicit
elimination has a smooth b direction.

This handles the entire h=0 axis with two zero roots, and its
intersection with the triple-zero point. No special choice excluding
three order3 critical roots is needed.

## 3. Nonzero cofactor roots: all critical configurations

Away from the cone vertex the cover (2) is etale over the ordered-root
cone. If all roots are nonzero, any two are local coordinates:
their sum is nonzero by uv+uw+vw=0, and
\(w=-uv/(u+v)\) has both partial derivatives nonzero.

At a point of the finite fiber let I be the indices for which
y_i=0 and u_i is critical for f. Every such critical point has
the same local ramification order e, equal to2 over0 or3 over beta.
The following cases exhaust the configurations.

If I is empty, the local branch factors on E come only from
unramified roots. There are at most three distinct branch lines,
with pairwise distinct tangents; h is a unit here.
The A1/D4 compound-Du-Val argument of KMFI91 applies even though
other roots of f in this fiber may be critical.

If |I|=1, the Jacobian of (4) has full rank. When all three y_i
vanish, use the other two, unramified root coordinates: their
derivatives together with db are independent. If fewer vanish,
the unit y_i derivatives supply the missing rows. Thus V is smooth.

If |I|=2 and the third y_i is a unit, eliminate it and b.
The resulting hypersurface is
\[
                  y_1^2-y_2^2=f(u_1)-f(u_2).
\]
The slice fixing u_2 at its critical value is A_(e-1), so this
threefold is cDV and canonical.

If |I|=2 and the third y_i also vanishes, the third root is
unramified. Eliminate its coordinate using that f' there is a unit,
and eliminate b. On the slice y_3=0 the remaining surface equation
has a nondegenerate quadratic part in y_1,y_2.
The formal Morse lemma with the remaining root coordinate x as
parameter puts it in the form
\[
                            Y_1Y_2=H(x).
\]
Here H is not identically zero. Otherwise the zero critical locus
would give a varying triple (possibly with repeated labels) of
roots of f(T)-b satisfying
\(1/u_1+1/u_2+1/u_3=0\).
Averaging its Galois conjugates over the transitive root group
would give
\[
              (3/11)\sum_i1/t_i=0,\qquad
                           \sum_i1/t_i=1/b,
\]
a contradiction. Thus H has finite order and the slice is A_k
for a finite k. The threefold is again cDV and canonical.

The formal Morse step can also be read directly: after the implicit
elimination, the equation is even in y_1,y_2, with nonzero opposite
quadratic coefficients. Its critical section in these variables
is y_1=y_2=0. If its value H vanished identically, the generic root
triple just excluded would result.

It remains to treat |I|=3, including repeated cofactor roots.
This is the purpose of the next explicit crepant calculation.

## 4. The three-critical-point complete intersection

Near three nonzero critical roots, choose x_1,x_2 as deviations
of the first two roots. The third deviation is
\[
 x_3=\alpha x_1+\beta_0 x_2+O((x_1,x_2)^2),
                    \qquad \alpha\beta_0\ne0 .
\]
Here beta_0 is a local slope, distinct from the branch value beta.
After eliminating b the leading equations of V are
\[
\begin{split}
 y_1^2-y_3^2&=A x_1^e-C(\alpha x_1+\beta_0x_2)^e,\\
 y_2^2-y_3^2&=B_0 x_2^e-C(\alpha x_1+\beta_0x_2)^e ,
\end{split}                                                   \tag{7}
\]
with A,B_0,C nonzero. All omitted terms have higher weight.

For e=2 use weights(1,1,1,1,1), equation weights(2,2).
For e=3 use weights(2,2,3,3,3), equation weights(6,6).
In both cases the residue exponent along the exceptional parameter is
\[
                \sum w_i-1-\sum\deg_w(F_i)=0.                \tag{8}
\]

We check every part of the exceptional weighted surface.
On its index cover, a singular point with a nonzero x coordinate
must have y_1=y_2=y_3=0.
Indeed a dependence (lambda,mu) of the two Jacobian rows first gives
lambda*y_1=mu*y_2=(lambda+mu)*y_3=0.
If one of lambda,mu,lambda+mu were zero, the x derivatives and
(7) would force every coordinate zero, impossible projectively.

Normalize x_1=1 and put t=x_2. A remaining singular point is
therefore a common zero of
\[
 F_1(t)=C(\alpha+\beta_0t)^e-A,\quad
 F_2(t)=C(\alpha+\beta_0t)^e-B_0t^e.
\]
At such a zero t and alpha+beta_0*t are nonzero. The slopes a,b
satisfy
\[
 a=eC\beta_0(\alpha+\beta_0t)^{e-1}\ne0,\quad
 b=-eC\alpha(\alpha+\beta_0t)^{e-1}/t\ne0,\quad
 a-b=eB_0t^{e-1}\ne0 .
\]
Eliminating t gives a nondegenerate quadratic surface equation
with coefficients proportional to b,-a,a-b.
Thus every such exceptional surface point is A1. The total strict
transform is a cDV deformation of that section and is canonical.

For e=2 these are ordinary projective charts.
For e=3 the x charts have finite index2 covers; the same argument
holds on the cover and the invariant pulled-back residue form
descends without acquiring a pole.
The remaining y charts have x_1=x_2=0 and all y_i nonzero.
Their index3 cover is smooth: two y ratios are eliminated by (7),
leaving the exceptional parameter and x_1,x_2 as coordinates.
The residual action has weights(1,1,1) modulo3.
The residue form in (8) is invariant and regular, so it extends
after taking a resolution of the quotient, by the finite-cover
pole test used above.

These charts cover the weighted exceptional surface. Equation(8)
and canonicality on the index covers prove extension of the residue
generator, including all further exceptional divisors.
Polynomial multiplication by B in (5) preserves regularity.
This finishes the q=1 chart at both finite branch fibers.

This calculation also corrects an overcautious exploratory restriction:
three order3 critical roots do not need to be excluded by a choice
of the tail normalized to0. Clearing the half-integral weights
includes the Jacobian of the additional quadratic parameter change;
selected fractional epsilon orders alone did not decide extension.

## 5. The remaining q=0 charts

If q=0, h*r!=0, the cofactor has two distinct finite roots u,-u
and one simple root at infinity. An etale ordering of the binary
cubic makes the two finite roots u,v independent local coordinates,
including at u+v=0. The infinite root contributes a unit to the
resultant: the factor q^11 cancels the pole of f(w).
After adjoining a square root of this unit the norm double cover is
\[
                    Z^2=(f(u)-b)(f(v)-b).
\]
Putting B=b-(f(u)+f(v))/2 gives
\[
                 (Z-B)(Z+B)=-(f(u)-f(v))^2/4 .              \tag{9}
\]
The right side is not identically zero in the independent coordinates
u,v. A general linear surface slice makes it a nonzero one-variable
power times a unit. Thus (9) is compound A and canonical.
The dualizing section (1) extends.

At q=h=0, use r=1. If b!=0 the cover is smooth.
For b=0 use the surface section q=0.
Write \(f(T)=E(T^2)+T O(T^2)\), with O(0)=1. Then exactly
\[
 R_h|_{q=0}=(b-E(-h))^2+h O(-h)^2 .
\]
Let H=h O(-h)^2, B=b-E(-h), W=O(-h)Z. The section becomes
\[
                         W^2=H B^2+H^2,
\]
or uv=-B^4/4 after completing the square. It is A3.
The total threefold is therefore cDV and canonical.
Finally at q=r=0, h!=0, the cofactor is constant and the resultant
is a nonzero constant; the cover is smooth.

The q=1 chart and these q=0 cases exhaust the entire projective
cofactor plane. Together with KMFI91 at infinity, they prove
KMFI92 globally.

## 6. Exact geometric genus and consequence for the proposed realization

The finite-fiber proof extends the whole dualizing sheaf there.
On q=1 the form (1) is a local residue generator, so any polynomial
multiplier in h,r remains regular after the finite cover.
The q=0 charts were shown canonical. Thus all global dualizing
sections pass the finite fibers.

The proper double cover is normal and Gorenstein, with covering
line bundle O(2,6) over P1 times P2. Projection of its dualizing sheaf
has summands O(0,3) and O(-2,-3). Therefore every regular top form on
a resolution is among the ten cubic-numerator sections
H(h,q,r)*db wedge dq wedge dr/Y. A regular form on a resolution
restricts to the smooth locus and extends reflexively to this
dualizing sheaf, so no additional forms are missed.

At infinity the scaling of KMFI91 gives a monomial q^a*r^d*h^j
the exact exceptional exponent3a+2d-7, with a+d+j=3.
Every nonzero weight piece remains nonzero in the function field
of the limiting double-cover surface. A negative exponent is
therefore an actual pole at the generic point of the new infinity
fiber, and cannot occur for a regular form.
The only surviving monomials are q^3,q^2*r,q*r^2, of exponents2,1,0.
Their transformed numerators are regular sections on the entire
projective fiber as in KMFI91. They all extend and are linearly
independent. The geometric genus is exactly3.

The form is nonzero in the generic field, so its extension is nonzero.
A unirational smooth projective model would admit a dominant rational
map from projective space, pulling it to a nonzero regular three-form
on a resolution of projective space. This is impossible in
characteristic zero. Independent variables cannot repair the defect:
a rational field containing E would itself give such a dominant map.
Hence E is neither unirational nor stably rational.

By KMFI90 every dominant base change selecting an eleven-sheet
marked component with full PSL2(11) action contains E.
It therefore cannot have a rational, ordinary-affine target.
This proves KMFI93 for both exceptional actions and every normalized
Klein input in the stated family.

The original unit obstruction concerned the literal affine component.
The new obstruction concerns the target field and survives arbitrary
birational source/target models and additional independent variables
within this dominant descent. It does not cover non-dominant
specializations or an unrelated extension with the same action.
An absolute exceptional-degree11 realization or obstruction, and
the full action spectrum, remain open.

## Verification and assurance

Run

    .venv/bin/python research/jacobian/monodromy-theory/verify_norm_finite_fibers.py

The written proof exhausts local configurations, both weighted
blowups and all projective coefficient charts. Exact controls replay
the finite cover, residue multiplier, origin pencil, crepant weights,
exceptional slope identities, axes eliminations, and A3 boundary
normal form. They supplement the universal arguments; they are not
a bounded search substituted for the case exhaustion.
The cDV-to-canonical and tame differential inputs are pinned retained
primary sources. Formal verification, independent whole-proof replay
and external review are not claimed.
