# Cubic two-pair Image positivity and the sharp balanced degree four

Work over a field of characteristic zero, with contraction pairs
\((\xi_1,z_1),(\xi_2,z_2)\), and
\(\mathcal E_2(\xi^\alpha z^\beta)=\partial_z^\alpha z^\beta\).

**Theorem (complete cubic stratum).** For every bihomogeneous polynomial \(F\) of
bidegree \((3,3)\),
\[
 \mathcal E_2(F^m)=0\quad(m\geq1)
 \quad\Longrightarrow\quad
 \mathcal E_2(QF^m)=0\quad(m>3\deg_z Q)
 \tag{1}
\]
for every polynomial multiplier \(Q\). The coefficient three is sharp.
Together with the [quartic counterexample](TWO_PAIR_IMAGE_MATHIEU_COUNTEREXAMPLE.md),
this proves that, among bihomogeneous balanced forms in two pairs, the exact
least failing bidegree \((d,d)\) has \(d=4\).

Over \(\mathbb C\), the proof classifies the entire sixteen-parameter sphere
family: after an orthogonal linear change of coordinates and nonzero scalar
multiplication, a pure phase is strictly one-sided in circular weight or is
\[
 p_*=y+xt-x^3/2.
\]
The zero phase is included in the one-sided case. The exceptional phase is
safe and lies outside the orthogonal nullcone.

The central step is a geometric reduction to a constant-profile frame. If no
such frame exists, explicit boundary blow-ups give a smooth model proper over
the phase line on which the spherical volume form is holomorphic. A proper
\(V\)-filtration argument rules out its nonzero constant exponential period.
An exact eleven-parameter certificate then classifies the remaining frame.
This proof retains the spherical functional and does not import the Gaussian
ternary-cubic classification.

The argument and exact replays are internal. External mathematical review and
formal verification of the complete theorem have not been performed.

## The full coefficient space and its spherical functional

Set

\[
 R=\xi_1z_1+\xi_2z_2,\quad X=\xi_1z_2,\quad
 Y=2\xi_2z_1,\quad T=\xi_1z_1-\xi_2z_2,
 \qquad T^2+2XY=R^2.
\]

Restriction to R=1 identifies the balanced space V_3 with

\[
 A_{\leq3},\qquad A=\mathbb C[x,y,t]/(q-1),\qquad q=t^2+2xy.
\]

Indeed, eliminate every xy using xy=(1-t^2)/2. The resulting basis is

\[
 1,t,t^2,t^3;\quad
 x,xt,xt^2,x^2,x^2t,x^3;\quad
 y,yt,yt^2,y^2,y^2t,y^3.
 \tag{2}
\]

To see independence, order monomials so that xy is the leading term of q-1;
the monomials in (2) are distinct unreduced monomials. Every member of (2)
has a balanced-degree-three lift, by homogenizing it with R. Conversely the
sixteen original monomials
xi_1^i xi_2^(3-i) z_1^j z_2^(3-j), 0<=i,j<=3, map to this basis by an
invertible rational matrix. In the displayed basis order and lexicographic
(i,j) order its determinant is 2^(-24). Thus this identification discards no
coefficient stratum.

Let L be the normalized sphere functional. Its exact algebraic definition is

\[
 L(x^a y^b t^c)=0\quad\text{unless }a=b\text{ and }c\text{ is even},
 \qquad
 L(x^a y^a t^{2k})=\frac{a!(2k-1)!!}{(2a+2k+1)!!}.
 \tag{3}
\]

Here (-1)!!=1. Equivalently, take the phase constant term after substituting
y=(1-t^2)/(2x), and then integrate t uniformly on [-1,1]. This formula
annihilates (q-1)A and is O_3(C)-invariant. For every balanced polynomial B of
bidegree (n,n), with sphere restriction b,

\[
 \mathcal E_2(B)=(n+1)!L(b).
 \tag{4}
\]

One can derive (4) without analytic integration: the formal beta identity
\(\int_0^1 s^i(1-s)^j ds=i!j!/(i+j+1)!\), followed by phase extraction,
gives precisely the monomial contraction. Thus (1) for balanced multipliers
is exactly the eventual mixed-moment assertion for L.

Antipodal parity gives a second useful full coordinate model:

\[
 p=H+E\pmod{q-1},\qquad H\text{ homogeneous cubic},\quad
 E\text{ homogeneous quadratic}.
 \tag{5}
\]

This representation is unique: homogenize odd linear terms by q, and even
constants by q; comparing degrees in (q-1)(a+ell) proves injectivity. The
parameter dimensions are ten and six. These are two linear subspaces of the
full space, not two components whose separate classification settles their
mixed sums. In particular the summands in L((H+E)^m) have different total
degrees. Different homogeneous summands have different Gaussian-to-sphere factors;
no common factor may be pulled out of this sum.

### The fixed compact volume period

In Cartesian coordinates with \(q=u_1^2+u_2^2+u_3^2\), let
\[
 \Omega=\frac1{2\pi}\operatorname{Res}_{q=1}
       \frac{du_1\wedge du_2\wedge du_3}{q-1}.
\]
Its restriction to the oriented real unit sphere \(\Sigma\) is the
normalized area form, because the unnormalized residue is half the area
form. Thus \(L(a)=\int_\Sigma a\Omega\) and \(\int_\Sigma\Omega=1\).
The form is algebraic and holomorphic on the affine quadric. In circular
coordinates it is a nonzero constant multiple of \(dy\wedge dt/(2y)\).
For a fixed phase \(p\), the compact cycle gives an entire period
\[
 I(\tau)=\int_\Sigma e^{\tau p}\Omega
        =\sum_{m\geq0}\frac{\tau^m}{m!}L(p^m).
\]
Compactness justifies the expansion and all parameter derivatives locally
uniformly in \(\tau\). All-order purity is exactly \(I(\tau)=1\).
The seed is the original fixed form \(\Omega\), without a parameter-dependent
normalization.

### Why balanced multipliers suffice

Suppose Q has bidegree (a,b). If a>b its product with F^m contracts to zero.
Otherwise r=b-a and P_m=E_2(QF^m) is homogeneous of coordinate degree r. For
each multi-index gamma with |gamma|=r,

\[
 \mathcal E_2(\xi^\gamma QF^m)=\partial_z^\gamma P_m.
 \tag{6}
\]

The multiplier xi^gamma Q has bidegree (b,b). Vanishing of all these scalar
moments forces every coefficient of P_m to vanish, in characteristic zero.
A balanced multiplier of degree b restricts to a sphere polynomial of degree
at most b. Hence a sphere bound m>C deg g implies the original Image bound
m>C deg_z Q. Split an arbitrary Q into its finitely many bidegrees to obtain
the same statement for arbitrary, unbalanced multipliers.

## An exact null-conic frame test

For the canonical lift p=H+E in (5), consider a nonzero null vector v. There
exists a cubic representative of p which is affine linear in the coordinate
in direction v if and only if

\[
 q(v)=0,\qquad E(v)=0,\qquad \nabla_q H(v)\in\mathbb C v.
 \tag{7}
\]

Here the gradient uses the bilinear form of q. Representatives of degree at
most three differ by (q-1)(a+ell), with a constant a and a linear form ell.
At a null vector, adding q ell changes the gradient of H by 2 ell(v)v and
adding a q does not change E(v). Consequently (7) permits a representative
whose cubic gradient vanishes at v and whose quadratic part vanishes there.
For a cubic polynomial P with homogeneous pieces P_3,P_2,

\[
 D_v^2P(z)=2B(\nabla_q P_3(v),z)+2P_2(v).
\]

This is zero exactly when P is affine linear in the v coordinate. The
converse proves necessity of (7) for such a representative.

This criterion is effective without choosing an orthogonal frame. Parametrize
the null conic by v(s,r)=(s^2,-r^2/2,sr) and form the binary sextic
f(s,r)=H(v(s,r)) and binary quartic g(s,r)=E(v(s,r)). Then (7) is equivalent
to a common projective zero of f_s,f_r,g. Euler's identity gives f=0 at such
a zero, and the derivative condition says that the zero of f on the conic
is repeated. If f is identically zero, retain both zero derivatives and
apply the same criterion; no division by f or by a discriminant is allowed.

## A proper holomorphic volume seed cannot have constant exponential period

**Lemma.** Let \(Y\) be a smooth irreducible complex algebraic variety,
let \(f:Y\to\mathbb A^1_t\) be proper and nonconstant, and let
\(\omega\) be a global algebraic holomorphic top form on \(Y\).
For a fixed compact middle-dimensional cycle \(\gamma\subset Y\),
the entire function
\[
 I(\tau)=\int_\gamma e^{-\tau f}\omega
\]
cannot be a nonzero constant.

The statement does not require isolated critical points. Properness is
relative to the affine target; \(Y\) need not be proper over a point.
Nonconstancy on irreducible \(Y\) is essential.

**Proof.** Use the decreasing Kashiwara--Malgrange filtration
\(V^\alpha\) for which \(\partial_t t-\alpha\) acts nilpotently on
\(\operatorname{gr}_V^\alpha\). Factor \(f\) through its graph
\[
 i:Y\hookrightarrow Y\times\mathbb A^1_t,\qquad
 \pi:Y\times\mathbb A^1_t\to\mathbb A^1_t.
\]
Put \(\mathcal B_f=i_+\mathcal O_Y\). Negativity of the local
Bernstein--Sato roots implies
\[
 \delta(t-f)\in V^{>0}\mathcal B_f.
 \tag{P1}
\]
This holds for arbitrary nonzero local \(f\), including non-isolated
critical loci. Multiplying by a regular coefficient of \(\omega\)
preserves the same positive filtration bound.

Let \(M=\mathcal H^0f_+\mathcal O_Y\) and let \(m\in M\) be the class
of the top form \(\omega\). Near \(t=0\), one may choose a single
\(\alpha>0\) for the relevant top-form section in (P1). Indeed the
fiber \(f^{-1}(0)\) is compact; finitely many local bounds suffice,
and their positive minimum is positive. Properness permits shrinking
the target so that no part of the source over that neighborhood is
missed by this finite cover. Off \(t=0\) the \(V\)-filtration imposes
no restriction. If the fiber is empty, the same conclusion is
immediate.

The projection \(\pi\) is proper on the support of \(\mathcal B_f\).
The proper push-forward theorem for \(V\)-filtrations identifies the
filtration on \(M\) with the filtration induced from the relative
de Rham complex. Its degree-zero top-form class therefore satisfies
\[
 m\in V^\alpha M\qquad\text{for some }\alpha>0.
 \tag{P2}
\]
In particular the actual form is tracked through the direct image;
no parameter-dependent gauge or change of seed has occurred.

Apply the algebraic Fourier transform, with
\(t\mapsto-\partial_\tau\), \(\partial_t\mapsto\tau\).
The seed belongs to the unlocalized Fourier module \(\widehat M\).
The exponential period is its evaluation on the compact cycle.
Stokes makes this evaluation respect the relative de Rham relations,
and differentiation in \(\tau\) inserts \(-f\). Thus the period
evaluation respects all Weyl-algebra relations on the cyclic
submodule generated by the Fourier seed.

Suppose \(I(\tau)=c\neq0\), and divide the seed evaluation by \(c\).
The cyclic submodule
\[
 \widehat N=\mathbb C[\tau]\langle\partial_\tau\rangle\,\widehat m
 \subset\widehat M
\]
then has a quotient \(\widehat N\to\mathcal O_{\mathbb A^1_\tau}\)
taking \(\widehat m\) to \(1\). This is well defined because every
relation of the seed is respected by the period evaluation; the
image is the cyclic Weyl module generated by \(1\), namely
\(\mathbb C[\tau]\).

Inverse Fourier transform gives a submodule \(N\subset M\) and a
quotient
\[
 N\longrightarrow\delta_0,\qquad m\longmapsto\delta_0,
 \qquad
 \delta_0=D_{\mathbb A^1_t}/D_{\mathbb A^1_t}t.
 \tag{P3}
\]
The \(V\)-filtration on a submodule is the induced filtration, so
(P2) gives \(m\in V^\alpha N\). Every morphism of modules with
\(V\)-filtrations preserves the filtration. But
\[
 V^\alpha\delta_0=0\qquad(\alpha>0),
\]
because \(\delta_0\) is supported on \(t=0\).
Equation (P3) would therefore send \(m\) to zero, a contradiction.

There is no localization issue in this proof. The Fourier seed and
its cyclic quotient are constructed before inverting \(\tau\).
Even after localization, an analytic constant generates
\(\mathcal O_\tau\), not a \(\tau\)-torsion module. \(\square\)

### Exact source statements used

Q. Chen, B. Dirks, and M. Mustaţă,
*An introduction to \(V\)-filtrations*,
[arXiv:2402.07816v1](https://arxiv.org/pdf/2402.07816v1):

- Example 5.19, printed page 27, explicitly states (P1), using
  negativity of Bernstein--Sato roots and Corollary 5.15.
- Theorem 8.1, printed pages 41--44, gives proper-on-support
  push-forward for a smooth projection: the filtered relative
  de Rham direct image injects into the ordinary direct image,
  and its image is the canonical \(V\)-filtration.
  Its four-step proof is included there. Theorem 8.3, page 46,
  gives the general proper push-forward formulation.
- Proposition 3.19(ii), page 14, identifies the filtration on
  submodules and quotients. Corollary 3.20, page 15, gives
  functoriality and strictness.
- Example 3.18, page 14, gives \(V^\alpha=0\) for \(\alpha>0\)
  on every coherent module supported on \(t=0\), including
  \(\delta_0\).

The [retained primary version](sources/chen_dirks_mustata_V_filtrations_2402_07816v1.pdf) has
SHA-256 c00c78a3106e9f2b4120e59207c917e961787385902820d1b5fdb4e0745b7fe3.

The proof uses the decreasing \(V^\alpha\) convention throughout.

## Every pure cubic admits a constant-profile frame

This section applies the proper holomorphic-seed lemma above to explicit
local graph calculations. It includes every boundary root and every
coefficient degeneration.

Let \(p=H+E\) be a sphere cubic on \(X=\{q=1\}\). If the boundary
restriction \(H|_{q=0}\) is identically zero, then \(H=q\ell\) and
\(p\) has sphere degree at most two. The complete quadratic
classification applies.

Suppose \(H|_{q=0}\neq0\). Its boundary roots are finite. Assume,
toward a contradiction, that no null frame has a constant profile
\[
 p=y(a+bt+ct^2)+G(x,t),\qquad b=c=0.
 \tag{C1}
\]
At every boundary root there are the following exhaustive cases:
the point is not a frame incidence; its frame profile has
\(b^2-4ac\neq0\); or its profile is repeated with \(c\neq0\).
Indeed an incidence profile with zero discriminant and \(c=0\)
also has \(b=0\), contrary to the assumption.

Use the actual closure of the graph in
\[
 \bar X\times\mathbb A^1_s,\qquad
 \bar X=\{T^2+2XY=R^2\}\subset\mathbb P^3.
\]
Its equation is
\[
 H(X,Y,T)+R E(X,Y,T)-sR^3=0.
\]
This closure is proper over the target. Since \(H|_{q=0}\ne0\), the equation
has no boundary-supported component; locally it is the saturated closure
from \(R\ne0\).

First take the chart
\((X,Y,T,R)=(u^2-z^2,-1/2,u,z)\). Its equation is exactly
\[
 F=f(u)+zg(u)+z^2h(u,z)-s z^3=0,
\]
where \(f=H(u^2,-1/2,u)\), \(g=E(u^2,-1/2,u)\), and \(h\) is polynomial.
At a non-incidence boundary root either \(f'\ne0\) or \(g\ne0\), by the
binary frame criterion. Hence the graph is smooth there. The chart
\(X\ne0\) covers the other null directions.

At an incidence point choose its null frame and move that point to
\([X:Y:T:R]=[0:1:0:0]\). Use instead
\[
 (X,Y,T,R)=((z^2-u^2)/2,1,u,z).
\]
For the affine-\(y\) representative in (C1), its local graph equation is
\[
 F=P(u,z)-s z^3=0,\qquad
 P(u,z)=c u^2+buz+a z^2+R_{\geq3}(u,z).
\]
Indeed, the term \(x^it^j\) in \(G\) contributes
\(2^{-i}(z^2-u^2)^i u^j z^{3-i-j}\), of total order at least \(3+i\).
This derives the full quadratic tangent cone for every coefficient choice.
In either chart the original residue volume, up to a nonzero constant and
orientation, is
\[
 \Omega=\frac{du\wedge dz}{2z^2}.
\]
Such constant factors do not affect regularity.

### Non-incidence points

At a non-incidence point the preceding expansion proves that the graph is
smooth. Let
\[
 \eta_F=\operatorname{Res}\frac{du\wedge dz\wedge ds}{F}.
\]
Since \(F_s=-z^3\), the volume on the smooth graph is
\[
 \Omega=-\frac z2\,\eta_F.
 \tag{C2}
\]
The dualizing residue \(\eta_F\) is a regular two-form on this smooth
hypersurface. Consequently \(\Omega\) extends across these boundary
points, and in fact vanishes on their boundary.

### Ordinary profile nodes

If \(b^2-4ac\neq0\), blow up the fixed boundary section. In the chart
\(u=zv\), the strict graph equation is \(F=z^2F_1\), with
\[
 (F_1)_s=-z,\qquad
 \Omega=\frac{dv\wedge dz}{2z}=-\tfrac12\eta_{F_1}.
 \tag{C3}
\]
At every exceptional intersection with finite \(v\), the profile
polynomial has a simple root, so the strict transform is smooth.
Thus the volume is regular there.

If \(c=0,b\neq0\), the branch tangent to \(z=0\) is covered by
\(z=uv\). In this chart \(F=u^2F_2\), with
\[
 (F_2)_s=-uv^3,\qquad
 \Omega=\frac{du\wedge dv}{2uv^2}=-\frac v2\,\eta_{F_2}.
 \tag{C4}
\]
The strict transform has nonzero \(v\)-derivative \(b\) at the
exceptional point, so it is smooth and the last expression is
regular. Both charts include their full exceptional boundary.
Because \(H|_{q=0}\neq0\), the original graph equation has no
boundary-supported component; the strict transform is the closure
from \(z\neq0\).

### Repeated profiles with nonzero quadratic coefficient

After a nonzero scalar normalization, the local profile is
\(y(t-r)^2\). Put \(w=u-rz\), and blow up the boundary section in
the chart \(w=zv\). This chart contains its complete exceptional
intersection, because the repeated tangent direction is \(w=0\).
The strict transform equation has the exact form
\[
 F_1=v^2+z\bigl(k(v)-s\bigr)+z^2h(v,z)=0.
 \tag{C5}
\]
It is smooth at the boundary \(v=z=0\) except possibly at
\(s=K:=k(0)\). At that point its quadratic part in
\((v,z,s-K)\) is nondegenerate: its Hessian determinant is \(-2\),
independently of the other coefficients. It is therefore an
ordinary \(A_1\) surface double point.

The volume is again the dualizing generator
\[
 \Omega=-\tfrac12\eta_{F_1},
\]
by (C3). Blow up the ordinary double point. Its strict transform is
smooth and its exceptional divisor is a smooth conic. The
resolution is crepant: in a smooth threefold the blow-up of a point
adds \(2E\) to the canonical divisor, while the double-point
hypersurface has total transform equal to its strict transform
plus \(2E\). Adjunction cancels these two contributions.
Consequently the pullback of \(\eta_{F_1}\), and hence of
\(\Omega\), is holomorphic on the resolution.

This argument permits non-isolated critical loci of the map to the
target. Only smoothness of the resolved total surface and
holomorphicity of its volume are required by the proper seed
lemma.

### Global contradiction

The modifications just described are proper, are isomorphisms over
the original affine quadric \(X\), and concern finitely many fixed
boundary sections and finitely many \(A_1\) points. Performing them
gives a smooth irreducible surface \(Y\), proper over the target
line, which contains \(X\) as a dense open subset. Equations
(C2)--(C5) prove that the original volume extends to a global
holomorphic algebraic two-form on \(Y\). Its map \(f:Y\to\mathbb A^1\)
is nonconstant because it extends the nonconstant phase \(p\).

The original real sphere remains a compact cycle in \(Y\), with
nonzero volume integral. If \(p\) were pure, its exponential period
on this fixed holomorphic form would be the nonzero constant one.
The proper holomorphic-seed lemma excludes this, without any
isolated-critical-point assumption.

Therefore the assumption preceding (C1) is false. Every pure sphere
cubic of actual degree three admits a null frame with constant
profile \(p=a y+G(x,t)\), \(\deg G\leq3\).
The complete eleven-parameter constant-profile classification can
then be applied. This proves a full-family reduction rather than a
generic statement in the repeated-root parameter \(r\); no
specialization exception in \(r\) is omitted.

## Complete classification and the original Image conclusion

First suppose the sphere phase has actual degree at most two. The full
nine-parameter quadratic proposition below gives one-sided support in an
orthogonal frame, including the zero polynomial.

Otherwise the preceding theorem provides a frame in which
\[
 p=ay+ht^3+(v_0+v_1x)t^2+(c+dx+ex^2)t+\sum_{j=0}^3r_jx^j.
\]
The eleven-parameter certificate below treats all coefficients, including
\(a=0\), and proves that purity forces \(h=0\) and precisely one of:

1. \(a=0\), \(v_0=c=r_0=0\). All surviving monomials have positive circular
   weight.
2. \(a\ne0\), \(v_0=v_1=e=0\), and
   \[
   p=ay+(c+dx)t-\frac{x(c+dx)^2}{2a}.
   \]

In the second case, the orthogonal null rotation (A4) with \(k=c/a\) gives
\[
 p=ay+dxt-\frac{d^2x^3}{2a}.
\]
If \(d=0\), this is one-sided. If \(d\ne0\), an overall scaling and the
orthogonal torus \(x\mapsto\lambda x,y\mapsto y/\lambda\) normalize it to
\(p_*\). For example choose nonzero \(\alpha,\lambda\) with
\(\alpha a/\lambda=\alpha d\lambda=1\). These exist over \(\mathbb C\).
Conversely, one-sided phases are pure by phase extraction, and \(p_*\) is
pure by the volume-preserving polynomial shear (8). This proves both
directions of the full classification.

For a sphere multiplier \(g\), one-sided phases satisfy mixed vanishing for
\(m>\deg g\). The shear and its inverse have degree at most three, so the
exceptional family satisfies mixed vanishing for \(m>3\deg g\). Orthogonal
linear substitutions preserve degree and the sphere functional, hence these
bounds hold before normalization. No lift of an arbitrary orthogonal
reflection to a linear substitution on the unbalanced original variables
is needed: the balanced-multiplier bridge (6) uses only the invariant sphere
functional.

A multiplier of original bidegree \((a,b)\), with \(b\geq a\), is tested by
the finitely many balanced multipliers \(\xi^\gamma Q\) of bidegree
\((b,b)\). Their sphere degrees are at most \(b\). Equations (4) and (6)
therefore prove (1), including arbitrary unbalanced multipliers. If
\(a>b\), contraction is automatically zero. Splitting a general multiplier
into bidegrees yields the bound \(m>3\deg_zQ\).

The proof over \(\mathbb C\) descends to every characteristic-zero field:
for each fixed \(F,Q\), embed the finitely generated field of their
coefficients into \(\mathbb C\). Every zero contraction identity for the
embedded coefficients descends by injectivity. The uniform bound depends
only on \(\deg_zQ\), so no choice of embedding affects it.

The same spherical argument covers each fixed balanced bidegree below three. At degree
four the rational witness in the linked counterexample has all pure
contractions zero and a fixed multiplier with every mixed contraction
nonzero. Thus positivity through degree three and failure at degree four
are both attained. The all-order mixed formula (10) proves sharpness of
the coefficient three in (1).

## Safe mixed cubics outside the orthogonal nullcone

Let delta=x partial_t-t partial_y. It is an infinitesimal orthogonal action,
delta(q)=0, and delta(x)=0. Sphere invariance gives L(delta h)=0; therefore
L(x delta h)=L(delta(xh))=0. The derivation x delta is locally nilpotent,
so its exponential preserves L by a finite polynomial expansion on every
input. In particular the automorphism phi=exp(-x delta) is

\[
 \phi(x)=x,\qquad\phi(t)=t-x^2,\qquad
 \phi(y)=y+xt-\tfrac12x^3=:p_*.
 \tag{8}
\]

Its inverse replaces x delta by -x delta in the exponential. Both maps
preserve q. Thus L(p_*^m)=L(y^m)=0 for every m>=1, and

\[
 L(gp_*^m)=L(\phi^{-1}(g)y^m)=0\qquad(m>3\deg g).
 \tag{9}
\]

The reason is that deg(phi^(-1)g)<=3 deg g, and its positive phase weight
cannot exceed its degree. Equation (6) gives the Image cutoff
m>3 deg_z Q for the balanced cubic

\[
 F_*=R^2Y+RXT-\tfrac12 X^3.
\]

This coefficient three is necessary even on this safe family. For D>=1,
only the highest positive-weight monomial in phi^(-1)(y^D) can balance
y^(3D), so

\[
 L(y^D p_*^{3D})=
 (-\tfrac12)^D L((xy)^{3D})
 =\frac{(-1)^D 2^{2D}((3D)!)^2}{(6D+1)!}\ne0.
 \tag{10}
\]

Consequently E_2(Y^D F_*^(3D)) is (10) times (10D+1)!, and is nonzero.
There is also an explicit invariant proving that F_* is outside the
orthogonal nullcone (the set whose orbit closure contains zero). Its odd
homogeneous component is H_*=qy-x^3/2, and

\[
 L(H_*^4)=L((y-x^3/2)^4)=-2L(x^3y^3)=-4/35\ne0.
 \tag{11}
\]

The map (H,E) to L(H^4) is an invariant polynomial of positive degree, so it
vanishes on that nullcone. Equation (11) also shows that purity of a mixed
sum does not force purity of its two parity summands. This does not make
F_* an Image counterexample; (9) proves its safety.

## The odd parity stratum as a consequence of the full classification

Every odd sphere polynomial of degree at most three has a unique homogeneous
cubic lift \(H\): replace its linear part by \(q\) times that part. The
odd stratum therefore has dimension ten.

Orthogonal substitutions preserve antipodal parity. The exceptional shear
\(p_*\) has the nonzero even part \(xt\), so no member of its orthogonal
and nonzero scalar orbit is odd. The full classification therefore makes
every pure odd phase one-sided. Its homogeneous cubic lift consequently has
only the four positive-weight cubic monomials:
\[
 H=x(ax^2+bxt+ct^2+dxy).
\]
Its original Image cutoff improves to \(m>\deg_zQ\) by the weight argument
or the balanced-multiplier bridge. This conclusion uses the full spherical
proof above and has no Gaussian classification dependency.

There remains a useful consistency identity: for a homogeneous polynomial
\(K\) of even degree \(n\), the Gaussian Wick functional satisfies
\(G(K)=(n+1)!!L(K)\). Thus \(G(H^{2r})=(6r+1)!!L(H^{2r})\), while odd
powers vanish by parity. The varying factors on mixed homogeneous degrees
are precisely why this identity alone does not settle the sixteen-parameter
family.

## The even parity stratum: a proof using three moments

Every even member of \(A_{\leq3}\), where
\(A=\mathbb C[x,y,t]/(t^2+2xy-1)\), is the restriction of a unique
homogeneous quadratic \(E\). A constant \(c\) is homogenized as \(cq\).
The injectivity follows by comparing degrees in a proposed identity
\(E=(q-1)c\). Hence this stratum has dimension six.

Let \(B\) be the bilinear form with \(q(v)=B(v,v)\), and write
\(E(v)=B(v,Av)\) for the unique \(B\)-selfadjoint operator \(A\). The
first three sphere moments are
\[
 \begin{split}
 L(E)&=\tfrac13\operatorname{tr}A,\\
 L(E^2)&=\tfrac1{15}
 \left((\operatorname{tr}A)^2+2\operatorname{tr}A^2\right),\\
 L(E^3)&=\tfrac1{105}
 \left((\operatorname{tr}A)^3+
 6\operatorname{tr}A\,\operatorname{tr}A^2+
 8\operatorname{tr}A^3\right).
 \end{split}
 \tag{E1}
\]
These identities hold for every selfadjoint \(A\), including defective
ones. One derivation first proves the Gaussian identities by Wick
pairings, then divides by \(3!!,5!!,7!!\), respectively. Alternatively
they are the first three coefficients of the formal identity
\[
 G(e^{sE})=\det(I-2sA)^{-1/2};
\]
this identity follows by diagonalizing a dense open set of symmetric
matrices, multiplying the one-variable formal Gaussian identity, and
then extending each polynomial coefficient to all symmetric matrices.
Thus no spectral diagonalization assumption is imposed on \(A\).

If the first three sphere moments vanish, equation (E1) successively
gives
\[
 \operatorname{tr}A=\operatorname{tr}A^2=
 \operatorname{tr}A^3=0.
\]
Newton's identities in dimension three now give the characteristic
polynomial \(\lambda^3\). Cayley--Hamilton proves \(A^3=0\).

Assume \(A\neq0\). If \(A^2\neq0\), choose \(v\) with \(u=A^2v\neq0\).
Then
\[
 Au=0,\qquad B(u,u)=B(v,A^4v)=0.
\]
If \(A^2=0\), choose \(v\) with \(u=Av\neq0\). Then again
\[
 Au=0,\qquad B(u,u)=B(v,A^2v)=0.
\]
In both cases \(\ker A\) contains a nonzero null vector \(u\). Complete
\(u\) to a Witt basis and identify it with the \(y\)-axis. Because
\(A\) is selfadjoint,
\[
 D_uE(v)=2B(u,Av)=2B(Au,v)=0.
\]
Thus \(E\) is independent of \(y\), and has the form
\[
 E=ax^2+bxt+ct^2.
\]
In these Witt coordinates \(\operatorname{tr}A=c\). Therefore \(c=0\).
The resulting form
\[
 E=ax^2+bxt
 \tag{E2}
\]
has only strict positive circular weights two and one. Its all-order
pure sphere moments vanish by phase extraction, and every fixed
sphere multiplier has eventual mixed vanishing. The case \(A=0\)
gives \(E=0\) and the same conclusion.

Consequently, for the full even parity stratum, the following
conditions are equivalent: the first three pure sphere moments
vanish; all positive pure sphere moments vanish; and the phase is
orthogonally equivalent to (E2). No Gaussian ternary-cubic
classification is used in this argument.

In original balanced bidegree \((3,3)\), the even phase lifts as
\(F=R(aX^2+bXT)\). The factor \(R\) has torus weight zero, so the
positive-weight argument together with (6) proves
\[
 \mathcal E_2(QF^m)=0\qquad(m>\deg_z Q)
\]
for every polynomial multiplier \(Q\).

## The full quadratic sphere family: a sharp six-moment criterion

Let \(q(v)=B(v,v)\) be nondegenerate in dimension three, and write a
sphere polynomial of degree at most two uniquely as
\[
 p(v)=E(v)+\ell(v),\qquad
 E(v)=B(v,Av),\quad \ell(v)=B(b,v),
\]
where \(A\) is \(B\)-selfadjoint. A constant in \(p\) has been absorbed
into \(E\) by multiplication by \(q\). This is the full
nine-dimensional family, with six quadratic and three linear
coefficients.

**Proposition.** The following are equivalent:

1. \(L(p^m)=0\) for \(1\leq m\leq6\).
2. \(L(p^m)=0\) for every \(m\geq1\).
3. After an orthogonal change of coordinates, with \(q=t^2+2xy\),
   \[
      p=ax^2+bxt+cx.
      \tag{Q1}
   \]

Consequently every pure member has eventual mixed vanishing, and its
balanced-degree-three lift satisfies
\[
 \mathcal E_2(QF^m)=0\qquad(m>\deg_zQ)
\]
for every original polynomial multiplier \(Q\). The prefix six in
condition 1 is sharp within this quadratic sphere family.

### Exact moment equations with the sphere weights retained

Set \(T_j=\operatorname{tr}(A^j)\) and \(B_j=B(b,A^jb)\).
The first moment is \(T_1/3\), so assume \(T_1=0\).
Cayley--Hamilton gives
\[
 A^3=\tfrac12T_2 A+\tfrac13T_3 I.
 \tag{Q2}
\]
Thus
\[
 T_n=\tfrac12T_2T_{n-2}+\tfrac13T_3T_{n-3}\quad(n\geq4),
 \qquad
 B_n=\tfrac12T_2B_{n-2}+\tfrac13T_3B_{n-3}\quad(n\geq3).
\]

Here is a reproducible derivation of the needed invariant moments.
For the Gaussian functional \(G\) with covariance \(B\), the formal
quadratic identity is
\[
 G(e^{sE+z\ell})
 =D(s)e^{z^2 C(s)},\quad
 D(s)=\exp\!\left(\sum_{n\geq1}\frac{2^{n-1}}nT_n s^n\right),
 \quad C(s)=\tfrac12\sum_{n\geq0}2^nB_n s^n.
 \tag{Q3}
\]
All inverses in this formal identity are expanded at \(s=0\), where
\(I-2sA\) is a unit. It therefore imposes no invertibility or
diagonalizability assumption on \(A\).

Crucially, homogeneous summands are converted to sphere moments
separately. If \(j+2k=m\), the term \(E^j\ell^{2k}\) has degree
\(2j+2k\). Therefore
\[
 M_m:=L(p^m)
 =m!\sum_{k=0}^{\lfloor m/2\rfloor}
   \frac{[s^{m-2k}]\,D(s)C(s)^k}
        {k!\,(2m-2k+1)!!}.
 \tag{Q4}
\]
This uses the Gaussian series only to compute explicit identities.
It does not identify the mixed-degree sphere functional with a
Gaussian functional.

Expanding (Q4) using (Q2) gives the following triangular system.
After the preceding equations have been set to zero, respectively,
\[
\begin{aligned}
 M_2&=(5B_0+2T_2)/15,\\
 M_3&=2(21B_1+4T_3)/105,\\
 M_4&=16(75B_2+4T_2^2)/2625.
\end{aligned}
\]
Consequently
\[
 B_0=-\frac25T_2,\qquad
 B_1=-\frac4{21}T_3,\qquad
 B_2=-\frac4{75}T_2^2.
 \tag{Q5}
\]
There is no division by any parameter. Under (Q5), the next moments
are exactly
\[
 M_5=-\frac{2048}{24255}T_2T_3,\qquad
 M_6=-\frac{8192}{165540375}(147T_2^3+500T_3^2).
 \tag{Q6}
\]
If both vanish, then \(T_2T_3=0\). In either branch the other
equation forces \(T_2=T_3=0\). Equation (Q5) also gives
\(B_0=B_1=B_2=0\), and (Q2) gives \(A^3=0\).

### The common null direction, including defective cases

If \(A^2\neq0\), the nilpotent operator \(A^2\) has rank one.
It is selfadjoint, so
\[
 A^2=\kappa\,u\otimes B(u,\cdot),\qquad \kappa\neq0,
\]
where \(u\neq0\) spans its image. The equations \(A^3=0\) and
selfadjointness show that \(Au=0\) and \(B(u,u)=0\).
Since
\[
 B_2=B(b,A^2b)=\kappa B(b,u)^2=0,
\]
we also have \(\ell(u)=0\).

If \(A^2=0\) but \(A\neq0\), its image is totally isotropic, since
\(B(Av,Aw)=B(v,A^2w)=0\). A totally isotropic subspace in dimension
three has dimension at most one. Hence \(A\) has rank one and
\[
 A=\kappa\,u\otimes B(u,\cdot)
\]
with \(u\neq0\) null and \(Au=0\). The equation \(B_1=0\) again gives
\(\ell(u)=0\).

In both cases complete \(u\) to a Witt basis, identifying it with
the \(y\)-axis. Selfadjointness and \(Au=0\) imply \(D_uE=0\);
thus \(E=ax^2+bxt+dt^2\). Its trace is \(d\), which is zero.
The equation \(\ell(u)=0\) says that
\(\ell=\alpha x+\gamma t\). Its squared norm is \(B_0=\gamma^2\),
so \(\gamma=0\). This proves (Q1).

If \(A=0\), the condition \(B_0=0\) makes \(b\) null. For
\(b\neq0\), choose \(u=b\) and the same Witt-coordinate argument
puts \(\ell\) in the \(x\)-line. If \(b=0\), the phase is zero.
These cases exhaust the full family.

Every nonzero monomial in (Q1) has positive circular weight.
Hence all positive pure moments vanish, and a fixed multiplier
cannot balance sufficiently large powers. This proves both the
converse and eventual mixed vanishing. In balanced degree three
the lift is \(F=R E(X,Y,T)+R^2\ell(X,Y,T)\); its nonzero original
torus weights are at least two. The arbitrary-multiplier weight
argument therefore gives the stated cutoff \(m>\deg_zQ\).

### Sharpness and exact checks

In Cartesian coordinates with \(q=u^2+v^2+w^2\), take
\[
 p_{\rm sharp}
 =u^2-v^2+\frac{2i}{15}
       \left(\sqrt6\,(u+v)+\sqrt{33}\,w\right).
 \tag{Q7}
\]
Here
\[
 A=\operatorname{diag}(1,-1,0),\quad
 T_2=2,\ T_3=0,\quad
 B_0=-4/5,\ B_1=0,\ B_2=-16/75.
\]
Equations (Q4)--(Q6) give
\[
 L(p_{\rm sharp}^m)=0\quad(1\leq m\leq5),\qquad
 L(p_{\rm sharp}^6)=-\frac{65536}{1126125}\neq0.
\]
Thus five moments do not suffice within the same quadratic family.

The retained structure checker implements (Q2)--(Q6) in exact
symbolic arithmetic. Independent Cartesian monomial evaluation
checks the generic first three moments, three explicit matrix/vector
controls through the sixth moment (including rank-two nilpotent
and zero matrices), and (Q7) directly through order six. Replacing
the factor four in the last equation of (Q5) by one is rejected
by the fourth-moment identity. These checks certify the displayed
algebraic moment equations. The preceding nilpotent linear algebra
and weight argument supply the classification and its all-order
consequence.

## A complete seven-parameter phase slice

Let \(L\) denote the normalized spherical moment functional, and consider
\[
 p=y(c_0+c_1x+c_2x^2)+d_0+d_1x+d_2x^2+d_3x^3.
 \tag{S1}
\]
This slice has dimension seven. Its complete classification is a proper
subfamily result and supplies no stratification of all sixteen parameters.

**Proposition.** For (S1), the following conditions are equivalent:

1. \(L(p^m)=0\) for \(1\le m\le6\).
2. \(p\) has strictly one-sided phase support in the displayed coordinates,
   including the zero polynomial.
3. \(L(p^m)=0\) for every \(m\ge1\).

More precisely, in the \(c_0\ne0\) branch, \(p=c_0y\). In the
\(c_0=0\) branch, \(c_1=d_0=0\), leaving only positive weights. The fixed
bound six is sharp for this slice.

**Proof.** The monomial rule
\[
 L(x^iy^j)=\delta_{i,j}\frac{i!}{(2i+1)!!}
 \tag{S2}
\]
follows by phase averaging and \(xy=(1-t^2)/2\):
\[
 L((xy)^i)=2^{-i}\int_0^1(1-t^2)^i\,dt
          =\frac{i!}{(2i+1)!!}.
\]
Thus for \(p=yC(x)+D(x)\), the exact identity is
\[
 L(p^m)=\sum_{k=0}^m\binom mk\frac{k!}{(2k+1)!!}
             [x^k]C(x)^kD(x)^{m-k}.
 \tag{S3}
\]

First suppose \(c_0\ne0\). Divide \(p\) by \(c_0\), write
\(C=1+ax+bx^2\), and use \(d_i\) for the rescaled coefficients of
\(D\). The first moment gives \(d_0=-a/3\). Successively setting
moments two, three, and four to zero gives
\[
\begin{aligned}
 d_1&=-\frac{a^2}{30}-\frac{2b}{5},\\
 d_2&=\frac{a(31a^2-108b)}{1890},\\
 d_3&=-\frac{289a^4}{28350}
       +\frac{47a^2b}{1050}-\frac{2b^2}{75}.
\end{aligned}
 \tag{S4}
\]
The coefficient of each new variable \(d_1,d_2,d_3\) in its equation is,
respectively, \(2/3,2/5,8/35\), so no division by a parameter is hidden
in (S4). After these substitutions, the next two moments are
\[
\begin{aligned}
 L(p^5)&=-\frac{4a}{9823275}
 (2219a^4-11874a^2b+13176b^2),\\
 L(p^6)&=-\frac8{13408770375}
 (482933a^6-1979487a^4b+218781a^2b^2+1571724b^3).
\end{aligned}
 \tag{S5}
\]
If \(a\ne0\), put \(z=b/a^2\). Vanishing in (S5) forces
\[
 P(z)=2219-11874z+13176z^2=0,
 \quad Q(z)=482933-1979487z+218781z^2+1571724z^3=0.
\]
These two equations are incompatible: the exact certificate is
\[
 A(z)P(z)+B(z)Q(z)=76274195239375,
 \tag{S6}
\]
where
\[
\begin{aligned}
 A(z)&=6676468225812z^2+2727590293062z-7662133538547,\\
 B(z)&=35364219296-55969842888z.
\end{aligned}
\]
It follows that \(a=0\). Then (S5) becomes
\(L(p^6)=-32b^3/34125\), so \(b=0\). Equations (S4) and the first
moment give \(d_0=d_1=d_2=d_3=0\). Undoing the scaling gives \(p=c_0y\).

If \(c_0=0\), the only phase-zero terms are \(c_1xy+d_0\); all other
terms have positive phase weight. They cannot contribute to any pure
moment. The first moment gives \(d_0=-c_1/3\); after this substitution,
the second moment is \(c_1^2/45\). Therefore \(c_1=d_0=0\), leaving
strictly positive phase weights. In either branch strict one-sidedness
proves all pure moments vanish, and it proves eventual mixed vanishing
for each fixed polynomial multiplier by the same weight argument.

Finally, take
\[
 p_*=y+x^2y-\frac25x-\frac2{75}x^3.
 \tag{S7}
\]
This is (S4) with \(a=0,b=1\), so
\[
 L(p_*^m)=0\quad(1\le m\le5),
 \qquad L(p_*^6)=-\frac{32}{34125}\ne0.
\]
Thus five moments do not suffice, even inside this slice. The proof above
shows that six do. This is an exact cutoff for (S1), not a cutoff for the
full cubic family. \(\square\)

The symbolic certificate is replayed by
`research/foundations/scripts/verify_two_pair_image_cubic_phase_slice.py`
with the existing `sympy==1.14.0` environment. It checks the literal
triangular identities, both exceptional branches, the nonzero Bezout
constant, rejection of an altered Bezout coefficient, and (S7). It uses
no floating-point arithmetic. The retained JSON records these equations
and their exact scope; the proof's support argument supplies the
all-order implication after the finite equations have classified (S1).



## Opposite quadratic height profiles

Consider the six-parameter family
\[
 p=yC(t)+xA(t),\qquad \deg C,\deg A\le2.
 \tag{OP1}
\]
It includes nonconstant negative-weight coefficients \(yt\) and \(yt^2\),
which are absent from the constant-\(y\) chart. It is still a proper
subfamily of the thirteen-parameter affine sphere family.

**Proposition.** The first ten spherical moments of (OP1) vanish if and
only if \(C A=0\). Consequently every pure point in this family has
strictly one-sided phase support and eventual mixed vanishing.

**Proof.** Odd powers have no phase-zero term. In the even powers the only
phase-zero selection uses equal numbers of the two summands. Therefore
\[
 L(p^{2r})=\binom{2r}{r}2^{-r}I(h^r),\qquad
 h(t)=(1-t^2)K(t),\quad K=CA,\quad
 I(f)=\frac12\int_{-1}^1f(t)\,dt.
 \tag{OP2}
\]
Write \(K=\sum_{j=0}^4k_jt^j\). For \(1\le r\le5\), let
\(M_r=I(h^r)\). These are explicit homogeneous polynomials over the
rationals. For completeness, their deterministic coefficient formula is
\[
 M_r=\sum_{n_0+\cdots+n_4=r\atop d=\sum jn_j\text{ even}}
 \frac{2^r(r!)^2(d-1)!!}{(2r+d+1)!!\prod n_j!}
 \prod_{j=0}^4k_j^{n_j}.
 \tag{OP3}
\]
In particular,
\[
 M_1=\frac23k_0+\frac2{15}k_2+\frac2{35}k_4,
 \qquad k_0=-\frac15k_2-\frac3{35}k_4.
 \tag{P4}
\]
The following finite exact coefficient classification closes every branch.
If \(K\ne0\), let \(j\) be its highest nonzero coefficient. Scale \(K\)
so that \(k_j=1\), and put all higher coefficients equal to zero. After
(P4), compute the ideal generated by the next \(j\) equations
\(M_2,\ldots,M_{j+1}\) in \(k_1,\ldots,k_{j-1}\). For each
\(j=4,3,2\), its reduced Groebner basis over \(\mathbb Q\), in graded
reverse lexicographic order, is the singleton \(\{1\}\). Thus none of
these equations has a simultaneous solution over \(\mathbb C\).
For \(j=1\), (P4) gives \(k_0=0\) and
\(M_2=8/105\ne0\). For \(j=0\), \(M_1=2k_0/3\ne0\).
This exhausts all nonzero \(K\), so \(M_1=\cdots=M_5=0\) forces
\(K=0\).

The exact checker computes (OP3), performs all four highest-coefficient
normalizations, and verifies these unit ideals in the pinned SymPy 1.14.0
rational polynomial domain. The literal generator polynomials and branch
outputs are retained in the phase-slice verification JSON. As an additional
internal check, rational unit-ideal certificates generated with Singular
were multiplied out independently in SymPy; no numerical approximation
is involved. The unit-ideal computations, rather than tests of sample
coefficients, establish the finite completeness argument.

Equation (OP2) now proves \(CA=0\). Since \(\mathbb C[t]\) is an integral
domain, either \(C=0\) or \(A=0\). This gives strictly positive or
strictly negative phase support, respectively, and proves all-order pure
vanishing and eventual mixed vanishing. Conversely these one-sided
families satisfy every pure moment equation. \(\square\)

The fifth univariate equation is substantive: in the normalized
\(k_4=1\) chart the first four equations generate a proper ideal. The
checker verifies this negative control; merely dropping the last equation
would not prove the proposition. No cutoff for arbitrary cubic phases is
inferred from this profile result.

## A complete ten-parameter affine slice

Consider

\[
 p=ay+(v_0+v_1x)t^2+(c+dx+ex^2)t+\sum_{j=0}^3r_jx^j.
 \tag{A1}
\]

The ten coefficients vary independently, including \(a=0\). This slice
omits \(yt,yt^2,t^3\) from the thirteen-parameter affine family.

**Proposition.** For (A1), vanishing of the first eight pure sphere moments
is equivalent to either of the following disjoint conditions:

1. \(a=0\) and \(v_0=c=r_0=0\).
2. \(a\ne0\), \(v_0=v_1=e=0\), and
   \[
   \sum_{j=0}^3r_jx^j=-\frac{x(c+dx)^2}{2a}.
   \tag{A2}
   \]

In both cases all pure moments vanish, and every fixed multiplier has eventual
mixed vanishing. The sphere cutoff is \(m>3\deg g\); the Image cutoff is
\(m>3\deg_z Q\). No optimality of the finite prefix eight is asserted.

**Proof.** First take \(a=0\). Phase averaging removes every term of positive
\(x\)-weight, leaving the univariate polynomial \(v_0t^2+ct+r_0\).
The first moment gives \(r_0=-v_0/3\). The second and then the third give

\[
 L(p^2)=\frac{4v_0^2}{45}+\frac{c^2}{3},\qquad
 L(p^3)\big|_{c^2=-4v_0^2/15}=-\frac{256v_0^3}{4725}.
 \tag{A3}
\]

Thus \(v_0=c=r_0=0\), leaving strict positive weights. This also proves the
converse and mixed vanishing in this branch.

If \(a\ne0\), divide by \(a\) and rename coefficients, so \(a=1\). The
orthogonal null rotation

\[
 x=x',\qquad t=t'+kx',\qquad y=y'-kt'-k^2x'/2
 \tag{A4}
\]

preserves the slice, changes the constant \(t'\)-coefficient from \(c\) to
\(c-k\), and preserves \(v_0,v_1\). Take \(k=c\). Hence it suffices to
classify the branch \(a=1,c=0\).

Formula (3) computes each pure moment exactly. The first four moments solve
successively for \(r_0,r_1,r_2,r_3\), with respective nonzero constant
coefficients \(1,2/3,2/5,8/35\). The solutions are

\[
\begin{aligned}
 r_0&=-v_0/3,\\
 r_1&=-(2v_0^2+3v_1)/15,\\
 r_2&=-4v_0(31v_0^2+54v_1)/945,\\
 r_3&=-(14175d^2+4624v_0^4+10152v_0^2v_1+3024v_1^2)/28350.
\end{aligned}
 \tag{A5}
\]

After substitution, write \(L(p^j)=\kappa_jN_j\), \(5\leq j\leq8\), where

\[
 (\kappa_5,\kappa_6,\kappa_7,\kappa_8)=
 \left(\frac8{9823275},-\frac8{13408770375},
 \frac{32}{273648375},-\frac{256}{1465387048125}\right).
 \tag{A6}
\]

The checker and its retained certificate give the four literal polynomials
\(N_j\). If \(v_0=0\), their specializations are

\[
 N_8=9144576v_1^4,\qquad N_6\big|_{v_1=0}=-58046625e^2.
 \tag{A7}
\]

Thus \(v_1=e=0\), and (A5) leaves only \(r_3=-d^2/2\), as required.

If \(v_0\ne0\), put

\[
 b=v_1/v_0^2,\qquad D=d^2/v_0^4,\qquad E=de/v_0^5,\qquad F=e^2/v_0^6.
\]

The four equations \(N_j/v_0^j=0\) are linear in \(D,E,F\). Their augmented
matrix, with columns \(D,E,F,1\), is

\[
 \begin{pmatrix}
 155925&155925&0&B_5\\
 77395500-116093250b&77395500&-58046625&B_6\\
 226800-297675b&226800+127575b&-212625&B_7\\
 4819500-216877500b&4819500-65063250b&-75907125&B_8
 \end{pmatrix},
 \tag{A8}
\]

where

\[
\begin{aligned}
 B_5&=35504+94992b+52704b^2,\\
 B_6&=30907712+63343584b+3500496b^2-12573792b^3,\\
 B_7&=53504+17328b-211968b^2-108864b^3,\\
 B_8&=7986496-27544704b-90351504b^2-16632864b^3+9144576b^4.
\end{aligned}
\]

The coefficient columns satisfy \(C_D-C_E-2bC_F=0\), so the coefficient
matrix has rank at most two. Consistency forces all three-by-three augmented
minors to vanish. Columns \(D,F,1\) and row triples \((5,6,7)\), \((5,6,8)\)
give, after removing nonzero rational factors,

\[
\begin{aligned}
 f(b)&=2485161b^3+6647562b^2+4908168b+1097984=0,\\
 g(b)&=10216206b^4+24356619b^3-25372278b^2-41158800b-11711872=0.
\end{aligned}
 \tag{A9}
\]

But

\[
 \operatorname{Res}_b(f,g)
 =-3179537325694882741702679455938048000000000\ne0.
 \tag{A10}
\]

The receipt also retains rational polynomials \(U,V\) with the checked identity
\(Uf+Vg=1\). This excludes \(v_0\ne0\), including all possible degeneracies
of the linear equations. Undoing (A4) and the nonzero scalar normalization
recovers (A2).

For sufficiency in the nonzero branch, let \(h(x)=(c+dx)/a\) and use

\[
 \phi(x)=x,\quad \phi(t)=t-xh(x),\quad
 \phi(y)=y+h(x)t-xh(x)^2/2.
\]

This is \(\exp(-h(x)\delta)\), with \(\delta=x\partial_t-t\partial_y\).
Since \(\delta h=0\), infinitesimal sphere invariance gives
\(L(h\delta g)=L(\delta(hg))=0\). Local nilpotence proves that the
polynomial automorphism preserves \(L\). Condition (A2) says \(p=a\phi(y)\).
Its inverse increases degrees by at most three, proving the stated all-order
and mixed conclusions exactly as in (9). This completes the proposition.

The [affine-slice checker](../scripts/verify_two_pair_image_cubic_affine_slice.py)
derives the moments directly from the sphere monomial rule and verifies the
normalizations, all exceptional branches, (A5)–(A10), and the Bezout identity.
This result does not prove that every pure cubic enters (A1).


## Excluding the additional height-cubic term: a complete eleven-dimensional slice

The preceding ten-dimensional classification extends to
\[
 p=ay+h t^3+(v_0+v_1x)t^2+(c+dx+ex^2)t+\sum_{j=0}^3r_jx^j.
\]
All eleven coefficients vary, including the branches \(a=0\) and
\(h=0\). This slice still excludes \(yt\) and \(yt^2\).

**Theorem.** The equations \(L(p^m)=0\), \(1\le m\le10\), hold if and
only if \(h=0\) and one of the following conditions holds:

1. \(a=0\) and \(v_0=c=r_0=0\).
2. \(a\ne0\), \(v_0=v_1=e=0\), and
   \[
   \sum_{j=0}^3r_jx^j=-\frac{x(c+dx)^2}{2a}.
   \]

Consequently every pure phase in this eleven-dimensional slice is safe
for every fixed sphere multiplier and hence, by the balanced-multiplier
reduction, for the original Image functional. The ten-dimensional slice
retains its stronger eight-moment certificate.

**Proof of the branch \(a=0\).** Phase averaging leaves only
\(h t^3+v_0t^2+ct+r_0\). The first moment gives \(r_0=-v_0/3\).
Write \(v=v_0\). For
\(u=h t^3+v(t^2-1/3)+ct\), the next three moments satisfy
\[
 M_2=\frac{105c^2+126ch+45h^2+28v^2}{315},
\]
\[
 M_3\bmod M_2=\frac{8v(56c^2+105ch+45h^2)}{2205},
\]
\[
 M_4\bmod M_2=-\frac{16R(c,h)}{735735},
\]
where the remainders are polynomial division in \(v\), and
\[
 R(c,h)=28028c^4+84084c^3h+93275c^2h^2+46410ch^3+8955h^4.
\]
If \(h=0\), the quadratic argument already gives \(v=c=0\).
For \(h\ne0\), scale \(h=1\) and put \(z=c/h\). Define
\[
 P(z)=56z^2+105z+45,\qquad B(z)=105z^2+126z+45,
 \qquad R(z)=R(z,1).
\]
If \(v\ne0\), the equations force \(P=R=0\), contradicting
\[
 (-1197196z^3-2954952z^2-2305927z-635895)P
 +(2392z+3213)R=157140.
\]
If \(v=0\), they force \(B=R=0\), contradicting
\[
 (-468468z^3-1845844z^2-2151045z-739960)B
 +(1755z+3756)R=336780.
\]
Thus \(h=v=c=r_0=0\); all remaining monomials have strict positive phase
weight. This proves both directions and safety in the \(a=0\) branch.

**Proof of the branch \(a\ne0\).** Divide by \(a\) and rename
coefficients, so that \(a=1\). The null rotation
\[
 x=x',\qquad t=t'+kx',\qquad y=y'-kt'-k^2x'/2
\]
keeps \(h\) unchanged and changes the constant \(t'\)-coefficient to
\(c-k\). It preserves the eleven-dimensional slice; in particular its
new \(x't'^2\)-coefficient is \(v_1+3hk\). Taking \(k=c\) gives
\(c=0\). If \(h=0\), the ten-dimensional theorem applies.

Suppose \(h\ne0\). The torus change \(x\mapsto\lambda x\),
\(y\mapsto y/\lambda\), followed by multiplication of \(p\) by
\(\lambda\), keeps its \(y\)-coefficient equal to one and changes
\(h\) to \(\lambda h\). Set \(\lambda=1/h\), so \(h=1\).
It remains to exclude the simultaneous vanishing of the first ten
moments in this normalized branch.

Here is an exact finite certificate for that exclusion. For the
normalized phase
\[
 y+h t^3+(v_0+v_1x)t^2+(dx+ex^2)t+\sum r_jx^j,
\]
associate the coefficient variables
\((h,v_0,v_1,d,e,r_0,r_1,r_2,r_3)\) with the triples
\[
\begin{aligned}
 \omega&=(1,1,2,2,3,1,2,3,4),\\
 \sigma&=(0,0,1,1,2,0,1,2,3),\\
 \tau&=(3,2,2,1,1,0,0,0,0).
\end{aligned}
\]
For every nonnegative exponent vector \(\alpha\) with
\(\omega\cdot\alpha=m\), let \(k=\sigma\cdot\alpha\) and
\(j=\tau\cdot\alpha\). Its coefficient in \(L(p^m)\) is zero for odd
\(j\), and otherwise equals
\[
 \frac{m!(j-1)!!}{(2k+j+1)!!\,\prod_i\alpha_i!}.
\]
Indeed, phase balance forces exactly \(k\) factors of \(y\); their
factorial in the multinomial coefficient cancels the \(k!\) in the
sphere monomial formula. This determines each moment as an exact rational
polynomial, without sampling coefficient values.

The first four moment equations have nonzero constant coefficients
\(1,2/3,2/5,8/35\) in \(r_0,r_1,r_2,r_3\), respectively. For \(h=1\)
their solutions are
\[
\begin{aligned}
 r_0&=-v_0/3,\\
 r_1&=-3/14-2v_0^2/15-v_1/5,\\
 r_2&=-13v_0/21-3d/7-124v_0^3/945-8v_0v_1/35,\\
 r_3&=-29/91-5318v_0^2/3465-173v_1/385-26v_0d/21-e/3\\
     &\hspace{1em}-2312v_0^4/14175-188v_0^2v_1/525-8v_1^2/75-d^2/2.
\end{aligned}
\]
Substitute these expressions into moments five through ten, and let
\(N_5,\ldots,N_{10}\) be their primitive integer numerators, retaining
the original signs. The retained certificate
[`two_pair_image_cubic_eleven_slice_certificate.json.gz`](two_pair_image_cubic_eleven_slice_certificate.json.gz)
contains six explicit integer polynomials \(A_5,\ldots,A_{10}\) and a
positive integer \(D\), together with the six \(N_i\), such that
\[
 \sum_{i=5}^{10}A_iN_i=D\ne0.
\]
All coefficients are explicit integers; compression merely encodes the
retained JSON. The checker derives the \(N_i\) independently from the
preceding coefficient formula, compares them with the certificate inputs,
and verifies this identity by exact polynomial multiplication over
\(\mathbb Q\). Hence the six residual moments cannot vanish together
in any characteristic-zero field. This excludes \(h\ne0\) and proves
the theorem.

**Replay and assurance.** The existing affine-slice checker now verifies
both the ten-dimensional eight-moment theorem and this eleven-dimensional
ten-moment theorem:

```bash
.venv/bin/python research/foundations/scripts/verify_two_pair_image_cubic_affine_slice.py
```

It checks every normalization and exceptional branch, the two univariate
Bézout identities, the first four moments against an independent direct
monomial expansion, and all six products in the retained multivariate
certificate. Adding one to a certificate multiplier is rejected. Neither
Singular nor the certificate-generation algorithm is needed for replay.
The check is internal exact replay, not external mathematical review or
formal verification. No exhaustion of the full thirteen-dimensional
frame slice or the sixteen-dimensional cubic family follows from this
result alone.

## The three fixed exceptional fields for repeated quadratic profiles

After a nonzero repeated quadratic profile has been scaled to
\((t-r)^2\), the orthogonal null rotation changes the coefficient of
\(t^3\) by \(-2k\). It therefore gives the normal form
\[
 p=y(t-r)^2+B_0+B_1t+B_2t^2
   +x(A_0+A_1t+A_2t^2)+x^2(D_0+D_1t)+Ex^3.
 \tag{RE1}
\]
The rotation preserves the displayed repeated profile and the degree
bound. For (RE1), the successive coefficients of \(B_0,A_0,D_0,E\)
in moments one through four are
\[
 1,\quad \frac{2P_2(r)}{15},\quad
 \frac{2P_3(r)}{105},\quad \frac{8P_4(r)}{15015},
 \tag{RE2}
\]
where
\[
\begin{aligned}
 P_2(r)&=5r^2+1,\\
 P_3(r)&=21r^4+18r^2+1,\\
 P_4(r)&=429r^6+715r^4+195r^2+5.
\end{aligned}
\]
These three polynomials are irreducible over \(\mathbb Q\) and pairwise
coprime. No inverse of a pivot in (RE2) is allowed on its zero locus.

**Proposition (fixed exceptional fields).** Suppose
\(rP_2(r)P_3(r)P_4(r)=0\). The first twelve spherical moments of (RE1)
vanish if and only if \(G=0\) in that normal form. In particular these
exceptional repeated-profile cases are one-sided after the stated linear
rotation, and satisfy eventual mixed vanishing.

The proof uses exact finite-field certificates with an explicit
characteristic-zero specialization argument. The finite fields are used
as a proof through positive weighted homogeneity; they are not a numerical
approximation or a rational reconstruction heuristic.

### Weighted projective specialization

Let \(K\) be a number field with a fixed finite prime, and let
\(f_1,\ldots,f_s\) be polynomials with coefficients integral at that
prime. Suppose all \(f_j\) are homogeneous for positive integer weights
\(w_1,\ldots,w_n\). If their reductions have no nonzero common zero over
the algebraic closure of the residue field, then the original polynomials
have no nonzero common zero over \(\overline K\).

Indeed, a nonzero common zero can be chosen algebraic over \(K\), by the
Nullstellensatz applied to an open chart where one coordinate is nonzero.
Extend the chosen valuation to that finite field and then, if needed, to
a further finite extension. For a solution \((z_i)\), choose \(\lambda\)
with valuation
\[
 v(\lambda)=-\min_{z_i\ne0}\frac{v(z_i)}{w_i}.
\]
Every \(\lambda^{w_i}z_i\) is integral and at least one is a unit.
Weighted homogeneity preserves every equation. Reduction therefore gives
a nonzero residue-field solution, a contradiction.

Equivalently, substituting \(z_i=Z_i^{w_i}\) turns the equations into
ordinary homogeneous equations with the same nonzero-solution condition
over algebraically closed fields. Their projective zero scheme is proper
over the localization of the number ring; a nonempty generic fibre must
specialize to a nonempty closed fibre. The projectivity-to-properness
input is [Stacks, Lemma 29.44.5, tag 01WC](https://stacks.math.columbia.edu/tag/01WC).
The valuation proof above supplies the needed implication directly.

### Complete exceptional branch split

The additional branch \(r=0\) is formed over \(\mathbb Q\), using residue characteristic 1000003. Eliminate \(B_0,A_0,D_0,E\) by the first four moments, then use moments five through twelve in the five active variables \(B_1,B_2,A_1,A_2,D_1\), with weights \(1,1,2,2,3\). Its exact residue basis has 219 elements, quotient dimension 797 and maximum standard weight sixteen. The pure powers \(B_1^{17},B_2^{17},A_1^9,A_2^9,D_1^6\) all reduce to zero. This proves this rational branch by the same specialization lemma.

All other moment equations and eliminations below are formed in
\(K_j=\mathbb Q[r]/(P_j)\), for \(j=2,3,4\). Since each \(P_j\) is
irreducible, this treats every one of its complex roots by conjugacy.
The exact residue fields and active variables are:

| Branch | Residue characteristic | Variables after elimination | Positive weights |
|---|---:|---|---|
| \(P_2=0\) | 1013 | \(B_1,B_2,A_0,A_2,D_1\) | \(1,1,2,2,3\) |
| \(P_3=0\) | 10007 | \(B_1,B_2,A_1,A_2,D_0\) | \(1,1,2,2,3\) |
| \(P_4=0\) | 1031 | \(B_1,B_2,A_1,A_2,D_1,E\) | \(1,1,2,2,3,4\) |

For \(P_2=0\), use moments one through four to eliminate
\(B_0,A_1,D_0,E\). The replacement second-moment pivot is
\(-4r/15\ne0\). For \(P_3=0\), eliminate \(B_0,A_0,D_1,E\);
the replacement third-moment pivot is
\(-8r(1+3r^2)/105\ne0\). For \(P_4=0\), eliminate only
\(B_0,A_0,D_0\) and retain \(E\) as a sixth variable. The fourth
moment is then a residual equation rather than an elimination.
Coprimality of the three pivot polynomials justifies all remaining
inversions. The verifier also checks the actual pivots in each residue
field, so none of these inversions silently reduces to zero.

Each residual moment of order \(m\) is exactly homogeneous of weight
\(m\). The residue ideal is generated by moments five through twelve in
the first two branches and by moments four through twelve in the last
branch. Direct exact Groebner computations give respectively 244, 244,
and 294 reduced basis elements; all three quotients have dimension 797
as residue-field vector spaces. Their standard monomials have maximum
weight sixteen. The verifier additionally reduces the following literal
pure powers to zero:
\[
\begin{array}{c|l}
P_2=0&B_1^{17},B_2^{17},A_0^9,A_2^9,D_1^6,\\
P_3=0&B_1^{17},B_2^{17},A_1^9,A_2^9,D_0^6,\\
P_4=0&B_1^{17},B_2^{17},A_1^9,A_2^9,D_1^6,E^5.
\end{array}
 \tag{RE3}
\]
Thus each residue ideal has only the origin as a common zero over the
algebraic closure. The verifier checks irreducibility of the same-degree
minimal polynomial over each finite field, all rational coefficient
denominators as units at the chosen prime, the positive weight identity
for every monomial, and every pure-power normal form in (RE3).
The weighted specialization lemma now proves that all active parameters
also vanish over \(K_j\). Every eliminated coefficient is homogeneous
of positive weight in those active parameters, so it vanishes as well.
This proves the proposition.

The exact checker is
`research/foundations/scripts/verify_two_pair_image_cubic_repeated_exceptions.py`.
It rederives the moments by the complete multinomial phase selection
\(n_y=n_A+2n_D+3n_E\), then the exact sphere beta functional. It uses
SymPy 1.14.0 and Singular 4.3.2 and writes the branch formulas, input
hashes, coefficient-domain checks and pure-power certificates to the
retained verification JSON. The negative control
\(y(t-r)^2+xt\) has second moment \(-4r/15\ne0\) on all three fields.

These fixed-frame certificates provide additional checks of the complete
classification. They prove the stronger statement \(G=0\) in precisely
the displayed parameter fields. The full Image theorem uses the existence
of an alternative constant-profile frame and requires no uniform
fixed-frame zero-coefficient assertion.

### An additional exact generic minor

The generic repeated-profile verifier works over \(\mathbb Q(s)\),
\(s=r^2\), with \(B_1=rb_1,A_1=ra_1,D_1=rd_1\). After the first four
pivots are inverted, residual moments through twelve are homogeneous of
weights \((1,1,2,2,3)\) in \((b_1,B_2,a_1,A_2,d_1)\). At weight eighteen,
form their integer-polynomial Macaulay coefficient matrix. It has 874 rows
and 1187 columns. The retained selector chooses an 874-square minor
\(D(s)\); the exact definition uses the rederived integer matrix, not a
floating-point determinant. Its value at \(s=1\) is 42 modulo 32003, so
\(D\) is a nonzero characteristic-zero polynomial. For every specialization
where \(sP_2P_3P_4D(s)\ne0\), all weight-eighteen monomials lie in the
moment ideal. This includes the pure powers
\(b_1^{18},B_2^{18},a_1^9,A_2^9,d_1^6\), so the active coefficients vanish.
The determinant certificate alone does not cover the roots of \(D\);
its precise scope is retained independently of the full geometric proof.

The [generic verifier](../scripts/verify_two_pair_image_cubic_repeated_generic.py),
[selector](two_pair_image_cubic_repeated_minor_selector.json), and
[receipt](two_pair_image_cubic_repeated_generic_verification.json) retain the
matrix generation, exact column choice, integer-input digest, degree bound,
and modular nonvanishing witness. The modular witness proves that this
specified integer determinant is nonzero; it is not used to identify its
characteristic-zero roots.

## Reproduction and assurance

Use the existing environment pinned in
[research/requirements.txt](../../requirements.txt), from the repository root:

```sh
.venv/bin/python research/foundations/scripts/verify_two_pair_image_cubic_structure.py
.venv/bin/python research/foundations/scripts/verify_two_pair_image_cubic_affine_slice.py
.venv/bin/python research/foundations/scripts/verify_two_pair_image_cubic_phase_slice.py
```

The [structure checker](../scripts/verify_two_pair_image_cubic_structure.py)
verifies the full coefficient maps, sphere identities, nine-parameter quadratic
criterion, graph and blow-up formulas, ordinary-double-point Hessian, residue
cancellation, and the exceptional shear identities. The geometric and
\(V\)-filtration arguments are written proofs; their analytic/algebraic
D-module source statements are not certified by this Python run. The
[primary-source audit](two_pair_image_cubic_primary_sources.json) records the
exact version, theorem locators, retained PDF and inspected hypotheses.

The [affine verifier](../scripts/verify_two_pair_image_cubic_affine_slice.py)
checks the complete ten- and eleven-parameter slices, all normalization and
zero-pivot branches, and the exact six-product integer Bezout certificate
in [compressed JSON](two_pair_image_cubic_eleven_slice_certificate.json.gz).
The [phase-slice verifier](../scripts/verify_two_pair_image_cubic_phase_slice.py)
checks the sharp seven-parameter six-moment theorem and the opposite quadratic
profiles. Their receipts distinguish exact identities from bounded regression.

Additional fixed-frame arithmetic checks use Singular 4.3.2 and the existing
SymPy 1.14.0 environment:

```sh
.venv/bin/python research/foundations/scripts/verify_two_pair_image_cubic_repeated_exceptions.py
.venv/bin/python research/foundations/scripts/verify_two_pair_image_cubic_repeated_generic.py
```

The full result does not assume a universal finite pure-moment prefix. Its
global geometric step uses all-order purity, while finite prefixes certify
their stated slices. The ordinary and exceptional normal forms give all-order
mixed vanishing with the displayed uniform degree bound.

Internal agents separately audited the proper-seed argument, all boundary
charts, complete case coverage, and the exact eleven-parameter certificate.
These are internal mathematical checks and replays. External mathematical
review and formal verification remain pending. The
[completion audit](two_pair_image_cubic_completion_audit.json) records the
full objective and the evidence for each requirement.
