# Finite images of quadratic radial Gaussian phases

Let (K) be a field of characteristic zero. Use a Witt Gaussian pair
\(Z,W\) and an independent standard Gaussian \(Y\), with
\(E[Z^iW^j]=\delta_{ij}i!\). Put \(S=ZW\). Let

\[
 P=P_0+R,\qquad P_0=A(Y)S^2+H(Y)S+G_0(Y),
\]

where \(A,H,G_0\) are arbitrary polynomials and the monomials of \(R\)
have strictly one sign of \(Z\)-minus-\(W\) weight. This includes \(R=0\).
Write \(\Phi_P[Q]=E[Qe^{tP}]\) for the periods of the original polynomial
markers, \(V_P\) for their linear image, and \(M_P=K[t]V_P\).

**Theorem.** The full actual image \(M_P\) is finite over \(K[t]\) exactly
when \(A=H=0\) and \(G_0=c+bY\). In this case the unit is primitive,
\(\Phi_P[1]=\exp(ct+b^2t^2/2)\), and \(M_P=\Phi_P[1]K[t]\).
For nonconstant \(P\), \(V_P=M_P\). For \(P=c\), \(V_P=K e^{ct}\).
Safety is equivalent to \(P_0=0\), so safety implies purity. For
\(G_0=0\), polynomiality of the original unit already forces \(P_0=0\),
and is equivalent to purity and safety.

The new step excludes the scalar identity
\(\Phi_{P_0}[1]=p(t)\exp(ct+b^2t^2/2)\) for \(\deg A\ge2\) and
\(G_0=c+bY\), with arbitrary \(p\in K[t]\), \(p(0)=1\). This is the
normalization forced by full finiteness. The unrestricted scalar exponent
for coefficients with higher root multiplicity remains open. The separate
[logarithmic-cover theorem](GAUSSIAN_QUADRATIC_RADIAL_LOGARITHMIC_COVER_FINITE_IMAGES.md)
supplies the arbitrary-exponent result for roots of order at most two.

## 1. The original relative source

The [angular theorem](GAUSSIAN_ANGULAR_CORE_FULL_IMAGE_INVARIANCE.md)
preserves both the entire module and the unit when passing to \(P_0\).
Balanced original monomials are exactly \(S^iY^j\), with factorial
moments in \(S\). Full finiteness of the original module forces affine
\(G_0=c+bY\), with unit \(pe^{ct+b^2t^2/2}\), by the
[boundary theorem](GAUSSIAN_INVARIANT_CORE_BOUNDARY_REDUCTION.md).
The [affine-coefficient theorem](GAUSSIAN_QUADRATIC_RADIAL_AFFINE_COEFFICIENT_FINITE_IMAGES.md)
handles \(A=0\) and \(\deg A\le1\), including their zero-boundary unit
criterion. We prove the new scalar exclusion over \(\mathbb C\).

Put \(r=1/t\), \(h=H-r\), and gauge the exponential by
\(-ct-b^2t^2/2\). On the normalization of \(V^2=A(Y)\), with a consistent
choice of \(\sqrt r\), write

\[
 G=-\tfrac12(Y-b/r)^2,\quad x=\frac{h}{2\sqrt r V},\quad F=G-x^2,
 \qquad \eta=\frac{\sqrt r}{V}\,dY.                       \tag{1}
\]

The exact relative Gamma source, including the time connection, is the
one in Sections 2–3 of the logarithmic-cover theorem. Its derivation
there uses arbitrary polynomial \(A\). In particular, with spatial
exponential derivatives \(D_S,D_Y\), the original quotient is by
\(D_Y\mathbb C(t)[S,Y]+D_S(S\mathbb C(t)[S,Y])\). On \(A\ne0\) the
relative \(S\)-cohomology has basis \(v=[dS]\),
\(e_B=[(1-tH-2tAS)dS]\). Scaling \(v\) to \(u=Vv/\sqrt r\) gives the
connection \(\nabla e_B=dG\,e_B\), \(\nabla u=dF\,u-dx\,e_B\).
Relative Gysin at the reduced support of \(A=0\) injects the original
source into the localized source. Its kernel term is relative \(H^0\)
of the union of \(S\)-lines with their endpoints, and vanishes. The
finite étale cover over \(A\ne0\) has a constant trace projector, also
when the normalization is disconnected. On the cover the boundary
line is invariant and \(x,u\) are odd. These operations give a
differential embedding of the original source. It needs no formal
evaluation of inverse spatial markers.

In what follows the exact sequence on the cover, followed by the
appropriate deck projector, is

\[
 0\longrightarrow H^1(G)\longrightarrow H_{\rm rel}^2
   \longrightarrow H^1(F)\longrightarrow0.               \tag{2}
\]

Both phases are spatially nonconstant, so their zeroth exponential
cohomologies vanish. The localized curve is affine, and the relative
\(S\)-complex is concentrated in cohomological degree one. This gives
(2) by the two-term curve complexes. The relative rapid-decay pairing
is the natural dual cone of the two absolute pairings for the total
exponential pair. The pairing and parameter comparison used here are
the ones in Section 3 of the logarithmic-cover theorem. No generic
semisimplicity of the extension is required.

## 2. Uniform bounds on the region of small phase

Use absolute values on a fixed parameter sector, taking compatible
ramifications for rational powers of \(r\). Constants may depend on the
fixed phase and sector. For sufficiently small rational \(\epsilon>0\)
we show that every descending inverse-phase arc, truncated to
\(|F|\le |r|^{-\epsilon}\), has \(\int |\eta|\) bounded by
\(C|r|^\delta(1+|\log |r||)\) for some \(\delta>0\). The same bound
holds on the bounded-phase compact and neck arcs in the marked curve
frame. The length property used in the estimates is that a uniformly
bounded number of arcs of fixed real algebraic degree has length
\(O(\rho)\) inside a disc of radius \(\rho\). It applies also after
a fixed algebraic coordinate change and bounded multiplicity pullback.
This follows from the [curve-length lemma](../../../unpublished/gaussian-boundary/analytic_core.tex)
by the uniform algebraic degree bound on projections. Inverse images
of straight or circular phase arcs have this property. For compact and
neck chains in the marked model, use finitely many algebraic charts of
the marked curve family, with straight or circular arcs in local
coordinates on the fibres. Concretely, in smooth charts choose a
spatial coordinate étale with the parameter, and lift finitely many
short straight segments between local sections. Use radial and circular
arcs in the algebraic node coordinates on a neck. After a finite
parameter cover all required end sections are algebraic. Projecting
these finitely many fixed algebraic constructions to \(Y\) preserves
bounded degree. Retracing on inverse-phase pieces is bounded by the
degree and the number of phase segments.

Only finite roots of \(A\) of order \(m>2\) need new local estimates.
At the other finite points and at infinity the fixed (dY/V) is
regular or logarithmic. At infinity this follows from the total
degree of \(A\), regardless of finite multiplicities. In a small fixed
disc about a higher root \(a\), put \(z=Y-a\) and write
\(A=z^m a_0(z)\) with \(a_0(0)\ne0\). Thus
\(|\eta|\asymp |r|^{1/2}|z|^{-m/2}|dz|\) sheetwise.

### Zero boundary slope

For \(b=0\), \(G\) is bounded in a fixed finite disc. In the region
\(|F|\le |r|^{-\epsilon}\), equation (1) yields
\[
 |h|^2\le C |r|^{1-\epsilon}|z|^m.                       \tag{3}
\]
Where \(|h|\ge |r|/2\), we have \(|z|\ge c |r|^{(1+\epsilon)/m}\).
Dyadic annuli and the length bound therefore give volume order at
least \(1/m-\epsilon(1/2-1/m)\). Where \(|h|<|r|/2\), one has
\(|H|\asymp|r|\). This region only occurs at a root of \(H\), of order
\(n\ge1\). Here \(|z|\asymp |r|^{1/n}\) and \(|H'|\asymp
|r|^{1-1/n}\). Use \(h\) as an algebraic coordinate on the branches
in this annulus. Equation (3) places its image in a disc of radius
\(C|r|^{(1+m/n-\epsilon)/2}\). The volume is bounded by the length in
that disc times \(C|r|^{1/2-m/(2n)-(1-1/n)}\), hence its order is at
least \(1/n-\epsilon/2\). Both bounds are positive for small
\(\epsilon\). When \(H=0\), \(h=-r\) and the first bound applies.

### Nonzero boundary slope

For \(b\ne0\) and bounded \(Y\), the leading term of \(G\) is
\(-b^2/(2r^2)\). For \(\epsilon<1\) the truncated region satisfies
\[
 |h|^2/|z|^m\asymp |r|^{-1}.                              \tag{4}
\]
If \(H(a)\ne0\), (4) gives the annulus \(|z|\asymp |r|^{1/m}\),
whose unrefined volume order is \(1/m\). Otherwise let \(n\ge1\) be
the order of \(H\) at \(a\). The following annuli and coordinate discs
cover the truncated region.

1. If \(|H|\ge2|r|\), then \(|h|\asymp|z|^n\), giving
   \(|z|^{m-2n}\asymp |r|\). Since \(|z|\) stays bounded and
   \(|z|^n\ge c|r|\), this forces \(m\ge3n\). For \(m>3n\),
   scale \(z=r^{1/(m-2n)}u\) with \(|u|\asymp1\). The leading
   function of \(r^2F\) is a nonzero multiple of \(u^{2n-m}\) plus
   \(-b^2/2\). Its nonzero roots are simple.
2. If \(|H|\le |r|/2\), then \(|h|\asymp |r|\), giving
   \(|z|\asymp |r|^{3/m}\) with \(m\le3n\). For \(m<3n\), scale
   \(z=r^{3/m}u\). The leading function is a nonzero multiple of
   \(u^{-m}\) plus \(-b^2/2\), again with simple roots. The same
   annulus and leading function apply if \(H=0\).
3. At equality \(m=3n\), or if \(|H|\asymp|r|\) with
   \(m=3n\), use \(z=r^{1/n}u\). Write \(w=H_n u^n\), with
   \(H_n\) the initial coefficient of \(H\) at \(a\). The limiting
   function of \(r^2F\) is a nonzero multiple of
   \((w-1)^2/w^3\) plus \(-b^2/2\). Its derivative with respect to
   \(w\) is proportional to \((w-1)(3-w)/w^4\). A zero is never
   \(w=1\), and a zero at \(w=3\) has multiplicity two. Thus every
   zero on the compact unit annulus has multiplicity at most two.
   The double case occurs, for example with \(m=3\), \(H=Y\),
   \(A=-2Y^3/27\), \(b=1\) and \(a=0\). There the limit is
   \(-(u-3)^2(4u-3)/(8u^3)\).
4. If \(|H|\asymp|r|\) and \(m\ne3n\), (4) gives
   \(|h|\asymp |r|^L\), where \(L=(m/n-1)/2\). Since \(|h|=O(|r|)\),
   we must have \(m>3n\) and \(L>1\). Use \(h=r^L u\) on branches
   near the algebraic roots of \(H(Y)=r\). Here \(|z|\asymp|r|^{1/n}\)
   and \(r^2F\) tends to a nonzero multiple of \(u^2\) plus
   \(-b^2/2\) on each branch. Its roots are simple.

These facts about the limit functions yield uniform small discs for
the actual truncated phase. On a scaled compact annulus, put
\(f_r=r^2F\). Near a zero of multiplicity \(j\) of the limit, factor
\(f_r\) as a nonvanishing bounded unit times the product of its
\(j\) moving linear factors. This follows by local Weierstrass
preparation after ramification. Outside these neighborhoods the
limit is bounded away from zero. Therefore
\(|f_r|\le |r|^{2-\epsilon}\) puts every point within
\(C|r|^{(2-\epsilon)/j}\) of at least one factor root. This includes
all subleading cancellations and splittings, without assuming a simple
critical point. The length bound applies inside the resulting discs.

The preliminary volume orders in cases 1–4 are respectively
\((1-n)/(m-2n),\ 3/m-1,\ 1/n-1,\ 1/n-1\). The first two are
strictly greater than \(-1\). The last uses
\(dY=r^L du/H'\). Add \(2-\epsilon\) in each simple case, and at
least \(1-\epsilon/2\) in case 3. All resulting orders are strictly
positive when \(\epsilon\) is sufficiently small. The proof treats a
whole truncated inverse image, rather than a list of stationary points.

## 3. Absolute chains and the complete tails

Use the complete marked curve comparison of Section 5 of the
logarithmic-cover theorem. Its description of the cycles uses the
phase and the punctures. The differential enters the subsequent
estimates. Mark all punctures and the zero, one and pole sections of
\(F\), and separate them after a finite time cover. The proper model
and nearby-cycle comparison retain every compact cell, finite/finite
annulus, descending attachment, and their monodromy. The regular
graded part has flat representatives with bounded phase on their
compact and neck pieces and descending tails of bounded imaginary
phase. The positive-genus and split-cover cases use the same
comparison on each component. The estimates in Section 2 apply to
these geometric chains and to incomplete inverse paths beginning
with bounded real part and bounded imaginary part of \(F\).

Here is the justification for truncating a descending path at
\(|F|=|r|^{-\epsilon}\). Its complement has \(\Re F\le
-c|r|^{-\epsilon}\), with the same conclusion if the path already
starts there on the descending side. Small detours of bounded radius
about exceptional phase values preserve the estimate. In fixed
coordinate discs use polynomial cutoffs \(|r|^M\) containing all
critical, starting, transition and marking scales. Off the deeper
tails the volume has at most polynomial growth in \(|r|^{-1}\), by
the dyadic length estimate. On the deeper tails the actual polar
term and its derivative dominate all other terms, including the
scalar offset of order two. For a finite \(A\)-root of order \(m\),
the leading phase is \(C(r)z^{-m}\) on the base and
\(C=-(H(a)-r)^2/(4r a_0(0))\). Thus
\[
 \left|\frac{\eta}{dF}\right|\le
 C_1|r|^{1/2}|C(r)|^{1/m-1/2}|F|^{-1/2-1/m}.             \tag{5}
\]
This holds on the ramified branches as well. At infinity the fixed
differential is regular or logarithmic and the usual inverse-pole
bound is at most polynomial in \(|r|^{-1},|F|\). The phase has a
true pole there because the term \(-Y^2/2\) in its (r^0) part
cannot identically cancel in the Laurent parameter; more explicitly,
the (r^0) part is \(-Y^2/2+H/(2A)\), and cancellation of its
top term would make \(H^2/A\) have still higher positive degree.
All finite punctures are phase poles over \(\mathbb C(r)\). The
deeper tails are integrated against the actual descending exponential
starting at the cutoff. The polynomial losses in (5) and the other
pieces are dominated by \(\exp(-c|r|^{-\epsilon})\).

At finite/finite necks the truncated bound controls the whole bridge,
including shrinking annuli. Logarithmic periods from monodromy
preserve a positive power. Hence every regular graded seed period of
the absolute quotient tends to zero. Each open descending filling
with the indicated initial phase bounds also has seed tending to
zero, regardless of whether its initial action is exponentially
negative. These estimates are on absolute integrals of the chain
pieces. They also give polynomial growth for every fixed algebraic
amplitude on these chains. Away from the deeper tails its additional
poles cost finite powers of \(|r|^{-1}\), and on the deeper tails an
algebraic inverse of the actual phase gives polynomial growth in
\(|F|\) and \(|r|^{-1}\), which the exponential absorbs.

## 4. Lifting the boundary cycles with controlled growth

On a generic ray at (r=0), the solutions with at most polynomial
growth form the zero step of the Stokes filtration. The strictly
smaller step on that ray consists of the exponentially decaying
solutions, and the quotient is the regular graded part. Use the
Stokes correspondence and strictness in Sabbah, Theorem 3.1 and
Proposition 5.12, in the [retained primary text](sources/sabbah_stokes_0912.2762v5.pdf), in the local version
for meromorphic connections on a disc. In particular the dual solution
sequence of (2) is exact on the moderate steps. The same fact applies
to the restriction to an embedded cyclic submodule. We use growth
of the full vector of periods to place a lift in the filtration.

The invariant boundary source on the punctured (Y\)-line has its
Gaussian contour and one loop per distinct (A\)-root. The loops pair
by residues with polar forms, and restricting from the unpunctured
line gives the one usual Gaussian class. For \(b=0\) choose fixed
contours and loops avoiding the (A\)-roots. Their Gaussian phases
have bounded imaginary part and real part bounded above. Every
boundary solution is moderate. For \(b\ne0\) choose a generic ray
where \(\Re(-b^2/(2r^2))>c|r|^{-2}\) and
\(|\Im(b/r)|>c|r|^{-1}\). This is an open set of directions. Every
puncture factor (e^{G(a)}\) then grows exponentially. The translated
Gaussian contour \(Y=b/r+\mathbb R\) gives the entire moderate boundary
line. Indeed the pole loops have evaluations (2\pi i e^{G(a)}\) on
the respective simple-pole form and zero on the other simple-pole
forms, whereas the translated contour is moderate on every form.

On the chosen paths lift the endpoint by one of the Gaussian tails
from Section 7 of the logarithmic-cover theorem. In the completed
coordinate \(U=VS/\sqrt r+x\), use the half-line \(U=x\pm is\),
\(s\ge0\), choosing the sign of \(\Im x\). Its exponential relative
to the endpoint has absolute value at most \(e^{-s^2}\). Hence the
endpoint density is \(\eta e^G J_\pm(x)\), with the uniform Gaussian
tail bound \(|J_\pm(x)|\le C/(1+|x|)\). On the fixed contours this
gives \(O(\sqrt{|r|})\) for the seed, and on the translated contour
also \(|Y|\ge c|r|^{-1}\). The same positive bound follows there,
using \(|1/V|=O(|Y|^{-\deg A/2})\) and the translated Gaussian weight.

When \(x\) crosses the real axis, the two half-lines differ by a
closed Gaussian thimble in relative homology. Attach that thimble
over a descending inverse-(F\) path to a phase pole. At the switch
\(\Re F=\Re G-x^2\le C\), and \(\Im F=\Im G\) is bounded. There
are uniformly finitely many switches from the fixed-degree algebraic
equations on the contour pieces, with identically real segments
handled by one choice throughout. The starts are algebraic along
the ray. Leftward phase rays with finite detours of bounded radius
can be chosen continuously on a generic sector, with the ramified
inverse lift used if a start is critical. They reach actual poles,
since the nonconstant meromorphic phase is finite as a map to the
projective line. Section 3 bounds the seed of every closing surface
by a positive power times finite logarithms, including starts with
small \(x\). In completed coordinates closing surfaces are products
of an absolute Gaussian chain with an open descending path. Their
ends have rapid decay in both variables. The glued chains vary flatly
because they vary continuously as relative cycles on the sector. Under
the deck change the half-density changes sign along with
\(J_-(-x)=-J_+(x)\), so they represent the original boundary lifts.

These lifts are moderate as solutions on the entire localized source.
In the two-element spatial basis their other periods are given by
fixed rational amplitudes times the Gaussian and endpoint densities.
On the chosen contour and loops those amplitudes have at most
polynomial growth, and Gaussian tails absorb the spatial powers at
infinity. On the closing surfaces Section 3 gives the corresponding
polynomial bound for arbitrary fixed rational amplitudes. This ensures
that their images under the strict moderate restriction map span the
moderate boundary space. Any moderate functional on the relative
source differs from their constant linear combination by a moderate
absolute solution. The latter has decaying seed using Section 3 for
its regular graded part, and its strictly decaying part tends to zero
on the seed as well. This proves decay for the whole moderate step of
the relative system.

## 5. Exclusion and consequences

A hypothetical actual unit (pe^{ct+b^2t^2/2}\) supplies the original
formal functional whose value on (v/p\) is one. Its restriction to the
cyclic differential submodule spanned by covariant iterates of (v/p\)
is rational, taking these iterates to (1,0,0,\ldots\). The horizontal
functional on this submodule is moderate on every generic sector.
Strictness of the restriction from the dual localized module extends
it to a moderate functional on the chosen sector. The ordinary-point
relative period pairing identifies this extension with a flat cycle.
On (v/p\) its value stays one. Sections 3–4 show instead that its
value tends to zero: division by (p(t)\) adds the nonnegative order
\(\deg p\). This proves the asserted scalar exclusion in the prescribed
boundary normalization. All coefficients of a hypothetical identity
over (K\) lie in a finitely generated field embeddable in
\(\mathbb C\), and the formal identity transfers coefficientwise.

The finite-image necessity and zero-boundary claim now follow from
Section 1. For (P_0=c+bY\), let (\gamma\) be the least absolute
weight of a monomial of (R\ne0\). For each original (Q\) of degree
\(D\), the exact finite formula is
\[
 \Phi_P[Q]=e^{ct+b^2t^2/2}\sum_{j=0}^{\lfloor D/\gamma\rfloor}
 \frac{t^j}{j!}E[Q(Z,W,Y+bt)R(Z,W,Y+bt)^j].              \tag{6}
\]
It follows from angular balance and Gaussian translation. If (R=0\)
only (j=0\) occurs. The constant marker gives the primitive vacuum,
and the [original fixed-marker theorem](GAUSSIAN_AMPLITUDE_IDEALS_AND_WEIGHT_INVARIANCE.md#7-every-finite-actual-image-and-its-original-fixed-marker-image)
gives the stated literal image. The [recomposition theorem](GAUSSIAN_SAFE_RECOMPOSITION_CLASSIFICATION.md)
gives the ideal ((U-c)\) for constant core (c\), and the zero ideal
otherwise. In the safe branch (6) gives the cutoff \(\lfloor D/\gamma\rfloor\).

On every fixed zero-boundary coefficient space with arbitrary finite
degree bounds for (A,H\), the field pure locus is now the origin.
The unit moment ideal has finite colength by the Nullstellensatz.
Noetherianity gives some finite prefix containing a power of the
coefficient maximal ideal. Its actual containment certificate gives
all homogeneous moment tails and the original-marker ring cutoff as
in Section 9 of the affine-coefficient theorem. The next radial-degree
case has a genuine (S^3\) term. Already in original degree six the
remaining support after boundary reduction and the quadratic-(Y\)
classification is \(a S^3+A_2(Y)S^2+H_4(Y)S+c+bY\), with
\(a\ne0\), \(\deg A_2\le2\), and \(3\le\deg H_4\le4\).

## Evidence

The [checker](verify_gaussian_quadratic_radial_unrestricted_finite_images.py)
verifies the exact scaled phase functions, the double resonant face,
the seed Jacobians and tail powers, and bounded original relative
Wick identities with repeated roots. The [receipt](gaussian_quadratic_radial_unrestricted_finite_images_verification.json)
pins these controls, the proof and its comparison inputs. The uniform
chain and filtration arguments are proved in the text. The confidence
scope is an internal written proof with exact bounded controls.
