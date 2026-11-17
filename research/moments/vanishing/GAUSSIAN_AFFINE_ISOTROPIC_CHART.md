# Gaussian H(2,1): affine charts, radial phases and obstructions

Internal written proof, 17 November 2026. This extends the algebraic cubic
subcase in [the boundary manuscript](../../../unpublished/gaussian-boundary/manuscript.tex),
subsection "A complete algebraic subcase", to arbitrary polynomial degrees.
The mixed-source extension below also uses the one-base manuscript's formal
leading-coefficient obstruction. Independent mathematical review,
formal verification and external review remain pending. It does not settle
the general constant-fibre question H(2,1).

## Statement

Let U1,U2,Y be independent standard Gaussians, and set
Z=(U1+i U2)/sqrt(2), W=(U1-i U2)/sqrt(2). Thus

\[
 E[Z^pW^q]=\delta_{pq}p!,\qquad E[Y^{2j}]=(2j-1)!!.
\]

All expectations below mean coefficientwise Gaussian contraction of formal
series. For arbitrary h,v,b,k in C[z], consider

\[
 P=W h(Z)+v(Z)Y^2+b(Z)Y+k(Z),\qquad
 \Phi_Q(t)=E[Qe^{tP}].
\]

Call P pure when Phi_1=1, and safe when Phi_Q is polynomial for every
fixed polynomial Q in Z,W,Y. Exactly the following disjoint forms are pure:

1. **Positive weight:** h is divisible by z^2, and v,b,k are divisible
   by z. These phases are safe, with deg_t Phi_Q <= deg Q.
2. **Shear:** h=c is a nonzero constant, v=0, b is arbitrary, and
   k=-z b^2/(2c). These phases are safe. With d=max(0,deg b), using
   d=0 for b=0, a sufficient bound is

   \[
   \deg_t\Phi_Q\le(2d+1)\deg Q.
   \]

3. **Rational marked period:** h=c+dz with c,d nonzero, and for an
   arbitrary polynomial ell,

   \[
   v=\frac{h-h^3/c^2}{2z}
     =-d-\frac{3d^2}{2c}z-\frac{d^3}{2c^2}z^2,\qquad
   b=h^2\ell,\qquad k=-\frac{c^2}{2}zh\ell^2.
   \]

   These phases are unsafe:

   \[
   \Phi_Z=\frac{ct}{1-dt},\qquad
   E[ZP^m]=c d^{m-1}m!\quad(m\ge1).
   \]

The third branch has no nonisotropic constant-quadratic direction: if a
constant vector e satisfies D_e^2 P in C, then e is a multiple of the
isotropic W direction. Consequently no affine complex orthogonal change
of variables puts a phase in this branch in H(2,1).

In particular, in the constant-fibre slice

\[
 P=aY^2+b(Z)Y+W h(Z)+k(Z),
\]

purity forces a=0 and one of the first two forms. This proves H(2,1)
for this all-degree subfamily. It imposes both B=b(Z) and an affine
dependence of g on W; arbitrary B(Z,W),g(Z,W) remain outside its scope.

The later sections prove further obstructions for the mixed family
P=W(h(Z)+c(Z)Y)+aY^2+u(Z)Y+k(Z). Every pure member is safe when
deg c<=8, with no degree bound on h,u,k. More generally, nonconstant c
with no root of multiplicity greater than two is pure exactly in the
positive-weight zero-root form (8). Nonconstant squarefree c of any
degree therefore admits no pure members. Arbitrary root multiplicities
are also excluded when a is nonzero or when c,h are coprime. A further
all-degree classification allows arbitrary multiplicities whenever h
vanishes at least as strongly as c at each common root. In particular
it classifies c|h, including h=0. The sharp four-level calculation
extends the classification to all common-root orders outside
ceil(2m/3)<=ord h<=m-3, where m=ord c. In particular every common
root of multiplicity at most eight is covered. Adding a
constant nonzero quadratic W term forces every pure member to be
quadratic, even with an arbitrary polynomial WY coefficient. The
nonconstant-coefficient exclusion uses a connected algebraic double
cover and positive differential order, retaining its positive-genus
component cycles. A relative endpoint argument proves
that every radial pure phase is zero in all degrees, and consequently
settles weakly one-sided angular support in all degrees. The polynomial
quadratic-form automorphism orbit of Long's quartic is excluded as a
source of H(2,1) counterexamples. A further all-degree exclusion allows
a polynomial coefficient ell(Z) of W^2 with deg ell>=2 and simple/double
roots in gcd(ell,c^2), subject to the explicit clean-source conditions
below. It includes gcd(ell,c)=1 without multiplicity restrictions and
the c=0 simple/double-root case with ell(0)!=0 or h(0)!=0. Acyclicity
on the localization divisor preserves the original seed.
The exceptional residual family
Z^(2nu+1)W^2+Z^nu WY-W/2+k(Z) has no pure member for any integer nu>=1
or polynomial k, by an elementary all-order coefficient-growth argument.
A six-parameter mixed quartic family is completely classified by an exact
algebraic certificate followed by an all-order weight proof. Polynomial
cubic-phase examples realize nonpositive normalized-volume order for
every integer pair in the remaining common-root strip. They show why
strict decay cannot be deduced from the W-affine shape alone.

## Exact circular elimination

Let r in t C[[t]] be the unique solution r=t h(r), and let
J=1-t h'(r), K=1-2t v(r). Both J and K have constant term one.
For any numerator A(z,t), polynomial in z at each t order,

\[
 E[A(Z,t)e^{tWh(Z)}]=\frac{A(r,t)}{J}.
\]

Indeed, the left side is

\[
 \sum_{n\ge0}t^n[z^n]\bigl(A(z,t)h(z)^n\bigr).
\]

Expand 1/(z-th(z)) as a geometric series in t. Its formal residue against
A(z,t) is that same sum. The factorization
z-th(z)=(z-r)u(z,t), where u is a t-adic unit and u(r,t)=J,
then gives the displayed residue. All operations at a fixed t order
involve only finitely many coefficients.

Eliminating Y next gives the exact identity

\[
 \Phi_1=\frac{\exp(R)}{J\sqrt K},\qquad
 R=t k(r)+\frac{t^2 b(r)^2}{2K},                         \tag{1}
\]

with the square root having constant term one. Likewise
Phi_Z=r Phi_1. No convergence of a real exponential integral is used.

We use an elementary rational-exponential fact. If rational functions
F,G over C, with G nonzero, satisfy exp(F)=G as formal germs at a
regular point, then F and G are constant. Differentiation gives
F'=G'/G. At a finite pole of F, its derivative has pole order at
least two; a logarithmic derivative has at most simple poles. Hence
F has no finite poles and is polynomial. At infinity G'/G=O(1/z),
so the polynomial F' vanishes. This argument proves the fact without
analytic continuation of the exponential.

## Necessity and completeness

Suppose first h(0) is nonzero. The formal substitution t=z/h(z) is
invertible near zero, with r=z. Set

\[
 A=h(z),\qquad C=h(z)-zh'(z),\qquad D=h(z)-2zv(z).
\]

Squaring (1) under purity gives an exponential of a rational function
equal to a rational function. The preceding fact, and the value at zero,
force R=0 and J^2 K=1. Equivalently,

\[
 C^2D=h^3,\qquad 2Dk+z b^2=0.                           \tag{2}
\]

These identities imply deg h<=1 in every degree. To prove this, write
the distinct roots of h as alpha_i, with multiplicities m_i. All alpha_i
are nonzero. If n=deg h>=2, reduce

\[
 \frac{C}{h}=1-z\sum_i\frac{m_i}{z-\alpha_i}
            =\frac{N(z)}{\prod_i(z-\alpha_i)}.
\]

The numerator has the same positive degree as the denominator, with
leading coefficient 1-n. It is nonzero at every alpha_i, since its
value there is -alpha_i m_i product_{j!=i}(alpha_i-alpha_j).
It therefore has a zero away from the zeros of h. This contradicts
C^2 dividing h^3 in (2). If h is constant c, (2) gives the second
branch. If h=c+dz, c,d nonzero, then C=c and D=h^3/c^2. Thus v
has the stated form. The second identity in (2) requires h^3 to
divide b^2, since h and z are coprime. The linear polynomial h must
therefore divide b at least twice. Write b=h^2 ell; the formula for k
follows. This is the third branch.

If h(0)=0, the unique root is r=0. Write eta=h'(0), and denote the
values of v,b,k at zero by v0,b0,k0. Applying the same rational fact
directly in t to the square of (1) yields

\[
 (1-\eta t)^2(1-2v_0t)=1,\qquad
 tk_0+\frac{t^2b_0^2}{2(1-2v_0t)}=0.
\]

The polynomial product identity gives eta=v0=0. The remaining identity
gives k0=b0=0. These conditions are exactly the first branch.
Conversely, substitution of each of the three branches into (1) gives
Phi_1=1, with the normalized square root; hence the list is complete.

## Safety and cutoffs in the first two branches

Give Z weight 1, W weight -1 and Y weight 0. Every monomial of P in
the first branch has weight at least 1. A nonzero Gaussian contraction
has weight zero. Each monomial of Q has weight at least -deg Q.
Thus E[QP^m]=0 for m>deg Q.

For the second branch, put

\[
 \mathcal D=Z\partial_Y-Y\partial_W,\qquad
 f(Z)=-b(Z)/c,\qquad T=\exp(f\mathcal D).
\]

The derivation f D is locally nilpotent: on each monomial each application
lowers 2 deg_W+deg_Y by one. Thus T and its inverse are polynomial
automorphisms. Explicitly,

\[
 TZ=Z,\qquad TY=Y-Zb/c,\qquad
 TW=W+bY/c-Zb^2/(2c^2),\qquad P=cTW.                  \tag{3}
\]

This automorphism preserves the Gaussian functional algebraically.
Gaussian integration by parts gives, for every polynomial S,

\[
 E[f\mathcal D S]
 =E[\partial_W(f\partial_Y S)-\partial_Y(f\partial_W S)]=0,
\]

because f depends only on Z. The finite exponential therefore satisfies
E[TS]=E[S]. Applying this coefficientwise to the exponential and using
circular translation gives the all-order polynomial identity

\[
 \Phi_Q(t)=E[(T^{-1}Q)(Z+ct,W,Y)].                    \tag{4}
\]

Here

\[
 T^{-1}Y=Y+Zb/c,\qquad
 T^{-1}W=W-bY/c-Zb^2/(2c^2).
\]

For Q=Z^p W^q Y^s, the Z degree of T^{-1}Q is at most
p+(2d+1)q+(d+1)s. This proves the stated bound, and the finer
monomial bound. For b nonzero of degree d with leading coefficient
beta, Q=W^q attains degree (2d+1)q: its leading coefficient in (4) is
(-beta^2/(2c^2))^q c^{(2d+1)q}. For b=0, Q=Z^q attains degree q.

## Unsafe branch and admissible fibre directions

For h=c+dz the root is r=ct/(1-dt). Since Phi_Z=r Phi_1, the fixed
marker Z proves the third branch unsafe at every positive order.

Let e=alpha partial_Z+beta partial_W+gamma partial_Y. If D_e^2P is
constant, its coefficient of Y^2 must vanish. That coefficient is
alpha^2 v''=-alpha^2 d^3/c^2, so alpha=0. With alpha=0, the second
directional derivative is 2 gamma^2 v(Z). Since v is nonconstant,
gamma=0. The only possible vectors are multiples of partial_W.
The Gaussian bilinear norm of this direction is zero. An ordinary
Gaussian fibre direction has nonzero norm. A complex orthogonal
transformation preserves this distinction, and a translation does
not change whether a polynomial second derivative is constant.
Thus the entire unsafe branch is excluded from H(2,1), in any affine
Gaussian frame. Long's quartic is c=d=1, ell=0 in this branch.

## Constant mixed source: purity forces a quadratic phase

There is a further all-degree H(2,1) subfamily, with genuine dependence
of the source on both circular base coordinates:

\[
 P=W(h(Z)+\lambda Y)+aY^2+u(Z)Y+k(Z),\qquad\lambda\ne0,
                                                               \tag{5}
\]

where h,u,k have arbitrary degrees. Every pure member of (5) is
quadratic, and hence safe. Combined with the constant-v slice above,
this settles H(2,1) for

\[
 B(Z,W)=\lambda W+u(Z),\qquad g(Z,W)=Wh(Z)+k(Z),
 \qquad\lambda\in\mathbb C.
\]

This extension depends on the **Leading-coefficient obstruction**
(label cor:leadingobstruction) in
[the one-base manuscript](../../../unpublished/gaussian-quadratic-fibre-jordan-cancellation/main.tex).
Its exact statement is: for a one-variable polynomial phase phi(t,z) of
degree D>=2 with leading coefficient a_D in C(t), and nonzero rational
normalization H, the strict inequality v_infinity(H^D/a_D)>0 excludes
a rational horizontal quotient of the differential seed submodule
taking the seed to one. That manuscript supplies a written formal
stationary-ball proof. Independent review of that imported argument
remains pending; the present application does not strengthen its
assurance. The earlier three-branch classification does not use it.

Here is a coefficientwise bridge retaining the original seed. Put

\[
 s(t,z)=\frac{z/t-h(z)}\lambda,\quad H=(\lambda t)^{-1},\quad
 \varphi=-s^2/2+t(as^2+u(z)s+k(z)).
\]

Let r=t(h(r)+lambda Y), J=1-t h'(r), and
F=t(aY^2+u(r)Y+k(r)). The formal circular identity, now over C[Y], gives

\[
 L(q)=E_Y[q(t,r)e^F/J],\qquad L(1)=\Phi_1.
\]

All coefficients are polynomials in Y, with finitely many negative t
powers allowed when q is in C(t)[z]. Since r_Y=lambda t/J and
F=phi(t,r)+Y^2/2, Gaussian integration by parts gives

\[
 L((\partial_z+\varphi_z)q)
 =H E_Y[(\partial_Y-Y)(q(t,r)e^F)]=0.
\]

The time connection is also the required one, including normalization:

\[
 \partial_t L(q)-L(q_t+(\varphi_t+H'/H)q)
 =H E_Y[(\partial_Y-Y)(r_t q(t,r)e^F)]=0.              \tag{6}
\]

To check (6), use 1/J=H r_Y and (r_Y)_t=(r_t)_Y when differentiating
L. Thus L is an actual horizontal formal functional on the spatial
twisted quotient, with its unit seed unchanged. If Phi_1=1, its values
on successive connection derivatives of the seed are 1,0,0,...;
restriction to their C(t)-span supplies exactly the rational quotient
forbidden by the imported criterion.

Expand the phase:

\[
 \varphi=-\frac{z^2}{2\lambda^2t^2}
 +\frac{zh+az^2}{\lambda^2t}
 -\frac{h^2}{2\lambda^2}-\frac{2azh}{\lambda^2}
 +\frac{zu}{\lambda}
 +t\left(\frac{ah^2}{\lambda^2}-\frac{uh}{\lambda}+k\right).
                                                               \tag{7}
\]

If its actual z degree D is at least three, its leading coefficient
is a nonzero linear combination of t^{-1},1,t; therefore its valuation
at infinity is at most one. But v_infinity(H)=1, so
v_infinity(H^D/a_D)>=D-1>0. This excludes purity. Consequently
deg_z phi<=2. Comparing successively the coefficients of t^{-1},1,t
in (7) gives deg h<=1, deg u<=1, deg k<=2. Hence P is quadratic.

For completeness, quadratic purity implies safety by a finite matrix
argument. Write P=X^T A X+l^T X+k in ordinary Gaussian coordinates,
with symmetric A. Its period is

\[
 \det(I-2tA)^{-1/2}
 \exp\!\left(tk+\frac{t^2}{2}l^T(I-2tA)^{-1}l\right).
\]

Purity and the rational-exponential fact force det(I-2tA)=1 and the
displayed exponent to vanish. Thus A is nilpotent, the inverse precision
is polynomial in t, and its mean t(I-2tA)^{-1}l is polynomial as well.
Every marked period is then a finite Gaussian Wick polynomial in this
mean and covariance, proving safety. This conclusion uses no bound on
the original degrees in (5).

## A nonconstant mixed source with zero root

One further elementary case allows a nonconstant coefficient of W in B.
Consider

\[
 P=W(h(Z)+c(Z)Y)+aY^2+u(Z)Y+k(Z),\qquad h(0)=c(0)=0.
\]

This phase is pure if and only if

\[
 a=0,\quad h,c\in z^2\mathbb C[z],\quad u,k\in z\mathbb C[z].
                                                               \tag{8}
\]

Every such phase is safe, with deg_t Phi_Q<=deg Q by positive weight.
To prove necessity, the Lagrange root is identically zero. Write
eta=h'(0), gamma=c'(0), u0=u(0), k0=k(0). Then

\[
 \Phi_1=E_Y\!\left[\frac{e^{t(aY^2+u_0Y+k_0)}}{1-t(\eta+\gamma Y)}\right].
\]

Set delta=1-2at, mu=t u0/delta,
D=1-t eta-t gamma mu, and q=t gamma/(sqrt(delta)D). Completing the
Gaussian square coefficientwise gives

\[
 \Phi_1=\frac{e^{tk_0+t^2u_0^2/(2\delta)}}{\sqrt\delta\,D}
             M(q^2),\qquad M(s)=\sum_{n\ge0}(2n-1)!!s^n,       \tag{9}
\]

where (-1)!!=1. The tilt identity follows first for an exponential
linear marker and hence by differentiation for every polynomial marker;
the geometric denominator contributes only finitely many terms to each
t coefficient, so (9) is a valid formal identity.

If gamma is nonzero, q is a convergent germ with q(0)=0 and
q'(0)=gamma. It has a convergent inverse t=t(q). Purity in (9) would
therefore make M(q^2) a convergent power series in q. Its even
coefficients (2n-1)!! have unbounded nth root, which contradicts a
positive convergence radius. This uses convergence only for the
explicit rational and exponential scalar germs, not for the original
Gaussian integral. Thus gamma=0. The remaining identity is a rational
exponential times (1-t eta)^{-1}(1-2at)^{-1/2}; the same rational
argument used above forces eta=a=u0=k0=0. This proves (8).

## Nonconstant squarefree mixed source: no pure phases

If c is any nonconstant squarefree polynomial in the preceding mixed
family, there are no pure phases, regardless of h,u,k,a. No upper
degree bound is imposed on c. Squarefree means that every complex root
of c is simple. Repeated roots are explicitly excluded from this result.

The change from one puncture to several requires a logarithmic bound
for a rational spatial phase. We give the extension, including the
global completeness and tail estimates, before applying it.

### Fixed-puncture logarithmic seed bound

Let S be a fixed finite subset of P1, let U=P1 minus S, and let
phi(t,z) be in C(t)(z), with all spatial poles in S. Let omega be a
fixed meromorphic differential on P1 with at most simple poles in S.
For the full exponential de Rham system on U with distinguished seed
[omega], the zero-exponential part of every scalar flat seed period
is O(1+log|t|) on a generic ray at infinity. A nonzero rational scalar
normalization H with v_infinity(H)>0 therefore excludes a rational
horizontal quotient of its differential seed submodule taking the
seed to one. The exclusion also applies to a unit horizontal formal
functional given at another parameter point.

Here "zero-exponential" refers to the summand with zero negative-power
exponential part in the formal irregular decomposition at infinity;
it does not mean that the entire connection is regular. The estimate
does not assert boundedness of arbitrary nonzero exponential factors.

This extends **Logarithmic seed bound** (label thm:laurentlog) in
[the boundary manuscript's Laurent proof](../../../unpublished/gaussian-boundary/laurent_logarithmic_bound.tex).
It uses the same three imported comparison facts:

1. Moderate nearby cycles recover the formal regular summand, and
   commute with proper direct image, retaining monodromy.
2. In the normal-crossing exponential model, a component with a
   vertical phase pole contributes no moderate nearby cycles. A
   finite-phase/pole node leaves the localized exponential connection
   on the finite component; a node between finite components retains
   its annulus complex and its gluing maps.
3. The exponential de Rham pairing with rapid-decay cycles is perfect,
   and horizontal functionals admit period representations away from
   finitely many exceptional parameter values.

The exact source statements are Roucairol, Section 1.2 and Lemma 1.4
([primary text](https://arxiv.org/html/math/0604134v2)); Sabbah,
Section 1.b, especially Proposition 1.5
([primary text](https://arxiv.org/html/0706.3570v2)); and Hien--Roucairol,
Theorems 2.7 and 3.5
([published primary text](https://www.numdam.org/item/10.24033/bsmf.2564.pdf)).
These are comparison statements, not logarithmic growth bounds. Their
retained source audit is the boundary manuscript's
[primary-source manifest](../../../unpublished/gaussian-boundary/sources/laurent_primary_manifest.json).
The proof below supplies the extra quantitative estimate. Independent
review of the imported and extended arguments remains pending.

**A complete marked model.** The case omega=0 is immediate. Otherwise
choose a projective coordinate
with infinity in S; this is possible since a nonzero logarithmic
differential on P1 has a pole. All the hypotheses are invariant under
that fixed coordinate change. In this lemma r denotes 1/t, rather
than the Lagrange root used in the Gaussian application. If phi is
nonconstant in z,
mark all points of S and all branches of the fibres phi=0,1,infinity.
Their convergent Puiseux expansions split after a finite ramification
of r. Separate the finitely many marked sections by blowing up their
collisions. Contact orders are finite and decrease, so the construction
terminates. A section cannot meet a node uv=r: the positive orders of
u and v would have to sum to one. Thus every blow-up centre is smooth
on a single central component; the central fibre stays reduced, with
node equations uv=r. The separated zero, pole and value-one sections
remove the indeterminacies of phi, exactly as in the Laurent proof:
a remaining zero/pole crossing would carry a value-one branch through
the same point. Only the original S sections are punctures of the
generic exponential connection; the auxiliary markings remove no points.

Every rational function z-alpha, alpha finite in S, has its horizontal
divisor in those marked sections. Locally on the resulting surface it
is u^a v^b V, with a,b integers and V a holomorphic unit. Since P1 has
no nonzero holomorphic differentials, partial fractions express omega
as a fixed complex linear combination of dlog(z-alpha). Consequently
its pullback is a fixed sum of logarithmic monomials and unit terms.
In particular no parameter-dependent residue introduces a power of r.

Apply comparison facts 1 and 2 to this proper model, localized at the
central fibre and all S sections. The regular nearby-cycle complex is
computed by finite-phase components localized at their irregular ends
and original punctures, together with the actual annulus complexes
at finite/finite nodes. Remove small discs at the nodes and punctures
and triangulate the remaining compact pieces. Attach descending phase
tails at irregular ends and circles at regular punctures; glue with
the full annulus complex at each surviving node. This computes the
nearby-cycle maps and monodromy, rather than just the Euler characteristic.
Proper comparison therefore supplies a complete finite basis of the
regular part from these chains. Its coefficients are fixed constants
on a small parameter sector. The same construction includes every
point of S; there is no restriction to two punctures or to simple
critical points. The de Rham system is finite dimensional: reductions
by the twisted derivative bound pole orders at each of the finitely
many regular or irregular punctures, and bound the polynomial part at
infinity. No particular rank formula is needed here.

**Uniform estimates.** On each compact finite-phase piece the phase
and differential coefficients are bounded. At a finite/finite node,
fixed r gives v=r/u and hence

\[
 d\log(z-\alpha)
 =\left(a-b+\frac{uV_u-vV_v}{V}\right)\frac{du}{u}.
\]

The coefficient is uniformly bounded. A bridge across
|r|/epsilon<=|u|<=epsilon has logarithmic length O(1+|log r|),
and circles have bounded logarithmic length. At a finite/pole node,
a local phase coordinate gives phi=w^{-p}, p>=1, and r=wv up to
a holomorphic unit. The logarithmic differential is bounded by a
constant times |dw|/|w|. Along a descending tail its exponential integral
is bounded by a constant times integral_{s0}^infinity e^{-s} ds/s.

The continuation outside a node chart also needs a uniform estimate.
Choose disjoint fixed discs about the finite points of S, an infinity
disc with coordinate x=1/z, and a compact complement. In each puncture
disc use x=z-alpha, or x=1/z at infinity. A common M contains all
transition scales in |r|^M<=|x|<=epsilon, since phi has a fixed rational
spatial denominator and finitely many rational parameter coefficients.
On these annuli a descending inverse-phase path lies on Im(phi)=0,
a real algebraic curve of uniformly bounded degree after clearing
the rational denominator. The
[curve-length estimate](../../../unpublished/gaussian-boundary/analytic_core.tex),
label lem:curvelength, bounds integral |dx|/|x| by O(1+|log r|)
using dyadic annuli. There are only finitely many inverse branches;
strict monotonicity of the real phase prevents retracing with unbounded
multiplicity. The compact complement has the same bounded-degree
length bound, with bounded differential coefficients.

Beyond the transition annuli the leading actual polar monomial
dominates uniformly. Since omega has only a logarithmic pole there,
|omega/dphi|<=C/|phi| on a descending inverse branch. Terms in omega
regular at that puncture satisfy an even smaller bound. Truncate a
local node tail at |w|=|r|^theta, 0<theta<1, before continuing it.
Its phase is already <=-|r|^{-theta p}; all subsequent contributions
are bounded by this exponential times the logarithmic annular bound,
or by the corresponding exponential tail integral. Choose the ray to
avoid unbounded critical values and parameter Stokes directions; the
finitely many critical values are convergent Puiseux functions, so
this excludes only finitely many limiting directions. This justifies
continuation to an irregular point of S. Regular puncture circles have
bounded phase and length. These estimates also allow a fixed rational
amplitude, with at most an additional algebraic power of |r|^{-1}.

Thus the complete regular-period basis has the asserted seed estimate.
The formal exponential decomposition separates its zero-exponential
part from nonzero factors. Comparison facts 1--3, or the same regular
part argument in the Laurent theorem, show that every regular flat
seed functional is a constant linear combination of this complete basis.
After multiplication by H all those seed values tend to zero.

To exclude a unit seed quotient, use the
[seed-functional lemma](../../../unpublished/gaussian-boundary/seed_period_comparison.tex),
label lem:seedfunctional. A unit formal functional first restricts
to a rational quotient on the differential span of the seed. At an
ordinary point extend its fibre value to the full dual connection;
the analytic horizontal extension has seed value one identically.
The total space (A1 minus finitely many parameter values) times U is
smooth affine; phi is regular there, and the projection and phi are
algebraically independent when phi is nonconstant in z. Thus comparison
fact 3 applies. Formal exponential separation would make the
zero-exponential seed value of this extension equal to one, contradicting
its decay after normalization. If phi is independent of z, all periods
are constant residue periods times exp(phi(t)); the rational-exponential
fact excludes H exp(phi(t)) being a nonzero constant when H has positive
infinity valuation. This proves the fixed-puncture lemma, including
the constant spatial phase.

### Application to the mixed Gaussian source

Put s=(z/t-h(z))/c(z) and phi=-s^2/2+t(as^2+u(z)s+k(z)). Its
spatial poles lie in S=zeros(c) union {infinity}. Squarefreeness gives

\[
 \frac{dz}{c(z)}
 =\sum_{c(\alpha)=0}\frac{1}{c'(\alpha)}d\log(z-\alpha).
\]

The normalized volume is (1/t) dz/c(z); its scalar normalization has
infinity valuation one. The fixed-puncture lemma excludes purity once
the original functional extends to all required localized amplitudes.

The actual functional extends to the required Laurent amplitudes.
If c(0) is nonzero, c(r) is a t-adic unit. If c(0)=0 and h(0) is
nonzero, then r=t h(0)+O(t^2), so c(r)=c'(0) h(0)t+O(t^2) is
invertible as a formal Laurent series with polynomial Y coefficients.
These cover all cases except the zero-root case already handled.
In that remaining case squarefreeness forces c'(0) nonzero, contradicting
(8). Localization is not asserted when c(r)=0.

For clarity, the spatial and time bridges for general c are as follows.
Set J=1-t(h'(r)+c'(r)Y), F=t(aY^2+u(r)Y+k(r)), and
L(q)=E_Y[e^F q(r)/J]. Since r_Y=t c(r)/J,

\[
 L(c(q_z+\varphi_zq))=\frac1t E_Y[(\partial_Y-Y)(e^Fq(r))]=0.
\]

This is the full spatial relation for volume dz/c, not a restriction to
one selected partial-fraction component. With A(t,z)=1/(t c(z)),
differentiation gives

\[
 \partial_tL(q)-L(q_t+(\varphi_t-1/t)q)
 =E_Y[(\partial_Y-Y)(e^F A(t,r)r_tq(r))]=0.
\]

This follows from 1/J=A(t,r)r_Y and equality of mixed derivatives,
as in (6). It preserves the original seed, so the fixed-puncture
logarithmic obstruction applies. Nonconstant affine c is the degree-one
special case. The exclusion of squarefree c does not extend across
repeated roots: for example c=z^2 in (8) permits positive-weight pure
phases.

## A double root: the remaining quadratic mixed coefficient

Suppose c(z)=kappa(z-alpha)^2, with kappa nonzero. Outside the case
alpha=h(0)=0 there are no pure phases in the mixed family, with
unrestricted degrees in h,u,k. In that exceptional case (8) applies.
Together with the preceding results, this classifies all pure members
when deg c<=2. The only pure nonconstant-c branch in that range is

\[
 c=\kappa z^2,\qquad a=0,\qquad
 h\in z^2\mathbb C[z],\quad u,k\in z\mathbb C[z],
\]

which is safe with cutoff m>deg Q. Thus H(2,1) holds for every
P=W(h(Z)+c(Z)Y)+aY^2+u(Z)Y+k(Z) with deg c<=2, even when the
other polynomial degrees are arbitrary. In particular it holds for
the entire total-degree-at-most-four stratum affine in circular W.

The proof again uses the complete marked-model comparison, but now
the differential has a higher-order pole before changing coordinates.
The logarithmic lemma alone would not justify this conclusion.

Put x=1/(z-alpha) and r=1/t. The transformed volume is
-(r/kappa) dx. The spatial phase is a Laurent polynomial in x. Its
positive polynomial part has degree exactly four, with leading coefficient

\[
 a_4(r)=-\frac{(1-2a/r)(\alpha r-h(\alpha))^2}{2\kappa^2}.
                                                               \tag{12}
\]

Outside alpha=h(alpha)=0 this coefficient is nonzero and has
r-valuation at most two. Every coefficient of a negative Laurent
power has r-valuation at least minus one. To check the latter bound,
expand the phase as

\[
 -\frac{z^2r^2}{2c^2}
 +r\left(\frac{zh+az^2}{c^2}\right)
 -\frac{h^2}{2c^2}-\frac{2azh}{c^2}+\frac{zu}{c}
 +r^{-1}\left(\frac{ah^2}{c^2}-\frac{uh}{c}+k\right),
\]

then substitute z=alpha+1/x. All coefficients are linear combinations
of r^{-1},1,r,r^2. Terms x^4 in (12) are unaffected by u and k,
because these contribute at most x^2 at this puncture. The condition
alpha=h(alpha)=0 is equivalent to alpha=h(0)=0 and is exactly the
exception treated by (8).

**Generic component charts.** In the marked P1 model used above, each
vertical component has, after a common ramification, an exact affine
coordinate

\[
 X=\frac{x-\xi(r)}{r^\rho}
\]

whose generic residue on that component is a degree-one rational
coordinate. Here xi is a finite Puiseux Laurent polynomial and rho
is rational. Consequently the relative differential dx has valuation
rho there. Normalize v(r)=1 throughout: if the model parameter is R
with r=R^N, actual component orders are N times these valuations.

For completeness this chart assertion follows by induction on the
smooth-point blowups, after making the common ramification r=R^N.
The original component has x=X, with xi=0 and integer scale m=0.
Suppose x=xi(R)+R^m X and the residue of X has degree one. At a
smooth finite residue point X=lambda, the pair (R,X-lambda) is a
regular coordinate system. The new exceptional chart X-lambda=R U
gives the exact expression x=(xi+lambda R^m)+R^{m+1}U, with
degree-one residue U. At a smooth infinity point of X use V=1/X;
the chart V=R U gives x=xi+R^{m-1}/U, with degree-one residue 1/U.
The old component keeps its generic chart. These operations preserve
a finite Laurent centre and an integer scale m. Marked sections do
not meet nodes, so these cases exhaust the centres in the marked
construction. Set rho=m/N. The definition of X is exact; no omitted
higher-order coordinate terms can change a top coefficient.
Characteristic zero ensures the residue derivative is nonzero.

Consider any finite-phase component. Its generic phi has nonnegative
valuation. If rho<=-1, then
v(x)=min(v(xi),rho)<=-1, since the residue X is nonconstant.
Every negative Laurent term b_q x^{-q}, q>=1, has valuation at least
-1-q v(x)>=q-1>=0. The entire negative part is therefore integral
there. On the other hand the polynomial part after x=xi+r^rho X has
exact top X^4 coefficient a4 r^{4rho}, of valuation at most
2+4rho<0. Its coefficients cannot cancel across distinct powers of
the nonconstant residue coordinate X. It is nonintegral, and cannot
be cancelled by the integral negative part. This contradicts finite
phase. Thus rho>-1 on every finite-phase component, and r dx has
strictly positive generic valuation on every such component.

**From generic valuations to a complete seed bound.** A generic check
alone would omit the necks and irregular ends. On a node chart R=uv,
relative differentiation at fixed R gives

\[
 r\,dx=f(u,v)\frac{du}{u},\qquad
 f=r(u x_u-v x_v),\qquad r=R^N.
\]

The meromorphic coefficient f has poles only on vertical divisors and
the marked x=infinity section. Since a4 is nonzero, this horizontal
section is an irregular phase pole; x=0 is not a pole of dx. No
marked section meets a node. At a
finite/finite node the positive generic orders on both components
therefore give f=u^A v^B V with A,B positive integers and V bounded
holomorphic. The relative integral on a neck is bounded by a constant
times |R|^{min(A,B)}(1+|log R|), which tends to zero. This also shows
why a large moving centre in x does not create a missing volume term:
the displayed relative derivative differentiates it to zero.

At a finite/pole node, the order A on the finite component is positive;
the other order B may be negative. Substitution u=R/v gives
f=R^A v^{B-A}V. Since the resolved phase is a morphism, it is a
pole monomial v^{-p} times a unit there, without an opposite zero
factor. Absorb the unit into the phase coordinate. The descending
integral is O(|R|^A), because e^{-s} integrates any fixed algebraic
power of s. If the node chart ends, truncate at a small power of |R|
and continue the descending inverse branch. The remaining exponential
decay dominates every algebraic growth of this rational differential;
the uniform transition-annulus argument above supplies that bound.
On compact finite pieces the positive component order is explicit;
regular puncture circles and any other irregular ends satisfy the
same fixed-chart estimate. There are finitely many pieces.

The complete regular-part basis supplied by proper nearby-cycle
comparison thus takes the normalized seed r dx to values tending to
zero. The same unit-seed restriction and analytic period extension
used above exclude a rational unit seed quotient. This establishes
the obstruction at (12), rather than merely checking some stationary
points. The actual Gaussian branch localizes at z-alpha when alpha
is nonzero, or when alpha=0 and h(0) is nonzero; its leading Laurent
coefficient is independent of Y in both cases. Those are exactly the
cases used in this proof. The remaining nonlocalizable branch was
handled without localization by (8).

### Several simple or double roots

The same local argument gives an all-degree classification without a
bound on deg c. Suppose c is nonconstant and every root has multiplicity
at most two. Then the mixed phase is pure if and only if

\[
 a=0,\qquad h,c\in z^2\mathbb C[z],\qquad
 u,k\in z\mathbb C[z].                              \tag{13}
\]

These are precisely the safe positive-weight phases already handled by
(8). In particular every pure nonconstant c in this range has a double
root at zero. Nonzero constant c is not included in this assertion.

It remains to exclude purity outside the zero-root case h(0)=c(0)=0.
At a double root alpha, write c=(z-alpha)^2 q(z), q(alpha) nonzero,
and use x=1/(z-alpha). The normalized volume is

\[
 -r\,\frac{dx}{q(\alpha+1/x)}.
\]

Its extra factor is a generic unit on every component reducing to
alpha. The phase now has a convergent Laurent expansion at x=infinity;
its positive part is quartic with leading coefficient (12), replacing
kappa^2 by q(alpha)^2. Outside alpha=h(alpha)=0 it is nonzero of
valuation at most two. Every coefficient of the negative Laurent
tail still has valuation at least minus one, since all spatial
denominators are independent of r and the parameter powers remain
r^{-1},1,r,r^2. When v(x)<0 this infinite tail converges in the
parameter valuation. For a component scale rho<=-1 its qth negative
term has valuation at least q-1, tending to infinity. It is therefore
integral and cannot cancel the nonintegral positive quartic top
coefficient in the exact chart. Thus rho>-1, and the normalized volume
has strictly positive order 1+rho, exactly as above.

Every other type of component also has positive normalized volume
order. At a simple root dz/c is a holomorphic unit times a logarithmic
differential, whose generic pullback order is nonnegative. At infinity
dz/c is logarithmic or regular because deg c>=1. Away from these
punctures it is regular; a component with nonconstant original-coordinate
residue has normalized order one. These cases exhaust the generic
centres of the marked P1 model. The only higher-order horizontal
amplitude poles are the double roots; each is a phase-irregular
puncture by its nonzero quartic coefficient. The positive-order
finite/finite and finite/pole neck estimates therefore apply to the
full normalized differential r dz/c, without a hidden higher-order
regular-puncture pole. Complete regular-part comparison and the
unit-seed argument exclude purity.

Localization is valid in every excluded case. If c(0) is nonzero it
is an ordinary t-adic unit; otherwise h(0) is nonzero and, for a root
of multiplicity m=1 or 2 at zero, c(r) has a nonzero leading term
proportional to h(0)^m t^m, independent of Y. The one excluded local
exception alpha=h(alpha)=0 is exactly the zero-root case already
classified by (8). This proves (13).

## Arbitrary multiplicities: nonzero quadratic fibre and unshared roots

Two further obstructions use the Gaussian-square structure of the phase
rather than a multiplicity bound. In the mixed family with nonconstant c:

1. If a is nonzero, there are no pure phases, for any root multiplicities
   of c and any h,u,k.
2. If every common root of c and h has multiplicity at most two in c,
   purity is equivalent to (8). In particular gcd(c,h)=1 excludes purity
   for every nonconstant c, without a degree or multiplicity bound.

Here h=0 has every root of c as a common root; this convention prevents
the second assertion from imposing an unintended coprimality condition.
The positive-weight zero-root form in (8) remains sufficient for arbitrary c.

Write r=1/t and s=(rz-h)/c. The transformed phase is

\[
 \varphi=\frac{as^2+us+k}{r}-\frac{s^2}{2}.
\]

Work on a finite-phase component whose generic centre reduces to a root
alpha of c, of multiplicity m. Let N=v(z-alpha)>0 and let rho be the
relative scale in an exact z chart, so rho>=N. All of u,k are regular
at alpha. If a is nonzero and v(s)<0, the term as^2/r has strictly
smaller valuation than every other term. Finite phase therefore forces
s integral. If h(alpha) is nonzero this is impossible, since
v(s)=-mN. At alpha=0 outside the zero-root case, h(0) is nonzero,
so no finite-phase component lies there either.

For alpha nonzero with h(alpha)=0, put ell=ord_alpha h, allowing
ell=infinity when h=0. If ell>=m then h/c is regular. Integrality of
s=rz/c-h/c gives mN<=1, and hence

\[
 v\!\left(r\frac{dz}{c}\right)
 =1+\rho-mN\ge1-(m-1)N\ge\frac1m>0.
\]

If 1<=ell<m, the identity h+cs=rz forces ell N=1: h has strictly
smaller valuation than cs, whereas rz has valuation one. In
D=r-h'-c's, the term -h' is then the unique smallest term, with
valuation (ell-1)N=1-N. The exact relative identity ds=(D/c)dz gives

\[
 v\!\left(r\frac{dz}{c}\right)
 =1+v(ds)-v(D)\ge N>0.
\]

Here an integral rational function s has relative differential of
nonnegative generic valuation in the exact component coordinate;
a constant leading residue only increases that valuation. These cases
cover every root for a nonzero a. Away from the roots and at infinity
the normalized volume has positive order as already proved. Every
horizontal volume pole is phase-irregular: its leading coefficient is
proportional to (a/r-1/2)(r alpha-h(alpha))^2, which is nonzero outside
alpha=h(alpha)=0. The complete positive-order neck and tail argument
therefore excludes a unit seed. The exceptional zero-root case itself
forces a=0 by (8), completing assertion 1.

For assertion 2 only unshared roots of higher multiplicity remain.
For a=0, finite phase forces v(s)>=-1: below -1 the term -s^2/2
uniquely dominates us/r and k/r. At an unshared root h(alpha) is
nonzero, so v(s)=-mN and again N<=1/m. The same volume inequality
gives a positive order at that root. Shared roots of multiplicity one
or two have already been treated. Complete comparison excludes purity
outside the zero-root case, and (8) supplies its exact classification.

## Common roots at which h vanishes at least as strongly as c

Let

    P=W(h(Z)+c(Z)Y)+aY^2+u(Z)Y+k(Z),

where a is constant, c is nonconstant, and h,c,u,k are complex
polynomials. Suppose every common root alpha of c,h satisfies

    ord_alpha h >= ord_alpha c,                         (H5)

with ord_alpha 0=infinity. Then P is pure if and only if

    a=0,    h,c in z^2 C[z],    u,k in z C[z].            (H6)

Every pure member is safe, with cutoff m>deg Q. There is no degree or
root-multiplicity bound. In particular this classifies the entire
structural subfamily c|h, including h=0. For that subfamily the condition
on h in (H6) follows from the condition on c.

The nonzero-a exclusion and the zero-root classification (8) are already
proved above. Thus it remains to exclude purity when a=0 and
(h(0),c(0))!=(0,0). Put r=1/t and

    s=(rz-h)/c,
    phi=-s^2/2+(us+k)/r,    eta=r dz/c.

### The original source and exact localization

The source is the actual mixed Gaussian functional, not an arbitrarily
chosen local solution. Denote its Lagrange root by zeta:

    zeta=t(h(zeta)+c(zeta)Y),    zeta(0,Y)=0,
    J=1-t(h'(zeta)+c'(zeta)Y),
    F=t[u(zeta)Y+k(zeta)],
    L_t(q)=E_Y[e^F q(zeta)/J].

Then L_t(1)=E[e^{tP}]. If c(0)!=0, c(zeta) is an ordinary t-adic unit.
If c(0)=0 and h(0)!=0, write c(z)=C z^m+O(z^(m+1)), C!=0. Then

    zeta=t h(0)+O(t^2),
    c(zeta)=C h(0)^m t^m+O(t^(m+1)).

Its leading coefficient is nonzero and independent of Y. Hence in both
cases every required amplitude localized at c has a coefficientwise
formal Laurent evaluation with polynomial Y coefficients. The retained
integration-by-parts identities give

    L_t(c(q_z+phi_z q))=0,
    partial_t L_t(q)=L_t(q_t+(phi_t-1/t)q).

These are the full spatial and parameter connections for the normalized
seed eta, and preserve its original value. Purity would give a horizontal
unit seed functional. No localization is asserted in the remaining case
c(0)=h(0)=0; that case is handled directly by (8).

The phase is spatially nonconstant. If c is not proportional to z this
follows from its r^2 coefficient -z^2/(2c^2). If c=kappa z, the present
case forces h(0)!=0, and its r coefficient zh/c^2 has a pole at zero.
On the fixed punctured affine curve obtained by deleting the roots of c
and infinity, phi is regular, and the base projection and phi are
algebraically independent. Thus the same ordinary-parameter exponential
period comparison and seed-extension lemma used above apply.

### Exact scales and all possible leading faces

Use the complete marked P1 model from the preceding rational-phase
argument, including all finite-phase components and annulus gluing maps.
Normalize v(r)=1 after finite ramification. At an unshared root alpha of
c of multiplicity m, the already proved square argument gives v(s)>=-1
on a finite-phase component. Since h(alpha)!=0, v(s)=-mN, where
N=v(z-alpha)>0. Consequently N<=1/m and

    v(eta)>=1-(m-1)N>=1/m>0.

It remains to treat a common root. It is nonzero in the present case.
Write m=ord_alpha c and h=c v locally, with v regular, as permitted by
(H5). With x=(z-alpha)^(-1), the phase has convergent Laurent groups

    phi=r^2 F2(x)+r F1(x)+F0(x)+r^(-1) Fm1(x),
    F2=-z^2/(2c^2),       F1=z v/c,
    F0=-v^2/2+zu/c,       Fm1=k-u v.

Their actual highest exponents, omitting an identically zero group, obey

    j2=2m,    j1<=m,    j0=A<=m,    jm1=B<=0.

The leading coefficient of F2 is nonzero because alpha!=0. These bounds
allow every vanishing or cancellation in the coefficient functions.

On a component centred at alpha put v(x)=-N<0. Its exact disc chart is

    x=xi(r)+r^rho X,    rho>=-N,
    epsilon=rho+N>=0.

Here xi is a Puiseux centre and X has nonconstant generic residue. Such
charts follow the finite sequence of point blowups defining a component
of the marked P1 model: absorb each constant leading residue into the
centre, and clear fractional exponents by finite parameter ramification.
Relative differentiation at fixed parameter gives

    v(eta)=1-(m-1)N+epsilon.                            (H1)

The moving centre contributes no differential. Positivity is immediate
if 1-(m-1)N>0, including every m=1 component. Otherwise m>=2 and

    N>=1/(m-1).                                        (H2)

In this range the level-one group lies strictly above level two:

    (1-j1 N)-(2-2mN)>=mN-1>0.

It therefore cannot enter a dangerous leading face. The remaining group
weights are

    L2=2-2mN,    L0=-AN,    Lm1=-1-BN.

Let L be their minimum. It is negative because L2<0. With tau=r^N x the
leading phase is r^L f(tau), where f has at most three monomials, one
from each level attaining L. Their exponents are distinct: equality of
exponents and weights would force equality of the levels. The convergent
negative Laurent tails have weights tending to infinity and do not add
leading monomials.

Finite phase forces the leading residue of tau to be a constant nonzero
root tau0 of f. Let its multiplicity be mu. A k-term Laurent polynomial
has nonzero-root multiplicity at most k-1: applying powers of
tau partial_tau gives the invertible Vandermonde matrix in its k
distinct exponents. Thus mu<=2, and mu=1 if only two levels occur. A
one-monomial face has no nonzero root and supplies no finite-phase
component with L<0.

The final exact chart necessarily satisfies

    epsilon>=-L/mu.                                    (H3)

Indeed its X^mu coefficient in phi has valuation exactly L+mu epsilon.
The mu-th derivative of f at the full Puiseux centre is a unit. Every
higher-weight group remains higher-weight after these derivatives at a
unit tau centre; negative Laurent powers have unit denominators there.
Positive-order shifts of the centre do not change this conclusion.
Consequently lower-X-power cancellations and further critical-point
refinements cannot cancel this coefficient. Finite generic phase
requires its valuation to be nonnegative, proving (H3).

Combining (H1) and (H3), it suffices to check

    E=1-(m-1)N-L/mu>0.                                 (H4)

There are precisely three cases.

* If level two occurs, L=2-2mN. For mu=1,
  E=-1+(m+1)N>=2/(m-1)>0. For mu=2, E=N>0.
* If the highest level is zero, only levels zero and minus one occur,
  so mu=1. Here A>0 and comparison with level two gives
  N<=2/(2m-A). Together with (H2) this forces A>=2. Now
  E=1+(A-m+1)N. It is at least one if A>=m-1; otherwise
  E>=(A+2)/(2m-A)>0. The denominator is positive since A<=m.
* If the highest level is minus one, the face is a single monomial and
  has no nonzero root, so this case cannot supply a finite-phase component.

These alternatives include all omitted groups, coefficient cancellations,
and successive exact scales. Therefore every finite-phase component
centred at a common root has positive normalized seed order, without a
multiplicity bound.

### All centres, punctures, and crossover estimates

At any ordinary finite point the normalized volume is r times a regular
differential and has positive generic order. A component with nonconstant
finite z residue has order at least one. At infinity, if n=deg c>=1 and
v(z)=-N<0, then

    v(eta)>=1+(n-1)N>=1.

These cases exhaust all generic centres. At a c-root alpha, the highest
spatial polar coefficient of phi is proportional to
(r alpha-h(alpha))^2 and is nonzero over C(r) outside
alpha=h(alpha)=0. Thus every higher-order horizontal seed pole is a
phase-irregular puncture. At infinity dz/c is logarithmic for n=1 and
holomorphic for n>=2; every phase-regular puncture has at most a
logarithmic seed pole.

Apply the complete compact-piece, neck, and descending-tail comparison
already proved above, with its stated primary-source inputs and
independent-review status. Here is why arbitrary root multiplicity causes
no missing crossover estimate. On a node R=uv write eta=f(u,v)du/u.
After extracting its two vertical divisorial orders, the remaining
coefficient is bounded holomorphic; no marked horizontal pole meets the
node. At a finite/finite node the orders A,B are positive integers, so
the whole crossover annulus contributes

    O(|R|^min(A,B)(1+|log R|)).

This includes transitions between different Newton faces and between
successive cancellation scales. At a finite/pole node the positive
finite-side factor R^A survives substitution u=R/v. The resolved phase
is a pole monomial times a unit without a competing zero factor, and
descending exponential integration absorbs every fixed extra pole order.
Although those orders depend on m, they are finite for the fixed phase;
no uniform upper bound on m is needed. Beyond node charts the fixed-degree
rational phase, differential, and algebraic inverse branches have only
algebraic parameter and phase growth. Truncation at phase size a negative
power of |R| supplies an exponential cutoff dominating that growth, as
in the preceding uniform tail argument. Regular puncture circles and
approach paths are at most logarithmic. The finite marked model supplies
only finitely many pieces.

Hence every zero-exponential normalized seed period tends to zero. A
hypothetical original unit functional transports through the actual
source bridge and, at an ordinary parameter, extends to a full analytic
curve period with seed value one. The retained perfect period comparison
and formal exponential separation contradict the complete decay bound.
This excludes purity outside the zero-root form (H6). Its converse and
safety follow from strict positive angular support, as already proved
for (8).

## A sharp common-root range, including multiplicity eight

The preceding common-root theorem can be extended using the relations
between all four parameter levels. In the same mixed family with
nonconstant c, suppose every nonzero common root alpha with

    m=ord_alpha c,    1<=ell=ord_alpha h<m

satisfies

    m-ell<=2  OR  3ell<2m.                             (SR1)

Then purity is equivalent to the positive-weight form (H6), and every
pure member is safe. Unshared roots and roots with ord h>=ord c have
no multiplicity restriction. In particular this proves the complete
classification when every nonzero common root has multiplicity at most
eight in c, without a bound on deg c,h,u,k.

More precisely, (SR1) is the sharp range for uniform strict positivity
of the normalized volume on all finite-phase components at a nonzero
common root. Its complement is exactly

    ceil(2m/3)<=ell<=m-3.                              (SR2)

The general cubic barrier below constructs a nonpositive-order chart for every
integer pair in (SR2). Those examples fail purity; this is a boundary of
the volume method, not a classification of possible pure phases there.

The nonzero-a exclusion, direct zero-root classification, actual source
localization, and complete comparison inputs are exactly those in the
preceding common-root theorem. It remains only to prove the new local
positivity assertion for a=0 and 1<=ell<m. Put p=m-ell>0 and
x=(z-alpha)^(-1). The phase has the four convergent Laurent groups

    phi=r^2 F2+r F1+F0+r^(-1) Fm1,
    F2=-z^2/(2c^2),      F1=zh/c^2,
    F0=-h^2/(2c^2)+zu/c, Fm1=k-uh/c.

Their actual highest x exponents, omitting a zero group, satisfy

    j2=2m,   j1=m+p,   j0=A<=max(2p,m),   jm1=B<=p.

The first two leading coefficients are nonzero, since alpha!=0. Also
A<=2m-2. These bounds already allow all coefficient cancellations.

Use the exact charts and normalized volume identity (H1): v(x)=-N<0,
with extra relative scale epsilon>=0, and

    v(eta)=1-(m-1)N+epsilon.

Only N>=1/(m-1) needs treatment. The group weights are

    L2=2-2mN,   L1=1-(m+p)N,   L0=-AN,   Lm1=-1-BN.

Their minimum L is negative. With tau=r^N x, the leading phase is
r^L f(tau), where f has at most four distinct monomials. The convergent
negative Laurent tails add no leading monomials. Finite phase forces
the generic tau residue to be a constant nonzero root of f, of
multiplicity mu. The same Vandermonde argument gives mu<=3, with
mu<=2 unless all four levels occur. The exact Taylor-coefficient
argument (H3) remains valid with mu=3: the mu-th derivative at a unit
tau centre is a unit, and all higher-weight Laurent terms stay higher
under those derivatives. Thus

    epsilon>=-L/mu,
    v(eta)>=E=1-(m-1)N-L/mu.                           (SR3)

First suppose mu<=2. Classify the face by its highest parameter level;
these alternatives exhaust all faces and omitted-group cases.

* Highest level 2: L=2-2mN. For mu=1,
  E=-1+(m+1)N>=2/(m-1)>0. For mu=2, E=N>0.
* Highest level 1: L=1-(m+p)N. Comparison with level two gives
  N<=1/ell. For mu=1, E=(p+1)N>0. For mu=2,
  E=1/2+(2-ell)N/2: this is at least 1/2 if ell<=2, and at least
  1/ell if ell>2.
* Highest level 0: only levels zero and minus one occur, so mu=1.
  Here A>0, and comparison with level two gives N<=2/(2m-A).
  Combining with N>=1/(m-1) forces A>=2. If A>=m-1 then
  E=1+(A-m+1)N>=1. Otherwise
  E>=(A+2)/(2m-A)>0. The denominator is positive because A<=2m-2.
* Highest level -1: a single monomial has no nonzero root, so cannot
  supply a finite-phase component with L<0.

It remains to consider mu=3. All four levels must occur, forcing

    N=1/ell,   A=2p,   B=2p-ell,   L=-2p/ell.          (SR4)

If ell<2p, then B>0. Because k is regular, the leading pole of Fm1
must come from -uh/c. Thus u is nonzero, and with q=ord_alpha u,

    p-q=2p-ell,    q=ell-p>=0.

If ell<p this is impossible. Otherwise let C,H,U be the nonzero
leading coefficients of c,h,u in z-alpha. Write T=tau^ell. The
four-level face factors exactly as

    f(tau)=tau^(2p-ell)(alpha T-H)
                 [-T(alpha T-H)/(2C^2)+U/C].           (SR5)

Indeed the four coefficients are -alpha^2/(2C^2), alpha H/C^2,
-H^2/(2C^2)+alpha U/C, and -UH/C. The bracket is quadratic and
has value U/C!=0 at the root T=H/alpha of the linear factor.
Since T=tau^ell is etale at every nonzero tau, no nonzero root
can have multiplicity three. Hence mu=3 requires ell>=2p.

In that remaining case, (SR3)--(SR4) give

    E=(3-p)/(3ell).                                   (SR6)

It is positive for p<=2. This proves (SR1), since ell<2p is equivalent
to 3ell<2m. Nonpositive order can therefore occur only for p>=3 and
ell>=2p, equivalently (SR2). In particular m=ell+p>=3p>=9 there;
the first possible pair is (m,ell)=(9,6).

This local calculation completes the global exclusion by the precise
source and period argument already proved for (H5). Unshared roots,
roots with ell>=m, ordinary centres and infinity retain their positive
orders. Every higher-order horizontal volume pole is phase-irregular
outside the directly handled zero-root case. At finite/finite nodes the
two positive orders give a positive power times logarithmic neck length;
at finite/pole nodes descending exponential integration absorbs every
fixed extra pole order. Arbitrarily large but finite multiplicity causes
no change in that integrability statement. The complete marked model,
including transitions between faces, successive refined scales and all
annulus gluing maps, is retained. Thus all zero-exponential seed periods
decay, contradicting the transported original unit seed. The converse
and safety are the already proved positive-weight implication in (H6).

## A constant quadratic W term

The following all-degree family has genuine higher W degree:

\[
 P=\lambda W^2+Wh(Z)+cWY+aY^2+u(Z)Y+k(Z),\qquad
 \lambda\ne0,
\]

where lambda,c,a are constants and h,u,k are arbitrary polynomials.
Every pure member is quadratic, and therefore safe.

In the full phase -ZW-Y^2/2+tP, eliminate W,Y over C(t). Their
quadratic precision and determinant are

\[
 M=\begin{pmatrix}-2\lambda t&-ct\\-ct&1-2at\end{pmatrix},\quad
 \Delta=-2\lambda t+(4\lambda a-c^2)t^2.
\]

The remaining phase is phi=N/(2Delta), with normalization
H=Delta^{-1/2}, where

\[
 \begin{aligned}
 N={}&z^2+t(-2zh-2az^2)\\
 &+t^2(h^2+4azh-2czu-4\lambda k)\\
 &+t^3(-2ah^2+2chu-2\lambda u^2+2(4\lambda a-c^2)k).
 \end{aligned}                                             \tag{14}
\]

This is a generic quadratic elimination of the spatial twisted module
and its constant volume seed, using
[the determinant-factor lemma](../../../unpublished/gaussian-boundary/repeated_line_exclusion.tex),
label lem:quadraticelimination. Its matrix is invertible over C(t).
Its singularity at t=0 does not justify a conditional Wick integral;
the original unit functional first supplies its rational seed quotient,
which is transported through the generic module identification. The
normalization is algebraic, so the full stationary-ball theorem is
used after finite ramification. The rational-normalization corollary
alone would not justify this application.

If the actual z degree D is at least three, its nonzero leading
numerator N_D(t) is divisible by t. For d=deg Delta in {1,2},

\[
 Dv_\infty(H)-v_\infty([z^D]\varphi)
 =(D/2-1)d+\deg N_D>0.
\]

The full leading-coefficient obstruction excludes purity. Thus phi
has degree at most two. In (14), the coefficient of t forces h affine.
The coefficient of t^2 then gives deg k<=max(2,deg u+1). If deg u=n>=2,
the term -2lambda u^2 in the coefficient of t^3 has degree 2n; every
other term there has degree at most n+1. This is impossible. Hence u
is affine and k quadratic, proving the claim by quadratic safety.
Allowing c or lambda to depend on z makes the determinant spatially
variable and is outside this elimination argument.

### An arbitrary WY coefficient with constant nonzero W^2 coefficient

For the same family with lambda nonzero constant, now allow c=c(Z)
to be an arbitrary nonconstant polynomial. There are no pure members,
with no degree bound on c,h,u,k. Together with the constant-c result,
every pure member of

\[
 \lambda W^2+Wh(Z)+c(Z)WY+aY^2+u(Z)Y+k(Z),\qquad\lambda\ne0,
\]

is quadratic. This conclusion uses a curve period argument; the preceding
polynomial leading-coefficient calculation alone does not prove it.

Let r=1/t, F(z)=c(z)^2-4lambda a, and let C_r be the smooth projective
completion of

\[
 V^2=q_r(z)=F(z)+2\lambda r.
\]

If n=deg c>=1, q_r has degree 2n and is squarefree over C(r): its
derivative is the fixed polynomial F', and its value at a root of F'
cannot vanish identically in r. Thus the curve is geometrically connected
and has genus n-1. Put Delta=-t^2 q_r. The effective phase is N/(2Delta),
where (14) remains the exact algebraic numerator with c=c(z).
The normalized spatial differential is a fixed nonzero constant times

\[
 \eta=r\,\frac{dz}{V}.                                  \tag{17}
\]

**The algebraic source bridge.** Here the fibre matrix depends on z,
so the constant-matrix elimination lemma requires an extension. Work over
A=C(t)[z,q_r^{-1}], and complete the square in W,Y. Polynomial fibre
contractions reduce the top spatial module to its constant fibre class.
Differentiation in either z or t produces the effective phase derivative
and the half-trace correction

\[
 -\tfrac12\operatorname{tr}(M^{-1}M_z)
 =-\tfrac12\partial_z\log\Delta,\qquad
 -\tfrac12\operatorname{tr}(M^{-1}M_t)
 =-\tfrac12\partial_t\log\Delta.
\]

Terms from the moving centre are spatial exact terms. This proves the
connection and seed identification, rather than merely the stationary
equations. The scalar-volume connection is exactly that of (17):
d log(r/V)=-d log t-(1/2)d log q_r. The original Wick functional
extends to the required localization, since

\[
 q_r^{-1}=\frac{t}{2\lambda}
 \sum_{j\ge0}\left(-\frac{tF(z)}{2\lambda}\right)^j
\]

has polynomial z coefficients at every t order. Fibre contractions and
generic reduction therefore do not discard the original unit seed.
There is no conditional integral in an isotropic W direction at t=0.

On the curve with V inverted, the exponential de Rham system splits
under V mapped to -V. Its odd forms are A dz/V, and its odd primitives
are A/V. Their twisted derivative is the half-density relation
partial_z+phi_z-q_r'/ (2q_r), exactly the Gaussian reduction. The time
lift commutes with this involution. A unit horizontal functional on the
odd system extends by zero on the even system to a unit functional on
the full curve system. This is the full system used for period comparison.

The affine total space is smooth: its equation solves
r=(V^2-F(z))/(2lambda), so it is an open part of the affine (z,V)
plane after the indicated parameter and spatial poles are removed.
The phase is necessarily nonconstant in z, since at t=0 it has
Laurent leading term -z^2/(4lambda t)+O(1). Thus the base projection
and phase are algebraically independent, as required by the exponential
period comparison. Additional branch punctures from the localization
may be regular for the phase; dz/V is holomorphic there, so they cause
no extra seed growth.

**A proper marked curve model.** Semistable reduction for smooth projective
curves gives a model after a finite parameter cover; the hypotheses just
checked include geometric connectedness and constant field C(r). The
precise imported theorem is
[The Stacks Project, Theorem 55.18.1, Tag 0CDN](https://stacks.math.columbia.edu/tag/0CDN).
[Section 55.14](https://stacks.math.columbia.edu/tag/0CDB), especially
Example 55.14.1 and Lemma 55.14.3, resolves its nodal surface singularities
while retaining a reduced nodal fibre. These are additional existence
inputs, not period-growth estimates.

Split all marked branches (the punctures and fibres phi=0,1,infinity)
after a further finite cover, and resolve any newly introduced equations
uv=R^m as in that section. The marked sections then miss the nodes;
separating their smooth-point collisions preserves reduced normal crossings.
The value-one marking removes the remaining phase indeterminacies by
the same zero/pole-crossing argument given above. The result is a proper
regular curve model with reduced nodal central fibre, a morphism phi
to P1, and separated horizontal sections. Its finite components may
have positive genus. Proper nearby-cycle comparison and the full gluing
complex retain their entire cohomology; they must not be replaced by a
count of critical points or by a genus-zero computation.

**Positive order on every component.** Normalize v(r)=1 and use the
relative differential at fixed parameter on a generic central component.
For any meromorphic function f, v(df)>=v(f): differentiating its leading
residue cannot lower the parameter valuation. If z has nonconstant
finite residue, or reduces to a finite alpha with F(alpha) nonzero,
then v(V)=0 and (17) has order at least one.

Suppose z reduces to a root alpha of F of multiplicity m, and put
N0=v(z-alpha)>0 and rho_z=v(dz)>=N0. If mN0<1, the F term dominates
and v(V)=mN0/2; hence

\[
 v(\eta)\ge1+N_0-mN_0/2>0.
\]

The strict inequality follows from mN0<1 (and is immediate for m=1).
If mN0>1, v(V)=1/2, and v(eta)>=1/2+N0>0. When mN0=1 and there
is no additional cancellation the same estimate applies. If cancellation
raises v(V) above 1/2, differentiate the curve relation at fixed r:

\[
 2V\,dV=F'(z)\,dz,\qquad v(F')=(m-1)N_0=1-N_0.
\]

It follows that v(eta)=1+v(dV)-v(F')>=v(V)+N0>0. This treats even
the components requiring a refined scale at a multiple branch collision.
At spatial infinity, put v(z)=-N0<0. Since deg c=n>=1,
v(V)=-nN0 and v(dz)>=-N0, giving

\[
 v(\eta)\ge1+(n-1)N_0>0.
\]

All generic centres have been covered. This order calculation holds
on every central component, independently of whether its phase is finite.
It uses the curve relation, not affine-disc coordinates on a rational curve.

The only horizontal poles of eta are at infinity. In x=1/z its local
form is a unit times r x^{n-2}dx: it is logarithmic for n=1 and
holomorphic for n>=2. It is holomorphic at every finite branch point,
as dz/V=2dV/F'(z) at a generic parameter. Thus the complete finite-piece,
finite/finite neck, and finite/pole tail estimates apply. Positive generic
orders on both sides give a positive factor times logarithmic neck length;
descending exponential tails integrate every fixed pole order. Beyond a
node chart the algebraic cover, differential and inverse-phase branches
have only algebraic parameter growth. One can use z away from branch
discs and V within them; the finitely many algebraic transition scales
and bounded-degree projected phase curves give that bound. The exponential
truncation dominates it. Regular puncture circles have at most logarithmic
amplitude, so the positive component order persists.

Proper comparison with these complete glued complexes therefore makes
every zero-exponential seed period tend to zero. Extending a hypothetical
unit seed quotient to an actual curve period would make that seed value
one, a contradiction. This proves the nonconstant-c exclusion.

## A varying quadratic W coefficient: clean localization

Let

    P=ell(Z)W^2+W(h(Z)+c(Z)Y)+aY^2+u(Z)Y+k(Z),

where all displayed coefficient functions are polynomials, a is constant,
and L=deg ell>=2. Choose the monic gcd and define

    g=gcd(ell,c^2),   e=ell/g,   C=c^2/g,
    f=4ae-C,         J0=-2e+t f,
    R=C[z f-2e h]+4e^2 c u.

Then gcd(e,C)=1. Assume:

1. Every root of g has multiplicity at most two.
2. At every root alpha of g, (alpha,h(alpha))!=(0,0).
   Equivalently, g(0)!=0 or h(0)!=0.
3. If either e or C is spatially nonconstant, then R is not the zero
   polynomial. No R condition is imposed when both are constant.

Under these assumptions P is not pure. There is no degree bound on any
coefficient and no multiplicity bound on ell outside the gcd g.

The moving-source condition in item 3 is automatic if e is neither a
nonzero constant nor a nonzero scalar multiple of z. In particular this
theorem contains the coprime case g=1, L>=2. For c=0 one has g proportional
to ell and e constant, recovering the fixed-cover simple/double-root branch.

Two useful consequences have simpler hypotheses:

- If gcd(ell,c)=1 and deg ell>=2, there are no pure members, without
  any root-multiplicity bound. For example ell=(z-1)^5 and c=z are covered,
  with arbitrary a,h,u,k.
- If c=0, deg ell>=2, all roots of ell have multiplicity at most two,
  and ell(0)!=0 or h(0)!=0, there are no pure members.

Shared factors are also allowed. For example ell=(z-1)^2(z-2)^3 and
c=z-1 have g=(z-1)^2 and e=(z-2)^3. The fixed-root source condition
is automatic because its only root is nonzero, and the moving-source
condition is automatic by the criterion above. Thus there is no pure
member for any a,h,u,k in this example. The proof of these statements
retains the original Gaussian seed; it does not define a Wick
expectation of 1/ell.

### Exact source criterion and localization of the original seed

Put K=C(t), d=1-2at. The W,Y precision matrix, source vector, and
determinant of the full phase Sfull=-zw-y^2/2+tP are

    M=[[-2t ell,-tc],[-tc,d]],   Lvec=(th-z,tu),
    Delta=-2t ell+t^2(4aell-c^2)=t g J0.

Set X=A3_K and use the twisted differential d+dSfull throughout.
Let D be the reduced divisor Delta=0. Its fixed and moving components
are disjoint over K.
Indeed at a root of g, e and C cannot both vanish, so J0 does not vanish
identically as an element of K at that root.

If either e or C is nonconstant, J0 is irreducible and separable in K[z].
For gcd(e,f)=1 makes this t-linear polynomial primitive irreducible in
C[z,t]; its positive spatial degree makes it primitive over C[t] as well.
The ratio e/f cannot be constant in this case. If both e and C are
constant, J0 is an invertible scalar in K and there is no moving spatial
divisor.

A null vector of M on Delta=0 is (d,tc). Its source is

    S=d(th-z)+t^2 c u.

For the moving divisor the following criterion is exact:

    J0 divides S in K[z]  iff  R=0 in C[z].

To see this, the parametrization t=2e/f of the irreducible moving pencil
has nonzero denominator f, and direct substitution gives

    f^2 S|_(t=2e/f)=C[z f-2e h]+4e^2 c u=R.

Thus hypothesis 3 makes S a unit in the finite field K[z]/(J0), and after
base change its value is nonzero at every geometric moving root. Reducing
R=0 modulo e gives -z C^2=0 mod e. Since gcd(e,C)=1, this forces e|z,
proving the automatic criterion stated above.

At a fixed root alpha of g, both ell(alpha) and c(alpha) vanish. The
restricted phase is

    (t h(alpha)-alpha)w-d y^2/2+t u(alpha)y+t k(alpha).

The W source is nonzero over K by hypothesis 2. On a moving geometric
root, set Y'=Y-tcW/d. The restricted phase is

    -d Y'^2/2+t u Y'+(S/d)W+t k,

whose W source is again nonzero. In either case partial_W+b is bijective
on polynomials for b!=0, with finite inverse

    f -> sum_(j=0)^deg f (-1)^j b^(-j-1) partial_W^j f.

Consequently every reduced determinant component has an acyclic restricted
exponential de Rham complex. The reduced divisor is a disjoint union of
smooth planes after a finite separable field extension.

Here is the full localization argument, including higher poles. The
untwisted logarithmic residue exact sequence is
[The Stacks Project, Section 50.15, Lemmas 50.15.2 and 50.15.7, Tag 0FMU](https://stacks.math.columbia.edu/tag/0FMU).
Its smooth-divisor hypotheses hold for D inside X. Since dS_full is
regular across D, the same sequence is a short exact sequence for the
twisted differentials:

\[
 0\longrightarrow\operatorname{DR}(X)
 \longrightarrow\operatorname{DR}_X(\log D)
 \longrightarrow\operatorname{DR}(D)[-1]\longrightarrow0.
                                                               \tag{C2}
\]

Locally take a normal equation x for D and filter meromorphic logarithmic
forms by x^{-m}Omega(log D), m>=0. On the quotient of consecutive pole
orders, m>=1, the normal connection residue is the nonzero scalar -m.
The regular twisting term contributes x partial_x S_full to that residue,
which vanishes on D. The quotient is therefore the tangent twisted
complex tensored with the two-term normal complex of differential -m;
it is contractible. Passing through the finite pole orders and their
filtered colimit proves that the logarithmic complex is quasi-isomorphic
to the full complex localized at Delta. No formal completion is used.

Since the last term of (C2) is acyclic, localization of the original
integration source at Delta is an isomorphism on all de Rham cohomology.
The canonical original-to-meromorphic map commutes with the total
integrable parameter connection. Although the moving divisor varies with t
and the intermediate logarithmic pole filtration need not be parameter-stable,
the resulting cohomology isomorphism is horizontal. Its inverse
transports the original horizontal functional and preserves its unit
volume seed. This is the reason localization is valid here; there is
no assertion that inverse powers of ell have polynomial Wick expansions.

### Normalized cover, phase, and horizontal poles

Quadratic elimination gives the effective phase

    phi=t k+(1/2)Lvec^T M^-1 Lvec

and the Gaussian volume, up to a nonzero constant, is

    eta=r dz/U,    U^2=G+2r ell,    r=1/t,
    G=c^2-4aell=g(C-4ae).

The connection retains both the spatial and parameter half-trace terms.
The source clean condition also forces phi spatially nonconstant. Indeed,
if N0=Lvec^T adj(M)Lvec, the Schur identity is

    d N0=S^2+Delta(tu)^2.

On every reduced determinant component its right side is nonzero; hence
the numerator of phi does not vanish there and phi has an actual spatial
pole. Such a component exists because Delta has spatial degree at least
L>=2. This includes the case with no moving residual divisor.

Write g=s^2 b with s and b squarefree and coprime, allowing b constant.
Then the generic normalization of the cover is

    V^2=b[C-4ae+2r e],    U=sV,
    eta=r dz/(sV).

The bracketed pencil is either a nonzero scalar in C(r), or separable of
positive degree. It is coprime to b over C(r). Consequently this normalized
equation is generically smooth. After deleting finitely many parameter
exceptions its affine total space is smooth as well: at a fixed b-root
the z derivative is b' times the nonzero residual pencil; at a moving
branch the r derivative is 2be, nonzero by gcd(e,C)=1. It is
geometrically connected whenever
the squarefree right side has a spatial root. Otherwise its geometric
normalization is the disjoint union of two P1s, and all period estimates
are applied componentwise, as in the fixed-cover branch. Passing to a
finite parameter cover if needed does not affect the order estimates.

At each generic fixed root of g of multiplicity one, eta is regular in
the branch coordinate. At each multiplicity-two root it is logarithmic.
At a moving simple root it is regular. Its infinity order is logarithmic
or holomorphic, because the generic raw degree

    E=max(deg ell,2 deg c)>=L>=2

and eta is dz divided by the square root of a raw polynomial of degree E.
Here c=0 is interpreted by omitting 2 deg c. Thus all finite horizontal
poles are logarithmic. The phase polar divisor is additionally deleted
before applying the exponential comparison theorem. On this smooth affine
open set the projection and phase are regular and algebraically independent:
the actual spatial phase pole proves nonconstancy on each connected
component of the geometric generic curve. These check the hypotheses of
the same exponential period comparison used above.

### Positive volume order over finite centres

Use the complete proper marked-curve model constructed above, with
semistable reduction, nodal resolution and separated phase and puncture
sections. Retain all positive-genus component cohomology and the full
annulus gluing complex. Normalize every divisorial
valuation by v(r)=1. Relative differentiation does not lower valuation.
The estimates below use U=sV as a rational function on the normalization;
U need not be a normal affine coordinate.

First suppose G is not identically zero. At a finite centre alpha put

    l=ord_alpha ell,   q=ord_alpha G,   N=v(z-alpha)>0.

At any point where G vanishes, one has

    min(q,l)=min(ord_alpha c^2,ord_alpha ell)<=2,

where the minimum is zero if alpha is not a shared root. This equality
still holds when the leading terms of c^2 and 4aell cancel.
Here the order of the zero polynomial is infinity, so c=0 is included.

If the G term dominates, qN<1+lN, then

    v(eta)>=1+(1-q/2)N.

For q<=2 this is at least one. If q>2, the minimum condition gives
l<=2<q and N<1/(q-l). The displayed bound is therefore strictly greater
than or equal to the positive endpoint bound

    (q/2+1-l)/(q-l)>=1/2.

If the r ell term dominates, then (q-l)N>1, forcing q>l and hence l<=2.
Thus

    v(eta)>=1/2+(1-l/2)N>=1/2.

At equal valuation without cancellation either estimate applies. With
cancellation, necessarily q>l, N=1/(q-l), and v(U)>qN/2. The leading
coefficient of G'+2r ell' is proportional to q-l, so

    v(G'+2r ell')=(q-1)N.

Using 2U dU=(G'+2r ell')dz yields

    v(eta)>=1+v(U)-(q-1)N
          >1+(1-q/2)N>=1/2.

The same estimates include q=0 or l=0 when needed. Finite components
with nonconstant z residue have v(eta)>=1. Thus every finite centre has
strictly positive normalized volume order.

If G=0 identically, coprimality of e and C forces e and C to be constants
(including c=0,a=0). The normalized curve is fixed and eta is a nonzero
constant times sqrt(r) dz/sqrt(ell). Since all roots of ell now belong
to g and have multiplicity at most two, the fixed-cover logarithmic
differential argument applies directly. To spell it out, on the fixed
normalization of v^2=ell(z), omega=dz/v is regular at simple roots and
logarithmic at double roots. At infinity, for L even use x=1/z to get
a unit times x^(L/2-2)dx; for L odd use z=x^-2 to get a unit times
x^(L-3)dx. These are logarithmic or holomorphic for L>=2.
On any marked model of this fixed curve, omega has nonnegative generic
relative order. A component dominating the original curve has order
zero. At a component centred at a point, a local coordinate x has
v(x)=N>0 and v(dx)>=N: regular forms have nonnegative order, and
dx/x has order at least zero. Multiplication by sqrt(r) makes the
order at least 1/2 on every component, including either component of
a disconnected normalization. This proves the required positivity
without a nonconstant G.

### Infinity and complete period comparison

When G is nonzero, at infinity put v(z)=-N<0, L0=deg ell>=2 and D0=deg G>=0.
If D0>=L0, G dominates and

\[
 v(\eta)\ge1+(D_0/2-1)N\ge1.
\]

If L0>D0, the transition is at N*=1/(L0-D0). For N<N*, G
dominates with the same first bound. For D0>=2 it is at least one.
For D0=0, N<1/L0 gives 1-N>1-1/L0>=1/2; for D0=1,
N<1/(L0-1) gives 1-N/2>1-1/(2(L0-1))>=1/2.
For N>N*, the r ell term dominates, and

\[
 v(\eta)\ge\tfrac12+(L_0/2-1)N\ge\tfrac12.
\]

At the transition without cancellation these bounds apply. With
cancellation v(U)>-D0 N/2. The cancelled leading terms of G and
2r ell imply that the leading coefficient of G'+2r ell' is proportional
to D0-L0, which is nonzero. Therefore

\[
 v(G'+2r\ell')=-(D_0-1)N,
\]

and the differential identity gives

\[
 v(\eta)\ge1+v(U)+(D_0-1)N
 >1+(D_0/2-1)N>0.
\]

For D0=0 the first derivative G' is zero, but the second term has
exactly the stated valuation and coefficient; this case is included.

When G=0 the fixed-cover infinity estimate applies instead. Therefore
every central component in a complete marked model has strictly positive
normalized seed volume order. The generic horizontal volume poles are
only logarithmic. The complete compact-piece, neck, and tail
estimates established in the preceding double-cover argument retain all
positive-genus cycles and gluing maps, and imply that
every zero-exponential seed period tends to zero.

The Gaussian system is the odd sector of this full curve system; extending
a hypothetical transported horizontal seed functional by zero on the even
sector is valid even in the disconnected case. Ordinary-parameter analytic
comparison then supplies a full curve period with the same unit seed,
contradicting the positive-order bound. Thus P is not pure.

### The moving-source exception

The degree threshold is necessary: ell=z, c=1, h=u=k=a=0 gives
P=ZW^2+WY, whose angular weights are both strictly negative. It is
pure and safe, although the clean-source conditions hold.

The condition R!=0 cannot be inferred just from the clean fixed roots.
For example, alpha!=0 and

    ell=z(z-alpha)^2, c=z-alpha, h=-1/2, a=u=0

give g=(z-alpha)^2, e=z, C=1, J0=-2z-t, and R=0. The moving null source
vanishes, so the actual localization map may have a source kernel. This
example is not claimed pure; it is excluded from the theorem. Likewise
e constant with R=0 is realizable. Those exceptional branches require a
separate actual-seed retention argument.

## A uniform obstruction on the moving-null-source-zero locus

For every integer nu>=1 and every polynomial k in C[Z], the phase

    P_nu = Z^(2nu+1) W^2 + Z^nu W Y - W/2 + k(Z)

fails Gaussian purity. Here E[Z^p W^q Y^r] is zero unless p=q and r
is even, in which case it equals p! (r-1)!!. All coefficients of k may
be arbitrary complex numbers; no degree bound is imposed.

This concerns the actual Wick functional. No source-localization
isomorphism, period-selection assertion, or finite-prefix sufficiency
is used. The result includes nu=1, ell=Z^3,c=Z,h=-1/2,a=u=0.
More generally ell=Z^(2nu+1),c=Z^nu have g=gcd(ell,c^2)=Z^(2nu),
e=Z,C=c^2/g=1,J0=-2Z-t, and moving null source S=-Z-t/2=J0/2.
Thus the residual moving source vanishes for every nu in the theorem.
The fixed root is clean, since h(0)=-1/2 is nonzero.

### Exact marked base periods

Put P0_nu=Z^(2nu+1)W^2+Z^nu WY-W/2 and
B_N^(nu)(t)=E[Z^N exp(tP0_nu)], N>=0. All series are formal.
First integrate the ordinary Gaussian Y:

    B_N^(nu)(t)=E_ZW[Z^N exp(t Z^(2nu+1)W^2
                                  +t^2 Z^(2nu)W^2/2-tW/2)].

Take i copies of the first quadratic-W term, j-i of the second,
and b copies of the linear-W term, where j>=0 and 0<=i<=j.
Angular balance forces

    b=N+(2nu-1)i+2(nu-1)(j-i)
     =N+2(nu-1)j+i.

The contracted W degree is N+2nu*j+i. The total t degree is N+2nu*j.
The sign from (-W/2)^b is (-1)^(N+i), and all powers of 2 in the
coefficient multiply to 2^(N+(2nu-1)j). Consequently

    B_N^(nu)(t)=(-t/2)^N sum_{j>=0} 2^(-(2nu-1)j)t^(2nu*j)
        sum_{i=0}^j (-1)^i (N+2nu*j+i)!
          / [i!(j-i)!(N+2(nu-1)j+i)!].

The factorial quotient in the inner sum is
(2j)! binom(N+2nu*j+i,2j). The finite-difference identity

    sum_{i=0}^j (-1)^i binom(j,i) binom(N+2nu*j+i,2j)
      =(-1)^j binom(N+2nu*j,j)                       (K1)

follows by expressing the binomial coefficient as [x^(2j)] of
(1+x)^(N+2nu*j+i). The sum becomes

    [x^(2j)] (1+x)^(N+2nu*j) (1-(1+x))^j
      =(-1)^j [x^j] (1+x)^(N+2nu*j).

It follows exactly that

    B_N^(nu)(t)=(-t/2)^N sum_{j>=0} (-1)^j (2j-1)!!
          2^(-2(nu-1)j) binom(N+2nu*j,j)t^(2nu*j).    (K2)

Use (-1)!!=1. This is a formal identity for every N,j and every
integer nu>=1; it is not inferred from finitely many moment values.

### Polynomial perturbations cannot cancel the leading factorial term

Purity first forces k(0)=0 because E[P_nu]=k(0). Write

    k(Z)=sum_(ell>=1) k_ell Z^ell,
    A=sum_(ell>=1)|k_ell| < infinity,
    k(Z)^m=sum_N a_(m,N) Z^N.

Then a_(m,N)=0 for N<m, and sum_N |a_(m,N)|<=A^m.
Define positive numbers

    L_j^(nu)=(2j-1)!! 2^(-2(nu-1)j) binom(2nu*j,j),
    L_0^(nu)=1.

For all N,j>=0, the ratio of the absolute coefficient in (K2) to
L_j^(nu) is

    2^(-N) binom(N+2nu*j,j)/binom(2nu*j,j)
      = product_(ell=1)^N (2nu*j+ell)/(2((2nu-1)j+ell)) <= 1.  (K3)

Every factor is at most one, since nu>=1. For j>=1 the exact ratio is

    L_j^(nu)/L_(j-1)^(nu)
      =2^(-(2nu-1)) (2j)(2j-1)/j^2
         * product_(s=0)^(2nu-1)(2nu*j-s)
           / product_(s=0)^(2nu-2)((2nu-1)j-s).        (K4)

Dividing this finite product by j and taking j to infinity gives

    (L_j^(nu)/L_(j-1)^(nu))/j
      -> c_nu=2^(3-2nu)(2nu)^(2nu)/(2nu-1)^(2nu-1)>0. (K5)

Hence L_j^(nu) is eventually increasing and tends to infinity. It
therefore eventually exceeds the finite set of earlier terms, so for
all sufficiently large q,

    max_(0<=j<q) L_j^(nu)=L_(q-1)^(nu).

The full actual seed expands by (K2) as

    Phi_(nu,k)=B_0^(nu)
       +sum_(m>=1,N>=m) (a_(m,N)/m!)t^m B_N^(nu).

At coefficient t^(2nu*q), a perturbation contributes only when
m+N+2nu*j=2nu*q. Thus m+N is a positive multiple of 2nu and j<=q-1.
All rearrangements are t-adically locally finite: m+N<=2nu*q, N>=m,
and j>=0. The coefficient estimate (K3) and the l1 norm bound yield,
for all sufficiently large q,

    [t^(2nu*q)]Phi_(nu,k)=(-1)^q L_q^(nu)+E_q,
    |E_q| <= (exp(A)-1)L_(q-1)^(nu).                 (K6)

The actual sums for a fixed coefficient are finite; extending their
nonnegative numerical bound to sum_(m>=1) A^m/m!=exp(A)-1 is harmless.
By (K4)-(K5), L_(q-1)^(nu)/L_q^(nu)->0. Therefore

    [t^(2nu*q)]Phi_(nu,k)/((-1)^q L_q^(nu)) -> 1.      (K7)

In particular these coefficients are eventually nonzero. This proves
that Phi_(nu,k) cannot be the unit series for any polynomial k. When
A=0 the unperturbed term is already nonzero for every q>=1.

No finite cutoff uniform in k or nu is claimed. The proof is an
arbitrary-degree obstruction derived from the exact all-order kernel.

### An explicit estimate when nu=1

For nu=1, L_j=(2j-1)!!binom(2j,j) is increasing from j=0, since

    L_j/L_(j-1)=2(2j-1)^2/j >= 2,   j>=1.

Thus (K6) holds for every integer q>=1, with exact relative error

    |[t^(2q)]Phi_(1,k)-(-1)^q L_q|/L_q
      <=(exp(A)-1)q/[2(2q-1)^2]
      <=(exp(A)-1)/(2q).                             (K8)

Every integer q>=1 with 2q>exp(A)-1 therefore gives a nonzero even
coefficient. This gives an explicit bound in the nu=1 case and covers
arbitrary complex k.

The exact checker supplies direct Wick controls for the marked formula,
finite-difference and ratio identities. The coefficient extraction and
uniform estimates above prove the assertion for every nu,N,j,q and
every complex polynomial k; those quantifiers do not come from the
finite controls. This theorem uses no analytic period comparison.

## All-degree radial obstruction and angular support

For arbitrary polynomials B,g in C[R], with R=ZW, the phase

\[
 P=aY^2+B(R)Y+g(R)
\]

is pure only when P=0. The endpoint of the factorial functional
E[R^n]=n! must be retained in proving this statement.

### A relative leading-coefficient bound

Let K=C(t), let phi in K[x] have degree D>=1 and leading coefficient
a_D, and let H be a nonzero algebraic scalar with H'/H in K. If

\[
 Dv_\infty(H)>v_\infty(a_D),                             \tag{15}
\]

the exponential de Rham system of the pair (A1,{0}), with normalized
volume H dx, has no rational horizontal quotient of its differential
seed submodule taking [1] to one. Algebraic normalization is handled
on a finite cover. This is a relative assertion, proved below; it is
not an application of the ordinary rank-(D-1) criterion.

Put A=K[x] and D_x=partial_x+phi_x. The relative complex has degree-zero
term xA and degree-one term A dx, so its quotient and connection are

\[
 M_{\rm rel}=A/D_x(xA),\qquad
 \nabla[q]=[q_t+(\varphi_t+H'/H)q].
\]

The two twisted derivatives commute and preserve the endpoint condition.
The leading degree of D_x(xq), for deg q=m, is m+D. Thus the classes
1,x,...,x^{D-1} are a basis and the rank is D. If N_abs=A/D_x A,
the exact sequence is

\[
 0\longrightarrow Kb\longrightarrow M_{\rm rel}
 \longrightarrow N_{\rm abs}\longrightarrow0,\qquad
 b=[\varphi_x].
\]

Indeed D_x is injective on polynomials when phi is nonconstant, so
D_x A/D_x(xA) has dimension one. The extra boundary class satisfies

\[
 \nabla b=(\varphi_t(t,0)+H'/H)b,                       \tag{16}
\]

because D_x(phi_t-phi_t(t,0)) is a relative relation. The boundary
summand is therefore retained with its actual parameter action.

At an ordinary parameter the closed rapid-decay cycles give the perfect
pairing for N_abs. Add one path from x=0 to a decay sector at infinity.
Its value on b is -H exp(phi(t,0)), by the fundamental theorem of
calculus, and is nonzero. The closed-cycle pairing, plus this path,
is consequently a perfect relative pairing of rank D. Equivalently,
the relative de Rham complex is the cone of the absolute complex to
point restriction, shifted by minus one; the corresponding chain
cone allows a boundary at 0. The ordinary perfect pairing and elementary
point evaluation identify the two exact sequences. Flat transport
makes the pairing horizontal. Every analytic horizontal functional on
M_rel is thus a relative rapid-decay period. The ordinary pairing used
here is the Hien--Roucairol comparison already specified above; the
relative extension follows from this explicit extra path and cone.

For the full regular part at infinity, use the same proper marked P1
model, now marking the endpoint 0, the puncture infinity and the
fibres phi=0,1,infinity. Only infinity is removed; 0 is an allowed
boundary vertex, not another puncture. The comparison of nearby cycles
commutes with the cone defining the relative complex and with proper
direct image. Thus the complete regular-part calculation retains
finite-phase component chains and annulus gluing, now relative to the
endpoint section wherever its phase is finite. If phi(t,0) has a
nonzero polar action, the point term has no zero-exponential part.
Otherwise the endpoint is a separated smooth boundary vertex of a
finite piece. This explicitly accounts for the extra rank-one term.

On an exact component chart x=xi+r^rho X, r=1/t, finite phase forces
the top coefficient a_D r^{D rho} to have nonnegative valuation.
Therefore rho>=-v(a_D)/D. The generic order of H dx on every finite
component is at least v(H)-v(a_D)/D, strictly positive by (15).
The only horizontal pole of dx is infinity, always phase-irregular.
In particular the volume is regular at the allowed endpoint 0.
The positive generic factor consequently extends uniformly to paths
ending at that endpoint. At finite/finite nodes the coefficient of
du/u has positive order on both adjacent components; its integral
is O(r^epsilon log r). At finite/pole nodes the positive order of
the finite component remains after descending exponential integration.
The uniform tail continuation has only algebraic volume growth and is
dominated by the same exponential cutoff proved above. After finite
ramification these are exactly the positive-order node calculations
already displayed, with H dx in place of r dx.

Every zero-exponential relative seed period therefore tends to zero.
A unit formal seed functional restricts to a rational quotient on its
differential seed span, extends analytically to the full module at an
ordinary point, and is represented by a flat relative cycle by the
complete relative pairing. Its zero-exponential seed value would be
one, contradicting this decay. This proves (15). The endpoint class
b can itself have a constant period; this does not give a unit period
on [1]. For example phi=t x^2 has a constant endpoint period on 2tx
while its seed periods scale as t^{-1/2}.

### Application to radial phases

Put d=1-2at and H=d^{-1/2}. Eliminating Y gives

\[
 \varphi=-R+tg(R)+\frac{t^2B(R)^2}{2d},\qquad
 L(q)=H E_R[q\exp(\varphi+R)].
\]

The factorial integration-by-parts identity is
E_R[f'-f]=-f(0). Coefficientwise it implies L(D_R f)=0 whenever
f belongs to R K[R]. The endpoint term is zero for these relations,
and time differentiation gives precisely the relative connection.
Thus the actual moment functional descends to M_rel with seed Phi_1.
No convergence of the original exponential integral is assumed.

If D=max(deg g,2deg B)>=2, it is the actual degree of phi. Its top
coefficient, with g_D and c_D from g and B^2 not both zero, is

\[
 a_D=tg_D+\frac{t^2c_D}{2d}.
\]

For a=0 its infinity valuation is -1 or -2, whereas v(H)=0.
For a nonzero,

\[
 a_D=\frac{g_D(r-2a)+c_D/2}{r(r-2a)},\qquad r=1/t.
\]

Its valuation is -1, or zero when the constant numerator cancels;
it cannot vanish identically. Since v(H)=1/2, (15) is strict in
every case. Purity is excluded. The remaining case has B=b constant
and g=g0+g1 R. Its exact period is

\[
 \Phi_1=\frac{\exp(tg_0+t^2b^2/(2d))}{\sqrt d(1-g_1t)}.
\]

Rational-exponential separation forces d(1-g1t)^2=1 and its exponent
zero. Hence a=g1=b=g0=0. This proves the all-degree radial assertion.

### An elementary quartic proof

Let R=ZW, so E[R^n]=n!. Consider the full radial H(2,1) family of
total degree at most four:

\[
 P_0=\alpha R^2+\beta RY+aY^2+\gamma R+bY+k.
\]

Its only pure member is P0=0. This is an all-order statement for every
complex choice of the six coefficients.

Indeed, integrating Y and then expanding the factorial R functional
gives

\[
 \Phi_1=\frac{e^{kt+b^2t^2/(2d)}}{\sqrt d\,(1-C)}
             F\!\left(\frac{A}{(1-C)^2}\right),\qquad
 d=1-2at,
\]

\[
 A=\alpha t+\frac{\beta^2t^2}{2d},\quad
 C=\gamma t+\frac{\beta b t^2}{d},\quad
 F(s)=\sum_{n\ge0}\frac{(2n)!}{n!}s^n.                \tag{10}
\]

For example E[R^{2n}e^{CR}]=(2n)!/(1-C)^{2n+1}; expansion of
exp(A R^2) proves (10) coefficientwise. All scalar factors preceding F
are convergent germs at zero. If A is not identically zero, then
q=A/(1-C)^2 is a nonzero analytic germ of order m=1 or 2. Write
q(t)=v(t)^m with v a convergent germ and v'(0) nonzero, choosing an
analytic mth root of its unit factor. A unit period would make F(q(t))
convergent. Inverting v then makes F(v^m) convergent, contradicting
the unbounded nth root of its nonzero coefficients (2n)!/n!.
Consequently A=0, forcing alpha=beta=0. The remaining period is

\[
 (1-\gamma t)^{-1}(1-2at)^{-1/2}
      \exp\!\left(kt+\frac{b^2t^2}{2(1-2at)}\right).
\]

The rational-exponential argument forces
(1-gamma t)^2(1-2at)=1 and its exponent zero. Thus gamma=a=b=k=0,
proving the claim without any finite-moment inference.

### Weakly one-sided angular support in every degree

The all-degree radial theorem gives a uniform support obstruction.
Suppose all monomials of an H(2,1) phase have nonnegative Z-minus-W
weight, with no degree restriction. Write P=P0+Pplus by angular weight.
In a pure moment,
any factor from Pplus leaves strictly positive weight and cannot
contract. Therefore Phi_1(P)=Phi_1(P0). The weight-zero part has exactly
the all-degree radial form aY^2+B0(R)Y+g0(R): the constant fibre
coefficient forbids a nonconstant coefficient of Y^2. Purity forces
P0=0. Conversely this condition gives strictly positive
weight, so P is pure and safe with cutoff m>deg Q. The same conclusion
holds with all weights nonpositive, after exchanging Z and W.

Thus every H(2,1) phase with weakly one-sided angular support is pure
exactly when its weight-zero part vanishes. A counterexample in any
degree must have monomials of both angular signs in every allowable
centred Gaussian frame retaining an ordinary fibre direction. In
degree at most four the same conclusion also follows from the elementary
quartic proof, independently of the relative leading-coefficient lemma.
The remaining mixed-sign phases are not classified here.

## An exact classification of a mixed quartic family

There is a complete algebraic classification beyond the W-affine and
quadratic-W families. For arbitrary complex alpha,beta,gamma,delta,e,eta, let

\[
 P=\alpha Z^3+\beta W^3+\gamma Z^2WY+\delta ZW^2Y
       +e((ZW)^2-2)+\eta(Y^2-1).
\]

Then P is pure if and only if e=eta=0 and either alpha=gamma=0 or
beta=delta=0. Every such phase is safe, with E[QP^m]=0 for m>deg Q.
In fact the six equations E[P^j]=0 for j=2,3,4,5,6,8 already force
these forms. This is a finite algebraic completeness certificate: after
those equations the remaining forms have an all-order weight proof.

Put A=alpha beta and D=gamma delta, and write Mj=E[P^j]. Exact Wick
contraction gives

\[
 M_2=2(6A+6D+10e^2+\eta^2),\qquad
 M_3=8(81e(A+D)+9\eta D+74e^3+\eta^3).
\]

The moments through order five depend on alpha,beta,gamma,delta only
through A,D. To see why this is valid even when some coefficients vanish,
let i,j,k,l be the numbers of the first four summands in a surviving
Wick monomial. Angular balance gives 3(i-j)+(k-l)=0, and Y parity makes
k+l even. Thus i-j is even. If i differs from j, then |i-j|>=2 and
|k-l|>=6, requiring at least eight factors. Below that order i=j and
k=l. No division by alpha or gamma is used.

Suppose first eta is nonzero. Scale the whole phase by 1/eta, so eta=1.
The equations M2=M3=0 give

\[
 A=-D-5e^2/3-1/6,\qquad D=(122e^3+27e-2)/18.
\]

After these substitutions the next two moments are M4=-4F(e) and
M5=-16G(e)/3, where

\[
\begin{split}
F(x)&=24400x^5+2364x^4+6864x^3-652x^2+468x-63,\\
G(x)&=4196800x^6+537276x^5+1192320x^4-13760x^3
       +63780x^2-1215x-722.
\end{split}
\]

Their leading coefficients stay nonzero modulo 7. The 11 by 11
Sylvester determinant, with six shifted F rows followed by five shifted
G rows, is 5 modulo 7. Its integer determinant is therefore nonzero,
so F and G have no common complex root. This excludes eta nonzero.
The checker retains both polynomials, forms this matrix, and checks its
determinant by modular row elimination; it does not infer coprimality
from numerical roots.

With eta=0, M2=0 gives A+D=-5e^2/3, and then M3=-488e^3 forces
e=0 and D=-A. For the four remaining terms,

\[
 M_6=9!\,(20A^3+180A^2D+540AD^2+300D^3).
\]

Substitution D=-A reduces this to 80 times 9! A^3. Thus A=D=0.
Modulo the ideal (alpha beta,gamma delta), the eighth moment is

\[
 M_8=420\,12!\,(\alpha^2\delta^6+\beta^2\gamma^6).
\]

The two displayed summands have product A^2 D^6=0. Their sum vanishes,
so each vanishes, giving alpha delta=beta gamma=0. Together with
alpha beta=gamma delta=0 this says that every product of a coefficient
from (alpha,gamma) and one from (beta,delta) is zero. Hence one entire
pair is zero, including all coefficient-zero cases.

Conversely, the surviving monomials have angular weights 3,1 or -3,-1,
all of one strict sign. Each factor of P contributes weight of absolute
value at least one, whereas a monomial of Q has weight of absolute
value at most deg Q. Therefore E[QP^m]=0 for m>deg Q. This proves
both the all-order converse and the marker cutoff.

## Long's quartic under polynomial quadratic-form automorphisms

The earlier directional test excludes the affine Gaussian orbit of the
unsafe comparison family. For Long's particular quartic there is a
stronger algebraic exclusion. Let T be any polynomial automorphism of
C[z,w,y] preserving Q=2zw+y^2. Then

\[
 T\!\left((1+z)\left(w-\frac{2+z}{2}y^2\right)\right)
\]

cannot have the form a y^2+B(z,w)y+g(z,w) with constant a.
This includes automorphisms of arbitrary degree. The assertion is
purely about membership in this polynomial orbit; no claim that every
Q-preserving polynomial automorphism preserves Gaussian expectation
is needed or used.

Write F=Tz, G=Ty, H=Tw, so 2FH+G^2=Q and F,G,H are nonconstant
polynomial coordinates. The transformed quartic is

\[
 L=(1+F)H-\frac{(1+F)(2+F)}2G^2
   =\frac{1+F}{2F}\left(Q-(1+F)^2G^2\right).          \tag{11}
\]

Work over C(z,w) and let f=deg_y F, g=deg_y G. If f>0 and f+g>1,
the bracket in (11) has unique top degree 2(f+g)>2. Its rational
prefactor has degree zero in y, so L has degree 2(f+g)>2. The only
remaining case with f>0 is f=1,g=0. Then the identity 2FH+G^2=Q
forces deg_y H=1. Write F=A y+B and H=D y+E. Comparing y^2 gives
2AD=1 in C[z,w], so A,D are nonzero constants. The y^2 coefficient
of L is 1/2-A^2G^2/2, which cannot be constant because G is a
nonconstant coordinate.

If f=0, F is a nonconstant polynomial in z,w. For g>=2, (11) gives
degree_y L=2g>2. If g=0, comparison of the y^2 coefficient in 2FH+G^2=Q
would make F a unit, a contradiction. Thus g=1. Write G=A y+B with
A in C[z,w] nonzero. The coefficient of y^2 in L is

\[
 \frac{(1+F)\left(1-(1+F)^2 A^2\right)}{2F}.
\]

If it were a constant a, the polynomial identity
(1+F)(1-(1+F)^2A^2)=2aF reduced modulo 1+F would force a=0.
It would then force (1+F)^2A^2=1, impossible because F is nonconstant.
This exhausts all degree cases. Any additional linear complex
orthogonal choice of ordinary fibre direction can be composed with T
and still preserves Q, so the obstruction covers every such frame.

## Sharp cubic barriers throughout the remaining common-root strip

This section proves sharpness of the uniform strict-positive-volume range.
It gives no pure phase or H(2,1) counterexample, and it does not assert a
rational or constant unit period.

Fix integers p>=3 and ell>=2p, put m=ell+p, and choose H!=0. Set

    c=(z-1)^m,
    h=H z^2 (z-1)^ell,
    u=-(H^2/6) z^3 (z-1)^(ell-p),
    k=-(H^3/54) z^5 (z-1)^(ell-2p),
    P=W(h(Z)+c(Z)Y)+u(Z)Y+k(Z).

All coefficients are polynomials, and Y remains the original ordinary
Gaussian. At alpha=1 they have ord c=m and ord h=ell. The integer
conditions are exactly

    ceil(2m/3)<=ell<=m-3.

Thus the construction covers every pair in the remaining local strip.

### The exact cube and zero-action critical branch

More generally, whenever the following expressions are polynomial, put

    u=-h^2/(6zc),    k=-h^3/(54zc^2).

For r=1/t and s=(rz-h)/c, direct expansion gives

    phi=-s^2/2+(us+k)/r
       =-(rz-2h/3)^3/(2rz c^2).                        (CB1)

Indeed, multiplying -2r c^2 phi out gives

    r(rz-h)^2-2uc(rz-h)-2kc^2
      =r^3 z^2-2r^2 zh+(4/3)rh^2-(8/27)h^3/z
      =(rz-2h/3)^3/z.

The displayed polynomial family satisfies these identities exactly.
The apparent pole at z=0 cancels because h contains a z^2 factor.

Set L=2H/3. A critical branch zeta tending to one is determined by

    r=L zeta(zeta-1)^ell.

After r=R^(3ell), choose A^ell=1/L. Then

    zeta=1+A R^3+O(R^6).

For F(r,z)=rz-(2/3)h(z), the centre value is exactly zero, and

    F_z(r,zeta)
      =-L zeta(zeta-1)^(ell-1)[(ell+1)zeta-1]
      =-ell A^(-1) R^(3ell-3)+higher order.

Thus phi and its first two z derivatives vanish exactly at the centre,
while the third derivative is nonzero for nonzero sufficiently small r.
This is a cubic critical point of exact zero action, rather than a
truncated cancellation of its critical value.

### The exact stationary chart and normalized seed order

Use

    z=zeta+R^(2p+3)X.

The displacement is smaller than zeta-1 by the positive factor R^(2p).
Consequently the linear term of F has leading term

    -ell A^(-1) R^(3ell+2p)X,

and all higher Taylor terms have larger order. The denominator in (CB1)
has leading term 2 A^(2m) R^(9ell+6p). Hence

    phi=ell^3/(2A^(2m+3)) X^3+O(R),
    eta=r dz/c=R^(3-p)[A^(-m)+O(R)]dX.                 (CB2)

At each fixed R order the coefficients of the perturbations are
polynomials in X. The chart defines a divisorial disc valuation, realized
by smooth point blowups after the displayed finite ramification. Its
nonconstant finite cubic phase is retained in the complete marked model.
Since v(r)=1, the normalized volume has exact order

    (3-p)/(3ell).

The order is zero for p=3 and negative for p>3, with nonzero leading
coefficient A^(-m). Thus no uniform strictly-positive-volume statement
can extend to any further integer pair in the local strip. The first
zero-order pair is (m,ell)=(9,6); the first negative-order pair is
(12,8), where the order is -1/24.

The nonzero local seed is also explicit. For the fixed cubic C X^3,
C=ell^3/(2A^(2m+3)), choose the dual functional ell_0 with
ell_0(1)=1 and ell_0(X)=0 on
C[X]/(partial_X+3C X^2)C[X]. Applying it coefficientwise to the
positive-order exponential perturbation and volume unit in (CB2)
gives a zero-exponential local seed value

    R^(3-p)[A^(-m)+O(R)].

Every coefficient uses finitely many polynomial terms. Spatial exactness
and parameter horizontality follow from the same moving-coordinate
product rule as in the retained formal stationary construction. A
constant rescaling can normalize the leading coefficient to one. This
does not make the resulting series rational or identically constant.

For the original multiplicity-nine example, take ell=6,p=3,m=9 and
q(z)=(z-1)^3. Then c=q^3, h=H z^2q^2,
u=-(H^2/6)z^3q and k=-(H^3/54)z^5. In the smaller ramification
r=T^6, zeta=1+A T+O(T^2), A^6=3/(2H), its chart is

    z=zeta+T^3 X=zeta+r^(1/2)X,
    phi=108 A^(-21)X^3+O(T),
    eta=[A^(-9)+O(T)]dX.

This recovers the exact order-zero barrier while placing it in the full
sharp family.

### Every displayed phase fails purity

Exact circular and ordinary Gaussian contraction gives, for the entire
family,

    E[P]=0,    E[P^2]=2 binom(2m,2)!=0.

The only surviving square term is (WcY)^2, whose expectation is

    2! [z^2] (z-1)^(2m)=2 binom(2m,2).

All other square and cross terms have odd Y degree or too much positive
z-minus-W weight: h contains z^2, u contains z^3 and k contains z^5.
The first moment vanishes for the same reason. At m=9 the second moment
is 306, independently of H.

These are exact obstructions to universal strict positivity, not pure
phases. Additional consequences of purity or a way to select the actual
Gaussian period could still exclude the local strip; the polynomial
shape and the Gaussian-square coefficient relations alone cannot do so
through a strict-positive-volume argument.

## Verification scope and remaining question

The all-degree classifications have precisely the coefficient and root
hypotheses stated in their sections. The affine-chart proof uses the
formal residue calculation and rational identities already present in
the boundary manuscript's cubic algebraic subcase; it adds the complete
all-degree classification, marker bounds and the exclusion of every
unsafe member from every nonisotropic constant-fibre frame. No novelty
claim against the external literature is made.

The companion [exact checker](verify_gaussian_affine_isotropic_chart.py)
checks symbolic identities for the arbitrary-ell branch, the shear,
the mixed-source bridges, the simple-puncture partial fractions,
the reciprocal double-root phase and the polynomial-orbit coefficients,
as well as finite Wick regression controls for the radial quartic formula.
It also checks the variable-precision determinant and source identities,
the sharp four-level common-root factor and scale identities, and
the representative cubic chart with its nonzero second moment,
and the marked coefficient identities in the all-degree residual-source
obstruction. The local common-root proof has symbolic controls for its
phase groups and exact scale inequalities.
For the six-parameter quartic family it verifies the exact finite
completeness certificate, including the nonzero modular Sylvester
determinant and the final cross-product equations. The proof explains
why this particular finite system forces all-order safety. The other
finite controls do not prove completeness or the all-degree theorems;
the preceding arguments supply those proofs. A run is recorded in
[the verification receipt](gaussian_affine_isotropic_chart_verification.json).

The additional mixed-source theorems use the explicitly identified
one-base leading-coefficient obstruction or the proved fixed-puncture
extension of the logarithmic seed bound, with the displayed spatial
and time bridges. The all-degree radial proof uses the explicitly
constructed relative endpoint module and its complete relative period
pairing. Normalization, boundary and phase identities have exact controls;
those checks do not certify the imported stationary-ball or marked-model
arguments. The radial quartic and polynomial-orbit obstructions also
have the elementary all-order proofs displayed above.

General H(2,1) remains open. In the W-affine mixed family, the remaining
common-root orders are precisely the strip

    m=ord_alpha c>=9,    ceil(2m/3)<=ord_alpha h<=m-3,

at a nonzero common root, with a=0 and outside the direct zero-root case.
All other root-order pairs are classified without a degree bound.
The strip is sharp for a uniform strictly-positive-volume proof: the
polynomial cubic barriers realize order zero or negative order at every
pair there, while failing purity. Additional consequences of purity or
control of the actual Gaussian period are still needed in this strip.

For polynomial ell(Z)W^2, the exclusion above leaves degree-one ell,
roots of multiplicity at least three in gcd(ell,c^2), and the fixed or
moving source exceptions. In the moving case the explicit failure
condition is R=0 and requires ell/gcd(ell,c^2) to be constant or a
scalar multiple of z. The monomial residual family in the uniform
obstruction above is excluded independently, even for high multiplicity.
The remaining cases also include the unclean zero-root source.
General W degree at least three in g, or
at least two in B, also remains open beyond the exact six-parameter
quartic classification. All statements continue to require the ordinary
nonisotropic Gaussian fibre Y. The scoped results and the method
obstruction do not settle arbitrary complex B(Z,W),g(Z,W).
