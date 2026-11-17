# All-degree Gaussian classifications in affine isotropic charts

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

## A radial quartic obstruction and weakly one-sided support

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

This implies a broader support obstruction. Suppose all monomials of
an H(2,1) phase have nonnegative Z-minus-W weight, and its weight-zero
part has total degree at most four. Positive-weight terms may have
arbitrary degree. Write P=P0+Pplus by angular weight. In a pure moment,
any factor from Pplus leaves strictly positive weight and cannot
contract. Therefore Phi_1(P)=Phi_1(P0). The weight-zero part has exactly
the radial form above: the constant fibre coefficient forbids RY^2.
Purity forces P0=0. Conversely this condition gives strictly positive
weight, so P is pure and safe with cutoff m>deg Q. The same conclusion
holds with all weights nonpositive, after exchanging Z and W.

In particular every total-degree-at-most-four H(2,1) phase with weakly
one-sided angular support is pure exactly when its weight-zero part
vanishes. A quartic counterexample must have monomials of both angular
signs in every allowable circular base frame. This does not classify
the remaining mixed-sign quartics.

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

## Verification scope and remaining question

The proof above is unrestricted in all polynomial degrees. It uses the
formal residue calculation and rational identities already present in
the boundary manuscript's cubic algebraic subcase; it adds the complete
all-degree classification, marker bounds and the exclusion of every
unsafe member from every nonisotropic constant-fibre frame. No novelty
claim against the external literature is made.

The companion [exact checker](verify_gaussian_affine_isotropic_chart.py)
checks symbolic identities for the arbitrary-ell branch and for the shear,
the directional-derivative obstruction, and finite Wick regression controls.
Those finite controls do not prove completeness or the all-degree theorem;
the preceding argument supplies that proof. A run is recorded in
[the verification receipt](gaussian_affine_isotropic_chart_verification.json).

The additional mixed-source theorems use the explicitly identified
one-base leading-coefficient obstruction or the proved fixed-puncture
extension of the logarithmic seed bound, with the displayed spatial
and time bridges. Their normalization and phase identities are checked
symbolically; those checks do not certify the imported stationary-ball
or marked-model arguments. The radial quartic and polynomial-orbit
obstructions have the elementary all-order proofs displayed above.

General H(2,1) remains open. The next cases beyond this result include
B=c(Z)W+u(Z) where deg c>=3 and c has repeated roots outside the zero-root case,
and higher W degree in B or g. A repeated root gives a higher-order
pole in dz/c; the fixed logarithmic bound no longer applies. Its
extension requires a bound retaining the special Gaussian-square
structure of the rational phase, or another argument. Arbitrary radial
B(R),g(R) beyond quartic degree also remain outside the radial proof:
the factorial functional has an endpoint, so the ordinary Gaussian
leading-coefficient obstruction cannot simply be substituted for it.
