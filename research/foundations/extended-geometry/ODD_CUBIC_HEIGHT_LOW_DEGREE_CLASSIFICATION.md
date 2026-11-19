# Odd cubic height: rational reconstruction and the sharp unsafe degree-six threshold

Date: 19 November 2026.

All fields below have characteristic zero. This result concerns
winding-one odd cubic-height phases. It does not classify arbitrary
quintics or septics.

Let `S=K[x,y,t]/(t²+2xy-1)` and use the exact polynomial functional

\[
\Lambda(f)=\operatorname{CT}_x\frac12\int_{-1}^1
 f\left(x,\frac{1-t^2}{2x},t\right)dt.
\]

Consider every regular polynomial phase of the form

\[
p=\frac{t}{x}(C(x)+t^2D(x)),\qquad C,D\in K[x].                (OC1)
\]

Regularity is exactly `C(0)+D(0)=0`. Call p pure if all positive pure
moments vanish, and safe if every fixed polynomial marker has eventual
mixed vanishing. Least degree means the polynomial filtration of S,
equivalently the radial-primitive balanced degree of its homogeneous
two-pair lift.

**Theorem OddCubicHeightLowDegree.** Every pure phase(OC1) of least
degree at most5 is safe. More precisely:

* If `C(0)=0`, purity is equivalent to both C,D being divisible by x².
  This is the strictly positive angular-weight safe branch, in every degree.
* If `C(0)=c0!=0` and the least degree is at most5, purity forces
  `C=c0,D=-c0`, so `p=2c0*t*y`, a strictly negative angular-weight safe
  phase.

The retained odd-height family, Theorem5.5 of
[TWO_PAIR_IMAGE_MATHIEU_COUNTEREXAMPLE.md](TWO_PAIR_IMAGE_MATHIEU_COUNTEREXAMPLE.md), has a radial-primitive
non-proper-power unsafe k=1 member of least degree6, with its original
fixed marker surviving at every even positive order. Combining that
antecedent with the low-degree obstruction proves that6 is the sharp
least unsafe degree IN THIS CLASS. No new degree6 construction is claimed.
In particular there is no primitive degree-five Image counterexample in
the entire class. Degree-seven phases and all other quintic shapes remain
outside the conclusion.

The proof includes an independent all-degree structural result.

**Theorem OddCubicRationalCenter.** If `C(0)!=0` and p is pure, put
`X=sqrt(-D/C)` with `X(0)=1`, `A=C/X`, and `u=x/A`. Then u and A are
rational functions in K(x), even though their initial definitions are
formal algebraic series. This rationality is not an extra ansatz.

## 1. The exact formal even-moment criterion

Put `L=1/X`. Then `C=A/L,D=-A/L³`, and under `v=Xt`,

\[
p=(A/x)v(1-v^2),\qquad dt=L\,dv.
\]

Odd pure moments vanish by t-parity. For even m let

\[
\kappa_m=\frac12\int_{-1}^1 v^m(1-v^2)^m\,dv\ne0.
\]

It is a nonzero rational number: over the reals its integrand is
nonnegative and nonzero. Characteristic zero transfers that nonvanishing
without an ordered-field hypothesis on K. A polynomial primitive of
`v^m(1-v²)^m` has order m+1 contact at each endpoint ±1. Since
`X=1+O(x)`, endpoint substitution gives

\[
\Lambda(p^m)=\kappa_m[x^m]A^m L \quad(m\text{ even}).           (OC2)
\]

The same endpoint congruence is valid below coefficient m as well, hence
for every nonnegative integer j,

\[
\Lambda(x^j p^m)=\kappa_m[x^{m-j}]A^m L
\quad(m\text{ even}).                                      (OC3)
\]

These identities hold over the formal coefficient ring: they are not
deduced from a finite moment prefix.

Let `v(u)` denote the formal inverse of `u=x/A`. Formal residue
substitution turns the necessary diagonal conditions into

\[
F(u)=L(v(u))\frac{u v'(u)}{v(u)}=1+uH(u^2).                  (OC4)
\]

Indeed `[u^m]F=[x^m]A^mL` and `F(0)=1`; purity removes exactly the
positive even coefficients. Formal change of variable follows first on
monomials: derivatives have residue zero, while `u'/u` has residue1.
It therefore holds coefficientwise for the Laurent series used here.

## 2. A Galois character projector forces u to be rational

The finite-group step is useful separately.

**Lemma CharacterSpanDescent.** Let E/F be a finite separable extension
in characteristic zero, let L/F be a normal closure, and write
`G=Gal(L/F), H=Gal(L/E)`. Let a belong to E. If the F-linear span of
the conjugates of a contains a nonzero vector b transforming by a
one-dimensional character `chi:G->F*`, then chi is trivial on H and
b belongs to E.

Indeed the projector `P_chi=|G|^-1 sum_g chi(g)^-1 g` fixes b. If an
h in H had chi(h)!=1, then `h(a)=a` and equivariance force
`P_chi(a)=chi(h)P_chi(a)`, hence this projector annihilates a and every
conjugate. It would annihilate b too, a contradiction. Thus H fixes
every chi-eigenvector, including b, proving the lemma.

**Corollary QuadraticInverseBranchDescent.** Let R,S belong to K(x),
with R nonconstant. If two inverse branches of R over z have S-values
differing by `1/sqrt(z)`, then sqrt(R(x)) belongs to K(x). The difference
places sqrt(z) in the normal closure; apply CharacterSpanDescent with
a=S(x), b=1/sqrt(z) and its quadratic character. Equality of the stated
formal local algebraic branches is enough.

Because `A=C L`, the function in(OC4) is `F(u)=v'(u)/C(v(u))`.
Let `V(x)` be the formal primitive of `1/C(x)` with V(0)=0. Integrating
(OC4) gives

\[
V(v(u))=u+E(u^2),\qquad E(z)\in zK[[z]].
\]

Now set

\[
R(x)=u(x)^2=-\frac{x^2D(x)}{C(x)^3}\in K(x),\qquad
S(x)=\frac1{C(x)R'(x)}\in K(x).
\]

Both are well defined: R has order2 at zero. Differentiating the last
formal identity gives

\[
S(x)=\frac1{2u(x)}+E'(R(x)).                                (OC5)
\]

The two inverse branches `x_+(u)=v(u)` and `x_-(u)=v(-u)` of `R(x)=u²`
therefore satisfy the EXACT algebraic-function identity

\[
S(x_+(u))-S(x_-(u))=1/u.                                   (OC6)
\]

Here an equality of formal expansions of algebraic functions is an
identity of those functions; no finite approximation is used.

Take a normal closure L of the finite separable extension `K(x)/K(z)`,
where `z=R(x)`, and write `G=Gal(L/K(z))`, `H=Gal(L/K(x))` for the
chosen embedding. Both inverse-branch values in(OC6) are conjugates of
S(x) in L, so their difference puts `sqrt(z)` in L. Let
`chi:G->{+1,-1}` be its quadratic character. Define the character
projector

\[
P_\chi(a)=\frac1{|G|}\sum_{g\in G}\chi(g)g(a).
\]

If some h in H had `chi(h)=-1`, then `h(S)=S` would imply
`P_chi(S)=P_chi(hS)=chi(h)P_chi(S)=0`. Equivariance of this projector
would annihilate every G-conjugate of S, and therefore their difference.
But `1/sqrt(z)` is itself a nonzero chi-eigenvector, so its projector is
itself. This contradicts(OC6). Thus chi is trivial on H, and

\[
\sqrt z\in L^H=K(x).
\]

The selected square root with leading term `x/C(0)` is precisely the
formal u. Therefore `u,A=x/u` lie in K(x). This argument already gives
descent to K, rather than only rationality after adjoining algebraic
constants. This proves OddCubicRationalCenter.

Rationality of u does NOT imply that F in(OC4) is a rational function of
u. An argument that assumes rational F would lose algebraic branches;
none of the low-degree elimination below makes that assumption.

## 3. Polynomial factorization and exhaustive degree branches

Scale the phase to set `c0=1`; scalar multiplication preserves all zero
moments and all eventual-vanishing properties. Rationality gives
`-C³/D=A²` in K(x). Over an algebraic closure, the multiplicities of
every irreducible factor of C and D consequently have the same parity.
Thus there are polynomials

\[
C=H J^2,\qquad D=-H K^2,\qquad H(0)=J(0)=K(0)=1,            (OC7)
\]

with H squarefree. H,J,K need not be pairwise coprime. It is harmless
that this factorization is made over an algebraic closure: impossibility
there implies impossibility over the original K, and the final constant
phase descends immediately.

Writing `h=deg H,j=deg J,k=deg K`, the exact least degree of p is

\[
d=\max(2,h+2j,h+2k+2).                                     (OC8)
\]

To see exactness, rewrite the regular phase as
`2ty+t sum_(i>=1)c_i x^(i-1)+t³ sum_(i>=1)d_i x^(i-1)`.
For d>=3 its nonzero top part lies in K[x,t], so it is not divisible by
the leading quadric `t²+2xy`; degree2 has a nonzero ty term. The sphere
relation cannot reduce that degree.

The hypothesis d<=5 leaves exactly these degree boxes:

| h | j | k |
|---|---|---|
|0|<=2|<=1|
|1|<=2|<=1|
|2 or3|<=1|0|

The normalized even-moment equations from(OC2) are

\[
f_m=[x^m]H^m J^{3m+1}K^{-m-1}=0\quad(m=2,4,6,8,\ldots).     (OC9)
\]

All coefficient identities below are over Q and therefore valid over
every characteristic-zero field. Only these necessary equations are
used after the all-order rationality theorem has supplied(OC7).

## 4. Nonconstant H: two exact unit-ideal certificates and one short calculation

If h=1, scale the angular coordinate to make `H=1+x`. Write
`J=1+b x+c x²,K=1+d x`. Then

\[
f_2=21b^2-21bd+14b+7c+6d^2-6d+1,
\]

so

\[
c_*=-3b^2+3bd-2b-6d^2/7+6d/7-1/7.
\]

After this substitution, the THREE equations `f4=f6=f8=0` generate
the unit ideal in Q[b,d]. The exact rational Bezout identity is supplied
in `odd_cubic_height_low_degree_certificates.json`, key `h1`.
It has degree bound16 and multipliers with71,38,24 terms. The file
specifies each primitive generator, its nonzero rational scaling relative
to the substituted f_m, and every multiplier coefficient. Direct
polynomial multiplication gives exactly1, so the equations have no
simultaneous solution over any characteristic-zero extension.

For h=2 or3 and j=1, scale to make `J=1+x`; write
`H=1+b x+c x²+d x³,K=1`. This includes degree2 by setting d=0.
Here

\[
f_2=b^2+14b+2c+21,
\quad c_*=-b^2/2-7b-21/2.
\]

Again the substituted equations `f4=f6=f8=0` generate the unit ideal.
The exact degree12 Bezout identity is the same file's key `h3_j1`,
with22,11,8 multiplier terms. Thus these branches are impossible.

For h=2 or3 and j=0, `J=K=1`; put `H=1+a x+b x²+c x³`.
The equations are simply `[x^m]H^m=0`. The second gives `b=-a²/2`,
and the fourth then is `a(24c-7a³)/2=0`. If a=0, the sixth is
`15c²=0`, forcing c=b=0. If a!=0, put `c=7a³/24`; the sixth becomes
`1397a^6/192=0`, impossible. Thus H is constant, contradicting h=2 or3.
This exhausts all nonconstant H cases.

## 5. Constant H: an explicit homogeneous two-variable elimination

Now H=1, `J=1+b x+c x²,K=1+d x`. Equation f2 gives

\[
c_*=-3b^2+3bd-6d^2/7.
\]

After substitution, define g4=f4 and g6=f6. They are

\[
g_4=-\frac{56693b^4-101920b^3d+67977b^2d^2-19474bd^3+1952d^4}{49},
\]

\[
g_6=\frac3{343}(4099193b^6-11965212b^5d+14438879b^4d^2
 -9211314b^3d^3+3271002b^2d^4-611800bd^5+47048d^6).
\]

If d=0, `g4=-1157b^4`, so b=0 and then c=0. If d!=0, homogeneity
allows d=1. The polynomials `g4(b,1),g6(b,1)` have a rational Bezout
identity `A(b)g4+B(b)g6=1`, retained in the certificate file under
`h0_d_nonzero`. Its multiplier degrees are5 and3. This excludes d!=0.
Therefore J=K=1. Restoring c0 gives exactly `C=c0,D=-c0`.

Angular coordinate scaling here is a genuine symmetry: substitute x->s x,
y->y/s in the sphere and multiply the phase by s to preserve C(0)=1.
Thus the normalizations of H, J or d used above lose no nonzero branch.

## 6. The zero-center branch

If C(0)=D(0)=0, p has only nonnegative angular powers. Its weight-zero
part is `ct+dt³`, with `(c,d)=(C'(0),D'(0))`. Every pure moment equals
the corresponding interval moment of this part. If d=0, the second
moment is c²/3 and forces c=0. If d!=0, normalize d=1 and set r=c/d.
The second and fourth moments are

\[
P_2(r)=r^2/3+2r/5+1/7,
\quad P_4(r)=r^4/5+4r^3/7+2r^2/3+4r/11+1/13.
\]

They have resultant

\[
7663616/86995605571875\ne0.
\]

An exact rational Bezout identity equal to1 is included in the
certificate file under `zero_center_core`; direct multiplication, rather
than a numerical root test, verifies it. Thus d!=0 is impossible too.
Consequently C,D are both divisible by x². Conversely that condition
makes every term of p strictly positive in angular x-weight, so all pure
moments and all sufficiently high mixed moments with any fixed marker
vanish. This proves the zero-center assertion in every degree.

The unit-center constant phase `2c0ty` is similarly safe: every positive
power has angular weight -m, and a fixed polynomial marker has bounded
positive angular support. Hence sufficiently high powers have no constant
term. This finishes OddCubicHeightLowDegree.

For completeness, angular safety transfers to every fixed original
polynomial Image marker, not only balanced markers. Split a marker into
bihomogeneous components Q of bidegree(a,b). If a>b, its contraction
against a balanced phase power is zero by degree. If a<=b, that output
is homogeneous in z of degree b-a. For each multi-index gamma of that
size, differentiation of the monomial contraction identity gives

\[
\partial_z^\gamma\mathcal E_2(QF^m)
       =\mathcal E_2(\xi^\gamma QF^m).
\]

The marker on the right is balanced, hence is a homogeneous polynomial
in the four bilinears R,Z,W,T. The factorial/angular normalization and
the safe sphere branch annihilate it eventually. There are finitely
many gamma and bihomogeneous components, so all top derivatives of
each homogeneous output eventually vanish, hence the output itself
vanishes. This proves safety for arbitrary original markers.

## 7. The retained k=1 witness makes the threshold sharp

Recall the existing odd-height family in
[TWO_PAIR_IMAGE_MATHIEU_COUNTEREXAMPLE.md](TWO_PAIR_IMAGE_MATHIEU_COUNTEREXAMPLE.md),
Theorem5.5 and equations(5.11)--(5.12). Its k=1 member is precisely the
degree6 phase used here. The following formulas are recalled for the
sharpness implication and independent finite controls; they are not a
new witness construction.

Take `C=(1+x)²,D=-(1+x)^4`. Then `X=1+x,A=1+x,L=1/(1+x)`, so
`A^mL=(1+x)^(m-1)`. Its coefficient at x^m is zero, while its coefficient
at x^(m-1) is1. Equations(OC2)--(OC3) therefore prove all-order purity
and fixed-marker failure, with odd moments zero by parity.

For the standard two-pair bilinears

\[
R=\xi_1z_1+\xi_2z_2,\quad Z=\xi_1z_2,\quad W=2\xi_2z_1,
\quad T=\xi_1z_1-\xi_2z_2,\quad T^2+2ZW=R^2,
\]

the explicit balanced-degree6 lift is

\[
F_6^{odd}=T(R+Z)^2\big(2R^2W-T^2(2R+Z)\big).
\]

The exact factorial normalization gives

\[
\mathcal E_2((F_6^{odd})^m)=0,
\qquad
\mathcal E_2(Z(F_6^{odd})^m)=(6m+2)!\kappa_m
\quad(m\text{ even}),
\]

where `kappa2=8/105,kappa4=128/15015` and every even kappa is nonzero.
This establishes failure of eventual mixed vanishing with the original
fixed marker Z. At `(xi1,xi2,z1,z2)=(1,1,1,-1)`, R=0 and F=8, so the
phase is radial-primitive. Its exact multiplicity along the absolutely
irreducible rank-two bilinear T is1: after dividing by T, the value at
the T=0 point `(1,1,1,1)` is144. Hence it is not a proper power, even up
to a scalar over extensions. All those witness facts are retained
antecedents; the new threshold is scoped to class(OC1).

## 8. Exact computational certificate and remaining scope

Verification package:

* `odd_cubic_height_low_degree_certificates.json`: all four exact
  unit-ideal/univariate Bezout certificates,115544bytes.
* `verify_odd_cubic_height_low_degree.py`: standard-library verifier
  independently derives the moment polynomials, eliminates c, checks each
  primitive scaling and the identities equal to1, and rejects an altered
  multiplier. No SymPy or Groebner implementation is imported for replay.
* `odd_cubic_height_low_degree_verification.json`: passed replay,
  including48 independent raw interval checks and original four-variable
  contractions of the sharp degree6 witness through m8.

The certificates are exact finite implications inside a complete algebraic
case split. The rational-center bridge uses ALL even moments. No claim
that the first eight original moments force rational reconstruction, or
form a universal purity cutoff, follows from this proof.

Primary composition results were inspected before the prime-degree work:

* Pakovich--Muzychuk, *Solution of the polynomial moment problem*, Proc.
  London Math. Soc.99(2009),633--657, actual Theorem1.1, author PDF p3:
  <https://www.math.bgu.ac.il/~pakovich/Publications/s.pdf>.
  This concerns a FIXED univariate polynomial P and interval moments
  `integral P^i q`, classified as sums of compositional solutions.
* Pakovich, *On the polynomial moment problem*, Math. Res. Lett.10(2003),
  401--410, actual Theorem1 p402:
  <https://www.math.bgu.ac.il/~pakovich/Publications/opmp.pdf>.
  For indecomposable P under its endpoint hypotheses, Q factors through P.

These results are important antecedents for inverse-branch and composition
methods. They do not identify a prime balanced degree with the degree of
that fixed P, and they are not imported to infer a global prime-degree
obstruction. The Galois projector and all low-degree algebraic implications
needed here are supplied explicitly. No worldwide priority, formal proof
assistant verification, or external human review is asserted.


The [source ledger](odd_cubic_height_low_degree_source_ledger.json) pins
the retained normalization and sharpness input and records the exact
primary-source antecedent scopes. The [checker](verify_odd_cubic_height_low_degree.py)
and [receipt](odd_cubic_height_low_degree_verification.json) separate
the exact finite algebraic implications from the unrestricted Galois
reconstruction.

```bash
python3 research/foundations/extended-geometry/verify_odd_cubic_height_low_degree.py \
  --write-receipt research/foundations/extended-geometry/odd_cubic_height_low_degree_verification.json
```
