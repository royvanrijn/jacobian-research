# Complete fibres of Davenport–Zannier endpoint coordinates

Fix a characteristic-zero field \(K\) and a rational function
\(u\in K(x)\) with \(u(0)=0\), \(u'(0)=1\). Assume that

\[
P=\frac{u}{xu'}\in K[x].
\tag{S1}
\]

This finite algebraic condition characterizes the endpoint coordinates
that admit a polynomial constant-endpoint profile at some starting
height. The nontrivial coordinates are exactly the marked
Davenport–Zannier coordinates described below. We classify their
**entire polynomial pure fibres**, including nonconstant endpoint
profiles, at every integer starting height. There is no coefficient
degree bound.

The key conclusion is stronger than an infinite-fibre criterion:
weighted Davenport–Zannier trees permit **no** nonconstant-endpoint
deformation at their fixed coordinate. Their admissible fibres are
singletons. Ordinary trees permit exactly one explicit sequence of
deformations, and that sequence is exhaustive.

The theorem concerns the exact two-height sphere phase

\[
p=t^a\frac{C(x)+t^2D(x)}x,
\qquad C,D\in K[x],\quad C(0)=-D(0)=1,
\quad t^2+2xy=1,
\tag{S2}
\]

with fixed endpoint coordinate \(x(-D/C)^{a/2}/C=u\); the selected
square root has value 1 at zero. Purity means that every positive pure
angular moment vanishes. The statement does not classify coordinates
for which (S1) fails, arbitrary sphere phases, or unrestricted
primitive degrees 5 and 7.

## 1. The all-height singleton/infinite theorem

If \(P=1\), then \(u=x\), and its entire fibre at every integer
\(a\geq0\) is \(C=1,D=-1\).

Suppose henceforth that \(r=\deg P\geq1\). Then \(P\) is squarefree,
\(u(\infty)=c\in K^*\), and there are coprime monic polynomials of the
same degree \(N\) such that

\[
\frac uc=\frac AB,\qquad A(0)=0\text{ simply},\quad B(0)\ne0,
\qquad \deg(A-B)=N-r.
\tag{S3}
\]

Let \(h\) be the largest multiplicity of an unmarked geometric root of
\(A\), with \(h=0\) if there is no such root. In particular \(N\geq r\).

**Theorem DavenportZannierCoordinateFibres.** Fix any integer
\(a\geq0\). The full polynomial pure fibre at \(u\) is as follows.

| Condition | Entire fibre |
|---|---|
| \(a<h\) | Empty |
| \(a\geq h\), \(a\) even | The one constant-endpoint pair below |
| \(a\geq h\), \(a\) odd, \(N>r\) | The same singleton |
| \(a\geq h\), \(a\) odd, \(N=r\) | Exactly the sequence (S5), \(k\geq0\) |

The constant-endpoint pair is

\[
C_0=\frac{xP^a}{u},\qquad D_0=-C_0P^2.
\tag{S4}
\]

In the last row define

\[
R=\frac{c+u}{c-u}=\frac{A+B}{B-A}\in K[x],\qquad
L_k=\frac{1+R^k}{2},
\qquad
C_k=C_0L_k^a,\quad D_k=D_0L_k^{a+2},\quad k\geq0.
\tag{S5}
\]

Here \(R(0)=1\), \(\deg R=r=N\), and these are distinct pairs.
The selected square root is \(X_k=PL_k\), so all divisions in (S4)
are already polynomial before deformation. The exact least sphere
degrees are

\[
\deg_S p_0=(a+2)(r+1),\qquad
\deg_S p_k=(a+2)\bigl(r(k+1)+1\bigr).
\tag{S6}
\]

Only \(k=0\) has constant endpoint. Thus the finite nonempty fibres in
this whole coordinate class are **singletons**, not unspecified finite
sets. All statements hold over the original field and after arbitrary
characteristic-zero extensions.

## 2. The logarithmic geometry of the coordinate

Equation (S1) gives

\[
\frac{du}{u}=\frac{dx}{xP(x)}.
\tag{S7}
\]

A logarithmic differential of a nonzero rational function has a simple
pole at each zero and pole, with residue its nonzero integral order.
Consequently the roots of \(P\) are simple and are exactly the finite
zeros and poles of \(u\) other than zero. The normalization gives
\(P(0)=1\).

If \(P\) is constant, it is 1, and (S7) gives \(u=x\). If \(P\) is
nonconstant, \(u\) cannot have a zero or pole at infinity: such an
order would give \(u'/u\) a nonzero \(x^{-1}\) leading term, whereas
the right side of (S7) has order \(r+1\) there. Thus
\(u(\infty)=c\ne0,\infty\). Comparing the first nonconstant term of
\(u/c\) at infinity gives

\[
\operatorname{ord}_\infty(u-c)=r.
\tag{S8}
\]

Reducing \(u/c=A/B\) with monic polynomials gives (S3). In particular
\(u^{-1}(c)=\{\infty\}\) exactly when \(N=r\). As a differential,
(S7) has no finite zeros and has its only possible zero at infinity,
of order \(r-1\). That last fact supplies the rigidity.

Conversely the retained marked Davenport–Zannier theorem proves that
(S3), with \(r\) the total number of distinct roots of \(AB\) minus 1,
gives (S7). We use that correspondence and its tree interpretation;
the new fibre argument starts from (S7).

## 3. No larger intermediate field supports the marked reflection

Let \(v(s)\) be the inverse branch of \(u\) at zero, so rational
functions of \(x\) embed in \(K((s))\) by \(x=v(s)\).

**Lemma LogarithmicStableFieldRigidity.** Assume (S1) with
\(r=\deg P\geq1\), as in Section 2. Suppose
\(K(u)\subseteq E\subseteq K(x)\) and replacement \(s\mapsto-s\)
preserves \(E\). Then \(E=K(u)\).

By Lüroth, \(E=K(w)\); write \(u=U(w)\). The formal replacement
induces a genuine involution \(\sigma\) of this intermediate rational
field, with \(\sigma(U)=-U\). In particular the rational differential

\[
\theta=\frac{dU}{U}
\quad\text{satisfies}\quad \sigma^*\theta=\theta,
\qquad w^*\theta=\frac{du}{u}.
\tag{S9}
\]

We may inspect its geometric orders over an algebraic closure; the
rational maps and the degree conclusion remain over \(K\). At a
point of ramification index \(e\), a differential of order \(m\)
pulls back to order

\[
e(m+1)-1.
\tag{S10}
\]

Indeed in local parameters the map is \(z=\lambda t^e+\cdots\),
\(\lambda\ne0\), and \(z^m dz\) has the stated leading order.
Characteristic zero makes \(e\lambda\) nonzero.

If \(\theta\) had a zero at a point \(Q\), every point over \(Q\)
would be a zero of \(du/u\). By Section 2, all such points would have
to be \(x=\infty\). Therefore \(\theta\) could have only one zero,
namely \(Q=w(\infty)\). Invariance in (S9) would make \(\sigma(Q)=Q\).
But

\[
U(Q)=u(\infty)=c\in K^*,\qquad
U(\sigma Q)=-U(Q),
\]

which would give \(c=-c\). This is impossible. Hence \(\theta\) has
no zeros.

A rational differential on \(\mathbb P^1\) has divisor degree \(-2\):
the divisor of its rational coefficient has degree zero and \(dw\)
has a double pole at infinity. All poles of \(dU/U\) are simple and
are exactly the distinct zeros and poles of \(U\). With no zeros,
there must be exactly two such points: one zero and one pole. Their
multiplicities in \(U\) are the same integer \(m=\deg U\).
The marked zero of \(u\) is simple, so the pullback multiplicity
\(m e\) there equals 1. Thus \(m=1\), and \(E=K(U)=K(u)\).

This proves the lemma without requiring a rational return on the
original \(x\)-line, a polynomial generator \(w\), or a classification
of rational-map decompositions.

## 4. Exhaustion of every odd-height fibre

Let \(a\) be odd and suppose (S2) is pure at the fixed coordinate.
Set \(X=\sqrt{-D/C}\), \(X(0)=1\). It is rational: since \(u\) is
rational and \(a\) is odd,
\(X=(uC/x)/(X^2)^{(a-1)/2}\in K(x)\).

The exact endpoint criterion in the retained all-height theorem is

\[
T(v(s))+T(v(-s))=2,\qquad
T=\frac{X^{a-1}}{Cu'}=\frac PX.
\tag{S11}
\]

Thus \(\eta=2/T-1=2X/P-1\) is normalized at zero and satisfies
\(\eta(v(s))\eta(v(-s))=1\). The field \(K(u,\eta)\) is preserved
by formal reflection, which sends \((u,\eta)\mapsto(-u,\eta^{-1})\).
The lemma therefore gives

\[
\eta=H(u),\qquad H(z)H(-z)=1,\qquad H(0)=1
\tag{S12}
\]

for a rational function \(H\in K(z)\).

Polynomiality of \(C=xX^a/u\) and \(D=-xX^{a+2}/u\) ensures that
every finite pole of \(\eta\) is a pole of \(u\). At a finite nonzero
zero of \(u\) of order \(s\), polynomiality requires
\(\operatorname{ord}X\geq\lceil s/a\rceil\geq1\); \(P\) has order
one there, so \(\eta\) is regular. At any other finite point that is
not a pole of \(u\), \(X\) and \(\eta\) are regular as well.
The same conclusion is the finite-pole part of the retained general
fixed-coordinate theorem.

If \(H\) has a pole at a finite value \(\beta\), every finite point
of \(u^{-1}(\beta)\) would then be a forbidden pole of \(\eta\).
Such a fibre can be entirely at infinity only if \(\beta=c\), and
only if \(N=r\). There can be no pole of \(H\) at infinity either:
the identity in (S12) gives
\(2\operatorname{ord}_\infty H=0\).

If \(N>r\), \(H\) has no poles anywhere and is constant. Its marked
value is 1. Thus \(\eta=1\) and \(X=P\).

If \(N=r\), the only possible pole of \(H\) is \(c\). Reflection in
(S12) puts every zero at \(-c\), and the normalization gives exactly

\[
H(z)=\left(\frac{c+z}{c-z}\right)^k,
\qquad k\geq0.
\tag{S13}
\]

Negative \(k\) would put a forbidden pole at \(-c\). The formula
also includes \(H=1\). Consequently (S5) exhausts all possibilities.
Conversely every (S5) has the same \(u\), satisfies (S11), and is pure
when its base pair is polynomial.

At every unmarked zero of \(u\), (S12) gives \(\eta=1\). Hence
\(X/P=(1+\eta)/2=1\) there, with no additional zero. The exact local
order of \(C\) is \(a-s\). Thus no deformation repairs a height
\(a<h\). At a pole of \(u\) of order \(m\), the base pair has orders
\(a+m\) and \(a+2+m\), both nonnegative; at the marked zero the
factor \(x\) cancels its simple zero. Therefore the base pair is
polynomial exactly when \(a\geq h\). This proves the odd rows of the
theorem, including the empty case.

For an even integer \(a\geq0\), the endpoint identity detects every
positive exponent and forces \(T=1\), not only its even part. Thus
\(X=P\) directly, and the same local orders give precisely the even
row. This includes \(a=0\), without requiring \(P^{-1}\) to be
polynomial. The trivial \(u=x\) case follows from this observation
for even heights and from its retained exact singleton fibre for odd
heights.

Finally \(\deg C_0=ar+1\), \(\deg D_0=(a+2)r+1\). For \(k\geq1\),
\(\deg L_k=kr\); for \(k=0\), \(L_0=1\). Substitution into the exact
least-degree formula
\(\max(\deg C+a-1,\deg D+a+1)\) proves (S6).

## 5. Shabat polynomials and ordinary trees

**Corollary ShabatFibreParameterization.** At a fixed odd height
\(a\geq1\), all infinite fibres satisfying (S1) are in bijection with
polynomials \(R\in K[x]\) such that

\[
R(0)=1,\qquad R'(0)\ne0,
\tag{S14}
\]

all finite critical values of \(R\) lie in \(\{-1,1\}\), and every
root of \(R-1\) other than the marked root zero has multiplicity at
most \(a\). Degree-one polynomials are included. The formulas are

\[
u=\frac{2}{R'(0)}\frac{R-1}{R+1},\qquad
P=\frac{R^2-1}{2xR'},
\tag{S15}
\]

followed by (S4)--(S5). These are actual normalized polynomials over
\(K\), not just geometric equivalence classes.

Indeed (S5) gives (S15). Since \(du/u\) has no finite zero, a critical
point of \(R\) can only lie above \(1\) or \(-1\). Conversely that
critical-value condition says \(R'\mid R^2-1\); the quotient has a
simple zero at each distinct root of \(R^2-1\). The marked root
cancels the extra \(x\), and Taylor expansion at zero gives
\(P(0)=1\). Thus \(P\) is squarefree of degree \(\deg R\), and (S15)
satisfies (S1). The stated multiplicity cap is exactly the height
condition already proved.

In the classical weighted-plane-tree correspondence for a DZ pair,
the total edge weight is \(N\) and the number of edges is \(r\).
All weights are positive integers. Hence

\[
N=r\quad\Longleftrightarrow\quad
\text{every edge has weight }1.
\tag{S16}
\]

These are the ordinary trees, equivalently Shabat polynomials. A marked
black leaf records the simple zero of \(A\); the other black degrees
are at most \(a\). Thus (S16) gives the exact geometric distinction
between singleton and infinite odd-height fibres. It is not an
assertion that every geometric tree has a model over a specified
nonclosed field. The retained marked Galois descent criterion applies.

## 6. Exact counts of rigid and deformable coordinate classes

Fix an odd \(a\geq1\) and \(r\geq1\). Count geometric classes modulo
the specified source scaling \(x\mapsto sx\), \(s\ne0\). The retained
constant-endpoint count is

\[
A_{a,r}=\frac1r\binom{ar}{r-1}.
\tag{S17}
\]

**Theorem SingletonAndInfiniteFibreCounts.** The exact numbers of
infinite and singleton coordinate classes are respectively

\[
I_{a,r}=\frac1r[z^{r-1}](1+z+\cdots+z^a)^r,
\qquad A_{a,r}-I_{a,r}.
\tag{S18}
\]

An explicit integer formula is

\[
I_{a,r}=\frac1r\sum_{j=0}^{\lfloor(r-1)/(a+1)\rfloor}
(-1)^j\binom rj
\binom{2r-2-(a+1)j}{r-1}.
\tag{S19}
\]

To prove the new count, remove the marked black leaf and its edge.
Its white neighbor becomes a rooted white vertex with an indicated
corner. At each white vertex, the black branches are ordered. An
unweighted black branch with \(j-1\) white children has \(j\) edges
and is allowed exactly when \(1\leq j\leq a\). If \(T\) counts these
residual rooted trees by edges, then

\[
T=\frac1{1-\sum_{j=1}^a z^jT^{j-1}},\qquad
F=zT=z(1+F+\cdots+F^a).
\tag{S20}
\]

The removed edge makes the desired coefficient \([z^r]F\). Formal
Lagrange inversion gives exactly (S18); alternatively apply the
residue substitution proof already used for (S17). Expanding
\((1-z^{a+1})^r(1-z)^{-r}\) gives (S19). The marked edge eliminates
plane-tree automorphisms, so there is no further symmetry factor.

For example, at height \(a=3\):

| \(r\) | All coordinate classes | Infinite fibres | Singleton fibres |
|---:|---:|---:|---:|
| 1 | 1 | 1 | 0 |
| 2 | 3 | 1 | 2 |
| 3 | 12 | 2 | 10 |
| 4 | 55 | 5 | 50 |
| 5 | 273 | 13 | 260 |
| 6 | 1428 | 36 | 1392 |

For \(a=1\), both \(A_{1,r}\) and \(I_{1,r}\) equal 1: every
nontrivial coordinate in this class deforms. These are source-scaling
classes, not orbits under all sphere or original-variable automorphisms.

There is also an exact count after allowing **all** the deformations.
At exact least degree \(D=(a+2)(s+1)\), \(s\geq1\), the number of pure
phase classes satisfying (S1) is

\[
A_{a,s}+\sum_{\substack{r\mid s\\r<s}}I_{a,r}.
\tag{S21}
\]

No other nontrivial degrees occur in this coordinate class. A singleton
has \(r=s\); an infinite fibre contributes exactly once when
\(r(k+1)=s\). The \(r=s,k=0\) terms are already included in
\(A_{a,s}\), proving (S21) without double counting. The coordinate
is recovered from its phase, so distinct coordinate classes do not
acquire new identifications under this specified source scaling.

## 7. Examples, exclusions and finite controls

The star polynomial \(R=2(1+\lambda x)^r-1\), \(\lambda\ne0\),
has a marked simple \(R=1\) point and only the critical value \(-1\).
It is admissible at every height \(a\geq1\), and (S5) gives its full
odd-height fibre. At \(k=0\) it is the retained constant-endpoint
star family.

A non-star example is

\[
R=4(1+x)^3-3(1+x)=1+9x+12x^2+4x^3.
\tag{S22}
\]

Its unmarked \(R=1\) root has multiplicity 2. Thus its fibre is empty
at \(a=1\), is a singleton at even \(a\geq2\), and is infinite at
odd \(a\geq3\). At \(a=3\) its full sequence has degrees
\(20,35,50,\ldots\).

For a rigid weighted example put \(U=1+x\), \(V=1+2x/3\), and

\[
u=\frac{xV^2}{U^3},\qquad P=UV.
\tag{S23}
\]

Here \(r=2,N=3,h=2\). Its fibre is empty at heights 0 and 1, and
is a singleton at **every** height \(a\geq2\):
\(C=V^{a-2}U^{a+3}\), \(D=-V^aU^{a+5}\). The lack of odd-height
deformations is a conclusion of the stable-field lemma, not a failed
search over bounded coefficients.

The polynomial condition (S1) is essential. Put \(z=1+x\),
\(w=(z-1)/(z+1)\), and

\[
u=\frac23(w^3+3w)=\frac83\frac{z^3-1}{(z+1)^3}.
\tag{S24}
\]

Then \(u'(0)=1\), and the displayed rational function of \(z\)
changes sign under \(z\mapsto1/z\), but

\[
\frac{u}{xu'}=
\frac{(z^2+z+1)(z+1)}{3(z^2+1)}
\]

is not polynomial. Its only finite pole is the inversion-fixed point
\(z=-1\), so the permitted reciprocal units are \(z^k\), \(k\geq0\).
The two finite critical points \(z=\pm i\) give the exact congruence
\(k\equiv2\pmod4\). Thus this degree-three reciprocal quotient has
an infinite fibre at every odd height; it is outside the rigidity
theorem. Its first cubic-height example is a retained one-pole
fixed-zero example, not a new witness claim.

The checker also uses the degree-two quotient
\(u=-(z^2-1)/(z^2-4z+1)\). Its critical roots satisfy
\(z^2-z+1=0\), and \(\eta=z^{3+6j}\), \(j\geq0\), supplies polynomial
pure pairs. These controls prevent extrapolating the Shabat conclusion
to arbitrary rational endpoint coordinates.

## 8. Sources and assurance

The exact all-height endpoint criterion and least-degree interface are
retained in [ALL_START_HEIGHT_MARKER_SAFETY.md](ALL_START_HEIGHT_MARKER_SAFETY.md)
and [ALL_START_HEIGHT_SHARP_UNSAFE_DEGREE.md](ALL_START_HEIGHT_SHARP_UNSAFE_DEGREE.md).
The marked DZ correspondence, total geometric count and marked-field
descent are imported from
[CONSTANT_ENDPOINT_DAVENPORT_ZANNIER_CLASSIFICATION.md](CONSTANT_ENDPOINT_DAVENPORT_ZANNIER_CLASSIFICATION.md).
The general fixed-coordinate reciprocal-unit interface is in
[ODD_START_FIXED_COORDINATE_FIBRE_CLASSIFICATION.md](ODD_START_FIXED_COORDINATE_FIBRE_CLASSIFICATION.md).
The singleton proof here does not assume its ineffective finite-base
enumeration has been carried out.

The field import is Milne, *Fields and Galois Theory*, version 5.10,
Theorem 9.19 and Lemma 9.24, printed pages 117--118, in the retained
[primary PDF](sources/milne_fields_galois_v510.pdf). The local differential
calculation (S10) is proved above; it agrees with the characteristic-zero
calculation in the Stacks Project, [tag 0C1B](https://stacks.math.columbia.edu/tag/0C1B),
Section 53.12, preceding Lemma 53.12.1. No global Riemann--Hurwitz
classification is required for the stable-field lemma.

For the tree interpretation, the actual primary statements inspected
are Pakovich--Zvonkin, *Minimum Degree of the Difference of Two Polynomials
over Q, and Weighted Plane Trees*, arXiv:1306.4141v1,
Definition 2.19 (page 9), Theorem 2.24 (page 10), and Theorem 3.1
(page 11), in the retained
[primary PDF](sources/pakovich-zvonkin-1306.4141v1.pdf).
Definition 2.19 explicitly identifies all-weight-one trees with Shabat
polynomials. The ordinary-tree counting argument here is elementary;
no novelty is claimed for rooted-tree enumeration or the classical
DZ/dessin correspondence.

The [exact checker](verify_davenport_zannier_coordinate_fibres.py),
[source ledger](davenport_zannier_coordinate_fibres_source_ledger.json),
and [receipt](davenport_zannier_coordinate_fibres_verification.json)
check explicit star and non-star fibres, the height cutoff, rigid
weighted examples, raw moments, the two outside-scope controls, and
the counting formulas against independent bounded tree enumeration.
They do not replace the unrestricted stable-field proof, the
all-height exhaustion, or the classical geometric correspondence.
No external human review or formal verification is claimed.

```bash
python3 research/foundations/extended-geometry/verify_davenport_zannier_coordinate_fibres.py \
  --repo-root . --output /tmp/davenport_zannier_coordinate_fibres_review.json
```
