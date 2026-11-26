# Pure quintic normal forms on the Hopf quadric

Let \(K\) be any characteristic-zero field and put
\[
 A=K[x,y,t]/(t^2+2xy-1),\qquad
 \Lambda(f)=\operatorname{CT}_x\frac12\int_{-1}^{1}
                 f\bigl(x,(1-t^2)/(2x),t\bigr)\,dt .
\]
The integral assigns rational values to monomials. Write \(\deg f\) for the
least total degree of a representative. A phase is pure if
\(\Lambda(f^m)=0\) for every \(m\ge1\).

## Theorem (Quintic purity and orbit equivalence)

A phase of least degree five is pure if and only if it is
carried by a polynomial \(K\)-automorphism of \(A\) to
\[
                         xP(x,t),\qquad \deg P\le4.                 \tag{1}
\]
Two nonzero representatives \(xP(x,t)\) and \(xP'(x,t)\), both of degree at
most five, are equivalent under polynomial \(K\)-automorphisms if and only
if, for some \(a\in K^*\), \(\varepsilon\in\{1,-1\}\), and \(Q\in K[x]\)
with \(\deg Q\le3\), one has the polynomial identity
\[
             P'(x,t)=aP(ax,\varepsilon t+xQ(x)).                   \tag{2}
\]
The identity means equality in \(K[x,t]\). The scalar \(a\) ranges over the
roots in \(K\) of an explicit equation of degree at most five. Each choice
of \(a,\varepsilon\) determines \(Q\) when \(P\) depends on \(t\). Thus
the criterion gives a finite test over fields with effective arithmetic
and effective univariate root finding.

By the [rank-one trace and marker formulas](RANK_ONE_TRACE_AUTOMORPHISM_INVARIANCE_ALL_IMAGE_PAIRS.md),
this gives normal forms for pure radially primitive balanced quintics in two
contraction pairs, using their restrictions to \(A\). Every such restriction
has least degree five. The equivalence is on the quadric, where all
polynomial automorphisms preserve \(\Lambda\).

## Existence and descent of the representatives

The [stable leading-form period obstruction](STABLE_LEADING_SPHERE_PHASE_PERIOD_OBSTRUCTION.md)
excludes quintics for which every root of the leading binary decic has
multiplicity at most four. The [exact-half-root argument](SPHERICAL_LEADING_BINARY_FORM_MOMENT_OBSTRUCTIONS.md#4-proof-of-the-elementary-obstructions)
excludes a root of multiplicity exactly five. A pure quintic therefore has a
unique geometric root of multiplicity at least six. Its uniqueness makes it
fixed by every automorphism of an algebraic closure over \(K\).
Characteristic zero makes every algebraic coefficient separable, so the
root is a point of \(\mathbb P^1(K)\). Choose a \(K\)-vector representing it,
and complete it to the first column of a determinant-one matrix. The
[symmetric-square coordinate formula](SPHERICAL_LEADING_BINARY_FORM_MOMENT_OBSTRUCTIONS.md#2-the-full-sphere-filtration-and-trace-preserving-changes)
identifies the linear isometry of \(t^2+2xy\) induced by this matrix. Applying
the corresponding substitution moves the leading root to \(v=0\) in the
binary coordinates \((x,y,t)=(u^2,-v^2/2,uv)\). The substitution preserves
the trace and degree and is defined over \(K\).

Use the [unstable quintic face exhaustion](UNSTABLE_QUINTIC_SPHERE_PURITY_EXCLUSION.md).
Its face tests, including the cyclic-cover exclusions, leave precisely the
following possibilities for a pure phase in the normalized coordinates.

* If there are no positive angular terms, the constant-angular part is a
  polynomial \(p(t)\), whose positive moments equal those of the whole
  phase. Müger and Tuset's [one-variable theorem](https://arxiv.org/pdf/1910.12148v1),
  Theorem 2.1, implies that \(p=0\). For coefficients in an arbitrary
  characteristic-zero field, apply it after embedding the finitely
  generated coefficient field into \(\mathbb C\) and rescaling the interval
  from \([-1,1]\) to \([0,1]\). The moments are rational polynomial
  expressions in the coefficients. The normal basis then gives
  \(f\in yK[y,t]\), and interchanging \(x\) and \(y\) gives (1).
* At the quadratic-height line \(\rho=H=2\), exchange \(x,y\). The phase
  now has the form \(x^{-1}(C(x)+tB(x)+t^2D(x))\) of the
  [quadratic-height classification](QUADRATIC_HEIGHT_HOPF_CLASSIFICATION.md#main-theorem).
  If \(C(0)=0\), that classification forces \(C,B,D\) to be divisible by
  \(x^2\), giving (1). Otherwise the nonconstant rational center has even
  least degree and is excluded in degree five. For the constant center
  \(a_0\in K^*\), write \(B=xb\). The phase is
  \(2a_0y+b(x)t-xb(x)^2/(4a_0)\). The substitution
  \[
      x\mapsto x,\quad t\mapsto t+xh(x),\quad
      y\mapsto y-th(x)-xh(x)^2/2,\qquad h=-b/(2a_0),
  \]
  carries \(2a_0y\) to this phase. Invert it and interchange \(x,y\) to
  obtain a linear representative of the form (1). The center is computed
  from \(C,B,D\), and the shear uses only elements of \(K\).
* For the remaining cases, the [order-five collision analysis](UNSTABLE_QUINTIC_SPHERE_PURITY_EXCLUSION.md#order-five-and-the-shears)
  uses the maps
  \[
  \theta_{c,r}(y)=y,\quad
  \theta_{c,r}(t)=t+cy^{r+1},\quad
  \theta_{c,r}(x)=x-cy^rt-\tfrac12c^2y^{2r+1}.
  \]
  In the \(H=5/3\) case the face is \(Bx+Ay^2t+Cy^5\),
  \(B\ne0\). Repetition gives \(A^2+2BC=0\). The inverse of
  \(\theta_{-A/B,2}\) leaves \(Bx+P_4(y)+(a_1+b_1y)t\) with
  \(\deg P_4\le4\). The next supporting lines eliminate the \(y^4\) term.
  The \(y^3,yt\) terms are either both absent or form a repeated face
  with \(Bx\). In the repeated case, if \(d_1\) is the \(y^3\)
  coefficient, then \(b_1^2+2Bd_1=0\); the inverse of
  \(\theta_{-b_1/B,1}\) cancels that face. The remaining terms below
  weight three for \((x,y,t)=(3,1,2)\) belong to the span of
  \(1,y,y^2,t\), in addition to the retained \(Bx\).
  In the \(H=5/2\) case the repeated face must be a cube. It equals
  \(B(z+c/4)^3\), where \(B\ne0\); its \(z^2\) coefficient determines
  \(c\) by division by \(3B/4\). The inverse of \(\theta_{c,1}\) leaves
  \(x(Bt+a_1)+P_4(y)+tP_2(y)+d_1t^2\), with
  \(\deg P_i\le i\). The successive two-term face tests remove
  \(y^4,y^2t,y^3\). In both cases the output has degree at most two.

In the last case, the [full quadratic criterion](TWO_PAIR_IMAGE_CUBIC_STRATIFICATION.md#the-full-quadratic-sphere-family-a-sharp-six-moment-criterion),
Proposition (Q1), turns a pure quadratic into \(ax^2+bxt+cx\) by a linear
isometry. Its common null direction is defined over \(K\). The proof uses
the rank-one image of the square of the nilpotent self-adjoint operator, the
rank-one image of the operator itself, or the line of the null linear vector
in the zero-operator case. A nonzero null vector \(u\) over \(K\) extends
to a hyperbolic pair over \(K\) by choosing \(v\) with \(B(u,v)=1\) and
replacing \(v\) by \(v-B(v,v)u/2\), where \(B\) is the symmetric form
associated to \(t^2+2xy\). The norm of a generator of the orthogonal
complement is a nonzero square since the full determinant and the
hyperbolic-pair determinant both have squareclass \(-1\). This gives the
required Witt coordinates over \(K\). All shears above preserve the trace.
The face exhaustion allows arbitrary lower coefficients.

The resulting representatives have degree at most five. Conversely, the
monomials of (1) have strictly positive angular weights, so every positive
pure moment vanishes. A fixed marker has finitely many angular weights, so
its product with sufficiently high powers of (1) has zero constant term.
Transporting the marker by the same automorphism also recovers quintic
safety.

## Corollary (Sharp mixed moment cutoff)

For every pure phase \(f\) of least degree five and every
nonzero marker \(g\in A\),
\[
                       \Lambda(g f^m)=0\qquad(m>5\deg g).          \tag{3}
\]
For each integer \(D\ge1\) there is a pure quintic and a degree-\(D\)
marker with nonzero mixed moment at \(m=5D\). In the original two-pair
polynomial ring the corresponding cutoff is
\(\mathcal E_2(QF^m)=0\) for \(m>5\deg_z Q\), again sharp for each
positive marker degree.

Indeed the normalizing automorphism \(\alpha\) taking \(f\) to (1)
can be chosen with all coordinate images of degree at most five. In
the constant-center branch, \(b\ne0\) and \(\deg B=3\) by the
[exact degree formula](QUADRATIC_HEIGHT_HOPF_CLASSIFICATION.md#5-exact-least-degree-parity-and-the-image-normalization),
so \(\deg h=2\). In the \(H=5/3\) branch both inverse shears fix \(y\).
They act on \(t\) by adding polynomials \(yQ(y)\) and compose by
adding their \(Q\) parameters. Their composition is one shear with
\(\deg Q\le2\). The \(H=5/2\) branch uses a single shear with linear
\(Q\). These are the only nonlinear normalizations, and compositions
with the linear isometries before or after a shear preserve the bound.
The inverse has the same degree bound. Thus \(\alpha(g)\) has degree
at most \(5\deg g\). In the normal basis its least angular weight is
at least \(-5\deg g\). Every power \((xP)^m\) has least angular
weight at least \(m\), which gives (3) by trace invariance.

For sharpness, choose any \(c\in K^*\) and take
\(f=\theta_{c,2}(x)\), \(g=x^D\). Its degree-five leading term is a
nonzero multiple of \(y^5\). After applying \(\theta_{-c,2}\), the
unique term of weight \(-5D\) in the marker is
\((-c^2/2)^D y^{5D}\). At \(m=5D\) the mixed moment is
\[
                 (-c^2/2)^D\Lambda((xy)^{5D})\ne0
\]
since \(\Lambda((xy)^n)=2^n(n!)^2/(2n+1)!\).

For the original contraction, separate the bihomogeneous components
of \(Q\). A component of bidegree \((u,v)\) with \(u>v\) contracts to
zero by degree. For \(u\le v\), its output against \(F^m\) is
homogeneous of \(z\)-degree \(v-u\). Every top derivative is the
contraction of \(F^m\) against the balanced marker \(\xi^\eta Q\)
with \(|\eta|=v-u\), since contraction intertwines multiplication by
\(\xi_i\) with \(\partial_{z_i}\). The restriction of this balanced
marker has degree at most \(v\); the factorial factor in the
[balanced contraction formula](RANK_ONE_TRACE_AUTOMORPHISM_INVARIANCE_ALL_IMAGE_PAIRS.md#3-explicit-coinvariants-give-the-normalized-image-trace)
is nonzero. At \(m>5v\) all top derivatives vanish by (3), so the
homogeneous output is zero. The largest \(v\) is \(\deg_z Q\).
The balanced degree-\(D\) lift of \(x^D\) and the primitive balanced
degree-five lift of \(\theta_{c,2}(x)\) transfer the sharpness example.

## Equivalence of the representatives

We use the [Blanc–Dubouloz generator theorem](sources/blanc-dubouloz-a1-fibered-surfaces-2011.pdf),
Theorem 5.4.5(7),(8), as specialized in
[Coordinate degrees](HOPF_COORDINATE_DEGREES_AND_FULL_HYPERBOLIC_SPECTRA.md#2-the-precise-classical-input-and-a-useful-word-form).
For this split quadratic surface over \(K\) it gives generation by the
linear isometry group \(L\) and by the group \(J\) of substitutions
\[
 x\mapsto ax,\quad t\mapsto\varepsilon t+xQ(x),\quad
 y\mapsto a^{-1}\bigl(y-\varepsilon tQ(x)-xQ(x)^2/2\bigr) .          \tag{4}
\]
Here \(a\in K^*\), \(\varepsilon^2=1\), and \(Q\in K[x]\). This follows
from their degree-two case with two \(K\)-rational roots under
\((u,v,w)=(x,-2y,t)\). The fibered generators in their proof give (4) by
preserving \(uv=w^2-1\). In particular the original \(v\) summand of a
fibered shear is retained, as required by that equation.

Let \(B=L\cap J\), the subgroup in (4) with constant \(Q\). Let \(s\)
exchange \(x,y\) and fix \(t\). If \(\ell\in L\setminus B\), write
\(\ell(x)=ux+vy+wt\). It is isotropic, so \(2uv+w^2=0\). If \(v=0\), then
\(\ell(x)\in K^*x\) and the isometry belongs to \(B\). Hence \(v\ne0\).
Choosing \(a=1/v,\varepsilon=1,Q=-w/v\) in (4) gives \(b\in B\) with
\(b(y)=\ell(x)\). The isometry \((bs)^{-1}\ell\) fixes \(x\) and belongs
to \(B\). This proves \(\ell\in BsB\).

Merge consecutive factors in each subgroup of an alternating word in
\(L,J\), and absorb any interior factor in \(B\) into its neighbors.
Replace each factor in \(L\setminus B\) by an element of \(BsB\). The
result has form \(\sigma=j_l s j_{l-1}s\cdots s j_0\), with \(j_i\in J\)
and each interior \(j_i\) nonlinear. These reductions need only the
generation assertion. We compose as substitutions on elements of \(A\).

The associated graded ring for the degree filtration is
\(K[x,y,t]/(t^2+2xy)\). Indeed, the leading homogeneous part of any
multiple of \(t^2+2xy-1\) is a multiple of \(t^2+2xy\). The latter form
has rank three and is absolutely irreducible, since a product of two
linear forms has rank at most two. The graded quotient is a domain, so
\(\deg(fg)=\deg f+\deg g\) for nonzero \(f,g\). Also
\(\deg R(f)=(\deg R)\deg f\) when \(R\) is a nonzero univariate polynomial
and \(f\) is nonconstant. The monomials in \(K[x,t]\) have ordinary total
degree in this filtration, by the normal basis.

Write \(q_i=\max(1,1+\deg Q_i)\) for the parameter of \(j_i\), taking
\(q_i=1\) if \(Q_i=0\). For \(l\ge1\) the word gives
\[
 \deg\sigma(x)=\prod_{i=1}^l(2q_i-1),\qquad
 \deg\sigma^{-1}(x)=\prod_{i=0}^{l-1}(2q_i-1).                     \tag{5}
\]
To see this, the coordinate images of the leftmost \(j_l\) have degrees
\((1,2q_l-1,q_l)\). Appending a factor on the right substitutes the
current images \((X,Y,T)\) into the coordinate formula of that factor.
Immediately after appending a swap, the first image has degree \(D\)
equal to the previous second degree, which was largest, and the third
degree is at most \(D\). The next nonlinear \(j_i\) leaves the first
degree equal to \(D\). Its \(XQ_i(X)\) term has degree \(q_iD>D\), so
dominates \(\varepsilon_iT\) strictly. The resulting second degree
is \((2q_i-1)D\) because the positive degrees of the two terms in
\(T'^2+2X'Y'=1\) must agree. This proves the first product, including
the case that the leftmost \(j_l\) is linear. The final \(j_0\) merely
scales the first coordinate. Apply the same argument to the inverse
word, whose inverse fibered factors have the same \(Q\)-degrees.

If \(g=\sigma(f)\) with \(f,g\in xK[x,t]\) and both degrees at most
five, then \(\sigma(x)\) divides \(g\) and \(\sigma^{-1}(x)\) divides
\(f\). Additivity bounds both degrees in (5) by five. Every interior
factor contributes an odd factor at least three to both products.
There can therefore be at most one interior factor, and \(l\le2\).

Each \(j\in J\) preserves \(xK[x,t]\) in both directions, so the
endpoint factors can be removed when testing membership in that space.
With one swap, the middle sends \(xK[x,t]\) into \(yK[y,t]\); their
intersection is zero by their strict opposite angular weights. With
two swaps, the middle is \(sj_1s\). Its scaling and reflection can
also be removed since their conjugates by \(s\) preserve \(xK[x,t]\).
The remaining map fixes \(y\) and sends \(t\mapsto t+yQ(y)\) with
\(k=\deg Q\ge1\).

If its image of \(f=xP\) belongs to \(xK[x,t]\), it vanishes on both
lines \(x=0,t=\delta\) for \(\delta=1,-1\). Thus \(f\) itself vanishes
on the images of these two lines under the substitution. Their
projections to the \((X,T)\)-plane are
\[
                 X=-\delta Q(y)-yQ(y)^2/2,\qquad
                 T=\delta+yQ(y).
\]
Writing \(Q(U)=\sum_{i=0}^k q_i U^i\) temporarily, their equations are
\[
 R_\delta(X,T)=2X^{k+1}+(T+\delta)
       \sum_{i=0}^k q_i\bigl((1-T^2)/2\bigr)^i X^{k-i}=0.          \tag{6}
\]
For irreducibility, invert \(X\) and \(1-T^2\). The relation
\(y=(1-T^2)/(2X)\) identifies this open plane with the open
\((y,T)\)-chart where \(y\) and \(1-T^2\) are invertible.
Equation (6) is then a unit times \(T-\delta-yQ(y)=0\), which has
prime ideal because its quotient is a localization of \(K[y]\).
This quotient is nonzero since \(\delta+yQ(y)\) is nonconstant.
The discarded linear factors \(X,T-1,T+1\) do not divide
\(R_\delta\). At \(X=0\) the polynomial is
\((T+\delta)q_k((1-T^2)/2)^k\) and at \(T=\pm1\) it has nonzero
leading term \(2X^{k+1}\). Thus \(R_\delta\) is irreducible in
\(K[X,T]\). The two prime ideals are distinct on the open chart because
their linear equations in \(T\) differ by two. Each polynomial has
degree \(2k+1\), from the unique top \(T^{2k+1}\) term.

On the open chart the images have \(X\ne0\), so the vanishing of
\(f=XP(X,T)\) forces \(P\) to vanish on both curves. Equivalently,
each graph substitution into \(P\) vanishes as a polynomial in \(y\).
Both irreducible polynomials therefore divide \(P\) in \(K[X,T]\).
Their product has degree \(4k+2\ge6\), contradicting \(\deg P\le4\).
This excludes the two-swap case. Hence every equivalence between the
representatives lies in \(J\).

Finally, to bound \(Q\) in an equivalence in \(J\), write
\(f=\sum_{j=0}^m a_j(x)t^j\) with \(m=\deg_t f\). If \(m=0\), take
\(Q=0\). Otherwise suppose \(Q\ne0\). The \(t^{m-1}\) coefficient
after substitution is
\[
 \varepsilon^{m-1}\bigl(a_{m-1}(ax)+m a_m(ax)xQ(x)\bigr).
\]
Both this coefficient and \(a_{m-1}(ax)\) have degree at most
\(6-m\). Subtracting gives the same bound on their difference,
which has degree \(\deg a_m+1+\deg Q\). Here \(a_m\ne0\) is divisible
by \(x\). Therefore \(\deg Q\le4-m\le3\). Formula (4) now gives
exactly (2), with the converse realized by that substitution.

For the finite test write \(g=\sum_{j} b_j(x)t^j\). If its degree in
\(t\) differs from \(m\), the representatives are inequivalent.
Otherwise, the leading coefficient identity is
\(b_m(x)=\varepsilon^m a_m(ax)\). Choose a nonzero coefficient
\([x^r]a_m\), where \(1\le r\le5-m\). For each sign, compare that
coefficient to get
\[
                     a^r=\frac{[x^r]b_m}{\varepsilon^m[x^r]a_m}.
\]
Try the nonzero roots in \(K\), of which there are at most \(r\). For
\(m>0\), the next coefficient equation fixes
\[
 Q(x)=\frac{\varepsilon^{1-m} b_{m-1}(x)-a_{m-1}(ax)}
                 {m x a_m(ax)} .
\]
Require this rational expression to be a polynomial of degree at most
three, then check all coefficients of the proposed equality. For
\(m=0\) the only remaining test is \(g(x)=f(ax)\).

For example, \(y-x^2t-x^5/2\) is a degree-five mixed-angular phase
carried to the linear representative \(x\) by the inverse of the
shear with \(Q=x^2\) followed by interchanging \(x,y\).
Among degree-five representatives, \(x^5+xt\) and \(2x^5+xt\)
are equivalent by \(t\mapsto t+x^4\). Any triangular equivalence
between them has \(a\varepsilon=1\), hence \(a=\pm1\), and its
\(x^5\) coefficient forces a nonzero cubic term in \(Q\). Thus the
uniform cubic bound in (2) is attained. The representatives \(x^5\)
and \(x^4t\) are inequivalent because (4) preserves degree in \(t\).
