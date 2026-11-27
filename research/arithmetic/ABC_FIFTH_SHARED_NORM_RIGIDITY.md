# Rigidity of small primitive norms shared by a quadratic twist

The [small shared norm construction](QUADRATIC_TWIST_SHARED_SMALL_TORSION_NORMS.md)
gives a sharp count for primitive records with bounded quadratic
coefficients at minimum norm size. The radicand typically determines
the whole record there. A common-value surface compares the fixed
binary forms from the real twist field. Every rational line on such
a surface with nonzero form value comes from a projective symmetry,
and primitivity removes the apparent freedom to rescale the record.

**Theorem SN3.** Fix a positive squarefree integer \(t>1\), an odd
integer \(\ell\ge3\), and constants \(C,L>0\). For signed squarefree
\(d\ne1\) with \(1\le |d|\le X\), consider the integer records
\[
 \begin{split}
 x^2-dy^2&=z^2-tdw^2=m^\ell,\qquad 1\le m\le C X^{1/\ell},\\
 0<|y|,|w|&\le L,\qquad \gcd(m,xz)=1.
 \end{split}                                                     \tag{1}
\]
Except for \(O_{t,\ell,C,L,\epsilon}(X^{\beta_\ell+\epsilon})\)
radicands, there is at most one such record up to independent
sign changes of \(x,z,y,w\). Here
\[
                 \beta_3=\frac{1}{2\sqrt6},\qquad
                 \beta_\ell=\frac1{2\ell}\quad(\ell\ge5).
 \tag{2}
\]
For \(\ell=3,5\) and sufficiently large fixed \(C,L\), there are
\(\gg_{t,\ell}X^{1/\ell}\) negative squarefree \(d\) with a
unique such record and exact order-\(\ell\) ideal classes in both
quadratic fields. At \(t=\ell=5\) one may take \(L=1\).

## Root sets over the fixed real quadratic field

Take \(y,w>0\) by choosing signs. In \(F=\mathbf Q(\sqrt t)\) put
\(\theta=zy+xw\sqrt t=A+B\sqrt t\) and \(c=y^2-tw^2\ne0\).
Then
\[
        N_F(\theta)=c m^\ell,\qquad
        d=\frac{B^2/w^2-A^2/y^2}{c}.                             \tag{3}
\]
By the ideal factorization in the
[shared norm upper bound](QUADRATIC_TWIST_SHARED_SMALL_TORSION_NORMS.md#upper-bound-at-the-same-shared-norm-scale),
for each coefficient pair there are finitely many fixed
\(\gamma\in F^\times\) and fractional ideal lattices with
\(\theta=\gamma\lambda^\ell\). The lattice coordinates of
\(\lambda\) in fixed bases have size \(O_{t,\ell,C,L}(X^{1/(2\ell)})\).
Thus the second expression in (3) gives finitely many rational
binary forms of degree \(2\ell\) representing \(d\). Their root
sets, in the projective coordinate \(q=\lambda/\bar\lambda\), are
\[
      q^\ell=(\bar\gamma/\gamma)R^{\pm1},\qquad
                   R=\frac{y+w\sqrt t}{y-w\sqrt t}.               \tag{4}
\]
Both right sides are nonzero, and \(|R|>1\), so the \(2\ell\)
roots are distinct. The ratio of the larger circle radius to the
smaller in the \(q\)-coordinate is \(|R|^{2/\ell}\).

Each root set has projective rotations of order \(\ell\), with
fixed pair \(0,\infty\). The full stabilizer over \(\mathbf C\)
is finite since it acts faithfully on the root set. The classical
finite subgroups of \(\mathrm{PGL}_2(\mathbf C)\) are cyclic,
dihedral, \(\mathfrak A_4,\mathfrak S_4,\mathfrak A_5\).
A point stabilizer is cyclic, since a finite subgroup fixing
infinity injects into the group of affine multipliers.
[Beauville, introduction](https://arxiv.org/pdf/0909.3942)

For \(\ell\ge7\) only the cyclic and dihedral types have an
element of order \(\ell\). For \(\ell=5\), the additional
possibility \(\mathfrak A_5\) is excluded since every point
orbit has size at least \(60/5=12\), and the root set has size ten.

For \(\ell=3\), an \(\mathfrak A_5\)-orbit is again too large.
A six-point set invariant under \(\mathfrak A_4\) or
\(\mathfrak S_4\) has to be the octahedral orbit, up to projective
coordinates. Indeed \(\mathfrak A_4\), in its tetrahedral rotation
action, has six points fixed by its three distinct involutions,
and all other orbits have size four or at least twelve. A six-point
invariant set can therefore only be that orbit. In the
\(\mathfrak S_4\) octahedral action, a six-point orbit is the
fixed pairs of the three order-four rotation axes, since other
point stabilizers have order at most three. The two sets in these
models are the same six points \(\{\pm e_1,\pm e_2,\pm e_3\}\) on
the unit sphere. Any order-three axis in these actions is a body
diagonal. With the diagonal along the vertical axis, the six
points have signed heights \(\pm1/\sqrt3\). Stereographic
projection sending the two poles to \(0,\infty\) gives the
ratio of radii \((1+1/\sqrt3)/(1-1/\sqrt3)=2+\sqrt3\).
Every other projective coordinate with the same two poles
preserves this ratio between larger and smaller circles. Hence
(4) would require \(R^2=(2+\sqrt3)^3\). This requires \(t=3\)
by the distinct quadratic fields, and would make \(2+\sqrt3\)
a square in \(\mathbf Q(\sqrt3)\). In an expression
\((u+v\sqrt3)^2=2+\sqrt3\), \(2uv=1\) and
\(u^2\in\{1/2,3/2\}\), impossible for rational \(u,v\).

Thus the stabilizer of every root set in (4) is cyclic or
dihedral. In either case its cyclic subgroup of order \(\ell\)
is unique, as \(\ell\) is odd. Any projective equivalence
between the forms conjugates these subgroups and preserves
their fixed pair \(0,\infty\) in the respective coordinates.

## Equivalent forms give the same primitive record

A rational invertible linear equivalence between two forms
from (3) is consequently diagonal or antidiagonal in the two
embedding coordinates of \(\lambda\). Its two nonzero entries
are Galois conjugate, since the original map is rational.
Since \(\theta=\gamma\lambda^\ell\), on the \(\theta\)
coordinates this gives either \(\theta'=\rho\theta\) or
\(\theta'=\rho\bar\theta\) for some \(\rho\in F^\times\).
Write \(y',w'>0,c',R'\) for the corresponding data.

The root pairs in the \(\theta/\bar\theta\) coordinate are
\(\{R,R^{-1}\}\) and \(\{R',R'^{-1}\}\). Scaling by
\(\rho/\bar\rho\) and optionally inverting sends one pair
to the other. The ratio of moduli yields \(|R'|=|R|\), and
the image of the member of larger modulus gives
\(\rho/\bar\rho=R'/R=\pm1\).

If \(R'=R\), then \(y'=k y,\ w'=k w\) for some positive
rational \(k\) and \(\rho\) is rational. Equality of the
forms for \(d\) in (3) gives \(\rho^2=k^4\), so
\[
             |x'|=k|x|,\quad |z'|=k|z|,\quad m'^\ell=k^2m^\ell .
\]
If any \(v_p(k)>0\), then \(p\) divides both \(m'\) and the
two output traces, violating primitivity. If \(v_p(k)<0\),
interchange the two records. Hence \(k=1\), which gives the
same record up to signs.

If \(R'=-R\), then \(y'/w'=t w/y\), so write
\(y'=k t w,\ w'=k y\) with \(k>0\) rational. Now
\(c'=-k^2 t c\) and \(\rho=h\sqrt t\) for \(h\ne0\)
rational. The transformation swaps the absolute traces and
multiplies them by \(|h/k|\). Equality in (3) would
require \(h^2=k^4t\), impossible for nonsquare \(t\).
Therefore every rational equivalence of two of the forms
produces only sign changes of the record.

## Counting exceptions over varying radicands

For each pair of binary forms \(P,Q\) from (3), clear their
denominators by a common integer and consider the surface
\(P(a,b)-Q(a',b')=0\) in \(\mathbf P^3\).
It is smooth: the two partial derivatives of each binary
form vanish simultaneously only at its zero vector, by
homogeneity and separability. A rational line on the surface
with nonzero common form value must be the graph of an
invertible linear equivalence. Parametrize the line by two
linear coordinate pairs. If one projection has rank two
and the other lower rank, a separable degree-\(2\ell\)
binary form would equal a single power or zero. If both
projections have rank one and the common form is nonzero,
the proportional powers force the same linear parameter,
contradicting that the image is a line. Both projections
therefore have rank two, giving the graph. Such a line
contributes only sign changes by the preceding section.

In degree \(D=2\ell\), the smooth-surface bound off lines
is \(O_\epsilon(H^{3/\sqrt D+\epsilon}+H^{1+\epsilon})\)
at height \(H=O_{t,\ell,C,L}(X^{1/(2\ell)})\). This gives
the exponents (2). The bound also counts integer lifts of
that size by summing over their common divisor with fixed
positive slack in the exponent. A nonrational line accounts
for finitely many rational projective points, since each
also lies on a distinct conjugate line. A smooth surface
of degree \(D\ge6\) has finitely many lines: a line has
normal bundle of degree \(2-D<0\) in the surface and is
isolated in the projective scheme of lines. Each of the
finitely many resulting points contributes \(O(H)\)
integer lifts. This proves the bound on the number of
radicands with inequivalent records.
[Stewart–Xiao, equation (2.5)](https://arxiv.org/pdf/1605.03427)

## A growing coefficient window

**Corollary SN4.** Fix \(t\) and odd \(\ell\ge3\) as above.
Allow \(1\le L\le X^K\) for any fixed \(K\) and replace the
norm cutoff in (1) by \(m^\ell\le C X L^2\).
Then the number of radicands admitting inequivalent records is
\[
 \ll_{t,\ell,C,K,\epsilon}
      (XL)^\epsilon X^{\sigma_\ell/(2\ell)}
          L^{4+2\sigma_\ell/\ell},
 \qquad
       \sigma_\ell=\max(1,3/\sqrt{2\ell}).                         \tag{6}
\]
In particular, for \(t=\ell=5\) and
\(0<\eta<1/110\), a power-saving exceptional bound holds
in the imaginary window \(X\le |d|<2X,\ m\le X^{1/5+\eta}\)
even with arbitrary integer nonzero quadratic coefficients.
It is \(O_\epsilon(X^{1/10+11\eta+\epsilon})\).
At least order \(X^{1/5}\) fifth-twisted imaginary pairs at
this window have genuine class witnesses and remain unique
at this larger norm cutoff.

For each \(y,w\), the remainder ideal \(\mathfrak c\) in
the extraction of \(\theta=\gamma\lambda^\ell\) has
exponents less than \(\ell\) and is supported on the primes
over \(2ty w(y^2-tw^2)\). The number of choices is
\(O_{t,\ell,\epsilon}(L^\epsilon)\) since
\(2tyw|y^2-tw^2|\ll_t L^4\) and \(t,\ell\) are fixed.
The integral class representatives \(\mathfrak j\) of
\(\mathfrak a\) are fixed in \(F\). For each principal
\(\mathfrak c\mathfrak j^\ell\), choose a generator whose
two absolute embeddings are comparable, using a unit of
the fixed field. Its absolute norm is a positive integer
at least one. After multiplication by one of the fixed
unit representatives modulo \(\ell\)-th powers, both
absolute embeddings of \(\gamma\) are still bounded below
uniformly. Consequently the lattice coordinates of
\(\lambda\) satisfy
\[
                        H\ll_{t,\ell,C} X^{1/(2\ell)} L^{2/\ell}.
\]
Indeed the norm and coefficient bounds give
\(|x|,|z|\ll_{t,C}\sqrt X L\), so both embeddings of
\(zy+xw\sqrt t\) have size \(O_{t,C}(\sqrt X L^2)\).
The fractional ideal lattices for \(\lambda\) belong to
the same fixed list from the class representatives.
There are \(O_{t,\ell,\epsilon}(L^{2+\epsilon})\) forms
over the coefficient window. Apply the smooth-surface
estimate to every ordered pair, using the bound uniform
in the surface of fixed degree and the degree-dependent
bound on the number of lines there. The rational graph
equivalences remain exactly as in the fixed-coefficient
proof. This gives (6).

In the imaginary window the norm equation itself gives
\(|y|,|w|\ll X^{5\eta/2}\) for \(\ell=5\).
Using \(L=O(X^{5\eta/2})\) in (6) gives the stated
exponent, strictly below \(1/5\) when \(\eta<1/110\).
For the lower count in this dyadic window, the small
shared norm construction can use a sufficiently small
rectangle in its sign region around fixed scaled
parameters with the form value strictly between \(-2\)
and \(-1\). The squarefree and representation arguments
apply to this rectangle, by homogeneity and differences
of comparable rectangles. The coefficient-one variant
does the same. Thus \(\gg X^{1/5}\) fields in this window
have a small primitive shared witness, and the bound
(6) removes only \(o(X^{1/5})\) of them. No shared prime
factorization of the norm is imposed.

## Lower counts with uniqueness

The [sharp lower construction](QUADRATIC_TWIST_SHARED_SMALL_TORSION_NORMS.md)
at \(\ell=3,5\) supplies \(\gg_{t,\ell}X^{1/\ell}\) negative
squarefree \(d\) satisfying (1) with fixed \(C,L\) and both
classes of exact order \(\ell\). The exceptional exponents
in (2) are smaller than \(1/\ell\). Changing signs in
(1) only conjugates the ideal roots in each field.

For the lower count with \(L=1,t=\ell=5\), take
\(A+B\sqrt5=(a+b\sqrt5)^5\) with \(a/b>0\) fixed in a
sufficiently small interval. The form \(d=A^2-B^2\) is
negative of size comparable to \(b^{10}\).
As in (4), \(A-B\) and \(A+B\) are degree-five forms with
simple, disjoint projective roots. Impose \((a,b)\equiv(1,0)\)
modulo the product of \(2,5\), all primes up to ten, and all
primes dividing the content or discriminant of this form.
It follows that \(d\) and \(m=5b^2-a^2\) are units there.
The squarefree sieve and essential representation count
give \(\gg X^{1/5}\) distinct squarefree \(d\) in a fixed
rectangle with \(b\asymp X^{1/10}\) and \(|d|\le X\), using
the same local density check as in the lower construction.

Now \(x=2B,z=2A\), and
\(x^2-d=z^2-5d=5B^2-A^2=m^5\).
At \(p\mid m\) the reduction of \(a+b\sqrt5\) vanishes in
one component only, since the squarefree value of \(d\)
forces \(\gcd(a,b)=1\). Hence \(A,B\) are both units and
(1) is primitive. The irrational coefficients are exactly
one in the two squarefree radicands \(d,5d\).
The ideals of norm \(m\) are nonprincipal by the comparison
\(4m<\min(|d|,|5d|)\) and have classes of exact order five.
Deleting the exceptional radicands preserves the lower count.

For example, \((a,b)=(1,10)\) gives \(d=-4927727492499\),
\(m=499\), \(x=5100100\), and \(z=2510002\).
The factorization of \(-d\) is
\(3\cdot7\cdot23\cdot83\cdot743\cdot165437\).
The two square tests in (1) for \(1\le m\le1000\) and
\(y=w=1\) give only this \(m\) and the corresponding
absolute traces.

Rigidity here concerns shared norms in the minimum-height
range with bounded coefficients. It constrains attempts to
get independent torsion classes by repeating that construction
in the same pair of fields. Witnesses with larger coefficients
or different norms in the two fields remain possible.
