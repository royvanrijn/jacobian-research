# The target cylinder and the outer branch of the order-\((12,42)\) row

The second resultant supplies a target coordinate that
has a finite limit along the residual position cover.
Near position infinity this cover has nontrivial local
monodromy. In the split outer profile, neither affine
ramification curve approaches that finite target limit.
The two facts are incompatible. The remaining outer chart
also gives an exact conductor count for its image, and
the polar curve gives an additional resultant identity.

Work with the planar transfer from the
[two Euler faces and outer profile](POLYNOMIAL_WEYL_ORDER_TWELVE_FORTY_TWO_OUTER_FIBRES.md).
Write \(H\) for the number of high roots,
\(D=\deg h/3=5H\) and \(M=D/2\).
All the geometric roots of \(h\) have multiplicity fifteen,
and \(H\) is even and at least two.
The planar resultants are
\[
 \begin{split}
 W&=Q^2-P A(P)^2+2cA(P),\\
 Z&=W^3-\delta^3P^5Q-e WP^2Q ,
 \end{split}
\]
where \(A\) is monic cubic and the leading coefficients
in \(X\) of \(W,Z\) are
\(\delta h^{17/3}\) and \(K B h^{49/3}\),
with \(K=-\kappa\delta^2/2\ne0\).
Here \(B'-(h'/6h)B=1\).
These are the
[shared planar resultants](POLYNOMIAL_WEYL_ORDER_TWELVE_FORTY_TWO_PLANAR_FIBRES.md#the-shared-planar-resultant).

## A target chart with the residue variable

**Theorem PlanarFortyTwoUnsplitOuterProfile.**
The outer profile of a hypothetical reduced \((12,42)\) Weyl
pair after the two-Euler planar transfer is unsplit.
Its terminal scale is \(k=M=5H/2\), and the genus of the
generic \(P\)-fibre is \(10H-11\).
For \(H=2\) this leaves geometric degree \(16\),
twenty punctures and genus \(9\).

We may work over \(\mathbb C\) by extending the finitely
generated coefficient field of the transferred data.
Write \(t=P^{-1/12}\) as a target chart parameter.
The second parameter \(u\) will be \(t^{99}Z\).
Here \(W,Z\) also denote their explicit polynomials in
the two target coordinates. Put \(V=t^{34} W\) and
\(a(t)=t^{36} A(t^{-12})\), so \(a(0)=1\).
There is a holomorphic branch \(\widetilde q(t,V)\) near
\(t=0,V=\delta\) specified by
\[
     \widetilde q(t,V)^2
        =a(t)^2-2ct^{48}a(t)+t^{50}V,\qquad
                   \widetilde q(0,V)=1 .
\]
The target values in these coordinates are
\(P=t^{-12}\), \(Q=t^{-42}\widetilde q\).
The second resultant gives
\[
         t^3 u=
         V^3-\delta^3\widetilde q(t,V)
                  -e t^2 V\widetilde q(t,V).
\]
At \(t=0,V=\delta\) the derivative of the right-hand
side with respect to \(V\) is \(3\delta^2\ne0\).
Thus near \((0,u_0)\) for each finite \(u_0\) it fixes
\(V=V(t,u)\) holomorphically. In particular for small
nonzero \(t\) a bounded \(u\)-disc gives an analytic
disc in the affine target line \(P=t^{-12}\).
The two identities ensure \(t^{99}Z=u\) there.

Let \(\pi\) be the finite cover of the affine target
obtained by normalization in \(\mathbb C(x,X)\).
This is the usual finite part of a Stein factorization
of a resolved map, as in
[Borisov, Section 3](https://people.math.binghamton.edu/borisov/documents/papers/ampleramification.pdf).
For the argument it suffices that the branch locus is
algebraic and that its curve components can be read at
the generic \(P\)-level. In fact normalization is
finite, and localization at the generic \(P\)-value gives
the finite normalization of the affine \(Q\)-line in
the normalized generic fibre's function field.
On that curve, the Keller condition gives no ramification
at affine source points. At finite high roots the
finite-\(Q\) points are unramified. The \(Q\)-poles do
not lie over the affine \(Q\)-line.
The [outer point calculation](POLYNOMIAL_WEYL_ORDER_TWELVE_FORTY_TWO_OUTER_FIBRES.md#the-outer-newton-path)
therefore accounts for every horizontal affine branch
component of \(\pi\). Its image is one of the parametrized
curves \((P_0(v),Q_0(v))\) from the outer charts.
Any vertical branch components, and any isolated points
of the non-étale locus, project to finitely many \(P\)-values.

Suppose the outer profile is split.
Both outer curves have \(P_0\)-degree six and their
\(Z_0=Z(P_0,Q_0)\) is polynomial in \(v\).
As \(P_0(v)=t^{-12}\) goes to infinity, \(v\) must go to
infinity. For \(Z_0\ne0\), the magnitude of
\(u=t^{99}Z_0(v)\) is asymptotic to a positive
constant times
\[
                         |v|^{\deg Z_0-99/2}.
\]
Since the exponent is nonzero, the image tends either
to \(u=0\) or to infinity. An identically zero \(Z_0\)
also gives \(u=0\). Therefore every closed disc around
a finite nonzero \(u_0\), sufficiently small to miss
zero, has no branch point of the pulled-back cover
for small nonzero \(t\).
Vertical branch values and isolated exceptions are also
avoided since \(P=t^{-12}\) is large.

We now exhibit nontrivial monodromy on such a disc.
Choose an analytic branch \(\lambda=h^{1/6}\) near
\(x=\infty\), meromorphic there with pole order \(M\).
It exists since the leading power of \(h\) at infinity
has degree \(6M\). At \(w=X^{-1}=0\) with \(x\) finite
and sufficiently large, the coordinate \(t\) of the
actual pair is
\[
                    t=\frac{w}{\lambda(x)}(1+O(w)).
\]
It is invertible as a function of \(w\) near zero.
On substituting \(w=w(t,x)\) into the actual polynomials,
one has \(t^{42}Q\to1\), \(t^{34}W\to\delta\), and
\[
                     t^{99}Z\ \longrightarrow\
                               u_{\rm res}(x)=K B(x)/\lambda(x).
\]
The choices of signs and of \(V\) agree with the target
chart just constructed. The function \(u_{\rm res}\)
has finite nonzero value \(u_\infty\) at infinity,
since \(\deg B=M\). It satisfies
\(u_{\rm res}'=K/\lambda\).
Consequently \(u_{\rm res}-u_\infty\) has order \(M-1\)
in \(1/x\). Here \(M-1\ge4\).

A small circle around \(u_\infty\) thus has \(M-1\)
local inverse lifts near \(x=\infty\), permuted by
an \((M-1)\)-cycle. Choose the circle so that these
lifts stay in a compact annulus of large finite
\(x\)-values and have nonzero derivative.
The equation for \(t^{99}Z\) after solving for \(w\)
is a holomorphic perturbation of \(u_{\rm res}\) on
this annulus. The implicit function theorem, uniformly
on a finite cover of the lifted loops, preserves that
permutation for small nonzero \(t\).
These lifts are affine source points. They give
distinct sheets of the finite normalization over the
branch-free locus. Indeed the affine source maps
regularly and birationally into the normalization:
an element integral over the target is in the normal
source coordinate ring. Over a branch-free neighbourhood
two different affine inverse branches cannot map to one
normalization sheet. If they did, the inverse function
theorem would give two local inverse branches of
the birational map onto an open subset of one sheet.

For split outer charts the disc bounded by this circle is
branch-free in the pulled-back finite cover, on which
every loop in the disc has trivial monodromy.
This contradicts the \((M-1)\)-cycle.
Thus the outer chart is unsplit. Substituting \(k=M\)
in the [exact genus formula](POLYNOMIAL_WEYL_ORDER_TWELVE_FORTY_TWO_OUTER_FIBRES.md#exact-geometric-degree-and-genus)
gives the result.

## The plane curve of outer branch values

**Corollary PlanarFortyTwoOuterBranchSemigroup.**
In the remaining unsplit profile the outer chart's
image is a birational polynomial parametrization with
pole-order semigroup \(\langle12,42,34,99\rangle\).
The image has total affine \(\delta\)-invariant
\(99\). The outer chart gives a generic transverse
ramification index \(M-1\) over this image.
Each regular high-root chart at infinity has image
semigroup \(\langle2,7\rangle\) and total affine
\(\delta\)-invariant three.

Let \(P_0,Q_0,W_0,Z_0\) denote restrictions to the
unsplit outer boundary \(x=\infty\), with parameter
\(v\). Their respective degrees in \(v\) are
\(12,42,34,99\). For the last degree, the leading
\(X\)-coefficient \(K B h^{49/3}\) has \(x\)-degree
\(M+49D=99M\); Laurent shifts of \(X\) preserve that
coefficient, and \(Z\) is regular in the terminal
chart because it is polynomial in \(P,Q\).
The gcd of the degrees is one, so the field
\(k(P_0,Q_0)\) is \(k(v)\): the intermediate degree
would divide each degree by the tower law.

The equation of \(Q_0\) over \(k(P_0)\) is monic of
degree twelve with coefficients in \(k[P_0]\),
by integrality in \(k[v]\) over \(k[P_0]\).
Thus the monomials
\[
       P_0^i Q_0^j W_0^a Z_0^b,\qquad
       i\ge0,\quad 0\le j<2,\quad 0\le a<3,\quad 0\le b<2,
\]
span \(k[P_0,Q_0]\).
Indeed \(j+2a+6b\) enumerates the twelve reduced
\(Q_0\)-degrees, and \(W_0,Z_0\) are monic of
degrees two and six in \(Q_0\) over \(k[P_0]\).
The top \(v\)-degrees \(42j+34a+99b\) are distinct
modulo twelve, so no two of the displayed monomials
have the same degree.
The Apéry representatives modulo twelve are
\[
      \{0,34,42,68,76,110\}
         \ \cup\ \bigl(99+\{0,34,42,68,76,110\}\bigr).
\]
They give conductor \(198\) and \(99\) gaps.
Equivalently each missing degree contributes one to
the dimension of \(k[v]/k[P_0,Q_0]\), giving total
affine \(\delta=99\) via normalization.
The transverse index \(M-1\) is the outer chart
calculation. In a nonzero simple characteristic chart
the degrees are \(2,7\), and the same monic-basis
argument with \(0\le j<2\) gives the semigroup
\(\langle2,7\rangle\) and its three gaps.

## The polar norm

**Corollary PlanarFortyTwoPolarNorm.**
Let \(f=\prod_{h(b)=0}(x-b)\) be monic and squarefree.
In the no-small-root profile,
\[
            \operatorname{Res}_X(P_X,Q_X)=C f^{1253},
                         \qquad C\in k^\times .
\]
On the normalized polar curve \(P_X=0\), counting
components with their multiplicities, the pole divisor
of \(P_x\) over each finite high position has degree
\(23\). In the unsplit outer chart the total intersection
multiplicity of \(P_v,Q_v\) on the outer boundary is
\(9H/2\).

On the polar, \(P_x Q_X=\kappa\), so there are no
finite affine poles or zeros of \(P_x\).
At a high position the two simple characteristic charts
give units \(P_x\), by their nonzero critical slopes.
The two other roots of \(r_b'\) off \(r_b=0\),
counted with multiplicity, have critical \(P\)-values
asymptotic to nonzero multiples of \((x-b)^{-6}\).
They contribute \(14\) to the pole degree of \(P_x\):
its restriction is the derivative of the critical
\(P\)-value in \(x\).
At the seven critical points in the scalar chart the
specialization of \(w^2P\), with \(x-b=w^7\),
has nonzero values. The seven roots of its derivative
form one deck orbit, from the
[generic position monodromy calculation](POLYNOMIAL_WEYL_ORDER_TWELVE_FORTY_TWO_PLANAR_FIBRES.md#monodromy-of-the-generic-position-projection).
The critical \(P\)-value has pole order two in \(w\),
so \(P_x\) has pole order nine there.
This gives \(23=14+9\) and accounts for all polar sheets
at the finite high positions.
At other finite positions the leading coefficient
of \(P_X\) is a unit, so no pole or zero is lost there.

The norm of \(P_x\) over \(k(x)\) on the polar thus
equals \(C_0 f^{-23}\), \(C_0\ne0\).
The norm of \(Q_X\) is its inverse up to a constant.
The leading coefficient of \(P_X\) in \(X\) is
a scalar times \(f^{30}\) and \(\deg_X Q_X=41\).
The resultant is that leading coefficient to power
\(41\) times the norm of \(Q_X\), giving
\(30\cdot41+23=1253\).

At the unsplit outer boundary, with \(t=1/x\),
the polar has eleven sheets, counted with multiplicity,
from \(P_v(0,v)\) of degree eleven.
On it \(Q_X=t^{-M}Q_v\), so \(P_x=\kappa t^M/Q_v\).
The total outer order of \(P_x\) is \(23H\) by
the polar norm. Subtracting from \(11M\) gives
\(\sum\operatorname{ord}Q_v=11M-23H=9H/2\).
The sum on the normalized polar, with component
multiplicities, is the intersection multiplicity
of the two derivatives at the boundary.

The [residue-monodromy theorem](POLYNOMIAL_WEYL_ORDER_TWELVE_FORTY_TWO_RESIDUE_MONODROMY.md)
classifies the square–cube residue factors and proves that this unsplit
profile has full symmetric or alternating geometric monodromy,
according to the parity of the number of pairs of high roots.

The [position double locus](POLYNOMIAL_WEYL_ORDER_TWELVE_FORTY_TWO_POSITION_DOUBLE_LOCUS.md)
uses the polar resultant at zero momentum difference. Its two resultant
leading coefficients determine the norm of that difference and the
saturation obstruction for the graph lift.
