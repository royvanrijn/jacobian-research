# Rational slices for signed weights in the plane

The leading-symbol obstructions for a Weyl pair and for a plane
Jacobian pair ask the same question for a weighted homogeneous
polynomial \(p\): does its Hamiltonian derivation have a
rational slice? For strictly positive weights the geometric
normal forms below are the classification of Abhyankar and
Oka. For mixed weights the signed toric chart gives a pole
bound for every multiplicity, and a single root-interpolation
test when the weighted degree is positive.

[Oka, Theorems 3.2 and 4.1 and §4.3](https://math.ac.vn/uploads/files/0702123.pdf)
gives the plane leading-face obstruction and the positive-weight
classification for homogeneous rational slices. His Assertion 4.2
and Lemmas 4.4–4.5 study local multiplicities, the necessary
crossing condition and a mixed case with equal multiplicities.
The [weighted Weyl slice theorem](POLYNOMIAL_WEYL_LEADING_FACTOR_RESIDUE.md)
supplies the same hypothesis for every nonconstant principal
symbol at weights of positive sum. The proof here uses the
same cyclic-cover exactness mechanism as the
[generic-fibre differential obstruction for Keller maps](../../jacobian/monodromy-theory/KELLER_LEADING_PAIR_DE_RHAM_OBSTRUCTION.md).

## A finite test at arbitrary signed weights

**Theorem WeightedPlaneHamiltonianSlicePoleTest.**
Let \(k\) have characteristic zero, let \((r,s)\) be primitive
integer weights with \(w=r+s>0\), and let \(p\in k[x,y]\)
be a nonconstant weighted homogeneous polynomial of degree \(d\).
A rational \(F\) with \(J(p,F)=1\) always exists if \(d=0\).
For \(d\ne0\), choose integers \(a,b\) with \(ar+bs=1\)
and put
\[
  z=x^a y^b,\quad t=y^r/x^s,\quad x=z^r t^{-b},\quad y=z^s t^a,
  \qquad p=z^d P(t).
\]
Write \(P(t)=\gamma t^\rho\prod_i(t-\lambda_i)^{m_i}\)
over the algebraic closure, with distinct \(\lambda_i\ne0\).
Put \(\beta=a-b,\ \nu=w/d,\ N=\sum_i m_i\).
Let \(q_0,q_\infty,q_i\) be the positive-integer parts of
the potential pole orders, in the following precise sense:
\[
\begin{aligned}
 q_0&=\beta-\nu\rho
                 &&\text{if this is a positive integer, and zero otherwise},\\
 q_\infty&=\nu(\rho+N)-\beta
                 &&\text{if this is a positive integer, and zero otherwise},\\
 q_i&=-\nu m_i   &&\text{if this is a positive integer, and zero otherwise}.
\end{aligned}
\]
Define \(B(t)=t^{q_0}\prod_i(t-\lambda_i)^{q_i}\).
Then a rational slice exists exactly when there is a polynomial
\(A(t)\) of degree at most \(\deg B+q_\infty\) satisfying
\[
 t\left(\frac AB\right)'
           +\left(\beta-\nu\,\frac{tP'}{P}\right)\frac AB=1.
 \tag{1}
\]
After clearing denominators this is one finite linear system.
In particular, \(d>0\) forces \(q_i=0\) and a simple zero of
\(A/B\) at every \(\lambda_i\). If \(\nu m_i=1\) for any \(i\)
there is no slice.

To prove this, the toric chart gives
\[
                 dx\wedge dy=z^{w-1}t^{\beta-1}\,dz\wedge dt .
\]
For \(d=0\) take
\(F=-z^w t^{\beta-1}/(wP'(t))\) and check \(dp\wedge dF=dx\wedge dy\).
For \(d\ne0\) use \(p,t\) as independent parameters. Over
\(k(p,t)\) the function field of the torus is the cyclic radical
extension \(z^d=p/P(t)\). Its degree is \(|d|\), since
\(Z^{|d|}-p/P(t)\) for positive \(d\), or
\(Z^{|d|}-P(t)/p\) for negative \(d\), is Eisenstein at \(p\)
or \(p^{-1}\), respectively. Extend constants first to an
algebraic closure to project onto the character of \(z^w\).
The differential on the generic level is
\[
                      \frac{z^w}{d p}t^{\beta-1}\,dt .
\]
Differentiation at fixed \(p\) commutes with the cyclic action.
Any rational primitive therefore has a character component
\(z^w t^\beta h(p,t)/(d p)\) whose derivative is that
differential. Substitution gives
\[
                     t h'+(\beta-\nu\,tP'/P)h=1 .
\]
The coefficients are independent of the level \(p\). A rational
solution \(h(p,t)\) therefore gives a solution in the rational
function field in \(t\) by specializing \(p\) away from the
finitely many values making a denominator identically zero.
Conversely any solution in \(t\) gives a slice by the same
formula for \(F\), homogeneous of weight \(w-d\).
The argument over an algebraic closure suffices for \(k\):
the finitely many algebraic coefficients of a slice lie in
a finite extension of \(k\), and the normalized field trace
over \(k(x,y)\) preserves \(J(p,F)=1\). This also allows
the linear system to be tested over a splitting field.

At any finite nonzero point outside the roots of \(P\), a pole
of \(h\) would give an uncancellable higher pole in \(t h'\).
At \(t=\lambda_i\) a nonzero leading term of order \(j\) gives
order \(j-1\) on the left unless \(j=\nu m_i\).
Thus a pole requires the indicated negative integral resonance.
If \(\nu>0\), the only possible order at this root is one:
orders at most zero give a pole, orders above one cannot
produce a nonzero constant, and at order one \(\nu m_i=1\)
would cancel that constant. At zero the indicial coefficient
is \(j+\beta-\nu\rho\) for leading order \(j\). At infinity,
for a leading term of degree \(j\) it is
\(j+\beta-\nu(\rho+N)\). A pole at either point therefore
requires the corresponding integral resonance, giving the
other two bounds. A rational \(h\) with these pole bounds
has exactly the form \(A/B\) of (1).
These bounds and the differential equation form a necessary
and sufficient test for the weighted symbol itself. They
supply necessary conditions face by face for a polynomial
or Weyl commutator partner.

## Strictly mixed weights and positive degree

**Corollary MixedWeightHamiltonianSliceRootInterpolation.**
Use the theorem's notation with \(r<0<s,\ r+s>0,\ d>0\),
and assume \(p\) has at least two monomials.
Write \(G(t)=t^{-\rho}P(t)\), so \(G(0)\ne0\).
A slice exists exactly when \(m=\beta-\nu\rho\) is a
positive integer and the numbers
\[
                    \frac{T_m(\lambda_i)}{\lambda_i^m},
         \qquad T_m(t)=\left[\bigl(G(t)/G(0)\bigr)^\nu\right]_{\le m}
 \tag{2}
\]
are all equal to the same nonzero constant. Here the bracket
truncates the binomial expansion at zero through degree \(m\).
In particular \(m\) is at least the number of distinct roots.
Exchanging the variables gives the symmetric test.

Indeed an exponent point \((i,j)\) of weight \(d\) contributes
\(t^{-bi+aj}\) to \(P\), and its endpoint exponent \(E\)
satisfies
\[
              \nu E-\beta=(j-i)/d .
\]
At the maximal exponent of \(t\) this quantity is less than
one, since \(s\ge2\) and for \(j-i>0\) we have
\(d=s(j-i)+(r+s)i\ge s(j-i)\).
Thus \(q_\infty=0\), and since \(\nu>0\) the \(q_i\) vanish.
A constant \(h\) cannot vanish at a root, so \(q_0=m>0\)
is necessary.

Substitute \(h=t^{-m}A\) into (1). Then
\[
                 G A'-\nu G'A=t^{m-1}G,\qquad \deg A\le m.
 \tag{3}
\]
The Taylor expansion at zero through order \(m\)
determines the candidate from its constant term \(A_0\):
\[
                         A=A_0 T_m+t^m/m .
\]
This follows on multiplying by \((G/G(0))^{-\nu}\)
and integrating the formal equation at zero.
The candidate works exactly when it vanishes at each
root of \(G\). Necessity follows from the pole test.
For sufficiency the residual in (3) vanishes at zero
to order at least \(m\) by construction. Its vanishing
at the roots has the full multiplicities of \(G\),
since \(G'\) loses at most one order at each root and
\(A\) vanishes there. Its degree is at most \(m+\deg G-1\),
so divisibility by \(t^m G\) makes it zero.
The common \(A_0\) exists exactly under (2).

For instance, in Oka's Example 4.1 take
\(p=y^{2n}(1+(x^3 y)^n)^{1+3n}\) at weights \((-1,3)\),
with \(t=(x^3y)^{-1}\).
Then \(m=n,\ \nu=1/(3n)\),
\(G=(1+t^n)^{1+3n}\), and the quotient in (2) is \(1/(3n)\)
at all \(n\) roots. This verifies an arbitrary-root family.
More generally, when \(m\) equals the number of distinct
roots, put \(D=\prod_i(t-\lambda_i)\). Equation (3) forces
\(A=\eta D\) and becomes the explicit condition
\[
                  \eta(1-\nu m_i)D'(\lambda_i)=\lambda_i^{m-1}
                      \quad\text{for every }i
\]
with one nonzero scalar \(\eta\). Conversely these equations
give (3) by interpolation in degree at most \(m-1\).
For an unequal-multiplicity example at weights \((-1,2)\),
take \(p=y(1-x^2y)(1-2x^2y)^4\).
Here \(m=2,\ \nu=1/2\), and \(\lambda_1=1,\lambda_2=2\)
have multiplicities one and four. The common \(\eta=-2\)
gives the slice
\[
                   F=-xy(1-x^2y)(1-2x^2y)/p.
\]
The equal-multiplicity polynomial
\(x^2y^2(1-(x^2y)^2)^3\) at the same weights fails:
\(m=3,\ T_m=1-3t^2/2\), and the values in (2) disagree
at the two roots \(1,-1\).


## Positive weights and compatibility of exposed edges

**Corollary WeightedHamiltonianSlicePositiveEdgeCompatibility.**
For primitive strictly positive weights, the weighted homogeneous
polynomials admitting rational slices are, up to a nonzero scalar,
the following forms over the algebraic closure:
\[
\begin{array}{ll}
 x^i y^j, & i\ne j,\\
 x^i(y-\lambda x^s)^j, & r=1,\ i\ne j,\\
 y^i(x-\lambda y^r)^j, & s=1,\ i\ne j,\\
 (y-\lambda x)^i(y-\mu x)^j, & r=s=1,\ i\ne j .
\end{array}
\]
The monomial exponents are nonnegative, the displayed
non-coordinate factors have nonzero slope parameters,
and in the last row \(\lambda\ne\mu\).
Consequently, a polynomial whose leading symbols at all
strictly positive weights admit rational slices has at most
one edge exposed by a strictly positive weight. This applies
to each member of a scalar-commutator Weyl pair and to each
member of a plane Keller pair.

For the classification, write the endpoint exponent points
as \((i_0,j_0),(i_1,j_1)\) corresponding to the minimum and
maximum powers of \(t\) in \(P\). The quantities
\((j_0-i_0)/d,(j_1-i_1)/d\) lie in \([-1,1]\)
because \(d=r i+s j\) and all exponents are nonnegative.
A negative integral resonance at zero is therefore only
\(-1\), requiring \(r=1,j_0=0\).
A positive integral resonance at infinity is only one,
requiring \(s=1,i_1=0\).
The candidate \(h\) is a Laurent polynomial with powers
among \(-1,0,1\), including each outer power only under
its resonance. It has at least one zero for each distinct
non-coordinate factor, and therefore at most two such zeros.
One zero requires either \(r=1,j_0=0\) or \(s=1,i_1=0\),
which gives the displayed monomial times triangular factor.
Two zeros require both resonances, hence \(r=s=1\) and
neither coordinate factor. With no zeros the original
form is monomial. Each listed form becomes a monomial
in polynomial coordinates of constant nonzero Jacobian.
For the monomial \(X^iY^j\), the expression
\(X^{1-i}Y^{1-j}/(i-j)\) is a slice when \(i\ne j\);
when \(i=j>0\), the generic-level form is a nonzero
constant times \(dX/X\), which has no rational primitive.
This recovers the Abhyankar–Oka positive-weight list.

Order the positive exposed edges by increasing \(x\)-exponent.
Their negative slopes strictly decrease. Any shallow edge
(with slope between minus one and zero) must reach \(y=0\)
at its right endpoint by the second row of the classification.
It cannot have a following positive exposed edge.
Nor can it have a preceding one: that edge would be
shallow too and would have to reach \(y=0\) first.
By symmetry a steep edge has no other positive exposed edge.
If there is neither a shallow nor a steep edge, the only
possible slope is minus one. This gives the compatibility
statement. The weighted Weyl slice theorem applies to each
member of a scalar-commutator pair.

For plane Keller pairs the analogue is Oka's Theorem 3.2.
The following formal proof also works at the signed weights
of positive sum used here. If \(J(F_1,F_2)=\kappa\ne0\),
rescale the two polynomials by their weighted degrees to
obtain \(F_{1,\epsilon}=p+O(\epsilon)\) and
\(J(F_{1,\epsilon},F_{2,\epsilon})=\kappa\epsilon^n\),
where \(n=\deg F_1+\deg F_2-r-s\ge0\) in the chosen weights.
If a coefficient below order \(n\) of the second series is
the first remaining coefficient, it Poisson commutes with
\(p\) and is algebraic over \(k(p)\) by the Jacobian criterion.
Lift its simple root of the minimal polynomial to the
formal series field, replacing \(p\) in the coefficients by
\(F_{1,\epsilon}\). The lifted series commutes with
\(F_{1,\epsilon}\) by applying its Hamiltonian derivation
to that polynomial equation and dividing by the invertible
derivative. Subtract the lift at the corresponding order.
After finitely many subtractions, the coefficient at order
\(n\) has bracket \(\kappa\) with \(p\). Thus the signed tests
apply to both settings, including weights supporting a
vertex and symbols of weighted degree zero or negative.


## Gluing the tests around a Newton polygon

Call an edge tested if a primitive normal of positive sum
exposes it. Call a vertex tested if it is the unique maximizer
in some direction of positive sum. The Newton polygon here
is the convex hull of the exponent support, with no added
origin. In a segment the same definition applies to its
endpoints and to the exposed segment.

**Theorem WeightedHamiltonianSliceBoundaryGluing.**
Over an algebraically closed field of characteristic zero,
fix a lattice Newton polygon in the nonnegative quadrant
and prescribe nonzero coefficients at its vertices.
It can be realized by a polynomial whose every nonconstant
leading symbol at weights of positive sum admits a rational
Hamiltonian slice exactly under the following conditions.

- Every tested vertex other than the origin has distinct
  coordinates.
- On every tested edge with the given endpoints there exists
  a weighted homogeneous face polynomial admitting a slice,
  with nonzero coefficients at the endpoints.

Each edge polynomial that meets the second condition can
be used after diagonal scaling of its variables and scalar
multiplication, independently for each edge. Its end
coefficients can be made equal to the prescribed vertex
coefficients. Coefficients away from the tested faces are
unrestricted within the polygon.

Indeed a nonconstant monomial has a rational slice precisely
when its two exponents are distinct, as proved above for the
positive list, and the same elementary formula or logarithmic
obstruction is independent of the weights. Thus the vertex
condition is necessary. The edge condition is also necessary
by taking actual face polynomials.

An edge normal is unique up to positive scaling on the
tested side. In its toric chart write the face polynomial
as \(z^d P(t)\), with the first and last powers of \(t\)
being \(\rho\) and \(\rho+N\), \(N>0\). Diagonal scaling in
\(x,y\) can multiply \(t\) by an arbitrary member of
\(k^\times\), since its exponent vector is primitive in
the unimodular toric change. In particular the ratio of
the last to the first coefficient of \(P(t)\) can be
multiplied by any \(N\)th power. Every desired ratio can
be obtained over the algebraic closure. A further scalar
multiplication sets the first coefficient, and hence also
the last coefficient, to its prescribed vertex value.
These operations preserve rational slice existence: the
Jacobian scales by a nonzero constant, which is absorbed
into the primitive. This argument does not require \(d\)
to be nonzero.

In a polygon two distinct edges meet only at a vertex, so
after this normalization they specify consistent coefficients
on all tested edges. Set all other vertex coefficients to
their prescribed values and take arbitrary remaining
coefficients supported in the polygon. The nonzero vertex
coefficients give exactly that polygon. Every face exposed
by a weight of positive sum is either one of the tested
edges or a tested vertex. Consequently the two conditions
are sufficient, and give the stated compatibility
description. For a fixed polynomial the finite pole test
and root interpolation above decide each actual edge,
while this gluing result describes existence with prescribed
vertex data.

**Corollary WeightedHamiltonianSliceArbitrarilyManyMixedEdges.**
For every \(M\ge1\), the polynomial
\[
                  f_M=\sum_{j=2}^{M+2}x^{j(j+1)/2}y^j
\]
has exactly \(M\) edges exposed by weights of positive sum,
all strictly mixed. Every leading symbol at such weights
admits a rational Hamiltonian slice.

Put \(I_j=j(j+1)/2\). The successive upper-edge vectors are
\((j+1,1)\), whose slopes strictly decrease. Each outward
normal \((-1,j+1)\) has positive sum. When there is a
two-dimensional polygon, the closing chord has positive
slope less than one and its outward normal on the lower
side has negative sum. In the segment case the opposite
normal likewise has negative sum. Hence the only edge
conditions are at these \(M\) upper edges, and every vertex
is off the diagonal, since \(I_j>j\) for \(j\ge2\).

On the edge from \(j\) to \(j+1\) the face polynomial and a
slice are
\[
     p_j=x^{I_j}y^j(1+x^{j+1}y),\qquad
           F_j=\frac{x^{1-I_j}y^{1-j}}{I_j-j}.
\]
The first monomial of \(p_j\) contributes one to
\(J(p_j,F_j)\). The second contributes zero: its
determinant factor is \(2I_j-j(j+1)=0\).
Thus all the edge and vertex tests pass. In particular
the number of mixed edges cannot be bounded using slice
existence on the faces alone.

For one root, the mixed interpolation test has a simpler
equivalent condition. In the notation of (3), if
\(G=\gamma(t-\lambda)^N\) then a slice exists precisely when
\(m\) is a positive integer and
\(\nu N\notin\{1,\ldots,m\}\).
Indeed (3) becomes
\[
           (t-\lambda)A'-\nu N A=(t-\lambda)t^{m-1}.
\]
The operator on the left is diagonal on powers of
\(t-\lambda\), and the right side has nonzero coefficients
at all powers one through \(m\).
For the displayed family \(m=1\) and \(\nu N=2/(j+1)\),
which also verifies the slice condition by interpolation.


## The differential of the full polynomial

**Corollary WeightedHamiltonianSliceMixedChainGlobalObstruction.**
For \(M\ge2\), the polynomial \(f_M\) has no rational
Hamiltonian slice, although every nonconstant leading symbol
at weights of positive sum has one.

Write \(N=M+2\) and let \(p=f_M\), \(\eta\) a transcendental
level, over the algebraic closure of \(k(\eta)\).
The affine fibre \(p=\eta\) is smooth by generic smoothness
in characteristic zero. It lies in the torus because \(p\)
has \(x^3y^2\) as a factor. The relative differential
\(\omega\) on the fibre is regular there, given up to sign by
\(dx/p_y\) or \(dy/p_x\) where the respective denominator
is nonzero.

For the punctures of its complete normalization, the convex
hull for \(p-\eta\) has vertices the origin and
\((I_j,j)\) for \(2\le j\le N\). Each puncture lies over
an exposed edge of this hull. To see this, at any puncture
at least one of \(x,y\) has a pole, since otherwise the place
would have a center on the affine fibre.
Take the primitive direction proportional to the negatives
of their two orders. The relation \(p-\eta=0\) cannot have
a unique term of minimal order, so the direction exposes
an edge. Use the toric variables \(z,t\) of the pole test
for this direction and set \(q=1/z\).
Then \(q\) vanishes and \(t\) is a unit.
The equation is
\[
                       q^d(p-\eta)=\Phi(q,t)=0
\]
where \(d\) is the height of that edge and
\(\Phi(0,t)\) is its face polynomial in the toric chart.
Every root of this face polynomial in the torus is simple:
all the edges are binomial, including those from the origin.
Each root therefore gives a unique local branch with
uniformizer \(q\), and all punctures arise in this way.
Indeed the implicit function theorem gives the completed
local ring at a root as the power series ring in \(q\).

Writing \(w=r+s\) in this chart, one has
\[
              \operatorname{ord}(\omega)=d-w-1,\qquad
              \operatorname{ord}(x\omega)=d-w-1-r .
\]
The formulas follow from the volume
\(z^{w-1}t^{\beta-1}\,dz\wedge dt\),
\(p_t=z^d\Phi_t\), and the fact that \(\Phi_t\) is a unit
on the branch. On each middle edge
\((r,s)=(-1,j+1)\), \(d=I_j>w=j\), so both forms are regular.
On the edge from the origin to \((3,2)\) the direction
is \((-2,3)\), with \(d=0\).
The vector \((3,2)\) is primitive, so this edge gives
one puncture, where \(\omega\) has order minus two and
\(x\omega\) is regular.
On the last edge from \((I_N,N)\) to the origin the
direction is \((N/g,-I_N/g)\), \(g=\gcd(I_N,N)\), again
with \(d=0\). The order of \(x\omega\) is
\[
                          (I_N-2N)/g-1\ \ge\ 0
\]
because \(N\ge4\) and the difference divided by \(g\) is a
positive integer. The form \(\omega\) is regular there as
well. Thus \(x\omega\) is a nonzero regular differential
on the component containing the sole double pole of
\(\omega\). This component has positive genus.

If a rational primitive of \(\omega\) existed on this
component, it would have just one pole, of order one,
since differentiation raises pole order by one in
characteristic zero. Such a function is a degree-one map
to the projective line, impossible in positive genus.
A rational Hamiltonian slice would restrict to just such
a primitive. This proves the obstruction for the full
polynomial.

For \(M=1\), the global rational slice is \(1/(x^2y)\), by the
monomial bracket formula.
