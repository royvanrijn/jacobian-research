# Uniform sphere moment tails and effective safety loci

Let \(r\ge2\), let \(X\subset\mathbb A^r\) be the smooth affine quadric
\(q=\sum_i x_i^2-1=0\), and let \(\Lambda\) be the normalized mean on
its real sphere. For any field of characteristic zero, define the
same functional by its rational monomial moments. A monomial with an
odd exponent has mean zero, and
\[
 \Lambda(x^{2\alpha})=
 \frac{\prod_i(2\alpha_i-1)!!}{\prod_{j=0}^{|\alpha|-1}(r+2j)}.
\]
For a phase \(f\) and a marker \(g\), write
\(M_m(g,f)=\Lambda(gf^m)\). Call \(f\) **safe** if for every polynomial
marker \(g\) the sequence \(M_m(g,f)\) is eventually zero. This condition
includes eventual vanishing for the marker \(1\).

**Theorem.** For fixed \(r,d\) there are computable integers \(A\ge1\)
and \(B\ge0\) such that for every phase of degree at most \(d\) and
marker of degree at most \(e\), eventual vanishing is equivalent to
\(M_m(g,f)=0\) for all \(m\ge A+Be\). For each \(r,d,e\) there is a
computable \(U\ge A+Be\) such that it suffices to test
\(A+Be\le m\le U\). The safety locus in the phase coefficient space
is Zariski closed. An algorithm given \(r,d\) produces a finite list of
single-moment polynomial equations defining this locus. All these
assertions apply over arbitrary fields of characteristic zero.

The onset bound comes from cancelling vertical poles on a proper
model. The finite window and the finite equations use generic-field
moment algorithms and strict enlargement of radical ideals.

## 1. Clearing poles over the zero phase value

Work first over \(\mathbb C\). The Poincaré residue form
\(\omega=\operatorname{Res}_X(dx_1\wedge\cdots\wedge dx_r/q)\) restricts
to a nonzero constant multiple of the standard oriented volume form.
For nonconstant \(f\), start with the smooth projective quadric
compactifying \(X\). Resolve the rational map to \(\mathbb P^1\) given
by \(f\) by principalizing its base ideal. Then resolve the reduced
union of the preimage of the boundary and the zero fiber. On restricting
to the preimage of \(\mathbb A^1\), this gives a proper morphism
\(f:Y\to\mathbb A^1\) with \(Y\) smooth and a morphism to the original
projective quadric. Write \(X'\) for the preimage of \(X\). The reduced
boundary \(Y\setminus X'\) and the reduced zero fiber have simple normal
crossings together. These constructions use principalization on a
smooth variety over characteristic zero. [Bierstone, Grigoriev, Milman
and Włodarczyk, Theorems 1.0.1–1.0.2](https://arxiv.org/pdf/1206.3090).

Let \(H\) be the union of the irreducible boundary components that
are not contained in the zero fiber. Components mapping to a value
other than zero are allowed in \(H\). At any point over zero the
branches of \(H\) are distinct from all the local factors of \(f\).
Choose \(a\ge1\) so that \(f^a\omega\) has no poles on divisors over
zero and \(b\ge0\) so that every \(f^b x_i\) has no poles on those
divisors. Each such pole order is finite and multiplication by \(f\)
strictly increases valuation along every component of the zero fiber.
All other poles of these forms and functions can only occur on \(H\).
If \(\deg g\le e\), then \(f^{be}g\) is a polynomial in \(f^b\) and
\(f^b x_i\), by using \((f^b)^{e-|\alpha|}\prod_i(f^b x_i)^{\alpha_i}\)
for each monomial. Hence \(f^{a+be}g\omega\) is a top form with poles
only on \(H\).

The following version of the
[proper-volume obstruction](../../../unpublished/primitive-balanced-image-spectrum/proper-volume.tex)
applies to these forms. If a top form \(\eta\) on \(Y\) has poles only
on \(H\), every polynomial exponential period of \(\eta\) on a
compact closed cycle in \(Y\setminus H\) vanishes. The argument uses
the following properties of the meromorphic graph module.

- Use the left graph module \(i_{f+}\mathcal O_Y(*H)\) along the
  target coordinate \(s=0\), with the convention that
  \(\partial_s s-\alpha\) is nilpotent on \(\operatorname{gr}_V^\alpha\).
  Locally over zero one can write \(f=u\prod_i z_i^{a_i}\) with
  \(a_i>0\), and \(H=\{\prod_j v_j=0\}\), using disjoint sets of
  normal-crossing coordinates. Rescaling \(f\) by the unit \(u^{-1}\)
  transports the graph filtration and changes the graph generator by
  multiplication by a unit. For the monomial phase and local coefficient
  \(v^{-\beta}=\prod_j v_j^{-\beta_j}\), the operator
  \(\prod_i\partial_{z_i}^{a_i}\) on \(v^{-\beta}f^{s+1}\) gives
  \(v^{-\beta}f^s\) times
  \(\prod_i a_i^{a_i}\prod_{j=1}^{a_i}(s+j/a_i)\). Its roots are negative.
  The coefficient \(\prod_j v_j^{-1}\) generates \(\mathcal O_Y(*H)\)
  here under differential operators, so the rational roots of its
  b-function give existence of the local filtration. The negativity
  for each \(\beta\) puts \(v^{-\beta}\delta(s-f)\) in \(V^{>0}\).
  Off the zero fiber the filtration is the entire graph module, and
  uniqueness glues the local filtrations.
  Regular coefficients preserve the filtration. A finite cover of the
  zero fiber gives a common positive rational index for
  \(\eta\delta(s-f)\). [Chen, Dirks and Mustaţă, Remarks 5.9–5.10,
  Theorem 5.11 and Corollary 5.15](https://arxiv.org/pdf/2402.07816v1).
- Apply filtered direct image along the smooth projection of the graph
  module. The projection is proper on its support because \(f\) is
  proper. The image class of the top form has the same positive
  filtration bound. With Fourier convention \(\tau=\partial_s\) and
  \(\partial_\tau=-s\), the direct image is represented by the twisted
  de Rham complex with forms meromorphic along \(H\). Its restriction
  to \(Y\setminus H\) computes the period pairing with a compact cycle
  there. In particular the restriction sends any relation on the cyclic
  class of \(\eta\) to a differential relation on its exponential
  period. This is the same specialization and Čech de Rham pairing as
  for a regular form, applied after restriction. [Chen, Dirks and
  Mustaţă, Theorem 8.1](https://arxiv.org/pdf/2402.07816v1).
- A nonzero polynomial period would give a surjection from the Fourier
  transform of the cyclic class to the polynomial Weyl module. Inverse
  Fourier transform gives a quotient to the point module at zero sending
  the class nontrivially. Strictness for the V-filtration precludes this
  because the point module has zero filtration piece at every positive
  index. [Chen, Dirks and Mustaţă, Example 3.18, Proposition 3.19 and
  Corollary 3.20](https://arxiv.org/pdf/2402.07816v1).

The construction allows centers inside \(X\). The induced map
\(\pi:X'\to X\) is proper birational between smooth complex manifolds
and has degree one. Let \(\Gamma\) be the real sphere cycle. Pull back
its compactly supported Poincaré dual and take the inverse Poincaré dual
on \(X'\) to obtain a class \(\Gamma'\) in ordinary middle-dimensional
homology. The degree-one projection formula gives
\(\pi_*[\Gamma']=[\Gamma]\). Concretely, pairing with any ordinary
middle-degree cohomology class on \(X\) reduces this identity to
integration of the pullback of a compactly supported top-degree class.
The class \(\Gamma'\) has a compact cycle representative in \(X'\), and
pullback preserves the periods of the closed exponential top forms.

If \(M_m(g,f)\) is eventually zero, the exponential period of
\(g\omega\) is a polynomial. Its derivative of order \(a+be\) is, up to
sign, the period for \(f^{a+be}g\omega\). The preceding obstruction
forces that derivative to vanish. Equivalently, all moments from
\(a+be\) onward vanish. For constant phases, onset \(1\) suffices since
\(M_m(g,f)=f^m\Lambda(g)\).

## 2. Uniformity in the phase and an effective search

Consider the universal family on the phase coefficient space
\(S=\mathbb A^{\binom{r+d}{d}}_{\mathbb Q}\). The phases constant on the
quadric form a linear subspace, obtained by reduction modulo \(q\).
On each remaining irreducible component of the parameter space, perform
the construction of §1 on the generic fiber. The finite blowup sequence,
the morphism to \(\mathbb P^1\), and the pole orders \(a,b\) are defined
over the parameter function field, allowing a finite extension if
needed to split divisors and choose coordinates. Spread out to a
nonempty open set, using a finite dominant cover when a finite extension
was used and shrinking so it is surjective there. In particular the
following data spread out after shrinking.

- The centers are closed and relatively smooth. The blowup charts of
  a relatively smooth center in a relatively smooth ambient scheme
  commute with passage to the fibers. This can be checked in relative
  coordinates in which the center is cut out by coordinate functions.
- The rational pencil extends to the stated morphism on the model. A
  finite collection of charts covers the preimage of the zero value.
  On these charts its equation is a unit times the indicated powers
  of vertical factors. The pulled-back equation for the original
  boundary is a unit times powers of the indicated boundary factors,
  including any vertical factors which belong to the boundary. The
  combined distinct factors have independent relative differentials
  whenever they vanish together.
- On the same charts \(f^a\omega\) and all \(f^b x_i\) have regular
  expressions after multiplication by powers of the boundary factors
  other than the vertical factors. These are equalities of finitely
  many rational forms and functions. Projectivity, or constructibility
  of the image of the excluded loci, allows the open to be shrunk so
  the charts still cover the entire zero fiber in every fiber.

The remaining parameter locus has smaller dimension on each irreducible
component. Noetherian induction gives finitely many such pieces and
their models. Take \(A,B\) to be the maxima of the respective orders
on all pieces, with \(A\ge1\). Replacing \(a,b\) by these maxima
preserves the pole conditions. The construction uses only the coordinate
functions and unweighted volume form, so the same pair \(A,B\) works
for every marker degree. For an individual characteristic-zero field
point all its coefficient data lie in a finitely generated subfield,
which embeds into \(\mathbb C\). The rational monomial moment formula
then transfers the conclusion back to the field of definition.

Here is a terminating algorithm for \(A,B\) using the same existence
argument. Enumerate finite lists of finite-type parameter charts over
\(\mathbb Q\) mapping to \(S\), each with a nonnegative pair of pole
orders and the finite data just described. The charts can present
finite covers or ordinary locally closed parameter subsets. Each model
is specified as a finite sequence of blowups of explicitly given closed
centers in the original projective quadric over its parameter chart.
Use relative local Jacobian presentations for the smooth center and
ambient charts, with the center ideal generated by the indicated
members of a relative coordinate system. The usual blowup chart
equations then identify the fiberwise model. The last charts specify
the morphism to \(\mathbb P^1\), the unit and factors for \(f\), and the
two kinds of pole-cleared expressions. For a constant phase chart the
alternative certificate is the identity of the phase with its constant
restriction modulo \(q\).

The tests are finite algebraic ones. Ideal identities, including
gluing and equality of rational expressions, can be supplied as
polynomial identities after localization, with denominators explicitly
nonzero on the relevant charts. Relative smoothness, the independence
of the differentials of any subcollection of simultaneously vanishing
factors, agreement with the pencil on the chart open where it is
defined, and coverage of the model and zero fiber can all be expressed
by polynomial equalities, inequations and quantifiers in fiber
coordinates. The last test asks whether every parameter point has a
preimage in the union of the listed parameter charts. The entire test
is decidable by characteristic-zero algebraically closed field
elimination. These presentations can be enumerated with equations,
inequations and all identities encoded over \(\mathbb Q\). The finite
spread-out constructions above admit such presentations, including
finite extensions of generic fields through parameter covers. Therefore
the search terminates at an accepted list. Output the maxima of its
orders. This gives computability without a numerical complexity bound.
[Seidenberg, A New Decision Method for Elementary Algebra, §6(c)](https://www.jstor.org/stable/1969640).

## 3. Finding a finite test window for one marker box

For a finitely generated characteristic-zero coefficient field \(K\),
let \(T\) be a compact tube around the real sphere in
\(\mathbb C^r\setminus\{q=0\}\). The normalized period of the rational
form with coefficient \(g/[q(1-zf)]\) on \(T\) has Taylor series
\(J(z)=\sum_{m\ge0} M_m(g,f)z^m\). The regular-certificate
creative-telescoping algorithm works over \(K\) and computes a nonzero
operator \(L\in K[z]\langle D_z\rangle\) annihilating \(J\). The tube is
closed and the certificate has no pole on it for small \(z\), apart from
isolated values of the parameter where the conclusion follows after
clearing denominators. The normalization is the unweighted volume,
which is nonzero. [Bostan, Lairez and Salvy, Proposition 4 and Theorem
12](https://arxiv.org/pdf/1301.4313v2).

For any prescribed \(N\), subtract from \(J\) its Taylor polynomial
\(P\) of degree less than \(N\). For \(h=-L(P)\), the residual is
annihilated by \(L\) if \(h=0\), and by \((hD_z-h')L\) otherwise.
Given an annihilator \(\sum p_jD_z^j\), put
\(\rho=\min_j(\operatorname{ord}_0 p_j-j)\) over its nonzero terms.
The order of a nonzero power-series solution is a nonnegative integer
root of
\[
 \sum_{\operatorname{ord}_0 p_j-j=\rho}
           \operatorname{lc}_0(p_j)(n)_j.
\]
Here \(\operatorname{lc}_0\) is the coefficient of the lowest nonzero
power and \((n)_j\) is the falling factorial. Comparison of the first
possible term gives the formula, and the term with largest \(j\) makes
the displayed polynomial nonzero. Over a prime quotient function field,
integer roots can be found by clearing denominators and taking a nonzero
rational coordinate polynomial in standard monomials of the term
coefficients. Its finite set of integer roots is then tested in the
full polynomial. Computing moments through the largest possible order
decides whether \(J-P=0\) and, if it is nonzero, returns an index
\(m\ge N\) and a nonzero moment. This is an exact decision for a
prescribed tail. [Global sphere purity recognition, §§2–4 and
7](GLOBAL_SPHERE_PURITY_RECOGNITION.md).

On the joint coefficient box for \(f,g\), every \(M_m\) is a polynomial
over \(\mathbb Q\). Set \(N=A+Be\) and begin with the zero ideal in
the coordinate ring. At each iteration compute the minimal primes
of the current radical ideal. Apply the prescribed-tail decision in
the fraction field of each quotient. If a tail is nonzero, append an
index \(m\ge N\) returned there, by adjoining the universal moment
polynomial \(M_m\) and taking the radical. The polynomial is nonzero
modulo that prime, so the radical increases strictly. The algorithm
stops by the ascending chain condition. At termination all tail moments
vanish modulo every minimal prime and hence throughout the final
zero locus. The maximum index adjoined, or \(N\) if none was adjoined,
is a valid endpoint \(U\). A test of every moment from \(N\) to \(U\)
therefore decides the complete tail. By §2 this is also an eventual
vanishing test.

## 4. Deciding safety over a coefficient field

Fix \(f\) over a finitely generated characteristic-zero field \(K\).
The full integration module supplies finitely many monomial markers
whose periods generate those for all markers. Write \(D_{x,\tau}\)
for the polynomial Weyl algebra over \(K\), and put
\[
 M=\bigl(K[x,1/q]/K[x]\bigr)[\tau]e^{\tau f},\qquad
 \delta=[1/q]e^{\tau f}.
\]
The notation for the exponential denotes the usual connection twist.
The cyclic generator \(\delta\) has exact left ideal
\[
 I=(q,\ L_{ij}-\tau L_{ij}(f)\ (i<j),\ \partial_\tau-f),
 \qquad L_{ij}=x_i\partial_{x_j}-x_j\partial_{x_i}.
\]
To see exactness, first untwist by conjugating the Weyl action. On
\(q=0\) the charts \(x_j\ne0\) cover the support. In such a chart
\(q,(x_i)_{i\ne j}\) are étale coordinates and the \(L_{ij}\) span
the tangent derivatives at fixed \(q\). The quotient by \(q\), tangent
derivatives, and \(\partial_\tau\) has the usual local normal-derivative
filtration, with one copy of the structure sheaf of the hypersurface
in each order. Mapping to \([1/q]\) identifies the successive quotients
by \(\partial_q^j(1/q)=(-1)^j j!/q^{j+1}\). This proves local exactness
on the support. Off the support both modules vanish. Twisting back by
the polynomial exponential gives the displayed ideal. The module is
holonomic by the smooth hypersurface presentation and the twist.

The algebraic integration algorithm computes generators of the full
\(D_\tau\)-module \(M/\sum_i\partial_{x_i}M\). For this cyclic input
one may take the residue classes of \(x^\alpha\delta\) with
\(|\alpha|\le k_1\), where \(k_1\) is the largest integral root of the
integration b-function, or no classes if \(k_1<0\) or no root exists.
The required weight Gröbner and elimination calculations use the
coefficient field operations, and the integral roots are computed as
in §3. These computations commute with scalar extension, or can be
checked after a complex embedding by faithful flatness. The resulting
presentation is the full integration quotient.
[Oaku, §4, immediately before Definition 4 and Algorithm 3, steps 1–4](https://www.lab.twcu.ac.jp/oaku/Oaku2013.pdf).

Over a complex embedding of \(K\), the compact tube defines a
\(D_\tau\)-linear map from the integration quotient to exponential
period series. The period of \([v]\,e^{\tau f}\) uses the representative
\(v\,e^{\tau f}\,dx_1\wedge\cdots\wedge dx_r\) on the tube. If \(v\)
is regular across \(q=0\), the period vanishes by Stokes on the disk
tube. The \(\partial_{x_i}\) terms are exact on the boundary tube and
have zero period there. Differentiation in \(\tau\) commutes with
integration and realizes the twisted \(\partial_\tau\) action. After
the common nonzero normalization the image of \(g\delta\) is
\(\sum_{m\ge0}M_m(g,f)\tau^m/m!\). In particular these monomial period
series and their relations are defined over \(K\).

A marker has eventually vanishing moments precisely when its exponential
period is a polynomial. Polynomials are stable under \(D_\tau\), so
all markers have polynomial periods if and only if each of the
calculated monomial generators does. The test in §3 with the onset
from §2 decides this for each generator in finitely many exact
operations. If a generator fails, the same test returns that actual
monomial and an index at or beyond its onset with nonzero moment.
The theorem about the onset implies the marker has infinitely many
nonzero moments. This procedure applies to every phase, including
phases which already fail for the marker \(1\).

## 5. Computing equations for the safety locus

Testing all monomial markers is equivalent to testing all polynomial
markers, because a polynomial has finitely many terms. By §2, safety
is exactly the common zero locus of the rational-coefficient moment
polynomials
\[
 \{\,M_m(x^\alpha,f)\mid \alpha\in\mathbb N^r,
                           \ m\ge A+B|\alpha|\,\}.
\]
This proves closedness. To compute a finite list of these equations,
begin with the zero ideal in the phase coefficient ring. Take the
minimal primes of the current radical and apply §4 to each function
field phase. If one is unsafe, adjoin the single-moment equation
returned for it and take the radical. This strictly increases the
radical because the moment is nonzero in that quotient field. Every
safe point satisfies the adjoined equation. The ascending chain
condition gives termination, and when all the generic phases are
safe every moment in the displayed set belongs to every minimal
prime. Their intersection is the current radical. Thus every point
of the final zero locus is safe. This specialization uses polynomial
identities only. It works over all fields of characteristic zero.

In particular the algorithm computes a sufficient monomial marker
degree for a given phase box, by taking the largest degree in the
output equations, and also computes a finite set of moment indices to
test for those markers. The algorithm does not require finding a
rational point on a coefficient variety.

## 6. Bihomogeneous balanced Image

For a fixed number of paired variables and a fixed bidegree, put the
usual sphere phase on the coefficient space of balanced original
polynomials by \((z,\xi)=(x+iy,x-iy)\). The radial factorial identity
transfers eventual mixed vanishing for each balanced bihomogeneous
original marker to eventual sphere vanishing for its associated marker.
Every polynomial in the sphere coordinates decomposes into bidegrees
in \(z,\xi\). Terms with unequal degrees vanish against powers of the
balanced phase by circle invariance. Testing all balanced bihomogeneous
markers is therefore equivalent to testing all sphere markers.

For an original marker of unequal bidegree, terms with more \(\xi\)
than \(z\) contract to zero. For each remaining homogeneous output
degree, all scalar coefficients are tested by using
\(\partial_z^\gamma E(H)=E(\xi^\gamma H)\) with \(|\gamma|\) equal
to that degree. The resulting scalar tests are finite sums of balanced
bihomogeneous tests. This proves the equivalence for all original
polynomial markers in the definition of balanced Image safety. For the
coefficient equations, a sphere monomial marker in \(x,y\) becomes a
polynomial in \(z,\xi\) with rational coefficients times a common power
of \(i\), by \(x=(z+\xi)/2\) and \(y=(z-\xi)/(2i)\). It has at most one
balanced bidegree piece, since the total degree is fixed. The zero
test for its moment is therefore a rational-coefficient contraction
equation, by dropping the common power of \(i\) and the nonzero radial
factor. Thus the safety locus for balanced Image in a fixed paired
bidegree is closed and effectively defined over \(\mathbb Q\) by mixed
moments. [Global sphere purity recognition, §§6–7](GLOBAL_SPHERE_PURITY_RECOGNITION.md).

The proof gives an existence and termination result for every fixed
degree box. It supplies no numerical onset or safety-locus equations
for a specific box without running the corresponding algorithms.
