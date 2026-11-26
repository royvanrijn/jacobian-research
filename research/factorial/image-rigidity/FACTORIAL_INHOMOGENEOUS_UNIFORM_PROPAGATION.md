# Uniform propagation of zero factorial windows for polynomial phases

Let \(\Lambda(X^\alpha)=\prod_i\alpha_i!\) and write
\(M_j(F)=\Lambda(F^j)\).

**Theorem.** For each dimension \(m\ge1\) and total degree bound \(D\ge0\)
there are effective integers \(B\ge1\) and \(H\ge0\) with the following
properties for every phase \(F\) of degree at most \(D\) over a
characteristic-zero field.

- A block of \(B\) consecutive zero factorial moments starting at any
  positive order \(n\) forces \(M_j(F)=0\) for all \(j\ge n\).
- If the factorial moments are eventually zero, then \(M_j(F)=0\)
  for every \(j>H\).

The bounds are independent of the coefficients, support and starting
order. A family in which each nonzero phase has infinitely many nonzero
moments therefore has a uniform finite nonvanishing window. Its length
is bounded in terms of dimension and total degree.

The ordinary generating series in this argument is formal. Integration
by parts identifies it with a formal functional on an exponentially
twisted relative de Rham group. An integral structure for the Stokes
filtration of the top cohomology of each face supplies a countable set
of possible indicial roots. Specialization then gives a finite collection
of endpoint polynomials for the recurrences.

## The formal factorial functional

Take a coefficient field \(K\) of characteristic zero and let
\(\theta=t\partial_t\). Work over \(K(t)\). For
\(I\subseteq\{1,\ldots,m\}\) set the coordinates indexed by \(I\) equal
to zero, write \(F_I\) for the resulting polynomial and put
\[
 U_I=\operatorname{Spec}K(t)[X_j:j\notin I][(1-tF_I)^{-1}],
 \qquad
 S_I=\sum_{j\notin I}X_j .
\]
In each column of the relative affine de Rham complex use the
twisted differential
\[
                 d_I=d_X-d_XS_I\wedge .
\]
At column \(p\) take the direct sum of the complexes for \(|I|=p\).
For \(J=I\cup\{i\}\) with \(i\notin I\), the horizontal differential
\(\delta\) restricts from \(U_I\) to \(U_J\) with sign
\((-1)^{\#\{j\in I:j<i\}}\). The total differential is
\(d_{\mathrm{tot}}=\delta+(-1)^p d_I\) on column \(p\).
The restriction maps commute with the twisted differentials and
\(\theta\). The action of \(\theta\) on coefficients thus defines the
connection on total cohomology. At each affine face it is the
relative de Rham connection for \(d-dS_I\), whose exponential twist
is independent of \(t\).

The highest total degree is \(m\). Put
\(dX_{I^c}=\bigwedge_{j\notin I}dX_j\) in increasing index order. On
total degree \(m\) define
\[
 \mathcal A\bigl((h_I\,dX_{I^c})_I\bigr)=
   \sum_I (-1)^{\binom{|I|}{2}+\sum_{i\in I}(i-1)}\Lambda_I(h_I)
   \quad\in K((t)).
 \tag{1}
\]
Each \(h_I\) is expanded at \(t=0\), using the geometric expansion
of \((1-tF_I)^{-1}\) and Laurent expansions of scalars in \(K(t)\).
The Laurent coefficients are polynomials, so \(\Lambda_I\) acts
coefficientwise. It is the factorial functional on the surviving
coordinates, with the empty product equal to one.

The basic integration-by-parts identity, strictly on polynomials, is
\[
 \Lambda_I((\partial_i-1)h)
          =-\Lambda_{I\cup\{i\}}(h|_{X_i=0}).
 \tag{2}
\]
For \(h=X_i^k\) times a monomial in the other coordinates, the claim
follows from \(k(k-1)!=k!\) if \(k\ge1\) and from direct evaluation
if \(k=0\). The Laurent expansion carries it over coefficientwise.
For a degree-\((m-1)\) term \(h\,dX_{(I\cup\{i\})^c}\) in column
\(p=|I|\), the vertical contribution to (1), prior to applying (2),
has extra sign \((-1)^{p+i-1-q}\), where
\(q=\#\{j\in I:j<i\}\). The horizontal contribution has extra sign
\((-1)^q\). The weight in (1) changes by \((-1)^{p+i-1}\)
from \(I\) to \(I\cup\{i\}\), so the two contributions cancel by (2).
These terms span total degree \(m-1\). Thus \(\mathcal A\) kills
exact forms and still commutes with \(\theta\).

Let \(\omega\) be the class with sole component
\[
                \omega_{\varnothing}
               =\frac{dX_1\wedge\cdots\wedge dX_m}{1-tF}.
\]
Then \(\mathcal A(\omega)=\sum_{j\ge0} M_j(F)t^j\), whose constant
term is one. In particular \(\omega\) is a nonzero cohomology
class. These identities require only formal series and apply
even when the moment series has zero radius of convergence.

## Formal monodromy of the face complexes

Here is the geometric input over \(K=\mathbb C\). For each \(I\),
the top relative cohomology of \(U_I\) with twist \(d_X-d_XS_I\wedge\)
is finite-dimensional over \(\mathbb C(t)\). Its formal local
monodromy on the summand with exponential factor zero has
algebraic eigenvalues at \(t=0\) and \(t=\infty\).
We give the direct image and the hypotheses for this assertion.

Realize \(U_I\) over an open subset of the \(t\)-line by taking the
open set \(1-tF_I\ne0\) in
\(\mathbb A^{m-|I|}\times\mathbb A^1_t\). It is smooth. Compactify
into a smooth projective variety \(X\) over \(\mathbb C\) with
projective map \(f:X\to\mathbb P^1_t\), a strict normal crossing
boundary \(D_X\), and with \(g=-S_I\) extending to a morphism
\(X\to\mathbb P^1\). For this, first take the closure of the graph
of \(g\) in the product compactification and resolve the pair,
preserving the given smooth open set. The zero and pole divisors
of \(g\) are then disjoint and its poles are contained in the
boundary. For a zero-dimensional face the twist is trivial and
the assertion is immediate.

Choose a dense open \(V\subseteq\mathbb P^1_t\) over which \(f\) and
all horizontal boundary strata are smooth, excluding \(0\) and
\(\infty\). Include the fibers over the excluded points in the
boundary, resolving the divisor over those excluded points to
keep it strict normal crossing. These adjustments preserve the
generic fiber and achieve conditions (i)–(iii) in Sabbah's
Section 4.2. The regular meromorphic connection used there is
\((\mathcal O_X(*D_X),d)\), a rank-one torsion connection with
trivial local system. Use \(D_1=\varnothing\).

The differential module in question is the generic fiber of
\(H^0 f_+(\mathcal O_X(*D_X),d+dg)\). Indeed the localized
connection is the open direct image of the twist on \(X\setminus
D_X\). Over the generic base the morphism to the \(t\)-line is
smooth of relative dimension \(m-|I|\). In the left
\(\mathcal D\)-module direct image, the relative de Rham shift by
that dimension identifies \(H^0 f_+\) with the top relative
cohomology. The generic fiber of the open set is affine, so
the coordinate-ring de Rham complex computes it. Holonomic direct
image gives finite rank over the generic base.
[Roucairol, Introduction and Theorem 0.3](https://arxiv.org/pdf/math/0604134)

Sabbah's Lemma 4.10 applies to the right-hand direct image in his
display (3.1*) with exactly these geometric data. The lemma
requires a torsion input regular connection and conditions
(i)–(iii), independently of rigidity of a target connection.
It gives a Stokes-filtered local system defined over a cyclotomic
number field. By Definition 4.2 and Remark 4.3 the graded for
each exponential factor is a local system over that field after
the Stokes ramification. In the irregular Riemann–Hilbert
correspondence this graded is the local system of the associated
regular formal factor. Thus its monodromy has algebraic
eigenvalues. For the zero exponential factor, a positive power
of the unramified monodromy is the ramified monodromy, giving
algebraic eigenvalues also before the base change.
[Sabbah, Remark 1.2, Definition 4.2, Remark 4.3, and Lemma 4.10](https://ems.press/content/serial-article-files/31137)

For total degree \(m\), the column spectral sequence has associated
graded groups that are subquotients of
\(E_1^{p,m-p}\), the direct sums of these **top** face cohomologies.
It is a spectral sequence of differential modules.
Formal completion is exact, and the zero-exponential summand
of the formal decomposition is preserved by morphisms and exact
on subquotients. This follows by taking a common ramification
in the formal decomposition theorem and using the absence of
maps between distinct exponential factors. In a logarithmic
basis for the internal Hom of two regular factors, such a map
would be a nonzero Laurent vector \(v\) satisfying
\[
 \bigl(u\partial_u+A(u)+u\partial_u(\phi_2-\phi_1)\bigr)v=0
\]
in a ramified local parameter \(u\). Here \(A(u)\) has power
series entries and the nonzero polar difference
\(\phi_2-\phi_1\) has pole order \(h\ge1\).
The last term has valuation \(\operatorname{ord}_u(v)-h\),
strictly below the possible valuations of the other terms.
Thus the maps are block diagonal and exactness holds within
each block. The monodromy eigenvalues on the regular factor
of a subquotient or extension come from those on its
constituents. This proves the same algebraicity assertion for
the total cohomology in degree \(m\) and for every cyclic
differential submodule of it.
[Sabbah, Definition 1.1 and the following Levelt–Turrittin statement](https://ems.press/content/serial-article-files/31137)

## Endpoint roots of an annihilator

Let \(K\) be a finitely generated extension of \(\overline{\mathbb Q}\).
The affine face complexes and their total cohomology commute with
extension from \(K(t)\) to \(\mathbb C(t)\) under any embedding
\(K\hookrightarrow\mathbb C\) fixing \(\overline{\mathbb Q}\), by
flatness of field extension. In particular the cyclic rank of
\(\omega\) is finite. Take the first dependence of
\(\omega,\theta\omega,\theta^2\omega,\ldots\), clear its scalar
denominators and remove any common power of \(t\). This gives
\[
             L=\sum_{k=0}^b t^k P_k(\theta),
             \qquad P_0P_b\ne0,\qquad P_k\in K[T].
 \tag{3}
\]
The dependency and cyclic rank persist under the stated field
extension. Equation (1) shows that \(L\) annihilates the factorial
moment series.

To examine \(P_0\) after an embedding, choose a root \(\lambda\)
with no positive integer translate among its roots. The formal
series \(t^\lambda(1+\sum_{j>0}c_j t^j)\) solves (3), since its
recursion after the first coefficient divides by
\(P_0(\lambda+j)\) at order \(j\). The formal series and its
successive Euler derivatives define a nonzero horizontal map
from the formal cyclic module to the rank-one regular
module with Euler exponent \(\lambda\). In fact the initial
cyclic derivatives are a basis, their last relation is (3),
and the scalar solution makes that relation hold in the target.
Thus the map is a quotient and its monodromy eigenvalue
\(\exp(2\pi i\lambda)\), or its inverse under the dual convention,
is algebraic by the preceding section. Other roots of \(P_0\)
differ by an integer from a root with no further positive
translate, so their exponentials are also algebraic. At infinity
the change of variable \(u=1/t\) changes \(t^{-b}L\) to an operator
with leading endpoint polynomial \(P_b(-T)\). The same argument
there shows the algebraicity of \(\exp(2\pi i\lambda)\) for
every root of \(P_b\).

The set
\[
 \mathcal R=\{\lambda\in\mathbb C:
                         \exp(2\pi i\lambda)\in\overline{\mathbb Q}\}
\]
is countable. Every coefficient of each monic endpoint polynomial
maps under every embedding \(K\hookrightarrow\mathbb C\) to an
elementary symmetric polynomial of some finite tuple in
\(\mathcal R\). Its possible images belong to a fixed countable
set. If such a coefficient were transcendental over
\(\overline{\mathbb Q}\), it could be included in a transcendence
basis of \(K\), and it could be sent to any transcendental complex
number by extending the embedding of the rational function field
to the algebraic extension \(K\). Thus the two monic endpoint
polynomials already belong to \(\overline{\mathbb Q}[T]\).
Algebraicity of these coefficients is sufficient for the bounds.

## Finite specialization and propagation

Take an irreducible closed subvariety \(A\) of the space of
coefficients of the degree-\(\le D\) polynomials over
\(\overline{\mathbb Q}\), and take \(K=\overline{\mathbb Q}(A)\).
The generic relation \(L\omega=0\) is witnessed by a finite tuple
of rational forms in the total complex, with
\(L\omega=d_{\mathrm{tot}}\eta\) at the level of forms.
Restrict to a nonempty open subset of \(A\) on which the rational
functions in the parameters occurring in the relation and in
\(\eta\) are defined as elements of the specialized field in
\(t\). Also keep the scalar factors of the two monic endpoint
polynomials nonzero. This is possible because the coefficient
denominators need only stay nonzero as polynomials in \(t\), and
each endpoint polynomial is a scalar multiple of its monic
polynomial over \(\overline{\mathbb Q}\). The remaining
denominators in the coordinates are powers of \(1-tF_I\) as
in the affine complex. The specialized identity therefore holds
in that complex for all points of the open subset. Apply (1) to
get (3) for all the corresponding factorial moment sequences with
the same monic endpoints. Induction on the dimension of the
closed remainder gives a finite stratification with such
operators.

Let \(b_*\) be the largest span \(b\) in these operators. Choose
\(N_*\ge1\) strictly larger than every nonnegative integer root
of all the initial endpoint polynomials, and choose \(H\ge0\)
at least as large as every nonnegative integer root of the
terminal endpoints. When an endpoint has no nonnegative integer
root it puts no constraint on the corresponding choice.
The coefficients of the moment series satisfy
\[
 P_0(j) M_j(F)
       =-\sum_{k=1}^b P_k(j-k)M_{j-k}(F),\qquad M_i(F)=0\quad(i<0).
 \tag{4}
\]
For \(j\ge N_*\) the recurrence determines the next term.
A zero block of length \(B=N_*+b_*+1\) at any positive start
forces every subsequent term to vanish. If the series is a
polynomial of degree \(h\), its leading coefficient is nonzero
and the coefficient of \(t^{h+b}\) in its annihilator equation is
\(P_b(h)M_h(F)\), hence \(h\le H\).
For an arbitrary characteristic-zero field, embed the subfield
generated by the finitely many coefficients into \(\mathbb C\).
The factorial moment identities survive that embedding.

The bounds are effective. At the generic point of each stratum,
enumerate operators (3) and exactness witnesses
\(\eta\) over the computable function field \(K\).
The searches use rational form arithmetic and enumerate the two
monic endpoints separately over \(\overline{\mathbb Q}\).
Their existence proves termination. The nonempty open subset
can be determined from the finite denominators and scalar factors.
The complement decomposes effectively into irreducible components
of smaller dimension. The resulting finite set of algebraic
endpoint polynomials permits exact root and integer tests.

In particular eventual vanishing at fixed \(m,D\) is equivalent to
\[
             M_{H+1}(F)=\cdots=M_{H+B}(F)=0 .
 \tag{5}
\]
The equations are homogeneous polynomials in the coefficient
variables. Thus an exact projective ideal calculation decides
whether there is a nonzero eventually vanishing phase at any
fixed dimension and degree. The analogous test for vanishing
at all positive orders uses the first \(B\) equations. The proof
gives no small complexity estimate for either test.

## Nonvanishing consequences and scope

The window applies wherever eventual nonvanishing has an
independent proof. In particular Long's nonflat theorem gives
a coefficient-uniform all-start window for every binary phase
of degree at most \(D\) whose top homogeneous part restricts
to a nonconstant polynomial on \(X+Y=1\). The full binary conclusion using
Long's corollary also uses the separate exponential-integral
nonvanishing theorem specified in that corollary.
[Long, Theorem 1.1 and Corollary 1.3, preprint 2611.00813v1](https://hexagonmath.org/2611.00813v1)

The propagation proof uses the Stokes integral structure and
formal integration by parts. Its window depends on dimension and
degree. The support-size window and unrestricted factorial
nonvanishing remain open.
