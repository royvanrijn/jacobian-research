# A geometric lower bound of sixteen for ternary cubic moment prefixes

Let the original Gaussian Wick functional over a characteristic-zero field
have coordinates \(Z,W,Y\), with \(E[ZW]=E[Y^2]=1\) and all other basic
pairings zero. Let \(V=K[Z,W,Y]_{\leq3}\), of dimension twenty, and put
\(\mu_j(P)=E[P^j]\). A phase is pure if \(\mu_j(P)=0\) for every
\(j\geq1\). The coefficient space includes all polynomials of total
degree at most three; it is not the space of homogeneous ternary cubics.

**Theorem.** There is a ternary polynomial with algebraic coefficients,
of total degree at most three, such that

\[
 \mu_1(P)=\cdots=\mu_{15}(P)=0,
 \qquad \mu_m(P)\ne0\quad\hbox{for some }m>15.
 \tag{1}
\]

Consequently the least universal geometric pure-moment prefix is at
least sixteen, over every algebraically closed field of characteristic
zero. In particular neither ten nor fifteen moments suffice. This is an
existence theorem: it does not specify the nonzero later moment, assert
\(\mu_{16}(P)\ne0\), or give an upper bound for a universal prefix.

The argument uses the retained ternary-cubic pure-phase classification.
It is not a consequence of finite moment experiments. The theorem is an
internal deduction with the same mathematical-review limits as that
classification; independent whole-proof review and formal verification
are not asserted.

## 1. The precise classification input

The [Gaussian boundary manuscript](../../../results/gaussian-boundary/manuscript.tex),
Proposition `prop:frameextract`, Proposition `prop:wlinear`, and their
assembly in Section `sec:completion`, prove the following statement over
\(\mathbb C\): every pure ternary phase of total degree at most three
is, after a complex orthogonal linear coordinate change and nonzero
scalar multiplication, either strictly one-sided in circular weight or
equal to

\[
 P_*=W+ZY-\tfrac12 Z^3.                                      \tag{2}
\]

The lower-degree part uses the manuscript's quadratic classification.
Here strictly one-sided means that all monomials have strictly positive,
or all strictly negative, weight for
\(\operatorname{wt}(Z,W,Y)=(1,-1,0)\). No translation of the spatial
coordinates or addition of a constant is among these equivalences.

The exceptional phase is itself pure. One direct all-order justification
is the Gaussian-preserving polynomial shear

\[
 \phi(Z)=Z,\qquad \phi(Y)=Y-Z^2,\qquad
 \phi(W)=W+ZY-Z^3/2.
\]

Its infinitesimal derivation
\(Z(-Z\partial_Y+Y\partial_W)\) is locally nilpotent, has zero
ordinary divergence, and annihilates \(2ZW+Y^2\). Gaussian integration
by parts therefore gives \(E[\phi(A)]=E[A]\) for every polynomial
\(A\). Since \(\phi(W)=P_*\) and \(E[W^j]=0\) for \(j\geq1\),
this proves purity without testing a moment prefix. The same proof is
retained in
[the mixed-cutoff and moment-jet note, Section 3](MIXED_CUTOFF_JETS_AND_INFINITE_TAILS.md#3-the-exceptional-shear-exact-original-images-and-sharp-jet-order).

## 2. A polynomial separates the exceptional and one-sided cases

Write a phase in its ordinary homogeneous spatial pieces as
\(P=F_3+F_2+F_1+F_0\). Define the polynomial in its coefficients

\[
 \mathcal J(P)=E[F_1^3 F_3].                                \tag{3}
\]

This definition retains the original Gaussian functional. Linear
orthogonal substitutions preserve both spatial homogeneous degree and
the Wick functional, so \(\mathcal J\) is invariant under those
substitutions. Also \(\mathcal J(cP)=c^4\mathcal J(P)\).

If \(P\) is strictly one-sided, then every monomial of \(F_1^3F_3\)
has the same strict sign of circular weight. Every such monomial has
Gaussian mean zero. This includes a zero homogeneous piece, and proves
\(\mathcal J(P)=0\) for every phase orthogonally equivalent, up to a
scalar, to a one-sided phase. Thus even the Zariski closure of that
one-sided family lies in \(\{\mathcal J=0\}\). No invariant-theory
nullcone theorem or closure assumption is needed.

For (2), the first and third homogeneous pieces are \(W\) and
\(-Z^3/2\), respectively. Hence

\[
 \mathcal J(P_*)=-\tfrac12 E[W^3Z^3]=-\tfrac12 3!=-3.        \tag{4}
\]

Put \(U=\{\mathcal J\ne0\}\), a principal open neighbourhood of
\(P_*\). If \(\mathcal S\) denotes the pure locus, the classification
and (3)--(4) give an equality of sets of complex points

\[
 \mathcal S\cap U=
 (\mathrm O(2ZW+Y^2)\times\mathbb G_m)\cdot P_*.
 \tag{5}
\]

Every point of the orbit on the right is pure and remains in \(U\).
Conversely the one-sided alternative is excluded throughout \(U\).

## 3. The orbit has dimension four

The complex orthogonal group of a nonsingular ternary form has
dimension three. For an elementary dimension count, an orthonormal
first vector varies on a smooth two-dimensional quadric, a second
orthonormal vector in its nonsingular orthogonal plane varies on a
one-dimensional quadric, and the final unit vector has two choices.
Thus the orbit in (5), and its Zariski closure, have dimension at most
\(3+1=4\): the dimension of an algebraic image cannot exceed that of
its source.

For completeness the upper bound is attained. The three infinitesimal
orthogonal operators

\[
 D_0=Z\partial_Z-W\partial_W,\qquad
 D_+=Z\partial_Y-Y\partial_W,\qquad
 D_-=W\partial_Y-Y\partial_Z
\]

annihilate \(2ZW+Y^2\). Together with scalar multiplication their
orbit tangents at \(P_*\) are

\[
\begin{aligned}
 P_*&=W+ZY-Z^3/2,\\
 D_0P_*&=-W+ZY-3Z^3/2,\\
 D_+P_*&=Z^2-Y,\\
 D_-P_*&=ZW-Y^2+3Z^2Y/2.
\end{aligned}                                               \tag{6}
\]

They are linearly independent. The third and fourth have unique
\(Z^2\) and \(Z^2Y\) coefficients; the first two are independent
already on their \(W,ZY\) coefficients. Therefore (5) is exactly
four-dimensional. Only the upper bound is needed for (1).

The set \(\mathcal S\) is Zariski closed: it is the intersection of
the closed moment zero sets, or equivalently the zero set of the ideal
generated by all \(\mu_j\). Thus (5) also identifies the reduced pure
locus on this open set with the exceptional orbit.

## 4. Fifteen equations leave too much dimension

Work over \(\mathbb C\) first. At the point \(P_*\), the local ring
\(R=\mathcal O_{V,P_*}\) is a regular local ring of dimension twenty.
All the \(\mu_j\) lie in its maximal ideal by purity of \(P_*\).
For every \(N\geq0\), successive application of the local dimension
inequality gives

\[
 \dim R/(\mu_1,\ldots,\mu_N)\ \geq\ 20-N.                 \tag{7}
\]

The precise classical input is
[Stacks Project, Lemma 10.60.13](https://stacks.math.columbia.edu/tag/00KW):
if \(A\) is Noetherian local and \(x\) is in its maximal ideal, then
\(\dim A\leq\dim A/xA+1\). Iterating this statement proves (7).
The equivalent height bound for minimal primes over an ideal with
\(N\) generators is
[Lemma 10.60.12](https://stacks.math.columbia.edu/tag/0BBZ).
The actual statements and proofs of both lemmas were inspected for
this application. Each successive quotient is Noetherian local and
each next moment is still in its maximal ideal.

For \(N=15\), equation (7) gives local dimension at least five. The
open condition \(\mathcal J\ne0\) does not change this local ring.
If every complex point with these fifteen vanishing moments in that
open neighbourhood were pure, its reduced zero set would be contained
in the orbit closure of (5), of dimension four. That contradicts (7).
More explicitly, choose an irreducible component of the fifteen-moment
zero set through \(P_*\) with dimension at least five. Its intersection
with \(U\) minus the four-dimensional orbit closure is nonempty and
contains a closed complex point. That point has its first fifteen
moments zero, is outside the exceptional orbit, and cannot be one-sided
because \(\mathcal J\ne0\). The classification makes it nonpure.
This proves the complex form of (1). No assertion that the sixteenth
moment is nonzero follows from the dimension argument.

To obtain algebraic coefficients, choose such a complex point and one
index \(m>15\) with \(\mu_m\ne0\). The finite polynomial system

\[
 \mu_1=\cdots=\mu_{15}=0,
 \qquad a\mathcal J-1=0,
 \qquad b\mu_m-1=0                                      \tag{8}
\]

has rational coefficients and a complex solution. Its defining ideal
is proper over \(\overline{\mathbb Q}\), since a polynomial identity
expressing one as a combination of its generators would remain such
an identity after extending to \(\mathbb C\). The weak Nullstellensatz
therefore supplies a solution over \(\overline{\mathbb Q}\).
Any algebraically closed characteristic-zero field contains a copy of
\(\overline{\mathbb Q}\), so the same algebraic witness proves the
existence claim there. This descent does not presume an extension of
the analytic classification to arbitrary enormous coefficient fields.

The Nullstellensatz input is
[Stacks Project, Theorem 10.34.1](https://stacks.math.columbia.edu/tag/00FS).

## 5. A related scheme-theoretic consequence

This paragraph also uses the separately proved exact Jacobian statement
in the mixed-cutoff note, Section 3: at \(P_*\), the differentials of
\(\mu_1,\ldots,\mu_{10}\) are independent and every subsequent moment
differential vanishes. It concerns all twenty original coefficient
directions. Hence the Zariski tangent space of the full moment scheme
at \(P_*\) has dimension ten.

By (5)--(6) its reduced locus near \(P_*\) is the smooth exceptional
orbit of dimension four. Algebraic group orbits over \(\mathbb C\)
are smooth; here the orbit dimension and tangent dimension are both
four directly from (6). Thus the local full moment ideal is nonradical:
otherwise its tangent space would be that of this reduced orbit and
would have dimension four, contradicting the proved dimension ten.
The same conclusion holds at every point of the exceptional orbit by
orthogonal and scalar transport. This does not determine the full
nonreduced structure or its nilpotence index.

## Evidence and remaining scope

[The exact checker](verify_gaussian_cubic_prefix_geometry.py) uses only
the Python standard library. It reconstructs (4), all four tangents
in (6), their rank and the quadratic-form preservation identities with
rational polynomial arithmetic. It also enumerates every possible
one-sided linear/cubic monomial product used in the weight-separation
argument and checks that each has zero original Wick value. A symmetric
negative-weight check and explicit zero-piece controls are included.
These finite controls validate the displayed algebra. They do not
verify the imported all-pure classification, the dimension theorems,
Nullstellensatz, the existence of the unspecified witness, or a
universal upper prefix. Those conclusions use the written proof.

Run from the repository root:

```text
.venv/bin/python research/moments/vanishing/verify_gaussian_cubic_prefix_geometry.py
```

The accompanying
[receipt](gaussian_cubic_prefix_geometry_verification.json) records the
checked algebra, source hashes and exact assurance limits. The theorem
depends on the retained written classification; a revision to that input
must trigger a review of this deduction. Publication, external review,
formal verification and worldwide priority are separate questions.
