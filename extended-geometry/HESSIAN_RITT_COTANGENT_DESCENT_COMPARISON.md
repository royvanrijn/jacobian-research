# Cotangent descent for the Hessian--Ritt diagram

> **Status.** Cotangent descent for the full simplicial bar presentation of
> the actual derived Hessian intersection is proved in Section 2.  The
> face-poset bar compresses canonically to cellular chains; exact
> subdivision maps are certified for the filled braid and relative
> half-braid.  Thus the remaining comparison is precisely the presentation
> of the actual bar coefficients as one coherent, homotopy-cofinal Ritt
> face diagram.  The required local coefficients are identified on the
> degree-thirty conormal fibers and on the certified completed
> degree-forty-two ideal flag.  All six degree-forty-two first conormal
> flags and the fourth jets of the two remaining labelled sectors are now
> computed.  Thin/boundary equality and quadratic-overlap vanishing are
> also proved after completion for both new sectors.  Thus all three
> degree-forty-two half-braids satisfy the first-Postnikov comparison;
> coherent completed extension transport remains open.  Full algebraic
> \(H^2\)-descent requires genuine three-cells and remains open beyond the
> displayed Coxeter regression.

This note separates two statements that were previously bundled into the
phrase “cellular cotangent descent”:

1. the actual derived intersection has a canonical bar presentation whose
   cotangent complex satisfies descent; and
2. the much smaller coefficient-decorated Ritt cell complex is an effective
   cellular compression of that bar presentation.

The first statement is formal.  The second contains the Ritt geometry.

## 1. The actual derived intersection

Fix a Hessian coefficient ambient algebra \(P\).  For every requested
normalized decomposition word \(v\), let \(A_v\) be the coordinate algebra
of its factor chart, regarded as a \(P\)-algebra through normalized
composition followed by Hessian projection.  For a finite collection
\(D\) of words, define

\[
 A_D^{\mathrm{der}}
 =
 \mathop{\bigotimes\nolimits^{\mathbf L}}_{v\in D,P}A_v.
                                                                    \tag{1.1}
\]

Equivalently, \(A_D^{\mathrm{der}}\) is the homotopy colimit in
commutative \(P\)-algebras of the star diagram

\[
 P\longrightarrow A_v,\qquad v\in D.                        \tag{1.2}
\]

Its derived affine spectrum is the derived fiber product of the factor
charts over the Hessian coefficient space.  Its classical truncation is the
scheme-theoretic Hessian intersection because the normalized triangular
reconstruction makes every factor chart a closed incidence subscheme.

Thus (1.1), rather than a bare graph of factor words, is the canonical
meaning of the actual derived Hessian intersection.

## 2. Cotangent complexes commute with this colimit

Let \(F:\mathcal I\to\operatorname{CAlg}_P\) be any small diagram of
connective commutative \(P\)-algebras and put

\[
 A=\operatorname*{colim}_{i\in\mathcal I}F(i).               \tag{2.1}
\]

There is a canonical equivalence in \(\operatorname{Mod}_A\)

\[
 \boxed{
 L_{A/P}
 \simeq
 \operatorname*{colim}_{i\in\mathcal I}
 \left(A\otimes_{F(i)}^{\mathbf L}L_{F(i)/P}\right).}         \tag{2.2}
\]

To prove (2.2), let \(M\) be an \(A\)-module.  The universal property of
the cotangent complex and the colimit give

\[
\begin{aligned}
 \operatorname{Map}_A(L_{A/P},M)
 &\simeq \operatorname{Der}_P(A,M)\\
 &\simeq
 \lim_{i\in\mathcal I}\operatorname{Der}_P(F(i),M)\\
 &\simeq
 \lim_{i\in\mathcal I}
 \operatorname{Map}_A
 \left(A\otimes_{F(i)}^{\mathbf L}L_{F(i)/P},M\right).
\end{aligned}                                                \tag{2.3}
\]

The last line corepresents the colimit on the right of (2.2), so Yoneda
proves the equivalence.

Applying (2.2) to (1.2) proves:

> **Bar cotangent-descent theorem.** The actual derived Hessian
> intersection always satisfies cotangent descent over its full simplicial
> bar construction.

No Ritt classification or overlap calculation is needed for this theorem.

## 3. From the bar construction to cellular chains

The homotopy colimit in (2.2) is computed by the simplicial replacement, or
two-sided bar construction, of \(F\).  Its normalized chains contain one
summand for every nondegenerate composable string in \(\mathcal I\).

Let \(K\) be a finite regular CW complex and suppose that the bar diagram
descends to a constructible cellular cosheaf \(\mathcal L\) on its face
category.  Then its homotopy colimit is computed by cellular chains

\[
 C_p^{\mathrm{cell}}(K;\mathcal L)
 =
 \bigoplus_{\sigma\in K^{(p)}}\mathcal L_\sigma,             \tag{3.1}
\]

with the signed incidence and corestriction differential.  If all
\(\mathcal L_\sigma\) are perfect over the completed base algebra, finite
duality converts (3.1) into the cellular cochain totalization used for
tangent and deformation complexes:

\[
 \mathbf R\!\operatorname{Hom}
 \left(
 \operatorname*{hocolim}\mathcal L,A
 \right)
 \simeq
 \operatorname*{holim}
 \mathbf R\!\operatorname{Hom}(\mathcal L_\sigma,A).         \tag{3.2}
\]

Therefore the cochain orientation of the executable prototype is not in
conflict with the algebra colimit in (2.2): it appears after perfect
duality.

There are now two distinct comparisons:

\[
\begin{array}{ccc}
 \operatorname{Bar}_D(L)^\vee
 &\xrightarrow{\ \mathrm{presentation}\ }&
 \operatorname{Bar}_{\operatorname{Fc}(K_D)}(\mathcal L^\vee)\\
 &&\mathrel{\Big\downarrow}\scriptstyle{\mathrm{cellular\ subdivision}}\\
 &&C^\bullet_{\rm cell}(K_D;\mathcal L^\vee).
\end{array}                                                    \tag{3.3}
\]

Here \(\operatorname{Fc}(K_D)\) is the face category.  The right vertical
map is formal once \(\mathcal L\) is a genuine constructible cellular
coefficient functor.  The top horizontal map is the Ritt-specific
presentation assertion.

More precisely, suppose there is a homotopy-coherent functor

\[
 u:\mathcal B_D\longrightarrow\operatorname{Fc}(K_D)         \tag{3.4}
\]

from the indexing category of the actual normalized bar construction and a
perfect face-category coefficient system \(\mathcal L\) such that the bar
cotangent diagram is equivalent to \(u^*\mathcal L\).  If the relevant
comma category of \(u\) is contractible over every face, then homotopy
cofinality identifies the two bar totalizations in (3.3).  Cellular
subdivision then identifies the face-category bar with (3.1), naturally in
\(\mathcal L\).  This gives the following useful criterion.

> **Cellular-descent criterion.** The actual Hessian--Ritt cotangent
> diagram satisfies the cellular comparison provided that:
>
> 1. its completed factor, move, and relation maps form a coherent perfect
>    functor on the Ritt face category; and
> 2. the bar-to-face functor (3.4) is homotopy cofinal.
>
> No separate chain-level exactness calculation is then required.

The point is not merely terminological.  Condition 1 contains the
coefficient transport and condition 2 is combinatorial.  Neither should be
hidden inside a single degree-specific Gröbner calculation.

For reference, the use of cellular coefficient functors and their cellular
cochain complexes is standard in
[Curry's thesis](https://arxiv.org/abs/1303.3255); the compatibility between
regular CW face posets and cellular chain models also appears in
[Clark--Tchernev](https://arxiv.org/abs/1310.2315).

### 3.1 Exact subdivision regression

The checker constructs the normalized face-poset bar explicitly in the two
local shapes used by the prototype.

* For a three-edge half-braid relative to its endpoints, the cellular chain
  dimensions are \((2,3)\), the normalized face-bar dimensions are
  \((5,6)\), and the canonical subdivision map has acyclic mapping cone.
  Both complexes have homology \((0,1)\).
* For the filled braid disk, the corresponding dimensions are
  \((6,6,1)\) and \((13,24,12)\).  Again the mapping cone is acyclic, and
  both homologies are \((1,0,0)\).

The subdivision maps are written with exact rational matrices.  Tensor
tests in coefficient dimensions \(2,4,6\) remain quasi-isomorphisms.  Since
tensoring over \(\mathbb Q\) is exact, this proves the same statement for
every finite perfect constant block, including the underlying vector spaces
of the certified non-split degree-forty-two conormal modules.  Their module
extension is retained because all cellular maps are equivariant identity
maps on the coefficient block; the calculation does not split the
coefficient module.

## 4. The skeletal range

For a connective cellular coefficient system, the cochain complex begins

\[
 C^0\overset{d_0}{\longrightarrow}
 C^1\overset{d_1}{\longrightarrow}
 C^2\overset{d_2}{\longrightarrow}
 C^3\longrightarrow\cdots.                                  \tag{4.1}
\]

The complete two-skeleton determines

\[
\begin{aligned}
 H^0&=\ker d_0,\\
 H^1&=\ker d_1/\operatorname{im}d_0.
\end{aligned}                                                \tag{4.2}
\]

No cell of dimension at least three can change (4.2).  By contrast,

\[
 H^2=\ker d_2/\operatorname{im}d_1                           \tag{4.3}
\]

cannot be read from a two-skeleton unless the next differential is known.

This distinction is visible in the Coxeter regressions.

* The filled three-factor braid has cellular dimensions
  \((6,6,1,0)\), ranks \((5,1,0)\), and cohomology
  \((1,0,0,0)\).  Its two-cell already fills the hexagon.
* The four-factor Coxeter two-skeleton has dimensions
  \((24,36,14)\), ranks \((23,13)\), and apparent cohomology
  \((1,0,1)\).
* The genuine permutohedron three-cell has oriented boundary

\[
 (1,-1,-1,1,-1,-1,1,1,1,1,-1,-1,1,1),                     \tag{4.4}
\]

raises the last rank to one, and changes the cohomology to
\((1,0,0,0)\).

Tensoring by a two-dimensional perfect constant coefficient duplicates the
top relation and gives cohomology \((2,0,0,0)\).  Thus the top-cell
phenomenon is functorial in perfect coefficients.

## 5. The coefficient-effectivity condition

Let \(K_D\) be the Ritt cell complex generated by factor words, elementary
moves, commuting squares, braid hexagons, and the required higher Coxeter
cells.  Let \(\operatorname{Bar}_D\) denote the full bar diagram computing
(1.1).  For every cell \(\sigma\), let \(L_\sigma\) be its completed local
cotangent or relative-defect complex.  The remaining geometric condition is:

> **Coefficient effectivity.** The local decoration \(L_\sigma\) is the
> homotopy colimit of the bar objects over the latching boundary of
> \(\sigma\), and the maps are compatible under face incidence.

Equivalently, the comparison from cellular chains to normalized bar chains
must be a quasi-isomorphism in the required Postnikov range.  By the
criterion in Section 3, this reduces further to coherent coefficient
transport plus homotopy cofinality of (3.4); the face-bar-to-cellular step
itself is no longer an additional conjecture.  Under coefficient
effectivity and local perfectness, (2.2)--(3.2) prove

\[
 L_{A_D^{\mathrm{der}}/P}^{\vee}
 \simeq
 \operatorname{Tot}
 C^\bullet(K_D;L_\sigma^\vee).                              \tag{5.1}
\]

This is the desired Hessian--Ritt cellular-descent comparison.  The proof is
now formal once coefficient effectivity is supplied.

## 6. What is currently effective

On every tame factor-order block modeled by permutations, the unlabelled
combinatorial two-skeleton is complete: adjacent transpositions are the Ritt
moves, while commuting and braid relations present all relations among
factor-order words.  The remaining issue there is the coefficient
decoration, not Coxeter connectivity.

### 6.1 Degree thirty

The exact degree-thirty checker transports both polynomial and Hessian path
ideals across all six factor-order charts.  For every labelled half-braid,
one path is the reduced boundary and the other has its certified conormal
fiber at the monomial point.  The relative three-edge cellular block

\[
 D^2\longrightarrow D^3
\]

is a free cellular resolution with \(H^1=D\) and no \(H^0\) or \(H^2\).
Consequently the comparison is effective on the first conormal fiber of
every certified degree-thirty braid sector.  This proves the fiberwise
degree-thirty comparison on \(H^0\) and \(H^1\).

It does not yet identify the completed restriction maps as one cellular
coefficient system, nor reconstruct the nonlinear Artin algebras from a
Chevalley--Eilenberg differential.  Those are higher filtered statements.

### 6.2 Degree forty-two

On the normalized \(2\mathbin\circ7\mathbin\circ3\) chart, HRCELL4--HRCELL6
identify the actual conormal tower, prove its first homology sequence short
exact, retain its non-split extension, and distinguish finite-base-change
Tor.  Thus coefficient effectivity holds locally through \(H_1\) for this
flag chart.

The completed flag and comparison maps have not yet been transported to the
other five factor-order charts.  Therefore a completed global
degree-forty-two braid comparison is not claimed.  The first-order and
fourth-jet advances below sharply constrain that remaining transport.

There is now, however, an all-six-chart necessary first-order calculation.
Put every normalized factor-chart tangent image into the same Hessian
coefficient space

\[
 T\mathcal K_{42}
 =\langle W^2,\ldots,W^{41}\rangle.                          \tag{6.1}
\]

Formula (1.2) computes each image without elimination.  At
\((t,z)=(1,0)\), all six vertex images have rank nine.  The literal
intersections for the six adjacent Ritt moves have ranks

\[
\begin{array}{c|cccccc}
\text{edge}
 &237\!-\!327&237\!-\!273&273\!-\!723
 &327\!-\!372&372\!-\!732&723\!-\!732\\ \hline
\dim(T_v\cap T_w)&8&5&6&6&5&8.
\end{array}                                                   \tag{6.2}
\]

Thus the actual move coefficients are not a constant rank-two Dickson
system.  The intersection of all six vertex tangent images has rank three.
It contains the rank-two tangent plane of the reduced Dickson component,
and its quotient is generated in Hessian coefficient space by

\[
 \pi_H\bigl((W+1)^{36}-1\bigr).                              \tag{6.3}
\]

In particular, the extra common direction is Hessian-visible; forgetting
the linear coefficient does not put it in the Dickson tangent plane.  This
is the global first-order shadow expected from the common spectator layer.
Because (6.2) is computed as literal subspace intersection in one ambient
space, all incidence inclusions commute.  It supplies the vertex/move
incidence lattice needed by a degree-forty-two face diagram, but not yet
the covariant restriction maps and completed two-cell latching comparison.
Those require transporting the non-split conormal flag, not merely its
tangent support.

Intersecting the four vertex images on each half-braid gives more.  The
three opposite-pair representatives and their omitted cuts are

\[
\begin{array}{c|ccc}
\text{representative}&273&237&327\\ \hline
\text{thick composite omission}&6&14&21\\
\text{thin prime omission}&7&3&2.
\end{array}                                                   \tag{6.4}
\]

In every column the thick-path, thin-path, full-boundary, and reduced
Dickson tangent dimensions are

\[
 (4,3,3,2),                                                   \tag{6.5}
\]

so their conormal ranks inside the nine-dimensional factor chart are

\[
 \boxed{(5,6,6,7).}                                          \tag{6.6}
\]

Reversing the paths supplies the opposite three charts.  Hence the entire
degree-forty-two braid now has the same first conormal flag
\(5<6<7\) previously certified only on \(273\).  This proves all-chart
first-order transport of the sector and spectator lines.  It does not prove
that the new completed sector modules have the same annihilators or
non-split extension class: those are nonlinear invariants, and the
exploratory rotated-chart checker is computing them separately.

The first nonlinear comparison is also complete through the fourth
maximal-adic jet for the two new representatives.  In both cases the thin
path equals the boundary modulo \(\mathfrak m^4\), the ideal flag has the
expected inclusions, and the conormal ranks remain \((5,6,6,7)\).  The
quotient-ring lengths are

\[
\begin{array}{c|ccc|ccc}
 &\multicolumn{3}{c|}{S/(I_{\rm thick}+\mathfrak m^q)}
 &\multicolumn{3}{c}{S/(I_\partial+\mathfrak m^q)}\\
\text{omission}&q=2&q=3&q=4&q=2&q=3&q=4\\ \hline
6&5&14&29&4&9&16\\
14&5&13&25&4&9&16\\
21&5&13&26&4&9&16.
\end{array}                                                   \tag{6.7}
\]

Here the cut-\(6\) row is the existing certified jet profile, while the
other rows are the new exact computations.  Since the Dickson graph lengths
are \((3,6,10)\), the sector and spectator layer dimensions are

\[
\begin{array}{c|c|c}
\text{composite omission}
&\dim(I_\partial/I_{\rm thick})_{q=2,3,4}
&\dim(K/I_\partial)_{q=2,3,4}\\ \hline
6&(1,5,13)&(1,3,6)\\
14&(1,4,9)&(1,3,6)\\
21&(1,4,10)&(1,3,6).
\end{array}                                                   \tag{6.8}
\]

Thus the spectator layer transports uniformly through order four, whereas
the sector layer is genuinely labelled by the omitted composite cut.
This is positive evidence for cellular descent but negative evidence for a
constant defect coefficient.

Modulo \(\mathfrak m^4\), \(z^3\) kills each new sector quotient and \(z\)
kills the spectator quotient.  The exponent three is only a truncated
statement: it cannot be compared directly with the completed \(z^8\)
annihilator on the cut-\(6\) sector, because every sufficiently high
\(z\)-power vanishes in a fourth jet.  Upgrading (6.7)--(6.8) to completed
modules requires a new Artin--Rees/Nakayama cutoff for the cut-\(14\) and
cut-\(21\) flags.

That upgrade is now complete for both new cuts.  Begin with cut \(14\).
Write \(I\) for the thick
cut-\(14\) path ideal, \(L\) for the thin cut-\(3\) path ideal, \(J\) for
the full-boundary ideal, \(K\) for the Dickson graph ideal, and
\(\mathfrak m=(K,\tau,\zeta)\).  Exact reduction gives

\[
 J\subseteq L+\mathfrak mJ,\qquad L\subseteq J.              \tag{6.9}
\]

Nakayama therefore proves \(\widehat J=\widehat L\); this is completed
thin-path/boundary equality, not fourth-jet stabilization.  For the sector
source

\[
 S_{14}=J/(I+KJ)
\]

modulo \((\tau,\zeta)^2\), the explicit Artin--Rees cutoff
\(J\cap\mathfrak m^4\subset I+KJ+(\tau,\zeta)^2J\) has zero remainder.
The corresponding denominator and numerator colengths are \(19\) and
\(16\), so this base-square quotient has dimension \(3\).

For the quadratic overlap

\[
 \mathcal O_{14}
 =\frac{J\cap(I+K^2)}{I+KJ},                                 \tag{6.10}
\]

the numerator has \(49\) standard generators.  The exact cutoff

\[
 \bigl(J\cap(I+K^2)\bigr)\cap\mathfrak m^5
 \subset I+KJ+(\tau,\zeta)^2\bigl(J\cap(I+K^2)\bigr)          \tag{6.11}
\]

has zero remainder, and the two finite colengths are both \(34\).
Consequently
\(\mathcal O_{14}/(\tau,\zeta)^2\mathcal O_{14}=0\), and
Nakayama proves

\[
 \boxed{\widehat{\mathcal O}_{14}=0.}                         \tag{6.12}
\]

The cut-\(21\)/thin-cut-\(2\) chart satisfies the same exact containments.
The cut-specific standard-basis data are

\[
\begin{array}{c|cc|ccc}
\text{omission}
&\#(J\cap\mathfrak m^4)&
 \dim S_c/(\tau,\zeta)^2S_c&
\#\operatorname{gens}(J\cap(I+K^2))&
\#(\mathcal O\text{-cutoff intersection})&
\text{overlap colengths}\\ \hline
14&490&3&49&1278&(34,34)\\
21&495&3&61&1278&(34,34).
\end{array}                                                   \tag{6.13}
\]

Every displayed cutoff remainder is zero: order four for the sector source
and order five for the overlap.  Hence Nakayama also gives

\[
 \widehat J_{21}=\widehat L_{21},\qquad
 \boxed{\widehat{\mathcal O}_{21}=0}.                         \tag{6.14}
\]

Thus the completed first-Postnikov conormal sequence is short exact on all
three degree-forty-two half-braids (cuts \(6,14,21\)).  For each new sector
the source base-square quotient has dimension three, whereas the \(q=2\)
quotient-length difference in (6.7) is one; the two-dimensional loss after
finite base change is Tor, not a completed quadratic overlap.  The
extension classes and braid restriction coherence remain separate
obligations.

The first tensor-presentation transport test changes the expected picture.
For cut \(14\), tensor the presented conormal projection with
\(B/(\tau,\zeta)^4\).  Its exact dimensions and cocycle ranks are

\[
 0\longrightarrow\mathbf Q^9\longrightarrow\mathbf Q^{13}
 \longrightarrow\mathbf Q^4\longrightarrow0,\qquad
 \operatorname{rank}\delta=\operatorname{rank}[\delta\mid c]=32. \tag{6.15}
\]

An explicit simultaneous section for the \(\tau,\zeta\) actions exists;
the seven normal actions vanish.  Hence this module extension splits
through base order four, and therefore through every lower order.  The
cut-\(21\) calculation independently splits at orders two and three, with
dimensions \(3\to5\to2\) and \(6\to9\to3\).  These facts do **not** prove a
completed split: one must still lift compatible sections to every order.
They do prove that the cut-\(6\) base-square non-splitting obstruction does
not transport uniformly to the rotated sectors.

The inverse-limit lifting problem is now closed for both rotated sectors.
Rather than choosing new sections independently, take the exact order-seven
tensor presentation and quotient that single module extension by
\((\tau,\zeta)^q\), for \(q=5,6,7\).  The resulting compatible towers have

\[
\begin{array}{c|ccc|cc}
q&\dim K_q&\dim T_q&\dim S_q&
\operatorname{rank}\delta_q&
\operatorname{rank}[\delta_q\mid c_q]\\ \hline
5&12&17&5&55&55\\
6&15&21&6&84&84\\
7&18&25&7&119&119.
\end{array}                                                   \tag{6.16}
\]

This table is identical for omissions \(14\) and \(21\).  The chosen
order-seven section reduces exactly to the displayed order-six and
order-five sections.  The vector spaces of section differences have
dimensions \(5,6,7\).  Their restriction maps have ranks \(3\) and \(4\),
so each has a two-dimensional cokernel.  Thus independent finite splitting
would not by itself justify an inverse limit: a coherence obstruction really
has room to occur.

The completed presentations remove that ambiguity.  Put
\(u=1+\tau\) and \(B=\mathbf Q[[\tau,\zeta]]\).  In both sectors the
spectator module is

\[
 S_c\cong\mathbf Q[[\tau]]e_4,\qquad \zeta e_4=0.             \tag{6.17}
\]

For cut \(14\), the total presentation contains

\[
 -e_5+3u e_6=0,\qquad
 -\zeta e_4+u\zeta e_5-2\zeta^2e_6=0.
\]

Consequently

\[
 s_{14}(e_4)=e_4+(-3u^2+2\zeta)e_6,\qquad
 \zeta s_{14}(e_4)=0.                                       \tag{6.18}
\]

For cut \(21\), the relevant relations are

\[
\begin{aligned}
 -2e_5+3u e_6-8\zeta e_7&=0,\\
 -e_6+4u e_7&=0,\\
 -3\zeta e_4+2u\zeta e_5-4\zeta^2e_6&=0.
\end{aligned}
\]

They give

\[
 s_{21}(e_4)=e_4+(-4u^3+8u\zeta)e_7,\qquad
 \zeta s_{21}(e_4)=0.                                       \tag{6.19}
\]

Both displayed elements project to \(e_4\).  Exact reduction of every
spectator relation under these maps has zero remainder in the corresponding
total presentation, and the two retraction residuals have zero remainder in
the spectator presentations.  Hence

\[
\boxed{0\longrightarrow K_c\longrightarrow T_c
       \longrightarrow S_c\longrightarrow0
       \text{ splits over }B,\quad c=14,21.}                 \tag{6.20}
\]

Equivalently, the completed extension obstruction

\[
 [\operatorname{id}_{S_c}]
 \in\operatorname{coker}\!\left(
 \operatorname{Hom}_B(S_c,T_c)\longrightarrow
 \operatorname{End}_B(S_c)\right)                           \tag{6.21}
\]

is zero.  The polynomial sections (6.18)--(6.19) give a point of the full
inverse limit of section torsors, so the corresponding
\(\varprojlim^1\)-obstruction class is zero despite the nonsurjective finite
restriction maps.

This also explains the sector asymmetry in (6.8).  Before
presentation-first tensoring, the cut-\(21\) quotient-ring fourth jet has
one extra sector dimension.  After tensoring the presented first-conormal
modules, both sectors have

\[
 \dim K_q=3q-3,\qquad \dim T_q=4q-3,\qquad \dim S_q=q
 \quad(2\leq q\leq7).                                       \tag{6.22}
\]

Thus the extra cut-\(21\) quotient-ring dimension lies in the non-flat
base-change/Tor discrepancy; it is not a different completed extension
rank or a split/non-split asymmetry.  A genuine labelled asymmetry remains:
the completed corrections are respectively
\((-3u^2+2\zeta)e_6\) and \((-4u^3+8u\zeta)e_7\).  These different
polynomials record the omitted composite cut.  Comparing them under the
actual restriction maps between factor charts is still the separate braid
coherence problem.

### 6.1. The embedded-support restriction test

The characteristic-zero sixth-jet class \(c_6\) from
[the support-saturation compiler](SUPPORT_SATURATION_COMPILER.md) is a new
input to this cellular programme.  It is nonzero on the specialized
cut-\(6\) jet and is killed by the full boundary ideal.  Thus generic
vanishing away from the boundary does not justify a global saturation
claim, even when every previously inspected generic stratum vanishes.

At present there is no proved comparison map identifying this
local-cohomology class with the cut-\(6\) cotangent non-splitting class of
HRCELL2--HRCELL4.  Constructing that map is the first task.  If an image is
defined, the completed splittings (6.18)--(6.20) make the rotated
cut-\(14\) and cut-\(21\) restrictions the next finite calculations.  One
must then transport all three representatives through the labelled edge
maps and evaluate the closed braid cochain.

There are two mathematically distinct outcomes: the rotated restrictions
and braid homotopy may kill the class, or they may exhibit a coherent
embedded class that persists beyond the cut-\(6\) chart.  Local splitting
on the rotated sectors does not decide between them.  Likewise, the
surjective order-seven-to-order-six saturation transition presently known
modulo \(32003\) is evidence for finite-jet lifting, not a
characteristic-zero persistence theorem.  This restriction/coherence test,
rather than another global colon, is the immediate `HRCELL` interface with
support saturation.

## 7. Result and remaining proof obligation

The comparison problem now has the following exact status.

1. The actual derived Hessian intersection satisfies full bar cotangent
   descent by (2.2).
2. A genuine perfect coefficient functor on the Ritt face category has a
   canonical face-bar-to-cellular quasi-isomorphism.
3. The Ritt Coxeter cells give the complete unlabelled two-skeleton.
4. A complete coefficient-decorated two-skeleton determines \(H^0,H^1\),
   but not \(H^2\) without the next cellular differential.
5. Coefficient effectivity through \(H_1\) is proved on the degree-thirty
   conormal fibers and on the certified completed local degree-forty-two
   flag.  In degree forty-two the first conormal flag now transports over
   all six charts, and thin/boundary equality plus first-Postnikov overlap
   vanishing are proved after completion on all three half-braids.
6. The two rotated first-conormal extensions split after completion, with
   explicit polynomial sections compatible with every base truncation.
7. The all-degree finite Ritt presentation remains equivalent to proving
   coherent coefficient transport and homotopy cofinality for the universal
   power move, Dickson move, commuting square, and labelled braid, followed
   by composition base change.

The next concrete calculation is first to construct the comparison map for
\(c_6\) described in Section 6.1.  Its rotated restrictions then use the
same factor-chart transition maps needed to transport (6.18)--(6.19) and
verify the braid and commuting-cell restrictions.  Completed local
splitting is no longer the obstruction; identification of the embedded
class and labelled restriction coherence are.  Closing them would promote
the local \(H_1\) comparisons to the full degree-forty-two cellular
coefficient diagram.

## Reproduction

Run

```bash
.venv/bin/python scripts/verify_hessian_ritt_cotangent_descent.py
.venv/bin/python scripts/verify_degree42_ritt_cut14_postnikov_overlap.py
.venv/bin/python scripts/verify_degree42_ritt_cut21_postnikov_overlap.py
.venv/bin/python scripts/verify_degree42_ritt_cut14_tensor_split_q4.py
.venv/bin/python scripts/verify_degree42_ritt_inverse_limit_sections.py
.venv/bin/python scripts/verify_degree42_ritt_completed_splits.py
```

The commands write
`artifacts/generated-results/hessian_ritt_cotangent_descent.json`,
`artifacts/generated-results/degree42_ritt_cut14_postnikov_overlap.json`,
and
`artifacts/generated-results/degree42_ritt_cut21_postnikov_overlap.json`,
and
`artifacts/generated-results/degree42_ritt_cut14_tensor_split_q4.json`,
`artifacts/generated-results/degree42_ritt_inverse_limit_sections_q5_q7.json`,
and
`artifacts/generated-results/degree42_ritt_completed_splits.json`.
