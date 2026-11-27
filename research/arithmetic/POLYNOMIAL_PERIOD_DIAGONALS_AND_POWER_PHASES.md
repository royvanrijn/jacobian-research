# Period diagonals, power phases, and the degree-eight cubic split

Grouping several ramification points into one branch value imposes
equations inside a fixed-residue differential space. It is not the
operation of deleting a critical-value collision divisor used for
the [single-cycle-tail theorem](POLYNOMIAL_SINGLE_CYCLE_TAIL_COMPONENTS.md).
Those equations can disconnect a connected ambient space. Their
components can retain genuine finite-group lifting obstructions even
after all cyclic power-phase data agree.

We give an exact period description for every minimal polynomial
passport, a finite rational-cell reduction with a distinguished branch,
and a canonical power-phase invariant. The degree-eight example
admits a complete algebraic calculation:
\[
 f_\pm(z)=\left[z^2\left(z^2+z+
                   \frac{9\pm3\sqrt3}{16}\right)\right]^2.
                                                               \tag{1}
\]
These are representatives, up to source affine changes and target
scaling, for its two same-group, same-class cubic-lifting components.
They have the same power-phase pattern. The third representative is
\[
                         f_0(z)=(z^4+z^2)^2.                 \tag{2}
\]
The general period reduction does not by itself classify the
components of all its sliced cell graphs. That remains part of the
full arbitrary-passport problem.

## Primary inputs and their boundaries

Chen--Mullane--Salter--Tahar,
[The residual map as a stratified fibration over the resonance arrangement](sources/chen-mullane-salter-tahar-residual-stratified-2026.pdf),
arXiv:2611.07499v1, November 2026, supply the relevant primary period
model. We use the actual statements and construction:

- Proposition 8.2, printed p45, identifies integral-residue
  differentials with logarithmic derivatives of rational maps.
  Branch values are allowed to coincide; a prescribed grouping
  of them is an additional condition.
- Theorem 5.3, p24, and Section 6 identify the fixed-residue fiber
  with its glued period domains, for vertically generic residues.
- Definitions 5.10 and 5.16--5.18, pp26--30, give distinguished
  saddles, integer cycle equations and their positive period cones.
- Section 6 gives the face identifications; Lemma 7.8 and its
  proof, p44, supply finiteness of graph types for fixed orders.

Their residue convention is the usual coefficient residue,
unlike the period-residue convention in
[Chen--Gendron--Prado--Tahar](sources/chen-gendron-prado-tahar-isoresidual-curves-2026.pdf).
All normalizations below specify this difference. The new input is a
preprint; its theorems are imported, not independently verified here.
The residual-fibration homeomorphism alone does not identify
arbitrary prescribed critical-value diagonals.

## Theorem G1. Every polynomial passport is a colored period diagonal

Distinguish a finite branch of a nonidentity passport and translate
its value to zero. Write its cycle partition as
\(\lambda=(\ell_1,\ldots,\ell_k)\), with sum \(N\).
List the nontrivial cycles at every other branch as lengths
\(a_j+1\), \(1\le j\le n\), and color each \(j\) by its branch.
There are \(c=r-1\) nonempty colors, where \(r\ge2\) is the number
of finite branches. Minimality says
\[
                         \sum_{j=1}^n a_j=k-1.               \tag{3}
\]
On the fixed-residue fiber use
\[
 \eta=\frac1{2\pi i}\frac{df}{f},\qquad
 \operatorname{res}(\eta)=
      \frac1{2\pi i}(\ell_1,\ldots,\ell_k,-N).                 \tag{4}
\]
The coefficient residues in (4) are vertically generic and
nonresonant. If \(W_{ij}=\int_{z_i}^{z_j}\eta\), then the exact
prescribed passport locus is
\[
 \begin{split}
 W_{ij}&\in\mathbb Z&&\text{for zeros of the same color},\\
 W_{ij}&\notin\mathbb Z&&\text{for zeros of different colors}.
 \end{split}                                                  \tag{5}
\]
The conditions are independent of the integration paths.
With zeros and roots labeled, this is a smooth complex locus of
dimension \(c-1=r-2\), whenever nonempty.

**Proof.** Every absolute period of \(\eta\) is an integer
combination of the \(\ell_j\), so changing paths changes
\(W_{ij}\) by an integer. Proposition 8.2, or the direct product
formula in F1, recovers the polynomial up to a nonzero scalar.
Its critical-value ratio is
\[
                   \frac{f(z_j)}{f(z_i)}
                              =\exp(2\pi i W_{ij}).           \tag{6}
\]
This proves (5). These are exactly the nonzero branch groups;
none meets the root fiber. Positive finite residues and their
negative sum give nonresonance, as in F1.
Moreover \(2\pi i\operatorname{res}(\eta)\) has real, nonzero
proper partial sums. This is precisely vertical genericity
in the new primary period model.

Relative periods from one zero to the others are local
coordinates of dimension \(n-1\). Choose one representative
of each color. The same-color conditions are \(n-c\)
independent complex affine equations with integer right
sides. Different-color coincidences are excluded. This
gives the dimension and smoothness assertion.
\(\square\)

Root labels of equal multiplicity and zero labels of equal
order within a color are then forgotten. Equal-type branch
colors may also be forgotten. If the distinguished type
occurs more than once, forgetting which branch was
distinguished includes the operation
\(f\mapsto f-f(z_j)\), followed by a new root fiber.
This last operation is not an affine change of period
coordinates on one fixed chart.
Thus the full braid convention is the finite label-forgetting
quotient of the distinguished-branch model; an individually
marked branch convention must retain its markings.

## Theorem G2. A finite rational-cell reduction for the distinguished branch

For fixed labeled data in G1, the locus (5) is represented
by finitely many rational polyhedral cells with explicit
face identifications. Its connected components are those
of the finite incidence graph of the retained cells.
The construction uses only period domains of smooth
differentials with the specified distinct zeros; points
in the metric completion with collided zeros or arboreal
sources are excluded.

More explicitly, in every critical-graph chart choose a
spanning tree of the distinguished saddle graph. Its path
periods between zeros satisfy
\[
                    |\operatorname{Re}W_{ij}|
                                 \le (n-1)N.                \tag{7}
\]
Hence only the integers in this interval need occur in (5).

**Proof of the bound.** In the flat coordinate
\((2\pi i)^{-1}\log f\), a horizontal level corresponds to
\(|f|=R\). For a regular \(R\) this is a union of smooth
closed curves, each mapping to the target circle by a
positive covering degree. The sum of those degrees is
\(N\), so their total horizontal flat length is \(N\).
The sublevel set has no holes: a bounded hole on which
\(|f|>R\) would violate the maximum principle.

Each strip of the vertically generic critical graph
contributes its width to such a horizontal section.
Thus there are no half-space faces and every strip
width is at most \(N\). A distinguished transverse
saddle crosses one strip and has real period equal
to its width; a vertical saddle has real period zero.
A simple tree path has at most \(n-1\) edges, proving (7).

**Proof of the reduction.** The primary Theorem 5.3 gives
the period atlas and exact gluing. Its cycle equations
have coefficients \(0,\pm1\) and right sides given by
absolute residue periods. For (4) these right sides
are integers. Its positivity inequalities and face
identifications are therefore rational. There are
finitely many critical graphs for the fixed orders.

In each chart subdivide by
\(\operatorname{Im}W_{ij}=0\) and
\(\operatorname{Re}W_{ij}=m\) for all pairs and all
integers \(|m|\le(n-1)N\). Also refine by the finite
polyhedral gluing subdomains in Section 6, so that
identifications are along whole refined cells.
Coordinate changes express paths as integer combinations
of saddle paths plus integer absolute periods.
Consequently all subdivisions and identifications
are exact rational operations.

Retain precisely those relatively open cells satisfying
(5). Keeping the whole subdivision is necessary:
deleting an unwanted equality can cut a real face,
even when it has complex codimension one in the
ambient differential locus. Every retained cell is
convex. Join a retained cell to a retained face
where the source period model identifies them.
A segment entering that face supplies a path.
Conversely, a finite polyhedral decomposition has
the same path components as this incidence graph,
by subdivision and piecewise linear paths.
No connection through an excluded boundary face
is introduced. This proves the statement.
\(\square\)

G2 is a finite reduction, not an implementation of all
critical-graph enumeration and gluing. The ordinary
equal-label quotients act by finite combinatorial
renaming. Exchanging the distinguished branch when its
type repeats requires the additional action described
after G1. Computing these graphs, their full-branch
quotients, and all monodromy and lifting labels for
arbitrary passports is still an open obligation.

## Theorem G3. Canonical cyclic power extraction and phase data

Let \(D=\gcd(\ell_1,\ldots,\ell_k)\). For every polynomial
in G1 there is a polynomial \(h\) of degree \(N/D\) such
that, after target scaling,
\[
                             f=h^D.                          \tag{8}
\]
The root divisor determines \(h\) up to a nonzero scalar.
If \(D>1\), the generated group preserves the standard
residue partition with \(D\) blocks. Only the distinguished
factor acts on that block quotient, by a full \(D\)-cycle;
all other factors fix the blocks.

For two labeled critical points in one color,
\[
        \chi_{ij}=\frac{h(z_j)}{h(z_i)}
              =\exp\!\left(\frac{2\pi i}{D}W_{ij}\right)
                                   \in\mu_D                 \tag{9}
\]
is well-defined, path-independent and locally constant.
Forgetting point labels gives, for each color, the
cycle-length data in its \(D\) phase bins modulo cyclic
rotation. These are the cycle types of the corresponding
block return permutations.

**Proof.** Divide the root multiplicities by \(D\) in the
product formula for \(f\); this gives (8). Equivalently,
this is the power correspondence of the primary Example
8.5, here with its polynomial divisor specified.
The factorization through \(w\mapsto w^D\) gives the
asserted block system and quotient action.
Every absolute period in (4) is divisible by \(D\).
Thus the integer \(W_{ij}\) in a same-color relation
has a well-defined class modulo \(D\), giving (9).
The ratio is independent of the scalar choice of \(h\).

If the common \(f\)-value is \(v\), its \(D\) preimages
under \(w^D\) are cyclically indexed. The critical
points of \(h\) lying over each preimage supply exactly
one return permutation in that block. Changing the
choice of the first \(D\)-th root rotates all bins
for that color. This proves the return-type description.
\(\square\)

There is at most one finite branch with all cycle parts
divisible by some integer greater than one. Such a factor
is fixed-point-free and has index at least \(N/2\);
two would contradict total index \(N-1\). Thus the
power branch, when present, is intrinsically distinguished.
Its phase data do not assert that these are all block
systems or that phases classify every component.

## Theorem G4. Exact algebraic form of the cubic degree-eight pair

For the passport \((4,2,2),(2,2,1^4)\), every polynomial,
up to source affine change and target affine change,
is one of (1)--(2). The power-phase profiles of (1) are
both \((T,T)\); that of (2) is \((\mathrm{DT},1)\).
The two polynomials (1) have trivial deck group and
represent the two literal components of size \(16\)
with the same embedded group of order \(576\).
The polynomial (2) has deck order two and represents
the literal component of size \(8\), of group order \(64\).
The two components (1) have identical intrinsic branch
classes and opposite nontrivial cubic values, as proved
in [D5](POLYNOMIAL_INTRINSIC_CLASSES_AND_LIFTING.md#theorem-d5-a-cubic-area-invariant-explains-the-eight-sheet-coupling).

**Proof.** The unique multiplicity-four root can be put
at zero. G3 gives \(f=H^2\) with
\[
                    H=z^2(z^2+az+b).                        \tag{10}
\]
The requirements on its root fiber are \(b\ne0\) and
\(a^2-4b\ne0\). Its two remaining critical points
\(r,t\) satisfy
\[
                    4z^2+3az+2b=0,
\]
so they are distinct when \(9a^2-32b\ne0\).
Reduce \(H\) modulo this quadratic:
\[
 H(z)=Az+B,\quad
 A=\frac{a(9a^2-32b)}{64},\quad
 B=\frac{3a^2b}{32}-\frac{b^2}{4}.
                                                               \tag{11}
\]
The grouped branch condition is \(H(r)^2=H(t)^2\),
with both values nonzero. There are two cases.

If \(H(r)=H(t)\), (11) and distinctness of \(r,t\)
give \(a=0\). Source and target scaling set \(b=1\),
giving (2). Both critical points lie in the same
phase bin for the squaring map.

If \(H(r)=-H(t)\), their sum is
\[
 H(r)+H(t)=
         \frac{-27a^4+144a^2b-128b^2}{256}.                  \tag{12}
\]
Here \(a\ne0\), since \(a=0\) would force \(b=0\).
Normalize \(a=1\) by source scaling. Equation (12)
has exactly the two solutions
\[
                         b=\frac{9\pm3\sqrt3}{16}.
                                                               \tag{13}
\]
Each satisfies all root, critical-point and nonzero-value
conditions. The two points occupy opposite phase bins.

Any source affine equivalence preserves the unique
multiplicity-four root and hence fixes zero.
The remaining source scaling preserves \(b/a^2\).
Thus the two solutions (13) are inequivalent.
If \(a\ne0\), comparing the degree-eight and
degree-seven coefficients of \(f(\gamma z)=f(z)\)
forces \(\gamma=1\). If \(a=0\), the degree-eight
and degree-six coefficients give exactly
\(\gamma=\pm1\). These are the stated deck groups.

There are therefore exactly three rigid absolute
classes with the distinguished branch type first.
Allowing both type placements, a full two-strand
orbit has \(2N/|\mathrm{Deck}|\) literal tuples.
This gives \(16,16,8\). The independently retained
L8/D5 permutation certificates identify the two
deck-trivial classes with the same group \(H_{576}\)
and the even polynomial with the group of order
\(64\). D5 proves the intrinsic-class agreement and
genuine reduced cubic separation. This completes
the algebraic identification. \(\square\)

## Connected ambient fibers can have split diagonals

In G4 the two zeros of the logarithmic derivative
have order one. The ambient fixed-residue fiber is
connected by the primary isoresidual component
theorem. Its integral diagonal nevertheless has
the three rigid classes just computed. Even after
fixing the phase \((T,T)\), it still has the two
genuinely cubic-separated classes.

The phenomenon also occurs in positive dimension.
For every \(b\ge2\), the
[cyclic cubic family D6](POLYNOMIAL_INTRINSIC_CLASSES_AND_LIFTING.md#theorem-d6-an-infinite-family-with-identical-intrinsic-classes)
has passport
\[
           (2^{2b}1^{3b})^3,\quad(b1^{6b}),\qquad N=7b.
\]
Distinguish one of the first three branches.
The ambient differential has \(4b\) zeros of order
one and one of order \(b-1\), so its labeled
fixed-residue fiber is connected. Grouping the
other two copied branches imposes \(4b-2\)
independent complex period equations. The resulting
diagonal has dimension two, as G1 predicts.
Within the specified cyclic block-size-\(b\) sector
it contains two same-group, same-intrinsic-class
components distinguished by the cubic invariant.
Their root gcd is one, so G3 gives no nontrivial
phase label. No ambient-connectivity assertion can
replace classification of these diagonals.

## Evidence and open work

The [checker](verify_polynomial_period_diagonals.py),
[receipt](polynomial_period_diagonals_verification.json)
and [source ledger](polynomial_period_diagonals_source_ledger.json)
retain exact quadratic-field calculations, the quartic
remainder and sum identities, a direct dihedral polynomial
control, power-phase return profiles, finite orbit and
cubic witnesses, rational slice controls and source pins.
They do not implement the global graph construction G2.
That construction depends on the explicitly imported
primary period model and the written bound and slicing
argument. All current audit and replay are internal.
Independent replay, external review and formal
verification remain outstanding.
