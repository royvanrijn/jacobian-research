# Global low-degree support census below ((7,6,4))

## Status and exact result

This note introduces a support-first, valuation-aware pipeline for polynomial
Keller collisions in dimension three.  It makes no weighted, equivariant, or
linearity-in-one-variable assumption.

The cardinality-unbounded census below the invariant degree profile
((7,6,4)) is **not complete**.  The full normalized support space already has
a 105-digit number of collision-admissible labelled points.  The current exact
theorem is the complete sparse stratum:

> **Sparse degree-seven collision exclusion.**  Let (k) be a
> characteristic-zero field and let
> 
> 
> \[
> F:\mathbb A_k^3\longrightarrow\mathbb A_k^3
> \]
> 
> 
> satisfy (F(0)=0), (JF(0)=I), and (F(e_1)=0).  If every coordinate has
> ordinary degree at most seven, then the nonlinear part of (F), counted
> with coordinate multiplicity, contains at least seven monomials.

Equivalently, every normalized collision with invariant degree profile
lexicographically below ((7,6,4)) has at least seven nonlinear monomial
occurrences.  The proof exhausts exact supports of sizes at most six, not
coefficient boxes.  This is a lower bound: no support of size seven is
constructed or claimed to exist.

A separate dense calculation proves that the complete degree-at-most-two
collision ideal is the unit ideal.  Consequently exactly the following four
profiles in the requested range are currently eliminated without a support
bound:

\[
(1,1,1),\quad(2,1,1),\quad(2,2,1),\quad(2,2,2).
\]

No other full degree profile is declared empty.  No smaller counterexample,
unrestricted Pareto frontier, or cardinality-unbounded no-go theorem is
claimed.

## 1. Collision normalization and the invariant profile

Let (G) be Keller, let (a\ne b) satisfy (G(a)=G(b)), put
(A=JG(a)), and choose (L\in GL_3(k)) with (Le_1=b-a).  Then

\[
F(u)=(AL)^{-1}\bigl(G(a+Lu)-G(a)\bigr)
\tag{1.1}
\]

satisfies

\[
F(0)=0,\qquad JF(0)=I,\qquad F(e_1)=0,
\qquad\det JF=1.
\tag{1.2}
\]

Thus choosing (e_1) loses no collision; (e_2,e_3) are simultaneous
coordinate permutations.

The sorted degrees of the three displayed coordinates are not invariant
under (1.1), because its target factor can mix rows of different degrees.
The pipeline instead uses the filtered degree flag of

\[
V_F=\operatorname{span}_k\{F_1,F_2,F_3\}.
\]

It is the unique nonincreasing triple (d_1\ge d_2\ge d_3) such that, for
every (q),

\[
\operatorname{codim}_{V_F}
\bigl(V_F\cap k[x_1,x_2,x_3]_{\le q}\bigr)
=\#\{i:d_i>q\}.
\tag{1.3}
\]

After recentering at the normalized source and target origins, affine source
changes act by a filtered linear isomorphism on (V_F), while the linear part
of an invertible affine target change only chooses another basis of (V_F).
Thus (1.3) is affine-invariant.
There are exactly 74 nonincreasing positive triples lexicographically below
((7,6,4)).  They are listed in the stage-one artifact.

There is an important implementation consequence.  A normalized collision
frame need not be a degree-adapted basis of (V_F).  Support-only data cannot
decide whether top-degree coefficient rows cancel.  For a proposed profile
(d), the later coefficient algebra must impose

\[
\operatorname{rank} C_{>q}=\#\{i:d_i>q\},
\tag{1.4}
\]

where (C_{>q}) is the three-row coefficient matrix on all monomials of
degree greater than (q).  Vanishing minors impose the upper rank bound; a
Rabinowitsch equation for the ideal of required maximal minors imposes the
lower bound.  Treating the sorted coordinate degrees in (1.2) as an invariant
would silently omit maps.

## 2. The finite support universe

There are

\[
\binom{10}{3}-4=116
\]

nonlinear monomials of total degree two through seven in three variables.
They may occur independently in each output, giving 348 optional
support atoms before collision constraints.

Write (a_i) for the number of pure (x_1)-powers in the nonlinear support
of (F_i).  Exact-support feasibility of (F(e_1)=0) gives

\[
a_1\ge1,\qquad a_2\ne1,\qquad a_3\ne1.
\tag{2.1}
\]

Indeed,

\[
1+\sum_{x_1^m\in\operatorname{Supp}(F_1)}c_{1m}=0,
\qquad
\sum_{x_1^m\in\operatorname{Supp}(F_i)}c_{im}=0\quad(i=2,3).
\]

One nonzero coefficient cannot make either of the latter sums vanish.

The exact number of labelled supports satisfying only (2.1) is

\[
463548420604713283156113141997491731202027903828060298076586416157674159026322138160605594632849843027968.
\tag{2.2}
\]

After fixing (e_1), the only finite coordinate permutation retained is the
simultaneous swap (x_2\leftrightarrow x_3),
(F_2\leftrightarrow F_3).  Burnside's lemma gives

\[
231774210302356641578056570998745865601013951914035617533123602568962510567884904265612311840597244641280
\tag{2.3}
\]

orbits.  We do not quotient by the continuous affine stabilizer of (e_1):
it changes Newton supports by creating and cancelling monomials, so a naive
support quotient would not preserve completeness.

Equations (2.2)--(2.3) explain why writing every unrestricted support to a
flat JSON list is not a meaningful near-term computation.  The pipeline
stores the full Boolean universe and rank gates compactly, and exhausts
cardinality strata by closure.

## 3. Determinant buckets

For fixed exact supports (S_1,S_2,S_3), write

\[
F_i=\sum_{\alpha\in S_i}c_{i,\alpha}x^\alpha,
\]

including the fixed linear exponent (e_i) with coefficient one.  Choose
((\alpha,\beta,\gamma)\in S_1\times S_2\times S_3).  Summing the six column
permutations in the determinant shows that this triple contributes

\[
\det\!\begin{pmatrix}\alpha\\\beta\\\gamma\end{pmatrix}
c_{1,\alpha}c_{2,\beta}c_{3,\gamma}
x^{\alpha+\beta+\gamma-(1,1,1)}.
\tag{3.1}
\]

Thus the exponent in the request is correct, but the exponent determinant
in (3.1) is essential: triples with zero determinant contribute nothing.
The full degree-seven universe has 1,404,292 nonzero triples in 1,330
buckets; the largest bucket has 4,332 triples.

On an exact coefficient torus every displayed coefficient is nonzero.
Therefore a nonconstant bucket with one contribution cannot vanish.  The
support closure begins with each possible pure (x_1)-power in (F_1).  At
a singleton bucket it branches over every other determinant triple in the
same bucket.  If a support becomes balanced, it branches over every optional
next monomial, so balanced supersets are not missed.  A temporary single
pure-axis term in (F_2) or (F_3) is branched over every possible second
axis term before (2.1) is applied.

This gives the complete coefficient-free ledger:

| nonlinear support | collision-admissible supports | no-singleton supports | residual-symmetry orbits |
|---:|---:|---:|---:|
| 1 | 6 | 0 | 0 |
| 2 | 1,995 | 0 | 0 |
| 3 | 330,860 | 0 | 0 |
| 4 | 36,491,940 | 30 | 15 |
| 5 | 3,011,315,766 | 85 | 47 |
| 6 | 198,314,182,399 | 1,694 | 851 |

The closure visited 88,142 partial supports.  This is an exhaustive sparse
search, not a sample of the much larger rows in the second column.

The singleton rule and exact-support saturation generalize the plane method
in [the certified sparse JC(2) census](../plane-jc/CERTIFIED_SPARSE_SUPPORT_EXCLUSIONS.md).
The new ingredient is the three-row determinant bucket (3.1).

## 4. Infinity valuations are face classes

There are infinitely many (w\in\mathbb Z^3), so no finite program can emit
them literally.  For a fixed finite support only finitely many exposed-face
patterns occur.  The pipeline enumerates these equivalence classes exactly.

For each of the seven nonempty coordinate strata and each active integer
weight, put

\[
m_i(w)=\max_{\alpha\in S_i}\langle w,\alpha\rangle.
\]

If an escaping Puiseux arc has bounded image and (m_i(w)>0), its leading
coefficient in (F_i) must cancel.  Hence the exposed face of (S_i) must
contain at least two monomials.  If (m_i(w)\le0), no such multiplicity
condition is valid.  In particular, it would be incorrect to reject a
support merely because some generic weight exposes a vertex.

Z3 enumerates the face-selection Boolean patterns with unbounded integer
weights, then minimizes an integral representative in (L^1)-norm.  A
coordinate that is identically zero on the arc is handled by its coordinate
stratum rather than by a fictitious finite weight.  Among the 913 residual
support representatives, 687 have one admissible face class and 226 have
three.  None is eliminated at this necessary gate.

This polyhedral organization is consistent with the broader use of Newton
faces in algorithms for nonproperness; see El Hilany--Tsigaridas,
[Computing efficiently the non-properness set of polynomial maps in the
plane](https://arxiv.org/abs/2101.05245), and El Hilany,
[The tropical non-properness set of a polynomial
map](https://arxiv.org/abs/2207.00989).  Their hypotheses and conclusions are
not imported here.  Our multiple-face test is only a necessary valued-curve
gate, not a sufficient nonproperness theorem.

## 5. SMT, modular algebra, and exact lifting

For each fixed support, a Boolean SMT problem asks whether coefficient signs
can give both signs in every determinant bucket and every nontrivial collision
sum.  All 913 representatives are satisfiable.  This test is necessary over
ordered coefficient fields and has no exclusion force over (mathbb C).

Attach a coefficient variable to every nonlinear atom and put

\[
C=\prod c_{i,\alpha}.
\]

The exact coefficient-torus ideal is

\[
I_S=\bigl(
[x^\delta](\det JF-1),\ F(e_1),\ \rho C-1
\bigr).
\tag{5.1}
\]

No coefficient is normalized to one.  Singular computes (5.1) for every
residual support representative over

\[
\mathbb F_{11},\quad\mathbb F_{13},\quad\mathbb F_{17}.
\]

All (3\cdot913) modular ideals are unit ideals.  These are routing results,
not characteristic-zero proofs: a fixed special fiber can be empty even when
the generic fiber is nonempty.

The same 913 ideals are then computed over (mathbb Q).  Singular returns
((1)) in every case, and an independent SymPy exact Gröbner replay agrees
on every row.  Consequently no coefficient system survives exact lifting,
which proves the sparse theorem stated above.

For the dense quadratic row, start with all 18 quadratic coefficients in

\[
F_i=x_i+Q_i.
\]

The collision fixes the three (x_1^2)-coefficients to (-1,0,0).  The 19
nonzero determinant coefficient equations in the remaining 15 variables
have exact Gröbner basis ((1)).  This calculation includes every quadratic
support simultaneously and proves the four unrestricted profile exclusions.

## 6. Boundary accounting

Two different boundaries must not be conflated.

1. **Source infinity.**  The seven coordinate strata and all exposed-face
   equivalence classes are compiled before coefficient algebra in stage four.
2. **Coefficient support boundary.**  The Rabinowitsch equation in (5.1)
   checks the entire exact coefficient torus.  Setting a coefficient to zero
   gives a smaller exact support, and every smaller support through size six
   is a separate census row.

No affine coefficient chart was chosen and no coefficient was set to one.
Therefore a projective coefficient chart cannot resurrect a point after
(I_S=(1)): that equality already proves that the exact-support affine torus
is empty.  Projective coefficient charts are needed when a gauge normalization
must be covered or when compactifying a surviving component; neither occurs
in this closed sparse stratum.

Since exact lifting has no survivor, zero projective component charts remain
to audit.  “Nothing survives boundary analysis” in the artifacts refers only
to support at most six and to the dense quadratic row.

## 7. Reproduction and stage artifacts

Generate the pinned ledgers with

```bash
.venv/bin/python scripts/compile_global_low_degree_census.py
```

Recompute and compare every pinned decision with

```bash
.venv/bin/python scripts/verify_global_low_degree_census.py
```

The manifest is
[`global_low_degree_census_manifest.json`](../artifacts/generated-results/global_low_degree_census_manifest.json).
It pins eight JSON artifacts:

1. invariant degree profiles and their filtered-rank gates;
2. exact supports and symmetry orbits;
3. every determinant-bucket exponent and contributing exponent triple;
4. infinity face valuations;
5. sign SMT;
6. modular coefficient algebra;
7. exact rational Gröbner results;
8. source/coefficient boundary status.

The reference run used the versions recorded inside the artifacts.  The
verifier shares the support and valuation compiler with the generator, so it
is a deterministic recomputation rather than an independent implementation
of the full census.  Stage six and the primary stage-seven replay use
Singular; SymPy supplies an independent backend for the exact rational
Gröbner calculation.

## 8. Open global frontier

The present work does not justify replacing the weighted, equivariant, or
other template-specific searches.  Those searches probe supports far above
the completed sparse layer and retain their own exact conclusions.

The next global tasks are:

1. enumerate exact determinant-balanced supports of size seven, with
   restartable symmetry shards;
2. add the invariant rank conditions (1.4) before claiming profile-specific
   elimination;
3. strengthen the face gate from cardinality to common leading-coefficient
   torus solvability;
4. continue modular and exact algebra only on those survivors;
5. compactify any surviving positive-dimensional coefficient component and
   audit all of its genuine gauge charts.

Until these tasks terminate for every support cardinality, the honest global
answer is: four profiles are completely eliminated, support at most six is
globally eliminated in raw degree at most seven, and the unrestricted census
below ((7,6,4)) remains open.
