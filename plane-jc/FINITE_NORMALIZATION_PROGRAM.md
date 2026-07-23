# Finite normalization in dimension two

> **Status.**  The finite-flatness theorem in Section 1 is unconditional.
> It applies to every plane Keller map over an algebraically closed field of
> characteristic zero and removes the closed-point nonflatness obstruction
> that occurs for three-dimensional finite normalizations.  The divisor and
> log-surface programme in Sections 3--6 is a research programme, not a proof
> of JC(2).

This is a second plane-Jacobian programme, parallel to the
[degree frontier](DEGREE_FRONTIER_125.md).  The degree programme starts from
published minimal-standard-pair and Newton-polygon reductions.  The present
programme starts directly from an arbitrary hypothetical counterexample and
uses its canonical Zariski--Main finite normalization.

Work over an algebraically closed field \(k\) of characteristic zero.  Let

\[
 F=(P,Q):U=\mathbb A^2_{x,y}\longrightarrow
 Y=\mathbb A^2_{u,v}
\]

satisfy \(J(P,Q)\in k^\times\).  Put

\[
 A=k[u,v]\simeq k[P,Q],\qquad K=k(x,y),\qquad
 B=\operatorname{Norm}_A(K),
\tag{0.1}
\]

and write

\[
 U\mathop{\lhook\joinrel\longrightarrow}^{j}
 \bar X=\operatorname{Spec}B
 \mathop{\longrightarrow}^{\pi}Y
\tag{0.2}
\]

for the canonical Zariski--Main factorization.  The first arrow is an open
immersion and the second is finite.  The degree
\(d=[K:\operatorname{Frac}(A)]\) is the geometric degree of \(F\).

## 1. The surface finite-flatness theorem

### Theorem 1.1 -- the canonical cover is finite free

The \(A\)-algebra \(B\) is finite locally free of rank \(d\).  Since
\(A=k[u,v]\), it is in fact free as an \(A\)-module:

\[
 B\simeq A^{\oplus d}.
\tag{1.1}
\]

Thus every scheme-theoretic fiber of \(\pi\) has length exactly \(d\).
There is no zero-dimensional Fitting defect:

\[
 \operatorname{Fitt}^{A}_{d}(B)=A.
\tag{1.2}
\]

#### Proof

The ring \(B\) is a finite normal \(A\)-algebra and has dimension two.
Every normal Noetherian surface is Cohen--Macaulay: normality gives
\(S_2\), and in dimension two this is the full Cohen--Macaulay depth
condition.

Let \(\mathfrak p\subset A\), and let \(\mathfrak q\subset B\) lie above it.
The integral finite extension has
\(\dim B_{\mathfrak q}=\dim A_{\mathfrak p}\).  A regular system of
parameters of \(A_{\mathfrak p}\) is a system of parameters in every
\(B_{\mathfrak q}\).  Since these local rings are Cohen--Macaulay, that
sequence is regular on every localization of the finite semilocal
\(A_{\mathfrak p}\)-module \(B_{\mathfrak p}\).  Hence

\[
 \operatorname{depth}_{A_{\mathfrak p}}B_{\mathfrak p}
 =\dim A_{\mathfrak p}.
\]

Auslander--Buchsbaum over the regular local ring \(A_{\mathfrak p}\)
(equivalently, miracle flatness) makes \(B_{\mathfrak p}\) free.  This holds
for every \(\mathfrak p\), so \(B\) is finite locally free.  Its rank is its
generic rank \(d\).  Quillen--Suslin then makes the projective \(A\)-module
\(B\) free.  The fiber-length and Fitting assertions follow.

The argument uses only normality, dimension two, and the regularity of the
target.  The Keller condition is used to construct the distinguished
étale open \(U\), not to prove the flatness of \(\pi\).

### Dimension comparison

For a finite normalization over a regular threefold, normality gives only
\(S_2\), hence depth at least two; maximal Cohen--Macaulay depth would be
three.  A closed-point nonflatness scheme can therefore remain.  Over the
regular surface \(A\), \(S_2\) already is maximal depth.  The obstruction
studied in the
[cubic normalization frontend](../cancellation/CUBIC_NORMALIZATION_FRONTEND.md)
is consequently specific to dimension at least three and disappears
automatically here.

This is stronger than generic flatness and stronger than flatness away from
finitely many points.  No curvilinearity, Cartier-boundary, or collision
fiber hypothesis is needed in dimension two.

## 2. Canonical consequences of the distinguished \(\mathbb A^2\)

Let

\[
 E=\bar X\setminus U,\qquad
 E^{(1)}=E_1\cup\cdots\cup E_r
\tag{2.1}
\]

be the divisorial part of the missing boundary.

### Proposition 2.1 -- the boundary freely generates the class group

Restriction along \(j\) gives

\[
 B^\times=k^\times,\qquad
 \operatorname{Cl}(B)\simeq
 \bigoplus_{i=1}^{r}\mathbb Z[E_i].
\tag{2.2}
\]

In particular:

1. the classes of the missing-boundary primes are independent;
2. every divisor class on \(\bar X\) has a unique boundary expression;
3. no nonzero divisor supported on \(E^{(1)}\) is principal; and
4. if \(E^{(1)}\) is empty, then \(\operatorname{Cl}(B)=0\).

#### Proof

The localization sequence for the normal scheme \(\bar X\) and its open
subscheme \(U\) is

\[
 B^\times\longrightarrow
 \Gamma(U,\mathcal O_U)^\times\longrightarrow
 \bigoplus_i\mathbb Z[E_i]\longrightarrow
 \operatorname{Cl}(B)\longrightarrow
 \operatorname{Cl}(U)\longrightarrow0.
\tag{2.3}
\]

Restriction embeds \(B^\times\) into
\(\Gamma(U,\mathcal O_U)^\times=k^\times\), and constants give the reverse
inclusion.  Also \(\operatorname{Cl}(U)=\operatorname{Cl}(\mathbb A^2)=0\).
Thus the middle divisor map in (2.3) is an isomorphism.

The boundary also reconstructs the polynomial source ring inside \(K\):

\[
 k[x,y]
 =
 \{h\in K:v_D(h)\ge0
 \text{ for every prime divisor }D\notin\{E_1,\ldots,E_r\}\}.
\tag{2.4}
\]

Equivalently, it is the union of the fractional ideals allowing arbitrary
effective pole order along the \(E_i\).  This valuation description is
canonical and does not require a choice of source coordinates.

### Proposition 2.2 -- ramification is missing boundary

The map \(F=\pi|_U\) is étale.  Hence every codimension-one ramification
prime of the finite flat cover \(\pi\) is one of the \(E_i\).  Write

\[
 e_i=e(E_i/\pi(E_i)),\qquad
 f_i=[k(E_i):k(\pi(E_i))].
\tag{2.5}
\]

Then:

- \(E_i\) is a **branch boundary** when \(e_i>1\);
- \(E_i\) is a **missing unramified boundary** when \(e_i=1\);
- each \(\pi(E_i)\) is an irreducible affine plane curve; and
- over its generic point the full degree ledger satisfies
  \[
  \sum_{D\mid \pi(E_i)}e(D/\pi(E_i))f(D/\pi(E_i))=d,
  \tag{2.6}
  \]
  where the sum includes both boundary primes and affine primes meeting
  \(U\).

Zariski--Nagata purity implies that a nontrivial cover cannot hide all of
its branching in closed points.  If \(\pi\) were étale everywhere, the
geometric simple connectedness of \(\mathbb A^2\) would force the connected
cover to have degree one.  Thus a hypothetical counterexample of degree
\(d>1\) has at least one branch-boundary divisor.

### Proposition 2.3 -- exact nonproperness and the affine-sheet budget

Let \(S_F\subset Y\) be the set of points at which \(F\) is not proper.
Then

\[
 S_F=\pi(E).
\tag{2.7}
\]

For an irreducible component \(C\subset S_F\), put

\[
 b_C=\sum_{E_i\mapsto C}e_i f_i
\tag{2.8}
\]

and let

\[
 a_C=\sum_{\substack{D\mapsto C\\D\cap U\ne\varnothing}}
 f(D/C)
\tag{2.9}
\]

be the degree carried by affine sheets.  Since \(F\) is étale, every affine
row has ramification index one.  The generic finite-flat degree equation is

\[
 d=a_C+b_C.
\tag{2.10}
\]

Moreover,

\[
 a_C\ge1,\qquad 1\le b_C\le d-1.
\tag{2.11}
\]

#### Proof

If \(y\notin\pi(E)\), the finite fiber \(\pi^{-1}(y)\) lies in the open
set \(U\).  Since \(\pi(E)\) is closed, the same holds over a neighborhood
of \(y\), and there \(F\) is the restriction of the finite map \(\pi\).
Thus \(F\) is proper at \(y\).

Conversely, suppose that \(F\) were finite over a neighborhood \(V\) of
\(y\).  The ring of \(F^{-1}(V)\) is normal, finite over
\(\mathcal O_Y(V)\), and has fraction field \(K\).  By uniqueness of
normalization it equals the restriction of \(B\) to \(V\).  Hence
\(\pi^{-1}(V)\subset U\), contradicting the existence of a point of \(E\)
above \(y\).  This proves (2.7).

Let \(g\in k[u,v]\) be a prime equation of \(C\).  The polynomial
\(g(P,Q)\in k[x,y]\) is nonconstant and therefore has a height-one prime
factor.  Its zero curve maps dominantly to \(C\): the quasi-finite map
cannot contract it.  This supplies an affine row and proves \(a_C\ge1\).
The Keller condition makes every affine row unramified.  Applying the
fundamental equality for the finite flat cover over the generic point of
\(C\) gives (2.10), and (2.11) follows.

### Corollary 2.4 -- the first geometric-degree exclusions

A hypothetical plane Keller counterexample has geometric degree

\[
 d\ge3.
\tag{2.12}
\]

If \(d=3\), every target curve carrying ramification has exactly:

1. one boundary row \((e,f)=(2,1)\); and
2. one affine row \((e,f)=(1,1)\).

There is no room over that target curve for another boundary or affine row.

#### Proof

A nontrivial connected finite cover \(\pi\) must have a ramified boundary
prime, by purity and the absence of nontrivial connected finite étale covers
of \(\mathbb A^2\).  Such a row contributes at least two to \(b_C\), while
Proposition 2.3 gives \(a_C\ge1\).  Thus \(d\ge3\).

When \(d=3\), equality forces \(b_C=2\) and \(a_C=1\).  The only ramified
boundary row of contribution two is \((2,1)\), and the only affine ledger of
total degree one is the single row \((1,1)\).

This is a geometric-degree statement, not the total-coordinate degree
frontier at \(125\).

### Proposition 2.5 -- missing-boundary curves are rational

For every \(E_i\), the normalization of its projective closure in \(X^c\)
is \(\mathbb P^1\).  Consequently the normalization of the affine curve
\(E_i\) is

\[
 \mathbb P^1\setminus\{p_{i1},\ldots,p_{is_i}\}
\tag{2.13}
\]

for a nonempty finite set of punctures over the target line at infinity.

#### Proof

Resolve the pair \((X^c,X^c\setminus U)\).  The strict transform of \(E_i\)
is a component of the SNC boundary of the smooth completion
\(U=\mathbb A^2\).  Every component of such a boundary is a rational curve.
The strict transform is the normalization of the projective closure of
\(E_i\), proving the assertion.  Its affine part removes exactly the points
lying outside \(\bar X=\Pi^{-1}(\mathbb A^2)\).

Thus a boundary row is not an arbitrary finite cover of its target curve:
its source normalization has genus zero, and its puncture number is part of
the canonical log ledger.

## 3. The branch and missing-boundary census

The first classification object should be the decorated finite-flat cover

\[
 \mathcal N(F)=
 \bigl(B/A,\ U\subset\bar X,\ 
 \{E_i,\pi(E_i),e_i,f_i\}_{i=1}^{r}\bigr).
\tag{3.1}
\]

For every target curve \(C=(g=0)\subset Y\), factor the Cartier pullback

\[
 \operatorname{div}_{\bar X}(g)
 =
 \sum_{D\mid C}e(D/C)D.
\tag{3.2}
\]

The affine terms are closures of components of
\((g(P,Q)=0)\subset U\); the remaining terms are missing-boundary primes.
Equation (3.2), together with (2.2), is the basic divisor ledger.  It keeps
three distinctions that a branch curve alone loses:

1. ramified boundary sheets versus unramified missing sheets;
2. boundary sheets versus affine sheets over the same target curve; and
3. normal degree \(f_i\) versus ramification index \(e_i\).

The first concrete objective is an exhaustive list of possible ledgers
compatible with finite flatness, the Keller open, and degree \(d\).  Useful
finite-algebra invariants are the trace pairing, discriminant ideal,
Kähler different, conductor of any intermediate order, and the Fitting
ideals of \(\Omega_{B/A}\).  Unlike the threefold frontend, no Fitting ideal
of the underlying \(A\)-module \(B\) remains to be controlled.

### External nonproperness input

For \(k=\mathbb C\), the Jelonek--Lasoń theorem says that every irreducible
component of \(S_F\) is a parametric curve.  If

\[
 \delta=\max\{\deg P,\deg Q\},
\tag{3.3}
\]

it admits a polynomial parametrization of degree at most \(\delta-1\).
Their argument extends to every algebraically closed characteristic-zero
field by the Lefschetz principle.  In particular, the normalization of every
target component \(C\subset S_F\) is \(\mathbb A^1\), and its projective
closure has one place at infinity.

This creates a direct interface between the two JC(2) programmes:

\[
 \begin{array}{c}
 \text{Newton frontier: }\delta\ge125\\[2mm]
 \Downarrow\\[-1mm]
 \deg(\text{a polynomial parametrization of }C)\le\delta-1\\[2mm]
 \text{finite normalization: }
 \widetilde E_i=\mathbb P^1\setminus\{s_i\text{ points}\}
 \longrightarrow \widetilde C=\mathbb A^1 .
 \end{array}
\tag{3.4}
\]

After projective completion, the last arrow is a degree-\(f_i\) map
\(\mathbb P^1\to\mathbb P^1\).  If its points over infinity have local
degrees \(m_{i1},\ldots,m_{is_i}\), then

\[
 \sum_{j=1}^{s_i}m_{ij}=f_i,
 \qquad
 \sum_{p\in\mathbb P^1}(\operatorname{ram}_p-1)=2f_i-2.
\tag{3.5}
\]

Equations (2.10), (3.5), and boundary adjunction should be compiled
together.  They couple the finite-flat sheet budget to the puncture and
intersection matrix rather than treating the nonproperness curve only as a
reduced target equation.

### Proposition 3.1 -- arbitrary-puncture residue rigidity

For a boundary prime \(E_i\), let \(s_i\) be its puncture number and \(f_i\)
its residue degree over its target nonproperness curve.  The ramification
forced on the affine normalization of \(E_i\) is

\[
 R_i^{\mathrm{aff}}
 =
 (2f_i-2)-\sum_{j=1}^{s_i}(m_{ij}-1)
 =
 f_i+s_i-2.
\tag{3.6}
\]

Consequently, if the normalized residue map is immersive at every affine
point, then

\[
 f_i=1,\qquad s_i=1.
\tag{3.7}
\]

#### Proof

The puncture contribution to Riemann--Hurwitz is

\[
 \sum_j(m_{ij}-1)=\sum_jm_{ij}-s_i=f_i-s_i.
\]

Subtracting it from the total ramification \(2f_i-2\) gives (3.6).
Immersion on the affine normalization makes this number zero.  Since
\(f_i,s_i\ge1\), the only solution is (3.7).

This removes the former one- and two-puncture restriction from the residue
budget.  Once a boundary theorem supplies affine residue immersion, every
number of punctures is handled uniformly.  The exact regression is
`cas/puncture_profile_budgets` in
[`cas/plane_boundary_exclusion.py`](cas/plane_boundary_exclusion.py).

Reference: Z. Jelonek and M. Lasoń,
*Quantitative properties of the non-properness set of a polynomial map*,
arXiv:1411.5011, Theorems 1.2 and 3.2.

## 4. The projective log-surface package

Let \(X^c\) be the normalization of \(\mathbb P^2\) in \(K\).  The finite
map extends canonically:

\[
 \Pi:X^c\longrightarrow\mathbb P^2.
\tag{4.1}
\]

The affine normalization is
\(\bar X=\Pi^{-1}(\mathbb A^2)\).  Let \(H\) be the reduced divisor over the
line at infinity and put

\[
 D^c=H\cup E^{(1)}=X^c\setminus U
\tag{4.2}
\]

up to possible codimension-two subsets.  In characteristic zero the
codimension-one canonical formula on the normal surface is

\[
 K_{X^c}=\Pi^*K_{\mathbb P^2}
       +\sum_G(e_G-1)G.
\tag{4.3}
\]

Consequently the components over the target line at infinity cancel in the
log correction, while every internal missing component contributes its full
normal multiplicity:

\[
 K_{X^c}+D^c
 \sim
 \Pi^*(K_{\mathbb P^2}+L_\infty)
 +\sum_{i=1}^{r}e_iE_i
 =
 -2\Pi^*L_\infty+\sum_{i=1}^{r}e_iE_i.
\tag{4.4}
\]

Formula (4.4) is a central numerical constraint.  On a log resolution
\(\rho:S\to X^c\), all exceptional curves lie in the boundary because
\(U=\mathbb A^2\) is smooth.  The reduced SNC boundary
\(D=S\setminus U\) then has the standard completion properties of the
affine plane:

- its components are rational and its dual graph is a tree;
- their classes form a \(\mathbb Z\)-basis of \(\operatorname{Pic}(S)\);
- their intersection matrix is unimodular; and
- \(K_S+D\), the pullback of (4.4), and the resolution discrepancies must
  agree coefficient by coefficient.

The basis assertion follows from the same localization sequence as (2.3),
now using \(\mathcal O(U)^\times=k^\times\) and
\(\operatorname{Pic}(U)=0\).  It turns every divisor identity into an
integral lattice identity rather than merely a numerical equivalence.

## 5. Four structural contradiction engines

The programme should test the following engines in parallel on the same
decorated cover.

### 5.1 Class group and factorial-open constraints

Use the canonical isomorphism
\(\operatorname{Cl}(B)\simeq\mathbb Z^r\) to compare:

- relations forced by pullbacks of branch curves;
- classes forced by the different and the dualizing module
  \(\omega_B\simeq\operatorname{Hom}_A(B,A)\);
- constraints on \(\operatorname{Cl}(B)\) coming from a finite free
  \(A\)-algebra; and
- any proposed contraction or identification of boundary components.

Any proposed boundary atlas that makes a nonzero boundary combination
principal is impossible.

### 5.2 Intersection matrices

Resolve the normal projective cover and express the pullback of a target line,
the internal boundary, the infinity boundary, the ramification divisor, and
the canonical class in the boundary basis of \(\operatorname{Pic}(S)\).
Reject graphs that fail unimodularity, adjunction, the Hodge index theorem,
or the required degree identity

\[
 (\Pi^*L)^2=d.
\tag{5.1}
\]

### 5.3 Units and punctures

Normalize every boundary curve and record its punctures at intersections
with the rest of \(D\).  The global unit rank of \(U\) is zero.  Local units,
norms, and residues on the punctured boundary curves must therefore cancel
in the global divisor lattice.  This is especially restrictive for
\(\mathbb G_m\)-type or multiply punctured boundary components.

### 5.4 Log canonical and ramification inequalities

Combine (4.4) with adjunction on each resolved boundary component.  The
effective ramification weights \(e_i\), the negative term
\(-2\Pi^*L_\infty\), and the rational-tree boundary must produce the log
canonical divisor of an \(\mathbb A^2\) completion.  Test log Kodaira
dimension, negativity of contractible subtrees, discrepancy bounds, and
Riemann--Hurwitz on every noncontracted boundary component.

These engines are complementary.  A dual graph that passes intersection
tests may still violate the class-group basis; a divisor ledger that passes
the class group may still have impossible puncture units or log
discrepancies.

## 6. Work programme and interface with the degree frontier

### FN2-A -- canonical algebra and branch divisor

Develop degree-independent presentations of the finite free algebra \(B/A\),
its trace pairing, different, and discriminant.  The goal is intrinsic
control, not a chosen free basis: Quillen--Suslin proves a basis exists but
does not make one canonical.

### FN2-B -- missing-boundary classification

Classify the possible collections
\((E_i,\pi(E_i),e_i,f_i)\), including affine sheets above the same target
curves.  Begin with the exact constraints

\[
 a_C+b_C=d,\quad a_C\ge1,\quad
 \widetilde E_i\subset\mathbb P^1,\quad
 \widetilde{\pi(E_i)}=\mathbb A^1,
\tag{6.2}
\]

and the Riemann--Hurwitz ledger (3.5).  Prove bounds on the number of
boundary components, their normal degrees, and their punctures.

### FN2-C -- log completion

Normalize the projective target compactification, resolve the source pair,
and compile the boundary tree, intersection matrix, canonical divisor,
different, and pullback line class.  The existing
[log-boundary compiler](LOG_BOUNDARY_COMPILER.md) provides a model for the
required ledger, but the input here is the canonical finite cover rather
than a published Newton polygon.

### FN2-D -- \(\mathbb A^2\) rigidity

Exploit simultaneously:

\[
 \mathcal O(U)^\times=k^\times,\quad
 \operatorname{Pic}(U)=0,\quad
 U\simeq\mathbb A^2,\quad
 \det(D_i\cdot D_j)=\pm1,
\tag{6.1}
\]

together with the log-canonical identity (4.4).  The desired theorem is that
no decorated finite-flat cover of degree \(d>1\) can contain such a
distinguished open.

### Interface with Newton polygons

The two programmes should exchange certificates without becoming logically
dependent:

- a Newton/Laurent candidate supplies explicit valuations and boundary
  centers to the finite-normalization compiler;
- the structural programme can reject a candidate before coefficient
  elimination;
- a surviving structural ledger tells the Newton programme exactly which
  missing sheets, ramification indices, and global boundary clusters its
  local polygons must realize.

The current degree statement remains:

\[
 \max\{\deg P,\deg Q\}\ge125
\]

conditional on the published standard-pair enumeration and the locally
reproduced final exclusions.  The finite-normalization theorem is
unconditional but, by itself, gives no new degree bound and does not prove
JC(2).  Its value is that it starts from every hypothetical plane
counterexample and replaces a potentially unbounded coefficient search by a
finite-flat surface and log-boundary classification problem.

## 7. Importing the lower-face and sparse-minimality methods

The recent Gaussian results contain two different proof mechanisms, and they
should be imported separately.

First, the lower-face proof of \(\operatorname{GMC}(2)\) selects an exposed
valuation face, proves that the face has a nonzero invariant, and uses prime
dilation to prevent higher terms from cancelling it.  The literal prime
step has no direct analogue for a fixed finite cover: there is no Frobenius
family of sheet ledgers.  The transferable principle is instead:

> choose an extremal boundary component by a linear functional, isolate its
> lowest divisor contribution, and prove that every possible cancellation
> term has strictly larger intersection or ramification cost.

Second, the sparse three-real-variable theorem proves a Pareto statement,
not an unsupported absolute minimum.  It fixes a bounded complexity
rectangle, enumerates finitely many supports, saturates away coefficient
deletions, and requires a structural all-order certificate for every
survivor.  The normalization analogue uses the local complexity vector

\[
 \kappa_C=
 \bigl(
 d,\ \#\{\text{boundary rows over }C\},\
 \#\{\text{affine rows over }C\},\
 \sum_i s_i,\
 \sum_i(f_i+s_i-2)
 \bigr).
\tag{7.1}
\]

At fixed \(d\), the possible coarse signatures are finite:

\[
 d=
 \sum_i e_if_i+\sum_j a_j,\qquad
 a_j\ge1,\qquad
 1\le s_i\le f_i.
\tag{7.2}
\]

The exact atlas and Pareto filter are implemented in
[`cas/finite_normalization_signatures.py`](cas/finite_normalization_signatures.py).
Through geometric degree eight the numbers of all, transversely ramified,
and ramified-residue-immersive signatures are:

| \(d\) | all | ramified | ramified and residue-immersive |
| ---: | ---: | ---: | ---: |
| 2 | 1 | 0 | 0 |
| 3 | 6 | 1 | 1 |
| 4 | 19 | 4 | 4 |
| 5 | 54 | 15 | 11 |
| 6 | 132 | 42 | 25 |
| 7 | 312 | 115 | 51 |
| 8 | 691 | 282 | 97 |

These are combinatorial possibilities, not candidate maps.  As in the GMC
minimality proof, a surviving finite support is only the start: it must be
promoted by an identity valid on the full geometry.

### Proposition 7.1 -- the residual-different identity

Here is the first such identity.  Work on a clean projective finite surface
chart where the target curve \(C\), a boundary curve \(E\) above it, and the
source and target surfaces are smooth.  Write

\[
 \pi^*C=eE+M,\qquad
 R_\pi=(e-1)E+R',
\tag{7.3}
\]

where neither \(M\) nor \(R'\) contains \(E\).  Let

\[
 g:E\longrightarrow C
\]

be the residue map of degree \(f\).  Then its ramification divisor satisfies

\[
 R_g\sim (R'-M)|_E,
\qquad
 \deg R_g=(R'-M)\cdot E.
\tag{7.4}
\]

If \(E\simeq C\simeq\mathbb P^1\), this becomes

\[
 \boxed{R'\cdot E=M\cdot E+2f-2.}
\tag{7.5}
\]

#### Proof

Adjunction and the finite canonical formula give

\[
\begin{aligned}
 K_E
 &=(K_X+E)|_E\\
 &=(\pi^*K_Y+(e-1)E+R'+E)|_E\\
 &=(\pi^*K_Y+eE+R')|_E.
\end{aligned}
\]

On the other hand,

\[
 g^*K_C
 =\pi^*(K_Y+C)|_E
 =(\pi^*K_Y+eE+M)|_E.
\]

Subtracting gives (7.4).  Riemann--Hurwitz gives
\(\deg R_g=2f-2\) in the rational case, proving (7.5).

### Corollary 7.2 -- an exposed ramification leaf

If \(R'\cdot E=0\), then

\[
 f=1,\qquad M\cdot E=0.
\tag{7.6}
\]

Thus a boundary component disjoint from all residual ramification can survive
only as a residue-degree-one sheet which is also disjoint from every other
sheet above the same target curve.  Any connectedness or conductor argument
forcing \(M\cdot E>0\) eliminates that extremal component immediately.

This is the surface counterpart of a noncancellable exposed face.  The next
structural task is to prove that a Pareto-minimal counterexample has such an
extremal ramification leaf, or else that every leaf pays a positive
residual-different cost whose sum exceeds the global ramification vector
\(K_S+3\Pi^*L\).
