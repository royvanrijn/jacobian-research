# Finite valuation criteria for polynomial algebraization

This note tests the proposed principle that polynomiality in a
controlled-boundary incidence diagram should be decidable from finitely many
valuation, residue, and Rees-module conditions.

The conclusion is three-part.

1. For a **fixed finite reconstruction ansatz with a complete polar ledger**,
   the principle is true and elementary at its core. Normality reduces
   regularity to height-one valuations, and a fixed pole bound leaves only
   finitely many principal-part residues.
2. The unrestricted statement with only a normal affine source and an SNC
   compactification is false. SNC does not imply finite generation of the
   boundary-value semigroup, the associated graded algebra, or a divisorial
   Rees algebra. The needed finite-generation and strictness properties must
   be hypotheses or compiler certificates.
3. There is a structural positive class between these two extremes.  On a
   toroidal chart whose reconstruction module splits into saturated
   polyhedral contact-character modules, finite Khovanskii data and Rees
   strictness follow from the monoid presentation.  Polynomiality is then
   decided by finitely many facet inequalities and character-residue
   equations.

Thus the useful general theorem is a **finite valuation/Hilbert-basis
criterion with a strict Rees presentation**, not an automatic
finite-generation theorem for every SNC boundary.  In compiler language:

\[
\boxed{\text{complete polar ledger}
 +\text{ finite pole range}
 +\text{ Hilbert-basis degree data}
 +\text{ finite initial presentation}
 +\text{ Rees strictness}}
\]

turns the experimental residue tables into a proof.

## 1. The finite pole theorem

Let

\[
 X=\operatorname{Spec}A
\]

be a normal affine variety with function field \(K\). Let
\(U\subset X\) be an open chart, and let

\[
 E_1,\ldots,E_s
\]

be the codimension-one components of \(X\setminus U\), with normalized
valuations \(\nu_i=\operatorname{ord}_{E_i}\).

Fix rational reconstruction coordinates

\[
 f_1,\ldots,f_d\in K
\]

which are regular on \(U\).

### Proposition 1.1 -- complete polar ledger

The reconstruction is polynomial on \(X\) if and only if

\[
 \nu_i(f_j)\ge 0
 \qquad(1\le i\le s,\ 1\le j\le d).
\tag{1.1}
\]

#### Proof

A normal noetherian domain is the intersection of its height-one local
rings inside its fraction field:

\[
 A=\bigcap_{\operatorname{ht}\mathfrak p=1}A_{\mathfrak p}.
\tag{1.2}
\]

Every height-one prime not among the \(E_i\) meets \(U\), where each \(f_j\)
is already regular. At \(E_i\), the local ring is a DVR and membership is
equivalent to \(\nu_i(f_j)\ge0\). Hence (1.1) is equivalent to
\(f_j\in A\) for every \(j\). QED

No higher-rank or nondivisorial valuation is needed in Proposition 1.1.
Such valuations become relevant when the recorded model does not have a
complete polar ledger, or when one replaces element membership by integral
closure of an ideal or module.

The proposition also shows what the boundary compiler must certify: every
possible pole of every reconstructed coordinate must either occur on its
finite divisor list or be ruled out by regularity on the retained chart.
An SNC graph by itself is not that certificate.

### Proposition 1.2 -- denominators give a finite complete ledger

Suppose the reconstruction coordinates are presented as

\[
 f_j=\frac{a_j}{d_j},\qquad a_j,d_j\in A,\quad d_j\ne0.
\]

Let \(\mathfrak p_1,\ldots,\mathfrak p_s\) be the height-one primes occurring
in the Weil divisors of the \(d_j\).  Then

\[
 f_1,\ldots,f_d\in A
 \quad\Longleftrightarrow\quad
 \operatorname{ord}_{\mathfrak p_i}(f_j)\ge0
 \quad\text{for every }i,j.
\tag{1.3}
\]

If \(\pi:\widetilde X\to X\) is a proper birational model on which the
denominator divisor is resolved, it is enough—and is still
necessary—to check the strict transforms of these primes and every
exceptional divisor lying over their support.

#### Proof

At a height-one prime \(\mathfrak q\) outside the divisor of every \(d_j\),
each denominator is a unit, so every \(f_j\) lies in \(A_{\mathfrak q}\).
The remaining height-one primes are precisely the finite list
\(\mathfrak p_i\), and Proposition 1.1 proves (1.3).

For the resolved model, a regular function on \(X\) pulls back to a regular
function on \(\widetilde X\), proving necessity.  Conversely,
nonnegativity on the strict transform of each \(\mathfrak p_i\) gives
nonnegativity at \(\mathfrak p_i\), while primes outside the denominator
divisor were already harmless.  Proposition 1.1 again gives sufficiency.
The exceptional inequalities record cancellation at infinitely near
boundary centers and are redundant for membership in \(A\), but necessary
for regularity of the resolved reconstruction chart. QED

Thus “polar completeness” is not an additional finiteness conjecture for a
fixed rational formula.  It is certified by declaring its denominators and
resolving their finite Weil-divisor support.  The genuine extra hypotheses
enter only when the compiler is expected to cover an entire symbolic ansatz
without expanding its formulas.

## 2. Why leading residues are finite

Suppose now that the \(f_j\) vary in a finite ansatz. More precisely, let
\(R\) be a noetherian coefficient ring, let \(V\) be a finite free
\(R\)-module, and let

\[
 \phi:V\longrightarrow K^d
\]

be an \(R\)-linear reconstruction map.  Assume that the boundary valuations
are trivial on the coefficient field (as they are for ordinary family
parameters), so taking \(R\)-linear combinations does not introduce new
poles. Assume fixed pole bounds

\[
 \nu_i(\phi(v)_j)\ge-b_i
\tag{2.1}
\]

for all ansatz generators and coordinates, where \(b_i\ge0\).  Taking the
maximum over a finite generating set gives such a bound whenever the
reconstruction formulas are fixed.

At the generic point of \(E_i\), choose a uniformizer \(\pi_i\). The polar
part lies in the finite-length DVR module

\[
 \pi_i^{-b_i}A_{E_i}/A_{E_i}.
\tag{2.2}
\]

Consequently \(\phi(v)\) is regular at \(E_i\) exactly when its successive
coefficients at

\[
 \pi_i^{-b_i},\ldots,\pi_i^{-1}
\tag{2.3}
\]

vanish. If several ansatz monomials attain the same least valuation, the
first condition is precisely the vanishing of the sum of their initial
forms in the residue field \(k(E_i)\). If it vanishes, one passes to the
next graded piece. This process stops after at most \(b_i\) steps.

Thus a fixed ansatz has only

\[
 d\sum_{i=1}^{s}b_i
\tag{2.4}
\]

possible negative residue slots (some may vanish identically). They are
identities in residue fields, not merely numerical tests at sample points.
In a finitely presented coefficient family, clearing denominators and
reducing the finitely supported numerators in a presentation of \(k(E_i)\)
turns them into finitely many coefficient equations.

This observation is stronger and cleaner than saying that one checks only
the first tie. A cancellation at the lowest value may expose another
negative value, so all negative graded pieces down to value \(-1\) must be
checked. Finiteness comes from the pole bound.

## 3. The corrected semigroup--Rees criterion

A Hilbert basis of value vectors does not by itself control addition.
Distinct functions can have the same value and different initial residues,
and initial forms of chosen algebra generators need not generate the full
associated graded algebra. The appropriate strengthening is a finite
Khovanskii/SAGBI presentation together with filtered-module strictness.

Let \(\nu=(\nu_1,\ldots,\nu_s)\). Let \(C\subset K\) be the finitely
generated chart algebra used by the reconstruction formulas, and let
\(M\subset K^d\) be the finitely generated module of allowed reconstructed
expressions. Filter both by the divisorial valuations.  Write

\[
 \Gamma=\{\nu(g):0\ne g\in C\}\subset\mathbb Z^s
\]

for the relevant value semigroup.  In the saturated case a **Hilbert
basis** means the minimal finite generating set of
\(\operatorname{cone}(\Gamma)\cap\operatorname{gp}(\Gamma)\); in the
nonsaturated case the compiler must instead supply generators and relations
for \(\Gamma\) itself; if it uses the saturation for enumeration, it must
also list the holes that meet the ansatz's finite degree set.  This
distinction prevents the normalization of the semigroup from silently
replacing the actual chart algebra.

### Theorem 3.1 -- finite filtered-presentation criterion

Assume:

1. **polar completeness:** Proposition 1.1 applies to every coordinate in
   the ansatz;
2. **bounded support:** the ansatz has fixed finite pole bounds \(b_i\);
3. **finite degree data:** the value semigroup \(\Gamma\) is affine, given
   by a finite Hilbert basis in the saturated case or by a finite
   presentation of the actual semigroup, and the bounded ansatz uses only a
   finite declared set of value degrees;
4. **finite initial presentation:** chosen algebra and module generators
   form a finite Khovanskii/SAGBI basis, with a finite presentation of the
   associated graded algebra and module;
5. **Rees strictness:** the homogenized presentation maps onto the actual
   multi-Rees algebra and module with the declared kernel.  Saturation of
   the relation module with respect to the product of the boundary
   homogenizing variables is a sufficient finite certificate in the
   compiler presentations used here.

Then polynomiality of every reconstructed coordinate is equivalent to a
finite certificate consisting of:

- the finitely many semigroup degree inequalities occurring in the bounded
  ansatz;
- vanishing of the finitely many negative-degree initial residues, including
  every equal-valuation cancellation;
- the finite presentation and saturation/strictness certificate for the
  Rees module.

#### Proof

The Hilbert-basis and Khovanskii hypotheses reduce every ansatz expression
to a finite multigraded normal form and make every residue map computable
from finitely many generators and relations.  Fix \(E_i\).  For each
\(m=b_i,b_i-1,\ldots,1\), reduction modulo the next filtration step gives
the coefficient of \(\pi_i^{-m}\) in

\[
 \pi_i^{-b_i}A_{E_i}/A_{E_i}.
\]

Rees strictness says that a zero computed in the presented associated
graded module is the initial form of a genuine filtered cancellation, not
an artifact of a missing relation or homogenizing-variable torsion.
Successively killing these \(b_i\) residues is therefore equivalent to
\(\nu_i(f_j)\ge0\).  Repeating this independently for the finitely many
pairs \((i,j)\) gives all the inequalities (1.1).  Proposition 1.1 then
identifies them with \(f_j\in A\) for every \(j\).

Conversely, if every \(f_j\) lies in \(A\), its image in each negative
quotient (2.2) is zero.  The strict presentation records those zero initial
forms and hence satisfies the listed certificate. QED

The role of the clauses is now precise:

- the Hilbert basis controls which valuation degrees can occur and makes
  their multiplicative bookkeeping finite;
- the associated graded ring controls cancellation inside a degree;
- Rees strictness controls whether graded cancellations lift to the
  filtered module;
- the complete polar ledger converts the resulting inequalities into global
  regularity.

The Hilbert basis is therefore neither the regularity theorem nor the
cancellation theorem.  It is the finite indexing device between them.

For a fixed explicit rational expression, Proposition 1.1 and its finite
principal parts are usually enough; no Hilbert-basis theorem is needed.
The Rees formulation becomes useful when one wants one certificate for an
entire symbolic ansatz or wants the compiler to manipulate generators and
relations rather than expanded coordinates.

### Corollary 3.2 -- operational compiler criterion

For a fixed finite reconstruction ansatz over a normal affine source, the
following procedure is necessary and sufficient:

1. declare a common denominator for every reconstructed coordinate;
2. compile every height-one component of its divisor and, on the chosen
   resolved chart, every exceptional divisor above it;
3. compute the finite pole bound at each compiled divisor;
4. set every negative principal-part residue equal to zero, grouping all
   terms of equal order before taking the residue.

The reconstruction is polynomial exactly on the coefficient locus cut out
by those finitely many residue equations.

If the compiler replaces expanded formulas by semigroup normal forms, add
the Hilbert-basis, finite-initial-presentation, and Rees-strictness
certificates of Theorem 3.1.  These certificates justify the compression;
they do not add new geometric valuations.

#### Proof

Proposition 1.2 makes the divisor list complete.  Section 2 proves that the
negative principal parts form a finite list and that their vanishing is
equivalent to nonnegative order at each divisor.  Proposition 1.1 converts
those orders into membership in \(A\).  The last assertion is exactly
Theorem 3.1. QED

## 4. Smallest test: the three polynomiality theorems

The three established families satisfy Theorem 3.1 in a particularly simple
form: their relevant filtrations are principal, their pole orders are
bounded explicitly, and the filtered modules are free. Hence the Rees
modules are finitely generated and saturated without a separate global
finite-generation theorem.

### 4.1 Weighted admissibility

Use the notation

\[
 W=u\gamma,\qquad
 p=H',\qquad
 q=(Wp-H)/c,
\]

and

\[
 C=x\gamma,\qquad
 B=\frac{c+p(W)/\gamma}{x},\qquad
 A=\frac{u+q(W)/\gamma^2}{x^2}.
\]

There are two principal filtrations.

First, at \(\gamma=0\),

\[
 W\mid p(W),\qquad W^2\mid q(W)
\]

are exactly the nonnegative-valuation conditions which remove the
\(\gamma^{-1}\) and \(\gamma^{-2}\) terms. They are equivalent to

\[
 H'(0)=0,\qquad H(0)=H'(0)=0.
\]

Second, at \(x=0\), one has \((u,\gamma,W)=(1,1,1)\). The only negative
graded pieces are:

\[
\begin{array}{c|c|c}
\text{coordinate}&\text{required order}&\text{residue conditions}\\ \hline
B&1&c+p(1)=0,\\
A&2&q(1)=-1,\quad
1+\kappa+a_0(2+\kappa)=0,
\end{array}
\tag{4.1}
\]

where \(\kappa=p'(1)/c\). These are precisely

\[
 H(1)=0,\qquad H'(1)=-c,\qquad
 a_0=-\frac{1+\kappa}{2+\kappa}.
\]

Thus weighted admissibility is the complete finite principal-part
certificate. There are no further residues because the pole bounds are one
and two.

### 4.2 The cancellation jet

For the \((m,r)\) cancellation family, every positive \(z\)-degree term is
already regular. The entire polar part is

\[
 R|_{z=0}=\frac{Cx}{A^{r+1}}
 \Phi_{m,r}(A,h(A)).
\tag{4.2}
\]

Since \(A=1+xy^m\) is prime and \(A\nmid x\), Proposition 1.1 gives

\[
 R\in K[x,y,z]
\quad\Longleftrightarrow\quad
\Phi_{m,r}(A,h(A))\equiv0\pmod {A^{r+1}}.
\tag{4.3}
\]

The \(r+1\) coefficients of \(A^0,\ldots,A^r\) are exactly the successive
negative associated-graded residues. This is the operator
\(\mathcal L_{m,r}(h)\).

Its degree-zero equation is \(M_{m,r}(h(0))=0\). Since this polynomial is
squarefree, the derivative of the first residue with respect to \(h(0)\) is
nonzero at a root. The remaining residues therefore determine
\(h_1,\ldots,h_r\) successively. The “unique cancellation jet” is simply
the triangular lifting of the vanishing initial form through the principal
\(A\)-adic Rees filtration.

### 4.3 The quadratic-gauge jump

In the reciprocal chart

\[
 P=tq,\qquad S=x/t,
\]

a prospective monomial \(cP^\alpha S^k\) contributes \(t\)-orders

\[
\alpha-k+2\quad\text{and}\quad\alpha-k
\tag{4.4}
\]

to the slope and intercept coordinates. For an isolated term with
\(k\ge4\), the second value is the smaller one, so nonnegativity gives

\[
\boxed{\alpha\ge k.}
\tag{4.5}
\]

The minimal lift is \(\alpha=k\), recovering the coefficient-weight jump

\[
(w_1,w_2,w_3,w_4,w_5,\ldots)=(0,1,1,4,5,\ldots).
\]

Degrees one and three occupy the same coupled low-order filtered block.
Their initial residues cancel exactly when

\[
 h(t)\equiv\frac{g_1}{g_3}(1+3t)\pmod {t^2}.
\tag{4.6}
\]

Degree two has zero intercept initial form. Thus the apparent exceptions are
not exceptions to the finite criterion: they are its equal-valuation
residue block. For a more general lower-weight multi-term lift, one must
compute the finite initial-residue matrix rather than apply (4.5) term by
term.

These three arguments use one statement: a bounded rational ansatz over a
normal source is polynomial exactly when every negative principal part
vanishes. Their different-looking gates are respectively a two-prime jet,
an \(A\)-adic jet of length \(r+1\), and a \(t\)-weight initial-form
calculation.

## 5. Falsification tests for the broad form

### 5.1 An incomplete boundary ledger

Compactify \(X=\mathbb A^2=\operatorname{Spec}k[x,y]\) by
\(\mathbb P^2\), whose boundary is the SNC divisor \(L_\infty\). Then

\[
 f=\frac1x
\]

has positive order along \(L_\infty\), but it is not regular on \(X\)
because it has a pole along the affine divisor \(x=0\).

Thus nonnegative values on every *recorded compactification-boundary prime*
do not prove polynomiality. The ledger must also prove that no affine polar
prime was omitted.

### 5.2 Valuations detect integral closure, not module membership

In \(k[x,y]\), let

\[
 I=(x^2,y^2).
\]

The element \(xy\) satisfies the divisorial inequalities defining the
integral closure of \(I\), but

\[
 xy\notin I,\qquad xy\in\overline I.
\]

Hence even a complete Rees-valuation test can prove only integral
dependence, not membership in a nonsaturated presentation module. This is
exactly the gap addressed by the strictness/saturation clause in Theorem
3.1.

### 5.3 SNC does not imply finite generation

Finite generation of symbolic or saturated multi-Rees algebras is a genuine
extra theorem, not a formal consequence of regularity or SNC boundary.
There are prime ideals in polynomial rings with non-finitely-generated
symbolic Rees algebra, and valuation semigroups on smooth projective
varieties need not be finitely generated.

Accordingly, failure of the proposed general statement need not reveal a
mysterious higher-rank valuation. The first failure may already be:

1. an omitted ordinary prime divisor;
2. a value semigroup with no finite Hilbert basis;
3. failure of a chosen generating set to be a Khovanskii basis;
4. Rees torsion or nonsaturation, as in ideal membership versus integral
   closure.

## 6. Consequence for the boundary compiler

The plane boundary compiler can become theorem-producing once it emits five
checkable certificates:

1. **polar completeness:** the compiled model dominates the normalized graph
   of every reconstruction denominator in the declared ansatz;
2. **boundedness:** a finite Newton/pole box containing every candidate;
3. **finite degree data:** a Hilbert basis for the saturated value semigroup,
   or a presentation of the actual affine semigroup together with every
   hole relevant to the bounded ansatz;
4. **initial completeness:** a finite Khovanskii/SAGBI or Gröbner
   presentation for the boundary filtration;
5. **strictness:** saturation of the reconstructed-coordinate Rees module.

After these certificates, chart enumeration is no longer the proof.
Theorem 3.1 turns the compiled finite data into a uniform polynomiality
theorem.

Without them, the compiler remains a correct evaluator of supplied charts
but cannot infer that an unlisted divisor, initial generator, or Rees
relation does not exist.

## 7. A structural toroidal class

Theorem 3.1 is a certificate theorem: it explains what finite data suffice,
but does not say when the initial and Rees certificates exist.  There is,
however, a natural class in which both certificates are forced by the
geometry.  The essential hypothesis is stronger than saturation of the
contact monoid.  The reconstruction module must split into contact-character
pieces.

Let \(R\) be a noetherian normal domain, let \(P\) be a fine, saturated,
sharp, torsion-free affine monoid, and put

\[
 L=P^{\mathrm{gp}},\qquad \sigma=\operatorname{cone}(P).
\]

For a rational polyhedron \(Q\subset L_{\mathbb R}\) with recession cone
\(\sigma\), write

\[
 \Delta=Q\cap L.
\]

Then \(R[\Delta]\) denotes the free \(R\)-module with basis
\(\{\chi^u:u\in\Delta\}\), with its natural \(R[P]\)-module structure.
Gordan's lemma for modules makes \(R[\Delta]\) a finite \(R[P]\)-module.
This is the precise meaning used below for a saturated polyhedral character
module.  Bare saturation of an abstract module, without the character
splitting, is not enough.

### Theorem 7.1 -- toroidal character-module algebraization

Let \(X=\operatorname{Spec}A\) be normal affine, and fix a finite
reconstruction ansatz with a complete polar ledger
\(E_1,\ldots,E_s\) and fixed pole bounds.  Suppose that, at every boundary
stratum meeting the polar ledger, after strict henselization and completion:

1. the boundary chart, up to a formal power-series factor, is
   \[
      R[[P]];
   \]
2. the valuations of the boundary components are the primitive integral
   facet functionals
   \[
      \ell_i:L\longrightarrow\mathbb Z
   \]
   of \(\sigma\);
3. on the character factor, the chart algebra used by reconstruction is
   \(R[P]\), up to localization in degree-zero coefficients; and
4. the reconstruction module embeds homogeneously into the Laurent
   character module and is a finite direct sum
   \[
      \bigoplus_{j=1}^d R[\Delta_j]e_j,
      \qquad
      \Delta_j=Q_j\cap L,
   \tag{7.1}
   \]
   for rational \(P\)-polyhedra \(Q_j\), with flat coefficient modules if
   \(R\) is replaced by a finite module on the stratum.

Then polynomial algebraization of the ansatz is decided by finitely many
facet inequalities and finitely many negative character-residue equations.
Moreover, a Hilbert basis of \(P\), together with finite generators of the
\(\Delta_j\), gives a finite algebra-and-module Khovanskii basis, and its
homogenized module presentation is Rees-strict.  Thus the
initial-completeness and strictness clauses of Theorem 3.1 are automatic in
this class.

#### Proof

Since \(P\) is saturated,

\[
 P=\sigma\cap L
   =\{u\in L:\ell_i(u)\ge0\text{ for every facet }i\}.
\tag{7.2}
\]

The finitely many facet functionals therefore decide whether a character is
regular on the toroidal chart.  The Laurent character module has the direct
decomposition

\[
 R[L]^d=\bigoplus_{j=1}^d\ \bigoplus_{u\in L}R\chi^u e_j.
\tag{7.3}
\]

After terms with the same full character have been collected, cancellation
can occur only inside one summand of (7.3).  In particular, equality of the
numerical divisorial value vectors is not by itself enough: distinct
tangential characters remain distinct initial forms.  For the fixed pole
bounds, only finitely many negative degrees occur, so regularity is
equivalent to the vanishing of their finitely many coefficients.  If a
coefficient has its own boundary expansion, the same argument in successive
graded pieces gives the finite principal-part residue equations of Section
2.

Each \(\Delta_j\) is a finitely generated \(P\)-module by the rational
polyhedral form of Gordan's lemma.  Hence (7.1) has a finite homogeneous
semigroup-module presentation.  Its filtration is split by (7.3).
Homogenizing the complete monomial and binomial relations therefore presents
the actual multi-Rees module.  Equivalently, its kernel is saturated with
respect to the product of the homogenizing variables: localization embeds
the presentation into the free Laurent character module, where those
variables are nonzerodivisors.  Thus a graded cancellation is exactly the
initial form of a filtered cancellation, proving Rees strictness.

The finite inequalities and residues now give nonnegative order at every
prime in the polar ledger.  Proposition 1.1 converts these local conditions
into membership in \(A\). QED

### 7.2 What the contact-monoid theorem already supplies

The
[saturated contact-monoid theorem](../extended-geometry/COMBINATORIAL_COMPLETED_LOCAL_RINGS.md)
gives, on every selected phase branch and up to formally smooth factors,

\[
 \widehat{\mathcal O}^{\,\mathrm{norm}}_{x,\beta}
 \simeq k'[[Q_{x,\beta}^{\mathrm{sat}}]].
\tag{7.4}
\]

Thus it already proves hypotheses 1 and 2 of Theorem 7.1 for the normalized
completed atlas.  The remaining structural problem is narrower:

> prove that the universal reconstruction coordinates are finite sums of
> contact characters and that their relation module is a homogeneous,
> stratum-flat module of the form (7.1).

This is the direct bridge from the contact-monoid theorem to polynomiality.
The contact theorem constructs \(P\); the new calculation must construct the
exponent modules \(\Delta_j\).

The complete polar ledger remains necessary.  Formula (7.4) controls the
chosen boundary chart, not an affine polar divisor omitted from that chart.

### 7.3 Hierarchy of candidate structural classes

The proposed classes do not all supply the same part of Theorem 7.1.

| class | finite structure supplied | remaining obstruction |
|---|---|---|
| toric | saturated monoid, facet valuations, character basis, binomial relations | polar completeness and character form of the reconstruction |
| toroidal | the toric package étale-locally | descent, chart monodromy, and flatness over boundary strata |
| complexity-one \(T\)-variety | polyhedral-divisor and weight-space descriptions | vertical valuations on the base curve, round-down in section spaces, and module strictness |
| spherical | multiplicity-free weight data and available toric degenerations | compatibility of the chosen boundary valuation with a finite basis and with the reconstruction module |
| Mori dream/Cox chart | finite Cox generators and finitely many GIT chambers | a Cox generating set need not be a Khovanskii basis for the boundary valuation |
| finite-type cluster algebra | finitely many clusters and cluster variables | Laurentness is not regularity on the affine reconstruction chart; upper-cluster and polar-completeness issues remain |
| finitely generated Khovanskii degeneration | the algebra-side finite initial presentation | a module Khovanskii basis and Rees-strict lifting still have to be proved |

Consequently the classification target should concern filtered pairs, not
varieties alone:

> classify controlled-boundary reconstruction pairs \((C,M,\nu)\) for
> which the contact valuation has a finite algebra-and-module Khovanskii
> basis and the module degeneration is strict.

Theorem 7.1 gives the first nontrivial positive class.  Complexity-one
\(T\)-varieties are the natural next class because their weight spaces are
controlled by polyhedral divisors on a curve; the extra vertical valuations
make them the smallest test of whether the toric proof survives a
nontrivial coefficient base.

### 7.4 Minimal falsification suite

Any proposed extension of Theorem 7.1 should first survive five tests:

1. **omitted pole:** add an affine denominator divisor not represented by
   the toroidal boundary;
2. **same value, different character:** use two Laurent characters having
   the same selected divisorial values but different tangential weights;
3. **semigroup hole:** replace \(P\) by a nonsaturated affine semigroup and
   test an exponent in its saturation but not in \(P\);
4. **nonstrict module:** use \((x^2,y^2)\subset k[x,y]\), where valuation
   inequalities admit \(xy\) but module membership does not;
5. **monodromy:** glue toroidal charts with monodromy permuting character
   generators and test whether the local split modules descend.

These distinguish polar completeness, full initial-form completeness,
normality of the monoid, Rees strictness, and toroidal descent.  Passing
them would justify moving from the split toroidal theorem to a genuine
toroidal or complexity-one structural theorem.

## 8. External framework

The relevant external facts are:

- normal-domain membership is codimension-one membership; see the
  [Stacks Project, Lemma 15.24.18](https://stacks.math.columbia.edu/tag/0AVB);
- strict normal crossings gives local regular parameters and clean
  intersections, but no global finite-generation statement; see the
  [Stacks Project section on normal crossings divisors](https://stacks.math.columbia.edu/tag/0CBN);
- affine modifications are naturally expressed through Rees-type algebras;
  see Kaliman--Zaidenberg,
  [*Affine modifications and affine hypersurfaces with a very transitive automorphism group*](https://arxiv.org/abs/math/9801076);
- finite generation of saturated multi-Rees algebras requires substantive
  hypotheses; see Das--Roy,
  [*Finitely generated saturated multi-Rees algebras*](https://arxiv.org/abs/2112.14587);
- non-finitely-generated symbolic Rees algebras already occur for prime
  ideals in polynomial rings; see Sannai--Tanaka,
  [*Infinitely generated symbolic Rees algebras over finite fields*](https://arxiv.org/abs/1703.09121);
- valuation semigroups need not be finitely generated even in smooth toric
  geometry; see Altmann--Haase--Küronya--Schaller--Walter,
  [*On the finite generation of valuation semigroups on toric surfaces*](https://arxiv.org/abs/2209.06044).
- finite Khovanskii bases connect higher-rank valuations, tropical data,
  and toric degenerations; see Kaveh--Manon,
  [*Khovanskii bases, higher rank valuations and tropical geometry*](https://arxiv.org/abs/1610.00298);
- normal affine \(T\)-varieties are described by proper polyhedral divisors;
  see Altmann--Hausen,
  [*Polyhedral Divisors and Algebraic Torus Actions*](https://arxiv.org/abs/math/0306285);
- affine and polarized projective spherical varieties admit toric
  degenerations; see Alexeev--Brion,
  [*Toric degenerations of spherical varieties*](https://arxiv.org/abs/math/0403379);
- Mori dream spaces are controlled by finitely generated Cox/GIT data, but
  that result concerns the Cox grading rather than an arbitrary selected
  boundary valuation; see Hu--Keel,
  [*Mori Dream Spaces and GIT*](https://arxiv.org/abs/math/0004017);
- finite-type cluster algebras have finitely many clusters and are
  classified by finite root systems; see Fomin--Zelevinsky,
  [*Cluster Algebras II: Finite Type Classification*](https://arxiv.org/abs/math/0208229).

The resulting research target is therefore not “prove Rees finite
generation from SNC.”  Theorem 7.1 handles the split toroidal case.  Beyond
that case, the target is:

> identify the controlled-boundary diagrams for which the reconstruction
> algebra has a finite Khovanskii presentation and the reconstruction module
> is Rees-strict, then make those two properties part of the compiler
> certificate.
