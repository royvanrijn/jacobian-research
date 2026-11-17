# Filtered principal-pair transfer for GVC

15 November 2026. Internal written theorem; independent review outstanding;
no novelty claim. This extracts an all-degree principle behind the
[quadratic-phase proof](QUADRATIC_PHASES_FOR_ALL_OPERATOR_ORDERS.md).
The [registry](../../MATH_STATUS.json) is authoritative.

## 1. Exact scope

Choose positive integer weights w=(w_1,...,w_n) and a positive integer r.
Polynomial coefficients may lie in any commutative Q-algebra, including a
nonreduced one. Give a constant derivative monomial partial^alpha the same
weight w dot alpha as x^alpha. Assume
\[
 \deg_w P\le r,\qquad
 \text{every monomial of }\Lambda\text{ has derivative weight }\ge r.
\]
The operator has constant coefficients and may be a polynomial or a formal
power series in derivatives; it acts locally finitely on each polynomial.
Write
\[
 P=f+g,\quad \deg_w f=r,\ \deg_w g\le r-\delta_P;
 \qquad\Lambda=B+R,
\]
where B is the derivative-weight-r part and every term of R has weight at
least r+delta_Lambda. The gaps are positive integers. Zero parts are allowed;
for a zero remainder choose any positive gap, for example1. A declared degree
bound for a zero polynomial can likewise be zero. Multipliers are polynomials.

The conclusion is
\[
 \boxed{\ (\Lambda,P)\text{ is pure and GVC-safe}
 \iff (B,f)\text{ is pure and GVC-safe}.\ }                      \tag{1}
\]
Purity and safety actually transfer separately. Here safety means that every
fixed multiplier has an eventually zero marked sequence. No assumption that
every principal pair is safe is made: the retained counterexamples are
balanced principal pairs and remain counterexamples after these perturbations.

For context, [van den Essen--Willems--Zhao0903.1478v2,Lemma4.9](https://arxiv.org/html/0903.1478v2)
proves pure/mixed vanishing under a strict Newton-polytope separation.
Its statement and complete proof were read. The present argument addresses the
balanced boundary by retaining its principal pair, rather than declaring it
safe from support alone. The proof below uses only grading, binomial expansion
and the definition of eventuality; it imports no new external theorem.

## 2. Pure identities and ordinary pure ideals

Every term in \((B+R)^m(f+g)^m\) containing k factors g and l factors R has
input weight at most mr-k*delta_P and derivative weight at least
mr+l*delta_Lambda. Unless k=l=0 it is zero. Consequently the exact identity is
\[
 \Lambda^mP^m=B^mf^m\quad(m\ge1).                              \tag{2}
\]
The right side has weight zero and hence is a scalar in the coefficient ring.
Equation (2) is an identity before specialization, not just equality of zero
sets. Therefore every finite-prefix ordinary pure coefficient ideal, and the
full ordinary pure ideal, is the extension/pullback of its principal ideal.
For freely adjoined lower-phase and higher-operator coefficients the whole
pure scheme is a product, with nilpotents and embedded components retained.
Faithfully flat polynomial extension preserves the least generating prefix;
an additional specialization can shorten it. This does not identify ordinary
mixed cutoff ideals, which can change under these perturbations.

## 3. The principal marked component and the reverse implication

Let Q be weight-homogeneous of weight e. The full marked contraction has
weight at most e, and its weight-e component is exactly
\[
 [\Lambda^m(QP^m)]_e=B^m(Qf^m).                                \tag{3}
\]
All other contributions have strictly smaller weight, so cannot cancel it.
If the full pair is safe, (3) makes the principal pair safe for homogeneous Q;
decomposing any polynomial multiplier into finitely many weighted components
gives safety for every Q. Conversely, if a principal marker is nonzero
infinitely often, at least one of its homogeneous components is nonzero
infinitely often. That component gives the same obstruction for every allowed
perturbation. No division, field argument, or reducedness assumption is used.

## 4. The finite-defect argument for the forward implication

For an arbitrary fixed Q, put e=deg_w Q. Expand both powers. A term indexed
by k lower-phase insertions and l higher-operator insertions can contribute
only if
\[
 k\delta_P+l\delta_\Lambda\le e.                              \tag{4}
\]
There are finitely many such pairs, independently of m. For each pair set
\[
 s=m-\max(k,l),\quad
 M_{k,l}=Qg^k f^{\max(l-k,0)}.
\]
The corresponding differentiated expression, apart from its binomial scalar,
is exactly
\[
 R^l B^{\max(k-l,0)}
       \bigl[B^s(M_{k,l}f^s)\bigr].                            \tag{5}
\]
The extra powers commute because all operators have constant coefficients.
If (B,f) is safe, each fixed marker M has an eventually zero sequence. For
large enough m, all the finitely many brackets in (5) vanish. Terms outside
(4) already vanish by weight. This proves forward safety in (1).
The proof applies over the stated coefficient ring and to formal operator
symbols, because only finitely many derivatives act on each polynomial.

If deg_w P<r and there is no principal phase, strict degree counting already
gives vanishing for m>e/(r-deg_w P). If the principal operator vanishes, the
analogous strict gap does the same. These zero-principal cases also fit the
identities above; they are included on coefficient-family boundaries.

## 5. A transferred linear cutoff and coefficient-family reduction

Suppose the principal pure-safe pair has a bound
\[
 B^j(Mf^j)=0\quad\text{for }j>C\deg M,
\]
with C>0. This is a conditional input to this quantitative statement, not a
claim of a coefficient-uniform C for arbitrary principal families. Let d=deg Q,
let D_g,D_f bound ordinary degrees of g,f, and set
\[
 \delta=\min(\delta_P,\delta_\Lambda),\qquad
 A_* =\max(D_g/\delta_P,D_f/\delta_\Lambda).
\]
In (5), deg M<=d+kD_g+max(l-k,0)D_f. Equation (4) bounds this by d+eA_*,
and max(k,l)<=e/delta. Thus a sufficient full cutoff is
\[
 m>C d+e(1/\delta+C A_*)\quad\text{and hence}\quad
 \boxed{m>C'\deg Q,\quad
 C'=C+w_{\max}(1/\delta+C A_*).}                              \tag{6}
\]
This is generally nonoptimal. Constant multipliers use (2). The finite maximum
of max(k,l)+C*deg M over the budget (4) can give a smaller marker-specific
bound. The earlier quadratic-phase proof used stronger shifted-power data and
can therefore retain its better bound; (6) does not replace it.

For a free coefficient family of lower phases and higher operator jets,
the pure-safe field locus is exactly the inverse image of the principal
pure-safe locus under coefficient projection. A zero-perturbation section
gives the converse implications. Therefore closedness and constructibility
are each equivalent for the two loci. With fixed weights,r and positive gaps,
ordinary degrees D_f,D_g are bounded, so a common principal C gives a common
full C' by (6). Conversely a full family bound restricts to its principal
section. Thus field uniformity and closed-locus questions in this balanced
architecture reduce exactly to principal data. This does not solve them for
arbitrary principal pairs or pairs admitting no such positive-weight balance.

## 6. Positive classes and rigid failure families

If the principal operator has a split homogeneous symbol, retained GVCSS
makes every pure principal phase safe, so (1) closes its whole filtered family.
Another positive class uses the retained binary theorem GVC2ENV: if the
principal operator differentiates in at most two independent directions,
treat passive variables as coefficient parameters over their fraction field.
Its arbitrary-degree principal phases are safe, and the resulting identities
descend to the original polynomial ring. The filtered theorem then allows
higher jets involving additional variables. These positive applications are
over characteristic-zero fields and retain the assurance limits of the
imported repository theorems; they are not an assertion of those positive
theorems over arbitrary nonreduced rings.

Two negative applications require no search:

* DVG1 has a homogeneous cubic phase and homogeneous order-three operator.
  Add any polynomial of degree<=2 to the phase and any derivative terms of
  order>=4 to the operator, including formal higher jets. All pure identities
  remain zero. With its homogeneous marker -c, the degree-one marked component
  is still \((-1)^{m+1}(m+1)!(m!)^2t\) for every m>=2. Every member fails.
  Arbitrary nonzero Q-algebras, including nilpotent perturbation coefficients,
  cannot erase this unit-coefficient top component.
* The retained ternary GVC3HC phase has degree12 and operator Delta^6 of
  order12. Every phase perturbation of degree<=11 and operator perturbation
  of order>=13 preserves its pure nullity and nonzero degree-two marker
  component. Taking the operator perturbation zero gives arbitrary lower-degree
  deformations of the fixed Delta^6 witness. No new Delta power is produced.

These are full infinite coefficient families, not finite-prefix inferences.
They leave protected minimum-degree/dimension records unchanged. In particular
lower terms or higher jets cannot repair a fixed bad principal pair; changing
the principal data, the balance hypotheses, or the class of operators is needed.

## Verification and boundaries

The [checker](verify_filtered_principal_transfer.py),
[receipt](filtered_principal_transfer_verification.json) and
[source audit](filtered_principal_transfer_source_audit.json) retain finite
controls of literal pure equality, marked top components, the defect budget,
ordinary-ideal independence, a nilpotent coefficient specialization and the
quantitative inequality. A small bounded DVG1 perturbation check is labelled
as such; its all-order proof remains the retained canonical result. No new
ternary witness producer or general principal-classification solver is run.

Independent review of (5)--(6) and the ring-level scheme statement is
outstanding. Zero/negative weights, variable-coefficient operators and formal
nonpolynomial multipliers are outside the proof. General Gaussian moment
uniformity, unbalanced higher-degree GVC families, full mixed ordinary ideals
and optimal numerical frontiers remain open.
