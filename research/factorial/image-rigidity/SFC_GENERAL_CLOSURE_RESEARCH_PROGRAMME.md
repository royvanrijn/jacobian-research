# A research programme for complete Strong Factorial closure

**Status: 21 November 2026.** This is a ranked programme with precise proof
targets, verified entry points, and falsification tests. It is not a proof of
unrestricted SFC. The exact caustic/composition certificate and scalar-operator
obstructions are proved in the linked notes. The resultant positivity evidence
is a finite collection of exact symbolic identities. All universal closing
targets below remain conjectural.

The inspected GitHub main revision is
`19eb9912fb60142f8b620dac24896931d4ff3e46`. The current internal baseline is
[SFC_STRUCTURAL_CLOSURE.md](SFC_STRUCTURAL_CLOSURE.md), including its latest
algebraic-window, initial-resultant and pullback proofs. No external review,
formal verification, or literature-priority claim is made.

## 1. Two separate closure problems

For a fixed ordinary support \(a_1,\ldots,a_s\in\mathbb N^d\), write

\[
 M_m(c)=\Lambda\!\left(\left(\sum_{j=1}^s c_jX^{a_j}\right)^m\right)
 =\sum_{|k|=m}\binom{m}{k}\prod_{i=1}^d(a_i^{\rm row}\cdot k)!\,c^k.
\]

The target is that \(M_n,\ldots,M_{n+s-1}\) have no common nonzero coefficient
vector, for every \(n\ge1\). Coefficients may lie in any characteristic-zero
field; passing to an algebraic closure supplies no positivity assumption.

There are two milestones:

1. **Sharp monomial-linear closure:** for
   \(F=(\prod_iX_i^{r_i})\sum_i c_iX_i\), all nonnegative integer profiles,
   all starting orders, and all coefficient vectors.
2. **Unrestricted ordinary-support closure:** arbitrary exponent matrices.
   This requires a separate transfer theorem or a direct proof.

The fixed-window theorem handles all starting orders in the first class,
but with a window larger than \(s\). The initial support-size window is closed,
and some infinite structural classes have the sharp window at all orders.
Neither statement settles milestone 1; milestone 1 does not imply milestone 2.

| Programme | Closing target | Scope if proved | Assessment |
|---|---|---|---|
| Shifted resultant combinatorics | Positive certificates for every normalized shifted resultant | All monomial-linear profiles, including zero exponents | Best defined global algebraic target |
| Pointed pencils and composition | A forbidden sharp gap forces a common right factor | Entire polynomial-\(H\) monomial-linear class | Strongest new geometric representation |
| Factorial transfer on finite ideals | Preserve the ordinary Macaulay map's surjectivity under support moves | Arbitrary ordinary supports, with an explicit seed construction | Largest reach; substantial missing bridge |
| Adversarial audit | Challenge period-span, positivity, genericity and limiting arguments | Protects all three programmes | Run beside the proof searches |

The first two should share their combinatorics: inverse-tree expansions supply
identities; resultants or finite jet ideals turn identities into nonvanishing
certificates. Expanding more moments alone does not resolve complex cancellation.

## 2. Shifted resultants: positivity in the exponent parameters

Put

\[
 g_m(c;r)=[z^m]\prod_{i=1}^s(1-c_i z)^{-mr_i-1},\qquad
 {\cal R}_{s,n}(r)=\operatorname{Res}_{c}(g_n,\ldots,g_{n+s-1}).
\]

Each \(g_m\) is homogeneous of degree \(m\) in \(c\), with coefficients
polynomial over \(\mathbb Q\) in the formal \(r_i\). Nonvanishing of this one
resultant excludes every nonzero complex coefficient vector simultaneously.

Normalize \({\cal R}_{s,n}(0)\) to \(1\). At \(r=0\), the \(g_m\) are the
complete homogeneous functions \(h_m(c)\). A block of \(s\) zeros in
\(\sum h_m t^m=1/\prod_i(1-c_i t)\) propagates by its denominator recurrence.
A polynomial reciprocal then forces every \(c_i=0\). This works over every
field, so the integral resultant at zero is a unit, \(+1\) or \(-1\).

The exact initial theorem is

\[
 {\cal R}_{s,1}(r)=
 \prod_{\varnothing\ne I\subseteq[s]}
 \left(1+\sum_{i\in I}r_i\right)^{(|I|-1)!(s-|I|)!}.
\]

See [the initial-resultant proof](SFC_INITIAL_RESULTANT_AND_SHIFTED_WINDOWS.md).

### Sufficient all-order conjecture

The weakest proposed positivity target uses the product Newton basis:

\[
 {\cal R}_{s,n}(r)=\sum_{\alpha\in\mathbb N^s}
 A_{s,n,\alpha}\prod_i\binom{r_i}{\alpha_i},\qquad
 A_{s,n,\alpha}=\Delta^\alpha{\cal R}_{s,n}(0)\ge0.
\]

The sum is finite, its constant coefficient is one, and all binomial factors
are nonnegative at the required nonnegative integer exponents. This would
close monomial-linear SFC. The coefficients are integers: at integer exponent
tuples the moment forms, their resultant, and hence its finite differences
are integral. This makes an interpretation as counts of decorated forests
particularly natural, although no such interpretation is yet proved.

Ordinary monomial-coefficient positivity is a stronger target and would also
give nonvanishing for all nonnegative real exponent parameters. Schur positivity
is stronger still. The eliminated coefficients \(c_i\) remain arbitrary in
all three formulations. A negative monomial or Schur coefficient need not
refute the weaker Newton-basis target; a negative Newton coefficient would.
None of those failures alone would refute SFC.

The [positivity audit](SFC_SHIFTED_RESULTANT_POSITIVITY_AUDIT.md) checks binary
starting orders \(1,\ldots,5\), symbolically in both exponent parameters.
Every normalized resultant has nonnegative coefficients; after the displayed
elementary factors are removed, the tested residuals are also Schur-positive.
At starting order \(5\), the full resultant has degree \(45\) and \(841\)
nonzero positive coefficients; the residual has \(295\) positive Schur terms.

These are parameter identities, not integer samples. They still cover only
five starting orders and support size two, not an all-order or higher-support
theorem. An important negative control is the order-two residual at one zero
exponent, proportional to

\[
 17r^3+37r^2+22r+4,
\]

with discriminant \(-428\). The hoped-for positivity is not explained by
universal real-rootedness or a positive-semidefinite affine determinant pencil.

### Where trees should enter

Use the Koszul/Macaulay determinant presentation of the resultant, insert an
exact inverse-tree expansion, and seek a sign-reversing involution leaving
positively weighted forests or tableaux. Account for determinant signs, not
only the positive moment coefficients.

The first proof target is an identity for arbitrary starting order, even in
the binary symbolic family. In parallel, a fully symbolic ternary window is
a stronger falsification test than additional binary samples. A successful
formula must recover both the initial subset product and shifted nonlinear
residuals.

Varchenko's arrangement forms are a useful comparison: the initial factors
resemble edge weights. But the standard full contravariant determinant is a
product of linear weight factors, while the shifted residual is already
nonlinear. No direct identification is established. A subquotient, minor, or
different weight map would have to be constructed.

## 3. Pointed pencils: retain the branch and detect composition

In the polynomial-\(H\) class, put

\[
 a=z\prod_i(1-c_i z)^{r_i},\quad
 H(0)=0,\quad H'=\prod_i(1-c_i z)^{r_i-1},\quad D=a'/H'.
\]

Repeated coefficient groups are permitted when each group's total exponent
is at least its variable multiplicity. Then \(S=\deg D=\deg a-\deg H\) is
the actual support size. The sharp target is

\[
 \operatorname{ord}_0\{H-b(a)\}\le\deg b+S.
\]

The pencil

\[
 a_\varepsilon=a+\varepsilon H,\qquad
 a_\varepsilon'=H'(D+\varepsilon)
\]

has \(S\) moving critical points, counted with multiplicity. Its other critical
points are fixed, with
critical-value lines \(t=\varepsilon H(\alpha)\) through \((0,0)\), because
\(H'(\alpha)=0\) implies \(a(\alpha)=0\). The marked source point \(z=0\)
lies over \((\varepsilon,t)=(-1,0)\) on the moving curve, not its concurrence
point. It must remain marked on the normalization.

### A proved exact certificate

Set \(T=a-DH\) and

\[
 {\cal C}(t,\varepsilon)
 =\operatorname{Res}_z(a+\varepsilon H-t,D+\varepsilon).
\]

Over any characteristic-zero field,

\[
 K(D,T)=K(a,H),\qquad
 {\cal C}=c\,{\cal F}^{\,e},\qquad
 e=[K(z):K(a,H)],\qquad \deg_t{\cal F}=S/e,
\]

where \({\cal F}\) is monic in \(t\) and absolutely irreducible and \(c\ne0\).
The identities \(T'=-D'H\) and \(D=da/dH\) prove the field equality; a field
norm gives the exact multiplicity.

This multiplicity is the maximal common polynomial right-factor degree.
The existing [descent theorem](SFC_RAMIFIED_PULLBACK_CLOSURE.md) then gives
another ordinary factorial frame of support \(S/e\), preserving every
normalized moment. It covers all common right factors, not only cyclic ones.
The new certificate recognizes exactly when that established reduction applies.

This is an application of plane-curve duality, not a novelty claim for duality.
Full proof, base-field details and examples are in
[the pencil note](SFC_POINTED_PENCIL_AND_CAUSTIC_DESCENT.md).

### Proposed closing lemma

> For an actual root-aligned frame, a block of \(S\) zero moments forces
> \({\cal C}\) to be a proper power.

Equivalently, prove the sharp contact bound for primitive pairs
\(K(a,H)=K(z)\). A smallest-support counterexample would be primitive; the
proposed implication would force smaller-support descent, a contradiction.

This remains an unproved form of the sharp problem. The useful additional
structure is a moving curve retaining the pair, the distinguished branch,
and the exact descent obstruction. One concrete attack is to derive from
the vanishing jet a subresultant relation forcing
\(\gcd({\cal C},\partial_t{\cal C})\ne1\).

### A precise bridge to Jacobian inversion and trees

For \(Y(t,\varepsilon)=H(a_\varepsilon^{-1}(t))\),

\[
 Y_\varepsilon+Y Y_t=0,\qquad
 Y_t=(D(a_\varepsilon^{-1}(t))+\varepsilon)^{-1}.
\]

With \(B=H\circ a^{-1}\),

\[
 Y=\sum_{k\ge0}\frac{(-\varepsilon)^k}{(k+1)!}
       \partial_t^k(B^{k+1}).
\]

The explicit change of variables in the pencil note converts this into
Zhao's inverse/Burgers deformation. His Theorems 2.4 and 5.2 supply inverse
identities and a binary rooted-tree expansion over a \(\mathbb Q\)-algebra.
Our series is usually algebraic, not polynomial: no nilpotence or termination
conclusion follows.

The challenge is a uniform contact-to-composition identity using root alignment
and the marked branch. Burgers and moving degree alone are insufficient: the
note gives a degree-two moving curve without alignment whose contact exceeds
\(\deg b+S\).

Singularity is also insufficient. The primitive degrees-\((5,2)\) control has
a cusp at its marked caustic point, but the caustic is irreducible. The proposed
conclusion requires a repeated global factor, not merely vanishing partial
derivatives at one point.

Another bridge is the zero-cycle center theorem of Álvarez--Bravo--Christopher--
Mardešić: identically zero full displacement implies composition. It could
apply if a forbidden finite gap forced a nonzero integral cycle with persistent
vanishing. Constructing that cycle and proving the full two-variable identity
are both missing. Tangential vanishing is insufficient.

The profile \((1,1,2)\), for arbitrary coefficients and starting order, is a
useful stress test: degrees \(5\) and \(2\) rule out composition. Use it to test
a structural mechanism, not to launch another fixed-support catalogue.

## 4. Full-support transfer must act on ordinary moment ideals

For arbitrary support let

\[
 d_j=n+j\ (0\le j<s),\quad
 N=1+\sum_{j=0}^{s-1}(d_j-1),\quad R=K[c_1,\ldots,c_s].
\]

The sharp window is equivalent to surjectivity of the finite map

\[
 \Phi:\bigoplus_{j=0}^{s-1}R_{N-d_j}\longrightarrow R_N,\qquad
 (q_j)\longmapsto\sum_jq_jM_{d_j}.
\]

If the forms have only the origin as a common zero, they are a homogeneous
regular sequence and their complete-intersection Hilbert series vanishes in
degree \(N\). Conversely surjectivity puts every degree-\(N\) monomial in the
ideal and excludes a projective zero. This uses ordinary multiplication.

The support moves have exact coefficient multipliers:

- Adjoining a fresh unit row in column \(j\): multiply \(c^k\) by \(k_j!\).
- Merging rows \(u,v\): multiply \(c^k\) by
  \(\binom{(u+v)\cdot k}{u\cdot k}\).

A uniform theorem transporting surjectivity along such factorial moves would
be a route to arbitrary supports, with a specified seed construction and
collision-free intermediate supports. Removing private factors from the
proved amplification endpoint is another sufficient transfer target.

The useful certificate is a finite chain map or explicit factorization between
the two Macaulay maps, controlling every denominator and exceptional parameter.
An invertible coefficient multiplier does not transport ordinary ideals.
The [de-amplification](SFC_DEAMPLIFICATION_BRIDGES.md) and
[row-merger](SFC_ROW_MERGING_OBSTRUCTIONS.md) audits establish that boundary.

### Why general hypergeometric rank theory is insufficient

The natural formal exponential integral for exponent matrix \(A\) uses
\([I_d\mid A]\). Its semigroup is already \(\mathbb N^d\), and its semigroup
algebra a polynomial ring, for every ordinary support. Cohen--Macaulayness
and the corresponding absence of GKZ rank jumps do not distinguish the
unresolved factorial examples.

Euler--Koszul contiguity still supplies precise module maps and exceptional
parameter sets. Positive Gamma-shape increments have appropriate one-step
isomorphisms, with GKZ Euler parameter the negative of the Gamma shape.
This does not say that all actual parameters avoid the full strongly resonant
locus.

The missing step is preservation of the selected factorial functional and
the finite ordinary map \(\Phi\) after restriction. The
[positivity/GKZ audit](SFC_SHIFTED_RESULTANT_POSITIVITY_AUDIT.md) states the
sign convention and an exact ideal-level countercontrol. A \(D\)-module
isomorphism cannot replace this commutative statement.

## 5. Exact obstructions eliminate tempting shortcuts

A rational homogeneous differential operator annihilating the actual algebraic
series also annihilates all conjugates and their constant linear combinations.
Some combinations have forbidden short zero blocks when the actual branch does
not.

The [differential audit](SFC_DIFFERENTIAL_BRANCH_OBSTRUCTIONS.md) proves this
inside the ordinary profile \((1,2)\), at two exact quadratic coefficient
parameters. Every support-width operator must acquire a positive-integer
resonance. At the second parameter it must also annihilate a particular
nonconstant linear polynomial.

Thus neither universal operator premise survives:

- Support-size width with no late positive-integer resonance.
- Support-size width with only constants as polynomial solutions.

For polynomial \(H\), trace zero also rules out a nonzero rational inhomogeneous
source as a branch selector. Higher order and apparent desingularization remain
computational tools, but cannot remove these structural obstructions.

Do not infer sharp windows from period-module dimension, arbitrary positive
Gram nondegeneracy, limits of nonzero determinants, or generic regularity.
The retained controls show why the information lost at those steps matters.

## 6. Concrete next deliverables

1. **Resultant team:** an all-\(n\) binary Schur/forest identity; in parallel,
   challenge positivity on a full symbolic ternary window. Return an identity
   or a negative coefficient, not just numerical samples.
2. **Geometry/tree team:** derive finite-jet-to-composition while retaining
   the aligned pencil and marked normalization point. Test primitive,
   genuine-descent and misaligned controls. Degree in \(t\) alone is ruled out.
3. **Full-support team:** construct a map between ordinary Macaulay
   presentations for a support move. Test multiplication compatibility;
   generic GKZ rank and invertible Hadamard scaling are insufficient.
4. **Audit team:** exact witnesses; coefficient collisions and coordinate
   faces; active attempts to falsify every strengthened intermediate claim.

The main recommendation is programme 1, with programme 2 supplying the strongest
structure for tree identities and composition. Keep the separate full-support
transfer explicit: monomial-linear success must not be called unrestricted SFC.

## 7. Primary research inspected

- Zhao, *Inversion Problem, Legendre Transform and Inviscid Burgers'
  Equations*, Theorems 2.4 and 5.2:
  <https://arxiv.org/pdf/math/0403020>.
- Álvarez--Bravo--Christopher--Mardešić, *Infinitesimal Center Problem on Zero
  Cycles and the Composition Conjecture*, Theorem A:
  <https://arxiv.org/pdf/2006.07600>.
- Bravo--Mardešić--Novikov--Pontigo-Herrera, *Infinitesimal and tangential
  16-th Hilbert problem on zero-cycles*, Bull. Sci. Math. 202 (2025), 103634.
  Theorem D refutes a broad Chebyshev expectation even for simple cycles:
  <https://arxiv.org/pdf/2312.03081>.
- Schulze--Walther, *Hypergeometric D-modules and twisted Gauss--Manin systems*:
  <https://arxiv.org/pdf/0712.2021>. One-step contiguity is Remark 3.6 in
  arXiv v3; numbering differs in the author-hosted version.
- Matusevich--Miller--Walther, *Homological Methods for Hypergeometric
  Families*, Corollary 9.2: <https://arxiv.org/pdf/math/0406383>.
- Varchenko, *Quantum Integrable Model of an Arrangement of Hyperplanes*:
  <https://sigma-journal.com/2011/032/sigma11-032.pdf>.

The audit notes state hypotheses and limitations. No source is represented
as already proving a sharp factorial window.
