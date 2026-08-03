# Affine support versus Newton and boundary geometry

## Status and purpose

The support-six theorem proves that a noninvertible plane Keller map \(F\)
must satisfy

\[
\sigma_{\mathrm{aff}}(F)\geq7.
\]

The tempting next step is to force a minimal standard counterexample to have
\(\sigma_{\mathrm{aff}}\leq6\) from its Newton polygon or finite-normalization
ledger.  This note tests that bridge and records two exact results:

1. geometric degree, reduced nonproperness data, boundary-row counts, and
   Newton vertex count cannot upper-bound \(\sigma_{\mathrm{aff}}\), even for
   Keller automorphisms;
2. a sparse monomial-Jacobian Newton block reaches the support theorem only
   when it passes an exact Kummer-character descent gate.

The live \((75,125)\) F2 terminal block fails that gate.  Thus the direct
bridge does not close JC(2), but its failure identifies the missing datum:
the full lower-band character packet, not another support-cardinality layer.

The minimal-pair input is the lower-edge framework of
[Guccione--Guccione--Valqui](https://arxiv.org/abs/1605.09430), where a
standard minimal counterexample produces weighted-homogeneous edge data and
Poisson equations such as \([G,R]=R^2\).  No external theorem is used in the
two propositions below.

## 1. Coarse geometry cannot control affine support

For a plane map \(F\), recall

\[
N_{p,L}(z)
 =(JF(p)L)^{-1}\bigl(F(p+Lz)-F(p)\bigr),
\qquad
\sigma_{\mathrm{aff}}(F)=\min_{p,L}\sigma(N_{p,L}).
\]

### Proposition 1.1 -- unbounded affine support with fixed coarse geometry

For every \(d\geq4\), there is a Zariski-open set of polynomials

\[
p(t)=a_2t^2+\cdots+a_dt^d,\qquad a_2\cdots a_d\ne0,
\]

such that the triangular Keller automorphism

\[
F_p(x,y)=(x+p(y),y)
\]

satisfies

\[
\sigma_{\mathrm{aff}}(F_p)\geq d-2. \tag{1}
\]

Meanwhile:

- \(F_p\) has geometric degree one;
- its nonproperness set and missing-boundary ledger are empty;
- the Newton polygon of \(x+p(y)\) has the same three vertices
  \((1,0),(0,2),(0,d)\), independently of the \(d-3\) interior coefficients.

#### Proof

Choose \(p\) so that the polynomials

\[
p^{(2)},p^{(3)},\ldots,p^{(d-1)}
\]

are pairwise root-disjoint.  Such \(p\) form a nonempty Zariski-open set.
Indeed, for fixed \(2\leq i<j\leq d-1\), the incidence conditions

\[
p^{(i)}(c)=p^{(j)}(c)=0
\]

are two independent linear equations on \((a_2,\ldots,a_d)\) for each
fixed \(c\): the coefficient \(a_i\) occurs in the first functional and not
the second.  The incidence variety therefore has dimension at most \(d-2\),
strictly below the \(d-1\)-dimensional coefficient space.  Removing the
finitely many projected incidence closures and the coordinate hyperplanes
leaves the required open set.

Fix a normalization center \((r,c)\) and \(L\in\operatorname{GL}_2(k)\).
If \(\ell(z)\) is the nonzero linear form given by the second row of \(L\),
then the nonlinear part of \(N_{(r,c),L}\) is

\[
v\left(
 p(c+\ell(z))-p(c)-p'(c)\ell(z)
\right), \tag{2}
\]

for a nonzero vector \(v\).  The coefficient of total degree \(j\) in (2)
is \(v\,p^{(j)}(c)\ell(z)^j/j!\).  At a fixed \(c\), pairwise
root-disjointness allows at most one of
\(p^{(2)}(c),\ldots,p^{(d-1)}(c)\) to vanish, while \(p^{(d)}\) never
vanishes.  Hence at least \(d-2\) distinct homogeneous degrees occur.
Each contributes at least one coordinate-monomial occurrence, proving (1).
\(\square\)

This rules out any bridge depending only on geometric degree, the reduced
nonproperness curve, missing-boundary row counts, or the number of Newton
vertices.  A counterexample-specific bridge must use coefficient or
lower-band information absent from those coarse invariants.

The checker verifies the pairwise derivative resultants for the explicit
regression family

\[
p_d(t)=\sum_{j=2}^{d}\frac{t^j}{j!},
\qquad 4\leq d\leq12.
\]

The finite regression is not the proof of the generic proposition.

## 2. The Kummer-character descent gate

Sparse terminal Newton blocks generally have monomial rather than constant
Jacobian.  This changes the problem.

### Proposition 2.1 -- Kummer descent

Let \(r\geq0\), put

\[
u=\frac{x^{r+1}}{r+1},
\]

and suppose \(P,Q\in k[u,y]\subset k[x,y]\).  Then

\[
[P,Q]_{x,y}=x^r[P,Q]_{u,y}. \tag{3}
\]

Consequently,

\[
[P,Q]_{x,y}=c x^r
\quad\Longleftrightarrow\quad
[P,Q]_{u,y}=c. \tag{4}
\]

If the descended constant-Jacobian map has affine-normalized nonlinear
support at most six, the certified support-six theorem makes it a polynomial
automorphism of \(k[u,y]\).  The original monomial-Jacobian block is then a
Kummer pullback of that automorphism; it contains no additional hidden
constant-Jacobian component.

The membership condition \(P,Q\in k[u,y]\) is equivalent to every
\(x\)-exponent in both supports being divisible by \(r+1\).  Equivalently,
both components lie in the trivial character of the
\(\mu_{r+1}\)-action \(x\mapsto\zeta x\).

Equation (3) is the chain rule, so this gate is exact in every
characteristic zero field.

### Proposition 2.2 -- character-resolved bracket equations

Decompose a monomial-Jacobian block into \(\mu_{r+1}\)-characters,

\[
P=\sum_{a\in\mathbb Z/(r+1)}P_a,\qquad
Q=\sum_{b\in\mathbb Z/(r+1)}Q_b.
\]

Differentiation with respect to \(x\) lowers the character by one, while
differentiation with respect to \(y\) preserves it.  Hence

\[
[P_a,Q_b]\quad\text{has character}\quad a+b-1. \tag{5}
\]

Since \(x^r\) has character \(r=-1\), the equation
\([P,Q]=c x^r\) splits into \(r+1\) independent equations:

\[
\sum_{a+b-1=\chi}[P_a,Q_b]
=
\begin{cases}
c x^r,&\chi=-1,\\
0,&\chi\ne-1.
\end{cases}
\tag{6}
\]

Thus only complementary pairs \(a+b=0\) contribute to the target character;
all other character sectors must cancel internally.  This is the
character-resolved system that missing Laurent bands must satisfy.

## 3. Application to the live F2 terminal block

The audited \((75,125)\) F2 derivation uses the integral Kummer coordinate
\(X=x^{1/5}\).  Its normalized terminal type-I block is

\[
\begin{aligned}
P_{\mathrm I}
 &=X^4y\left(1+X^{17}y^5\right),\\
Q_{\mathrm I}
 &=-X\left(1+3X^{17}y^5+\frac95X^{34}y^{10}\right),
\end{aligned}
\tag{7}
\]

and direct differentiation gives

\[
[P_{\mathrm I},Q_{\mathrm I}]_{X,y}=X^4. \tag{8}
\]

Thus (7) has only five monomial occurrences, but it is not a five-support
Keller map.  For the Kummer descent \(u=X^5/5\), its \(X\)-characters are

\[
\begin{array}{c|c|c}
&X\text{-exponents}&\text{residues modulo }5\\ \hline
P_{\mathrm I}&4,21&4,1\\
Q_{\mathrm I}&1,18,35&1,3,0.
\end{array}
\]

Neither component lies in \(k[u,y]\).  The terminal block therefore fails
the exact hypothesis of Proposition 2.1.  Applying the constant-Jacobian
support-six theorem to (7) would be invalid.

This is not merely a coordinate nuisance.  The different Kummer characters
encode which bracket pairs can contribute to the \(X^4\) right-hand side.
The exact sector calculation gives

\[
\begin{array}{c|c}
\text{character}&\text{sum of terminal brackets}\\ \hline
0&0\\
1&6X^{21}y^5-6X^{21}y^5=0\\
3&9X^{38}y^{10}-9X^{38}y^{10}=0\\
4&X^4.
\end{array}
\]

They must be followed through the missing lower bands because every new band
adds terms to one of these independent equations.

## 4. Bridge outcome

The attempted implication

\[
\text{minimal Newton/boundary data}
\quad\Longrightarrow\quad
\sigma_{\mathrm{aff}}\leq6
\]

does not follow from the current coarse invariants:

- Proposition 1.1 rules out a general bound from Newton vertices,
  geometric degree, or coarse boundary ledgers.
- Proposition 2.1 shows that sparse monomial-Jacobian edge blocks reach the
  support theorem only through a Kummer-character descent.
- The first live \((75,125)\) block fails that descent exactly.

The viable replacement is a **character-resolved lower-band bridge**:

> For every minimal standard-pair branch, compile each Laurent band by its
> Kummer character, prove the band list exhaustive, and show that either the
> packet descends to a support-at-most-six Keller map or its nontrivial
> characters force an impossible log-boundary/ramification ledger.

For F2, the corrected chart has \([t,z]_{X,y}=-z\).  The common-power band is
layer 40, the \(X^4\) target is layer 4, and the 35 intervening zero layers
are exactly 39 through 5.  The finite B0 degree/halfspace envelope and all
five character sectors of those layers are now compiled exactly in
[`F2_75_125_DERIVATION.md`](F2_75_125_DERIVATION.md).  The same replay now
includes every lower B0 band through bracket layer `-200`: `2,418`
jet-reduced parameters and `240` zero layers.  Inconsistency of this
over-envelope would exclude the row directly; a surviving point would still
need an exhaustive B1 polygon mask before it represented a plane candidate.
The corrected top tangent is now exact and shows why this route does not
triangularize so early.  The former substitution by (C_0^2) and (C_0^4)
assumed unproved divisibility.  Every exact P-band direction instead has the
Q follower `q=-3*C0^2*p`; the first-five kernel dimensions are
`6,6,7,7,10`.  The extra layer-35 term is the commuting (C_0^4), while the
formal `lambda*C0^(-1)` resonance is at layer 10 and is not an independent
source-band mode.  The next optional coefficient test must retain all P-band
variables and impose the nonlinear forcing by the exact cokernel/Fitting row.
The first such calculation is now complete: a resultant `1701*a^8` recovers
a source root through descent 7 and leaves the local descent-8 ratio
`27*y^2-9*y+1=0`.  The Q-band-one normalization excludes its four fixed
Kummer supports.  At this earliest spacing only a nonzero double root of `R` remains; it passes the
first local target jet, and its fifth multiple is a layer-zero lower-tail
Fitting condition rather than the primitive equation `E5=0`.  The full
target cokernel and the layer-zero Artinian quotient both have rank `14`.
Old B0 generators span them before the earlier equations are imposed, and
the forced lowest-`u` edge correction has an exact Bezout witness.  That
witness extends to an all-`r` formal shear, but the shear is necessarily
infinite: its literal polynomial truncation fails at transverse order two.
The polynomial order-two repair is an explicit two-parameter family but
never terminates quadratically; exact `r=3` calculations also exclude cubic
and quartic termination.  The exact first Kummer-return packet at `v^5` has
since been compiled in
[`F2_75_125_DERIVATION.md`](F2_75_125_DERIVATION.md): its maximal-minor ideal
is the unit ideal, both descent-eight branches survive, and the isolated
edge recursion remains surjective through `v^10`.  The live obstruction is
therefore in the global Hermite quotient, not in another transverse
termination test.  This is now a support theorem rather than a coefficient-
box observation: the all-`r` unit minors select `18*r-1` original-polynomial
source combinations strictly below both certified edges (`53` at `r=3`).
Quotienting their triangular `w=0` control by exact confluent CRT leaves a
rank-`24` global Hermite module over the rank-two descent-eight algebra.  The
fixed `w=1` summand is now one normalized identity plus a ten-variable block
of determinant `75000`; eliminating it leaves thirteen cokernel coordinates.
The substitution has now been carried, and the endpoint-disjoint power block
is eliminated through layer `29`.  It leaves `927` active source coordinates,
with the coupled Schur/Fitting system beginning at layer `28`.  Only that
coupled image or a genuinely new lower Newton edge can close this branch.
Later first-defect spacings remain in the full system.

The first boundary pivot is complete at the contact-only level justified by
the common edge.  Its quadratic factor has four exhaustive contact
partitions, but those contacts alone do not determine normal scales or
finite-normalization rows.  Even the strongest naive contact-to-row surrogate
survives the finite-flat packet budget.  The retained audit is
[`F2_BOUNDARY_HANDOFF.md`](F2_BOUNDARY_HANDOFF.md).  Subsequent exact work in
[`F2_KUMMER_ORBIT_TRANSFER.md`](F2_KUMMER_ORBIT_TRANSFER.md) and
[`F2_TERMINAL_RESIDUE_COVER.md`](F2_TERMINAL_RESIDUE_COVER.md) bypasses that
surrogate: it reduces the live row to one principal chain or two copies and
  computes its degree-six target residue cover.  The 35 nonlinear layers still
should not be continued sequentially; the remaining problem is global
source/target gluing.  Its first processed consequences are geometric degree
at least six, or at least twelve for two distinct packets over one target
divisor; no affine-sheet increment applies because the certified target
valuation is centered at infinity.  Purity instead requires a separate
ramified row over an affine nonproperness curve.  The target-node fibers also
fix three interior source-boundary attachment points, so the global route now
starts from an incidence skeleton rather than an undecorated `(1,6)` row.
<!-- status-consumer: PF2GC1 33dbc5ff48b5d064 -->

## 5. Reproduction and claim boundary

Run:

```bash
.venv/bin/python plane-jc/cas/verify_affine_support_newton_bridge.py
.venv/bin/python plane-jc/cas/classify_f2_75_125_layers.py
```

Intentional artifact regeneration uses `--refresh`.  The pinned artifact is
[`artifacts/generated-results/jc2_affine_support_newton_bridge.json`](../artifacts/generated-results/jc2_affine_support_newton_bridge.json).

What is proved:

- geometric degree, reduced nonproperness data, boundary-row counts, and
  Newton vertex count do not upper-bound \(\sigma_{\mathrm{aff}}\);
- the Kummer descent criterion (3)--(4);
- the character-resolved bracket equations and the exact nontrivial
  character profile of the F2 terminal block.

What remains open:

- a support bound for minimal counterexamples;
- the exhaustive F2 lower-band classification;
- any exclusion of \((75,125)\) or improvement of the degree frontier.
