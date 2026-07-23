# Reproduction of the \((72,108)\) exclusion

> **Classification of content.**  Sections marked **external reduction** restate
> Guccione--Guccione--Horruitiner--Valqui (GGHV) or Helali.  Sections marked
> **local reproduction** report exact reruns and line-by-line algebra checks.
> **Repository interpretation** is explanatory, not part of the external
> proof.  The final section lists open conceptual work and makes no theorem
> claim.

See [PROVENANCE.md](PROVENANCE.md) for the exact PDF, archive, source commit,
hashes, licenses, and software versions.

## 1. Published Newton reduction — external reduction

GGHV Proposition 4.3 begins with the last admissible chain

\[
 A_0=(8,28),\qquad A_1=(11/4,7),\qquad (m,n)=(3,2).
\]

Its original standard-pair degree orientation is

\[
 (\deg P,\deg Q)=(3(8+28),2(8+28))=(108,72).
\]

Coordinate reversal gives \((72,108)\).  After the paper's sequence of a
coordinate flip, Laurent triangular translations, edge reduction, and the
monomial morphism \(x\mapsto x^{-1},y\mapsto x^4y\), a counterexample produces
\(P,Q\in L^{(1)}=k[x,x^{-1},y]\), \([P,Q]=x^2\), in one of:

\[
\begin{aligned}
\text{Case 1: }N(P)&=\operatorname{conv}\{(0,0),(1,0),(8,14),(8,16),(0,8)\},\\
N(Q)&=\operatorname{conv}\{(0,0),(2,1),(12,21),(12,24),(0,12)\};\\[2mm]
\text{Case 2: }N(P)&=\operatorname{conv}\{(0,0),(1,0),(8,14),(8,16)\},\\
N(Q)&=\operatorname{conv}\{(0,0),(2,1),(12,21),(12,24)\}.
\end{aligned}
\]

The proposition's proof is not merely a choice of orientation.  It uses
GGV's earlier Corollary 7.4, Proposition 3.12, Proposition 2.5, Proposition
8.2, and the classification of predecessor/successor directions.  In the
intermediate coordinates it obtains the possible lower corners
\((-2,0)\), \((-3,0)\), and, in the split-edge case, \((16,4)\); the reduced
edge has \((24,7)\), and the opposite candidates are
\((24,7),(17,5),(10,3),(3,1)\).  Proposition 8.2 leaves \(k\in\{1,2\}\);
parallelism excludes \(k=2\), and \(k=1\) supplies endpoints \((-1,0)\) and
\((2,1)\).  These are exactly the vertices visible after the final monomial
map.

## 2. Coordinate and support transcription — local reproduction

Set

\[
 t=xy^2,\qquad z=y^{-1}.
\]

Then

\[
 [t,z]_{x,y}=-1,\qquad x^2=t^2z^4,
\]

and a lattice monomial transforms as

\[
 x^iy^j=t^iz^{2i-j}.
\]

The local support enumerator checks every lattice point, not only the
vertices:

| Case | \(\#(N(P)\cap\mathbb Z^2)\) | \(P\)-bands | \(\#(N(Q)\cap\mathbb Z^2)\) | \(Q\)-bands |
| --- | ---: | --- | ---: | --- |
| 1 | 61 | \(z^2,z^1,z^0,\ldots,z^{-8}\) | 125 | \(z^3,z^2,z^1,\ldots,z^{-12}\) |
| 2 | 25 | \(z^2,z^1,z^0\) | 47 | \(z^3,z^2,z^1,z^0\) |

This transcription is executable in
[`cas/laurent_band_frontend.py`](cas/laurent_band_frontend.py).  Its typed
normal-form certificate records that Proposition 4.3 is exhaustive and maps
the single printed chain to these two distinct Laurent cases; the compiler
then enumerates the supports and derives the bracket layers below.

For the common upper part, after subtracting irrelevant constants,

\[
 P=A(t)z^2+B(t)z+C(t),\qquad
 Q=D(t)z^3+E(t)z^2+F(t)z+G(t).
\]

The generic band formula

\[
 [P_i(t)z^i,Q_j(t)z^j]_{x,y}
 =z^{i+j-1}\bigl(iP_iQ_j'-jP_i'Q_j\bigr)
\]

gives, coefficient by coefficient,

\[
\begin{aligned}
2AD'-3A'D&=t^2 &(J4)\\
2(AE'-A'E)+(BD'-3B'D)&=0 &(J3)\\
(2AF'-A'F)+(BE'-2B'E)-3C'D&=0 &(J2)\\
2AG'+(BF'-B'F)-2C'E&=0 &(J1)\\
BG'-C'F&=0. &(J0)
\end{aligned}
\]

These signs were derived independently from the chain rule.  Case 2 uses all
four bands and therefore the full coefficient system.  Case 1 uses a
necessary upper truncation: the descent adds \(P\)-bands through \(z^{-5}\)
and \(Q\)-bands through \(z^{-4}\), producing all compatibility equations at
bracket layers \(z^1,z^0,z^{-1},z^{-2},z^{-3}\).  Lower bands cannot repair a
compatibility failure at those already determined layers.  The archive's
phrase “complete coefficient system” is accurate for Case 2 but should be
read as “complete necessary truncated system” for Case 1.

## 3. Scaling normalization — local reproduction

Let \(a_1,a_8,c_8\ne0\) be the coefficients of the three nonzero vertices
\((1,0),(8,14),(8,16)\) of \(P\).  Their nonvanishing is part of the assertion
that these points are vertices of the exact Newton polygon.  For

\[
 \widetilde P=\rho P(\lambda x,\mu y),\qquad
 \widetilde Q=\sigma Q(\lambda x,\mu y),
\]

choose in the algebraically closed ground field

\[
 \mu^2={a_8\over c_8},\quad
 \lambda^7={a_1\over a_8}\mu^{-14},\quad
 \rho=(a_1\lambda)^{-1},\quad
 \sigma=(\rho\lambda^3\mu)^{-1}.
\]

Direct substitution gives the three normalized coefficients equal to one and
keeps \([\widetilde P,\widetilde Q]=x^2\).  No zero coefficient is divided by.
Thus

\[
 A=t+a_2t^2+\cdots+a_7t^7+t^8,qquad
 D=d_2t^2+\cdots+d_{12}t^{12}.
\]

## 4. First block and exact coefficient field — local reproduction

Triangular coefficient solving in \((J4)\) determines all eleven
\(d_2,\ldots,d_{12}\) and leaves six equations in
\(a_2,\ldots,a_7\).  All pivots here are nonzero rational integers.  A fresh
Singular 4.4.1 characteristic-zero run of the archived input
`firstblock_Q_exact.sing` returned:

```text
DP_SIZE=56
LEX_SIZE=6
```

The reduced lexicographic ideal has one irreducible polynomial \(H(a_7)\) of
degree 35 and five relations linear in \(a_2,\ldots,a_6\).  Hence every
solution defines an embedding of

\[
 K_0=\mathbb Q[u]/(H(u)),\qquad [K_0:\mathbb Q]=35.
\]

The run performs `modStd` followed by exact rational reconstruction and FGLM;
the local symbolic verifier then substitutes the five linear relations back
into all six original equations and checks zero modulo \(H\).  Irreducibility
of \(H\) over \(\mathbb Q\) was checked exactly.

## 5. Case 2 ideal — external certificate, locally reproduced

Exact linear solution of \((J3)\) and \((J2)\) over \(K_0\) leaves three free
coefficient variables, named \(r,s,h\).  The remaining \((J1),(J0)\)
coefficients give 25 residual polynomials.  Four, indexed
\(R_0,R_1,R_7,R_9\), already generate the unit ideal in

\[
 K_0[r,s,h]
\]

with degree-reverse-lexicographic (`dp`) order in the Singular input.  The
preserved machine-readable certificate proves

\[
 1=T_0R_0+T_1R_1+T_7R_7+T_9R_9.
\]

The certificate SHA-256 is
`cfbc3c39d7a28013671144f43ef76f0498542eaf6d562dd624bba3311194e4aa`.
The local verifier regenerates the 25 equations, checks the four saved
generators byte-for-byte at the algebraic-coefficient level, and evaluates
the displayed identity exactly.

## 6. Case 1 descent and exhaustive branches — local reproduction

The five added bracket layers solve linear band coefficients and leave 13
compatibility polynomials in

\[
 K_0[r,s,h,u_1,u_2,u_3].
\]

The numbers of new compatibility equations at layers
\(z^1,z^0,z^{-1},z^{-2},z^{-3}\) are respectively

\[
 0,2,2,4,5.
\]

The first nonzero compatibility polynomial factors exactly as

\[
 \kappa(s-c)^2(s+c)^2,qquad \kappa,c\in K_0^\times.
\]

The factorization file records both linear factors with multiplicity two.
Because \(K_0\) is a field and \(\kappa,c\ne0\), the split

\[
 s=c\quad\text{or}\quad s=-c
\]

is exhaustive and the branches are disjoint.

On each branch two exact equations obey

\[
 E_6-\lambda E_2=S\Lambda^2,qquad S\in K_0^\times,
\]

where \(\lambda=13/9\) or \(-13/9\) and

\[
 \Lambda=r+\alpha h+\tfrac13u_2+\gamma,qquad
 \alpha=-\tfrac23\ \text{or}\ \tfrac23.
\]

Thus every solution has \(\Lambda=0\).  Substituting the resulting affine
expression for \(r\) loses no branch.  Three further equations are discarded
only after exact verification that each is a nonzero field multiple of a kept
equation.

An invertible affine translation of \(h\) makes the coefficient of \(u_3\) in
one kept equation equal to \(b_hh\) with \(b_h\in K_0^\times\).  This is a
coordinate change, not a localization.

### The \(h=0\) strata

Before dividing by \(h\), the seven branch equations are specialized at
\(h=0\).  Separate saved multipliers prove

\[
 1=\sum_i A_iE_i\big|_{h=0}
\]

for both \(s=c\) and \(s=-c\).  Their hashes are respectively

```text
664de005e99bc6a0e61ba479ba64ed57ddbfa9c5399b955cbb12237ac70f8186
d7844c2f8edea62be4e3a7fe8a160dc4b7b70efd1262993ec6d2ed7e78722a1a
```

### The localized \(h\ne0\) strata

The equation linear in \(u_3\) is written

\[
 b_hhu_3-N(h,u_1,u_2)=0.
\]

On \(h\ne0\), substitute \(u_3=N/(b_hh)\) into the other equations and clear
powers of \(h\).  The implementation homogenizes this substitution and
removes only common powers of \(h\); this is equivalent precisely on the
already isolated open stratum \(h\ne0\).

After a checked weighted rescaling, all coefficients lie in the degree-five
subfield

\[
 L=\mathbb Q[w]/(w^5-w^4+3w^3+3w^2+26).
\]

The embedding is checked by evaluating the degree-35 polynomial and by
verifying that the rescaled coefficients occupy only powers divisible by
seven before reduction to \(L\).  Four resulting polynomials
\(F_1,\ldots,F_4\in L[h,u_1,u_2]\) have the exact identity

\[
 \boxed{h=T_1F_1+T_2F_2+T_3F_3+T_4F_4.}
\]

It follows that any localized common zero has \(h=0\), contradicting the
definition of the stratum.  The human-readable multiplier file is 89,105,967
bytes and has SHA-256
`0e48ffab32469ef8405a6945b16cf1521ddeb3c592ae4e5051968110a4dc656a`.
The involution

\[
 (h,u_1,u_2)\longmapsto(h,-u_1,-u_2)
\]

and exact row scalars transport the first branch system and identity to the
second; the transport is evaluated polynomial-by-polynomial.

### Direct unit-certificate form

The two-stage construction composes to one ordinary Nullstellensatz
certificate in each pre-division branch ideal.  Write the seven equations
after the invertible affine changes as

\[
 G_0=b_h(hu_3-N),\qquad G_j=a_ju_3+b_j \quad (j>0).
\]

Every archived equation has \(u_3\)-degree one, and the elimination removes no
common power of \(h\).  Consequently its generators satisfy the polynomial
identities

\[
 F_j=hG_j-\frac{a_j}{b_h}G_0.
\]

Thus \(h=\sum_jT_jF_j\) is already an identity \(h=\sum_iC_iG_i\) in
the pre-division ideal; it is not merely membership after localization.  If

\[
 D_i=\frac{G_i-G_i|_{h=0}}h,
 \qquad B=-\sum_iA_iD_i,
\]

the specialized certificate becomes \(1=\sum_iA_iG_i+Bh\).  Substitution
gives the single branch certificate

\[
 \boxed{1=\sum_i(A_i+BC_i)G_i.}
\]

The second branch is obtained by the checked sign involution.  The combined
multipliers need not be expanded: these formulas are a deterministic
straight-line representation of the same unit certificate.  The optional
checker `cas/verify_case1_unit_composition.py` verifies all six elimination
lifts, the degree-five descent isomorphism, and the polynomial lift of the
specialized identity against the pinned archive.

This certificate also lifts across the earlier substitution \(\Lambda=0\).
Let \(\widetilde G_i\) denote the branch equations before solving for \(r\),
and divide exactly by the monic affine polynomial \(\Lambda\) to write
\(\widetilde G_i=G_i+\Lambda Q_i\).  If
\(1=\sum_iU_iG_i\), set

\[
 X=\sum_iU_i\widetilde G_i,\qquad R=-\sum_iU_iQ_i.
\]

Then \(1=X+R\Lambda\).  The already checked relation
\(S\Lambda^2=\widetilde E_6-\lambda\widetilde E_2\), with
\(S\in K_0^\times\), gives the direct identity

\[
 1=X(1+R\Lambda)
   +\frac{R^2}{S}
      \bigl(\widetilde E_6-\lambda\widetilde E_2\bigr).
\]

Hence “original branch ideal” may mean either the seven-generator
pre-division presentation or the branch ideal before eliminating \(r\); both
have a direct unit certificate.  The factored presentation is preferable
because flattening these products only increases the archive size.

## 7. Division and saturation audit

| Operation | Potential excluded locus | Why it is safe / where the complement is handled |
| --- | --- | --- |
| Normalize by \(a_1,a_8,c_8\) | a named vertex coefficient vanishes | Exact Newton polygon means every vertex coefficient is nonzero |
| Choose square/seventh roots for \(\mu,\lambda\) | roots absent | Ground field is algebraically closed; all radicands are nonzero |
| Set \(\sigma=(\rho\lambda^3\mu)^{-1}\) | scaling determinant zero | all scaling parameters were constructed nonzero |
| Triangularly solve \((J4)\) for \(d_i\) | pivot zero | pivots are fixed nonzero rationals in characteristic zero |
| Gaussian band solves | pivot zero | pivots are explicitly tested nonzero elements of the irreducible number field; changing a pivot does not localize the coefficient variety |
| Replace the six first-block equations by the lex basis | loss under FGLM | fresh exact Singular Gröbner/FGLM run; same ideal, not a saturation |
| Split \(\kappa(s-c)^2(s+c)^2=0\) | \(\kappa=0\) or missed factor | \(\kappa,c\) are verified nonzero field elements; both factors are taken |
| Infer \(\Lambda=0\) from \(S\Lambda^2=0\) | \(S=0\) | \(S\) is a verified nonzero field scalar |
| Solve affine \(r\)-relation | coefficient of \(r\) zero | normalized square has coefficient one in \(\Lambda\); only field scalars are inverted |
| Remove proportional equations | proportionality scalar zero | each ratio is checked as a nonzero field element and the full polynomial identity is evaluated |
| Translate \(h\) | noninvertible coordinate change | coefficient \(b_h\) is a verified nonzero field scalar; inverse affine map exists |
| Substitute \(u_3=N/(b_hh)\) | \(h=0\) | both complete pre-division \(h=0\) branch systems have independent unit certificates |
| Divide out common powers of \(h\) after clearing denominators | components supported at \(h=0\) | same two pre-division certificates eliminate the complementary strata |
| Pass from degree 35 to degree 5 coefficients | wrong conjugate / lost embedding | exact subfield identity and weighted divisibility checks; a unit/membership identity over \(L\) remains one after extension to \(K_0\) |
| Sparse Macaulay linear solves | a chosen determinant might vanish | matrices are over fixed exact fields, not parameter rings; selected pivots/minor are exact nonzero constants |
| Select 1925 of 2010 scalar rows | ignored equations | selected minor has determinant \(10\bmod71\), proving nonsingularity over \(\mathbb Q\); recovered solution is then checked against all 2010 rows |

There is no unexamined saturation factor.  The only parameter localization is
at \(h\), and its complement is eliminated before localization.

## 8. Certificate formats and independent checking

The primary computations use:

| Stage | Base ring | Variables | order / method | certificate |
| --- | --- | --- | --- | --- |
| First block | \(\mathbb Q\) | \(a_2,\ldots,a_7\) | `dp` Gröbner, then `lp` FGLM | 56-element degree-order basis and 6-element lex basis |
| Case 2 | \(K_0\) | \(r,s,h\) | Singular `dp`; Macaulay degree 8 | four-generator unit identity |
| Case 1, \(h=0\) | \(K_0\) | \(u_1,u_2,u_3\) after branch/affine changes | exact lifted Singular identities | two unit identities |
| Case 1, direct branch certificate | \(K_0\), with the hard solve descended to \(L\) | \(h,u_1,u_2,u_3\) | factored composition of the specialized identity and the 2010-by-1925 membership system | \(1=\sum U_iG_i\) |

The 89 MB hard multiplier text is not logically necessary as an authoritative
stored artifact.  A compact reconstructible form consists of:

1. the exact four-generator system, coefficient-field polynomial and ordered
   385-monomial multiplier support;
2. the deterministic rule expanding the 402 field rows to 2010 rational
   scalar rows and 1925 scalar columns;
3. the 1925 selected scalar-row indices and the exact target vector;
4. the modular determinant check for the selected minor;
5. an exact solve followed by verification against all 2010 scalar rows; and
6. hashes pinning every input and, optionally, the reconstructed solution and
   human-readable expansion.

The determinant \(10\bmod71\) proves that the selected rational minor is
invertible because its denominators are units modulo 71.  It does not by
itself prove consistency of the 85 omitted scalar rows: the final all-row
exact replay is essential.  Hashes provide integrity and provenance, not a
replacement for either algebraic check.

The primary verifier imports the generation modules and reevaluates serialized
identities.  The independent checker in [cas/verify_h_certificate_independent.py](cas/verify_h_certificate_independent.py)
does not import FLINT, SymPy, Singular, or the generating code.  It parses the
text certificate, implements arithmetic in \(L\) directly with `gmpy2`,
expands 13,410 field products to 335,250 rational scalar products, and checks
all 2010 scalar rows.  The smaller composition checker verifies that this
membership identity and the \(h=0\) certificate constitute a direct unit
certificate for the seven pre-division equations.

The complete replay command and hashes are in [cas/README.md](cas/README.md).

## 9. Final contradiction

Case 2 is empty by a unit identity.  In Case 1, every solution belongs to one
of \(s=\pm c\), and each branch ideal has the direct unit identity constructed
above.  Therefore both Newton polygons from GGHV Proposition 4.3 are
impossible.
Together with GGHV Theorem 2.1, this gives the larger-coordinate degree bound
125 stated in [DEGREE_FRONTIER_125.md](DEGREE_FRONTIER_125.md).

## 10. Conceptual replacement audit — repository-original, presently negative

Several compact structures are visible:

- the first eliminant is sparse in \(a_7^7\), exposing a cyclic order-seven
  grading and the degree-five subfield;
- the two branches are a sign involution rather than unrelated systems;
- the forced relation is a square \(S\Lambda^2\);
- the hard calculation is linear algebra for ideal membership once the four
  generators and multiplier support are fixed.

These observations compress the computation but do not yet replace it by an
intersection, conductor, semigroup, or lattice contradiction.  In particular,
the admissible chain already passes the known determinant and divisibility
tests; the surviving obstruction depends on coefficient moduli over a
degree-five field.

The log-boundary audit now narrows the missing conceptual step further.  One
must correct the final transformed bracket by
`R_actual=K+3H+div(X^2)`; the two terminal dicriticals then have normal
indices three and five, not zero ramification.  Their exact degree-twelve
residues have a priori normalization-cover degree `1,2,4`.  General
polynomial-decomposition remainder ideals exclude degrees `2,4` in Case 2
without using `(J0)` or the `(J1)` compatibility equations, so that residue
is birational at that intermediate stage.  The remaining degree-twelve row
is now empty as well: the seven residual `(J1)` compatibility cubics,
localized at the forced endpoint `G_12 != 0`, generate the unit ideal over
the exact degree-35 field without `(J0)`.  Separately, `(J0)` and `(J1)`
reduce through
`H=gcd(C',G')`, and the initial coefficients force `t|H`; the degree-one gcd
row has only one possible origin-order pattern.  At the opposite end, three
low coefficients of `remainder(G',C')` and the terminal `(J0)` coefficient
generate an exact unit ideal, excluding gcd degree seven without a residual
`(J1)` compatibility equation.  These reductions are proved and checked in
[FRONTIER_LOG_SCALE_AUDIT.md](FRONTIER_LOG_SCALE_AUDIT.md).  The Case-2
endpoint certificate is a smaller coefficient contradiction, not yet a
pure ramification/residue-degree replacement.

The analogous Case-1 decomposition sieve requires data not present in the
archived necessary truncation: the full polygons continue from the recorded
`P:z^-5,Q:z^-4` bands down to `P:z^-8,Q:z^-12`, and those terms enter the
alternate-chart residue.  They must be derived, not silently set to zero.

No Smith-normal-form obstruction of the boundary intersection matrix and no
conductor inequality recovering the hard identity has been found.

The promising repository-original question is whether the order-seven grading
and the two approximate-root filtrations define incompatible value semigroups.
Until that is proved, it is a conjectural continuation, not a strengthening of
the marked-point theorem or a family theorem.
