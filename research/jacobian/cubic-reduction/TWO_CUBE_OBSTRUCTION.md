# Two arbitrary new cubes cannot shorten the saved89-form identity

The saved89-form padded companion identity cannot be shortened by
retaining any subset of the old forms, adding at most **two arbitrary
new linear cubes**, and re-solving all retained output coefficient
columns. This holds over every characteristic-zero extension field.
The new forms are not confined to a finite pool or to the saved scalar
cubic span S.

This is a replacement obstruction for the one identity in
[KCUB5](POWER89.md). It is not a global minimum89 for the tensor or for
cubic-linear Keller maps. Three or more new forms remain possible.

## Output-column rank separates the cases

Write P=sum beta_j c_j for the89 independent old cubes. A shorter
representation with two new cubes removes at least three old forms.
If its two new output coefficient columns are independent, passing to
the quotient by S forces both new cubes into S. The
[fixed-span theorem](SPAN_OPTIMALITY.md) then forbids shortening.

Otherwise the two output columns are proportional. Their common
projective direction must equal that of every deleted old column.
Over an algebraic closure, absorb their scalar coefficients into cube
roots of the new linear forms. The scalar contribution is

g=ell^3+m^3 in S,

with a common output direction d. If beta_j=s_j d, deleting column j
requires a_j(g)=s_j in the unique old-cube coordinate expansion.
It suffices to test a triple of deleted columns; there are214 triples
with proportional original coefficient columns.

## A mixed-derivative support lemma

If ell and m are proportional, they combine to one cube, already
excluded by the fixed-span theorem. Suppose they are independent.
For a missing edge(i,k) in the old source co-occurrence graph, every
cubic in S has mixed second derivative zero. Therefore

\[
0=\partial_i\partial_k g
 =6\ell_i\ell_k\ell+6m_i m_k m.
\]

Independence forces ell_i ell_k=m_i m_k=0. Each new form's support is
thus a clique, even though its cube need not lie in S individually.
The15 maximal cliques cover all225 nonempty cliques. There are120
unordered pairs of maximal cliques, including equal pairs.

On each support pair, introduce disjoint coefficient variables for
ell and m. The membership equations for g in S must use a **common
linear projection** of the ambient cubic space; separately normalized
residual equations from different supports cannot simply be added.
Build the combined residual equations and the89 coordinate cubics.
For every eligible triple j,k,l impose

\[
a_j(g)=s_j,\qquad a_k(g)=s_k,\qquad a_l(g)=s_l.
\]

Seven support pairs exclude every triple because some required
coordinate vanishes identically. The other113 pairs give607 systems.
Every system has an explicit rational unit identity. The largest
multiplier degree is12; the607 identities contain29,023 multiplier
monomials in total. These are exact exclusions over arbitrary
characteristic-zero extensions, without coefficient-height bounds.

## Independent verification

[Primary construction](power_secant_two.py) uses a FLINT cube-space
projection and Singular unit-ideal calculations.
[Certificate export](certify_power_secant_two.py) retains each multiplier
identity. [The manifest](power-secant-two/certification.json) pins all
113 case certificates and their source inputs.

[The independent checker](verify_power_secant_two.py) uses only standard-
library Fraction arithmetic and the retained sparse polynomial engine.
It independently reconstructs the89-cube coordinates, all combined
membership row spaces, every coordinate map modulo membership, all214
coefficient triples, all120 support pairs and their exhaustive clique
coverage. It checks that all102 missing source edges have zero mixed
second derivatives on S, then multiplies out all607 unit identities.
[Validation](power-secant-two/validation.json).

The initial18 equal-support systems are retained separately as a
weaker discovery calculation, not a premise of the full proof.

## Subsequent construction attempts

The fixed930-form pool was also tested for three-new-form exchanges:

- Rank-two residual circuits:2,061 direction circuits and50,138 actual
  triples; at most three old forms cancel, so no shorter identity is
  found. [Primary checkpoint](power-three-rank2/checkpoint.json).
- Rank-one residual exchanges:123,379 tasks, including two outside forms
  with a genuinely new inside form. No candidate passes the necessary
  scalar-projection test for four canceled columns.
  [Primary checkpoint](power-three-rank1/checkpoint.json).
- Tangent cubics p^2 q, represented by three cubes, with930 fixed p
  directions and all q solved linearly from p^2 q in S. Of these,
  455 tangent spaces are nonzero. None has four potentially nonzero
  cancellation coordinates in one proportional output direction.
  New q forms were not restricted to the old pool.
  [Primary checkpoint](power-tangent-search/checkpoint.json).

These are retained discovery searches, not an unrestricted three-form
obstruction. No smaller map was found. A next construction test must
introduce at least three new forms, enlarge the current search classes,
or change the underlying companion. The unrestricted fixed-tensor
shared-cube interval remains[34,89].
