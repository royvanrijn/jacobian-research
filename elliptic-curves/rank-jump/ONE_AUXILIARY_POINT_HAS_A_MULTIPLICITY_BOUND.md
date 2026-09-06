# One auxiliary point is limited by elliptic multiplicity

A single rational point on an auxiliary curve cannot supply arbitrarily
many independent specialization directions through prescribed maps.
Modulo the original generic subgroup, the number is bounded by
\[
\boxed{\operatorname{rank}_{\mathbb Z}\operatorname{Hom}_K(J_C,E).}
\]
When \(\operatorname{End}_K(E)=\mathbb Z\), this is the multiplicity of
\(E\) in \(J_C\), and is at most the genus of \(C\).

For both retained genus-five auxiliary constructions, exact finite-field
certificates now show that this multiplicity over \(\mathbb Q\) is
**two**, not five. Their genus-three quotient contributes no rational
elliptic factor isogenous to either control curve. This is an exact
bound on independent maps, not a new Mordell--Weil rank or a proof of
rational solubility.

## The map bound, without an auxiliary base point

Let \(K\) be a characteristic-zero field, \(C/K\) a smooth projective
geometrically connected curve of genus \(g\), and \(E/K\) an elliptic
curve. Every morphism \(\phi:C\to E\) induces
\[
\phi_*:J_C\longrightarrow E.
\]
The kernel of this assignment consists precisely of constant morphisms,
identified with \(E(K)\). This holds even when \(C(K)\) is empty:
over an algebraic closure the Jacobian universal property makes a map
with zero induced homomorphism constant; descent makes that constant
point \(K\)-rational.

Moreover, the image has finite index in
\(\operatorname{Hom}_K(J_C,E)\). Choose any \(K\)-rational divisor \(D\)
of positive degree \(d\), for example a closed point as a divisor.
The morphism
\[
C\longrightarrow J_C,\qquad P\longmapsto[dP-D]
\]
induces multiplication by \(d\) on \(J_C\). Composing with any
\(J_C\to E\) shows that the image contains
\(d\operatorname{Hom}_K(J_C,E)\). Therefore
\[
\operatorname{rank}\bigl(\operatorname{Mor}_K(C,E)/E(K)\bigr)
=\operatorname{rank}\operatorname{Hom}_K(J_C,E)=h.
\]
The needed universal properties, including the formulation without a
rational base point, are in
[Milne, *Jacobian Varieties*, Propositions 6.1 and 6.4](https://www.jmilne.org/math/xnotes/JVs.pdf).

If \(\operatorname{End}_K(E)=\mathbb Z\), complete reducibility gives
\[
J_C\sim_K E^h\times B,\qquad \operatorname{Hom}_K(B,E)=0,
\]
so \(h\le g\). The decomposition up to isogeny and the finiteness of
endomorphism groups are standard facts in
[Milne, *Abelian Varieties*](https://www.jmilne.org/math/xnotes/AVs.pdf).

In particular, \(\operatorname{End}_K(E)=\mathbb Z\) for
\(K=\mathbb Q\) or \(\mathbb Q(u)\). The differential of an endomorphism
embeds its ring into \(K\) in characteristic zero. Endomorphisms are
integral over \(\mathbb Z\); an element of these fields algebraic over
\(\mathbb Q\) lies in \(\mathbb Q\), hence an integral such element lies
in \(\mathbb Z\). This argument concerns endomorphisms defined over \(K\);
it makes no geometric non-CM assumption.

## Specialization and the original generic subgroup

Now take \(K=\mathbb Q(u)\). Suppose \(E\) is the original family and
\(G\subset E(K)\) spans its full generic Mordell--Weil group over
\(\mathbb Q\). Let \(\phi_1,\ldots,\phi_n:C\to E\) be prescribed over \(K\).
At a specialization \(u_0\) where the families, maps, and relevant
sections specialize as required, take **one**
\(P\in C_{u_0}(\mathbb Q)\). Then
\[
\dim_{\mathbb Q}
\left\langle[\phi_{1,u_0}(P)],\ldots,[\phi_{n,u_0}(P)]\right\rangle
\le h\le g,
\]
where the classes lie in
\((E_{u_0}(\mathbb Q)/\operatorname{sp}_{u_0}G)\otimes\mathbb Q\).

To prove the bound, consider the linear map sending the \(i\)-th
coordinate vector to \(\phi_{i,*}\). Every integral relation in its
kernel gives a constant morphism
\[
\sum_i n_i\phi_i=R,\qquad R\in E(K).
\]
After specialization and evaluation at \(P\), the same relation holds
modulo the specialized generic rational span. These relations leave at
most \(\dim\langle\phi_{i,*}\rangle\le h\) dimensions.

Thus a one-point explanation of \(+8,+10,+14\) by this type of auxiliary
curve requires elliptic multiplicity at least \(8,10,14\), respectively,
and genus at least those values. Explaining a recent \(+9\) this way
requires multiplicity at least nine. A genus-two auxiliary curve with
one rational point cannot by itself supply nine independent directions
through fixed original-family maps.

These are necessary conditions, not sufficient ones. Having \(h\) maps
does not make their evaluated points independent. Their specialized
images may collapse, or the common curve may have no rational point.

### Scope that must stay explicit

- With \(\ell\) auxiliary rational points, the same argument gives the
  weaker upper bound \(\ell h\). A single elliptic auxiliary curve can
  carry many independent rational points; this theorem does not bound
  their number by its genus.
- Translations by original generic sections vanish in the stated
  quotient. If the proposed generic subgroup is smaller than the full
  generic rational span, omitted generic directions must be counted
  separately.
- A finite cover of the **parameter line** is different. Its generic
  fibre over \(\mathbb Q(u)\) is zero-dimensional, not the positive-
  dimensional auxiliary curve \(C/K\) in this theorem. A quadratic
  parameter cover can create several sections, as the earlier
  split-cubic example proves; this is not contradicted by genus zero
  of its total parameter curve.
- New maps appearing only on a special auxiliary fibre are not
  prescribed generic maps. Extending the theorem to such an event
  requires its own specialization argument.
- On a fixed curve over \(\mathbb Q\), all rational constant translations
  belong to \(E(\mathbb Q)\). The fixed-control calculation below counts
  independent maps modulo those translations. It does not pretend
  that this quotient is the original production exceptional quotient.

## Exact multiplicity on the two genus-five controls

The [preceding construction](SYNCHRONIZING_COVERS_CAN_ADD_AN_OBSTRUCTION.md)
gives, for each alignment \(\phi={\rm id},{\rm swap}\),
\[
J_{C_{+,\phi}}\sim_{\mathbb Q}E_{\rm rat}^2\times J_{D_\phi},\qquad
J_{C_{-,\phi}}\sim_{\mathbb Q}E_{\rm Sha}^2\times J_{D_\phi}.
\]
The genus-three factor was previously unconstrained as an elliptic
quotient source. We now prove
\[
\operatorname{Hom}_{\mathbb Q}(J_{D_\phi},E_{\rm rat})
=\operatorname{Hom}_{\mathbb Q}(J_{D_\phi},E_{\rm Sha})=0
\]
for both alignments.

The [initial protocol](AUXILIARY_ELLIPTIC_MULTIPLICITY_PROTOCOL.json)
fixes primes \(3,7,11\), counting only over the first three extensions
of an eligible prime. The identity model has good displayed reduction
at 11, where the genus-three curve counts are
\[
\#D(\mathbb F_{11})=16,\quad
\#D(\mathbb F_{11^2})=130,\quad
\#D(\mathbb F_{11^3})=1360.
\]
Newton identities give
\[
W_{\rm id}(T)=T^6+4T^5+12T^4+36T^3+132T^2+484T+1331.
\]
The elliptic Frobenius polynomials and division remainders are

| Elliptic target | Frobenius polynomial at 11 | Remainder of \(W_{\rm id}\) |
|---|---|---|
| \(E_{\rm rat}\) | \(T^2+5T+11\) | \(392T+1540\) |
| \(E_{\rm Sha}\) | \(T^2-5T+11\) | \(952T-3740\) |

Both remainders are nonzero. The distinct elliptic traces also prove
that the two elliptic curves are not isogenous over \(\mathbb Q\).

The swapped model has bad displayed reduction at all three original
primes. That result is preserved as **UNKNOWN** in the first artifact.
A [separate bounded completion](AUXILIARY_MULTIPLICITY_COMPLETION_PROTOCOL.json)
allows only \(13,17,19\); the first, 13, settles it. Counts \(8,188,2168\)
give
\[
W_{\rm swap}(T)=T^6-6T^5+27T^4-100T^3+351T^2-1014T+2197.
\]
Both elliptic polynomials at 13 equal \(T^2+2T+13\), and the remainder
is \(-432T+1248\ne0\).

### Why one reduction proves the exclusion

A nonzero rational homomorphism from \(J_D\) to an elliptic curve is
surjective and, up to isogeny, makes that elliptic curve a factor of
\(J_D\). At a common good reduction prime its Frobenius polynomial must
divide that of \(J_D\). The nonzero remainders contradict this necessary
condition. No converse, geometric isogeny assertion, or unproved Tate
conjecture is used.

Consequently
\[
\operatorname{rank}\operatorname{Hom}_{\mathbb Q}(J_{C_{+,\phi}},E_{\rm rat})
=\operatorname{rank}\operatorname{Hom}_{\mathbb Q}(J_{C_{-,\phi}},E_{\rm Sha})
=2.
\]
The cross-target Hom spaces are zero. Every rational map from a signed
genus-five curve to its corresponding elliptic control is, after
multiplying by an integer, a linear combination of its two natural maps
plus a rational constant. The genus-three factor cannot supply a third
independent map to that target over \(\mathbb Q\).

This does not determine any of these Jacobians' full Mordell--Weil
ranks or their geometric Hom spaces after extending the ground field.
Nor does it make the negative-sign curves rationally soluble: their
projections to the known Sha torsors still exclude rational points.

## Ranked implications for the rank-jump question

1. **Incidence mechanism, necessary capacity bound:** count the
   multiplicity of the target elliptic curve in an auxiliary Jacobian,
   rather than treating all its genus or rank as available directions.
   Here a genus-five object has only two independent maps to the target.
2. **Still viable solubility mechanisms:** several independent elliptic
   quotients evaluated at a common rational point, or sections over a
   finite parameter cover. Both require actual rational solubility and
   independence modulo the original generic subgroup.
3. **Weak explanations:** a large auxiliary Jacobian rank, a large genus,
   or many formulas obtained by translating existing maps. None implies
   additional independent directions from one point.
4. **Missing production implication:** exhibit an original-family
   construction with enough target-elliptic multiplicity (or a genuine
   finite-cover section block), prove its rational specialization
   condition, and verify the surviving quotient independence. The
   present small controls supply no \(+8\) or larger construction.
5. **For Agent 1:** this is a structural **incidence** upper bound on a
   proposed mechanism, not a high-rank selector. Auxiliary rational-point
   tests remain **solubility**; chart exposure remains **visibility**.
   No scoring or search protocol change is justified.

## Reproducibility

- [Initial bounded finite-field result, including the preserved UNKNOWN](../../artifacts/generated-results/elliptic-curves/rank_jump_auxiliary_elliptic_multiplicity_v1.json)
- [Swapped-case completion at 13](../../artifacts/generated-results/elliptic-curves/rank_jump_auxiliary_multiplicity_completion_v1.json)
- [Independent finite-field verification](../../artifacts/generated-results/elliptic-curves/rank_jump_auxiliary_multiplicity_verification_v1.json)

The producer uses Sage finite fields. The independent verifier constructs
its own fields in pure Python with irreducible quadratic/cubic moduli,
enumerates square sets, counts all affine and infinite points, repeats
the Newton identities, and divides the Frobenius polynomials exactly.
Its largest field actually used has 2,197 elements. No rational point
search, production sweep, or Jacobian descent is run.

    sage -python elliptic-curves/rank-jump/auxiliary_elliptic_multiplicity.py check
    sage -python elliptic-curves/rank-jump/auxiliary_multiplicity_completion.py check
    python3 elliptic-curves/rank-jump/verify_auxiliary_multiplicity.py check

Both captures have exclusive checkpoints and 30-second limits. Changes
are confined to new rank-jump analysis files and the analysis index.
