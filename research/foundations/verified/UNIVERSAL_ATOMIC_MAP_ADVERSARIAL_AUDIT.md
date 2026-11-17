# Adversarial audit of the universal atomic map

This note attacks the
[universal absolute atomic-map theorem](UNIVERSAL_RELATIVE_KELLER_MAP.md#31-one-absolute-atomic-map-in-every-rank)
at every logical gate.  It does not replace that canonical proof.  Its purpose
is to distinguish consequences of exact algebra from imported geometric
theorems, exhibit the failure boundary of every hypothesis, and pin explicit
connected, split, and disconnected witness fibers.

The verdict is:

\[
\boxed{\text{No obstruction was found.  UAM1 survives for every }N\ge3
\text{ over characteristic-zero fields.}}
\]

The qualifications are genuine: the representing target is noncanonical,
the claim is not made in rank two or in positive characteristic, and neither
ambient-dimension minimality nor descent from `BS_N` is asserted.

## 1. The fixed map and its inverse equation

For `N>=3`, retain

\[
 {\bf u}=(u_4,\ldots,u_N)
\]

as unchanged coordinates and apply the three-variable quadratic-gauge map
with seed

\[
 G_{\bf u}(S)=S+S^3+\sum_{j=4}^Nu_jS^j.
\]

This gives

\[
 {\cal U}_N({\bf u},x,y,z)
 =({\bf u},{\cal K}_{N,{\bf u}}(x,y,z)):
 \mathbb A^N\longrightarrow\mathbb A^N.
\]

At target `({\bf u},pi,b,c)`, its inverse polynomial is

\[
 E(S)=S+bS^2+\pi S^3+\sum_{j=4}^Nu_j\pi^jS^j-\frac c2. \tag{1}
\]

All attacks below concern one of the implications from these two formulas.

## 2. Gate-by-gate attack

| Gate | Adversarial question | Resolution |
|---|---|---|
| Polynomiality | Does promotion introduce denominators or require `u_N` to be invertible in the map itself? | No.  Every `u_j` occurs polynomially in the displayed coordinates.  The map extends across all parameter hyperplanes.  Inverting `u_N`, `pi`, and the discriminant is needed only to select the degree-`N` finite-etale inverse cover. |
| Full Jacobian | Can derivatives of the three moving coordinates with respect to the promoted parameters spoil the determinant? | No.  The upper-right block is zero because the promoted outputs are exactly the promoted inputs.  The full matrix is `[[I,0],[*,D_v K]]`, so its determinant is the already-proved vertical determinant one. |
| Generic degree | Could reconstruction lose a generic root or create extra sheets? | No.  Work over `k(u_4,...,u_N)`.  The existing three-variable function-field comparison applies to the single seed `G_u`; (1) has leading coefficient `u_N*pi^N` (`pi` for `N=3`) and the reconstruction identifies its root extension with the source function field.  Thus the degree is exactly `N`. |
| Parameter compiler | Can an arbitrary normalized polynomial fail to have the form (1)? | Only on `h_3=0`.  On `h_3!=0`, the inverse formulas `pi=h_3`, `b=h_2`, `c=-2h_0`, and `u_j=h_j/h_3^j` recover every coefficient exactly. |
| Admissible translation | Might every translation of a squarefree presentation have vanishing first or third jet? | No in characteristic zero.  `P'` and `P^[3]` are nonzero polynomials for `deg(P)>=3`; their product cannot vanish at every point of an infinite field.  Lean already formalizes this as `exists_admissible_translation`. |
| Abstract algebra | Might a disconnected finite-etale algebra lack a one-polynomial presentation? | No over an infinite field.  The primitive elements form a nonempty open; the Lean development proves the stronger product-of-fields monogenicity statement and constructs a power basis. |
| Complete fiber | Does coefficient equality prove only a bijection of geometric roots? | No.  The represented-fiber theorem is natural over every commutative test algebra.  Translation and multiplication by the nonzero scalar `P'(a)^(-1)` identify the quotient algebra with `k[T]/(P)`, so the result is scheme-theoretic. |
| Etaleness | Can the compiled quotient acquire nilpotents? | Not when `P` is squarefree.  Translation and nonzero scaling preserve separability.  A repeated-root input still satisfies the coefficient identity, but its discriminant is zero and it is correctly excluded from the finite-etale open. |
| Monodromy | Does promoting the parameters reduce `S_N` to a subgroup? | No.  The coefficient change from `({bf u},pi,b,c)` to `(h_0,h_2,...,h_N)` is birational on `pi*u_N!=0`.  The presentation base of Section 1 is exactly `V_N times A1`, and its root cover is the universal monic-polynomial root cover restricted to a nonempty open.  Deleting the purely transcendental origin coordinate does not change geometric monodromy. |
| Absolute atomicity | Is `S_N` enough, or is integer primality being used implicitly? | `S_N` in its natural action is primitive for every `N>=2`.  The primitive-monodromy theorem rules out proper intermediate function fields, forces a degree-one factor in every polynomial decomposition, and turns that factor into an automorphism.  No primality assumption on `N` occurs. |
| Stable atomicity | Could adjoining identity coordinates create an intermediate field? | No.  Stabilization replaces `K subset L` by `K(t_1,...,t_r) subset L(t_1,...,t_r)` with the same Galois closure group and point stabilizer. |
| Constant extension | Can geometric monodromy or primitivity drop after extending the characteristic-zero constants? | No.  Both are computed after passage to algebraically closed constants; the same geometric generic cover is obtained. |

The monodromy and atomicity rows are written theorem arguments, not outputs of
the symbolic checker.  The checker certifies the algebraic rows and guards the
failure boundaries below.

### 2.1 External theorem cross-check

The noncomputational inputs were checked against the following standard
sources.

- [Stacks Project, Section 10.143](https://stacks.math.columbia.edu/tag/00U0)
  identifies finite etale algebras over a field with finite products of
  finite separable field extensions.  The field primitive-element input is
  [Stacks Project, Lemma 9.19.1](https://stacks.math.columbia.edu/tag/030N);
  the product-algebra strengthening used here is proved in
  `AbstractFiniteEtale.lean`.
- [Serre, *Topics in Galois Theory*, Theorem
  4.4.5](https://indico.math.cnrs.fr/event/11410/attachments/4760/7351/SerreTopicsGaloisTheory.pdf)
  supplies the symmetric monodromy of a Morse polynomial.  The repository's
  generic-slice argument reduces the two-parameter pencil to that theorem.
- [Osada, *The Galois groups of the polynomials
  \(X^n+aX+b\)*](https://doi.org/10.1016/0022-314X(87)90029-1)
  supplies the all-rank arithmetic `S_N` witness family `T^N-T-1` used
  below.
- [Zieve--Mueller, *On Ritt's polynomial decomposition
  theorems*](https://arxiv.org/abs/0807.3578), Section 2, records the
  intermediate-field/monodromy-block correspondence used to attack
  compositional decompositions.
- [Stacks Project, Zariski's Main
  Theorem](https://stacks.math.columbia.edu/tag/00Q9) is the external
  geometric input in the degree-one Keller lemma used by the atomicity
  theorem.

These references support the theorem-level gates; they do not turn the
symbolic witness search into a proof of monodromy or atomicity.

## 3. Counterexample attempts that hit the stated boundary

The following attacks fail for precise reasons.

1. **Set `pi=0`.**  Equation (1) loses every term of degree at least three.
   This is a genuine degree-drop target, but it lies outside the universal
   incidence open.
2. **Set `u_N=0`.**  For `N>=4`, the top term vanishes.  The ambient map
   remains polynomial and Keller; only this target chart ceases to describe a
   degree-`N` fiber.  A compiled monic polynomial always has `u_N!=0`.
3. **Use the translation `a=0` for `P=T^N-2`.**  Then `P'(0)=0`.  Translation
   is not optional, but the admissibility theorem supplies `a=1`.
4. **Use `P=T^4+T` at `a=0`.**  Here `P'(0)=1` but `P^[3](0)=0`.  This shows
   why avoiding only the first-jet divisor is insufficient.
5. **Use `P=(T-1)^2(T^2+1)`.**  At `a=-3` both required jets are nonzero and
   the compiler identity holds, but the discriminant is zero.  The selected
   fiber is not finite etale, exactly as the theorem states.
6. **Demand `S_N` for every specialized fiber.**  This is false and was never
   claimed.  For example `T^4-2` has a proper arithmetic Galois group, while
   the generic inverse cover of the fixed map still has geometric monodromy
   `S_4`.
7. **Identify two presentations of the same abstract algebra.**  They can
   compile to different targets and different surrounding stable map data.
   This defeats canonical descent from `BS_N`, not existence inside one fixed
   map.

Thus every attempted counterexample either violates a displayed open
condition or attacks a stronger statement than UAM1.

## 4. Exact witness cards

Targets are ordered as

\[
 (u_4,\ldots,u_N;\pi,b,c).
\]

There is first an all-rank connected stress family.  Osada's theorem gives

\[
 \operatorname{Gal}(T^N-T-1/\mathbb Q)=S_N
 \qquad(N\ge3).
\]

At `a=1`, its normalized translate is

\[
 \frac{(1+S)^N-(1+S)-1}{N-1},
\]

and its target has

\[
 \pi=\frac{N(N-2)}6,\qquad b=\frac N2,\qquad
 c=\frac2{N-1},
\]

\[
 u_j=\frac{\binom Nj}{N-1}
      \left(\frac6{N(N-2)}\right)^j
 \quad(4\le j\le N).                                  \tag{2}
\]

Thus every fixed map `U_N` has an explicit connected fiber whose arithmetic
Galois group is already the full symmetric group.  This does not prove
generic *geometric* monodromy by computation; it is an independent
specialization stress test.  The checker verifies (2), the inverse identity,
and nonzero discriminant through rank twelve.

| Algebra type | `P(T)` | `a` | Compiled target |
|---|---|---:|---|
| connected rank 3 | `T^3-2` | 1 | `(1/3,1,2/3)` |
| connected rank 4 | `T^4-2` | 1 | `(1/4;1,3/2,1/2)` |
| split rank 4 | `(T-1)(T-2)(T-3)(T-4)` | 0 | `(-25/2;1/5,-7/10,24/25)` |
| product rank 4 | `(T^2-2)(T^2-3)` | 1 | `(-27/32;-2/3,-1/6,2/3)` |
| connected rank 5 | `T^5-2` | 1 | `(1/16,1/160;2,2,2/5)` |
| product rank 5 | `(T^2-2)(T^3-3)` | 1 | `(-1715/4096,2401/32768;-8/7,-1/7,4/7)` |
| connected rank 6 | `T^6-2` | 1 | `(81/4000,243/100000,243/2000000;10/3,5/2,1/3)` |
| product rank 6 | `(T^2-2)(T^4-3)` | 1 | `(-26/81,8/81,-8/729;-3/2,0,1/2)` |

For the connected family, `T^N-2` is Eisenstein at two.  The split witness
has discriminant `144`.  The product factors are irreducible Eisenstein
polynomials and are pairwise coprime.  The checker verifies every target,
the normalized inverse identity, and nonzero discriminant.

For example, the connected sextic target selects

\[
 \frac{(1+S)^6-2}{6}
 =\frac{S^6}{6}+S^5+\frac52S^4+\frac{10}{3}S^3
   +\frac52S^2+S-\frac16.
\]

The disconnected sextic target selects

\[
 -\frac18S^6-\frac34S^5-\frac{13}{8}S^4
 -\frac32S^3+S-\frac14,
\]

which is `P(1+S)/P'(1)` for
`P=(T^2-2)(T^4-3)`.

## 5. Formalization ledger

The following parts are now checked in Lean without `sorry`:

1. existence of an admissible translation:
   `FiniteEtaleKeller/Admissibility.lean`;
2. monogenicity of characteristic-zero finite-etale algebras:
   `FiniteEtaleKeller/AbstractFiniteEtale.lean`;
3. natural scheme-level represented-fiber equivalence:
   `FiniteEtaleKeller/GeneralGaugeNormalization.lean` and
   `FiniteEtaleKeller/GeneralGaugeRealization.lean`;
4. the determinant of an abstract unchanged-coordinate block matrix
   `[[I,0],[C,D]]`:
   `FiniteEtaleKeller/UniversalPromotedBlock.lean`;
5. the literal promoted polynomial map on the coordinate type
   `Fin 3 ⊕ {j // 4 ≤ j ≤ N}`, the identification of its actual Jacobian
   with `[[D_x,D_u],[0,I]]`, the cardinality `N` of that type, and the
   determinant-one specialization of the universal gauge:
   `FiniteEtaleKeller/UniversalPromotedMap.lean` and
   `FiniteEtaleKeller/UniversalPromotedGauge.lean`;
6. the unchanged-parameter coefficient compiler, its exact normalized
   translated-polynomial identity, full selected inverse degree, and nonzero
   compiled top parameter:
   `FiniteEtaleKeller/UniversalParameterCompiler.lean`;
7. nonzero scalar invariance of polynomial quotient algebras and the
   canonical equivalence
   `K[S]/(P(a+S)/P'(a))` with `K[T]/(P(T))`:
   `FiniteEtaleKeller/UniversalParameterQuotient.lean`;
8. connected, split, and product quartic target identities:
   `FiniteEtaleKeller/UniversalParameterWitnesses.lean`.

Not formalized in Lean:

1. the explicit identification of a selected full fiber of the literal
   promoted map with the quotient produced by the compiler;
2. geometric degree `N` for the literal promoted map;
3. generic `S_N` monodromy of the universal root cover;
4. primitive-monodromy absolute and stable atomicity;
5. the external degree-two exclusion.

Accordingly UAM1 remains `proof_type = hybrid` and
`formal_verification = false`.  The Lean development certifies its algebraic
compiler and fiber core, not the complete monodromy-to-atomicity chain.

## 6. Reproduction

Run

```bash
.venv/bin/python scripts/verify_universal_relative_keller_map.py
cd formal/finite-etale-keller
lake build FiniteEtaleKeller.UniversalParameterCompiler
lake build FiniteEtaleKeller.UniversalParameterQuotient
lake build FiniteEtaleKeller.UniversalParameterWitnesses
lake build FiniteEtaleKeller.UniversalPromotedBlock
lake build FiniteEtaleKeller.UniversalPromotedMap
lake build FiniteEtaleKeller.UniversalPromotedGauge
```

The Python checker also tests the parameter-boundary, bad-translation, and
repeated-root attacks.  The Lean modules print their axiom dependencies;
none contains `sorryAx`.
