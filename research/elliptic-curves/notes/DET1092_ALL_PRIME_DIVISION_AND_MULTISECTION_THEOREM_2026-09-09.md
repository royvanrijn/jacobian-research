# Every genuine bisection is generically independent: the all-prime barrier

## Results and precise scope

**New deductions/applications of established modular and covering theory.**
Let `E/Q(t)` be the determinant1092 parent and let
`M=E(Qbar(t))=M17`, the already certified full, torsion-free geometric
Mordell--Weil group. All its generators are defined over `Q(t)`.

For a prime `ell` and `Z in M \ ell*M`, the smooth projective curve
`T_(ell,Z)` of solutions to `[ell]Q=Z` is geometrically connected and has

\[
\boxed{\deg(T_{\ell,Z}/\mathbf P^1)=\ell^2,
\qquad g(T_{\ell,Z})=(\ell-1)(11\ell-13).}
\tag{1}
\]

The connected curve of nonzero `ell`-torsion points has degree `ell^2-1`
and genus one larger. Thus for `ell=2,3,5` the division genera are
`9,40,168`; the torsion genera are `10,41,169`.

**Constructive consequence.** For any nonconstant base change
`C -> P1`, either of the following hypotheses suffices:

- `g(C)<9`, with no degree or ramification restriction;
- `deg(C/P1)=2`, with no genus or ramification restriction.

Under either hypothesis `E(Qbar(C))` has no torsion and the inherited
group is fully saturated, at **every** integer:

\[
\boxed{nQ\in M\Longrightarrow Q\in M\quad(n\ge1).}
\tag{2}
\]

Consequently every genuine geometrically irreducible bisection, and every
genuine multisection of degree greater than1 whose normalization has genus
below9, supplies a point outside `M tensor Q` over its function field.
If the curve and map are over `Q`, this certifies a rank-at-least18 subgroup
over `Q(C)`. No new height or finite-group rank calculation is needed for
that function-field conclusion. The map must be birational onto the actual
multisection: merely parametrizing an inherited section with redundant
degree does not qualify.

**Boundary.** This is not a guarantee at every rational specialization.
It does not select a302 member, prove a rational point source on an arbitrary
cover, or distinguish the eight null fibres. The known dependent conic
specialization is entirely compatible with this theorem. The central
remaining task is rational split incidence with specialization independence,
not generic independence of a genuine bisection.

## 1. Audit: what this adds

**Verified prior applications.** The
[trace-parity theorem](DET1092_TRACE_PARITY_DESCENT_2026-09-09.md)
already proves the genus-nine halving cover and2-saturation under genus
below9. The existing
[multisection-height theorem](DET1092_TWO_FIBRATION_SEED_CONSTRUCTION_2026-09-08.md)
proves independence for smooth multisections unramified over singular fibres.
The explicit rational, genus-one and genus-two covers already have their
own rank18 certificates.

**New scope.** The proof below extends division exclusion to every prime,
excludes all torsion over low-genus bases, and removes the old
singular-fibre-unramified hypothesis in the stated degree/genus ranges.
In particular it applies to the normalization of a singular multisection,
provided the normalization map and irreducibility are verified. None of the
old concrete rank18 constructions is presented as newly discovered.

## 2. Simple j-poles exclude all prime-degree isogenies

**Verified equation-side input.** On the global short model
`y^2=x^3+A(t)x+B(t)`, the checker verifies `deg(A)=8`, `deg(B)=12`,
squarefree discriminant of degree24, and `gcd(c4,Delta)=1`. The scaled
infinity fibre is smooth. Therefore the parent's nonconstant `j`-map has
exactly24 poles, all simple.

**Established literature.** `X0(ell)` classifies an elliptic curve with
a cyclic subgroup of order `ell`. Cusp widths give the pole orders of
its forgetful `j`-map. See
[Milne, *Modular Functions and Modular Forms*, sections2 and8](https://www.jmilne.org/math/CourseNotes/MF.pdf).

**New application; proof.** A Galois-stable line in `E[ell]` over
`Qbar(t)` would give a nonconstant map `P1 -> X0(ell)` factoring the
parent's `j`-map. At cusp0 the width is `ell`: conjugating
`[[1,m],[0,1]]` by `[[0,-1],[1,0]]` gives
`[[1,0],[-m,1]]`, which lies in `Gamma0(ell)` exactly when `ell | m`.
Hence `j` has a pole of order `ell` there. A nonconstant map of proper
curves is surjective, so its pullback would have a pole of order at least
`ell`. This contradicts simple poles. Thus no such line, and no
prime-degree isogeny, exists over `Qbar(t)`.

This argument is uniform in `ell`; it is not a finite modular-polynomial
test. It remains valid after every constant number-field extension.

## 3. The geometric mod-ell monodromy is full SL2

**Established local geometry.** An `I1` fibre acts on `E[ell]` by a
nontrivial transvection, conjugate to `[[1,1],[0,1]]`. This is the usual
nodal/Tate-curve monodromy; the semistable setup is reviewed in
[Schuett--Shioda, *Elliptic Surfaces*](https://arxiv.org/abs/0907.0298).
The Weil pairing puts geometric monodromy in `SL2(F_ell)`.

**New elementary application.** By section2 its action is irreducible.
Conjugate one transvection by an element moving its fixed line. In the
basis of the two distinct fixed lines these give a nontrivial upper and
a nontrivial lower unipotent. Their powers give every upper and lower
unipotent, because the coefficient field is the prime field `F_ell`.
These generate `SL2(F_ell)` by elementary elimination. Therefore

\[
\rho_{E,\ell}(G_{\overline{\mathbf Q}(t)})
       =\operatorname{SL}_2(\mathbf F_\ell)
\quad\text{for every prime }\ell.
\tag{3}
\]

The checker verifies the symbolic unipotent/Bruhat identities, not just
a few group orders. No assertion about the full adelic image is needed.

## 4. Primitive division torsors are connected

**New deduction using standard Kummer torsors.** The Galois group acting
on the `ell^2` solutions to `[ell]Q=Z` is an affine subgroup
`G <= V semidirect SL2(F_ell)`, where `V=F_ell^2`, projecting onto the
full linear group. Its translation kernel `G intersection V` is an
invariant vector subspace. By irreducibility it is zero or all of `V`.

For odd `ell`, if the kernel were zero, `G` would be the graph of a
1-cocycle `c` on `SL2(F_ell)`. The central element `z=-I` gives

\[
c(g)+gc(z)=c(z)+zc(g),\qquad
2c(g)=(I-g)c(z).
\]

Thus `c(g)=(I-g)v`, `v=c(z)/2`, and every affine transformation fixes
`v`. There would be a division point over `Qbar(t)`, contrary to
`Z notin ell*M`. Therefore the translation kernel is full; the action is
transitive on all `ell^2` points.

For `ell=2`, the previous halving argument is elementary: an orbit of
size1 is a generic half; an orbit `{P,P'}` of size2 gives the nonzero
rational2-torsion point `P+P'-Z=P'-P`. Both are impossible. Every proper
partition of4 has a part of size1 or2. Hence the four-point torsor is
connected. This handles the characteristic-two coefficient case without
dividing by2 in the cocycle argument.

For nonzero torsion, `SL2(F_ell)` is transitive on `V minus {0}`.
Its coordinate curve is therefore connected of degree `ell^2-1`.

## 5. Exact genera and sharpness

**New all-prime application of the local halving argument.** At each
nodal fibre, the section `Z` lies in the smooth Neron locus. Over a
strictly henselian characteristic-zero trait its reduction lies in `G_m`.
The `ell`-power map is surjective on the algebraically closed residue
field and etale, so an `ell`-division point lifts locally. Translating by
this division point identifies the local affine torsor action with the linear action
on `E[ell]`.

The transvection `(x,y) -> (x+y,y)` has `ell` fixed points and `ell-1`
cycles of length `ell`. Its contribution to ramification is therefore
`(ell-1)^2`. There is no ramification over smooth elliptic fibres,
including infinity, because multiplication is etale there. At the24
nodal fibres the total contribution is `24*(ell-1)^2`.

**Established formula, new application.**
[Riemann--Hurwitz](https://stacks.math.columbia.edu/tag/0C1B) gives

\[
g(T_{\ell,Z})=1-\ell^2+12(\ell-1)^2
             =(\ell-1)(11\ell-13).
\]

Removing the zero torsion point removes one unramified sheet and leaves
the same local ramification contribution; its degree is one smaller,
so the nonzero-torsion genus is one larger.

**Sharp bound.** The minimum division genus is9, since
`g(T_(ell,Z))-9=(ell-2)*(11ell-2)`. Genus9 really permits dependent
multisections: for any primitive displayed section `Z=(cx,cy)`, the
degree4 curve `[2]Q=Z` has genus9 and its universal point is a half of
an inherited point. One explicit affine model is

\[
m^4-6c_xm^2-8c_ym-3c_x^2-4A=0,
\quad x_Q=(m^2-c_x)/2,\quad y_Q=m(x_Q-c_x)-c_y.
\]

Take `Z` to be the first displayed generic basis section, for example.
The universal polynomial and doubling identities were already checked in
the trace-parity packet. This sharpness example uses no exceptional point.
It is ramified above nodal fibres, so it does not satisfy the older
unramified multisection-height theorem's hypotheses.

## 6. Saturation and generic rank18 without a height calculation

**Proof of (2).** If a base curve has genus below9, it cannot map
nonconstantly to a primitive division curve or a nonzero-torsion curve:
their genera are at least9 and10, respectively. This follows directly
from Riemann--Hurwitz. Thus no new prime-order torsion is possible.

Suppose `nQ in M`. For a prime `ell | n`, put `R=(n/ell)Q` and `Z=nQ`.
If `Z notin ell*M`, the coordinates of `R` induce a map to `T_(ell,Z)`,
impossible. Otherwise write `Z=ell*S` in `M`. Then `R-S` is
`ell`-torsion and therefore zero. Replace `n` by `n/ell` and repeat.
Eventually `Q in M`. This treats all composite denominators, not only
primitive points or squareclasses.

For a degree2 base change, the same proof works by degrees: a map over
`P1` to a primitive division curve would require `ell^2 | 2`, and a
nonzero-torsion point would require `ell^2-1 | 2`. Neither is possible.
More generally these exclusions work for any degree `d` divisible by
neither `ell^2` nor `ell^2-1`, for every prime `ell`.

For the normalization `C` of a genuine multisection of degree greater
than1, its tautological point cannot already lie in `E(Qbar(t))`:
otherwise the image would be a section. Full saturation then implies
that point is outside `M tensor Q`. If the data are rational, the17
rational inherited sections and this rational function-field point give
the asserted rank18 subgroup over `Q(C)`.

This also shows why a generic trace-parity class is not itself an extra
direction. Independence is attached to a genuine irreducible multisection
and its marked generic point, not merely to the inherited trace label.

## 7. What the controls now tell us

**Verified prior applications; no new specialization.**

- The orbit8044 conic, the norm8 genus-one carrier and the historical
  genus-two carriers are genuine bisections. Their function-field points
  are independent by the same theorem, without separate ramification gates.
- The calibrated302 member and the three productive conic specializations
  have independently certified new rational directions. Those pointwise
  certificates are still necessary for the specified parameters.
- The same conic at `-528/3635` specializes to inherited points, with exact
  halving cycles. Generic independence does not prevent this finite-height
  specialization failure.
- The nine blinded recoveries are degree1 generic sections when referred
  to the full MW17 group. They do not meet the genuine-degree-greater-than1
  hypothesis; independence over an MW16 core is a different assertion.
- Nonsplitting on any null-control parameter is an incidence failure of
  that selected cover, not failure of its generic independence theorem.

**Open endpoint.** A generic-input procedure selecting a useful member
through302 remains missing. This result removes generic dependence and
hidden division as explanations for failure of genuine bisection
constructions. It does not remove the rational splitting or specialized
independence requirements, and it supplies no new rational point by itself.

## 8. Checkpoint and proof status

**Verified application.** The checker completed in0.59s under a25s cap.
It verifies the24 simple `j`-poles, smooth infinity, pinned full geometric
MW17 certificate, cusp matrix, unipotent identities, central cocycle
identity and exact genus polynomial. Tiny permutation regressions at
`ell=2,3,5,7,11` check the local cycle counts; they are not the proof for
all primes.

The all-prime group/curve argument is a written mathematical proof with
symbolic and source checks, not formal verification or an independently
duplicated proof. No prime-degree division polynomials were generated,
and no point search, new parameter, multisection enumeration, class group,
unit group, generic Selmer dimension, production change or detached job
was used.

```bash
timeout 25s /home/royvanrijn/.local/bin/sage -python research/elliptic-curves/cas/verify_det1092_prime_division_barrier.sage
```

- [Protocol](../../artifacts/generated-results/elliptic-curves/det1092_prime_division_barrier_v1/protocol.json).
- [Symbolic and pinned-input replay](../../artifacts/generated-results/elliptic-curves/det1092_prime_division_barrier_v1/replay.json).
- [Existing38 exact split classifications](DET1092_SPLIT_SEED_HALVING_DICHOTOMY_2026-09-08.md).
