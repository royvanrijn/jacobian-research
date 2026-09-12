# Close the carrier bank; test one complete arithmetic dependency

The fixed-word continuation branch is now closed by the stronger
[finite-Selmer-specialization theorem](FIXED_WORD_HAS_FINITE_SELMER_SPECIALIZATIONS_2026-09-12.md).
The present note retains the carrier diagnostic, translation/pairing lemmas,
explicit continuation and earlier genus bound used in that proof.

The low-shell carrier experiment is closed at its original stopping condition.
Its 182 completed bisections represent **182 distinct quadratic extensions of
the fixed parameter line**, all nonsplit at `t0=3/17`. No class-specific label
test was reached. Neither translation by a generic section nor pairing can
repair these particular incidence failures. The incomplete shells, 73 split
residuals and interrupted candidate remain separate and unchanged.

The first arithmetic compatibility check also has a precise endpoint. The
complete 1,676-atom column-6 word has an explicit coefficientwise continuation
in the MW16-05 cubic algebra. It specializes to the certified word at `3/17`,
but has odd valuations at twelve distinct geometric good-fibre parameters.
Consequently it is not a rational-point Kummer class over `Q(t)`. Moreover,
any finite parameter cover making **this particular continued word** a
rational-point class has genus at least five. This obstructs one specified
continuation, not the soluble class at the control or alternative formulas.

## The final nonsplitting table

The [CSV table](../../artifacts/generated-results/elliptic-curves/marked_carrier_closure_v1/nonsplitting.csv)
contains every exact branch polynomial, its constant, its value at `3/17`,
the homogeneous value and a small local witness. The
[certificate](../../artifacts/generated-results/elliptic-curves/marked_carrier_closure_v1/result.json)
also records the complete local squareclass at every prime at most 1009.
Only the [previously completed maps](../../artifacts/generated-results/elliptic-curves/marked_two_class_propagation_v1/carriers.json)
are read. No Riemann–Roch, enumeration, point search or interrupted-candidate
calculation is repeated.

| Diagnostic | Exact result |
|---|---:|
| Completed branch polynomials, all degree two | 182 |
| Distinct monic squarefree branch polynomials | 182 |
| Distinct extensions over the fixed `t`-line | 182 |
| Distinct rational squareclasses of `d(3/17)` | 182 |
| Negative at the real place | 21 |
| Nonsquare in `Q_2`, `Q_3`, `Q_5`, `Q_17` | 15, 20, 14, 23 |
| Largest single-prime exclusion among primes at most 1009 | 109, at 739 |
| One tested place excluding the whole bank | None |

The local counts overlap. Distinct squareclasses are not asserted to be
linearly independent. No claim is made about a common obstruction at an
untested larger prime.

For a squarefree degree-one or degree-two polynomial `d=c*m`, with `m` monic,
the extension is determined by `m` and the class of `c` in `Q*/Q*^2`.
Two such polynomials have the same extension precisely when their monic
polynomials agree and their constant ratio is a rational square. This is the
degree-two specialization of the repository's
[existing extension protocol](../../elkies-k3/scripts/hash_bisection_extensions.py).
Exact integer square roots suffice for constant-ratio comparisons; large
integer factorization is unnecessary. Here even the monic polynomials differ.
There are no independently chosen changes of the base coordinate.

For each nonzero rational value `q`, the local test strips numerator and
denominator powers of `p` exactly. At odd `p`, a square requires an even
valuation and a square unit modulo `p`; at 2 it requires an even valuation
and unit 1 modulo 8. Negative valuations are retained. The homogeneous model
is `S^2=D(T,U)=sum_i c_i*T^i*U^(2-i)`. At `(T,U)=(3,17)`,
`D(3,17)=17^2*q`; the reciprocal chart has value `D(3,17)/9` at `U/T=17/3`.
Both chart values give exactly the same local square decisions. Thus the
denominator prime 17 is handled without an invalid affine reduction.

## Translation and pairing cannot supply the missing control point

Let `C` be a normalized bisection and `T` a rational section. Translation
by `T` is invertible over the generic fibre and fixes `t`, so it identifies
the function fields of `C` and its image over `Q(t)`. The induced isomorphism
of smooth projective normalizations commutes with their finite maps to
`P1`. It therefore identifies the fibres over `3/17`, including their
rational points. Nonsplitting survives **every** generic-section translation,
even if an unnormalized image has a different equation or divisor class.

A rational point on the normalized component of a fibre product projects
to rational points of both constituent normalized covers. If either has
no rational point above `3/17`, that component cannot have the required
rational point there. This needs no independence or class-label test.

These deductions exclude translates and pairings of the completed misses.
They do not complete either shell, settle candidate 171, or constrain a
different carrier. In particular, soluble additional strict classes do not
imply representation on low-intersection rational bisections. The
[conditional propagation theorem](MARKED_TWO_CLASS_PROPAGATION_2026-09-12.md)
remains valid; the sampled existence hypothesis received no positive support.

## A defined continuation of the full column-6 word

Use the frozen compact model `Y^2=X^3+A(t)*X+B(t)`, with exact coefficient
arrays and all sixteen sections in the
[geometry packet](../../artifacts/generated-results/elliptic-curves/marked_two_class_propagation_v1/geometry.json).
Set `u=289/2` and define

```
c1(t) = (u^4*A(t)+27)/81,
c0(t) = (u^6*B(t)+3*u^4*A(t)+27)/729,
f_t(Z) = Z^3+Z^2+c1(t)*Z+c0(t),
K_t = Q(t)[theta_t]/(f_t).
```

The exact coordinate transport is
`z=(u^2*X-3)/9`, `w=u^3*Y/27`, giving `w^2=f_t(z)`.
At the control these are `z=4*x_original`,
`w=8*y_original+4*x_original`. The specialized cubic is exactly

```
Z^3+Z^2-2919231625641258502793755607986240*Z
 +45440201616242830029801770634418828098464819545088.
```

It is irreducible modulo 23. Monicity then proves that `f_t` is irreducible
over `Q(t)`, since a factorization would specialize to one at the control.
The [definition packet](../../artifacts/generated-results/elliptic-curves/frozen_dependency_continuation_v1/definition.json)
contains all expanded coefficients, every atom and the complete factor list.

Keep every frozen `alpha_i=a_i+b_i*theta` in the first certified dependency,
including its full 1,676 indices, and define

```
alpha_i(t) = a_i+b_i*theta_t,
n_i(t) = Norm(alpha_i(t))
       = a_i^3-a_i^2*b_i+c1(t)*a_i*b_i^2-c0(t)*b_i^3,
pi_i(t) = n_i(t)*alpha_i(t).
```

For the transported generic sections `(z_j,w_j)`, define
`gamma_j(t)=s_j*(z_j(t)-theta_t)`, where each positive rational square `s_j`
is determined by the original generic representative. All sixteen formulas
and their specialization identities are checked. Thus the specified word is

```
beta_6(t) = product_(i in I6) pi_i(t)
            * gamma_3(t)*gamma_5(t)*gamma_6(t)
            * gamma_7(t)*gamma_9(t)*gamma_10(t).
```

Its norm is a square identically: `Norm(pi_i)=n_i^4` and
`Norm(gamma_j)=s_j^3*w_j^2`. Every atom norm and generic representative
specializes exactly to the original certificate. All these factors are
units in the algebra over `Q[t]_(t-3/17)`, which supplies the specialization
homomorphism `theta_t -> theta`. This is not an evaluation homomorphism on
the entire rational function field, and it identifies no number-field
prime ideals across parameters. The existing square-equivalence circuit
to compact column 6 is retained as inherited provenance, not recomputed.
Column 7's complete 1,572-atom dependency remains frozen and untested here.

## Exact obstruction to this continuation

Choose atom 2 by the frozen rule “least index in the complete word”:

```
alpha_2(t) = -6324860115183736256 + 286*theta_t.
N(t) = -54/10793861 * n_2(t).
```

Here `N` is a primitive integer polynomial of degree 12. No factorization
or search for another atom is required. At the single predeclared prime
`p=1000003`, the [certificate](../../artifacts/generated-results/elliptic-curves/frozen_dependency_continuation_v1/result.json)
supplies a Bezout identity proving

```
gcd(N, N' * disc_Z(f_t) * product_(i in I6, i != 2) n_i
                       * product_(j=0..15) w_j) = 1  in F_p[t].
```

All rational coefficient denominators are units at `p`, and reduction
preserves the degree of `N`. Therefore the same coprimalities hold in
`Q[t]`. In particular `N` is squarefree, all its twelve geometric zeros
are good-fibre parameters, and the other 1,675 projected atoms and all
sixteen generic representatives are units above them.

At any irreducible factor `q` of `N`, the cubic is etale. Its reduction has
the simple root `r=-a_2/b_2`, and two other geometric roots. The valuations
of `alpha_2` on those three sheets are `(1,0,0)`: its norm has a simple zero
and only the first factor vanishes. The base scalar `n_2` has valuation one
on each sheet. Hence the valuations of `pi_2`, and of the complete word
`beta_6(t)`, are

```
(2,1,1).
```

If the residual quadratic is irreducible, its single degree-two prime has
valuation one; the displayed pattern is over geometric sheets. The other
atoms were checked in full, so atom 2 is an obstruction witness, not a
substitute for the original dependency.

For completeness, rational-point Kummer classes have even valuations at
these good divisors. After an unramified splitting of the monic cubic,
if a point has integral `z`, at most one factor `z-theta_k` has positive
valuation, and their product is `w^2`; that valuation is even. If `z` has
negative valuation, the term `z^3` determines `f_t(z)`, so
`3*v(z)=2*v(w)` forces `v(z)` even; every `z-theta_k` has that same valuation.
The identity point has trivial Kummer class. If a base change introduces
rational 2-torsion, its Kummer entries are root differences and the
derivative at the matching root; all are units because the discriminant
is a unit. This proves the needed
horizontal, residue-characteristic-zero assertion directly. The standard
number-field unramified descent framework is discussed in
[Poonen–Schaefer, Sections 12.2–12.5](https://math.mit.edu/~poonen/papers/descent.pdf).

Thus the odd valuations prohibit a rational-point Kummer interpretation of
this continued word over `Q(t)`. Squares and generic Kummer corrections
cannot change the obstruction. A norm-square identity alone did not
continue the original even-ideal-valuation dependency.

There is also a useful exact base-change consequence. Let a smooth
geometrically integral curve `B` map to the parameter line with degree `d`.
If the pullback of this word becomes a rational-point Kummer class, each
ramification index above each of the twelve zeros of `N` must be even:
an odd index would retain an odd valuation on the two residual sheets.
Each target point therefore contributes at least `d/2` to ramification.
Also `d` is even and at least two. Riemann–Hurwitz gives

```
2*g(B)-2 >= -2*d + 12*d/2 = 4*d,
g(B) >= 2*d+1 >= 5.
```

In particular no rational or genus-one parameter base repairs this fixed
coefficientwise continuation. This is a necessary condition, not a promise
of solubility on any higher-genus cover. No such cover is searched.
Other point-producing recipes can specialize to the same marked class
while defining different classes over `Q(t)`; they are not excluded.
This is not a genus bound for a carrier whose point map merely specializes
to column 6.

## Replay, provenance and next gate

The independent replay uses SymPy and imports no producer polynomial code.
It verifies the field transport, all sixteen generic sections, all 1,676
atom specializations, an independent resultant for atom 2's norm, the full
modular product, and the supplied Bezout identity. It passes in 0.682 CPU
seconds; the plain-Python producer used 0.204 CPU seconds. These are component
times, not an end-to-end research benchmark. Both are below 162 MB peak RSS.
An initial Sage library mapping failure under a 1 GiB address-space limit
occurred before the modular test; its source and failure are retained.
The Python implementation uses the same atom, prime, word and memory bound.
There is no alternate-witness attempt. Input hashes are checked; no claim
of operating-system blinding is made.

```sh
python3 research/elliptic-curves/rank-jump/diagnose_marked_carrier_nonsplitting.py --check
sage -python research/elliptic-curves/rank-jump/verify_frozen_dependency_continuation.py
```

Both commands read the retained evidence without extending a search or
overwriting receipts. The complete original two-class arithmetic and
compaction provenance remains in the earlier archive.

The research gate returns to the
[successful adaptive principal-relation constructor](TWO_CONSTRUCTED_STRICT_CLASSES_AND_302.md):
reuse the equation-and-generic-data procedure to construct a fresh complete
dependency for each future input. The
[finiteness theorem](FIXED_WORD_HAS_FINITE_SELMER_SPECIALIZATIONS_2026-09-12.md)
closes this fixed output as an infinite template; it does not impose that
word on a parameter-dependent constructor. Every new output still needs
its full provenance and separate solubility, independence and ideal audits.
No new collector, carrier shell or completion of candidate 171 is scheduled.

The desired infinite marked rank-18 construction is still open. Finite
congruence refinements can preserve a supplied independence certificate,
but strictness and ordinary ideal-class independence at new bad primes
remain separate endpoints requiring their own certificates or uniform proof.
