# One native quadratic cover cannot supply the whole observed +8 jump

For the fixed native cover `orbit-1795d`, exact Frobenius moments and the
functional-equation sign prove
\[
\boxed{1\le\operatorname{rank}E^{(q)}(\mathbf Q(t))\le7},\qquad
q(t)=409689-1439214t+328441t^2.
\]
This bounds **all** arithmetic generic sections of this twist, including
undisplayed sections of arbitrarily large height. It is not a bounded
section-search conclusion. The earlier bound 20 was geometric Hodge
capacity; the new bound 7 is arithmetic and must not be extended to all
algebraic constant fields.

On the observed +8 fibre `08234-009`, at \(t_0=-4112/1937\),
\[
q(t_0)=(4307419/1937)^2.
\]
The original family has exact generic rank 17. Character decomposition
therefore bounds the **full** arithmetic generic group over this quadratic
cover between 18 and 24. The retained fibre subgroup has rank 25. Thus
**at least one of its directions lies outside the rational span of the
specialized full generic pullback group**, at either rational lift of t0.

This excludes the hypothesis that hidden sections on this particular cover
explain the entire observed +8 jump. It does not determine which retained
point lies outside that span, the exact twist rank, or the original fibre's
full rank. Other cover directions and additional specialization events
remain possible.

## Fixed experiment and good reduction

The [moment protocol](NATIVE_TWIST_FROBENIUS_MOMENTS_PROTOCOL.json) fixes the
smallest-coefficient member of the previously selected positive quartet.
This is retrospective selection; no replacement twist or new rational
parameter is tested. The curve is
\[
E^{(q)}:\quad Y^2=X^3+A(t)q(t)^2X+B(t)q(t)^3,
\]
with the exact published-R17 A,B retained in the frozen input.

The known section is rederived from the stored lift P by
\(T=P-\sigma(P)\). Writing \(y(T)=\sqrt q\,y_1(T)\), its twist coordinates
are \(X=qx(T)\), \(Y=q^2y_1(T)\). The exact polynomial identity is checked;
the degrees are at most 6 and 9. Its twist height is 6, from the canonical
height-12 pullback calculation in
[Theorem F2](../../elkies-k3/RANK_MUTATION_AND_LIFT_THEOREMS.md).
This supplies the non-torsion lower bound one.

At p=131, exact polynomial calculations verify:

- the original discriminant has degree 24 and is squarefree;
- q remains a separable degree-two polynomial, coprime to that discriminant;
- the fibre at infinity is smooth.

Since p>3, the original 24 nodal fibres and the twist's two additional
I0* fibres retain their tame, disjoint configurations. The corresponding
minimal elliptic surfaces have good surface reduction: the nodal total
spaces are smooth, and the I0* resolutions spread over the unramified
branch data. The geometric trivial lattices have unchanged ranks 2 and 10.

Specialization of divisor classes into the smooth reduction preserves
intersection pairings and is injective. The rational Mordell–Weil space is
the complement of the trivial lattice; its image remains in that complement.
Classes defined over Q specialize to Frobenius-fixed divisor classes after
Tate twist. Hence the arithmetic generic rank is at most the multiplicity
of the eigenvalue p in the nontrivial surface cohomology. Nonintegral
coefficients of an individual section do not undermine this divisor-class
argument, which does not assume good affine reduction of every point.

## Exact counts, not prime averages

The conductor degrees are 24 for the original and 28 for the twist.
Over the rational parameter line their nontrivial L-polynomial degrees
are consequently 20 and 24. The inverse roots have complex absolute value
p, and rank over the finite function field is at most the central
multiplicity; neither equality in BSD nor equality in the Tate conjecture
is needed. These facts are stated in
[Ulmer, Lecture 1, Theorems 9.3 and 12.1](https://arxiv.org/pdf/1101.1939).

For each n=1,2, all fibres over F_(p^n), including infinity, are counted.
The trace on the nontrivial cohomology is
\[
T_n=-\sum_{t\in\mathbf P^1(\mathbf F_{p^n})}a_{p^n}(t).
\]
For smooth fibres the usual point count gives a. Nodal fibres contribute
their exact local Euler trace +1 or -1; they are not discarded. At a
branch of q the twisted local trace is zero. Elsewhere the twist trace
is the original trace multiplied by the quadratic character of q(t).
The smooth infinity fibre is treated with the correct leading coefficients.

| Surface | L-polynomial degree N | T1 | T2 |
|---|---:|---:|---:|
| Original R17 | 20 | 1884 | 319520 |
| Native twist 1795d | 24 | 122 | 33710 |

One C++ worker performs the count, with a 60-second bound. It uses
\(\mathbf F_{131^2}=\mathbf F_{131}[\alpha]/(\alpha^2-2)\) and weights
conjugate parameter pairs. A separate NumPy/Sage implementation independently
recounts **all 8,779 fibre orbits**, building its character table from actual
finite-field squares. There are also 148 independent Sage elliptic-trace
cross-checks. Both complete recounts agree exactly.

## A two-moment upper bound

Write the normalized inverse roots as \(z_j=e^{i\theta_j}\), and let m
be the multiplicity of z=1. Put
\[
S_1=T_1/p,\quad S_2=T_2/p^2,\quad x_j=\cos\theta_j.
\]
Conjugate pairing gives
\(\sum x_j=S_1\) and \(\sum x_j^2=(N+S_2)/2\).
For any rational c other than one,
\[
m(1-c)^2\le\sum_j(x_j-c)^2
=\frac{N+S_2}{2}-2cS_1+Nc^2.
\]
Every omitted summand is nonnegative. This proves an upper bound from just
two **exact** moments; a full Frobenius polynomial is unnecessary.

Taking c=-1/2 already gives
\[
m\le\frac{3N+4S_1+2S_2}{9}.
\]
For the twist this is \(1366940/154449<9\), so m<=8. Optimizing c gives
\(c=-206805/395882\) and the slightly sharper rational bound
\(5332004/602687\), with the same integer consequence. For the original
surface, the optimized bound is \(1538972/90491<18\). Its known rank-17
group makes this an exact rank-17 calibration over F_131(t), as well as
recovering the original arithmetic generic rank.

## The functional equation improves eight to seven

The [separate sign protocol](NATIVE_TWIST_MOMENT_PARITY_PROTOCOL.json)
uses the same fixed surface and performs no additional point counts.
The original discriminant factors over F_131 into irreducibles of degrees
1,2,3,4,7,7. At each nodal place the split character is computed from -c6;
the local root number is minus this character. The computation uses
residue-field norms of 864B and 864Bq^3, and an independent verifier checks
the characters directly in each residue field.

| Factor degree | Original local sign | Twisted local sign |
|---:|---:|---:|
| 1 | -1 | +1 |
| 2 | +1 | +1 |
| 3 | -1 | +1 |
| 4 | +1 | -1 |
| 7, first factor | +1 | -1 |
| 7, second factor | -1 | -1 |

The I0* places are tame ramified quadratic twists of good reduction.
For such a two-dimensional unramified representation, the square of the
quadratic Gauss-sum sign gives the local factor chi(-1). The residue
degrees of the branch places sum to two, so their product is
\(\chi_{131}(-1)^2=1\). Infinity contributes +1. This is the tame
local epsilon-factor calculation underlying the familiar I0* formula;
compare [Desjardins, Proposition 3.1](https://jtnb.centre-mersenne.org/item/10.5802/jtnb.1112.pdf).

Both global functional-equation signs are therefore -1. Their central
multiplicities must be odd. The twist's integer bound m<=8 improves to
m<=7, and arithmetic rank is at most m.

This uses parity of a known polynomial's central zero, **not** a parity
conjecture for the unknown Mordell–Weil rank. The conclusion is
rank<=7; it does not assert that the actual rank itself is odd.

## The remaining specialization direction is unavoidable

By [the exact character decomposition, Theorem F4](../../elkies-k3/RANK_MUTATION_AND_LIFT_THEOREMS.md),
\[
\operatorname{rank}E(\mathbf Q(t,\sqrt q))
=17+\operatorname{rank}E^{(q)}(\mathbf Q(t))\in[18,24].
\]
The retained nonzero square at t0 gives two rational places of this cover.
At either place, the specialization homomorphism's image has rank at most
24. The independent witness subgroup W from
[the frozen +8 certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_solubility_first_v1.json)
has rank 25. If H is the rational span of the full specialized generic
pullback group, then
\[
\dim\bigl((W\otimes\mathbf Q)/(H\cap(W\otimes\mathbf Q))\bigr)\ge1.
\]
No injectivity of specialization in the t-direction is assumed for this
dimension argument. Losing generic directions could only increase the gap.

The curve's rank is not bounded above by 24: its known rank-25 subgroup
already demonstrates why a generic bound cannot be read as a fibre bound.
The certified statement concerns the source of its directions, not their
visibility or the full rank of the specialization.

## Ranked consequences

1. **Incidence, proved restriction on the one-cover mechanism:** no more
   than seven generic arithmetic directions occur on this native twist.
   A full +8 block on the chosen fibre cannot consist solely of hidden
   generic sections on this cover. Large-height sections are covered by
   the bound too.
2. **Solubility, still necessary:** the successful square condition does
   make every rational generic section on the cover specialize rationally.
   It does not trivialize other parameter carriers or Selmer torsors.
   Their simultaneous global lifting remains the central unresolved event.
3. **Still possible:** additional directions from other singleton covers,
   product-character twists, or classes that arise only at specialization.
   The earlier four-cover construction proves a three-dimensional subblock,
   and this result does not close the uncovered part of the observed +8.
4. **Weak explanation:** a finite Nagao average or absence of a second
   polynomial section. Neither supplied this rank bound; exact moments,
   purity, good reduction and the functional-equation sign are its inputs.
5. **Missing computation:** an exact rank for this twist, or a corresponding
   bound on the relevant product-character twists. The cached base-fibre
   traces can be reused for a small predeclared retrospective comparison
   without recounting or selecting new original-family parameters.
6. **Visibility:** this changes no chart, score, worker or search budget.
   Agent 1 could eventually use such bounds to reject an impossible generic
   block explanation, but cannot use them to discard high-rank fibres.

Theorem F2 already proves injectivity of all 39,120 bisection squareclasses.
The preceding 39,119-equation common-field audit is an independent replay
on our frozen finite-chart subset, not a strengthening of that theorem.
The new ingredient here is the full arithmetic generic rank upper bound
for this fixed native twist.

## Reproducibility and scope

```sh
sage -python elliptic-curves/rank-jump/verify_native_twist_frobenius.py check
sage -python elliptic-curves/rank-jump/native_twist_moment_parity.py check
python3 elliptic-curves/rank-jump/native_twist_jump_boundary.py check
```

The first command repeats the complete independent finite-field recount
and the exact section/geometry verification. The second repeats the finite
factorization and local signs. The third checks the fixed cover square and
subgroup accounting against the frozen witness certificate. Inputs, source,
compiler invocation, finite-field ledger and proof dependencies are pinned.
No rational specialization or new rational point was searched.

- [Fixed equations and lift input](../../artifacts/generated-results/elliptic-curves/rank_jump_native_twist_frobenius_inputs_v1.json)
- [Exact moments, known section and complete orbit ledger](../../artifacts/generated-results/elliptic-curves/rank_jump_native_twist_frobenius_v1.json)
- [Functional-equation sign and refined bound](../../artifacts/generated-results/elliptic-curves/rank_jump_native_twist_moment_parity_v1.json)
- [Independent recount verification](../../artifacts/generated-results/elliptic-curves/rank_jump_native_twist_frobenius_verification_v2.json)
- [Full-pullback versus witnessed-fibre dimension certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_native_twist_jump_boundary_v1.json)

The initial successful recount is preserved as verification v1. Its
integer-versus-string JSON key mismatch was corrected in v2; the
`rank_jump_native_twist_verifier_format_v1.json` record retains the original
producer source. All arithmetic values and conclusions agree. Existing
search outputs and the mathematical-status registry were left untouched.
