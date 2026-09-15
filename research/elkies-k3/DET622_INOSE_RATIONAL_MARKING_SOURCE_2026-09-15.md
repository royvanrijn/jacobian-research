# Determinant 622: a rational marking via a quadratic Q-curve

The subsequent [global root obstruction](DET622_GLOBAL_ROOT_OBSTRUCTION_2026-09-15.md)
excludes every MW17 fibration on this NS. The rational-marking existence
proved here remains valid. Its explicit twist and section remain uncomputed,
but are no longer prerequisites for the different-NS MW17 objective.
The UNKNOWN geometric endpoint below records this arithmetic gate alone.

There exists a projective K3 over Q with full rational saturated NS

    U + E8(-1) + E8(-1) + <-622>.

It arises from the published exceptional point on X0+(311), through an
explicit rational Inose equation followed by a constant quadratic twist.
This is a new source for this research queue, outside its retained827-row
catalogue; no literature-newness claim is made. The twisting squareclass
has not yet been computed. Rootless-frame existence and an MW17 equation
remain UNKNOWN.

Unlike a bare rational point on a period quotient, the argument below
constructs a rational surface and proves that a single constant twist
makes **all** of its geometric divisor classes rational. It does not
assume an arithmetic identification from the complex period group alone.

## Exact non-CM input and rational equation

[Adzaga et al., Theorem1.1 and Section5.1](https://arxiv.org/pdf/2105.04811)
identify the point [6:8:-1:-2] on their canonical model of X0+(311) as
noncuspidal and non-CM. It corresponds to a cyclic311-isogeny between
elliptic curves with conjugate j-invariants A+B sqrt(D), where

```text
A = 31244183594433270730990985793058589729152601677824000000
B = 1565810538998051715397339689492195035077551267840000
D = 39816211853 = 11*17*9011*23629.
```

These are published arithmetic theorem inputs, not a new isogeny
calculation. The checker verifies that the displayed point lies smoothly
on the two canonical equations. The retrieved paper and its hash are
retained. It does not rerun quadratic Chabauty or evaluate a degree311
modular polynomial.

Set

    I = (A²-B²D)/1728²,
    J = ((A-1728)²-B²D)/1728².

The following equation has rational coefficients:

    X0: y² = x³ - 3IJ t⁴ x + t⁵(t² - 2IJ² t + I²J³).

Over Qbar this is the Inose surface of the displayed pair. Indeed, for
an equation y²=x³+a t⁴x+t⁵(b t²+c t+d), its normalized Inose invariants
are -a³/(27bd) and c²/(4bd). Here they are exactly I and J, the symmetric
j-invariants required by the standard Inose normal form; see
[Kumar–Kuwata, Section2.1](https://arxiv.org/pdf/1409.2931).
No rational cube root of I or square root of J is required.

The checker verifies the normalization identities and that the residual
quartic discriminant is squarefree. The minimal elliptic K3 model has
II* fibres at0 and infinity and four geometric I1 fibres. In particular
both reducible fibres are at rational base points.

## Full geometric NS, not a section subgroup

The curves are non-CM and nonisomorphic. Their geometric Hom group has
rank one. The cyclic isogeny of prime degree311 is primitive: writing it
as m times another isogeny would make m² divide311, so m=1.

Shioda's full Inose Mordell–Weil lattice theorem, stated in
[Kumar–Kuwata, Proposition3.1](https://arxiv.org/pdf/1409.2931) and
[Utsumi, Proposition3.1](https://arxiv.org/pdf/2209.02463), identifies this
torsion-free group with Hom(E1,E2), with heights scaled so that the
primitive311-isogeny section has height622. Thus geometric MW rank is
one and geometric Picard rank is19.

The trivial lattice U+E8(-1)+E8(-1) is unimodular. Shioda–Tate therefore
gives the full NS displayed above, not merely a finite-index section
sublattice. Its determinant is622 and its transcendental complement is
U+<622>. The checker records both full Gram matrices.

## The rationalizing twist

The zero section and fibre are defined over Q. Every component of each
II* fibre is defined over Q: Galois preserves its incidence graph and
multiplicities, and the affine E8 graph has no nontrivial automorphism.
Thus the entire rank18 trivial lattice is represented by rational divisors.

Let P generate the full geometric Mordell–Weil group, which is Z and has
no torsion. Galois acts through a character

    chi: Gal(Qbar/Q) -> {+1,-1},   sigma(P)=chi(sigma)P.

This character is continuous: P uses finitely many algebraic coefficients.
It is consequently the character of a constant quadratic extension
Q(sqrt(delta)), including the trivial squareclass delta=1. Twist the
elliptic equation by that character:

    X_delta: delta*y² = x³ - 3IJ t⁴ x + t⁵(t² - 2IJ² t + I²J³).

Under the geometric isomorphism to X0, its new Galois action on P is the
product chi²=1. Hence P is a rational section of X_delta. The fibre
configuration and rational E8 components are unchanged, so the primitive
section together with the rational trivial lattice generates every
geometric NS class over Q. This proves full rational saturated marking,
including actual divisor representatives, without a residual Brauer issue.

The proof determines **existence** of the rationalizing squareclass. It
does not yet identify delta, display P, or claim that the untwisted X0
already has rational MW rank one. The explicit equation is the untwisted
model; a fully marked explicit equation requires the remaining twist step.

## Next obligations

Identify the twisting squareclass and a certified rational generator, then
seek a primitive rootless U in this NS. The determinant388 global
obstruction does not apply to determinant622. Its auxiliary root-energy
method can be tried after the full source gate, with literal discriminant
forms rather than determinant guesses. An MW17 endpoint still requires
its nef fibration, equation and17 saturated rational sections.

```sh
sage -python research/elkies-k3/scripts/certify_det622_inose_source.py --check
```

The [certificate](../artifacts/generated-results/elkies-k3-det622-inose-source-v1/certificate.json)
separates exact model arithmetic, published theorem inputs and the written
twisting proof. The twist value, primitive section, rootless frame and
MW17 equation are explicitly uncomputed. Independent implementation,
formal verification, external review and public novelty are unclaimed.
