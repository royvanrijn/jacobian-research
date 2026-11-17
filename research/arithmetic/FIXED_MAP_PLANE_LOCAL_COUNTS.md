# Local and Hasse counts on every fixed rational target plane

The two coordinate planes with rational sections are also the only
fixed rational affine planes with order `B^3` locally soluble target
counts. Every other fixed plane has at most order `B^2` locally
soluble targets. Hasse-failing targets have the latter upper bound
on every plane. The reason is that a transitive quintic root action
has no small orbit on its quadratic factors: each orbit has size
five or ten.

These are bounds on each fixed plane, with constants allowed to
depend on it. They do not improve the full ambient Hasse upper bound
by summing over moving planes.

## 1. The fixed map and the plane factor-selection cover

Retain the map `Phi`, target `y=(u,b,c)`, and `pi=-u`. For
`pi!=0` its inverse equation is
\[
 E_y(S)=\pi^5S^5-\frac32\pi^4S^4+\frac32\pi S^3
  +\left(-\frac54\pi-\frac9{32}b\right)S^2
  +\frac9{16}S-\frac9{32}c.
\]
Write `U={pi Disc(E_y)!=0}`. By FMA1--FMA3 in
[the ambient classification](FIXED_MAP_AMBIENT_LOCAL_CLASSIFICATION.md),
every Hasse-failing target lies in `U` and has an irreducible
quadratic factor times an irreducible cubic factor with matching
resolvent and the specified local tests.

Let `P` be a rational affine target plane different from `u=0`
and `c=0`. Theorem PS2 of
[the plane-section classification](FIXED_MAP_PLANE_SECTIONS.md)
proves that `E_y` is geometrically irreducible of degree five over
the function field of `P`. In particular `U_P=U intersection P`
is a dense open subset of the plane.

Here the generic irreducibility is an actual geometric assertion,
not an inference from specializations. PS2 proves it by the
primitive linear-parameter criterion. On a plane
`c=alpha pi+beta b+delta`, its inverse equation is
`N_pi(S)-(9/32)b(S^2+beta)`. Over `Qbar(pi)`, the two polynomials
`N_pi` and `S^2+beta` are coprime except on `c=0`: a nonzero
constant root `rho` of the latter leaves the nonzero leading term
`pi^5 rho^5` in the former, and `rho=0` gives precisely that
exception. The other plane forms `b=alpha pi+beta` and
`pi=a!=0` have a nonzero scalar coefficient of the free parameter
`c`, giving the same irreducibility. The remaining plane is `u=0`.

**Lemma PLC1 (degrees of the plane quadratic-factor cover).**
Over `U_P`, let `X_(2,P)` select a monic quadratic factor of the
monic polynomial `E_y/pi^5`. It is a finite etale cover of rank
ten. Each of its integral components has degree five or ten over
`U_P`; thus it has at most two such components. In particular no
such plane has a quadratic-times-cubic factorization over its
rational function field.

**Proof.** The ten geometric points of a generic fibre correspond
to the unordered two-element subsets of the five roots. Generic
irreducibility makes the arithmetic Galois group `G` transitive on
the roots. Hence `5` divides `|G|`, so Cauchy's theorem gives a
five-cycle in `G`. Its cyclic subgroup acts freely on the
two-element subsets: a nontrivial subgroup of a group of prime
order is the whole group, and a subset invariant under a five-cycle
can have only zero or five elements.

Every `G`-orbit on the two-element subsets is therefore a union
of orbits of size five. Since there are ten subsets, its size is
five or ten. The generic field degrees of the integral components
of the finite etale cover are these orbit sizes. The smooth
integral base is normal, so these generic components extend to
disjoint integral components of the cover, each dominating the
base. No degree-one orbit exists, which proves the last assertion.
\(\square\)

Only arithmetic transitivity is required for this lemma. Geometric
irreducibility from PS2 is stronger than needed; the argument does
not assert that every component of `X_(2,P)` is geometrically
integral.

## 2. The complete fixed-plane counting dichotomy

Use the actual target height `H(y)=H([1:u:b:c])`. For each
fixed plane `P`, let `A_P(B),G_P(B),L_P(B)` count respectively
Hasse-failing, rationally soluble and everywhere locally soluble
targets in `P(Q)` with height at most `B`. The exact partition is
`L_P(B)=G_P(B)+A_P(B)`.

**Theorem PLC2 (all fixed-plane local and Hasse bounds).** For
every fixed rational affine plane,
\[
                         A_P(B)\ll_P B^2.          \tag{2.1}
\]
More precisely:

- If `P` is `u=0` or `c=0`, then `A_P(B)=0`, and
  \[
   G_P(B)=L_P(B)=\frac4{\zeta(3)}B^3+O(B^2).
                                                               \tag{2.2}
  \]
- On every other fixed rational affine plane,
  \[
                     G_P(B)\ll_P B^2,\qquad
                     L_P(B)\ll_P B^2.              \tag{2.3}
  \]

Thus the two coordinate planes are exactly the rational affine
planes with order `B^3` locally soluble target counts. No fixed
rational affine plane has an order `B^3` Hasse-failing count.

**Proof.** On either coordinate plane, PS1 provides a rational
source point at every rational target, including its boundary.
The primitive projective-pair count in PS3 is
`4B^3/zeta(3)+O(B^2)`, proving (2.2).

On any other plane, every Hasse target belongs to `U_P` and has
a rational point of `X_(2,P)`: select its irreducible quadratic
factor. By PLC1 this point lies in one of at most two integral
components of degree five or ten. For each component apply the
following exact primary interface: Buggenhout, Cluckers, Salberger,
Santens and Vermeulen, *Serre's question on thin sets in projective
space*, arXiv:2506.13471v2, 20 January 2026, Theorem 1.1,
PDF p.4. For an integral source and a quasi-finite dominant
`Q`-morphism to `P^n` of degree at least five, the number of image
points of standard height at most `B` is `O_f(B^n)`:
[primary PDF](https://arxiv.org/pdf/2506.13471v2),
[retained source](sources/buggenhout-et-al-thin-projective-2506.13471v2.pdf).
The statement was read directly at the applied theorem.

Here the morphism is the component followed by
`U_P -> Pbar`, with `Pbar` the projective plane closure of `P`.
It is quasi-finite, dominant and has degree five or ten. Its
source is integral. The theorem does not require a proper source
or geometric integrality. A fixed rational linear isomorphism
`P^2 -> Pbar subset P^3` compares the two standard heights by
fixed constants. The theorem with `n=2`, summed over the at most
two components, proves (2.1). No deleted-boundary estimate is
needed for Hasse targets, since none lies outside `U_P`.

The bound for `G_P` in (2.3) is PS2, including its boundary
estimate. The exact partition then gives the bound for `L_P`.
\(\square\)

On each fixed plane, the proportion of Hasse targets among all
rational targets is consequently `O_P(1/B)`. On a plane other
than the two sections, the same holds for locally soluble targets,
since the full target count on that plane has order `B^3`.

## 3. Scope and proof dependencies

The result rules out a generic irreducible quadratic-times-cubic
family on a rational target plane, and rules out plane-sized
`B^3` Hasse accumulation. It does not assert that the Hasse set on
another plane is finite or empty, or that it fails to be Zariski
dense there. The bound has no constant uniform over all rational
planes, so summing it over moving planes is not an ambient count.
Nonlinear surfaces and the genus-two resolvent fibration remain
separate sources of possible arithmetic accumulation.

The source-level dependencies are the exact ambient fibre criterion,
PS1--PS3's generic plane and rational-section proofs, the elementary
five-cycle orbit argument above, and the stated projective thin-set
theorem. The affine theorem's integer-weight and absolutely
irreducible leading-form hypotheses are not assumed or used. The
imported quantitative theorem is not reproved here. This is a
written mathematical deduction, so no artificial finite checker is
added; no formal verification, external review or independent
whole-proof replay is claimed.
