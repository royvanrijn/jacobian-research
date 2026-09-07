# A full governing block at fixed Selmer incidence

The two cross-cochains complete the governing data for the retained small
curve's **entire three-dimensional 2-Selmer space**. Among the predeclared
odd primes at most 199, the twists by 41 and 113 preserve this full space
and have zero full Cassels--Tate pairing. The twist by 97 preserves the
same space but has CT rank two, so its Mordell--Weil rank is at most one.

The ranks at 41 and 113 remain UNKNOWN. No point search or elliptic descent
was performed. This experiment closes the restricted-versus-full-pairing
gap in this control; it does not close the rational-versus-higher-Sha gap.

## Fixed classes and the two new norm witnesses

Use the already certified control

\[
E_0:y^2=x^3-11x^2-14x-1,
\quad K=\mathbf Q(\theta),\quad \theta^3-11\theta^2-14\theta-1=0.
\]

The [exact quotient comparison](NORM_LIFTS_CAN_BE_ENTIRELY_SHA.md) proves

\[
S=\operatorname{Sel}_2(E_0)=\langle\alpha,\beta,\eta\rangle,
\quad \dim S=3,
\]

where

\[
\alpha=\theta^2-10\theta+1,\quad
\beta=\theta^2-13\theta+12,\quad \eta=-1-\theta.
\]

The first two classes were selected from the cubic class group before
their elliptic solubility was tested. They are strict at {2,163,infinity}.
The third is the retained non-strict class of the fixed point (-1,1).
This is a retrospective small control, not a masked high-rank selector.

In this ordered basis the original full pairing is

\[
B_0=\begin{pmatrix}0&1&0\\1&0&0\\0&0&0\end{pmatrix}.
\]

The upper-left block comes from the independently verified
[dyadic switch](UNPOINTED_NORM_COCHAINS_AND_THE_DYADIC_BLOCK_SWITCH.md).
The last row vanishes because eta is rational on E0. The full dimension
three comes from the earlier exact Selmer argument, not from the sizes
of the matrices computed here.

The [protocol](FULL_SMALL_GOVERNING_BLOCK_PROTOCOL.json) fixes two norm
equations, each with a 60-second worker cap. Their exact witnesses are

\[
\begin{aligned}
X_{02}&=-119/25-(163/50)\theta+(13/25)\theta^2,\\
Y_{02}&=-929/250-(402/125)\theta+(83/250)\theta^2,\\
X_{12}&=1,\qquad Y_{12}=(2+\theta)/5,
\end{aligned}
\]

with X_02^2-alpha Y_02^2=eta and X_12^2-beta Y_12^2=eta.
No elliptic point is used to construct these witnesses. The norm
solubility gate opens because one argument is strict and the other is
S-unramified, as proved in the preceding note.

Apply its norm-cochain formula to these two witnesses and reuse the
existing alpha,beta witness. Exact changes of generators give three
small octics, ordered by the slots (01,02,12):

\[
\begin{aligned}
g_{01}(T)&=T^8-18T^6+92T^4-112T^2+16,\\
g_{02}(T)&=T^8-16T^6+90T^4-216T^2+841,\\
g_{12}(T)&=T^8+18T^4-200T^2+625.
\end{aligned}
\tag{1}
\]

Their splitting groups have orders 24,96,96. Every norm identity and
degree-eight change of generator is independently checked using rational
polynomial arithmetic. The 24-group's special cyclic relation was proved
in the preceding note; no S3 independence theorem is invoked for this
cyclic cubic.

## Ramification: distinguish an order index from a field obstruction

The polynomial discriminants of g02 and g12 contain primes 5 and 29.
Those factors must not be treated automatically as governing ramification.
PARI supplied candidate integral bases. The independent verifier checks
closure under multiplication for all 64 basis products in each algebra,
checks that the bases contain one, and obtains the order discriminant

\[
2^{18}163^4
\]

in both cases. A finite free unital ring proves its elements integral;
the determinant formula computes its discriminant. Thus the field
discriminant divides this integer, with no claim that the candidate order
is globally maximal. Together with disc(g01)=2^36*163^4, this proves that
the full governing compositum is unramified outside {2,163}.

The table still excludes non-squarefree reductions of the displayed
polynomials at 5 and 29. This is a limitation of those factorization
presentations, not a claim of field ramification there.

## Three independent central obstruction bits

The automorphism tau(theta)=theta^2-12theta-2 cycles the alpha conjugates,
with tau(alpha)=beta and the third conjugate equal to alpha beta times
a square. Their strict squareclasses span a two-dimensional space over K.
The conjugates of eta have product one; eta and tau(eta) have independent
real sign vectors and hence add two independent squareclasses. The joint
class field L therefore has

\[
\operatorname{Gal}(L/\mathbf Q)=(V\times V)\rtimes C_3,
\quad [L:\mathbf Q]=48.
\]

Here V is the even-mask subspace of F2^3, and its dot product is the
additive Weil pairing. If e is the alpha sign mask and c the eta sign
mask, the beta mask is lambda e, where lambda cyclically shifts coordinates.
The three cochains have coboundaries

\[
\omega(\sigma,\tau)=
\bigl(e\cdot g\lambda e',\ e\cdot gc',\ \lambda e\cdot gc'\bigr).
\tag{2}
\]

They yield an extension with central subgroup F2^3. This subgroup is
actually full, not merely an upper bound for a compositum:

* Commuting two suitable strict translations produces the central vector
  (1,0,0).
* Commuting a strict translation with two appropriate eta translations
  produces (0,1,0) and (0,0,1), since e and lambda e are independent.

Surjectivity onto the order-48 class group supplies lifts of these
translations in the actual Galois image. Their commutators force all
eight central elements into that image. Consequently

\[
\boxed{[N:\mathbf Q]=48\cdot8=384.}
\]

The verifier realizes this group faithfully on the three octics' 24
labelled roots and checks all 384^2 composition identities. Its derived
subgroup has order 128, so its abelianization is C3, with maximal abelian
subfield exactly K.

For a fixed-point-free Frobenius g, define the vector

\[
\psi_{ij}=e_2((g-1)^{-1}a_i,a_j)+\gamma_{ij}.
\]

All eight vectors occur, each on exactly 32 of the 256 elements lying
above the two three-cycles. Although the classes share Galois structure,
their **three obstruction bits are independent** in this governing group.
This control supplies no reduction of the number of full-matrix conditions.

## Exact prime conditions and full CT matrices

Let ell be a positive prime outside {2,163} satisfying

\[
\ell\equiv1\pmod8,\qquad
\left(\frac{\ell}{163}\right)=1,\qquad
f\bmod\ell\text{ irreducible}.
\tag{3}
\]

These conditions make the twist character locally trivial at 2,163 and
infinity, with its only new ramification at a fixed-point-free prime.
[Morgan, Lemma 2.8 and Proposition 3.3](https://arxiv.org/pdf/2309.02374v2)
therefore preserve every local Kummer image and the **entire** Selmer group,
and compute the full CT difference from the three cochains.

No End(V)=F2 hypothesis is needed for these comparison results. We computed
the governing group directly in the cyclic case rather than applying the
paper's later maximal-image theorem.

At a squarefree reduction of an octic in (1), factor degrees (1,1,3,3)
give its bit zero, and (2,6) give its bit one. Write (u,v,w) for the three
bits. The full pairing on S is

\[
B_\ell=\begin{pmatrix}
0&1+u&v\\1+u&0&w\\v&w&0
\end{pmatrix}.
\tag{4}
\]

In particular, it is zero precisely when (u,v,w)=(1,0,0). That condition
concerns all three entries; it is stronger than merely making the strict
two-class block isotropic.

The fixed prime range gives:

| ell | Governing bits (01,02,12) | dim Sel2 | Full CT rank | Full radical dimension | MW rank bound |
|---:|---|---:|---:|---:|---:|
| 41 | (1,0,0) | 3 | 0 | 3 | <=3 |
| 97 | (0,0,1) | 3 | 2 | 1 | <=1 |
| 113 | (1,0,0) | 3 | 0 | 3 | <=3 |

These are exactly the eligible primes in the declared table; the range
was not enlarged to obtain a zero matrix. The models, if needed, are
y^2=x^3-11 ell x^2-14 ell^2 x-ell^3. None was sent to an elliptic
point-search or rank-descent routine.

For an independent arithmetic check, at each eligible prime we evaluate
N(X+sqrt(beta)) in the cubic finite field, selecting the square root with
the prescribed norm. Its Legendre symbol reproduces the octic bit. Across
all 29 joint inert-prime rows, 86 of 87 such norm values are nonzero and
agree. The remaining value is degenerate modulo 3 for the old alpha,beta
norm witness; its bit is supplied by the reduced octic factorization.
All nine norm evaluations for the three eligible twists are nondegenerate.

The rank bound at 97 is unconditional because rational Kummer classes lie
in the full CT radical. It is not an exact-rank statement. At 41 and 113,
zero full CT gives second-descent lifting and does **not** prove rank three.
Higher-divisible Sha remains possible even if Sha is finite.

## An arithmetic baseline, not a rank probability

Condition (3)'s local-square restrictions mean splitting completely in
D=Q(zeta8,sqrt(-163)), a degree-eight extension. The governing group has
abelianization C3, so N intersects D trivially. Chebotarev in ND therefore
preserves the equal counts of all eight governing vectors among primes
satisfying (3). Each full CT matrix has conditional density 1/8.

This is a theorem-derived obstruction baseline, not a frequency inferred
from three twists and not a rank distribution. In particular, the two
zero matrices in the small table do not suggest an elevated rank-three
rate. Also, requiring ell=1 modulo 163 would force splitting in the
cubic subfield and destroy the desired inertness; the quadratic-residue
condition in (3) is the relevant local condition.

## What this changes, and what remains

1. **Incidence:** the three tested twists have exactly the same full Selmer
   space and local Kummer images. This is stronger than equal dimensions
   or a shared strict subspace.
2. **Solubility obstruction:** explicit field factorizations distinguish
   rank-two CT from zero full CT before supplying any points on the twists.
   The exact full-radical computation replaces an earlier restricted-block
   diagnostic.
3. **Weak mechanism:** sharing a cyclic class-group block reduces the class
   field, but does not force several obstruction entries to vanish together.
   The three central bits remain independent. Neither governing-field
   existence nor a small class field is by itself a positive rank feature.
4. **Next decisive question:** determine whether the full-radical classes
   at 41 and 113 are rational or survive into higher Sha. These two fixed
   controls now give precise targets for a further descent or a global
   rational construction, without another candidate sweep.
5. **Production gap:** this experiment concerns a small cyclic-cubic twist
   family, not t in R17/A1. Production cubic classes selected independently
   of exceptional points and an applicable family-comparison theorem are
   still needed. No active search scoring or protocol change follows.

There is no visibility measurement in this experiment.

## Reproducibility

The [input artifact](../../artifacts/generated-results/elliptic-curves/rank_jump_full_small_governing_block_inputs_v1.json)
contains the two deterministic norm outputs and exact reduction maps.
The [verification artifact](../../artifacts/generated-results/elliptic-curves/rank_jump_full_small_governing_block_v1.json)
retains the integral orders, group counts, every prime row and full matrices.
Candidate order bases came from PARI nfbasis; the stored multiplication and
discriminant checks establish the needed ramification bound independently.

```bash
python3 elliptic-curves/rank-jump/verify_full_small_governing_block.py check
```

Both norm workers completed in seconds, and the independent standard-library
replay takes under a second. Inputs and sources are hash-bound; outputs are
immutable. The active high-rank search, shared navigation, MATH_STATUS and
STATUS remain untouched.
