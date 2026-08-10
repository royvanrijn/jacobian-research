# Hessian-pencils as an optimized source of `DC_2` symbols

> **Status and scope.** This note deliberately stops extending the normalized
> degree-five symbol and starts a different classical-symbol generator.  It
> gives an exact finite search over 1,540 sparse cubic/quartic Hessian pencils,
> scores 238 noncommuting two-pencil words, and quantizes the top eight rows in
> a declared native-support lattice through `hbar^5`.  Every one-pencil
> survivor is a square-zero Moyal-flat shear.  The hostile two-pencil rows are
> known polynomial automorphism controls, and their factorwise Weyl lifts show
> exactly where the native-support score gives a false positive.  This is a
> search experiment, not a restricted nonquantization theorem and not a result
> about `DC_2`.

The pencil identities are required over each full natural coefficient base
`Q[a]` or `Q[a,b]`.  The experiment does not classify exceptional coefficient
subloci inside a support whose generic two-parameter pencil fails a gate.
The two-pencil ranks and scores are evaluated at the clean point `(a,b)=(1,1)`.

The point is to replace

\[
 \text{continue one exceptional symbol to a larger support}
\]

by

\[
 \text{generate classical families}
 \longrightarrow \text{hard geometric gates}
 \longrightarrow \text{score}
 \longrightarrow \text{quantize only a shortlist}.
\tag{1}
\]

The parent family calculation already says that the old quintic survivor is
globally killed at order seven; see the
[relative family package](QUANTUM_RESIDUE_OBSTRUCTION.md#12-relative-family-package).
The marked-root degree-six, degree-seven, and degree-eight survivors are also
killed at order seven.  Those closures are controls for the selection
procedure, not reasons to increase their correction supports.

## 1. The Hessian-to-symplectic bridge

Use Darboux coordinates

\[
 z=(q_1,q_2,p_1,p_2),\qquad \{p_i,q_j\}=\delta_{ij},
\]

and let `Pi` be the corresponding constant Poisson matrix.  For a polynomial
potential `A`, put

\[
 N_A=\Pi\operatorname{Hess}(A).
\tag{2}
\]

This is the skew-metric analogue of the relative nilpotent endomorphisms in
the [HC4 master reduction](../HC4_RELATIVE_NILPOTENT_MASTER_REDUCTION.md).
Pointwise, `N_A` lies in the symplectic Lie algebra.  If

\[
 N_A^4=0,
\tag{3}
\]

then its Cayley transform is the polynomial matrix

\[
 C_A=(I-N_A/2)^{-1}(I+N_A/2)
    =I+N_A+\frac12N_A^2+\frac14N_A^3.
\tag{4}
\]

It is pointwise symplectic.  Pointwise symplecticity is not enough: a
Jacobian must also be integrable.  The search therefore imposes the exact
row-closure equations

\[
 \partial_k(C_A)_{ij}=\partial_j(C_A)_{ik}
 \quad\text{for all }i,j,k.
\tag{5}
\]

When (5) holds, radial homotopy integrates the rows to a polynomial four-tuple
`F_A`, after which the checker independently verifies

\[
 JF_A\Pi JF_A^{\mathsf T}=\Pi,
 \qquad \det JF_A=1.
\tag{6}
\]

Thus (2)--(6) are a genuine generator of rank-two classical symplectic Keller
symbols.  They do **not** turn an arbitrary HC4 Hessian pencil into a
symplectic map.  The extra Hamiltonian and row-integrability equations are the
bridge, and they are restrictive.

## 2. Exact sparse census

For every one- or two-monomial support in homogeneous degrees three and four,
and every mixed cubic--quartic two-monomial support, use its natural affine
coefficient base `Q[a]` or `Q[a,b]`.  The complete counts are:

| degrees | supports | `N^4=0` | Cayley-integrable | `N^2=0` | Moyal-flat |
|---|---:|---:|---:|---:|---:|
| `3` | 210 | 56 | 40 | 40 | 40 |
| `3+4` | 700 | 124 | 84 | 84 | 84 |
| `4` | 630 | 84 | 60 | 60 | 60 |
| total | 1,540 | 264 | 184 | 184 | 184 |

Here Moyal-flat means that every possible higher odd term vanishes; for these
degrees it is enough to check `Pi^3` and `Pi^5`.  The important negative
result of the experiment is

\[
 \boxed{
 \text{every sparse Cayley-integrable row found here collapses to }N_A^2=0.
 }
\tag{7}
\]

Consequently

\[
 C_A=I+N_A,
 \qquad F_A=z+\Pi\nabla A,
\tag{8}
\]

and all 184 rows are exact one-pencil Weyl controls.  The smallest direct
HC4 cross-link therefore does not yet produce a hostile symbol.

## 3. Optimization data

The search keeps hard admission gates separate from the heuristic score.
Every scored row must first be a polynomial four-tuple satisfying (6).  A
future `DC_2` candidate must additionally have a proved nonautomorphic or
noninjective classical status; the two-pencil rows below deliberately fail
that final gate and are calibration controls.

For a common correction lattice, take the componentwise monomial support at
`(a,b)=(1,1)`.  The exact finite complex records:

- localization pole-width;
- dimension and relative size of the restricted obstruction cokernel;
- rank jump after adjoining the Moyal section;
- Hamiltonian gauge-image rank for Hamiltonians of degree at most four; and
- parameter-base complexity, counted as parameters plus relations plus
  localization divisors.

The scalar score rewards relative cokernel size, raw odd-Moyal support, and a
section rank jump, and penalizes gauge rank, poles, and base complexity.  It
is only a queueing device inside this one native lattice.  The component
vector, not the scalar, must be used when comparing a root-boundary lattice,
a reciprocal localization, or another ordering.

## 4. Noncommuting two-pencil calibration

There are 28 one-monomial primitive pencils.  Up to reversal/inversion and
parameter sign, 238 pairs do not Poisson commute.  Twenty-four resulting
two-step words have a nonzero raw odd Moyal defect.  All maps have pole-width
zero and base `A^2`.

The top eight split into two symmetry classes:

| representative `(A,B)` | score | `C^2` | image | cokernel | gauge | native result |
|---|---:|---:|---:|---:|---:|---|
| `(p1*q2^3, p2^2*q1^2)` | 133 | 39 | 14 | 25 | 3 | rank `14 -> 15` at `hbar^3` |
| `(p1*q2^3, p2^3*q1)` | 89 | 30 | 12 | 18 | 3 | lifts through `hbar^5` |

Each row represents four polarization symmetries in the checked shortlist.
The first class is exactly the sort of row a naive hostility score would
select: zero poles, a two-parameter base, a large restricted cokernel, small
gauge image, and a nonzero obstruction section.

It is nevertheless an automorphism false positive.  For the displayed top
row the factorwise Moyal exponential gives

\[
 U_2=(0,0,-3q_2,0),\qquad U_4=0.
\tag{9}
\]

The single monomial `q2` is outside the native support of the third component.
After adjoining it, the exact Moyal canonical relations hold through
`hbar^5`; factorwise composition gives an exact Weyl automorphism in every
order.  The other three top rows have the same one-monomial escape after the
corresponding polarization symmetry.

Thus the optimization needs two safeguards:

1. known automorphism controls cannot enter the `DC_2` candidate queue, no
   matter how hostile their restricted coordinates look; and
2. a high score should trigger one deliberate support-saturation check before
   an expensive higher-order calculation.

This is useful calibration rather than a failed search: it distinguishes a
structurally narrow lattice from a genuinely hostile classical symbol.

## 5. Next candidate queue

The experiment leaves two genuinely different next gates.

1. **Reciprocal `R21` admission.**  The
   [classical-symbol census](DC2_CLASSICAL_SYMBOL_CENSUS.md) already supplies
   a simple quadratic coefficient base and controlled localization exponent
   one, but its polynomial rank-two completion is open.  Proving or refuting
   that admission should precede PBW work.  If admitted, it is the smallest
   currently recorded non-marked-root family on which to run the score.
2. **Higher-nilpotence Cayley integrability.**  Search beyond two-monomial
   supports, or solve the coefficient scheme defined by (3) and (5), with the
   explicit objective `N^2 != 0`.  A survivor would be a genuinely new
   Hessian-derived symplectic family.  Only after checking its classical orbit
   status should its boundary lattice and restricted PBW complex be built.

The first part of this gate has now been reached by the
[regular-index-four frontier](DC2_HIGHER_NILPOTENCE_R21_FRONTIER.md): an
all-degree triangular family has `N^4=0` and `N^2!=0` on `c!=0`, and its
Cayley matrix is integrable. Exact inversion and the Moyal calculation show
that it is a polynomial automorphism control. Moreover, row closure in the
whole surrounding triangular ansatz forces `N^2` to be constant. The next
search must therefore move the nilpotent flag rather than merely add support
inside that chain. The same frontier verifies stable unimodular charts for
the reciprocal powers in `R21`, excludes its affine-contact shortcut, and
uses eleven exact tame shears to remove the graph defect through degree four.
Its explicit degree-five remainder cannot be finitely resummed inside the
same fiber-preserving subgroup, because `R=x*S` is reducible. The remaining
graph--Darboux search must genuinely mix the stable variable into the base.
The exponent-two chart now does so: twenty-six exact shears kill its defect
through degree six, leaving a nonzero degree-seven row. Euler homotopy gives
the all-order recurrence and formally kills every later homogeneous row.
A constant-Pfaffian dilation has a polynomial vector field but a
nonpolynomial time-one inverse, while the exact `b=0` kernel excludes every
normalizer whose target-`b` component has a polynomial Poisson mate on its
zero fiber. In particular, no elementary target-`b` coordinate can work.
The tame family `F_k=a+c^k*(b+a*d)` shows that no-slice coordinates do exist,
starting in degree three. Its curve pole order and transverse weight are tied
as `(k,k)`, whereas the R21 kernel has split signature `(2,3)`. The complete
constant rings are Danielewski surfaces with intrinsic exponents `k` and
three, so the constant ring forces `k=3` while the generic time divisor
forces `k=2`; every row in this family is excluded. For `k=2` an exact
conjugacy exists on `I!=0`, with Jacobian `-I/5184`. Any surviving resummation
must correct that affine modification over `I=0` or change the target
polarization.

The two-pencil words remain permanent zero-pole calibration rows.  Continuing
the old quintic support, or automatically extending the marked-root degree
ladder, is not selected by this experiment.

## Reproduction

Run:

```bash
.venv/bin/python scripts/search_dc2_hessian_symbol_candidates.py \
  --shortlist 8 \
  --output \
  artifacts/generated-results/dc2_hessian_symbol_optimization.json
```

The command verifies the nilpotence, Cayley, integrability, symplectic,
Jacobian, Moyal, correction-rank, gauge-rank, and factorwise quantization
claims exactly over `Q`.
