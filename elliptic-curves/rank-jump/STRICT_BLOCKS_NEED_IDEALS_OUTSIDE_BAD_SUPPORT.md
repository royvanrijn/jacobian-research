# Strict blocks need ideal classes surviving outside the bad support

The latest constructor test produced no new class. On the fixed 103b2
high/low pair, none of the 60 squared bad-prime ideals yielded a generator
in the frozen reduced-lattice test. All six generic positive controls
were recovered. These misses do **not** prove nonprincipality or small
class rank.

A separate retrospective deduction explains a limitation of a broader
bad-support strategy. On the earlier completely measured R17 high and
low controls, **no nonzero class in the known strict block has an S-unit
representative**. On the earlier MW16 high, at most one dimension does.
Thus bad-support units alone cannot carry those observed blocks.

This does not identify the mechanism creating the fresh +10/+11 or
historic +12/+14 classes. It narrows a possible constructor and keeps
the direct arithmetic attempt distinct from a diagnostic using known
points.

## Equation-only experiment on the fresh matched pair

The [protocol](BAD_PRIME_PRINCIPALIZATION_PROTOCOL.json) fixes the already
retained 103b2 fibres at 3726/881 and -1049/2296. The arithmetic consumes
their equations, certified maximal orders, complete bad-prime sets, and
three generic classes as positive controls. Exceptional points and their
derived classes are excluded from these inputs.

For every prime ideal P above S={2,bad primes}, form I=P². PARI ideal
reduction supplies J and an exact multiplier a with I=(a)J. Trace-LLL
reduction of J then supplies a basis in which all 49 primitive vectors
in [-2,2]³ up to sign are tested. A vector z with |N(z)|=N(J) proves
(z)=J, hence (az)=P². The norm-square generator is still subject to
strict local squareness and independence; finding one would not alone
prove class creation or rational solubility.

This differs from the preceding good-prime test by using the bad-prime
population and reducing the ideal within its class before looking for a
generator. It is a finite principal-generator detector, not a class-group
algorithm.

| Frozen fibre | Rational primes in S | Squared prime ideals tested | Candidate generators | Generic controls recovered |
|---|---:|---:|---:|---:|
| 103b2, 3726/881 | 12 | 29 | 0 | 3/3 |
| 103b2, -1049/2296 | 13 | 31 | 0 | 3/3 |

Both workers completed within the unchanged sixty-second, 256 MiB PARI
stack bounds. No pair of candidate ideals had identical retained reduced
HNF either. Such an equality would have given an exact principalization
of their quotient squares; different reduced HNFs do not prove distinct
ideal classes, since this reduction is not a canonical class invariant.

The [certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_bad_prime_principalization_v1.json)
retains every target, reduced ideal, multiplier, basis transformation,
norm and generator. The
[independent replay](../../artifacts/generated-results/elliptic-curves/rank_jump_bad_prime_principalization_verification_v1.json)
checks all **3234 norms** with rational multiplication matrices, and
verifies each ideal transport by an integral unimodular change of lattice.
It also independently reconstructs the complete prime-ideal population.
The verifier's first attempt encountered a Sage/Fraction conversion
error; it was corrected before any verifier artifact was written.

## A bound on the entire strict S-unit carrier

This section uses only the **earlier oracle-derived Artin matrices** as
a retrospective diagnostic. They do not enter the new experiment above.
They concern different, earlier high/low controls; no result here is
silently transferred to the fresh pair.

Let C=Cl(O_K,S) and U be the strict Kummer space. For beta in U, define

    Phi_S(beta) = [J_beta] in C[2],       (beta)=J_beta².

Strictness makes every valuation even. The fundamental exact criterion is

\[
 \boxed{\Phi_S(\beta)=0
 \iff [\beta]\text{ has a representative in }\mathcal O_{K,S}^{\times}.}
\tag{1}
\]

Indeed, if the half ideal is principal after localization, write it as
(a) there; then beta/a² is an S-unit. Conversely such a representative
has trivial half ideal in the localized class group. This argument allows
all fractional ideals and every square multiple; it is not limited to a
finite list of principal generators or to small units.

For the known strict subspace V with basis beta_i, the existing matrix is

    M_ij = chi_beta_i(Phi_S(beta_j)).

By (1), any class of V with an S-unit representative lies in the **right**
kernel of M. Therefore

\[
 \dim(V\cap\mathrm{image}(\mathcal O_{K,S}^{\times}))
 \le \dim V-\operatorname{rank}M.
\tag{2}
\]

The existing ideal identities and Artin evaluations are certified in
[the elementary S-class factor theorem](SOLUBLE_ELEMENTARY_S_CLASS_BLOCKS.md).
The new
[bound artifact](../../artifacts/generated-results/elliptic-curves/rank_jump_strict_sunit_carrier_bound_v1.json)
enumerates the complete binary right kernels independently of the old
rank routine:

| Earlier control | Known strict dimension | Artin rank | Strict S-unit dimensions in that space, at most |
|---|---:|---:|---:|
| MW16-05 high, 307/206 | 10 | 9 | 1 |
| Published R17 high, -2300/843 | 8 | 8 | 0 |
| Published R17 low, -1561/3133 | 6 | 6 | 0 |

The one possible MW16 direction is the previously retained **right-kernel
mask 637**. It is not proved to be an S-unit: other strict characters
could detect its half ideal, or that ideal could be a nonzero element of
2C. A zero column pairing against the known characters is only a bound.

The low control also has no strict S-unit directions in its known space.
Thus this property is not a high-jump discriminator. Its use is to reject
S-unit-only constructions of the already measured blocks. Their half
ideals must survive localization, and most even survive in C/2 as an
elementary direct factor.

There is a crucial limitation on generic corrections. Equation (2) does
not exclude arbitrary products of non-strict S-units with non-strict
generic classes whose local obstructions cancel. To bound that larger
constructor one must additionally measure the full generic half-ideal
image in C[2], not just its strict subspace. No such enlargement or
unproved inference is made here.

## What this changes

* **Incidence:** a mechanism based solely on S-unit squareclasses is too
  small for the earlier measured soluble blocks. The missing class
  supply must involve half ideals that remain nonprincipal after the
  bad primes are inverted, or a separately certified generic correction
  construction beyond the scope of (2).
* **Weak constructor:** isolated bad-prime-square principalizations,
  even after ideal reduction, supplied nothing on the fresh matched
  pair. The successful positive controls certify only generator
  recovery for those controls. This is not evidence of absent classes
  on either fibre.
* **Missing computation:** obtain independent mixed ideal relations and
  their square principalizations outside the inherited class span, with
  complete strict local tests. A basis or sufficient independent
  characters of the cubic S-class group remains unknown on the fresh
  pair independently of exceptional points.
* **Solubility:** no new CT entry or rational elliptic cover was computed.
  The fixed-incidence six-direction switch remains a separate control.
  None of these outputs is a point-visibility feature or a candidate
  score for Agent 1.

The new files preserve the null experiment and the exact diagnostic
without changing active searches, their parameters, or MATH_STATUS.

```sh
timeout 60 sage -python elliptic-curves/rank-jump/verify_bad_prime_principalization.py check
python3 elliptic-curves/rank-jump/strict_sunit_carrier_bound.py check
```
