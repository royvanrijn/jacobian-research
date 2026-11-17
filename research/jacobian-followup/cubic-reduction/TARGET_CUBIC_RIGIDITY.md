# Cubic rigidity of target changes for the saved twelve-variable map

KCUB26, 15 November 2026. Fix the integral twelve-variable cubic K of
[KCUB3](PROOF12.md). Two complete target-change obstructions hold over every
characteristic-zero coefficient field k.

1. If R in k[Y1,...,Y12] has degree<=5 and R(K) has degree<=3, then R is
   affine. Consequently a nonlinear polynomial target map of degree<=5
   cannot turn K into another cubic map. This statement also holds after
   arbitrary invertible affine source and target frame changes.
2. For each i, the intersection of k(K_j:j!=i) with the degree<=3 source
   polynomials is exactly span_k(1,K_j:j!=i). Thus if a polynomial or
   rational function q of the other eleven displayed outputs has polynomial
   cubic pullback, q is affine in those outputs. There is no degree bound
   on q. In particular every single-coordinate polynomial target
   automorphism preserving cubic degree in the displayed output frame
   is affine.

No smaller positive map is constructed. These conclusions do not prove a
full target-orbit, arbitrary affine sparsity or global dimension/term
optimum. They do not exclude coordinated nonlinear source/target changes,
temporary higher degrees, or general target maps of degree>=6 involving
all outputs. Arbitrary affine mixing of target coordinates is covered by
statement1, not by the unbounded rational-function assertion in statement2.
All current positive bounds and the strict62/12 profile partition remain
unchanged.

## Nonlinear targets of degree at most five

There are6175 nonlinear target monomials of degrees2 through5. Their
higher-than-cubic source coefficient columns have a6175-by6175 minor
nonzero modulo p=1000000007, recorded in the
[higher target certificate](target-pullback12-higher/checkpoint.json).
The earlier degree2/3 control has442 columns and7276 source coefficient
rows, with its own nonzero minor retained in the
[target pullback certificate](target-pullback12/checkpoint.json).

The independent checker reconstructs every selected pullback coefficient
modulo p and the entire selected minor, and performs row-based modular elimination. The primary producer
used column elimination; no primary rank assertion is trusted. A nonzero
integer minor modulo p is a nonzero rational minor, hence remains nonzero
in every characteristic-zero extension. Thus the6175 columns are linearly
independent. The affine part of R contributes nothing in source degree>3,
so all nonlinear coefficients of R must vanish.

For affine source A and affine target B, degree is preserved by A and by
composition of a scalar target polynomial with B. If R(B(K(A(X)))) is
cubic and degR<=5, then (R composed B)(K) is cubic and has target degree<=5.
The fixed-frame result makes R composed B, hence R, affine. This proves
the stated affine-frame extension without sampling frames.

## Polynomial derivations from the exact inverse Jacobian

Let J=JK and let V=J^-1. The explicit integer-polynomial
[inverse witness](target-derivation12/inverse.json) has4620 monomial
occurrences and degree15. The independent checker verifies both JV=I
and VJ=I as full polynomial matrix identities. Define

```
Delta_i = sum_a V_(a,i)(X) * d/dX_a.
```

The chain rule gives Delta_i(K_j)=delta_ij. Therefore every rational
function of the other eleven outputs is killed by Delta_i. No denominator
bound, polynomial target-degree truncation or generic specialization is
used in this necessary condition.

Apply Delta_i to all455 source monomials of degree<=3, including1.
For each of the twelve i, a443-by443 minor of this coefficient matrix has
nonzero determinant modulo p. The twelve selected row/column sets are
recorded in the [derivation certificate](target-derivation12/checkpoint.json).
Thus each operator on the455-dimensional cubic source space has kernel
dimension at most12. Its kernel already contains the12 independent
polynomials1,K_j for j!=i. Their independence follows from the identity
linear part of K. Hence, for every i,

```
ker(Delta_i on k[X]_(degree<=3)) = span_k(1,K_j:j!=i).
```

The complete dimension check is

| Outputs omitted one at a time | Source columns | Certified rank | Kernel dimension |
| --- | ---: | ---: | ---: |
|i=1,...,12, separately|455|443|12|

The cubic bound here is a bound on the desired resulting source polynomial,
not on the unknown target q. It is what makes the unbounded-degree target
claim finite to verify.

## From the kernel equality to arbitrary rational target functions

Suppose q in k(Y_j:j!=i) and f=q(K_j:j!=i) is a polynomial of degree<=3.
Then Delta_i f=0, so f=a0+sum_(j!=i) a_j*K_j. Pullback by K is injective
on polynomials: because K=X+terms of order>=2, the lowest nonzero homogeneous
part of any polynomial P is unchanged in P(K). It therefore also embeds
the corresponding rational function fields. The equality for f consequently
forces q=a0+sum a_j*Y_j. This proves statement2, including rational q with
arbitrary degree and denominator.

For a single-coordinate polynomial target automorphism, the changed row
has form a*Y_i+q(Y_j:j!=i), a a nonzero scalar. If its pullback remains
cubic, q(K_j:j!=i) is cubic and q must be affine. The special shear case
Y_i+q is included. This does not classify a sequence whose intermediate
maps exceed cubic degree or a map changing several rows simultaneously.

## Exact reproduction and trust boundary

Run

```
.venv/bin/python research/jacobian-followup/cubic-reduction/verify_target_rigidity12.py
```

The checker uses standard-library integer/Fraction arithmetic only. It:

* reconstructs the6175 target coefficient columns through degree5 modulo p,
  verifies the selected integer minor and the exact442-column lower-degree control;
* differentiates K and verifies both polynomial inverse identities;
* reconstructs all455 source monomials and all twelve derivation operators;
* verifies a nonzero443 minor for each operator by independent row elimination;
* checks the exact kernel generators via Delta_i(K_j)=delta_ij and their
  independent constant/linear jets.

The modulus is checked prime by trial division. All certificate matrix
entries are integers, so the rank conclusions are characteristic-zero
statements, not heuristic finite-field evidence. The chain-rule,
injectivity and kernel-dimension arguments above are written mathematical
proofs; no proof-assistant formalization or external referee is asserted.
Initial primary screen receipts and producer snapshots are preserved.

- [Independent validation](target-derivation12/validation.json)
- [Degree2/3 control producer](target_pullback12.py)
- [Degree4/5 primary extension](target_pullback12_higher.py)
- [Inverse/derivation producer](target_derivation12.py)
- [Trusted cubic input and full determinant/collision replay](integer-validation12.json)

## Consequences for construction work

A further cubic simplification cannot start with a nonlinear target map of degree<=5, in any affine frame of K. It also cannot start with a
nonlinear single-coordinate target change of any degree in the displayed
output axes. KCUB19 separately closes same-axis nonlinear source shears;
KCUB25 closes identity-factor reduction by stable conjugacy.

Useful remaining routes include coordinated unequal source/target changes,
temporary degree growth, higher-degree target changes after a genuinely
changed output frame, different graph directions with joint repair, or
another construction. The ordinary cubic upper bound remains12 variables
with36 terms; its global dimension/term minima are open. The other bounds
remain cubic-linear89, raw3D13, normalized[10,48], global3D degree[5,7],
and exact cubic-component degree9. No full simplicity completion follows.
