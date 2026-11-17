# The ten associated primes of the full three-factor moment ideal

Claim: `G3PRODUCTASS20261114`. Characteristic zero. This is an internal
written proof supported by exact SymPy and Singular calculations, without
independent replay, external review or formal verification.

Let R=Q[a,b,c,d,e,f] be the universal Gram ring, with covariance matrix
[[a,d,e],[d,b,f],[e,f,c]]. Write I=I_infinity=I10. The
[global ideal theorem](THREE_FACTOR_GLOBAL_MOMENT_IDEAL.md) supplies its five
moment generators. The [colon theorem](THREE_FACTOR_COLON_GEOMETRY.md)
supplies J^13 contained in I, where J=(af,be,cd,abc,def).

The associated primes of R/I are exactly:

| Type | Primes |
|---|---|
| Minimal rows | (a,d,e), (b,d,f), (c,e,f) |
| Minimal stars | (b,c,f), (a,c,e), (a,b,d) |
| Embedded diagonal axes | (b,c,d,e,f), (a,c,d,e,f), (a,b,d,e,f) |
| Embedded origin | (a,b,c,d,e,f) |

There are no embedded planes, off-diagonal axes, or noncoordinate primes.
The six isolated primary components are also determined by saturation.
The [subsequent constructive theorem](THREE_FACTOR_PRIMARY_DECOMPOSITION.md)
now gives the four embedded primary ideals and intrinsic origin torsion77.
The earlier four-component decomposition of I8:M10 is a different result.

## 1. A full row family, including the elimination pole

On b*c nonzero, adjoining square roots of b and c is faithfully flat and
allows form rescalings to b=c=1, f=z. Consider the resulting ideal in
Q[z][a,d,e], with weights wt(a)=2, wt(d)=wt(e)=1 and wt(z)=0.
Each M_(2m) has weighted degree 2m.

For degree k, form every monomial cofactor of each of the five generators
whose weighted degree is k. Their coefficient columns give the complete
relation matrix on the weighted monomials of degree k. Exact Smith normal
forms over the PID Q[z] have only nonzero constant invariant factors for
k=0,...,9. The free ranks are:

| k | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
|---|---|---|---|---|---|---|---|---|---|---|
| rank | 1 | 2 | 3 | 4 | 4 | 4 | 3 | 2 | 1 | 0 |

Degree 10 vanishes by a short identity, avoiding the expensive Smith step.
Every degree-10 monomial containing d or e is a multiple of degree 9.
The only other monomial is a^5. The coefficients of a in M2 and a^2 in M4
are t=1+2z^2 and h=27+216z^2+72z^4. They are coprime. Multiplying these
moments by a^4 and a^3, then using the recorded polynomial Bezout identity,
kills a^5. Degrees 9 and 10 zero imply every later degree zero by induction
on the generators a,d,e.

Thus the entire normalized quotient is Q[z]-free of rank24, with no
excluded parameter value. In particular, z^2=-1/2 is included: the earlier
formula solving M2 for a had a pole there, but the family does not.

The positive-weight ideal (a,d,e) is nilpotent. Over the fraction field the
quotient is an Artinian local algebra. If g has nonzero constant term in
Q[z], its multiplication determinant on the free module is that constant
term to the 24th power, hence is nonzero. Therefore g is a nonzerodivisor.
The normalized ideal is (a,d,e)-primary. The same argument applies after
the faithfully flat normalization above and descends. Consequently

    Q_row,A = I:(bc)^infinity

is the isolated (a,d,e)-primary component. Its two label permutations give
the other rows. The family on b*c nonzero is finite locally free of rank24
over Q[b,b^-1,c,c^-1,f].

## 2. A full star family without excluding a=0

On d*e nonzero, form rescalings normalize d=e=1 and leave a free. The
radical is (b,c,f). Since J becomes (b,c,f), the known J^13 containment
makes this quotient a finite Q[a]-module.

On a nonzero it is the normalized a=d=e=1 star algebra after an invertible
rescaling; the retained exact star calculation gives a free rank6 family
there. At a=0, M2 eliminates f=-(b+c)/4. The next two moments become a
homogeneous binary quadratic and cubic. Their exact nonzero resultant
and the nonzero leading coefficient exclude a common projective zero.
They form a regular sequence of degrees2 and3, giving length6 and no
degrees at least4. Every later moment lies in those degrees.

A finite module over Q[a] is a free module plus torsion. Any torsion here
is supported at a=0, while equality of the generic and a=0 fiber lengths
forces that torsion to vanish. Hence the full star family is Q[a]-free
of rank6. The same determinant argument proves primaryness. Therefore

    Q_star,A = I:(de)^infinity = I:(ade)^infinity

is the isolated (b,c,f)-primary component. Its permutations give the other
stars. No central-variance exception remains on this open set.

## 3. Why only finitely many coordinate candidates remain

The six open sets bc,ac,ab,de,df,ef nonzero now carry primary quotients,
so no embedded associated prime meets them. Every prime contains a
minimal prime from the known six-plane radical. Inside those coordinate
planes, removing these open sets leaves only the six coordinate planes
with free coordinates ad,ae,bd,bf,ce,cf, their coordinate axes and origin.

I is homogeneous for the three-dimensional torus rescaling the three
linear factors. Its finite set of associated primes is torus invariant;
connectedness makes each individual prime invariant. On each remaining
two-coordinate face the two covariance weights are linearly independent.
The torus acts transitively on its open coordinate torus after algebraic
closure. Thus no proper noncoordinate invariant closed subset remains.
Every possible embedded prime is coordinate.

Label symmetry reduces the candidates to four tests: free coordinates
ad, free coordinate a, free coordinate d, and no free coordinate. For each,
pass the free variables to the coefficient field and compute I:m, where
m is generated by the remaining variables. Its quotient by I is the
socle at the candidate coordinate origin. Flat localization and the
coordinate-candidate argument justify recovering the original associated
prime from this test; this is not a numerical specialization.

| Free variables | Socle dimension | Conclusion |
|---|---:|---|
| a,d | 0 | All six candidate embedded planes excluded |
| a | 7 | All three diagonal axes occur |
| d | 0 | All three off-diagonal axes excluded |
| none | 8 | Origin occurs |

The axis-a socle has filtered Hilbert polynomial 2t^13+t^14+4t^17.
Its displayed 13 nonzero reduced quotient generators are redundant; the
socle dimension is7. The origin socle is graded, with Hilbert polynomial
5t^21+3t^30. One origin witness is a^5*b*e^2*f^13: it survives modulo I,
and every variable times it lies in I. Exact quotient calculations, both
forward membership tests and complete Hilbert numerators are retained.
The completeness assertion uses Singular's exact ideal-quotient algorithm;
the membership tests alone would not prove it.

## 4. The intrinsic embedded defect

Let U be the intersection of the six isolated primary components.
The exact saturations computed for row bc and star ade stabilize at
powers13 and11 respectively. Their symmetric intersection contains I
strictly: 103 reduced standard-basis generators of U survive modulo I.
These counts depend on the chosen standard bases.

D=U/I is the canonical embedded defect. It is supported on the three
diagonal axes and origin, and has dimension1 and multiplicity1617.
Its complete Hilbert numerator over denominator (1-t) is retained in the
certificate. By symmetry, the generic torsion length at each diagonal-axis
prime is 1617/3=539. This is the intrinsic local finite-length torsion of
R/I at that prime, not the length of a chosen embedded primary component.
The subsequent [graded-truncation theorem](THREE_FACTOR_PRIMARY_DECOMPOSITION.md)
computes intrinsic origin torsion77 and supplies the full primary
decomposition; neither follows from multiplicity1617 alone.

## Evidence and retained failures

The [checker](verify_three_factor_associated_geometry.py),
[certificate](three_factor_associated_geometry_certificate.json), and
[receipt](three_factor_associated_geometry_verification.json) bind the
moment inputs, exact Smith ranks, star fiber, saturations, intersection,
socle tests and full Hilbert numerators. Replay with:

```sh
.venv/bin/python research/jacobian-followup/moment-vanishing-audit/verify_three_factor_associated_geometry.py
```

Use `--audit-existing` to audit retained Singular output and run the small
symbolic star-fiber checks.
Replay uses SymPy1.14.0 and Singular4.3.2. Each Singular run is capped at
50 seconds and 768MiB; the accepted row Smith checker has a35-second CPU
and 768MiB bound and completes using the written degree-10 identity.

The direct full `primdecGTZ` call timed out after50 seconds without a
component. The first full row Smith run stopped at its degree-10 CPU
bound, retaining degrees0..9. Both are terminal UNKNOWN computations,
not live jobs or evidence for a decomposition. Earlier row-resultant
explorations preserve an invalid first input and a second input with an
incorrectly copied resultant coefficient; neither is used. The corrected
third input is an exploratory fiber check, not the freeness proof.
The final all-parameter Smith proof supersedes the need to interpret that
fiber check as a generic or flatness certificate.

This result concerns centered products of three linear forms in the
universal Gram ring. Affine factors, constrained coefficient pullbacks,
positive-characteristic extensions and arbitrary ternary cubics retain
their separate obligations.
