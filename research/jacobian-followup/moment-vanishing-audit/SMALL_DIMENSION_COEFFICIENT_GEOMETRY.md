# Three-factor coefficient geometry in dimensions three and four

Claim `G34MINIMAL20261114`, characteristic zero over C. Let C have rows
A,B,C, let pi(C)=C*C^T, and let I_N be the full pure-moment ideal of their
centered Gaussian product. The [all-order theorem](THREE_FACTOR_GLOBAL_MOMENT_IDEAL.md)
gives I_N=I10,N. Write I=I10 in the six Gram entries

    G = [[a,d,e],[d,b,f],[e,f,c]],   R=C[a,b,c,d,e,f].

The [universal primary theorem](THREE_FACTOR_PRIMARY_DECOMPOSITION.md)
and [coefficient transfer](THREE_FACTOR_COEFFICIENT_PULLBACK.md) are inputs.
The following closes the minimal-component and generic-thickness questions
in N=3,4. It also identifies some embedded structure in N=4. Full primary
decompositions and global nilradical indices in these dimensions remain open.

| Dimension | Minimal component types | Number | Component degree | Generic moment length |
|---|---|---:|---:|---:|
| 3 | One factor zero | 3 | 1 | 48 |
| 3 | Nonzero isotropic factor, orthogonal to the other two; closure | 3 | 6 | 24 |
| 3 | Two factors parallel and isotropic | 3 | 4 | 12 |
| 4 | Isotropic factor orthogonal to the other two | 3 | 8 | 24 |
| 4 | Two factors spanning a maximal isotropic plane, either family; closure | 6 | 4 | 6 |

All components have codimension3. The total scheme degree is720 in both
dimensions; the reduced degrees are33 and48. These are degrees of the
actual coefficient schemes, rather than degrees inferred by nonflat transfer.

## 1. Exact flat loci

For a fixed Gram matrix of rank r, put k=3-r and d=N-r. The rank-s stratum
of k null vectors in a nondegenerate d-space has dimension

    D(d,k,s)=s*(d+k)-s*(3*s+1)/2,
    0<=s<=min(k,floor(d/2)).

Adding the frame dimension r*N-r*(r+1)/2 and taking the maximum gives:

| N | Gram rank0 | Rank1 | Rank2 | Rank3 | Expected fiber dimension |
|---|---:|---:|---:|---:|---:|
| 3 | 4 | 4 | 3 | 3 | 3 |
| 4 | 7 | 6 | 6 | 6 | 6 |

The strata formula is the same incidence calculation used in the earlier
coefficient note. The rank list is exhaustive. The source and target are
smooth, and equal fiber dimensions on the indicated open sets give flatness
by the earlier miracle-flatness argument. The larger remaining fibers
exclude flatness there. Thus pi_3 is flat exactly over rank>=2, and pi_4
is flat exactly away from the zero Gram matrix. Neither map is globally flat.

## 2. The nine ternary minimal components

The universal radical has row ideals (a,d,e) and permutations, and star
ideals (b,c,f) and permutations. Their inverse images exhaust the null set.

For a row ideal in N=3, the equations are

    A*A=A*B=A*C=0.

There are exactly two components: Z_A=(A=0), and R_A, the closure with
A nonzero isotropic and B,C in A-perp. Both have dimension6. The latter
is irreducible: the nonzero isotropic cone is irreducible and its two
orthogonal-vector fibers form a vector bundle. The three equations are
a complete intersection of quadrics, of degree8. At generic Z_A, B,C span
a nondegenerate plane; two equations eliminate two coordinates of A,
and the remaining equation is a nonzero multiple of the square of its
normal coordinate. Hence Z_A has multiplicity2 in this row pullback.
At generic R_A the three differentials are independent, so its multiplicity
is1. Consequently deg R_A=8-2=6.

For a star ideal the equations are B*B=B*C=C*C=0. An isotropic subspace
in dimension three has dimension at most one. Thus the support S_A has
B,C on one isotropic line, with A arbitrary; it is irreducible of dimension6.
At a nonzero pair, use a Witt basis Z,T,W with Z*W=T*T=1. The transverse
pair B=Z+xW, C=qZ+zT+yW has Gram ideal

    (2x,y+q*x,2q*y+z^2)=(x,y,z^2).

This complete intersection has generic length2. The original three
quadrics have degree8, so deg S_A=4. The nine components Z_A,R_A,S_A
and their permutations are distinct and exhaust the null set.

At generic Z_A the Gram matrix has rank2, where pi_3 is flat. The universal
row algebra has length24, and its local closed fiber has length2 as above.
The moment algebra therefore has length48. In normal A-coordinates its
moments have degrees2,4,6,... . For B=X,C=Y the first three have no
nonzero projective common zero, checked exactly after putting
u=x^2,v=y^2,w=z^2. This complete-intersection condition holds on a nonempty
open parameter set. A (2,4,6) complete intersection has length48; hence
the remaining moment equations cannot shorten the generic algebra. Its
maximal-ideal Hilbert function is

    [1,3,5,7,8,8,7,5,3,1],  nilpotence index10.

Generic R_A and S_A are dense orbits of O(3) times the three independent
factor scalings. Exact transverse slices are, in the same Witt basis,

    R_A: A=Z+xW, B=T+yW, C=Z+T+zW;
    S_A: A=Z/2+W, B=Z+xW, C=Z+zT+yW.

The six group tangents have rank6, and adjoining the three slice tangents
gives rank9, checked symbolically. The action-and-slice maps are therefore
etale locally at these points. Since the moment ideal is preserved by the
group, the slice algebras are the generic transverse algebras on the dense
orbits; this argument supplies genericity beyond a numeric sample.

The R_A slice is saturated by C*C=1+2z, and the S_A slice by 2*A*C=y+2.
These are units at the selected point. They remove other components of
the global slices. Exact standard bases and powers of (x,y,z) give

    R_A: length24, [1,2,3,4,4,4,3,2,1], index9;
    S_A: length12, [1,2,2,2,2,2,1], index7.

The terminal zero remainder of the indicated maximal-ideal power proves
that the saturated finite algebra is supported at the origin. In particular,
the unlocalized S_A slice length60 is not its generic local length12.
Summing the nine minimal contributions gives

    deg(I_3)=3*(48*1+24*6+12*4)=720,
    deg(rad I_3)=3*(1+6+4)=33.

## 3. Quaternary minimal components and the nonzero-Gram open set

In N=4 a row incidence has dimension9 over the nonzero isotropic cone;
the A=0 locus has smaller dimension8. Thus the three row equations form
a prime complete intersection of degree8. Their generic Gram rank is2.

A generic null pair B,C spans a maximal isotropic two-plane. Such planes
have two irreducible families. The two incidence closures are distinct,
each of dimension9, with their dependent-pair boundary of dimension8.
The three star quadrics form a complete intersection, generically smooth
on each family. The two families are exchanged by an orthogonal isometry,
so each component has degree4. Generic Gram rank is2. Flatness there and
reduced generic closed fibers transfer the universal lengths24 and6.
This proves the table and

    deg(I_4)=3*8*24+6*4*6=720,
    deg(rad I_4)=3*8+6*4=48.

The complete universal associated-prime list can be transferred on the
nonzero-Gram open set, where pi_4 is flat. Its six minimal primes give
the nine components above. Each of its three diagonal axes has one
irreducible inverse-image support on this open set: A is nonisotropic,
and B,C lie on an isotropic line in the nondegenerate three-space A-perp.
The local closed fiber again has length2, by (x,y,z^2). Therefore the
intrinsic generic torsion at each of these three horizontal axis components
has length2*539=1078. There are exactly twelve associated primes on this
open set. This statement concerns intrinsic torsion, not the colength of
an arbitrarily chosen embedded primary component.

## 4. The generic rank-two zero-Gram section

Use the four-dimensional Witt form with two paired coordinate planes,
and write the coefficient matrix as [U,V]. On a dense chart, a rational
orthogonal change puts

    U=[[1,0],[0,1],[s,t]].

The five entries a,b,d,e,f and one antisymmetric coordinate freely
parametrize V. The remaining Gram entry satisfies exactly

    L=c+s^2*a+t^2*b+2*s*t*d-2*s*e-2*t*f=0.

The original chart also has a freely varying invertible two-by-two matrix
from the first two rows of U. Thus its ring is a smooth faithfully flat
extension of B=R[s,t]/(L). This is an exact chart, not a tangent approximation.

The row relation is lambda=(-s,-t,1). When every entry of lambda is nonzero,
independent row scalings identify B with a Laurent extension of the fixed
five-variable section

    H=I modulo ell,  ell=a+b+c+2d+2e+2f.

The moment generators are semi-invariant under these scalings. Flatness
of pi_4 away from Gram0 descends through the smooth chart to B. For each
universal minimal prime its extension to B on s*t nonzero is prime:
the single relation solves one remaining Gram coordinate with unit
coefficient. Each diagonal axis has empty inverse image away from Gram0,
because L on that axis is a unit times its remaining diagonal entry.
The flat associated-prime formula therefore shows that H has exactly
six associated primes away from its origin: its six restricted minimal
Gram primes. No other off-origin prime is silently omitted here.

Let M=R/I and T=H^0_m(M), with m=(a,b,c,d,e,f). The universal theorem
gives length(T)=77. Since ell avoids all nine non-origin associated
primes of M, the entire kernel K=(0:_M ell) is contained in T.
Using the retained exact saturation standard basis, a monomial-difference
basis of T constructs all77 multiplication columns over Q. Each column
is reduced in the original full moment ideal. Grading-wise exact ranks give

    H_K=10*t^20+5*t^21+3*t^30,  length(K)=18.

The basis is complete because its leading monomials are precisely the
finite difference between the two retained leading ideals. Triangular
normal forms give unique coordinates; there is no rank estimate from
floating point or sampling. Also

    H_(T/ell*T)=9*t^18+6*t^19+3*t^20.

For a stronger torsion bound let U0 be the intersection of the six
isolated universal primary ideals and put D=U0/I. This module has dimension1,
multiplicity3*539=1617 and finite origin torsion T. Multiplication by ell
has the same kernel K on D as on M. Since ell avoids the six associated
primes of R/U0, it is a nonzerodivisor there. Hence

    D intersect ell*M = ell*D,
    D/ell*D injects into M/ell*M,
    length(D/ell*D)=1617+18=1635.

The injected module has finite length. Thus the origin is an associated
prime of H, and its intrinsic torsion has length **at least1635**. Combined
with the off-origin classification this proves exactly seven associated
primes for the generic section. Both families of rank-two zero-Gram
quaternary matrices have such a dense chart; their generic row relations
have no zero entry. Each family is therefore an embedded associated
component of I_4. Together with the twelve nonzero-Gram primes this proves
at least fourteen associated primes globally. It does not exclude further
primes on the remaining zero-Gram strata.

Finally, the exact multiplication sequence determines the section Hilbert
series. If N_I is the retained numerator of R/I with denominator(1-t)^6,
then the five-variable section numerator is

    N_H=N_I+t*(1-t)^5*(10*t^20+5*t^21+3*t^30).

Its full coefficient list is retained in the certificate. The correction
is essential: ell is a zero divisor. We have not proved depth(R/U0)>=2;
therefore1635 is not promoted to the exact origin torsion length. Possible
additional local cohomology of the isolated-component intersection remains
a separate obligation.

## Evidence and remaining gates

The [checker](verify_small_dimension_minimal_geometry.py),
[certificate](small_dimension_minimal_geometry_certificate.json) and
[receipt](small_dimension_minimal_geometry_verification.json) bind the exact
inputs, logs, group Jacobians, finite rank list, slice Wick substitutions,
77-dimensional multiplication matrix, section numerator and degree sums.
The chart and associated-prime arguments above are written proofs; the
checker does not claim an independent verification of algebraic geometry.

```sh
.venv/bin/python research/jacobian-followup/moment-vanishing-audit/verify_small_dimension_minimal_geometry.py
```

Accepted Singular runs are the primitive-integer localized ternary
maximal-ideal calculation and the finite origin multiplication map.
The initial rational-expression slice input produced an exponent-precedence
error and is retained as invalid. The valid unlocalized pair slice and
the subsequent localization are both retained, with their distinct scopes.
Direct five-variable section basis, full six-variable principal colon,
and Hilbert-driven section basis calls each reached their50-second limit
without producing the requested basis. They are terminal computational
UNKNOWN, not evidence for a decomposition or for the exact torsion length.
A separate minimal-resolution/depth test on the retained isolated-component
intersection reached its20-second limit before producing a resolution.
It is also terminal UNKNOWN; it supplies no depth bound.

The [subsequent generic-section theorem](QUATERNARY_NULL_SECTION_PRIMARY.md) now gives all seven primary ideals constructively, with origin powers35, regularity<=34 and intrinsic torsion in[1635,2576]. This uses a regularity bound and rigorous modular upper-bound lifting; the depth and saturation timeouts above remain UNKNOWN.

Next gates are the full ternary embedded structure; exact quaternary
zero-Gram torsion and the remaining associated primes where the rank-two
row relation has a zero entry, or coefficient rank is at most one; full
N=3,4 primary decompositions and global nilradical indices. No new external
theorem was imported in this step. No independent replay, external review
or formal verification is asserted. Affine factors, arbitrary cubics and
lower-power GVC remain separate open problems.
