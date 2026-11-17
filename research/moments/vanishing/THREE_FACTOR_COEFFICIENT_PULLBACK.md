# Three-factor moment schemes in actual Gaussian coefficients

Claim: `G3COEFFICIENT20261114`. Work over C. This is an internal written
geometric and invariant-theoretic proof with exact symbolic transcription
checks, without independent replay, external review or formal verification.

Let C be a generic 3-by-N matrix, with rows the coefficient vectors of
three centered linear Gaussian forms. Put

    S_N = C[entries of C], R = C[a,b,c,d,e,f],
    pi_N(C) = C*C^T = [[a,d,e],[d,b,f],[e,f,c]].

Write I_k=(M2,...,Mk) for the universal Gram ideals when k is even, and
I=I_infinity=I10. An ideal with a subscript N denotes extension to S_N.
The [universal primary theorem](THREE_FACTOR_PRIMARY_DECOMPOSITION.md)
gives ten primary ideals Q_1,...,Q_10 whose intersection is I.

## Results

The least initial moment orders generating the full coefficient ideal,
and the least initial pure prefixes detecting its zero set, are:

| Gaussian dimension | Ideal-generation cutoff | Set-theoretic pure cutoff |
|---|---:|---:|
| N=1 | 2 | 2 |
| N=2 | 6 | 6 |
| N>=3 | 10 | 8 |

The Gram map is faithfully flat exactly for N>=5. In this range the full
coefficient moment scheme has six minimal primes, degree720, reduced
degree48, and exact nilradical indices14 for I8_N and13 for I_N.
Its full primary decomposition transfers from R, with one dimension-six
splitting:

| Dimension | Primary components / associated primes | Generic torsion over Gram zero |
|---|---:|---|
| N=5 | 10 | One component, intrinsic length154 |
| N=6 | 11 | Two components, intrinsic length77 on each |
| N>=7 | 10 | One component, intrinsic length77 |

The three embedded axis components still have intrinsic generic torsion
length539. The six minimal components retain generic lengths24 for rows
and6 for stars. These torsion lengths are not chosen primary-component
quotient lengths.

“Gram zero” here is the positive-dimensional locus C*C^T=0, with its
reduced support when necessary. It is not the coefficient origin C=0.
For N>=5 the coefficient origin is not an associated prime, and its
finite-length local-cohomology submodule H^0_(entries C)(S_N/I_N) is zero.

The [subsequent binary factor theorem](BINARY_CUBIC_FACTOR_SCHEMES.md) now
completes N=2 with45 primary components, degree36, nilradical index9 and
intrinsic coefficient-origin torsion13. Full N=3,4 coefficient schemes and
further constraints remain open. The [N=3,4 refinement](SMALL_DIMENSION_COEFFICIENT_GEOMETRY.md) closes the minimal components and generic thickness, with nine components and degree720 in both dimensions, and identifies embedded quaternary null families with torsion at least1635. Full exceptional embedded structure remains open. N=1 is elementary and complete.
These results do not cover affine factors or change arbitrary cubic GMC.

## 1. Fiber dimensions and the flatness boundary

For k vectors in a nondegenerate d-dimensional quadratic space, the
rank-s stratum of the zero-Gram locus has dimension

    D(d,k,s) = s(d+k) - s(3s+1)/2,
    0 <= s <= min(k,floor(d/2)).                         (1)

It parameterizes an isotropic s-plane and a spanning ordered k-tuple in
it. It is irreducible unless d=2s>0, in which case it has two irreducible
components. For k>=2 these facts follow by setting p=0 in
Casabella--Sammartano's [Theorem3.5, v2](https://arxiv.org/html/2512.25058v2).
The k=0,1 cases are immediate. We use these rank-stratum facts, not the
paper's primality theorem for the different ideal containing only
off-diagonal orthogonality equations.

Fix a Gram matrix of rank r. Congruence in the three row indices puts it
in the form diag(I_r,0). Choose the first r orthonormal vectors, a Stiefel
variety of dimension rN-r(r+1)/2. The other k=3-r vectors lie in its
nondegenerate orthogonal complement of dimension d=N-r and have zero
Gram matrix. Thus fiber rank strata have dimensions

    rN-r(r+1)/2 + D(N-r,3-r,s).                          (2)

For N>=7, the maximum occurs at s=3-r and is3N-6. The difference from
each lower s is a positive linear function of N, already positive at7.
For N=5,6 the exhaustive four values r=0,1,2,3 give the same maximum
3N-6. Every irreducible component of a fiber, cut out by six equations
in affine3N-space, has dimension at least3N-6. Hence all fibers have pure
dimension3N-6 and are complete intersections.

The source is Cohen--Macaulay and the target is regular. Applying
[miracle flatness, Stacks00R4](https://stacks.math.columbia.edu/tag/00R4)
gives flatness. Every complex symmetric3-by-3 matrix is a Gram matrix
in dimension at least3, so the map is surjective and faithfully flat.

For N=3, the zero fiber has dimension4, exceeding the generic dimension3.
For N=4 it has dimension7, exceeding6. These jumps preclude flatness.
For N=1,2 the homomorphism R to S_N has a nonzero kernel given by rank
relations, which also precludes flatness over the domain R.

## 2. The nine nonorigin primes remain prime

It suffices to treat one row, star and axis. Flatness makes their
extensions complete intersections of codimensions3,3,5 respectively.
Such rings are unmixed, so a unique irreducible generic component and
generic reducedness suffice for primality.

* Row A imposes A.A=A.B=A.C=0. On A nonzero its isotropic base is
  irreducible and B,C each range freely in A-perp. This open set is smooth
  and irreducible of dimension3N-3. The complement A=0 has dimension2N,
  strictly smaller for N>=5, and cannot be a component.
* Star A imposes B.B=C.C=B.C=0, with A free. When B,C span an isotropic
  two-plane, the incidence space is irreducible for N>=5 and has dimension
  3N-3. The rank-one and rank-zero complements have dimensions2N and N.
  The three equations have independent differentials at a spanning pair.
* Axis A imposes B.B=C.C=B.C=A.B=A.C=0. For independent B,C, choose an
  isotropic two-plane W, a basis B,C in W, and A in W-perp. This irreducible
  vector-bundle construction has dimension3N-5. The other ranks have
  dimensions2N-1 and N, both strictly smaller. Generically choose A
  anisotropic outside W, possible since dim(W-perp/W)=N-4>=1. Then the
  three rows are independent, the full Gram differential is surjective,
  and the chosen five equations have independent differentials.

Permuting rows gives all nine primes. In particular the extension J_N of
the universal radical J is the intersection of six prime complete
intersections and is itself radical. Each minimal component has degree8.

In dimension five, a fixed nonzero rank-one Gram matrix does have a
geometrically split fiber. The total axis preimage is nevertheless
irreducible by the incidence construction above. The argument does not
assume that every fiber of the flat Gram map is irreducible.

## 3. The zero-Gram fiber in dimensions five and six

The ideal m*S_N, for m=(a,b,c,d,e,f), is a complete intersection of six
quadrics when N>=5. Its rank strata are given by (1) with k=3,d=N.

For N>=7, rank3 is the unique irreducible top stratum and has surjective
Gram differential. Lower strata have smaller dimension, so m*S_N is
prime. This also agrees with Kraft--Schwarz's
[Theorem9.1(4)](https://arxiv.org/html/1112.3634v1).

For N=6, rank3 has two components, corresponding to the two families of
maximal isotropic three-planes. Both are generically reduced, while all
lower strata have smaller dimension. Consequently m*S_6 is radical with
two minimal primes and no embedded primes.

For N=5, the single top stratum has rank2. Use Witt coordinates with
bilinear form pairing the first two coordinates with the next two, plus
one square. On the dense chart where the first two coordinate columns
have rank2, row operations put the matrix into

    [1 0 v11 v12 z1]
    [0 1 v21 v22 z2]
    [0 0 v31 v32 z3].

Five zero-Gram equations solve

    v11=-z1^2/2, v22=-z2^2/2, v21=-v12-z1*z2,
    v31=-z1*z3, v32=-z2*z3.

The last equation is z3^2=0. Thus m*S_5 has generic length2. Its support
is irreducible, and its complete-intersection quotient is unmixed;
therefore it is primary but nonreduced. Flatness did not imply a reduced
fiber. The direct chart proves the multiplicity, beyond a reducedness
criterion alone.

## 4. The transferred primary decompositions

Flat extension preserves finite intersections. The
[flat associated-prime formula, Stacks05C2](https://stacks.math.columbia.edu/tag/05C2)
applied to R/Q_i shows that its extended associated primes are exactly
those of the generic fiber over its one base associated prime. The nine
prime extensions just proved each give one primary extension. For the
origin component Q_0, the preceding zero-fiber analysis gives one
associated prime when N=5 or N>=7, and two when N=6.

Therefore for N=5 and N>=7, the ten ideals Q_i*S_N are already an
irredundant primary decomposition of I_N. The selected origin extension
is primary even in dimension five, where its radical is the reduced
zero-Gram support rather than the nonreduced ideal m*S_5.

For N=6 let p123 and p456 be the determinants of the first and last
three columns of C in Euclidean coordinates, and put

    h_plus=p123+i*p456, h_minus=p123-i*p456.

On the chart p123 nonzero, write C=U*[I3,M]. The zero-Gram equation gives
M=i*O with O orthogonal, so p456=-i*det(O)*p123. The two signs det(O)=+1
and-1 give dense charts in the two top components. Thus h_plus vanishes
on one generic component and h_minus on the other. Replace Q_0*S_6 by

    (Q_0*S_6):h_plus^infinity,
    (Q_0*S_6):h_minus^infinity.

They are the two isolated primary components of this extension, and
their intersection equals it. Along with the other nine ideals this is
an explicit irredundant eleven-component decomposition. All saturations
in the original ten component formulas commute with flat base change.

## 5. Lengths, Hilbert series and nilpotence

Let N_k(t) be the retained full Hilbert numerator of R/I_k over
denominator(1-t)^6, where each Gram entry has degree1. Gram entries have
coefficient degree2. Tensoring a finite graded free resolution with S_N
therefore gives, for N>=5,

    H(S_N/I_k*S_N;t) = N_k(t^2)/(1-t)^(3N).              (3)

The numerator has a third-order zero at1. Its leading coefficient gives
degree2^3*90=720 for both k=8 and k=10. The radical has six degree-eight
components and reduced degree48.

Since extension is faithfully flat and J_N is radical,

    (J_N)^j contained in I_k*S_N
        iff J^j contained in I_k.

Thus the sharp nilradical indices remain14 and13. In particular the
explicit universal lower witnesses still survive.

At each row, star or axis prime, the extended base prime is prime, so
the local closed fiber is a field. Flatness preserves finite local
lengths24,6 and539 there. For the intrinsic base-origin torsion T of
length77, flat extension and the associated-prime list show that its
localization is exactly the corresponding zero-Gram local torsion. The
local closed-fiber length is2 in dimension five,1 at each of the two
dimension-six primes, and1 in dimensions at least seven. A composition
series of T gives lengths154,77/77,77 as displayed above.

No associated prime in this range is the coefficient origin: the largest
codimension is6, whereas S_N has dimension3N>=15. Hence there is no
finite-length submodule supported at that point.

## 6. Sharp generation even without flatness

For N>=3 the six Gram entries are algebraically independent and generate
the orthogonal invariant ring S_N^O(N), by the orthogonal first
fundamental theorem, recalled in Kraft--Schwarz Section9. Averaging over
the compact real orthogonal group gives an R-linear Reynolds projection
S_N to R. Consequently every extended ideal contracts to itself.
This argument does not require flatness.

The universal M10 outside I8 therefore remains outside I8*S_N for every
N>=3. The five universal generators pull back in every dimension, so10
is the sharp ideal cutoff for all N>=3, including the nonflat N=3,4
cases. The earlier rank-three Gram witness gives sharp set cutoff8 in
every N>=3 by padding its coefficient matrix with zero columns.

For N=2, the rank relation is

    det(Gram)=A+2E-s=0,
    A=abc, E=def, s=af^2+be^2+cd^2.

Use the symmetric invariant coordinates from the
[global ideal proof](THREE_FACTOR_GLOBAL_MOMENT_IDEAL.md), including
q=(af^2)(be^2)+(af^2)(cd^2)+(be^2)(cd^2). Modulo M2,M4 and this determinant,

    A=-4E, s=-2E, q=7E^2,
    M6=-4082400E^3,
    M8=11735539200E^4,
    M10=-9807557760000E^5.

Hence M8,M10 belong to (M2,M4,M6,det). The universal all-order generation
theorem then gives I_infinity,N=I6,N in dimension two.

For sharpness, take all off-diagonal Gram entries1 and let its three
diagonal entries be the roots of z^3+2z^2+7z+4. Their elementary symmetric
values are s=-2,q=7,A=-4; thus the determinant is zero, M2=M4=0 and
M6=-4082400. The cubic discriminant is-728, so the diagonal entries are
distinct. This Gram matrix has rank2: rank1 with all off-diagonal entries1
would force all diagonal entries1. Every complex symmetric rank-two
matrix is realized by three linear forms in two Gaussian coordinates.
All odd moments vanish by homogeneity. Thus the first five pure moments
can vanish with the sixth nonzero, proving both sharp cutoffs6.

For N=1 write P=(uvw)*X^3. Then M2=15*(uvw)^2 generates every positive
moment, and the coefficient ideal has the explicit primary decomposition

    ((uvw)^2) = (u^2) intersect (v^2) intersect (w^2).

Its degree is6 and its nilradical index2. This closes dimension one.

## Reproduction and remaining scope

The [checker](verify_three_factor_coefficient_pullback.py),
[certificate](three_factor_coefficient_pullback_certificate.json), and
[receipt](three_factor_coefficient_pullback_verification.json) bind the
earlier exact inputs. They verify the symbolic stratum identities and
all-N inequalities, exceptional Witt charts, separator signs, rank-two
identities, witness discriminant and Hilbert degree. They do not replace
the written geometric, invariant-theoretic and flat-transfer arguments.

```sh
.venv/bin/python research/moments/vanishing/verify_three_factor_coefficient_pullback.py
```

The [source audit](three_factor_coefficient_source_audit.json) records the
actual selected sections read, the October30 v2 identity of the frames
paper, the older reduced-null-cone source, and the two Stacks lemmas.
No producer code or formalization was executed. Full coefficient
primary/thickness descriptions in N=3,4, further constrained coefficient
subfamilies and positive characteristic remain separate obligations. The
subsequent [N=2 theorem](BINARY_CUBIC_FACTOR_SCHEMES.md) closes the binary
case by a different, nonflat constructive argument.
The [N=3,4 refinement](SMALL_DIMENSION_COEFFICIENT_GEOMETRY.md) closes their minimal components and generic thickness, identifies quaternary embedded null families and derives exact generic-section Hilbert data. Full exceptional primary structure and nilradical indices remain open. Nothing here extends the centered-product theorem to general cubics.
