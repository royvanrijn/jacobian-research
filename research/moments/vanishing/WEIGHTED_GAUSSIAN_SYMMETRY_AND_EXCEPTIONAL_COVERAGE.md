# Weighted Gaussian symmetry transfer and exact exceptional coverage

16 November 2026. Internal written proof. Independent replay, formal
verification and external review are not claimed; external review remains
deferred. BOTH full exceptional Mathieu--Zhao cases remain **OPEN**.
All earlier goals/results, capacity and automorphism constraints remain.

The user's Gaussian-to-Mathieu symmetry bridge is valid. The stronger
no-linear-term corollary is valid as well, in the supplied nonnegative
covariant-weight-two Gaussian architecture. Coverage can be determined
exactly for the canonical homogeneous exceptional lifts: it recovers a
STRICT SUBSET of the already closed torus-safe families and no surviving
mixed phase. Allowing an arbitrary invertible LINEAR Gaussian coordinate
change does not enlarge that coverage beyond the criterion below.
Nonlinear preprocessing requires separate hypotheses and is not implicit
in the weighted theorem.

## 1. The completed Gaussian input and its actual witness

We use Sections1,3,4 of the completed
[nonnegative weighted Gaussian theorem](SEMIPOSITIVE_WEIGHTED_GAUSSIAN_REDUCTION_AND_SAFE_NORMAL_FORMS.md),
read in full. Its hypotheses are a fixed nonsingular centered Gaussian
quadratic Q0 and supplied nonnegative rational LINEAR coordinate weights
for which BOTH Q0 and the phase P have weight2. Coefficients have weight0.
In adapted linear coordinates its normal form is

    Q0=Z dot W+(1/2)T^T B0 T,
    P=W dot H(Z)+F(Z,T),
    g=t H(g), g in tR[[t]]^a.

The residual weights satisfy 0<w(T)<2, and F has neither T-constant
nor T-linear terms. Under ALL-POSITIVE PURITY, the original marked
periods of the fixed linear position markers are

    Phi[Z_i]=g_i(t),

and ALL-original-marker safety is equivalent to every g_i being a
polynomial. This input holds over EVERY commutative Q-algebra R with
the specified nonsingular Gaussian and weighted coordinates defined.
The Gaussian construction, finite Morse proof, producer and receipts
are completed inputs; none is changed or replayed here.

## 2. Symmetry kills the required unsafe witness

Let a finite linear symmetry group preserve the centered Gaussian and
P. If its representation on linear functions has no invariant vector,
then for EVERY linear function ell and EVERY m>=0,

    E[ell P^m]=E[(Reynolds ell) P^m]=0.                (1)

This includes linear functions in ANY adapted Gaussian coordinates,
since they remain linear combinations of original coordinates.
In the weighted class, purity now gives g_i=0 identically, so every
fixed linear safety test is polynomial. The Gaussian theorem proves
FULL safety for ALL original polynomial markers.

For the binary exceptional lifts, the symmetry -I alone proves (1):
P is even and ell P^m is odd. There is no invariant linear marker.
Consequently every PURE exceptional lift that admits the supplied
weighting is Gaussian safe. The exact
[radial/marker transfer](GAUSSIAN_RADIAL_TRANSFER_AND_EXCEPTIONAL_UNIFORMITY.md)
then makes the original quotient phase safe for every fixed original
regular invariant marker. This is an all-order, all-degree implication,
not a finite moment inference. With the data defined over R it also
holds over every Q-algebra R by the same symmetry and exact unit radial
factors. No new normalization or ordinary-ideal calculation is needed.

## 3. The simpler no-linear-term corollary

In the displayed normal form, the ENTIRE ordinary linear part of P is

    sum_i W_i H_i(0).

There are no other linear terms: F has no T-constant/T-linear term,
and a lone zero-weight Z term cannot have weight2. Thus if P has no
ordinary linear term, H(0)=0. A linear coordinate change preserves this
hypothesis. The fixed point is then exactly g=0: zero is a solution and
the t-adic solution is unique. Directly, a nonzero solution's first
nonzero t-order would increase in tH(g), a contradiction.

It follows that EVERY PURE phase in this weighted class with no linear
term is safe, in EVERY dimension/degree and over EVERY Q-algebra R.
No symmetry is needed for this corollary. Purity remains a premise;
a nonzero quadratic Gaussian phase is not made safe by omitting linear
terms. This is a corollary of the completed Gaussian theorem, not a
universal statement about Gaussian phases outside its architecture.

## 4. Conformal linear Euler operators on matrix space

For coverage, work over C with V=Mat_2, Q=det(M). Quaternion coordinates
identify Q with the nonsingular Gaussian radial quadratic up to a
nonzero constant. A supplied weighting after any invertible linear
coordinate change gives a semisimple linear vector field E satisfying

    E(Q)=2Q, E(Ptilde)=2Ptilde.                       (2)

Conversely, a diagonalizable E with nonnegative rational eigenvalues
and (2) supplies exactly such coordinates and weights. No affine
translation or nonlinear change is included.

Every linear conformal field with E(Q)=2Q has a UNIQUE decomposition

    E=Eu+L_A+R_B, A,B in sl2,
    Eu(M)=M, L_A(M)=AM, R_B(M)=MB.                   (3)

Indeed Eu(Q)=2Q and the infinitesimal determinant-preserving fields
are precisely sl2_left direct_sum sl2_right: the six displayed fields
are independent and the orthogonal Lie algebra in four variables has
dimension6. This is ordinary linear algebra, also checked symbolically.

Let Gamma=2O or2I act on the RIGHT and let Ptilde be invariant. Averaging
conjugates of E over Gamma leaves its action on Ptilde and Q unchanged.
It kills the right traceless part: an invariant traceless 2x2 matrix
commutes with the irreducible standard Gamma representation and is0.
Thus

    Ebar=Eu+L_A,
    L_A(Ptilde_j)=(2-j)Ptilde_j                       (4)

for every ordinary homogeneous degree-j component.

For these two exceptional groups there is also rigidity of E itself.
The adjoint Gamma representation on sl2 is irreducible. To see this
without a character table, an invariant line spanned by a nonzero
nilpotent matrix would give a Gamma-invariant kernel line in C2 and
force a finite unitary Gamma cyclic. A semisimple matrix line would
put Gamma in a torus normalizer, hence make its projective group cyclic
or dihedral. Both contradict the octahedral/icosahedral groups. Complete
reducibility then rules out an invariant plane as well.

The right infinitesimal stabilizer of Ptilde is a Gamma-stable subspace
of sl2. If nonzero, it is all sl2, so Ptilde is a polynomial in Q.
Here right-SL2 invariants are C[Q]: on Q!=0, matrices of a fixed nonzero
determinant form one right orbit, and regularity excludes finite poles
in Q. But (2) makes a nonzero radial homogeneous phase possible only
in degree2. Therefore for a nonzero homogeneous Ptilde of degree D>2,
its right infinitesimal stabilizer is zero. Comparing E and Ebar in
(2),(4) forces B=0 in (3).

## 5. Exact coverage of the canonical homogeneous lift

Let f be an original exceptional quotient phase and Ptilde its canonical
homogeneous determinant lift of prescribed EVEN degree D>2. The lift
is unique, right-Gamma invariant, and restricts to f on Q=1. Assume this lift is nonzero. We claim:

**Ptilde admits the nonnegative rational weighting of Section1, after
ANY invertible linear Gaussian coordinate change, if and only if some
left SL2 basis makes it a SINGLE left-torus weight k with**

    |k|=D or |k|=D-2.                                (5)

Necessity: by Section4, E=Eu+L_A. Because E is semisimple, A is
semisimple; and because 2-D!=0, A!=0. Conjugate it on the left to
A=diag(a,-a). Equation (2) then says

    a k=2-D                                         (6)

on EVERY nonzero monomial of Ptilde, so there is just one left weight k.
The four Gaussian coordinate weights of E are 1+a,1+a,1-a,1-a.
Nonnegativity forces a real and |a|<=1. Hence |k|>=D-2. A degree-D
matrix polynomial has weights -D,-D+2,...,D. This proves (5).

Sufficiency: if (5) holds, let H be that left-torus Euler field and set

    E=Eu+((2-D)/k)H.

It gives Q weight2 and Ptilde weight2. Its coordinate weights 1+/-
(2-D)/k are nonnegative rational numbers, so it is an allowed supplied
weighting. At |k|=D they are 2/D and 2-2/D, strictly positive. At
|k|=D-2 they are0 and2, the zero-weight boundary of the architecture.
This proves both directions, in ALL degrees. No degree enumeration
or weight-recognition producer is used.

The degree2 invariant space is C Q, by the absence of invariant binary
quadratics; its nonzero members restrict to constants and are not pure.
Degree0 contributes only the zero phase to weight2. A nonzero pure
exceptional phase has higher degree, where (5) applies. In fact a
nonzero weighted homogeneous lift with D>2 is automatically pure and
safe, without separately assuming purity: k!=0 and left-torus invariance
kills every unit moment and every sufficiently high marked moment.

## 6. Comparison with the closed families and surviving mixed phases

Condition (5) lies in the ORIGINAL left nullcone. It therefore recovers
already proved torus safety. The |k|=D case is a one-row binary invariant;
the |k|=D-2 case has row bidegrees (D-1,1) or (1,D-1). Apparent mixed
weights in an unadapted basis can occur after a left SL2 rotation, but
they remain in this already closed family. There is NO surviving mixed
phase, with closed left orbit outside the nullcone, admitting this
canonical-lift weighting. The symmetry bridge adds no such cases.

The criterion is for the PRESCRIBED lift degree D. Padding to a higher
homogeneous degree by powers of Q can destroy a weighting; it cannot
rescue a surviving mixed phase.

It does not even cover every already-safe one-sided phase. For example,
the retained primitive PURE SAFE controls f_O=F8+G12 and f_I=F12+G20 have
canonical lifts

    Q^2 F8+G12 (D=12), Q^4 F12+G20 (D=20).

Their two distinct Peter--Weyl components are left highest-weight
vectors with weights (8,12) or(12,20). If a traceless A had the same
nonzero eigenvalue on both, its lowering part would vanish on each
highest vector, and its diagonal part would have to give equal unequal
positive multiples. This is impossible. Equation (4) therefore rules
out ANY conformal linear Euler field satisfying (2), even without the
nonnegative requirement. These are exact existing all-order controls,
not new finite degree closures.

More generally, for an INHOMOGENEOUS invariant Gaussian polynomial
admitting the weighting, (4) gives a scalar multiple of Q in degree2
and strictly one-sided left weights in every higher nonzero component.
There is no invariant linear part or other invariant quadratic, and no
constant weight2 part. If there is a nonradial component, the same
stabilizer argument forces B=0 and a nonzero real left eigenvalue a;
all higher charges (2-j)/a have the same sign. The quadratic term is
cQ. On Q=1 the phase is consequently c plus a one-sided phase. Original
Haar purity forces c=0. This extension also adds no genuinely mixed
family; it must not be confused with the radial moment identity, which
uses a FIXED homogeneous lift.

If one first applies an original Haar-preserving algebraic automorphism
and THEN tests the canonical lift, the covered phase belongs to the
already-safe automorphism saturation of the nullcone. Degree may change,
but safety is preserved. No theorem says every pure phase lies in that
saturation. Arbitrary nonlinear Gaussian preprocessing is not licensed
by the weighted theorem or by this coverage calculation.

## 7. What remains

The user's symmetry implication and no-linear-term corollary are proved
at their exact supplied-weight scope. The missing coverage question has
a negative answer for surviving mixed exceptional phases under that
LINEAR Gaussian-weight architecture. This is a coverage theorem, not a
new proof of the full pure-to-safe inclusion.

The full target stays intrinsic boundary splitting on EVERY pure original
phase, or an ACTUAL pure exceptional phase with ONE fixed original bad
marker. Both cases remain OPEN. Capacity/support/aut/primitive inputs
remain, including at least6/8 weights for a genuine failure after every
allowed automorphism/torus. The completed family transfer, previous
Long/projectable-lift obstructions and their receipts are preserved.
External review remains deferred; no degree/prefix campaign or additional
decidability result is substituted for the full obligation.
