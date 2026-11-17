# Time-zero saturation and reconstruction by one principal-part map

16 November2026. Internal all-dimension, all-degree proof over EVERY
commutative Q-algebra. The classical Koszul argument is elementary; no
novelty or independent/formal/external assurance is claimed. The
[checker](verify_gaussian_origin_and_global_lattices.py),
[receipt](gaussian_origin_global_lattice_verification.json) and
[source audit](gaussian_origin_global_lattice_source_audit.json) retain
finite identities, not a bounded substitute for this proof.

## 1. Corrected exact Gaussian source theorem

Let R be any commutative Q-algebra, P in R[X_1,...,X_n], and set

    V=R[t,X], nabla_i=partial_i-X_i+t P_i,
    H=V/sum_i nabla_i V, D_t[A]=[partial_t A+P A].

All statements also hold for a fixed nonsingular symmetric rational
precision H0, replacing X by H0 X and covariance by H0^-1.
The actual functional is Phi[A]=E[A exp(tP)] in R[[t]].

Multiplication by t is INJECTIVE on the ENTIRE exact source H. This
strengthens the retained completion theorem and corrects its old sentence
that the actual kernel might contain t-torsion. That possibility applies
to general holonomic sources, but not to this exact Gaussian source at
its normalization origin. Punctual source torsion at OTHER points remains
possible and must be retained; section4 gives a genuine Gaussian example.

Here is the proof. Form the homological Koszul complex C of the commuting
R[t]-linear operators nabla_i on V, with H_0(C)=H. Its terms have
injective multiplication by t, and C/tC is the unperturbed Gaussian
Koszul complex. The polynomial automorphism U=exp(-Delta/2), finite on
each polynomial, satisfies

    partial_i-X_i = -U X_i U^-1.

Thus C/tC is isomorphic to the Koszul complex of multiplication by
-X_1,...,-X_n on R[X]. These variables form a regular sequence over ANY
R, so H_j(C/tC)=0 for j>0. For the exact sequence of complexes
0->C --t-> C -> C/tC ->0, the homology segment

    H_1(C/tC) -> H_0(C) --t-> H_0(C)

proves injectivity. There is no Noetherian or field hypothesis.
For general precision, use U=exp(-Delta_(H0^-1)/2) and an invertible
constant linear change among the Koszul operators.

The exact classical inputs are
[Stacks Lemma062F](https://stacks.math.columbia.edu/tag/062F)
(regular implies Koszul-regular over an arbitrary ring) and the
[long homology sequence](https://stacks.math.columbia.edu/tag/010V).
The explicit first-syzygy contraction below also proves the needed
positive Koszul degree directly without a general solver.

The rest of the corrected completion theorem is retained in full.
Modulo t^k, multiplication by exp(tP) is finite and invertible and
conjugates nabla_i to partial_i-X_i. Gaussian reduction gives

    H/t^k H -> R[t]/t^k, [A] -> Phi[A] mod t^k,            (1)

an isomorphism for every k>=1. Class1 maps to the unit E exp(tP), not
necessarily to constant1. These are structure-sheaf quotients, not
Weyl-module quotients. Completing gives R[[t]] with ordinary time
differentiation. A D_t-map to R[[t]] is determined by its t=0 Gaussian
functional and hence is a scalar multiple of Phi, exactly as in the
[sealed original proof](GAUSSIAN_FORMAL_COMPLETION_AND_POLYNOMIAL_HOM.md).

In particular

    K=ker Phi=intersection_k t^k H,       tK=K.

The new injectivity shows t acts BIJECTIVELY on K. The largest
t-surjective submodule is still exactly K, now with no t-torsion.
This does not force K=0: the source is finite over the Weyl algebra,
not necessarily finite over R[t]. The old nonfaithful-completion and
rank controls remain valid. The source/actual generic ranks need not be1.

## 2. Finite cancellation of a putative t-torsion certificate

Suppose finite polynomial amplitudes satisfy

    t f=sum_i nabla_i g_i.                                (2)

Reducing at t=0 gives a first Gaussian Koszul syzygy. It can be removed
EXPLICITLY. In standard precision put w_i=U^-1 g_i(0,X). Then
sum X_i w_i=0. For its homogeneous X-degree-d part define

    h_ij^(d)=(partial_i w_j^(d)-partial_j w_i^(d))/(d+1).

Differentiating sum X_i w_i^(d)=0 and using Euler homogeneity gives
w_i^(d)=-sum_j X_j h_ij^(d). Consequently the skew matrix

    h_ij=U sum_d h_ij^(d)

satisfies g_i(0,X)=sum_j nabla_j(0)h_ij. Replace

    g_i by g_i-sum_j nabla_j(t)h_ij.

The new g_i have zero t-constant coefficient. Their total divergence
is unchanged because the nabla commute and h_ij is skew. Divide each
by t in the polynomial ring. Equation(2) now gives a finite original
identity f=sum_i nabla_i(g_i/t). Every division by d+1 is by a rational
unit, so this cancellation works with arbitrary zero divisors in R.

For rational precision H0 with covariance C0, first replace the syzygy
vector g by H0*g, use U=exp(-Delta_C0/2), then carry its skew matrix
back by h -> C0*h*C0^T. The focused implementation checks these finite
identities in standard, split and nilpotent coefficient controls.

## 3. One universal point quotient restores the unlocalized source

Let delta_(0,R)=R[t,t^-1]/R[t], with the usual derivative. The target
R((t))/R[[t]] is canonically this module of FINITE negative principal
parts. Since t is injective on H, regard H as a submodule of H[t^-1].
The formal period extends to R((t)); take its principal part:

    pp_P:H[t^-1] -> delta_(0,R).

Then the EXACT sequence is

    0 -> H -> H[t^-1] --pp_P-> delta_(0,R) ->0.            (3)

Proof of the kernel: write an element as t^-k h. Zero principal part
means Phi(h) is divisible by t^k. The injectivity in(1) implies h in
t^k H, so t^-k h belongs to H. For surjectivity the actual unit period
Phi[1] has constant coefficient1. Its multiples by finite Laurent
polynomials give every desired finite principal part by triangular
coefficient inversion. If R!=0, (3) does not split even as R[t]-modules:
the middle has no t-torsion, so every map from the t-torsion target back
to it is zero.

The ACTUAL image M=Phi(H) also satisfies M/t^kM=R[t]/t^k. Indeed
K lies in every t^kH. Therefore the same proof gives

    M=M[t^-1] intersect R[[t]],
    0 -> M -> M[t^-1] --pp-> delta_(0,R) ->0.              (4)

These statements hold for nonpure phases, arbitrary nilpotents and
non-Noetherian targets. They are at the NORMALIZATION ORIGIN only.
The polynomial constant phase, for example, retains its nontrivial time
connection; principal parts are taken after its actual formal gauge.

The target map is FINITE DATA on any given finite set of localized
source generators. For g_i=t^-k_i[A_i(t,X)],

    pp_P(g_i)=sum_(j=1..k_i) [t^(k_i-j)]Phi[A_i] t^-j.   (5)

These are finitely many ORIGINAL Wick coefficients. In the point-module
basis generated by t^-1, replace t^-j by
(-1)^(j-1)D_t^(j-1)(t^-1)/(j-1)!.
Thus a correct full time-localized presentation, retaining its original
amplitude maps, can be restored by the kernel of ONE explicit finite
map to delta_0. This is not the full actual-period image kernel.

Time localization is therefore permissible WITH this restoration.
It is still invalid to discard amplitudes, replace the source by a
larger annihilating equation, or invert every nonzero time polynomial:
the last operation loses source point modules away from0 that (3)
cannot recreate. No general restored-source CAS backend is executed here.

## 4. Why other punctual source kernels must survive

For the genuine unary Gaussian phase P=X^2/2,

    nabla=partial_X+(t-1)X,
    (t-1)[X]=0,       D_t[X]=[X^3/2].

The class [X] is NONZERO. If X=nabla g for polynomial g(t,X), its
highest X-degree forces g to be independent of X, after which its
coefficient would be1/(t-1), not polynomial. The nonzero map from the
simple point module delta_1 sending its generator to [X] is injective.
All these odd amplitude periods vanish, so delta_1 lies in ker Phi.
Multiplication by t is invertible on delta_1, as required by the theorem.

This example explains the exact correction: origin t-torsion is absent,
but arbitrary point-supported source kernels are not. It preserves the
need for the retained actual-kernel algorithm on general sources and
for full polynomial-Weyl presentations before safety decisions.

## 5. The origin boundary is a specific hyperplane, not arbitrary data

Work over a characteristic-zero field, so H and M are holonomic. Let
L=M[t^-1], and let J be its intermediate extension across0. The retained
[boundary-subspace theorem](RATIONAL_CORES_AND_BOUNDARY_SUBSPACES.md)
gives J subset M subset L and

    L/J=delta_0 tensor E_0.

The principal-part map annihilates J, since J has no point quotient.
It therefore induces a NONZERO scalar functional ell:E_0->K. Equation(4)
identifies the selected boundary subspace EXACTLY as

    M/J=delta_0 tensor ker ell.                            (6)

It is a codimension-one hyperplane, in every dimension and phase degree.
The same assertion holds for the exact source H and L=H[t^-1]. Over
general coefficient rings the exact sequences(3)-(4) remain valid, but
no field Grassmannian/semisimplicity assertion is imported.

If boundary vectors are lifted as g_i=t^-k_i[A_i], ell(g_i) is just the
single coefficient [t^(k_i-1)]Phi[A_i], since a vector killed by t in
delta_0 is a multiple of t^-1. Thus FINITE actual jets determine this
hyperplane once the localized module and its embedded boundary basis
are known. No actual generic-rank assumption or convergence premise enters.

This is a universal structural restriction on Gaussian origin lattices,
not a classification of their generic connections. It does not force
M=J: dim E_0 can exceed1. Nor does it make the unit essential in every
pure actual image. The
[global applications](GLOBAL_EULER_AND_GAUSSIAN_PULLBACK_LATTICES.md)
compute both possibilities and a genuine counterexample to essentiality.

## 6. Authority and remaining scope

GMCFORMALCOMP20261115 is amended only where necessary: its origin
t-torsion warning is withdrawn and its canonical current source is this
corrected statement. Its original sealed proof/checker artifacts remain
unchanged as historical evidence. Completion, unique formal Hom and all
later field/all-ring safety results are retained; no assurance is promoted.

This closes time-zero t-saturation and the canonical one-point restoration
and boundary-hyperplane questions. It does not compute the full actual
kernel, classify every generic Gaussian connection, or constrain boundary
points away from0. Both old capped jobs remain terminal UNKNOWN. The
two-O(2)-block closure and all Oaku/formal/ordinary-ideal packets remain
at their credited scopes. Sabbah/Maisonobe retain convergent rational-slope
credit; arbitrary nonconvergent formal/all-ring extension priority remains
unassessed. The broad goal remains active.
