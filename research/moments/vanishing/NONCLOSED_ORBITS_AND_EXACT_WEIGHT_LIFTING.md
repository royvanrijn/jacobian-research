# Nonclosed pure orbits and exact lifting from torus-fixed periods

16 November 2026. Internal written proofs, with selected exact controls.
Independent replay, external review and formal verification are not claimed.
External review is deferred under the user's closure-first priority: the
main finish line remains a COMPLETE exceptional Mathieu--Zhao theorem.
See the [checker](verify_orbit_weight_lifting.py),
[receipt](orbit_weight_lifting_verification.json) and
[source audit](orbit_weight_lifting_source_audit.json).

This makes progress on genuinely mixed phases. On BOTH full exceptional
quotients, a nonzero eventually pure phase has a nonclosed left orbit
EXACTLY when it lies in the left nullcone. It is then safe for ALL original
markers, with sufficient cutoff2m>deg_matrix(Q), in every degree.

The remaining non-null pure phases have CLOSED left orbits and finite
cyclic stabilizers. Such an orbit is not a failure certificate: the known
safe Haar-shear witnesses already have precisely this orbit type. A full
failure must stay in this closed-orbit sector after EVERY algebraic
Haar-preserving automorphism, and must have at least THREE weights for
every maximal left torus. The two-weight closure below includes arbitrary
mixed isotypic degrees, all right multiplicities and arbitrary coefficient
rings. Its general form is an all-rank theorem for affinely independent
torus supports with an exact reduction of ALL original marked sequences.

Both FULL exceptional MZ decisions remain **OPEN**. The precise unresolved
step is controlling cancellations between MULTIPLE balanced monomials in
the remaining automorphically closed, mixed, image-interior sector.
No external reviewer is the next research step, and no growing coefficient
or moment-prefix search is introduced.

## 1. A curve and a two-by-two Smith factorization

Put A=C[SL2/Gamma], Gamma=2O or2I, and let F be nonzero and eventually
pure. A finite matrix-degree filtered space V contains F and is a finite
rational left-SL2 module. All eventual equations Haar(F^m)=0 are invariant
polynomial equations on V. Assume its affine left orbit is not closed.

Choose a boundary point and an algebraic curve in the orbit closure,
with punctured curve in the orbit. Normalizing the curve and completing
at that point gives a regular vector arc v(t) in V[[t]], with v(0)
outside the orbit. Over an algebraic closure of C((t)), its generic point
is g(t)F for some g in SL2. The finitely many matrix coordinates of g
are algebraic, so one finite extension suffices. Replacing the parameter
by a finite ramification puts g in SL2(C((t))) while keeping the vector
arc regular and its boundary specialization unchanged.

Here curve selection only uses that the orbit is locally closed and
dense open in its closure; a general curve through a boundary point
avoids the boundary generically. A lift of its generic geometric point
requires only a finite extension. No lifting over the entire curve or
assumption that the orbit map is a trivial torsor is made.

Smith normal form over the discrete valuation ring C[[t]] gives

    g(t)=U(t) diag(t^a,t^-a) V(t),
    U,V in SL2(C[[t]]), a>=0.                               (1)

The determinants of the two invertible factors are units. A diagonal
unit commuting with the middle factor adjusts both determinants to1;
a Weyl matrix interchanges the two exponents if necessary. If a=0,
g is regular and invertible at0, contrary to the boundary specialization.
Thus a>0.

Since U^-1 is regular, diag(t^a,t^-a)V(t)F is regular. Decompose V(t)F
into torus weights. A negative-weight component of V(0)F would have
a nonzero constant coefficient multiplied by a negative power of t.
Higher-order terms in that SAME weight cannot cancel its first pole,
and different weights are independent. Therefore V(0)F has no negative
weight. Its actual torus limit is its zero-weight part F0.

The phase V(0)F is eventually pure, and its torus limit is eventually
pure by the invariant polynomial equations. But F0 belongs to the
ENTIRE original torus-fixed algebra, whose eventual pure radical is0
by the retained [torus closure](VALUATIVE_BOUNDARY_VOLUMES_AND_TORUS_INVARIANT_CLOSURE.md).
Hence F0=0. All weights of V(0)F are strictly positive, so an actual
left cocharacter contracts F to0. This proves that F is unstable.

Conversely, a nonzero unstable vector has0 in its orbit closure and
cannot have a closed orbit. Thus, ON the eventual pure locus,

    nonclosed left orbit iff left-nullcone membership.      (2)

Every original marker Q has left weights of absolute value at most its
matrix degree. All function weights are even because -I acts trivially.
After the degree-preserving linear left conjugation above, the minimum
phase weight is at least2. Therefore

    2m>deg_matrix(Q) => Haar(QF^m)=0.                       (3)

This is a full mixed-component result, not a statement about the leading
component alone. The existence of a nonclosed orbit is essential; a
closed orbit does not supply the arc with a boundary point used in(1).
The proof gives its own rank-one orbit-closure argument and does not
silently import a more general orbit theorem without its hypotheses.

## 2. What a non-null pure orbit must look like

Suppose F is eventually pure and outside the left nullcone. Equation(2)
says that its left orbit is closed. Its stabilizer is finite as well.
Indeed a nonzero element of its stabilizer Lie algebra in sl2 is either
semisimple or nilpotent. A semisimple element makes F torus-fixed and
therefore0, impossible. A nilpotent element makes F fixed by a root
unipotent subgroup. Its invariant algebra is C[z,w]^Gamma after left
conjugation. Its nonconstant weights all have one strict sign, and its
constant term must vanish by eventual purity. Such an F is unstable,
again impossible. A characteristic-zero zero-dimensional algebraic
stabilizer is finite.

The [strict leading theorem](STRICT_LEADING_NULLCONES_AND_FULL_THIN_CLOSURE.md)
gives a UNIQUE common leading left root of multiplicity>D/2. The full
stabilizer fixes this root, so it lies in its Borel subgroup. A finite
subgroup splits the fixed line by averaging and is conjugate into the
diagonal torus; it is cyclic. Consequently every non-null eventually
pure phase has a closed orbit with finite cyclic stabilizer, whether
it is safe or failing.

Intrinsic Haar invariance transports purity and all-marker safety under
every algebraic automorphism Phi. If ANY Phi(F) has a nonclosed left
orbit, equation(3), applied to Phi(Q), proves safety of F. The degree
after a nonlinear automorphism must be retained:

    2m>deg_matrix(Phi(Q)) is sufficient in that certificate. (4)

The original marker cutoff is not assumed unchanged by Phi.
Thus every FULL failure has closed left orbit with finite cyclic
stabilizer for EVERY Phi(F), not just for one selected coordinate system.

For a known pure F, the original orbit type also has a finite exact
test. Its strict leading root is unique. In the full Peter--Weyl
decomposition F=sum F_d, contraction to0 is possible exactly when this
same root has multiplicity>d/2 in EVERY nonzero tuple F_d, with zero
constant. This is the common-root one-parameter test for the entire
finite module. If it fails, (2) makes the pure orbit closed. It is not
a purity decision for an arbitrary input.

The known [Haar-shear family](HAAR_SHEARS_AND_GENERAL_NULLCONE_CLASSIFICATION.md)
is a necessary countercontrol. For n=8 or12 and t!=0,

    F_t=f_n(x+t*q*z,y+t*q*w), q=f_n(z,w),

is pure, safe and outside the left nullcone by its retained nonzero
invariant. Hence its ORIGINAL orbit is closed with finite cyclic
stabilizer by the theorem above. Undoing the shear sends it to the
one-sided highest vector f_n(x,y), whose orbit is nonclosed. Closedness,
finite cyclic stabilizer and two-sided support alone do NOT imply failure.

## 3. Exact all-marker lifting for an affinely independent weight support

Here is a general theorem independent of the period and orbit arguments.
Let R be a commutative Q-algebra, let A be a commutative R-algebra with a
rational torus grading A=direct_sum_(w in Z^r) A_w, and let h:A->R be
invariant and R-linear. Write h0 for its restriction to A_0. Suppose
ker(h0) is a Mathieu--Zhao subspace in A_0: every pure B in A_0 has
h0(CB^ell)=0 eventually for every fixed C in A_0.

Let

    P=P_1+...+P_s,  P_i in A_(w_i),

where the DISTINCT weights w_i are AFFINELY INDEPENDENT. Then every
pure P is safe for ALL original markers in A, in arbitrary degrees.
The same assertion holds with eventual purity in place of purity.

If0 lies outside conv{w_i}, a strict rational separator makes every
power one-sided and every fixed marked sequence eventually zero.
Suppose instead that0 lies in that convex hull. The augmented weight
matrix with columns(w_i,1) has rank s, while the weight matrix W has
rank s-1. Its integer kernel is generated by a unique primitive vector

    a=(a_1,...,a_s)>=0,  gcd(a_i)=1,
    M=sum_i a_i>0,  B=product_i P_i^(a_i) in A_0.             (5)

Some a_i can be zero when0 lies on a proper face. The complete
nonnegative relation monoid is N*a. Therefore the exact identities are

    h(P^m)=0 if M does not divide m,
    h(P^(M*ell))=multinomial(M*ell;a_i*ell)*h0(B^ell).        (6)

All multinomial factors are nonzero rational units in R. Thus purity,
and separately eventual purity, of P is equivalent to the corresponding
condition on B. This is equality of the ORDINARY pure coefficient ideals
under any finite parameterization, not merely equality of radicals.

Now let Q be a fixed homogeneous marker of weight tau. Consider

    W*n=-tau,  n in N^s.                                   (7)

If there is no solution, its entire marked sequence is zero. Otherwise
the integer solutions differ by Z*a. The nonnegative solutions form
ONE ray n=n0+ell*a, ell>=0, for its minimal nonnegative representative
n0. This remains true when some a_i=0; the corresponding coordinates
are fixed by(7). Put m0=sum n0_i and C=Q*product P_i^(n0_i) in A_0.
The entire original marked sequence is then EXACTLY

    h(QP^m)=0 outside m=m0+M*ell,
    h(QP^(m0+M*ell))
      =multinomial(m0+M*ell;n0_i+ell*a_i)*h0(CB^ell).        (8)

The assumed MZ property on A_0 proves eventual vanishing in(8).
An arbitrary Q has only finitely many homogeneous components, so one
finite maximum of their cutoffs proves safety for that original Q.
No marker is required to lie in a smaller phase-generated algebra.

For eventual rather than full purity, choose an integer N for which
B^N is pure and apply the MZ property to the finitely many markers
CB^j, 0<=j<N. This covers every residue class of ell and proves the
eventual version without changing the hypothesis on A_0.

This proof explicitly handles zero weights, faces, primitive divisibility,
all original markers, and ordinary coefficient equations. It does not
extend from a single ray to a relation monoid of rank at least2.

## 4. Complete two-weight closure on both FULL exceptional quotients

Apply Section3 to the original left torus and use the retained ALL-ring
torus theorem: on (A_Gamma tensor R)_0 the eventual pure radical is
exactly the nilpotents for EVERY commutative C-algebra R. A set of at
most two distinct weights is affinely independent. This proves full
MZ for EVERY such phase, with arbitrary mixed isotypic degrees and
all right multiplicities, over every R.

There is a more explicit complete eventual-purity criterion:

* If the support is strictly on one side of0, every phase is pure and
  safe. The weight bound is m*min|w_i|>deg_matrix(Q).
* For support{0,a}, a!=0, write P=P0+Pa. Eventual purity is equivalent
  to NILPOTENCE of P0 in the original algebra. If P0^nu=0, all marked
  moments vanish whenever m>deg_matrix(Q)/|a|+nu-1.
* For support{a,-c}, a,c>0, put g=gcd(a,c), p=c/g,q=a/g,M=p+q and

      B=P_a^p*P_(-c)^q.

  Eventual purity is equivalent to NILPOTENCE of this actual T-fixed
  product B. If B^nu=0, a sufficient original all-marker bound is

      m>=M*nu+ceil(deg_matrix(Q)/min(a,c)).                  (9)

  Indeed a contributing term has counts
  i=(c*m-tau)/(a+c), j=(a*m+tau)/(a+c), |tau|<=deg Q;
  under(9), i>=p*nu and j>=q*nu, so the term itself contains B^nu.
* For support{0}, eventual purity is exactly nilpotence, with the
  immediate cutoff m>=nu.

The sufficiency statements for nilpotents concern EVENTUAL moments.
They do not identify nilpotence with vanishing of every first pure
moment. For opposite weights, the literal pure ideal is precisely the
pullback of the ordinary pure ideal for B through(6). Over fields the
original algebra is a domain, so B=0 forces one charged component to
vanish. For support{0,a}, P0=0 over fields. Thus a nonzero pure FIELD
phase with at most two weights always has one STRICT sign. Its uniform
all-marker cutoff is2m>deg Q, since all weights are even.

Over rings the two signs may live on different components; do not demand
one global sign after base change. For example in R=C[e,f]/(ef),
P=e*f_n(x,y)+f*f_n(z,w) is pure and safe although both charge blocks
are nonzero. Their balanced product is0. Ordinary nilpotent information
is retained throughout and no arbitrary coefficient ideal is replaced
by its radical.

In every fixed finite coefficient family in these two weights, purity
also bounds the nilpotence exponent in(9) uniformly over all target
rings. In its finite ordinary pure quotient, the coefficients of B
(or P0) vanish at every complex point, hence generate a finitely
generated nilpotent ideal. Noetherianity makes one exponent work in
that universal quotient. This is a finite-family existence conclusion;
no numerical universal prefix or exponent is computed here.

Consequently every non-null pure field phase has at least THREE
weights for every maximal left torus, and any FULL failure has at
least three weights after EVERY algebraic automorphism as well.
This closes a full mixed-degree class, not only a sparse subset of
one Peter--Weyl component.

## 5. The exact global implication that remains

The common strict leading root determines the only possible contracting
line, but opposite lower weights can still prevent contraction. Section1
now resolves the orbit-closure alternative completely on the pure locus.
The unresolved candidates are pure mixed phases with image interior,
strict leading nullcone and CLOSED finite-cyclic left orbits after every
algebraic automorphism, with at least three weights for every left torus.
All original fixed markers remain in the required conclusion.

One SUFFICIENT global route to a complete positive theorem would be:

    For every nonzero pure F, some algebraic automorphism Phi makes
    the left orbit of Phi(F) nonclosed.                     (10)

Then(2)--(4) prove the full MZ assertion. Equivalently, this route asks
whether every pure phase lies in the automorphism saturation of the
left nullcone. Its sufficiency is proved; its universal hypothesis is
OPEN and is not asserted necessary for MZ. A minimal-degree representative
in an automorphism orbit exists by well-ordering, but no degree-reducing
replica for every remaining stable pure representative has been proved.
The mere existence of the unique strict leading root does not supply it.

There is an equally precise obstruction to extending Section3 blindly.
For three or more weights in one angular dimension, different balanced
monomials can contribute to the same pure moment. Torus-fixed MZ alone
does not prevent their scalar period sums from cancelling. The standard
interval example, on C[U,T,T^-1] with h=integral_0^1 CT_T, has

    F=(1-T^-1)((1-U)+U*T)
      =(U-1)T^-1+(1-2U)+U*T,
    h(F^m)=0,
    h(T^-1 F^m)=(-1)^(m+1)/(m+1), m>=1.                   (11)

Indeed integration of ((1-U)+U*T)^m gives
(1+T+...+T^m)/(m+1), by the elementary beta integral for each
binomial coefficient. Multiplying by (1-T^-1)^m makes the pure
constant term the complete alternating binomial sum, zero for m>=1.
The marker T^-1 deletes exactly its last term and gives(11).

Its zero-weight algebra C[U] has zero eventual pure radical by the
retained fixed-curve period theorem, yet its three-weight extension
fails MZ. The affine dependence is essential. This Laurent algebra
does NOT embed into the original exceptional algebra preserving its
units and periods: T is a nonconstant unit, whereas the exceptional
algebras have only constant units. Equation(11) is a countercontrol to
an invalid inference, not an exceptional counterexample.

The fresh [six-variable cotangent classification](COMPLETE_COTANGENT_PHASE_CLASSIFICATION.md)
uses a proven prime-separating marked period on a specific two-generator
algebra. Its prime certificate cannot be assumed for the actual
four-variable Gaussian lift or exceptional Haar law. The safe two-sided
Haar-shear witnesses above must pass any proposed global test. A valid
arithmetic route must construct and prove its actual marker/period
identities and handle multiple balanced monomials; the six-variable
example supplies neither automatically.

The next main task is this remaining global implication or an actual
full failure, preferably closing one whole exceptional quotient if the
two geometries diverge. Qualitative thin and single-component closures
are retained evidence, not an expanding research queue. External review
is deferred; the mathematical assurance flags remain unchanged and
internal scrutiny of all dependencies continues.
