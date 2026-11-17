# Strict leading nullcones, full thin closure and ordinary isotypic safety

16 November 2026. Internal written proofs; independent specialist review
and formal verification remain pending. No novelty is claimed. The
[checker](verify_strict_leading_thin_closure.py),
[receipt](strict_leading_thin_verification.json) and
[source audit](strict_leading_thin_source_audit.json) distinguish selected
exact controls from the universal arguments below.

The previous [half-fibre and single-isotypic theorem](SEMISTABLE_LOG_PERIODS_AND_SINGLE_ISOTYPIC_SAFETY.md)
left an equality case for genuinely mixed phases. That equality case is
now excluded. On BOTH exceptional quotients, every nonzero eventually
pure FULL phase has a STRICTLY unstable leading tuple: a common left
root of multiplicity greater than half the matrix degree.

Consequently the eventual pure radical on the ENTIRE thin-image sector
is zero. Every polynomial in ANY nonconstant actual regular real seed
is excluded from nonzero eventual purity, with no isotypic restriction.
The old full thin existence gates, including degree36/60 and bad-marker
degreeD-12/D-20, are subsumed; they are retained conditional results,
not remaining existence questions. The separate symmetrized Hopf-height
law is not a pure-phase law and remains separate.

The full MZ questions remain **OPEN** for mixed-isotypic phases whose
compact images have interior. A strict leading nullcone need not make
the FULL phase unstable; opposite lower weights are retained.

## 1. The central quotient changes the boundary pole from three to two

Both binary groups contain -I. Work with the lifted invariant functions
on PGL2=SL2/{+I,-I}, before taking the remaining finite right quotient.
The compact core is SO3, with a nonzero normalized Haar volume period.
The finite central cover is etale, so it preserves the moment equations
up to a fixed nonzero normalization.

Compactify PGL2 as X=P3, with boundary B=(Delta=0), Delta=xw-yz.
A regular function of minimal matrix degree D=2k has the unique form

    F=H/Delta^k,   H homogeneous of degree2k,
    H|B!=0.                                                   (1)

Here H is the even determinant homogenization already retained in the
[radial bridge](GAUSSIAN_RADIAL_TRANSFER_AND_EXCEPTIONAL_UNIFORMITY.md).
The boundary is P1_left times P1_right, and H|B has bidegree(2k,2k).
The invariant top form has divisor -2B: K_P3+2B is trivial. More
concretely, on x=1 put y=y,z=z,w=yz+v, so Delta=v. Up to a nonzero
constant the ORIGINAL invariant volume is

    Omega=dy wedge dz wedge dv / v^2.                         (2)

For a divisor E over this chart put a=ord_E(v). Its volume order is

    ord_E(Omega)=A_X(E)-1-2a.                                 (3)

A boundary divisor above a FINITE value of F must have
ord_E(H)>=k*a. No arbitrary Jacobian weight replaces (2).

## 2. A local critical-pole lemma

Let K be a characteristic-zero field, and work at a smooth codimension2
centre with completed transverse ring K[[v,z]]. Suppose

    H(0,z)=z^k times a unit, k>=1,
    Omega=unit*dv wedge dz/v^2.

For a divisorial valuation E centered there, assume

    ord_E(H)>=k*ord_E(v),
    A(E)>=ord_E(v)+(1/k)*ord_E(H).                            (4)

The second inequality is the log-canonical inequality for
(v=0)+(1/k)div(H). If mult_(v,z)(H)<k, then

    A(E)>2*ord_E(v);                                         (5)

in particular the original volume has no pole on E.

Proof. Inequality(4) first gives A(E)>=2a. If a pole persists,
integrality in(3) forces A(E)=2a, ord_E(H)=ka and volume order -1.
Weierstrass preparation replaces H by a unit times a monic polynomial
of degree k in z with reduction z^k. Over an algebraic closure of K,
take one finite ramification v=t^e splitting that polynomial:

    H(t^e,z)=unit*product_(i=1)^k (z-phi_i(t)),
    phi_i(0)=0,   phi_i in Kbar[[t]].                         (6)

Only formal root factorization is needed, not analytic continuation.
For completeness, finite ramification follows by the usual degree
induction. Depress a monic polynomial to remove its penultimate
coefficient. If it is not z^k, set q=min_i ord(a_i)/i, ramify so e*q
is integral, substitute z=t^(eq)w and divide by t^(eqk). The reduction
is a nonmonomial monic polynomial with zero penultimate coefficient.
In characteristic zero it has at least two distinct roots. Hensel
factorization splits it into smaller-degree monic factors. Apply the
induction to each. The polynomial z^k is the immediate repeated-root
case. This is the standard formal Weierstrass/Hensel argument.

Lift E to a divisorial valuation E' on this finite cover; write
b=ord_E'(t)>0. A simple logarithmic volume pole pulls back to a simple
logarithmic pole: in characteristic zero the local order formula is
rho*(-1)+(rho-1)=-1, where rho is the ramification index at E'.
But the actual pulled-back form in these transverse coordinates is

    e*unit*dt wedge dz/t^(e+1).

Thus A(E')=(e+1)b. For each i the pair (t=0)+(z-phi_i(t)=0)
is a simple-normal-crossing pair in regular coordinates. Its elementary
log-discrepancy inequality gives

    b_i:=ord_E'(z-phi_i(t)) <= A(E')-b=e*b.                  (7)

Valuations of(6), and the equality ord_E(H)=ka, give
sum_i b_i=k*e*b. Therefore EVERY b_i=e*b. It follows that
ord_E'(phi_i-phi_j)>=e*b, so all roots agree modulo t^e.
Their mean is the negative penultimate Weierstrass coefficient divided
by k. It belongs to K[[v]] and has zero constant term, hence has order
at least e in t. Therefore every phi_i has order at least e in t.
Every i-th Weierstrass coefficient consequently has v-order>=i.
The prepared polynomial has (v,z)-multiplicity>=k, and multiplication
by a unit preserves this multiplicity. This contradicts the hypothesis.

The argument works after passing to the algebraic closure of the generic
centre field. No valuation of that field is introduced. It does not
assert that differential poles on arbitrary affine models are harmless.

## 3. Excluding the half-degree leading case for FULL phases

Let F be a nonzero eventually pure phase in the full original exceptional
algebra. Constants are excluded immediately. Write D=2k>0 and(1).
The previous half-fibre theorem already gives k_max>=k. Suppose k_max=k.
Remove the common left factor of H|B. On every left fibre the residual
section restricts to a nonzero invariant right binary form of degree2k.
Its roots have multiplicity at most2k/r, where r=6 or12. In particular
their coefficient at scale1/k is at most2/r<1.

Exactly as in the earlier adjunction proof, the boundary pair

    (B,(1/k)div(H|B))                                       (8)

is log canonical. Its only log-canonical centres are the one or two
left fibres L of multiplicity k. Near each such L it is plt, with L
the coefficient1 divisor. Everywhere else it is klt. This description
includes scheme multiplicities and all right-invariant multiplicities.

Kollar7.5.2 (the same [primary source](https://arxiv.org/pdf/alg-geom/9601026v2),
printed43, with its proof) gives that

    (X,B+(1/k)div(H)) is log canonical near B.                (9)

All its hypotheses hold: X and B are smooth, div(H) is Cartier, and
H does not contain B. Thus every finite boundary divisor has A(E)>=2a,
so (3) allows at worst a simple pole. Such a pole requires equality
in both inequalities, hence is a log-canonical place of(9).

There are no such places centered at a closed point of L. To see this
without assuming an lc-centre correspondence theorem, choose a local
smooth divisor R_B through that point, transverse to L. A sufficiently
small positive rational epsilon keeps the restricted pair on L klt:
the old restricted coefficients are strictly below1. Adjunction gives
that (B,(1/k)div(H|B)+epsilon R_B) is still lc locally. Extend R_B
to a Cartier divisor R on X and apply7.5.2 again. A valuation centered
at the point has ord_E(R)>0, so its ORIGINAL log discrepancy is at
least epsilon*ord_E(R)>0. It cannot be the required lc place. Off the
fibres L, plt from7.5.1 excludes exceptional lc places as before.
The divisor B itself maps to infinity and is not a finite boundary.

It remains to consider divisors whose centre is the generic point of L.
Choose a left torus fixing that root and use the chart in Section1,
with L=(v=z=0). Under diag(t,t^-1), after projective normalization x=1,

    (y,z,v)->(y,t^-2*z,t^-2*v).

If F_ell is its weight-ell component, the numerator Delta^k F_ell
has (v,z)-homogeneous order k-ell/2. Consequently

    mult_L(H)=k-(maximum occurring left weight)/2.           (10)

If no positive weight occurs, the actual left cocharacter has a
nonzero zero-weight limit F0. The leading tuple has multiplicity
EXACTLY k, so its zero-weight coefficient is nonzero; lower isotypic
components cannot cancel this top component. All eventual moment
equations pass to this limit as invariant polynomials in the finite
filtered module. F0 is a nonzero original torus-fixed eventually pure
phase, contrary to the complete retained torus-invariant closure.

Thus a positive weight occurs, and(10) gives mult_L(H)<k for EVERY
critical fibre L. Apply the local critical-pole lemma to each generic
centre. It removes the remaining finite volume poles.

Resolve the projective graph as an isomorphism on PGL2 and restrict
over A1. Every finite boundary divisor now has nonnegative volume
order, including vertical exceptional divisors. The original volume
extends holomorphically by codimension1 regularity on the smooth
resolution. The map is PROJECTIVE. The retained projective log-period
theorem applies to the compact SO3 cycle with nonzero volume period
and contradicts eventual purity. This proves the strict condition

    EVERY nonzero eventually pure FULL phase has k_max>D/2. (11)

This proof genuinely treats lower isotypic components. It uses them
in(10); it never discards them in a leading degeneration. It also does
not imply that the full phase is unstable: different lower components
can carry opposite weights. Apply(11) separately to every algebraic
automorphic image, with its own degree, using intrinsic Haar invariance.

There is a further stabilizer consequence. The strict common root is
unique: two roots each of multiplicity> D/2 cannot occur in a nonzero
degree-D tuple. The full left stabilizer preserves the top component,
so it fixes this projective root and lies in its Borel subgroup. If it
is finite, averaging a projection onto its invariant line gives an
invariant complementary line. It is conjugate into the diagonal torus
and is cyclic. The retained finite-stabilizer theorem therefore improves
for full failures to a FINITE CYCLIC stabilizer, for every automorphism-
conjugated left action. This necessary condition does not certify failure.

## 4. Every full thin-image phase is excluded

Let h be ANY nonconstant regular real-valued seed on the compact core,
with minimum matrix degree H. Compact dagger preserves the degree
filtration and acts on the rank-one boundary by the fixed-point-free
antipodal involution on each P1 factor. Its common left roots therefore
come in distinct equal-multiplicity pairs. Their total left degree is H,
so each has multiplicity at most H/2.

For a polynomial p of degree d>0, the associated-graded domain gives
deg(p(h))=dH and leading section lc(p)*(h_top)^d. Every common left root
of this leading section has multiplicity at most dH/2. It fails(11).
Thus NO nonzero polynomial in h is eventually pure. Nonzero constants
are excluded separately. This includes mixed seeds and all polynomial
reparametrization degrees.

The retained curve-normalization theorem expresses EVERY thin-image
regular phase as p(h) for an ACTUAL regular real seed h. Therefore

    eventual pure radical on the full thin sector = {0}.     (12)

Here thin means empty interior of the compact image, equivalently
real differential rank at most1. Both full quotient models and all
matrix degrees are covered. Every nonzero pure phase, and hence any
full MZ failure, now has compact image with interior/rank2 somewhere.

The former thin degree bounds36/60 and bad-marker boundsD-12/D-20,
Jordan-loop/inner-branch existence questions for actual Haar laws, and
the normalized real mixed-seed existence question are SUBSUMED.
Keep their proofs and the abstract positive loop law: that law was
never an actual Haar-seed realization. Finite seed-power bad-marker
theorems remain valid general conditional statements. The separate
even/symmetrized Hopf-height law is not asserted to be a pure-phase
identity and is NOT settled here.

Fixing any such real h also closes the restricted coefficient-ring
problem. In C[c_0,...,c_d], every tail ideal

    I_M=(Haar(p(h)^m):m>=M), M>=1,

has zero locus only the coefficient origin by(12). Its radical is the
maximal coefficient ideal. Over ANY commutative C-algebra R, an
eventually pure p(h) therefore has all coefficients nilpotent. Their
finitely generated ideal is nilpotent, so p(h) itself is nilpotent.
The converse is immediate for eventual moments. Hence the restricted
eventual radical is exactly the nilpotent polynomials in R[h], and
EVERY original marker Q has Q*p(h)^m=0 eventually. Nilpotents need not
be pure in their first moments.

For each fixed d, Noetherianity also gives integers N,L with
(c_0,...,c_d)^L contained in the first-N ordinary pure ideal. That
finite prefix forces p(h)^L=0 over every R. This is an existence
theorem; no general numerical N,L are computed. It uses universal
coefficient ideals over C, so no impossible embedding of an arbitrary
extension of C back into C is assumed.

## 5. Closing the ordinary nonreduced single-isotypic question

The previous all-degree single-isotypic field theorem has a separate
ordinary consequence, using the already retained
[derivative/symbolic-power method](NILPOTENT_LIFTING_AND_FORMAL_SAFE_COMPLETION.md)
and [positive-family stabilization](POSITIVE_FIELD_THEOREMS_AND_ORDINARY_CUTOFF_STABILIZATION.md).
Its proof here uses the explicit field cutoff2m>deg Q. It does NOT
depend on general Gaussian holonomic uniformity.

Fix W_D=V_D tensor(V_Dstar)^Gamma, D=2k>0, with right multiplicity r>0,
and omit the constant coordinate, which the first moment kills over
every coefficient ring. Put n=(D+1)r and work in the REGULAR ambient
coefficient ring S=C[c_1,...,c_n]. Let P be the universal phase and

    I=(Haar(P^m):m>=1),  I_N=(Haar(P^m):1<=m<=N),
    J=the radical nullcone ideal=sqrt(I).

The nullcone is the proper image of the rank-rk vector bundle over
P1 parametrizing tuples divisible by a common linear factor^(k+1).
It is irreducible, and its common root is unique for every nonzero
tuple. Hence its dimension is rk+1 and its codimension is

    eta=r(k+1)-1.                                           (13)

Thus J is prime; eta>=1. Choose N,nu>=1 with J^nu contained in I_N.
They exist because sqrt(I)=J and S is Noetherian. They are in principle
obtainable by enumerating pairs and exact ideal-membership certificates
over the fixed algebraic field of definition; no such campaign is run.

For a marker Q of degree e and a coefficient derivative of order s,

    partial_c^alpha Haar(QP^m)
      =(m)_s Haar(Q*B_alpha*P^(m-s)),
    deg(B_alpha)<=Ds.                                       (14)

At every nullcone point this vanishes when

    2m>e+(D+2)s.                                            (15)

Derivatives of order greater than m vanish identically. Consequently
all derivatives of orders<K vanish on the nullcone if
2m>e+(D+2)(K-1). Radical Zariski--Nagata identifies this with symbolic
power J^(K). Ein--Lazarsfeld--Smith, TheoremA/2.2, in the SMOOTH
coefficient space gives J^(eta*nu*L) contained in J^(nu*L).
The exact radical differential statement is Sidman--Sullivant2.9.
These are the same primary imports already audited in the retained
nilpotent proof, applied BEFORE the possibly nonreduced pure quotient.
We obtain the ORDINARY containment

    2m>e+(D+2)(eta*nu*L-1)
       => Haar(QP^m) in J^(nu*L) subset I_N^L.               (16)

In the symbolic term above J^(eta*nu*L) denotes a symbolic power;
J^(nu*L) and I_N^L in(16) are ordinary powers. No equality between
symbolic and ordinary powers of the same exponent is asserted.

Set E_D=(D+2)(eta*nu-1). After any coefficient map to ANY commutative
C-algebra, a pure phase in C+W_D therefore satisfies

    2m>deg(Q)+E_D => Haar(QP^m)=0                            (17)

for every ORIGINAL marker. The target ring may be non-Noetherian and
have arbitrarily deep nilpotents. The additive intercept depends on the
fixed universal component, not on that target ring or its thickness.
It is not asserted zero, and the field cutoff is not silently copied
unchanged to nonreduced rings.

Formula(16) also gives a finite ordinary pure prefix explicitly from
the certificate pair: K=max(N,floor(E_D/2)) generates I, since all
remaining constant-marker moments already lie in I_N. For an integer
T>=ceil((E_D+1)/2), the full ordinary ideal of marked moments with
m>T*deg(Q), including constant markers, equals I. Thus the restricted
ordinary cutoff ideals stabilize, retaining their nilpotent structure.
No numerical sufficient V12 prefix or least cutoff is computed. The
ordinary V12 invariant-quartic deficit is preserved; equality with the
invariant ideal or reducedness of I is NOT claimed.

## 6. Remaining all-degree objective

Both FULL MZ cases remain OPEN on mixed-isotypic phases with compact
image interior, strictly unstable leading tuple, and finite cyclic
stabilizers/automorphic-torus exclusions. All original markers and matrix
filtrations remain in force. Leading instability alone does not eliminate
opposite lower weights or supply a marked-safety proof.

No degree/support/density/prefix search is expanded. The main remaining
geometric question is how the global periods and opposite lower weights
interact in this strict-leading, image-interior architecture. The whole
proof, particularly the projective period input and the critical local
root/valuation argument, still requires independent specialist review.
