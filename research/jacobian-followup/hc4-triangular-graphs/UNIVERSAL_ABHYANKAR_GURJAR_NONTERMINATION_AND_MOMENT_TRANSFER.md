# Every plane Keller extension admits all-order inversion-term nontermination

HC4MYMM0. 16 November 2026. ALL actual complex polynomial plane Keller
extensions, every sheet degree, and EVERY nonzero polynomial amplitude.

The target framing of retained HC4MYMK0 preserves the literal embedded
target B, its full finite normalization S and its trace-dual Lambda. In
that framing every Abhyankar–Gurjar inversion term of every nonzero
polynomial amplitude is nonzero, with an exact degree in every order.
Nevertheless the formal sum for amplitude one is ONE; for a target
pullback it is a polynomial, and for an automorphism every polynomial
amplitude has a polynomial sum. This is a uniform obstruction to using
termwise moment vanishing as a universal plane inversion argument.
It is NOT an HC4 counterexample or a new failure of a moment conjecture.
The full actual-plane trace-pairing/marked-energy gate remains OPEN.

## 1. Retained framing and the classical inverse formula

Use MK0's construction for ANY actual plane Keller K. After affine source
coordinates and invertible target changes, for any integers M,N>=2 it
gives a plane Keller map F satisfying

    F(0)=0,       DF(0)=I,       det DF=1,
    deg F=D=M*N*q >=4,
    F_D=e2*phi,    phi homogeneous of degree D, phi(e2)!=0.   (1)

Here q=deg K and the target changes are two explicit tangent-identity
shears; their inverses and the first-jet issue are already proved in MK0.
The notation F is local to this note, not MK0's four-variable doubled map.
This target framing is NOT dynamical conjugacy of the original K. It
keeps B=C[F]=C[K], S and Lambda literally unchanged in the intrinsic
source ring. Source coordinates are only being renamed.

Put H=z-F(z), so ord(H)>=2 and H_D=-e2*phi. For a polynomial a(z) and
an integer m>=0 define the classical grouped inversion term

    T_m(a)=sum_(alpha1+alpha2=m) 1/alpha! *
                           partial^alpha(a*H^alpha).        (2)

The Abhyankar–Gurjar formula, reproduced coefficientwise in retained GS1
Section4, gives

    sum_(m>=0) T_m(a)=a composed with F^(-1),                (3)

as a formal series at zero. The usual Jacobian factor in that formula
is one here. Indeed ord(T_m(a))>=m+ord(a) for nonzero a, because each
H^alpha has order at least 2m and differentiation has total order m.
Thus every coefficient of (3) receives only finitely many contributions.
No analytic convergence or substitution at a nonzero numerical point
is asserted. The inverse on the right is the unique FORMAL inverse.

GS1 supplies a direct formal-residue proof of the exact formula needed
here. Retained GL1 supplies the Wick conventions, including moving
constant terms. These retained checkers are NOT replayed. The inversion
formula is classical; no priority claim is made for it.

## 2. Every amplitude, every order: exact nonzero top term

Let a!=0 have degree k and highest homogeneous part a_k. In (2), every
term has degree at most k+(D-1)m. Only alpha=(0,m) can contribute to
that highest degree: the degree-D part of H has ONLY its second component
nonzero. Consequently the complete top term is

    (T_m(a))_top = (-1)^m/m! * partial_2^m(a_k*phi^m).       (4)

This polynomial is NONZERO for EVERY m>=0 and EVERY nonzero a. To see
this without any genericity assumption on a, regard a_k and phi as
polynomials in z2 over C[z1]. By phi(e2)!=0, deg_(z2) phi=D. If
deg_(z2) a_k=j, then deg_(z2)(a_k*phi^m)=j+D*m. Its leading coefficient
is nonzero in the domain C[z1]. In characteristic zero an m-th derivative
of such a polynomial is nonzero, since j+D*m>=m. Hence

    T_m(a)!=0,
    deg T_m(a)=k+(D-1)*m,        ALL m>=0, ALL a!=0.         (5)

This is the unbounded written proof. Bounded symbolic controls only check
its transcription. It also covers amplitudes whose top part vanishes at
e2, such as a=z1; evaluation at e2 alone would NOT prove that case.

For amplitude one a useful scalar check of (4) is

    (T_m(1))_top(e2)
       =(-1)^m/m! * (D*m)!/((D-1)*m)! * phi(e2)^m !=0.      (6)

Equations (1)--(6) hold for every original extension, not just a supplied
inverse example. They do not impose any bound on its original generic
sheet degree or monodromy. The full original trace gap survives unchanged.

## 3. Polynomial sums with no zero individual terms

For a=1, (3) says sum T_m(1)=1. Every positive-order term is nonzero by
(5), and the infinite tail cancels coefficientwise. More generally, for
ANY nonzero target polynomial b and a=b composed with F, (3) says

    sum_(m>=0) T_m(b composed with F)=b.                    (7)

This is true even before deciding whether the original K is invertible.
Thus the target algebra already supplies amplitudes with polynomial
sums and nonzero terms in EVERY order in every framed extension.

Starting from the identity map gives actual polynomial automorphisms

    F=(z1+z2^M, z2+(z1+z2^M)^N),
    F^(-1)(u,v)=(u-(v-u^N)^M, v-u^N).                     (8)

For them (3) is a polynomial for EVERY polynomial amplitude, while (5)
still holds for every amplitude and order. In particular polynomiality
of the inverse does NOT imply eventual vanishing of the grouped terms
in this mixed-degree presentation. These are whole-plane inverse controls,
not constructed noninverses or unmarked finite-algebra realizations.

The quantities T_m are grouped ordinary derivatives in a classical inverse
formula. They are NOT the Nambu derivatives in the retained finite-energy
criterion, and NOT the Pascal composition differences of MK0. No failure
of Nambu termination or marked-energy integrality follows.

## 4. What a Gaussian moment producer actually sees

Let Z=(Z1,Z2), W=(W1,W2) have the algebraic Wick functional

    E(W^alpha*Z^beta)=delta_(alpha,beta)*alpha!.

This is the circular normalization from GS1/GL1, or a nonsingular formal
Gaussian in four spatial variables. For independent translation parameters
c=(c1,c2), define

    P_c(W,Z)=W dot H(c+Z).

Multinomial expansion and Wick contraction give EXACTLY

    E(a(c+Z)*P_c^m)=m!*T_m(a)(c).                          (9)

By (5), these are nonzero polynomials in c for EVERY m and every nonzero
a. In particular over C(c1,c2) the phase has NO zero pure moments at all.
The formal exponential sum at scalar coefficient one nevertheless equals
a(F^(-1)(c)), interpreted c-adically through (3). Replacing the sum by
termwise vanishing changes the hypothesis.

At c=0 there is a different, automatic phenomenon. For ANY polynomial
H of order>=2, with NO Keller hypothesis, P_0^m has W-degree m and
Z-degree at least 2m. For a monomial marker W^alpha Z^beta a nonzero
Wick contraction requires

    m+|alpha| >= 2m+|beta|.

Thus E(W^alpha Z^beta P_0^m)=0 whenever m>|alpha|-|beta|; in particular
every fixed polynomial marker has zero moments for m above its total
degree. Pure moments are all zero. Unshifted Gaussian safety here detects
the degree imbalance, BEFORE any Keller condition or inversion theorem.
It cannot be promoted to the required polynomial identities in c.

These two facts are compatible with GS1's translation argument: that
argument requires pure contraction vanishing as a polynomial BEFORE
passing to EVERY translated scalar phase. The needed hypothesis fails
here at every positive order. No Gaussian or Image conjecture is newly
refuted by (9).

## 5. Low-dimensional transfer and the homogeneous hypothesis

For F_tau=z-tau*H, a plane determinant calculation gives

    det DF_tau=1+tau*(tau-1)*div(H),                       (10)

because det(I-JH)=1 implies det(JH)=tr(JH)=div(H). By (4) at m=1,
div(H) has nonzero top -partial_2 phi. Thus this linear homotopy is NOT
Keller for general tau. It differs from the scaling family
tau^(-1)F(tau*z), which IS Keller but has tau-dependent H. Formal-parameter
specialization and vanishing of each scalar-power coefficient cannot
be interchanged.

For a homogeneous correction H of degree r>=2, in contrast, T_m(1) is
homogeneous of degree (r-1)m. Equation (3) with a=1 then forces each
positive term zero separately. This homogeneous separation is the classical
step in GS1, and cannot be used for an arbitrary mixed-degree plane map
while retaining its dimension and actual target. No reduction to cubic
homogeneity in the SAME plane dimension is supplied.

The new shared FIRFC2FULL20261116 factorial theorem and
GMCORDINARYCUTOFF20261116 Gaussian producer were screened at their CURRENT
ledger scopes, not imported as new theorems into this proof. The latter
decides moment ideals for supplied phases; it does not assert pure
vanishing for an arbitrary Keller phase. The former concerns the
functional L(y^alpha)=alpha! on TWO polynomial variables. Passing from
paired Wick variables to that functional requires diagonal extraction,
which is linear and NOT multiplicative. For example, if

    p=W1*Z2+W2*Z1,

then its diagonal projection is zero but the projection of p^2 is
2*y1*y2. Consequently L applied to powers of one projected phase does
not reproduce the original powers in general. No compatible algebra
homomorphism or automatic binary factorial reduction has been proved.

Primary context: Zhao,
[*New Proofs for the Abhyankar–Gurjar Inversion Formula...*, 0907.3991v2](https://arxiv.org/html/0907.3991v2),
Section3 Theorems3.1--3.5/Corollary3.6 and Section4
Lemma4.1/Propositions4.2--4.5 pertinent proofs were read. The paper
explicitly distinguishes det(I-JH)=1 from nilpotence and the entire linear
homotopy. Its Weyl-symbol and formal-logarithm foundations are NOT newly
reproduced here; the needed coefficient formula uses the retained GS1
residue proof instead. No new external theorem import or primary cache.
The older KPC provenance note already detects the FIRST nonzero contraction
for a particular mixed-degree three-dimensional normalization. That
distinction is prior; (5) proves ALL orders/ALL amplitudes for EVERY plane
extension under the existing preserving-B framing.

## 6. Verification and remaining objective

The new checker provides finite exact controls for (2), (4), (6), (8)--(10),
formal sums modulo a small source-degree cutoff, circular Wick
normalization, and diagonal-projection failure. It does not replay any
registered mathematical checker or infer (5) from a finite prefix.
The universal conclusion follows from the polynomial-in-z2 argument.

FULL THEOREM TARGET: prove or break full HC4/JC2.
CURRENT LAST UNSOLVED IMPLICATION remains

    ACTUAL polynomial plane source
       ==> Tr(Lambda_L*Lambda_L) subset B.                 OPEN

Full marked-energy integrality is OPEN. The new result screens a direct
low-dimensional moment transfer; it supplies no automatic regularity,
monic inverse equation, same-dimension homogeneous reduction, or actual
noninverse. All prior positive subclasses remain retained. Published
generic-degree-four PLANE exclusion remains. No independent replay,
formal verification or external review; review/outreach remain DEFERRED.
