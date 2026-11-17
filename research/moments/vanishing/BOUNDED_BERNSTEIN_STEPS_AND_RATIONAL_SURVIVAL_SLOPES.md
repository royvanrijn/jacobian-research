# Bounded Bernstein steps and rational first-survival slopes

16 November 2026. Internal written theorem in every dimension and degree,
over every characteristic-zero field, for arbitrary formal symbols.
Independent replay, formal verification and external review remain false.
No novelty is claimed. The [checker](verify_rational_survival_slopes.py),
[receipt](rational_survival_slopes_verification.json) and
[source audit](rational_survival_slopes_source_audit.json) separate exact
transcription controls from the universal proof.

The candidate bounded-step and lattice arguments came from the user's
parallel research thread, preserved in its
[unregistered draft](../../jacobian/cross-programme-connections/BOUNDED_BERNSTEIN_STEPS_AND_RATIONAL_APOLAR_SLOPES_2026_11_16.tex)
and [audit](../../jacobian/cross-programme-connections/bounded-bernstein-rational-slopes-audit-2026-11-16.json).
Here they are audited against the actual module and primary Bernstein
theorem. The finite parallel-line conclusion below is an additional
deduction. The earlier field and all-ring formal closures remain separate,
unchanged proof packets.

## 1. The full first-survival theorem

Let K have characteristic zero, let P be a nonzero polynomial of degree
D>0, and let lambda be a nonzero member of (xi)K[[xi_1,...,xi_n]]. Put
T=lambda(partial_x), r=ord(lambda)>0, and

    tau_k=min{j>=0:T^k(P^j)!=0}, min(empty)=infinity.

Exactly one of the following holds.

* Some tau_N is infinite, and every tau_k with k>=N is infinite.
* Every tau_k is finite. There are an integer J, a rational rho>0 and
  a constant B with, for all sufficiently large k,

      0<=tau_(k+1)-tau_k<=J,
      |tau_k-rho*k|<=B,             r/D<=rho<=J.          (1)

More precisely write rho=a/b in lowest terms, a,b positive integers.
There is a finite set C of integers such that, for every sufficiently
large k,

      b*tau_k-a*k belongs to C.                          (2)

Each retained line b*Y-a*X-c=0, c in C, is a factor of one nonzero
rational polynomial satisfied by the actual integer pairs (k,tau_k).
Constants used only finitely often can be omitted from C. This does
not assert an eventually periodic selection among those parallel lines.
No purity or eventual marker-safety hypothesis is used in (1)-(2).

## 2. An ordinary Bernstein equation on the actual module

Use the actual holonomic image from the
[generic Bernstein proof](GENERIC_BERNSTEIN_EQUATIONS_AND_FORMAL_GAP_CLOSURE.md):

    R=K[[xi_1,...,xi_n,eta]], D_R=R<partial_xi,partial_eta>,
    E=K[x][[eta]], M=D_R*1 subset E,
    xi_i=partial_x_i+eta*P_i,
    partial_xi_i=-x_i, partial_eta=partial_eta(actual)+P.

The retained proof establishes that these formal coefficient actions
are coefficientwise well-defined, M is holonomic, and

    v_k=lambda(xi)^k*1=e^(-eta*P) T^k(e^(eta*P)),
    ord_eta(v_k)=tau_k.                                  (3)

It also establishes a nonzero p in K[X,Y] with
p(k,tau_k)=0 at every sufficiently large finite tau_k.
These inputs concern the actual cyclic image M, not all completed
polynomial multipliers or a substitute scalar time source.

Apply Mebkhout--Narvaez-Macarro, Theorem3.1.1, to M, its element1,
and the coefficient function lambda in R. Formal power series rings
satisfy its hypotheses, as explicitly covered in Section4.2. The result is

    B(s)*lambda^s*1=Q(s)*lambda^(s+1)*1,
    0!=B in K[s], Q in D_R[s],                            (4)

in M_lambda[s]*lambda^s. Both the differential order and parameter
degree of Q are finite. The primary source is
[Mebkhout--Narvaez-Macarro, 1991](https://doi.org/10.24033/asens.1627),
printed236-237 and241. This is its ordinary equation for lambda,
not the second, eta-exponent equation producing p.

Specialization initially gives (4) only in M_lambda. We need an
identity in M on the ORIGINAL powers. Expand the finite operator by
Leibniz and write the difference of (4) as lambda^(s-A0)*u(s),
with one A0>=0 and u in M[s]. There are finitely many coefficients.
Their vanishing after localization gives one A1>=0 such that
lambda^A1*u(s)=0 in M[s]. For k>=A0+A1 the actual difference is
lambda^(k-A0)*u(k)=0. Consequently

    B(k)*v_k=Q(k)*v_(k+1) in M                           (5)

for all sufficiently large integers k. The existing positive exponent
k-A0 kills the localization torsion. There is no replacement of v_k
by v_(k+A1), and no shift of the desired sequence.

Write Q in left-coefficient PBW order and let J be its largest
partial_eta exponent. Every coefficient in R preserves eta-order.
Each partial_xi_i acts by -x_i and preserves eta-order. Each
partial_eta acts by actual differentiation plus multiplication by P
and lowers eta-order by at most1. Therefore

    ord_eta(Q(k)*v_(k+1))>=tau_(k+1)-J.

After excluding the finitely many integer roots of B, (5) gives the
upper step bound in (1). Monotonicity gives its lower bound:
T^(k+1)P^j!=0 implies T^kP^j!=0. If some v_N is zero, all later
ones are zero instead; otherwise the above argument applies to all
late indices. In the latter case J cannot be0, since the positive
operator order gives tau_k>=r*k/D.

This argument uses finite operators on the actual module. It assumes
neither a scalar time differential equation nor a product law for
the two-parameter apolar ideals.

## 3. A polynomial curve with bounded integer steps

We prove the abstract statement used here. Suppose a_k are integers,
0<=a_(k+1)-a_k<=J eventually, and p(k,a_k)=0 eventually for some
nonzero p in Q[X,Y]. Then

    a_k=rho*k+O(1), rho rational and nonnegative.          (6)

Bounded steps imply a_k=O(k). Remove repeated factors and factors
depending only on X, and pass beyond their finitely many exceptional
real X values and the discriminant zeros. Over this last real
interval the real curve branches are graphs Y=f_i(X). The classical
convergent Newton--Puiseux expansions at infinity give the following
facts. Branches growing faster than linearly cannot contain a_k.
The remaining branches have bounded derivatives. For each pair,
their difference is either bounded or tends in absolute value to
infinity. No branch expansion or analytic assertion about the formal
symbol lambda itself is involved.

Group the remaining branches by bounded difference. This is an
equivalence relation. Distinct groups cannot exchange the selected
branch at consecutive sufficiently large integer indices: their
separation tends to infinity, whereas the bounded derivatives and
the step bound keep such a change uniformly bounded. Thus one fixed
group supplies a_k for EVERY sufficiently large integer k.

A nonlinear algebraic graph has, eventually, a nonzero second
derivative of fixed sign. Indeed its second derivative is algebraic;
a nonzero algebraic function has only finitely many zeros on this
last nonsingular branch interval. Such a branch is strictly convex
or concave and, by its bounded first derivative, has arc length O(N)
over a fixed starting point through X=N.

Here is the needed elementary lattice estimate. If L lattice points
lie on one such arc in increasing X order, their L-1 consecutive
chord vectors are distinct nonzero integer vectors: strict convexity
or concavity makes their slopes strictly monotone. At most O(R^2)
integer vectors have Euclidean norm at most R. Hence the total norm
of L-1 distinct vectors is at least c*(L-1)^(3/2), for an absolute
c>0. The total chord length is at most the arc length O(N).
Therefore the branch has only O(1+N^(2/3)) lattice points up to N.
This proof applies to these open graph arcs without any closed-curve
or curvature-radius hypothesis.

The classical provenance is
[Jarnik, 1926, Section1 and beginning of Section2](https://www.dml.cz/handle/10338.dmlcz/500690),
printed500-502. His first theorem concerns arcs, and his second concerns
closed curves. The weaker bound needed here is proved in the preceding
paragraph; his sharp constants and full asymptotic are not imported.

An affine line not having both rational slope and rational intercept
contains at most one lattice point: two such points determine both.
If our fixed group contained no rational affine line, finitely many
nonlinear branches and nonrational lines would supply O(N^(2/3))
points. They must supply N-O(1) distinct points (k,a_k). This is
impossible. Some rational line rho*X+c belongs to the group, and
every selected branch differs from it by O(1), proving (6).
Monotonicity forces rho>=0.

For the constraint over arbitrary K, expand its finitely many
coefficients in a basis of their Q-linear span. Evaluation at integer
pairs forces each rational component to vanish. Select a nonzero
component and apply the preceding proof. Thus no ordering or embedding
of K into the real numbers is required.

Finally, write rho=a/b in lowest terms, with b>0. Bounded error and
integrality put b*a_k-a*k in a FINITE set C of integers. For c in C,
the polynomial p(X,(aX+c)/b) is either identically zero or has only
finitely many integer roots. Pass beyond all roots of the finitely
many nonzero such polynomials. Every remaining selected c therefore
gives a line factor b*Y-a*X-c of p. This proves the stronger finite
parallel-line statement. Combining it with Section2 and tau_k>=rk/D
proves Sections1's assertions.

## 4. Running apolar heights, passive suspension and safety

Define the ACTUAL heights and their running maximum by

    h_m=max{k>=0:T^kP^m!=0}, H_m=max_(0<=j<=m)h_j.

Every height is finite since r>0. We do not assume that T^kP^j is
nonzero for all j>=tau_k. What holds exactly is

    tau_k<=m iff some j<=m has h_j>=k,
    H_m=max{k:tau_k<=m}.                                 (7)

In the all-finite case, inversion of (1) gives

    H_m=m/rho+O(1),
    limsup_(m->infinity) h_m/m=1/rho,
    liminf_(m->infinity)(m-h_m)/m=1-1/rho.                 (8)

For the limsup, H_m>=h_m gives one inequality. If eventually
h_j<=(L+epsilon)j, then H_m<=max(C,(L+epsilon)m) for a fixed early
maximum C, giving the other. The running maximum is also exactly
the height of the passive suspension P+u: in (P+u)^m the different
powers u^(m-j) prevent cancellation between T^kP^j. No bounded-error
law for the ORIGINAL unsuspended h_m is claimed.

In fact the qualitative all-marker safety alternative, without
assuming the first pure equations, is

    all polynomial Q have T^m(QP^m)=0 eventually
       iff some tau_N=infinity OR rho>1.                 (9)

For the nilpotent alternative, h_m<N for every m; the retained
[formal commutator lemma](FORMAL_SHIFTED_POWERS_AND_APOLAR_NILPOTENCE_GAPS.md)
kills all marked outputs for m>=N+deg Q. For rho>1 use
h_m<=m/rho+B and delta=1-1/rho>0. The same lemma kills the output
whenever delta*m>deg Q+B. This is one affine degree cutoff.

Conversely Q=1 first gives T^jP^j=0 for j>=M. For k>=j>=M it
follows that T^kP^j=0. The finitely many powers j<M are also killed
for large k by positive operator order. Testing the other markers
Q=P^d for the finitely many shifts 0<=d<=d0 now proves
tau_k>k+d0 eventually. Since d0 is arbitrary, tau_k-k tends to
infinity; this is the retained shifted-powers equivalence with its
quantifiers intact. The law (1) then forces rho>1. This proves (9).

Under ALL positive pure identities, tau_k>k and hence rho>=1.
If rho=1, (2) says tau_k=k+c_k for one finite set of POSITIVE
integer offsets. At least one fixed d>=1 occurs infinitely often.
For this original polynomial marker Q=P^d,

    T^k(QP^k)=T^kP^(k+d)!=0 infinitely often.             (10)

Thus critical rational slope1 produces an actual fixed bad marker.
If rho>1, or the nilpotent alternative holds, safety follows and
the retained full formal field theorem gives a common LINEAR degree
cutoff on the pure locus. This strengthens its structure theorem;
it does not prove that every pure input is safe. The known unsafe
pure seeds remain in the critical-slope alternative.

Zero symbol, zero phase and constant phase cases are immediate.
If lambda has nonzero constant term, T is invertible on every finite
polynomial-degree space, tau_k=0 for nonzero P, and all-marker safety
fails already for Q=1. The positive-slope assertion requires the
stated lambda(0)=0. No rational-slope statement over arbitrary
coefficient rings is made; the earlier all-Q-algebra affine-bound
theorem remains unchanged.

## 5. Controls and exact limits of the conclusion

For P=xy and T=partial_x^N+partial_y^M, N,M>=1, the derivative terms
have distinct exponent vectors and positive integer coefficients.
Consequently

    tau_k=min_(0<=i<=k) max(N*i,M*(k-i)),
    rho=N*M/(N+M).                                      (11)

Rounding i=M*k/(N+M) proves bounded error directly. Adding one
derivative summand proves a step bound max(N,M). In the pure safe
examples M=N+1,N>=2, the reduced denominator is 2N+1. These are
unbounded at FIXED phase degree2. No denominator bound depending
only on phase degree is possible. The ordinary unsafe pure DVG seed
is also retained; selected small first-survival controls have
tau_k=k+1 and prevent mistaking rationality for pure-to-safe GVC.

The abstract sequence a_k=k+indicator(k is a square) is monotone,
has jumps at most2 and satisfies (a_k-k)(a_k-k-1)=0. Its intercept
is not eventually periodic. It satisfies even the finite parallel-
line conclusion, so no periodicity can follow from our abstract
hypotheses alone. It is not claimed to be an actual apolar sequence.
Without bounded steps, the polynomial-curve sequence a_k=k^2
shows that monotonicity and a polynomial relation do not imply a
linear law. The earlier nonmonotone square-index control is also
preserved as an obstruction to omitting the actual step hypothesis.

The proof gives neither a general algorithm producing Bernstein
certificates for arbitrary infinite formal input, nor a least
numerical marker cutoff or a universal denominator. It does not
revive raw-time ODEs, product Rees models or quasi-superadditivity.
Both complete exceptional Haar MZ cases remain open. Unrestricted
pure-to-safe GVC remains false with its existing exact counterexamples;
no Mathieu conclusion follows automatically.
