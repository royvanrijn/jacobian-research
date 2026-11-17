# Hopf purity, homogeneous section weights and the real-height obstruction

15 November 2026. Internal written proofs with exact coefficient controls;
no novelty, independent specialist replay, external review or formal
verification is claimed. See the [checker](verify_hopf_height_criterion.py),
[receipt](hopf_height_verification.json) and
[source audit](hopf_height_source_audit.json).

The [regular Bezout family](REGULAR_BEZOUT_MAPS_AND_WEIGHTED_HAAR_TAILS.md)
already supplies a fixed infinite marked tail. This note gives structural
information about its remaining **unweighted** purity condition:

1. For any norm-one Hopf triple with opposite torus weights, purity is
   exactly the uniform-interval even-moment condition on one height.
   If it holds, the unweighted marker u itself has the positive bad tail.
2. For homogeneous exceptional section parameters, all but five possible
   weights are rejected by the same known nonzero first pure moment,
   without any bound on the degree of the parameter.
3. The entire torus-equivariant branch with a real compact height is
   excluded by a local density and global sphere argument.

The five remaining weights, mixed-weight sections and complex heights
are not classified. Both exceptional MZ questions remain **OPEN**.

## 1. A complete all-order criterion for a norm-one torus Hopf triple

Let A be a commutative C-algebra with a rational C^* action and a
normalized invariant functional h. Suppose

    u v+tau^2=1,
    wt(u)=a, wt(v)=-a, wt(tau)=0, a!=0.

Put

    P=(1+u)(v-(2+u)tau^2).

All products below are identities in A. Introducing a Laurent symbol
u and using uv=1-tau^2 is legitimate for computing the weight-zero
coefficient: every balanced u^k v^k becomes (1-tau^2)^k. Equivalently,
the identities follow in the universal graded quadric algebra and
then specialize to A. Write

    C_m(T)=sum_(j=0)^m (-1)^j binom(m,j) binom(m+2j,m) T^(2j),
    B_m(T)=sum_(j=0)^m (-1)^j binom(m,j) binom(m+2j,m-1) T^(2j).

The weight-zero parts of P^m and u P^m are respectively C_m(tau)
and B_m(tau). Indeed

    P=(1+u)/u * (1-(1+u)^2 tau^2),

so each coefficient is an ordinary finite binomial coefficient.
The degree of C_m is2m and its leading coefficient is
(-1)^m binom(3m,m), which never vanishes.

Define ell_0(q)=(1/2) integral_(-1)^1 q(T)dT. For every m>=1,

    ell_0(C_m)=0,
    ell_0(B_m)=c_m=4^m(m!)^2/(2m+1)!>0.                     (1)

Here is a direct all-order proof. Set

    A_m(q)=integral_0^q (1-s^2)^m ds.

The derivative of A_m has a zero of order m at q=1, so
A_m(1+u)=A_m(1)+O(u^(m+1)). Integrating the Laurent expression
for P^m in T from0 to1 gives

    u^(-m)(1+u)^(m-1) A_m(1+u).

Its constant coefficient is zero: the constant term A_m(1) is
multiplied by a polynomial of degree m-1, and the remainder starts
at u^(m+1). Multiplication by u instead selects coefficient m-1,
which is A_m(1)=c_m. This proves both identities in (1).

Since1,C_1,C_2,... form a triangular basis of the even polynomials,
we obtain the exact equivalence

    h(P^m)=0 for all m>=1
      iff h(tau^(2j))=1/(2j+1) for all j>=0.                 (2)

No reality hypothesis is used in (2). Whenever it holds, every B_m
has its uniform value, and hence

    h(u P^m)=c_m>0 for every m>=1.                           (3)

Thus this architecture's purity already forces an unweighted fixed
bad multiplier. No Jacobian weight is needed in (3) **after** the
height condition has been proved. The prior weighted-transfer theorem
remains a different statement: it supplies its marked tail for all
section parameters before purity is known.

The height condition can also be written as the formal-series identity

    h(1/(1-t tau^2))=sum_(j>=0) t^j/(2j+1)
                         =artanh(sqrt(t))/sqrt(t).           (4)

Equation (4) is a prescribed logarithmic period, not a finite prefix
test or a claim that the actual height has real values.

## 2. Application to every torus-equivariant regular section

Use the original left torus on SL_2/Gamma. In the Bezout notation,

    wt(U)=m, wt(V)=n, wt(F)=-n, wt(G)=-m, L=n+m.

For a homogeneous b of weight L, the regular map

    Phi_b=[[U+Fb,V+Gb],[F,G]]

is equivariant for the target torus of weights(m,n,-n,-m). Its
Hopf triple

    u_b=-2(U+Fb)(V+Gb), v=2FG,
    tau_b=(U+Fb)G+(V+Gb)F=tau_0+2FG b

has weights L,-L,0 and satisfies u_b v+tau_b^2=1. Therefore its
Hopf pullback is pure exactly when

    h_Gamma(tau_b^(2j))=1/(2j+1) for every j>=0.              (5)

If (5) holds, its unweighted regular marker u_b already gives the
tail c_m in (3). This is an exact all-degree criterion for this
entire equivariant section family, not a proof of existence of b.

## 3. Five possible section weights, with no degree restriction

Let b be any nonzero homogeneous source function of weight beta;
it need not have weight L. On the target let delta=F partial_X+
G partial_Y. Decompose the Hopf phase as

    P_H=P_(-1)+P_0+P_1+P_2,
    P_(-1)=v,
    P_0=uv-2tau^2=1-3tau^2,
    P_1=-3u tau^2,
    P_2=-u^2 tau^2.

These have weights kL and delta-depths0,2,4,6 respectively. The
canonical Phi_0 is torus-equivariant. Finite Taylor substitution gives

    P_b=sum_(k=-1)^2 sum_j b^j/j! *Phi_0^*(delta^j P_k).

The weight of each summand is

    kL+j(beta-L).                                           (6)

For j>=1, it can be zero only if

    beta/L in {-1,0,1/3,1/2,3/5,2/3,3/4,1}.                (7)

The centre of SL_2 acts trivially on both exceptional quotients, so
every source weight is even. Intersecting (7) with2Z leaves exactly

    2O: beta in {-20,0,10,12,20},
    2I: beta in {-32,0,16,24,32}.                            (8)

For every other homogeneous b, all additional summands have nonzero
weight. Consequently h(P_b)=h(P_H o Phi_0), whose exact retained values are

    -1758938/118915225 on2O,
     97777523/2181764186 on2I.

They are nonzero, so the entire class is excluded from purity. This
uses one structural first-moment calculation for all polynomial
degrees. The five weights in (8) are only necessary candidates;
none is certified pure. For mixed-weight b, products between its
different components can have the weights in (6), so (8) cannot be
applied componentwise as an exclusion theorem.

## 4. A bounded scalar density on a sphere forces two critical points

**Lemma.** Let f:S^2->R be smooth, and let mu be a measure with a
smooth strictly positive density. If f_*mu has an essentially bounded
density with respect to Lebesgue measure, f has exactly two critical
points, one nondegenerate minimum and one nondegenerate maximum.

Proof. Bounded density would give, uniformly in v and small epsilon,

    mu{|f-v|<=epsilon}<=C epsilon.                           (9)

At a degenerate critical point, diagonalize the Hessian. If its rank
is zero, Taylor's theorem shows that a coordinate square of side
c epsilon^(1/3) lies in this band. Positive source density gives a
lower bound c' epsilon^(2/3), contradicting (9). If its rank is one,
use a rectangle of sides c epsilon^(1/2) in the nondegenerate
direction and c epsilon^(1/3) in the null direction. The quadratic
term and all cubic remainders are O(epsilon), so its mass is at
least c' epsilon^(5/6), again a contradiction. This also excludes
nonisolated critical sets, since the argument applies pointwise.

At a nondegenerate saddle, the Morse lemma gives f-v=xy in local
coordinates. The area of |xy|<=epsilon in a fixed small square is
bounded below by c epsilon log(1/epsilon), contradicting (9).
Thus every critical point is a nondegenerate extremum. The critical
set is finite by compactness and nondegeneracy. The Morse Euler
identity on S^2 is

    #minima-#saddles+#maxima=2.

There are no saddles and there is at least one maximum and minimum,
so there is exactly one of each. This finishes the proof.

Only the Morse lemma and Euler critical-point identity are classical
imports; see [Milnor, Lemma2.2 and Sections3,5](https://www.maths.ed.ac.uk/~v1ranick/papers/milnmors.pdf).
The local lower bounds above give the required density obstruction;
they do not assume a generic or Morse input function.

**Rotation-group consequence.** If a noncyclic finite subgroup of
SO(3) preserves f, its pushforward cannot have bounded density. If
it did, uniqueness would make the group fix each of the two extrema.
A finite rotation group fixing a point is contained in the rotations
about its axis, and is cyclic. This is a contradiction.

## 5. The real compact-height branch is closed in every degree

Realize Gamma unitarily and use the compact left torus in SU_2.
A left-torus-invariant function on SU_2/Gamma lifts to a function
on T\SU_2=S^2 invariant under the rotation group Gamma/{+1,-1}.
Haar descends to the usual smooth positive sphere measure. For2O
and2I those rotation groups are noncyclic.

Suppose tau_b in (5) is **real-valued on this compact core** and
satisfies (5). Let nu be its positive pushforward probability measure.
The symmetric measure (nu+(-1)_*nu)/2 has all even moments1/(2j+1)
and all odd moments0. Compact real moment determinacy, by polynomial
approximation on a containing interval, identifies it with uniform
measure on[-1,1]. Positivity then implies

    nu<=2*Uniform[-1,1],

so nu has a bounded density. The rotation-group consequence of the
lemma excludes this. Thus no torus-equivariant regular section with
a real compact tau_b can produce a pure Hopf phase, in any degree.

This also proves the broader geometric statement: no smooth real
polyhedral-invariant height on S^2 has the uniform even moment
sequence. It applies to all noncyclic finite rotation groups.

The real condition concerns the auxiliary tau_b, not the generally
complex Hopf phase. It is still a genuine restriction. Complex
heights can have the same analytic moments as real interval measure
without the same positive pushforward measure, so the density proof
does not apply to them. Nothing here proves that purity forces
tau_b real or conjugate to a real function. The canonical frozen
heights are not even assumed real; their existing first-moment
exclusion remains a separate exact result.

## 6. Postcomposition preserves a weighted gate, not ordinary purity

A current related source is [Fasel--Hornslien, Section2.6](https://arxiv.org/html/2604.12541v2#S2.SS6),
revised5 November2026. It supplies explicit polynomial quadric
self-maps representing prescribed motivic Brouwer degrees. No
exceptional MZ result or Haar-preservation assertion is imported.

For an ordinary regular self-map psi of SL_2, let d be its actual
scalar action on H^3 and write psi^*Omega_*=J_psi Omega_*.
Pullback of top-period classes proves, for all regular a,

    h_*(J_psi*(a o psi))=d*h_*(a).                           (10)

Thus for d!=0, postcomposing Phi_b supplies a normalized weighted
transfer with weight W_b*(J_psi o Phi_b)/d. This adds the new
Jacobian, rather than removing it. The actual H^3 scalar, generic
field degree and a motivic degree label must not be substituted for
one another without the relevant comparison. No such postcomposition
was asserted to make an unweighted pure phase here.

The remaining construction space is explicit: the five homogeneous
weights in (8), complex heights in the equivariant branch, and genuine
mixed-weight sections. These remain sufficient-family targets within
the broader full admissible pure-null classification, not replacements
for the original exceptional MZ objective.
