# Full reductive Mathieu classification and recomposition failures

17 November 2026. Written deductions from the
[fixed original exceptional counterexamples](EXPLICIT_EXCEPTIONAL_PURE_PHASES_WITH_FIXED_BAD_MARKERS.md).
The acting group may now be disconnected. Every polynomial recomposition
of the exceptional phases has an explicit marker determined by its degree;
its ACTUAL ordinary marked resolvent is transcendental. The same subalgebra
has an exact purity/safety classification over every coefficient algebra.
External review remains deferred; independent replay and formal verification
are not claimed. No new degree campaign or finite-prefix inference is used.

## 1. Complete classification, including disconnected acting groups

Let G be a complex linear algebraic group with reductive identity component
K=G^0, let H be closed with reductive identity component, and put

    X=G/H,       D=[K,K].

For normalized Reynolds h on C[X], the complete criterion is

    ker h is Mathieu--Zhao  iff  D is contained in H.             (1)

There is no spherical or component-group restriction. The negative case
has an actual pure phase and ONE fixed original marker with infinitely
many nonzero marked moments.

**Proof.** The connected components are the finitely many open-and-closed
K-orbits X_i=K/H_i, with H_i=K intersect g_i H g_i^-1. The component
stabilizers have reductive identity components and are reductive in
characteristic zero. All components are affine. Their regular function
algebra is a finite product, and the normalized functional is

    h((f_i))=(1/r) sum_i h_i(f_i),                               (2)

where each h_i is normalized Reynolds. These are the component and
normalization arguments in Sections1 and3 of the retained
[component theorem](SPHERICAL_PAIRS_AND_COMPONENT_NULLCONES.md); neither
argument uses sphericity.

If D is not contained in H, the connected pair (K,K intersect H) is
non-toral. The now proved connected classification supplies a pure p
and a fixed bad q on that component. Extend both by zero on the other
components. This is an ORIGINAL regular function operation because the
component is open and closed. Equation (2) scales every moment by1/r,
so purity and the SAME infinitely many nonzero marked orders persist.

If D is contained in H, it is contained in every H_i because D is
characteristic in K and hence normal in G. Each component is a torus
T_i=K/H_i. Embed its character lattice in M=X^*(K/D). The retained
finite-trace theorem GVC2FTS gives

    sum_i CT(f_i^m)=0 for all m>=1
       implies CT(f_i^m)=0 for every i and m>=1.                  (3)

Duistermaat--van der Kallen then puts 0 outside the convex hull of each
nonzero component support. Separate each support strictly. For any fixed
original marker q=(q_i), sufficiently large m makes every CT(q_i f_i^m)
zero. There are only finitely many components, so one maximum cutoff
works in (2). This proves the Mathieu property. The proof also covers
eventual purity: apply the all-positive argument to f^N, then use the
finitely many markers q f^j for 0<=j<N. QED.

The positive proof is the existing nonspherical toral half of the
component theorem. The new exceptional failures remove the only
restriction in its negative half. No stronger EVERY-order nonzero
exceptional witness is required.

For clarity, the already proved moment-nullcone criterion remains separate:

    Z_X=N_G  iff  D subset H AND
    theta: direct_sum_i X^*(K/H_i) -> X^*(K/D) is injective.      (4)

Equation (4) is retained from the
[full nullcone classification](HAAR_SHEARS_AND_GENERAL_NULLCONE_CLASSIFICATION.md).
Here theta is defined after imposing D subset H.
The extra simultaneous-character condition is not part of (1).

## 2. The exact arithmetic input for the exceptional phases

Fix either original pair P,Q from the counterexample note. Use

    O: c=28, a=3;          I: c=66, a=5,
    b_j=h(Q P^j).

The preceding proof gives two uniform facts:

    b_0=0;
    denominator(b_j) divides (c*j-a)! for every j>=1;            (5)
    v_p(b_j)=-1 whenever p=c*j-a is prime.                      (6)

Here (5) follows directly from the integral zero-weight polynomial of
row degree c*j-a-1 and the original beta-integral monomial formula.
It is a bound for every order, not a numerical observation. In (6),
all relevant primes are greater than5. Also gcd(c,a)=1 in both cases.

## 3. A marker chosen by degree for EVERY polynomial recomposition

Let k be any characteristic-zero field containing the coefficients in
question, and let H(T) in k[T] have degree d>=1. Choose

    0<=e<d,       gcd(c*d,c*e-a)=1,       Q_d=Q P^e.             (7)

Such e always exists. For each prime ell dividing d but not c, choose
e different from a/c modulo ell. If ell divides c, it cannot divide a,
so no condition is needed. The Chinese remainder theorem combines these
choices modulo rad(d), and a representative can be taken below d.
Moreover gcd(c,c*e-a)=gcd(c,a)=1, giving (7). Taking the smallest such e
is a deterministic rule depending ONLY on d and the quotient.

We claim that the SAME original marker Q_d detects infinitely many
nonzero moments of H(P), for every H of degree d. No coefficient of H
is used to choose Q_d. The nonzero orders may depend on H.

First suppose H has rational coefficients, with leading coefficient h_d.
For any sufficiently large prime

    p=c*(d*m+e)-a,                                              (8)

the coefficients of H are p-integral and h_d is a unit. The leading
term of H(T)^m is h_d^m T^(dm). Expanding the marked moment gives

    h(Q_d H(P)^m)=h_d^m b_(dm+e)
                  + sum_{j<dm} [T^j]H(T)^m * b_(j+e).           (9)

The leading summand has valuation -1 by (6). Every other summand is
p-integral by (5), since c*(j+e)-a<p; b_0 is zero. Thus the whole
marked moment has valuation -1. The progression (8) is coprime by (7),
so Dirichlet gives infinitely many such primes and arbitrarily large m.

The passage to arbitrary complex, or other characteristic-zero, coefficients
does not require a split-prime assumption. The following elementary lemma
handles it.

**Coefficient valuation lemma.** If L is a finitely generated field over
Q and S is a finite subset of L^*, then for all but finitely many rational
primes p there is a valuation w of L extending v_p with every member of
S a unit. It can be normalized so that w(p)=1.

**Proof.** Choose a transcendence basis t with L finite algebraic over
Q(t). Give Q(t) the p-adic Gauss valuation, the minimum coefficient
valuation on polynomials, extended by subtraction to rational functions.
For each alpha in S choose monic equations for alpha and alpha^-1 over
Q(t). For all but finitely many p every coefficient of these finitely
many equations is Gauss-integral. Extend this valuation to L. Such an
extension exists by the domination lemma for valuation rings; see
[Stacks, Lemma10.50.2](https://stacks.math.columbia.edu/tag/00I8).
The monic equations force w(alpha)>=0 and w(alpha^-1)>=0. Hence
w(alpha)=0. An algebraic extension of this discrete valuation has
rational value group and admits the stated normalization. QED.

Apply the lemma to the field generated over Q by the finitely many
coefficients of H. At every sufficiently large p in (8), use such a w
with all nonzero coefficients units. Equation (9) again has a unique
summand of negative valuation, namely -1. Thus the original marked
moment cannot vanish. The polynomials P and Q_d have not changed.

This proves the recomposition theorem over EVERY characteristic-zero
coefficient field. In particular, over C,

    H(P) is eventually pure iff H(0)=0;
    every NONZERO such H(P) is unsafe, with fixed marker Q_d.   (10)

Indeed h(H(P)^m)=H(0)^m by all-order purity of P. The zero polynomial
is safe. Nonzero constant recompositions are not eventually pure.

For example, the powers P^3 on2O and P^5 on2I, not detected by the
original prime progression with e=0, are both detected with e=1 and
the fixed marker QP. This is an application of (7), not a degree search.

## 4. ACTUAL ordinary marked resolvents are transcendental

For every complex H of positive degree d, the preceding Q_d gives

    C_(Q_d)^(H(P))(z)=sum_{m>=0} h(Q_d H(P)^m)/z^(m+1)          (11)

transcendental over C(z). This is the actual convergent Cauchy germ
outside a sufficiently large disk, because H(P) is bounded on the
compact Haar cycle. We now prove the transcendence assertion rather
than inferring it just from an infinite tail.

We use a finite-generation form of the Eisenstein argument. If L has
characteristic zero and A(t) in L[[t]] is algebraic over L(t), then all
coefficients of A belong to a finitely generated Z-subalgebra of L.
Here is a proof sufficient for this application.

Take F(t,Y) in L[t,Y] with F(t,A)=0 and F_Y(t,A)!=0, and let
r=ord_t F_Y(t,A). Choose N>r, put a(t)=A modulo t^N, and write
A=a+t^N B. Taylor expansion shows that

    G(t,Y)=t^(-N-r) F(t,a(t)+t^N Y)

is a polynomial in L[t,Y]. Its reduction modulo t is linear, with
nonzero Y coefficient u. To see divisibility of its constant term,
use F(t,A)=0: F_Y(t,a) has order r and every Taylor term of order
at least2 has order at least2N>N+r. Every coefficient of Y^i for
i>=2 vanishes modulo t after the division. Therefore recursive
coefficient comparison in G(t,B)=0 expresses all coefficients of B
in the ring generated by the finitely many coefficients of G, u^-1,
and the coefficients of a. This proves the assertion. It is the
finite-generation consequence of
[Eisenstein's theorem](https://arxiv.org/abs/2412.11738); no analytic
convergence theorem or bounded degree calculation is needed.

Let L be generated over Q by the coefficients of H. If the series in
(11) were algebraic over C(z), its formal series would be algebraic
over L(z). Indeed, for any fixed degree bounds a polynomial relation
is a homogeneous linear system over L in finitely many unknown
coefficients. Its rank is unchanged by extending L to C. An infinite
system with finitely many columns has the rank of a finite subsystem,
so a nonzero complex solution gives a nonzero L solution.

Apply the finite-generation assertion to

    A(t)=sum_{m>=0} h(Q_d H(P)^m)t^m.

Choose a finite set S_A of nonzero elements of L whose Z-algebra contains
all coefficients of A, as above, and combine it with the nonzero
coefficients of H. For every sufficiently large prime p in (8), the
coefficient valuation lemma supplies ONE w making all these generators
units. Then all coefficients of A are w-integral. But (9) proves that
its m-th coefficient has valuation -1. This contradiction proves (11).

In particular the original unit germ is1/z for H(0)=0, whereas this
specific ordinary marked germ is transcendental. Algebraic exponential
generating functions cannot replace it. This establishes an actual
failure of the proposed ordinary-algebraicity route on these families.

## 5. Exact purity and safety on R[P] over EVERY coefficient algebra

Let R be any commutative unital C-algebra, and extend h R-linearly to
R tensor C[X]. For H(T)=sum_{j=0}^d c_j T^j in R[T], one has

    h_R(H(P)^m)=c_0^m,
    H(P) is all-positive pure iff c_0=0,
    H(P) is eventually pure iff c_0 is nilpotent,
    H(P) is safe for ALL original R-markers iff EVERY c_j
       is nilpotent iff H(P) is nilpotent.                     (12)

The unit formula proves the first three statements. If every c_j is
nilpotent, their finitely generated coefficient ideal is nilpotent:
if c_j^(r_j)=0, its power 1+sum_j(r_j-1) is zero. Thus H(P) is
nilpotent, which proves safety for all markers at one common cutoff.

Conversely, if some coefficient is not nilpotent, choose a prime ideal
of R where it survives and pass to the fraction field of R/prime.
If the resulting H has positive degree, Section3 supplies the explicit
original marker QP^e with infinitely many nonzero moments in that
field. The same moments in R must be nonzero. If the resulting H is
a nonzero constant, the unit marker already gives nonzero moments
in every order. In either case H(P) is not safe. The zero ring causes
no exception. This proves all of (12), including the final equivalence.

This classifies the complete polynomial subalgebra R[P] in all degrees
and with nilpotents retained. It uses the existing
[family-transfer framework](FULL_MIXED_HAAR_FAMILIES_AND_TORSION_NORMALIZATION.md)
without asserting a separate nonreduced failure gate or replacing the
full pure ideal by its radical. Over C its nonzero pure elements are
ALL unsafe; over R the safe elements are exactly the nilpotent ones.

## 6. Status and scope

The full exceptional decisions and the complete reductive classification
are now recorded as internal written theorems. Earlier moment-nullcone,
capacity, automorphism, primitive-fibre and family-transfer results remain
valid. The new counterexamples make the proposed universal exceptional
pure-to-safe and universal corrected-coupling-vanishing implications false.

No minimal counterexample degree, EVERY-order nonzero exceptional tail,
explicit finite-fibre graph, or numerical corrected cycle pairing is
claimed. These are separate follow-up questions, not missing conditions
for the full Mathieu decision or classification (1). Other Gaussian and
Image programmes remain under their existing goals. External review is
still deferred at the user's request.
