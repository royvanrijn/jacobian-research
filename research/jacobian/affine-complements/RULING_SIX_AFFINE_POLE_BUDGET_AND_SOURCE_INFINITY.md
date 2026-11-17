# Affine responses bound the complete ruling-six pole degree and source infinity

PF2RULED164. 17 November 2026. Actual complex polynomial plane Keller maps.

For the ENTIRE complete quadratic ruling-six class, affine first
parameter response improves the local pole bound by a strict inequality.
The resulting complete pole ledger forces every intermediate point over
u=infinity to have u-index at most two. On the original source curve,
every point over t=infinity is consequently unramified for the projection
t. If none of the pole divisors has a quartic face, the intermediate
projection is itself unramified over infinity.

The theorem supplies necessary constraints on arbitrary coefficients.
The later [simple-escape theorem](RULING_SIX_SIMPLE_ESCAPE_GRAPHS_AND_QUARTIC_FACE_REDUCTION.md)
excludes B and C. The [strict infinity ledger](RULING_SIX_RETAINED_SOURCE_PACKETS_AND_STRICT_INFINITY_BOUNDS.md)
strengthens the bounds using a retained source prime over each finite
image. Face A and JC2 remain open; the inequalities below remain valid.

## 1. Actual setup and complete-end inputs

Use the complete marked quadratic source image and notation of
[PF162](RULING_SIX_FINITE_END_COVERAGE_AND_QUADRATIC_BOUNDARY_ALGEBRA.md):

    V: w^2=N+M*z, N=u*Sigma^2, M=u*A,
    u=t^2, w=t*(Sigma+A*x), z=2Sigma*x+A*x^2,
    m=deg M>=2.                                             (1)

The actual target B=C[P,Q] lies in O(V), E=C(u,w), L=C(x,t),
F=C(P,Q), and t^2=u. Normalize the actual source volume as

    dP wedge dQ=du wedge dw/M.

Let P be ANY actual target coordinate of ruling degree six:

    P=b3*z^3+w*(a2*z^2+a1*z+a0)+b2*z^2+b1*z+b0, b3!=0.      (2)

All coefficient polynomials have arbitrary degree. Put

    D=deg b3, delta_c=ord_c b3, b=[E:F].                     (3)

The quantity b is the INTERMEDIATE generic map degree, not ruling
degree. The original generic map degree is 2b.

For an algebraic closure k of C(p), take all smooth complete
geometric generic components P=p. Their nonzero coframes are

    omega=dQ=-partial_p(w)*du/M.                            (4)

The sum of the degrees of Q on those components is b. Indeed
E/F is finite of degree b and base change from C(P) to k preserves
the dimension of its finite generic algebra; if geometric components
split, their degrees add. No geometric-integrality assumption is used.

The complete-end audit in PF162 gives the following indispensable
inputs. Every affine point and every M-root end is holomorphic.
The leading coefficient b3 is nonconstant, so D>=1, and every
u-infinity end is holomorphic. A pole can only occur at a finite
c off M with specialization degree two. There are exactly two
finite generic roots over that c and four escaping roots, counting
ramification. Any pole lies in one of the following weighted faces:

    A: ell=4, H0=2h, h=alpha*L, delta_c=4L, L>=3.
    B: ell=3, H0=2h, h=alpha*L, delta_c>=4L, L>=4.
    C: ell=3, H0=(3L-B0)*alpha,
       delta_c=B0+3L, 4<=B0<L.                              (5)

Here u-c=s^alpha on a normalized branch, w=c0*s^-h+...,
and ell is the nonzero initial root multiplicity. L is an integer.
For B, delta_c=4L when the face has degree six; delta_c>4L
when its last factor is constant and the sixth coefficient lies
strictly above the face.

There is at most ONE pole-carrying cluster at each c: each such
cluster needs at least three of the four escaping roots. Within
that cluster the sum of the indices alpha of all pole points is
at most ell. This is a count of normalized branches with their
indices, not a count of Puiseux conjugates as separate points.

The all-degree polynomial-mate argument
[PF159, Section 1](POLYNOMIAL_MATE_RESPONSES_AND_QUADRATIC_RULING_FIVE_ESCAPE_COVERAGE.md#1-a-polynomial-mate-forces-an-affine-first-parameter-response-at-a-pole)
also gives an indispensable condition: at EVERY actual pole of
omega the later first parameter response is AFFINE.

## 2. Affine response makes the multiplicity inequality strict

Write the first nonconstant coefficient of a branch as

    w=w0+s^r*(lambda+O(s)),  w0 in C[s,s^-1],
    lambda not in C.                                      (6)

The entire polynomial P* after substitution has initial polynomial
H(Y), with H(lambda)=p. Every coefficient has nonnegative s order,
and its reduction has degree one at an actual pole, because H is
affine. This is the polynomial-mate conclusion, rather than a
simple-root assertion about the earlier weighted face.

The coefficient of Y^ell in P*(w0+s^r*Y) has EXACT order

    ell*(r+h)-H0.                                         (7)

The leading ell-th derivative of the weighted initial polynomial
at c0 is nonzero. Every later term of the p-independent truncation
has higher order, so it cannot change that leading derivative.

Since ell>=3 and H is affine, coefficient (7) cannot have order
zero. Negative order would also contradict the entire initial
equation P*=p. Therefore

    ell*(r+h)-H0>0, r>H0/ell-h.                           (8)

The weaker non-strict inequality does not use the actual mate
fully. This strictness is essential for the bounds below.

The constants of the extended partial_p derivation on k are C.
Hence partial_p(lambda)!=0, and the index r in (6) is the EXACT
order of partial_p(w), not merely a lower estimate.

At a pole put q=-r-alpha>=1. Then ord omega=-q-1 and the
primitive Q has a pole of EXACT order q. From (5),(8),

    A: q<alpha*(L/2-1);
    B: q<alpha*(L/3-1);
    C: q<alpha*(B0/3-1).                                  (9)

## 3. Local and global complete pole budgets

Let beta_c be the sum of the pole orders of Q at all points over
c, over all geometric components. Sum (9) over the pole cluster.
Its total alpha is at most ell; all displayed multipliers are
positive whenever a pole is possible. Since beta_c is an integer,
strictness gives

    A: beta_c<=2L-5=delta_c/2-5;
    B: beta_c<=L-4<=floor(delta_c/4)-4;
    C: beta_c<=B0-4<=floor(delta_c/4)-4.                    (10)

Consequently B with L=4 and C with B0=4 carry NO actual pole.
The necessary pole ranges sharpen to L>=5 for B and B0>=5
for C. The A range L>=3 remains.

Let the three sets of pole-carrying divisors be Aset, Bset, Cset.
The exact necessary ledger is

    1<=b=sum_c beta_c
      <= sum_Aset(2L_c-5)+sum_Bset(L_c-4)
                                      +sum_Cset(B0_c-4).  (11)

There is at least one such divisor because a nonzero rational exact
differential on a complete curve has a pole. Every other finite
divisor, every M-root and every infinity point is retained in the
complete-end argument; none contributes a hidden pole.

Since sum_c delta_c<=D, equation (10) implies the simpler uniform
bound

    1<=b<=floor(D/2)-5.                                   (12)

For A the comparison is equality in its local upper bound. For B
and C, delta_c>=20 after the strengthened ranges, so their bound
floor(delta_c/4)-4 is at most delta_c/2-5. Sum over the nonempty
pole-divisor set and then use the entire leading divisor degree D.

If Aset is EMPTY, the stronger bound is

    1<=b<=floor(D/4)-4.                                   (13)

These inequalities keep leading zeros which do not support poles:
they consume degree D and cannot increase the right side of (11).
No coprimality of b3 and M or restriction on the number, position
or order of its zeros is imposed.

## 4. Infinity local degrees force alpha<=2

At an arbitrary normalized intermediate point over u=infinity,
write u=s^-alpha. The degree-six monic parameter equation has

    h(u)=M(u)^3/b3(u), ord_s h=(D-3m)*alpha.

The full monic-parameter lemma, including every leading cancellation,
gives

    ord_s omega
      >=ceil(alpha*(D+3m-6)/6)-1.                         (14)

This is the PF149 infinity calculation with n=3. Its proof uses the
first nonconstant coefficient of w and does not impose a bound on
lower coefficient degrees or on the response depth.

The point is holomorphic. A nonconstant primitive Q is finite there;
its local multiplicity over its value is ord omega+1. It is bounded
above by the degree of Q on that component, and hence by b. Thus
every infinity point satisfies the new global/local constraint

    ceil(alpha*(D+3m-6)/6)<=b.                             (15)

Combining (12),(15), alpha>=3 would give

    b>=D/2+3*(m-2)/2>D/2-5>=b,

a contradiction. Therefore EVERY intermediate infinity point has

    alpha in {1,2}.                                      (16)

If Aset is empty, alpha>=2 would give b>=D/3+(m-2),
which is strictly larger than D/4-4. Consequently

    Aset empty => alpha=1 at EVERY infinity point.          (17)

In particular the exact integer inequalities

    ceil((D+3m-6)/6)<=b<=floor(D/2)-5
    [or, if Aset empty, b<=floor(D/4)-4]                   (18)

must hold. With m>=2 they imply D>=16 universally, and D>=48
if Aset is empty. For the first assertion, the real inequalities
give D>=15; D=15 would require the integer b to be both at least
three and at most two. The second follows already from D/6<=D/4-4.
Larger m is retained in the sharper exact constraints (18).

More generally (11),(15), using every pole divisor, can exclude
coefficient configurations whose actual available pole budget
does not reach any of their infinity local degrees.

## 5. All original points over t=infinity are unramified

Return to the ORIGINAL quadratic source extension L/E, t^2=u.
At a point in (16), the residue field after geometric base change
is algebraically closed.

If alpha=1, the local quadratic extension is ramified of index
two. Above it t has pole order (2*alpha)/2=1.
If alpha=2, the leading unit has a square root in the complete
local ring. There are two split points and at each t has pole
order alpha/2=1. Thus every original point over t=infinity is
unramified for the source projection t.

The total degree of that projection on the original geometric
generic source curve is six, because P has degree six in x.
Consequently there are EXACTLY six points over t=infinity,
all with pole order one. This total remains valid if generic
components split. On E the possible infinity index partitions are
made only of ones and twos; if Aset is empty they are six ones.

Unramified t projection does not force a degree-one residue field
over C(P). For a precise local guard, the unrelated polynomial
P0=x+(t^2*x)^6 is a submersion on the whole plane. At t=infinity,
put s=1/t and x=s^2*Y. Its exact equation becomes

    Y^6+s^2*Y=p.

All six geometric branches are unramified in s and their coframe
dt/(P0)_x is holomorphic, but the original infinity residue field
is C(lambda), lambda^6=p, of degree six over C(p). This is not an
actual Keller pair or a quadratic-image realization. It prevents
promoting source unramification or holomorphicness alone to
degree-one residue or a forbidden graph image.

This conclusion is about the original source curve and its
projection, not normality of a translated open or replacement of
the full finite target normalization. All original target-boundary
primes, including P-vertical ones, remain in the original S.
No full-normalization first-derivative integrality or termination
of HC4 normal tails follows here.

## 6. Scope, controls and remaining obstruction

The argument applies to arbitrary polynomial coefficients in (2),
every complete quadratic source model (1), all companion degrees,
and all original ramification. It improves the entire ruling-six
necessary ledger rather than treating another selected family.

Equations (11),(18) alone do not contradict every large-D numerical
ledger. PF165 separately closes the B/C alternatives, and PF166 sharpens
the remaining A ledger. A pole faces with intermediate infinity indices
one or two remain unresolved. No polynomial Keller mate in that class
is constructed here.

The checker verifies the strict multiplicity calculation,
all finite ramification-index partitions allowed by ell=3,4
through complete period reductions, the local/global coefficient
comparisons, integer lower-bound controls and the original quadratic
infinity lifting table. The checker passes 230 exact controls. The unbounded proof is
(7)--(18), not a finite coefficient prefix. Replay with

    .venv/bin/python research/jacobian/affine-complements/verify_ruling_six_pole_budget.py

No new external classification is imported. HC4/JC2, entire
ruling-six closure, derivative integrality and normal-tail
termination remain open. Independent whole-proof review,
formal verification and external human review remain unclaimed.
