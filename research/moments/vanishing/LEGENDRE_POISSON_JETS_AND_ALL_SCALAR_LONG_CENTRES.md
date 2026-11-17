# Legendre--Poisson jets and all scalar Long centres

16 November 2026. Internal written proof, using the retained proper-volume
unit-period theorem. Independent replay, formal verification and external
review are not claimed. Both COMPLETE exceptional MZ cases remain **OPEN**.

Every scalar shift of the primitive opposite-weight Long architecture is
excluded from eventual purity on both exceptional quotients, unless one
of its two weight functions is zero. In that degenerate case the unique
eventually pure centre is 1, and the resulting one-weight phase is safe.
This is an all-degree classification of this architecture, not of the
whole exceptional algebra. General nonlinear composites remain unresolved.

See the [checker](verify_haar_legendre.py),
[receipt](haar_legendre_verification.json) and
[source audit](haar_legendre_source_audit.json).

## 1. Opposite weights give Legendre moments

Let A be a commutative complex algebra with a rational Cstar action and
normalized invariant linear functional h. For a nonzero integer a, take
u in A_a, v in A_-a and put

    P=(1+u)(1+v+uv), tau=1+2uv,
    ell(R)=h(R(tau)).                                  (1)

No unit or domain assumption is needed here. Weight projection gives

    h(P^m)=sum_j binom(m,j)^2 (uv)^j(1+uv)^(m-j)
            evaluated by h
          =ell(L_m),                                  (2)

where L_m is the usual degree-m Legendre polynomial. One can verify the
polynomial identity by its Rodrigues formula

    L_m(t)=1/(2^m m!) D_t^m[(t^2-1)^m].                (3)

Expansion and the binomial theorem give (2); equivalently expand the
constant term of (1+u)^m(1+uv+v)^m with uv held fixed.
Thus L_0,L_1,... are a basis of C[t].

Fix ANY complex c. If P-c is eventually pure, write
mu_j=h((P-c)^j), so mu_0=1 and mu_j=0 for j>=M. For EVERY m,

    ell(L_m)=sum_(j<M) mu_j binom(m,j)c^(m-j)
            = D_mu c^m,
    D_mu=sum_(j<M) mu_j/j! * partial_c^j.              (4)

The first expression uses only j<=m and is polynomial also at c=0.
Define ell_c uniquely by ell_c(L_m)=c^m. Equation (4) says, on ALL
polynomials, ell=D_mu ell_c at the chosen c. This is an exact consequence
of the entire vanishing tail, not of finitely many computed moments.

## 2. The Cauchy function of this finite Poisson jet

For |c|<1 the elementary Legendre Poisson density is

    q_c(t)=(1-c^2)/(2(1-2ct+c^2)^(3/2)), -1<=t<=1,
    integral q_c(t)L_m(t)dt=c^m.                     (5)

For clarity, the generating identity from (3) is
sum_m L_m(t)c^m=(1-2ct+c^2)^(-1/2). Applying
(2c partial_c+1)/2 gives q_c=sum_m (2m+1)L_m c^m/2.
Integration by parts in (3) gives orthogonality to lower degrees and
integral L_m^2=2/(2m+1), proving (5). The last norm follows by applying
the same integration by parts to the leading coefficient
2^(-m) binom(2m,m) and the beta integral of (1-t^2)^m.

Let A(s,c)=ell_c(1/(s-t)), interpreted first as its convergent expansion
at infinity. Put

    B^2=1+c^2-2cs,
    R=((1+c-B)(1-c+B))/((1+c+B)(1-c-B)), L=log R.

Direct integration of (5), with z=sqrt(1-2ct+c^2), gives

    A(s,c)=-2c/B^2 + (1-c^2)/(2B^3) L.              (6)

Indeed q_c dt=-(1-c^2)dz/(2c z^2) and
s-t=(z^2-B^2)/(2c). Integrate from z=1-c to 1+c and use partial
fractions. The branch in (6) is the one agreeing with the Cauchy germ.
The logarithm has the exact derivatives

    partial_c L=-2/B,
    partial_s L=2(cs-1)/((s^2-1)B).                  (7)

These identities follow by differentiation of R and B; no branch
constant enters them. At c=0 the removable expression is
L=log((s+1)/(s-1)), B=1, A=L/2.

There is no restriction to small c in the resulting finite-jet identity.
Each ell_c(t^k) is a polynomial in c, since it is a finite Legendre
expansion. Their Cauchy series converges near infinity uniformly on
compact parameter sets: the coefficients of the finite Legendre
expansion of t^k, or the three-term recurrence, give an exponential
bound there. Continue (6) as this germ. For every c!=0 and generic s,
including c=1 and c=-1, its square root and logarithm admit local
analytic choices in c. Hence all finite c derivatives agree with the
continued polynomial moments. The value c=0 may instead be handled
directly by the polynomial density

    q(t)=sum_(j<M) mu_j (2j+1)L_j(t)/2.              (8)

For c!=0, apply D_mu to (6). Each derivative of B is odd under
B -> -B; each positive c derivative of L is also odd, by (7).
More precisely differentiation preserves the parity of a rational
function of B. Therefore the logarithm coefficient is odd and the
remaining terms are even:

    A_tau(s)=C(s)L(s,c)+E(s),
    C=D_mu[(1-c^2)/(2B^3)], E in C(s).              (9)

Only powers of B^2 occur in the denominators of E. In particular its
only possible finite pole is

    d=(1+c^2)/(2c).                                 (10)

This claim includes c=1,-1: the coefficient 1-c^2 may vanish, but its
derivatives must be retained. Terms where derivatives hit L contain an
even product of odd functions of B, and are rational with the stated
denominators. At c=0, division in (8) instead gives

    A_tau(s)=q(s)/2 * log((s+1)/(s-1))+E(s),
    E polynomial, q nonzero since integral q=1.     (11)

## 3. Proper original-volume height models forbid both alternatives

Assume the height moment functional ell in (1) is represented by the
normalized compact period of a NONCONSTANT regular height on a smooth
affine quadric Q_N, N>=2, with a projective smooth extension carrying
the ORIGINAL holomorphic volume. These are precisely the substantive
hypotheses of the retained
[proper unit-period theorem](PROPER_VOLUME_UNIT_PERIODS_AND_COMPLEX_HOPF_CLOSURE.md).
It supplies a local meromorphic differential unit retraction for the
actual module generated by A_tau and literal 1, at EVERY puncture.
The retraction persists after a finite cover. It is used before any
quotient by polynomial functions and involves the ACTUAL period relations.

If C in (9) is nonzero, pass to the quadratic cover adjoining B and
divide A_tau-E by C. The actual module then contains L. At least one
of the two endpoints is unramified on this cover and has a nonzero
residue for partial_s L:

    at s=1, if c!=1: B=+/- (c-1), residue=+/-1;
    at s=-1, if c!=-1: B=+/- (c+1), residue=+/-1.   (12)

Every c has at least one available endpoint; at c=1 use -1, and at
c=-1 use 1. A unit retraction would send L to a meromorphic primitive
of its rational derivative in the local cover field. This contradicts
the nonzero residue in (12). Zeros of C at the endpoint cause no
problem: meromorphic division is allowed. Formula (11) gives the
same contradiction at c=0. In particular endpoint cancellations of
the apparent logarithm in A_tau itself do not suffice.

If C=0 and c!=0, (9) is rational with at most the ONE finite pole d
in (10). As the actual Cauchy germ it decays at infinity with leading
coefficient 1/s. Hence

    A_tau(s+d)=sum_(k=1)^K b_k/s^k.

Its expansion says h((tau-d)^m)=0 for all sufficiently large m.
But the retained
[projective original-volume period theorem](SEMISTABLE_LOG_PERIODS_AND_SINGLE_ISOTYPIC_SAFETY.md)
excludes eventual purity of every nonconstant polynomial in this height.
Thus the rational alternative is impossible as well. At c=0 the
nonzero polynomial q in (11) already excludes this alternative.

We have proved the general implication:

    a nonconstant opposite-product height with such a proper
    ORIGINAL-volume moment model excludes eventual purity of
    P-c for EVERY c in C.                            (13)

The rational alternative is essential. For example at c=1, ell_1 is
evaluation at t=1, its Cauchy function is 1/(s-1), and its logarithm
coefficient vanishes. It would be incorrect to claim C never vanishes.
No proper-volume model is assumed for P or P-c themselves.

## 4. Complete affine-centre classification on both exceptional quotients

For Gamma=2O or 2I let A_Gamma=C[SL2]^Gamma with its original Haar
functional and any maximal left torus. Every neutral height tau lies
in C[Q_2]^barGamma. The retained
[torus volume construction](VALUATIVE_BOUNDARY_VOLUMES_AND_TORUS_INVARIANT_CLOSURE.md),
Sections 2--3, gives the proper ORIGINAL-volume moment model required
in Section 3 for EVERY nonconstant such tau, with all complex
coefficients and lower terms. The quadric here is the auxiliary Q_2;
its normalized compact period is exactly the neutral Haar functional.
We do not assert that a general mixed-weight phase on SL2/Gamma has
a proper-volume model.

The algebra A_Gamma is a domain with only constant units. One way to
check the unit assertion is to pull a unit back to SL2. The degree
filtration has associated graded C[x,y,z,w]/(xw-yz), a graded domain;
therefore the degree of a product is the sum of degrees. A product
equal to 1 has both degrees zero. Thus for nonzero u and v of opposite
nonzero weights, uv cannot be constant: a nonzero constant would make
u and v units. Consequently tau=1+2uv is nonconstant and (13) applies.

If uv=0 the domain property gives u=0 or v=0. Then P is 1 plus a
single weight function (possibly zero). For every m, weight projection
gives h((P-c)^m)=(1-c)^m. This vanishes eventually exactly when c=1.
For that centre P-1 is a one-weight phase or zero, safe against EVERY
fixed original marker: its growing weight eventually misses the
finitely many opposite weights of that marker. We obtain the exact law

    P-c eventually pure  iff  uv=0 and c=1,
    and every phase in this degenerate pure locus is SAFE.       (14)

This is valid in every degree on BOTH exceptional quotients and under
the retained algebraic automorphism transport of Haar and the torus.
The domain and unit argument in this application is over C; (14) is
not promoted to arbitrary coefficient rings.

## 5. Every torus-equivariant Bezout Long pullback and every scalar centre

Use the retained exceptional Bezout data, degrees (n,m)=(8,12) or
(12,20), and put L=n+m. For ANY regular b of left weight L,

    Y_b=V+G_0 b, wt(Y_b)=n, wt(F_0)=-n,
    g_b=(1+Y_b)(1+F_0+Y_b F_0),
    tau_b=1+2F_0Y_b=U_bG_0+Y_bF_0.

Here U_bG_0-Y_bF_0=1. Neither F_0 nor Y_b is zero. If Y_b=0,
the determinant identity would make G_0 a unit, which it is not.
Apply (14) with u=Y_b and v=F_0:

    for EVERY b of weight n+m and EVERY c in C,
    g_b-c is NOT EVEN EVENTUALLY PURE, in ANY degree.             (15)

This phase is the primitive Long pullback in the
[cohomological pullback theorem](COHOMOLOGICAL_PULLBACK_AND_A_UNIFORM_ORIGINAL_MARKER.md).
It differs from the Hopf phase closed in the proper-unit proof. Both
equivariant architectures are now excluded, but for different moment laws.
The all-b nonzero-component and fixed-marker results remain valid.

No implication from eventual purity of H(g_b) to that of g_b-c has
been proved. Nonlinear H, other homogeneous b weights, genuinely mixed
b and the ENTIRE exceptional pure locus remain separate. The FULL
target is still purity forcing every corrected actual local torus
valuation zero, or an actual invariant pure failure. Known one-weight
safe phases are expressly allowed by (14); no blanket proper-volume
assumption would allow that. Review remains DEFERRED, not a next gate.
