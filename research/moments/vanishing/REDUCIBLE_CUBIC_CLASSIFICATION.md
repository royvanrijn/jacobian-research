# Complete homogeneous reducible ternary cubic moment classification

Claim `G3REDUCIBLE20261114`, over C in characteristic zero. Let x have three
independent standard real Gaussian coordinates, with complex polynomial
coefficients. For a nonzero homogeneous reducible cubic P, the following
are equivalent:

1. Its first ten positive pure Gaussian moments vanish.
2. Every positive pure Gaussian moment vanishes.
3. It has a factorization P=L*Q2, with L=l^T*x, Q2=x^T*A*x, A symmetric,
   such that l^T*l=0 and A*l is proportional to l.

For every such pure-null P and every multiplier R of total degree D,

    E(R*P^m)=0 for m>D.                                    (1)

Both the uniform pure gate10 and mixed cutoff D are sharp in this entire
class. The complete tangent line–smooth-conic subclass has the sharper
optimal pure gate8. A transverse line–smooth-conic cubic always has a
nonzero pure moment by10, and10 is sharp on that subclass too.

This closes all homogeneous reducible ternary cubics for GMC. Together
with the existing smooth-leading and
[cuspidal exclusions](CUSPIDAL_CUBIC_SHARP_THRESHOLD.md) and the subsequent
[nodal Newton exclusion](NODAL_NEWTON_EXCLUSION.md), every homogeneous
ternary cubic is now GMC-safe. Lower-degree perturbations and the separate
coefficient-scheme questions remain open.

The subsequent [affine triangle theorem](AFFINE_TRIANGLE_EXCLUSION.md)
closes every lower-degree perturbation of three linearly independent
leading factors. Concurrent leading lines, line-smooth-conic curves and
repeated lines remain outside that new theorem.

The subsequent [concurrent-line theorem](CONCURRENT_LINE_EXCLUSION.md)
also covers every affine perturbation of three distinct concurrent lines.
Only line-smooth-conic and repeated-line leading curves remain possible
cubic GMC failure classes.

The subsequent [affine line-conic theorem and squarefree-leading
corollary](LINE_CONIC_AFFINE_EXCLUSION.md) close every reduced leading
ternary cubic. Only repeated leading lines remain as possible cubic
GMC failures; the earlier arguments below retain their own scopes.

## 1. A small universal invariant calculation

First assume Q2 is nondegenerate, so det A!=0. Set

    D(t)=det(I-2tA)=1+d1*t+d2*t^2+d3*t^3,
    H(t)=l^T*adj(I-2tA)*l=h0+h1*t+h2*t^2.

Then d3=-8det A!=0. The joint Gaussian generating identity is

    E exp(t*Q2+s*L)=D(t)^(-1/2)*exp(s^2*H(t)/(2D(t))).        (2)

This is the homogeneous, uncentered quadratic identity. There is no
trace-centering factor. It follows from Gaussian quadratic completion
or coefficientwise Wick pairing, and is valid as a formal identity for
complex A and l. Extracting s^(2k)*t^(2k) gives

    M_(2k)=E(P^(2k))=(2k)!^2/(2^k*k!)*C_k,
    C_k=[t^(2k)] H(t)^k*D(t)^(-k-1/2).                      (3)

Odd moments vanish by total-degree parity. In particular

    C_1=h2-(3/2)d1*h1+((15/8)d1^2-(3/2)d2)*h0.              (4)

The exact finite assertion used below is

    d3!=0 and C_1=...=C_5=0  iff  h0=h1=h2=0.               (5)

It holds for arbitrary six complex coefficients of D,H, without a
matrix-realizability assumption. Thus it also describes the radical
of these five coefficient equations after inverting d3, rather than
asserting that their nonreduced ideals are equal.

Here is an exhaustive proof of (5). Independently rescale L by r and
Q2 by s; this changes P by a nonzero scalar. On the invariants,

    d_i -> s^i*d_i,   h_j -> r^2*s^j*h_j.

The C_k all scale by r^(2k)*s^(2k). Nonzero complex rescalings therefore
allow the following complete branch list. Equation (4) eliminates h2
with coefficient1 in each normalized branch.

| Branch | Normalization / remaining identity | Contradiction when d3!=0 |
|---|---|---|
| h0*d1!=0 | h0=d1=1 | The ideal(C2,C3,C4,C5) contains1. |
| h0!=0, d1=0, h1!=0 | h0=h1=1 | The ideal(C2,C3,C4) contains d3. |
| h0!=0, d1=h1=0, d2!=0 | h0=d2=1 | C2=-7/8. |
| h0!=0, d1=h1=d2=0 | h0=1 and h2=0 | C3=63*d3^2/8. |
| h0=0, h1!=0 | h1=1 | The ideal(C2,C3,C4) contains d3^2. |
| h0=h1=0 | C1=h2 | H=0, the permitted branch. |

The three displayed ideal memberships have explicit rational Bezout
coefficients in the retained artifact. The checker clears their common
denominators and verifies the full integer polynomial identities. In
particular, the generic unit is a checked characteristic-zero identity;
there is no modular lifting or finite-prefix extrapolation in (5).
The converse follows immediately because H=0 kills every positive C_k.

## 2. Intrinsic classification and the mixed bound

Since

    H(t)/D(t)=l^T*(I-2tA)^(-1)*l,

H=0 implies l^T*A^j*l=0 for every j>=0. The cyclic subspace
K=span{l,Al,A^2*l,...} is totally isotropic: inner products of any two
generators are of this form because A is symmetric. A totally isotropic
subspace of a nondegenerate complex three-space has dimension at most1,
as K is contained in K-perp and their dimensions sum to3. Thus

    l^T*l=0,   A*l=lambda*l.

Conversely these equations force H/D=0. Since A is invertible here,
lambda!=0. The point [l] lies on the conic Q2=0 and its tangent line is
L=0, because grad Q2(l)=2lambda*l. Therefore every pure-null smooth-conic
product is tangent; no transverse product can be pure-null.

In Gaussian Witt coordinates Z,W,T with E(ZW)=E(T^2)=1 and all other
basic pairings zero, put L=Z up to a scalar. The eigenvector condition
removes the W^2 and WT terms of Q2, leaving

    P=Z*(alpha*ZW+beta*T^2+gamma*ZT+delta*Z^2).              (6)

For a smooth conic alpha*beta!=0. Give Z charge1, W charge-1 and T
charge0. Every monomial of P has charge at least1, while a multiplier
of total degree D has charge at least-D. Wick pairing can only retain
charge0. This proves (1) for every multiplier, without a sampled tail.

If the quadratic factor is degenerate, it splits over C, and the
[completed centered three-factor theorem](NODAL_CUBIC_REDUCTION.md)
applies: its pure gate is8 and its mixed cutoff is D. It also supplies
the factorization in item3. In its isotropic-row case choose that factor
L; the other two factors are orthogonal to l, so A*l=0. In its
isotropic-pair case the two factors are proportional in dimension three;
choose either as L, and again A*l is proportional to l. Conversely any
factorization in item3 still has form(6), allowing degenerate coefficients,
and the same charge argument proves pure and mixed vanishing. This
completes the stated equivalence for every reducible homogeneous cubic.

The uniform mixed bound is sharp for every D>=1. Take

    P=Z*(ZW+T^2),    R=W^D.

Its conic is smooth and tangent, P is pure-null, and

    E(R*P^D)=sum_(j=0)^D binomial(D,j)*(D+j)!*(2(D-j)-1)!! >0. (7)

Here(-1)!!=1. Every summand is positive, so the boundary m=D cannot be
removed from a uniform statement.

## 3. The optimal tangent gate is eight

A line tangent to a smooth conic has an elementary linear normal form

    F=X*(XZ+Y^2).

Indeed, put the line at X=0 and its restricted conic equation at Y^2.
The full conic is Y^2+X(A*X+B*Y+C*Z), with C!=0 by nondegeneracy;
replace Z by A*X+B*Y+C*Z. Retain the full invertible Gaussian covariance
G=[[a,d,e],[d,b,f],[e,f,c]] after this change.

The automorphism Y->Y+sX, Z->Z-2sY-s^2X preserves F. It sends d to
d+sa and b to b+2sd+s^2a. Together with the diagonal transformations
(X,Y,Z)->(lambda X,mu Y,mu^2/lambda Z), which scale F, this gives:

| Branch | Normalization | Finite nonvanishing gate |
|---|---|---:|
| a!=0, transformed b!=0 | a=b=1, d=0 | 8 |
| a!=0, transformed b=0 | a=1, b=d=0 | 6 |
| a=0, d!=0 | a=b=0, d=1 | 6 |
| a=d=0 | No further normalization needed | All positive moments vanish |

In each of the first three charts M2 eliminates one variable by a nonzero
constant coefficient. Exact two-variable ideal memberships force det G=0
if moments through the displayed gate vanish. They are checked by native
rational computation and exact symbolic division. In the fourth branch,
each X in F^n must pair with Z, but F^n has strictly more X than Z; hence
all moments vanish. Invertibility there forces b*e!=0. The Gaussian dual
of the line is then exactly the tangency point, in agreement with item3.

For sharpness choose a positive real root e of

    p(e)=3381e^7+8358e^6+8808e^5+5165e^4+1730e^3+263e^2-11e-6.

Such a root exists because p(0)<0<p(1). Put

    f^2=(28e^4+21e^3+9e^2-e-1)/(7e+3),
    c=-4e^2-2e-1,
    G=[[1,0,e],[0,1,f],[e,f,c]].

Its determinant is

    -2*(14e^4+28e^3+19e^2+6e+1)/(7e+3),

which is nonzero at the selected positive root. The first six moments
vanish; parity supplies the seventh zero. Exact reduction modulo p
checks these equations, while the reduced eighth-moment numerator is
coprime to p and is nonzero. This covariance is realized explicitly by
X=xi, Y=eta, Z=e*xi+f*eta+sqrt(det G)*tau for independent real Gaussians
xi,eta,tau. Thus the tangent gate8 is sharp.

## 4. The transverse gate ten is also sharp

The upper bound follows already from (5), since a transverse line cannot
be the tangent at[l]. For sharpness use the invariant chart h0=d1=1,
with h2=3*(4d2+4h1-5)/8 from C1=0. The three equations C2=C3=C4=0
give a finite rational algebra of length28. Exact saturation retains
that length after inverting

    d3*disc(1+t+d2*t^2+d3*t^3)*(4d2+4h1-5).

The ideal obtained by additionally setting this product to zero is the
unit ideal. The checker independently counts the28 standard monomials.
Thus this is a nonempty algebraic sharpness certificate on the required
open set;28 distinct points are not asserted.

Every point in this algebra has a Gaussian realization. Factor
D(t)=product_i(1-2lambda_i*t), with distinct nonzero lambda_i, and write

    H(t)/D(t)=sum_i r_i/(1-2lambda_i*t).

Choose A=diag(lambda_i) and l_i^2=r_i. Partial fractions give exactly
the specified D,H; sum r_i=h0=1, so L is nonisotropic. Also

    h2=4det(A)*l^T*A^(-1)*l !=0,

which makes the line transverse to the smooth conic. Its moments through8
are zero, and parity gives the ninth zero. Its tenth moment is nonzero
because the ideal(C2,C3,C4,C5) in this chart is the unit ideal. This proves
sharpness10 for transverse line–conic cubics and hence for all reducible
homogeneous ternary cubics.

## Evidence and remaining scope

The [checker](verify_reducible_cubic_moments.py),
[certificate](reducible_cubic_moments_certificate.json) and
[receipt](reducible_cubic_moments_verification.json) bind the exact inputs.
The checker derives C_k by a separate multinomial formula, verifies all
three rational Bezout identities as integer polynomial identities,
compares(3) with the full six-covariance Wick moments of both conic types
through10, verifies the tangent charts and algebraic witness, and checks
the length28 sharpness algebra and partial-fraction construction.

```sh
.venv/bin/python research/moments/vanishing/verify_reducible_cubic_moments.py
```

The direct three-variable coupled transverse covariance open chart
reached its30-second limit before producing a basis. Its separate
denominator-zero chart completed, but it is not the proof of the full
class. The invariant calculation closes that class by a different exact
argument. The initial sharpness input had a rejected polynomial expression;
the corrected version uses explicit multiplication and distinct names.
Both inputs are retained. No timeout or invalid input is promoted.

The [source audit](reducible_cubic_source_audit.json) separates the
uncentered joint quadratic identity from the earlier centered-quadratic
formula, records the targeted nonexhaustive literature search, and states
that no new external theorem was imported. All arguments and calculations
are internal; independent replay, external review and formal verification
are not asserted.

The subsequent [nodal Newton theorem](NODAL_NEWTON_EXCLUSION.md) excludes
the entire homogeneous irreducible nodal class and completes all
homogeneous ternary cubics. General nonhomogeneous nodal,
cuspidal-leading and reducible-leading polynomials remain open. No
arbitrary-cubic upper prefix is supplied; its lower bound remains10.
The N=3,4 coefficient schemes, exact generic quaternary null torsion and
lower-power GVC remain separate obligations.
