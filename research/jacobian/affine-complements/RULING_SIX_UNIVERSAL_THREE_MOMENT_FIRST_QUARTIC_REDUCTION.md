# Three adapted moments force a square first quartic

PF2RULED187. 17 November 2026. Actual complex polynomial plane Keller maps,
in the complete marked quadratic obstruction class.

At EVERY actual A pole support, the first separated quartic initial must be

    R0(Z)=(Z^2+j0)^2, j0 in C*.

Equivalently, writing the complete depressed quartic as

    G=X^4+2J*X^2+q*X+J^2+Delta,

and t=ord q, ell=ord J, delta=ord Delta, every survivor satisfies

    ell finite, ell<L, 2ell<delta, 3ell<2t.              (1)

Thus pure cubic/linear, dominant-Delta, square/Delta boundary, and every
other NONSQUARE first quartic are impossible for arbitrary slope L,
arbitrary leading degree D, arbitrary source coefficient degrees and
arbitrary actual companion degree. Repeated roots of the first initial
are fully included.

This theorem alone does not exclude the remaining first-square region (1).
The subsequent [square-deformation theorem](RULING_SIX_SQUARE_DEFORMATION_AND_MARKED_QUADRATIC_EXCLUSION.md)
excludes that region and closes the retained marked quadratic ruling-six
model. Global model/import coverage, target integrality, HC4 and JC2
remain open.

## 1. Universal source and companion inputs

Retain the original fields, full normalizations and actual volume form:

    w^2=N+M*z, dP wedge dQ=du wedge dw/M,
    P=H^2+R, H=H3, R=d*z+e*w+f.                         (2)

At an actual A support c, put v=u-c and Y=v^L*w. The prior unrestricted
source and balanced-remainder results give

    L>=3, M(c)!=0, ord_c h=2L, ord_c K=L,
    H=a3*w^3+a2*w^2+a1*w+a0, a3=h/M, a2=K/M,
    v^L*H(v,v^-L*Y) mod v = C*Y*(Y-gamma)^2,
    C,gamma in C*,
    ord_c e=b, ord_c d=L+b, 0<=b<=L-2,
    t=L+b, L<=t<=2L-2.                                  (3)

All coefficients of H are regular at c. In particular H specializes
to a LINEAR polynomial with nonzero derivative. P specializes to a
nonconstant QUADRATIC in w: precisely two roots stay finite, and the
other four escape. No point or prime is removed.

PF176 and the inverse prefix used in PF185 give, for the globally
reduced ACTUAL companion,

    Q-sum_(n=0)^k c_n*P^(n/2)-c_-1*P^(-1/2)=O(w^-5),
    k>=3 odd, c_n,c_-1 in C, c_k!=0.                     (4)

This holds before the degree-thirty restriction in those later proofs.
It uses the full Jacobian and the noncube Kummer equation proved in
PF173. No target change by a half-power or degree cutoff is made.

## 2. Three exact zero identities and the universal bounded center

Define the rational coefficient and three polynomial weights

    a=-a2/(2a3)=-K/(2h),
    B0(w)=1, B1(w)=w-a,
    B2(w)=w*(w-a)^2/a.                                  (5)

Here a has pole order L and is nonzero in C(u).
The inverse-half expansion begins

    P^(-1/2)=a3^-1*w^-3-a2*a3^-2*w^-4+O(w^-5).

Consequently

    [w^-1] B0*P^(-1/2)=[w^-1] B1*P^(-1/2)=0,
    [w^-1] B2*P^(-1/2)
          =a3^-1*(-a2/(a3*a)-2)=0.                     (6)

The cubic degree of B2 also makes the weighted residue of the
O(w^-5) remainder zero. Polynomials Q and even half-powers contribute
zero. Thus, for ALL positive odd n, put

    W_(j,n)=[w^-1] B_j(w)*P^(n/2), j=0,1,2.

The SAME constants in (4) satisfy THREE local identities

    sum_(1<=n<=k,n odd) c_n*W_(j,n)=0, j=0,1,2.          (7)

The adapted center is boundedly close to the separated quartic midpoint
for every slope, not only L=7. Indeed factor the scaled source cubic as

    v^L*H(v,v^-L*Y)
       =C(v)*(Y-r(v))*[(Y-m_H(v))^2+J_H(v)].

The regular linear specialization of H has an analytic finite simple
root, so r(v)=O(v^L). The exact coefficient root sum gives

    a=v^-L*(m_H+r/2).

The whole scaled remainder first deforms the square at order t>=L.
The coprime quadratic/quartic Hensel factorization therefore gives
m_G-m_H=O(v^t), with m_G the quartic midpoint at p=0. With
b_G=v^-L*m_G,

    a-b_G=v^-L*(m_H-m_G+r/2)=O(1).                       (8)

Only a constant nonzero resultant is inverted in this factor recursion.
The center identity uses all source coefficients and does not put a
source function into the target.

## 3. A quartic residue vector can vanish only for a square

Let R0(Z) be ANY monic quartic and let n be a positive odd integer.
Define

    r_(j,n)=[Z^-1] Z^j*R0(Z)^(n/2), j=0,1,2,

with the monic branch at infinity.

If all three coefficients vanished, let S(Z) be the polynomial part.
Then

    R0^(n/2)=S+O(Z^-4).

Apply the differential operator

    L_n=2R0*d/dZ-n*R0'.

It annihilates R0^(n/2), and sends O(Z^-4) to O(Z^-1).
Thus L_n(S) is both a polynomial and O(Z^-1), so L_n(S)=0.
It follows that the derivative of S^2/R0^n is zero. Leading coefficients
give S^2=R0^n. Since n is odd, every root multiplicity of R0 is even:
R0 is a square.

Therefore, for EVERY nonsquare monic quartic and EVERY positive odd n,

    (r_(0,n),r_(1,n),r_(2,n)) != (0,0,0).                (9)

No generic-root or discriminant assumption is needed. Triple/simple,
double/simple/simple and other nonsquare multiplicity patterns are covered
by the same argument.
A monic quartic with exactly two double roots is itself a square and
is correctly left outside (9).

## 4. The full contour comparison at any subcritical first scale

Suppose the first separated quartic scale is X=v^rho*Z, rho<L/2,
and its monic initial R0(Z) is nonsquare. After clearing a finite
Puiseux power, all rescaled coefficients are analytic, and

    P(v,w) ~ v^(4rho-2L)*tau^2*R0(Z),
    w=b_G+v^(rho-L)*Z, tau in C*.

Both root clusters have even size, two and four. The selected square-root
branch is single-valued outside discs around these clusters. Its coefficient
residue splits into their two contour integrals. The other separated factor
is an analytic unit on the whole escape disc, with nonzero constant limit
tau^2. These are the same full contour facts as PF184; no finite cluster is
discarded.

Set lambda=L-2rho>0 and epsilon=v^(rho-L). By (8),

    B_j(b_G+epsilon*Z)=epsilon^j*(Z^j+o(1)),
    j=0,1,2.                                           (10)

For B2 this uses w/a tending to one and a-b_G=O(1).
Define the scaled vector

    V_n=(W_(0,n), v^(L-rho)*W_(1,n),
                       v^(2(L-rho))*W_(2,n)).

All three zero identities (7) still hold componentwise after this
scaling. Its ESCAPE contribution has the common asymptotic

    V_n=tau^n*v^(-p_n)*(r_(0,n),r_(1,n),r_(2,n))
                                          +o(v^(-p_n)),
    p_n=lambda*n+L-rho.                                 (11)

On a FIXED finite circle, the three unscaled weights are bounded by

    B0=O(1), B1=O(v^-L), B2=O(v^-L).

Hence the FINITE contributions to V_n have respective upper pole bounds

    0, rho, 2rho-L=-lambda.                              (12)

The top common escape order is always strictly larger than ALL these
finite bounds, since

    p_n-rho=lambda*(n+1)>0.                              (13)

By (9), its leading vector is nonzero. The exact whole leading vector
therefore equals (11). Since p_n grows strictly with n, the top
character c_k*V_k has a strictly larger pole than every smaller odd
character. It cannot cancel in (7). This excludes EVERY nonsquare first
quartic with rho<L/2.

These comparisons take place on one common Puiseux cover. They do not
identify ramification indices with residue-field degrees, or assume
that the four geometric sheets are distinct normalized points.

## 5. The parameter cap excludes the critical first scale

Over C(p), write the actual depressed quartic as

    G=X^4+2J*X^2+q*X+D0, D0=J^2+Delta.

PF182 proves that BOTH D0 and Delta have a nonzero affine-p coefficient
at order v^(2L), and earlier coefficients are p-independent. In particular

    ord_v D0<=2L, ord_v Delta<=2L.

Its first Newton scale is therefore

    rho=min(ord J/2, ord q/3, ord D0/4)<=L/2.            (14)

If rho<L/2, every coefficient of its first quartic occurs strictly
before order 2L and is independent of p. Thus the first initial used
in Section4 at p=0 is exactly the generic first initial.

If rho=L/2, the whole first initial is a monic quartic whose ONLY
p dependence is a nonzero affine term in its constant coefficient.
Its other coefficients are constants. Over C(p), that quartic is
separable, with no zero root. For example a repeated root would also
be a root of its p-independent derivative, forcing p constant.
On every geometric branch,

    ord_v G_X=3L/2.

The whole F0_Y has the same order because the separated quadratic
is a unit near the escape cluster. Therefore

    ord_s omega=alpha*(L+1-3L/2)-1
               =-alpha*(L/2-1)-1.

Since L>=3, EVERY such branch is an actual Q pole. The strict
affine-response bound of PF164 for an A pole requires
ord_v G_X<3L/2, a contradiction. Equivalently, the first response is
a degree-four equation in the parameter and is not affine.

Thus rho=L/2 is impossible. The first scale is strictly smaller,
and Section4 forces its first quartic to be a square.

## 6. Exact necessary first-square region and remaining scope

The quartic is depressed and monic. A square first initial must be
(Z^2+j0)^2; its linear square-root coefficient is zero because the
cubic coefficient is zero. Also j0!=0: the first scale (14) is attained
by at least one coefficient other than the leading one, so its first
initial cannot be Z^4 alone.

It follows that J has finite order ell=2rho<L, its linear term q
is strictly later than weight 3rho, and the constant coefficient has
exactly the leading term of J^2. Hence

    3ell<2t, delta>2ell,

which is exactly (1). Conversely these strict order conditions do
give a square first initial. No converse about an actual companion
is asserted.

The theorem excludes all competing first initials, including the
entire cubic/linear region, the dominant-Delta region, the even
square/Delta boundary, and cubic/constant or triple/simple boundaries.
The deeper two-square-center deformation and its first parameter
responses are handled by the subsequent square-deformation theorem.
The present note proves the first-initial reduction only; no target
integrality or JC2 closure is claimed here.

## 7. Exact evidence and replay

The checker verifies the inverse-prefix cancellation, the exact root-sum
center, the quartic differential-operator identity, nonzero moment vectors
on square and nonsquare controls including repeated roots, the three full
escape/finite scalings, whole sextic moment controls with bounded center
mismatch, and the parameter-cap edge. Its finite controls test the formulas;
the differential and contour arguments above prove the unrestricted theorem.

Replay:

    .venv/bin/python research/jacobian/affine-complements/verify_ruling_six_universal_first_quartic.py

Essential proof inputs are PF164, PF171, PF173, PF176, PF178, PF179,
PF182, PF184, PF185 and PF186. All original source, normalization and
actual-Keller hypotheses are retained. Internal exact replay is not formal
verification or external human review.
