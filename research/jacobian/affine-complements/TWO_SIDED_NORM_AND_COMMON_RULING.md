# Two-sided norms with arbitrary unramified boundary order

15 November 2026. PF2RULED49 removes the simple-zero hypothesis on the
unramified prime from a substantial all-degree sector of PF2RULED48.
Both complete principal parts are retained. This is a written polynomial
and valuation proof, not a bounded Laurent search. A separate unbounded
family of affine exact-volume controls shows why common source rulings
cannot be excluded from those geometric hypotheses alone.

The full quartic two-prime and higher-degree problems remain open.
Complementary coordinates and the inequalities below remain hypotheses.

## 1. Actual charts and the two characteristic defects

Let normal affine N contain the marked source U=A2_(x,y) and disjoint
omitted affine lines R,D. The source volume dx wedge dy is globally exact,
with orders e-1 and zero on R,D respectively. Assume smooth affine
individual opens U union R and U union D. The SAME source coordinate
pair satisfies y|R=0, v_R(y)=r and x|D=0, v_D(x)=s>=1.
Other omitted primes are allowed and retained in N.

First suppose r>=max(2,e). Use the entire finite-root charts established
in PF2RULED48 and PF2RULED46, without independently shearing the source:

    y=z^r, x=h(z)+z^(e-r)*t,
    h=sum c_j*z^(-j), r-e+1<=j<=q, c_q!=0, c_r=0;
    x=w^s, y=g(w)+w^(1-s)*v,
    g=sum d_j*w^(-j), s+1<=j<=p, d_p!=0.               (1)

The support gcds with r and s are one. At D the initial lower bound s
from PF2RULED46 improves to s+1 because exactness kills d_s. At s=1 this
is the usual whole chart y=g(x)+v with no x^-1 residue. The arbitrary-index
R chart still requires the explicit smoothness and regular-coordinate
hypotheses of PF2RULED48; these are not deduced for every higher-index prime.

For h define

    j_R(a)=max{j:c_j!=0, r does not divide a*j},
    S_R=sum_(a=1)^(r-1)(2q-j_R(a)), Delta_R=r*q-S_R.

Define S_D,Delta_D similarly from s,p and g. Empty sums at s=1 give
S_D=0,Delta_D=p. These are full-support invariants. The sufficient bounds are

    q>=r+1, q>=s, p>=r,
    Delta_R>=e, Delta_D>=s+1.                          (2)

Then N cannot exist with these source charts. In particular the conclusion
is independent of finite Keller field degree and target normal positivity.
For a quartic normalization e=2 and the required smooth component opens
are already available. If gcd(q,r)=gcd(p,s)=1, the defects are exactly q,p.
Thus every such pair with q>=max(r+1,s), p>=max(r,s+1) is excluded, with
ALL lower coefficients allowed by the whole-chart hypotheses. Neither
leading coprimality nor the sufficient defect is required to classify an
individual affine extension; they are only sufficient for this exclusion.

## 2. Cut only the incompatible norm band

Define the genuine full norm polynomials

    P_R=y^q*product_(z^r=y)(x-h(z)),
    P_D=x^p*product_(w^s=x)(y-g(w)).                    (3)

Each has a nonzero constant at its zero base coordinate. P_R has top term
x^r*y^q; its other terms have x degree <=r-1 and y degree <=q-1.
Its x^(r-1)*y^(q-1) coefficient is -r*c_r=0. These assertions follow
directly by multiplying the Laurent roots, or by the quotient-algebra norm.

Put P0=x^(p-r)*P_R. Describe each of its nonleading monomials as

    a_(k,l)*x^(p-k)*y^(q-l), 1<=k<=r, l>=1.

Delete precisely those with l<max(2,k), calling their sum E. The missing
(k,l)=(1,1) coefficient is essential: every deleted term has k>=2 and
1<=l<=k-1. This is a prescribed finite band of the FULL norm polynomial,
not a truncation of the input principal part. Define

    F=P0-E+y^(q-s)*P_D-x^p*y^q.                       (4)

All exponents are nonnegative by (2). The leading monomial is x^p*y^q;
every other monomial is in 0<=i<=p-1,0<=j<=q-1. In the D norm's
nonleading terms, if b is the y-degree loss and a the x-degree loss,

    1<=b<=s, a>=b+1.                                 (5)

Indeed each selected Laurent factor has exponent at least s+1; a
nonzero conjugate-symmetric coefficient has total exponent a*s.
This gives a*s>=b*(s+1), hence a>=b+1. No leading-symbol approximation
or cancellation assumption is involved.

## 3. Both whole chart equations have finite simple roots

At the reference R branch, the derivative of P0 has z order

    B_R=S_R+2q-pq; M_R=B_R+e-r.                       (6)

This is the full conjugate-difference identity of PF2RULED48. A deleted
monomial has order minus M_R equal to

    Delta_R-e+q*(k-2)-r*(l-1)
      >=Delta_R-e+(q-r)*(k-2)>=0.                     (7)

The D norm correction has the same difference, with k=a,l=b. By (5),

    Delta_R-e+q*(a-2)-r*(b-1)
      >=Delta_R-e+(q-r)*(b-1)>=0.                     (8)

All correction dependence on t has strictly larger order, since
q+e-r>=1. The main norm's higher t powers also have larger order:
z^(e-r)*t has strictly larger order than every conjugate difference of h.
Thus the entire substitution, after factoring z^M_R, is

    alpha*t+beta+O(z), alpha!=0.                      (9)

At D, take y^(q-s)*P_D as the reference polynomial. Its y derivative
has w order B_D=S_D+2p-pq, and M_D=B_D+1-s. The remaining R correction
consists precisely of terms with l>=max(2,k). For such a term its order
minus M_D is

    Delta_D-1+p*(l-2)-s*(k-1).                       (10)

For k=1 this is >=Delta_D-1. For k>=2 it is at least

    Delta_D-1-s+(p-s)*(k-2)>=0.                       (11)

Here Delta_D<=p, so (2) implies p>=s+1. As at R, all correction
dependence on v is of higher order because p+1-s>0. The whole D
substitution, after w^M_D, is gamma*v+epsilon+O(w), gamma!=0.
The implicit roots are finite points of the actual root charts.

All branches at infinity are accounted for. The corrections are strictly
above the reference R infinity face: M_R minus the face order is
2q-Delta_R+e-r>=q+e-r>0. At D the analogous difference is
2p-Delta_D+1-s>=p+1-s>0. Thus the r conjugate R roots and s conjugate
D roots exhaust the two faces. Primitivity of each full Laurent support
descends each group to one branch on the source curve. Both coordinates
tending to infinity, or one tending to infinity with the other tending to
a nonzero finite value, are impossible because the top monomial strictly
dominates the entire lower rectangle. Hence there are no further ends.

The projective normalization of any irreducible component of F=0 now maps
nonconstantly into U union R union D, and hence into the original affine N.
Regular functions on a proper integral complex curve are constant, whereas
an affine variety separates distinct image points. This contradiction
proves the theorem. F itself need not be globally regular on N.

## 4. Simple index-two R and arbitrary order on D

There is a stronger separate row when e=2,r=1. Retain the simultaneous
constant in x=y*t+A+h(y), with h=sum_(j=2)^q c_j*y^-j and C(y)=y^q*h.
Use the arbitrary-s D norm P_D from (3). Suppose only

    q>=max(2,s), Delta_D>=s+1.                         (12)

In particular p>=s+1. The polynomial is

    F=x^p*y^q-x^(p-1)*C(y)-A*c_q*x^(p-2)
                         +y^(q-s)*P_D-x^p*y^q.       (13)

The first three terms have R leading chart order -pq+2q+1 and a nonzero
affine-linear t coefficient, as proved in PF2RULED48. The D norm
correction has order at least -pq+3q-1 there, using (5), and its t
dependence is higher. It therefore cannot precede that coefficient
because q>=2. At D the two remaining corrections have deficits (1,l>=2)
and (2,q>=2). Equations (10),(11) apply unchanged. Both implicit roots
are finite and the same two-face argument proves the contradiction.

Thus a simple R zero does NOT require a simple D zero for exclusion.
For coprime p,s, every p>=s+1, q>=max(2,s) is excluded, with arbitrary
lower principal parts and A retained. At s=1 this recovers the entire
simple/simple sector of PF2RULED48. Failure of (2) or (12), including
p<r, q<s or insufficient defects, is not a surviving affine realization.

## 5. Why a common ruling is a real geometric escape

Here is an actual unbounded affine control family, with no finite Keller
map asserted. Choose q>=2 and distinct nonzero a,b. Glue three entire
planes U_(x,y), W_R_(t,y), W_D_(v,y) over y!=0 by

    x=y*t+a*y^-q on R; x=v+b*y^-q on D.               (14)

Every pair of directional translations has a pole. Between R and D they
are (a-b)*y^-q and (b-a)*y^(-q-1). Between U and either added plane this
is immediate from (14). These are the whole line-bundle torsor charts
of Dubouloz Section2.8, with multi-index (1,0,1) in the order U,R,D.
Thus the linear transition exponent is mu_i-mu_j, as in Section2.7.
Its exact iff criterion
therefore makes the entire glued surface affine. It is smooth, its
marked complement is exactly two disjoint affine lines, and the source
volume has divisor R. Both individual opens are affine by the same test.

There is even an explicit global polynomial primitive on the marked plane:

    theta=(y*dx+q*x*dy)/(q-1), dtheta=dx wedge dy.

In the R and D planes it is respectively

    (y^2*dt+(q+1)*y*t*dy)/(q-1),
    (y*dv+q*v*dy)/(q-1),                              (15)

so it extends globally. Affineness, exactness, disjoint boundary lines and
orders 1,0 alone do not exclude the common-coordinate configuration.
In particular they cannot replace the global finite target algebra.
This does not establish any finite Keller realization or claim that every
common-coordinate surface has this form. No uniqueness theorem for the
componentwise coordinate pencils is asserted.

## 6. Reproducibility and a resultant sign pitfall

The [checker](verify_two_sided_norm.py) writes
[the receipt](two_sided_norm_receipt.json). Its primary inputs and exact
assurance boundary are in [the manifest](sources/two-sided-norm/manifest.json).
It verifies both complete substitutions, supports, contact defects and
Newton faces in representative single and multiple-characteristic data,
including arbitrary coefficients in a small symbolic row, all indices
e=1,2,3,4 and nonsimple orders on D. Such finite controls support the
written unbounded proof; they do not prove it.

The norm is computed as the determinant of multiplication in the quotient
algebra. A second check first reduces the resultant's second input modulo
the monic root equation. This avoids the degree-swap sign behavior of the
installed SymPy1.14.0, discussed in its primary
[issue25406](https://github.com/sympy/sympy/issues/25406). For example its
raw resultant of z^3-1 and 10z^7-3-2z^3 is -875; the root-product norm,
direct multiplication determinant and Sage give +875. The older sealed
PF2RULED48 controls did not include this odd/odd case and remain unchanged.
Its mathematical norm definition and signed resultant identity are not
invalidated by this software behavior; its raw resultant helper is not
safe for arbitrary new input signs. Reuse the new determinant helper.

Full compatible-coordinate classification, small-pole and failed-defect
sectors, the finite algebra on common-ruling controls and the unrestricted
quartic/higher-degree problem remain open. Formal verification, independent
full reproof, external review and novelty priority remain false.
