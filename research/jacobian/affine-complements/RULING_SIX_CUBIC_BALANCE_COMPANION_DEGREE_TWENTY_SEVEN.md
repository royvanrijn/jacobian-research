# Cubic balance forces companion ruling degree twenty-seven

PF2RULED185. 17 November 2026. Actual complex polynomial plane Keller maps,
in the complete marked quadratic obstruction class.

In the remaining degree-D=30 class, EVERY globally reduced actual
companion has ruling degree exactly TWENTY-SEVEN:

    deg_w Q=27.                                         (1)

All source-coefficient degrees remain unrestricted. This is a necessary
degree reduction, not an assumed cutoff or a construction of a companion.
The subsequent [adapted cubic moment theorem](RULING_SIX_ADAPTED_CUBIC_MOMENT_AND_DEGREE30_EXCLUSION.md)
uses the unused weighted prefix equation to exclude this last class
and hence all of D=30. The later
[square-deformation theorem](RULING_SIX_SQUARE_DEFORMATION_AND_MARKED_QUADRATIC_EXCLUSION.md)
closes every leading degree in the retained marked quadratic ruling-six
model. Global model/import coverage, higher ruling degrees, trace
perfection and HC4/JC2 remain open.

The proof reduces the unbounded set of possible top characters to
two indices using the finite cluster AND the first three positive
centered moments. It leaves indices k=3 or9. An exact polynomial-part
formula then excludes k=3 by its actual pole degree.

## 1. Original source and the complete residual class

Keep the original fields, full normalizations and actual volume form
of PF171--PF184. PF184 proves that D=30 forces

    P=H^2+R, H=H3,
    R=d*z+e*w+f, w^2=N+M*z,
    M=constant*u^2, h=constant*(u-c)^14, c!=0,
    L=7, b_field=[C(u,w):C(P,Q)]=8,
    ord_c d=10, ord_c e=3.                              (2)

Every Q pole belongs to the single A cluster, and its FOUR
geometric sheets are all pole-carrying. The depressed quartic is

    G=X^4+2J*X^2+q*X+J^2+Delta,
    ord q=10, ord J>=7, ord Delta=14.                   (3)

The parameter cap and the affine coefficient at order fourteen
are retained over C(p). No point or prime is deleted.

At the cubic/linear scale X=v^(10/3)*Z, v=u-c, the leading
quartic is

    R0(Z)=Z^4+q0*Z, q0 in C*.                            (4)

The three nonzero roots and the remaining inner root have
X order at least 10/3. The source polynomial has leading
factor v^(-2/3)*tau^2*R0(Z), tau in C*, and

    dw=v^(-11/3)*dZ.

Here a local p-independent center a(v) is included in
w=a(v)+v^(-11/3)*Z. It has pole order seven and integer
Laurent exponents. It need not be a target function.

## 2. One more constant coefficient in the inverse expansion

PF165 gives a globally reduced companion degree 3k, k>=3 odd.
Use PF176's formal inverse variable U=P^(1/6). The actual
Jacobian implies that its coefficients q_j are complex constants
for j>-5, and the cubic Kummer symmetry makes q_j zero
unless three divides j in that range.

The q_-3 coefficient is consequently a COMPLEX constant,
while q_-4=0. In the original w Laurent variable this gives

    Q-sum_(n=0)^k c_n*P^(n/2)
          -c_-1*P^(-1/2)=O(w^-5),
    c_n,c_-1 in C, c_k!=0.                              (5)

This is an exact formal expansion of the ACTUAL companion.
No odd half-power is applied as a target change.

For j=0,1,2 and positive odd n define coefficient residues

    T_(j,n)=[w^-1](w-a(v))^j*P^(n/2).                   (6)

The polynomial Q and even half-powers have zero residues.
The two identities of PF184 and the next coefficient of (5) give

    sum_(n<=k,n odd) c_n*T_(0,n)=0,
    sum_(n<=k,n odd) c_n*T_(1,n)=0,

    sum_(n<=k,n odd) c_n*T_(2,n)=-c_-1*M/h.              (7)

The last right side is computed ALGEBRAICALLY, not by discarding
a finite contour. Indeed P^(-1/2) starts with (M/h)*w^-3.
Its residues after multiplication by (w-a)^j are zero for
j=0,1, and exactly M/h for j=2. The O(w^-5) remainder
cannot contribute to any of the three residues.

Since ord_c h=14 and M(c)!=0, the right side of the last
identity has pole order at most FOURTEEN. It may be zero.

## 3. Bounds and exact coefficients for every positive odd character

Use PF184's two-cluster contour decomposition. The finite cluster
has two bounded roots and P converges to a nonconstant quadratic
on a fixed surrounding circle. For each fixed positive odd n,

    finite contribution to T_(j,n)=O(v^(-7j)).           (8)

This bound is essential for j=2: it can dominate a small-index
escape term. It is retained throughout the proof.

At the four-root escape cluster, (4) gives the general bound

    escape contribution to T_(j,n)
          =O(v^(-(n+11(j+1))/3)).                       (9)

The coefficient at that first possible order is

    tau^n*[Z^-1]Z^j*(Z^4+q0*Z)^(n/2).                  (10)

It is NONZERO whenever

    d_(j,n)=(2n+j+1)/3 is an integer,

because then its exact value is

    tau^n*binom(n/2,d_(j,n))*q0^d_(j,n).                 (11)

To verify it, expand
Z^(2n+j)*(1+q0*Z^-3)^(n/2) at infinity.
For positive odd n the binomial top is a noninteger half-integer,
so the displayed integer coefficient never vanishes.

All actual T_(j,n) have INTEGER Laurent exponents in v:
the uncentered residues are rational functions of u and the
center a(v) has integer exponents. Thus (8),(9) imply the
whole-moment pole bound

    pole_order T_(j,n)
       <=max(7j, floor((n+11(j+1))/3)).                  (12)

When (11) is nonzero AND its escape order exceeds 7j,
that escape term determines the EXACT whole leading pole.
Neither the finite contribution nor any later rescaled term
can cancel it.

## 4. The top index is three or nine

Let k>=3 be the largest odd index in (5), with c_k!=0.

If k is congruent to ONE modulo six, then k>=7.
Choose j=0. The exact top pole order is

    (k+11)/3>0.

Every lower positive odd n is at most k-2. Its bound (12)
is strictly smaller, because its escape bound is at least
2/3 lower, and its finite bound is zero. The top term
cannot cancel in the first identity (7). This excludes
ALL such k, with no upper cutoff.

If k is congruent to FIVE modulo six, choose j=1.
Now the exact top pole is

    (k+22)/3>=9>7.

Every lower n has a strictly smaller bound in (12).
The finite contribution is only of order seven.
The second identity (7) is impossible.

The remaining congruence is k=3 modulo six. If k>=15,
choose j=2. Its exact pole order is

    (k+33)/3>=16>14.

All lower characters have strictly smaller pole bounds.
The finite contribution and the explicit right side in
(7) have poles of order at most fourteen. Thus the top
term cannot cancel here either.

The ONLY indices not excluded by this arbitrary-index
argument are

    k in {3,9}.                                        (13)

For k=3 the third-character escape bound is twelve,
below its finite bound fourteen. For k=9 it is fourteen,
equal to that bound and to the possible right-side pole.
These two cases were NOT excluded by an unjustified
negative-character or finite-cluster comparison.

## 5. A cubic companion index has too little actual pole degree

Suppose k=3, so deg_w Q=9. The polynomial part of sqrt(P)
is H. The complete positive polynomial part of P^(3/2) is

    [P^(3/2)]_+
       =H^3+(3/2)*H*R+(3/8)*T1,
    T1=[R^2/H]_+, deg_w T1<=1.                          (14)

All later binomial terms have negative w degree.
Equation (5) therefore gives

    Q=beta*[H^3+(3/2)*H*R+(3/8)*T1]+alpha*H+E(P),
    beta in C*, alpha in C, E in C[P].                  (15)

Even polynomial terms do not create poles on P=p.

Write

    H=a3*w^3+a2*w^2+a1*w+a0,
    R=r2*w^2+r1*w+r0.

Polynomial division gives the exact formula

    T1=(r2^2/a3)*w
             +(2r2*r1/a3-r2^2*a2/a3^2).                 (16)

At the compulsory support, (2) and the A square initial give

    ord a3=14, ord a2=7, ord r2=10, ord r1=3.
    T1 has slope order six and intercept order at least -1.

All four sheets have w=O(v^-7). Consequently

    T1=O(v^-1) on EVERY sheet.                          (17)

We also need a bound on the REGULAR source cubic H.
Factor its full scaled polynomial at the simple finite root:

    v^7*H(v,v^-7*Y)
      =A(v)*(Y-r(v))*[(Y-m_H(v))^2+J_H(v)],
    A(0)!=0, r(0)=0, m_H(0)=gamma!=0.                   (18)

Before deformation the quartic is the square of that monic
quadratic. The remainder first deforms it at order t=10.
The Hensel-factor and depression calculation of PF179 gives

    J-J_H=O(v^10), m_G-m_H=O(v^10),                     (19)

where m_G is the depressed quartic midpoint.
By (3), ord J_H>=7. Every cubic-balance sheet has
ord(Y-m_G)>=10/3, and (19) gives the same lower bound
at m_H. Thus the quadratic bracket in (18) has order
at least min(20/3,7)=20/3. The other factor is a unit, so

    H=O(v^-1/3) on ALL FOUR sheets.                      (20)

This estimate does not put H in the target field.
On the actual generic fibre P=p, replace R by p-H^2
in (15). Then

    H^3+(3/2)*H*R=-(1/2)*H^3+(3/2)*p*H.

Equations (17),(20) imply Q=O(v^-1) on each sheet.
An index-alpha point therefore has Q pole order at most
alpha. The total quartic sheet mass is four, and EVERY
actual Q pole lies there. Hence its total degree would be

    b_field<=4,

contradicting b_field=8 in (2). This excludes k=3.

Combining this exclusion with (13) proves k=9 and (1).

## 6. Exact scope and replay

The conclusion is a proved necessary degree, not a bounded
search domain chosen for convenience. The remaining degree-27
actual companion must still satisfy the full Jacobian, all
constant residue equations, all retained finite-normalization
conditions and polynomial source regularity. Their incompatibility is
proved in the subsequent adapted cubic moment theorem, by the exact
zero weighted residue and its unique top pole.

The checker verifies the inverse-prefix coefficient orders,
the algebraic third-residue right side, all three exact
congruence coefficients and pole comparisons, whole sextic
moment controls including the small-index finite-cluster
exceptions, and the exact polynomial-part/pole estimate
excluding k=3.

Replay:

    .venv/bin/python research/jacobian/affine-complements/verify_ruling_six_companion_degree27.py

Essential sources are PF164, PF165, PF171, PF176, PF179,
PF182 and
[PF184](RULING_SIX_TWO_MOMENT_RESONANCE_AND_DEGREE30_CUBIC_REDUCTION.md).
No target integrality, trace perfection, formal verification,
external human review or full ruling-six/HC4/JC2 closure is claimed.
