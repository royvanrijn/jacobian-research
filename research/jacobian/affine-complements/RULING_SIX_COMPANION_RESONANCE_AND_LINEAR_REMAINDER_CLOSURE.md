# Companion resonance excludes every linear square remainder

PF2RULED176. 17 November 2026. Actual complex polynomial plane Keller maps
in the complete marked quadratic obstruction class.

In the regular normal form P=H3^2+d*z+e*w+f, the coefficient d
cannot vanish. This excludes the ENTIRE linear-remainder branch,
with unrestricted coefficient degrees and unrestricted actual companion
degree. The ruling-six problem is reduced to a genuinely quadratic
remainder d*z+e*w+f with d nonzero. [PF188](RULING_SIX_SQUARE_DEFORMATION_AND_MARKED_QUADRATIC_EXCLUSION.md) excludes the full retained marked
quadratic class. Automatic trace perfection and HC4/JC2 remain open.

The proof uses a stronger consequence of the actual companion than
the rational-exactness hierarchy alone: a nontrivial CONSTANT linear
relation among odd half-power residues. At an actual A pole in the
linear branch, those residues have strictly increasing pole orders.
The highest term cannot cancel.

## 1. Exact hypotheses and retained inputs

Retain the original source and full normalizations of PF165--PF175:

    E=C(u,w), F=C(P,Q), L=C(x,t), u=t^2,
    w^2=N+M*z, N=u*Sigma^2, M=u*A,
    dP wedge dQ=du wedge dw/M.                            (1)

The actual companion can be globally reduced by subtracting polynomials
in P. Relabel that reduced companion Q. PF165 gives

    deg_w Q=3k, k>=3 odd,
    leading_w Q=beta*(h/M)^k, beta in C*.                 (2)

PF171 supplies regular source functions, without placing them in F,

    P=H3^2+R2,
    H3=(h*w+K)*z+B0*w+e0, leading_w H3=h/M,
    R2=d(u)z+e(u)w+f(u), d,e,f in C[u].                  (3)

PF165 supplies at least one actual A pole cluster, at an off-M
base value c. Its integer slope is L>=1. With v=u-c and
Y=v^L*w, its square initial has

    v^L*H3(v,v^-L*Y) modulo v = C*Y*(Y-c0)^2,
    C,c0 in C*.                                         (4)

Here L in (4) denotes the LOCAL SLOPE, distinct from the original
source field L in (1). The slope is integral by the A classification.
All coefficients of the scaled cubic on the left of (4) belong to
C[[v]]. This also follows from the A initial of P and the lower
weighted degree of the linear remainder used below.

PF173 proves that G=M^2*h is NOT a cube in C(u): it has at least
two distinct roots, and the actual leading primitive equation
excludes every such cube. Consequently

    a^3=M/h, K0=C(u,a)

is a cyclic extension of degree exactly three. Indeed M/h=M^3/G
has the inverse cube class of G. No unproved Kummer irreducibility
assumption is added.

Finally, the actual-pole support theorem of PF175 says, if d=0,

    a_e=ord_c e < L                                     (5)

at every actual A supporting base value. The zero polynomial has
infinite order, so (5) also makes e nonzero. PF172 already excludes
the constant-remainder case. These are hypotheses and conclusions
inside the COMPLETE obstruction class; no claim about arbitrary
polynomials of w degree six outside that class is made.

## 2. The actual companion is a constant half-power polynomial to high order

Let X be a new formal variable and select the inverse at P infinity

    w=w(u,X)=a*X+c0(u)+c1(u)*X^-1+...,
    P(u,w(u,X))=X^6,                                     (6)

in K0((X^-1)). To construct it, set w=X*V and use the equation
X^-6 P(u,X*V)=1. At X^-1=0, V=a is a simple root, since
the V derivative is 6/a. Constant-unit recursion gives the
unique series with this leading term, with all coefficients in K0.

Write the actual polynomial companion evaluated on (6) as

    Q(u,w(u,X))=sum_j q_j(u)*X^j, q_j in K0,
    j<=3k, q_(3k)=beta.                                  (7)

Extend d/du uniquely to the separable algebraic field K0 and
differentiate while holding X fixed. The full actual Jacobian gives

    partial_u Q(u,w(u,X))=-1/(M*P_w(u,w(u,X))).            (8)

The right side starts at X^-5, with leading coefficient -a/(6M).
Hence q_j'=0 for EVERY j>-5.

The constants of K0 under this derivation are exactly C. For completeness,
if b'=0 and b is algebraic over C(u), differentiate its monic minimal
polynomial. The resulting polynomial in b has smaller degree; minimality
forces every differentiated coefficient to vanish. Thus all coefficients
belong to C, and algebraic closedness of C gives b in C. Therefore

    q_j in C for j>-5.                                  (9)

Let zeta be a primitive cube root of unity and let sigma(a)=zeta*a.
Both sigma(w(u,X)) and w(u,zeta*X) solve (6) with leading term
zeta*a*X. Uniqueness gives

    sigma(w(u,X))=w(u,zeta*X),
    sigma(q_j)=zeta^j*q_j.                               (10)

Combining (9) and (10),

    j>-5 and 3 does not divide j => q_j=0.               (11)

In particular the X^-1, X^-2 and X^-4 coefficients vanish.
The square root of P chosen to have leading_w (h/M)*w^3 maps
under (6) to X^3, because (h/M)*a^3=1. Consequently there are
CONSTANTS c_n=q_(3n), 0<=n<=k, with c_k=beta nonzero, such that

    Q-sum_(n=0)^k c_n*P^(n/2)=O(X^-3)=O(w^-3).           (12)

This is a formal Laurent statement, not polynomial termination,
a target change by an odd half-power, or an integrality assertion.
Order is preserved on extending the coefficient field and substituting
the invertible formal parameter X^-1; thus (12) is also valid in the
original Laurent series field C(u)((w^-1)).

## 3. A constant relation among odd characters

Use the coefficient residue convention of PF175,

    res_w A=[w^-1]A, m_n=res_w P^(n/2).                   (13)

The polynomial Q has zero coefficient residue. So do constants and
the even half-powers P^(n/2). Taking res_w of (12) yields the exact
RATIONAL FUNCTION IDENTITY

    sum_(1<=n<=k, n odd) c_n*m_n=0,
    c_n in C, c_k=beta nonzero, k>=3 odd.                (14)

This is stronger than the separate rational-exactness statements
m_n*du/M exact. No converse from that hierarchy is used.

In the linear branch d=0, m_1=0 exactly, as follows either from
PF175 or the w degree of (e*w+f)/(2H3). Thus the top character
in (14) has index at least three.

## 4. Exact pole growth at any actual A support

Suppose d=0. Fix one actual A supporting base value, with the
notation of (4),(5), and put

    a_e=ord_c e>=0, q=L-a_e>0,
    e=e_a*v^a_e+..., e_a in C*.

The full scaled functions have coefficients in C[[v]] and satisfy

    H3(v,v^-L*Y)=v^-L*(H0(Y)+v*H1(v,Y)),
    H0(Y)=C*Y*(Y-c0)^2,
    R2(v,v^-L*Y)=v^-q*(e_a*Y+v*R1(v,Y)).                (15)

The second formula includes f, since f is regular and q>0.
Here H1 and R1 are polynomials in Y with power-series coefficients.
The leading Y^3 coefficient of H0 is a unit. Therefore Laurent
expansion of every inverse power of the scaled H3 at Y infinity
has v-integral coefficients. Taking a fixed coefficient commutes
with reduction modulo v.

For a positive odd n, the complete binomial formula of PF175 is

    m_n=sum_(j=(n+1)/2)^floor((3n+1)/4)
             binom(n/2,j)*res_w(H3^(n-2j)*R2^j).          (16)

The upper bound is allowed to be loose for a linear R2.
All omitted terms with nonnegative H3 exponent are polynomial;
terms beyond that bound have w degree below -1.

Under w=v^-L*Y, coefficient residue scales as

    res_w A(w)=v^-L*res_Y A(v^-L*Y).                     (17)

Thus the j term of (16) has v order AT LEAST

    -L*(n+1)+(L+a_e)*j.                                 (18)

Write n=2r+1>=3, so r>=1. The first index is j=r+1,
and (18) there is -(r+1)*q. Every later j has strictly
larger lower bound, with margin L+a_e>0.

The leading coefficient at this FIRST index is exactly

    binom(r+1/2,r+1)*(e_a^(r+1)/C)
                 *res_Y[Y^r/(Y-c0)^2]
      =binom(r+1/2,r+1)*(e_a^(r+1)/C)*r*c0^(r-1),        (19)

which is NONZERO. The last equality follows by expanding
Y^(r-2)*(1-c0/Y)^-2 at infinity. Characteristic zero,
C,c0,e_a nonzero and r>=1 are all essential here.

Equations (18),(19) prove, with no genericity exception,

    ord_c m_(2r+1)=-(r+1)*(L-a_e), EVERY r>=1.            (20)

The positive odd characters m_3,m_5,m_7,... therefore have strictly
increasing pole orders. This proof retains the entire Laurent
series and all original normalized primes. It does not identify
the degree of a first coefficient response with a residue-field
degree or discard a finite branch.

## 5. Exclusion and exact remaining scope

In the constant identity (14), the top nonzero term c_k*m_k
has a pole at c of strictly larger order than every smaller odd
character. The character m_1 is zero. The leading pole coefficient
of c_k*m_k cannot cancel, contradicting (14).

Therefore

    P=H3^2+d*z+e*w+f in this complete class => d!=0.      (21)

This excludes BOTH possible remainders of w degree zero or one.
A cubic function H3 and the quadratic remainder remain source
functions; neither is promoted to the actual target field or
assumed integral over it. Equation (21) does not exclude the
remaining case d!=0, even if its local order is large. [PF188](RULING_SIX_SQUARE_DEFORMATION_AND_MARKED_QUADRATIC_EXCLUSION.md) gives the
global exclusion of that retained marked quadratic case. No bounded companion degree,
bounded coefficient degree, finite experimental prefix, or
trace/perfection assumption is used to prove (21).

## 6. Exact controls and replay

The checker records exact internal controls for the inverse-chain
rule, the cyclic coefficient characters, coefficient-residue
scaling, the first nonzero binomial coefficient, and complete
finite moment sums on full cubic-plus-linear coefficient models.
It also shows why the noncube hypothesis cannot be omitted:
P=w^6+u,Q=w is an actual polynomial Jacobian pair outside the
complete obstruction class, with a nonresonant constant X term.

The former positive D=30 coefficient control M=u^2,
h=(u-1)^14,e=(u-1)^6 passes the stated leading primitive equation,
the first/third rational-character test, central units, and A
support. It nevertheless lies in the now excluded linear class.
No actual companion for that control was ever asserted. Its
moments at u=1 have orders -2,-3,-4,... by (20), so it cannot
satisfy any actual companion identity (14).

Replay with

    .venv/bin/python research/jacobian/affine-complements/verify_ruling_six_companion_resonance_linear_closure.py

The arbitrary-degree theorem is the written argument (6)--(21).
The finite controls do not certify all mathematical dependencies.
No external human review, formal proof, target landing, or full
ruling-six/HC4/JC2 closure is claimed.

Essential proof sources:
[PF165](RULING_SIX_SIMPLE_ESCAPE_GRAPHS_AND_QUARTIC_FACE_REDUCTION.md),
[PF171](RULING_SIX_REGULAR_CUBIC_APPROXIMATE_ROOT_AND_QUADRATIC_REMAINDER.md),
[PF173](RULING_SIX_LEADING_PRIMITIVE_AND_POLYNOMIAL_ODE.md),
[PF175](RULING_SIX_ODD_HALF_POWER_RESIDUE_HIERARCHY.md).
