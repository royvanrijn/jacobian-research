# F7 boundary fibres, ten occurring cardinalities, and the nonproper locus

14 November 2026. Claims **GAOF7BND1**, **GAOF7CARD1**, **GAOF7NP1**.
All varieties and fibres are over C unless a finite field is specified.

## Source boundary and remaining open questions

Use exactly Gao's F7 in [arXiv:2608.00222v2](https://arxiv.org/html/2608.00222v2),
Section 4.6 and Section 10. Gao already proves the curve-type simple-root
criterion, triangular reconstruction, generic twelve-point fibres, seven-point
axis fibres, and an isolated quintuple-root stratum. Those are prior results.
The [earlier exact-fibres note](https://ideosphere.ai/papers/keller-fibers)
records the F7 transverse formulas as known data; it does not give the
boundary classification or the occurrence certificates below.

The new results count every boundary fibre, prove that all cardinalities
in {1,2,4,5,6,7,8,9,10,12} occur, and determine the complete nonproper-value
hypersurface. Eleven is impossible. The [separate GAOF7SURJ1 proof](../gao-f7-surjectivity/PROOF.md)
now excludes zero; whether three occurs remains open. The full nonzero-chart multiplicity stratification is NOT claimed;
in particular the possibility of one/two points on that chart is not decided.
No conclusion about JC2 or unrestricted HC4 follows.

## Map and exact reconstruction

For v=(xy,x^2*z1,x^3*z2,x^4*z3), define

    gamma=1-61v1+9012v1^2+(238/19)v1*v2-(35823670/19)v1^3+8v3+v4,
    u1=1+59v1+v2-5178v1^2+19v3,
    u2=-7v1-27v2, u3=-1-234v1-5v2,
    t=gamma*u1, a=gamma^3*u2, b=gamma^4*u3, C=gamma*x.

The sweep data are

    G2=a^2+a+t*b+t^3,
    G3=2t*G2+a+t^3+t^2*b,
    G4=3t^2*G2+b+t^4,
    X1=-det d(G2,G3,G4)/d(t,a,b),
    X=(X1,G2+tX1,G3+t^2*X1,G4+t^3*X1),
    S=X+gamma*(1,t,t^2,t^3).

Then F7=(C,S1/C,S2/C^2,S3/C^3,S4/C^4), using polynomial extensions.
The checker proves polynomiality by the weighted-jet argument used for
F6: E_i=S_i(gamma,gamma*u1,gamma^3*u2,gamma^4*u3)/gamma^i is polynomial;
substituting v_j=e^j*(source coordinate) has no coefficients below weight i.
Positive weights make these finite tests sufficient for full polynomiality.

The stage determinant is 119377, the sweep determinant in (gamma,t,a,b)
is gamma, and the other determinant factors are x^10, gamma^8, gamma
from padding, and C^-10. Consequently det J(F7)=119377 everywhere.
Fibres are zero-dimensional and reduced.

The stage inverse, valid polynomially for arbitrary (gamma,u1,u2,u3), is

    v1=(5u2-27(u3+1))/6283,
    v2=(7(u3+1)-234u2)/6283,
    v3=(u1-1-59v1-v2+5178v1^2)/19,
    v4=gamma-1+61v1-9012v1^2-(238/19)v1*v2
       +(35823670/19)v1^3-8v3.

For x!=0, dividing v_j by x^j recovers the original source uniquely.

## GAOF7BND1: every target-boundary fibre

Write an original target as (c,A,B,D,Z). The x=0 restriction is
(0,H1,H2,H3,H4), where exact weighted jets give

    H1=y,
    H2=-722y^2+61z1,
    H3=(72851307y^3-24938y*z1-1957z2)/19,
    H4=(130524671288y^4+1724330y^2*z1-111055y*z2
        -22572z1^2-361z3)/19.

The triangular diagonal is (1,61,-103,-19), so this is an isomorphism
onto c=0 and contributes exactly one point to every boundary fibre.

Since C=x*gamma and gamma=1 at x=0, the rest of that preimage is the
disjoint piece gamma=0,x!=0. Here the scaled sweep formulas reduce to

    E1=1, E2=u1, E3=u1^3+u1^2+u2, E4=u1^4+u1^3+u3.

Thus A=1/x. There is no such source when A=0. When A!=0 there is exactly
one, obtained by

    x=1/A, u1=B*x^2,
    u2=D*x^3-u1^3-u1^2,
    u3=Z*x^4-u1^4-u1^3,

then the stage inverse with gamma=0. All divisions are nonzero and these
formulas satisfy every target equation. Therefore

| Boundary target | Fibre size |
| --- | ---: |
| c=0, A=0 | 1 |
| c=0, A!=0 | 2 |

This is the whole boundary fibre, not just its x=0 restriction. In particular
the origin has one point; (0,1,0,0,0) has two rational source points.

## Nonzero chart: the prior simple-root criterion, replayed

For c!=0 set Y=(cA,c^2 B,c^3 D,c^4 Z), and put

    b(t)=Y4-3Y2*t^2+2Y1*t^3-t^4,
    a(t)=Y3-2Y2*t+Y1*t^2-t^3-t^2*b(t),
    R(t)=a(t)^2+a(t)+t*b(t)+t^3+Y1*t-Y2,
    gamma(t)=Y1-X1(t,a(t),b(t)).

R is monic of degree twelve. All sweep residuals are divisible by R, and
gamma-R' is divisible by R. The latter is Gao's general curve-type lemma;
the checker verifies the specific polynomial identity without exceptions.
Exactly the simple roots lift, using x=c/gamma, u1=t/gamma,
u2=a/gamma^3,u3=b/gamma^4, followed by the stage inverse.

For any given target let d1=deg gcd(R,R') and d2=deg gcd(R,R',R''). Then

    #F7^(-1)(target)=12-2d1+d2.

This gives an exact pointwise counting algorithm together with the boundary
table. It does not decide which of the parameter-dependent multiplicity
strata are nonempty, which remains a separate classification task.

## GAOF7CARD1: ten certified cardinalities

Twelve and seven are prior results, replayed at Y=(1,2,-1,3) and Y=0.
Ten occurs at Y=(-1,0,0,1), where gcd(R,R')=t, giving one ordinary double
root and ten simple roots. Nine occurs at Y=(0,0,-1,0), where the gcd is
t^2 and R has a triple root at zero and nine simple roots.

For eight, let alpha be any root of the irreducible quintic

    f(z)=16z^5+48z^4+48z^3+10z^2-13z-5.

Set beta=alpha*(6alpha^3+12alpha^2+7alpha+1)/(2alpha+1) and

    Y=(4alpha^3+6alpha^2+2alpha-beta, alpha^2+alpha, alpha, beta).

The denominator is nonzero. An exact gcd over Q(alpha) is t^3. The root
zero is therefore quadruple and all other roots are simple, giving eight
source points. The checker uses exact algebraic-number arithmetic, not an
approximate root. Irreducibility transports the identity to every conjugate.

The cardinalities four, five and six use the following positive lifting
certificates. For a list [q0,...,qd], the polynomial is sum(qj*t^j).

| Simple roots | p | Y mod p | Multiplicity factors mod p |
| --- | ---: | --- | --- |
| 4 | 17 | (11,10,15,10) | Q1=[8,1,5,13,1], Q2=[13,12,15,14,1] |
| 5 | 13 | (11,11,8,8) | Q1=[6,4,9,8,7,1], Q2=[5,10,1], Q3=[11,1] |
| 6 | 13 | (0,9,12,10) | Q1=[3,11,7,4,9,6,1], Q2=[6,10,10,1] |

In each row R=product(Qm^m) modulo p. Each Qm is monic and squarefree,
and the factors for different m are coprime. Form the twelve integer
coefficient equations of R-product(Qm^m)=0 in degrees zero through eleven.
Unknowns are Y1..Y4 followed by all nonleading Qm coefficients in increasing
m and increasing polynomial degree. The Jacobian minors in the checker are:

- four points: 12 by 12, determinant 3 modulo 17;
- five points: 12 by 12, determinant 10 modulo 13;
- six points: 12 by 13, first twelve columns have determinant 7 modulo 13.
  Fix the remaining coefficient, the t^2 coefficient of Q2, to the integer 10.

Multivariable Hensel lifting now gives a p-adic solution of twelve equations
in twelve variables in each row. These solutions are algebraic over Q:
in the finitely generated field of their coordinates, differentiating the
equations gives an invertible Jacobian times the vector of differentials
equal to zero. All coordinate differentials vanish, so in characteristic
zero that field has transcendence degree zero. It is a number field and
admits an embedding into C. The lifted discriminants and pairwise resultants
remain nonzero, since they are units modulo p. Thus the embedded complex
polynomials have exactly the indicated multiplicities and exactly four,
five and six simple roots. Taking c=1 produces actual F7 targets of those
cardinalities.

The modular Y-values in this table are residues of algebraic targets;
they are NOT asserted to be the exact rational targets. The finite-field
search is not itself a lifting proof: the separate checker reconstructs
the factorization equations, differentiates them in SymPy, and checks the
nonzero minors. Its default run does not repeat the search.

The boundary table supplies one and two. Consequently

    {1,2,4,5,6,7,8,9,10,12}

all occur globally. Eleven is impossible for any degree-twelve simple-root
count: one nonsimple root consumes at least two of the twelve roots, while
a squarefree polynomial gives twelve. It is also absent from the boundary.
The separate GAOF7SURJ1 proof excludes zero, so the only global cardinality
still undecided is three, requiring c!=0. This does not mean that the strata of the known sizes have
been exhausted.

There is a further exact exclusion. Complete the square:

    R=(a(t)+1/2)^2+L(t),
    L=t*b(t)+t^3+Y1*t-Y2-1/4.

Here a+1/2 is monic of degree six, while L has degree five and leading
coefficient -1 for every Y. If R were a square of a monic degree-six H,
then (H-a-1/2)*(H+a+1/2)=L. The second factor has degree six and the first
is either zero or has nonnegative degree, contradicting deg L=5 and L!=0.
Hence R is never a square; every all-even multiplicity partition is excluded.
In particular a hypothetical omitted target requires odd multiplicities
at least three. The complete empty-fibre exclusion is in GAOF7SURJ1.

The nonsquare reduction gives exactly ten hypothetical empty-fibre types,
all excluded in GAOF7SURJ1:

    9+3, 7+5, 7+3+2, 6+3+3, 5+5+2, 5+4+3,
    5+3+2+2, 4+3+3+2, 3+3+3+3, 3+3+2+2+2.

Each has deg gcd(R,R') at least seven. A uniform bound at most six would
prove surjectivity, but has not been established. GAOF7SURJ1 instead excludes
the two power-factorization families covering these ten types.
The eight possible three-point types consist of three simple roots plus

    9, 7+2, 6+3, 5+4, 5+2+2, 4+3+2, 3+3+3, 3+2+2+2.

They require gcd degree at least five. A uniform bound at most four would
exclude three-point fibres as well. Neither uniform gcd bound is proved here.

## GAOF7NP1: normalization and the full nonproper hypersurface

The finite root incidence R(t;Y)=0 is polynomially isomorphic to A4 with
coordinates (gamma,t,a,b). The inverse formulas above are polynomial even
when gamma=0, and the sweep residual identities give both compositions.
It is finite over Y because R is monic and gamma,a,b are polynomial in
(t,Y). Its critical divisor is the irreducible A3 given by gamma=0, by
the determinant det dS=gamma. The discriminant D(Y)=disc_t(R) therefore
has irreducible zero set, the finite image X(A3).

At Y=(-1,0,0,1) there is only one ordinary double root, zero, and
R_Y2(0)=-1. The other ten roots are simple. This gives a first-order
discriminant zero in the Y2 direction, so D is reduced as well as
irreducible. The single-double-root locus is a nonempty open of D;
there is a unique critical source above it. Thus X:A3->D is finite and
birational. Since A3 is normal, it is the normalization of D.

For the original target coordinates define

    E(c,A,B,D,Z)=disc_t R(t;cA,c^2 B,c^3 D,c^4 Z)/c^5.

The fourth target-coordinate symbol D in this formula is distinct from
the discriminant D(Y). The quotient E is polynomial and

    E(0,A,B,D,Z)=-182528*A^5.                         (*)

Here is a proof without expanding the large discriminant. At c=0,

    R0=t^5*H7(t),
    H7=(t-1)*(t^6+t^5+t^4-t^3-t^2-t+1),
    H7(0)=-1, disc(H7)=713.

The monic factors t^5 and H7 are coprime. Over the formal field
K=C(A,B,D,Z), Hensel factorization gives monic factors P5(t,c),Q7(t,c)
reducing to t^5,H7 respectively. Put c=rho^4,t=rho*h. The checker verifies
that R(rho*h;rho^4*A,rho^8*B,rho^12*D,rho^16*Z)/rho^5 is polynomial in
rho, with constant term -h^5+A*h. Hence P5(rho*h,rho^4)/rho^5 reduces to
h^5-A*h. Its discriminant is -256*A^5, nonzero over K. The five cluster
roots consequently have pairwise differences of order rho. Their
discriminant contributes rho^20=c^5 with leading coefficient -256*A^5.
The other factor contributes 713, and the squared resultant between the
two factors tends to H7(0)^10=1. This proves exact c-valuation five and
leading term -256*713*A^5=-182528*A^5. Since the original discriminant
is a polynomial in c,A,B,D,Z, the generic-field identity proves polynomial
divisibility and (*) identically, including A=0.

The complete nonproper-value locus is V(c*E). On c!=0, off the discriminant
all twelve inverse branches remain bounded locally. At a multiple root,
approach by simple-root points in the smooth irreducible root incidence
A4. Then gamma=R' tends to zero and x=c/gamma diverges. This establishes
exactly the discriminant on that chart.

Every c=0 target is also nonproper. Approach it with the last four original
target coordinates fixed. Choose a simple root alpha of H7 other than 1.
Such roots exist, are nonzero, and have R0'(alpha)!=0. A nearby analytic
root gives gamma->R0'(alpha)!=0 and x->0. Bounded source coordinates
would force all vi=x^i*(source coordinate) to zero, hence gamma->1,u1->1
and t->1, contradicting alpha!=1. Thus a source branch escapes at every
boundary target, including those with two finite source points.

Finally E is irreducible: on c!=0 the weighted target change is an
isomorphism and E is the irreducible discriminant up to a unit. Its
nonzero restriction (*) excludes a factor c. The full nonproper locus
has exactly two irreducible hypersurface components. Their intersection
has support c=A=0 and, as a restriction of E to c=0, multiplicity five
along that hyperplane. Its support is precisely the boundary one-point
locus, but the fifth-power intersection is additional information.

## Incomplete attempts and verification limits

A bounded thirty-second Singular calculation of a sixfold-root ideal
returned only `halt 1`. It is retained in sixfold-input.json/.sing,
sixfold-groebner.log and sixfold-attempt.json. It proves no unit ideal
or multiplicity bound. The first exploratory finite-field histogram
miscounted some gcd-degree-two cases; its receipt is preserved as
lift-preflight-v1.json and explicitly superseded by the corrected run.
Neither exploratory histogram is used as a characteristic-zero exclusion.

Run from the repository root:

```sh
timeout 60s .venv/bin/python research/jacobian/gao-f7-fibres/verify.py
```

The exact checker covers algebra, weighted jets, rational and number-field
witnesses, independent modular Jacobians and source hashes. Hensel lifting,
normalization, formal discriminant degeneration and analytic escape are
written geometric arguments. No independent end-to-end review, formal
verification, external review, author endorsement or novelty priority is
claimed. MATH_STATUS.json remains authoritative.

## Subsequent three-point exclusion

The eight lists above are the exhaustive combinatorial possibilities before
coefficient obstructions. The later
[GAOF7CUBE1 proof](../gao-f7-three-point-certificate/PROOF.md) excludes
(9), (6,3), (3,3,3) by ruling out every cubed cubic divisor, including
repeated/shared roots. The later
[GAOF7MULT1 proof](../gao-f7-linear-cube/HIGH_ROOT_PROOF.md) proves the
sharp maximum root multiplicity five and excludes(7,2). The subsequent
[GAOF7FIVE1 proof](../gao-f7-quintuple/PROOF.md) classifies all118 fivefold
targets and leaves only(4,3,2),(3,2,2,2) and global cardinality three open.
The partition enumeration in this note and its checker is unchanged.
