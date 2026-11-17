# Positivity on the minimal resolution, with a direct ruled-surface exclusion

14 November 2026. Work over C. This follows the positivity statement of the
[7 November campaign note](sources/minimal-resolution-ramification-positivity-swarmHQ-root-sol56-20261107T122700Z.md).
The new contribution is a direct exclusion of its ruled endpoint, replacing
the inherited smooth-finite genus and cubic-scroll donor interfaces. The
consequences include a general-line genus test and an exact contraction
procedure for complete numerical boundary data. No novelty priority,
formal verification or independent external review is claimed.

## 1. Statement and exact model

Let g:Y -> P2 be a generically finite morphism of degree N>1, with Y smooth
irreducible projective. Let U be an open subset isomorphic to the whole A2,
assume g is etale on U, and assume g(U) misses a line L. Suppose no smooth
rational (-1)-curve is g-contracted. Write

    T = g*L as an actual effective divisor,
    H = [T], K = K_Y, R = div(det(dg)), Q = 2K+3H.

R is the full determinant divisor, including contracted components. Then:

1. Q and R are nef and big.
2. Their common null curves are exactly the g-contracted smooth rational
   (-2)-curves.
3. p_a(R)>=2 and h0(Y,O(Q))>=2.
4. If gamma is the genus of a smooth general member of |H|, then gamma>=2,
   rho(Y)>=2, and

       p_a(R)=2+9*gamma-rho(Y),
       Q^2=16+24*gamma-4*rho(Y)-3*N>0,
       2 <= rho(Y) <= 9*gamma.

For every irreducible curve C of arithmetic genus a and square s, put
m=H.C. The boundary test is

    Q.C = 4*a-4-2*s+3*m >= 0,
    R.C = 2*a-2-s+3*m >= 0.

Each zero occurs precisely when (a,s,m)=(0,-2,0). In particular, a smooth
rational curve with m>0 must satisfy

    s <= floor((3*m-5)/2).

This last formula excludes a rational square-zero curve of H-degree one.
It applies after smooth H-null (-1)-curves have been removed. It is not a
test on the singular Stein surface, nor on an arbitrary unminimalized SNC
resolution. Passing a boundary numerical test does not certify existence
or all-curve nefness.

For a hypothetical plane Keller map, normalize target P2 in its source
function field and take the minimal resolution of that finite normal
surface. Zariski's Main Theorem embeds the same whole affine source into
the normalization. The resolution is an isomorphism there. Its contracted
curves are precisely those of g, and none is a smooth rational (-1)-curve.
Thus this is the actual model to which the theorem applies.

## 2. Whole-source boundary facts

Y is rational. If B_i are its divisorial boundary components, localization
and O(U)^*=O(Y)^*=C*, Cl(U)=0 give

    direct_sum Z[B_i] = Pic(Y).

In particular the number of boundary curves is rho(Y); their classes are
an integral basis. Two divisors supported on the boundary and linearly
equivalent are equal as actual divisors: the rational function relating
them restricts to a unit on U, hence is constant.

There are no isolated boundary points outside these curves. Indeed, if p
were such a point, choose a normal affine neighborhood W disjoint from
the boundary curves and other isolated boundary points. The two coordinate
functions of U=A2 on W minus p extend across p by normality (Hartogs).
They define W -> A2. Composing with A2=U -> Y agrees with the inclusion
W -> Y on the dense open W minus p; separatedness makes them equal on W.
This would put p in U, a contradiction. This purity step is needed for the
Euler calculation below.

H is nef, H^2=N, R is effective and R~K+3H. The supports of T and R are
in the boundary. If a prime C is not g-contracted, its coefficient in R
is e_C-1 by the characteristic-zero DVR different formula. If it maps to
L, its coefficient in T is e_C. These statements are not applied to a
contracted prime.

Also Y cannot be P2. Its unique boundary curve must be a line C, because
its class is an integral Picard generator. Then T=mC and R=(3m-3)C.
H is ample, so g is finite and C maps onto L. The preceding different
formula gives 3m-3=m-1, hence m=1 and N=1, a contradiction.

## 3. Direct ruled endpoint: no genus or cubic-scroll interface

Suppose Y=F_n and a ruling fibre F has H.F=1. Write S^2=-n, S.F=1,
H=S+bF, b>=n, and N=2b-n>1. The canonical formula gives R.F=1.

### 3.1 The ample case b>n

H is ample, hence g is finite. R has exactly one horizontal component C,
of coefficient one and C.F=1. Likewise T has exactly one horizontal
component B, of coefficient one and B.F=1. These are smooth sections:
their maps of degree one to the ruling base are finite birational maps
to the normal curve P1.

C and B are distinct. If C=B, its coefficient one in T would give e_C=1,
whereas its coefficient one in R would give e_C=2. Since rho(F_n)=2,
these are all the boundary curves. Thus there are no vertical components
of either divisor and

    T=B, R=C, U=F_n minus (B union C).

Let t be the number of distinct intersection points of B and C. Compactly
supported Euler characteristic gives

    1=chi_c(U)=4-(2+2-t)=t.

Tangency multiplicities do not affect this count. Thus B meets C at just
one point p. The map g_B:B=P1 -> L=P1 has degree H.B=H^2=N.
Because g*L=B is smooth and reduced, the ramification divisor of g_B is
exactly R restricted to B. Locally take a target equation u for L; its
pullback is a local source coordinate normal to B. The determinant of dg
restricted to B is then the derivative of the other target coordinate
along B. This proves the divisor identity, including at p.

Riemann--Hurwitz now puts all 2N-2 ramification at p. In characteristic
zero its contribution is e_p-1<=N-1. Therefore 2N-2<=N-1, impossible
for N>1. This directly replaces the earlier cubic-scroll donor dependency.

### 3.2 The nonample case b=n

Now N=n>=2 and S is g-contracted. Hence S occurs in R with positive
coefficient. Since R.F=1, that coefficient is one and all other primes
of R are fibres. The canonical formula gives

    R=S+sum r_j F_j,  sum r_j=2N-2.

Every F_j has H-degree one, hence is not g-contracted. The finite DVR
formula gives r_j=e_j-1<=N-1. At least two distinct F_j are required.
Together with S this gives three boundary curves on a surface with
Picard rank two, a contradiction. Thus no such ruled endpoint exists.

## 4. Q is nef, nonzero and big

If Q were not nef, the surface cone theorem gives a Q-negative
K-negative extremal ray, because Q is nonnegative on the K-nonnegative
part of the cone. The surface contraction theorem gives three cases.

A birational contraction has a smooth rational (-1)-curve E. Minimality
excludes H.E=0; otherwise Q.E=-2+3H.E>=1. A contraction to a point gives
Y=P2, already excluded. A contraction to a curve has general fibre F
with K.F=-2. Hodge index and H^2>0 imply H.F>0, so Q.F<0 forces H.F=1.

This last surface really is a Hirzebruch surface. All contracted fibre
components lie on one numerical ray, hence have square zero. A connected
reducible fibre would have a component of negative square, a contradiction.
A multiple fibre is excluded by H.F=1. Flatness over the smooth curve
and arithmetic genus zero make every fibre smooth P1; rationality makes
the base P1. Section 3 excludes this endpoint. Hence Q is nef.

Since R=(Q+3H)/2, R is nef and R^2>0. If Q is numerically zero, K^2=9N/4.
Rationality gives K^2=10-rho<=9; integrality forces N=4 and rho=1. A
smooth rational projective surface with rho=1 is P2, already excluded.

If Q^2=0, nonzero nefness gives Q.H>0, so K.Q=-3H.Q/2<0. Riemann--Roch
and H2(nQ)=0 for large n give arbitrarily many sections of nQ. Remove
the fixed part, leaving a nonzero movable class M. Q.M=0, M^2>=0, and
Hodge index imply M^2=0 and M is a positive numerical multiple of Q.
Two general movable members have no common component and zero intersection,
so have no base points. Stein factorization gives a smooth general connected
fibre F with F^2=0 and K.F<0. Adjunction makes it P1 with K.F=-2. But

    0=Q.F=-4+3H.F

contradicts integrality. Thus Q^2>0, and both Q and R are nef and big.

## 5. Null curves, sections and identities

If Q.C=0, bigness and Hodge index give C^2<0. The relation
2K.C+3H.C=0 gives K.C<=0; adjunction gives K.C>=-1. If H.C>0 the relation
instead gives K.C<=-2. Thus H.C=K.C=0, and adjunction yields p_a(C)=0,
C^2=-2. Such an integral curve is smooth rational and g-contracted.
The converse is immediate. R has the same null curves by R=(Q+3H)/2.

Q.R>0, and adjunction gives Q.R=2p_a(R)-2. Thus p_a(R)>=2.
The sequence O(K) -> O(K+R) -> omega_R, and H0(K)=H1(K)=0, identify
H0(Q) with H0(omega_R). Cartier-divisor duality gives

    h0(omega_R)=p_a(R)+h0(O_R)-1>=2.

No reducedness or connectedness of R is needed. A smooth general H-member
is connected (H is nef and big) and hence irreducible; adjunction gives
K.H=2gamma-2-N. Together with K^2=10-rho this proves the identities in
Section 1 and rho<=9gamma.

## 6. Excluding general-line genus one

Genus zero is already excluded by rho<=9gamma. Suppose gamma=1. For a
smooth general H-member C, the adjunction sequence and rationality give

    h0(K+H)=h0(omega_C)=1.

Let D>=0 represent K+H. Since H.D=2gamma-2=0, every component of D is
g-contracted. A contracted irreducible curve E has E^2<0. Minimality
and adjunction imply K.E>=0: the only possible negative case would be a
smooth rational (-1)-curve. Thus K.D>=0. If D is nonzero, Hodge index gives

    D^2<0, but D^2=(K+H).D=K.D>=0,

a contradiction. Hence D=0 and K+H~0. Consequently R~2T, and the whole-A2
boundary-basis fact strengthens this to equality of actual divisors R=2T.
At least one component E of T is not contracted, since H.T=N>0. Its
coefficient in T is e>=1, so its coefficient in R is e-1, whereas R=2T
would give 2e. This is impossible. Therefore gamma>=2.

This is a scoped rederivation; no novelty is asserted. Classical work on
low-genus Keller fibres must be compared with the precise general-target-line
statement before making a priority or stronger polynomial-component claim.

## 7. Transport from a complete smooth boundary

For beta:Y' -> Y the blowup of a point, g'=g beta satisfies

    H'=beta*H, K'=beta*K+E,
    R'=beta*R+E, Q'=beta*Q+2E.

In particular Q'.E=-2, R'.E=-1 and p_a(R')=p_a(R)-1. A raw negative
intersection here is a blowup artifact, not an exclusion.

Conversely remove only a smooth rational E with E^2=-1 and H.E=0. If
v_i=C_i.E for surviving curves, the exact smooth blowdown rules are

    M_new=M_surv+v*v^T,
    (K_new.C_i_new)=(K_old.C_i)-v_i,
    H_new.C_i_new=H_old.C_i,
    p_a(C_i_new)=p_a(C_i)+v_i*(v_i-1)/2.

The last rule is essential when a boundary ceases to be SNC. All contracted
curves are in the boundary because g is etale on U, so a complete boundary
can be reduced to the required relatively minimal smooth model this way.
Contraction of every H-null curve at once computes a singular Stein model;
it is a different operation and cannot be substituted here.

The implementation checks these formulas and the changes K^2:+1,
R^2:+1, Q^2:+4, p_a(R):+1 at every smooth blowdown. The degree and gamma
are invariant. It tests boundary intersections and numerical identities,
and explicitly supplies no realization certificate.

## 8. Current application boundary

The quartic Cox packets do not yet supply a complete projective boundary
matrix and pullback on the actual minimal resolution. Neither does the
surviving F2 (75,125) programme. Both applications remain UNKNOWN at this
interface; their local packets cannot be completed by inventing curves.

The retained 23-component (72,108) packages are useful numerical regressions
but concern a transformed Laurent pair with Jacobian X^2, not a Keller map
on the stated affine chart. The source's explicit transformed-Poisson
correction must be retained. Those packages cannot provide a new Keller
exclusion under the theorem. The replay records this failed eligibility
test, rather than ignoring its affine ramification.

## Sources and assurance

- The pinned campaign positivity note supplies the target theorem and its
  Q-bigness/null-curve strategy. Its genus and smooth-finite dependencies
  are retained verbatim for provenance; Sections 2--4 replace that chain.
- [Hassett, Rational surfaces over nonclosed fields](https://www.claymath.org/wp-content/uploads/2022/03/cmip08c.pdf),
  Sections 1.3 and 2.2, especially Theorem 2.14: ruled surfaces and the
  surface cone/contraction theorem. These are named classical imports.
- [Stacks, Riemann--Hurwitz, Tag 0C1B](https://stacks.math.columbia.edu/tag/0C1B):
  characteristic-zero ramification degree and local exponent e-1.
- Other named standard imports: divisor localization, normal Hartogs
  extension, Bertini, Hodge index, rational-surface Riemann--Roch and
  Serre duality, and compactly supported Euler additivity over C.

The proof is written geometry with exact algebraic regressions, not a
formal proof. The checker does not certify the cone theorem, the boundary
purity/Euler argument, or existence of a global Keller compactification.
No surviving quartic/F2 case or unrestricted JC2/HC4 claim is closed.
