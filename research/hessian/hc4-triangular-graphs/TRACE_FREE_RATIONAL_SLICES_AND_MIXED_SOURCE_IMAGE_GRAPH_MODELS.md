# Trace-free rational slices and an all-degree mixed source-image exclusion

HC4MYMV0. 16 November 2026. Characteristic zero over C. This gives a
GLOBAL exactness criterion for a rational-slice source image, and excludes
an explicit all-degree family of surjective etale plane-source images
from having an etale plane TARGET. It does not exclude every rational
slice for the 2+3 model. Full Wright/HC4/JC2 remains OPEN.
No novelty, independent replay, formal verification or review is asserted.

The previous goal turn made PROGRESS with MU0's exact first mixed-fibre
order and marked quadratic trace pole. The present argument uses a regular
primitive, rather than more moments or a monodromy-only obstruction.

## 1. A global criterion that retains the rational slice

Let R=C[x,y], L=Frac R, h(y) monic of degree m>=2, and

    u=h(y), E=C(u,s), s=x/h'(y)+beta(y), beta in C(y).

Assume the ACTUAL normal affine source-image algebra C=R intersect E
is smooth and X=Spec R -> V=Spec C is SURJECTIVE ETALE. This is MR0's
conclusion for an actual Keller factor; it is also verified directly
for the family below. The rational two-form

    Omega=ds wedge du

pulls back to dx wedge dy. By etale local descent it is a regular unit
volume on V. Any unit on V is scalar, since it pulls back to a unit of
the WHOLE source plane. Let Pi=(1/m)*Tr_(L/E). Its restriction to C(y) is the
normalized trace for C(y)/C(h), and put

    gamma(u)=Pi(beta), beta0=beta-gamma(u).

Then Omega has a GLOBAL REGULAR Kahler primitive on V if and only if
there exist G in C(y), F0 in C[x,y] satisfying BOTH

    Tr_(C(y)/C(u))(G)=0,
    dG=beta0*du,
    F0-G belongs to E.                                    (1)

The last condition is essential: G must have a POLYNOMIAL SOURCE
representative modulo E. Rational exactness alone is not claimed to
make a primitive regular on V.

### Necessity by actual field trace

Let theta be a global regular primitive. Its source pullback differs
from x*dy by a closed polynomial one-form, so

    theta=x*dy+dF0, F0 in R.                              (2)

The polynomial Poincare assertion can be proved by integrating one
coefficient in x and then the remaining coefficient in y. Trace commutes
with the differential in a finite separable extension, as in the retained
KO2 proof: its embeddings commute with the uniquely extended derivations.
Since theta is a form over E and x*dy=(s-beta)*du, normalized trace gives

    theta=(s-gamma(u))*du+dPi(F0).

Subtract (2), obtaining

    d(F0-Pi(F0))=(beta-gamma(u))*du.                       (3)

The right side has no dx coefficient. Its rational primitive in L
therefore belongs to C(y). Set G=F0-Pi(F0). It has trace zero, and
F0-G=Pi(F0) belongs to E. This proves ALL parts of (1) globally.

### Sufficiency, including source-image regularity

Given (1), write H=F0-G in E and define the rational E-form

    theta=(s-gamma(u))*du+dH.

On the actual source it is exactly x*dy+dF0, a regular polynomial
one-form. Surjective etale descent makes theta regular at EVERY point
of V, not just at its generic divisors: the pullback of its locally free
one-form module is the source one-form module. Its derivative is Omega.
This proves sufficiency, without a chart atlas or a formal-to-regular
extension assumption.

The first condition in (1) has a cheap necessary residue formulation:

    beta0(y)*h'(y)*dy must have zero finite residues.       (4)

Conversely zero residues supply a rational G; its constant can be chosen
to make its trace zero, since dTr(G)=0. They do NOT supply the remaining
polynomial representative F0 in (1). The obstruction is invariant under
s->s+q(u), because beta0 is unchanged. A polynomial source shear changes
the differential in (4) by an exact rational differential, so cannot
repair a nonzero residue.

An actual etale target V->A2 supplies theta=P*dQ (up to a scalar), and
must therefore pass this criterion. This does not say that an arbitrary
regular primitive supplies two global target functions with unit Jacobian.
The permitted rational slice has not been set to zero.

## 2. Explicit source-image algebras in every polynomial degree

Let h be ANY monic polynomial of degree m>=2. Choose a polynomial A(u)
whose roots include EVERY finite critical value of h, with EVERY root
of A having multiplicity at least TWO. Additional roots are permitted.
Set

    S_(h,A): A(u)*w=h(v)-u.                               (5)

This is a smooth integral AFFINE surface. Off A=0 use the w derivative
of its defining equation. Over A=0, A'=0, so its u derivative is1,
proving smoothness there as well. Integrality follows by primitive
linear-in-w irreducibility: A(u) and h(v)-u are relatively prime.

Its regular unit volume is

    Omega=dv wedge du/A(u).                               (6)

For example near A=0 the differential of (5) gives
(1+A'*w)du=h'(v)dv-A*dw, so (6) equals
-dv wedge dw/(1+A'*w), a regular unit there. Off A=0 it is already
regular. This is a hypersurface volume, not an asserted target pair.

Define the polynomial source map explicitly by

    u=h(y), c(y)=A(h(y))/h'(y),
    v=y+c(y)*x,
    w=[h(y+c(y)*x)-h(y)]/A(h(y)).                         (7)

These are SOURCE POLYNOMIALS. Indeed h'^2 divides A(h): at a critical
point of local degree n the denominator order is2(n-1), while the
numerator has order at least2n. Taylor's formula gives

    w=x+sum_(j=2)^m [h^(j)(y)/j!]
                *[A(h)^(j-1)/h'^j]*x^j,
    A(h)^(j-1)/h'^j=[A(h)/h'^2]^(j-1)*h'^(j-2).           (8)

Thus (8) proves polynomiality uniformly, without expanding a growing
source-degree ansatz. Pullback of (6) is dx wedge dy, so the map of
smooth surfaces is ETALE everywhere.

It is SURJECTIVE. At a root b of A and any root c of h(c)=b,
the line y=c maps to

    (u,v,w)=(b,c,x),

isomorphically onto that WHOLE fibre component. All the higher terms
in (8) vanish there. Off the roots of A, every h-fibre has m simple
roots, and each source line maps by a nonconstant affine-linear change
of x onto the target A1-fibre. Its generic degree is m.

The source field realizes the permitted rational slice

    E=C(u,v)=C(u,s), s=v/A(u)=x/h'+y/A(h),
    beta(y)=y/A(h(y)).                                    (9)

Regular functions descend along this SURJECTIVE ETALE map, so the
EXACT source-image algebra is

    R intersect E = O(S_(h,A)).                           (10)

One may check this also at every height-one prime and then use normality.
Units on S are scalar: a unit pulls back to a unit of the whole plane,
and dominance makes the same scalar on S. These are genuine explicit
source-image algebras, not merely abstract fibre passports. They are NOT
Keller maps to A2: their target is the surface (5).

For the user's first model choose

    h=y^2*(y-1)^3, a=-108/3125,
    A(u)=[u*(u-a)]^2.                                    (11)

Its u-fibres are EXACTLY 2C2+3C3 at0 and one double plus THREE distinct
unit components at a. The source map is a five-sheeted non-Galois
surjective etale affine pseudo-covering. Its relative monodromy is S5,
as in MU0. Thus the mixed source-cover geometry itself is realizable.
Its global etale plane TARGET is ruled out next.

## 3. An all-degree rational-map obstruction to exact volume

Suppose A has r>=2 distinct roots b_j, of multiplicities k_j>=2,
and N=deg A=sum k_j. Let mu be the average of the m roots of h(Y)-u;
it is constant, by the coefficient of Y^(m-1). For (9),

    gamma(u)=mu/A(u),
    beta0*du=(y-mu)*h'(y)*dy/A(h(y)).                     (12)

We prove this rational differential has NO rational primitive, in EVERY
degree. This immediately fails the necessary part of (1).

Assume G is a rational primitive. At a point c with h(c)=b_j of local
degree n_c, put epsilon_c=1 if c=mu and0 otherwise. Its derivative has
pole order

    n_c*(k_j-1)+1-epsilon_c.

A simple pole is already impossible for a rational derivative. Otherwise
G has pole order n_c*(k_j-1)-epsilon_c. No other finite pole is possible.
The total pole degree of the rational map G:P1->P1 is consequently

    deg G=m*(N-r)-epsilon,
    epsilon=1 if A(h(mu))=0, and0 otherwise.               (13)

At infinity its derivative has leading term a nonzero constant times
y^(-m*(N-1)). Hence G has a finite limit there and

    ord_infinity(G-G(infinity))=m*(N-1)-1.                 (14)

A local degree cannot exceed the total degree of a rational map.
Equations(13)--(14) give

    m*(r-1)<=1-epsilon,

impossible for m>=2,r>=2. This adapts PF137's retained elementary
pole-degree argument to a trace-free ALGEBRAIC base pullback. It is
not a bounded residue computation or a transfer of its separated-variable
hypersurface hypothesis.

Thus Omega is NOT EXACT on S_(h,A) whenever A has at least two distinct
roots. Since its units are scalar, ANY etale morphism S_(h,A)->A2 would
pull back the plane volume to a scalar multiple of Omega and make it
exact. Therefore

    NO etale plane target exists for this ENTIRE family.  (15)

In particular (11) is excluded for arbitrary proposed target degree.
More generally every h with two or more finite critical values, and
EVERY A as in Section2, is excluded. This includes the whole 2+odd-N
sequence, but does not exclude other rational beta for those h.

The root-count condition is real. For h=y^m, A=u^k with m,k>=2,
the surface u^k*w=v^m-u has the regular Euler field of weights
(u,v,w)=(m,1,m*(1-k)). The volume (6) has nonzero weight
1-m*(k-1). Cartan's formula gives the GLOBAL regular primitive
i_Euler(Omega)/(1-m*(k-1)). These are credited one-root pseudo-plane
controls, not newly realized Keller targets. They prevent turning
(15) into a claim excluding every source-image cover.

## 4. What this closes and what still has to be proved

The natural graph models (5)--(11) realize all the source-cover features
of the first mixed fibre and retain a genuinely rational slice, yet fail
the GLOBAL regular primitive required by a target coframe. This excludes
an explicit all-degree realization architecture. It also shows why the
source-image cover alone cannot supply exactness or a Keller target.

For an arbitrary permitted beta, the correct next constraint is the full
trace-free criterion (1): a rational primitive G and a polynomial source
representative modulo E. The first is now a precise residue obstruction;
the second is a GLOBAL algebra condition. Even passing both supplies
exact volume, not a target pair or full finite-B trace-lattice regularity.
PF138 already cautions against promoting exactness of a given atlas to
target realization; its original scope and completed tests are retained.

MU0's forced marked quadratic trace pole remains valid on any actual
mixed factor. The family here has a real surjective plane-source cover
but no original target B, so it is not a counterexample to HC4/JC2 or
automatic trace regularity in the actual Keller setting. No field-trace
functional is substituted for a Mathieu ordinary marked-period source.

MR0/MU0, PF137/PF138 and KO2 are credited at their exact scopes. The
source primitive, rational trace, explicit polynomial map and rational-map
degree proofs are reproduced above. Primary pseudo-covering literature
was screened for context, not adopted as an additional classification.
No new external theorem, primary cache/download or registered mathematical
checker replay is used. New controls verify transcriptions, not independent
whole-proof review. Full HC4/JC2 and its original full trace/marked-energy
implication remain OPEN; published quartic PLANE exclusion and deferred
external review/outreach remain intact.
