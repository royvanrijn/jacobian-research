# Logarithmic descent and cyclic exact source images with arbitrary Picard rank

HC4MYMX0. 16 November 2026. Written all-degree arguments over C.
Full HC4/JC2 and ACTUAL source => Tr(Lambda_L*Lambda_L) subset B remain OPEN.
No novelty, independent replay, formal verification or external review is asserted.

MW0 completed additive primitive descent. This note tests its proposed
multiplicative upgrade and identifies the obstruction exactly. An explicit
family retains a SURJECTIVE ETALE whole-plane source, cyclic relative field,
exact unit volume and arbitrarily large free Picard group. Its plane TARGET
is not constructed. Thus the obstruction survives on the cyclic geometry
left by PF2RULED142; this is not a claim about the actual finite trace lattice.

## 1. The multiplicative differential quotient is exactly Pic(V)

Let R=C[x,y], L=C(x,y). Assume q:A2->V=Spec(C) is SURJECTIVE ETALE,
V smooth integral affine, with function field E subset L. The map need
not be finite. Define the multiplicative subgroup

    M={g in L* : dg/g belongs to Omega_R^1+Omega_E^1}.

Then there is a natural isomorphism, up to the harmless choice of sign
for the Picard class,

    M/E* = Pic(V).                                           (1)

All forms and fields are embedded in Omega_L^1 as in MW0. We use Pic(A2)=0
and R*=C*, so the denominator E* already contains every regular source unit.

For g in M let z=g2/g1 on the ACTUAL smooth relation T=X times_V X.
Its logarithmic differential is regular. At the generic point of a divisor,
write z=t^a*b with b a unit. If a!=0 then dz/z=a*dt/t+db/b has a pole in
characteristic zero. Thus z has no zero or pole. Normality makes both z
and its inverse regular at ALL points of T. Moreover z23*z12=z13 on the
triple relation. It is a unit cocycle defining a line bundle on V trivial
on X. This uses effective faithfully flat descent for line bundles, the
standard input already retained at
[Stacks tag023R](https://stacks.math.columbia.edu/tag/023R).

The class vanishes exactly when z=r2/r1 for some r in R*. Then g/r has
equal rational pullbacks and belongs to E*, proving injectivity of (1).

Conversely, take any line bundle on V. Its pullback to X is trivial.
Choose a rational section over E and a regular trivialization after
pullback; their ratio is a g in L*. The corresponding z=g2/g1 is a regular
unit cocycle on T. The regular differential dz/z is an ADDITIVE cocycle
with values in the pullback of Omega_V^1. Here etaleness identifies
Omega_X^1=R tensor_C Omega_V^1 and the analogous relation modules.
Amitsur exactness for that C-module supplies a regular source one-form
alpha with alpha2-alpha1=dz/z. Therefore dg/g-alpha has equal rational
pullbacks, hence is an E-form. Thus g belongs to M, proving surjectivity.
This is the same all-module faithfully flat exactness proved and credited
in MW0; it does not require a finite map or a Galois cover.

The contrast is exact: MW0 splits regular ADDITIVE cocycles using affine
coherent cohomology. Unit cocycles instead retain Pic(V), even on an affine
surface. Equation (1) is an obstruction to this logarithmic upgrade, NOT
an identification of Pic(V) with the finite trace-pairing defect.

## 2. Requiring a closed source form leaves only torsion/Kummer classes

Let M_cl consist of g for which

    dg/g=alpha+eta, alpha a CLOSED regular source one-form,
    eta a rational E-form.

Then, with m=[L:E],

    M_cl={g in L* : g^m in E*},
    M_cl/E* = Pic(V)_tors.                                (2)

In particular every torsion order divides m. This is not an assertion
that every g in M admits a closed alpha.

Polynomial Poincare gives alpha=dF0, F0 in R. Let N be the field norm
and Pi the normalized field trace. Taking trace of the logarithmic
identity gives

    dlog(g^m/N(g))=d(m*(F0-Pi(F0))).                       (3)

The elementary rational exponential lemma says that dlog(H)=dA in a
characteristic-zero function field forces A,H constant. To prove it,
take a smooth projective model: a pole of nonconstant rational A gives
dA a pole of order at least two, while a rational logarithmic differential
has at most a simple pole at a generic divisor. With no poles A is
constant on the projective model, and dH=0 makes H constant. Applied to
(3), this gives g^m/N(g) in C*, proving g^m in E*. The converse is immediate:
dlog g=(1/m)dlog(g^m) is already an E-form.

For the Picard assertion, a class represented by g is torsion of order
n precisely when z^n=r2/r1 for a source unit r. Such an r is constant,
so this is equivalent to g^n in E*. That condition implies membership
in M_cl and therefore also g^m in E*. Hence (2).

If n is the least positive exponent with g^n in E*, E(g)/E is a cyclic
Kummer extension of degree n: C contains all roots of unity, its Galois
group acts by constant root ratios, and minimality of n determines its
order. This extension is unramified over V: equivalently the torsion
line bundle and its trivialized nth power define a finite etale mu_n
torsor, which is connected by this minimality. It is an intermediate
CYCLIC QUOTIENT, not necessarily a deck subgroup of L/E when L/E is
non-Galois. No monodromy/deck conflation is made.

## 3. The actual ruled source-image Picard group

For an MR0 source image, write the finitely many nontrivial u-fibres as

    rho^-1(a)=sum_(j=1)^r_a e_(a,j)*C_(a,j).

Enlarge this finite list to include the zeros of a local-slice denominator.
Outside it the algebra is C[u,A(u)^-1][s], a UFD with units generated
by constants and the u-a. Smooth divisor localization consequently gives

    Pic(V)=direct_sum_a Z^(r_a)/Z*(e_(a,1),...,e_(a,r_a)). (4)

Indeed the omitted component classes generate, and the displayed full
fibre divisors are EXACTLY the relations coming from units of that UFD.
Thus the free rank is sum_a(r_a-1), and the torsion summand at a is
Z/gcd_j(e_(a,j)). No surface Picard classification is imported.

These classes even have polynomial source representatives for (1):
if I_(a,j) is the set of roots c of h(c)=a assigned to C_(a,j), put

    g_(a,j)(y)=product_(c in I_(a,j))(y-c).

Source etaleness makes div(g_(a,j)) exactly q^*C_(a,j), with order one
on each source line. Hence its unit relation cocycle represents that
component class. Formula (1) supplies a differential splitting for each
such SOURCE POLYNOMIAL, although its regular source term usually is
not closed.

For MU0's first mixed 2+3 image, every fibre gcd is one. Thus Pic(V) is
torsion-free; at zero its rank contribution is one, and at the other
critical value its double plus q unit components contribute q>=1 more.
There are nontrivial free logarithmic obstruction classes, but M_cl=E*.
The mixed family is now a CONDITIONAL arbitrary-degree test only:
PF2RULED142 gives a comparison Keller map of degree 2b instead of 5b,
so it cannot occur as such a factor in a MINIMAL-degree counterexample.

## 4. Exact cyclic plane-source images of unbounded Picard rank

The free obstruction in (1) is compatible even with EXACT unit volume
and the cyclic source-cover geometry that remains after PF142.
Fix ANY g>=2, k>=2, and ANY square-free polynomial B(u) of degree r>=1
with B(0)!=0. Put u=t^g and

    a(t)=t^2*B(t^g)^k,
    R0(t)=t^-1*B(t^g)^(1-k),
    sigma(t)=-B(t^g)+g*(1-k)*t^g*B'(t^g).

The uniform identity

    R0'=sigma/a                                          (5)

is immediate by differentiation. Sigma is nonzero at EVERY root of a:
at t=0 use B(0)!=0; at a root b!=0 of B use k>=2 and B'(b)!=0.

Consider the smooth integral affine surface

    W: a(t)*z=v^g-sigma(t)^g.                            (6)

Off a=0 its z derivative proves smoothness; over a=0 its v derivative
is nonzero. Irreducibility follows from primitive linearity in z.
For each epsilon in mu_g there is a WHOLE plane chart

    v=epsilon*sigma+a*x_epsilon,
    z=sum_(l=1)^g binom(g,l)*(epsilon*sigma)^(g-l)
                           *a^(l-1)*x_epsilon^l.       (7)

It selects exactly the epsilon component at every root of a. The inverse
x_epsilon=(v-epsilon*sigma)/a extends there using z divided by the other
g-1 nonzero root factors. These charts cover W and intersect in
D(a) times A1. This is a direct open embedding and covering argument,
not just a generic parametrization.

The cyclic action

    (t,v,z) -> (epsilon*t,epsilon*v,epsilon^-2*z)         (8)

preserves (6) because a has weight two and sigma weight zero. It is FREE:
away from t=0 no nonidentity element fixes t; at t=0, v^g=sigma(0)^g!=0
and none fixes v. Its translates of the epsilon=1 chart are exactly (7).
Thus V=W/mu_g is smooth affine and the chart restriction gives a
SURJECTIVE ETALE map A2->V of generic degree g. With E=C(V),

    C=R intersect E=O(V)

by regular-function descent along this actual surjective etale map.
Units on V are scalar by pullback to the whole plane.

The invariant regular unit volume

    Omega=dv wedge dt/a

pulls back to dx wedge dt. Its exactness follows already from (5) and
PF138's actual-chart cocycle argument, and can be certified directly in
the quotient using MW0. Indeed

    E=C(u,s), h=t^g,
    s=v/[g*t^(g+1)*B(u)^k]
      =x/[g*t^(g-1)]+beta(t),
    beta=R0'/[g*t^(g-1)], Pi(beta)=0, G=R0, Pi(G)=0.     (9)

Here s is invariant under (8), and beta*du=dR0. The rational slice is
retained, with poles at all the selected base roots.

There is even an explicit polynomial source representative. Since
gcd(sigma,u*B^k)=1 in C[u], choose a polynomial C0(u) with

    C0*sigma = -B modulo u*B^k,
    D0=(B+C0*sigma)/(u*B^k) in C[u].

Then the SOURCE POLYNOMIAL

    F0=t*C0(u)*x+t^(g-1)*D0(u)

satisfies F0-G=g*u*C0(u)*s in E. Thus (s*du)+d(F0-G) is a GLOBAL regular
primitive on V by surjective etale descent. This is one Bezout identity
for every g,k,B, not an increasing polynomial-degree search.

The u-ruling of V has one irreducible g-multiple fibre at zero and, at
each root of B, exactly g disjoint REDUCED A1 components. To see the
latter, the g values of t and g labels of v in W give g^2 components;
(8) acts freely on t and simultaneously rotates the labels, giving g
orbits, each quotient an A1. All other fibres are reduced A1s. Formula
(4) therefore gives exactly

    Pic(V)=Z/g direct_sum Z^(r*(g-1)),
    chi(V)=1+r*(g-1).                                  (10)

Euler follows by adding the special fibres to the ordinary ruling.
In (1) there are thus arbitrarily many independent FREE classes despite
the whole-plane surjective source, actual source-image equality, cyclic
relative field, scalar units and exact unit volume. The torsion/closed
subgroup is exactly Z/g; the class at the multiple fibre is represented
by t, with t^g=u. No plane target is supplied by (5)--(10).

## 5. Effect on the full target and source audits

MW0 is retained: the polynomial representative for an additive rational
primitive is automatic. The logarithmic upgrade instead has exactly the
Picard obstruction (1). Requiring a closed regular source form reduces
it to torsion (2), but the new exact cyclic family shows why such
closedness must not be assumed. Its free Picard classes cannot admit
that stronger splitting. Exact volume alone supplies no flat splitting
for these classes.

PF138's arbitrary-root exact moving-label controls and LA0's one-fibre
cyclic pseudo-plane controls are credited. Equations (5)--(10) combine
these mechanisms into an explicit cyclic SURJECTIVE plane-source family
with arbitrary free Picard rank; they do not realize an actual intermediate
Keller TARGET or contradict MV0's different constant-label graph exclusion.
PF142's non-power degree descent is now the current minimal-degree
reduction, independently of volume exactness and the rational slice.

The logarithmic Picard classes have NOT been identified with the full
finite trace-dual pairing or its defect. A determinant/line-bundle test
alone cannot replace its multiplication data. The remaining implication
still requires the ACTUAL target coframe together with the original full
finite order, all omitted ramified primes and marked energy. No radical
field formula is asserted to regularize that order; no Mathieu ordinary
marked-period/B-lattice comparison is supplied. Full HC4/JC2 remains OPEN.

Primary effective descent is reused at its stated scope; MW0's Amitsur
proof is reused for the module of one-forms. No compact-Kahler flat-bundle
classification or solvable-monodromy inversion theorem is imported.
New exact controls support signs and the explicit identities. The
unbounded results follow from the written arguments, not a bounded sample.
All registered predecessor checkers remain sealed and unexecuted.
