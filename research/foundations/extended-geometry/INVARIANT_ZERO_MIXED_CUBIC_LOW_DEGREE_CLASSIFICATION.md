# Invariant-zero mixed cubics: low-degree safety and the central classification

Work over any characteristic-zero field K and its coefficient extensions in
A=K[x,y,t]/(t^2+2xy-1), with
Lambda=CT_x (1/2)integral[-1,1] after y=(1-t^2)/(2x).
A phase is pure when all positive pure moments vanish, and safe when
every fixed polynomial marker eventually has zero mixed moments.

Let

    f=x^-1 P,   P=C+tB+t^2D+t^3E,   C,B,D,E in K[x],
    E0!=0,  C0=D0=0,  B0=-E0,
    27CE^2-9BDE+2D^3=0.                               (A1)

These are the CENTRAL unit-E depressed-cubic hypotheses, not the
zero-center condition that all four constant coefficients vanish.
They are not asserted necessary for general mixed-cubic purity.
The central classification is proved first. Section8 then gives
degree5/7 safety in the ENTIRE invariant-zero mixed cubic branch,
including E0=0 and the off-center unit cases. Section9 gives a
separate automatic nonpower lemma for exact height3 phases.

**Theorem CentralDepressedCubicThroughDegreeSeven.**
Every pure phase(A1) of least degree at most7 belongs to exactly one
of the following types.

1. For b in K*, and a polynomial c in xK[x] of degree at most2,

       E=-b,   P=b[(t-c)-(t-c)^3].                     (A2)

   These phases are safe. Their exact least degrees are2 or5.

2. For b,a in K*, U=1+a*x, and constants p0,p1 in K,

       c=p0+p1*x-p0/U,
       E=-b U^4,
       P=b U^2(t-c)-b U^4(t-c)^3.                     (A3)

   These phases have exact least degree6 and are unsafe. The SAME
   fixed sphere marker x survives at every positive EVEN order.

Every displayed choice is regular and pure. Thus there is no pure
central-unit phase of exact least degree7; all pure central-unit
quintics are safe. This theorem neither classifies the generic
nonzero-invariant branch nor settles the global degree5/7 question.

## 1. Literal degree and polynomial depression

For every regular mixed cubic presentation, the exact least degree N is

    N=max(degC-1,degB,degD+1,degE+2),                   (A4)

with deg0=-infinity. Indeed the regular polynomial representative is

    2C0*y+2B0*yt+(C-C0)/x+t(B-B0)/x
       +t^2(D+C0)/x+t^3(E+B0)/x.

It has no xy monomial and is reduced under a degree-compatible order
with leading relation monomial xy. A nonzero reduced top homogeneous
part cannot be divisible by t^2+2xy; its leading monomial would then
have an xy factor. This proves(A4) without a filtration guess.

Define

    c=-D/(3E),   H=B-D^2/(3E).                         (A5)

The rational function c is a formal power series at zero with c0=0.
The crucial fact is that H is a POLYNOMIAL.

At a geometric E-root of multiplicity e, suppose D has order d with
2d<e. In(A1), the term D^3 then has order3d, strictly smaller than
both ord(BDE)>=d+e and ord(CE^2)>=2e. Its nonzero coefficient2
prevents the three terms summing to zero. Therefore2d>=e at every
E-root, and E divides D^2. The argument includes D=0 trivially.
Hence H lies in K[x], with H0=B0=-E0 nonzero.

The invariant equation is exactly

    P=E(t-c)^3+H(t-c),   C=-E c^3-Hc.                  (A6)

No constant term remains after depressing the cubic. Put

    g=t(H+Et^2)/x.                                    (A7)

It is a regular unit-center odd cubic.

## 2. Depression cannot increase the least phase degree

Let N=deg_A f, e=degE and delta=degD. Formula(A4) gives

    e<=N-2,  delta<=N-1,  degB<=N,  degC<=N+1.

If degH>N, its highest term must come from -D^2/(3E), and

    degH=2delta-e>N.

This forces j=delta-e>=1, since otherwise2delta-e<=delta<=N-1.
In the invariant expression

    C=BD/(3E)-2D^3/(27E^2),

the second term has degree3delta-2e=degH+j. The first has degree
at mostN+j, strictly smaller. Thus C has degree degH+j>=N+2,
contrary to degC<=N+1. Therefore

    degH<=N,    deg_A g=max(degH,degE+2)<=N.             (A8)

This is a least-degree bound, not only a bound on coefficient sums.

## 3. Exact pure and polynomial-x marked moment preservation

Let G(x,t)=Et^3+Ht. By(A6), P(x,t)=G(x,t-c(x)).
The coefficients are formal power series at zero, while P and G
are polynomials in x,t. At x=0,

    G(0,t)=E0(t^3-t),

which vanishes at t=1 and t=-1.

For every m>=1, translate the numerator interval by c=O(x):

    integral[-1,1] P(x,t)^m dt
       =integral[-1-c,1-c] G(x,z)^m dz.

The difference from the interval[-1,1] is a sum of two endpoint
integrals of lengthO(x). Near either endpoint z=epsilon+O(x),
epsilon=+1 or-1, G(x,z)=O(x). Equivalently G belongs to the ideal
(x,z-epsilon) in the corresponding formal local ring.
Its mth power integrated over that endpoint displacement is
O(x^(m+1)). Thus the formal polynomial averages differ by
a multiple of x^(m+1).

Extracting the coefficient at x^(m-ell) proves the EXACT identity

    Lambda(x^ell f^m)=Lambda(x^ell g^m),
       m>=1, ell>=0.                                 (A9)

In particular f is pure iff g is pure. The assertion extends by
linearity to every polynomial in x as marker.

This is NOT a safety-transfer theorem for a general rational shift.
An arbitrary fixed marker can have negative angular support, and
there is no claim that its rationally translated expression is a
polynomial marker in A. Safety will be transferred only when c/x
is a polynomial.

## 4. Exact retained odd-cubic inputs

The following scopes are the actual internally retained written
theorems, applied to g in(A7), with old coefficients named H,E.

- ODD_CUBIC_HEIGHT_LOW_DEGREE_CLASSIFICATION.md, main theorem:
  a pure regular unit-center odd cubic of least degree at most5
  has H=b,E=-b constant. Its zero-center branch is irrelevant here,
  since H0=-E0 is nonzero.
- The same source, OddCubicRationalCenter, gives the rational map

      u=x sqrt(-E/H)/H,  u=x/b+O(x^2),
      u^2=-x^2 E/H^3.                                  (A10)

- ODD_CUBIC_DEGREE_SEVEN_SAFETY.md, Section1:
  no pure unit-center odd cubic has exact least degree7.
  Its Sections2--3 give the all-degree tools used below:
  squarefree nonconstant H is impossible; exactly two geometric
  H-roots are impossible; for h distinct roots and c_H=degH>=2,
  the rational part R of a primitive of dx/H has degree c_H-h,
  and

      deg u<=2(c_H-h).                                 (A11)

  If u has a zero at infinity of orders>0, then its simple zero
  at zero gives deg u>=s+1.
  These statements use exact/log separation and the induced
  intermediate-field involution; it need not extend to K(x).
- ONE_POLE_ODD_CUBIC_CLASSIFICATION.md, Theorem
  OnePoleOddCubicPurity:
  constant H forces E constant; for one geometric H-root and
  c_H>=2, the possible pure least phase degrees are3c_H or5c_H-2.
  Its explicit branches, not just these degree bounds, will
  determine the c_H=2 case.

Whenever a cited source normalizes H0=1, apply it to g/b,
b=H0. Nonzero scalar rescaling preserves purity, safety, coefficient
degrees, root counts, map degrees and infinity orders.

These are all-order results at their stated hypotheses, not
conclusions from a finite vanishing prefix. No prime-existence
principal or coarse old degree family is used as a premise.

## 4A. The exact/log separation on a complete common curve

**Lemma CentralOddCubicExactLogSeparation.** For a pure unit-center
odd cubic with rational coordinateu as in(A10), write

    dx/H=dR+sum_i a_i dlog(x-r_i),

with R rational, normalized at zero. Then

    R(v(u))-R(v(-u))=2u.                               (A10a)

Here v is the specified inverse branch of u at zero. This supplies
the exact/log separation used in the retained degree tools.

Indeed the all-order endpoint criterion gives
V(v(u))-V(v(-u))=2u, where V'=1/H. Take a finite common algebraic
function field over the algebraic closure k of K containing both
inverse branches, and its nonsingular complete curve. Differentiation
of the formal identity gives

    dF=-sum_i a_i dlog((v(u)-r_i)/(v(-u)-r_i)),
    F=R(v(u))-R(v(-u))-2u.

Every logarithmic differential has at most simple poles. If the
rational function F had a pole of positive orderr, dF would have a
pole of orderr+1>=2, with nonzero leading coefficient in characteristic
zero. Therefore F has no poles on the complete curve and is constant.
Its expansion at the specified base point gives F(0)=0. This proves
(A10a). Equivalently, a nonconstant function has a pole at a place
extending the infinity place of its rational function field, so the
no-pole assertion needs no analytic comparison or period theorem.

For precision the root tools can now be recovered algebraically.
If H is squarefree and nonconstant, R is constant, contradicting
(A10a). If H has c_H>=2 and h distinct roots, its repeated-root
principal parts give degR=c_H-h and no infinity pole. On
k(u,R), formal u->-u induces the automorphism

    u->-u,   R->R-2u.

It need not extend to k(x). The degree tower and this automorphism
give deg(R-2u)=degR. Thus
deg u<=degR+deg(R-2u)=2(c_H-h), as used in(A11).
A zero of orders>0 at infinity, together with the simple zero at
zero, gives deg u>=s+1.

Exactly two roots of H are impossible. Their residues sum to zero
because degH>=2; they cannot both vanish, since then R'=1/H would
make degR=c_H-2 while R-R(infinity) has an infinity zero of order
c_H-1, impossible for a rational map of that degree. Hence the
logarithmic part is a nonzero scalar times
dlog((x-r1)/(x-r2)). Subtracting the exact-part identity forces
the ratio of that Mobius function at v(u),v(-u) to be constant1.
Mobius injectivity would give v(u)=v(-u), contrary to their opposite
nonzero linear coefficients. These arguments justify the exact/log
and root exclusions without importing an unstated separation lemma.

## 5. Unique pure unit-center odd cubic of degree six

Suppose g is pure and has exact least degree6. Write
c_H=degH,e=degE. Rationality(A10) implies c_H and e have the
same parity. The exact degree max(c_H,e+2)=6 therefore leaves

    c_H=0,2,4 with e=4,
    or c_H=6 with e=0,2,4.                            (A12)

If c_H=0, the constant-H case of OnePoleOddCubicPurity makes E
constant, so g has degree2, a contradiction.

For c_H=4 or6, a single geometric H-root is excluded by its possible
degrees3c_H and5c_H-2, neither of which is6. Exactly two roots are
excluded by the retained two-root argument. If h>=3, use(A10):
u vanishes at infinity with order

    s=(3c_H-e-2)/2.

At c_H=4,e=4 one has s=3 and
deg u>=4, while(A11) gives deg u<=2(4-h)<=2.
At c_H=6,e<=4 one has s>=6 and
deg u>=7, while(A11) gives deg u<=2(6-h)<=6.
Both are contradictions. Squarefree cases, where R is constant,
are already excluded before using the nonconstant-map bound.

Thus c_H=2,e=4. Two distinct H-roots would be squarefree,
so H has a single double root. Restore H0=b and write

    H=b U^2,  U=1+a*x,  a!=0.                          (A13)

The one root and coefficient a belong to K; explicitly
a=H'(0)/(2b). After dividing the phase by b,
OnePoleOddCubicPurity has N=1. Its inversion branch gives

    E/b=-U^4(U+1)^2/4,

of degree6, hence phase degree8, and is excluded.
The unique fixed-zero branch has divisor d=1 and quotientq=1:

    ((2U-1)-1)/(U-1)=2,
    E=-b U^4.                                        (A14)

This proves the claimed degree-six uniqueness, over K itself.

For clarity its sufficiency and unsafe marker can also be seen
directly. Under v=Ut, even numerator moments have the endpoint
congruence b^m kappa_m U^(m-1), while odd moments vanish by parity.
Consequently

    Lambda(g^m)=0 for all m>=1,
    Lambda(xg^m)=b^m a^(m-1) kappa_m !=0 for every
      positive EVEN m,                               (A15)

where kappa_m=integral[0,1]v^m(1-v^2)^m dv.
Both b,a are nonzero and every even kappa_m is a nonzero rational.

## 6. Returning to the mixed phase

By(A8), a pure f with degree at most7 gives a pure unit-center
odd g with degree at most7. The low-degree theorem, exact-degree7
exclusion and Section5 leave only g of degree2 or6.

If g has degree2, H=b,E=-b are constants. Then c=-D/(3E) is a
polynomial with c0=0. The original f is the genuine polynomial
Hopf shear t->t-c(x) of g=2bty. Its coefficient h=-c/x lies in
K[x], so direct finite-exponential moment invariance transports
every fixed marker. Thus f is safe.

For completeness, invariance follows without an all-Aut theorem:
with J=x partial_t-t partial_y and h in K[x],

    Lambda(hJ a)=CT_x xh
       [a(x,0,1)-a(x,0,-1)]/2=0.

The boundary bracket is polynomial. Exponentiating the locally
nilpotent derivation hJ preserves Lambda, and its inverse transports
arbitrary fixed polynomial markers.

If c=0 or degc=1, formula(A4) gives degree2. If degc=j>=2,
then C=b(c^3-c), B=b(1-3c^2), D=3bc, E=-b give exact
least degree3j-1. The bound at most7 thus forces degc<=2,
yielding exactly type(A2), of degree2 or5.

If g has degree6, use(A13)--(A14). The only possible finite
pole of c=-D/(3E) is the root of U. If its order isr>=3,
polynomiality of B=H+3Ec^2 fails, since that term has order4-2r<0
while H has order2. Thus r<=2. If r=2, polynomiality of
C=-Ec^3-Hc fails: the first term has order-2 and the second
has order0, so cancellation is impossible. Therefore r<=1 and

    c=p+q/U,  p in K[x],q in K.                        (A16)

If p is bounded, all of C,B,D have degree at most4, and E has
degree4, giving exact phase degree6. If degp=j>=1, the degrees
are exactly

    degC=4+3j,  degB=4+2j,  degD=4+j,  degE=4.

Thus(A4) gives exact least degree

    max(3j+3,2j+4,j+5,6)=3j+3.

The value is6 for j1 and at least9 for j>=2.
Hence degree at most7 forces degp<=1, and c0=0 makes
q=-p0. This is exactly type(A3), always of degree6.

Conversely every type(A3) has polynomial C,B,D,E:
the denominators U,U^2,U^3 are removed by the displayed U^4
and U^2 factors. It is pure by(A9) and(A15), and unsafe by the
same fixed-x even moments(A9). Every type(A2) is the safe polynomial
shear just proved. This finishes the exhaustive classification
through least degree7.

## 7. Scope of the reduction

Combining this central theorem with the separately proved off-center
I0 branch supplies safety at least degrees5 and7 throughout the
unit-E invariant-zero mixed cubic class. The off-center source is
MIXED_CUBIC_ZERO_CENTER_AND_CONDITIONAL_OFFCENTER.md, Sections3--4,
at its explicit C0=+-3B0 hypotheses; no generic constant-shape theorem
is needed for this combination. Algebraically the three cases are
exhaustive because(A1) at x0 gives
2C0(9B0^2-C0^2)=0 and B0=-E0!=0.

This does not classify generic I!=0 or arbitrary higher J-height.
Section8 treats the E0=0 invariant-zero case; an E0=0 nonzero-center
phase necessarily lies outside I0. No unrestricted global degree5/7
Image conclusion follows.

The old odd-cubic theorems are used only at their literal scopes
listed in Section4. Full-field rational depression, the degree
bound, exact polynomial-x marked moments and the mixed low-degree
conclusion are proved explicitly above. No external review,
formal verification or historical priority is asserted.


## 8. The entire invariant-zero branch at degrees five and seven

Now remove the unit-E and central assumptions. Retain the regular
mixed cubic presentation, the two regularity conditions, and I=0.

**Lemma VanishingEConstantInvariantZeroIsZeroCenter.**
If E0=0, regularity gives B0=0. Evaluating I at zero then gives
2D0^3=0. Characteristic zero forces D0=0 and then C0=0.
Thus all four center coefficients vanish.

For a pure such phase, the algebraic first-four-moment theorem in
MIXED_CUBIC_ZERO_CENTER_AND_CONDITIONAL_OFFCENTER.md, Section2,
forces its weight-zero polynomial, of degree at most3, to vanish.
Consequently C,B,D,E are all divisible by x^2. The phase has
strictly positive angular weight and is safe for every fixed marker.
This uses ONLY that finite algebraic theorem and bounded Laurent
support, not the classical all-height interval-rigidity theorem.

If E0 is nonzero, B0=-E0 is nonzero. The invariant at zero is

    I(0)=2C0(9B0^2-C0^2).

Thus exactly the three cases C0=0, C0=3B0, C0=-3B0 occur.
The first is the central theorem above. In the other two cases,
the literal off-center invariant-zero result of
MIXED_CUBIC_ZERO_CENTER_AND_CONDITIONAL_OFFCENTER.md, Sections3--4,
gives safety of every pure phase of least degree5 or7.
Its extra-invariant and unit/offcenter hypotheses have just been
checked; its classical Section1 interval theorem is not imported.

**Corollary InvariantZeroMixedCubicDegreeFiveSevenSafety.**
Every pure regular mixed cubic with I=0, of exact least degree5 or7,
is safe over every characteristic-zero field. The central unit
degree7 case is actually impossible; the assertion above needs only
safety in the remaining cases. No necessity of I=0 is asserted.

## 9. Automatic scalar nonpower in exact height three

The following is independent of the invariant-zero hypothesis.

**Lemma HeightThreeNegativeCenterAbsoluteNonpower.**
Let f be a regular mixed cubic of exact J-height3 and suppose
(C0,B0) is not(0,0). Then f cannot equal alpha*h^r in A_L for ANY
coefficient-field extension L/K, r>=2 and alpha in L*.

The Laurent chart embeds A_L in L[x,x^-1,t], with J=x partial_t,
so derivation height is precisely degree in t. Such a factorization
would give3=r*deg_t(h). Hence r3 and deg_t(h)=1.

Every regular height-one h lies in L[x,t]. To see this, write
x^N h=a(x)t+b(x) with N a minimal nonnegative common pole order,
a,b polynomials and at least one of a0,b0 nonzero if N>0.
For N>0, reduction in
A_L/(x)=L[y,t]/(t^2-1) would force a0=b0=0 by independence of1,t,
a contradiction. Thus N0 and h is polynomial in x,t.

It would follow that f=alpha*h^3 has no negative angular x power.
But f's x^-1 coefficient is exactly
(1-t^2)(C0+B0t), nonzero by the hypothesis. This contradiction
proves the lemma. It is an absolute scalar nonpower statement,
with no field-embedding or safety-after-extension assumption.

**Corollary PureUnsafeHeightThreeAbsoluteNonpower.**
Every PURE UNSAFE regular exact-height3 mixed cubic is an absolute
scalar nonpower. Indeed if C0=B0=0, regularity also gives D0=E0=0,
and the algebraic first-four theorem above makes the pure phase
strictly positive in angular weight and safe. Thus an unsafe phase
has the nonzero negative center required by the lemma.

There is also a direct original-ring consequence. If a balanced
homogeneous original polynomial F restricts to a phase satisfying
the lemma, then F is not a scalar proper power in the FULL original
ring over extensions. For if F=alpha*H^r, the smallest/largest dual
degrees in H multiply by r; homogeneity of F forces H homogeneous
in the dual variables. Repeating independently in the coordinate
variables makes H balanced. Its restriction would supply the
excluded power of f in A_L. This is the forward grading argument,
not a quotient-only reverse assertion. No starting unsafe witness
or original degree5/7 existence is supplied by the lemma.

The exact central controls and source-guarded replay for this note
are supporting algebra, not a proof of the old rational-u theorem
from a finite moment prefix. All imported odd-cubic statements,
their exact finite certificate and the two narrow mixed425 imports
are authenticated at their stated scopes.

## 10. Source and exact-control boundary

The source ledger authenticates the three retained odd-cubic notes
listed in Section4, the ACTUAL
odd_cubic_height_low_degree_certificates.json used by the old
degree-at-most-five theorem, and the mixed425 note only at its
algebraic first-four/offcenter locators. The old certificate's
primitive generators, moment scalings and four exact Bezout
identities are replayed using rational arithmetic.

For the finite cyclic intermediate-field step in the imported
one-pole classification, the retained primary is J.S. Milne,
Fields and Galois Theory v5.10 (September2022), Theorem3.17,
printedp39, in sources/milne_fields_galois_v510.pdf. That is the
finite Galois correspondence, not a new automorphism theorem.
The norm-product Remark5.47 printedp84 is context for other old
branches; it is not required by the degree-six one-pole cases here.

The checker verify_invariant_zero_mixed_cubic.py authenticates all
these reviewed raw bodies, replays the retained exact certificate,
and checks the new universal depression identities, the six
degree-six boxes, raw polynomial-x moment comparisons, and
necessary original-grading/pole boundary controls. Its finite
controls do not infer the old all-order rational-u theorem or the
new classification from a vanishing prefix. The proof above uses
no full constant-affine-shape theorem, no FPYZ interval theorem,
and no general all-automorphism theorem.

```bash
python3 research/foundations/extended-geometry/verify_invariant_zero_mixed_cubic.py \
  --repo-root . --output /tmp/invariant-zero-mixed-cubic.json
```
