# An unramified sheet bound and a pairwise complete-curve obstruction

15 November 2026. PF2RULED48 gives two all-degree reductions for the
unramified omitted prime left in the quartic packet. The first uses the
FULL finite map. The second excludes whole source-chart configurations
by constructing a complete curve; it needs neither a positive target
normal bundle nor a finite field-degree computation.

The full quartic two-prime packet and unrestricted higher degrees remain
open. Complementary source coordinates, simple vanishing on the
unramified prime, and the explicit higher-order defect below are real
hypotheses, not consequences asserted for every normalization.

## 1. A marked unramified sheet replaces the meridian centralizer

Let pi:N->A2 be an actual plane Keller finite normalization of degree d.
Let D be an unramified omitted prime. It is A1, pi is etale at every
finite point of D, and its map onto H=pi(D) is a birational polynomial
immersion, by PF2BATT1/PF2OMIT2. Keep every other omitted prime and the
actual transitive finite monodromy group G.

Choose a smooth projective target X resolving and separating the full
branch union together with H at infinity. The immersed normalization
h:P1->Hbar is smooth and transverse to the removed boundary E. Put

    eta_H=deg N_h=-K_X.Hbar-2=Hbar^2-2*delta_H,
    delta_H=sum of finite delta invariants of H.

If eta_H>0, then

    d<=Hbar^2/eta_H=1+2*delta_H/eta_H,
    (d-1)*eta_H<=2*delta_H.                            (1)

No global symmetric-monodromy or zero-excess condition on other images
is needed. H may itself be a branch support because another omitted
prime ramifies over that same image; D's own marked sheet is unramified.

Proof. Add H to the removed target divisor if it is not already a branch
support. Its additional meridian then has identity monodromy, so the
actual quotient G is unchanged. Take the compact NORMALIZATION
neighborhood V of h as in PF2RULED47 and apply Nori WLT(A),(C):

    K=image(pi1(V')->Pi), [Pi:K]<=Hbar^2/eta_H.           (2)

This time the image L in G fixes ONE marked sheet. On the punctured
normal-disk bundle away from finite special parameters, the local etale
factor through D supplies a globally distinguished sheet. It is preserved
by the base loops because D->H^nu is an isomorphism. Normal meridians
also fix it, even if other sheets ramify over H. At any finite image
singularity or intersection, the full local cover has a degree-one factor
at the selected point of D; the whole local group fixes its letter.
Other factors can have arbitrary nonabelian monodromy on other letters.

At infinity the local complement of Hbar and E is a two-torus. The
normal H meridian fixes the sheet as above. A loop along H around its
one missing parameter also fixes it, by the marked D=A1 section on the
punctured normalization; this is the other torus generator. Thus the
whole infinity local image fixes the same letter. The same connected
collar/van-Kampen argument as PF2RULED47 proves L subset G_a, the
stabilizer of that sheet. Transitivity gives

    d=[G:G_a]<=[G:L]<=[Pi:K].                           (3)

Combining (2),(3) and adjunction proves (1). This does not restrict to
the nonfinite source open U union D. For the quartic packet, a positive
eta_H therefore requires 3*eta_H<=2*delta_H, independently of the normal
degree of the ramified image C. For eta_H<=0 no inequality is asserted.

The same argument works for any closed unramified affine-line curve in
N with birational immersive image, including a retained source coordinate
line. Such a line maps birationally: its polynomial immersion factors
through its A1 normalization, and a higher-degree polynomial residue map
would have a finite critical point, contradicting immersion. Positivity
must still be checked on the simultaneous target model.

## 2. The source-side problem must retain simultaneous coordinates

Now work just with a normal affine N containing U=A2_(x,y). Suppose two
distinct disjoint omitted affine lines R,D have volume orders e-1 and
zero respectively, e>=1, and the source volume is globally exact. Assume
the actual source coordinates x,y have these roles:

    y regular on U union R, y|R=0, v_R(y)=r>=1;
    x regular on U union D, x|D=0, v_D(x)=1.            (4)

In particular x,y are one polynomial automorphism coordinate pair in
the SAME source plane. Assume the individual opens U union R and
U union D are smooth affine. Other omitted primes may be present; they
play no role in the curve constructed below. For an actual finite Keller
normalization their affineness is the retained BGV component-open theorem;
smoothness is automatic at indices2 and1. For higher index on R it is
an explicit hypothesis. The index-two and unramified whole charts are
PF2RULED38-46; Section4 derives the needed arbitrary-index extension.

The unramified chart, retaining x rather than changing it, is

    y=v+h_D(x), h_D(x)=sum_(j=1)^p d_j*x^(-j),
    d_p!=0, d_1=0, p>=2.                              (5)

The nonnegative part is absorbed into v. Exactness gives d_1=0. No
source shear by a polynomial in y is used, since it could destroy the
role of x on D. This simultaneous-coordinate point matters in Section3.

## 3. The entire simple/simple complementary index-two sector is impossible

Suppose r=1 and e=2. The index-two chart on R has the WHOLE form

    x=y*t+A+h_R(y), h_R=sum_(j=1)^q c_j*y^(-j),
    c_q!=0, c_1=0, q>=2.                              (6)

The constant A must be retained. Gauge removes terms of degree at least
one, but removing A by a source translation would change x|D=0. Put

    C(y)=y^q*h_R(y), D(x)=x^p*h_D(x).

They have degrees at most q-2 and p-2, with C(0)=c_q and D(0)=d_p.
Define the actual nonzero source polynomial

    F=x^p*y^q-x^(p-1)*C(y)-y^(q-1)*D(x)
                         -A*c_q*x^(p-2).              (7)

Its monomials lie in the rectangle 0<=i<=p-1,0<=j<=q-1, apart from its
single top monomial x^p*y^q. At R, substituting x=A+h_R(y)+y*t gives

    y^(-q*(p-2)+1) * (c_q^(p-1)*t+constant+O(y)).       (8)

All coefficients here are regular polynomial expressions in y,t after
the displayed Laurent power is factored off. The constant correction
in (7) is what cancels the order-zero source-coordinate term A before
the t term. The coefficient of t is nonzero. At D, substituting
y=h_D(x)+v gives

    x^(-p*(q-2)) * (d_p^(q-1)*v+constant+O(x)).         (9)

Again the coefficient of v is nonzero. The implicit function theorem
therefore supplies finite boundary points on R and D for the two curve
branches at infinity. Their functions x,y are the SAME source functions,
so this is a statement about the actual embedding, not an abstract atlas.

There are no other branches at infinity. If x,y both tend to infinity,
the top monomial dominates every other monomial. If x tends to infinity
and y tends to a nonzero finite limit, it also dominates, and similarly
after exchanging x,y. The two remaining Newton faces are primitive:
their direction vectors are (1,q) and (p,1). They give exactly one
x=infinity,y=0 branch and one x=0,y=infinity branch. Equations (8),(9)
identify their finite points in the respective actual charts.

Normalize any irreducible component of F=0 in its smooth projective
completion. Its map into U extends at every point at infinity by these
charts, yielding a nonconstant morphism from a projective curve to N.
Every regular function on a projective integral complex curve is constant
(Stacks Lemma33.26.2), whereas an affine N separates distinct image
points. This is impossible. Hence (4) with r=1 cannot occur in an affine N.
The conclusion holds for ALL p,q and ALL lower Laurent coefficients,
with A arbitrary, without any finite Keller-degree assumption.
F is a polynomial in the marked SOURCE plane; it need not extend as a
regular function on N. The obstruction is the curve's complete extension,
not a claimed global principal divisor or a formal boundary identity.

## 4. Higher vanishing order and any index: an explicit characteristic defect

More generally let U union R be smooth, let e>=1, and suppose
r>=max(2,e). Keep the simultaneous source x,y. The root charts at R are

    y=z^r, x=z^(e-r)*t+h(z),
    h=sum_(j=r-e+1)^q c_j*z^(-j), c_q!=0,
    gcd(r,{j:c_j!=0})=1, c_r=0.                       (10)

Here is the whole-chart justification for arbitrary e. Normalize the
componentwise affine open in y=z^r. It is etale over R and has r boundary
copies of A1. Delete the retained line z=0 and all but one boundary copy.
As in PF2RULED40 this is an affine reduced A1 fibration over A1_z, hence
an entire plane. In the overlap x=z^k*t+h the source volume order e-1
gives r-1+k=e-1, so k=e-r. Gauge removes every exponent at least e-r,
including the constant; no source-coordinate translation is needed.
Separatedness gives the support gcd. Comparing the polynomial primitives
r*z^(r-1)*x dz and r*z^(e-1)*t dz gives residue r*c_r, hence c_r=0.
This does not infer existence of a regular coordinate or smoothness for
arbitrary higher-index primes; those hypotheses are explicit here.
Define

    j(a)=max{j:c_j!=0 and r does not divide a*j},
    S=sum_(a=1)^(r-1)(2q-j(a)), Delta=r*q-S.            (11)

This is the same full-support contact sum S as PF2RULED43, not a
truncation or a generic-coefficient estimate. Suppose additionally

    p>=r and Delta>=e.                                (12)

Then the configuration (4) is impossible in affine N.

Here is the complete polynomial construction and proof. Let

    P(x,y)=y^q*product_(z^r=y)(x-h(z))
     =(-1)^((r+1)q)*Res_z(z^r-y,z^q*x-z^q*h(z)).        (13)

It is a genuine polynomial. Its x^r coefficient is y^q; all other
coefficients have y degree at most q-1. Its value at y=0 is a nonzero
constant, since c_q!=0. Its root reference branch is primitive by the
support gcd. The coefficient of x^(r-1)*y^(q-1) is -r*c_r=0.
Thus, with

    P0=x^(p-r)*P, A0(x)=[y^(q-1)]P0,

we have deg A0<=p-2. The polynomial

    F=P0-y^(q-1)*(A0(x)+D(x))                          (14)

has leading monomial x^p*y^q, all other monomials in the same lower
rectangle as Section3, and its entire y^(q-1) column is -D(x).

At the R reference branch, the derivative of P0 has z order

    B=S+2q-pq.                                        (15)

Indeed the product of its other r-1 conjugate differences has order
-sum j(a); the norm denominator contributes rq and the prefactor
contributes -q(p-r). This is (15) by (11).
Put M=B+e-r, the order of the first chart-varying term. The correction
in (14), because deg(A0+D)<=p-2, has order at least

    r*(q-1)-q*(p-2)=M+(Delta-e).                       (16)

By (12) it cannot precede the t term. Every conjugate difference of h
has order at most -(r-e+1), whereas the variation z^(e-r)*t has strictly
higher order. Higher powers of t therefore have higher z order. The
correction's t dependence also has higher order, since q+e-r>=1.
Consequently the WHOLE substitution in (14) is

    F(h+z^(e-r)*t,z^r)=z^M*(alpha*t+beta+O(z)),
    alpha!=0.                                        (17)

The r conjugate roots supplied by the implicit function theorem descend
to the one primitive R branch. At D, the top column -D(x) gives (9),
so its one remaining branch also reaches a finite point of D.

For completeness there are no additional infinity branches. On the R
Newton face the correction in (14) has strictly smaller pole weight,
because 2q-r>0 under (10),(12). Indeed S>=(r-1)q implies Delta<=q,
so q>=e, while q>=r-e+1; adding these gives 2q>=r+1. The face of
x^(p-r)P is unchanged and
has total x-root multiplicity r. The r constructed conjugate roots
exhaust it. On the D face the y-degree jump is one and D(0)!=0, giving
one branch. No other infinity face is possible: both growing coordinates
or one growing coordinate with a nonzero finite other coordinate would
make the top monomial dominate, as in Section3. Any irreducible component
therefore extends completely inside U union R union D. The same
projective-curve contradiction proves the assertion.

If gcd(q,r)=1, then every j(a)=q, so S=(r-1)q and Delta=q. Thus ALL
coprime-pole data satisfying q>=e and p>=r are excluded in this configuration,
with no restriction on lower coefficients. This includes the whole
index-two q=r-1 family when r>=3 and p>=r. Multi-characteristic data
satisfying Delta>=e are covered too. In the index-two specialization,
the r=2,q=1 borderline and data with Delta<=1 or p<r are not uniformly
excluded by this construction.

## 5. What the pairwise theorem does and does not close

For every actual finite Keller normalization, global exactness supplies
the zero residues c_r,d_1 used here. In the remaining quartic R,D
packet, choose their intrinsic regular source coordinates from
PF2RULED45/46. If they can be chosen as the complementary pair (x,y)
and the D coordinate has simple zero, then:

* simple zero on R is entirely excluded, at all pole heights;
* for r>=2 and p>=r, the remaining data must have r*q-S<=1;
* in particular coprime q,r with q>=2 and p>=r are excluded.

The same obstruction applies pairwise in higher field degrees even if
other primes are present, since the complete curve maps into the ORIGINAL
affine N and uses only the two displayed added lines. We do not need an
affineness or finiteness theorem for the open obtained by deleting the
other primes. No common-coordinate or complementary-coordinate theorem
for arbitrary primes has been proved here. Coordinates may be the same
up to affine change, may form a noninvertible coordinate pair, or may
have higher boundary vanishing on D. Those are real remaining cases.

The polynomial in (7) retains A; independently normalizing the two source
coordinates would miss this constraint. Formula (14) retains every
Laurent coefficient via the full norm resultant and its entire top column.
Neither result is a bounded coefficient search. Failed defect controls
are preserved as failures of this sufficient construction, not asserted
affine realizations or proofs that those cases survive other obstructions.

## 6. Sources and assurance

The marked-sheet argument imports exactly Nori WLT(A),(C) and the
normalization-neighborhood construction already pinned and page-audited
in PF2RULED47. Its fixed-sheet proof is separate from the centralizer
proof: an unramified marked sheet, not the identity of the whole local
monodromy, is essential. The whole source charts and their exact residue
tests are PF2RULED38/40/43/46, with their original primary sources and
assurance limits unchanged. The complete-curve proof is a written norm,
valuation and two-face argument; the only general proper-curve input is
the constant-global-function fact linked in the source manifest.

The [source manifest](sources/unramified-sheet-pair/manifest.json) retains
the precise inputs. [verify_unramified_sheet_pair.py](verify_unramified_sheet_pair.py)
writes [unramified_sheet_pair_receipt.json](unramified_sheet_pair_receipt.json).
The checker verifies general degree identities, arbitrary-coefficient
simple-chart cancellations, the full norm constructions in representative
single/multiple-characteristic controls, the exact full-support defect,
all infinity-face bounds and failed-defect/cusp controls. Finite controls
do not prove the universal result. Formal verification, independent full
reproof, external review, novelty priority and full-goal completion remain
false. The next gate is genuinely simultaneous source-coordinate and
finite-algebra compatibility, together with the target normal-degree
constraints, not another growing Laurent-coefficient search.
