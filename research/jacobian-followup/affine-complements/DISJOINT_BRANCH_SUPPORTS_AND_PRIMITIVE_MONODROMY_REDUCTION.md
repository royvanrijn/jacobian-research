# Disjoint branch supports exclude subset and binary classical actions

PF2RULED121, 16 November 2026. For every hypothetical nonautomorphic
complex PLANE Keller map, two geometric branch meridians have nonempty,
DISJOINT moved-sheet supports. They may come from two branches of the
same irreducible target curve. This applies to the ORIGINAL finite
normalization, all omitted primes, arbitrary indices and singular N.

Consequently every faithful action induced on k-subsets of an m-set,
2<=k<=m-2, is excluded, including transitive subgroups. Combining this
with PF120 and an explicitly imported Guralnick--Magaard theorem also
excludes ALL THREE remaining binary classical primitive action families.
Primitive monodromy is therefore necessarily NATURAL A_d/S_d, or binary
affine. These are necessary containers, not realizations or an inverse
theorem. Imprimitive actions and the global finite algebra remain OPEN.

## 1. Distinct local target branches give disjoint supports

Let pi:N->A2 be the original finite normalization, U=A2 the WHOLE
Keller open, and B the reduced FULL finite branch curve. Retain every
omitted prime, including unramified primes. PF105/PF2OMIT6 give pairwise
disjoint, unibranch omitted curves, each with normalization A1. The
normalization of each irreducible component C of B is also A1. The
critical locus of pi is contained in the ramified omitted curves;
normality/purity add no isolated branch component elsewhere.

Fix y and two DISTINCT analytic branches beta_1,beta_2 of B at y.
Proper finiteness splits the inverse image of a sufficiently small ball
at y into the separate finite analytic germs at points p over y. At p
there is at most ONE omitted-curve germ. Its image, if ramified, is
ONE irreducible analytic branch of B. This remains true at a cyclic
quotient point: an irreducible critical axis has irreducible image.

Choose based local meridians gamma_1,gamma_2 around generic smooth
points of beta_1,beta_2, respectively, in the same local complement.
In each finite local factor at p, a meridian about any branch OTHER
than that factor's branch image acts trivially. Indeed that local
cover is unramified across the other branch; its meridian bounds a
transverse disc there, and basing only conjugates the identity.
Thus each local sheet block is moved by at most one gamma_i. Both
meridians act nontrivially somewhere, because both branches belong
to the full branch locus. After identifying with global generic sheets,

    supp(rho(gamma_1)) intersect supp(rho(gamma_2)) = empty.       (1)

This uses neither smoothness of N nor immersion/zero excess. It does
not say that these two elements generate G, that one chosen meridian
per component generates G, or that source fibres are nonempty.

## 2. The full branch union must have a finite multibranch point

Suppose instead every finite point of B has just one analytic branch.
Then distinct irreducible components of B are disjoint, and the finite
normalization A1->C of each component is injective. A finite bijection
of complex curves is a proper homeomorphism, so each C is contractible.

If B has one component, this contradicts PF113: an irreducible
contractible FULL branch curve is impossible, even with additional
unramified nonproper images outside B.

If B has at least two components, Chau's published no-A1-component
Remark excludes a smooth component C: every branch image belongs to
S_F, the whole nonproper-value set. Therefore pick a singular C.
The Lin--Zaidenberg theorem, in the exact Palka Theorem A form already
imported in PF113, supplies ONE polynomial target coordinate system
in which

    C={x^a-y^b=0},   a,b>=2, gcd(a,b)=1.

Let C' be another component, parametrized nonconstantly by A1. Since
C and C' are disjoint, the polynomial x(t)^a-y(t)^b has no finite
zero; it equals a constant c!=0. The smooth affine curve

    H_c={x^a-y^b=c}

has smooth projective completion of genus (a-1)(b-1)/2>=1. For
completeness, projection to x has degree b, a finite branch values
and one at infinity, each fully ramified (gcd(a,b)=1), so
Riemann--Hurwitz gives 2g-2=-2b+(a+1)(b-1). A nonconstant map
A1->H_c extends to P1 into this smooth projective curve. A second
application of Riemann--Hurwitz would give -2>=0, impossible.
Hence such C' cannot be disjoint from C.

Only C was rectified; no simultaneous coordinate rectification or
positive weighting of the entire reducible union was assumed. We
have proved that B has a finite point with at least two branches.
Equation (1) now proves the asserted universal pair of disjoint,
nonidentity BRANCH meridians. This is stronger than merely producing
two arbitrary disjoint-support elements of G.

The deep geometric inputs (Chau, Palka and the actual boundary results)
are imports/retained proofs. The deductions here are written for all
degrees, not inferred from sampled groups or finite local models.

## 3. Every interior subset action has intersecting supports

Let Sym(m) act on Omega, the k-subsets of {1,...,m}, with
2<=k<=m-2. We prove that ANY two nonidentity elements have overlapping
supports on Omega. This does not assume primitivity or G=Sym(m).

Suppose tau(i)=j!=i and put

    F_ij={A in Omega : i in A, j not in A}.

Every A in F_ij is moved by tau, since tau(A) contains j. If sigma
had disjoint support from tau, sigma would preserve EACH A in F_ij.
But the membership signatures of the m underlying letters across
F_ij are all distinct. Letter i belongs to every A, j to none, and
each other letter belongs to some but not all. Two other letters are
separated by a (k-1)-subset of the remaining m-2 letters, possible
because 1<=k-1<=m-3. Preserving every A therefore fixes every letter,
so sigma=1, a contradiction.

Together with Section2, this excludes EVERY faithful transitive
subgroup in this subset action. The endpoints k=1,m-1 are natural
actions and are NOT excluded. The exclusion needs no CFSG and no
list of possible branch-cycle types. In particular ALL interior
alternating/symmetric subset actions in PF120 are now CLOSED.

## 4. The remaining binary classical branch cycles are transvections

PF120 proves that if actual primitive monodromy has nonalternating
socle, EVERY branch meridian fixes at least d/2 sheets. This is an
all-meridian result, not just the existence of one small-support
element in PF118.

We explicitly IMPORT the following narrow part of Guralnick--Magaard,
*On the minimal degree of a primitive permutation group*, J. Algebra
207 (1998), 127--145, Theorem1, printed128--129: in its nonalternating
binary affine/classical cases an element fixing AT LEAST half the
points is a transvection, hence an involution (including the triality
equivalent orthogonal action). PF120 already excludes product actions.

Source access boundary: the primary author-uploaded original-paper
full text was read on Academia, including the theorem and introduction.
The original publisher PDF could not be retrieved and was NOT visually
audited. The source manifest pins our narrow review record and links
the primary text. No suspect OCR subset-parameter bound is used; the
primitive containers themselves come from PF120's audited modern
Burness--Guralnick classification. The prime-order-only modern bound
is NOT substituted for this ALL-ELEMENT, inclusive-half statement.
We do not prove the imported CFSG classification independently.

Thus every branch meridian in the three remaining classical actions
is a transvection. As a separate consequence, every ramified prime
has e=2,f=1, and PF106 forces the entire original N to be smooth;
this observation alone would not exclude those actions.

## 5. Binary classical transvections never have disjoint supports

Write V for a symplectic space of dimension 2a over F2, with polar
form B. Its transvection is t_v(x)=x+B(x,v)v, v!=0. In orthogonal
groups require q(v)=1. Distinct transvections can commute only when
B(v,w)=0. If B(v,w)=1 their product has order3 and they do not
commute. Disjoint-support permutations always commute. Thus only
distinct orthogonal v,w need consideration; equal transvections have
the same nonempty support.

(i) Sp_(2a)(2) on the minus-type quadratic forms, a>=3. Here
t_v fixes a quadratic form q precisely when q(v)=1, so support
intersection means q(v)=q(w)=0. Extend v,w to a symplectic basis,
with first two hyperbolic pairs v,f_1 and w,f_2. The two values
q(f_i) are arbitrary. Arf(q)=1 exactly when the remaining dimension
2a-4 quadratic form has Arf1. Consequently

    |supp(t_v) intersect supp(t_w)|
        =4*2^(a-3)*(2^(a-2)-1)
        =2^(a-1)*(2^(a-2)-1)>0.                     (2)

The standard count of minus forms in dimension 2s is
2^(s-1)(2^s-1); it follows by counting the symplectic-basis values
with sum of pairwise products1. We use it only for s=a-2>=1.

(ii) O_(2a)^+(2) on nonsingular vectors, a>=4, and O_(2a)^-(2)
on nonzero singular vectors, a>=3. Put sigma=+1 or -1 for the type,
and epsilon=1 or0 for the required q-value. A vector is moved by
both t_v,t_w exactly when B(x,v)=B(x,w)=1 and q(x)=epsilon.
In a symplectic basis extending the orthogonal nonsingular v,w,
these constraints fix the coefficients of f_1,f_2 to1. The v,w
coefficients are free and cancel out of q(x), since q(v)=q(w)=1.
Counting the remaining quadratic values, or expanding their Gauss
sum, gives

    |supp(t_v) intersect supp(t_w)|
       =2^(2a-3)+sigma*(-1)^epsilon*2^(a-1)
       =2^(a-1)*(2^(a-2)-1)>0.                     (3)

In both retained rows sigma*(-1)^epsilon=-1. The constraint already
excludes the zero vector in the singular action. The quadratic count
uses sum_x(-1)^q(x)=sigma*2^a, proved by multiplying the two-variable
Gauss sums of symplectic pairs. Triality just permutes the equivalent
permutation description, so preserves support intersections. The
O6minus/U4(2).2 degree27 identification is retained; the small
alternating-socle plus case is handled by Section3/natural actions.

Equations (2)--(3) forbid the disjoint branch pair from Section2.
ALL THREE remaining binary classical primitive families are therefore
excluded in EVERY degree. These overlap proofs are unbounded finite
linear algebra; small controls are verification aids only.

## 6. Current primitive reduction and the surviving affine profile

Combining PF120's complete imported primitive containers with the
two exclusions above leaves ONLY

    G=A_d or S_d in the NATURAL d-letter action; OR
    G=V:H, V=F2^a, d=2^a, H<=GL_a(2) irreducible.

In the affine case inertia is a transvection. Affine inertia normally
generates G, so its image normally generates H; in particular H is
generated by linear transvections. This uses a NORMAL quotient, not
an assertion about a selected set of ordinary meridian generators.

PF118--119 then give exactly ONE full branch image C, e2/f1 on its
d/4 ramified primes, smooth N, and

    t0+(d/4)M=d/2-1,    M=0 or1.

M=0 means an immersive normalized image, t0=d/2-1; M=1 means one
ordinary cusp, t0=d/4-1. C is noncontractible, hence has double
normalization values. At each such value its two branches are smooth,
have arbitrary contact, and the original Keller open has empty fibre.
The disjoint pair there has complementary half-supports. This is
compatible with binary affine transvections, so it is NOT an exclusion.
Degrees2/3 and the published degree4 retirement remain intact.

No general induction from primitive to arbitrary monodromy is valid
without recognizing the intermediate source; it need not be A2.
Imprimitive cases, the natural alternating/symmetric cases, binary
affine cases, actual embedded multiplication, global trace moments and
the full equality B=C[x,y] remain OPEN. Earlier sealed proofs, all
primes, PF117's closed local trace gate, shared analytic/multiplication
priorities and published quartic retirement are preserved. No historical
priority, formal verification, external review or full goal completion
is claimed.
