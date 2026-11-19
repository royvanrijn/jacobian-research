# Degree-nineteen safety in the full odd-cubic class

Date: 19 November 2026.

Let K have characteristic zero. The statements concern precisely regular
winding-one odd-cubic phases p=t(C+t²D)/x, not arbitrary degree19
Image polynomials or wider height profiles.

**Theorem OddCubicDegreeNineteenSafety.** Every all-pure phase of exact
least sphere degree19 in this class is eventually safe for every fixed
sphere and original polynomial marker. There is no pure unit-center
phase of that exact degree. Zero-center pure phases are safe.

The proof first reduces every possible degree pattern by exact field,
degree and local-order arguments. The final rational identities have
literal integer certificates; no finite-prefix purity inference is used.

## 1. Setup

The zero-center pure branch is safe by the retained all-degree theorem.
Normalize the unit center C(0)=-D(0)=1. Constant C gives C1,D-1 and
degree2, so assume C nonconstant. The canonical polynomial-descent
theorem gives a pure quotient of degrees cbar,ebar, with full stable
field and genuine Mobius return j, and a normalized polynomial of
degreeg such that

    c=g(cbar-1)+1, e=g(ebar-1)+1.

The quotient center has degree cbar>=2. Its selected rational coordinate
is ubar=z+O(z²), with ubar²=-Dbar*z²/Cbar³. The exact least degree is
max(c,e+2); rational selected square roots force equal parity of c,e,
and likewise of cbar,ebar. Dominance is preserved under descent.

We use the retained one-root classification (pure degrees3cbar or
5cbar-2), two-root exclusion, squarefree-center exclusion, and the
already proved unit-center degree7 exclusion. The rational primitive
Rbar is regular at infinity, Omega is logarithmic and j-invariant,
and j fixes the base0. Put Tbar=1/(Cbar*ubar'); its purity identity
is Tbar(z)+Tbar(jz)=2.

An infinity pole of ubar of orderr makes A=j(infinity) a finite
Cbar-root of multiplicityr+1: Rbar(jz)=Rbar(z)-2ubar(z) has a pole
at infinity while Rbar is regular there. For cbar>=3 a moving infinity
zero has order1. Otherwise both dz/Cbar and dubar have positive order
at infinity, so j*(dz/Cbar)=dz/Cbar-2dubar has a positive zero there;
this would give a finite zero of dz/Cbar, impossible. Centers of
degree2 are handled directly by the one/two-root exclusions. A fixed infinity zero makes j=-z and its order is odd.
In that case the Cbar roots occur in opposite pairs, so their number
Hbar is even and at least4. The reflection inequalities are

    cbar >= (5/6)deg(ubar)+Hbar/2,
    5ebar >= 3cbar+6Hbar.

For clarity, the capacity bounds follow from the same local orders
used below. On a paired root orbit of multiplicities m<=n and common
pole order k, one has k<=n-1 and 2k<=3m. Thus m+n>=5k/3+1.
Summing over paired roots gives the first inequality, since the full
pole degree is twice the sum of k. Combining it with
degubar>=s+1 and s=(3cbar-ebar-2)/2 gives the second inequality.

A finite nonzero original u(infinity) would make its exact degree a
multiple of3 by the retained finite-limit theorem, hence cannot occur
at prime degree19. Polynomial descent preserves that infinity type.

## 2. D dominance

Here e17,c<=17, and g(ebar-1)=16. Thusg divides16.

### Infinity pole

If ubar has pole orderr>=1, then ebar=3cbar+2r-2, so

    g(3cbar+2r-3)=16.

For integersg>=1,cbar>=2,r>=1 the only possibilities are

    (g,cbar,r)=(1,3,5),(1,5,2),(2,3,1).

In the first case the required Cbar-root multiplicity6 exceeds degree3.
In the third, the required multiplicity2 in degree3 forces two distinct
roots, already excluded. In the middle, a triple root in degree5 leaves
only pattern3+2 (two roots, excluded) or3+1+1. In the latter the logarithmic
residue at A is zero by j-invariance and regularity at infinity; the
other two roots must be exchanged by j. Their residues are equal and
sum0, hence are zero, contradicting their simple-root nonzero residues.
Therefore no infinity-pole case remains.

### Moving infinity zero

Its order is1, so ebar=3cbar-4 and

    g(3cbar-5)=16.

The only possibilities are

    (g,cbar,ebar)=(1,7,17),(4,3,5),(16,2,2).

The second has quotient degree7, whose unit-center branch is already
excluded. The last has degree4 and center degree2: a one-root center
would have pure degree6 or8, and two roots are excluded. Thus only
cbar7,ebar17,g1 remains for the ten-cell analysis below.

### Fixed infinity zero

Since cbar<=ebar, the growth of Tbar at infinity is(cbar-ebar)/2.
Strict inequality would make Tbar(infinity)=0 at a fixed point, whereas
sigma(Tbar)=2-Tbar requires value1. Hence cbar=ebar=1+16/g, yielding

    cbar=17,9,5,3,2.

The first four have even zero order cbar-1 at a fixed point, contrary
to anti-invariance of ubar. The final center-degree2 case is excluded
as above. This completes D dominance except for c7,e17,g1 with a
moving simple infinity zero.

## 3. C dominance

Here c19,e<=17 and g(cbar-1)=18. The possible cbar are

    19,10,7,4,3,2.

The cbar2 case is excluded by the one/two-root results. For cbar>=3,
ebar<=cbar-2 makes the infinity zero order at least cbar, so the return
fixes infinity. Reflection capacity gives

    cbar>=3Hbar+5>=17.

Henceg1,cbar19. The inequality5ebar>=57+6Hbar, together with odd
ebar<=17 and evenHbar>=4, forces ebar17,Hbar4. Its zero order is19.
Therefore degubar>=20; reflection capacity bounds it by20, and it is
exactly20.

There are two root pairs with common pole orders k,l, k+l=10.
For a pair of multiplicities m<=n one needs

    n>=k+1, m>=ceil(2k/3).

An order-zero pair would cost at least two roots of multiplicity1 and
make total Cbar degree at least20. For positive orders, up to exchange,
(k,l)=(1,9),(2,8),(3,7),(4,6),(5,5). The second and fifth cost at
least20. The other three cost exactly19 and force respectively

    (m,n) pairs (1,2),(6,10);
                (2,4),(5,8);
                (3,5),(4,7).

The three exact rational identities in Section5 below give
nonzero literal constants2944,640,520. Thus every C-dominant case
is excluded.

## 4. The moving c7,e17 root table is exhaustive by hand

Now g1, the infinity zero is simple, and DR=7-Hbar. Since infinity
and its moving partner and the fixed base are three distinct zeros,
degubar>=3; meanwhile degubar<=2DR. Thus3<=Hbar<=5, after the
one/two/squarefree exclusions. The finite Cbar roots are permuted by j. Repeated roots are detected
by poles of Rbar(jz), which cannot occur outside Cbar roots here;
there is no infinity pole because ubar vanishes there. Simple roots
are detected by nonzero Omega residues, and Omega is regular at
infinity. These two observations prove the asserted root permutation.
Besides the base0 there is one other fixed point B. If Hbar is odd,
B is a Cbar-root; if Hbar is even it is not. If B is not a pole of
ubar it is another zero, requiring degubar>=4.

At a paired orbit m<n, the exact ubar pole order is n-1 (the two
primitive-part pole orders cannot cancel); if m=n it lies between0
andn-1. In both cases2k<=3m. At a fixed root of even multiplicitym,
the exact pole order ism-1; at oddm it is0 or a positive odd order
at mostm-2.

For Hbar3 the partitions are5+1+1,4+2+1,3+3+1,3+2+2. Applying
those local rules and the zero lower bound gives, respectively, one,
two, one and two possibilities: exactly rows1--6 of the moving-cell table in Section6. For Hbar4 the partitions are4+1+1+1,3+2+1+1,2+2+2+1.
The first has an impossible(1,4) pair, since pole order3 violates
2k<=3. The second permits only(1,1),(2,3); the third only
(1,2),(2,2) with both pair pole orders1. These are rows7--8.
For Hbar5, pattern3+1+1+1+1 is impossible: fixing the triple root
gives degree at most1, while pairing it with a simple root violates
2k<=3. Pattern2+2+1+1+1 permits either a fixed simple root with
two(1,2) pairs, or a fixed double root with(1,1),(1,2). These are
rows9--10. No other positive partitions occur.

Section6 below reconstructs the exact coordinate change and
the forced numerator of each of those ten rational functions. Its
stdlib checker verifies the entire cleared purity identity and the
106 literal integer multiplier terms, yielding nonzero constants in
all ten cases. This excludes the only remaining D-dominant branch.

Consequently there is no pure unit-center phase of exact degree19
within this winding-one odd-cubic class. The zero-center branch is
safe, so every pure phase of that exact degree is safe. This is not
an obstruction for arbitrary degree19 Image polynomials or for wider
height profiles.

## 5. The three fixed-infinity identities

For each remaining root-pair pattern, choose the signed roots r,-r,s,-s
with the lower multiplicities at r,s, and put t=r/s. The pole orders
(k,l) are(1,9),(3,7),(4,6). Since the simple base zero and order19
infinity zero exhaust degree20, u=Kz/[(z²-r²)^k(z²-s²)^l].
Its derivative has common-denominator numerator

    K[-19z4+((2l-1)r²+(2k-1)s²)z²+r²s²].

Against the same denominator the reciprocal C numerator is a scalar
times(z-r)^A(z-s)^B, with(A,B)=(1,4) in the first pattern and(2,3)
in the other two. Taking its even part and comparing constant, z² and
z4 coefficients gives respectively

    11t²-4t+1=0, t²+4t+19=0;
    5t²-3t+2=0, t+11=0;
    4t²-3t+3=0, t+11=0.

The actual product-rule and reciprocal numerators are reconstructed by
the stdlib checker. The equations are incompatible by the literal
identities

    (3t-1)(11t²-4t+1)+(155-33t)(t²+4t+19)=2944;
    (5t²-3t+2)-(5t-58)(t+11)=640;
    (4t²-3t+3)-(4t-47)(t+11)=520.

The sign t+11 means r=-11s. Thus all C-dominant cases are excluded.

## 6. The ten moving-zero normal forms and their certificates

All remaining cases have c7,e17,g1. The complete root/pole table is

| C multiplicities | fixed multiplicity | paired multiplicities | pair pole orders | fixed pole order | degu |
|---|---|---|---|---|---|
|5,1,1|5|(1,1)|0|3|3|
|4,2,1|1|(2,4)|3|0|6|
|4,2,1|4|(1,2)|1|3|5|
|3,3,1|1|(3,3)|2|0|4|
|3,2,2|2|(2,3)|2|1|5|
|3,2,2|3|(2,2)|1|1|3|
|3,2,1,1|none|(1,1),(2,3)|0,2|0|4|
|2,2,2,1|none|(1,2),(2,2)|1,1|0|4|
|2,2,1,1,1|1|(1,2),(1,2)|1,1|0|4|
|2,2,1,1,1|2|(1,1),(1,2)|0,1|1|3|


Here g1,c7,e17 and s1. Let A=j(infinity), which is finite and nonzero;
j fixes0 with derivative-1. Here T(infinity)=0, so T(A)=2. The
simple u-zero at A has nonzero derivative; T=1/(Cu') therefore makes
C(A) nonzero. Rescale the original coordinate by A,
and rescale u by the same factor, so j(x)=x/(x-1). This preserves
the normalized center, polynomial degrees and purity identity.
Introduce

    x=2w/(1+w),  w=x/(2-x).

Then j sends w to-w. The original infinity is w=-1, its moving
partner x1 is w1, the base is w0, and the other fixed point x2 is
w=infinity. Write

    Ch(w)=(1+w)^7 C(2w/(1+w)).

This is a polynomial with Ch(0)=1. Its degree is7-M, where M is the
multiplicity of the possible fixed root x2 (zero if it is not a root).
Its finite roots are exactly the paired roots in the ten-cell table.
Consequently, choosing inverse root locations a,b with arbitrary signs,

    Ch=product (1-a_i*w)^m_i(1+a_i*w)^n_i.

The parameters satisfy a_i!=0,1,-1, and a_i²!=a_j². These inequalities
hold because the roots are distinct and neither moving infinity point
is a C-root. They are not needed to make the eventual polynomial
certificates nonzero: the certificate ideals below are unsaturated.

Put U(w)=u(2w/(1+w)). It is odd, U'(0)=2, and has simple zeros at
w0 and w=1,-1. Its paired pole orders are the k_i in the table and
its possible pole at infinity has the table's fixed-root order K.
Its numerator and denominator are therefore forced to have the form

    U=2w(1-w²)*(1+q*w²)^epsilon /
                      product (1-a_i²*w²)^k_i.        (N1)

Here epsilon1 occurs exactly for table rows2,3,5 (one-based), and
epsilon0 for the other seven. This is exhaustive, not an ansatz:
in those three rows the numerator is odd of degree at most5, with
the three specified zeros, while in the other rows it has degree3.
In rows3,5 its degree must equal5 and hence q!=0. In row2 degree3
is also allowed by taking q0. Valid phases also require q!=-1 and
no cancellation of the specified poles. The unsaturated certificates
exclude all parameter values, so weakening these guards loses no
valid phase. The coefficient2 follows from U'(0)=2.

The differential is exactly

    dx/C = 2(1+w)^5 dw/Ch(w).

Since U' is even, T+Tj=2 becomes the FULL rational identity

    U'(w)=(1+w)^5/Ch(w)+(1-w)^5/Ch(-w).                (N2)

The rational functions in(N2) have explicitly known denominators.
After multiplying by their least common symmetric denominator, every
coefficient must vanish. The stdlib verifier constructs this numerator
directly from(N1),(N2), without a symbolic-algebra package. Its nonzero
coefficients, up to factors plus/minus2, are the following integer
polynomials, in exactly the order of the ten-cell table:

1. 3a²-5; a²+13.
2. a(aq+5a-2); 3a²q-13a²+20a-5q-5;
   4a²+10a+3q-13.
3. a(3aq+1); a²q-a²-10a+5q+5; a²+5a+3q-13.
4. a²+5; 3a²-13.
5. a(aq+1); a²q-a²+10a-5q-5; 3a²+5a+3q-13.
6. a²-5; a²-13.
7. b(a²b+1); 3a²b²-3a²+b²-10b+5;
   a²-3b²-5b+13.
8. a(ab²+1); 3a²b²-a²-10a-b²+5;
   a²+5a+b²-13.
9. a²b²-5ab+a+b; 3a²b²-a²+10ab-10a-b²-10b+5;
   a²-ab+5a+b²+5b-13.
10. b(a²b-1); a²b²-3a²-b²-10b+5;
    a²-b²-5b+13.

For each row, the retained exact integer multiplier polynomials in
odd_cubic_degree_nineteen_safety_certificates.json sum these inputs to a nonzero
integer. The ten constants, in order, are

    44,3280,33920,28,106240,8,9152,128,11648,2816.

There are106 multiplier monomials in total and maximum coefficient
identity degree7. The portable verifier multiplies and adds every
literal term, reconstructs every input from the actual full rational
identity, and rejects a corrupted constant. Thus none of the ten
normal forms satisfies(N2) over any characteristic-zero field.
Sympy was used only to discover the integer multipliers. A Groebner
basis output is not the retained proof or the verifier's input.

Combined with the exact global case reduction and the three fixed
infinity certificates, this excludes every pure unit-center phase of
degree19. The zero-center pure branch is safe by the retained all-degree
theorem. This proves OddCubicDegreeNineteenSafety within its stated class.


## 7. Sources and exact reproduction

The proof imports the canonical polynomial-descent theorem from
[ODD_CUBIC_CANONICAL_POLYNOMIAL_DESCENT.md](ODD_CUBIC_CANONICAL_POLYNOMIAL_DESCENT.md),
the one-root classification from
[ONE_POLE_ODD_CUBIC_CLASSIFICATION.md](ONE_POLE_ODD_CUBIC_CLASSIFICATION.md),
and the unit-center degree7 exclusion and two-root lemma from
[ODD_CUBIC_DEGREE_SEVEN_SAFETY.md](ODD_CUBIC_DEGREE_SEVEN_SAFETY.md).
The all-degree zero-center safety and finite-nonzero-infinity obstruction
are in
[ODD_START_HEIGHT_RATIONALITY_AND_QUINTIC_SAFETY.md](ODD_START_HEIGHT_RATIONALITY_AND_QUINTIC_SAFETY.md).
The underlying rational inverse coordinate and exact/log separation are
retained dependencies of the canonical descent theorem. Every necessary
source edition and input is pinned in the source ledger.

The certificate input contains all ten normal-form records, every
literal multiplier monomial and the three reflection certificates. Run

    python3 research/foundations/extended-geometry/verify_odd_cubic_degree_nineteen_safety.py

The checker uses only Python's standard library, rejects optimized
assertion-disabled execution, reconstructs the full rational identities,
checks all106 literal moving-case multiplier terms, and checks the
three complete reflection identities. It also verifies the exact global
integer case reductions and corrupted-input controls. Its receipt is
odd_cubic_degree_nineteen_safety_verification.json. These controls audit
the written proof and exact algebra; they are not external review,
formal verification or a general finite-prefix purity test. No global
classification beyond this odd-cubic class is asserted.

The [source ledger](odd_cubic_degree_nineteen_safety_source_ledger.json), [literal certificate inputs](odd_cubic_degree_nineteen_safety_certificates.json), [portable checker](verify_odd_cubic_degree_nineteen_safety.py) and [receipt](odd_cubic_degree_nineteen_safety_verification.json) keep the complete global reduction, root-orbit coverage and exact algebra together. The certificate reconstruction is independent of the symbolic-algebra discovery code.
