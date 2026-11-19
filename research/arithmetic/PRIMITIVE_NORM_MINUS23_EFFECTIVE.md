# Effective discriminant -23 primitive norm sums

Let Q=x²+xy+6y², P its positive primitive numerical values, and
U={n inP:23 does not dividen}. Write R_A(N) for the number of
ordered numerical pairs u,v inA withu+v=N. All intermediate shell
counts below are signed ordered integer quadruple counts.

Write R=Z[omega], omega=(1+sqrt(-23))/2, and O=R+Rj,
with j²=-1 and ja=bar(a)j. Its norm is Q plus Q. As proved in
the split-target note, O is maximal in the definite quaternion algebra
of discriminant23, with units plus/minus1, plus/minusj, total
projective mass11/6 and principal probabilityw=3/11.

**Theorem CN3EffectiveAll.** For every N>=2^249,

    N belongs to P+P if and only if529 does not divideN.

More quantitatively, at these admissible targets,

    R_P(N)>=N^(5/7)/15000000000.

The onset is explicit but deliberately coarse. No conclusion about
the unresolved interval between the retained finite scan and this
onset is asserted.

## E1. The joint small-prime baseline

First take23 coprime toN and write N=2^a3^bm, gcd(m,138)=1.
Set d_p(0)=1 and d_p(t)=(p+1)p^(t-1) for t>=1. Let

    B_N=4w d2(a)d3(b)sigma(m),        w=3/11.

The positive all-four-coordinate primitive operators at two and three
are D_(p,0)=1 and D_(p,t)=T_(p^t)-T_(p^(t-2)). On the
complement of norm characters,

    ||D_(p,t)||<=2(t+1)p^(t/2)        (includingt=0).

The exact shell with neither binary block divisible by two or three
has count G23(N). CN3JointSmallPrimes gives

    G23(N)>= (5/12)B_N-144 tau(N)sqrt(N).                (E1)

Indeed its full error144(a+1)(b+1)tau(m)sqrt(N) equals the
displayed error, since(a+1)(b+1)tau(m)=tau(N). This step uses
the actual mixed normalizer/Iwahori operators, including the two
principal axis points in the doubly divisible branch.

## E2. Bad primes not dividing the target

Let p>=5, p!=23 and p not dividingN. Modulo p, a quaternion
with one binary block zero belongs to the normalizer of R_p's
Cartan subgroup: split ifp splits inR, nonsplit ifp is inert.
Conversely this normalizer is precisely the union of R_p-units
and jR_p-units, hence precisely this bad locus. Its maximal-compact
index is c_p=p(p-1)/2 in the inert case and p(p+1)/2 in the
split case. The unit norm image is all Zp-units in either case:
use the finite-field norm in the nonsplit Cartan and diag(u,1) in
the split Cartan, with principal-unit lifts. Thus the fine quotient
has only the constant global norm character.

The principal point of this normalizer level has actual stabilizer
{plus/minus1,plus/minusj} and probability mass w/c_p. Counting
bad shell elements which remain four-coordinate primitive at two
and three therefore gives the exact return

    b_p(N)=4(D2 D3 T_m f_p)(x_p).

The main term is B_N/c_p. Its error is at most

    4||D2||||D3||tau(m)sqrt(m)<=16tau(N)sqrt(N).          (E2)

This is an ordinary congruence return at a prime not dividing the
index. It makes no unramified spectral assertion at the level primep.

## E3. Bad primes dividing the target

Write m=p^t q, t>=1, gcd(q,138p)=1. Ifp splits inR, separate
the rank-one branch outsidepO and the full-content branchpO.
The former is exactly the retained positive endpoint correspondence
with the commuting D2,D3 factors retained:

    rank-one bad=8(V_(p,t)D2 D3 T_q f_A)(x_A).

Its main term is16w d2(a)d3(b)p^(t-1)sigma(q)/(p+1).
The complementary norm estimate gives error at most

    64t(a+1)(b+1)tau(q)sqrt(N)
       <=64tau(N)sqrt(N).                               (E3)

Here t tau(q)<=tau(m). The content branch is in exact bijection,
by division of all four coordinates byp, with the two-and-three
primitive shell at N/p². Its main term is
4w d2(a)d3(b)sigma(p^(t-2))sigma(q), where the sigma factor
is zero fort=1. Its error is at most16tau(N)sqrt(N)/p.
Thus the complete bad count satisfies

    b_p(N)<=delta_(p,t) B_N+80tau(N)sqrt(N),             (E4)

where

    delta_(p,t)=
      [4p^(t-1)/(p+1)+sigma(p^(t-2))]/sigma(p^t)
        <= (5p-3)/(p²(p+1))=:delta_p.                   (E5)

The last bound also exceeds the unit-target split value1/c_p.
It can be checked by the elementary geometric sums, or by the
monotone Haar ratio in the separate mass-margin note. In particular
the full four-coordinate content branch has NOT been omitted.

Ifp is inert and dividesN, a zero binary block modulo p forces
the other block to vanish, by anisotropy of the quadratic field norm.
There is no rank-one bad branch. The content branch alone gives
the same weakened error(E4), with its main ratio at mostp^(-2),
and hence at most the uniform inert bound

    delta_p=2/[p(p-1)].                                 (E6)

Consequently(E4), with80 as a common error constant, holds for
everyp>=5,p!=23, whether or notp dividesN. All auxiliary compact
norm images here are full local unit groups, so only constants were
removed in(E2)--(E4).

## E4. A strong ramified selection at23

We force both binary norms to be prime to23. Put pi=sqrt(-23).
At this odd ramified place the given maximal order has the explicit
description Z23[j] plus Z23[j]pi, withpi²=-23 andpi j=-j pi.
Its quotient modulo pi isF23[j]=F_(23²), sincej²=-1 is
irreducible modulo23, and the ramified quadratic subring R maps
toF23. For a norm-unit quaternion, either binary norm is zero
modulo23 precisely when the residue lies in

    F23^* union jF23^*.

This subgroup has order44 in F_(23²)^* of order528, so its
inverse-image compact C23 has index12. The principal point has
massnu=w/12 and actual stabilizer of order four. Its compact norm
image consists exactly of units square modulo23: both F23 scalars
and j have square norm, while scalar squares supply every element of
1+23Z23 by Hensel lifting. Thus its global norm characters are exactly1 and
chi23=(./23) on positive rational indices. This character is
understood as the quadratic Hecke character, including its infinite
component; it need not be an even Dirichlet character.

Here is the full norm-character quotient argument. The projective
center forces a norm character to be quadratic. By the Hasse--Schilling
norm theorem (Voight, Main Theorem14.7.4), the norms of invertible
rational quaternions are exactly the positive rationals. Thus the relevant
finite norm quotient is the finite ideles modulo positive rationals,
idele squares, and the compact unit norm images. Multiplying by a
positive rational clears all finite valuations. At a level whose compact
norm images are all local units this leaves the trivial quotient, as
used above and in the joint small-prime argument. At C23 it leaves
exactly F23^*/(F23^*)², of order two. The quotient is attained:
local reduced norms are surjective, including at23 by the explicit
residue-field norm and the norm-valuation-one elementpi. Its nontrivial
character has the stated Legendre values on positive indices. The odd
infinite component makes it trivial on all rational principal ideles;
no evenness restriction is imposed on its finite Dirichlet part.
Consequently1 andchi23 are orthogonal unit norm characters, and
the principal point has value one under each, so its indicator has
coefficientnu on each.

The principal point indicator has the same Fourier coefficientnu
on both norm characters. Moreoverchi23(2)=chi23(3)=1. Hence
the count b23ram(N) of bad elements retaining the small-prime
four-coordinate primitive conditions has main term

    B_N*(1+chi23(N))/12 <=B_N/6.

After both characters are removed, the coprime Hecke norm bound
and the D2,D3 norm bounds give complementary error at most
16tau(N)sqrt(N), as in(E2). Thus

    b23ram(N)<=B_N/6+16tau(N)sqrt(N).                    (E7)

This gate excludes norm-divisibility by23, which is stronger than
excluding coordinate-gcd23. It is exactly the selection needed
for U+U and the subsequent single23 scaling.

## E5. Strict infinite local margin

The exact rational mass packet proves

    sum_(p>=5,p!=23)delta_p<6/25.                       (E8)

For completeness, enumerate primes toR=65536, split according
to(p/23), use(E5),(E6), and round each term UP to a multiple of
2^(-80). The rounded numerator is285050268339193439840907.
Every omitted mass is less than5/p². Summing over all odd
integers greater than the evenR gives tail at most5/(2R)+5/R².
The resulting upper bound is exactly

    285096386606752597273227/1208925819614629174706176
       <6/25.

Thus the main term after(E1),(E7), and any finite subset of the
other bad primes is at least

    (5/12-1/6-6/25)B_N=B_N/100.                         (E9)

## E6. Elementary large-prime tail

For a positive k the signed binary representation count satisfies
r_Q(k)<=2tau(k): each generator gives a principal ideal of normk,
each principal ideal has two generators, and the total number of
quadratic ideals of normk is at mosttau(k). Also r_Q(0)=1.

For every X>=1 the full binary lattice count S_Q(X), including zero,
satisfies

    S_Q(X)<=6X.                                        (E10)

Indeed Q=(x+y/2)²+23y²/4, so the elementary row bound is
(4sqrt(X/23)+1)(2sqrt(X)+1), which is at most
2X+3sqrt(X)+1<=6X.

The exact Euler-maximization packet below proves

    tau(k)<=C k^(1/7),         C=1424.                  (E11)

Suppose Y>=1. For a primep>Y, p<=sqrt(N), any binary block
divisible byp is p times a point counted by S_Q(N/p²). The
remaining binary block has at most2C N^(1/7) possible signed
representations for its required norm, also when that norm is zero.
Taking either block gives an unrestricted upper bound

    4C N^(1/7) S_Q(N/p²)<=24C N^(8/7)/p².

The selected G23 shell already excludes a zero binary block, since
zero is divisible bytwo. Therefore a prime greater thansqrt(N)
cannot be a surviving coordinate-gcd obstruction, and no additional
tail beyond that range is omitted. Summing all integers above the
integerY and using sum_(n>Y)n^(-2)<=1/Y gives tail at most

    24C N^(8/7)/Y.                                    (E12)

## E7. Uniform positive coordinate count and exact onset

Take the integerY=floor(N^(1/4)), withN>=16. Subtract(E7) and
all(E4) forp<=Y from(E1), then subtract(E12). The number of
finite primes is at mostY. The summed spectral error is at most

    (160+80Y)tau(N)sqrt(N)
       <=240Y tau(N)sqrt(N)<=240C N^(25/28).

SinceY>=N^(1/4)/2, the tail in(E12) is at most
48C N^(25/28). Also d2(a)d3(b)sigma(m)>=N, so B_N>=12N/11.
Thus the signed shell S_U(N) whose two binary norms belong toU
satisfies

    S_U(N)>=(3/275)N-288C N^(25/28).                    (E13)

In particular whenever

    N^(3/28)>=52800C,                                  (E14)

we have S_U(N)>=3N/550. Each ordered pair of numerical values
has at most4C²N^(2/7) signed representations, by(E11), hence

    R_U(N)>=3N^(5/7)/(2200C²).                          (E15)

The exact constant certificate proves(E11) without floating point:
forp>=128, a+1<=2^a<=p^(a/7). Forp<128 enumerate the
maximum of(a+1)^7/p^a at0<=a<=28. This is complete because
the successive ratio decreases witha and(30/29)^7<2. Multiply
the finitely many rational maxima. Their product is at most1424^7.

The same certificate verifies the exact integer inequalities

    2^(245*3)>=(52800C)^28,
    2^(249*3)>=23³(52800C)^28.                          (E16)

Thus(E15) applies to every23-unitN>=2^245, and also to every
23-unitN>=2^249/23 directly through(E14). The latter rational
threshold is used for scaling; it does not rely on rounding to the
larger power-of-two unit cutoff.

## E8. One ramified scaling and the numerical-value bound

ForN>=2^249 with529 not dividingN, either23 does not divideN,
and(E15) already supplies a U+U decomposition, orN=23M with
23 coprime toM. ThenM>=2^249/23 satisfies(E14), and multiplying
both U values by23 preserves primitivity. Concretely multiplication
of x+yomega bysqrt(-23)=2omega-1 sends its coordinate column
to(-x-12y,2x+y), with determinant23. Away from23 this is
invertible. If both output coordinates vanished modulo23, then
x+12y would vanish too, forcingQ(x,y)=0 modulo23, contrary toU.
Thus its output coordinates are primitive. The ordered numerical
pairs lift injectively. Hence

    R_P(N)>=3(N/23)^(5/7)/(2200C²)
             >N^(5/7)/15000000000,

since23^5<10^7 and22000C²<45000000000.

For a primitive binary norm, the23-adic valuation is at most one,
with square unit quotient after removing that possible factor23.
If a sum is divisible by23, nonsquareness of-1 modulo23 forces
both summands divisible by23. Dividing by23 then gives two unit
square norms whose sum cannot be zero modulo23. Thus529 dividesN
is an exact obstruction at every size, proving the other direction.

This gives an explicit all-large theorem, while
the complete exact40-exception list still requires a finite bridge
or a further structural reduction below2^249.

## Joint two-and-three operator proof

Use Q=x^2+xy+6y^2, O=R+Rj, w=3/11 and the exact local
correspondences from the retained split-target note. Both two and three
split in R. Let N=2^a 3^b m, a,b>=1, gcd(m,138)=1.
Let b23(N) be the signed ordered O-norm-N elements which are outside
2O and 3O, and have a binary coordinate block divisible by two and
a (possibly different) binary coordinate block divisible by three.
Zero binary blocks are included in these bad conditions.

**Proposition J23.** Put d=2^(a-1)3^(b-1). Then

    |b23(N) - (16w/3)d sigma(m)| <=64ab tau(m) sqrt(N).

When a=b=1 the sharper error 16 tau(m)sqrt(m) holds.
The main term is exactly 1/9 of the full two-and-three-primitive
shell baseline 48w d sigma(m).

## Two points, rather than a single axis point

At the product Iwahori level, the principal fibre has the two axis
choices at two and the two axis choices at three. The global projective
unit j exchanges both pairs simultaneously. Thus these four flag pairs
give TWO distinct principal points x,y. Each has trivial projective
stabilizer and probability mass

    nu=2w/((2+1)(3+1))=w/6=1/22.

At an individual prime there are four equally large image/kernel axis
cases. At the pair of primes, the 16 cases fall into four orbits under
left and right multiplication by j. The row and column labels x,y
record the two invariant relative image-axis and kernel-axis choices.
Within each orbit, choose one image/kernel pair at each prime; its
principal integral return is divided by the two generators plus/minus1.
Consequently, with h=1_x+1_y and

    A=V_(2,a) V_(3,b) T_m,

the exact signed-count identity is

    b23(N)=8 sum_(i,j in{x,y}) (A 1_j)(i)
           =8< h,A h >/nu.                              (J1)

Each V_(p,t)=B_p^(t-1)F_p is the positive nonbacktracking
endpoint operator of degree p^(t-1). The local factors commute,
as they have different prime supports; they are self-adjoint.
The only norm character at the product Iwahori level is constant:
the Iwahori compact unit norms are all local units. There is no
additional quadratic contribution to (J1).

## Constant term and full complementary norm

Since mean(h)=2nu and ||h-2nu||^2=2nu-4nu^2,
the constant term in (J1) is

    (8/nu)*(2nu)^2*d*sigma(m)
        =32nu*d*sigma(m)=(16w/3)d*sigma(m).

On the orthogonal complement of constants, the retained full-space
bound for each V (including its nonspherical joint-projection kernel)
and the coprime Hecke bound give

    ||A|| <=4ab sqrt(2^a3^b) tau(m)sqrt(m).

Thus the error in (J1) is at most

    (8/nu)*(2nu-4nu^2)*||A|| <=16||A||,

proving the stated error. For a=b=1 both local operators are unitary
flips, so ||A||<=tau(m)sqrt(m). The primitive shell degree at two
is3*2^(a-1), at three4*3^(b-1); multiplying4w gives the displayed
baseline. The ratio 1/9 equals (4/9)*(1/4), the product of the
two rank-one bad densities.

## Consequence for joint small-prime selection

Let G23(N) count the signed ordered O-norm-N elements with neither
binary block divisible by two or by three. The same conditions already
exclude 2O and3O. Inclusion-exclusion on the two-and-three-primitive
shell gives its main factor

    1-4/9-1/4+1/9=5/12.

For t>=1 the positive primitive maximal-level operator is
D_(p,t)=T_(p^t)-T_(p^(t-2)), with T_(p^negative)=0. Its degree is
(p+1)p^(t-1) and its complementary norm is at most
2(t+1)p^(t/2). Apply the split-target return argument with the
other prime's commuting D operator. The base-shell error is at most
16(a+1)(b+1)tau(m)sqrt(N), the two individual bad errors are at
most32a(b+1)tau(m)sqrt(N) and32b(a+1)tau(m)sqrt(N), and the
joint error is bounded above. Hence the rigorously weakened estimate

    |G23(N)-20w d sigma(m)|
      <=16(9ab+3a+3b+1)tau(m)sqrt(N)                    (J2)

follows. The estimate is about coordinate gcds at two and three only.
It does not assert that either binary norm is globally primitive.

## Independent 22-state check

The exact checker constructs the complete maximal class list
O,I_(3,1),I_(2,0) from the split-target packet. It enumerates every
right-ideal norm-two and norm-three transition J=alpha I_j subset I_i.
Completeness follows because n_j O subset I_j implies
alpha in n_j^(-1)I_i, so all required alpha occur among the bounded
integer norm vectors beta=n_j alpha, of norm p n_i n_j.
Testing containment on all four ideal basis vectors is sufficient.
The p+1 image lines modulo p give all neighbours.

For local trivializations I_i,p=g_i M_2(Z_p), the matrix
g_i^(-1)alpha g_j maps the new vertex to the old. At the step prime,
its image and kernel are the forward and reverse flags; at the other
prime its inverse transports the retained flag. All reductions use
exact rational coordinates and Hensel roots of omega^2-omega+6.
Quotienting pairs of flags by the actual left-order units produces
22 states. Every state has trivial projective stabilizer and uniform
mass1/22. The two principal axis states have labels0,1.

The computed F2,F3 are commuting involutions. The B2,B3 have row
degrees2,3, commute with each other and the other-prime flip, and
Fp Bp Fp=Bp* under the uniform measure. For every1<=a<=6,
1<=b<=4, the exact operator return in(J1) equals a separate complete
signed-coordinate enumeration, through the largest target5184.

The joint-model section of the exact checker supplies these matrices and coordinate controls. These are finite checks of
the actual combined correspondence; they do not replace the written
unrestricted correspondence or its classical spectral input.

## All valuation combinations, including unit targets

At a unit target, the binary-block bad locus in GL2(Fp) is exactly
the normalizer of the split diagonal torus: diagonal or antidiagonal
matrices. Equivalently it preserves the unordered pair of coordinate
lines. Let C_p be the inverse image of this normalizer in the maximal
compact. Its index is p(p+1)/2, hence c2=3,c3=6. All compact
norm images remain Zp-units, since every diagonal matrix diag(u,1)
belongs to C_p. In particular these normalizer levels add NO global
norm character.

The principal normalizer point has all four actual O units as its
stabilizer, since j preserves the unordered axis pair. Its mass is
w/c_p. A unit-index Hecke return at this point counts the corresponding
binary-block bad coordinates divided by FOUR. Thus the unit-target
single-prime bad density is1/c_p: 1/3 at two,1/6 at three.

The joint returns are as follows. In this table Vp is V_(p,t),
and every displayed T_m is coprime to138.

| Target valuations | Principal axis points | Point mass | Exact joint bad count |
|---|---:|---:|---|
| a=b=0 | one | w/18 | 4(T_m f)(x) |
| a>=1,b=0 | one | 2w/(3*6)=w/9 | 8(V2 T_m f)(x) |
| a=0,b>=1 | one | 2w/(3*4)=w/6 | 8(V3 T_m f)(x) |
| a,b>=1 | two | w/6 each | 8 sum_(i,j)(V2 V3 T_m1_j)(i) |

For a mixed row, j exchanges the two oriented axis flags at the
target-dividing prime, while preserving the unordered axis pair at
the unit prime; thus there is one two-element orbit and its stabilizer
has actual order two. Multiplying the four target-dividing matrix-axis
cases by two stabilizer units gives its factor eight. For the first
row j already stabilizes both unordered pairs, so the actual stabilizer
has order four. This proves the counts, rather than applying the same
factor at all levels. Their main bad fractions against the appropriate
two-and-three-primitive shell are respectively1/18,2/27,1/12,1/9.

**Theorem CN3JointSmallPrimes.** Let a,b>=0, gcd(m,138)=1, N=2^a3^bm.
Let G23(N) count all signed ordered norm-N elements for which neither
binary block is divisible by two or three. Define

    d_p(0)=1,       d_p(t)=(p+1)p^(t-1) for t>=1;
    rho2(0)=2/3,    rho2(t)=5/9 for t>=1;
    rho3(0)=5/6,    rho3(t)=3/4 for t>=1.

Then

    |G23(N)-4w d2(a)d3(b)rho2(a)rho3(b)sigma(m)|
       <=144(a+1)(b+1)tau(m)sqrt(N).                    (J3)

In particular its main selected fraction is at least5/12 of the
full two-and-three-primitive baseline, uniformly in target valuations.

To verify the error, let D_(p,0)=1 and use D_(p,t) for t>=1.
The full primitive-shell count is4(D2 D3 T_m f0)(x0), with massw.
Its complementary error is at most4||D2||||D3||tau(m)sqrt(m),
at most16(a+1)(b+1)tau(m)sqrt(N). At a target-dividing prime the
single bad count is8 times the Vp return, with the other primitive
D factor; at a unit prime it is4 times the normalizer-level return
with the other D factor. Each single-bad error is bounded above by
32(a+1)(b+1)tau(m)sqrt(N). The joint-return table and its one- or
two-point coefficient estimate give an error at most
64(a+1)(b+1)tau(m)sqrt(N). Adding these four errors proves(J3).

The main factors use exact inclusion-exclusion. The individual bad
fractions are (1/3 or4/9) at two and (1/6 or1/4) at three; each
joint fraction in the table is their product. Therefore the surviving
products are5/9,25/54,1/2,5/12 in the four valuation rows. No
independence assumption replaces a return calculation here.

The new proof uses only the retained split-target theorem's tree and
classical transfer setup. The operators D and the coprime T_m remain
at maximal index primes; the V estimates are explicitly on the
complement of constants, with the two Iwahori nonspherical sectors
treated before taking their commuting product. The mathematical argument and the finite computations have been independently reconstructed internally. No external human review, formal verification or historical priority is asserted.


## Exact finite bridge for the 40-exception question

Let E be the explicit 40-element set in Proposition CN3Finite of
[the class-number-three note](PRIMITIVE_NORM_CLASS_NUMBER_THREE.md).
Every member is at most3818, hence below2^24.

**Lemma CN3FiniteBridge.** Under the proved bounded classification and
Theorem CN3EffectiveAll, the following statements are equivalent:

1. For every N>=2, N belongs to P+P exactly when529 does not divide N
   and N is not in E.
2. Every N with2^24<N<2^249 and529 not dividing N belongs to P+P.

The forward implication restricts the global formula to the displayed
interval, which contains no member of E. Conversely partition N>=2 into
N<=2^24, the displayed middle interval, and N>=2^249. Proposition
CN3Finite proves the first range, statement2 proves the admissible
middle targets, and CN3EffectiveAll proves the last range.
The local529 obstruction in E8 excludes every forbidden middle target;
the bounded result supplies every actual exception in E. This proves
both directions. The middle assertion has not been proved, so this
equivalence does not establish global exception completeness.

## Sources and exact finite support

The infinite input is the retained
[classical all-ideal Hecke transfer](../../results/eisenstein-restricted-prime-sums/sections/hecke-transfer.tex),
Lemma `lem:hecke-transfer`, whose complementary bound removes ALL
global reduced-norm characters, applied at maximal index primes, and the
[split-target endpoint theorem](PRIMITIVE_NORM_MINUS23_SPLIT_TARGETS.md).
At Iwahori level the full nonspherical sectors are bounded in that
theorem; no unramified assertion is made at a nonmaximal index prime.
The normalizer levels in this note have full local unit norms; the
explicit ramified23 compact has exactly the two stated norm characters.

The [portable exact checker](verify_primitive_norm_minus23_effective.py)
and [receipt](primitive_norm_minus23_effective_checks.json) include
the actual joint and mixed2/3 Brandt operators, the28-state fine23
axis compact with both norm characters, complete signed coordinate
comparisons, and all rational Euler/divisor/tail/onset comparisons.
The actual fine23 T5 operator also checks the negative-character branch:
T5 chi23=-6chi23, and returns at5,25,125,625 agree with the
complete signed-coordinate counts0,4,0,148.
The checker pins its split-target support source. Finite computations
support the written infinite correspondences and normalizations; they
do not replace the unrestricted spectral and lattice arguments.

Structural primary sources are Voight, Quaternion Algebras, author
edition1.0.7u: Main Theorem14.7.4 for the positive rational norm image,
Theorem25.1.1 for maximal projective mass,
Section23.5 for lattice trees and edge correspondences, and
Propositions41.3.1/41.3.6 for all-integral-ideal Hecke operators.
The shared transfer imports its pinned Jacquet--Langlands and Deligne
theorems. No fresh full revalidation of those deep proofs is claimed.
The written arguments and exact replay have internal review only;
external human review, formal verification and historical priority
remain separate.
