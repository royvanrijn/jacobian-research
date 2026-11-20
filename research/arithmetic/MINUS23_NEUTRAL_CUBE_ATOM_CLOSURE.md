# Every nonprincipal neutral cube atom is covered

20 November 2026. Internal written proof and exact finite certificate.

Let Q(x,y)=x²+xy+6y², P its positive primitive numerical values, and
U the values in P prime to23. A split prime is principal when represented
by Q; the other split primes are called nonprincipal. The two statements
below use the already retained all-ideal Hecke transfer, including the
local endpoint operators, exactly as in PRIMITIVE_NORM_MINUS23_EFFECTIVE.
No new spectral theorem is being asserted.

**Theorem N23CubeAtoms.** For every nonprincipal split prime p and every
c in{5,10,15}, one has cp³ in U+U.

**Proposition N23SparseAtomsEffective.** If p is any split prime other
than2,3 and N=cp³, c in{5,10,15}, then N>=2^43 implies N in U+U.
If p is principal and N=cp, the stronger sufficient onset is N>=2^40.
The latter does not settle all principal-prime atoms.

## 1. Classical comparison and exact imported estimates

Rouse–Thompson arXiv2206.00412v1 Theorem1 is an explicit unrestricted
quaternary representation bound for PRIME discriminant. The Gram
discriminant of Q+Q is23², and its coordinate-block gcd restrictions
are additional: that theorem does not give the present conclusion.
Its introduction explains the classical Eisenstein/Deligne decomposition.
Fuchs–Hsu–Rickards–Schindler–Stange, 27April2025, Theorem1.1 counts
primes represented by a fixed shifted primitive binary form; its
varying-prime average does not imply every cp³ is a primitive norm sum.
These actual primary statements were read before this deduction.

The new explicit prime bound used below is Rosser–Schoenfeld,
Approximate formulas for some functions of prime numbers (1962),
Section3 Corollary1(3.6), printedp69/PDFp6:

    pi(x)<1.25506 x/logx <4x/(3logx), x>1.             (A1)

Partial summation, discarding its negative lower-end term, gives

    sum_(prime ell>Y)1/ell² <=8/(3YlogY), Y>1.         (A2)

Indeed the tail is -pi(Y)/Y²+2 integral_Y^infty pi(t)/t³ dt,
which is at most(8/3) integral_Y^infty dt/(t²logt).
The latter integral is at most1/(YlogY).

The retained exact positive quaternion shell has O=R+Rj, massw=3/11.
For N=2^a3^b m, gcd(m,138)=1, put

    B=4w d2(a)d3(b)sigma(m), d_p(0)=1,d_p(1)=p+1.

The shell G with NEITHER binary block divisible by2 or3 has

    G >=rho2(a)rho3(b)B-144tau(N)sqrtN,                (A3)

where rho2(0)=2/3,rho2(1)=5/9,rho3(0)=5/6,rho3(1)=3/4.
This is CN3JointSmallPrimes/J3, including all unit/divisible rows.
For ell>=5,ell!=23 not dividing N, the bad block-gcd count is a
normalizer-level Hecke return with main B/c_ell and error

    <=4||D2||||D3||tau(m)sqrtm.                        (A4)

Here c_ell=ell(ell-1)/2 for inert ell, ell(ell+1)/2 for split ell;
D_(p,0)=1,D_(p,1)=T_p, and ||T_p||<=2sqrtp on the complement
of norm characters. The auxiliary compact has full local unit norm
image, hence only the constant norm character. Formula(A4) is the
EXACT error expression preceding the weakened constant in E2.
For a target-dividing split p, E4 gives

    b_p <=delta_p B+80tau(N)sqrtN,
    delta_p=(5p-3)/(p²(p+1)).                         (A5)

For every prime ell>=5,ell!=23, define delta_ell PIECEWISE: it is
(5ell-3)/[ell²(ell+1)] when ell splits and2/[ell(ell-1)] when ell
is inert. The latter is the uniform inert majorant. The retained
exact rational mass certificate E8 for this piecewise sequence proves

    sum_(ell>=5,ell!=23)delta_ell<6/25.                (A6)

All these are imported with precisely their current hypotheses;
in particular(A5) includes the full four-coordinate-content branch.

## 2. Two local exclusions, and the actual main coefficients

The analytic argument in Sections2--4 excludes p=2,3; the complete
finite bridge in Section5 includes and proves both of those atom cases.
In this analytic range p>=13, since all smaller odd primes are inert
in Q(sqrt(-23)). We have v5(N)=1. If a binary block were
divisible by5, its norm would vanish modulo5; anisotropy at this
inert prime forces the other block divisible by5 too. Then25|N,
a contradiction. Thus the entire5-bad event is ZERO, and its
majorant delta5=1/10 can be removed from(A6). The remaining mass
sum is less than7/50. This is an exact local exclusion, not a
primitivity inference from squarefree target content at other primes.

Second chi23(N)=(N/23)=-1: split p,2,3 have character+1, while
5 has character-1. Modulo23 the binary norm is(x+12y)². If one
binary norm vanished, the other would force N to be a square or
zero modulo23. Hence EVERY norm-N shell element already has both
binary norms in the23-unit locus. No ramified23 bad count or error
needs subtraction.

For t=1 or3 and N=cp^t, sigma(5p^t)=6sigma(p^t), and
sigma(p^t)>=p^t. The selected main term after all remaining bad
prime masses is therefore STRICTLY greater than N times the
following constants:

    c5 : (72/55)(5/9-7/50)=748/1375;
    c10: (108/55)(25/54-7/50)=872/1375;
    c15: (96/55)(1/2-7/50)=864/1375.                   (A7)

All three exceed1/2. These use the actual d2,d3 and selected
fractions, rather than the weaker all-target B>=12N/11 bound.

For t=3, tau(N)<=16 and tau(m)=8. Exactly one of a,b may be1.
Consequently(A3) costs at most2304sqrtN. Formula(A4) costs at
most64sqrtN for each unit-index bad prime: with a=b=0 its
constant is32, otherwise4*(2sqrtq)*8sqrt(N/q)=64sqrtN.
We always subtract the target-prime event p using(A5), whether
or not p<=Y; it costs at most1280sqrtN. Unit-index primes up
toY number at most pi(Y). Subtracting p again in a later tail
is harmless, but is not needed for correctness of the finite set.

Thus the total finite spectral error is at most

    3584sqrtN+64pi(Y)sqrtN.                           (A8)

There are no other target prime factors outside2,3,5,p. Unlike
full-coordinate content, binary-block gcds at primes NOT dividing N
are still present and are explicitly retained in(A4) and the tail.

## 3. An elementary divisor constant and binary lattice bound

For every positive k,

    tau(k)<=29k^(1/5).                               (A9)

For primes ell>=32, a+1<=2^a<=ell^(a/5). At each prime below32
maximize(a+1)^5/ell^a. The maximizing exponents, in increasing
prime order2,3,5,7,11,13,17,19,23,29,31, are

    6,4,2,2,1,1,1,1,1,1,1.

Completeness follows because the successive ratio decreases and
(34/33)^5<2; enumeration through32 includes every maximum.
Their exact product is16780370116608000/955049953<29^5.

Let S*(X) count NONZERO integer points with Q<=X. For every X>=1,

    S*(X)<=2X.                                      (A10)

For fixed y, the row count is at most2sqrt((X-23y²/4)_+)+1.
The even function f(y)=sqrt((X-23y²/4)_+) is decreasing for y>=0,
so sum_(integer y)f(y)<=f(0)+integral_R f. The number of nonempty
rows is at most4sqrt(X/23)+1. Subtracting the origin yields

    S*(X)<=2piX/sqrt23+(2+4/sqrt23)sqrtX
         <(176/133)X+(54/19)sqrtX.

Here pi<22/7 and sqrt23>19/4. The latter bound is <=2X when
sqrtX>=21/5, in particular X>=25. The exact finite table for
integer X1..24 is

    2,2,2,4,4,8,8,12,14,14,14,18,18,18,18,20,
    20,24,24,24,24,24,26,30,

which satisfies S*(X)<=2X. Taking the integer part proves the
remaining real X range. No asymptotic lattice estimate is used.

Since r_Q(k)<=2tau(k), a binary block divisible by ell>Y contributes
at most S*(N/ell²) choices for its divided NONZERO block and
2*29*N^(1/5) choices for the other signed block. Taking both block
positions gives at most232N^(6/5)/ell². The zero block is excluded
already by the2/3 selection. No prime>sqrtN is possible. By(A2)
the complete large-prime tail is at most

    (1856/3)N^(6/5)/(YlogY).                         (A11)

## 4. Literal uniform onset for cube atoms

Combining(A7)--(A11) gives the positive primitive23-unit shell count

    S_U(N)>N/2-3584sqrtN
       -(256/3)YsqrtN/logY
       -(1856/3)N^(6/5)/(YlogY).                     (A12)

Set z=N^(7/20), Y=3floor(z). If N>=2^43, then z>=2^15,
Y>=(29/10)z and logY>11. The latter follows from
Y>=3(2^15-1)>(11/4)^11>e^11. The upper bound Y<=3z is immediate.
The two Y-errors in(A12), divided by N^(17/20), are at most

    256/11+640/33=1408/33.

The remaining spectral error is at most3584/z<=7/64.
The integer inequality87^20<=2^129 gives N^(3/20)>=87,
so the main N/2 is at least(87/2)N^(17/20). Finally

    87/2-1408/33-7/64=139/192>0.                    (A13)

Thus S_U(N)>0 for all these targets. This proves the cube onset.

## 5. Complete finite nonprincipal-prime bridge

If p>=2^14 then cp³>=5*2^42>2^43, so Section4 applies.
For p<2^14, the exact input file contains one primitive positive
23-unit coordinate decomposition for each c and each nonprincipal
split prime. There are633 such primes and1899 decompositions.

The standard-library checker independently sieves ALL primes below
2^14 and computes the entire primitive-Q value set up to2^14 using
4Q=(2x+y)²+23y² and complete integer coordinate bounds. A split
prime not in this principal value set is precisely nonprincipal,
since the class group has order three. It checks equality of the
resulting domain with the input domain, then checks each literal
Q(x,y)+Q(u,v)=cp³, gcd(x,y)=gcd(u,v)=1, positivity, and23-unit
condition. PARI qfbsolve was used ONLY to FIND these coordinates;
its algorithm and classification are not trusted by the verifier.

This and Section4 prove N23CubeAtoms at EVERY prime. The certificate
is a finite bridge for a proved prime-parameter onset, not evidence
for an unrestricted extrapolation.

## 6. The remaining principal atoms and their smaller onset

For a principal split p and c notin P+P, division by p preserves
primitive individual membership. Hence a primitive pair at cp can
have neither norm divisible by p. We need the stronger coordinate
version BEFORE the other gcd sieves: in the2/3-selected shell at cp,
the entire p-bad locus is zero. If p divides one norm, it divides
both. Divide each integral binary block by a generator of either
prime above p that divides it. The resulting norm-c pair remains
primitive at2,3, since that division is by local units there.
But the ordinary Q values through15 are{0,1,4,6,8,9,12}; the only
norm pairs for5,10,15 are respectively(1,4),(1,9),(4,6),(6,9)
up to order, and each has a2- or3-divisible block. Zero is itself
divisible by2,3. This contradiction proves p-avoidance without
assuming the original blocks globally primitive.

Now tau(N)<=8, tau(m)=4. The2/3 spectral error is at most1152sqrtN,
and every remaining unit-index prime costs at most32sqrtN. No p
error is needed. The same main>N/2 and tail give(A12) with3584
replaced by1152 and256/3 replaced by128/3. Set z=N^(7/20),
Y=4floor(z), N>=2^40. Then z>=2^14, Y>=(39/10)z, and logY>11.
For the last statement, e<87/32 and(87/32)^11<4(2^14-1);
the elementary factorial-series bound proves e<87/32.
The two Y-errors divided by N^(17/20) are at most

    512/33+18560/1287,

and1152/z<=9/128. Since N^(3/20)>=64, the exact surplus is

    32-512/33-18560/1287-9/128=328385/164736>0.

This proves the principal-prime onset. The remaining principal
atoms cp below2^40 have NOT been certified. In particular this
packet does not complete inert-core5 or the global40-exception list.

Together with the reviewed neutral-atom reduction, the NEW cube
theorem removes the entire nonprincipal-cube obligation: only
cp in U+U for principal split p and c in{5,10,15} remains to close
that particular inert core. The fixed finite-polynomial-template
barrier still applies to those principal atoms. Other inert cores
and the global middle interval remain separate obligations.

Internal proof/replay/primary-source review do not imply external
human review, formal verification, or historical priority.

Exact support: [source ledger](minus23_neutral_cube_atom_source_ledger.json), [literal coordinate inputs](minus23_neutral_cube_atom_witnesses.json), [standard-library checker](verify_minus23_neutral_cube_atom_closure.py), and [receipt](minus23_neutral_cube_atom_closure_verification.json). The analytic conclusion uses the written all-prime argument; the finite packet verifies its complete remaining prime domain.
