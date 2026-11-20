# Every principal neutral atom and the complete inert core five

20 November 2026. The proof combines an effective selected-block sieve with a complete finite principal-prime bridge.

Let Q=x²+xy+6y², P its positive primitive numerical values, and
U={n inP:23 does not dividen}. A principal split prime means a prime
represented byQ, excluding the ramified prime23.

**Theorem N23PrincipalAtoms.** For every principal split primep and
c in{5,10,15}, cp belongs toU+U.

**Theorem N23InertCoreFive.** Among all positive integers N whose
inert-prime part is exactly5 and whose23-adic valuation is at most one,
the complete P+P exceptions are

    5,10,15,115,230,345.

Here the inert-prime part is the product of ALL inert prime powers inN;
thus N=5q with q supported on split primes and an optional single23.
This theorem settles that entire multiplier cone. It does not classify
other inert cores or establish the global40-exception list.

## 1. Inputs and the already audited principal-atom sieve

Use the exact C3 individual primitive criterion CN3Individual and the
all-valuation2/3 selected quaternion shell, full local bad-prime return
bounds and rational mass certificate of the retained minus23 notes.
The current MINUS23_NEUTRAL_CUBE_ATOM_CLOSURE proves all the
nonprincipal cube atoms cp³. Its Sections2,6 prove the following
facts for a principal splitp, N=cp:

* Both binary blocks selected modulo2,3 give signed shellG. Its error
  from the actual selected main is at most1152sqrtN.
* The5-bad event is zero since5 is inert and v5(N)=1.
* The23-bad norm event is zero since(N/23)=-1.
* Thep-bad event is zero BEFORE global gcd sieving: divide both norm
  blocks by generators abovep, which are local units at2,3, to get
  an impossible2/3-selected ordinary norm-c pair.
* Every other primeell>=5 not dividingN has a bad coordinate-block
  count with error at most32sqrtN. All such primes, including primes
  NOT dividingN, must remain in the sieve.
* After the remaining infinite bad-prime main masses are subtracted,
  the actual selected main is greater than

      (68/125)N,(872/1375)N,(864/1375)N

  for c5,10,15 respectively. Their minimum is(68/125)N.

These statements are imported with their exact written hypotheses;
in particular the p-avoidance does not assume already primitive blocks.
The new improvement below is confined to the large-prime tail.

## 2. Count the already selected other block

Write r_6(k) for the number of signed integer representations ofk byQ
whose coordinate block is divisible by neither2 nor3. For allk>=1,

    r_6(k)<=12 k^(1/6).                              (P1)

Indeed each principal integral ideal has exactly two generators.
At a split primeell of exponenta in the ideal norm there area+1
possible orientation choices. At2 and3 the selected condition says
the ideal is not contained inellR; at most TWO choices remain, one
with all exponent on each of the two prime ideals. Thus, omitting the
principal-class restriction only enlarges the count, and r_6(k)/2 is
bounded by the multiplicative function with prime-power factors:

    2 for ell=2,3 and a>=1;
    a+1 for other split ell;
    1 for inert/ramified ell (if a representation is possible).

For split ell>=64, a+1<=2^a<=ell^(a/6). At2,3, the ratio
2/ell^(a/6) is maximized ata1. All other split primes below64 are
13,29,31,41,47,59; their maxima of(a+1)^6/ell^a also occur at1,
because the first successive ratio(3/2)^6/ell is already less than1.
Inert and ramified prime-power factors contribute a ratio at most1.
Consequently the sixth power of the multiplicative bound constant is

    2^48/(2*3*13*29*31*41*47*59)
       =140737488355328/3986190273 <6^6.              (P2)

This proves(P1). These are only ideal-count UPPER bounds; no claim
that every orientation is principal, or that ordinary norms are
primitive numerical values, is made.

Both Rosser–Schoenfeld1962 bounds from Section3 Corollary1 apply:

    Y/logY <pi(Y),Y>=17;
    pi(t)<(4/3)t/logt,t>1.

Keeping the negative endpoint term in partial summation now gives

    sum_(prime ell>Y)1/ell²
      =-pi(Y)/Y²+2 integral_Y^infty pi(t)/t³ dt
      <5/(3YlogY),Y>=17.                             (P3)

The nonzero binary lattice bound S*(X)<=2X is already proved in the
cube-atom note. If a selected coordinate block is divisible by a
primeell>Y, the other block is STILL selected at2,3. Therefore(P1),
not an unrestricted divisor bound, counts it. Taking the two possible
bad blocks gives a tail at most

    2*(12N^(1/6))*(2N/ell²)=48N^(7/6)/ell²,

and by(P3) the full tail is at most80N^(7/6)/(YlogY). Zero binary
blocks were excluded already by the2/3 selector; ell>sqrtN cannot
occur. No large-prime coordinate-gcd obstruction is omitted.

## 3. An explicit principal-atom onset at2^30

The final signed U+U shell consequently satisfies

    S_U(N)>(68/125)N-1152sqrtN
                 -(128/3)YsqrtN/logY
                 -80N^(7/6)/(YlogY).                 (P4)

Set t=N^(1/3), Y=floor(3t/2). If N>=2^30 then t>=1024,

    (149/100)t<=Y<=(3/2)t,    logY>73/10.             (P5)

The first bound follows from floor(3t/2)>=3t/2-1 and t>=100.
For the logarithm, Y>=1536 and e<87/32; the exact integer inequality
(87/32)^73<1536^10 proves log1536>73/10. The elementary factorial
series bound used in the cube note proves e<87/32.

Divide(P4) by N^(5/6). Its positive main is at least(68/125)*32;
the spectral origin term is at most1152/1024; the two Y terms are
at most(64+8000/149)/(73/10). Their exact difference is

    (68/125)*32-(64+8000/149)/(73/10)-1152/1024
      =1750191/10877000>0.                            (P6)

Hence every principal atom cp withN>=2^30 belongs toU+U. This bound
uses the same all-prime separate-block sieve as the cube proof; a
total quaternion theta series alone is insufficient.

## 4. A completed native finite bridge below2^30

The native source verify_minus23_principal_atoms_finite.cpp was
compiled and run for the INCLUSIVE bound B=2^30-1 and anchor limit5000.
It reports ZERO failures on the COMPLETE remaining principal-prime
domain. The enumeration argument is as follows.

**Primitive norm flags.** Complete the square:

    4Q(x,y)=(2x+y)²+23y².

Negative y adds no new norm values, by simultaneous sign reversal.
At y0, only x=+/-1 is primitive, giving the flag1. For EVERY integer
1<=y<=floor(sqrt(4B/23)), set s=floor(sqrt(4B-23y²)); then the
complete integer x range is

    -floor((s+y)/2)<=x<=floor((s-y)/2).

The implementation adjusts every approximate square root by exact
integer inequalities. On that row it starts with all x permitted,
then sieves precisely those divisible by each distinct prime divisor
ofy. Thus the retained x are EXACTLY gcd(x,y)=1. Each such point sets
the bit for its positive norm. No bit is ever set by any other action.
It follows that the final bitset marks exactlyP intersect[1,B].

**All principal primes and all target coefficients.** An ordinary
Eratosthenes bit sieve lists EVERY primep<=B/5. Such a prime is a
principal split prime precisely when p!=23 and its primitive-norm bit
is set. The source checks all three c5/10/15 separately whenever
cp<=B. Thus no prime beyond the list can contribute a missing target
below2^30, and no prime in the list is sampled or omitted.

**Complete finite representation checks.** The anchor list is EVERY
primitive Q value a<=5000 with23 not dividinga (834 anchors). For each
targetcp the program checks these anchors until it finds a<cp with
the bit for cp-a set and23 not dividing(cp-a). Both values then
belong toU by the exact coordinate construction. ALL targets in the
complete principal domain passed. The anchors are only a device for
finding actual representations; their success is directly verified
for every required target and is not extrapolated from smaller ranges.

An independent direct signed-square enumeration authenticates all
primitive flags through100000, including the ellipse endpoints and
the row-gcd sieve. The ordinary norm4 is explicitly rejected as
primitive. A separate complete2^24 control run agrees with the old
canonical primitive-flag prefix (after accounting for the inclusive
endpoint). These controls support the unrestricted enumeration logic;
the actual target coverage is the full2^30 run itself.

Completed run counts:

    ellipse coordinate candidates, y>0: 703340510;
    primitive y>0 coordinates:427599363;
    primitive numerical norm flags:132325868;
    all primes through214748364:11849083;
    principal split primes through that cap:1974507;
    c5 principal targets covered:1974507;
    c10 principal targets covered:1026431;
    c15 principal targets covered:700565;
    failures:0.

All3701503 finite target cases are covered. Maximum anchor attempts
were76,63,67 respectively. The raw flags and target/anchor assignments
are written in specified little-endian formats and cryptographically
hashed by the wrapper. The counts, source, compiler command and final
digests are retained in the scratch replay record. Runtime is not a
mathematical premise. This finite bridge plus Section3 proves
N23PrincipalAtoms.

## 5. The exact neutral reduction, with both directions audited

Let S0 consist of positive integers having no inert prime factor,
no23 factor, and every nonprincipal split-prime exponent divisible
by3. Its atoms are principal split primes and cubes of nonprincipal
split primes. Put T=S0 union23S0. By CN3Individual, multiplying by
S0 preserves primitive numerical membership in BOTH directions: it
changes no active nonprincipal exponent residue, inert condition or
ramified exponent. A further single23 preserves membership onU.

For arbitraryq supported on split primes and an optional23, compare
the two scaled pairs5q=q+4q=2q+3q. Set e=(v2q,v3q) modulo3 and
R=number of active nonprincipal prime powers outside{2,3}. Each
final summand is primitive exactly when its active count is not1.
For R>=2 both pairs work. For R1 at least one pair works at every
phase. For R0 the phases where BOTH pairs fail are exactly

    (0,0),(1,0),(0,1),(1,2).

Equivalently the missing multipliers are exactly T union2T union3T
union18T. This is a statement about these two scaled pairs, not a
nonmembership assertion about their targets. The whole90T branch
is covered by90=12+78=Q(2,1)+Q(8,1), then neutral multiplication
and optional23 lifting. The only residual families are5T,10T,15T.

If t inT has nontrivial23-unit part, that part has some neutral atom
a. Both N23PrincipalAtoms and the retained complete cube-atom theorem
provide ca inU+U. Multiplication by the remaining neutral factor and
optional single23 gives ct inP+P. This works for EVERY split support
and exponent; the multiplier itself need not be a principal norm.

If t is1 or23 the possible failures are5,10,15 and their23 multiples.
The primitive norms through15 are exactly{1,6,8,12}; none of5,10,15
is their sum. Finally23M inP+P iffM inU+U: modulo23 a vanishing sum
of two square norms forces both divisible by23, while a primitive
norm has valuation at most one; division bysqrt(-23) identifies such
norms with23U. The converse follows by multiplication bysqrt(-23).
This excludes115,230,345 as well and proves N23InertCoreFive.

Conversely an exact core5 classification with these six exceptions
would imply every neutral atom ca lies inU+U: it has inert part5,
is larger than15 and prime to23, and its character modulo23 is-1,
so neither summand can be ramified. This supplies the actual reverse
for the neutral-atom criterion; it does not assert anything about
other inert cores or global exception completeness.

Only internal mathematical review and exact computation are claimed.
Independent replay, external human review, formal verification and
historical priority must retain their separate scopes.
