# The neutral-atom reduction and exact inert core seven

20 November 2026. Internal written proof with complete exact finite bridges.

Let Q(x,y)=x^2+xy+6y^2, let P be its positive primitive numerical values,
and let U be the values of P prime to23. A split prime is principal when
represented by Q. Put S0={positive integers supported on split primes,
with every nonprincipal split-prime exponent divisible by3}, and
T=S0 union23*S0.

**Theorem CoreSeven.** Among positive N whose inert-prime part is exactly7
and v23(N)<=1, the complete P+P exception list is {21,483}.

**Theorem AtomTwentyOne.** For every principal split prime p,21p belongs
to U+U. For every nonprincipal split prime p,21p^3 belongs to U+U.

The second theorem is equivalent to the first, under the exact phase and
ramified reductions proved below. These conclusions concern one entire
inert core. They do not prove the global forty-exception formula or change
the remaining generic numerical onset.

## 1. Exact finite phase reduction

Use [CN3Individual](PRIMITIVE_NORM_CLASS_NUMBER_THREE.md): a positive
integer n is in P exactly when it has no inert prime factor, v23(n)<=1,
and its number r(n) of distinct nonprincipal split primes with whole
exponent nonzero modulo3 is not1. Multiplication by S0 preserves this
criterion in both directions, and multiplication by23 preserves it on U.
The primes2,3,13,29 are split and nonprincipal. The prime7 is inert
(5=-23 mod7 is a nonsquare).

Write N=7q with q supported on split primes and optional23. Compare

    7q=q+6q=3q+4q.

Let e=(v2(q),v3(q)) mod3, and R count active nonprincipal primes outside
{2,3}. For a summand u supported on {2,3},

    r(qu)=R+#{i:e_i+v_i(u) is nonzero mod3}.

For R>=2 both pairs work. For R=1 at least one works in every phase.
For R=0 both fail precisely at e=(0,1),(2,0),(1,2). Thus the failures
of these TWO DISPLAYED PAIRS are exactly q in3T union4T union18T;
this is not yet a nonmembership assertion about their sums.

The latter two branches have the neutral seeds

    28=1+27=Q(1,0)+Q(1,2),
    126=39+87=Q(3,2)+Q(7,2).

All four values are positive primitive23-unit norms. Multiplication by
T handles every q in4T or18T. Only21T remains. If t in T has a
nontrivial23-unit part, choose a neutral atom a dividing that part.
Every such atom is a principal split prime or a cube of a nonprincipal
split prime. AtomTwentyOne gives21a in U+U; multiplication by the remaining
neutral factor and optional23 proves21t in P+P.

For t=1, the positive primitive norms <=21 are exactly{1,6,8,12,18};
none sum to21. For23M in P+P, reduction modulo23 forces both summands
divisible by23, because -1 is a nonsquare; primitive norms have23-adic
valuation at most1. Division by sqrt(-23) identifies such norms with23U.
Thus23M in P+P iff M in U+U. This excludes483=23*21 too.

Conversely CoreSeven implies AtomTwentyOne: each target21a has inert
part7, has no factor23, and is not21 or483. Its character modulo23 is-1,
since the neutral atom has character+1. If either norm in a P+P pair
were divisible by23, the other would make the target a square modulo23.
Hence both summands are in U. This proves the stated equivalence.

## 2. Imported all-prime shell bounds, with their exact scopes

Use the retained [quaternion bounds(A3)--(A6)](MINUS23_NEUTRAL_CUBE_ATOM_CLOSURE.md)
and [selected-tail bounds(P1)--(P3)](MINUS23_PRINCIPAL_ATOMS_AND_INERT_CORE_FIVE.md). Their logical inputs are
CN3JointSmallPrimes, the maximal-level transfer including ALL norm characters,
the split target-dividing endpoint bound including full pO content, and the
complete rational local-mass certificate. No new spectral import is made.

For N=21*p^t with t=1 or3 and split p outside{2,3}, its2/3 valuations
are0,1. The selected main has fraction(2/3)*(3/4)=1/2 of

    B=(12/11)*4*8*sigma(p^t),   B/N >=128/77.

The 7-bad block-gcd event is zero: v7(N)=1, and the norm is anisotropic
at7. The full23-bad norm event is zero since chi23(N)=-1. Every other
prime must remain in the separate coordinate-block gcd sieve, including
primes not dividing N. The retained piecewise local-mass majorants are

    delta_l=(5l-3)/(l^2(l+1)) for split l,
    delta_l=2/[l(l-1)] for inert l;
    sum_(l>=5,l!=23) delta_l <6/25.

Removing only delta7=1/21 therefore leaves a selected main strictly above

    M*N,  M=(128/77)*(1/2-6/25+1/21)=20672/40425.       (C1)

The selected other-block bound r6(k)<=12*k^(1/6), the nonzero binary
lattice bound S*(X)<=2X, and the actual Rosser--Schoenfeld lower AND upper
prime inequalities give the full large-prime block-gcd tail

    80*N^(7/6)/(Y*logY),  Y>=17.                        (C2)

This includes both possible bad blocks and all large primes. The zero block
is excluded by the2/3 selection. The scale is (C1), not an unjustified
generic O(N) normalization of a sigma(N)-dependent all-target shell.

## 3. Principal atoms: p-avoidance and an onset at2^31

For principal p the entire p-bad event vanishes BEFORE the global gcd sieve.
If p divides one binary norm at target21p, it divides both. Divide each
block by a generator of a prime above p that divides it. These generators
are local units at2,3, so the divided blocks still pass the2/3 selector and
their positive ordinary norms sum to21. The complete2/3-selected ordinary
Q values <=21 are{1,6,8,12,18}, which contain no such pair. To check
completeness,4Q=(2x+y)^2+23y^2 implies |y|<=1; at y=0 selection permits
only x=+/-1, and at |y|=1 the four possible values are6,8,12,18.

The joint selector costs at most1152*sqrtN and each remaining unit-index
bad prime costs at most32*sqrtN, exactly as the retained principal proof:
tau(N)=8, and after stripping the3 factor tau(m)=4. Together with (C1),(C2)
this gives

    S_U(N)>M*N-1152*sqrtN
       -(128/3)*Y*sqrtN/logY-80*N^(7/6)/(Y*logY).       (C3)

Set z=N^(1/3),Y=floor(3z/2). For N>=2^31 the retained inequalities
Y>=149z/100,Y<=3z/2,logY>73/10 hold already from z>=1024.
Also N^(1/6)>896/25, since (28/25)^6<2. Divide(C3) by N^(5/6):

    M*(896/25)-(64+8000/149)/(73/10)-9/8
       =13571626381/12562935000>0.                    (C4)

Thus every principal21p with N>=2^31 lies in U+U.

## 4. Cube atoms: a selected-tail onset at2^34

For t=3 the retained full target-p bound (A5) includes the four-coordinate
content branch. The joint selector error is at most2304*sqrtN, the p error
at most1280*sqrtN, and each unit-index bad prime at most64*sqrtN.
These follow from tau(N)=16,tau(m)=8 and the actual T3 norm at the
small prime. The same (C1),(C2) therefore give

    S_U(N)>M*N-3584*sqrtN
       -(256/3)*Y*sqrtN/logY-80*N^(7/6)/(Y*logY).       (C5)

Set z=N^(1/3),Y=floor(z). For N>=2^34, z>=2048,
999z/1000<=Y<=z, and logY>38/5: the exact rational inequality
(87/32)^76<2048^10 and e<87/32 prove the latter. Also N^(1/6)>50
because50^6<2^34. Dividing(C5) by N^(5/6) gives the positive surplus

    50*M-(256/3+80000/999)/(38/5)-7/4
       =84031919/40923036>0.                         (C6)

This analytic argument excludes p=2,3. They are included in the finite
witness bridge below. Every nonprincipal p>=1024 satisfies
21p^3>=21*2^30>2^34, so no larger prime is left to certify.

## 5. The ACTUALLY completed two finite bridges

For principal atoms the native source verify_core_seven_principal_atoms.cpp
was compiled with g++ -O3 -std=c++17 and run through the INCLUSIVE norm
bound2^31-1. Its coordinate/gcd-row/ellipse enumeration is the same fully
bounded construction as the retained core-five native proof. Every norm
bit is set by an actual primitive Q coordinate; the full list is complete.
An independent signed-coordinate rectangle authenticates all flags through
100000, and ordinary-only values2,3,4 are rejected. The complete prime
sieve through floor((2^31-1)/21) tests EVERY principal split prime and
checks21p against ALL834 primitive23-unit anchors <=5000 until it finds
a valid positive primitive23-unit complement.

The completed run reports:

    inclusive norm bound2147483647; prime cap102261126;
    all primes5884179; principal split primes980245;
    ALL980245 targets covered; failures0; max anchor trials70;
    primitive norm flags261095104;
    positive-y ellipse candidates1406700140;
    primitive positive-y coordinates855199379.

The actual flag stream is256MiB. Its SHA256 is
e3610c2abafb3a5ab3b8597247b1c76b0bb9a2950b964b39033d256ea62e4f7c;
the little-endian (prime,anchor) assignment stream SHA256 is
3000e52e53ac94833c8e938dea8649052dc58f09f124f33a4944517cd9ae3e29.
The native source SHA256 is
dfcc09b690a1b0a26e60158b0603de41ddb7580e01b9fab44a4f9caf9e5504c3.
The complete counts and stream digests are authenticated by the
[replay receipt](core_seven_completion_verification.json). The primary checker
regenerates the full streams in a temporary directory; the large flag stream
need not be retained. Runtime is not a mathematical premise. This is complete target coverage,
not a sample or anchor extrapolation.

For cubes, core_seven_cube_atom_witnesses.json contains exact coordinate
witnesses for ALL57 nonprincipal split primes below1024, including2,3.
PARI qfbsolve was used only to find these coordinates. An independent
standard-library replay reconstructs the complete prime/class domain with
coordinate bounds and trial division, then checks BOTH positive norms,
both gcds1, sum21p^3 and both23-unit conditions. The verifier does not trust
PARI's answer or its enumeration. These completed finite bridges plus
(C4),(C6) prove AtomTwentyOne, and Section1 proves CoreSeven.

## 6. Global boundary

Neutral multiplication preserves the entire inert-prime part. Consequently
closing finitely many inert cores cannot by itself reduce all large
targets to a bounded core list: the distinct cores11^k alone already give
unbounded families. The full forty-exception classification and the generic
averaged large-prime tail remain separate obligations. The false uniform
N/(YlogY) tail is not used here; neither is the already available weak
N^(19/21) error repackaged as a new theorem.

Assurance is internal written mathematical review and exact completed replay.
No external human review, formal verification or historical priority is asserted.

Proof support: [primary replay](verify_minus23_inert_core_seven.py),
[native source](verify_core_seven_principal_atoms.cpp),
[literal cube witnesses](core_seven_cube_atom_witnesses.json),
[exact inputs](core_seven_completion_inputs.json),
[source ledger](core_seven_completion_source_ledger.json), and
[receipt](core_seven_completion_verification.json).

From the repository root run
`python3 research/arithmetic/verify_minus23_inert_core_seven.py --output /tmp/core_seven_review.json`.
