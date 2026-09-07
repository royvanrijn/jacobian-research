# Generic chord blocks have private good-prime ramification

Every retained generic chord and tangent on the fixed 103b2 high/low pair
and historic ICARM356/398 has a private good-prime obstruction. Consequently
**no nonzero combination** of their norm-projected relation roots becomes
unramified outside S, even after arbitrary generic corrections. This is
an exact exclusion of the whole supplied constructor, not an individual
point-search miss.

Local three-root cancellations do occur in the dictionary. Their existence
does not create an unramified block because independent private obstructions
survive elsewhere. This directly tests a simultaneous-cancellation mechanism
at the incidence stage, before rational solubility or CT.

## A constructor outside simple multiplication of generic classes

On E:y²=f(x)=x³+Ax+B, let L(x)=s*x+h be the chord through P,Q, or the tangent
at P. If R is its third intersection, then

    f(X)-L(X)^2=(X-x(P))(X-x(Q))(X-x(R)).

At a cubic root theta this gives

    L(theta)^2=(x(P)-theta)(x(Q)-theta)(x(R)-theta).

Taking this explicit relation root is not multiplication in the F2 span
of the generic Kummer classes. Normalize L to the primitive integral
polynomial a+bX with b>0, and set

    alpha=a+b*theta,
    n=Norm(alpha)=a^3+A*a*b^2-B*b^3,
    eta=n*alpha,                 Norm(eta)=n^4.

Changing the linear form by a rational scalar multiplies eta by a rational
fourth power, so this normalization does not change its squareclass.
Eta is a norm-square H1 candidate, not automatically a Selmer class.
This differs from testing whether the x-axis intercept of a secant is
a rational elliptic point: no new point search is performed here.

The [frozen protocol](GENERIC_CHORD_RAMIFICATION_PROTOCOL.json) includes
every unordered generic pair with both relative signs, plus each tangent.
For m generic points this gives m² forms. On these four cases every form
is nonconstant and distinct. The input contains the equation, generic
points and certified bad support only; ranks are joined below solely
for interpretation.

## Exact paired outcome

| Frozen fibre | Retained rank / generic rank | Gain lower bound | Forms | Outside-S parity rank | Private witnesses | Possible extra Selmer dimension in this span |
|---|---:|---:|---:|---:|---:|---:|
| 103b2, 3726/881 | 27 / 17 | +10 | 289 | 289 | 289 | 0 |
| 103b2, -1049/2296 | 17 / 17 | 0 observed | 289 | 289 | 289 | 0 |
| ICARM356 | 29 / 17 | +12 | 289 | 289 | 289 | 0 |
| ICARM398 | 30 / 16 | +14 | 256 | 256 | 256 | 0 |

The rows are case-02,03,13,15 in the frozen panel. Full ranks remain
UNKNOWN, and the low row is a censored search control. All 1123 forms
are checked; no adaptive enlargement follows the zero kernels.

Writing D for the span of the displayed eta classes, the certificates prove

    dim D = m²,
    D intersect G = 0,
    Sel2(E) intersect (G+D) = G.

The constructor readily creates independent *ramified* H1 classes. It
creates none of the unramified directions required by the jump theorem.

## Why root supports, rather than norm parity alone, matter

Let p be an odd prime outside S, so f is integral and etale at p. If
p divides n for a primitive linear form, p does not divide b. The form
vanishes only at the degree-one prime corresponding to

    theta = -a/b mod p.

Its valuation there is v_p(n), and all its other prime-ideal valuations
are zero. In the three geometric root coordinates, eta therefore has
parity v_p(n)*(1+e_r), where e_r selects that root. For nonsplit p the
remaining conjugate coordinates have the same value.

For three distinct roots at a split prime,

    (0,1,1)+(1,0,1)+(1,1,0)=0.

Thus three candidates with odd norm valuations can cancel their
prime-ideal ramification even though their product norm has odd valuation.
Rational norm parity alone is not a valid exclusion for a product block.
The exact local regression f=X³-2 at p=31 uses roots 4,7,20; each linear
norm has valuation one, and the displayed triple cancellation occurs.
This is a local regression, not a global Selmer or rationality claim.

The worker builds a pairwise-coprime basis of the norms after removing S.
For each atom it groups the linear forms by their root residues, refining
by gcds of pair resultants if necessary. Root differences are checked to
be zero or units, so every prime in an atom has the same partition.
Perfect-square atoms impose no parity constraint; a nonsquare atom
contains at least one odd-valuation prime and supplies the displayed
constraint without identifying or factoring that prime.

The four rows contain 7,9,12,8 atoms respectively with three active root
groups. These local cancellations survive the grouping calculation,
but they cannot defeat the following private witnesses.

## A simpler independent proof of the zero kernels

For every form i, the verifier finds an integer c_i>1 with:

1. c_i is nonsquare and coprime to S and disc(f);
2. n_i=c_i^e*u with e odd and gcd(c_i,u)=1;
3. gcd(c_i,n_j)=1 for every other form j;
4. its root -a_i/b_i modulo c_i is simple.

Choose any prime p dividing c_i to odd order; one exists because c_i
is nonsquare. Then v_p(n_i) is odd. Eta_i has nonzero parity at the
other cubic root factor(s), whereas every eta_j for j!=i is a unit
there. Generic Kummer classes are unramified at this good odd place.
Any product selecting form i therefore fails the Selmer condition,
irrespective of all other selections and generic corrections. Applying
this argument to each i proves every asserted independence and exclusion.
The good-prime unramified condition is the standard descent condition;
see [Schaefer–Stoll](https://mathe2.uni-bayreuth.de/stoll/papers/p-descent-long.pdf).

No primality assertion is made about c_i. The private witnesses have bit
lengths 66–736,65–779,84–898,42–863 on the four rows. Exact coprimality,
integer division and nonsquareness suffice; no large integer factorization
or prime-ideal enumeration is required.

## Certificates and implications

The [arithmetic artifact](../../artifacts/generated-results/elliptic-curves/rank_jump_generic_chord_ramification_v2.json)
retains every form, source pair, norm, gcd atom, root partition and parity
row. The [independent verification](../../artifacts/generated-results/elliptic-curves/rank_jump_generic_chord_ramification_verification_v2.json)
reconstructs the short model from the original Weierstrass coefficients,
checks complete chord/tangent coverage, and computes norms as determinants
of cubic multiplication matrices. It then verifies all 1123 private
witnesses without using the producer's parity elimination or root partitions.
An initial verifier typo in the tangent identity was corrected before its
first valid artifact; the producer and arithmetic artifact were unchanged.

A separate regression exposed an overlapping-power bug in the initial
gcd-refinement routine: inserting 343 after 49 must split the existing
atom after stripping its first power. The corrected v2 producer replays
exactly the frozen inputs and independently checks pairwise coprimality.
All four arithmetic rows are identical to v1, and the independent private
witness verification passes again. Both versions are retained for provenance;
v2 is the reusable route. Three regressions cover this repair, the genuine
three-root cancellation, and its destruction by private obstructions.

The experiment rules out this specific generic relation-root dictionary
on successful fibres as strongly as on the matched low. It does not
exclude other relation roots or uncomputed ideal/unit constructions.

For a surviving block, all private odd good-prime supports must disappear
or be shared by compatible root supports. Three-root collisions at a few
primes are insufficient. A mechanism stated in terms of t must force
that *global* cancellation and then prove independence and local S
admissibility. Here the actual successful specializations do not do so
for these forms, so this dictionary cannot explain their jump.

No solubility or visibility feature is inferred, and no candidate-selection
change is supplied to Agent 1. The next independent constructor needs
shared even-valuation ideal relations from the outset; simply adding more
generic chord roots is not justified by these results. The required
unramified block and a sufficient condition on t remain unconstructed.

```sh
timeout 60 sage -python elliptic-curves/rank-jump/verify_generic_chord_ramification_v2.py check
python3 -m unittest discover -s elliptic-curves/rank-jump -p test_generic_chord_ramification.py
```
