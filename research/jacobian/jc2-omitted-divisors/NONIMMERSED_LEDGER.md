# Complete first-image losses across immersed and nonimmersed r7 charts

14 November 2026. PF2OMIT13 removes the immersion restriction from the
first-prime matrix accounting and adds the missing oriented residual row.
The global assumptions are those of PF2R7RAM2: a degree-six plane Keller
candidate containing the specified index-two, degree-(6,10), gap-21 prime,
with the named geometry and single-prime exclusion inputs retained.
This does not derive those assumptions for every terminal passport.

## 1. The exact cusp values in all coefficient charts

PF2R7FAM1/2 give two complete coefficient charts. The nonimmersion
parameters have precisely these normalized target cusp coordinates:

| Description | Parameter equation | Critical t* | Cusp U* | Cusp V* | det(gamma'',gamma''') at t* |
|---|---|---|---|---|---|
| L, class I | 3d^3+2=0 | -3d^2/2 | -1/4 | -49d^2/5832 | 6048d |
| K, class I | 3x^3+2=0 | -3x^2/2 | -1/4 | -49x^2/5832 | 6048x |
| K, class II | 15x^3-4=0 | -15x^2/4 | 1/48 | 686963x^2/58320 | -6300x |

The first two descriptions are the same coefficient points at L cap K,
with d=x. There are two classes, not three. Coordinates are in the
PF2R7FAM1 normalized target; the transformation from an unconstructed
global target remains unspecified.

These formulas are checked by literal polynomial reduction modulo the
displayed parameter equations. The checker divides p',q' by t-t* and
checks their residual resultant is a unit in each parameter algebra.
Hence this is the only common derivative root at every parameter in
each row. The displayed jet determinant is a unit, so it is an ordinary
A2 cusp. Dividing p-U*,q-V* by (t-t*)^2 gives a second unit-resultant
test: the cusp has no other normalization parameter over its target.
These are uniform identities over the squarefree cubic algebras, not
samples of algebraic roots. Their exact residual resultants are in the
[ledger](nonimmersed_ledger.json).

## 2. Correct the image norm matrix at a cusp

Let m(y)=6-rank M(y), with the explicit PF2OMIT11 matrix of q-V on the
monic p-U algebra. At every parameter, m is the length of the
parametrization fibre. It becomes a branch count only on the immersed
domain. For a nonimmersed coefficient point define epsilon(y)=1 at
its listed cusp value y*, and zero elsewhere. For an immersed
coefficient point put epsilon=0 everywhere.

The unique critical parameter is the cusp parameter, whose fibre is
exactly C[t]/(t-t*)^2. Every other parametrization fibre is reduced:
a multiple gcd root would make both p' and q' vanish there. Therefore
the normalization point count and the first-prime omitted ambient
length are, on **all** the coefficient charts,

    c(y)=m(y)-epsilon(y),
    L_R(y)=2c(y)+epsilon(y)=2m(y)-epsilon(y).          (C1)

PF2OMIT7 supplies the ambient length3 at the cusp and length2 at every
other R point. At the cusp the matrix has corank2 but there is one
normalization point and ambient length3. At a node it also has corank2,
but there are two points, each of ambient length2. Counting both cases
as two omitted folds would incorrectly lose one extra sheet at the cusp.

This completes the first-prime image/loss formula beyond the previously
immersed domain. It does not supply the other target image.

## 3. Complete oriented ramification/omission rows

For either nonimmersed class, PF2R7CUSP1 gives one ordinary A2 cusp
and no other nonimmersion: class I is A2+A5+8A1, class II A2+11A1.
PF2OMIT7 thus gives E_R=deg gcd(p',q')=1. The additional ramification
is forced by PF2R7RAM2. The entire remaining Euler budget is

    5=(e_R+E_R)+(e_S+E_S)+t=3+(e_S+E_S)+t.

Since e_S>=2, necessarily e_S=2,E_S=0,t=0. Its f is one and its
image B is immersive by PF2OMIT7. In particular B differs from the
cuspal C. There is no unramified omitted prime in this row, and both
ramified primes are omitted. The affine normalization is smooth since
its source open is smooth and every boundary point has index2.

Together with PF2OMIT6 this gives the full four-row ledger relative to
the distinguished gap-21 prime R:

| Oriented row | R over C: (e,f,E) | S over B: (e,f,E) | Unramified omitted prime | Nonimmersive image |
|---|---|---|---|---|
| A | (2,1,0) | (2,1,0) | One (1,1,0) over H | None |
| B | (2,1,0) | (2,1,1) | None | B, cusp value still unknown |
| C | (2,1,0) | (3,1,0) | None | None |
| N | (2,1,1) | (2,1,0) | None | C, cusp value listed above |

Rows B,N are one unordered ramification pattern with different locations
of the distinguished gap-21 image. They must remain distinct when
claiming which image equation or exceptional value is actually known.
The first three rows retain all their prior image-sharing constraints.
No other oriented row is permitted under these hypotheses.

## 4. Every exceptional source length in the new row

Let b(y) be the normalization branch count of the still unknown B.
The complete pointwise formula is

    #F^-1(y)=6-2c(y)-2b(y)-epsilon(y).               (C2)

Retained factors have length1. The classified nonimmersed C has only
ordinary cusp and double-image values. Its full conditional table is:

| Location on C | c | epsilon | Allowed b | Retained lengths |
|---|---:|---:|---|---|
| Outside C | 0 | 0 | 0,1,2,3 | 6,4,2,0 |
| Smooth C value | 1 | 0 | 0,1,2 | 4,2,0 |
| Double C value, including A5 | 2 | 0 | 0,1 | 2,0 |
| Its unique cusp y* | 1 | 1 | 0,1 | 3,1 |

At an A5 value there are two separate omitted fold points, not an
additional length-two excess on one point. The only R excess occurs
at its cusp. If B passes through y* it must have exactly one smooth
branch there; if it avoids y*, the source fibre there has three points.
All other source lengths are even. Consequently the **only** possible
singleton value is the explicit y*, and it is a singleton exactly when
B contains that point. Equivalently its reduced equation must satisfy
H_B(y*)=0, with nonzero target gradient there. The target branches need
not be transverse. A proposed singleton-zero map must send its zero
to the listed point under its actual target normalization; that
normalizing affine transformation is still required.

If the singleton occurs, its omitted factors are C[w]/w^3 at R and
C[w]/w^2 at S, plus the retained C. The PF2OMIT7 radial packet is
(3,1)+(2,1), with per-direction lengths3,2 and a retained (1,1) packet.
The generic original indices are2,2, not3,2.

## 5. Inverse and monodromy obligations

For an admissible primitive source coordinate T from PF2OMIT3, with
pole orders p_R,p_S>=1, the complete leading coefficient and raw
discriminant have forms

    a6=constant*H_C^p_R*H_B^p_S,
    Disc(R6)=constant*H_C*H_B*J^2.                   (C3)

At the explicit cusp with B absent, the retained degree is3, so a3
is a unit there. The three reciprocal infinity roots have valuation
2p_R/3 on a generic target line by PF2OMIT10. Hence

    ord a_(3+j)>=ceil(2j*p_R/3), j=1,2,3,
    ord a6=2p_R,
    ord Disc(R6)>=4p_R, ord J>=2p_R-1.              (C4)

If B contains the cusp, the remaining fold has reciprocal valuation
p_S/2 twice. The PF2OMIT10 singleton row B now applies with the roles
of the distinguished images specified: the minimum coefficient row
is(1,1,2,3,3), raw discriminant order>=11 and index order>=4. With
the actual poles the raw lower bound is
4p_R+p_S+12 min(2p_R/3,p_S/2), normal different order3.

Finally the new row also has monodromy S6. Repeat PF2OMIT8's argument:
the connected finite cover over A2 minus C union B is transitive,
branch meridians normally generate because the ambient A2 is simply
connected, and every generic inertia is a single transposition. The
graph of all conjugate transpositions is connected and generates S6.
The S5 source stabilizer is maximal, so there is no proper intermediate
source field. This argument uses generic inertia and is unchanged by
the cusp in the image; it does not assert etaleness at C intersect B.

The [checker](verify_nonimmersed_ledger.py) replays all cusp identities,
cofactor-resultant units, row budgets and exceptional length tables.
Global finite-flat conservation, the single-prime exclusion, boundary
regularity and topology are written/imported steps. Actual B, the
original target normalization and source realization remain unknown.
No whole-r7, unrestricted singleton-zero properness or JC2 conclusion
is claimed; formalization and external review are absent.
