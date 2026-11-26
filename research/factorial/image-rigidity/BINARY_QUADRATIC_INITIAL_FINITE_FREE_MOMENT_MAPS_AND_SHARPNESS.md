# Initial moment maps and sharp support windows for binary quadratics

Corollaries of the separately certified
[initial-window theorem](BINARY_QUADRATIC_INITIAL_ACTUAL_SUPPORT_WINDOW.md)
for every binary-quadratic monomial support. These are initial-window
claims only. The standard commutative-algebra inputs are identified below.

Fix a nonempty subset S of {1,X,Y,X^2,XY,Y^2}, let s=|S|, and use its
ACTUAL monomial coefficients a_i as independent affine coordinates.
Write f_a=sum_(i in S)a_i m_i and H_j(a)=Lambda(f_a^j).

**Theorem 1 (initial finite-free coefficient moment maps).** Over every characteristic-zero field K, the initial moment map

          Phi_S:A^s_K -> A^s_K,
          a |-> (H_1(a),...,H_s(a))

is finite free of rank s!: its coordinate ring is free over the ABSTRACT
target polynomial ring under Y_j -> H_j. In particular it is finite flat
and surjective, and every geometric fiber has scheme length s!.
The zero fiber has the origin as its only underlying point; if s>=2 it
is nonreduced. No reducedness or uniqueness of arbitrary fibers is asserted.
For the full six-coefficient family the rank is720.

The family of finiteness assertions for ALL nonempty S is equivalent to
the initial actual-support theorem. This reverse uses only homogeneity
and full linear coefficient coverage, not Cohen--Macaulayness.

**Corollary 2 (exact-support sharpness and projective length).** On every support S the initial actual-s window is optimal:
there is a phase with EXACT support S over Qbar with moments1,...,s-1
zero and moment s nonzero. It can be chosen over a number field of degree
at most (s-1)!. There are (s-1)! such projective phase classes counted
with scheme multiplicity. This is a scheme-length count, not a claim
of that many distinct phases. For s6 it gives projective length120.

## 1. Zero-only common fibers and graded freeness

Each H_j is homogeneous of coefficient degree j. The certified initial
theorem applied to every actual support contained in S implies that
H_1=...=H_s=0 has only the origin over an algebraic closure: a nonzero
phase with t<=s terms would already contradict its first t zero moments.

Put R=K[a_1,...,a_s] and I=(H_1,...,H_s). The zero-only assertion gives
finite colength. In the Cohen--Macaulay polynomial ring, a homogeneous
system of s parameters is a regular sequence. Its exact multiplication
sequences give

       Hilb_(R/I)(t)=product_(j=1)^s(1-t^j)/(1-t)^s,
       dim_K(R/I)=product_(j=1)^s j=s!.

For completeness the finite-FREE conclusion does not assume algebraic
independence. Let T=K[Y_1,...,Y_s] be an abstract polynomial ring graded
by degY_j=j and acting via H_j. Choose homogeneous lifts b_l of a K-basis
of R/I. Induction on coefficient degree makes

                  direct_sum_l T(-deg b_l) -> R

surjective: subtract the class modulo I, then divide the remainder into
H_j times terms of smaller degree. Source and target have equal Hilbert
series by the displayed quotient series, so this surjection is injective
on every finite-dimensional graded component. Thus R is free of rank s!
over T. Its faithful T-action also proves algebraic independence of H_j.

Tensoring with any geometric target gives a nonzero finite algebra of
dimension s!, hence a nonempty fiber of that scheme length. For the zero
target the support is just the origin. Length s!>1 then forces
nonreducedness if s>=2.

The polynomial-CM and local regularity inputs are precisely
Lemma10.104.7 in the
[retained Cohen--Macaulay section](../../arithmetic/sources/stacks_00N7_cohen_macaulay_2026_11_18.html)
and [Lemma10.104.2](../../arithmetic/sources/stacks_02JN_regular_sequence_2026_11_18.html):
https://stacks.math.columbia.edu/tag/00N7 and
https://stacks.math.columbia.edu/tag/02JN .
The homogeneous graded-basis proof is reproduced above.
The earlier [matching inverse-map application](QUADRATIC_MATCHING_FINITE_FLAT_INVERSE_MOMENT_MAPS.md)
is context only. Neither its fixed-correlation upper theorem nor its
application-specific moment-map conclusions is a premise here. The
zero-fiber input is the separately certified initial-window theorem.

## 2. The genuine restricted converse

Suppose each Phi_S is finite. A hypothetical initial actual-s zero block
of a nonzero phase with support S gives a nonzero point a in its zero
fiber. Homogeneity gives Phi_S(lambda*a)=0 for every lambda after scalar
extension to an algebraic closure. This entire affine line cannot lie
in a finite geometric fiber. Thus initial actual-s nonvanishing follows.
Conversely Section1 derives finite freeness, hence finiteness, from
initial nonvanishing. A single full-six map is not substituted for the
full support-indexed family: actual smaller supports need their shorter
windows. No implication to later starting indices is supplied.

## 3. Sharpness, exact support, and projective length

For s1 take any nonzero coefficient: the first-zero block is empty and
H_1 is nonzero. Suppose s>=2. The s-1 homogeneous positive-degree forms
H_1,...,H_(s-1) in s variables have a nonzero common zero over Qbar,
by the homogeneous common-zero/height theorem (Stacks tag03RB).
Every such nonzero point has EXACT support S: if some coefficient
vanished, its actual support t<=s-1 and its first t zero moments would
contradict the certified initial theorem on that smaller face.
The same initial theorem forces H_s nonzero.

The already established regular sequence H_1,...,H_s has the regular
prefix H_1,...,H_(s-1). The homogeneous quotient by this prefix has
dimension1 and Hilbert series

            product_(j=1)^(s-1)(1+t+...+t^(j-1))/(1-t).

Its Hilbert function stabilizes at (s-1)!. Its projective scheme is
zero dimensional, has no irrelevant torsion (the quotient is a
one-dimensional complete intersection), and has length (s-1)!.
Equivalently this is the usual degree of that projective complete
intersection, derived here from the regular-sequence Hilbert series.

All these forms have rational coefficients. Any closed projective point
therefore has number-field residue degree at most the total projective
length (s-1)!. Normalize a nonzero coordinate to1 to obtain an actual
phase over that field. Exact support and the next nonzero moment have
already been proved for every point. Multiplicities or several Galois
orbits can occur; neither reducedness nor a unique field is claimed.

The imported existence theorem is the exact homogeneous statement in
[retained Stacks tag03RB](sources/stacks-03rb-homogeneous-common-zero.html):
https://stacks.math.columbia.edu/tag/03RB . It is represented by
EXT-STACKS-03RB-HOMOGENEOUS-COMMON-ZERO. The Hilbert-series and
regular-sequence arguments above are standard commutative algebra;
this note claims a consequence for the certified initial quadratic map,
not a new general algebraic theorem.

The mathematical graph inputs are the initial actual-support theorem,
EXT-RECIPROCAL-CM-REGULAR-SEQUENCE, and the exact tag03RB statement.
There is no imported matching upper bound and no later-start assertion.
