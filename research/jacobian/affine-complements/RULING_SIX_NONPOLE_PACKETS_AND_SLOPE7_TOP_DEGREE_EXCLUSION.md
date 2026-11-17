# Nonpole packets and the top slope-seven degree-30 exclusion

PF2RULED183. 17 November 2026. Actual complex polynomial plane Keller maps,
in the complete marked quadratic obstruction class.

A single nonpole geometric sheet cannot remain alone in an actual A
quartic cluster. Applied to the remaining degree-D=30 class, this forces
ALL FOUR quartic sheets to be poles. The intermediate degree nine is
then impossible by the exact quartic Newton exponents and parameter
precision.

Together with PF182, every remaining D=30 survivor must satisfy

    deg M=2, deg h=14, one actual A support of slope L=7,
    M=constant*u^2, h=constant*(u-c)^14, c!=0,
    [C(u,w):C(P,Q)] in {6,7,8},
    ord_c e in {2,3,4,5}.                               (1)

These are necessary constraints. Explicit formal quartics for ALL THREE
remaining degrees pass the valuation and affine-response tests used here;
none is asserted to have an actual Keller mate. The subsequent [two-moment theorem](RULING_SIX_TWO_MOMENT_RESONANCE_AND_DEGREE30_CUBIC_REDUCTION.md)
excludes degrees six and seven. Its remaining degree-eight class has
companion ruling degree twenty-seven by the further prefix theorem;
that class, HC4 and JC2 stay open.

## 1. A nonpole packet contains at least two geometric sheets

Retain the whole source, actual polynomial companion Q and original
full normalizations. At an off-M A support put v=u-c and work over
an algebraic closure k of C(p). On a geometric generic P branch,
choose a fixed Puiseux parameter s with v=s^alpha. The source function
w has a pole and a finite first truncation

    w=w0+s^r*(lambda+higher terms),
    w0 in C[s,s^-1], lambda not in C.

The full polynomial-leading-coefficient lemma of PF159 gives

    H(lambda)=p, H in C[T], deg H=h>=1.                  (2)

The polynomial H(T)-p is irreducible over C(p): the extension
C(T)/C(H(T)) has degree h. It is separable. Each of its h distinct
conjugates lambda_i lifts uniquely by the ordinary implicit recursion,
since H_T(lambda_i) is nonzero. The resulting Puiseux series share
w0 and the index r and are distinct roots of the same entire equation
P=p, inside the SAME A quartic cluster.

This counts geometric ROOTS, not normalized points or residue-field
generators. Some roots can be related by s -> zeta*s and belong to
one ramified point; they still count as distinct geometric sheets.
No assertion that alpha*h is a residue degree is used.

The extended partial_p derivation commutes with these conjugations.
Since partial_p lambda_i=1/H_T(lambda_i) is nonzero, every conjugate
has the same exact order r of partial_p w. Hence all have the same
normalized coframe pole/holomorphic behavior. In particular, if the
chosen branch is a nonpole of Q, then all h conjugate sheets are
nonpoles. Here dQ is holomorphic exactly when the rational Q is finite.

If h=1, the first response is affine in p and its derivative is a
nonzero COMPLEX constant. The entire implicit recursion is then over
C[p] after the fixed Puiseux substitution, with no division by a
nonconstant polynomial in p. The actual Q evaluates in C[p]((s));
holomorphic dQ makes its residue a polynomial q0(p). The original
source escapes because w is polar and both t lifts over c are units.
The original full-normalization argument of PF165 therefore gives
the WHOLE closed graph Q=q0(P), forbidden by the retained Chau
component theorem. Thus h=1 is impossible for such a nonpole.

Consequently any nonpole packet contains at least TWO geometric
sheets. In particular the total nonpole sheet mass of an A quartic
cannot equal one. This lemma does not require all nonpole packets
to have even size or identify h with a closed-prime residue degree.

## 2. All remaining degree-30 sheets are poles

PF182 excludes all slope-six D=30 strata and leaves exactly

    L=7, b_field=[C(u,w):C(P,Q)] in {6,7,8,9}.

There is one actual supporting base, and all Q poles occur in its
quartic cluster. PF164 gives the strict pole-order bound

    q_Q<(5/2)*alpha.

Pole sheet mass at most two would give total pole degree below
five, contradicting b_field>=6. Pole mass three would leave precisely
one nonpole sheet, forbidden by Section1. Thus ALL FOUR sheets are
poles in every remaining D=30 stratum.

Write the critical depressed quartic as

    G=X^4+2J*X^2+q*X+J^2+Delta,
    t=ord q=7+ord_c e,
    ell=ord J in positive integers or infinity,
    delta=ord Delta>=t+1.

PF182 proves the exact parameter cap delta<=14. On each sheet put
mu=ord_v G_X. The whole coframe formula and strict affine-response
bound give

    8<mu<21/2,
    sum_(four sheets) mu=32+b_field.                    (3)

The derivative trace sum G_X(root)=-8q gives at least one mu<=t.
Since all four mu exceed eight, t>=9. PF178 gives ord_c e<=5.
Therefore every remaining D=30 case satisfies

    2<=ord_c e<=5.                                      (4)

## 3. Degree nine permits only two exact pole patterns

Suppose b_field=9. Positive integral pole orders satisfy
q_Q<(5/2)*alpha and their indices sum to four. Exhausting the
integer partitions of four gives ONLY

    one index-four point with Q pole order nine,
    OR indices (3,1) with Q pole orders (7,2).

Their exact derivative-valuation multisets are respectively

    [41/4,41/4,41/4,41/4],
    [31/3,31/3,31/3,10].                                (5)

The derivative trace now gives t>=10, so the complete parameter range is

    t in {10,11,12}, t+1<=delta<=14.                    (6)

The coefficient orders ell and delta are integral before taking roots.

## 4. The Newton mechanisms exclude both patterns

Use the arbitrary-coefficient quartic facts of PF180. They are exhaustive
in range (6). The possible equality delta=4t/3 is nonintegral for
t=10,11 and equals sixteen, outside the cap, for t=12. The possible
equality ell=2t/3 is likewise nonintegral for t=10,11. At t=12,
ell=8 has delta<=14<16, so the dominant-Delta mechanism applies.
All infinite-J cases are included in these comparisons.

A dominant Delta gives four roots with mu=3delta/4; the square/Delta
boundary gives AT LEAST TWO with that value, retaining its possible
double-zero initial. For integral delta this number is neither
41/4 nor 31/3 nor 10: the required delta values would be 41/3,
124/9 and 40/3. Neither multiset (5) can occur there.

Cubic/linear balance gives all four mu=t, also incompatible with (5).

It remains the square centers. Here ell<2t/3 and delta>2ell.
Put m0=min(delta,t+ell/2). Then

    m0>2ell, m0<2t-ell.

The full square-center calculation gives at least TWO roots with

    mu=(m0+ell)/2.                                      (7)

When delta=t+ell/2, cancellation at one sign can leave only the
other pair; the two guaranteed roots are enough below.

If m0=delta, the value (7) is a half-integer. It cannot be 41/4
or 31/3. If it equals ten, it occurs on at least two sheets, whereas
the second multiset (5) contains only ONE ten. Thus both patterns
are excluded.

If m0=t+ell/2<delta, equation (7) becomes

    mu=t/2+3ell/4.

It belongs to one quarter of the integers and cannot be 31/3.
Again a value ten on at least two sheets is incompatible with the
second pattern. For the first pattern it would have to satisfy

    2t+3ell=41.

Among t=10,11,12 the only integral candidate is t=10,ell=7.
It violates the required square-center inequality ell<2t/3.
Hence the first pattern is impossible as well.

This covers every coefficient order and all one-sign cancellations.
The degree-nine stratum is therefore excluded, proving (1).
No degree cutoff on Q or on later source coefficients was used.

## 5. Formal controls for all three remaining degrees

The following are FORMAL depressed quartics over C(p)((v)). All have
L=7, the correct nonzero affine-p coefficient at order v^14, and every
sheet satisfies 8<mu<21/2. The indices and pole orders below are the
ones forced by the coframe test, not assertions that the whole
coframe has an exact rational primitive or an actual Keller mate.

### Degree six

Take

    G=(X-v^3)^3*(X+3v^3)+v^13-p*v^14.

Its exact coefficients have

    J=-3v^6, q=8v^9,
    Delta=-12v^12+v^13-p*v^14,
    ell=6, t=9, delta=12.

The simple root X~ -3v^3 has mu=9, index one and Q pole order one.
Indeed substitution X=v^3*Z gives initial (Z-1)^3*(Z+3), with
nonzero derivative -64 at Z=-3. Its first p response is affine
at order v^5.

At the triple center write X=v^3+U. Put v=s^3,U=s^10*Z.
The ENTIRE equation after division by s^39 is

    4Z^3+1+s*Z^4-p*s^3.

Its three constant nonzero leading roots are simple. They lift over
C[p][[s]], have mu=29/3, and form one index-three point with
Q pole order five. Their first p response in U is at s^13=v^(13/3)
and is affine. Thus the total pole degree is 5+1=6, and the
four-sheet derivative sum is 3*(29/3)+9=38=32+6.

### Degree seven

Take

    G=(X-v^3)^3*(X+3v^3)+v^10*(X-v^3)-p*v^14.

Here

    J=-3v^6, q=8v^9+v^10,
    Delta=-12v^12-v^13-p*v^14,
    ell=6, t=9, delta=12.

The simple root X~ -3v^3 again has mu=9, index one, pole order
one and affine first p response at v^5.
At X=v^3+U the entire equation is

    4v^3*U^3+U^4+v^10*U-p*v^14.

The outer substitution v=s^2,U=s^7*Z, divided by s^27, gives

    4Z^3+Z+s*Z^4-p*s.

The two nonzero simple roots are Z=+i/2,-i/2. They have mu=10
and form one index-two point with Q pole order four. The first
p response is affine in U at s^8=v^4.
The inner substitution U=v^4*Z gives

    v^-14*G=Z-p+4v*Z^3+v^2*Z^4.

Its simple affine first response Z=p has mu=10, index one and
Q pole order two. Thus the total pole degree is 4+2+1=7 and
the four-sheet derivative sum is 3*10+9=39=32+7.

### Degree eight

Take

    G=X^4+2v^7*X^2+v^10*X+(1-p)*v^14.

It has t=10,ell=7,Delta=-p*v^14 and delta=14.
Its outer substitution v=s^3,X=s^10*Z, divided by s^40, gives

    Z^4+Z+2s*Z^2+s^2*(1-p).

The three nonzero roots of Z^3=-1 are simple and constant.
They lift over C[p][[s]] and have mu=10. Their first p response
in X is at s^12=v^4 and is affine. Together they form one
index-three point with Q pole order six.
The remaining substitution X=v^4*Z gives

    v^-14*G=Z+1-p+2v*Z^2+v^2*Z^4.

It has a simple affine first response Z=p-1, mu=10 and index one,
with Q pole order two. Thus the total pole degree is 6+2=8 and
the four-sheet derivative sum is 4*10=40=32+8.

None of these controls is asserted to be a sextic source realization,
a full-companion resonance solution, an exact rational primitive
for the entire coframe, or a constructed Keller map. They show that
all THREE remaining degrees require additional global or coefficient
constraints, rather than a continuation of just the tests used here.
The additional full-resonance theorem PF184 excludes degrees six and
seven. Degree eight remains open, with companion ruling degree
twenty-seven forced by PF185. These controls remain guards for the
strictly weaker local tests verified here.

## 6. Exact controls and proof boundary

The checker verifies conjugate counting on representative full first
responses, irreducibility over C(p), exact derivative/coframe orders,
the complete index-four pole partitions, every Newton order comparison
in (6), the two target multisets and the formal degree-six, seven and eight
quartic substitutions, their normalized pole indices and affine responses. The written irreducibility, Hensel and graph
arguments prove the nonpole-packet lemma for every response degree.
The Newton comparisons prove the degree-nine exclusion for every
later coefficient choice; finite controls do not replace them.

Replay:

    .venv/bin/python research/jacobian/affine-complements/verify_ruling_six_nonpole_packets_slope7_top_degree.py

Essential inputs are PF159, PF164, PF165, PF178, PF180, PF181 and
[PF182](RULING_SIX_PARAMETER_PRECISION_AND_SLOPE6_DEGREE30_EXCLUSION.md).
All original full-normalization primes remain retained. No target
landing, trace integrality, normal-tail termination, HC4 or JC2
closure follows. Independent whole-proof replay, formal verification
and external human review are not claimed.
