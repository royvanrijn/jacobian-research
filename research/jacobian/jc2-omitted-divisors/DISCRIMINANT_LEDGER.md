# Ramification omission and normalization index are separate ledgers

14 November 2026. PF2OMIT4 applies to the source-finite primitive inverse
equation R(T)=sum a_i*T^i of [PF2OMIT3](ADMISSIBLE_INVERSE.md). Its degree-d
projective root cover is X, its finite normalization is N, and its source
open is exactly the finite-T chart of N. Characteristic is zero. Write
D_R for the discriminant of the homogeneous degree-d binary form, rather
than re-taking a lower-degree discriminant after specialization.

## 1. The trace identity separates true ramification from the inverse order

Fix a target prime H. Over its DVR let O be the finite flat order of X
and B its integral closure. Both are free of rank d, with the same
generic algebra. Define the normalization index

    I_H=length(B/O).

If M includes an O-basis into a B-basis, their trace Gram matrices satisfy
G_O=M^t*G_B*M. Smith form over the DVR gives v_H(det M)=I_H. Thus

\[
 \boxed{v_H(D_R)=\sum_{D\mid H}(e_D-1)f_D+2I_H.}          \tag{D1}
\]

The normalized discriminant exponent is the first sum. Indeed, after
unramified residue-field splitting and completion, each tame branch has
uniformizer map h=t^e. Its different exponent is e−1; residue degree f
counts the copies. Trace-pairing foundations are recorded in
[Stacks 0BVH](https://stacks.math.columbia.edu/tag/0BVH) and
[Stacks 0BW0](https://stacks.math.columbia.edu/tag/0BW0).

To identify the raw binary discriminant with the order discriminant
valuation, choose a projective coordinate whose infinity avoids the
finite special fibre. Its leading coefficient is a DVR unit, so the
order is monogenic in that chart. A projective linear change with unit
determinant changes the binary discriminant by a unit. Hence (D1) also
applies when some original T-values are infinite.

Omission has different formulas:

\[
 r_H=d-\deg(R\bmod H)=\sum_{D\text{ omitted over }H}e_Df_D,
 \qquad
 v_H(a_d)=\sum_{D\text{ omitted over }H}p_Df_D,             \tag{D2}
\]

where p_D=-v_D(T)>0. The first equality follows from the source Smith-form
argument over the residue field of H, just as at a closed target. It
counts omitted length. The second counts pole order. Neither is the
normalized different or the normalization index in (D1).

There is a derivative form of this separation. Let c_D be the divisorial
conductor exponent of O in B. In a monogenic local chart S,

\[
 c_D=v_D(F'_S(S))-(e_D-1),\qquad
             \sum_{D\mid H} f_Dc_D=2I_H.                 \tag{D3}
\]

Here F is the locally monic chart equation. The trace dual of O is
F'(S)^(-1)O. An element q in the trace dual of B satisfies q*B contained
in the trace dual of O, and conversely. Therefore the conductor is
F'(S) times the trace dual of B, giving the first formula. Taking the
norm valuation, or comparing the two trace determinants, gives the
second. The different/derivative identification is also supported by
[Stacks 0BWB](https://stacks.math.columbia.edu/tag/0BWB).

On the original finite-T chart this is c_D=v_D(R_T)−(e_D−1). At an
omitted prime use z=1/T and Q(z)=z^d R(1/z). Along the root cover,
Q'(z)=−z^(d−2)*R_T, whence

\[
 \boxed{c_D=v_D(R_T)+(d-2)p_D-(e_D-1).}                  \tag{D4}
\]

The ledger now has three independent divisorial flags: ramification,
omission, and nonnormality of the chosen inverse order. Conductor
exponents here are codimension-one presentation data; they do not replace
individual ambient Artin lengths at exceptional target points.

## 2. What discriminant parity can and cannot identify

For a plane Keller source, if H does not divide a_d, its whole finite
normalization fibre lies in the source. The normalization is etale
there, so v_H(D_R)=2I_H is even. Thus the product of irreducible factors
of D_R having odd exponent is supported on V(a_d).

Under the stronger entire single-ramified-row hypothesis e=2,f=1 with
image C=(H_C=0), equation (D1) gives

\[
 \boxed{D_R=\text{constant}\,H_C\,J^2}                   \tag{D5}
\]

as divisors on the affine target. The odd-exponent part is exactly H_C.
Other factors are normalization-index factors. Unramified omitted
images are detected by a_d and its degree drops; they may be absent
from D_R entirely.

This parity test requires the declared indices and residue degrees.
For F=(x^2,y^2), the two actual ramification primes have (e,f)=(2,2),
and both are retained. The normalization trace discriminant in the
basis (1,x,y,xy) is 256*u^2*v^2. For T=x+y the primitive equation is

    T^4-2*(u+v)*T^2+(u-v)^2,

and its discriminant is 4096*u^2*v^2*(u-v)^2. The inclusion determinant
of its power basis into the normal basis is −4*(u−v). Thus u=v is a
pure index factor and the odd-exponent part is one, despite genuine
ramification over u=0 and v=0. No omitted prime occurs in this finite
map. This non-Keller control prevents dropping residue-degree weights.

## 3. Five generic target groupings for immersed r7

Retain the entire single e2,f1 ramified-prime hypothesis, degree six,
and an immersed degree-(6,10), gap21 image C. PF2OMIT2 gives exactly
three unramified omitted primes, each e=f=1. At most one can share
image C. Let q be that number, and let g_j be the number of remaining
primes sharing each distinct additional image H_j. The exhaustive
groupings are:

| q | Unramified primes grouped by other images | deg(R mod H_C) | Degrees modulo the other images | Lower bounds for raw discriminant orders at the other images |
|---:|---|---:|---|---|
| 0 | 3 | 4 | 3 | 6 |
| 0 | 2+1 | 4 | 4,5 | 2,0 |
| 0 | 1+1+1 | 4 | 5,5,5 | 0,0,0 |
| 1 | 2 | 3 | 4 | 2 |
| 1 | 1+1 | 3 | 5,5 | 0,0 |

In every row v_(H_C)(D_R) is odd and at least 1+2q. Leading-coefficient
orders are at least 1+q at C and g_j at H_j. Exact orders require pole
and conductor data, not only the row counts.

For the extra-image discriminant bounds, split the etale normalization
over a strict complete DVR. The g omitted reciprocal roots z_i all
vanish along H. Each pair contributes 2*v_H(z_i−z_j)>=2 to the
discriminant. Hence v_H(D_R)>=g*(g−1). More precisely the contribution
is at least 2*sum_(i<j) min(p_i,p_j). Other primitive-coordinate
collisions can increase the order. At C with q=1, the local quadratic
ramified factor and the extra linear infinity factor have resultant
divisible by H_C; the product discriminant gains at least two beyond
the normalized different exponent one.

The table supplies the number, images by grouping, e,f and generic
omitted lengths in every numerical type. It does not manufacture the
unknown H_j or a global inverse polynomial R_6.

## 4. Singleton values outside the finite branch locus force high jets

Let the reduced source fibre at zero be the origin, and assume the
finite normalization is etale over zero. In the completed target ring
Ahat=C[[u,v]], it splits into d copies. The retained branch has
T=t0 in the maximal ideal m. Each of the k=d−1 omitted branches has
reciprocal T-value z_j in m. Primitivity makes these functions distinct.
Since a_1(0) is nonzero, the inverse factors formally as

\[
 R(T)=g(T-t0)\prod_{j=1}^{k}(1-z_jT),\qquad g\in Ahat^*.
                                                               \tag{D6}
\]

Let w_j=ord_m(z_j)>0, sorted increasingly. For a Keller boundary, each
omitted point is etale on a smooth boundary curve; its image branch is
smooth and z_j is a unit times its equation to power p_D. Hence w_j
is that omitted prime's pole order, repeated when the prime has several
points over zero.

Elementary symmetric coefficients in (D6) prove

\[
 \boxed{\operatorname{ord}_0(a_i)\ge\sum_{j=1}^{i-1}w_j
       \ge i-1\ (2\le i\le d),\qquad
        \operatorname{ord}_0(a_d)=\sum_{j=1}^{k}w_j.}     \tag{D7}
\]

Indeed the coefficient from T times the product uses i−1 reciprocal
roots, while the term multiplied by t0 uses i roots and has strictly
higher order. For the leading coefficient there is just one product,
so its order is exact. The binary discriminant of (D6) is a unit times

\[
 \prod_{i<j}(z_i-z_j)^2\prod_j(1-t0z_j)^2.
\]

The second product is a unit. Consequently

\[
 \boxed{\operatorname{ord}_0(D_R)
   =2\sum_{i<j}\operatorname{ord}_0(z_i-z_j)
   \ge2\sum_{i<j}\min(w_i,w_j)\ge(d-1)(d-2).}            \tag{D8}
\]

For degree six this requires a_2 in m, a_3 in m^2, a_4 in m^3,
a_5 in m^4, a_6 in m^5, and D_R in m^20. The actual normalized
discriminant is a unit at zero. These raw zeros measure primitive
collisions of five omitted etale points, not ramification.

Under the single e2,f1 ramified-image hypothesis of (D5), zero outside
C makes H_C a local unit. Thus the polynomial index factor J must also
lie in m^10. This is a constraint on an actual inverse equation, not
an assertion that such an equation or global cover exists.

Formal split models z_j=(u+j*v)^w_j, t0=u+7v attain the coefficient
orders and the minimal pairwise bounds. They are local controls, not
connected global Keller covers or source-A2 realizations.

## 5. An off-branch singleton removes two r7 target groupings

In Section3, suppose the source also has a singleton fibre at zero
and zero is outside C. Every unramified prime mapping onto the same
image has the same number b_H(0) of normalization branches there.
The omitted length five therefore requires

\[
                   \sum_j g_j b_{H_j}(0)=5.             \tag{D9}
\]

The q=0,g=3 row is impossible, since 3 cannot divide5. The q=1,g=2
row is likewise impossible. Precisely these three generic groupings
remain compatible with the integer ledger:

| q | Other-image groups | Allowed branch counts over zero |
|---:|---|---|
| 0 | 2+1 | (b_H,b_K)=(0,5),(1,3),(2,1) |
| 0 | 1+1+1 | Unordered (0,0,5),(0,1,4),(0,2,3),(1,1,3),(1,2,2) |
| 1 | 1+1 | Unordered (0,5),(1,4),(2,3); the prime over C contributes0 |

Thus any such degree-six immersed-r7 map with a singleton fibre off
its branch image would need **three or four distinct nonproper target
components**, including C. In the q=1 row one other image must have
at least three branches at zero. A node on an image shared by two
unramified primes contributes two surplus units, rather than one;
prime and target-component counts must remain separate.

These are exclusions of two complete numerical grouping types under
the additional singleton/off-branch hypotheses. They do not exclude
all remaining r7 patterns or all singleton-zero Keller maps.

## 6. A complete marked-root line with 27 pure collision values

Restrict the degree-eight control to q=1,r=0. Its finite normalization
over A1_b is explicitly

\[
 B=C[b,t,a]/(t^4+t-b,\ a^2-4t^3-1).
\]

It is free with basis (1,t,t^2,t^3,a,a*t,a*t^2,a*t^3). Eliminating b
leaves the smooth affine elliptic curve a^2=4t^3+1, so B is normal.
Put Delta=−256*b^3−27. Direct trace Gram determinants give

\[
                    \operatorname{disc}(B)=256\Delta^3. \tag{D10}
\]

For the inverse R_8(T)=sum a_i*T^i, the binary-root order has basis
1 and zeta_j=a_8*T^j+a_7*T^(j−1)+...+a_(9−j)*T, 1<=j<=7.
This is the standard binary-form order basis and projective root
construction in [Wood, Sections2.1 and3](https://arxiv.org/html/1007.5508v1).
The checker expresses these eight elements in the displayed normal
basis. Every coefficient is in C[b], and the inclusion determinant is

\[
 J(b)=C_3(b)C_{24}(b),\qquad C_3=256b^3-b^2+6b+28.       \tag{D11}
\]

The exact C24 and the full inclusion matrix are retained in
[discriminant_equations.json](discriminant_equations.json). Literal
determinants, not only factor counts, verify

\[
                    D_{R_8}=256\Delta^3J(b)^2.           \tag{D12}
\]

Delta,C3,C24 are squarefree, pairwise coprime, and nonzero at b=0.
There are exactly three branch values and 27 pure index/collision values.
The complete affine b-value ledger is:

| Value | Normalization fibre lengths | Retained source lengths | Omitted lengths | Distinct finite inverse T-values |
|---|---|---|---|---:|
| Outside b*Delta*C3*C24=0 | eight1s | eight1s | none | 8 |
| b=0 | eight1s | seven1s | one1, unramified | 7 |
| Delta=0, three values | one4 and four1s | four1s | one4, ramified | 4 |
| C3*C24=0, 27 values | eight1s | eight1s | none | 7, one value shared by two points |

At an index value, the normal cover is etale and J has order one.
Thus exactly one pair of formal T-branches coincides to first order;
all other T-values remain distinct. This conclusion covers every one
of the 27 values, not only a sampled algebraic root.

The projective elliptic completion has one infinity point. Removing it,
the three a=0 points and the omitted t=0,a=−1 point leaves the source
curve with Euler characteristic −5. The complete deficit calculation is

\[
            \chi(U_{\rm line})=8-(3\cdot4+1)=-5.
\]

Counting only ramified omission would give −4 and miss the unramified
loss at b=0. The 27 primitive collisions contribute no source deficit.
This explicit line is a slice of the higher-dimensional marked-root
source, not a plane Keller counterexample.

## Verification and remaining inputs

[verify_discriminant.py](verify_discriminant.py) replays trace matrices,
the binary-order inclusion, exact discriminant factors and gcds, all
five r7 groupings and their singleton refinements, and formal weighted
jet controls. [discriminant_ledger.json](discriminant_ledger.json)
separates normalized ramification, omitted length, pole multiplicity
and index at each kind of row. Normalization, trace-dual and formal
etale-splitting arguments are written proof steps. No actual inverse
R_6 or target equations for an unconstructed remaining candidate are
provided; the full research goal remains open.
