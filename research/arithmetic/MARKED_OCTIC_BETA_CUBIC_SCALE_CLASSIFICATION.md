# The cubic scales of the fourth marked octic norm

The fourth norm function in
[MO3](QUADRATIC_MIRROR_CLASS_RANK_FOUR_FROM_MARKED_OCTIC.md)
forces a finite list of cubic parameter scales. Exactly
four rational cube classes allow it to be both virtual
and locally a cube at three. Reciprocity pairs them into
two rational isomorphism classes of genus-eleven covers.
Both support all four marked classes on explicit slices
where every bad prime of the three elliptic factors splits.

## Exact valuation restrictions

**Proposition BS1.** Use \(R^2=F(t)\) and \(\beta\) of
MO3, with \(t\in\mathbb Q^\times\) giving a quadratic
field \(K\). Put \(k_p=v_p(t)\).
The conditions for \(\beta\) to be a virtual cubic
element are precisely the following.
\[
\begin{array}{c|l}
 p & \text{condition}\\ \hline
 2 & (k_2\le0,\ k_2\equiv0)\ \text{or}\
      (k_2\ge1,\ k_2\equiv1)\\
 3 & (k_3\le0,\ k_3\equiv1)\ \text{or}\
      (k_3\ge1,\ k_3\equiv0)\\
 7 & k_7=0,\ \text{or}\ (k_7<0,\ k_7\equiv-1),\
              \text{or}\ (k_7>0,\ k_7\equiv1)\\
 11 & k_{11}\equiv0,\ \text{and if }k_{11}=0
                    \text{ then }v_{11}(T(t))\ne1\\
 p\notin\{2,3,7,11\} & k_p\equiv0 .
\end{array}
\tag{1}
\]
Here congruences are modulo three and
\[
    T(t)=49t^2/9+125/3+196/t^2=\beta+\bar\beta .
\]
Whenever (1) holds, \(\beta\) is a cube at both places
above three if and only if
\[
                   3^{-k_3}t\equiv\pm4\pmod9.
 \tag{2}
\]

The polynomial of the two conjugates is
\(Z^2-T(t)Z+11^3\). At \(p\ne11\), an integral trace
means two unit conjugates. A negative trace valuation
gives two distinct integral Newton slopes and the
two roots belong to \(\mathbb Q_p\), with valuations
\(v_p(T)\) and \(-v_p(T)\). At \(2\) the three trace
terms have valuations \(2k_2,0,2-2k_2\). At \(3\) they
have valuations \(2k_3-2,-1,-2k_3\), and at \(7\)
they have \(2+2k_7,0,2-2k_7\). Whenever the minimum
is negative it is unique at each of these primes.
At other primes with \(k_p\ne0\) the minimum is
\(-2|k_p|\), and for \(k_p=0\) the trace is integral.
These observations give all rows except \(11\).

At \(11\), the Newton polygon has distinct slopes
when \(v_{11}(T)<3/2\), with root valuations
\(v_{11}(T),3-v_{11}(T)\). If \(v_{11}(T)\ge2\),
there is a single slope and a ramified quadratic
completion. Its normalized root valuation is three.
For \(k_{11}\ne0\) the trace has valuation
\(-2|k_{11}|\), giving the stated \(11\)-adic criterion.

At \(3\) the trace always has negative valuation.
Under (1) its valuation is a multiple of three
and is at most \(-6\). The root of smaller valuation
divided by \(T(t)\) is congruent to one modulo nine,
as follows by dividing the root equation by \(T(t)^2\).
The other root has cube norm relative to it.
Write \(w=3^{-k_3}t\).
After removing its \(3\)-power the trace is congruent
to \(4w^2\bmod9\) for \(k_3<0\), and to
\(7w^{-2}\bmod9\) for \(k_3>0\).
Among the square units modulo nine, the cube
condition is equivalent in both cases to
\(w^2\equiv7\bmod9\), giving (2).

## The allowable cube classes and the cover isomorphisms

**Corollary BS2.** If \(\beta\) is virtually cubic,
the rational cube class of \(t\) belongs to
\[
  \{2^i3^j7^k: i,j\in\{0,1\},\ k\in\{0,1,2\}\}.
 \tag{3}
\]
All twelve classes occur on infinite quadratic
subfamilies. Requiring local cubes above three
leaves exactly
\[
                         49,\ 147,\ 14,\ 42.
 \tag{4}
\]
On the fixed local neighborhoods of MO3, virtuality
of \(\beta\) is equivalent to \(t\in(3/7)\mathbb Q^{\times3}\).

Formula (3) follows from (1) because every sign is
a rational cube. For a representative \(c\) in (3)
use \(t=cu^3\). At \(2\) and \(3\) choose a valuation
of \(u\) giving the required sign for \(k_p\) in
(1). At \(7\), for the three residues of \(k_7\)
choose \(0,1,-1\), respectively. At \(11\) choose
a unit \(u\) near a solution to \(cu^3=1\), which
exists because cubing is bijective on \(11\)-adic
units. The trace at \(t=1\) is a unit there.
At other primes the condition is automatic.
These choices define a nonempty open slice in the
parameter \(u\) by weak approximation. The nonzero
squarefree polynomial \(F(cu^3)\) of degree 24 gives
infinitely many quadratic squareclasses by the same
Stewart–Top argument as in MO4.

In (3), the unit \(w\) of (2) is \(2^i7^k\) times
a \(3\)-adic unit cube. Hence the remaining condition
is \(2^i7^k\equiv\pm4\bmod9\). It leaves \(k=2\)
for \(i=0\) and \(k=1\) for \(i=1\), giving (4).
Their realization with (2) follows from the same
local choices. On the MO3 slice the valuations at
\(2,3,7\) are \(-3,-2,-4\), respectively, so (3)
forces the cube class \(3/7\). At eleven that slice
is already ramified and satisfies (1). Thus the
cubic parameter change is forced by virtuality on
those neighborhoods.

**Proposition BS3.** For \(c,d\in\mathbb Q^\times\),
the curves \(C_c:R^2=F(cu^3)\) and
\(C_d:R^2=F(du^3)\) are isomorphic over \(\mathbb Q\)
if and only if, in \(\mathbb Q^\times/\mathbb Q^{\times3}\),
\[
                [c]=[d]\quad\text{or}\quad[cd]=[6].
 \tag{5}
\]
Thus the classes (3) give six isomorphism classes
and (4) gives two, with pairs \((49,42)\) and
\((147,14)\). In each reciprocal pair the isomorphism
carries the four marked cube classes to the same span.

The octic satisfies \(t^8F(6/t)=6^4F(t)\). The geometric
stabilizer of the branch set of \(C_c\) on the
\(u\)-line therefore contains both \(u\mapsto\zeta_6u\)
and \(u\mapsto b/u\), with \(b^3=6/c^2\).
It is a finite group, since the 24 branch points
are distinct. By the characteristic-zero list of
finite subgroups of \(\operatorname{PGL}_2\)
([Beauville, Introduction](sources/beauville-finite-pgl2-subgroups.pdf)),
it is dihedral. It contains a dihedral group of
order 12, and the three exceptional groups have
no element of order six. In a dihedral group with
rotation order divisible by six, the order-six
subgroup of rotations is unique. Its two fixed
points are \(0,\infty\). Any isomorphism of the
hyperelliptic curves descends through their unique
hyperelliptic quotients and must preserve this pair.
In fact the geometric reduced group has order 12,
since a scaling preserving \(F(cu^3)\) up to a
constant has sixth power one by the nonzero
coefficients of \(1,u^6\).

A rational isomorphism therefore acts on \(u\)
by \(u\mapsto a u\) or \(u\mapsto a/u\).
For a scaling, compare constant and degree-six
coefficients of the two branch polynomials to
obtain \(da^3=\pm c\). For inversion, use the
reciprocal polynomial identity to reduce comparison
to \(6/(da^3)=\pm c\). Conversely each of these
changes lifts to an isomorphism by the same identities.
They act on the original parameter \(t\) by sign
and reciprocity; the two curves in each pair of (4)
therefore realize the same four-dimensional function
space by the exact span in
[CK3](MARKED_CUBIC_BASE_CHANGE_RATIONAL_THREE_TORSION.md#consequence-for-cubic-norm-functions).
Under \((t,R)\mapsto(6/t,36R/t^4)\), the MO1
coordinates send \((x,y,w)\) to \((x,-y,w)\).
The \(Y_1,Y_3\) signs reverse and \(Y_2\) is unchanged.
Thus \(\alpha_1,\alpha_3\) invert modulo cubes and
\(\alpha_2\) stays fixed. The sign change of \(t\)
inverts \(\alpha_1,\alpha_2\) and preserves \(\alpha_3\).
The formula for \(P+R\) in MO3 shows that both
changes preserve \(\beta\).

## A split slice with four marked classes on both covers

**Theorem BS4.** For each class \(c\) in (4), there
are \(\gg_c X^{1/12}/(\log X)^2\) distinct real
quadratic discriminants \(D\le X\) from \(t=cu^3\)
for which all four marked classes are independent
virtual cubic classes and cubes in the completions
above three. The real and imaginary mirror class
\(3\)-ranks are at least three and four, respectively,
and \(2,3,5,7,11\) all split in the real field.
The three fixed elliptic twists of MO1 have positive
rank, with conductors \(1540D^2,13860D^2,6930D^2\).
These local conditions give fields disjoint from the
original MO3 slice, in which \(5,11\) ramify.

Here are local centers for the four scales in (4).
At each listed prime use \(u=p^e\), with exponents
\[
\begin{array}{c|rrrrr}
 c & p=2 & 3 & 5 & 7 & 11\\ \hline
 49  &-5&5&5&-5&5\\
 147 &-5&-5&5&-5&5\\
 14  &5&5&5&5&5\\
 42  &5&-5&5&5&5.
\end{array}
\]
At these centers the squareclasses split and the
valuations of both evaluations of each function
are multiples of three. At \(3\) all four normalized
unit pairs have residues \(1,-1\bmod9\), in some
order. These are open conditions and can be checked
by the finite-precision valuation and residue
calculations in the
[checker](verify_marked_octic_beta_scale_classification.py).
Away from the listed primes the three original
marks are virtual by the good-prime argument
of MO2, and \(\beta\) is virtual by (1).

For independence it suffices to give rows on the
representatives \(49,147\), using reciprocity for
their partners. For \(c=49\) the following table
gives residues of \(u,R\) and a cube root of unity
\(\zeta\), and exponents of the four character values.
\[
\begin{array}{c|rrr|c}
 q&u&R&\zeta&\text{exponents}\\ \hline
 13&4&2&3&(1,2,1,0)\\
 31&1&8&25&(1,1,2,0)\\
 37&2&13&26&(1,1,1,0)\\
 43&5&11&36&(0,2,0,1).
\end{array}
\]
For \(c=147\) use
\[
\begin{array}{c|rrr|c}
 q&u&R&\zeta&\text{exponents}\\ \hline
 13&1&1&3&(0,2,2,1)\\
 19&4&2&7&(2,2,1,2)\\
 31&4&12&25&(2,1,2,2)\\
 37&1&10&26&(2,1,2,1).
\end{array}
\]
Both matrices have rank four. The conjugate branch
gives the negative row because each norm is a cube.
CRT gives simultaneous open slices with all the
required local tests. Equivalently one may transport
the whole chosen slice for each partner by the
reciprocal isomorphism. The proof of MO3 applies
to the four classes, including the real unit loss,
the mirror descent and the positive ranks from the
three noncubic elliptic markings.

Every bad prime of all three original elliptic curves
now splits in \(K\), so the discriminant twist leaves
their local conductors there unchanged. At every
new ramified prime of \(K\) the twist of good
reduction has conductor exponent two. The conductors
follow using the original conductors of MO2.
The degree-24 squarefree polynomial and the fixed
local slices give the claimed count by the proof of
MO4, with \(c\) fixed.

CK2 gives exact rational \(3\)-torsion rank four on
the genus-eleven curve for every \(c\in\mathbb Q^\times\),
using Chebotarev transport of its two finite-field
certificates. This includes all twelve virtual scales
in BS2. Changing the constant scale cannot supply a
fifth generic direction. A common quadratic field obtained
from the two inequivalent covers could still carry
additional classes from the two evaluations.
No infinitude of such common fields is asserted.
