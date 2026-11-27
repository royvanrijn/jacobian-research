# Four virtual cubic classes after a cubic change of parameter

On the first genus-three gluing in
[HG1–HG3](GENUS_THREE_THREE_MARKED_ELLIPTIC_QUOTIENTS.md),
two different primes obstruct two of the three marked
norm directions. A second gluing admits a local parameter
slice on which all three are virtual units and cubes at
\(3\). Its octic also has a fourth norm identity. A cubic
change of parameter makes the valuations of this fourth
direction divisible by three at moving primes. The four
directions then give class \(3\)-rank at least three in
the real field and at least four in its imaginary mirror.

## A second marked octic

**Theorem MO1.** Put
\[
 F(t)=49t^8+750t^6-2403t^4+27000t^2+63504.
 \tag{1}
\]
The curve \(C:R^2=F(t)\) has genus three and its Jacobian
is isogenous over \(\mathbb Q\) to the product of
\[
 E_i:\quad Y^2=X^3+(h_i X+h_i+1)^2,\qquad
 (h_1,h_2,h_3)=(-16/11,24/11,3/22).
 \tag{2}
\]
The pulled-back cubic norm functions
\(\alpha_i=Y_i+h_iX_i+1+h_i\) are geometrically
independent. Each of the three fixed elliptic twists
by \(F(t)\) has generic rank one.

Use
\[
 x=(t+6/t)/2,\quad y=(6/t-t)/6,\quad
 w=\frac{108R}{121(6t)^2}.
\]
Then \(x^2-9y^2=6\), and
\[
 w^2=\frac{5778x^4+103518y^4-116883}{14641}.
\]
The three genus-one quotients under the sign involutions
in \(x,y\) have the same shape as in HG1, with data
\(v_i^2=a_i z_i^4+b_i z_i^2+c_i\) given by
\[
 \begin{array}{c|r|r|r|r|r}
 i & z_i & v_i & a_i & b_i & c_i\\ \hline
 1 & y & -11w/54 & 196/121 & 214/121 & 125/484\\
 2 & x & 11w/6 & 196/121 & -426/121 & -7875/484\\
 3 & x/y & 11w/(9y^2) &125/484&1443/242&-7875/484.
 \end{array}
\]
The quartic-to-elliptic maps are exactly (5) of HG1,
with \(X_i=U_i-1,\ Y_i=V_i\). The identities can also be
checked directly in (1). All three maps have degree
eight and lie in the three nontrivial sign eigenspaces.
The joint morphism on Jacobians is an isogeny.
The curves have traces \((-4,2,-4)\) at \(13\) and
\((6,6,0)\) at \(23\), giving pairwise rational
nonisogeny. This proves the generic ranks and independence
of the divided divisors as in HG1.
The octic is squarefree modulo \(13\).
It is positive on the real line, since \(x^2\ge6\)
and \(5778\cdot36>116883\).

## Three virtual classes on the original parameter line

**Theorem MO2.** There are
\(\gg X^{1/4}/(\log X)^2\) distinct real quadratic
fields \(K\) of discriminant \(D\le X\) from (1) such that
\[
 \operatorname{rk}_3\operatorname{Cl}(K)\ge2,\qquad
 \operatorname{rk}_3\operatorname{Cl}(\mathbb Q(\sqrt{-3D}))\ge3.
\]
The three \(D\)-twists of (2) have positive rank and
conductors \(44D^2,396D^2,198D^2\), respectively. The
bounds hold on one fixed open parameter slice in the
local conditions described below.

Take \(t\) in sufficiently small neighborhoods of the
following centers. The first row of data gives the
local squareclass of \(F(t)\), expressed by the valuation
parity and normalized unit where needed.
\[
 \begin{array}{c|r|l}
 p&t_p&\text{squareclass condition}\\ \hline
 2&-1/2&\text{split}\\
 3&-2/3&\text{split}\\
 5&-1/2&\text{odd valuation}\\
 7&-1&\text{odd valuation}\\
 11&-1&\text{split}.
 \end{array}
\]
Write \(L_i=h_iX_i+h_i+1\) and \(S_i=Y_i^2\).
In a split completion, the valuations of both conjugates
of \(\alpha_i\) are \(0\bmod3\) if the smaller of the
two is \(0\bmod3\), by the cube norm. When the valuations
of \(Y_i,L_i\) differ, this smaller value is the minimum
of those valuations. If the valuations agree, it is
their common value for odd primes, and this common
value plus one at \(2\). These claims follow by
subtracting the two conjugate expressions, or at
\(2\) by reducing the two unit factors modulo four.

At the \(2\)-adic center the pairs
\((v_2(Y_i),v_2(L_i))\) are
\((-6,1),(-6,-1),(1,0)\).
At the \(11\)-adic center they are \((-3,-3)\) for
each row. These give virtuality at both split primes.
At \(5,7\) the field has one prime above the rational
prime, where the cube norm forces the valuation to be
divisible by three. At \(3\) multiply the three
evaluations by \(3^6,3^3,1\), respectively.
The resulting pairs \((S_i,L_i)\) after scaling
by the square and first power of these factors are
\((1,0),(1,0),(7,4)\bmod9\).
The normalized sums are cubes on at least one branch
of \(\mathbb Q_3\), and on the other by the cube norm.
Indeed the units congruent to \(\pm1\bmod9\) are cubes,
and in the last pair the choice \(Y_i\equiv4\bmod9\)
gives sum \(8\). Off these five fixed primes the constants
of the linear forms are units and the pole or unit
argument of BT2 applies.

For independence impose \(t\equiv10\bmod13,\ t\equiv
3\bmod19,\ t\equiv2\bmod31\). At the respective split
branches \(w=2,1,14\), the residue triples are
\((10,4,7),(13,15,2),(10,12,5)\).
With cube roots of unity \(3,7,5\) in the residue fields,
their cubic character rows are
\((1,2,2),(2,2,1),(1,2,1)\), of rank three.
The local constraints are open and simultaneously
attainable by weak approximation.

In a real quadratic field the unit contribution to
the space of virtual units modulo cubes has dimension
one. The three independent virtual classes therefore
give real ideal class rank at least two. Adjoin
\(\zeta_3\) to \(K\), and use the descent of
[MT2](QUADRATIC_MIRROR_RANK_TWO_FROM_MARKED_MESTRE_TWISTS.md#the-common-local-conditions-and-the-mirror-classes)
on all three classes. The argument works with these
signs as well. The cube norms make the simultaneous
conjugation fixing the imaginary mirror act trivially
on the Kummer group. Virtual valuations give
unramifiedness away from \(3\), and the local cubes at
\(3\) give complete splitting there. Thus the imaginary
mirror has an elementary abelian unramified cubic
extension of degree \(27\), split over \(3\).

The base conductors are \(1540,13860,6930\).
All three cubic factors after removing \(X\) from
the \(3\)-division polynomials are irreducible.
As in HG2, the Kummer evaluations of the three twist
points give positive rank on all three curves. The
local twists at \(2,3,11\) are trivial and at \(5,7\)
each ramified twist gives conductor exponent two.
The same exponent holds at the other odd primes
dividing the radicand, where the base curves have good
reduction. Since \(2\) splits, the squarefree radicand
equals \(D\) and the stated law follows.
Finally apply
[Stewart–Top, Theorem 2](sources/stewart-top-ranks-of-twists.pdf)
to the squarefree binary octic in fixed numerator and
denominator congruences inside the stated neighborhoods.
The fields give distinct imaginary mirrors of
discriminant \(-3D\).

## The fourth norm identity and a cubic change of parameter

**Theorem MO3.** On a fixed local parameter slice of
\[
 t=3u^3/7,\qquad u\in\mathbb Q^\times,
 \tag{3}
\]
the real field \(K=\mathbb Q(\sqrt{F(t)})\) has four
independent virtual classes modulo cubes, all local
cubes at \(3\). Its real class \(3\)-rank is at least
three and the class \(3\)-rank of its imaginary mirror
is at least four. The latter has at least forty complex
cubic fields of fundamental discriminant \(-3D\) and
an unramified \((\mathbb Z/3)^4\) extension split at \(3\).
The prime classes above \(19,31,37,43\) are independent
modulo cubes in the imaginary mirror class group.
The three \(D\)-twists in (2) all have positive rank, with
conductors \(28D^2,252D^2,126D^2\), respectively.

Use the additional identity
\[
 P=7t^4+\frac{375}{7}t^2+252,\qquad
 F(t)=P^2-11\left(\frac{198}{7}t^2\right)^2,\qquad
 \beta=\frac{7(P+R)}{18t^2}.
 \tag{4}
\]
Then \(\operatorname N(\beta)=11^3\).
At a prime \(p\nmid2\cdot3\cdot5\cdot7\cdot11\) where
\(t\) is a unit, both conjugates of \(\beta\) are units.
If \(t\) has nonzero valuation there, the completion
splits: \(F/P^2\) is a square unit congruent to one
modulo \(p\). One branch has \(P+R\) of valuation
\(v_p(P)\), which equals \(0\) for \(v_p(t)>0\)
and \(4v_p(t)\) for \(v_p(t)<0\).
Thus the two valuations of \(\beta\) are
\(\pm2v_p(t)\). Condition (3) makes both multiples of
three, even though the primes moving with the parameter
are unrestricted. As before, the \(\alpha_i\) are
virtual at all these primes.

For the five fixed primes take the following centers
in the \(t\)-coordinate:
\[
 \begin{array}{c|r|l}
 p&t_p&\text{squareclass condition}\\ \hline
 2&3/56&\text{split}\\
 3&1/63&\text{split}\\
 5&-1/2&\text{odd valuation}\\
 7&3/2401&\text{split}\\
 11&19&\text{odd valuation}.
 \end{array}
\]
They lie in the image of (3) over the respective local
fields. For \(2,3,7\) the \(u\)-centers are \(1/2,1/3,1/7\).
For \(5,11\) the desired \(u^3\) is a unit, and cubing
is bijective on those local unit groups. At \(5,11\)
the cube norms take care of all four virtual-unit
conditions. At \(2\), the three \((v_2(Y_i),v_2(L_i))\)
pairs are \((-12,-4),(-12,-5),(1,0)\);
the corresponding pair for \(\beta\), separating
its multiple of \(R\) and its rational term in (4),
is \((-7,-7)\).
Thus the smaller branch valuations are
\(-12,-12,0,-6\), all multiples of three.
At \(7\) those four pairs are
\((-9,-6),(-9,-6),(1,0),(-6,-6)\), with minima
\(-9,-9,0,-6\).

At \(3\), multiply \(\alpha_1,\alpha_2,\alpha_3,\beta\)
by \(3^9,3^6,1,3^6\), respectively, all rational cubes.
The normalized pairs consisting of the squared irrational
term and the rational term are
\[
 (1,0),\quad(1,0),\quad(7,4),\quad(7,5)\pmod9 .
\]
The first two give unit cubes on both branches.
For the last two select the root \(4\) of \(7\bmod9\):
the sums are \(8,0\bmod9\), respectively. For the fourth
use the other root, \(-4\), to obtain \(1\bmod9\).
In each case one cube evaluation and the cube norm
give both local cubes. These calculations hold in
sufficiently small neighborhoods, since the valuations
and the finite unit residues are locally constant.

The following four residue rows enforce independence.
Each value \(w\) selects a split branch of the octic
in the displayed reduction. The last column records
the exponents after taking \((q-1)/3\) powers with
respect to \(\zeta\).
\[
 \begin{array}{c|r|r|r|r|c|c}
 q&u&t&w&\zeta&(\alpha_1,\alpha_2,\alpha_3,\beta)&
 \text{exponents}\\ \hline
 19&4&3&1&7 &(13,15,2,4)&(2,2,1,2)\\
 31&3&16&14&5 &(14,22,11,21)&(2,1,1,1)\\
 37&3&1&3&10 &(9,4,7,35)&(2,1,1,2)\\
 43&21&13&13&6 &(3,26,39,35)&(2,1,0,0).
 \end{array}
\]
The matrix is invertible modulo three.
All the congruences are conditions on \(u\) and can
be imposed together. The unit kernel consumes at
most one real direction, giving the real class bound.
The mirror descent now gives degree \(81\) and its
\((81-1)/2=40\) distinct fundamental cubic fields.
At the four indicated primes the splitting field
\(K(\zeta_3)\) is completely split and the four rows
are the Frobenius evaluations in the unramified
mirror extension. Thus the corresponding four prime
classes span the degree-\(81\) quotient.

The rank proof of MO2 applies to the three original
directions. The fixed bad primes now split at \(2,3,7\)
and ramify at \(5,11\). At the latter two primes the
ramified twists of the multiplicative base reduction
have exponent two. This gives the stated three
conductors.

**Corollary MO4.** The fields in MO3 have at least
\(\gg X^{1/12}/(\log X)^2\) distinct real discriminants
at most \(X\), and distinct imaginary mirrors of
absolute discriminant at most \(3X\).
The corresponding conductor triples have at least
\(\gg Z^{1/24}/(\log Z)^2\) values with largest conductor
at most \(Z\).

Indeed the degree-\(24\) polynomial \(F(3u^3/7)\) is
squarefree since \(F\) is squarefree and \(F(0)\ne0\).
Clear the fixed denominator by a square and apply
Stewart–Top, Theorem 2, on numerator and denominator
congruences for \(u\). The
[exact checker](verify_marked_octic_mirror_class_rank_four.py)
replays the maps, the new norm identity, both sets of
local residue data, and the binary form hypotheses.

## Rigidity of the two boundary fibers

**Theorem MO5.** On \(C\) let \(P_0^\pm\) and \(P_\infty^\pm\)
be the two points over \(0,\infty\), labeled by the signs
of \(R\) at zero and of \(R/t^4\) at infinity. Set
\[
 \delta_0=P_0^- -P_0^+,\qquad
 \delta_\infty=P_\infty^- -P_\infty^+,\qquad
 B_0=\delta_0+\delta_\infty .
\]
Then \([B_0]\) has exact order two, and the subgroup of
principal divisors in \(\mathbb Z\delta_0+\mathbb Z\delta_\infty\)
is generated by \(2B_0=\operatorname{div}(\beta)\).
For an odd prime \(\ell\) and any \(c\in\mathbb Q^\times\),
let \(C_{\ell,c}\) be the smooth model of
\(R^2=F(c u^\ell)\). The corresponding four boundary
points give a class \([B_{\ell,c}]\) of order exactly
\(2\ell\). The corresponding subgroup of principal
divisors on \(C_{\ell,c}\) is generated by
\(2\ell B_{\ell,c}\).
In particular the divided divisor of \(\beta\) after
pullback gives a rational class of exact order \(\ell\).

The polynomial in (4) shows
\(R(0)=\pm252\) and \(R/t^4=\pm7\) at the two infinities.
The rational function \(\beta\) has a pole of order two
on the two plus branches and no other poles.
Its quadratic norm is constant, so the matching
minus branches are its two zeros of order two.
Thus its divided divisor \(B_0\) has order at most two.
If \(B_0\) were principal, the poles
\(P_0^++P_\infty^+\) would give a degree-two map to
\(\mathbb P^1\). The unique hyperelliptic degree-two
pencil has invariant fibers under the hyperelliptic
involution. These two points do not form such a fiber.

The image of \(P_0^+\) on \(E_3\) in the map of MO1 is
\((75/121,1750/1331)\), and \(P_0^-\) maps to its
negative. The point has infinite order. Its order if
torsion would divide six, from the reductions of \(E_3\)
at \(13,23\), whose group orders are \(18,24\); in fact its
sixfold multiple is nonzero. Thus \([\delta_0]\) is
nontorsion. As \([B_0]\) has order two, any principal
combination \(m\delta_0+n\delta_\infty\) must have
\(m=n\) and \(n\) even. This proves the first subgroup
statement even over the algebraic closure.

The map \(C_{\ell,c}\to C\) is totally ramified over the
four boundary points, so the pullback of \(\beta\) has
divisor \(2\ell B_{\ell,c}\). Pushforward sends
\(B_{\ell,c}\) to \(B_0\) and similarly on both
differences, so the order of \([B_{\ell,c}]\) is even
and any principal combination upstairs has equal
coefficients. It remains to rule out order two.
Over the algebraic closure the function-field
extension is obtained by adjoining \(t^{1/\ell}\).
If \(\beta\) became an \(\ell\)-th power there, its original
class modulo \(\ell\)-th powers would belong to the line
of \(t\), by cyclic Kummer theory. The valuations of \(t\)
at \(P_0^+,P_0^-\) are both one, while those of \(\beta\)
are \(-2,2\). These are unequal modulo any odd prime
\(\ell\). Thus \(2B_{\ell,c}\) cannot be principal,
giving order \(2\ell\) and the asserted subgroup.
This proves that the two indicated boundary fibers
supply a single additional norm direction by this
mechanism. The arithmetic specialization in MO3
uses the \(\ell=3\) instance and also needs its finite
prime conditions.
