# Three marked cubic norm directions on a genus-three hyperelliptic curve

The curve below has three rational \(3\)-torsion directions
on its Jacobian. The three associated cubic norm functions
remain independent on a congruence family of quadratic
fields. Two separate valuations then remove two directions
from the virtual-unit subspace. The surviving direction can
come from a real quadratic unit.

## A split genus-three hyperelliptic Jacobian

**Theorem HG1.** The smooth projective model of
\[
 C:\quad R^2=t^8-6t^7-t^6+18t^5+93t^4-36t^3-4t^2+48t+16
          =F(t)
 \tag{1}
\]
has genus three and
\(\operatorname{Jac}(C)(\mathbb Q)[3]\simeq(\mathbb Z/3)^3\).
Its Jacobian is \(\mathbb Q\)-isogenous to the product of
the three curves
\[
 E_i:\quad Y^2=X^3+\bigl(h_iX+h_i+1\bigr)^2,\qquad
 (h_1,h_2,h_3)=(3/2,-3/4,-7/6).
 \tag{2}
\]
For the quadratic twist by \(F(t)\) each of these fixed
curves has generic rank one over \(\mathbb Q(t)\).
Their marked cubic norm functions on \(C\) are independent
over \(\overline{\mathbb Q}(C)\).

Here are explicit maps and the calculation behind the
construction. Parametrize \(6x^2+3y^2=2\) by
\[
 d_0=3(t^2+2),\qquad
 x=\frac{t^2-4t-2}{d_0},\qquad
 y=\frac{-2t^2-4t+4}{d_0},\qquad w=\frac{12R}{d_0^2}.
\]
The conic parameter comes from the line
\(y=2/3+t(x-1/3)\), and direct substitution gives
\[
 w^2=-\frac92 x^4-\frac{117}{16}y^4+\frac{59}{18}.
 \tag{3}
\]
The polynomial \(F\) is squarefree, for example modulo
\(7\). In (3) view \(x,y\) as the affine coordinates on the
conic and \(w\) as the function cutting out its double cover.
The involutions \(x\mapsto-x\) and \(y\mapsto-y\) fix \(w\).
Their three nonidentity combinations give genus-one
quotients, with respective equations
\[
 v_i^2=a_i z_i^4+b_i z_i^2+c_i,\qquad
 \begin{array}{c|l|l|r|r|r}
 i&z_i&v_i&a_i&b_i&c_i\\ \hline
 1&y&w&-135/16&3/2&25/9\\
 2&x&w/2&-135/16&39/8&1/144\\
 3&x/y&-w/(3y^2)&25/9&59/18&1/144 .
 \end{array}
 \tag{4}
\]
These identities follow by eliminating one squared conic
coordinate. Each has a rational point at \(z_i=0\),
since \(c_i\) is a nonzero square. The full group quotient
of \(C\) has genus zero: on eliminating \(x^2\), it is
a conic in \(y^2,w\). Thus the three genus-one quotients give
the three distinct nontrivial eigenspaces of the sign
group on the Jacobian.

For any even quartic in (4) the functions
\[
 U=\frac{v^2}{4z^2},\qquad
 V=\frac{v}{8z}\left(a z^2-\frac{c}{z^2}\right)
 \tag{5}
\]
map it to \(V^2=U(U^2-\frac b2 U+(b^2-4ac)/16)\).
Indeed \(V^2/U=\frac1{16}(a z^2-c/z^2)^2\), which equals
\(U^2-bU/2+(b^2-4ac)/16\).
This map factors through the two degree-two maps given
by \((a z^2,a z v)\) and the quotient by the marked
\(2\)-torsion on \(Y^2=X(X^2+bX+ac)\).
Equivalently, its abscissa is an even rational function
of degree four in \(z\), which gives degree four after
accounting for the two hyperelliptic covers.
For each row in (4), \(X_i=U_i-1,\ Y_i=V_i\) give a map
from \(C\) to the corresponding \(E_i\) of degree eight.
The factors are pairwise nonisogenous over \(\mathbb Q\):
the traces at \(11\) are \(6,-6,0\).
The sign eigenspaces show that the product map on
Jacobians is an isogeny. The three rational Hom spaces
to the targets have rank one, since an elliptic curve
over \(\mathbb Q\) has rational endomorphism ring
\(\mathbb Z\). The points in (5) are anti-invariant under
\(R\mapsto-R\). This gives the exact generic twist ranks.

On \(E_i\), the function \(Y+h_iX+h_i+1\) has divisor
\(3((0,-h_i-1)-O)\). Pullback along a degree-eight map
keeps the divided divisor nonzero by pushforward.
The three pulled-back classes are independent: each is
fixed by the involution defining its quotient, and the
other two involutions act by \(-1\) on it. This follows
on the quotient in (4), where the other involutions act
by \(z_i\mapsto-z_i\) and send (5) to its negative.
The Weil pairing bounds the rational \(3\)-rank by
three, since \(\mathbb Q\) does not contain a primitive
cube root of unity.

This construction is a specialization of the genus-three
hyperelliptic gluing in
[Howe–Leprévost–Poonen, Proposition 14](sources/howe-leprevost-poonen-split-jacobians.pdf).
The direct quotient calculation above establishes the
specified claims on (1). The question of genus three
with rational Jacobian \(3\)-rank three is recorded after
Theorem 9 in
[Levin–Yan–Wiljanen](sources/levin-yan-wiljanen-large-three-rank.pdf).
The explicit curve (1) satisfies that property.

## Three specialized Kummer directions and a common conductor

**Theorem HG2.** Impose the rational parameter congruences
\[
 t\equiv1\bmod 9,\qquad t\equiv1\bmod13,\qquad
 t\equiv4\bmod19,\qquad t\equiv1\bmod31 .
 \tag{6}
\]
The field \(K=\mathbb Q(\sqrt{F(t)})\) is real quadratic,
with discriminant \(D\).
The three elements
\(\alpha_i=Y_i+h_i X_i+h_i+1\) are independent in
\(K^\times/K^{\times3}\).
The three \(D\)-twists of (2) are pairwise nonisogenous
elliptic curves over \(\mathbb Q\), each of positive rank
and conductor \(80D^2\).

First \(F(t)>0\) for all real \(t\). In (3), substitute
\(y^2=2/3-2x^2\) and take \(0\le x^2\le 1/3\).
The right hand side is concave as a function of \(x^2\),
and its endpoint values are \(1/36\) and \(25/9\).
Thus the specialized fields are real. Also
\(F(1)=129\equiv3\bmod9\), so the field is ramified
at \(3\), with local squareclass \(3\).

At \(13,19,31\), (3) has the following split evaluations.
The second column specifies \(w\) rather than \(R\).
Write each cube residue \(\alpha_i^{(q-1)/3}\) as a power
of the indicated primitive cube root in \(\mathbb F_q\).
\[
 \begin{array}{c|r|r|r|c|c}
 q&t&w&\zeta&
 (\alpha_1,\alpha_2,\alpha_3)\bmod q&
 \text{exponents}\\ \hline
 13&1&6&3&(9,10,4)&(2,1,2)\\
 19&4&7&7&(15,3,17)&(2,1,1)\\
 31&1&6&5&(20,10,8)&(1,1,0).
 \end{array}
\]
The rows have nonzero determinant modulo three.
The quadratic norms of the \(\alpha_i\) are
\(-X_i^3\), so the other split branch negates each
corresponding row. Independence therefore holds on (6).
The cubic factor after removing \(X\) from the \(3\)-division
polynomial of each \(E_i\) is irreducible over \(\mathbb Q\).
Hence there is no nonzero anti-invariant \(3\)-torsion
point in \(E_i(K)\).
The [Cohen–Pazuki evaluation](sources/cohen-pazuki-2009-three-isogeny.pdf)
on the quadratic twist differs from
\(\alpha_i\) by a cube or inversion, as in BT2. The
noncubic classes give infinite-order twist points.

The base conductors of the three curves are \(90,90,30\).
As shown in HG3 below, \(2\) and \(5\) both split in \(K\),
so \(D\) is the positive squarefree radicand and is
coprime to \(10\). At \(3\), each twist is locally
isomorphic to the twist by \(3\), and those three fixed
twists all have conductor \(720\). At other primes
dividing \(D\) the base curves have good reduction and
the twist conductor exponent is two. Thus all three
conductors are \(720(D/3)^2=80D^2\). Twisting preserves
rational isogeny and nonisogeny.

## Two independent valuation obstructions

**Proposition HG3.** For any nonsplit quadratic
specialization of (1) away from poles and zeros of the
displayed norm functions, \(2\) and \(5\) split in \(K\).
At \(2\) the classes of \(\alpha_1,\alpha_3\) have zero
valuation modulo three, and \(\alpha_2\) has nonzero
valuation. At \(5\), the classes of \(\alpha_2,\alpha_3\)
have zero valuation modulo three, and \(\alpha_1\) has
nonzero valuation. Under (6), the virtual-unit subspace
of the marked span is exactly the line of \(\alpha_3\).
At \(t=1\), the class of \(\alpha_3\) is nontrivial and
comes from a unit of \(\mathbb Q(\sqrt{129})\).

At \(2\) the valuations of the two conic terms \(6x^2\)
and \(3y^2\) have opposite parity. Their minimum must be
one, so \(x\) is odd and \(y\) is even. Reduction
modulo eight shows \(v_2(y)=1\). Put \(y=2u\) with \(u\)
odd, so \(3x^2+6u^2=1\). Reducing modulo sixteen gives
\(x^2\equiv9\bmod16\). From the second row of (4),
\[
 X_2=\frac{-1215 x^4+126x^2+1}{576x^2}.
\]
Write \(x^2=1+8j\) with \(j\) odd. Then
\[
 X_2=-\frac{17+288j+1215j^2}{9x^2}
\]
is even. Its associated linear term is
\(L_2=(1-3X_2)/4\) of valuation \(-2\).
The ratio \(Y_2^2/L_2^2=1+X_2^3/L_2^2\) is
\(1\bmod 2^7\), so \(2\) splits in \(K\).
On one branch \(Y_2/L_2\equiv1\bmod4\), which gives
\(v_2(\alpha_2)=v_2(2L_2)=-1\). The other branch has
opposite valuation modulo three by the cube norm.
As \(U_2=X_2+1=w^2/(16x^2)\) is odd, \(v_2(w^2)=4\).
In each of the first and third rows of (4) this makes
\(U_i\) odd, hence \(X_i\) even and the respective linear
terms \(L_i\) have valuation \(-1\).
On one branch \(v_2(\alpha_i)=v_2(2L_i)=0\), and the norm
again gives a multiple of three on the other.

At \(5\), the conic forces \(v_5(y)=0\). A positive
valuation would require \(x^2\equiv2\bmod5\).
A negative valuation would require \(x\) to have the same
valuation and \((y/x)^2\equiv-2\bmod5\).
In the first row of (4),
\[
 X_1=\frac{-135 y^2/16-5/2+25/(9y^2)}{4}.
\]
If \(y^2\equiv1\bmod5\) then \(X_1/5\equiv0\bmod5\),
and if \(y^2\equiv4\bmod5\) then \(X_1/5\equiv1\bmod5\).
In either case \(L_1=(5+3X_1)/2\) has valuation one.
Hence \(Y_1^2/L_1^2\equiv1\bmod5\), which gives
splitting and \(v_5(\alpha_1)=1\) on one branch.
For \(\alpha_2,\alpha_3\) the constant terms of the
linear forms are units at \(5\). The cube norm then
gives the virtual-unit condition by the same local
argument as in BT2.

For \(t\) satisfying (6), \(3\) is ramified and there is
only one prime above it. Any element with rational cube
norm has valuation divisible by three at a rational
prime with just one prime above it. At every prime
other than \(2,3,5\), the coefficients of all three
linear terms are integral and their constants are units.
The pole or unit argument of BT2 applies. Thus the
two independent valuation obstructions cut out exactly
the specified line.

At \(t=1\), write \(r=\sqrt{129}\) with the sign giving
\(w=12r/81\). The field has trivial class group, and
\[
 \alpha_3=-\frac{251}{50}-\frac{217r}{250}
  =(16855-1484r)\left(-\frac{397}{10}-\frac{7r}{2}\right)^3 .
\]
Here \(16855-1484r\) is a fundamental unit, up to
inversion and sign. This is an explicit specialization
where a noncubic virtual unit yields the trivial ideal
class. The class group and fundamental unit in this
example are checked with proof enabled.

## An unbounded family with the same three ranks and conductor law

**Corollary HG4.** At least
\(\gg X^{1/4}/(\log X)^2\) real quadratic fields with
discriminant at most \(X\) satisfy HG2 and HG3. There are
\(\gg Z^{1/8}/(\log Z)^2\) such common conductor values
at most \(Z\).

The integral homogenization of \(F\) is a squarefree
binary octic. Apply
[Stewart–Top, Theorem 2](sources/stewart-top-ranks-of-twists.pdf)
with \(k=2\) and the fixed congruences (6) on the
numerator and denominator, the latter congruent to one.
The resulting distinct squarefree parts are all positive
and give the required discriminants by HG2–HG3.
The conductor formula is injective in \(D>0\), which
gives the second lower bound. The
[exact checker](verify_genus_three_marked_elliptic_quotients.py)
replays the three maps and the witnesses for the local
and conductor calculations.
