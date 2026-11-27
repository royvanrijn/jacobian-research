# Coherent prime triangles in a squarefree field window

For unit normalized coefficients, the Heron sextic of a prime
triangle is irreducible on every viable fixed linear slice of the
prime directions. A change of coefficients to \((5,1,4)\) gives a
linear and quintic factor on a prime progression slice. The
quintic large-square sieve then gives genuine prime triangles
with the field discriminant comparable to the unreduced norm
discriminant.

## The obstruction on the unit-coefficient slice

**Proposition SW1.** Let
\[
 H_0(P_1,P_2,P_3)=
 2\!\sum_{i<j}P_i^3P_j^3-\sum_i P_i^6 .
\]
On the rational line \(aP_1+bP_2+cP_3=0\) with primitive integer
coefficients, assume \(abc\ne0\) and \(a+b+c\) even, and that no
pair of coefficients sums to zero. The restriction of \(H_0\) is
an irreducible binary sextic over \(\mathbb Q\).
If a pair sums to zero, the line cannot support a primitive
unit-normalized prime triangle with all primes larger than
\(\max(|a|,|b|,|c|)\). A coefficient zero likewise cannot support
such a triangle at those scales. A nonhomogeneous fixed affine
linear relation also leaves an irreducible polynomial of degree
six.

The Gauss map of the Fermat cubic
\(C:X_1^3+X_2^3+X_3^3=0\) is
\(\nu:X\mapsto[X_1^2:X_2^2:X_3^2]\). Its image is the dual
sextic \(H_0=0\), as also follows from
[van Luijk, Section 2 and Example 2.13](https://arxiv.org/pdf/0708.2752).
The map is injective on geometric points: after rescaling two
points with the same image, they differ by coordinate signs.
Up to global sign at most one sign differs, and changing just
the sign of a nonzero cube cannot preserve the Fermat equation.
Thus a point and its image have the same residue-field degree.
The line intersects the dual at the images of the six-point
intersection divisor on \(C\) with the conic
\(aX_1^2+bX_2^2+cX_3^2=0\).

Exactly two of \(a,b,c\) are odd. Relabel so these are \(a,b\).
Every geometric intersection point reduces at \(2\) to
\(O=[1:-1:0]\), since the conic modulo \(2\) is the double line
\(X_1+X_2=0\). Use the rational flex \(O\) as origin on \(C\),
which has good reduction at \(2\). The rational group consists
of its three rational flexes by Fermat for cubes, and those
three points have distinct reductions. The sum of any rational
effective divisor supported in the intersection is therefore
\(O\). Indeed the sum reduces to \(O\), and is rational by
Galois invariance.

There is no rational point in the intersection, since any
rational Fermat point would require a pair among \(a,b,c\)
to sum to zero. A degree-two closed point would consist of
conjugates \(P,-P\). Negation on \(C\) with origin \(O\) is
interchange of \(X_1,X_2\), since that nontrivial involution
fixes \(O\). Their squared first two coordinates are different.
The equality \(X_1=-X_2\) gives \(O\), and \(X_1=X_2\) gives
a point of degree three from \(X_3^3=-2 X_1^3\).
So the conic equations for the two points force \(a=b\).
Normalize \(X_1+X_2=1\). The value \(z=X_3\) is now a
nonzero rational fixed by conjugation. The equations of the
Fermat cubic and conic give
\[
                     \frac{c}{a}=\frac{2z^3-1}{3z^2}.
\]
The 2-adic valuation on the right is nonpositive. It is
\(-2v_2(z)\) for \(v_2(z)\ge0\), and \(1+v_2(z)\) for
\(v_2(z)<0\). This contradicts \(a\) odd and \(c\) even.

A closed point of degree three gives a divisor \(D\) of class
\(3O\), since its sum is \(O\). Its complementary divisor in the
conic intersection likewise has class \(3O\), because a conic
section has class \(6O\). Both degree-three divisors are
therefore line sections in the given embedding, as the flex
tangent has divisor \(3O\). Their product cuts the same divisor
as the conic. The conic equals that product up to scalar, since
two quadratic forms with identical zero divisors on \(C\)
differ after scaling by a form vanishing on the whole cubic.
This contradicts \(abc\ne0\). A degree-six divisor with no
closed points of degree at most three is a single degree-six
point of multiplicity one. The dual restriction is thus
irreducible.

If \(a+b=0\), a prime triple in the indicated range has
\(p_1\equiv p_2\bmod p_3\). Hence
\(H_0(p_1,p_2,p_3)\equiv0\bmod p_3^2\), by its expression
\(- (p_1^3-p_2^3)^2+2p_3^3(p_1^3+p_2^3)-p_3^6\).
This violates the incident-prime condition in the
[exact reconstruction](ABC_COHERENT_NORM_PENCILS.md#exact-reconstruction-from-three-signed-coefficients).
With one coefficient zero, the two remaining primes cannot
satisfy the fixed linear relation once they exceed its
coefficients. For a nonhomogeneous affine relation, normalizing
by its nonzero constant identifies the affine plane with
a projective chart in the prime directions. The intersection
with the irreducible projective dual sextic remains irreducible.
Its affine degree is six because the line at infinity is
not a component of the dual sextic.

The parity hypothesis is essential. On \(3(P_1+P_2)=5P_3\),
substituting \((P_1,P_2,P_3)=(5U,5V,3(U+V))\) gives, up to
a nonzero constant, the product
\[
 (49U^2-2UV+49V^2)
 (196U^4-316U^3V-399U^2V^2-316UV^3+196V^4).
\]
This line has no triples of all odd primes. Thus changing only
the fixed linear relation in the unit-normalized construction
cannot obtain the large-square estimate for forms of factor
degree at most five.

## A squarefree progression with three fixed coefficients

**Theorem SW2.** For every sufficiently large \(X\), there are
\(\gg X^{1/3}/(\log X)^3\) distinct squarefree \(d\in[X,2X)\)
with \(d\equiv3\bmod4\), each admitting a genuine coherent
prime triangle with \((v_{12},v_{23},v_{13})=(45,9,36)\).
Its three primes lie at a common scale \(X^{1/6}\).
The two consecutive ideal-class differences are independent
in the three-torsion of \(\mathbb Q(\sqrt{-d})\).
Its reflected cubic field has fundamental discriminant \(d/3\),
with \(3\) completely split.

Take \(x=5,y=1,z=4\), set
\[
         h=270k,\qquad (p_1,p_2,p_3)=(p+h,p,p-h),
\]
and set \(A=p_1^3,\ B=16p_2^3,\ C=25p_3^3\).
Because \(5=1+4\), the Heron expression vanishes at \(h=0\),
where the three side lengths are proportional to \(1,4,5\).
It has a simple factor \(h\), and expansion gives
\[
 \begin{split}
 E&=\frac{4AC-(A+C-B)^2}{400}=81F(p,k),\\
 F(p,k)&=3k R(p,k),\\
 R(p,k)&=4p^5-5283p^4k+2371680p^3k^2-505197000p^2k^3
             +55269864000 p k^4-2694405870000 k^5 .              \tag{1}
 \end{split}
\]
On a sufficiently small positive fixed rectangle with \(k/p
\le1/5400\) and \(p,k\asymp Z\), \(R\asymp p^5\) and
\(F\asymp Z^6\). All three prime forms are then in a common
dyadic interval. Require \(p,k\) odd and \(p\equiv k\equiv1\bmod3\).
The polynomial gives \(F\equiv3\bmod4\), \(v_3(F)=1\),
and \(F/3\equiv1\bmod3\).

The full trace congruences follow from \(h\equiv0\bmod5\)
and \(h\) even: \(A-B\) is divisible by \(5\), and \(A-C\)
by \(4\). All incident coefficients are units for primes
greater than \(23\). At \(p_1,p_2,p_3\) respectively, the
Heron expression reduces to the negative squares of
\(184h^3,26h^3,8h^3\), up to the sign inside the square.
For distinct primes greater than \(23\) these are units.
Consequently all the prime-unit conditions of reconstruction
hold. For squarefree \(d=F(p,k)\) this reconstruction gives
\(\delta=9\) and the stated edge coefficients.

The cube-root ideals on edges \(12,23\) have norms
\(p_1p_2,p_2p_3\) and one-sided supports. The product reduces,
after removing the rational factor \(p_2\), to the one-sided
ideal on \(13\) of norm \(p_1p_3\). The product of the
conjugate of the first ideal and the second ideal has norm
\(p_1p_2^2p_3\) and one-sided supports as well. None of these
four ideals can be principal: a rational generator would
violate the support condition, and a nonrational integral
generator would have norm at least \(d/4\gg Z^6\).
Thus the two edge classes are independent and the reflected
cubic is irreducible by RF1. In RT1, \(s=\delta(A-C)/z\)
with \(z=4\) a 3-adic unit, so \(9\mid s\).
Since \(d/3\equiv1\bmod3\), RT1 gives the claimed fundamental
discriminant and splitting.

### The squarefree sieve on the prime progression

The form \(R\) is irreducible by reduction modulo \(7\).
For its degree-five reduction \(f(T)=R(T,1)\bmod7\),
\(\gcd(f,T^{49}-T)=1\), ruling out the factors of degrees
one or two that every reducible quintic must have. The
discriminant of \(R(T,1)\) is
\(2^{34}3^{49}5^{16}13^3 23^3\), by the resultant.
There is no repeated factor in \(F\) since \(R(p,0)=4p^5\).
At \(5,13,23\) require \(p\equiv1,\ k\equiv q\bmod q^2\)
for the specific prime \(q\). The three linear forms are
then units and \(v_q(F)=1\). Together with the conditions
at \(2,3\) this deals with the fixed bad primes. At every
other prime \(q\), the same choice \(p=1,k=q\bmod q^2\)
verifies that it is possible to avoid \(q^2\mid F\) while
all three forms are units.

Fix for the moment a large constant \(Y_0\). For every prime
\(q\le Y_0\) discard the residue pairs modulo \(q^2\) with
\(q^2\mid F\), in addition to the fixed congruences. By the
[Green–Tao prime forms theorem](https://arxiv.org/pdf/math/0606088)
the prime count for these finite conditions is asymptotic
to a positive constant times \(Z^2/(\log Z)^3\). Its
leading constant is bounded below uniformly in \(Y_0\).
In fact the three forms at primes outside \(2,3,5\) are
nonparallel distinct lines modulo \(q\), so their prime
local density is \(1+O(q^{-2})\). At each good \(q\) the
proportion of residue pairs modulo \(q^2\) with
\(q^2\mid F\) is \(O(q^{-2})\): the simple projective roots
give at most six primitive projective lifts, and the
pairs with both coordinates divisible by \(q\) have this
same order of density. Conditioning on the three line
values being units changes this by a bounded factor.
For each fixed \(Y_0\), applying the prime forms theorem
in each allowed class and adding gives the product of
these local factors. Their positivity checks and the
summable defects give a uniform positive lower constant.

For \(Y_0<q\le(\log Z)^8\), use the upper-bound sieve for
the three prime forms in each bad residue pair modulo
\(q^2\). In that residue pair write the new variables
in a box of side lengths \(\asymp N=Z/q^2\). At each
sieving prime \(\ell\nmid 2\cdot3\cdot5q\), the forbidden
union has \(3\ell-2\) elements modulo \(\ell\). For
squarefree sieving moduli \(m\), CRT gives the product
density and a remainder
\(O(3^{\omega(m)}(N+m))\), by counting complete residue
pairs and the box boundaries. Apply
[Richert, Theorem 9.1](sources/richert-sieve-methods-1976.pdf)
with sieve parameter \(D=Z^{1/20}\), supporting only
the indicated primes. The total of absolute remainders
at least common multiples of two integers at most \(D\)
is \(O(Z^{6/5})\), uniformly in the residue pair and
\(q\) in this range. The Selberg denominator is
\(\gg(\log Z)^3\): its weights at primes are
\(3/\ell+O(\ell^{-2})\). Taking products of primes up to
\(D^\epsilon\) gives total weight \(\asymp(\log D)^3\);
a fixed small \(\epsilon\) makes their weighted mean of the
logarithm of the product at most \(\tfrac12\log D\).
This last bound follows from
\(\sum_{\ell\le T}(\log\ell)/\ell\ll\log T\), and Markov's
inequality retains at least half the weight below \(D\).
The constants are uniform when the one prime \(q\) is
removed. Since there are \(O(q^2)\) bad residue pairs,
the total excluded in this prime range is
\[
        O\!\left(\frac{Z^2}{Y_0(\log Z)^3}\right)
           +o\!\left(\frac{Z^2}{(\log Z)^3}\right).              \tag{2}
\]

For \((\log Z)^8<q\le C Z\) with \(C\) containing the
rectangle, ignore primality but retain \(\gcd(p,k)=1\),
automatic for the prime values since \(0<k<p\).
Each primitive projective root modulo \(q^2\) gives a
lattice of determinant \(q^2\) with at most
\(O(1+Z^2/q^2)\) primitive pairs in the box. If they span
one dimension there are at most two primitive directions
with signs. If they span two, the lattice count follows
from successive minima, since both minima are \(O(Z)\).
Summing in this range gives
\(O(Z+Z^2/(\log Z)^8)\).
For \(q>CZ\), the factor \(k\) cannot be divisible by \(q\).
Apply
[Helfgott, Proposition 4.15](https://arxiv.org/pdf/math/0309109)
to the primitive variables and irreducible quintic \(R\).
It gives \(O_\epsilon(Z^{(5+\sqrt{113})/8+\epsilon})\),
a power saving from \(Z^2\).
Now choose \(Y_0\) a fixed constant large enough so that
(2) is less than half the uniform main-term lower bound.
This leaves \(\gg Z^2/(\log Z)^3\) squarefree pairs.

### Distinct fields and a cubic rank-growth consequence

The binary form \(F\) has nonzero discriminant and a real
linear factor \(k\). By
[Stewart–Xiao, Lemma 2.4](https://arxiv.org/pdf/1605.03427),
the pairs with \(0<|F|\ll Z^6\) representing a value
outside an essential representation orbit are
\(O_\epsilon(Z^{6(3/(4\sqrt6+3)+\epsilon)})\).
The exponent is strictly below two for sufficiently small
\(\epsilon\). The rational automorphism group is finite,
by the six distinct projective roots. Removing the
exceptional pairs thus leaves \(\gg Z^2/(\log Z)^3\)
distinct positive squarefree values. Choose the fixed
rectangle small enough around a point in its positive
region and scale \(Z\) as a constant times \(X^{1/6}\)
so that \(F(p,k)\in[X,2X)\) throughout. This is possible
by continuity and homogeneity of the positive value.
It proves the theorem at the stated field scale.

By [RF2](ABC_COHERENT_NORM_CUBIC_REFLECTION.md#line-sections-and-a-global-restriction-on-the-labels)
these reflected fields give cubic points from rational lines
on a fixed diagonal cubic isomorphic to
\(U^3+4V^3+5W^3=0\), which has a rational point.
Its Jacobian gains rank over all but finitely many of the
constructed fields. Indeed identify the plane cubic with
its Jacobian using that rational point. If rank did not grow
in a cubic field and \(P\) were its line-section point,
some positive multiple of \(P\) would be rational. Then
the differences of \(P\) with its conjugates are torsion
in the normal closure, of degree at most six. There are
only finitely many such torsion points by Northcott, and
the sum of the three conjugate points is the fixed class
of a line section. Thus \(3P\) can take only finitely many
values, and so can \(P\). This excludes at most finitely
many fields. The individual edge classes in the imaginary
field are nonzero, so CR1 also gives positive ranks on
each of the three fixed-colour Mordell twists at that \(d\).

For \((p,k)=(5479,1)\), the three primes in order are
\(5749,5479,5209\) and
\[
        d=46092984239139350979
           =3\cdot59\cdot613\cdot3049\cdot139329687071 .
\]
These factors are prime, and \(d/3\equiv1\bmod3\).
The expression (1), the mod-\(7\) irreducibility check,
and this factorization can be replayed with exact polynomial
arithmetic and `isprime` in PARI/GP.

The normalized class profile and common scale in this construction
are fixed. Its dyadic field result therefore concerns a genuine
fixed-profile family and does not give a saving in the occupied
rank-one count over growing profiles.

## Two independent points on one Mordell curve at squarefree scale

**Theorem SW3.** For all sufficiently large \(X\), there are
\(\gg X^{1/3}/(\log X)^3\) squarefree \(d\in[X,2X)\) with a
genuine coherent prime triangle of signed coefficients
\((8,1,1)\). Its primes are comparable to \(X^{1/6}\) with a
fixed ratio allowance. Its two consecutive ideal classes
are independent, and the same curve \(y^2=x^3-d/4\) has
rank at least two from the triangle. Its reflected cubic
field has fundamental discriminant \(3d\).

Use \(x=8,y=z=1\) and
\[
        p_1=11p+4k,\qquad p_2=11p-4k,\qquad p_3=p .
\]
Here \(A=p_1^3,\ B=p_2^3,\ C=(4p_3)^3\). The scaled
directions \((p_1,p_2,4p_3)\) lie on \(2(P_1+P_2)=11P_3\),
a rational line that cuts the Fermat dual in a quadratic
and a quartic factor. The normalized expression is the
identity
\[
 \frac{4AC-(A+C-B)^2}{64}
     =(20p^2+k^2)(263p^4-6496p^2 k^2-256k^4)=F_2(p,k).          \tag{3}
\]
The congruence \(8\mid A-B\) follows from \(p_1-p_2=8k\);
the other full congruences are automatic. At \(p_1,p_2\),
the Heron expression reduces to the negative square of
\((22^3-4^3)p^3\), and at \(p_3\) it reduces to the
negative square of \(2(4k)^3\).
Since \(22^3-4^3=2^3\,3^3\,7^2\), these are units at their
respective primes when the primes exceed \(7\). The
reconstruction applies with \(\delta=1\) on squarefree
\(d=F_2(p,k)\). Taking \(p,k\) odd and \(k\equiv0\bmod3\)
makes \(d\equiv3\bmod4\) and \(d\equiv1\bmod3\).
Choose a positive rectangle with \(k/p\le1/20\);
then \(F_2\asymp p^6\) and the three primes are of that
fixed common power scale.

The squarefree sieve proof above applies with degrees two
and four now. The quartic is irreducible modulo \(19\),
as its degree-four reduction has constant gcd with
\(T^{19^2}-T\). The discriminant of \(F_2(T,1)\) is
\(2^{54}3^{24}5\cdot7^{12}\cdot263\).
Include \(11\) among the fixed primes for the prime-form
conditions. At \(5,7,263\) take \((p,k)=(1,0)\bmod q^2\),
so the three forms are units and \(F_2=5260\bmod q^2\)
has valuation at most one. At \(11\) take
\((p,k)=(1,1)\bmod11^2\), where \(F_2\) is a unit and
the three forms are units. The conditions at \(2,3\)
already make \(F_2\) a unit. Away from these fixed primes
the first local choice \((1,0)\bmod q^2\) proves positivity.
For \(q>CZ\) the quadratic factor has size less than
\(q^2\) if \(C\) is chosen large enough, and the quartic
large-square bound is
[Helfgott, Proposition 4.12](https://arxiv.org/pdf/math/0309109).
All smaller-prime estimates are unchanged. The quartic
has a real projective root, by its opposite signs at
\((1,0)\) and \((0,1)\), so the same Stewart–Xiao
multiplicity estimate proves the distinct-field bound.
Shrink the rectangle and scale as for SW2.

The norm obstruction of SW2 gives independent ideal
classes since their possible principal products still
have norms \(O(Z^4)\) and \(d\asymp Z^6\).
The edge coefficients \(8,1,1\) are all cubes, so
[CR1](QUADRATIC_NORM_COLOUR_CLASS_RANK.md#from-one-colour-to-one-class-subspace)
gives independent points on the same Mordell curve,
with \(w=1\) and \(b=2,1\) respectively on the consecutive
edges. The reflected polynomial is the Fermat line
section with scaled marks \(p_1,p_2,4p_3\).
These are units at \(3\), with the first two congruent
and the third of opposite sign. The mixed-residue Hensel
calculation in RF3 uses just these unit and sign
conditions, so the local splitting field there has degree
two and the relative extension is unramified.
Away from \(3\) use RF1. This gives fundamental
discriminant \(3d\).

For example \((p,k)=(61,3)\) gives the primes
\((683,659,61)\) and squarefree
\(d=254836776019387=263\cdot283\cdot3423890903\).
The cubic field has discriminant \(764510328058161=3d\)
by exact PARI/GP `nfdisc` on (1) of RF1.

## One colour with all three primes in the same dyadic interval

**Theorem SW4.** For sufficiently large \(X\) there are
\(\gg X^{1/3}/(\log X)^3\) squarefree \(d\in[X,2X)\)
with genuine coherent prime triangles of signed coefficients
\((410^3,329^3,361^3)\), all three vertex primes in the same
dyadic interval at scale \(X^{1/6}\). Their two consecutive
ideal classes and the corresponding points on \(y^2=x^3-d/4\)
are independent. The reflected cubic field has fundamental
discriminant \(3d\).

Use the scaled prime directions
\[
       (P_1,P_2,P_3)=(329^2p_1,361^2p_2,410^2p_3),
       \qquad P_1+P_2=\lambda P_3,\qquad
       \lambda=\frac{143693}{84050}.
\]
The prime relation has right-hand coefficient
\(Q=\lambda 410^2=287386\). The weights
\(410=2\cdot5\cdot41,\ 329=7\cdot47,\ 361=19^2\)
are pairwise coprime.
An interior dyadic direction is
\[
            (p_1,p_2,p_3)=
               \left(\frac65,\frac32,\frac{3253707}{2873860}\right).
\]
The Heron value is positive at this direction by substitution.
On this line the binary form \(H_0\) has a quadratic and
quartic factor. One exact factorization, up to scalar, results
from \((P_1,P_2,P_3)=(143693U,143693V,84050(U+V))\). On
\(V=1\) its factors are
\[
 \begin{split}
 &126954678249\,U^2+171318643502\,U+126954678249,\\
 &44361550922906456001\,U^4-126458847552718837998\,U^3\\
 &\qquad +84685820122977118003\,U^2
            -126458847552718837998\,U+44361550922906456001 .
 \end{split}
\]
The quadratic has negative discriminant and the quartic
is irreducible modulo \(13\). The two factors have nonzero
discriminants and are coprime. The choice of slope comes from
the conjugate Fermat points with \(X_1+X_2=1,\ X_3=41/100\):
their first two squares are interchanged, and
\(X_1^2+X_2^2=\lambda X_3^2\).

To supply the full norm congruences, set
\((x,y,z)=(410^3,329^3,361^3)\).
At primes dividing \(410\) take \(p_3=1\) and
\(P_1=P_2=Q/2\) locally, except at \(2\) take
\(P_1=Q/2+4,\ P_2=Q/2-4\).
The derived \(p_1,p_2\) are units, and \(P_1^3-P_2^3\)
is divisible by the local part of \(x\). The choice at
\(2\) makes \(v_2(P_1^3-P_2^3)=3\) and \(v_2(P_3^3)=6\),
so the normalized \(E=H_0/(xyz)^2\) is odd and
\(E\equiv3\bmod4\). At \(q=7,47\), the valuation of
\(\lambda^3-1\) is two, and at \(19\) it is four. This
follows also from
\[
     143693^3-84050^3=3^4 7^2 19^4 31\cdot47^2\cdot67 .
\]
At each such \(q\), let \(\zeta\in\mathbb Z_q\) be the cube
root of unity congruent to \(\lambda\bmod q\), from Hensel's
lemma. Then \(v_q(\lambda-\zeta)=v_q(\lambda^3-1)\).
At \(7,47\) take \(P_3=410^2,\ P_2=\zeta P_3,\ 
P_1=(\lambda-\zeta)P_3\). This solves the relation,
makes the reconstructed primes \(q\)-adic units, and
gives \(P_3^3-P_2^3=0\). At \(19\) interchange the
two roles so that \(P_1^3-P_3^3=0\).
The full local congruences hold on open neighbourhoods
of these choices. At \(5,41\) one has \(A=B\) and
\(H_0=C(4A-C)\), with \(v_q(C)=v_q(x^2)\) and
\(4A-C\) a unit. At \(7,47\) one has \(B=C\) and
\(H_0=A(4B-A)\), with \(v_q(A)=v_q(y^2)\) and
\(4B-A\) a unit. The analogous calculation at \(19\)
uses \(A=C\). Consequently \(E=H_0/(xyz)^2\) is a
unit at each of these primes.

The lattice of integer prime directions of this primitive
linear relation has a basis over \(\mathbb Z\) that also spans
its local solution space at each prime. Impose the local
choices just described on its two parameters by CRT.
At \(3\) choose unit solutions; the scaled \(P_1,P_2\)
then agree modulo \(3\) and \(P_3\) has the opposite sign.
This makes \(E\equiv1\bmod3\).
At every other prime \(q\) outside the weight support,
\(q\ge11\). The dual Heron sextic over \(\mathbb F_q\)
has no line component: the coordinate-square map from
the smooth Fermat cubic is again injective on geometric
points and separable, and its irreducible image has degree
six. Thus the restricted binary sextic is nonzero modulo
\(q\). Its at most six projective roots and the at most
three vanishing prime directions leave an available unit
direction modulo \(q\).
Include all primes at most \(11\) and all bad primes
of the binary form and the linear forms in the fixed
set. Impose one of the unit directions at each prime of
this set outside the weight support and \(3\).
This deals with all fixed primes without removing any
square factor of \(E\).

Away from the fixed set, one can avoid the roots of
\(H_0\) and the three vanishing prime forms simultaneously
modulo \(q\), since their projective roots use at most
nine lines. The small- and intermediate-square sieve
and the quartic large-square bound of SW3 apply in the
two lattice parameters. A prime triple is a primitive
pair of such parameters, since any common divisor
would divide all three distinct primes. The binary form
has a real linear factor over \(\mathbb R\): on the
line the Heron expression is positive near the displayed
direction and negative at \(P_1=0\). Thus the same
essential-representation argument yields
\(\gg X^{1/3}/(\log X)^3\) distinct squarefree values
\(d=E\) in a dyadic window, by shrinking the
positive rectangle and scaling it as before.
Primes are larger than the fixed coefficients and the
prime divisors of the numerator of \(\lambda^3-1\);
the incident prime-unit condition follows by reducing
the Heron expression on the line at each vertex.
For example at \(p_1\), the opposite cube difference
is \((\lambda^3-1)P_3^3\), a unit. At \(p_3\) it is
the difference of two opposite unit cubes. Hence
reconstruction gives exactly the asserted common scale
and coefficients. The norm argument and CR1 give the
independent classes and points. The mixed-sign test on
the scaled prime directions at \(3\) is the same as in
SW3, so the reflected field has discriminant \(3d\).

An exact prime-direction control for this last relation is
\[
 \begin{split}
 p_1&=120227328909461958584273,\\
 p_2&=154773990478025346060233,\\
 p_3&=115467794261292522889201 .
 \end{split}
\]
PARI/GP certifies primality. Substitution verifies the linear
relation, the three full trace congruences, \(E>0\),
\(E\equiv3\bmod4\), \(E\equiv1\bmod3\) and
\(\gcd(E,p_1p_2p_3)=1\). This control uses the reconstruction
before squarefreeness of \(E\) is imposed.
