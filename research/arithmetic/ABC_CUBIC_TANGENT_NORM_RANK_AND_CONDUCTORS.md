# A cubic tangent family with one norm class and two independent points

For a fixed coherent prime progression, a tangent square reduces
the squarefree field parameter to degree four. In this family
the consecutive cubic norm classes always coincide. The two
corresponding Mordell points can still be independent, as two
fixed tests modulo \(2\) show. This gives a lower count for
rank-two Mordell curves by conductor.

## The tangent identity and the lost cubic direction

Write \(p_1=p-h,\ p_2=p,\ p_3=p+h\) with \(0<h<p\), and put
\[
 \begin{split}
 R(p,h)&=3p^4-6p^2h^2-h^4,\qquad E=h^2R,\\
 u_-&=2p^3-3p^2h-h^3,\qquad
 u_+=2p^3+3p^2h+h^3 .
 \end{split}
\]
These are the trace and parameter formulas from
[three-coefficient reconstruction](ABC_COHERENT_NORM_PENCILS.md#exact-reconstruction-from-three-signed-coefficients)
with \((x,y,z)=(1,1,2)\), so the Heron inputs are
\(A=p_1^3,\ B=4p_2^3,\ C=p_3^3\).
The identity \(4AC-(A+C-B)^2=4h^2R\) holds for all \(p,h\).

**Proposition CT1.** If the three \(p_i\) are distinct
primes greater than \(3\), these formulas yield
a genuine coherent triangle over \(\mathbb Q(\sqrt{-R})\).
The two consecutive virtual units have the same cube class.
For \(R>0\) with squarefree part \(d>3\), the designated
reflected binary cubic is reducible.

The three primes exceed \(3\), so \(h\) is divisible by
\(6\) and \(v_3(R)=1\). Thus the quadratic field is
well defined in either sign. Also \(2\mid A-C\). At
either outer prime the Heron numerator reduces to the negative square of the difference
between \(4p^3\) and \(8p^3\), while at \(p\) it is the
negative square of the difference between opposite odd cubes.
Thus \(\gcd(E,p_1p_2p_3)=1\).
For \(Y^2=-R\), the generators are
\[
    \alpha_-=(u_-+hY)/2,\qquad
    \alpha_+=(u_++hY)/2,\qquad
    \alpha_{13}=p^3-3ph^2+hY .
\]
They have integral traces, norms \((p_1p_2)^3,(p_2p_3)^3,
(p_1p_3)^3\) and \(\alpha_-\alpha_+=p^3\alpha_{13}\).
The traces and incidence units give one-sided prime
factorization in either quadratic sign. A further identity is
\[
             \alpha_-\,\overline{\alpha}_+
                   =\left(\frac{Y-p^2-h^2}{2}\right)^3.          \tag{1}
\]
It follows by expanding both sides and substituting \(Y^2=-R\).
The rational norm of \(\alpha_+\) is already a cube, so the
two virtual-unit classes agree. The coefficients when \(R\)
is squarefree are \((h,h,2h)\) in the doubled trace convention.

The reflected cubic reduction also has a direct geometric
description. The projective line with normal
\((p-h,p,p+h)\) meets the diagonal cubic
\(U_1^3+U_2^3/4+U_3^3=0\) at the fixed rational point
\([1:-2:1]\). This is the line section in
[RF2](ABC_COHERENT_NORM_CUBIC_REFLECTION.md#line-sections-and-a-global-restriction-on-the-labels).
When the two consecutive class differences are nonzero,
the three vertex labels fill an affine line in their
three-torsion coset.

## Rank from a second descent

Consider
\[
         M_{p,h}:\ y^2=x^3-E/4,\qquad
         P_-=(p(p-h),u_-/2),\quad P_+=(p(p+h),u_+/2).
\]

**Proposition CT2.** Suppose \(h\equiv1\bmod55\),
\(p\equiv42\bmod55\), and \(M_{p,h}(\mathbb Q)[2]=0\).
Then \(P_-,P_+\) are independent in the free part of
\(M_{p,h}(\mathbb Q)\).

For \(q=5,11\) the reductions have good reduction and
a unique nonzero rational two-torsion point \((r_q,0)\),
where \(r_q^3=E/4\bmod q\). Use the homomorphism into
\(\mathbb F_q^\times/\mathbb F_q^{\times2}\) given on
points away from \(O,(r_q,0)\) by \(x-r_q\). For example,
the product of these differences along a nonvertical
line is the square of the ordinate of that line at \(r_q\).
This gives the homomorphism, with the value at the
two-torsion point obtained by the local leading term
of \(x-r_q\). The two evaluations here are
\[
\begin{array}{c|c|c|c|c}
 q&p\bmod q&r_q&x(P_-)-r_q&x(P_+)-r_q\\ \hline
 5&2&3&4&3\\
 11&-2&9&8&4 .
\end{array}
\]
The squareclass rows are \((0,1),(1,0)\). Odd torsion
vanishes in both homomorphisms. If there were a dependence
modulo torsion, divide its two integer coefficients by
their gcd to make at least one odd. Applying the two
homomorphisms contradicts the dependence. This detects
rank two even though the cubic direction in (1) agrees.

There is a general restriction on recovering the lost
direction by another unramified norm function on this
same genus-one parameter curve.

**Proposition CT3.** Let \(C/\mathbb Q\) be a smooth
genus-one curve with a double cover of \(\mathbb P^1\).
The subgroup of \(\mathbb Q(C)^\times/\mathbb Q(C)^{\times3}\)
represented by functions with divisor divisible by \(3\)
and norm a cube in \(\mathbb Q(t)\) has dimension at most
one over \(\mathbb F_3\).

Take the rational divisor class of
\(\operatorname{div}(f)/3\) in \(\operatorname{Jac}(C)[3]\).
If it vanishes, its rational divisor is principal over
\(\mathbb Q\), by descent of principal divisors. Hence
\(f=c g^3\) for some \(g\in\mathbb Q(C)^\times\) and
\(c\in\mathbb Q^\times\). The norm condition makes
\(c^2\) a cube in \(\mathbb Q(t)\), so \(c\) is a rational
cube. This proves injectivity. The rational three-torsion
of an elliptic curve has dimension at most one by the
Weil pairing and \(\mu_3\not\subset\mathbb Q\).
For \(Y^2=-R(t,1)\) the divisors of the two consecutive
functions have orders divisible by \(3\), since the
finite norm zeros are one-sided and their two poles at
infinity have order \(3\). Thus this obstruction applies
to their parameter curve. It concerns the unramified
norm functions and leaves open other sources of special
ideal classes in the quadratic fibres.

## Squarefree dyadic fields and rank-two conductors

**Theorem CT4.** For each sufficiently large \(X\) and each
choice of sign of \(R\), there are
\(\gg X^{1/2}/(\log X)^3\) distinct squarefree
\(d=|R|\in[X,2X)\) from the above genuine triangles with
\(h=6k\), where all three prime vertices and \(h\) are
within fixed positive multiples of \(X^{1/4}\).
The curves \(M_{p,h}\) have rank at least two. For \(R>0\)
all three vertices can be in one dyadic prime window,
and each edge determines a nonzero ideal class of order
three in \(\mathbb Q(\sqrt{-d})\). The consecutive ideal
classes coincide. With \(k\) also required to be prime,
the same assertion holds with log power \(4\).

Choose a fixed small rectangle of positive \((p,k)\) at
scale \(Z\), with \(h=6k<p\). For \(R>0\) use a ratio
\(0<h/p<1/3\), and for \(R<0\) use \(h/p\) sufficiently
close to \(1\). The quartic has both signs on these
rectangles, and \(|R|\asymp Z^4\). By shrinking each
rectangle and then fixing \(Z\) as a suitable constant
multiple of \(X^{1/4}\), homogeneity gives
\(|R|\in[X,2X)\) throughout. In the positive cone the
ratio of the largest prime to the smallest is less than
two, uniformly on the rectangle, which can be shrunk
to put every vertex in one common dyadic prime window.

Impose \(p\equiv42,\ k\equiv46\bmod55\), so that CT2
applies. Require \(p,k\) odd and prime to \(3\), which
gives \(R\equiv3\bmod4\) and \(v_3(R)=1\).
At \(5,11\) the rank tests already give unit \(R\) and
unit prime forms. The same unit conditions hold at \(7\) by
taking \(h=1,\ p=2\bmod7\). At any other prime
\(q\ge13\), fix a unit \(k\bmod q\). The three
linear forms in \(p\) and the quartic exclude at most
seven \(p\) residues, so there is a direction with
unit prime forms and unit \(R\). These choices also
work when \(k\) itself is a fourth prime form.
The quartic \(R\) is irreducible by Eisenstein at \(3\)
in \(R(1,h)\), and \(R(t,1)\) has discriminant
\(-2^{12}3^3\). The substitution \(h=6k\) does not
introduce any other bad reduction prime.

Use the [squarefree sieve of SW2 and SW3](ABC_COHERENT_PRIME_TRIANGLES_SQUAREFREE_WINDOW.md#distinct-fields-and-a-cubic-rank-growth-consequence)
with the fixed irreducible quartic \(R(p,6k)\).
Here are the bounds with \(s=3\) or \(4\) prime forms.
First impose \(q^2\nmid R\) for \(q\le Y_0\), where
\(Y_0\) is fixed. The linear forms have finite complexity
and the choices above are locally admissible.
The finite-complexity prime forms theorem of
[Green–Tao](https://arxiv.org/pdf/math/0606088), completed by
[Green–Tao–Ziegler](https://arxiv.org/pdf/1009.3998),
gives an asymptotic of order \(Z^2/(\log Z)^s\).
The arbitrary finite local exclusions still give a
uniform positive lower constant as \(Y_0\) grows.
At a good prime the prime-form local factors are
\(1+O(q^{-2})\), and an \(O(q^{-2})\) fraction of the
unit pairs are removed by the primitive simple roots
modulo \(q^2\).

For \(Y_0<q\le(\log Z)^8\), there are \(O(q^2)\)
bad residue pairs modulo \(q^2\). In each pair the
[Selberg upper bound, Richert Theorem 9.1](sources/richert-sieve-methods-1976.pdf)
as used in SW2 holds with \(s\) in place of \(3\):
the sieve weights at good primes are
\(s/\ell+O(\ell^{-2})\). The boundary errors averaged
over the two divisor variables with cutoff \(Z^{1/20}\)
are \(O_s(Z^{6/5})\) per residue pair, since the
least common moduli are at most \(Z^{1/10}\) and
the root weights \(s^{\omega(m)}\) are subpower.
Summing gives
\(O_s(Z^2/(Y_0(\log Z)^s))+o(Z^2/(\log Z)^s)\).
In the range \((\log Z)^8<q\le C Z\), the primitive
pairs in each simple root lattice of determinant \(q^2\)
number \(O(1+Z^2/q^2)\). Primitivity holds since \(p\) is
prime, \(1\le k<p\). This range contributes
\(O(Z+Z^2/(\log Z)^8)\).
For \(q>C Z\) use
[Helfgott, Proposition 4.12](https://arxiv.org/pdf/math/0309109)
for large square factors of an irreducible binary
quartic at primitive pairs. Its \(O(Z^{4/3}(\log Z)^A)\)
bound suffices. Take \(Y_0\) large against the positive
main-term constant. This leaves
\(\gg Z^2/(\log Z)^s\) squarefree pairs.

To distinguish \(d\), apply
[Stewart–Xiao, Lemma 2.4](https://arxiv.org/pdf/1605.03427)
to the quartic, of nonzero discriminant and with
a real projective root by the two sign regions.
Of pairs representing nonzero values \(O(Z^4)\),
at most \(O_\epsilon(Z^{12/7+\epsilon})\) are outside
the essential representation orbits. The rational
automorphism group is finite because the form has
four distinct projective roots. This yields the
asserted count of distinct \(d\) for a fixed sign.

For all large squarefree \(d\) in the family there is
a prime \(q\ge5\) dividing \(R\) and not \(6k\), because
\(\gcd(R,k)=1\) follows from
\(R\equiv3p^4\bmod k\) and \(\gcd(p,k)=1,\ 3\nmid k\).
In particular \(E/4=9k^2R\) cannot be a cube. This
checks the two-torsion hypothesis in CT2. In the
imaginary sign, principality of a supporting
prime-ideal quotient would make its virtual unit a
cube up to a unit. The units are cubes when \(d>3\).
A cube root would be integral with norm \(p_i p_j\).
A nonrational integral element has norm at least
\(d/4\), which exceeds \(p_i p_j\asymp Z^2\), and a
rational integral element has square norm. Both are
impossible for distinct primes. This gives exact
order three.

## Counting the conductors

**Corollary CT5.** In the prime-\(k\) family of CT4,
for each sign there are \(\gg Y^{1/5}/(\log Y)^4\)
distinct conductors at most \(Y\) of these
rank-at-least-two curves.

The model
\(y^2=x^3-9k^2 R\) is integral. At every prime \(q\ge5\)
dividing \(k d\) the coefficient of \(x^0\) has
valuation \(1\) or \(2\). The model is minimal and
[Tate's algorithm](https://wstein.org/Tables/antwerp/tate/tate.pdf)
gives types II and IV, respectively, and conductor
exponent \(2\). Other primes \(q>3\)
give good reduction. At \(2,3\) the conductor exponents
are bounded by the fixed valuations of the model
discriminant. Thus the part of the conductor coprime
to \(6\) is \((k d/3)^2\), and the full conductor is
\(\asymp Z^{10}\). For a given conductor there are
\(O(1)\) possibilities for \(k\), because it is a prime
\(\asymp Z\) dividing its determined large-prime
square root, of size \(\asymp Z^5\). Once \(k\) is
specified and the sign fixed, \(R\) is determined and
the quartic equation for \(p\) has at most four
solutions. The fibres of the conductor count are
therefore bounded. Set \(Z\) to a sufficiently small
constant multiple of \(Y^{1/10}\).

The imaginary sign has an exact prime-\(k\) example
\((p,k)=(4057,211)\) with
\[
 R=651870845595003=3\cdot2414219\cdot90004379.
\]
The real sign has \((p,k)=(3617,431)\) with
\[
 R=-56185743375717=-3\cdot290219\cdot64532581.
\]
The four prime values in each case and the displayed
factorizations are certified by exact PARI/GP
primality and factorization. The rank proof is the
pair of residue rows in CT2.

The Mordell coefficient in CT4 keeps the varying
factor \(h^2\), and the field count by \(d\) and
conductor count by \(k d\) have different scales.
These tangent triangles lie in the reducible
reflected sector. The construction leaves the
aggregate occupied rank-one count and the bounds
for irreducible reflected paths open.

## Independence from the full simple prime support

**Proposition CT6.** Fix \(\eta>0\). For coprime positive
integers \(p>h\) with \(h\) even, define
\[
                 D_1=\prod_{\substack{q\ge5\\v_q(R)=1}}q .
\]
If \(D_1\ge p^\eta\) and \(p\) is sufficiently large in
terms of \(\eta\), the two displayed Mordell points are
independent, with no auxiliary congruence restriction.
In particular this holds for all sufficiently large
primitive pairs of this parity with squarefree \(R\).

Remove sixth powers from the integer coefficient
\(b=-h^2R/4\) to obtain a sixth-power-free parameter \(b_0\).
If \(q\ge5\) divides \(R\), it cannot divide \(p\), \(h\)
or \(p\pm h\): substitution in the quartic would give
a nonzero fourth power times \(1,3\) or \(4\) modulo \(q\).
Hence \(D_1\mid b_0\), with simple exponent at each of
its primes. Also \(|b_0|\ll p^6\). In the rescaled model
the logarithmic \(x\)-heights of the points remain
\(O(\log p)\), since their original abscissas are
\(O(p^2)\) and the denominator introduced by removing
sixth powers is \(O(p^2)\).
[Voutier–Yabuta, Theorems 1.2 and 1.3](https://arxiv.org/pdf/1305.6560)
give \(\widehat h(P_\pm)=O(\log p)\) and
\(\widehat h(T)\gg_\eta\log p\) for every non-torsion
\(T\) on the same Mordell curve. The lower bound
uses \(|b_0|\ge D_1\).

If the two points are dependent in the free part,
there is a relation \(aP_-+b'P_+\) in the torsion
subgroup, with \((a,b')\ne(0,0)\) and both
coefficients \(O_\eta(1)\).
If one point is torsion this is immediate. Otherwise
their common saturated rank-one subgroup has a
generator \(T\) modulo torsion. The two coefficients
of \(T\) in \(P_-,P_+\) are bounded by the square roots
of the height ratios, which bounds the coefficients
of a primitive relation.

For \(q>20\) dividing \(D_1\), the model is minimal at
\(q\) and its reduction is the cusp \(y^2=x^3\).
Both points specialize to nonsingular points with
nonzero finite coordinates. The nonsingular group
is additive, with coordinate \(x/y\) extended by
zero at infinity. Rational torsion has order prime
to \(q\). Indeed at the fixed prime \(5\) the
prime-to-\(5\) torsion injects into the quotient
by the formal group. On the sixth-power-free
Mordell model at \(5\) the good reduction has
order \(6\), and in the additive case the
nonsingular group has order \(5\) and the component
group has order at most \(4\) by Tate's algorithm.
Thus a rational torsion point contributes zero to
the additive group at \(q\).
Clearing the unit denominators and cancelling \(2p\)
in the reduced relation yields
\[
 q\mid L_{a,b'}(p,h),\qquad
 L_{a,b'}=(a+b')(2p^4-3p^2h^2-h^4)
          +(a-b')ph(p^2+h^2).
\]
This quartic cannot be a multiple of \(R\) for
\((a,b')\ne(0,0)\): comparing odd monomials first
gives \(a=b'\), and the three remaining coefficients
are not proportional to those of \(R\).
Since \(R(t,1)\) is irreducible, the resultant
of \(R(t,1)\) and \(L_{a,b'}(t,1)\) is nonzero.
As \(q\nmid h\), every \(q>20\) in \(D_1\)
divides that integer. The bounded relation
coefficients leave only a finite list of such
resultants. This forces \(D_1=O_\eta(1)\) in any
dependent pair, contrary to its lower bound for
large \(p\).

For squarefree \(R\), \(D_1\ge |R|/6\). Away from a
fixed neighborhood of the unique real root
\(\rho>1\) of \(R(t,1)\), \(|R(p,h)|\gg p^4\).
Inside it \(h\asymp p\); since the root is simple,
\(|R|\gg h^4 |p/h-\rho|\).
Roth's theorem with approximation exponent \(5/2\)
gives \(|R|\gg p^{3/2}\) there. Thus the hypothesis
on \(D_1\) holds with \(\eta=1\) for all large
squarefree primitive pairs. For the dyadic boxes
of CT4, \(|R|\asymp p^4\) already gives this
without rational approximation.
