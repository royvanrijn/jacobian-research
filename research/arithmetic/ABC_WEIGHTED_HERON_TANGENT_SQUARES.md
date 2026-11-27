# Tangent squares on weighted Heron sections

A fixed positive prime line can remove at most one
polynomial square from a weighted Heron form. When a
tangent square occurs, the remaining form is typically
irreducible at the exponent prime itself. In degree
five there is also a rational-point obstruction to
a second rational branch point.

## Rational tangency and the incident-prime condition

Fix an odd prime exponent \(\ell\) and positive rational
weights \(c_1,c_2,c_3\). The weighted Heron branch curve
\(\mathcal H\) is the coordinate-square image of
\[
              \mathcal F:\quad
              c_1 X_1^\ell+c_2 X_2^\ell+c_3 X_3^\ell=0.
\]
Its homogeneous equation in square coordinates \(U_i\)
is \(2\sum_{i<j}A_i A_j-\sum_i A_i^2=0\), where
\(A_i=c_i^2 U_i^\ell\). The image map is injective on
geometric points by the sign argument in
[HP1](ABC_COHERENT_HIGHER_PRIME_PROGRESSION_FIELDS.md#the-hyperelliptic-cover-and-its-power-classes).
It is an immersion away from the coordinate axes.

**Proposition HT1.** Let \(L:\sum_i\lambda_i U_i=0\)
be a rational line with all coefficients nonzero,
meeting the positive real cone and avoiding the
coordinate-axis points of \(\mathcal H\).
The restriction of the weighted Heron form to
\(L\) is separable or is a rational linear factor
squared times a separable degree-\((2\ell-2)\) form,
with coprime factors in the latter case. There is
at most one tangency, and its preimage on
\(\mathcal F\) is a rational point. For a fixed
weighted slice with \(U_i=s_i^2p_i\) and fixed rational
nonzero \(s_i\), an axis intersection instead
violates the incident-prime unit condition for
every sufficiently large triple of distinct primes.

The diagonal conic \(\sum_i\lambda_i X_i^2=0\)
pulls back \(L\). At a tangent point \(x\) its
normal must satisfy \(\lambda_i\propto c_i x_i^{\ell-2}\).
Take \(m=\ell-2\) and \(r_i=x_i/x_3\) for
\(i=1,2\); then \(r_i^m\in\mathbb Q^\times\) and
\(\lambda_1 r_1^2+\lambda_2 r_2^2+\lambda_3=0\).
Over \(\mathbb Q(\mu_m)\), each \(r_i^2\) is an
eigenvector of a Kummer character in the joint
splitting field. Any nontrivial character in
the last sum must occur twice, and the constant
term also needs another trivial term. Thus both
characters are trivial. As \(m\) is odd,
\(r_i\in\mathbb Q(\mu_m)\), so each \(r_i\) is a
rational number times an \(m\)-th root of unity.
Indeed its real radical representative is
also in this abelian extension and all its
conjugates are real, so it is the unique real
root and is rational.

A relation among three such roots of unity with
nonzero rational coefficients either uses the
same root throughout or has three equal
coefficients on a rotated primitive cube-root
triple. To see this, normalize one root to \(1\).
For \(\alpha+\beta u+\gamma v=0\), taking squared
absolute values gives \(u+u^{-1}\in\mathbb Q\),
hence an integer. For a root of odd order this
forces \(u=1\) or \(u\) primitive of order \(3\).
The same holds for \(v\), and the only
nonconstant relation then has \(v=u^2\) and all
coefficients equal. In the tangent relation the
coefficients have mixed signs even after the
rational coordinate squares are absorbed into
them. Thus the roots of unity are all trivial,
and the tangent point is rational. There is at
most one such point because the rational
\(m\)-th roots of the two prescribed ratios
are unique.

Rescale the normal so that
\(\lambda_i=c_i x_i^{\ell-2}\). The Hessian of
the difference of the Fermat equation and
\(\ell/2\) times the conic equation is
\(\ell(\ell-2)\operatorname{diag}(\lambda_i)\)
at \(x\). The tangent plane to the cone is
\(x^\perp\) for this nondegenerate form, and
the quotient \(x^\perp/\mathbb Q x\) is
one-dimensional and nondegenerate. Thus the
contact at the tangent point has exactly order
two. All other intersections are transverse.
This proves the restriction statement by
Bézout and the square map.

If \(U_i=0\) is an axis mark on \(L\cap\mathcal H\),
then \(c_j^2(U_j/U_k)^\ell=c_k^2\)
at \(U_j/U_k=-\lambda_k/\lambda_j\).
Reducing the linear prime relation modulo \(p_i\)
forces this same ratio there for all sufficiently
large primes. The Heron expression at \(p_i\)
then reduces to the negative square of an
opposite-side difference which vanishes.
This fails the required incident-prime unit,
as in
[the reconstruction test](ABC_COHERENT_NORM_PENCILS.md#exact-reconstruction-from-three-signed-coefficients).
A line with one coefficient zero cannot contain
arbitrarily large distinct primes in the two
other weighted coordinates because it fixes
their ratio. Thus the axis exception cannot be
used to extract further squares in a genuine
fixed-profile prime line.

## The local degree restriction on the residual form

Rescale by the rational tangent coordinates
\(U_i=x_i^2 V_i\). The signed Fermat coefficients
become \(c_i x_i^\ell\) and add to zero.
After permuting and scaling, they are
\(a,b,-(a+b)\) with coprime positive integers
\(a,b\). The tangent line is
\(aV_1+bV_2=(a+b)V_3\).
Set \(V_1=p+b h,\ V_2=p-a h,\ V_3=p\), and define
\[
 \begin{split}
 A&=a^2(p+bh)^\ell,\qquad B=b^2(p-ah)^\ell,\qquad
 C=(a+b)^2p^\ell,\\
 Q_{a,b}(p,h)&=\frac{4AB-(A+B-C)^2}{a^2 b^2(a+b)^2 h^2}.
 \end{split}                                                   \tag{1}
\]
The numerator is divisible by \(a^2,b^2,(a+b)^2\)
even as a universal polynomial: modulo each
linear factor one of \(A,B,C\) vanishes to order
two and the other two coincide. It is divisible
by \(h^2\) by tangency at \(h=0\). These distinct
polynomial factors give the displayed integral
quotient. Its value at \((1,0)\) is
\(\ell(\ell-2)\). For this last coefficient expand
the formal square roots at \(h=0\): the second
order term of
\(a(1+bh)^{\ell/2}+b(1-ah)^{\ell/2}-(a+b)\)
is \(\ell(\ell-2)ab(a+b)h^2/8\).
The other factors of the Heron product
have constant product \(8ab(a+b)\).

**Proposition HT2.** If \(\ell\nmid ab(a+b)\),
the residual form \(Q_{a,b}\) is irreducible over
\(\mathbb Q\) by Eisenstein at \(\ell\). If
\(\ell\mid ab(a+b)\) and the valuation of its
unique divisible weight among \(a,b,a+b\) is
\(k\) with \(\ell\nmid k\), the residual form over
\(\mathbb Q_\ell\) has exactly two irreducible
factors, of degrees \(\ell-2\) and \(\ell\).
Any factorization over \(\mathbb Q\) must have
these same degrees.

In the first case the numerator in (1) is the
\(\ell\)-th power modulo \(\ell\) of the same
expression with exponent \(1\).
At exponent \(1\) the numerator is
\(-a^2b^2(a+b)^2 h^2\) on the tangent line.
After dividing the denominator in (1) this
gives \(Q_{a,b}\equiv-h^{2\ell-2}\bmod\ell\).
The constant coefficient of \(Q_{a,b}(1,h)\)
has valuation one, so Eisenstein applies.

In the second case permute the signed weights
to put the divisible weight at \(a\), and
absorb the unit \(b\) into \(h\) so it is enough
to take \(b=1\) with \(a\in\mathbb Z_\ell\).
Signs of the weights can change under this
normalization. These changes preserve factor
degrees, and the normalized form again has (1).
Set \(a=0\) after taking the polynomial
quotient. The expansion before division by
\(h^2\) gives
\[
  4p^\ell(p+h)^\ell
      -(2p^\ell+\ell h p^{\ell-1})^2.
\]
Consequently \(Q_{a,1}\equiv
4p^\ell h^{\ell-2}\bmod\ell\). In \(Q_{a,1}(1,h)\)
the roots reducing to \(h=0\) form one Newton
segment of length \(\ell-2\) and rise \(1\),
so give an irreducible factor of that degree.
For the \(p=0\) cluster put \(h=1\) and give
\(p\) weight \(2k/\ell\). The initial terms of
\(A,B,C\) are \(a^2,p^\ell,p^\ell\), all of
weight \(2k\). The initial terms of \(B-C\)
cancel. The initial numerator of (1) is
\(4a^2p^\ell-a^4\), so the second segment
has length \(\ell\) and rise \(2k\). Its
reduced slope denominator is \(\ell\) when
\(\ell\nmid k\), proving the other
irreducible local factor.

**Corollary HT3.** For \(\ell=5\), if the unique
even weight among \(a,b,a+b\) has \(2\)-adic
valuation \(1\) or \(2\), every rational factor
of \(Q_{a,b}\) has even degree. In conjunction
with HT2 this proves irreducibility whenever
the \(5\)-adic weight condition of either case
of HT2 holds.

Put the even weight at \(a\) and take \(b=1\)
by a \(2\)-adic unit rescaling. Write \(k=v_2(a)\)
and give \(p\) weight \(k/4\) in \(h=1\).
The leading terms of \(A\) and \(B-C\) are
\(a^2\) and \(-5ap^4\), respectively, each of
weight \(2k\). In the numerator the cross
term has extra weight \(1\), and \(4AB\)
has weight \(2+2k+5k/4>4k\) for \(k=1,2\).
The initial numerator is therefore
\(-a^4-25a^2p^8\). After normalization this
gives a single \(2\)-adic Newton segment of
length eight and slope denominator \(4\)
or \(2\), so all factor degrees are even.
In the exceptional case \(\ell\mid ab(a+b)\)
under HT2, the only possible rational
factorization has degrees \(3,5\), which
are incompatible with this local condition.

For example, at \(\ell=5\) and \((a,b)=(1,1)\),
\[
        Q_{1,1}(p,h)=15p^8-80p^6h^2-110p^4h^4-20p^2h^6-h^8.
\]
It is Eisenstein at \(5\) after setting \(p=1\).
At \((a,b)=(5,1)\), the degrees over
\(\mathbb Q_5\) are \(3,5\), and HT3 still
gives irreducibility over \(\mathbb Q\).
The two local factor degrees can also be
replayed by applying PARI/GP
`factorpadic` at \(5\) to (1) with \(p\)
as variable and \(h=1\).

## A rational branch point obstruction at exponent five

**Proposition HT4.** For \(\ell=5\) and positive
rational weights \(a,b\), the residual form on
the tangent line has no rational projective root.

The diagonal conic and quintic of the normalized
tangent line are
\[
 aX^2+bY^2=(a+b)Z^2,\qquad
 aX^5+bY^5=(a+b)Z^5.
\]
A rational point on their intersection cannot
have \(Z=0\), by positivity, and cannot have
\(X=0\) or \(Y=0\), by comparing a second and
fifth power of the same rational ratio.
With \(Z=1\), if either coordinate is \(1\)
or \(-1\), both equations force \(X=Y=1\).
The same conclusion holds if \(X=Y\).
At any other point,
\[
 \frac{X^5-1}{X^2-1}
   =\frac{Y^5-1}{Y^2-1}.
\]
Using \((t^5-1)/(t^2-1)=t^3+t+1/(t+1)\)
and dividing by \(X-Y\) gives
\[
          (X^2+XY+Y^2+1)(X+1)(Y+1)=1.                   \tag{2}
\]
Writing \(s=X+Y,\ z=XY\), this implies
\[
 \bigl(2z-s^2+s\bigr)^2=s(s+1)(s^2+s+4).
\]
For \(s\ne0\), the substitution
\(x=4/s,\ y=4(2z-s^2+s)/s^2\) is on
\[
                   E_0:\ y^2=x^3+5x^2+8x+16.
\]

The rational group of \(E_0\) consists of
\(O,(-4,0),(0,\pm4)\). Here is a short descent.
After translating \(x\) by \(4\) the model
is \(y^2=x(x^2-7x+16)\), with \(2\)-isogenous
model \(y^2=x(x-1)(x+15)\). Use the usual
\(x\)-coordinate squareclass map for the
two curves. The first image can only
contain \(1,2\), because the real quadratic
factor of the model is positive and the
squarefree classes divide \(16\).
Class \(2\) would require coprime integers
\(u,v\) with
\(n^2=2u^4-7u^2v^2+8v^4\).
If \(u\) is odd this is impossible modulo
\(8\). If \(u\) is even and \(v\) odd it is
impossible modulo \(16\).
The second image contains \(1,-15,-3,5\),
using the point of abscissa \(-3\) and
the nonzero isogeny kernel point.
All possible classes divide \(-15\) and
the remaining coset is excluded by class
\(-1\), which would require
\(n^2=-u^4+14u^2v^2+15v^4\).
Opposite parities fail modulo \(8\), and
both odd fail modulo \(16\).
The two isogeny indices have product \(4\),
so the rank is zero. Indeed the dual
kernel point on the second curve has no
rational preimage under the first isogeny,
since the other two kernel points of
multiplication by \(2\) on the first
curve are nonrational.
Reduction at \(7,11\) bounds the torsion
order by \(8\), and the group has a unique
rational point of order two. An order-eight
point would have to halve \((4,4)\) or
\((4,-4)\) on the translated model.
The duplication formula for this abscissa
requires a rational root of
\((x^2-12x+16)(x^2-4x+16)\), impossible
by the two discriminants. The torsion
therefore has order \(4\), as claimed.

The only possible finite nonzero \(s\)
from \(E_0(\mathbb Q)\) is \(s=-1\),
which gives \(z=1\) and has no rational
pair \(X,Y\). For \(s=0\), equation (2)
gives \(z=0\), which does not satisfy the
positive-weight conic. Thus the only
rational intersection is the tangent point.
Since coordinate squaring is injective
on geometric points, a rational
intersection in the square coordinates
would also have a rational preimage.
The residual form therefore has no
rational projective zero.

## The only possible rational factor degrees at exponent five

**Theorem HT5.** For coprime positive integers \(a,b\)
and \(\ell=5\), the residual octic \(Q_{a,b}\) is
irreducible over \(\mathbb Q\) or has exactly two
irreducible factors of degrees \(3,5\). In the
reducible case \(5\mid ab(a+b)\).

Only the \(5\)-adic case omitted by HT2 requires
a further local analysis. As there, put the
divisible weight at \(a\) and absorb \(b\) into
the variable. Suppose \(a=5^{5m} a_0\) with
\(m\ge1\) and \(a_0\) a \(5\)-adic unit.
The length-three segment at \(h=0\) is unchanged.
In the length-five cluster put \(h=1\) and
\(p=5^{2m}u\). Divide the numerator by \(5^{20m}\)
and by \(4a_0^2\). Its reduction to precision
\(5^{2m+1}\) is the polynomial \(u^5-a_0^2/4\).
In fact the three scaled inputs, after division
by \(5^{10m}\), reduce to \(a_0^2,u^5,u^5\).
The approximation holds for coefficients,
since \((1+5^{2m}u)^5\equiv1\bmod 5^{2m+1}\)
and the other changes have at least that valuation.

Write \(c=a_0^2/4\). If \(c\) is not a fifth
power in \(\mathbb Q_5\), translate \(u\) by the
unique root of unity of order dividing four
with residue \(c\bmod5\). The resulting polynomial
has constant coefficient of valuation one, since
the fifth powers of principal units are exactly
\(1+25\mathbb Z_5\). The length-five segment has
slope denominator five and yields an irreducible
factor of degree five. If \(c=r^5\) in \(\mathbb Q_5\),
the derivative of \(u^5-c\) at each \(r\zeta_5^i\)
has valuation one. Hensel's lemma in
\(\mathbb Q_5(\zeta_5)\) lifts those distinct
approximations to five roots in disjoint balls
of radius at most \(5^{-2}\).
The rational root and the conjugate four give
degrees \(1,4\). Thus the local degree lists are
\(3,5\) and \(1,3,4\).
By HT4 there can be no rational linear factor.
The only candidate rational splittings are
\(3+5\) and \(4+4\).

The latter is impossible. At any odd prime
\(q\mid ab(a+b)\), put the divisible weight at
\(a\) and write \(k=v_q(a)\).
The same leading numerator \(4a^2p^5-a^4\)
at weight \(2k/5\) gives a local irreducible
factor of degree five whenever \(5\nmid k\).
Here \(q\) need only be odd so that \(4\) is a
unit. Two rational factors of degree four
would therefore force every odd prime
valuation of \(ab(a+b)\) to be divisible by
five. The pairwise coprimality of the three
weights makes the two odd weights fifth
powers and the even weight a fifth power
times a power of \(2\).
The relation between the weights now has the
form \(r^5+s^5=2^j w^5\) in coprime signed
nonzero integers, where \(0\le j\le4\).
For \(j=0\) Fermat's theorem excludes it.
For \(j=2,3,4\) use
[Ribet, Theorem 3](https://arxiv.org/pdf/math/9508208).
For \(j=1\),
[Darmon–Merel, Main Theorem (1)](https://perso.imj-prg.fr/loic-merel/wp-content/uploads/merel-pub/winding.pdf)
leaves only weights \(1,1,2\), up to order.
Those weights give an irreducible octic by
HT2. This eliminates \(4+4\).

## Finiteness of exceptional fifth-power tangent weights

**Corollary HT6.** Only finitely many positive
rational ratios \(a/b\) have reducible residual
octics. The possible exceptional ratios lift to
rational points on a fixed curve of genus \(71\).

Take \(v=a/b\) and set \(b=1\).
The degree-eight polynomial \(Q_{v,1}(t,1)\)
has leading coefficient \(15\) and discriminant
\[
 -2^{16}5^7 v^8(v+1)^8
 (v^2-v+1)^5(v^2+v+1)
 (v^2+3v+3)^5(3v^2+3v+1)^5.                            \tag{3}
\]
The identity follows by taking the discriminant
of (1). At \(v=2\), the reductions at
\(67,47,37\) are squarefree with degree lists
\((8),(1,7),(2,3,3)\), respectively.
Thus the Galois group of this specialization
is transitive and contains a \(7\)-cycle
and a transposition. The \(7\)-cycle excludes
nontrivial blocks in degree eight, so the
group is \(S_8\): in a primitive permutation
group the normal closure of a transposition
generates a connected graph and hence all
transpositions. The good specialization
embeds as a permutation group into the
generic arithmetic group. Thus that group is
also \(S_8\). The geometric group is a normal
subgroup and by (3) is not in \(A_8\),
so it is \(S_8\).

At \(v=0,-1,\infty\) inertia is a \(5\)-cycle.
For \(v=0\), the Newton initial equation in
the \(t=0\) cluster is \(4t^5=v^2\),
after removing a unit. The complementary
cubic \(15t^3+40t^2+20t+4\) at \(v=0\)
is simple and nonzero at \(0\).
Permuting the three signed weights treats
\(v=-1,\infty\) the same way.
The other eight branch values are the roots
of the four quadratic factors in (3).
Each has exactly one double root
and all other roots simple. For the three axis
quadratics this follows by looking at the
intersection of the Fermat quintic with the
diagonal conic at its one coordinate-axis point.
The conic vanishes to order two there, because
the other-coordinate ratio differs from its
axis value first in order five.
At \(v^2+v+1=0\) there is one additional
tangency at nonzero coordinates: with \(Z=1\),
the other ratios are the two distinct primitive
cube roots, as follows from the normal equations.
Its contact is exactly two by the same Hessian
calculation as in HT1.
The odd discriminant exponents imply
transposition inertia at all eight values.
There is no other branch value by (3).

Take the smooth projective cover corresponding
to unordered three-element subsets of the eight
roots. Its degree is \(56\) and it is geometrically
connected by the \(S_8\) action. A \(5\)-cycle
fixes one such subset and has eleven cycles
of length five. A transposition fixes
\(26=\binom63+\binom61\) such subsets and has
fifteen cycles of length two. Riemann–Hurwitz
therefore gives
\[
         2g-2=-112+3\cdot44+8\cdot15=140 .
\]
By HT5 a reducible positive specialization
defines a rational three-subset in the
unramified fibre over \(v\). Its lift is a
rational point on the subset cover. Faltings's
theorem gives finiteness, with no classification
of its rational points here.

The discriminant and specialization controls
in (3) can be replayed in PARI/GP from
~~~gp
x='x; v='v;
A=v^2*(x+1)^5; B=(x-v)^5; C=(v+1)^2*x^5;
Q=(4*A*B-(A+B-C)^2)/(v^2*(v+1)^2);
D=v^8*(v+1)^8*(v^2-v+1)^5*(v^2+v+1)*(v^2+3*v+3)^5*(3*v^2+3*v+1)^5;
print(poldisc(Q,x)/D);
f=subst(Q,v,2);
for(i=1,3,q=[67,47,37][i]; F=factormod(f,q); print(vector(matsize(F)[1],j,[poldegree(F[j,1]),F[j,2]])));
~~~

HT1 caps the polynomial square removal in
the fixed weighted linear model. HT5 and HT6
leave a finite and unclassified set of ratios
that could reduce the residual octic to two
factors of degree at most six. Finding a positive
ratio in that set with admissible prime and
squarefree conditions is a separate problem.
The alternative is a large-prime-square
estimate that applies to an irreducible
residual octic.

## The generic tangent cover at every higher odd prime

**Theorem HT7.** For every fixed odd prime \(\ell\ge5\),
the arithmetic and geometric generic Galois groups of
\(Q_{v,1}(t,1)\) over the normalized weight-ratio line
are \(S_{2\ell-2}\). Only finitely many positive rational
weight ratios have reducible residual forms. The
finiteness is for each fixed \(\ell\).

Put \(m=\ell-2\), \(n=2\ell-2\).
The three degenerate weights occur at \(v=0,-1,\infty\).
At zero the residual polynomial specializes to
\[
 t^\ell\left(4(t+1)^\ell-(2t+\ell)^2t^{\ell-2}\right).
\]
The bracket has degree \(\ell-2\), is nonzero at
\(t=0,-1,-\ell/2\), and has simple roots.
Indeed equating the logarithmic derivatives of its
two terms at a repeated root would give
\((2t-\ell+2)(2t+\ell)=4t(t+1)\), impossible
since \(\ell(\ell-2)\ne0\).
The Newton initial equation at the zero cluster
is \(4t^\ell=v^2\). It gives one cycle of length
\(\ell\), since \(\ell\) is odd. Permuting the
three signed weights gives the same inertia at
\(-1,\infty\).

At a nondegenerate weight, pull back the line to
the conic as in HT1. The axis intersections give
the three families
\[
          v=u-1,\qquad v=1/(u-1),\qquad v=-u,
              \quad u\in\mu_m\setminus\{1\}.
\]
At the first family, set \(Z=1\) and write the
two equations as
\[
 v X^2+Y^2=v+1,\qquad v X^\ell+Y^\ell=v+1 .
\]
At \(X=0,\ Y^2=u\), the determinant in \(Y,v\)
is nonzero, and \(v-v_0\) has order two as
a function of \(X\). For the latter assertion
subtract \(\ell/2\) times the conic equation from
the Fermat equation. Its derivative in \(v\)
at the axis point is \(\ell/2-1\), its derivative
in \(Y\) is zero there, and the \(X^2\) term
has nonzero coefficient \(-\ell v_0/2\).
This gives transposition inertia; the other
axes follow by permutation.

For a secondary tangency all coordinates are
nonzero. The normal equations in \(Z=1\)
force \(X^m=Y^m=1\).
Writing \(u=X^2,\ w=Y^2\), the nontrivial
tangent requires \(u,w\in\mu_m\setminus\{1\}\)
distinct, and
\[
                         v=-\frac{w-1}{u-1}.
\]
The Hessian calculation of HT1 gives quadratic
contact. The derivative in \(v\) of the Fermat
equation minus \(\ell/2\) times the conic equation
is \((1-\ell/2)(X^2-1)\ne0\).
Thus each secondary tangency gives a simple
branch of the residual cover.

All these branch values are distinct.
For \(u=\exp(2\pi i j/m), w=\exp(2\pi i k/m)\)
with \(1\le j,k<m\), the argument of
\((w-1)/(u-1)\) is \((k-j)\pi/m\) in
\((-\pi,\pi)\). For a fixed nonzero difference
\(k-j\), its absolute value is
\(\sin(\pi k/m)/\sin(\pi j/m)\), which is
strictly monotone in \(j\) by the sine-addition
formula. This distinguishes the secondary values.
They cannot meet \(u-1\), whose argument, after
division by \(\pi/m\), is a half-integer modulo
\(2m\), as \(m\) is odd. The reciprocal axis
family is excluded by the same argument.
A secondary value of absolute value one
has \(k+j=m\) and then belongs to \(\mu_m\),
which excludes the third axis family \(-\mu_m\).
The families \(u-1\) and \(-u\) cannot meet
because \(|u-1|=1\) would force an element
of order six. The reciprocal and third
families are disjoint by inversion.
If \((u-1)(w-1)=1\), their arguments
would give \(j+k=m\). The absolute value
condition would again require order six.
This distinguishes the first two axis families
as well. Tangencies and axis marks exhaust
possible repeated roots by the pullback
description. The branch profile is therefore
three \(\ell\)-cycles and
\((m-1)(m-2)+3(m-1)=(\ell-1)(\ell-3)\)
transpositions.

The generic arithmetic group is transitive
because it contains the group of the
irreducible specialization at \(v=1\) from HT2.
It contains a cycle of prime length
\(\ell>n/2\), so it is primitive. Indeed any
nontrivial invariant block system has both
block size and block count less than \(\ell\),
so an element of order \(\ell\) cannot act
nontrivially. A transposition in a primitive
group generates all transpositions by
conjugation, as in HT6. This gives the
arithmetic \(S_n\). The geometric group is
normal in \(S_n\) and contains a transposition,
so it too is \(S_n\).

For each \(1\le j\le\ell-1\), take the
geometrically connected cover with fibre
the \(j\)-subsets of the roots. It has degree
\(N_j=\binom nj\). A transposition contributes
\(\binom{n-2}{j-1}\) to the ramification sum.
An \(\ell\)-cycle fixes exactly
\(\binom{\ell-2}{j}\) subsets, counted as
zero if \(j>\ell-2\). All its other cycles
have length \(\ell\). Thus
\[
 2g_j-2=-2N_j
       +3(1-1/\ell)\left(N_j-\binom{\ell-2}{j}\right)
       +(\ell-1)(\ell-3)\binom{n-2}{j-1}.
\]
For \(\ell\ge7\) the final term is at least
\(2N_j\), using \(j(n-j)\ge n-1\) and
\((\ell-1)(\ell-3)/n\ge2\).
The middle term is positive. For \(\ell=5\)
the four genera are \(3,27,71,95\).
Thus every subset cover has genus above one.
A reducible positive rational specialization
has a rational factor of degree at most
\(\ell-1\), so supplies a rational point
on one of these finitely many covers.
Faltings's theorem proves the finiteness.
