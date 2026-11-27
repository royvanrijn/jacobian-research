# Prime progression witnesses at higher odd exponents

For every odd prime exponent, two fixed finite-field power tests
can enforce independence of consecutive coherent norm classes
in prime progressions. The two test directions exist because
the corresponding functions on a fixed hyperelliptic curve are
geometrically independent power classes. Quadratic twist
multiplicity then gives a count of distinct base fields.

## The hyperelliptic cover and its power classes

Fix an odd prime \(\ell\), and define
\[
 \begin{split}
 H_\ell(U,V,W)&=2\bigl(U^\ell V^\ell+U^\ell W^\ell+
                              V^\ell W^\ell\bigr)
                      -U^{2\ell}-V^{2\ell}-W^{2\ell},\\
 F_\ell(p,h)&=H_\ell(p-h,p,p+h),\\
 \mathcal C_\ell&:\quad y^2=-F_\ell(t,1).
 \end{split}
\]
Use \(A=(t-1)^\ell,\ B=t^\ell,\ C=(t+1)^\ell\) on the last curve
and put \(a_1=(A+B-C+y)/2,\ a_2=(B+C-A+y)/2\).

**Proposition HP1.** The binary form \(F_\ell\) is separable.
The functions \(a_1,a_2\) have geometrically independent
classes modulo \(\ell\)-th powers. Their simultaneous Kummer
cover of \(\mathcal C_\ell\) is geometrically connected and
unramified, of degree \(\ell^2\) and genus
\(1+\ell^2(\ell-2)\), in characteristic zero or all sufficiently
large characteristics.

The projective curve \(H_\ell=0\) is the image of the
Fermat curve \(X^\ell+Y^\ell+Z^\ell=0\) under coordinate
squaring. The map is injective on geometric points in
characteristic zero: up to global sign, two points with
the same squares differ by a sign of at most one coordinate,
and a nonzero odd power cannot change sign alone in the
equation. The image has degree \(2\ell\) since the pullback of a
line has degree \(2\ell\).
The progression line \(U-2V+W=0\) pulls back to the conic
\(X^2-2Y^2+Z^2=0\). At a point on both curves all coordinates
are nonzero. If \(Y=0\), \(Z/X\) would be \(\pm i\)
with odd power \(-1\), which is impossible. With
\(X=0\) or \(Z=0\), the modulus of the ratio of the other
coordinates would have to be both \(\sqrt2\) and \(1\).
If the gradients were dependent, putting \(m=\ell-2\) would
give
\[
       (Z/X)^m=1,\qquad (Y/X)^m=-2,\qquad
       1+(Z/X)^2=2(Y/X)^2 .
\]
The first equation makes the modulus of the left side of
the third at most two. The second makes its right side
have modulus \(2\cdot 2^{2/m}>2\). Thus the intersections
are transverse, giving \(2\ell\) distinct roots on the
progression line. Infinity is unramified since
\(F_\ell(1,0)=3\). This gives genus \(\ell-1\) for
\(\mathcal C_\ell\).

On \(\mathcal C_\ell\), conjugating \(y\) to \(-y\) gives
\(a_1\bar a_1=AB,\ a_2\bar a_2=BC\).
Both functions have pole order \(\ell\) at each of
the two infinities, since their outer terms there have
coefficients \(1\pm\sqrt{-3}\). At each root of their
norm they vanish on exactly one sheet with multiplicity
\(\ell\). For example at \(t=1\), \(B=1,\ C=2^\ell\)
and the norm of \(a_2\) is nonzero, while \(a_1\)
vanishes on one of the two sheets. The analogous endpoint
for \(a_2\) is \(t=-1\). At \(t=0\) the two functions
also have one-sided zeros, since \(A=-C\) there.
In particular their divisors are multiples of \(\ell\).

For a proposed relation take exponents represented in
\([-(\ell-1)/2,(\ell-1)/2]\) and replace negative powers
by positive powers of the conjugate, modulo \(\ell\)-th
powers. Their absolute values sum to \(m_0\le\ell-1\).
If the resulting function is an \(\ell\)-th power, its
root is regular away from infinity with pole orders at
most \(m_0\) there. It must be a polynomial in \(t\), since
on the affine curve it is of the form \(P(t)+yQ(t)\),
and a nonzero \(yQ(t)\) has pole order at least \(\ell\)
at one infinity even after cancelling with \(P(t)\).
For a nonzero relation one endpoint \(t=1\) or \(t=-1\)
would give the root a zero on one sheet and none on
the other. This contradicts being a polynomial in \(t\).
The same argument holds in all characteristics outside
the fixed bad primes, where the endpoint fibres and the
pole orders have the same description. This proves geometric independence.
The valuation divisibilities and characteristic prime to
\(\ell\) give an unramified Kummer cover. Riemann–Hurwitz
gives its genus.

## Two independent classes in actual prime progressions

**Theorem HP2.** Fix an odd prime \(\ell\). For each
quadratic sign there is a \(B_\ell>0\) such that at least
\(\gg_\ell X^{1/\ell}/(\log X)^{B_\ell}\) distinct
quadratic fields of that sign and discriminant of
absolute value at most \(X\) admit coherent prime
progression triangles at a fixed common scale
\(X^{1/(2\ell)}\). The coefficient on every edge, with the
triangle directed along the chosen order, is the same
positive integer before any sign adjustment for totally
positive generators. The two consecutive virtual units
are independent modulo \(\ell\)-th powers. In the
imaginary case their ideal classes are independent.
In the real case their ideal-class span is nonzero.
Except for the at most one quadratic field inside
\(\mathbb Q(\zeta_\ell)\), the same-colour Jacobian
\(J:y^2=x^\ell-F_\ell(p,h)/4\) has rank at least two.
The unreduced parameter \(|F_\ell(p,h)|\) is comparable
to \(X\); its squarefree part need not be.

Choose two distinct sufficiently large rational primes
\(q_i\equiv1\bmod \ell^2\) of good reduction for the
cover of HP1. For the first, prescribe Frobenius
\((1,0)\) in its Kummer group, and for the second
prescribe \((0,1)\), using a generator of \(\mu_\ell\)
at each prime. Such residue directions on
\(\mathcal C_\ell(\mathbb F_{q_i})\) exist by
Hasse–Weil on the twisted covers. Indeed the relevant
twist is geometrically connected of the same fixed
genus. The number of its rational points tends to
infinity with \(q_i\). Deleting the finitely many
points above zeros, poles or branch points of the
displayed functions and of \(t\) leaves a point with
the desired Frobenius. At the resulting direction
\((t_i,y_i)\), \(a_1,a_2\) are units and their
\((q_i-1)/\ell\)-th powers give the prescribed row of
power-residue tests. Exclude \(t=0,\pm1\) as well.

Take prime progressions \((p_1,p_2,p_3)=(p-h,p,p+h)\)
with \((p,h)\equiv(t_i,1)\bmod q_i\), \(p\) odd
and \(h\equiv0\bmod6\), \(3\nmid p\).
These form compatible congruences. Modulo any other
prime, setting \(p=1,h=0\) makes all three forms units.
Let \(E=F_\ell(p,h)\), let \(D\) be the signed squarefree
part of \(-E\), and write \(-E=D\delta^2\) with \(\delta>0\).
For the imaginary sign use a fixed positive rectangle
with \(h/p\) sufficiently small, so \(E>0\). For the
real sign choose \(h/p<1\) sufficiently large that
\((p+h)^{\ell/2}>p^{\ell/2}+(p-h)^{\ell/2}\), so \(E<0\).
In both cases the rectangle has \(p,h\asymp Z\) and
\(|E|\asymp_\ell Z^{2\ell}\).
For \(\ell\ge5\) the negative Heron condition already
holds near \(h/p=3/10\), so the three prime values can
also lie in the same dyadic interval for both signs.

Put \(A=p_1^\ell,\ B=p_2^\ell,\ C=p_3^\ell\).
For \(u_{12}=A+B-C,\ u_{23}=B+C-A,\ u_{13}=B-A-C\)
the signed trace identities give integral generators
\(\alpha_{ij}=(u_{ij}+\delta\sqrt D)/2\),
norms \((p_i p_j)^\ell\), and
\(\alpha_{12}\alpha_{23}=p_2^\ell\alpha_{13}\).
At \(p_1,p_3\) the Heron numerator is a unit once
the vertex primes exceed the prime factors of
\(2^\ell-1\); at \(p_2\) it is a unit because
\(p_1^\ell=-p_3^\ell\bmod p_2\). Here \(p_i\nmid h\)
by distinct primality and the chosen sizes.
The incident units and odd prime powers give the
\(\ell\)-th powers of the specified oriented prime
ideals, one-sided at every supporting prime. Thus
these are genuine coherent generators, in either sign.

At the selected auxiliary primes, \(q_i\nmid E\) and
\(q_i\) splits in \(\mathbb Q(\sqrt D)\) under
\(\sqrt D\mapsto y_i/\delta\bmod q_i\). Reduction of
the pair \(\alpha_{12},\alpha_{23}\) gives the two fixed
power tests. Any dependence of their virtual-unit
classes would be an \(\ell\)-th power, and its local
root is a unit because the product being tested is a
unit. The two rows exclude this relation. For an
imaginary quadratic field, every global unit has
trivial reduction modulo \(\ell\)-th powers at those
auxiliary primes, including the fields with extra
roots of unity. Thus the same rows also exclude a
dependence of the ideal classes. In the real case
the kernel of the map to the ideal class group has
dimension one. The evaluation in CR1 gives the
Jacobian rank from the virtual-unit independence,
with the sole possible cyclotomic quadratic exception.

By the [Green–Tao prime forms theorem](https://arxiv.org/pdf/math/0606088)
there are \(\gg_\ell Z^2/(\log Z)^3\) such prime pairs.
The forms have nonparallel linear parts and the
fixed residue choices are locally admissible.
After fixing a signed squarefree core \(D\), the
primitive pairs \((p,h)\) give rational points on
\(D z^2=-F_\ell(t,1)\) with \(t=p/h\); these ratios
are distinct by \(\gcd(p,h)=1\).
The [PP2 multiplicity argument](ABC_COHERENT_PRIME_PROGRESSION_NORMS.md#a-uniform-bound-on-rational-progressions-per-field)
works in genus \(\ell-1\): the 2-torsion is rational
over the fixed branch splitting field, the rank is
\(O_\ell(1+\omega(D))\) by the unramified Kummer bound,
and uniform Mordell–Lang bounds the number of points
by \(C_\ell^{1+\omega(D)}\). Similarly the
[PP3 prime-factor cutoff](ABC_COHERENT_PRIME_PROGRESSION_NORMS.md#quantitative-count-of-genuine-fields)
applies in degree \(2\ell\): replace the root bound
six by \(2\ell\) and choose the constant in
\(k=\lceil C_\ell\log\log Z\rceil\) large enough.
Outside \(O_\ell(Z^2/(\log Z)^4)\) unweighted pairs,
\(\omega(D)=O_\ell(\log\log Z)\).
We lose at most a fixed power of \(\log Z\) for
multiplicities. Taking \(Z\) a sufficiently small
constant times \(X^{1/(2\ell)}\) accommodates the
fundamental discriminant bound \(4|D|\le X\).

For \(\ell=5\), the progression \((41,47,53)\) gives
\[
 \begin{split}
 E&=100955999274692403
    =3\cdot2111\cdot25841\cdot616897951,\\
 (u_{12},u_{23})&=(-72994285,531684299).
 \end{split}
\]
The factors of \(E\) are prime. At \(11\) with
\(\sqrt{-E}=4\bmod11\), the fifth-power residue values
are \((5,4)\) in the basis generated by \(4\).
At \(101\) with \(\sqrt{-E}=49\bmod101\) they are
\((36,36)\) in the basis generated by \(95\).
The exponent rows are \((2,1),(2,2)\), hence independent.
These exact tests give explicit fixed congruences for
the fifth-power imaginary family. Conditions modulo \(\ell\)
suffice for the auxiliary primes when \(\ell\ge5\), since
imaginary quadratic units then have order prime to \(\ell\). The broader lower
count of [Bartz–Levin–Thamminana, Theorem 1.3](https://arxiv.org/pdf/2502.00845)
has order \(X^{1/3}/(\log X)^2\) for imaginary fields
of 5-rank at least two. The count here also prescribes
the shared colour and coherent prime progression.

## A second irreducibility obstruction in the unit slice

**Corollary HP3.** Let \(\ell\ge5\). Suppose the Jacobian
of the degree-\(\ell\) Fermat curve has injective
reduction on its rational Mordell–Weil group at \(2\).
On any homogeneous rational line \(aU+bV+cW=0\)
with primitive integral coefficients, \(abc\ne0\),
\(a+b+c\) even, and no opposite pair, the restriction
of \(H_\ell\) is irreducible of degree \(2\ell\).
This applies unconditionally for \(\ell=5\).

As in SW1, the conic \(aX^2+bY^2+cZ^2=0\) meets the
Fermat curve only above one flex \(O\) modulo \(2\).
The reduction condition implies that any degree-\(n\)
closed point in the intersection gives a rational
divisor equivalent to \(nO\). Interchange coordinates
so that \(O=[1:-1:0]\). The functions
\(Z/(X+Y)\), \(X/(X+Y)\) have their only poles at \(O\)
with orders \(\ell-1,\ell\) respectively.
The semigroup generated by these two orders already
has \((\ell-1)(\ell-2)/2\) gaps, the genus of the
curve, so it is the full pole semigroup there.
Consequently a degree-\(n\) closed point of the conic
section with \(n\le\ell\) requires \(n=\ell-1\) or
\(n=\ell\). In the former case its divisor is a
line through \(O\) with the fixed point \(O\) removed,
by the basis \(1,Z/(X+Y)\). In the latter case it is
a line section, since \(nO\) is the flex tangent class.
Either line would contain at least four distinct
points of the nondegenerate conic, which is impossible.
No point of degree at most \(\ell\) occurs, so the
intersection is one closed point of degree \(2\ell\).
The square map preserves that degree by injectivity
on geometric points. For \(\ell=5\),
[Klassen–Tzermias, Theorem 2](https://www.impan.pl/shop/publication/transaction/download/product/109631)
give \(J(\mathbb Q)\simeq(\mathbb Z/5\mathbb Z)^2\),
whose reduction is injective at the good prime \(2\).

Thus linear slices in the unit-normalized quintic
case do not reduce the field-square problem to the
binary forms of factor degree at most six. The
progression result HP2 bypasses that obstacle for
a cumulative discriminant count. A dyadic
squarefree-field result with all three prime marks
would need control of the square part or a way to
change the profile and reduce the degrees.

**Corollary HP4.** In degree five the irreducibility conclusion
also holds on any rational line with no zero coefficient or
opposite pair, provided the line has a positive real direction.
The coefficient parity assumption is unnecessary in this case.

Use the classification of points of degree at most five on
the Fermat quintic from
[Klassen–Tzermias, Theorem 1](https://www.impan.pl/shop/publication/transaction/download/product/109631).
The three rational flexes are excluded by the assumption on
opposite pairs. Its only quadratic pair has squared images
\([1:\omega:\omega^2]\) and \([1:\omega^2:\omega]\),
up to permutation, where \(\omega\) is a primitive cube
root of unity. Their line has all coefficients equal,
so it has no positive direction. There are no cubic points.
For degree four or five the conjugates lie on a rational
line by that theorem. They cannot also lie on a
nondegenerate conic, which has at most two points on
any line. Hence the intersection of the Fermat quintic
with the diagonal conic has no closed points of degree
at most five. It is one closed point of degree ten,
as is its squared image. In particular fixed scalings
\(P_i=b_i^2p_i\) in the common fifth-power colour cannot
produce a reducible decic from a genuine positive prime
line. A line through a rational cusp violates the
incident prime-unit test as in SW1.

## Removing a tangent square from the field parameter

**Theorem HP5.** For each fixed odd prime \(\ell\ge5\),
there are \(\gg_\ell X^{1/(\ell-1)}/(\log X)^{B_\ell'}\)
quadratic fields of either sign and discriminant at most
\(X\) in absolute value with genuine prime progression
triangles whose normalized doubled-trace coefficients are
\((1,1,2)\) up to signs. The consecutive virtual units are independent;
their ordinary ideal classes are independent in the imaginary
case and span a nonzero space in the real case. They give
rank at least two on the corresponding common-colour
Jacobian, outside the cyclotomic quadratic exception.

Replace the middle value \(B=p^\ell\) of HP1 by
\(B=4p^\ell\), and keep \(A=(p-h)^\ell,\ C=(p+h)^\ell\).
Define the binary form
\[
              R_\ell(p,h)=\frac{4AC-(A+C-B)^2}{4h^2}.
\]
This is integral of degree \(2\ell-2\): the even
polynomial in \(h\) upstairs vanishes at \(h=0\)
and its coefficients are multiples of four.
Its leading coefficient in \(p\) is \(\ell(\ell-2)\),
by expanding to second order in \(h\).
The common discriminant parameter is \(E=h^2 R_\ell\),
so the squarefree core comes from the smaller-degree form.
The geometric source of the square is tangency at
\(h=0\). The new branch curve is the square image of
\(X^\ell+2Y^\ell+Z^\ell=0\), cut by
\(X^2-2Y^2+Z^2=0\). There is no intersection on a
coordinate axis: when \(Y=0\) this is as in HP1,
and otherwise an axis point would require a ratio
with modulus both \(\sqrt2\) and \(2^{1/\ell}\).
At a dependent-gradient intersection one has
\((Z/X)^{\ell-2}=1,\ (Y/X)^{\ell-2}=-1\), so
\(|Z/X|=|Y/X|=1\). The conic forces
\(|1+(Z/X)^2|=2\), whence \((Z/X)^2=1\).
As \(\ell-2\) is odd, this gives
\([X:Y:Z]=[1:-1:1]\), the point \(h=0\).
Its line contact is exactly two by the leading
coefficient calculation. Thus \(R_\ell\) is separable.

In the variable \(t=p/h\), on \(y^2=-R_\ell(t,1)\),
the first two traces are
\((A+B-C)/2,(B+C-A)/2\) with \(h=1\).
Define the two norm functions from them as before.
Their zeros above the distinct endpoint fibres
\(t=1,-1\) are still one-sided, their norms are still
\(AB/4=(t-1)^\ell t^\ell\) and
\(BC/4=t^\ell(t+1)^\ell\), and each has pole order
\(\ell\) at the two infinities. Here \(y\) has pole
order \(\ell-1\); the trace has order \(\ell\) with
nonzero leading coefficient. At the mid fibre \(t=0\),
\(y^2=1\) and the zeros on the two sheets again have
multiplicity \(\ell\).

For a nonzero proposed relation choose centered
exponents as in HP1. If their absolute values sum
to \(\ell-1\), both have absolute value
\((\ell-1)/2\). Multiplying the exponents by the inverse of
\((\ell-1)/2\) modulo \(\ell\) reduces the two absolute
values to \(1,1\). In all cases we thus have a
representative whose absolute exponent sum is
strictly less than \(\ell-1\), using \(\ell\ge5\).
A root of the product would be a polynomial in \(t\),
since a term involving \(y\) has pole order at least
\(\ell-1\) at one of the infinities. The one-sided
endpoint zeros exclude this. Therefore the unramified
Kummer cover works exactly as in HP1 and HP2, now
on the genus-\(\ell-2\) twist curve.

For actual prime progressions take
\(A=p_1^\ell,\ B=4p_2^\ell,\ C=p_3^\ell\).
The trace numerators below are even because
\(A,C\) are odd and \(B\) is even. The incident-prime units hold for large
primes exactly as before, now using the fixed nonzero
difference \(4-2^\ell\) at the outer primes.
The quantity \(E\) is \(h^2R_\ell(p,h)\); choose the
auxiliary congruences with \(h=1\) at both test primes,
so the function residues are those of the indicated
hyperelliptic cover. The consecutive generators are
\(\alpha_{12}=(u_{12}+\delta\sqrt D)/2\) and
\(\alpha_{23}=(u_{23}+\delta\sqrt D)/2\), with traces
\(u_{12}=(A+B-C)/2\) and \(u_{23}=(B+C-A)/2\).
Each generator is integral by its monic trace–norm
polynomial. The third generator is
\(\alpha_{13}=(B-A-C)/2+\delta\sqrt D\), since direct
multiplication gives
\(\alpha_{12}\alpha_{23}=(B/4)\alpha_{13}\).
Their norms are \((p_1p_2)^\ell,(p_2p_3)^\ell,
(p_1p_3)^\ell\). In the doubled trace convention the
coefficients are \((\delta,\delta,2\delta)\) if the
unreduced parameter is \(E\). This gives the
two ideal and rank assertions by the same argument as
HP2. For the positive Heron sign take \(h/p\) small
by the leading coefficient above; for the negative
Heron sign take \(h/p<1\) near one, since
\(2^\ell>4\). Then \(|R_\ell|\asymp Z^{2\ell-2}\) in
both rectangles. The same prime-form density, twist
multiplicity and prime-factor truncation yield the
improved count on taking \(Z\) a sufficiently small
fixed constant times \(X^{1/(2\ell-2)}\). This argument
does not require the binary form \(R_\ell\) itself to
take squarefree values at the prime tuples.

For instance with \(\ell=5\) and \((p,h)=(47,6)\)
the squarefree part is
\(E/h^2=325427449513839=3\cdot37\cdot177379\cdot16528331\).
The consecutive traces are \(307520368,609859660\).
At \(11\) with square root of \(-E\) equal to \(4\),
the fifth-power residue values are \((9,4)\) in the
basis generated by \(4\). At \(31\) with square root
of \(-E\) equal to \(3\) they are \((16,2)\) in the
basis generated by \(16\), giving exponent rows
\((3,1),(1,4)\) modulo five. The determinant is one.
