# Primitive power and powerful values of binary forms

Fix \(d\ge5\) and a prime \(p\). Let \(F\) be a
squarefree homogeneous binary form of degree \(d\) with
integer coefficients, and let \(n\ne0\) be an integer.
After a bounded unimodular change if necessary, write
\(a=F(1,0)\ne0\), set
\(R=p a n\operatorname{disc}(F(T,1))\), and let \(K\)
be the splitting field of \(F(T,1)\) with \(\zeta_p\)
adjoined. The number of primitive integer pairs
\((s,t)\) for which \(F(s,t)=n z^p\ne0\) is at most
\[
        \exp O_{d,p}\!\left(
                   1+\omega(|R|)+\log|\mathrm{Cl}(K)[p]|\right).
 \tag{BT1}
\]
This is independent of the heights of \(s,t,z\).
In particular it is \(H^{o(1)}\) for forms of coefficient
height \(H^{o(1)}\) and any polynomial-height \(n\).
The prime \(p\) need not divide \(d\). The argument first
recovers the power classes of the linear factors in the
splitting field. Balanced ratios of factors then supply
points on a cover of the projective parameter.

## An actual anchor for the norm factorization

Use a bounded integral unimodular substitution if needed
to give \(a=F(1,0)\ne0\). For example a shear with one
of \(d+1\) distinct integer slopes has nonzero evaluation
at its first column. This preserves primitivity and changes
the height cap only by a bounded factor.
Let \(K\) split \(F(T,1)\) with \(\zeta_p\) adjoined.
For its roots \(\theta_i\), set \(\alpha_i=a(s-\theta_i t)\).
These are integral since the \(a\theta_i\) satisfy the
monic polynomial \(a^{d-1}F(T/a,1)\). They satisfy
\[
                \prod_i\alpha_i=a^{d-1}n z^p .
\]
At a prime ideal of \(K\) over \(q\) not dividing
\(an\operatorname{disc}(F(T,1))\), at most one
\(\alpha_i\) is divisible by the ideal. If \(q\mid t\)
none is divisible by primitivity. Otherwise the
differences \(a(\theta_i-\theta_j)t\) are units.
Thus the product formula implies that the valuation of
each \(\alpha_i\) there is divisible by \(p\).

For every \(i\), write
\((\alpha_i)=\mathfrak a_i\mathfrak b_i^p\) with the
exponents of \(\mathfrak a_i\) in \([0,p-1]\).
The number of choices for the joint residual ideals
\(\mathfrak a_i\) is at most \(\exp O_{d,p}(1+\omega(|R|))\).
Each \([\mathfrak b_i]\) solves
\([\mathfrak b_i]^p=[\mathfrak a_i]^{-1}\), so there
are at most \(|\mathrm{Cl}(K)[p]|^d\) joint class
choices. On each occupied record choose an actual solution.
The quotient of each \(\mathfrak b_i\) by its reference
ideal is principal, so \(\alpha_i/\alpha_{*,i}\)
is a \(p\)-th power up to a unit. The classes of these
units modulo \(p\)-th powers cost a bounded number of
records by Dirichlet's unit theorem. On each occupied
refined record choose an actual new reference solution.
Then \(\alpha_i/\alpha_{0i}=\gamma_i^p\) in \(K^\times\)
for every \(i\) on that record. All the factors are
nonzero by \(F(s,t)\ne0\).

## Projective covers from conjugate ratios

Index the roots \(\theta_i\) in \(K\). For any vector
\((e_i)\) of integers with \(\sum_i e_i=0\), the products
\[
  \prod_i\left(\frac{s-\theta_i t}{s_0-\theta_i t_0}\right)^{e_i}
                 =\left(\prod_i\gamma_i^{e_i}\right)^p
\]
are \(p\)-th powers in \(K\). For \(t\ne0\), putting
\(U=s/t\) cancels
the \(t\) powers. The resulting function of \(U\) is a
nonzero constant times a ratio of linear factors, and
its constant is fixed by the anchor.

Use \(e_i=1\) on two roots and \(-1\) on two others
when \(p\ge3\). The resulting cyclic degree-\(p\) cover
has four branch points and genus \(p-1\ge2\). When
\(p=2\) and \(d\ge6\), use three roots on each side.
This cover has six branch points and genus two.

The rank bound for these cyclic covers is the same as for
a split polynomial. More explicitly, if the branch orders
are \(e_i\in\{1,-1\}\) with sum zero, the branch-divisor
relations in the proof of (FP2) now have residues congruent
to \(h e_i\) for a common \(h\bmod p\). Conversely the
corresponding eigenfunction
\(V^h\prod(U-\theta_i)^{(a_i-h e_i)/p}\) realizes the
relation and has order zero at infinity. This proves
rationality of the entire \((1-\sigma)\)-kernel
of order \(p^{b-2}\) for \(b\) branch points.
Outside the places above \(R\), the branch sections are
distinct and the anchor factors have valuations divisible
by \(p\), just as for any solution of the norm equation.
Thus the constant in the ratio has valuation divisible by
\(p\), and scaling \(V\) locally by a uniformizer power
makes it a unit. The resulting cover has a smooth proper
model. Use \(V^p=f(U)\) at the simple zeros, \(V^{-p}=f(U)^{-1}\)
at the simple poles, and a finite étale cover over the
complementary charts, including infinity.
The unramified Kummer argument of (FP2) bounds the rank by
\(O_{d,p}(1+\omega(|R|)+\log|\mathrm{Cl}(K)[p]|)\).
Uniform Mordell–Lang bounds the points exponentially in
this rank.
[Split cyclic-cover descent and uniform Mordell–Lang](ABC_PRIME_POWER_FERMAT_DESCENT.md#uniform-rational-point-bound-from-the-splitting-field)

For \(d=5,\ p=2\), use two square roots of ratios, with
branch sets \(\{1,2,3,4\}\) and \(\{1,2,3,5\}\).
Use two numerator and two denominator roots in each.
The independent square classes give a geometrically
connected Klein-four cover. Its union of branch sets
has five points with index two at each, so it has genus
two. The three quadratic quotients have branch counts
\(4,4,2\), respectively, and thus genera \(1,1,0\).
The pullbacks from the two genus-one quotients span
the two different nontrivial character spaces of
differentials on the cover. Their Jacobians therefore
give an isogeny from their product onto the cover's
Jacobian over \(K\). Each quotient has full rational
two-torsion by the split branch argument just given,
and the same Kummer rank bound holds, since the
rank argument allows genus one. The genus-two cover
has rank at most the sum of the two ranks. One can
choose a \(K\)-point as basepoint on the cover as soon
as any solution in the given record lifts to it. Uniform
Mordell–Lang gives the same exponential point bound. The case \(t=0\) has only \((s,t)=(\pm1,0)\).
For \(t\ne0\) each point determines the rational
ratio \(s/t\) and hence the primitive pair up to sign.
Summing over the anchor records proves (BT1).

If the coefficients of \(F\) have height \(H^{o(1)}\),
then \(K\) has discriminant \(H^{o(1)}\). The scaled
roots and \(\zeta_p\) are algebraic integers with all
conjugates of size \(H^{o(1)}\). A linear combination
with bounded integral coefficients can separate the
field embeddings: each equality of two embeddings
rules out a proper hyperplane, and a fixed grid larger
than the number of such pairs avoids their union.
The resulting primitive integer has order discriminant
\(H^{o(1)}\), bounding the field discriminant.
Minkowski's class number bound now gives
\(|\mathrm{Cl}(K)[p]|=H^{o(1)}\). This proves the subpower
specialization of (BT1). More generally, with coefficients
of polynomial height just the subpower \(p\)-class bound
suffices.

## Square-full and cube-full values

Take a squarefree form of fixed degree \(d\ge5\) with
coefficient height \(H^{o(1)}\), and let \(\mathcal S\)
be any set of primes whose product is \(H^{o(1)}\).
Consider primitive \(s,t\) with
\(\max(|s|,|t|)\le H\) and \(F(s,t)\ne0\). Require
every positive prime exponent of \(F(s,t)\) away
from \(\mathcal S\) to be at least two.
The number of pairs with this property is at most
\[
                             H^{d/5+o(1)} .
 \tag{BT2}
\]
For \(p=2,3\), let \(R_p\) be the product of primes
outside \(\mathcal S\) with \(p\nmid v_q(F(s,t))\).
For every \(e\ge2\),
\(3\,\mathbf 1_{2\nmid e}+2\,\mathbf 1_{3\nmid e}\le e\).
Therefore \(R_2^3R_3^2\le |F(s,t)|\), and at least one
\(R_p\) is at most \(H^{d/5+o(1)}\). After choosing \(p\)
and \(R_p\), the \(p\)-free coefficient including its
sign has only \(H^{o(1)}\) possibilities: each exponent
is in a fixed finite range and the extra support is in
\(\mathcal S\). Apply (BT1) on each coefficient record.

When every positive exponent outside \(\mathcal S\)
is at least three, the upper bound improves to \(H^{d/6+o(1)}\).
In fact \(2\,\mathbf 1_{2\nmid e}+3\,\mathbf 1_{3\nmid e}
+\mathbf 1_{5\nmid e}\le e\) for \(e\ge3\).
The same argument chooses among \(p=2,3,5\).
These bounds save a power in the primitive parameter
count for degrees \(5\) through \(9\) in the square-full
case, and through \(11\) in the cube-full case.


Darmon–Granville give qualitative finiteness for the
relevant fixed generalized power equations. Swaminathan's
work on odd-degree forms studies the frequency of
nonexistence as the form varies at fixed leading
coefficient. The count (BT1) instead bounds the fibre of
an individual form by the indicated support and class
group data. The square-full polynomial estimates of
Wongcharoenbhorn–Meemark concern quadratic polynomials
in one variable.
[Darmon–Granville](https://www.math.mcgill.ca/darmon/pub/Articles/Research/06.Granville/paper.pdf)
· [Swaminathan](https://doi.org/10.1112/S0010437X23007649)
· [Wongcharoenbhorn–Meemark](https://arxiv.org/abs/2405.06968v2)

## Near-cusp count for a moving square–cube difference

Chan–Koymans count all integer pairs in a dyadic cubic-root
range \(Z\le z<2Z\) with \(y'\ge0\) and
\(0<|z^3-(y')^2|\le1728 Y\). For \(Y\ge1,\ Z\ge Y^{2/3}\)
their bound is \(O_\varepsilon(Z^\varepsilon
(Y^{2/3}+Y^{1/6}Z^{1/4}))\). In particular, for positive
integer coefficients \(A,B\), \(K\ge1\), and \(H\ge1\),
consider
\[
                H\le A x^3<2H,\qquad
                       0<|A x^3-B y^2|\le K .
\]
If \(A^2B^3 K^2\le H\) and \(x,y>0\), the pair count is
\[
 \ll_\varepsilon (A^{2/3}B H^{1/3})^\varepsilon
       \left(A^{4/3}B^2K^{2/3}
                 +A^{1/2}B^{3/4}K^{1/6}H^{1/12}\right).
\]
Indeed \((z,y')=(ABx,AB^2y)\) is an injective integral
embedding, and the difference is multiplied by \(A^2B^3\).
Use \(Z=A^{2/3}BH^{1/3}\) and
\(Y=\max(1,A^2B^3K/1728)\). Then \(Z\ge Y^{2/3}\)
under the displayed condition, and the full \(z\)-interval
has ratio less than two. This transfer counts varying
differences and uses no primitivity hypothesis.

For fixed \(A,B\) and \(K=H^\tau\) with \(0\le\tau<1/4\),
this improves the ordinary ambient cusp exponent \(1/6\).
The smallest exponent provided by this estimate is \(1/12\)
in the complete-term height scale. A bounded difference
radical is a different restriction. A use at the small
\(abc\) budget would need to impose that restriction in
addition to numerical closeness, or to have root heights at
a much smaller scale. In the quartic proof the binary quartic
invariant retains the actual difference. At that stage one
can formulate its additional radical constraint, which the
cited point count does not use.
[Chan–Koymans, Theorem 1.2 and Section 2](https://arxiv.org/html/2611.31390v1#S1)
