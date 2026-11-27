# Excluding longitudinal degree eighteen at an exceptional polar fold

**Theorem (ExceptionalFoldMinimalWeightBoundaryExclusion).** Let a complex
polynomial plane Keller map have an exceptional polar fold in its full finite
normalization. Write \(a,b\) for the normalized polynomial source coordinates
in which the whole fold chart is \(a=b^{-2}+bt\). Then
\[
\max(\deg_a P,\deg_a Q)\ge 24.
\]
The bound uses the longitudinal \(a\)-degree and allows arbitrary \(b\)-degrees.

The [weighted Jacobian obstruction](EXCEPTIONAL_FOLD_WEIGHTED_JACOBIAN_OBSTRUCTION.md)
gives the bound eighteen and determines the two lowest blocks at equality.
We study the first departure from a square after the degree-eighteen
normalization. A finite Laurent expansion excludes the late departures.
A second choice of direction excludes the early ones and then rules out
degrees nineteen through twenty-three.

## Normalization and the forced Laurent coefficients

Assume the maximum is eighteen. Put \(u=ab^2=1+b^3t\) and
\(h=u^5(u-1)\). By the equality statement of the weighted obstruction,
an affine target change gives
\[
 P=b^{-6}R,\quad Q=b^{-9}S,\qquad R=h^2+O(b),\quad S=h^3+O(b).
\]
Here \(R,S\in\mathbb C[u,b]\). Denote their coefficients by \(R_i,S_i\),
and let \(J=\operatorname{Jac}_{a,b}(P,Q)\ne0\), using this target
normalization. The [simultaneous admission rule](EXCEPTIONAL_FOLD_POLYNOMIAL_ADMISSION_AND_MARKED_TRACE.md)
gives
\[
 \operatorname{ord}_{u=0} R_i\ge \max(0,\lceil(6-i)/2\rceil).
 \tag{1}
\]
We use order \(+\infty\) for a zero coefficient.

In fact \(\deg_a P=12\). To see this, write \(p=\deg_a P\) and let
\(A(b),B(b)\) be the leading \(a\)-coefficients of \(P,Q\). The leading
blocks imply \(p\ge12\), \(\deg_a Q=18\), and \(\operatorname{ord}_b B=27\).
The \(a^{p+17}\) coefficient of the Jacobian is
\(pAB'-18A'B=0\), giving
\(\operatorname{ord}_b A=27p/18\).
If \(p>12\), the \(u^p\) terms of \(P\) have weights at least \(-5\),
since its weight \(-6\) block has \(u\)-degree twelve.
Then \(\operatorname{ord}_b A\ge 2p-5>27p/18\).
Thus \(R\) has \(u\)-degree twelve.

The bracket in \(\mathbb C(u)((b))\) is
\(b^2\operatorname{Jac}_{u,b}\). Cancel the lowest commuting terms of \(Q\)
as in the weighted obstruction. At any remaining weight \(j<5\), its
coefficient \(c(u)\) satisfies \(c^6(h^2)^j\in\mathbb C^*\). Since \(h\)
has a simple root at \(1\), rationality requires \(3\mid j\).
The cancelling term is a constant times \(P^{-j/6}\), with
\(P^{1/2}=b^{-3}h+O(b^{-2})\). Thus the possible weights through four
are \(-9,-6,-3,0,3\). The first gives \(P^{3/2}\), and we remove
the terms at \(-6,0\) by subtracting a multiple of \(P\) and a constant
from \(Q\). At weight five the residual coefficient is \(f/h^2\), where
\[
 6f'-\frac{(h^2)'}{h^2}f=J,\qquad f=\frac{J}{4}u(u-1).
\]
The displayed solution is unique among rational functions since a
nonzero homogeneous solution would have a fractional order at \(1\).
After multiplying by \(b^9\), we obtain
\[
 S\equiv R^{3/2}+d b^6 R^{1/2}+e b^{12}R^{-1/2}
                  + b^{14}f/h^2\pmod {b^{15}}.
\]
Translate \(P\) by \(2d/3\), so \(R\) increases by \((2d/3)b^6\).
Binomial expansion absorbs the term at \(b^6\) and replaces \(e\)
by \(E=e-d^2/6\). In the resulting normalization,
\[
 S\equiv R^{3/2}+E b^{12}R^{-1/2}+b^{14}f/h^2
 \pmod {b^{15}}.
 \tag{2}
\]
These target changes preserve (1) and the two lowest blocks.

## The first departure from a square

Expand in \(u^{-1}\) at infinity and let \(L\) be the polynomial part of
\(R^{1/2}\), with the branch \(L_0=h\). These operations are well defined
over \(\mathbb C[[b]]\), since the leading \(u\)-coefficient of \(R\)
has constant term \(1\). Write
\[
 M=R-L^2,\qquad \deg_u L=6,\qquad \deg_u M\le5,\qquad M_0=0.
 \tag{3}
\]
Let \(n\) be the \(b\)-order of \(M\), using \(n=+\infty\) for \(M=0\).

Set \(G=S-L^3-\tfrac32 LM\), which has polynomial \(u\)-coefficients.
Binomial expansion of (2) starts with
\[
 G\equiv \frac38\,\frac{M^2}{L}
          -\frac1{16}\,\frac{M^3}{L^3}
          +\sum_{k\ge4}\binom{3/2}{k}\frac{M^k}{L^{2k-3}}
          + E b^{12}\bigl(L^{-1}+O(b^n)\bigr)
          + b^{14} f/h^2 \pmod {b^{15}} .
 \tag{4}
\]
For \(1\le n\le4\), the coefficient at \(b^{2n}\) gives \(h\mid M_n^2\).
Multiplying by \(L^3\) and reducing modulo \(L^2\) through order \(3n\)
gives \(M^3\equiv 0\bmod(L^2,b^{3n+1})\).
Indeed the quadratic term and any contribution of \(E b^{12}/L\) are
multiples of \(L^2\), and all subsequent terms have order greater than
\(3n\). Polynomial division by \(L^2\) is valid because its leading
coefficient is a unit. Its quotient in this congruence starts in
order at least \(3n\), so \(h^2\mid M_n^3\).
Since \(h=u^5(u-1)\) and \(\deg_u M_n\le5\), this forces
\[
 M_n=\gamma u^4(u-1),\qquad \gamma\ne0.
 \tag{5}
\]

If \(n\ge5\), the only terms of (4) through order fourteen, apart
from \(b^{14} f/h^2\), have a single denominator \(L\).
Multiplying by \(L\) gives
\[
 LG\equiv \tfrac38 M^2+E b^{12}+b^{14}f/h\pmod{b^{15}}.
\]
The first two terms have polynomial coefficients and \(f/h\) has
a pole at \(u=0\), a contradiction. This includes \(M=0\).

Suppose \(n=4\), and put \(N=b^{-4}M\).
Through order eleven, (4) after multiplication by \(L\) gives
\(N^2=LD\bmod b^4\) with polynomial \(u\)-coefficients in \(D\).
Through order fourteen, multiplication by \(L^3\) and reduction modulo
\(L^2\) gives \(M^3\equiv 16 b^{14}fL\bmod(L^2,b^{15})\).
Dividing by \(b^{12}\) gives
\[
 N^2=LD\bmod b^4,\qquad N^3=L^2C+16 b^2 fL\bmod b^3.
 \tag{6}
\]
Here the existence of polynomial coefficient quotients follows as
above from the unit leading coefficient of \(L\) and from taking the
first nonzero \(b\)-order. We have \(D_0=N_0^2/h\) and \(C_0=N_0^3/h^2\).
At \(u=0\) let \(H\) and \(m\) denote the nonzero leading coefficients
of \(h\) and \(N_0\), of orders five and four. Then \(D_0,C_0\)
have orders three and two and leading coefficients \(m^2/H,m^3/H^2\).

Write \(s=\operatorname{ord}_{u=0} L_1\). Since
\(R_2=2hL_2+L_1^2\) in this case, (1) gives \(s\ge1\).
The cubic equation at \(b^1\) reads
\[
 3N_0^2N_1=h^2C_1+2hL_1C_0.
 \tag{7}
\]
If \(s\ge3\), this forces \(\operatorname{ord}_0 N_1\ge2\).
At \(b^2\) in the cubic equation, all terms except \(16fh\)
then have \(u\)-order at least seven. In detail, the terms are
\(3N_0^2N_2,3N_0N_1^2,h^2 C_2,2hL_1C_1,L_1^2 C_0,2hL_2C_0\).
But \(fh\) has order six.

If \(s=1\) or \(2\), let \(l\) be the leading coefficient of \(L_1\).
By (7), \(N_1\) has order \(s-1\) and leading coefficient \(2ml/(3H)\).
For \(s=1\), the \(u^4\) coefficient at \(b^2\) in the cubic equation
is contributed only by \(3N_0N_1^2-L_1^2 C_0\). It equals
\(m^3l^2/(3H^2)\ne0\).
For \(s=2\), the \(b^1\) square equation
\(hD_1=2N_0N_1-L_1D_0\) gives \(D_1(0)=m^2l/(3H^2)\).
The \(u^2\) coefficient at \(b^2\) in the square equation is
contributed only by \(N_1^2-L_1D_1\).
It equals \(m^2l^2/(9H^2)\ne0\). Hence \(n=4\) is impossible.

## The second direction

We use the following directional form of the rational Jacobian
equation. Let \(p,q\) be coprime integers with \(q>0\) and \(p<q\).
Choose \(u_0\in\mathbb C\) and put \(u-u_0=s^p x,\ b=s^q\).
The substitution acts in the Laurent field and allows negative \(p\).
In these variables the bracket is
\[
 b^2\operatorname{Jac}_{u,b}=\frac{1}{q}s^{q-p+1}\operatorname{Jac}_{x,s}.
 \tag{8}
\]
Suppose the lowest \(s\)-term of \(P\) is \(s^{-\mu}\pi(x)\) with
\(\mu>0\) and \(\pi\) polynomial. Then there is a polynomial
\(F(x)\) satisfying
\[
 \mu F'-(q-p)\frac{\pi'}{\pi}F=qJ.
 \tag{9}
\]
To prove this, eliminate the commuting terms of \(Q\) in
\(\mathbb C(x)((s))\). If its lowest remaining term is \(s^j c(x)\)
with \(j<\mu-(q-p)\), its bracket with the lowest term of \(P\)
vanishes. Then \(c^\mu\pi^j\) is a nonzero constant. Writing
\(P=s^{-\mu}\pi\, U\) with \(U=1+O(s)\), the term
\(s^j c\,U^{-j/\mu}\) has this lowest term and commutes with \(P\).
After these subtractions, (8) puts the first remaining term at
\(j=\mu-(q-p)\) and gives (9) with \(F=\pi c\).
A rational \(F\) cannot have a pole away from the roots of \(\pi\).
At a root of \(\pi\) of multiplicity \(d\), if \(F\) has order \(r\),
the putative leading coefficient on the left has factor
\(\mu r-(q-p)d\) at order \(r-1\). For \(r\le0\) the coefficient
is nonzero and is polar, and for \(r\ge2\) the left side vanishes
at the root. Thus \(F\) has a simple zero there and \(d\ne\mu/(q-p)\).
If \(\pi\) has at least two distinct roots, then \(\deg F\ge2\).
The leading term at infinity in (9) gives
\[
 \deg F=\frac{(q-p)\deg\pi}{\mu}.
 \tag{10}
\]

To apply (9), start with the monomial \(u^{10}\) of \(R_0\),
whose coefficient is \(H^2\). For \(i>0\), a term \(b^i u^j\)
ties with it when \(q/p=(10-j)/i\). Use the largest such slope
among terms of \(R\). At this slope \(P\) starts at
\(s^{-\mu}\pi(x)\) with \(\mu=6q-10p\), and
\(\pi\) has degree ten. Whenever the slope exceeds two, its
lowest \(x\)-power in the cases below is positive, and
another term besides \(H^2x^{10}\) makes \(\pi\) have a
nonzero root as well as zero. For \(i\ge5\), (1) puts
\((10-j)/i<2\), so these terms do not enter any of these
first faces. For \(i=3,4\) the general upper bounds from
(1) are \(8/3,9/4\), respectively.

For \(n=1\), \(R_1=2hL_1+M_1\) has \(u\)-order four by (5),
and \(R_2\) has order at least two by (1). The first slope is
six and the face has lower \(x\)-power four. Formula (10)
would give \(\deg F=25/13\).

For \(n=2\) or \(n=3\), (1) and \(R_2=2hL_2+L_1^2+M_2\)
give \(\operatorname{ord}_0 L_1\ge1\). Here either \(M_2=0\)
or \(M_2\) has order four. If \(\operatorname{ord}_0 L_1=1\),
the first slope is four, from \(R_1=2hL_1\) of order six and
\(R_2\) of order two. The lower \(x\)-power is two, and
(10) gives \(15/7\). If \(\operatorname{ord}_0 L_1=2\),
\(R_1\) has order seven and \(R_2\) has order at least four.
The first slope is three, with lower \(x\)-power at least
four, and (10) gives \(5/2\).

We now have \(\operatorname{ord}_0 L_1\ge3\).
If \(n=2\), \(R_2\) has order four from \(M_2\).
The first slope is again three, giving \(5/2\).

It remains to consider \(n=3\). If \(L_2(0)\ne0\),
then \(R_2\) has order five, \(R_1\) has order at least eight,
and \(R_3=2hL_3+2L_1L_2+M_3\) has order at least three.
The first slope is \(5/2\), whose face is
\(\pi=H^2 x^{10}+2H L_2(0)x^5\).
It has six distinct roots, while (10) gives \(\deg F=3\).
Thus \(\operatorname{ord}_0 L_2\ge1\). We obtain
\[
 \operatorname{ord}_0 R_1\ge8,\qquad
 \operatorname{ord}_0 R_2\ge6,\qquad
 \operatorname{ord}_0 R_3\ge4.
\]
If \(R_4\) has order one, the first slope is \(9/4\), with
lower \(x\)-power one, and (10) gives \(25/7\).
Otherwise the first slope is two. Equality occurs at \(R_1\)
if \(L_1\) has order three, at \(R_2\) if \(L_1\) has higher
order and \(L_2\) has order one, and at \(R_3\) from \(M_3\)
when both have higher order.

At slope two, \((p,q)=(1,2)\), so \(P\) starts at \(s^{-2}\pi\).
The degree-ten polynomial \(\pi\) is even, since all terms on
its face satisfy \(j+2i=10\). The indices \(i\ge5\) lie off
the face, so \(x=0\) is a root. The tie with a term of positive
\(b\)-order ensures that \(\pi\) also has a nonzero root.
Write the lowest term of \(Q\) as \(s^{-\nu}\theta(x)\).
The order-fifteen term of \(S_0=h^3\) gives \(\nu\ge3\).
If \(\nu>3\), a term of \(\theta\) must come from \(S_i\)
with \(i\ge1\) and \(j+2i<15\), so \(\deg\theta\le12\).
Yet the two lowest terms have bracket weight
\(1-2-\nu<0\), and their zero bracket gives
\(\theta^2/\pi^\nu\in\mathbb C^*\). Hence
\(\deg\theta=5\nu\), ruling out \(\nu>3\).

Thus \(\nu=3\) and \(\pi\) is a square up to a constant. By (9)
with \(\mu=2\) and \(q-p=1\), its roots cannot have multiplicity
two. The zero root has multiplicity at least four. Evenness
gives a pair of nonzero roots \(x=\pm\alpha\), each of
multiplicity at least four, contrary to \(\deg\pi=10\).
This excludes maximum longitudinal degree eighteen.


## Degrees nineteen through twenty-three

Assume the maximum \(a\)-degree is at most twenty-three. As in the
original direction, the lowest block \(b^{-m}p_0(u)\) of a target
coordinate with \(m>0\) satisfies \(\deg p_0=m\deg f\), where the
companion polynomial \(f\) has simple zeros at \(0,1\).
Thus \(m\le11\). The fold-line restrictions still have degrees
at most three. Any coordinate in an affine target frame has
restriction of degree at least two, and some coordinate has degree
three. These facts follow from birationality on the fold line and
the exclusion of an affine-line nonproperness component, as in the
[weighted obstruction](EXCEPTIONAL_FOLD_WEIGHTED_JACOBIAN_OBSTRUCTION.md#2-applying-the-equation-on-the-whole-fold-chart).

First consider a coordinate \(Q\) of lowest weight \(-11\).
Its lowest coefficient has degree twenty-two, the only possible
multiple of eleven, and has exactly the roots \(0,1\) because its
companion has degree two. If the mate \(P\) has a block at
\(-11\), its coefficient is proportional to the lowest coefficient
of \(Q\) by the vanishing of the lowest bracket. Subtract this
multiple of \(Q\) from \(P\). The resulting mate has lowest weight
\(-m'\) with \(6\le m'\le10\). If \(a_0,c\) are the two
multiplicities of the lowest coefficient of \(Q\), the lowest bracket
forces \(11\mid m'a_0\) and \(11\mid m'c\). Since eleven is
prime, both multiplicities are divisible by eleven. Their sum is
twenty-two, so each is eleven, the forbidden multiplicity in the
original direction. Thus no coordinate has lowest weight \(-11\).

Now consider a coordinate \(Q\) of lowest weight \(-10\).
Its lowest coefficient has degree twenty, the only multiple of
ten between twenty and twenty-three. The two roots \(0,1\) exhaust
the roots since the companion polynomial has degree two. Write
their multiplicities as \(a_0,c\), with \(a_0+c=20\), \(a_0\ge5\)
and \(c\ge4\) by admission. Multiplicity ten is forbidden by the
original directional equation. If \(P\) has a block at \(-10\),
it is proportional to the lowest coefficient of \(Q\) by the
vanishing of the lowest bracket. Subtract the corresponding multiple
of \(Q\) from \(P\). The resulting \(P\) has lowest weight \(-m'\)
with \(6\le m'\le9\), since its restriction still has degree at
least two. The lowest bracket for the resulting pair gives
\(10\mid m'a_0\) and \(10\mid m'c\). If \(\gcd(10,a_0,c)\) were
one or two, this would make \(m'\) a multiple of ten or five,
respectively, outside the interval. It cannot be ten by the
multiplicity bounds and the exclusion of multiplicity ten.
Consequently \(\gcd(10,a_0,c)=5\), \(c\in\{5,15\}\), and \(m'\)
is six or eight.

The mate \(P\) can now have restriction degree at most two, so \(Q\)
has cubic restriction. Its block at \(-9\) has order exactly three
in the local variable \(v=u-1\). Apply (8) with \(u_0=1\),
\(p=1\) and \(q=c-3\), fixing \(Q\) as the first coordinate.
The monomials \(b^{-10}v^c\) and \(b^{-9}v^3\) have the same
lowest order \(-10q+c\). For \(c=15\), blocks of weights
\(-10+i\), \(i\ge2\), lie higher since \(i q\ge24>15\).
For \(c=5\), the block at \(-8\) has \(v\)-order at least
\(\lceil8/3\rceil=3\), so its increase over the weight
\(-10q\) is at least \(2q+3=7>5\). The weights \(-7\) and
higher give an increase at least \(3q=6>5\).
It follows that the angular polynomial is \(x^3(A x^{c-3}+B)\)
with \(AB\ne0\). It has at least two distinct roots. With
\(\mu=10(c-3)-c\), equation (10) gives
\[
 \deg F=\frac{(c-4)c}{10(c-3)-c}
  =\begin{cases}1/3 & c=5,\\ 11/7 & c=15.\end{cases}
\]
Neither value is possible. Thus the lowest weights are at least \(-9\).

After cancelling a common cubic restriction if necessary, take
restriction degrees two and three for \(P,Q\). The lowest weight of
\(Q\) is \(-9\), and its lowest coefficient has degree eighteen,
since it is a multiple of nine between eighteen and twenty-three.
Its only roots are \(0,1\), and the cubic restriction gives order
three at \(1\). Thus the normalized lowest block is \(b^{-9}h^3\).
The lowest bracket now gives the block \(b^{-6}h^2\) for \(P\), by
the order at \(1\) and its degree-two restriction.

Suppose \(Q\) has terms of \(u\)-degree greater than eighteen.
For each such monomial \(b^{-9+i}u^{18+j}\) we have \(i\ge1\)
and \(1\le j\le5\). Choose positive coprime \(d,q\) with
\[
 \frac{d}{q}=\min\left\{\frac{i}{j}:
          b^{-9+i}u^{18+j}\text{ occurs in }Q,\ j>0\right\}.
\]
In the substitution \(u=s^{-d}x,\ b=s^q\), the monomial
\(b^{-9}u^{18}\) ties with the terms realizing the minimum.
Every term of \(u\)-degree below eighteen is higher, since its
\(b\)-weight is at least \(-9\). Thus \(Q\) starts at
\(s^{-\mu}\pi(x)\) with \(\mu=18d+9q\), and \(\pi\) has a root
at zero, a nonzero root, and degree \(18+\Delta\), where
\(\Delta\in\{1,2,3,4,5\}\) is the largest excess on this face.
The monomial of excess \(\Delta\) gives \(d/q\ge1/\Delta\).
On the other hand, (10) with \(p=-d\) gives
\[
 2\le\deg F=\frac{(q+d)(18+\Delta)}{18d+9q},
 \qquad \frac{d}{q}\le \frac{\Delta}{18-\Delta}.
\]
These inequalities give \(18\le\Delta^2+\Delta\), impossible
for \(\Delta\le3\).
For \(\Delta=4,5\), the function
\((1+z)/(18z+9)\) decreases for \(z\ge0\). Using \(z=d/q\ge1/\Delta\)
we obtain
\[
 \deg F\le\frac{(18+\Delta)(\Delta+1)}{18+9\Delta}
  =\begin{cases}55/27 & \Delta=4,\\ 46/21 & \Delta=5.\end{cases}
\]
Each value is less than three. Thus \(\deg F=2\) and (10) yields
\(d/q=\Delta/(18-\Delta)\). A monomial of excess \(\Delta\) on
the face must then have integer height \(i=\Delta^2/(18-\Delta)\),
which would be \(8/7\) or \(25/13\). This excludes both excesses.

We conclude that \(\deg_a Q=18\). The comparison of leading
coefficients of the constant Jacobian then gives \(\deg_a P=12\)
as before, and the degree-eighteen exclusion completes the proof.
