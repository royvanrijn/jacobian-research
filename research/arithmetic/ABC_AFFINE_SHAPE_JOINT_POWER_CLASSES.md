# Affine shapes and joint power classes of primitive rows

A prescribed affine shape of the ratios \(B_i/C_i\) leaves both a
translation and a scale free. Joint power classes on the two
complementary ratios recover the two parameters in succession.
The first curve removes a common factor from the product of the
\(A_i/C_i\), and the second recovers that remaining factor. This
uses the primitive quotient of each row and permits cancellations
of prime valuations between rows.

## Two power-class records determine an affine shape

Fix \(B_0,B_1>0\). Take \(s\) distinct positive primitive rows
\(A_i+B_i=C_i\le X^{B_0}\). Put \(b_i=B_i/C_i\) and \(a_i=1-b_i\).
The \(b_i\) are distinct, since the reduced fraction determines
the entire primitive row. Prescribe the affine shape
\[
       r_i=\frac{b_i-b_1}{b_2-b_1},\qquad r_1=0,\quad r_2=1.
\]
Every occurring shape has polynomial height, uniformly in the
rows. Prescribe a pool \(\mathcal P\) of primes whose product is
at most \(X^{B_1}\), possibly depending on the shape. For a
nonzero rational \(t\), the outside defect support at order \(q\)
is the set of \(p\notin\mathcal P\) with \(q\nmid v_p(t)\).

**AS1.** For \((q,s)=(2,4)\) or \((3,3)\), assume that the
product of the distinct primes in the union of the outside
defect supports of \(\prod_i a_i\) and \(\prod_i b_i\) is at most
\(K\). The number of original ordered \(s\)-tuples is
\(O_{B_0,B_1,\varepsilon}(K X^\varepsilon)\) for each fixed
shape. The same proof gives \(O(X^\varepsilon)\) completions for
a fixed shape and prescribed exact pair of \(q\)-power classes.

To fix both classes, enumerate their outside union radical, at
cost at most \(K\), and at each prime in this support or the pool
select both valuations modulo \(q\). All occurring outside
supports have polynomial height from the input rows. The pool
has polynomial height by assumption. The assignments cost
\(X^\varepsilon\), since each prime has at most \(q^2\) choices
and a polynomial-height integer has \(o(\log X)\) distinct prime
factors, uniformly in that height range. Each product is positive,
so these choices determine its class, including when \(q=2\).
Write \(\Delta_A,\Delta_B>0\) for the resulting \(q\)-free
integer representatives. Their heights are polynomial.

The change of scale
\[
 u=\frac{b_2-b_1}{a_1},\qquad
 a_i=a_1(1-r_i u)
\]
gives \(u\ne0\) and \(1-r_i u\ne0\) for every input row.
Since \(q\mid s\), the first class implies that \(u\) lies under
a rational point of
\[
    v^q=\Delta_A^{-1}\prod_{i=2}^s(1-r_i u).
    \tag{AS2}
\]
The factors have distinct roots because the nonzero \(r_i\) are
distinct. This is a split cubic for \(q=2\) and a split quadratic
for \(q=3\). The other branch lies at infinity. Once \(u\) is
specified, put \(\lambda_i=1-r_i u\).
All \(\lambda_i\) are nonzero and they are pairwise distinct
because \(u\ne0\). The second class now gives a point at
\(t=a_1\) on
\[
    w^q=\Delta_B^{-1}\prod_{i=1}^s(1-\lambda_i t).
    \tag{AS3}
\]
This is a split quartic for \(q=2\) and a split cubic for \(q=3\).
The transformations and coefficients have bounded rational
height exponents at fixed \(B_0,B_1\). The heights of \(u,t\)
are polynomial in \(X\), and the rational roots \(v,w\) also
have polynomial height by their displayed equations.

Both stages have the required uniform subpower count. For
\(q=2\), a split quartic is reduced to a split cubic by moving
one rational branch \(r\) to infinity with \(x=1/(t-r)\) and
\(y=w/(t-r)^2\). Scale the nonzero leading coefficient \(d\) of
the split cubic by \((x,y)\mapsto (dx,dy)\) to get a monic
Weierstrass model. Its rational roots give points of exact
order two. For \(q=3\), move one rational branch of a split
cubic to infinity by \(x=1/(t-r)\) and \(y=w/(t-r)\).
In the resulting split quadratic model write
\(y^3=d x^2+e x+f\). The change
\[
 U=4d y,\qquad V=4d(2dx+e),\qquad
 V^2=U^3+(4d)^2(e^2-4df)
\]
gives a nonsingular Weierstrass cubic. Its rational points over
the two roots of the quadratic have \(U=0\) and \(V\ne0\).
Their tangent is horizontal and has intersection multiplicity
three, giving rational points of exact order three. In both
cases the model coefficients and the resulting abscissas have
polynomial height uniformly. Dujella's bound applies with
height cutoff a sufficiently large fixed power of \(X\) and
counts \(O(X^\varepsilon)\) points on each model. The rational
changes used to move a branch omit only that branch, which is
absent in the positive rows because every \(a_i,b_i\) is nonzero.
[Dujella, Theorem 1.1 and Corollary 1.3](https://arxiv.org/pdf/2312.03655#page=2)

Thus the first stage bounds the number of \(u\), and the second
bounds the number of \(a_1\) for each \(u\), at subpower cost.
The fractions \(b_i\) determine the original primitive rows.
Resizing the inner allowances and restoring the two class records
proves (AS1).

## Weighted classes at other prime orders

**AS4.** The same bound holds at any fixed prime \(q\) with
\(s\ge3\) fixed, for the classes of \(\prod_i a_i^{e_i}\) and
\(\prod_i b_i^{e_i}\), where \(1\le e_i<q\) are fixed and
\(\sum_i e_i\equiv0\bmod q\). In particular, three rows with weights
\(1,1,q-2\) work at every odd prime \(q\). Constants may also
depend on \(q,s\).

The proof uses (AS2) and (AS3) with the indicated exponents on the
factors. The common scale \(a_1^{\sum_i e_i}\) remains a \(q\)-th
power. The first cover has \(s-1\) finite branches and one at
infinity, of order \(e_1\bmod q\). The second has \(s\) finite
branches and none at infinity, since its total degree is divisible
by \(q\). Each branch is totally ramified, and both smooth curves
have genus \((q-1)(s-2)/2\). In genus at least two apply the split
cyclic-cover counting argument of (IC2), which uses cyclotomic
isogeny descent and uniform Mordell–Lang. To put the first cover
in that form, take \(u=u_0+1/Z\) for a bounded integer \(u_0\)
avoiding its finite branches and multiply its root variable by
\(Z^{(\sum_i e_i)/q}\). The branch at infinity is then the factor
\(Z^{e_1}\), and all coefficients retain polynomial height.
The one possible omitted value \(u=u_0\) can be counted directly
in the second stage. The two genus-one cases are those of (AS1),
up to inversion of the cube class when all weights are two.
[Uniform bound for split cyclic covers (IC2)](ABC_SPLIT_INDEX_POWER_CORRESPONDENCES.md#a-joint-count-from-the-branch-set)

## Extraction at the small-budget scale

Set \(n=X^L\) with \(L>0\) fixed before the height limit. If a
family of \(T\ge n^d\) distinct primitive rows has at least
\(n^{-\tau}T^s\) ordered distinct tuples with shapes in a
dictionary of size at most \(n^\alpha\) and union defect radical
at most \(n^\beta\), then
\[
                  sd-\tau\le \alpha+\beta
                  \tag{AS5}
\]
in the counting exponents. The conclusion follows by summing
(AS1) or (AS4), with the inner allowance smaller than any fixed
margin times \(L\). Thus a strict joint record price below \(s/3\)
on tuples of full endpoint exponent would exclude saturation.
All coordinates \(a_i,b_i\) and both joint classes are unchanged
by common multiplication of any original row. The same test
therefore applies to minimally repaired primitive quotients,
including families with distinct repaired targets.

An upper bound on the number of affine shapes does not by itself
supply the union defect cap. On primes private to one row and
outside the pool, a valuation of its \(A_i,B_i\), or \(C_i\) that
is nonzero modulo \(q\) contributes to the union. For \(C_i\) it
contributes through both products and is counted once. A proof
at arbitrary valuation profiles needs an extraction controlling
the shape dictionary and the surviving joint defect together.


## The single combined parity class

**AS6.** One cannot replace the two classes in (AS1) by the
single parity class of \(\prod_i a_i b_i\). For integers
\(7\le m\le X/2\), take the four fractions
\[
            (b_i)_{i=1}^4=
 \left(\frac12-\frac3m,\ \frac12-\frac1m,\
       \frac12+\frac1m,\ \frac12+\frac3m\right).
\]
They are distinct in \((0,1)\), have denominators at most
\(2m\le X\) after reduction, and give actual positive primitive
rows of height at most \(X\). Their affine shape is always
\((0,1,2,3)\). The complementary fractions \(a_i\) permute
the \(b_i\), so \(\prod_i a_i b_i=(\prod_i b_i)^2\).
The first fraction determines \(m\), giving \(\gg X\) different
ordered tuples with trivial combined parity defect. This example
requires no individual radical caps on the original rows.
