# Codegrees with moving right and target power classes

A pair of primitive rows can have different right power classes
and different target power classes. Normalize the right pair by
its exact gcd, then recover the reduced omitted-root ratio in
a projective residue class. Elimination of the common right gcd
puts every prime shared by the two targets on one known integer.
Thus only the target defect on primes exclusive to one target
needs a new record. The same argument includes repeated left
labels and counts all distinct original row pairs in a column.

## A joint coefficient record

Fix \(e\ge2\), two signs \(\sigma,\tau\in\{-1,1\}\) and
prescribed positive integers \(P_1,P_2\), which can coincide.
Consider ordered pairs of distinct positive primitive rows
\[
 A_i+\sigma B_i=\tau T_i,\qquad A_i=P_i y_i^e,\qquad
 \gcd(A_i,B_i)=1,\qquad i=1,2.
\]
All complete values are at most \(X^B\) for fixed \(B>0\).
For a prescribed polynomial-height prime pool \(\mathcal P\),
define the two paid prime sets
\[
 \begin{split}
 \mathcal E_B&=\{p\notin\mathcal P:\,
            v_p(B_1/B_2)\not\equiv0\bmod e\},\\
 \mathcal E_T&=\{p\notin\mathcal P:\,
            p\text{ divides exactly one of }T_1,T_2,\quad
            v_p(T_1/T_2)\not\equiv0\bmod e\}.
 \end{split}
\]
Let \(H\) bound the product of the primes in the union.
The pool can depend on the prescribed labels and modulus.

Impose a projective lock
\[
 A_1/B_1\equiv A_2/B_2\pmod D
\]
at a prescribed positive integer \(D\), with all \(A_i,B_i\)
units modulo \(D\). A sufficient hypothesis is
\(D\mid\gcd(T_1,T_2)\). With \(g=\gcd(y_1,y_2)\)
and \(G=\gcd(B_1,B_2)\), write
\[
 y_i=g u_i,\qquad B_i=G L_i x_i^e
\]
using the canonical \(e\)-free parts \(L_i\) of \(B_i/G\).
Let \(Y,Z\ge1\) be upper bounds on \(u_i,x_i\), respectively.
The bounds \(X^{B/e}\) suffice for both types of roots, and
smaller bounds may apply to the gcd-stripped roots.

**MT1.** For \(e=3\), the number of pairs under these hypotheses
is
\[
                O_{B,B',\varepsilon}\!\left(
                  H X^\varepsilon\left(1+\frac{(YZ)^2}{D}\right)
                                      \right),                \tag{MT1}
\]
where \(X^{B'}\) bounds the prime product of \(\mathcal P\).
For each fixed \(e\ge36\), there is an additional factor
\(R^{3/\sqrt e}\), where
\(R=2\max\{2,X^{B/e},ZX^{B/e}\}\). One can replace \(X^{B/e}\)
in the second term of this maximum by a known bound on
\(g\), and the \(X^{B/e}\) multiplying \(Z\) by a known
bound on the \(e\)-th roots of the two target power parts.
The constants are uniform over prescribed coefficients,
pools and moduli at those heights.

If the logarithm of the relative original odds
\[
                      \frac{A_1/B_1}{A_2/B_2}
\]
lies in a prescribed interval of length \(\ell>0\), one may
replace the bracket by
\(1+((YZ)^2/D)\min(1,\ell/e)\).
Finally the same argument interchanges \(B_i,T_i\) with a
change of the signs. This exchanges which side contributes
its relative defect and which side contributes only its
exclusive defect. For that alternative, use the corresponding
size bounds and projective unit lock for the new right side.

## Recovery and the role of shared target primes

Record the union \(\mathcal E_B\cup\mathcal E_T\) at cost
\(O(H)\). The normalized right terms \(B_1/G,B_2/G\) are
coprime. Outside the pool the combined support of their
\(e\)-free parts is exactly \(\mathcal E_B\). This follows
because at each prime at most one normalized term has positive
valuation, and the difference is the original relative
valuation. Thus fixing \(L_1,L_2\) costs an inner subpower
factor on the recorded and prescribed supports, uniformly
since all complete values have polynomial height.
[Supported power enumeration (EB5)](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)

The fraction \(r=u_1x_2/(u_2x_1)\) is reduced. A prime in
both numerator and denominator would divide both \(u_i\),
both \(x_i\), or one same-row pair \(u_i,x_i\). The first two
possibilities contradict the two gcd normalizations. The last
contradicts primitivity of the original row. All these roots
and \(L_i\) are units modulo \(D\), since they divide the
corresponding \(A_i\) or \(B_i\). The projective lock therefore
implies
\[
                         r^e\equiv\frac{P_2 L_1}{P_1 L_2}
                                           \pmod D.
\]
There are \(X^{o(1)}\) root residues for fixed \(e\), using
the structure of the unit group of each prime power and the
polynomial height of \(D\). This also covers \(D=1\).
In one residue distinct reduced fractions have nonzero
determinant divisible by \(D\). Both coordinates are at
most \(YZ\). After scaling to the unit square, this gives
angular separation at least \(D/(2(YZ)^2)\) and hence
\(O(1+(YZ)^2/D)\) fractions. The original relative odds
equal \(P_1L_2 r^e/(P_2L_1)\). On the optional real interval
the range of \(\log r\) has length \(\ell/e\), which reduces
the angular range by a factor \(O(\min(1,\ell/e))\).
The usual divisor bound splits the recovered numerator
and denominator into \(u_1,x_2\) and \(u_2,x_1\), respectively,
at inner subpower cost.
[Omitted-coordinate packing (MP2)](ABC_MOVING_TARGET_POWERFUL_CODEGREES.md#a-bound-for-a-fixed-pair-of-left-labels)

Now set
\[
                   E=P_1L_2(u_1x_2)^e-P_2L_1(u_2x_1)^e .
\]
Multiplying each row by the other normalized right term
and subtracting gives
\[
                 E g^e=\tau\bigl(L_2x_2^eT_1-L_1x_1^eT_2\bigr).
\]
Here \(E\ne0\). In fact \(E=0\) would give \(A_1/B_1=A_2/B_2\);
two primitive positive pairs with that ratio coincide, and
the signs then determine the same target. This uses distinct
original rows alone and does not require distinct left power
classes. The same identity shows \(\gcd(T_1,T_2)\mid E\), since
\(g\) divides each \(A_i\) and hence is coprime to each \(T_i\).

Write \(T_i=K_i z_i^e\) with \(K_i\) canonically \(e\)-free.
A prime of \(K_1 K_2\) outside the prescribed pool either lies
on exactly one target and belongs to \(\mathcal E_T\), or
divides the fixed nonzero integer \(E\). In particular a
prime shared by the targets incurs no new target record
even when their two valuations differ modulo \(e\). Fix the
two coefficients at subpower cost on these supports. The
resulting equation is
\[
              E g^e-\tau L_2K_1(x_2z_1)^e
                         +\tau L_1K_2(x_1z_2)^e=0 .
\]
Its roots are primitive since \(g\) is coprime to each \(x_i\)
by same-row coprimality, and to each \(z_i\) by target
coprimality. All three coefficients are nonzero. The
coefficient-uniform cubic count gives the asserted subpower
completion at \(e=3\). At \(e\ge36\) the primitive generalized
Fermat bound of Browning and Verzobio gives \(R^{3/\sqrt e}\)
up to an arbitrarily small height allowance in the equal
box of side \(R\). The roots determine the two targets and
the two left terms. The original signed equation then gives
the right terms. Multiplying the counts proves (MT1).
[Heath-Brown, Theorem 1.3](https://arxiv.org/html/2601.07817v1#S1)
· [Browning–Verzobio, Theorem 1.2](https://arxiv.org/html/2608.24512v1#S1)

## All original pairs in a column

**MT2.** Suppose a family has \(W\) distinct original rows,
an \(A\)-element dictionary of left kernels, and a partition
into \(B_1\) right columns. The columns can be defined by
canonical right \(e\)-power classes. More generally the
joint relative-right and exclusive-target defect requirement
can be imposed directly on the eligible pairs in a column.
Use a uniform bound \(\Delta\) from (MT1) for each prescribed
ordered pair of left labels, including repeated labels, with
the same prescribed data. Suppose at least a fraction
\(\rho>0\) of all ordered distinct row pairs within columns
satisfy those hypotheses. Then
\[
                  W\ll B_1+A\sqrt{B_1\Delta/\rho}.            \tag{MT2}
\]
The formula counts actual rows in each cell, including cells
with many representations. If \(D\) or the pool varies among
a prescribed dictionary of charts, one can sum the appropriate
(MT1) bounds over that dictionary in place of \(\Delta\).

For each left label pair, (MT1) bounds all eligible actual
pairs, so the total is at most \(A^2\Delta\). If \(d_j\)
is the number of original rows in column \(j\), the eligible
pairs number at least \(\rho(\sum_j d_j^2-W)\ge
\rho(W^2/B_1-W)\). This proves the result, including the
case in which every column has at most one row. The
diagonal that is removed in this calculation consists only
of identical original rows. There is no maximum cell
multiplicity in the conclusion. This replaces the separate
multiplicity bound in (MP6) when the eligibility hypothesis
is imposed also on pairs having the same left label.
[Moving-target graph with distinct left labels (MP6)](ABC_MOVING_TARGET_POWERFUL_CODEGREES.md#counting-all-original-rows)

As a specific endpoint comparison, let \(n=X^L\),
\(e=3\), \(A,B_1\le n^{1/6}\) and \(D\ge(YZ)^2\).
If \(\rho\ge n^{-\theta}\) and the union defect cap is
\(H\le n^\beta\), the original row exponent is at most
\(\max\{1/6,\,1/4+(\beta+\theta)/2\}\) before inner
allowances. This gives a strict saving at every fixed
budget if \(\beta+\theta<1/6\) with a fixed margin.
A right column whose members have the same cube class
has \(\mathcal E_B=\varnothing\). Within that column the
target condition can be phrased using only the defect
of target primes which are absent from the paired target;
arbitrary positive target valuations on common supports
are permitted.

If the eligible fraction instead holds among all distinct
row pairs, take one column. Then
\(W\ll 1+A\sqrt{\Delta/\rho}\), so with \(A\le n^{1/6}\)
and the same lock the strict sector is
\(\beta+\theta<1/3\). In both comparisons the
fraction, divisor and defect conditions concern actual
complete rows. Extracting them on a sufficiently large
family at broad unrestricted profiles remains an
arithmetic obligation.


## Cubic completion at large omitted-power orders

**MT3.** If the fixed omitted-power order \(e\) is a multiple
of three, (MT1) and (MT2) hold with no generalized Fermat
loss when the target prime set is instead
\[
 \{p\notin\mathcal P:\ p\text{ divides exactly one of }T_1,T_2,\
                                      3\nmid v_p(T_1/T_2)\}.
\]
The right prime set still uses the relative defect at order
\(e\). In particular this keeps the short projected-root
heights and modulus packing at a large omitted-power order
while charging only the cubic defect on exclusive target primes.

After the same packing at order \(e\), fix the canonical
cube-free parts \(K_i\) in \(T_i=K_i z_i^3\) on the union
record, the pool and the primes of the fixed nonzero \(E\).
The eliminated identity is the diagonal cubic
\[
        E(g^{e/3})^3-\tau L_2K_1(x_2^{e/3}z_1)^3
                       +\tau L_1K_2(x_1^{e/3}z_2)^3=0.
\]
Its roots have collective gcd one for the same reason as
in (MT1), and their heights are at most \(X^{2B/3}\)
using the original complete-height caps. Thus the cubic
completion gives the original bound without an additional
power of \(X\). The root-residue count is used at fixed
\(e\), which may be chosen after fixing the budget but before
the height limit.
[Heath-Brown, Theorem 1.3](https://arxiv.org/html/2601.07817v1#S1)


## Complete-label graphs at arbitrary valuation profiles

**MT4.** Take a family of actual positive primitive rows
\(U+V=T\) of polynomial height, with dictionaries of sizes
\(A,B_1\) for the complete summands \(U,V\), respectively.
In a pair with the same complete \(V\) and distinct left
labels, charge only the outside prime set on exactly one
target whose valuation fails to be divisible by a selected
order \(q\ge2\). The order may vary with the pair. For each
fixed ordered pair of left labels, the number of such
pairs with this exclusive defect radical at most \(H\) is
\(O_{B,B',\varepsilon}(H X^\varepsilon)\). The modulus in
(MT1) is unnecessary. The pool may depend on the fixed
left labels and has polynomial prime product.

Indeed, the actual target difference is the nonzero integer
\(T_2-T_1=U_2-U_1\). Every prime of both targets belongs to
this fixed difference. After selecting the exclusive defect
support, the two \(q\)-free coefficient values are therefore
recoverable at subpower cost on the enlarged pool, even
when their valuations differ on the common target support.
For a fixed pair of coefficients, the uniform weighted
power-difference bound counts the two target roots at
subpower cost at any fixed \(q\). Alternatively, if \(q\)
exceeds a fixed cutoff depending on the height exponent
and the inner allowance, the first target root already has
subpower range and fixes the second. This gives uniformity
in the order. To allow it to vary, only a logarithmic
number of orders need be listed, since above the maximal
possible positive valuation both roots are one and the
exclusive defect set is constant. Each completed pair of
targets determines the common complete summand.
[Uniform weighted differences, Lemma 3.1](https://hexagonmath.org/pdf/2611.00621v1#page=5)

If at least a fraction \(\rho>0\) of all ordered distinct
row pairs sharing \(V\) meet the defect cap, the original
row count is \(O(B_1+A\sqrt{B_1 H/\rho}\,X^\varepsilon)\),
by the same moment as (MT2). In particular the minimum
complete-label sizes \(A,B_1\le n^{1/6}\) yield a strict
one-third sector when \(H\le n^\beta,\ \rho\ge n^{-\theta}\)
with \(\beta+\theta<1/6\) by a fixed margin. This specialization
places no restriction on the valuation profiles of the
complete summands. On a full-exponent family in this
dictionary range, the corresponding partial-pair inequality
can instead be read as a necessary lower bound on the
exclusive target defect at the prescribed popularity.
