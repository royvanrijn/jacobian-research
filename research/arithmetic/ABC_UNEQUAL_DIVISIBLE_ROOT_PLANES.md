# Distinct divisible power orders on a radical plane

Two distinct power orders give an unequal valuation identity even
when both free letters have divisible corrections. The smaller
order correction divides the first quotient forced by the actual
third radical. After fixing that quotient, one more record gives
a subpower count without any lower bound on the third letter's
valuations.

## The smaller correction divides the quotient

Fix \(p>q\ge2\), nonzero signed integers \(a,b\), an integral
plane normal \((\ell_1,\ell_2,\ell_3)\) with every entry nonzero,
and a prime pool \(P\). All fixed data and the prime product of
\(P\) have height at most a fixed power of \(X\). Consider
positive integers \(x,y,u,v,z\) with \(u\mid x\), \(v\mid y\),
signed complete letters
\[
  U=a x^p u,\qquad V=b y^q v,\qquad T=U+V\ne0,
  \qquad z=\operatorname{rad}_{\overline P}(T),
\]
and the plane equation \(\ell_1x+\ell_2y+\ell_3z=0\).
The absolute values of the three complete letters are bounded
by a fixed power of \(X\) and are pairwise coprime. The signs
allow any ordering or orientation of an original positive
primitive sum. In particular \(x,y,z\) are pairwise coprime.

Enlarge \(P\) by the primes of \(ab\ell_1\ell_2\).
Fix the part of the old \(z\) on those primes and absorb it
into \(\ell_3\), so \(z\) again denotes the radical outside
the enlarged pool. This refinement costs a height-subpower
factor uniformly. Put \(d=p-q\) and take the primitive
integer pair \((c,c')\) proportional to
\((a(-\ell_2)^p,b\ell_1^p)\), dividing both by their
positive gcd. All the coefficient data have polynomial
height at fixed \(p,q\).

**DU1.** The two exact divisibilities are
\[
             k=\frac{c y^d u+c'v}{z}\in\mathbb Z,
             \qquad v\mid k .                              \tag{DU1}
\]
Indeed, reduce the plane modulo \(z\) to substitute
\(\ell_1x=-\ell_2y\) in \(\ell_1^pT\).
The factor \(y^q\) and the normalized coefficient content
are units modulo \(z\). This proves integrality of \(k\).
Since \(d\ge1\) and \(v\mid y\), the numerator is divisible
by \(v\). Coprimality of \(v,z\) proves the second assertion.
This uses the original third support and places no constraint
on its positive valuations.

## Subpower recovery after two records

**DU2.** For fixed nonzero \(k\), either of the following
additional records gives \(O_{\varepsilon}(X^\varepsilon)\)
tuples, uniformly in all the polynomial-height data at fixed
\(p,q\):

- The integer \(s\) defined below.
- The exact larger-order correction \(u\).

In either case first record \(v\mid k\) at divisor cost.
With \(x=uX_1\) and \(y=vY_1\), eliminating \(z\) between
the plane and (DU1) gives
\[
 u(k\ell_1 X_1+c\ell_3 y^d)+
 v(k\ell_2 Y_1+c'\ell_3)=0 .
\]
Coprimality of \(u,v\) gives the unique integer \(s\) with
\[
   k\ell_1 X_1+c\ell_3 y^d=s v,\qquad
   k\ell_2 Y_1+c'\ell_3=-s u .                            \tag{DU2}
\]

For the fixed \(s\) record, take \(u\) as the parameter.
Equations (DU2) and (DU1) express \(y,X_1,z\) as
polynomials in \(u\) over \(\mathbb Q\). If \(s\ne0\),
\(y\) is affine and nonconstant with nonzero constant
term, \(X_1\) has degree \(d\), and \(z\) has degree \(d+1\).
The complete value becomes
\[
                  W(u)=a u^{p+1}X_1(u)^p+b v y(u)^q .
\]
The summand polynomials are coprime. At the only root of
\(y\), both \(u\) and \(X_1=s v/(k\ell_1)\) are nonzero.
The larger summand has degree \(p(d+1)+1\).
If all geometric roots of \(W\) belonged to \(z\), the
number of distinct roots of the union of the two summands
and \(W\) would be at most \((d+1)+1+(d+1)\).
Mason–Stothers would give maximum degree at most \(2d+2\),
a contradiction. If \(s=0\) and a tuple occurs, \(y\) and
\(X_1\) are nonzero constants. Then \(z\) is nonconstant
affine and \(W\) is a binomial in \(u\) with \(p+1\)
distinct roots. At least one root is outside \(z\) here too.

For the fixed \(u\) record, take \(y\) as the parameter.
The two polynomials
\[
       z(y)=\frac{c u y^d+c'v}{k},\qquad
       x(y)=-\frac{\ell_2 y+\ell_3 z(y)}{\ell_1}
\]
determine the other roots. The summands of
\(W(y)=a u x(y)^p+b v y^q\) are coprime because
\(x(0)\ne0\). If \(d\ge2\), the degree of \(W\) is
\(pd\), while Mason–Stothers would give at most \(2d\)
if every root of \(W\) were a root of \(z\).
If \(d=1\) and \(x\) is nonconstant, these bounds become
\(p\) and \(2\). If \(d=1\) and \(x\) is constant,
they become \(q\) and \(1\). Thus every branch leaves a
nonconstant factor of \(W\) with no common root with \(z\).

For both records, clear the fixed denominators and take the
factor of \(W\) remaining after removing all irreducible
factors over \(\mathbb Q\) shared with \(z\), to their full
multiplicities in \(W\). Its primitive integer normalization
\(G\) has polynomial coefficient height since all degrees
are fixed. Include in the prescribed pool the primes in
the denominators, the normalized contents, and the
nonzero resultant of \(G\) with the cleared \(z\).
A prime outside this fixed pool that divides \(G\) at a
tuple would divide \(W=T\), so it would also divide the
value of \(z\). The resultant excludes it. Uniform
supported-value counting gives subpower many nonzero
possible values of \(G\), and each has a bounded number of
preimages since \(G\) is nonconstant. This recovers the
remaining roots and the complete original letters.
[Mason–Stothers in Browning–Verzobio, Lemma 2.7](https://arxiv.org/html/2608.24512v1#S2.SS1)
[Uniform supported-value counting](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)

## The zero quotient

**DU3.** The branch \(k=0\) has a subpower count without
specifying \(s\) or either correction. Dividing the first
equation in (DU1) by \(v\) yields
\[
                  c\,y^{d-1}Y_1 u=-c'.
\]
If \(d\ge2\), this fixes \(y,Y_1,u\), and then \(v=y/Y_1\),
up to a bounded product of divisor counts of \(|c'|\).
The plane now makes \(x\) an affine function of \(z\) with
nonzero slope. Hence \(T=a u x(z)^p+b v y^q\) is a
polynomial with \(p\) distinct roots. Removing its possible
zero root leaves a nonconstant polynomial whose actual
values are supported on a fixed pool after adjoining the
denominator and its nonzero constant term. The same count
as in (DU2) applies.

If \(d=1\), the relation fixes \(Y_1,u\) at divisor cost.
The original free letters \(a u x^p\) and
\(b Y_1^q v^{q+1}\) then have the common power order
\(p=q+1\ge3\) on the roots \(x,v\). These roots are
coprime, and the plane is
\(\ell_1 x+\ell_2 Y_1v+\ell_3 z=0\).
Apply the equal-order plane count. Its proof only needs the
first two roots to be coprime and the third coordinate to be
the exact outside radical. The expanded pool may contain
the now fixed coefficient supports.
[Equal-order plane count (SR5)](ABC_SHALLOW_REPAIR_RADICAL_PLANES.md#two-equal-residual-powers-on-a-radical-plane)

## Normalized record cost

**DU4.** Let \(n=X^L\) with \(L>0\) fixed. Suppose the
plane, coefficients and pool range over a dictionary of
size at most \(n^\rho\). For each such record let the
quotient range over at most \(n^\kappa\) values. On each
nonzero quotient record, suppose either \(s\) has at most
\(n^\sigma\) choices or \(u\) has at most \(n^\alpha\)
choices. If both bounds are available use the smaller
dictionary. With the respective additional exponent \(b\)
the resulting count is
\[
                      O_\varepsilon(n^{\rho+\kappa+b}
                                            X^\varepsilon).
\]
The zero branch costs only the plane, coefficient and pool
dictionary. The individual dictionaries may depend on
earlier records. Their sizes may come from the exact
inequalities
\[
       |k|\le\frac{|c|y^d u+|c'|v}{z},\qquad
       |s|\le\frac{|k\ell_2|Y_1+|c'\ell_3|}{u}.
\]
Outside-pool radical budgets on a record can replace
its height budget, by the same supported-value count.

For example, on the stripped record suppose \(p=q+1\),
\(\max(|c|,|c'|)\le n^\gamma\), \(y,v\le n^h\),
\(z\ge n^{h-\xi}\), and \(u\le n^\alpha\).
Then \(|k|\le2n^{\gamma+\xi+\alpha}\), and the
fixed \(u\) route gives exponent
\(\rho+\gamma+\xi+2\alpha\).
A fixed margin below \(1/3\) gives a strict coefficient
saving on this sector with no separate cost for \(v\).
For squarefree projected roots in exact adjacent orders
on each letter, \(u\) and \(v\) are the products of the
primes in the higher of the two orders on their own
letter. Thus the example permits that part of the
lower-order letter to use its entire root budget.
More generally a root projection with \(u\mid x\), \(v\mid y\)
can encode many orders at once. At base order \(p\) any
valuation \(e\ge p(p-1)\) can be written
\(e=p\lfloor e/p\rfloor+r\) with
\(0\le r\le\lfloor e/p\rfloor\), and likewise for \(q\).

All estimates take the height limit at fixed budget and
fixed orders. The orders may be selected as functions of
the fixed budget before that limit. The count requires
affordable plane and quotient records. A diffuse profile
by itself supplies neither of these, and polynomial height
of the normalized coefficients does not imply the
small normalized coefficient bound in the example.
[Equal-order divisible roots](ABC_DIVISIBLE_POWER_ROOT_PLANES.md)
[Unequal exact powers](ABC_UNEQUAL_POWER_RADICAL_PLANES.md)
