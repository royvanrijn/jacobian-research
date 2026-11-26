# Counts on planes of divisible power roots

A plane on three projected roots can be used when neither
free letter is an exact power of its root. Suppose the original
free letters have the form \(a x^p u\) and \(b y^p v\) with
\(u\mid x\), \(v\mid y\). The third support first gives a linear
congruence in \(u,v\). The divisibilities \(u\mid x\) and
\(v\mid y\) then separate the two unmarked root factors.
A second valuation layer of the original third letter recovers
the remaining affine parameter at a cost that depends on
three explicit integer records.

## Divisible roots and two integral records

Fix \(p\ge2\). Let \(a,b,\ell_1,\ell_2,\ell_3\) be
nonzero integers of polynomial height in \(X\), and
prescribe a prime pool \(P\) of polynomial prime product.
Count actual positive primitive rows, with the two
additive signs included in \(a,b\), and put
\[
 U=a x^p u,\qquad V=b y^p v,\qquad T=U+V\ne0,\qquad
 \ell_1 x+\ell_2 y+\ell_3 z=0.
\]
Here \(x,y,u,v,z\) are positive integers, \(u\mid x,\ v\mid y\),
and \(z=\operatorname{rad}_{\overline P}(T)\). All complete
values have magnitude at most \(X^{B_0}\) for a fixed \(B_0\).
In particular \(x,y,z\) are pairwise coprime by original
primitivity. One may take one fixed choice of the roots
for each row. The bounds below also count all admissible
tuples for a fixed collection of records.

Enlarge \(P\) by the factors of \(ab\ell_1\ell_2\).
Fix the resulting \(P\)-supported part of \(z\), and absorb
it into \(\ell_3\). This has height-subpower cost.
The remaining \(z\) again denotes the outside radical.
Set \((c,c')\) to be the primitive integer pair proportional to
\((a(-\ell_2)^p,b\ell_1^p)\), using the positive common
divisor to normalize. Modulo \(z\) the plane gives
\(\ell_1 x\equiv-\ell_2 y\). Multiply the original equation
by \(\ell_1^p\) and use that \(y\) is a unit modulo \(z\).
The common content of the two coefficients is also a unit
there. Thus
\[
                  k=\frac{c u+c'v}{z}\in\mathbb Z .          \tag{DR1}
\]
For \(k\ne0\), write \(x=uX_1,\ y=vY_1\). Eliminate \(z\)
between the two relations and use \(\gcd(u,v)=1\). There is
an integer \(s\) satisfying
\[
 k\ell_1 X_1+c\ell_3=s v,\qquad
 k\ell_2 Y_1+c'\ell_3=-s u .                                \tag{DR2}
\]
These two records involve actual exact divisibilities
in the free letters. The second identity would not follow
from bounds on the heights of their correction factors.

If \(k=0\), coprimality fixes \(u,v\) from \(c,c'\).
The two free letters then have a common exact power order
\(p\) on the original plane in \(x,y,z\). If \(s=0\) and
\(k\ne0\), (DR2) fixes \(X_1,Y_1\). The two free letters
then have a common exact power order \(p+1\) on the new
plane \(c u+c'v=kz\). The equal-order plane argument applies
to coprime roots whether or not they are squarefree.
It gives a subpower count in both branches for \(p\ge3\).
For \(p=2,\ k=0\), the only exception from that argument
has the third letter squarefree outside a fixed polynomial-
height pool. The same positive-valuation-band discount applies
to that exception at the balanced endpoint on cheap records.
[Equal-order completion and the quadratic exception](ABC_SHALLOW_REPAIR_RADICAL_PLANES.md#two-equal-residual-powers-on-a-radical-plane)

## Higher third-letter valuations give the third record

Fix nonzero \(k,s\). Let \(z_2\) be the product of
primes of \(z\) with valuation at least two on the original
\(T\), excluding primes of the fixed data
\(ab k c c'\ell_1\ell_2\ell_3\).
Put \(z_0=z/z_2\). Its part supported on those fixed
primes costs a height-subpower factor to specify. Any
other prime on \(z_0\) has valuation one on \(T\), and
its contribution can be charged by its radical. Put
\[
                        L(u)=s u+c'\ell_3.
\]
**DR3.** In the fixed \(k,s,z_0\) record, the integer
\[
                  j=\frac{s u+(p+1)c'\ell_3}{z_2}            \tag{DR3}
\]
exists for every qualifying row. The same numerator is
\(-k\ell_2 Y_1+p c'\ell_3\) by (DR2).

For the proof, regard \(u,z\) as rational indeterminates.
In (DR1) and (DR2) put
\[
 v=\frac{kz-cu}{c'},\qquad
 X_1=\frac{s v-c\ell_3}{k\ell_1},\qquad
 Y_1=\frac{-s u-c'\ell_3}{k\ell_2}.
\]
The original polynomial in \(u,z\) is
\[
                  W(u,z)=a u^{p+1}X_1^p+
                                      b v^{p+1}Y_1^p .
\]
At \(z=0\) the two terms cancel identically by the
definition of \(c,c'\). For generic \(u\), the logarithmic
derivatives with respect to \(z\) of the first and
second terms there differ by
\[
           -\frac{p s k}{c L(u)}+\frac{(p+1)k}{c u}
          =\frac{k\{s u+(p+1)c'\ell_3\}}{c u L(u)} .
\]
Multiplying by the first term at zero gives
\(\partial_z W(u,0)\). Therefore it is a fixed nonzero
rational multiple of
\[
                 u^p L(u)^{p-1}\{s u+(p+1)c'\ell_3\}.
\]
This equality holds as a polynomial, including the
values excluded when computing the logarithmic
derivatives. At each prime dividing \(z_2\), all fixed
nonzero factors in front of the three displayed factors
are units. So are \(u\) and \(L(u)\), since \(u\mid x\)
and \(L(u)=-k\ell_2Y_1\), with \(Y_1\mid y\).
The polynomial \(W(u,z)/z\) evaluated at the actual
row is divisible by each such prime, because it
divides the squarefree \(z\) and occurs on \(T\) to
order at least two. Reducing the polynomial modulo that
prime gives its value at \(z=0\), namely
\(\partial_z W(u,0)\). Hence each prime divides
\(s u+(p+1)c'\ell_3\). Their squarefree product yields
(DR3).

## One-parameter recovery and count

**DR4.** For fixed nonzero \(k,s\), and fixed \(z_0,j\),
the number of original tuples is
\(O_{p,B_0,\varepsilon}(X^\varepsilon)\), uniformly in
polynomial-height coefficients and pools. The implicit
constant allows the fixed height bounds of all records.

To prove it put \(t=z_2\). Equations (DR1)–(DR3) give
successively
\[
 u=\frac{jt-(p+1)c'\ell_3}{s},\qquad z=z_0t,\qquad
 v=\frac{kz-cu}{c'},\qquad
 X_1=\frac{s v-c\ell_3}{k\ell_1},\qquad
 Y_1=\frac{-s u-c'\ell_3}{k\ell_2}.
\]
All four root factors are affine in \(t\) and have
nonzero constant terms. In fact \(u(0)\ne0\) and
\(L(u(0))=-p c'\ell_3\ne0\); the other constant terms
then follow from the identities on \(z=0\).
The expression \(W\) above becomes a polynomial
\(W_0(t)\) with \(W_0(0)=0\). If a tuple is present,
\(W_0\) is nonzero. It cannot be a monomial in \(t\).
For if it were, the two summand polynomials would be
coprime: their common factors would divide that
monomial and neither vanishes at zero. At least one
of \(u,v\) is nonconstant because \(k z_0\ne0\).
If both are nonconstant, all four factors are
nonconstant since \(s\ne0\). The summand degrees
are \(2p+1\), while Mason–Stothers bounds each
degree by at most four from the four affine roots and
zero. If exactly one of \(u,v\) is nonconstant, so
is exactly one of \(X_1,Y_1\), and one summand
has degree \(p+1\ge3\) while the bound is two.
This excludes the monomial in each case.
[Mason–Stothers in the generalized Fermat argument, Browning–Verzobio, Lemma 2.7](https://arxiv.org/html/2608.24512v1#S2.SS1)

Clear the fixed denominators and remove the full factor
at zero. The remaining nonconstant integral polynomial
has a nonzero constant term. At an actual \(t=z_2\), its
value can have no prime outside the given pool, the
recorded \(z_0\) support, the cleared denominator, and
that constant term. Indeed such a prime would divide
\(T\), hence divide \(t\), and then divide the constant
term. All these fixed pool data have polynomial height,
and the remaining polynomial value at an actual root
does too. Uniform supported-value enumeration gives
subpower many values, and a nonconstant fixed-degree
polynomial has bounded fibers. This counts \(t\) and
recovers all the roots and original complete letters.
There is no summation over common multipliers.
[Uniform supported-value count](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)

## Charging a correction instead of the integer \(s\)

The prime exclusions in (DR3) do not require \(s\).
The fixed coefficient of its derivative expression,
for example, is
\(a(-c/(k\ell_1c'))^p k/c\).
Thus after fixing nonzero \(k\) one may fix \(z_0\) and
the quotient \(j\) before choosing \(s\). If \(u\) is
fixed in place of \(s\), (DR3) gives
\[
  s(t)=\frac{jt-(p+1)c'\ell_3}{u},\qquad
  v(t)=\frac{kz_0t-cu}{c'},\qquad t=z_2 .
\]
Recover \(X_1,Y_1\) by (DR2).
At zero all four original factors are nonzero, as above.
If \(j\ne0\), \(v,Y_1\) are nonconstant affine
polynomials, and \(X_1\) is quadratic. The two
summand degrees are \(2p\) and \(2p+1\), with at most
four distinct nonzero roots between them. If \(j=0\),
\(s\) is a nonzero constant, \(v,X_1\) are affine
and nonconstant, and \(Y_1\) is constant. The larger
summand degree is \(p+1\), with at most two distinct
nonzero roots between them. In either case the same
Mason–Stothers argument rules out a monomial in \(t\).
The supported-value count applies. This variant
includes actual rows where \(s(t)=0\), as (DR2)
and the polynomial substitutions still hold.

One can instead fix \(v\) and use the symmetric
integer
\[
       j'=\frac{s v-(p+1)c\ell_3}{z_2}.
\]
Indeed \(cz_2 j+c'z_2 j'=s k z\), and \(c'\) is
a unit at each prime on \(z_2\). When both \(j,j'\)
are cheap, the same identity recovers \(s\) at no
separate record cost. Thus at depth two the third
quotient can be combined with whichever of the
integer \(s\), one of the two actual correction
factors, or the other quotient has the cheapest
available dictionary. One uses \(j\) for the
\(u\)-version and \(j'\) for the \(v\)-version.

## Costs and arbitrary high valuations

Suppose a dictionary of multiplicative coefficients,
plane normals and pools costs \(n^\rho\), where \(n=X^L\)
and \(L>0\) is fixed first. The nonzero \(k,s\) branch
costs at most \(n^{\rho+\kappa+\sigma+\omega+\zeta}X^\varepsilon\)
if the respective records \(k,s,j\) have at most
\(n^\kappa,n^\sigma,n^\omega\) choices and the part of \(z_0\)
of valuation one off the already recorded pool has
radical at most \(n^\zeta\).
The remaining part of \(z_0\) is supported on a fixed pool
after \(k,s\) are specified. Each record may have its
own dictionary depending on the earlier records.
For example one can bound their sizes by counting integers
in the following height ranges:
\[
\begin{split}
 |k|&\le (|c|u+|c'|v)/z,\\
 |s|&\le \min\{(|k\ell_1|X_1+|c\ell_3|)/v,\,
                (|k\ell_2|Y_1+|c'\ell_3|)/u\},\\
 |j|&\le \min\{(|s|u+(p+1)|c'\ell_3|)/z_2,\,
                (|k\ell_2|Y_1+p|c'\ell_3|)/z_2\}.
\end{split}
\]
The thresholds can be set on fixed bins of logarithmic
heights. An outside radical cap on any one of the
integers can replace its height bound by supported
enumeration at the same cost. For \(p\ge3\) the two zero
branches above require only the earlier records and the
subpower completion. The exceptional square branch can
be paid by the separately stated positive-scale discount.

The hypothesis \(u\mid x,\ v\mid y\) allows many different
actual valuations. For any fixed \(p\), if every positive
unrecorded valuation \(e\) is at least \(p(p-1)\), put the
power \(\lfloor e/p\rfloor\) in its projected root \(x\)
and the remainder in \(u\). The remainder is less than \(p\)
and at most \(\lfloor e/p\rfloor\). Thus this prescription
gives \(a x^p u\) with \(u\mid x\) on the actual letter,
with no relation between the exact orders at its different
primes. The lower orders can be put in the coefficient
record at their actual cost. In particular this includes
arbitrary squareful residual parts at \(p=2\), and for
exact adjacent orders it includes the two varying
valuation groups on each free letter.

With such projected roots, the plane hypothesis has the
same simultaneous-residue source as for radical roots.
Write the three actual complete values in the form
\(T_i=R_i^p E_i\), taking \(R_i=x,y,z\) respectively.
The correction \(E_i\) is permitted to be rational, with
numerator and denominator supported on \(T_i\).
At any original \(D_i\mid T_i\), the original sum and the
ratio of complementary corrections specify the unit
\(p\)-th power of the complementary projected root ratio.
Its root residues cost subpower at fixed \(p\).
After fixing supported parts at the moduli, apply the
congruence-lattice determinant bound to the three
projected roots. It uses their numerical heights and
primitivity in the same way as for radicals. In
particular a plane follows when the product of the
three moduli exceeds six times the product of the root
height bounds. Removing the supported parts of \(x,y\)
also preserves the divisibility of the residual \(u,v\)
by the residual projected roots. This supplies an
arithmetic count on those residue branches when the
record inequalities above have a strict margin.
It leaves open how to find sufficiently inexpensive
planes and quotient records across all balanced families.
[Simultaneous residue lattice](ABC_SHALLOW_REPAIR_RADICAL_PLANES.md#rank-and-occupancy-of-simultaneous-residue-records)

## A deeper quotient for the divisible-root projection

The same argument uses all the valuations of the third
letter up to any fixed depth. Take \(J\ge2\) and nonzero
\(k,s\) as above. Let \(z_J\) be the product of primes on
\(z\) having valuation at least \(J\) on the original \(T\),
away from the fixed data in (DR3) and from \(s p(p+1)\).
Fix \(z_0=z/z_J\). Only the primes off the enlarged fixed
pool with valuation below \(J\) need new radical records
here. Put \(t=z_J\) and use the polynomial
\[
                   F(U,t)=\frac{W(U,z_0 t)}{z_0t}.
\]
The division is polynomial because \(W(U,0)=0\).
At \(t=0\), by the derivative calculation in (DR3),
its nonzero factors are a constant, powers of \(U,L(U)\),
and \(s U+(p+1)c'\ell_3\). Hence
\[
                       u_0=-\frac{(p+1)c'\ell_3}{s}
\]
is a simple root, away from the prescribed primes.
In particular \(u_0\ne0\) and \(L(u_0)=-pc'\ell_3\ne0\).

**DR5.** Let \(Q(t)\) be the polynomial of degree less
than \(J-1\) with constant coefficient \(u_0\) such that
\(F(Q(t),t)=0\bmod t^{J-1}\). Solve recursively by dividing
by the nonzero derivative \(F_U(u_0,0)\).
At fixed \(p,J\) the coefficients have polynomial
height and their denominators are units on \(z_J\).
For any positive common denominator \(D\) with
those prime supports the original data give an integer
\[
                     j_J=\frac{D u-DQ(z_J)}{z_J^{J-1}} .
\]
For each fixed \(j_J,k,s,z_0\), there are \(O(X^\varepsilon)\)
original tuples, uniformly at fixed polynomial height
bounds. The quotient count at depth two is (DR3)
up to using a reduced rational denominator.

Indeed let \(\pi\mid z_J\). The case of (DR3)
gives \(u\equiv u_0\bmod\pi\), and \(z_0\) is a unit
there. Since \(T\) has order at least \(J\) at
\(\pi\), \(F(u,z_J)\) vanishes modulo \(\pi^{J-1}\).
The same is true of \(F(Q(z_J),z_J)\) by construction.
Their divided difference with respect to \(u\)
is a unit, reducing to \(F_U(u_0,0)\).
This proves the asserted integer divisibility,
because \(z_J\) is squarefree.

Upon fixing the quotient, \(u=Q(t)+j_Jt^{J-1}/D\)
is a polynomial in \(t=z_J\), and \(v,X_1,Y_1\)
are obtained by the same linear expressions as in
(DR4). Their constant terms are nonzero as there.
If \(\deg u=e\ge2\), the relation \(c u+c'v=kz_0t\)
makes \(\deg v=e\), and both \(X_1,Y_1\) also have
degree \(e\). The two summands have degree
\((2p+1)e\), with at most \(4e\) distinct roots
between them. As before, if \(W_0\) were a monomial
they would be coprime and Mason–Stothers would give
\((2p+1)e\le4e\), impossible for \(p\ge2\).
When all factors are affine the two cases of
(DR4) apply. Thus for every occurring tuple
there is a nonconstant factor after removing the
factor at zero. The fixed-support argument proves
the uniform count, with the new denominator and
nonzero constant term included in the pool.

If \(D,\|DQ\|_\infty\le n^\Theta\),
\(u\le n^{r_u}\), and \(z_J\ge n^{r_J}\), the
number of depth-\(J\) quotient records is at most
\[
 O_J\!\left(
              n^{\max\{0,\Theta+r_u-(J-1)r_J,\Theta-r_J\}}
        \right).
\]
This follows by dividing the bound for \(Du-DQ(z_J)\)
by \(z_J^{J-1}\). Alternatively one can charge the
outside radical of the quotient at its own price.
Both bounds pay separately for the primes below
depth \(J\) that were moved to \(z_0\).

### Deeper recovery with the correction \(u\) fixed

One can use the deeper third valuations while still
pricing \(u\) in place of \(s\). With \(k,u\) fixed,
regard \(W(u,z)\) as a polynomial in \(s,z\) via
(DR1)–(DR2). Its quotient by \(z\) has at zero the
same expression in (DR3). At its nonzero root
\(s_0=-(p+1)c'\ell_3/u\), the divisor \(L(u)=s u+c'\ell_3\)
has value \(-p c'\ell_3\) and the derivative of
\(s u+(p+1)c'\ell_3\) is \(u\). Hence the root
in \(s\) is simple off the prescribed primes
and the support of the fixed \(u\). The latter
support is already disjoint from \(z\). Use the
same recursion and exact quotient at depth \(J-1\),
now for \(s\) rather than \(u\). The prime exclusions
depend on \(k,u\) and the fixed coefficients, so they
can be prescribed before choosing \(s\).

After fixing this quotient, \(s(t)\) is a polynomial
with constant term \(s_0\ne0\), and \(v(t)\) is
affine and nonconstant. If \(\deg s=e\ge1\),
then \(X_1\) has degree \(e+1\), and \(Y_1\)
has degree \(e\). The larger summand degree
is \(p(e+1)+1>2e+2\), while the two summands
have at most \(2e+2\) distinct roots between them.
If \(s\) is constant, use the affine case of
the \(u\)-version of (DR4). Nonzero constant
terms at zero and Mason–Stothers again exclude a
monomial output. Hence this deeper version has
subpower completion per quotient as well. Its
height price is obtained from (DR5) by replacing
\(r_u\) by the height cap for \(s\) and using
the coefficients of the reconstruction for \(s\).
