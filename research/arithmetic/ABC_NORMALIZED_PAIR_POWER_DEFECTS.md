# Pair recovery from normalized power defects

Removing the common factor of two opposite summands gives
an exact power difference. Its coefficient records need
only the joint power defects of the two rows and a reduced
shift. For a common moving target, the same normalization
gives a primitive ternary equation. For unequal targets
one can also remove both common summand factors after a
canonical lift.

## Lifting two primitive targets

Take two distinct primitive positive rows \(A_i+B_i=C_i\)
with \(C_i\le X\). Define
\[
 H=\operatorname{lcm}(C_1,C_2),\qquad
 A'_i=(H/C_i)A_i,\qquad B'_i=(H/C_i)B_i,
 \qquad g=\gcd(B'_1,B'_2),\qquad U_i=B'_i/g .
\]
Set \(\Delta=U_2-U_1\). The lifted labels are at most \(X^2\).
At a prime of \(H\), at least one \(C_i\) attains the maximum
target depth. That row has a unit multiplier and both its
original summands are units. At any other prime the lift
multipliers are units. Consequently
\[
 \begin{split}
 \gcd(A'_1,A'_2)&=\gcd(A_1,A_2),\\
 \gcd(B'_1,B'_2)&=\gcd(B_1,B_2),\\
 \gcd(A'_1,A'_2,H)&=\gcd(B'_1,B'_2,H)=1.
 \end{split}                                                    \tag{NP1}
\]
If \(U_1=U_2\), the fractions \(B_i/C_i\) coincide. Each is in
lowest terms by primitivity, so the original rows coincide,
contrary to the hypothesis. Thus \(\Delta\ne0\), and
\[
 U_2A'_1-U_1A'_2=\Delta H,\qquad
       \gcd(A'_1,A'_2)\mid \Delta .
                                                               \tag{NP2}
\]
The second assertion uses \(\gcd(\gcd(A'_1,A'_2),H)=1\).
A recovered lifted row gives its original primitive row
by dividing by \(\gcd(A'_i,H)=H/C_i\).

Let \(q\ge2\). For a nonzero rational \(W\) define
\[
 d_q(W;P)=\prod_{\substack{p\notin P\\q\nmid v_p(W)}}p ,
\]
where \(P\) is a fixed pool with prime product of polynomial
height in \(X\). For the two pair ratios and the lifted target put
\[
  D_A=d_q(A'_1/A'_2;P),\qquad
  D_B=d_q(B'_1/B'_2;P),\qquad D_H=d_q(H;P).
\]
The ratio \(B'_1/B'_2=U_1/U_2\) is reduced, so
\(D_B=d_q(U_1U_2;P)\). The least common multiple of the
defect radicals records their support union.

Write
\[
          A'_i=a_i x_i^q,\qquad U_i=b_i y_i^q,\qquad H=h z^q
\]
with each coefficient \(q\)-free. Outside \(P\), a prime of
one of \(a_1,a_2\) is either in \(D_A\) or divides both
\(A'_1,A'_2\). The latter possibility puts it in \(\Delta\)
by (NP2). Thus the support of both coefficients is contained
in \(P\cup\operatorname{Supp}(D_A\Delta)\). The two \(b_i\) have
support in \(P\cup\operatorname{Supp}(D_B)\), since the \(U_i\)
are coprime. The coefficient \(h\) has support in
\(P\cup\operatorname{Supp}(D_H)\). For fixed supports, these
coefficients have subpower many possible complete values.
Their heights are at most \(X^2\), and support enumeration
is uniform even for varying \(q\).
[Uniform fixed-support count](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)

The target record adds only the common-target primes absent
from the two ratio defects. In detail,
\[
 D_H\mid\operatorname{lcm}\bigl(D_A,D_B,
                              \operatorname{rad}(\gcd(C_1,C_2))\bigr).
                                                               \tag{NP3}
\]
For a prime \(p\mid H\) occurring in just one target, the
other original row has a unit in at least one summand
position. In that position the quotient of the two lifted
summands has valuation \(\pm v_p(H)\), because the row with
\(p\) on its target also has both summands units. A prime of
\(D_H\) of this kind therefore occurs in one of \(D_A,D_B\).
The union uses the ratio defects separately, since taking
the defect of a product of ratios can cancel them.

## Recovery from an exact reduced shift

**NP4.** Fix \(\Delta\ne0\). If the two original targets
are equal and may vary, then for \(q=3\) the pairs with
\(\operatorname{lcm}(D_A,D_B,D_H)\le K\) number
\(O_{\varepsilon}(KX^\varepsilon)\).
For each fixed \(q\ge36\) the bound is
\(O_{q,\varepsilon}(KX^{6q^{-3/2}+\varepsilon})\).
If the lifted target \(H\) is fixed, including when it
comes from unequal original targets, one can omit \(D_H\)
and use \(O_\varepsilon(KX^\varepsilon)\) for any \(q\ge2\),
uniformly even in the order. The constants are uniform in
the pool under a fixed polynomial height bound.

Enumerate the union radical. Then it suffices to enumerate
each coefficient on the union with \(P\) and
\(\operatorname{Supp}(\Delta)\) at subpower cost. The
actual equation is
\[
 (b_2a_1)(y_2x_1)^q-(b_1a_2)(y_1x_2)^q=\Delta h z^q .
                                                               \tag{NP4}
\]
For prescribed \(H\), use the nonzero integer \(\Delta H\)
on the right and apply the coefficient-uniform weighted
power-difference bound. Each product root splits at divisor
cost into its \(x\) and \(y\) factors. Recover \(A'_i,U_i\),
then \(B'_i=H-A'_i\), and finally the original rows by (NP1).
The weights and roots have polynomial height. Uniformity
for varying \(q\) follows by applying the fixed-order lemma
to finitely many \(q\le Q\), then choosing the first root
at cost \(X^{4/Q}\) for \(q>Q\). The second root is unique
after the first is fixed. Choose \(Q\) after the exponent
allowance. If searching over all orders, orders above
\(4\log_2 X+1\) have the same defect supports and unit
roots as the integer cutoff just above \(4\log_2 X+1\).
The smaller orders give \(O(\log X)\) choices.
[Zeng, Lemma 3.1, inherited weighted differences](https://hexagonmath.org/pdf/2611.00621v1#page=5)
· [Quadratic case and the kernel-GCD recovery for Erdős #943](ERDOS_943_PELL_RECOVERED_KERNEL_GCD_WEIGHT.md)

For pairs with \(C_1=C_2=C\), all original summands are
coprime to \(C\). Hence the two cross-product roots are
coprime to \(z\). Heath-Brown's uniform primitive
diagonal-cubic theorem gives \(X^\varepsilon\) choices
per coefficient tuple for \(q=3\). For \(q\ge36\)
all root bounds are at most \(X^{2/q}\). Browning–Verzobio,
Theorem 1.2, with the equal box sides \(2X^{2/q}\), gives
\(O_{q,\varepsilon}(X^{6q^{-3/2}+\varepsilon})\).
After the roots and factors are specified, recover \(H=h z^q\)
and the original rows as before. With unequal original targets,
the two cross products in (NP4) can share lift factors with
\(H\), so the moving-target argument requires a different
primitive normalization.
[Heath-Brown, Theorem 1.3](https://arxiv.org/html/2601.07817v1#S1)
· [Browning–Verzobio, Theorem 1.2](https://arxiv.org/html/2608.24512v1#S1)

One can also record the full reduced pair \(U_1,U_2\).
It fixes \(\Delta\), so only the \(D_A,D_H\) defect union
is charged, and \(D_H\) is free for fixed \(H\).
Given \(U_1,U_2\) one uses (NP2) directly. For example, an
outside radical cap \(V\) for \(U_1U_2\) and a defect
union cap \(K\) give a total record cost \(VK\) up to
subpower factors. A bound on the union of all paid primes
can replace \(VK\). The root counts are the same as above.

## Pricing the reduced shift on the original rows

For a moving \(\Delta\), multiply the bound in (NP4)
by its dictionary size \(J\). An outside radical cap on
the shifts also gives a dictionary of that size up to a
subpower factor, since the shifts have polynomial height.
Their relation to an original projective determinant is exact:
\[
 \Delta=
 \frac{B_2C_1-B_1C_2}{\gcd(C_1,C_2)\gcd(B_1,B_2)}
   =\frac{H}{g}\left(\frac{B_2}{C_2}-\frac{B_1}{C_1}\right).
                                                               \tag{NP5}
\]
A prescribed nonzero determinant gives divisor-many possible
\(\Delta\). A bound on the projective difference and a lower
bound on \(g\) bound the numerical height of the reduced shift.
If a lock also gives \(U_1\equiv U_2\pmod D\), the count of
shifts of magnitude at most \(T\) is at most \(2T/D\),
since every nonzero shift is a multiple of \(D\).
The numerical conditions must be strong enough for this
dictionary to be small in the radical scale.

## Two exact gcds and the quadratic boundary

There is an elementary recovery when both reduced pairs
can be recorded. Let \(g_A=\gcd(A'_1,A'_2)\) and
\(g_B=\gcd(B'_1,B'_2)\). Then \(\gcd(g_A,g_B)=1\), since
a common prime would divide \(H\) as well as every summand.
Set \(V_i=A'_i/g_A\) and keep \(U_i=B'_i/g_B\).
The nonzero differences give
\(g_A(V_1-V_2)=g_B(U_2-U_1)\). It follows that
\[
 g_A=\frac{|U_2-U_1|}{t},\qquad
 g_B=\frac{|V_1-V_2|}{t},\qquad
               t=\gcd(U_2-U_1,V_1-V_2).
                                                               \tag{NP6}
\]
The four reduced integers determine both gcds, \(H\),
and both original rows, even if the targets differ.
An outside radical cap \(K\) on \(U_1U_2V_1V_2\) gives
\(O_\varepsilon(KX^\varepsilon)\) pairs with no target charge.

For moving \(H\) the subpower root count fails at \(q=2\),
even for equal original targets and fixed \(U_1,U_2\).
For even integers \(r\ge4\) put
\(z=r^2+2r-1,\ x_1=r^2+1,\ x_2=r^2-2r-1\).
Then \(z^2+x_2^2=2x_1^2\) and \(z>x_1>x_2>0\).
The integers \(z,x_1,x_2\) are odd, and
\(\gcd(z,x_2)\mid4r\) gives \(\gcd(z,x_2)=1\), since
\(z\equiv-1\bmod r\). Also \(\gcd(z,x_1)\mid2(r-1)\)
and \(\gcd(r^2+1,r-1)\mid2\) give \(\gcd(z,x_1)=1\).
Use the common target \(z^2\) and the two summands
\(x_1^2,x_2^2\). Their complements are in the ratio \(1:2\).
The two rows are primitive and have \(D_A=D_H=1\)
at order two with the reduced pair \((1,2)\).
Taking \(P=\{2\}\) gives \(D_B=1\) as well.
There are \(\gg X^{1/4}\) such pairs at height at most \(X\).
The example asserts an obstruction to the unconditional
root count with these records. The small combined original
radical budget is a further restriction.

## A quadratic moving-target count from the common complementary divisor

**NP8.** Keep the common original target \(C_1=C_2=C\)
moving, fix \(\Delta\ne0\), and take \(q=2\). Impose
\(\operatorname{lcm}(D_A,D_B,D_H)\le K\) and an outside
radical cap \(G\) on \(g=\gcd(B_1,B_2)\). If the two
square classes \(h/a_1,h/a_2\) span a group of rank at most
one in \(\mathbb Q^\times/\mathbb Q^{\times2}\), the number
of ordered pairs is
\[
                              O_\varepsilon(KG^{1/2}X^\varepsilon).
                                                               \tag{NP8}
\]
The coefficients \(a_1,a_2,h\) are the genuine squarefree
parts of \(A_1,A_2,C\). The size allowance is uniform under
a fixed polynomial-height bound on the prime pool.

First record the union of defects and the five squarefree
coefficients as in (NP4). With \(b_1,b_2,\Delta\) fixed,
\[
                         b_2y_2^2-b_1y_1^2=\Delta
\]
determines \(U_1,U_2\) at subpower cost by the same
fixed nonzero quadratic norm count. This uses the actual
reduced difference instead of prescribing the original
target. We are left with
\[
                         U_2a_1x_1^2-U_1a_2x_2^2=\Delta h z^2 .
\]
This is a nonsingular projective conic, since all three
coefficients are nonzero. The positive root triple
\((x_1,x_2,z)\) from an actual pair is primitive by the
original within-row coprimality.

If any actual pair occurs, take one primitive root triple
on the conic as a reference. Make it the first column of
a nonsingular integral matrix \(M\), completing with two
standard coordinate vectors. All entries have polynomial
height. In these coordinates the conic equation is
\[
                         x L(y,w)+Q(y,w)=0 ,
\]
where \(L,Q\) are homogeneous integral forms of degrees
one and two. The nonzero form \(L\) and \(Q\) have no common
projective zero, since a common zero would make the conic
reducible. The parametrization in the original coordinates is
\[
 (P_1(s,t),P_2(s,t),P_3(s,t))
                      =M(-Q(s,t),sL(s,t),tL(s,t)).
\]
It covers the conic once projectively. Away from the
reference the inverse is the line direction \([y:w]\).
The tangent direction maps to the reference. The primitive
parameter \((s,t)\) thus has polynomial height, by the
inverse linear map of polynomial height and the bound for
the original roots. The tangent direction itself has
polynomial height since \(L\) does.

For primitive \((s,t)\), the content of the three \(P_i(s,t)\)
has prime support in \(\det M\operatorname{Res}(L,Q)\).
Indeed, at any prime outside that integer the reduction
of \(M\) is invertible. Common divisibility of the three
coordinates would then force \(L(s,t)=Q(s,t)=0\) modulo
the prime. This is impossible by the resultant and
primitivity of \((s,t)\). Dividing the three coordinates by
their content recovers the actual primitive roots up to a
common sign.

Consider the binary quartic
\[
                         F(s,t)=hP_3(s,t)^2-a_1P_1(s,t)^2 .
\]
Its four zeros are simple. On the conic, the corresponding
points are exactly
\[
                   [\,\pm\sqrt{h/a_1}:\pm\sqrt{h/a_2}:1\,] .
\]
They are distinct, and the parametrization is a projective
isomorphism from \(\mathbb P^1\). The quartic is nonzero
and has degree four, so these four distinct points account
for all its zeros with multiplicity one. In rank one the
Galois action has two orbits of size two. In rank zero
there are four rational zeros. Thus the primitive integer
factors of \(F\) have degrees \(2,2\) or \(1,1,1,1\).
Their coefficients have polynomial height. For example
dehomogenize away from a possible linear factor at
infinity. Root bounds for a bounded-degree polynomial,
and divisibility of leading coefficients by Gauss's
lemma, bound the coefficients of its factors by a fixed
power of the original coefficient height.

Enlarge the fixed pool for this branch by the primes of
the coordinate content bound, the scalar content of
\(F\), \(U_1\), and all pairwise resultants of the
primitive factors of \(F\). Its prime product is still
of polynomial height. At a primitive parameter, an
outside prime cannot divide two of the factor values.
Moreover the quartic value is the square of the coordinate
content times \(gU_1\), so every outside prime of a factor
divides \(g\). The product of the outside factor radicals
is at most \(G\). All factor values are nonzero because
\(g>0\).

In the quadratic case record the cheaper factor, whose
outside radical is at most \(G^{1/2}\). Its full value is
specified at subpower cost on each support, since its
height and the fixed pool product are polynomial. A fixed
nonzero value of an irreducible binary quadratic determines
subpower many bounded-height pairs. To see this, complete
the square to get a nonzero norm equation
\[
  (2as+bt)^2-(b^2-4ac)t^2=4am
\]
for the form \(as^2+bst+ct^2\) and prescribed value \(m\).
Here \(a\ne0\) and \(b^2-4ac\ne0\) by irreducibility.
The uniform quadratic count applies, and the coordinate
change is injective. In the linear case take the two
factors with the smallest outside radicals. Their product
is at most \(G^{1/2}\), and their nonzero specified values
determine \(s,t\) by independent linear equations. Thus
in either case the parameter costs \(G^{1/2}X^\varepsilon\).
The primitive roots, \(C\), and both original rows now
follow. Multiplying by the coefficient cost proves (NP8).
[Browning–Verzobio, Lemma 2.3, citing Vaughan–Wooley](https://arxiv.org/html/2608.24512v1#S2.SS1)

In rank two the four zeros form one Galois orbit and give
an irreducible quartic. The rank-two case occurs, for
example, for \(2+11=3+10=13\), whose classes are
\(13/2,13/3\). The method gives a split-factor discount
when those classes generate a smaller group. The square
example above shows why some moving record is needed
even in rank zero. In (NP8) it is a value from the
factorization of the actual common divisor.

## Least-repair multiplicity consequences

In a balanced least-target-repair band, take \(N=Y^{1/M}\)
with fixed \(M,\lambda,a_->0\), and \(H_M(C)\le Y\),
\(\operatorname{rad}(C)\ge C^{a_-\lambda}\). Every original
value is at most \(N^{1/(a_-\lambda)}\).
Choose a set of target radicals and one complete target
per radical, each of multiplicity at least \(N^{u-o(1)}\).
Suppose at least \(N^{2u-\tau-o(1)}\) ordered distinct
pairs in each chosen fibre satisfy (NP4) with a common
fixed pool, union cap \(K=N^\kappa\), and shift dictionary
of size at most \(N^\nu\). For \(q=3\), the number of
these target radicals is
\[
                  \ll_\varepsilon
                    N^{\kappa+\nu+\tau-2u+\varepsilon}.
                                                               \tag{NP7}
\]
Here (NP4) is applied simultaneously across the moving
complete targets. This gives a fixed saving on the
critical tail \(N^{1-u}\) if \(\kappa+\nu+\tau<1+u\)
with fixed margin. For the full reduced-pair version use
its total prime price in place of \(\kappa+\nu\).
For the quadratic rank condition of (NP8), with
\(G=N^\gamma\), the same argument gives (NP7) with
\(\gamma/2\) added to the exponent.

One may instead use (NP6) on a family of distinct primitive
quotients with different targets allowed. If an \(N^{-\tau}\)
fraction of the family's distinct ordered pairs have the
double cofactor cap \(K=N^\kappa\), the family has at most
\(O_\varepsilon(1+N^{(\kappa+\tau)/2+\varepsilon})\) rows.
This applies to families with one representation on each
target when the double cofactor record is economical.
Both counts exclude arbitrary common multipliers by the
recovery of the original primitive rows.

High valuation overlap supplies numerical divisibility
of a shift. The price in (NP4) also needs a bound on
its reduced-shift dictionary and joint power defects,
or a cap for the full reduced-pair alternative.
The double cofactor recovery needs both opposite gcds
to remove sufficient radical mass.
[High-valuation completion stars](ABC_POPULAR_COMPOSITE_PACKET_LCM.md#forced-high-valuation-tails-in-common-completion-stars)
· [Joint least-repair multiplicity range](ABC_FIXED_LETTER_RADICAL_CAPS.md#the-remaining-multiplicity-estimate)

## Recovery on a split elliptic curve with independently moving targets

**NP9.** Take \((q,s)=(2,3)\) or \((3,2)\). For an
ordered tuple of \(s\) distinct original primitive rows
define \(H=\operatorname{lcm}_{i} C_i\),
\(B'_i=HB_i/C_i\), \(g=\gcd_i B'_i\), and \(U_i=B'_i/g\).
Fix the \(s-1\) exact nonzero reduced shifts
\(\Delta_i=U_i-U_1\) for \(i>1\). If
\[
 \operatorname{lcm}\left(
      d_q(U_1/U_2;P),\ d_q\left(\prod_i A_i/C_i;P\right)\right)
                      \le K ,
                                                               \tag{NP9}
\]
the number of ordered tuples is \(O_\varepsilon(KX^\varepsilon)\).
The common lower or upper bounds on original radicals and
the least-repair constraints may be retained. The targets
can coincide or vary independently.

At each \(p\mid H\) some original target attains the
maximum depth. Both its summand and its lift multiplier
are units there. At \(p\nmid H\) every multiplier is a
unit. This proves \(g=\gcd_i B_i\). The \(U_i\) are
distinct since coincident \(B_i/C_i\) would recover the
same primitive row, so all shifts are nonzero. First
record the union radical in (NP9). Write \(U_i=b_i y_i^q\)
for \(i=1,2\). A prime of either \(q\)-free coefficient
outside the ratio defect and the pool must divide both
\(U_1,U_2\), hence divides \(\Delta_2\). Thus fix
\(b_1,b_2\) at subpower cost and solve
\(b_2y_2^q-b_1y_1^q=\Delta_2\) by the weighted
power-difference bound. All weights and roots have
polynomial height since \(H\le X^s\). This fixes every
\(U_i\) by the reduced shifts.

Fix the positive \(q\)-free integer representative \(D\)
of \(\prod_i A_i/C_i\). Its prime support is contained in
the recorded union and the pool, with exponents below
\(q\). Thus its choices are subpower and its height is
polynomial. The parameter \(t=g/H\) satisfies
\[
                     \prod_{i=1}^s(1-U_i t)=D v^q
\]
with rational \(v\), since \(1-U_i t=A_i/C_i\).
Both \(t\) and \(v\) have polynomial rational height.
For \(q=2,s=3\), let \(L=-D U_1U_2U_3\).
Then \((x,y)=(L t,L D v)\) is a point of
\[
                      y^2=\prod_{i=1}^{3}(x-L/U_i).
\]
The monic cubic has integral coefficients of polynomial
height and three distinct rational roots, so the curve
has rational points of order two.

For \(q=3,s=2\), use an affine rational coordinate \(w\)
taking \(t=1/U_1,1/U_2\) to \(w=0,1\). In this
coordinate the equation is \(v^3=\eta w(w-1)\) with
\(\eta\) nonzero rational of polynomial height. Put
\(\alpha=4/\eta\). Then
\((x,y)=(\alpha v,\alpha(2w-1))\) lies on
\(y^2=x^3+\alpha^2\). The point \((0,\alpha)\) has
exact order three by its horizontal flex tangent.
If needed, scaling by a denominator of \(\alpha\)
gives an integral monic model of polynomial height
and preserves polynomial height of the point.

In both cases Dujella's bounded-height theorem counts
\(X^\varepsilon\) points uniformly, since the logarithmic
curve height is \(O(\log X)\) and the point \(x\)-height
is polynomial in \(X\). A curve point determines \(t\)
via \(x\) in the square case and via \(y\) in the cubic
case. Each fraction \(U_i t=B_i/C_i\) gives its original
row in lowest terms. This completes the count. The
argument records a product power class of all original
quotients, so prime cancellation in that class is
permitted. The support union in (NP9) keeps this class
separate from the first two opposite summands' ratio
class.
[Dujella, Theorem 1.1 and Corollary 1.3](https://arxiv.org/html/2312.03655v1#S1)
· [Zeng, Lemma 3.1](https://hexagonmath.org/pdf/2611.00621v1#page=5)

If the shift tuple ranges over \(N^\nu\) records in a
least-repair family and \(K=N^\kappa\), the count of
qualifying tuples is \(O(N^{\kappa+\nu+\varepsilon})\).
A lower bound of \(N^{-\tau}\) on their fraction among
distinct ordered tuples bounds the whole family by
\(O(1+N^{(\kappa+\nu+\tau)/s+\varepsilon})\). This is a
sufficient extraction for the occurrence endpoint even
when each target has one representation. Within fibres
of multiplicity \(N^u\), the same comparison gives a
strict saving on the critical target tail if qualifying
tuples occur at fraction \(N^{-\tau}\) in each fibre
and \(\kappa+\nu+\tau<1+(s-1)u\). The still-needed
extraction is an economical dictionary for the reduced
shifts and the product and ratio defects on a sufficiently
large family of actual tuples.

## A cubic-defect version with no exact shift record

**NP10.** For two distinct original rows with independently
moving targets, keep the common lift and its coprime reduced
opposite pair \(U_1,U_2\) from (NP1). Write
\(\Delta=U_2-U_1\ne0\). If
\[
 \operatorname{lcm}\!\left(
  d_3(U_1U_2|\Delta|;P),
  d_3(A_1A_2/(C_1C_2);P)\right)\le K ,
\]
there are \(O_{\varepsilon}(K X^\varepsilon)\) ordered pairs,
uniformly in fixed pools of polynomial prime product.
The normalized shift may have large raw radical and height.

The three integers \(U_1,U_2,|\Delta|\) are pairwise
coprime. After enumerating the indicated union support,
their cube-free parts can all be fixed at subpower cost
on that support and the fixed pool. The cube roots
are primitive. Heath-Brown's coefficient-uniform
primitive diagonal cubic bound on the exact equation
\(U_2-U_1=\Delta\), with its fixed sign, therefore recovers
\(U_1,U_2\) at subpower cost. The root heights are
polynomial, since the common lift of two targets is at
most \(X^2\). With these reduced labels and the cube
class of \(A_1A_2/(C_1C_2)\) fixed, use the split cubic
cover in (NP9) to recover \(g/H\) and the two primitive
rows. This uses the defect of the actual shifted
difference and leaves no separate shift or target
height record.
[Uniform primitive diagonal cubic count](https://arxiv.org/html/2601.07817v1#S1)
