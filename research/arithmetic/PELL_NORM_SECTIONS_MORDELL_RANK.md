# Mordell ranks from Pell units and actual class directions

The two reflected norm sections in
[CP1](QUARTIC_PELL_TWO_CLASS_DIRECTIONS.md)
and the Pell prime path of
[PI3](PELL_PERIOD_RING_CLASS_INERTIA.md#a-simultaneous-family-with-genuine-prime-classes)
give points on the associated \(j=0\) elliptic curves.
The distinction between the unit direction and the ordinary
class directions becomes independence in \(3\)-isogeny descent.
The first family has generic rank three. A variant of the
moving-prime construction has generic rank two while keeping
the same prime ideal and affine Frobenius controls.

## The quartic family has exact generic rank three

**Theorem PM1.** Use \(P,Z_\pm,H_\pm\) from CP1 and put
\(F=P^2+1\). The curve
\[
 \mathcal E:\ y^2=x^3+F(t)
\]
has rank exactly three over \(\mathbb Q(t)\) and rank exactly six
over \(\mathbb Q(\zeta_3)(t)\). The three points
\[
 (-1,P),\qquad (H_+,P+Z_+),\qquad (H_-,P+Z_-)
 \tag{1}
\]
are independent and generate the image of rational points
in the descent quotient by the dual \(3\)-isogeny over
\(\mathbb Q(t)\). Every sufficiently large specialization
in CP1's progression has three independent specialized
points. Its \(3\)-isogenous partner \(y^2=x^3-27F\) has
the same ranks in these assertions.

For \(k\ne0\) write \(E_k:y^2=x^3+k\). The usual \(3\)-isogeny
goes to \(E_{-27k}\), with dual back to \(E_k\).
For a nonsquare \(k\), its dual descent map on affine points
of \(E_k\) takes the value \(y+\sqrt{k}\) modulo cubes in
\(\mathbb Q(\sqrt{k})\), up to conjugating the square root.
This is [Bhargava–Elkies–Shnidman, Proposition 24 and Remark 25](sources/bhargava-elkies-shnidman-mordell-isogeny-selmer.pdf),
with the sign of the square root chosen as above.
The statement is over a field of characteristic prime to six
and applies over \(\mathbb Q(t)\). The kernel is the dual
isogeny image.

The norm identity from CP1 gives the two points of (1).
The value on the first point is \(U=P+R\) for \(R=\sqrt F\).
The value on the other two points, divided by \(U\), is
\[
             \frac{P+Z_\pm+R}{P+R}
                =1+Z_\pm(R-P)=\overline{\alpha_\pm}.
\]
These three classes are independent after specialization
by CP1. For positive nonsquare \(k\) the curve \(E_k\)
has no rational three-torsion. The \(3\)-division polynomial
is \(3x(x^3+4k)\), the first root would require \(k\)
a square and the others would give \(y^2=-3k\).
The same reasoning applies over \(\mathbb Q(t)\) with the
real ordering for large \(t\).
Since \([3]\) factors through the \(3\)-isogenies, independence
in the descent quotient proves independence of the points.
Any relation between the generic points would specialize on
CP1's good fibers. Thus the generic rank is at least three.

For the upper bound, work over \(\mathbb Q(\zeta_3)\), choose
a prime above seven, and reduce the twist
\(E_F\) of \(E_1\). The cyclic cover \(v^6=F(t)\) has good
reduction there because the reduction of \(F\) is squarefree
of degree eight, with unit leading coefficient, and the
cover is tame. One way to see smoothness is to take its
normalization over the integral \(t\)-line. At simple branch
points \(v\) is a local parameter. In the infinity chart
\(s=1/t\) the exponent of \(s\) in the Kummer function is
\(-8\). After taking the tame local roots of the unit part
one obtains the degree-three ramified charts on the two
infinity branches. These are smooth. The same tame argument
applies below to a squarefree sextic, whose infinity is
unramified.

A rational point of the isotrivial curve gives an equivariant
map from the smooth completion of this cyclic cover to
\(E_1\) via \((x,y)\mapsto(x/v^2,y/v^3)\).
The map extends over the local ring by good reduction of
\(E_1\) and the Néron mapping property for the smooth proper
source. A nonconstant map stays nonconstant, since the degree
of the pullback of an ample line bundle is constant on the
fibers. An equivariant constant map is the origin. In fact
a generator of \(\mu_6\) fixes no finite point under its
action on \(E_1\). Hence reduction injects the section group
over \(\mathbb Q(\zeta_3)(t)\) into the section group over
\(\mathbb F_7(t)\). In particular it injects their free parts.

The finite-field curve is nonconstant as a curve over
\(\mathbb F_7(t)\), though isotrivial. It has eight simple
geometric additive fibers at the zeros of \(F\), and one more
at infinity. After scaling near infinity the constant
coefficient has valuation four. Inertia acts tamely through
a nontrivial subgroup of \(\mu_6\) at these nine fibers,
so their geometric conductor exponents are all two.
The finite-field \(L\)-polynomial has degree
\(18-4=14\), and its reciprocal roots have weight two.
We use [Ulmer, Definition 1.1.4 and Theorems 9.3 and 12.1](sources/ulmer-function-field-l-polynomial-rank-excerpt.pdf)
for the nonconstant terminology, polynomial, and
rank upper bound by order of vanishing.

If \(\lambda_i\) are the reciprocal roots, the power sum
\(\sum_i\lambda_i^j\) is the negative sum of the local
elliptic traces over \(t\in\mathbb P^1(\mathbb F_{7^j})\),
using zero at the additive places. This follows by
expanding the logarithm of the Euler product.
The first seven values are
\[
       31,\ 177,\ 1126,\ 15845,\ 143491,\ 1124502,\ 2706595 .
\]
The [checker](verify_pell_norm_sections_mordell_rank.py)
computes them by grouping the \(F(t)\) into sextic residue
classes and computing the trace once per class.
Newton's identities and reciprocity give
\[
 \begin{split}
 (X-7)^6(X^2+7X+49)\bigl(&X^6+4X^5+42X^4+686X^3\\
                         &+2058X^2+9604X+117649\bigr).
 \end{split}
 \tag{2}
\]
The middle coefficient is nonzero, which fixes the sign
of the reciprocity relation. The rank over
\(\mathbb F_7(t)\) is at most six by (2).

On \(\mathcal E\) over \(\mathbb Q(\zeta_3)(t)\), the
endomorphism \([\zeta_3]-[\zeta_3^{-1}]\) is a degree-three
isogeny and interchanges invariant and anti-invariant
free parts for complex conjugation. Thus the rank over
\(\mathbb Q(\zeta_3)(t)\) is twice the rank over
\(\mathbb Q(t)\), giving the upper bounds.
Over \(\mathbb Q(t)\) the group modulo three is of dimension
three, and the descent quotient already has the three
independent images (1), so these generate that quotient.

## Exact specialized conductors for the quartic family

**Corollary PM2.** On CP1's progression put
\[
 R_6(F)=\prod_{\substack{\ell\ge5\\v_\ell(F)\not\equiv0\bmod6}}\ell.
\]
The two nonisomorphic isogenous rank-at-least-three curves
\(E_F,E_{-27F}\) have exact common conductor
\[
                              2^6 3^2 R_6(F)^2 .
 \tag{3}
\]
There are \(\gg X^{1/16}/\log X\) distinct
\(\mathbb Q\)-isogeny classes from this family of conductor
at most \(X\), and \(\gg_\varepsilon X^{1/16-\varepsilon}\)
distinct conductors for every \(\varepsilon>0\).

At \(\ell\ge5\), unit twists become good unramifiedly.
For a nonzero valuation modulo six the inertia action
is a nontrivial subgroup of \(\mu_6\) on the Tate module.
This is tame and has no fixed directions, so the conductor
exponent is two. At three, \(F\equiv1\bmod9\), a local
sextic power, which gives the local curve \(E_1\) of
conductor exponent two. At two, \(v_2(F)=1\) and
\(F/2\equiv1\) or \(5\bmod8\), so \(E_F\) is locally
isomorphic to \(E_2\) or \(E_{10}\). The unit sixth
powers in \(\mathbb Q_2\) are exactly the unit squares,
and the conductor exponent of each of these two models
is six. These small local conductor values are checked
by the exact local reduction algorithm.
The isogeny preserves the conductor; the sign change
excludes an isomorphism since the ratio would have
to be a rational sixth power.

Up to parameter size \(T\), (3) is \(O(T^{16})\).
The same CP1 multiplicity argument gives
\(\gg T/\log T\) quadratic fields. A \(\mathbb Q\)-isomorphism
between \(j=0\) models identifies their constants up
to a sixth power. Within each \(\mathbb Q\)-isogeny
class, the distinct curves from these quadratic fields
are bounded uniformly as well: see the isogeny-class
argument below. For a fixed conductor \(N\), there
are at most \(72\cdot5^{\omega(N)}\) sextic-free constants
giving \(j=0\) curves. Each prime at least five that
divides the constant is recorded in the conductor,
and its exponent has five possibilities. There are
at most six possibilities each at two and three
and two signs. This proves the bound on the number
of distinct conductors.

For the isogeny-class claim, let \(E_k,E_{k'}\) be
isogenous over \(\mathbb Q\), and choose an isomorphism
\(i:E_k	o E_{k'}\) using a sixth root of \(k'/k\)
over an algebraic closure. Any rational isogeny has the
form \(i\lambda\), with \(\lambda\) a nonzero geometric
endomorphism of \(E_k\). All those endomorphisms are
defined over \(\mathbb Q(\zeta_3)\), because their algebra
is generated by the order-three automorphism. If
\(i\lambda\) is rational, then \(\sigma(i)\lambda=i\lambda\)
for all \(\sigma\) fixing \(\mathbb Q(\zeta_3)\).
Surjectivity of \(\lambda\) forces \(i\) to be defined
over \(\mathbb Q(\zeta_3)\). Thus \(k'/k\) is a sixth
power there. For a positive ratio the real sixth root
is also in \(\mathbb Q(\zeta_3)\) since all sixth roots
of unity are. It belongs to the real subfield
\(\mathbb Q\). This proves that the positive models
from distinct quadratic fields are pairwise nonisogenous.

## Raising the Pell coefficient while preserving the prime path

**Theorem PM3.** Fix a finite nonempty set \(S\) of odd
primes containing three, integers \(a_p\ge1\), and
positive multiples \(b_p\) of three. Write
\[
 B=\prod_{p\in S}p^{b_p},\qquad m=\prod_{p\in S}p^{a_p}.
\]
For each \(e\mid m\), there is a prime progression with a
compatible positive-density Frobenius slice which gives
\(\gg X^{1/(2m)}/(\log X)^2\) distinct real quadratic
fields \(K\) of discriminant at most \(X\) with the
following properties.

- Each \(p\in S\) ramifies in \(K\), with
  \(K_p\simeq\mathbb Q_p(\sqrt p)\) and \(v_p(U_D)=b_p\).
  The moving prime has ordinary class order \(m\) and
  class order \(Bm\) in the order of conductor \(B\).
  Its \(m\)th power there generates the whole relative
  inertia kernel.
- The exact affine degree, discriminants, and ray conductors
  of AW1–AW4 hold with the same norm element, and the
  moving prime has prescribed translation Frobenius order \(e\).
- The constructed \(E_d:y^2=x^3+d\) and \(E_{-27d}\)
  have rank at least two, with two explicit independent
  points on \(E_d\). They have exact conductor
  \(2^6 3^5R_6(d)^2\).
  For conductor at most \(Y\) this gives
  \(\gg Y^{1/(4m)}/(\log Y)^2\) distinct isogeny
  classes, and \(\gg_\varepsilon Y^{1/(4m)-\varepsilon}\)
  distinct conductors.

Choose \(c,s\) as in PI3 with \(Bc=2^l s+1\),
where \(l\) is even and divisible by all \(p-1\) for
\(p\in S\), prescribing \(c\equiv1+2p\bmod p^2\)
at each \(p\). For the larger coefficient \(B\) these
conditions prescribe \(2^ls+1\) modulo \(p^{b_p+2}\)
and still give \(s\equiv-1\bmod p\). Additionally
choose \(c\equiv25\bmod27\), which agrees with
\(c\equiv7\bmod9\). The progression for \(s\) remains
reduced. Put \(A=2(Bc-1)\). Define
\[
 n=\frac{q^m-c^2}{A},\qquad
 d=n(B^2 n+2)=r_0^2D,\qquad T=1+B^2n .
\]
At two prescribe \(n\bmod32\) with
\(n\equiv2\bmod8\) and \(d/8\equiv1\bmod8\).
For \(n=2+8j\), the latter residue is
\((B^2+1)/2+6j\bmod8\), so such a residue exists.
Choose \(q^m\equiv c^2+A n\bmod32A\)
using the unit power permutation, and \(q\equiv1\bmod
\operatorname{lcm}(9,m)\).
The auxiliary Chebyshev non-power residues of PI3
apply with this \(B,c,A\) and are compatible with
these prescriptions. The prime-to-\(m\) residue
conditions and the fixed \(s\)-valuation of the
rational constant \(c(Bc-2)/(Bc-1)\) allow the same
Chebotarev selection as AW3.

The two norm elements are
\[
 u=T+B\sqrt d,\quad N(u)=1,\qquad
 \alpha=Bn+c+\sqrt d,\quad N(\alpha)=q^m .
\]
The \(p\)-adic square class of \(D\) is unchanged from
PI3–AW3. In fact \(v_p(n)=1\) and
\(n/p\equiv2\bmod p\), so \(d/p\) is a square
unit at all \(p\in S\), and \(p\nmid r_0\).
The Chebyshev residues ensure that if \(u=\epsilon_D^h\)
then \(p\nmid h\) for every \(p\in S\).
Comparing its radical coefficient \(B r_0\) with that
of the fundamental unit proves \(v_p(U_D)=b_p\).

For these local squareclasses the scalar unit quotient
\[
 O_{K_p}^{\times}/\bigl(\mathbb Z_p^\times(1+p^j O_{K_p})\bigr)
\]
is cyclic of order \(p^j\), generated by the class of
\(1+\sqrt D\). Here is a check of the exponent of
this generator. Its \(p\)th power has radical
coefficient of valuation one. For \(p\ge5\) the first
binomial term has valuation one and all others
are deeper. For \(p=3\) the first and last terms
at valuation one do not cancel since \(D/3\equiv1\bmod3\).
Subsequent \(p\)th powers increase the radical
coefficient valuation by one, and powers prime to
\(p\) preserve it. This gives order \(p^j\), equal to
the local quotient size. Consequently at conductor
\(p^{b_p}\) the unit image is zero, and the
relative kernel is the whole cyclic group of
order \(p^{b_p}\). At higher conductors the global
unit image has precisely the extra order, so the
primary groups stabilize there. The scalar unit
quotients for conductor \(B\) are the product of
these cyclic groups.

At \(q\) there is exactly one branch on which
\(\alpha\) vanishes, since \(Bn+c\equiv
c(Bc-2)/A\bmod q\ne0\) for large \(q\).
At each \(p\in S\) the local scalar class of
\(\alpha\) generates the cyclic quotient of
conductor \(p^{b_p}\), since both its scalar and
radical coefficients are units. The class cannot
be a \(p\)th power times a global unit there.
The same argument as in PI3 proves exact ordinary
order \(m\) and exact lift order \(Bm\).
In particular a cyclic \(p\)-power quotient that
detects the full relative inertia of order \(p^{b_p}\)
has degree at least \(p^{a_p+b_p}\).
On the affine side \(\gamma=\alpha/\bar\alpha\)
has exactly the local depth in AW1 and exactly
the same moving-prime residue as AW3. This
proves the asserted affine conclusions.

Since \(B\) is a cube and \(3\mid m\), the two points
on \(E_d\) are
\[
        (B^{-2/3},T/B),\qquad (q^{m/3},Bn+c).
 \tag{4}
\]
Their descent classes are \(u\) and \(\alpha\).
The first is a noncube unit by the choice at three.
The second maps to the nonzero ideal class
\([\mathfrak q^{m/3}]\) of order three under the
virtual cube class map. They are independent in the
descent quotient, hence as rational points, since
\(d>0\) is nonsquare and there is no rational
three-torsion.

At two, \(d/8\equiv1\bmod8\) gives the local model
\(E_8\) of conductor exponent six. At three, \(B\)
is divisible by \(27\), \(q^m\equiv1\bmod27\), and
\(c^2\equiv4\bmod27\). Hence \(d\equiv3\bmod27\),
giving the local model \(E_3\) of conductor exponent
five. Away from two and three the argument of PM2
gives the claimed conductor. The fixed progression
and Frobenius condition give the count as in PI3.
Up to \(q\) of size \(H\), \(d\ll H^{2m}\) and
the conductor is \(O(H^{4m})\). Distinct \(K\)'s
give distinct positive isogeny classes by the
argument of PM2. Its fixed-conductor multiplicity
bound gives the last count.

[PF2–PF3](PELL_UNIT_PRIMITIVITY_FROM_TRACE_EQUATIONS.md) strengthen the
PM2 and PM3 isogeny-class counts by showing that almost all specialized
Pell units are fundamental.
[CP5](QUARTIC_PELL_TWO_CLASS_DIRECTIONS.md) also strengthens the
PM2 count of distinct conductors, using a large prime visible in each
selected conductor.

## An exact generic rank-two moving-prime example

**Corollary PM4.** In PM3 take \(S=\{3\}\), \(B=27\),
\(m=3\), \(c=511\), and \(s=3449\).
With \(q\) replaced by the rational parameter \(t\)
in the formulas for \(n,d\), the two points (4)
generate a finite-index subgroup of \(E_d(\mathbb Q(t))\),
and its rank is exactly two. Over
\(\mathbb Q(\zeta_3)(t)\) its rank is exactly four.

The choices satisfy \(Bc-1=4s\), \(c\equiv25\bmod27\).
A prime fixture is \(q=23247793\).
Here \(n\equiv18\bmod32\), \(q\equiv1\bmod9\), and the
half-trace \(T\) is \(3\bmod5\), outside the
\(T_3\) Chebyshev image over \(\mathbb F_5\).
The fixed rational residue \(C=c(Bc-2)/(Bc-1)\)
is a noncube at this \(q\), so the moving prime
has wild affine Frobenius order three.
Applying PM3 to the fixture gives a specialization
where the two points are independent, proving the
generic lower bound.

Use the same reduction argument as PM1 at the
split prime nineteen. The sextic \(d(t)\) has
squarefree reduction of degree six, with unit
leading coefficient. The six simple geometric zeros
give conductor degree twelve and the fiber at
infinity is good. The \(L\)-polynomial has degree
eight. The first four Frobenius trace sums,
including the good point at infinity, are
\[
                        38,\ 722,\ 54872,\ 260642 .
\]
The reciprocal characteristic polynomial is
\((X-19)^4(X^2+19X+361)^2\).
Thus the finite-field rank is at most four.
Reduction of the equivariant morphisms bounds
the characteristic-zero rank over \(\mathbb Q(\zeta_3)(t)\)
by four. The degree-three anti-invariant endomorphism
then gives the asserted rank two over \(\mathbb Q(t)\).

The same generic rank proof works for every fixed
\(B=3^{3b}\) with \(b\equiv1\bmod6\) and every admissible
coefficient \(c\equiv511\bmod19\). These models have the
same good reduction at nineteen. The extra congruence
on \(c\) is compatible with the auxiliary prime selection
in PM3, since \(Bc-1\equiv2\bmod19\) and \(l\) can be
taken to be two. Thus the exponent \(b_3=3b\) of the
relative ring class inertia can be arbitrarily large
with generic rank exactly two.

## The exact generic lattices

**Corollary PM5.** The three points of (1) form a
\(\mathbb Z\)-basis of \(E_F(\mathbb Q(t))\), with height
matrix
\[
             \frac13
              \begin{pmatrix}
                8&2&2\\
                2&8&-4\\
                2&-4&8
              \end{pmatrix}.
 \tag{5}
\]
The two points of (4) form a \(\mathbb Z\)-basis over
\(\mathbb Q(t)\) in PM4 and its fixed-congruence extensions,
with height matrix \(2I_2\). The same two bases are
\(\mathbb Z[\zeta_3]\)-bases over \(\mathbb Q(\zeta_3)(t)\).
Both generic groups are torsion-free.

A torsion point on either generic curve gives a constant
equivariant map on the sextic cover after base change,
since its image is finite. The only equivariant constant
is the origin, as in PM1. This proves torsion-freeness.

Use the [Shioda height formula, §§11.6–11.8 and Table 4](sources/schuett-shioda-elliptic-surfaces-0907.0298.pdf).
For a section on these minimal models, with
\(\chi=2\) in PM1 and \(\chi=1\) in PM4, twice the
intersection with the zero section is
\(\max(d_{\rm den},d_{\rm num}-2\chi)\), where
\(d_{\rm num},d_{\rm den}\) are the reduced polynomial
degrees of its \(x\)-coordinate.
Indeed the finite negative \(x\)-valuations count
with half their pole order, and at infinity the
minimal \(x\)-coordinate is \(s^{2\chi}x\).
For PM1 the only reducible fiber is \(IV^*\) at
infinity, with component group of order three.
Each tripled point has no correction there.
In the order of (1), the numerator and denominator
degrees of their tripled \(x\)-coordinates are
\((24,16),(24,20),(24,20)\).
For the tripled sums, in the order \((0,+),(0,-),(+,-)\),
they are \((60,56),(60,56),(24,20)\).
The self-pairing of a tripled point is \(2\chi\) plus twice
its zero-section intersection. Polarization gives (5)
and determinant \(32/3\).

Every pairing on the full generic lattice has denominator
dividing three, because Table 4 gives only thirds as
corrections at \(IV^*\). Thus if (1) spans a subgroup
of index \(I\), \(I^2\) divides \(288=27(32/3)\).
The index is prime to three by PM1's independent
descent images. To test two, specialize the three
sections at the listed small rational \(t\)'s and
choose the displayed simple root \(\theta\) of
\(x^3+F\) modulo the indicated prime. The three columns
are the square-residue exponents of \(x(P)-\theta\).

| Prime | \(t\) | \(\theta\) | \(P_0\) | \(P_+\) | \(P_-\) |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 5 | 0 | 2 | 1 | 1 | 1 |
| 11 | 2 | 8 | 1 | 1 | 0 |
| 11 | 3 | 3 | 1 | 0 | 0 |

The root lifts by Hensel, and the local
two-descent character \(x-\theta\) kills doubles.
For instance its homomorphism property follows
by evaluating the chord polynomial at \(\theta\)
and multiplying the three \(x-\theta\) values,
which gives a square. The special values at the
origin and at \((\theta,0)\) follow by cancelling
the corresponding even local order.
Any generic double would specialize to a double
at each good rational fiber. The matrix is
invertible over \(\mathbb F_2\), so \(I\) is odd.
Thus \(I=1\).

In PM4 all fibers are irreducible. The \(x\)-degree
pairs for the two points and their sum in (4) are
\((0,0),(1,0),(3,2)\) for every admissible parameter
choice specified above. This follows also from their
intersection: the two sections are disjoint at
finite places and meet to order one at the smooth
infinity after rescaling. The self-pairings are
two and two and the self-pairing of the sum is
four. The full lattice is even integral by the
height formula without corrections. If the displayed
sublattice of determinant four had index two, the
full even rank-two lattice would have determinant
one, impossible modulo four. Thus these two points
are a basis.

For the integral claim over the quadratic constant extension
let \(\delta=\zeta_3-\zeta_3^{-1}\) and
\(M=E(\mathbb Q(\zeta_3)(t))\). Complex conjugation
\(\sigma\) acts semilinearly and \(M^+=E(\mathbb Q(t))\)
is its invariant group. Write \(N=\mathbb Z[\zeta_3]M^+\).
For any \(w\in M\), both \(w+\sigma w\) and
\(\zeta_3w+\zeta_3^{-1}\sigma w\) lie in \(M^+\).
Their combination gives \(\delta w\in N\).
If \(P=\delta w\in M^+\), torsion-freeness implies
\(\sigma w=-w\). The anti-invariant points identify with
rational points on \(E_{-27k}\). Multiplication by
\(\delta\) on these points is, under this identification,
the dual three-isogeny up to sign. This follows by
factoring \(\delta\) on \(E_k\) through the isomorphism
from its three-isogenous twist over \(\mathbb Q(\zeta_3)\);
the two compositions are \(\delta^2=-3\) and \([3]\)
up to that sign. The saturated descent images just
proved therefore give
\[
                      \delta M\cap M^+=3M^+ .
\]
Modulo \(\delta N\) every class in \(N/\delta N\) is
represented by a point of \(M^+\), since
\(\zeta_3\equiv1\bmod\delta\). Given \(w\in M\),
subtract \(\delta\) times a point of \(N\) from
\(\delta w\) to land in \(M^+\). The result belongs
to \(3M^+\subseteq\delta N\), so \(\delta w\in\delta N\)
and \(w\in N\). Thus \(M=N\), with the stated
Eisenstein bases.

## Scope and verification

The finite trace computations, local conductor values and
actual prime fixture are in the
[Sage checker](verify_pell_norm_sections_mordell_rank.py).
The finite-field rank bounds and the generic lower bounds
give exact generic ranks. For individual fibers only the
displayed lower bounds are proved. The conductor counts
are counts of these specific high-rank families. They do
not bound the unrestricted moments of elliptic curves
ordered by conductor.
