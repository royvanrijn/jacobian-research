# Simultaneous shallow repair residues and radical planes

A corrected least-repair residue at one original divisor
specifies a ratio of two radicals modulo that divisor.
Using three such records at once can force the residual
radical vectors into a rational plane even when the
divisors are too shallow to reconstruct any one pair.
Large plane fibers force explicit divisibility patterns
on the normal. If two of the actual residual letters
have the same exact power, the radical plane can then be
combined with the original sum. The third letter can
have arbitrary residual valuations.

## Rank and occupancy of simultaneous residue records

Use original positive primitive triples \(a+b=c\le X\).
Write their ordered complete letters \(T_1,T_2,T_3\),
with fixed signed relation. Fix respective divisors
\(D_i\mid T_i\), which are pairwise coprime. Allow a
prescribed prime pool \(P\) of polynomial prime product
containing all divisors' prime supports. By uniform
fixed-support enumeration one may first fix the entire
\(P\)-supported parts of all three complete letters at
subpower cost. Let \(P_i\) be the radicals of those
fixed parts. Write \(R_i=\operatorname{rad}(T_i)=P_i x_i\),
where \(x_i\le H_i\) are the pairwise coprime squarefree
radicals outside \(P\). All \(H_i\) are at most \(X\).

At each \(D_i\) fix the residue of the corrected
repair multiplier ratio of the other two letters,
as in (TL6). Repair orders may vary between the
three congruences but are fixed in the height limit.
For cyclic \(i,j,k\) the ratio \(R_j/R_k\) is a
power-root residue modulo \(D_i\), with the sign of the
corresponding additive orientation accounted for in
(TL6). There are \(X^{o_X(1)}\) root choices together at
fixed orders. Since \(P_j,P_k\) are units modulo \(D_i\),
each branch fixes
\[
                   x_j\equiv \mu_i x_k\pmod {D_i}
                       \qquad(i,j,k\ {\rm cyclic}) .       \tag{SR1}
\]
The unit residues \(\mu_i\) are fixed on the branch.
[Corrected multiplier residue recovery (TL6)](ABC_TWO_TARGET_SUMMAND_LOCKS.md#a-direct-weighted-residue-count-for-primitive-quotients)

**SR2.** Consider the distinct support vectors
\((x_1,x_2,x_3)\) on one branch. If they span a
three-dimensional space, their number is
\[
                     O\left(\frac{H_1H_2H_3}{D_1D_2D_3}\right).
                                                               \tag{SR2}
\]
In particular the span has dimension at most two if
\(D_1D_2D_3>6H_1H_2H_3\). If the span has dimension
at most one, it contains at most one such positive
primitive vector. If it has dimension two with primitive
integral normal \(\boldsymbol n=(n_1,n_2,n_3)\), either
one \(n_i=0\) and the number of actual original rows is
\(X^{o_X(1)}\), or the number of support vectors is at most
\[
 O\left(1+\frac{1}{Q(\boldsymbol n)\,\|\boldsymbol n^*\|}\right),
 \qquad
 Q(\boldsymbol n)=\prod_{i=1}^3
                    \frac{D_i}{\gcd(D_i,|n_i|)},\quad
 \boldsymbol n^*=
       \left(\frac{n_i}{H_jH_k}\right)_{i=1}^3 .           \tag{SR3}
\]
Each support vector has only subpower many original
rows, uniformly. The assertions count actual supports
and do not require a residue class to contain any row.

The congruences (SR1) define a lattice of covolume
\(D_1D_2D_3\). The map to the product of the three
quotient rings is surjective by the Chinese remainder
theorem, since at each \(D_i\) its linear form has a
coefficient one. If three independent candidate vectors
exist, let their span over the integers be a sublattice
\(\Gamma\). Its fundamental parallelepiped lies in the
box with side bounds \(3H_i\) in absolute value. Within
each coset of \(\Gamma\), attach that parallelepiped to
each candidate. The translates are disjoint up to their
boundaries and lie in a box of volume \(O(H_1H_2H_3)\).
Multiply the per-coset volume bound by the index of
\(\Gamma\) in the congruence lattice to get (SR2).
Independently, the determinant of any three support
vectors is divisible by \(D_1D_2D_3\) and is bounded in
absolute value by \(6H_1H_2H_3\).

If \(n_i=0\) in a plane of positive triples, the other two
coefficients are nonzero and the plane fixes the ratio
\(x_j/x_k\). This fixes each radical separately by their
coprimality, including the known factors \(P_j,P_k\).
The two complete values then range over subpower-sized
fixed-support dictionaries and determine the third.
Otherwise two distinct support vectors are
noncollinear by primitivity. Their cross product is
\(b\boldsymbol n\) for a nonzero integer \(b\).
Congruence (SR1) gives \(D_i\mid b n_i\) for all \(i\).
The moduli are coprime, so \(Q(\boldsymbol n)\mid b\).
Rescaling the coordinates to \(x_i/H_i\) gives a
cross product of length at least
\(Q(\boldsymbol n)\|\boldsymbol n^*\|\), while each vector
has norm at most \(\sqrt3\). Packing their directions
inside their common rescaled plane proves (SR3).

### Pricing high-occupancy planes

Suppose \(H_i\le n^{a_i}\), \(D_i\ge n^{d_i}\) and
\(n=X^L\). Fix \(r>0\). On any plane branch with at least
\(n^r\) distinct support vectors and all \(n_i\ne0\),
(SR3) gives integers
\[
          w_i=\frac{Q(\boldsymbol n)n_i}{D_i},\qquad
          |w_i|\ll \frac{H_jH_k}{D_i n^r}.
\]
The normal is projectively \((D_1w_1,D_2w_2,D_3w_3)\).
Thus for each prescribed modulus triple the planes
arising from all such branches together can be listed
at normalized exponent at most
\[
                  \sum_i (a_j+a_k-d_i-r)_+.                \tag{SR4}
\]
This prices the normals without summing over the corrected
residue dictionaries. A full-rank branch cannot meet
the occupancy threshold if
\(r>\sum_i a_i-\sum_i d_i\) with a fixed margin.
Subpower branch counts and multiplicity per support vector
are allowed by making the threshold slightly smaller
than the exponent in actual original rows.

For example, with \(a_1=a_2=a_3=1/3\), shallow moduli
\(d_i\ge1/3-o(1)\) and a nearly full branch
\(r=1/3-o(1)\), the projective normal can be charged
at vanishing normalized cost, regardless of its
ordinary coefficient height. Formula (SR3) is stronger
when the normal fails to absorb much of its three
corresponding moduli. The coefficient trace is necessary
for the price (SR4), which concerns high-occupancy
branches.

## Two equal residual powers on a radical plane

Prescribe a polynomial-height pool of primes \(P\) and
a polynomial-height integral plane normal
\(\boldsymbol n\), all of whose entries are nonzero.
The latter height bound is automatic for an actual
rank-two plane of bounded support vectors by choosing
two noncollinear vectors and taking their primitive
cross product. Let \(x_i=\operatorname{rad}_{\overline
P}(T_i)\). Suppose that off \(P\) two complete letters
are exact powers of their *own residual radicals* at a
common integer order \(q\ge2\). That is, their signed
values are \(A x_1^q\) and \(B x_2^q\), with nonzero
signed \(P\)-supported coefficients \(A,B\) of polynomial
height. Fix \(A,B\) at subpower cost when only the pool
was given. The third signed value equals
\(F(x_1,x_2)=A x_1^q+B x_2^q\ne0\), and has residual
radical \(x_3\). The plane equation is
\[
                  n_1 x_1+n_2 x_2+n_3 x_3=0 .
\]

**SR5.** For fixed \(q\ge3\) the number of original
rows in this record is \(O_\varepsilon(X^\varepsilon)\),
uniformly at fixed polynomial height bounds. For
\(q=2\) the same bound holds unless
\[
                         A(-n_2)^q+B n_1^q=0.              \tag{SR5}
\]
In the exceptional case the third complete letter is
squarefree outside a prescribed enlarged pool of
polynomial prime product. Since the first two complete
values determine the third, there is no common scale
parameter in this count.

If the expression in (SR5) is nonzero, reduce the plane
equation at each prime dividing \(x_3\) outside
\(n_1n_2\). The other two radicals are units there by
primitivity. Thus \(x_1/x_2\equiv-n_2/n_1\), and the
original signed equation shows that this prime divides
the expression in (SR5). The third letter is supported
on a fixed polynomial-height pool. Fix its complete
value at subpower cost and use the fixed nonzero-level
binary Thue bound for \(F\) when \(q\ge3\), or the
uniform diagonal binary quadratic norm bound when \(q=2\).

If (SR5) is zero, put \(d=\gcd(n_1,n_2)>0\) and write
\[
          F(X_1,X_2)=L_0(X_1,X_2)G(X_1,X_2),\qquad
          L_0=(n_1X_1+n_2X_2)/d .
\]
The factor \(G\) is integral by Gauss's lemma and has
polynomial coefficient height at fixed \(q\). The
binomial \(F\) has distinct geometric roots in
characteristic zero, so \(G\) is separable and has
nonzero homogeneous resultant with \(L_0\).
Include the primes of \(dn_3\) and this resultant
in the given pool. If a prime outside that pool divides
\(G(x_1,x_2)\), the original equation puts it in the
third letter's residual radical \(x_3\). The plane then
makes it divide \(L_0(x_1,x_2)\). This is impossible
because the pair \(x_1,x_2\) is primitive and the
resultant is a unit at that prime. Thus the nonzero
value of \(G\) is supported on the enlarged fixed pool.
Fix it at subpower cost. When \(q\ge4\), use the
coefficient-uniform Thue–Mahler bound for its at least
three distinct projective factors. For \(q=3\), write
\(G=a_0 X_1^2+b_0 X_1X_2+c_0 X_2^2\). If \(a_0\ne0\),
complete the square to get a nondegenerate binary
quadratic equation, with coefficients and variables of
polynomial height. The same uniform norm bound applies.
If \(a_0=0\), use the nonzero level of
\(X_2(b_0X_1+c_0X_2)\) and the divisor bound.
For \(q=2\), the plane and the supported quotient give
\[
                 |F(x_1,x_2)|
                   =|n_3/d|\ x_3\,|G(x_1,x_2)|.
\]
Every prime outside the enlarged pool has valuation at
most one on the right, proving the asserted alternative.
These arguments keep all signs and do not assume any
archimedean comparison of the original three letters.
[Evertse, Theorem 6](https://ir.cwi.nl/pub/2463/2463D.pdf#page=14)
· [Browning–Verzobio, Lemma 2.3](https://arxiv.org/html/2608.24512v1#S2.SS1)
· [Binary norm and weighted repair recovery](ABC_REPAIRED_WEIGHTED_TARGET_DIFFERENCES.md#proof-of-recovery-and-enumeration)

One may allow the common order to grow with height. For
orders above a large cutoff depending on the desired
inner exponent allowance, the first two roots have at
most \(X^{2/q}\) choices since each is a full root of
an original letter. Orders only up to \(\log_2 X\)
are needed if at least one root exceeds one. If both
roots are one, the fixed coefficients already recover
the complete row.

### Consequence for exceptional quadratic records

At every balanced endpoint there is an absolute positive
radical mass in a fixed compact positive scaled-degree
band on each original letter. Choose a fixed
\(\gamma>0\) smaller than the resulting normalized
lower bound, allowing the profile refinements, and
consider original rows with at least \(\gamma L\log X\)
of logarithmic radical on valuation at least two in
the third letter. In the exceptional quadratic case
all this mass is carried by the enlarged fixed
prime pool. Pricing that polynomial-height pool at
zero in (EB5) gives exponent at most
\((1-\gamma)L/3+\eta L+\varepsilon\) in the original
shell with raw total radical at most \(X^L\).
The uniform small-budget threshold and the allowance
hold at the prescribed fixed polynomial height bounds.
Hence modulus, normal, and coefficient records of
combined normalized price less than \(\gamma/3\) by
a fixed margin yield a strict saving in this sector.
One can use (SR4) for the normal portion of this price
when the branch occupancies and shallow moduli satisfy
its inequalities. For \(q\ge3\), subpower completion
applies once the plane and pool record are fixed.
[Positive scaled degrees on original letters](ABC_SHORT_HEIGHT_AND_REAL_GAP.md#positive-scales-in-every-endpoint-letter)
· [Polynomial-height support pools and effective radical costs](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)

The exact power assumption in (SR5) concerns the
entire residual part of each of the two letters.
It also applies after other valuation layers of
those letters have been fixed and put in the known
pool, at their actual record cost. It does not follow
from a high least-repair order or from a divisible
valuation class with varying multiples of \(q\).
The arguments leave the general diffuse radical
planes open.

## Sharp absorption by the normal outside the small-budget sector

There are actual primitive rows with one common
triple of corrected residue records and with the normal
absorbing all three corresponding moduli. They attain
the angular bound (SR3) in order of magnitude. Their
residual parts are squarefree and their original radical
budget is larger than the numerical height.

Fix \(1/2<\kappa<1\). Let \(t\) be large and odd, and put
\[
       D_1=t^2,\quad D_2=(t+1)^2,\quad D_3=(t+2)^2,
       \qquad H=t^{1+\kappa},\quad B=t^\kappa .
\]
There are choices of \(t\to\infty\) for which the plane
\[
                     D_1x_1+D_2x_2=D_3x_3                 \tag{SR6}
\]
contains \(\gg B^2\) vectors with \(x_i\asymp H\)
positive, squarefree, pairwise coprime and each
coprime to \(D_1D_2D_3\). Here and below the constants
may depend on \(\kappa\).

To verify this, the two vectors
\[
           (-3t-4,\,4t,\,t),\qquad
           (4t+5,\,-4t+4,\,1)
\]
are in the plane and their cross product is
\(4(D_1,D_2,-D_3)\). Since the moduli are pairwise
coprime, the lattice of all integer vectors in the
plane has primitive normal \((D_1,D_2,-D_3)\).
The displayed vectors span a sublattice of index four.
There is a basis of the full lattice with lengths
\(O(t)\), since any enlargement of index four lies in
the lattice generated by the displayed vectors divided
by four and there are only finitely many corresponding
basis changes. Its cross product has length
\(\asymp t^2\). Thus both singular values of the basis
matrix have size \(\asymp t\). The positive part of the
plane contains scaled boxes bounded away from the
coordinate walls because \(D_i/t^2\to1\). In the two
basis parameters \(u,v\) choose a rectangle with sides
\(\asymp B\) inside a region \(|u|,|v|\ll B\) on which
\(x_i\asymp H\). The coefficients of each \(x_i\) in
\(u,v\) are coprime, since their gcd divides the other
two components of the primitive cross product. At
least one of the two coefficients has size \(\gg t\)
by the sizes of the minors.

For the sieve, \(t\) can be chosen in \([T,2T]\) so
that it is odd and \(t\equiv1\bmod p\) for all primes
\(5\le p\le c\log T\), with a sufficiently small
absolute constant \(c>0\). The Chinese remainder
modulus is \(o(T)\), by the usual bound on the product
of primes up to \(c\log T\). None of these primes
divides a \(D_i\). For each prime at most a sufficiently
large fixed \(P_0\), choose parameter residues with
all three \(x_i\) units. Such residues exist: at \(2,3\)
one plane coefficient is zero and the other two are
units. At the other primes the three coordinate lines
cannot cover the plane. Fix the corresponding
parameter congruences, leaving \(\gg B^2\) points.

For primes \(p>P_0\) dividing a \(D_i\), exclude
\(p\mid x_j\) for any \(j\). In the fixed parameter
class modulo its bounded modulus, each exclusion costs
\(O(B'^2/p+B')\), where \(B'\asymp B\) denotes the
new parameter length. Since these \(p\gg\log t\) and
\(\omega(D_1D_2D_3)=O(\log t/\log\log t)\), the total
is \(o(B'^2)\). At other primes two coordinates
divisible by \(p\) force \(p\mid u,v\), by the
nonzero unit minors. Excluding this costs
\(O(B'^2/P_0)+o(B'^2)\) in total, by the bound
\(O((B'/p+1)^2)\) for \(p\ll B\).
For \(P_0<p\le B'\), excluding \(p^2\mid x_i\) costs
\(O(B'^2/p^2+B')\) per prime, giving
\(O(B'^2/P_0)+o(B'^2)\) by the prime counting bound.
For \(p>B'\), each integer value of \(x_i\) has only
\(O(1)\) representations in the rectangle: the
nonzero step in a fiber of the primitive linear form
has size \(\gg t\), whereas \(B=o(t)\). Summing multiples
of \(p^2\) up to \(O(H)\) thus gives
\(O(H/B')=O(t)=o(B'^2)\), since \(\kappa>1/2\).
The leading sieve constants in the shortened parameter
box are absolute. Choosing \(P_0\) sufficiently large
leaves \(\gg B^2\) vectors, as required.

Now take the actual letters \(D_i x_i\). They are
positive, primitive and satisfy the original sum.
The factored radicals are exactly the \(x_i\).
Modulo \(D_i\), the complementary ratio \(x_j/x_k\)
is a fixed unit by (SR6). Therefore all three
corrected residues are fixed by (TL6), at any fixed
repair orders. The primitive normal is exactly
\((D_1,D_2,-D_3)\) up to sign, and the occupancy
\(B^2=H^2/t^2\) has the order allowed by (SR3) when
\(Q(\boldsymbol n)=1\). The original height is
\(\asymp t^2 H\) and the radical of the three complete
letters is at least a constant times \(H^3\).
Accordingly this example tests the simultaneous
residue packing and normal absorption. Ruling out
such absorption in the small-budget problem requires
the large valuations on the original residual primes.
