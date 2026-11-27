# Coherent prime triangles from arithmetic progressions

Prime progressions give genuine coherent triangles with the same
coefficient on all three edges. Two fixed split residue tests force
the ideal-class differences to span a two-dimensional subspace.
Counting quadratic twists of a fixed genus-two curve controls how
often the same field can arise from different progressions.

## A prime progression and its ideal classes

**Proposition PP1.** Put \(M=6\cdot73\cdot613\). Suppose
\[
          (p,h)\equiv(23,6)\pmod M,\qquad 0<h\le p/10
\]
and \(p_1=p-h,\ p_2=p,\ p_3=p+h\) are primes greater than \(7\).
Define
\[
 \begin{split}
 E=F(p,h)&=3p^6-24p^4h^2-24p^2h^4-4h^6,\\
 E&=d\delta^2,\qquad d>0\text{ squarefree},\quad \delta>0.
 \end{split}                                                     \tag{1}
\]
There is a coherent triangle of oriented split primes over \(p_i\)
in \(\mathbb Q(\sqrt{-d})\) with all three signed coefficients
\(v_{ij}=\delta\) for \(i<j\). The two consecutive edge differences
\([\overline{\mathfrak p}_1\mathfrak p_2]\) and
\([\overline{\mathfrak p}_2\mathfrak p_3]\) are independent in
\(\mathrm{Cl}(\mathbb Q(\sqrt{-d}))[3]\).
The two associated integral points on \(y^2=x^3-E/4\), in the
integral model \(Y^2+Y=x^3-(E+1)/4\), are independent modulo
torsion. Here the change of coordinate is \(Y=y-1/2\).

Write \(A=p_1^3,\ B=p_2^3,\ C=p_3^3\). The symmetric expression
\[
             2(AB+AC+BC)-(A^2+B^2+C^2)=4AC-(A+C-B)^2
\]
is (1). It is positive and \(\asymp p^6\) under the bound on \(h/p\).
It is prime to the three \(p_i\). At \(p_1\), \(C=8B\bmod p_1\),
giving \(E=-49B^2\bmod p_1\). At \(p_3\) the same argument gives
\(-49B^2\bmod p_3\). At \(p_2\), \(A=-C\bmod p_2\), giving
\(E=-4C^2\bmod p_2\). These residues are nonzero.
The [three-coefficient reconstruction](ABC_COHERENT_NORM_PENCILS.md#exact-reconstruction-from-three-signed-coefficients)
with normalized coefficients \((1,1,1)\) gives traces
\[
               u_{12}=A+B-C,\quad u_{23}=B+C-A,\quad u_{13}=B-A-C.
\]
In particular for \(\alpha_{ij}=(u_{ij}+\sqrt{-E})/2\) one has
\[
   (\alpha_{ij})=(\overline{\mathfrak p}_i\mathfrak p_j)^3,\qquad
   \alpha_{12}\alpha_{23}=B\alpha_{13}.
\]
The three \(u_{ij}\) are odd, and \(E\equiv3\bmod4\), so the
corresponding \(Y\)-coordinates are integral. The projective
identity underlying \(F\) is the symmetric three-representation
identity that [Craig, Section 2](https://www.i-repository.net/contents/osakacu/sugaku/111F0000002-01402-7.pdf)
uses to start his higher-rank construction. See also the exact
identities in [Gillibert–Levin, Section 2.2](https://arxiv.org/pdf/1803.05800).

For the independence test, use the following residues. In the table
\(r^2=-E\bmod q\) and the last two columns are
\(((u_{ij}+r)/2)^{(q-1)/3}\bmod q\).
The modulus in the proposition fixes all these residues by the
values for the prime progression \((17,23,29)\).
\[
 \begin{array}{c|r|r|r}
 q&r&12&23\\ \hline
 73&21&64&1\\
 613&86&1&65
 \end{array}                                                     \tag{2}
\]
For an integer check, this progression has \(E=185684403\) and
\((u_{12},u_{23},u_{13})=(-7309,31643,-17135)\).
Both \(q\) split in \(\mathbb Q(\sqrt{-d})\), because
\(q\nmid E\) and we can send \(\sqrt{-d}\) to \(r/\delta\bmod q\).
Each nontrivial entry in (2) has order three. If an ideal-class
relation with exponents \(e_1,e_2\bmod3\) held, then
\(\alpha_{12}^{e_1}\alpha_{23}^{e_2}\) would be the cube of a
field element times a unit. The product is a unit at the two
selected prime ideals, so the field element can be reduced in
their residue fields. Every imaginary quadratic unit has
cube reduction at these primes: their orders divide \(2,4\) or
\(6\), and \(q\equiv1\bmod9\). Thus (2) forces
\(e_1=e_2=0\bmod3\). This argument also excludes \(d=3\), whose
class group is trivial. The rank conclusion follows from
[CR1](QUADRATIC_NORM_COLOUR_CLASS_RANK.md#from-one-colour-to-one-class-subspace),
with \(D=-d,w=\delta,b=1\).

The same check applies for each progression meeting the finite
congruences even when \(E\) has a large square factor. The
progression used in (2) has \(h/p>1/10\) and serves just to fix
the local residues. In its case \(E>0\) as well, so the coefficient
identity gives a genuine witness.

## A uniform bound on rational progressions per field

**Proposition PP2.** For each positive squarefree \(d\), the number
of coprime positive pairs \((p,h)\) with \(F(p,h)=d\delta^2\) for some
positive integer \(\delta\) is at most \(C^{1+\omega(d)}\) with an
absolute constant \(C\). There is no primality or height restriction
in this multiplicity bound.

The discriminant of \(F(t,1)=3t^6-24t^4-24t^2-4\) is
\(2^{16}3^7 7^6\), by the resultant with its derivative. The pairs
give rational points \((t,z)=(p/h,\delta/h^3)\) on the smooth genus-two
model of \(d z^2=F(t,1)\). Distinct coprime positive pairs have distinct
\(t\). Let \(J_d\) be its Jacobian. Take the fixed splitting field
\(L\) of \(F(t,1)\). All of \(J_d[2]\) is rational over \(L\) for every
\(d\), since the six branch points give its two-torsion and a
quadratic twist acts trivially on it. Away from primes of \(L\)
above \(2,3,7\) and the divisors of \(d\), the curve and Jacobian
have good reduction. At these good odd places multiplication by
\(2\) on the abelian scheme is finite étale. A local rational point
extends over the valuation ring, and its division torsor is
therefore unramified. Consequently the global Kummer map puts
\(J_d(L)/2J_d(L)\) into four copies of the group of square classes
in \(L\) with even valuations away from a set \(S_d\) of
\(O_L(1+\omega(d))\) finite primes.

This last square-class group fits in an exact sequence with the
\(S_d\)-units modulo squares and the two-torsion of the
\(S_d\)-ideal class group. The units contribute \(O_L(1+|S_d|)\)
dimensions and the ideal class group is a quotient of the fixed
class group of \(L\). It follows that
\(\operatorname{rank}J_d(\mathbb Q)\le\operatorname{rank}J_d(L)
=O(1+\omega(d))\).
If the curve has one rational point, use it to embed the curve
in \(J_d\) as a translate of a theta divisor of degree \(2\) in
the principal polarization. This genus-two curve contains no
positive-dimensional abelian translate.
[Gao–Ge–Kühne, Theorem 1.1](https://pmihes.centre-mersenne.org/item/10.5802/pmihes.26.pdf)
then bounds its rational points by
\(c(2,2)^{1+\operatorname{rank}J_d(\mathbb Q)}\). This proves the
multiplicity assertion, which is also immediate when the curve
has no rational point.

## Quantitative count of genuine fields

**Theorem PP3.** There is an absolute \(B_0>0\) such that for large
\(X\), at least \(\gg X^{1/3}/(\log X)^{B_0}\) distinct imaginary
quadratic discriminants \(-d\) with \(d\le X\) admit the coherent
prime triangles of PP1. Their three primes occur in a common
dyadic window at scale \(X^{1/6}\), and their coefficient
\(E=d\delta^2\) is comparable to \(X\). In each field these triangles
carry two independent three-torsion ideal classes and two
independent points on the same Mordell curve. Since squarefreeness
of \(E\) is unnecessary for this count, \(d\) can be smaller than
the scale \(X\).

Use \(4Z\le p\le5Z\), \(Z/10\le h\le Z/5\) and the congruences
of PP1. All primes are in \([3Z,6Z)\). Writing \(p=Ma+23\),
\(h=Mb+6\) gives affine forms in \(a,b\) with nonparallel
linear parts \(M(1,-1),M(1,0),M(1,1)\). Their complexity is one.
At primes dividing \(M\), their three constant residues
\(17,23,29\) are units. At any remaining prime, we can set
\(p=1,h=0\) modulo that prime to make all three values units.
The positive real region has area \(\asymp Z^2\) and the local
density at every prime outside \(2M\) is \(1+O(q^{-2})\), since
the three excluded lines have distinct slopes. The
[Green–Tao Main Theorem and Corollary 1.7](https://arxiv.org/pdf/math/0606088)
give the asymptotic for complexity at most two unconditionally.
In this rectangle it gives \(\gg Z^2/(\log Z)^3\) prime
progressions. In each one \(\gcd(p,h)=1\), since \(p\) is prime
and \(0<h<p\).

To keep a logarithmic loss when applying PP2, discard the rare
pairs with too many prime factors of \(F\). Put
\(k=\lceil100\log\log Z\rceil\) and \(Q=Z^{1/(2k)}\).
Since \(F\ll Z^6\), at most \(12k+O(k/\log Z)\) distinct prime
factors can exceed \(Q\). For each prime \(q\), the fraction of
residue pairs modulo \(q\) on \(F=0\) is at most \(6/q+1/q^2\).
Indeed the binary form has at most six projective roots modulo
every \(q\), since the coefficients have no common prime
factor. For each product of \(k\) distinct primes at most \(Q\),
the modulus is at most \(\sqrt Z\). Every residue pair for that
modulus occurs \(O(Z^2/\text{modulus}^2)\) times in the
rectangle. By CRT, summing over those products bounds the number
of pairs with at least \(k\) such prime factors by
\[
          \ll Z^2 \frac{(6\log\log Z+O(1))^k}{k!}
          \ll \frac{Z^2}{(\log Z)^4}.
\]
The last inequality follows by \(k!\ge(k/e)^k\) and the choice
of \(100\). After discarding these pairs there are still
\(\gg Z^2/(\log Z)^3\) prime progressions, all with
\(\omega(d)\le \omega(F)=O(\log\log Z)\).
Proposition PP2 now gives \((\log Z)^{O(1)}\) as a bound on the
number of kept progressions in each field. Since
\(E\asymp Z^6\) and \(E\equiv3\bmod4\), take
\(Z=cX^{1/6}\) with fixed small \(c>0\). Then \(d\le E\le X\)
and \(-d\) is the fundamental discriminant. This proves the
assertion.

This family concerns the existence of true, same-colour witnesses,
including simultaneous variation of the fields. The normalized
coefficient triple is always \((1,1,1)\), whereas the occupied
rank-one counting problem concerns coefficients in growing
ranges. A sharper lower bound on how often \(d\asymp E\) here
would require information about the square part of the binary
sextic at prime progressions, beyond the squarefree sieve for
unrestricted pairs.
