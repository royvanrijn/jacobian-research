# Primitive secants and triangles with only two radical caps

There are families of actual positive primitive additive triples
with two balanced radical caps, diffuse valuations on both capped
letters, and prescribed target multiplicity up to the half-budget
scale, for which almost all primitive target secants have
unbounded normalized radical. There is also a fixed-target fibre
with two radical caps in which almost all primitive difference
triangles have this property on all three sides. In both examples
the omitted summand radical is large on almost every row.

Use \(n=X^L\), with \(L>0\) fixed before \(X\to\infty\), and
then let \(L\to0\). Subscripts on \(o_X(1)\) refer to this inner
limit. The constructions below exist for all sufficiently large
\(X\) at each sufficiently small fixed \(L\). In particular the
two limits do not determine a single joint choice of \(L\) and \(X\).

## Diffuse power dictionaries

Take \(e=2^t\), \(t=\lfloor \tfrac12\log_2(1/L)\rfloor\), and
\(m=\lfloor L^{-1/4}\rfloor\). Write \(r=t-1\), the number of
nonlinear cyclotomic factors of \(z^e-x^e\).
For sufficiently small \(L\), these parameters satisfy
\(r\to\infty\), \(m\to\infty\), and \(eL\to0\) in the budget limit.
Set
\[
 s_i=eL k_i,\qquad
 k_i=\left\lfloor\frac{2+(2i-1)/m}{eL}\right\rfloor,
 \qquad
 H_L=\frac{1}{3m}\sum_{i=1}^m s_i
 \qquad (1\le i\le m).
\]
The integers \(k_i\) are positive, \(H_L\le1\), and \(H_L\to1\).

Put \(R=n^{1/(3m)}\), \(Q=\lfloor R/(\log R)^2\rfloor\).
In each of \(m\) disjoint subintervals of \([R/8,R/4]\) choose
a palette of \(Q\) primes. Each interval can be chosen to have
length a positive fixed multiple of \(R\) at the fixed \(L\).
The prime number theorem gives the required palettes for large
\(X\). Make a second set of \(m\) disjoint palettes of the
same size within \([R/2,R]\). A word in either dictionary
is a positive integer
\[
                       v=\prod_{i=1}^m p_i^{e k_i}
\]
with one prime from each palette of its dictionary. Denote
the two full dictionaries by \(\mathcal A_0,\mathcal C_0\).
The disjointness of the palettes gives unique words. Each
dictionary has size \(Q^m=n^{1/3+o_X(1)}\), every word has
radical at most \(n^{1/3}\), and all heights are
\(\asymp_L X^{H_L}\). The heights are at most \(X\).
Every \(c\in\mathcal C_0\) exceeds \(2a\) for every
\(a\in\mathcal A_0\) because each \(c\)-prime is at least twice
the corresponding \(a\)-prime. Their supports are disjoint.

For every word the scaled-degree measure is
\[
       \sum_{i=1}^m \log_n p_i\,\delta_{s_i}
       =\frac{1}{3m}\sum_{i=1}^m\delta_{s_i}+o_X(1)
\]
in total variation at fixed \(L\). As \(L\to0\) this tends
weakly to the wholly nonatomic measure of density \(1/6\)
on \([2,4]\). All statements about height, radical, and degree
measure continue to hold on subdictionaries.

### Few ordinary differences of bounded radical

For either full dictionary \(\mathcal D_0\), and any fixed
\(F>0\) at the inner height limit, the number of distinct
ordered pairs \(v,v'\) with
\(\operatorname{rad}(|v-v'|)\le n^F\) is
\[
                         \ll_{L,F} n^{1/3+F/(2r)+o_X(1)} .
 \tag{DO1}
\]
The inner constants may also depend on \(F\).
To see this, write the two words as \(v=x^e,v'=z^e\),
put \(G=\gcd(x,z)\), and set \(x_1=x/G,\ z_1=z/G\).
The \(r\) factors
\[
                   x_1^{2^j}+z_1^{2^j},
                   \qquad 1\le j\le r,
\]
all divide \(|v-v'|\) and their prime supports are disjoint
outside 2. At an odd prime dividing one of the factors,
\(x_1/z_1\) has \(2^j\)-th power \(-1\). Raising to a higher
power of two shows that this prime cannot divide another
factor. The radical cap therefore makes at least one factor
have radical at most \(2n^{F/r}\).

For each such factor there are at most \(n^{F/r+o_X(1)}\)
possible values, by enumeration of the radical followed by
uniform enumeration on its fixed prime support. Every fixed
value has only \(X^{o_X(1)}\) primitive pairs \(x_1,z_1\) of
polynomial height, by cyclotomic norm recovery. Indeed the
value is the norm of \(x_1-\zeta z_1\) for a primitive
\(2^{j+1}\)-th root of unity \(\zeta\) in a fixed cyclotomic
field at this \(L\). The ideal of the generator has
divisor-many possibilities. The absolute values at all
embeddings are between \(X^{-O_L(1)}\) and \(X^{O_L(1)}\),
using the upper bound for each conjugate and the nonzero
integral norm. Generators of a fixed ideal differ by units
whose logarithmic vectors in this range have polylogarithmically
many possibilities up to the finite kernel. The generator
recovers \(x_1,z_1\) since \(\zeta\) is irrational.
[Cyclotomic factor values and uniform fixed-field recovery](ABC_MINIMAL_REPAIR_CYCLOTOMIC_SPLIT.md#primitive-cyclotomic-factors-and-their-fixed-value-counts)
· [Uniform fixed-support enumeration](ABC_EXCEPTIONAL_SET_AND_POWERFUL_TRIPLES.md#a-radical-count-contradiction-to-the-claimed-smooth-number-lower-bound)

Stratify the original ordered pairs by the number \(s\) of
slots where their prime choices coincide. The slot supports
are disjoint, so the missing slots in \(x_1,z_1\) are precisely
these \(s\) slots. The original words are recovered from
\(x_1,z_1\) and the shared primes, giving at most \(Q^s\) lifts.
The trivial number of pairs of this type is at most
\(\binom{m}{s}Q^{2m-s}\). Combining the two counts bounds this
class by \(n^{o_X(1)}\) times
\[
             \min\{n^{2/3-s/(3m)},\
                           n^{F/r+s/(3m)}\}
                         \le n^{1/3+F/(2r)}.
\]
Summing the fixed number of classes proves (DO1).

### Thin dictionaries with small unit classes

For every modulus \(n\le d\le X\), every unit class contains
at most \(n^{1/6+o_X(1)}\) words of either full dictionary,
uniformly in \(d\) and the class. In each slot \(i\), restrict
to unit primes and form the character sum
\(S_i(\chi)=\sum_{p_i}\chi(p_i)^{e k_i}\).
In the \(2m\)-th moment, fix one \(m\)-tuple of primes.
The product of the other \(m\)-tuple lies in one coset
of the kernel of the \(e k_i\)-th power map on the units
modulo \(d\). This kernel has size \(X^{o_X(1)}\), since
\(e k_i\) is fixed at the inner limit and the unit groups
of prime powers give a bound \((e k_i)^{2\omega(d)}\).
Every allowed integer product has size at most
\(R^m=n^{1/3}<d\). Thus there are \(X^{o_X(1)}\) possible
products and \(X^{o_X(1)}\) factorizations per product.
The normalized character moment is at most
\(Q^m n^{o_X(1)}\). Hölder on the \(m\) slots bounds
\(\mathbb E_\chi|\prod_i S_i|^2\) by \(Q^m n^{o_X(1)}\).
Character orthogonality and Cauchy--Schwarz give the
claimed unit-class bound.
[Single-letter multiplicative residue count](ABC_BOUNDED_DEGREE_CURVES_AND_SECANT_GCD.md#one-complete-letter-in-a-unit-residue)

Choose \(B_L=\sqrt r\), \(F_L=B_L+4/3\), and
\(\epsilon_L=\min(1/100,r^{-1/4})\). There exists a
subdictionary \(\mathcal A\subseteq\mathcal A_0\) with the
following properties:
\[
\begin{split}
 |\mathcal A|&=n^{1/6+o_X(1)},\\
 |\mathcal A\cap (z+d\mathbb Z)|&\le 2n^{\epsilon_L}
        \quad(n\le d\le X,\ (z,d)=1),\\
 \#\{(a,a')\in\mathcal A^2:a\ne a',\
             \operatorname{rad}(|a-a'|)\le n^{F_L}\}
        &\le n^{F_L/(2r)+\epsilon_L}.
\end{split}                                                    \tag{DO2}
\]
Retain the words of \(\mathcal A_0\) independently with probability
\(n^{-1/6}\). The size has the stated exponent with probability
tending to one, by the variance bound. Each unit cell has
expected size \(n^{o_X(1)}\), uniformly. With
\(j=\lceil n^{\epsilon_L}\rceil\), its probability of having at
least \(j\) selected words is at most
\((\mathrm e\mu/j)^j\), where \(\mu\) is its expectation.
This tends to zero faster than any inverse power of \(X\),
even after taking the union over the at most \(X^2\)
modulus and class pairs. For the last claim the expectation
in (DO1) is multiplied by \(n^{-1/3}\). Markov's inequality
gives the displayed bound with probability tending to one,
since \(F_L,\epsilon_L\) are fixed in the inner limit and
the inner exponent error can be taken below \(\epsilon_L/2\).
These choices prove simultaneous existence.

## Target secants at prescribed multiplicity

Fix \(0\le u\le1/2\), including both endpoints, and retain
each target of \(\mathcal C_0\) independently with probability
\(n^{-u/3}\). There exists a resulting set \(\mathcal C\)
of size \(n^{(1-u)/3+o_X(1)}\) in which the number of distinct
ordinary target pairs of radical at most \(n^{F_L}\) is at most
\[
                 n^{1/3-2u/3+F_L/(2r)+\epsilon_L}.
 \tag{DO3}
\]
Indeed the size has the stated exponent with probability
tending to one, and (DO1) bounds the expectation in (DO3)
with the factor \(n^{-2u/3}\). Use Markov with the displayed
exponent allowance. When \(u=0\) the full dictionary also
works directly. Restrict to a single dyadic \(c\)-shell
containing a maximum number of chosen targets. Their original
heights are \(\asymp_L X^{H_L}\), so there are only a bounded
number of such shells at fixed \(L\). This preserves the
size exponent and (DO3). Its endpoint \(Y\) has
\(Y\asymp_L X^{H_L}\), and \(n\asymp_L Y^{L/H_L}\) with
\(L/H_L\to0\).

For each \(c\in\mathcal C\), choose independently a uniform
subset of \(\mathcal A\) of size
\[
                 J=\min\{|\mathcal A|,\lfloor n^{u/3}\rfloor\}.
\]
Include \((a,c-a,c)\) for its selected \(a\)'s.
The primes supporting \(a\) and \(c\) are disjoint, so every
row is positive and primitive. Each complete target has
exactly \(J=n^{u/3+o_X(1)}\) rows. There are
\(n^{1/3+o_X(1)}\) distinct rows in total. The actual
summand dictionary on \(a\) equals \(\mathcal A\) with
probability tending to one: the probability of missing
a prescribed \(a\) is at most
\(\exp(-J|\mathcal C|/|\mathcal A|)=\exp(-n^{1/6+o_X(1)})\).
These assertions allow \(J=1\) and \(J=|\mathcal A|\).
All heights, including \(c-a\), are
\(\asymp_L X^{H_L}\); the two tracked radicals are bounded
by \(n^{1/3}\) and both tracked degree laws are diffuse.

For a row pair with different targets define
\[
             g=\gcd(|a-a'|,|c-c'|),\qquad
             u_c=(c-c')/g .
\]
This is the common gcd of the affine secant and its
primitive target component, since the third coordinate is
the difference of the first two in these \((a,c)\) coordinates.
There is a choice of the random subsets for which
\[
 \#\{ \text{ordered row pairs with equal targets or }
                 (c\ne c'\text{ and }\operatorname{rad}|u_c|\le n^{B_L})\}
                         \le n^{1/2+o_{L\to0}(1)+o_X(1)} .
 \tag{DO4}
\]
The equal-target pairs are at most
\(|\mathcal C|J^2=n^{(1+u)/3+o_X(1)}\le n^{1/2+o_X(1)}\).
If ordinary target differences have radical at most \(n^{F_L}\),
(DO3) gives at most
\(n^{1/3+F_L/(2r)+\epsilon_L+o_X(1)}\) lifted row pairs.

For the remaining pairs it suffices to count those with
\(\operatorname{rad}(g)>n^{4/3}\). Fix distinct \(c,c'\) and
fix \(a'\in\mathcal A\). Because
\(\operatorname{rad}(a')\le n^{1/3}\), this event implies
there is a divisor \(d\mid |c-c'|\) with \(d>n\),
\((d,a')=1\), and \(a\equiv a'\bmod d\). One can take the
squarefree part of \(g\) supported outside \(a'\). The number
of possible divisors is \(n^{o_X(1)}\), and each unit class
has at most \(2n^{\epsilon_L}\) possibilities for \(a\) by
(DO2). For these two distinct targets any prescribed pair
\(a,a'\) is selected with probability
\((J/|\mathcal A|)^2\). The total expected number of such
bad lifted pairs is therefore bounded by
\[
      |\mathcal C|^2 J^2\, n^{\epsilon_L+o_X(1)}
                           /|\mathcal A|
                       =n^{1/2+\epsilon_L+o_X(1)} .
\]
Markov bounds it by \(n^{1/2+2\epsilon_L+o_X(1)}\) with
probability tending to one, simultaneously with activation of
the actual \(a\) dictionary. For all other pairs,
\(\operatorname{rad}(|c-c'|)>n^{B_L+4/3}\) and
\(\operatorname{rad}(g)\le n^{4/3}\) give
\(\operatorname{rad}|u_c|>n^{B_L}\). Since
\(B_L\to\infty\), \(\epsilon_L\to0\), and \(F_L/r\to0\),
this proves (DO4).

## Difference triangles in one fixed-target fibre

Choose \(P\) prime in \([\tfrac12 n^{1/3},n^{1/3}]\)
and disjoint from all palettes of \(\mathcal A\). For
\(m>1\) this disjointness holds automatically for large \(X\).
Set
\[
                T=P^{e(\lceil 3/(eL)\rceil+1)} .
\]
Then \(T\) is an exact \(e\)-th power,
\(\operatorname{rad}(T)\le n^{1/3}\), and
\(\log_X T=1+O(eL)+o_X(1)\). Moreover \(T>2X\) for
sufficiently large \(X\) at the fixed \(L\). Take the rows
\((x,T-x,T)\) for all \(x\in\mathcal A\). They are positive
and primitive, and form a single complete-target fibre
of size \(n^{1/6+o_X(1)}\). The free letter \(x\) has the
same bounded radical and diffuse degree law as before.
The free letter \(T-x\) is tall and the fixed target and
\(x\) have numerical height exponents tending to one. In the
shell with endpoint \(T\), the relative budget unit is
\(L/H_{T,L}\), where \(H_{T,L}=L e(\lceil 3/(eL)\rceil+1)/3\).
Indeed \(T\asymp_L X^{H_{T,L}}\) and \(H_{T,L}\to1\).

For \(x_1<x_2<x_3\) from this fibre put
\[
 g=\gcd(x_2-x_1,x_3-x_2),\qquad D_{ij}=(x_j-x_i)/g .
\]
Outside at most \(n^{1/3+o_{L\to0}(1)+o_X(1)}\)
triangles, one has
\[
            \operatorname{rad}(D_{12}),\
            \operatorname{rad}(D_{23}),\
            \operatorname{rad}(D_{13})>n^{B_L}.
 \tag{DO5}
\]
A triangle with an ordinary difference of radical at most
\(n^{F_L}\) is counted by (DO2) after choosing its remaining
vertex. This gives at most
\(O(n^{1/6+F_L/(2r)+\epsilon_L+o_X(1)})\) triangles.
If \(\operatorname{rad}(g)>n^{4/3}\), fix \(x_1,x_2\).
There is a divisor \(d>n\) of \(x_2-x_1\) prime to \(x_2\)
with \(x_3\equiv x_2\bmod d\), by removing the primes
of \(x_2\) from the radical of \(g\). The divisor bound and
(DO2) give at most \(n^{\epsilon_L+o_X(1)}\) choices of
\(x_3\) per fixed pair. These two exceptional counts prove
(DO5) by the same division inequality as for secants.
They are power smaller than the full \(n^{1/2+o_X(1)}\)
triangle scale.

In both constructions the untracked summand is itself
a difference of two coprime exact \(e\)-th powers.
If its radical is at most \(n^{B_L}\), at least one of
its \(r\) nonlinear cyclotomic factors has radical at
most \(2n^{B_L/r}\). Fixed-value norm recovery as in
(DO1), now without a shared-slot lift because the two
power bases are coprime, bounds the number of such
original rows by \(n^{B_L/r+o_X(1)}\). Thus adding even
a fixed normalized bound to the omitted radical discards
almost all original rows in these examples. They address
the necessity of that third constraint in any deduction of
affordable primitive secants or triangles. The question with
all three radical caps remains open.
[Primitive target-coordinate secants](ABC_UNRECORDED_COFACTOR_SECANT_TRANSFER.md#primitive-secant-coordinates-with-an-arbitrary-scaling-gcd)
· [Primitive triangles inside fixed fibres](ABC_FIXED_FIBRE_TRIANGLE_RADICALS.md#primitive-triangles-of-three-values-on-the-fibre)
