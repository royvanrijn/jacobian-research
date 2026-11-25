# A finite digit reduction for first \(p\)-adic coprimality certificates

At any fixed terminal \(N\) and prime \(p\), all factorial coefficient
valuations and first units relevant to Newton and residual certificates
are determined by a finite carry graph. Moreover, every profile of
arbitrary dimension has exactly the first \(p\)-adic data of a single
exponent row, after separate constant normalizations of the two
polynomials and a common \(p\)-adic scaling of their variable. The
single row can depend on the prime.

The carry graph also gives a finite integer-arithmetic formulation
of which first-layer certificates can fail. The formulation includes
all ordinary exponent profiles simultaneously, without bounding the
dimension or the exponents.

## Digit weights

Let \(j=N-1,N\), and index the moment coefficients by \(0\le k\le j\).
For an exponent row \((a_i,b_i)\) write
\(x_{i,j,k}=(j-k)a_i+kb_i\). Let \(s_p(x)\) denote the base-\(p\)
digit sum and \(\Phi_p(x)\) the product of the factorials of its
base-\(p\) digits in \(\mathbb F_p^\times\).
Write \(\Delta=\sum_i(b_i-a_i)\), and put
\[
 D_{j,k}=\sum_i\big(s_p(x_{i,j,0})-s_p(x_{i,j,k})\big),\qquad
 U_{j,k}=\prod_i\Phi_p(x_{i,j,k})/\Phi_p(x_{i,j,0}).
\]
With \(\delta\) the representative of \(\Delta\bmod (p-1)\)
in \(\{0,\ldots,p-2\}\), put \(q=(\Delta-\delta)/(p-1)\).
Normalize the ordinary moment polynomial by its constant term
and substitute \((-p)^{-q}t\). The resulting coefficient
at \(k\) has valuation and unit
\[
 \begin{aligned}
 V_{j,k}&=v_p\binom jk+W_{j,k},\\
 E_{j,k}&=\operatorname{unit}_p\binom jk\;(-1)^{W_{j,k}} U_{j,k},
 \qquad
 W_{j,k}=\frac{k\delta+D_{j,k}}{p-1}.
 \end{aligned}                                                       \tag{1}
\]
Here \(W_{j,k}\in\mathbb Z\), and \(\operatorname{unit}_p\) means
strip the full power of \(p\) and reduce modulo \(p\).
Legendre's formula gives the valuation difference in each row,
and Wilson's product formula gives
\(\operatorname{unit}_p(x!)=(-1)^{v_p(x!)}\Phi_p(x)\).
The common substitution removes \(kq\) from the valuation and
removes the associated sign. These identities give (1).

## Carry graph and reduction to one row

Read the pairs of digits \(A,B\in\{0,\ldots,p-1\}\)
from least to most significant. At one step, and for each
\(j,k\), update the carry \(c_{j,k}\) and output the digit
\(z_{j,k}\) by
\[
 (j-k)A+kB+c_{j,k}=p\,c'_{j,k}+z_{j,k},\qquad 0\le z_{j,k}<p.
 \tag{2}
\]
All carries start at zero. They stay in \(\{0,\ldots,j-1\}\),
since \(j(p-1)+j-1=jp-1\).
Append enough zero input digits to make all carries zero again.
For each digit step record the integer differences
\(z_{j,0}-z_{j,k}\), the multiplicative residues
\(z_{j,k}!/z_{j,0}!\), and the additive residue
\(B-A\bmod(p-1)\). Over an entire row these accumulate to its
contributions to \(D,U,\delta\).

Concatenate the padded digit words of any finite list of rows,
starting with the least-significant block. The carry vector
is zero at every boundary, so the resulting single pair of
integers has exactly the same \(D,U,\delta\) as the list.
Digit positions do not affect the last group coordinate,
because \(p\equiv1\bmod(p-1)\). If the original columns are
distinct, at least one input digit pair differs and the two
packed integers differ. Formula (1) proves the equivalence.
Packing at a different prime gives a different univariate pair
in general, so this reduction concerns one prime at a time.

## Exact recognition of the Newton and residue data

At fixed \(N,p\), the carry graph is finite, with directed
edges for each allowed input digit pair at each state. A vector
of nonnegative integer edge multiplicities occurs on a walk
starting and ending at the zero carry exactly when the indegree
and outdegree multiplicities agree at every vertex, and the
vertices of positive degree form a connected undirected graph
containing zero. The all-zero vector gives the empty walk.
For sufficiency, equal in- and outdegrees make every weakly
connected nontrivial component a single strongly connected
component. Otherwise a source in the directed graph of
components would have more outgoing than incoming edges.
Splicing directed cycles in this connected component yields a
walk using precisely the given multiplicities. To represent
distinct exponent columns, require positive total
multiplicity on edges with \(A\ne B\).

These conditions are expressible in Presburger arithmetic.
For example, enumerate the finitely many candidate positive
edge supports to impose connectivity, and on each support use
linear equations for the indegree and outdegree constraints.
The vector \(D\) is an integer linear combination of the edge
counts. The coordinates of \(\delta,U\) are linear congruences
after choosing discrete logarithms in the finite cyclic group
\(\mathbb F_p^\times\). The exponents at \(p=2\) are instead
trivial. Therefore (1) gives exactly the possible first \(p\)-adic
coefficient data, by finite formulas. With a fixed list of
permitted row types, the same construction can use independent
nonnegative multiplicities of those rows directly.

For clarity, the associated common-factor screen also has a
finite Presburger description. The lower Newton polygon is
determined by finitely many comparisons of rational slopes,
each a linear inequality with fixed integer coefficients
after clearing a fixed endpoint-width denominator. Fix a shared
slope \(-a/b\) in lowest terms. On its faces the supported
indices are congruent modulo \(b\). The residual polynomial
for each face uses the units \(E_{j,k}\) at those indices,
in powers of degree steps \(b\), after removing the monomial
at the left endpoint of the face. To justify this, take
\(\pi^b=p\) and substitute \(t=\pi^a z\) over the extension
of \(\mathbb Q_p\) with uniformizer \(\pi\). Divide by the
smallest occurring power of \(\pi\) and reduce its coefficients.
The residue polynomial is \(z^k\) times that face's residual
polynomial in \(z^b\), where \(k\) is the left endpoint index.

For a factor defined over \(\mathbb Q_p\), the same substitution
and reduction commute with multiplication after scalar
normalization. Its stripped residue in \(z^b\) must divide
the gcd of the two face residues, and its width at that slope
is \(b\) times the residual divisor degree. These widths sum
to the degree of a common rational factor across the shared
slopes. This gives a necessary common-degree list by taking
all sums of \(b\) times divisor degrees of the shared residual
gcds. A rational common factor must survive this list at every
prime used. At terminal six its degree must be \(2\) or \(4\),
since the order-six ordinary polynomial is positive on the
real axis and the lower degree is five.

One enumerates the finitely many edge endpoint choices,
on-face index subsets and residues of valuation differences
modulo each candidate face width. These residues determine
the reduced denominator \(b\) and hence the residual indices.
The units range over a finite set. The gcd and divisor-degree
lists can then be computed for each enumerated choice.
Consequently survival of this first screen at one prime,
including existential tests over arbitrary supports, is a
decision problem in Presburger arithmetic. For a fixed row
alphabet one can combine any finite list of primes by using
the same row multiplicities in all their formulas.

## A simultaneous obstruction at finitely many primes

**Theorem.** Fix \(N\ge2\), a finite nonempty set \(\mathcal P\) of
primes, positive integers \(L_p\) for \(p\in\mathcal P\), and any
nonempty open interval \(I\subset(0,\infty)\).
There is a fully shared ordinary profile with distinct exponent
columns and total-degree difference zero such that
\[
 \frac{M_j(t)}{M_j(0)}\equiv(1+t)^j\pmod{p^{L_p}}
 \qquad(j=N-1,N,\ p\in\mathcal P).
 \tag{3}
\]
It can be chosen with a common nonzero absolute gap in every row and
with curvature \(\sum_i(b_i-a_i)^2/\max(a_i,b_i)\in I\).
The increasing and decreasing rows need not occur in mirror pairs.

Start with pairs of mirrored rows \((A,A+D),(A+D,A)\)
using a common positive integer gap \(D\). For a prime \(p\),
order \(j\), and \(1\le k\le\lfloor j/2\rfloor\), the difference
from the \(k=0\) factorial valuation of this pair is
\[
 K_{p,j,k}(A)=
 \sum_{h=1}^{kD}v_p(jA+h)+\sum_{h=1}^{(j-k)D}v_p(jA+h)
                  -\sum_{h=1}^{jD}v_p(jA+h).
 \tag{4}
\]
The values at the remaining indices follow by reflection.
The fractions \(k/j\) in (4) for the two orders are distinct:
a nonzero equality across orders would imply that \(N-1\)
divides the smaller order's index. List them as
\(0<\tau_1<\cdots<\tau_m\le1/2\), and put \(\tau_0=0\).
Choose \(D\) large enough that for every \(i\) there is an integer
\(u_i\) with \(u_i/D\in(\tau_{i-1},\tau_i)\).

Fix one target prime \(p\in\mathcal P\), one index \(i\), and take
\(u=u_i\) or \(D-u_i\). For a large parameter \(H\), impose
\[
 A+u\equiv p^H\pmod{p^{H+1}} .
\]
For each other prime \(q\in\mathcal P\), choose a residue of \(A\)
modulo \(q^{B_q}\) for which no \(jA+h\) is zero modulo \(q^{B_q}\)
for \(j=N-1,N\) and \(1\le h\le jD\).
Such a residue exists for sufficiently large fixed \(B_q\):
each prohibited congruence excludes at most \(q^{v_q(j)}\)
residues, and there are only \((2N-1)D\) congruences.
The Chinese remainder theorem gives a positive \(A\) with
all these properties.

At the target prime, the term at \(h=ju\) in each of the two
orders has valuation \(H+v_p(j)\). Every other term in (4)
has a valuation bounded independently of \(H\), because
\(jA+h=j(A+u)+(h-ju)\) and the second summand is then nonzero.
At all other primes the valuations are bounded by construction.
As \(H\to\infty\), the vector of all \(K_{q,j,k}(A)/H\)
therefore converges, in the target-prime block, to the vector
with coordinates
\[
 {\bf1}_{\{u/D<k/j\}}+{\bf1}_{\{u/D<1-k/j\}}-1,
 \tag{5}
\]
and to zero in the other blocks. For \(u=u_i\), the vector in
(5) has ones precisely from position \(i\) onward in the
sorted fractions. Taking \(u=D-u_i\) negates the vector.
These \(m\) nested vectors form a basis in each target-prime
block.

Let \(P_H\) be the square matrix whose columns are the
unscaled vectors of (4) for all the positive choices, and let
\(b_H\) be the sum of the vectors for the negative choices.
Then \(P_H/H\) tends to the invertible block basis matrix,
and \(-P_H^{-1}b_H\) tends to the all-ones vector. For large
integer \(H\) it is a strictly positive rational vector,
since \(P_H,b_H\) have integer entries. Use its coordinates
as rational multiplicities on the positive choices, with
one on every negative choice. Multiplying all multiplicities
by a common denominator gives a finite multiset of pairs
for which every quantity in (4) sums to zero. The factorial
coefficient ratios to the constant term thus have valuation
zero at every specified prime and order.

Let \(R\) be a common multiple of
\(|(\mathbb Z/p^{L_p}\mathbb Z)^\times|\) for the chosen primes.
Repeating this whole profile \(R\) times raises each factorial
coefficient ratio to the \(R\)-th power. All of those ratios
are then 1 modulo \(p^{L_p}\). After multiplication by the
binomial coefficients, this gives (3).

Finally one can vary the bases without changing the argument.
For the finite constructed multiset, choose a positive integer
\(B\) divisible, at each \(p\in\mathcal P\), by a power of \(p\)
strictly larger than all valuations of all sampled terms
\(jA+h\), \(1\le h\le jD\). Adding any nonnegative multiple of
\(B\) to both exponents of any one row preserves all those
individual valuations. It therefore preserves the zero sums,
including if the two members of each pair receive different
shifts. Choose distinct lower bases across the two orientations
by varying those shifts. This preserves the signed total gap
and yields a profile with no reversed unequal rows. Repeating
the shifted profile by \(R\) still neutralizes all units.
Next add a further common shift \(TB\) to the bases for a
large integer \(T\). The resulting curvature \(S_T>0\) tends
to zero with \(T\), since the absolute gaps are fixed. Choose
a bounded open subinterval \((s_-,s_+)\) of \(I\), then take
\(T\) with \(S_T<s_+-s_-\). Repeating the shifted neutralized
profile \(\lfloor s_-/S_T\rfloor+1\) times puts the curvature
in that subinterval and retains (3). This completes the proof.

For terminal six, the first valuation-and-unit screen at any
fixed finite set of primes can therefore allow both candidate
degrees \(2\) and \(4\) simultaneously, even in a prescribed
positive curvature interval. The congruence (3) also applies
at any prescribed finite precision at those primes. A
precision bound that depends on the exponent profile is a
separate question.
