# The trace cofactor at split real-cubic primes

Let \(B/K\) be an abelian threefold whose geometric rational
endomorphism algebra is a totally real cubic field \(F\), with all
endomorphisms defined over \(K\). Write \(a_v\in\mathcal O_F\)
for the relative Frobenius trace. Consider good absolute
degree-one primes \(v\) with residue characteristic \(p\)
split completely in \(F\). Define
\[
                         b_v=\frac{N_{F/\mathbb Q}(a_v)}{p}.
\]
The reduction is nonordinary exactly when \(b_v\in\mathbb Z\).
In fact a split height-two factor fails exactly when its
relative trace is divisible by the corresponding coefficient
place, and the rational norm detects at least one such place.
The Weil bound gives \(|b_v|\le 8\sqrt p\). The bounds below
use the trace cofactor rather than the
[squareclass of the discriminant cofactor](CM_NEWTON_FROBENIUS_DISCRIMINANT_SIEVE.md).

## Counting small trace cofactors

**Theorem.** Count the above primes with \(p\le X\), \(b_v\) integral
and \(|b_v|\le H\). Unconditionally their number is
\[
 O_{B,K}\!\left((H+1)\frac{X(\log\log X)^{1/6}}{(\log X)^{7/6}}\right)
 \qquad(X\ge 3).
\]
Assuming GRH for the Dedekind zeta functions of the torsion
fields \(K(B[\ell])\) and their subfields, the count with
\(|b_v|\le p^\delta\) is
\[
 O_{B,K,\delta}\!\left(X^{16(1+\delta)/17}(\log X)^2\right)
 \qquad(0\le \delta<1/16).
\]
If GRH also holds for the splitting fields over \(\mathbb Q\)
of the rational Frobenius polynomials of \(B_v\) and their
composita with the normal closure of \(F\), the count is
\[
 O_{B,K,\delta}\!\left(X^{10(1+\delta)/11}(\log X)^2\right)
 \qquad(0\le \delta<1/10).
\]

**Residual trace sets.** At a sufficiently large auxiliary
rational prime \(\ell\) split in \(F\), the image \(G\) on
\(B[\ell]\) consists of triples of invertible two-by-two matrices
over \(\mathbb F_\ell\) with common determinant \(d\). This is
the full rank-two real-multiplication image of
[Lombardo, Theorem 1.4](https://arxiv.org/pdf/1411.1703).
For \(p\ne\ell\),
\[
 b_v\bmod\ell=d^{-1} t_1t_2t_3,\qquad t_i=\operatorname{Tr}g_i,
\]
whenever \(b_v\) is integral. Its class function on \(G\)
will be denoted \(b\).

Every element of \(G\) is conjugate in \(G\) into one of a
bounded number of abelian subgroups \(A\) of order
\(\asymp\ell^4\). Take the product with common determinant
of three groups of one of the following types: split torus,
nonsplit torus, or matrices
\(\{\left(\begin{smallmatrix}u&w\\0&u\end{smallmatrix}\right):
u\ne0\}\). The three types cover the characteristic polynomial
cases for each matrix, using a split torus for a scalar.
Each determinant image contains the squares and each
determinant kernel has order at least \(\ell-1\).
The individual conjugators can be made to have determinant
one by making their determinants squares and rescaling them.
For the first adjustment allow each listed group to be
conjugated by a fixed diagonal matrix of nonsquare
determinant. This leaves at most \(6^3\) subgroup choices.

On each \(A\), the fraction of elements with \(b\) in a set
of \(h\) residue classes is \(O((h+1)/\ell)\). For \(b\ne0\),
scaling all three matrices by the same unit scales \(b\)
by that unit. Thus all nonzero fibers have equal size.
For \(b=0\), some \(t_i=0\). At a given \(d\) there are
at most two trace-zero elements in either torus, and none
in the repeated-eigenvalue group. Dividing by the respective
determinant kernel proves the zero bound. Let \(C\) be the
subset of \(G\) where \(b\) lies in the specified residues.
The same fraction bound holds on \(G\), by common scaling
and the fact that at fixed determinant the trace-zero
probability in a full matrix factor is \(O(1/\ell)\).

**Unconditional bound.** Put \(L_\ell=K(B[\ell])\).
For each conjugacy class in \(C\), choose one of the
abelian subgroups \(A\) meeting it and set \(M=L_\ell^A\).
The bounds \(|G|\ll\ell^{10}\) and \(|A|\asymp\ell^4\)
give \(n_M=[M:\mathbb Q]\ll_{B,K}\ell^6\).
The possible ramified rational primes are in a fixed
finite set and \(\{\ell\}\). The standard local different
bound \(d\le e-1+e\,v(e)\), with \(v\) normalized on the
base local field, yields \(\log D_M\ll n_M\log\ell\).
For \(Q\) the largest absolute conductor norm of a
character of \(L_\ell/M\), it also gives
\(\log Q\ll n_M\log\ell\). Indeed its local Artin conductor
exponent is at most \(2d/e\), since in the lower ramification
formula \(|G_i|\le 2(|G_i|-1)\) at every nontrivial group
and a character has codimension of invariants at most one.
Summing with the local residue degrees gives the asserted
absolute conductor bound.

Take \(\ell\) split in \(F\) and of order
\(c(\log X/\log\log X)^{1/6}\), with \(c>0\) small enough
depending on \(B,K\). Such a prime in a fixed-ratio
interval exists by Chebotarev. Then
\(\log(D_M Q n_M^{n_M})\ll \ell^6\log\ell\) is small enough
to apply Thorner–Zaman Theorem 1.4. Its exceptional-zero
term always gives an upper bound of at most twice the
unmodified main term. Consequently each class contributes
\(\ll |C_{\rm class}|\operatorname{Li}(X)/|G|\).
Summing over \(C\), with \(h\le 2H+1\), gives the
logarithmic estimate. The at most \([K:\mathbb Q]\) places
over the auxiliary prime itself can be added directly.
[Thorner–Zaman, Theorem 1.4](https://msp.org/ant/2019/13-5/ant-v13-n5-p02-p.pdf#page=6).

**GRH bounds.** For an abelian quotient of order
\(\asymp\ell^4\), let \(D\) be the subset of elements
with \(b\) in the indicated residue classes. Its size
is \(O((h+1)\ell^3)\) and the sum of absolute values
of the Fourier coefficients of its indicator is at most
\(\sqrt{|D|}\), by Parseval and Cauchy–Schwarz. The
smoothed prime-ideal explicit formula for each abelian
character, on \(X/2<Nw\le X\), then gives the upper bound
\[
 \ll \frac{(h+1)X}{\ell\log X}
       +\sqrt X\,\sqrt{h+1}\,\ell^{3/2}n_M(\log X)^2 .
\]
For the explicit formula take a fixed nonnegative smooth
function supported at positive arguments and at least one
on \([1/2,1]\). Shifting its Mellin integral against the
Hecke logarithmic derivative gives the trivial-character
main term and a zero term
\(O(\sqrt X(\log(D_MQ)+n_M\log X))\) by GRH,
rapid Mellin decay, and zero counting. The prime powers
and ramified primes are absorbed in this bound with the
stated logarithmic allowance. GRH for the Hecke factors
follows from GRH for the extension's Dedekind zeta function,
using its abelian factorization. Every prime of \(K\) whose
class meets \(A\) contributes to this count through at least
one residue-degree-one prime of \(M\) above it with the
required class. This follows from the action on cosets \(G/A\).

Use the bounded collection of abelian subgroups with
\(n_M\ll\ell^6\). Take \(h\asymp X^\delta+1\) and
\(\ell\asymp X^{(1+\delta)/17}\). The two terms have
equal power \(16(1+\delta)/17\), giving the second bound
after summing the dyadic intervals.

For the stronger GRH assumptions, the rational Frobenius
splitting field of each \(v\), together with the normal
closure of \(F\), has bounded degree and discriminant
bounded by a fixed power of \(p\). For the latter bound,
adjoin the three roots of the squareclasses of the
relative discriminants, each of norm \(O(p^3)\), and
apply the discriminant tower formula. Under GRH and effective
Chebotarev, uniformly for \(p\le X\), there are
\(\gg Y/\log Y\) auxiliary primes \(\ell\in[Y,2Y]\)
split completely in this field, for \(Y=X^\theta\)
with any fixed \(\theta>0\). Discard the bounded bad
auxiliary set and \(p\) itself. All three residual
matrices at each of these \(\ell\) can be triangularized.
The common-determinant Borel \(P\) meets these classes.
Quotient it by the three unipotent radicals, taking
the diagonal abelian quotient of order \(\asymp\ell^4\).
The trace function descends to that quotient, with
the same size bound for \(D\). The fixed field
\(M=L_\ell^P\) now has degree \(O(\ell^3)\).
The same Hecke bound, summed over \(\ell\in[Y,2Y]\)
split in \(F\) and divided by the minimum
\(\gg Y/\log Y\) suitable for each counted \(v\),
gives the previous two-term bound with \(n_M\ll Y^3\).
Taking \(\theta=(1+\delta)/11\) gives the third
estimate. \(\square\)
[Lagarias–Odlyzko, Theorem 1.1 and smoothed explicit formula](https://aareyanmanzoor.github.io/assets/articles/lagarias-odlyzko.pdf).

## An unconditional trace bound in arbitrary dimension

**Proposition.** Let \(B/K\) have dimension \(g\) and geometric
rational endomorphism algebra a totally real field of degree \(g\),
with all its endomorphisms defined over \(K\). At good absolute
degree-one primes of characteristic \(p\) unramified in the
coefficient field and prime to the endomorphism-order index, put
\(r_v=g-f_p(B_v)\) and \(b_v=N_{F/\mathbb Q}(a_v)/p^{r_v}\).
Then \(b_v\) is an integer. The count with \(p\le X\) and
\(|b_v|\le H\) is unconditionally
\[
 O_{B,K}\!\left(
   (H+1)\frac{X(\log\log X)^{1/(2g)}}{(\log X)^{1+1/(2g)}}\right).
\]

**Proof.** At a coefficient place over \(p\) the trace is a nonunit
precisely when the local rank-two factor fails to be ordinary.
Such a place contributes its residue degree to \(r_v\) and at least
that much to the valuation of the trace norm. This proves integrality
also at nonsplit \(p\). For the counting argument fix a value
\(r\in[0,g]\) and a split auxiliary prime \(\ell\). On the
common-determinant group in \(g\) factors the test for \(b_v\) is
\(d^{-r}\prod_i t_i\). Conditional on the determinant, the traces
are independent, with probability \(O(1/\ell)\) for any given
trace by the fixed-determinant count. The probability of zero
product is \(O(g/\ell)\). Conditional on nonzero first \(g-1\)
traces, a given nonzero target imposes one value of the last trace.
Thus a set of \(h\) possible residues still has probability
\(O_g((h+1)/\ell)\).

Use the same abelian subgroups as in the cubic proof, now with \(g\)
factors. Their orders are \(\asymp_g\ell^{g+1}\) inside a group of
order \(\asymp_g\ell^{3g+1}\). Every class meets one of them,
and the fixed field and conductor bounds apply with
\(n_M=O_{B,K}(\ell^{2g})\). The class-wise upper bound therefore
holds for an auxiliary split prime of order
\((\log X/\log\log X)^{1/(2g)}\) times a sufficiently small
constant. Sum over the indicated class set and the \(g+1\)
possible values of \(r\). \(\square\)

## A separate sparsity condition without GRH

Fix a quadratic number field \(J\) not contained in the normal
closure of \(F\). The nonordinary split primes for which
\(b_v\) is a norm from \(J\) have density zero. This even
allows norms of nonintegral elements of \(J\). To prove it,
use auxiliary rational primes \(\ell\) split in the normal
closure of \(F\) and inert in \(J\). They have positive
density by disjointness. A nonzero rational norm from \(J\) has
even valuation at each such \(\ell\), and zero also avoids valuation
one. The condition
\(v_\ell(N_{F/\mathbb Q}(a_v))=1\) is detected modulo
\(\ell^2\), and its probability is \(\gg 1/\ell\).
For example require one trace to vanish modulo \(\ell\)
and have a nonzero trace digit at the next level, and
the other two traces to be units. At trace zero and
unit determinant the characteristic polynomial is
separable in odd characteristic. Lifts through
\(\ell^2\) with fixed determinant give each next trace
digit uniformly, since the trace and determinant
differentials are independent at a nonscalar matrix.
The conditions in the other two factors have probability
bounded away from zero.

For a finite set of sufficiently large \(\ell\), the events
are independent by the full product image and Chebotarev.
The product image follows by passing to the independence
extension and using the individual large images there.
Outside these finitely many auxiliary residue
characteristics, \(v_\ell(b_v)=v_\ell(N(a_v))\).
Thus the upper density of the norm condition is bounded
by products of \(1-c/\ell\) over arbitrarily many of the
eligible primes. These products tend to zero. A finite
union over such fixed quadratic fields still has density
zero.
[Serre, Theorem 1 and §3.1](https://www.college-de-france.fr/media/jean-pierre-serre/UPL8570957518086251784_independance.pdf).

The small-cofactor and fixed-norm conditions each cut out a
zero-density part of the split nonordinary locus. The
one-failure Newton polygon only yields the bound
\(|b_v|\le8\sqrt p\) used at the start. A density theorem
for the full split locus needs information about the
actual diagonal \(p\mid N_{F/\mathbb Q}(a_v)\), beyond
these conditions on the integer cofactor.
