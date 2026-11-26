# A uniform three-adic restriction on the shifted ternary cubic gate

For \(j\ge1\), set
\[
 F_j(A,B,C)=[z^j]\exp\!\left(j(Az+Bz^2/2+Cz^3/6)\right).
\]
These are the [noncollinear ternary limit forms](SFC_HOMOGENEOUS_TERNARY_SHIFTED_GATES.md#noncollinear-offsets).
Normalize \(v=v_3\) by \(v(3)=1\). If \(C\ne0\) and the three forms
at orders \(m,m+1,m+2\) vanish, let \(n\) be the unique multiple of
three in the window. Then \(A,B\ne0\) and, at each extension of the
three-adic valuation to the coefficient field,
\[
 v(A)-\tfrac12v(B)=\tfrac14,\qquad
 v(C)-\tfrac32v(B)\ge\tfrac12(1+v(n)).                            \tag{1}
\]
In particular the centered case \(A=0\) is excluded at every start, even
over an arbitrary characteristic-zero field. Rational triples are also
excluded, since \(v(A^2/B)=1/2\).
The centered assertion applies over arbitrary
characteristic-zero fields by the weak Nullstellensatz applied to the
rational system with \(A=0\) and an auxiliary inverse for \(C\).
The inequalities leave other cubic parameters open.

## Initial valuation restrictions

Use weighted scaling \((A,B,C)\mapsto (\lambda A,\lambda^2B,\lambda^3C)\).
It preserves the zeros. Take \(C=1\) in the algebraic closure.
Write \(P_j(X,Y,Z)=j![z^j]\exp(Xz+Yz^2/2+Zz^3/6)\).
Its coefficients count labelled partitions into singletons, pairs,
and triples. The three-cycle on three marked labels gives
\[
 P_{j+3}(X,Y,Z)\equiv(X^3+Z)P_j(X,Y,Z)\pmod3 .
\]
Indeed, a fixed partition leaves the three labels as singletons or as
a single triple. If three blocks are permuted, each contains just one
of the marked labels and no fixed label, and if a block contains the
three marked labels it is already full. This is also the specialization
to maximum block size three of [Carlitz's congruence as stated by Gessel, (5.5)](https://people.brandeis.edu/~gessel/homepage/papers/combinatorial-congruences.pdf).

At \(n=3k\), divide \(n!F_n\) by \(n^k\). The pure triple term has
three-adic valuation zero, as \(v(n!)=k+v(k!)\). All other terms
are integral partition coefficients times powers
\((n^{2/3}A)^a(n^{1/3}B)^b\) with \(a+b>0\). They therefore
have positive valuation if \(v(A)>-2v(n)/3\) and
\(v(B)>-v(n)/3\). In particular, vanishing at \(n\) implies
\(\min(v(A),v(B)/2,0)<0\).

Scale so that this minimum becomes zero. The rescaled \(C\) now has
positive valuation. At the order \(j\equiv1\bmod3\), the partition
congruence for \(P_j(jA,jB,jC)\) forces the rescaled \(A\) to
have positive valuation. Thus, in the normalization \(C=1\),
\(B\ne0\) and \(v(A)>v(B)/2\), allowing \(v(A)=+\infty\).
If \(v(B)>-v(n)/3\), then \(v(A)>-2v(n)/3\) too, which was
excluded by the pure triple term. Undoing the weighted scaling proves
the preliminary bounds \(v(A)-v(B)/2>0\) and
\(v(C)-3v(B)/2\ge v(n)/2\).

## Reduction by replacing singletons

Normalize \(B=1\), write \(\alpha=v(A)\), \(\gamma=v(C)\), and
suppose \(\alpha\ge1/3\). The case \(A=0\) is included by taking
\(\alpha=+\infty\). At start one the equations \(F_1=F_2=0\)
force \(A=B=0\), so assume \(m\ge2\).
At any order \(j\), let \(e=v(j)\). For \(a+2b+3c=j\), the
valuation of the corresponding term in \(F_j\) is
\[
 a(\alpha+e)+b e+c(\gamma+e-1)-v(a!)-v(b!)-v(c!).            \tag{2}
\]
Replace \(2r\) singletons by \(r\) pairs, where \(r=\lfloor
a/2\rfloor\). The excess valuation of the original term is at least
\[
 2r\alpha+re-v(a!)+v(r!).
\]
For \(a\ge2\), this is positive except possibly at unit orders
with \(a=3,9\). To check the bound, count multiples of \(3^i\)
in \((r,2r+1]\). With \(L=\lfloor\log_3(2r+1)\rfloor\), their
total for \(r\ge18\) is at most
\((r/2)(1-3^{-L})+L<2r/3\), since \(L\le r/6\).
At \(1\le r\le17\) the same count is
\[
 1,1,1,3,3,3,4,4,4,5,5,5,8,8,8,9,9 ,
\]
and for \(a=2,8\) the count ends at \(2r\) and is respectively
zero and one. This proves \(v(a!)-v(r!)<2r/3\) except for
\(a=3,9\). When \(e>0\), both exceptions are strictly lowered
by the same replacement since \(v(a!)\le r\).

We now show that \(\gamma>1\) is impossible. If an even order
\(j\) in the window is prime to three, compare with the pure pair
term. The number of removed pairs is \(R=(a+3c)/2\), so the
pair factorial contributes at least \(v(R!)\). For \(a=0,1\)
this offsets \(v(c!)\). For \(a=3,9\) the parity implies \(c\ge1\).
At \(a=3\) the excess valuation is at least
\(3\alpha-1+c(\gamma-1)>0\). At \(a=9\), \(R-c\ge5\)
gives \(v(R!)\ge v(c!)+1\), so the excess is at least
\(9\alpha-3+c(\gamma-1)>0\). Thus the pure pair term is
uniquely least.

Otherwise there are two odd unit orders around the even multiple
of three. At the one with \(j\equiv1\bmod3\), the replacement
for \(a=3,9\) increases the pair count to \(b'\equiv0\bmod3\).
The removed strip from \(b'!\) then contributes at least one
for length one and two for length four. The two exceptional
terms are therefore strictly above their replacements by (2).
At fixed \(a\in\{0,1\}\), the least admissible number of triples
is \(c_*=1-a\). Increasing it by \(2r\), \(r\ge1\), removes
\(3r\) pairs. This raises (2) since \(\gamma>1\) and
\(v((c_*+2r)!/c_*!)\le v((3r)!)\).
For finite \(\alpha\) the only possible tie at this \(j\) requires
\[
 \gamma=\alpha+1-v((j-1)/2)\le\alpha .
\]
This excludes \(\gamma>1\) for \(\alpha\le1\). If \(\alpha>1\)
then pair replacement lowers all terms with \(a\ge2\), including
at the odd order \(j'\equiv2\bmod3\). The same comparison there
requires \(\gamma=\alpha+1\), incompatible with the first tie.
When \(A=0\) only the \(a=0\) branch exists and its minimum is
unique. Thus \(\gamma\le1\).

For \(\gamma\le1\), replace three or nine singletons by one or
three triples, respectively. The respective excess valuations
in (2) are at least \(3\alpha+2e-\gamma\) and
\(9\alpha+6e-3\gamma\). The latter bound uses that the three new
consecutive triple indices contribute at least one to their factorial
valuation. The replacement is strict unless \(\alpha=1/3\),
\(\gamma=1\), and \(e=0\). At \(j\equiv0,2\bmod3\), every
\(a=1\) term has \(b\) a three-adic unit. Replacing its singleton
and a pair by a triple is strict, with excess in (2) equal to
\(\alpha+e+1-\gamma+v(c+1)-v(b)>0\).

## Opposite bounds on the centered Newton slopes

For \(j\ge2\), choose \(b_0\in\{0,1,2\}\) with \(2b_0\equiv
j\bmod3\), and set \(c_0=(j-2b_0)/3\), \(K=\lfloor c_0/2\rfloor\),
and \(s=c_0-2K\). After a nonzero monomial factor, the centered
polynomial is a polynomial in \(W=B^3/C^2\), of degree \(K\).
For \(b=b_0+3k,\ c=c_0-2k\), its coefficient valuation relative
to its constant coefficient is
\[
 U_k=(1+v(j))k+v\bigl(c_0!/(c_0-2k)!\bigr)-v(k!).
\]
This follows by taking the consecutive coefficient ratio and using
\(v((b_0+3k)!/b_0!)=k+v(k!)\).

At \(3\mid j\), all slopes of the lower coefficient polygon are
at least two. Indeed \(v(c_0!/(c_0-2k)!)\ge v((2k)!)\ge v(k!)\).
Equality in \(U_k\ge2k\) at \(k>0\) requires \(v(j)=1\), \(k=1\),
and \(c_0\equiv2\bmod3\), because \(v((2k)!/k!)>0\) when \(k\ge2\).
In particular slope two can occur only under these necessary conditions.

At \(3\nmid j\), all slopes are at most two. For \(1\le r\le K\),
\[
 U_K-U_{K-r}
 =r+v\bigl((s+2r)!/s!\bigr)-v\bigl(K!/(K-r)!\bigr)\le2r .
\]
Here \(v((s+2r)!)\le r\). Equality requires \(s=1\),
\(2r+1\) a power of three, and a length-\(r\) factorial strip
consisting of units. The strip condition gives \(r\le2\), so
\(r=1\) and \(3\nmid K\). In particular slope two requires
\(c_0\) odd and \(3\nmid K\).

A sum can vanish only if its least valuation is attained at least
twice. When only centered terms can minimize (2), their polygon
must therefore have slope \(-v(W)=2\gamma\). At the multiple of
three only centered terms can minimize even when \(\alpha=1/3\).
The lower slope bound and \(\gamma\le1\) force \(\gamma=1\),
slope two, and \(n/3\equiv2\bmod3\).

First let \(\alpha>1/3\). At the order equal to two modulo
three, only centered terms can minimize (2), so it too must
have slope two. Write \(n=3l\). If the other order is \(3l-1\),
its \(c_0=l-1\) must be odd, but then \(K=(l-2)/2\) is
divisible by three. The only remaining placement starts at
\(3l\) and contains \(3l+2\), whose \(c_0=l\) is odd.
The remaining \(3l+1\) is an even unit order. There the
\(a=1\) terms are strictly above pure pairs by the even-order
comparison at \(\gamma=1\), and its centered \(c_0=l-1\) is
even. Its centered polygon has no slope two.

Finally let \(\alpha=1/3\). At any even unit order, the
least centered term is unique at \(c=0\), since the centered
\(c_0\) is even and all slopes are strictly less than two.
An exceptional replacement of three or nine singletons by
triples cannot produce this term. Terms with \(a=1\) are
strictly above pure pairs as before. Thus there can be no
even unit order. The multiple of three \(n=3l\) must be the
even middle order, and \(j=n+1\) is an odd unit order.
Its centered \(c_0=l-1\) is odd, but \(K=(l-2)/2\) is
divisible by three. Its uniquely least centered term
therefore has \(c=1\). The \(a=1,c=0\) term lies strictly
below it, since \(v((j-1)/2)=v(3l/2)=1\), and the excess
valuation over the \(a=0,c=1\) term is \(\alpha-1\).
Within the \(a=1\) branch, increasing \(c\) from zero to
\(2r\) removes a strip of \(3r\) pairs. That strip has
valuation at least \(r>v((2r)!)\), so those terms are higher.
The exceptional \(a=3,9\) terms are at least as high as
centered terms, and all other replacements were strict.
This leaves a unique least term and excludes the boundary.

## The lower bound from the Gaussian pair at the unit orders

Continue with \(B=1\), \(\alpha=v(A)\in(0,1/3)\), \(\gamma=v(C)\),
and \(e=v(n)\ge1\). At order \(n\), replacing all but zero or one
singleton by pairs is strict for every \(\alpha>0\). Indeed for
\(a\ge2\), \(r=\lfloor a/2\rfloor\), the excess is at least
\(2r\alpha+re-v(a!)>0\). If \(\gamma<(1+e)/2\), the \(a=1\)
terms are also strictly lowered by replacing a singleton and one
pair by a triple. Here \(b\equiv1\bmod3\), and the excess is
\(\alpha+e+1-\gamma+v(c+1)>0\). Only centered terms could
then minimize (2), whose smallest Newton slope is at least \(1+e\)
by the centered coefficient formula. This contradicts the assumed
\(2\gamma<1+e\), proving the stronger lower bound for \(\gamma\).

At either order \(j\) prime to three, every term with \(c>0\)
has higher valuation than the term formed by replacing all \(c\)
triples by \(3c\) singletons. In (2) the excess is
\[
 c(\gamma-1-3\alpha)
 +v\bigl((a+3c)!/a!\bigr)-v(c!)\ \ge\ c(\gamma-3\alpha)>0.
\]
The inequality holds because a length-\(3c\) factorial strip has
\(c\) multiples of three and their quotients by three give a
length-\(c\) factorial strip. Thus a zero at either unit order
requires cancellation in the initial form of its pure Gaussian
part. We show the two Gaussian initial forms cannot both vanish
when \(0<2\alpha<1/2\).

For \(K\ge0\), \(\epsilon\in\{0,1\}\), and a three-adic unit
\(J\), write
\[
 P_{K,\epsilon,J}(x)=
   \sum_{r=0}^K \frac{J^r K!}{(K-r)!(2r+\epsilon)!}x^r .
\]
After removing a nonzero monomial factor, the Gaussian part of
\(F_j\) is \(P_{K,\epsilon,2j}(A^2)\), where \(2K+\epsilon=j\).
At \(x\ne0\) set \(\lambda=v(x)\). The initial form at \(x\)
means the reduction of the sum of terms of least valuation after
dividing by any one of those terms. Its vanishing is independent
of the divisor. We compare these initial forms by ternary division.

Write \(K=3K_0+k\) and \(r=3t+s\), with \(k,s\in\{0,1,2\}\).
Put \(\eta=1\) if \(s>k\), and zero otherwise, and put
\(\ell=\lfloor(2s+\epsilon)/3\rfloor\). Factoring the
multiples of three in the factorials gives the valuation
\[
 \begin{split}
 &(3\lambda-1)t+s\lambda+\eta-\ell\\
 &\qquad
  +v\bigl(K_0!/(K_0-t-\eta)!\bigr)-v\bigl((2t+\ell)!\bigr).
 \end{split}                                                       \tag{3}
\]
A residue class with \(K_0-\eta<0\) is omitted. Let \(M(D,E)\)
denote the minimum valuation of the terms of
\(P_{D,E,J}\) at an argument of valuation
\(\lambda'=3\lambda-1\). It depends on \(D,E,\lambda'\)
and not on the unit \(J\). Write \(M_0=M(K_0,0)\),
\(M_1=M(K_0,1)\), and \(M'_i=v(K_0)+M(K_0-1,i)\)
when \(K_0>0\). We have \(M_1\le M_0\) termwise.
Also \(M'_i\ge M_0-\lambda'\) for \(i=0,1\): compare the
term at \(t\) in \(M'_i\) to the term at \(t+1\) in
\(M_0\), whose factorial denominator is \((2t+2)!\).

Suppose \(0<\lambda<1/2\). The exact least-valuation residue
class in (3) is
\[
\begin{array}{c|ccc}
 & k=0 & k=1 & k=2\\ \hline
 \epsilon=0 &s=0&s=0&s=2\\
 \epsilon=1 &s=0&s=1&s=1 .
\end{array}
\]
Indeed, \(s=0\) always gives \(M_0\). For \(\epsilon=0\),
\(s=1\) has an additional \(\lambda\) if \(\eta=0\), and
otherwise gives \(M'_0+1+\lambda>M_0\).
The class \(s=2\) at \(\eta=0\) gives \(M_1+2\lambda-1<M_0\),
and at \(\eta=1\) gives \(M'_1+2\lambda>M_0\).
For \(\epsilon=1\), \(s=1\) at \(\eta=0\) gives
\(M_1+\lambda-1<M_0\), and at \(\eta=1\) gives
\(M'_1+\lambda>M_0\). The class \(s=2\) is higher than
\(s=1\) when they both have \(\eta=0\), and is higher than
\(M_0\) when \(\eta=1\). These strict comparisons use the
displayed bounds on \(M'_i\) and \(\lambda<1/2\).
The selected inner degree is
\[
             2(K_0-\eta)+\ell=\lfloor(2K+\epsilon)/3\rfloor .
\]

The initial form on that residue class equals the initial form
of \(P_{K_0-\eta,\ell,-J}(x^3/3)\), up to a nonzero residue
factor. To verify the sign, write \(N=3N_0+d\) and reduce
the unit part of \(N!\) modulo three. It is
\((-1)^{N_0}d!\) times the unit part of \(N_0!\), by grouping
the nonmultiples of three in blocks. In a fixed residue class
\(s\), the three factorials in the coefficient of \(x^r\)
therefore introduce the sign
\((-1)^{t+\eta-(2t+\ell)}\) and a digit-factorial ratio independent
of \(t\). The factor \(J^{3t+s}\) reduces to \(J^{t+s}\), so
the inner parameter has residue \(-J\). Its argument is
\(x^3/3\) by (3). This proves preservation of the initial-form
condition under the degree, parameter, and valuation changes.

Iteration from \(0<\lambda<1/2\) eventually reaches
\(\lambda'\le0\), since the distance from \(1/2\) triples
at each step. At a negative valuation only the highest-degree
term can minimize. At zero valuation, compare the coefficient
at \(K-r\) to the coefficient at \(K\): its excess valuation is
\[
 v\bigl((2K+\epsilon)!/(2K+\epsilon-2r)!\bigr)-v(r!)
 \ \ge\ v((2r)!)-v(r!) .
\]
For \(r\ge2\) this is positive. A tie with \(r=1\) is possible
only when the degree \(2K+\epsilon\) is two modulo three.
The terminal initial form in that case is a monomial times a
linear polynomial with nonzero constant and nonzero linear term.

The two original unit degrees differ by at most two. After at
least one ternary division, their reduced degrees differ by at
most one. If they are distinct, they cannot both be two modulo
three, so they cannot both have initial-form cancellation. If
they are equal, the two reduced parameters still have opposite
residues, because the original \(2j\)'s have opposite residues
and each division reverses both signs. For equal reduced degrees
the possible width-one initial forms therefore have distinct
nonzero roots. Hence the two unit orders cannot both vanish
when \(0<2\alpha<1/2\), which proves the lower bound in (1).

## The upper Gaussian clock

Let \(E_K(\lambda)\) and \(O_K(\lambda)\) be the minimum
valuations of the terms of \(P_{K,0,J}\) and \(P_{K,1,J}\),
respectively, at clock \(\lambda\). The unit \(J\) does not
affect either minimum. For \(K\ge1\) define
\[
 \delta_K=E_K-O_K,\qquad
 \Delta_K=E_K-O_{K-1}-v(K).
\]
Use \(\delta_0=0,\ \Delta_0=-\infty\) to omit unavailable
residue classes. Always \(\delta_K\ge0\) by the termwise comparison
of the two factorial denominators.

A clock \(\lambda>1\) gives a unique least constant term for all
states. At clock \(1\) the least terms consist of the constant
alone or, in the odd state with \(K\) prime to three, the constant
and the degree-one term. In fact \(v((2r+\epsilon)!)\le r\),
and equality for \(r>0\) requires \(\epsilon=1\) and
\(2r+1\) a power of three. The numerator factorial strip must
also be a unit, so \(r\le2\) and only \(r=1\) is possible.
In particular both minimum valuations are zero at all clocks
\(\lambda\ge1\).

Now let \(1/2<\lambda<1\). Write \(K=3K_0+k\), and evaluate
the primed deficits at \(K_0\) and \(\mu=3\lambda-1\).
The residue comparisons of (3) give
\[
\begin{array}{c|ccc}
 & k=0 & k=1 & k=2\\ \hline
 \delta_K &\max(0,\Delta'-\lambda)
          &\delta'+1-\lambda
          &\min(\delta'+1-\lambda,\lambda)\\[2pt]
 \Delta_K &\Delta'-\lambda
          &\max(0,\Delta'-\lambda)
          &\min(\delta'+1-\lambda,\lambda).
\end{array}                                                       \tag{4}
\]
For clarity, in an even state the minimizing class is \(s=0\)
if \(k=0,1\), while at \(k=2\) the minimum is
\(\min(E_{K_0},O_{K_0}+2\lambda-1)\) from \(s=0,2\).
In an odd state with \(k=1,2\) the class is \(s=1\), with
minimum \(O_{K_0}+\lambda-1\). For \(k=0\) the minimum is
\(\min(E_{K_0},O_{K_0-1}+v(K_0)+\lambda)\) from \(s=0,1\).
These give the first row of (4), and the second follows by
comparing the even state to the odd state with \(K-1\).
In the latter comparison with \(k=0\) one subtracts
\(v(K)=1+v(K_0)\). All unlisted residue classes are strictly
higher, using \(M'_i\ge M_0-\mu\), \(M_1\le M_0\), and
\(0<\lambda<1\), as in (3). Under this recursion
\[
                        \Delta_K(\lambda)\le \lambda
\]
for \(\lambda>1/2\). At \(\lambda\ge1\) this is immediate,
and in (4) it follows from \(\Delta'\le\mu\) and
\(\mu-\lambda=2\lambda-1<\lambda\).

There are just two possibilities for a tie between different
residue classes in (3). In an even state with \(k=2\) it requires
\(\delta'=2\lambda-1\). In an odd state with \(k=0\) it requires
\(\Delta'=\lambda\). To analyze both, put \(h=\lambda-1/2\)
at the current clock and consider target values
\[
                 \delta_K=a h
             \quad\text{or}\quad
                 \Delta_K=\tfrac12+b h,\qquad
                 a,b\in[1/3,2/3].                              \tag{5}
\]
The first tie gives the inner target of the first kind with
\(a=2/3\), and the second gives the inner target of the second
kind with \(b=1/3\).

For a \(\delta_K\) target, \(k=0\) in (4) descends to the
\(\Delta'\) target with \(b=(a+1)/3\). At \(k=1,2\), the
\(\lambda\) cap cannot be active, since \(ah<\lambda\).
Thus \(\delta'=(a+1)h-1/2\ge0\), which requires \(h\ge3/10\).
The inner clock is then \(1/2+3h>1\). Consequently
\(\delta'=0\) and both inner states have a uniquely least
constant term. At the current level the even and odd states
in this terminal case have uniquely least indices \(0,1\),
respectively.

For a \(\Delta_K\) target, \(k=0,1\) would require
\(\Delta'=1+(b+1)h\). This contradicts \(\Delta'\le0\) if the
inner clock is at least one. Otherwise \(h<1/6\), and the
required value exceeds \(1/2+3h\), contradicting the bound on
\(\Delta'\). Thus \(k=2\) is necessary, and the cap again
cannot be active since \(b<1\). The new target is
\(\delta'=(b+1)h\), of the first kind with \(a=(b+1)/3\).

Every nonterminal step alternates the kind of target and triples
\(h\), so a solution of (5) must end in the terminal case above.
For a \(\delta_K\) target the two states compared are \(E_K,O_K\),
and for a \(\Delta_K\) target they are \(E_K,O_{K-1}\). In each
nonterminal step their strictly minimizing residue classes are
\(s=0,1\), respectively. For \(\delta_K\) this is the case
\(k=0\) with a positive deficit, and for \(\Delta_K\) it is
\(k=2\) with deficit strictly below the cap in (4).
It follows by reversing the steps that each state has a unique
least index. After \(q\) nonterminal steps those indices are
\(0\) and \((3^{q+1}-1)/2\). The parity of \(q\) is even
for a \(\delta_K\) target and odd for a \(\Delta_K\) target.

We can now compare two Gaussian states of degrees differing by
at most one and with opposite unit parameter residues. At an even
state with \(k=2\) and a tie between classes, the two inner
states satisfy the first target in (5). At an odd state with \(k=0\)
and a tie, they satisfy the second target. Each tie therefore
gives a binomial initial form with odd width. In the even case
the two exponents are \(0,3r+2\), with \(r=(3^{q+1}-1)/2\)
odd. In the odd case they are \(0,3r+1\), with \(r\) even.
Equal degrees with this tie have opposite relative signs on the
two binomial terms, so their nonzero roots are distinct in
characteristic three. At two adjacent degrees, the neighboring
state of a tied even state descends to its odd inner state.
Likewise, the neighboring state of a tied odd state descends to
its even inner state. The neighboring initial form then has a
unique least term by the target analysis.

If there is no tie between residue classes, both states descend
by (3). Their degrees still differ by at most one. Indeed the
only choices other than \(\lfloor(2K+\epsilon)/3\rfloor\)
are the lower of the two branches for an even degree congruent
to four modulo six or an odd degree congruent to one modulo six.
Both neighboring degrees always descend to the higher branch.
Equal degrees have identical choices. The common clock increases
under descent, the two unit parameter residues remain opposite,
and the pair eventually reaches \(\lambda\ge1\), where the
terminal initial forms give no simultaneous cancellation.
For distinct adjacent degrees only one can be odd, and for equal
degrees the possible linear initial forms have opposite signs.
This proves the comparison for all clocks greater than \(1/2\).

The two unit orders of the cubic window are either adjacent or
differ by two. In the latter case they surround the multiple of
three. For an even middle order \(6a\) the unit degrees are
\(6a-1,6a+1\). The lower one descends to the odd state of
degree \(2a-1\), one of the two inner states of the upper one.
For an odd middle order \(6a+3\) the unit degrees are
\(6a+2,6a+4\). The lower one descends to the even state
of degree \(2a\), one of the two inner states of the upper one.
If the upper state has a residue-class tie, the lower one thus
has only a unique least term. Otherwise descent reduces to degrees
differing by at most one. The Gaussian comparison rules out both
placements at \(1/2<\lambda<2/3\), which proves
\(v(A)-v(B)/2=1/4\) in (1).

## Whole families at the remaining clock

The noncollinear gate is regular at every window containing either
of the following unit orders:

- An even order \(j=2K\) whose \(K\) has only digits zero and one in
  base three.
- An odd order \(j=2\cdot3^a+1\), where \(a\ge1\).

It is also regular at all windows whose unique multiple of three is
\(n=3^a\) or \(n=2\cdot3^a\), for \(a\ge1\).

For the first criterion, use the unique residue \(s=0\) at
\(\lambda=1/2\) for an even state with \(k=0,1\). Each step
decreases \(K\) to \(\lfloor K/3\rfloor\), leaves it in an
even state, and preserves the clock \(1/2\). All digits allow
this descent to \(K=0\), so the initial form cannot vanish.
For \(j=2\cdot3^a+1\), take \(K=3^a\), \(\epsilon=1\).
At the clock \(1/2\), the term at \(r\) has valuation
\[
 \frac{s_3(K-r)+s_3(2r+1)-2}{2},
\]
where \(s_3\) is ternary digit sum. This is positive at \(r=K\),
since \(s_3(2K+1)=3\). For \(0<r<K\) the two digit sums are
positive. Equality would require \(K-r=3^u\) and
\(2r+1=3^b\) with \(b>0\), but \(2(K-3^u)+1\) is one
modulo three if \(u>0\) and two modulo three if \(u=0\).
The constant term is therefore uniquely least.

If \(n=3^a\) is not the first order, \(j=n-1\) is an even
unit order whose half has all digits one. If \(n=2\cdot3^a\)
is not the last order, \(j=n+1\) is the odd order of the second
criterion. For \(n=3^a\) first, use \(j=n+2\) with
\(K=(3^a+1)/2\), whose digits are a last digit two preceded
by ones. For \(n=2\cdot3^a\) last, use \(j=n-1\) with
\(K=3^a-1\), whose digits are all two. Both are odd states
with no zero digit in \(K\). At clock \(1/2\) an odd state
with \(k=1,2\) always selects \(s=1\), descends to
\(K_0\), and remains an odd state. Neither can have an initial
zero, since the descent reaches \(K=0\).

## Exact quarter-clock residues

At the remaining clock \(1/2\), choose \(\pi^2=3\) and write
\(x=\pi y\) in \(P_{K,\epsilon,J}(x)\). If the least term
valuation is \(m\), reduce \(\pi^{-2m}P_{K,\epsilon,J}(\pi y)\)
in the residue field. Write that polynomial as \(E_K(J,y)\)
when \(\epsilon=0\) and \(O_K(J,y)\) when \(\epsilon=1\).
The notation in this section refers to residue polynomials, not the
minimum valuations in (4). Start with \(E_0=O_0=1\).
For \(K=3K_0+k>0\), the complete recursion is
\[
\begin{split}
 E_{3K_0+k}(J,y)
  &=E_{K_0}(-J,y^3) &&(k=0,1),\\
 E_{3K_0+2}(J,y)
  &=y^2 O_{K_0}(-J,y^3)
    +\mathbf 1_{K_0\equiv0\ (3)} E_{K_0}(-J,y^3),\\
 O_{3K_0+k}(J,y)
  &=-kJ y\,O_{K_0}(-J,y^3) &&(k=1,2),\\
 O_{3K_0}(J,y)
  &=E_{K_0}(-J,y^3)
    +\mathbf 1_{K_0\equiv2\ (3)} J y\,O_{K_0-1}(-J,y^3).
\end{split}                                                       \tag{6}
\]
Arguments \(J\) on the right are taken as residue units.
In particular the last term is omitted when \(K_0=0\).

To check the branching conditions, evaluate (4) at \(\lambda=\mu=1/2\).
The same comparisons apply at this endpoint, with equality now
allowed in the two indicated classes. The bound \(\Delta_K\le1/2\)
follows by induction on \(K\) from (4). Thus \(\delta_K=0\)
exactly when \(K\equiv0\bmod3\), and \(\Delta_K=1/2\)
exactly when \(K\equiv2\bmod3\). Indeed, in the first row the
three entries are zero, \(\delta'+1/2>0\), and \(1/2\).
In the second row the first two entries are nonpositive and the
third equals \(1/2\). The tie conditions \(\delta'=0\) and
\(\Delta'=1/2\) give the two indicators in (6).

For the units in (6), factor the multiples of three and use the
factorial unit rule in the preceding ternary descent. For \(r=3t+s\)
the fixed residue factor, in addition to the inner parameter \(-J\),
is
\[
 J^s (-1)^{\eta-\ell}
    \frac{k!}{(k-s+3\eta)!(2s+\epsilon-3\ell)!}
 \ \operatorname{unit}_3\!\left(\frac{K_0!}{(K_0-\eta)!}\right).
\]
Here \(\eta,\ell\) are as in (3) and
\(\operatorname{unit}_3(a)=3^{-v(a)}a\).
For an even branch \(s=2\) and \(k=2\), the factorial ratio is
\(2!=-1\) in the residue field, which cancels the sign
\((-1)^{-\ell}\); the factor \(J^2\) equals one. For an odd
branch \(s=1\) and \(k=1,2\), the ratio is \(k\) and the
sign is negative. For the odd additional branch \(s=1,k=0\),
the digit ratio is \(-1\), canceled by the unit of \(K_0\)
since \(K_0\equiv2\bmod3\). This gives (6).

As a consequence, whether a single unit order has a nonzero
quarter-clock residue root is decided entirely by its ternary
digits. Process the digits of \(K\) from least to most significant,
with current state \(E\) or \(O\). In state \(E\), a digit two
switches the state to \(O\), and gives a branch if the next digit
is zero. In state \(O\), a digit zero switches the state to \(E\),
and gives a branch if the next digit is two. A digit one leaves
the state unchanged. Stop when the remaining integer is zero,
treating a missing next digit as zero for the branch test. If no
branch occurs, the residue is a monomial and has no nonzero root.
At the first branch, the two summands in (6) occupy distinct
exponent residues modulo three, so neither can cancel the other
and the result is not a monomial. It has a nonzero root over the
algebraic closure of the residue field. This gives an exact
criterion for one unit order and extends the individual families
above without a test at successive orders.

## Simultaneous first-level cancellation at unbounded starts

First-level three-adic initial forms alone do not exclude all cubic
windows. For every \(t\ge1\), take the window starting at
\(n=12\cdot3^t\). These windows admit simultaneous cancellation
of all three initial forms for suitable valuations and nonzero
residues. This statement concerns initial forms only.

Set \(B=1\), \(v(A)=1/4\), and use \(A^2=\pi y\). The two unit
orders \(n+1,n+2\) are respectively the odd state at
\(K=6\cdot3^t\) and the even state at \(K=6\cdot3^t+1\).
Both descend by (6) to \(E_{2\cdot3^t}\), with opposite unit
parameters. Its trailing zero digits give further descents to
\(E_2(J,y)=1+y^2\), independent of \(J\). Hence both initial
forms are \((1+y^2)^{3^{t+1}}\), and both vanish when
\(\bar y^2=-1\).

At order \(n\), write \(e=v(n)\). The centered polynomial in
\(W=1/C^2\) has degree \(2\cdot3^t\), and its Newton slopes are
in \([e+1,e+2]\). The lower endpoint follows from the centered
coefficient formula, and the backward endpoint bound for a unit
order applies here with \(e\) added to the slopes. Fix any such
slope \(2\gamma\), and take \(v(C)=\gamma\).
Every term with \(a\ge2\) is strictly lowered by the pair
replacement at the multiple of three, and the \(a=1\) terms
are strictly lowered by replacing a singleton and a pair by a
triple, since \(\gamma\le1+e/2\). Consequently the initial
form is the centered edge polynomial. After removing its initial
monomial, it has positive degree and nonzero constant, so it has
a nonzero root in the algebraic closure of the residue field.
Choosing the residue of \(C\) to give that root is independent
of the choice of \(\bar y\). Indeed the two unit orders have
\(\gamma-3/4>0\), so their initial forms are independent of \(C\).
The residue choices and rational valuations can all be realized
in a finite extension of \(\mathbb Q_3\) by adjoining a suitable
root of three and lifting the finite residue field. They yield
simultaneous initial-form zeros for the chosen window, without
asserting any lift to common zeros of the original forms.

## The cubic coefficient clock at the order divisible by three

In the normalization \(B=1\), put \(\alpha=1/4\),
\(\gamma=v(C)\), and \(e=v(n)\) as in (1). If \(n\) is even,
any zero has \(\gamma<1+e/2\). If \(n\) is odd, a zero
requires either \(\gamma\le1+e/2\), or the separate resonance
\[
       \gamma=\alpha+e+1,\qquad
       \overline{\frac{(n-1)C}{6nA}}=-1.                       \tag{7}
\]
In every case in the ordinary clock band
\([(1+e)/2,1+e/2]\), the value \(2\gamma\) must be a slope of
the centered coefficient polygon at \(n\). The lower endpoint
of that band is possible only if \(n/3\equiv2\bmod3\).

At \(n\), pair replacement leaves only \(a=0,1\) as candidates
for least valuation. If \(n\) is even, compare to the pure pair
term. In (2) the excess valuation is at least
\[
       a(\alpha+e/2)+c(\gamma-1-e/2).
\]
Here the removed pair factorial strip has length
\((a+3c)/2\ge c\), and hence contributes at least \(v(c!)\).
For \(\gamma>1+e/2\) the pure pair term is uniquely least.

Whenever \(\gamma<\alpha+e+1\), replacing a singleton and a
pair by a triple strictly lowers an \(a=1\) term, since \(b\)
is a three-adic unit. Thus the centered polygon must have slope
\(2\gamma\) in this range. Its upper bound \(e+2\) is strict
when \(n\) is even, because the backward strip bound is strict
for \(s=(n/3)\bmod2=0\). Its lower slope bound \(e+1\) can be
attained only at the first edge through \(k=1\), where the
centered \(c_0=n/3\) must satisfy \(c_0\equiv2\bmod3\).

It remains to consider odd \(n\) with \(\gamma\ge\alpha+e+1\).
Compare to the term \(a=1,c=0\), whose pair count is
\(b_0=(n-1)/2\). Within \(a=1\), adding \(c>0\) triples
raises (2) by at least \(c(\gamma-1-e/2)>0\).
Within \(a=0\), the excess over the same reference is at least
\[
       c(\gamma-1-e/2)-(\alpha+e/2).
\]
Here \(c\) is odd and the removed pair strip has length
\((3c-1)/2\ge c\). This excess is positive unless
\(\gamma=\alpha+e+1\) and \(c=1\), when \(b_0\) is a
three-adic unit and the two terms have equal valuation.
The ratio of the \(a=0,c=1\) term to the \(a=1,c=0\) term
is \(b_0 C/(3nA)\). Their residue cancellation gives (7).
