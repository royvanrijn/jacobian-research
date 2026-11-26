# All-degree closure for monotone gap counts in an upper cluster

Fix an ordinary two-monomial profile and remove its equal rows.
Select a nonempty group of unequal rows. Let \(D\) be its
largest absolute gap, and let \(M\) be its largest exponent.
Let \(J\ge4\) be any integer at least the
[absolute angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md).
At each upper exponent \(M-e\), and separately for each
orientation present at that exponent, let \(c_d\) be the
multiplicity of gap \(d\). Write \(D'\) for the largest gap
in that height and orientation. Assume
\[
       0\le e\le D/(24J),\qquad
       D/2\le D'\le D,\qquad
                 0\le c_1\le c_2\le\cdots\le c_{D'} .
                                                               \tag{1}
\]
Zeros at the start of the sequence are allowed. There is no
bound on the number of rows or distinct gap values.

For the remaining rows put
\[
 \mathcal L=\sum_{\substack{i\ {\rm retained}\\JM_i>2J^2D}}
                     d_i\log(JM_i),
       \qquad M_i=\max(a_i,b_i),\quad d_i=|b_i-a_i|.
\]

**Theorem.** There is a computable \(D_0(J)\) such that
\[
        D\ge D_0(J),\qquad
        \mathcal L\le\frac{D\log M}{1024J^4}
                                                               \tag{2}
\]
imply physical SFC(2) in every degree over every
characteristic-zero field. The selected group may have many
different upper exponents, a least gap of one, and
independent monotone multiplicities in the two orientations.

For a fixed ordinary core, the log-height condition in (2)
holds automatically for all sufficiently large \(D\),
independently of the selected rows and counts.

## The exact literal count

The [curvature twelve theorem](MULTIVARIATE_ABSOLUTE_CURVATURE_12.md)
allows us to assume \(S<12\), and hence \(M>D^2/12\) by the
selected row of gap \(D\). Take a terminal
\(4\le m<J\), and put \(m'=m-1\).
Strip each row's smaller-end factorial at each order.

For integers
\[
                  D/8\le u\le D/4,\qquad u\equiv1\pmod m ,
                                                               \tag{3}
\]
consider the literal factor \(x_u=mM-u\).
At a height \(M-e\) write \(u_e=u-me\), so
\(D/12\le u_e\le D/4\le D'/2\).
For the given orientation define
\[
 C(v)=\sum_{1\le d\le v} c_d,\qquad
 v_j=\lfloor u_e/j\rfloor\quad(1\le j\le m),\qquad
                          v_0=D'.
\]
The occurrence condition for an increasing row in the
coefficient of index \(k\) is
\((m-k)d\le u_e<md\). Thus its orientation contributes
\(C(v_{m-k})-C(v_m)\) copies. For a decreasing row the
condition is \(kd\le u_e<md\), and that orientation
contributes \(C(v_k)-C(v_m)\) copies. These formulas
also hold at both endpoints, where one bound on \(d\)
is absent and \(v_0=D'\) counts all available gaps.

The sequence \(v_j\) is decreasing and strictly discretely
convex when \(D\) is sufficiently large in terms of \(J\).
At \(j=1\), its second difference is at least
\[
            D'-\tfrac32u_e-1\ge D/8-1.
\]
For \(2\le j<m\) it is bounded below by
\[
         \frac{2u_e}{j(j^2-1)}-2\ \ge\ \frac{D}{6J^3}-2 .
\]
Both estimates use the exact endpoint definition and the
floor error at the two outer arguments.

As the increments of \(C\) are nonnegative and
nondecreasing, \(C(v_j)\) is discretely convex.
In fact, if \(a=v_{j-1}-v_j\) and \(b=v_j-v_{j+1}\),
then \(a>b\ge0\). The left increment of \(C\) is a sum
of \(a\) counts at larger indices than the \(b\) counts
in the right increment. The first sum is at least
the second. It is strictly larger if the counts are
positive on that part of the grid. Reversal preserves
discrete convexity.

Suppose some nonempty height and orientation has first
positive gap at most \(D/(12J)\). Its counts are positive
from that gap to its ceiling. Since
\(v_m\ge\lfloor D/(12J)\rfloor\), this component gives
strict convexity at every internal index. Thus the total
literal count \(A_u(k)\) of \(x_u\) in the selected group
is strictly discretely convex. Its minimizing indices
are either one index or two adjacent indices. Denote
their interval by \(I_u\); it has width at most one.
For any index outside \(I_u\), the count rises by at
least one.

If every selected gap is larger than \(D/(12J)\), apply
the [large factorial mark theorem](MULTIVARIATE_LARGE_FACTORIAL_MARK_PRIME_BUDGET.md)
with minimum gap \(d_*>D/(12J)\), \(A=24J\), and \(B=1\).
Its cutoff \(J^2(A+B)d_*\) exceeds \(2J^2D\), and
the second bound in (2) implies its retained log
budget \(d_*\log M/(64J^3)\). The height condition
in that theorem holds for all sufficiently large \(D\).
This handles that case.

## The prime budget for the complete interval

Let \(P=2J^2D\). For a prime \(p>P\) dividing \(x_u\),
no other selected factor at order \(m\) can be
divisible by \(p\): its difference from \(x_u\) is a
nonzero integer of absolute value less than \(P\).
At the adjacent order, writing a rising factor in top
coordinates as \(r=m'(M-e)-h\), one has
\[
                mr-m'x_u=m'u-mm'e-mh .
\]
This is nonzero because \(u\equiv1\pmod m\) and
\(\gcd(m,m')=1\). Its absolute value is less than \(P\)
by the strip and gap bounds. Hence all selected
factors at order \(m'\) are units at \(p\).
The binomial coefficients are also units since \(p>J\).

For \(j=m,m'\) put \(B_j\) equal to \(j!\) times the
product of the retained stripped ramps at order \(j\).
Let \(E=v_p(x_u)\), \(U=v_p(B_m)\), and \(L'=v_p(B_{m'})\).
The marked ordinates are \(EA_u(k)\) plus a term
in \([0,U]\), and all adjacent slopes have magnitude
at most \(L'\). If \(E>U+mL'\), subtracting any
available slope times its index makes every index
outside \(I_u\) strictly larger than its nearest
index in \(I_u\). Its excess before the retained
and slope contributions is at least \(E\).
Consequently every shared slope at the marked order
has total width at most one. The
[Newton degree comparison](MULTIVARIATE_CENTRAL_DEEP_PRIME_DIVISOR_GATE.md#marks-in-arbitrary-blocks)
then gives common degree at most one, impossible
for a nonconstant common factor because the even
parent has no real root.

A possible common root therefore implies
\[
            (x_u)_{>P}\mid (B_mB_{m'}^m)_{>P}
                 \quad\text{at each tested }u .
\]
A prime greater than \(P\) divides at most one
mark, so the large-prime part of the product of all
marks divides \((B_mB_{m'}^m)_{>P}\).
For each prime at most \(P\), remove one mark of
maximal valuation if such a mark exists. The marks
are coprime to \(m\), and the progression index of a
mark divisible by \(p^b\) is in a single class
modulo \(p^b\). After the removal of a maximal mark
there are at most \(\lfloor(k_0-1)/p^b\rfloor\) terms
at that depth, where \(k_0\) is the original number
of marks. Thus the small-prime part of the remaining
product divides \((k_0-1)!\). This is the product
argument of the
[large factorial mark theorem](MULTIVARIATE_LARGE_FACTORIAL_MARK_PRIME_BUDGET.md#a-product-bound-at-a-possible-common-root),
applied to the new strictly convex interval.

For \(D\ge16J\), (3) gives \(D/(16J)\le k_0\le D\).
The elementary estimate
\(\pi(P)\le9P/\log P\) gives
\(\pi(P)\le k_0/4\) when \(\log D\ge1152J^3\).
All marks are greater than \(M\), and for \(D>144\)
we have \(\log k_0\le(2/3)\log M\).
The count of retained factors contributing at
primes greater than \(P\) gives
\[
   \log (B_mB_{m'}^m)_{>P}
      \le J^2\mathcal L \le \frac{k_0\log M}{32}.
\]
The remaining product thus has logarithm at least
\((3/4)k_0\log M\), but its divisibility upper bound
has logarithm at most
\((2/3+1/32)k_0\log M\), a contradiction.

The first two adjacent pairs use
[ordinary low-order coprimality](MULTIVARIATE_TWO_MONOMIAL_NONDIVISIBILITY.md),
and the absolute cutoff handles higher orders.
Rational Bézout identities give the characteristic-zero
assertion. For a fixed core, (2) holds with
\(\mathcal L=0\) for all sufficiently large \(D\).

## A nonmirrored family with a full interval of gaps

Take heights \(M-e\) for
\(0\le e\le\lfloor D/(24J)\rfloor\). At each height and
every \(1\le d\le D\), take \(2D+1\) increasing rows of gap
\(d\) and \(3d\) decreasing rows of gap \(d\). Both count
sequences obey (1). The signed gap of each height bucket
is zero since
\[
       (2D+1)\sum_{d=1}^D d
                 =3\sum_{d=1}^D d^2 .
\]
Let \(A_D\) be the sum of squared gaps of all these rows and
take \(M=\lfloor A_D/3\rfloor\). The total curvature tends
to three, the largest relative gap tends to zero, and the
number of different heights grows without bound.
The theorem gives all-degree SFC(2) for these original
profiles at all sufficiently large \(D\).
