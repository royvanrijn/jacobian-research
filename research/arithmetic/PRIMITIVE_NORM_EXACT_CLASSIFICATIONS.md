# Exact bounded primitive norm sumsets and global exceptions

For an integral positive-definite binary form, write
\[
 P_Q=\{Q(x,y)>0:(x,y)\in\mathbb Z^2,\ \gcd(x,y)=1\}.
\]
Zero coordinates are permitted; in particular \(1\in P_Q\) for all forms
below. The two summands are positive, so targets start at 2.

## What is established

The discriminant \(-11\) identity is now proved exactly:
\[
 P_{-11}+P_{-11}=\{N\ge2:2\mid N,\ 121\nmid N\}.
\]
The [complete quaternionic and finite proof](PRIMITIVE_NORM_MINUS11_EFFECTIVE.md)
first proves the stronger core identity \(C+C=\{N\ge2:N\text{ even},11\nmid N\}\).
It combines a fully explicit large-core sieve, 35 certified intervals
from \(2^{25}\) to \(2^{60}\), and exhaustive coverage of all 15252015
eligible targets through \(2^{25}\). The value-pair lower bound
\(N^{2/3}/100000\) remains explicit for eligible \(N\ge2^{71}\).
This is an internally reviewed deduction from classical inputs; priority
and external review have not been established.

The discriminant \(-8\) identity is also proved exactly:
\[
 P_{-8}+P_{-8}=\{N\ge2:16\nmid N\}.
\]
The [complete -8 proof](PRIMITIVE_NORM_MINUS8_EXACT.md) proves the two
stronger core identities \(C+C=\{N\ge2:N\text{ even},8\nmid N\}\) and
\(C+2C=\{N\ge3:N\text{ odd}\}\). It combines 107 complete rational
interval certificates with finite coverage through \(2^{27}\), followed
by split-prime descent including its separately proved prime boundary.
Its explicit value-pair bound \(N^{2/3}/4800\) holds for eligible
\(N\ge2^{78}\). Both exact binary identities have empty exceptional sets.

The general [quaternary deduction](PRIMITIVE_NORM_ANALYTIC_BRIDGE.md)
also proves eventual local-global existence for arbitrary fixed binary
form pairs. This broader eventual theorem and the three exact classifications
have different scopes.

The discriminant \(-7\) classification is now also proved exactly:
\[
 P_{-7}+P_{-7}=\{N\ge2:49\nmid N\}\setminus\{7,77,203\}.
\]
The [complete -7 proof](PRIMITIVE_NORM_MINUS7_EXACT.md) proves the
stronger core identity \(C+C=\{n\ge2:7\nmid n\}\setminus\{11,29\}\).
It joins full finite coverage through \(2^{24}\), 257 complete rational
intervals, 17 higher dyadic bins and the infinite core/prime bounds.
Exceptional-quotient descent and ramified scaling close every target.
Its value-pair bound \(N^{2/3}/20000\) holds for eligible \(N\ge2^{69}\).

Exact local-global equality also fails for \(D=-19\): 4 is locally
admissible for \(x^2+xy+5y^2\), but its only primitive positive value
at most 4 is 1. The remaining global exception lists below are open.

## Reproducible bounded classification

For \(D=-d\) with \(d\in\{7,11,19,43,67,163\}\), use
\(Q_D=x^2+xy+(d+1)y^2/4\), and for \(D=-8\) use
\(Q_D=x^2+2y^2\). Define \(\mathcal L_D\) by
\[
\begin{array}{c|c}
 D&\text{pair-primitive local criterion}\\\hline
 -8&16\nmid N\\
 -7&49\nmid N\\
 -11,-19,-43,-67,-163&2\mid N,\quad D^2\nmid N.
\end{array}
\]
The local proofs are in [PL3 of the local analysis](PRIMITIVE_NORM_LOCAL_STRUCTURE.md); the following
computation tests global representation against those explicitly stated
sets.

**Exact finite result.** For every \(D\) in the table and every integer
\(2\le N\le10^6\),
\[
 N\in P_D+P_D\quad\Longleftrightarrow\quad
 N\in\mathcal L_D\setminus E_D,
\]
where the complete bounded exception sets are retained in
[primitive_norm_exact_sumset_checks.json](primitive_norm_exact_sumset_checks.json).
Their sizes and largest elements are:

| Discriminant | Number of exceptions | Largest exception | Complete set when short |
|---|---:|---:|---|
| -8 | 0 | none | empty |
| -11 | 0 | none | empty |
| -7 | 3 | 203 | 7, 77, 203 |
| -19 | 3 | 380 | 4, 76, 380 |
| -43 | 18 | 3182 | see retained list |
| -67 | 38 | 9246 | see retained list |
| -163 | 129 | 161696 | see retained list |

The assertions about individual exceptions remain true without the
upper-bound qualification: e.g. \(380\notin P_{-19}+P_{-19}\).
For \(-19,-43,-67,-163\), the assertion that the lists have no further
members is bounded by \(10^6\); their maxima are not proved global cutoffs.
The \(-7\) exception list is globally complete by the separate proof above.

The checker is
[verify_primitive_norm_exact_sumsets.py](verify_primitive_norm_exact_sumsets.py).
Run from the repository root with Python 3 (standard library only):

```sh
python3 research/arithmetic/verify_primitive_norm_exact_sumsets.py \
  --bound 1000000 \
  --output research/arithmetic/primitive_norm_exact_sumset_checks.json
```

The recorded internal execution completed successfully under Python
3.12.3. It supplies exact finite computation, without formal verification,
independent external review, or a proof beyond the stated range.

## Coverage argument for the computation

For \(D=-8\), independent coordinate signs permit \(x,y\ge0\);
the enumeration uses \(0\le y\le\lfloor\sqrt{B/2}\rfloor\) and
\(0\le x\le\lfloor\sqrt{B-2y^2}\rfloor\), retaining gcd 1.

For odd \(d\), the identity
\[
 4Q_D(x,y)=(2x+y)^2+dy^2
\]
gives complete bounds. Simultaneous sign reversal makes \(y\ge0\),
and reflection \(x\mapsto-x-y\) makes \(2x+y\ge0\).
Both preserve the gcd and the represented value. Thus it is enough to
use
\[
 0\le y\le\lfloor\sqrt{4B/d}\rfloor,\qquad
 -\lfloor y/2\rfloor\le x\le
 \left\lfloor\frac{\lfloor\sqrt{4B-dy^2}\rfloor-y}{2}\right\rfloor.
\]
All square roots are computed by exact integer square root, and all
values and gcd tests use integers. Consequently the resulting flag
array is exactly \(P_D\cap[1,B]\).

A packed integer bitset marks all sums \(a+b\le B\) with
\(a\in P_D\cap[1,5000]\) and \(b\in P_D\cap[1,B]\).
Every marked target therefore has a genuine positive representation.
For each still-unmarked locally admissible target, the checker tests
all \(a\in P_D\cap[1,N/2]\) against the complete value flag array.
If none works, it records an exception. This last pass establishes
nonrepresentation rather than assuming that a short summand must exist.
The checker also verifies that the marked sums violate none of the
stated local criteria, and checks the primitive norm value arrays
through 16 against hand lists.

## Remaining exact obligations

The -7, -8 and -11 classifications are proved in their linked canonical sources.
For -19,-43,-67,-163 the bounded exception lists below remain globally
unclassified: an effective cutoff and complete finite coverage are still
needed. Increasing a scan bound without a completion argument does not
prove that those lists have no further members.


## An exact reduction with one coordinate equal to 1

Let \(Q(x,y)=x^2+2y^2\). The following elementary statement explains
what an ordinary ternary local-global argument supplies.

**Proposition.** For every integer \(N\ge2\), the following are equivalent:

1. \(16\nmid N\).
2. \(N=Q(x,y)+Q(u,v)\) with both norms positive and at least one of
   \(x,y,u,v\) equal to \(1\) or \(-1\).

In particular every target in the proposed \(-8\) sumset has a
representation with at least one binary pair primitive. The proposition
does not assert primitivity of the other pair.

**Proof.** Write
\[
 A(x,u,v)=x^2+u^2+2v^2,\qquad
 B(u,y,v)=u^2+2y^2+2v^2.
\]
The Gauss--Legendre three-square theorem says that a positive integer
is a sum of three integral squares exactly when it is not of the form
\(4^a(8b+7)\), with \(a,b\ge0\). A proof and precise statement are in
Paul Pollack, [*Not Always Buried Deep*, Theorem 4.9 and Lemma 4.14,
printed pages 134 and 138](https://www.pollack-math.net/NABDofficial.pdf).

The identity
\[
 2A(x,u,v)=(x+u)^2+(x-u)^2+(2v)^2
\]
also gives the converse: in any three-square representation of \(2n\),
choose an even coordinate for \(2v\). Such a coordinate exists; when
\(n\) is even all three coordinates are even, while when \(n\) is odd
exactly one is even. The other two coordinates have the same parity,
so their half-sum and half-difference recover \(x,u\). Therefore
\[
 n\in A(\mathbb Z^3)
 \quad\Longleftrightarrow\quad
 n\ne4^a(16b+14)\quad(a,b\ge0)
 \tag{E1}
\]
for positive \(n\).

Similarly,
\[
 B(u,y,v)=u^2+(y+v)^2+(y-v)^2.
\]
In the converse direction choose two coordinates of a three-square
representation having the same parity, and take their half-sum and
half-difference. Hence
\[
 n\in B(\mathbb Z^3)
 \quad\Longleftrightarrow\quad
 n\ne4^a(8b+7)\quad(a,b\ge0).
 \tag{E2}
\]

For \(N>2\), first attempt
\[
 N=Q(x,1)+Q(u,v),\qquad A(x,u,v)=N-2.
\]
If this fails, (E1) gives \(N-2=4^a(16b+14)\). Under \(16\nmid N\),
one must have \(a\ge1\), so \(N\equiv2\pmod4\). Consequently
\(N-1\equiv1\pmod4\), which is excluded from every exceptional family
in (E2). We may then write
\[
 N=Q(1,y)+Q(u,v),\qquad B(u,y,v)=N-1.
\]
In the first construction interchange \(x,u\), if necessary, to ensure
\((u,v)\ne(0,0)\). In the second interchange \(y,v\), if necessary,
to achieve the same. The relevant ternary value is positive, so these
interchanges always work. For \(N=2\), use \(Q(1,0)+Q(1,0)\).

Conversely, if \(16\mid N\), then \(N-2=16b+14\) is excluded by
(E1) and \(N-1=8c+7\) is excluded by (E2). By sign changes and
interchanging the two binary blocks, any coordinate equal to \(\pm1\)
would give one of these two impossible representations. This proves
the equivalence. \(\square\)

The missing coprimality does not follow from the ternary theorem. For
example, \(29=A(x,u,v)\) has only \(v=0\) and
\(\{|x|,|u|\}=\{2,5\}\). To verify this directly, test
\(|v|=0,1,2,3\); the corresponding sums of two squares are
\(29,27,21,11\), and only the first is represented. Thus
\[
 29\notin\{A(x,u,v):\gcd(u,v)=1\},
\]
even though 29 is represented by the ordinary ternary form. In
particular the construction with \(y=1\) fails to give two primitive
blocks at \(N=31\). The alternative unit-coordinate construction does
work there: \(31=Q(1,2)+Q(2,3)=9+22\).

This identifies the precise remaining step in this elementary route:
control the primitive second block in at least one of the two ternary
representations, rather than infer it from ordinary ternary regularity.


### The primitive ternary proposal contains an original binary subproblem

There is a precise limitation to replacing the binary problem by an exact
classification of \(A\) with a primitive block. Put
\[
 C=\{n\ge1:p\mid n\Longrightarrow p\equiv1,3\pmod8\},\qquad
 \mathcal A=\{x^2+u^2+2v^2:\gcd(u,v)=1\}.
\]
As throughout, a primitive binary block has positive norm; zero norms
are excluded from this definition. Then
\[
 P_{-8}=C\cup2C.
 \tag{E3}
\]
Here is the class-number-one input in (E3). In
\(\mathbb Z[\sqrt{-2}]\), rounding the two coefficients of a field
element to the nearest integers leaves norm at most
\(1/4+2/4=3/4<1\). Thus the norm is Euclidean and the ring is a unique
factorization domain. An odd rational prime splits precisely in classes
1 and 3 modulo 8, is inert in classes 5 and 7, and
\(2=-(\sqrt{-2})^2\) ramifies. An inert prime dividing an element norm
forces both rational coordinates divisible by that prime. At 2 a
primitive norm has valuation zero or one. Conversely, at every split
prime choose only one of its two conjugate prime factors, to the full
required exponent, and use the factor \(\sqrt{-2}\) at most once. Their
product has the prescribed norm and cannot be divisible by a rational
prime. This proves (E3). In particular \(C\) is closed under taking
positive divisors and under multiplication.

**Proposition.** For every positive integer \(t\),
\[
 t^2\in\mathcal A\quad\Longleftrightarrow\quad 2t\in C+C.
 \tag{E4}
\]
If \(t\) is odd, this is also equivalent to
\(2t\in P_{-8}+P_{-8}\).

**Proof.** A representation in \(\mathcal A\) gives
\[
 0<t^2-x^2=(t-x)(t+x)\in P_{-8}.
\]
Both factors are positive and have the same parity. They cannot both
be even, because a primitive norm in (E3) is never divisible by 4.
They are therefore odd, and their product lies in \(C\). Divisor
closure gives \(t-x,t+x\in C\), and their sum is \(2t\).
Conversely, if \(a+b=2t\) with \(a,b\in C\), both are odd. Put
\(x=(a-b)/2\). The identity \(t^2-x^2=ab\in C\subset P_{-8}\)
gives the desired primitive ternary representation.

Finally, for odd \(t\), the target \(2t\) is 2 modulo 4. Two summands
in \(P_{-8}\) with even sum have the same parity. They cannot both
be even: each even primitive norm is 2 modulo 4. Thus both are odd,
so (E3) makes them elements of \(C\). \(\square\)

Consequently an exact theorem asserting that all odd squares belong to
\(\mathcal A\) already proves the original \(-8\) identity for every
target \(N\equiv2\pmod4\). This is an exact reduction for this particular
primitive ternary proposal, not a limitation on arbitrary quaternary
methods. Ordinary ternary regularity leaves this binary sumset
subproblem intact. The companion checker independently tests (E3)
through 10000 and (E4) for \(1\le t\le100\), using trial prime
factorization on one side and gcd-coordinate norm enumeration on the
other. The proof above, rather than those finite tests, establishes the
unrestricted equivalence.
