# Radical compression for exact powers on the two sides of a difference

If the target and one summand are exact powers, the difference
splits into factors on disjoint prime supports, apart from
bounded primes. Several pairs of factors occur together in
primitive additive rows that recover the original two roots.
The half-price estimate on two selected letters then gives
a bound in terms of the radical of the original difference.

## Statement and counting input

Fix an integer \(m\ge1\), and choose one of the orders
\[
                     M=2^m \quad\text{or}\quad M=3\cdot2^m .
\]
Count positive primitive triples
\[
               y^M+h=x^M=c\le X,\qquad x>y\ge1 ,
\]
with \(\operatorname{rad}(xyh)<c^\lambda\). Interchanging
the original summands accounts for at most a factor of two.
Allow an additional cap \(\operatorname{rad}(h)\le c^{s\lambda}\)
with fixed \(0\le s\le1\). Define
\[
       k_M=\begin{cases}
                    m+1,&M=2^m,\\
                    2m+2,&M=3\cdot2^m.
                   \end{cases}
\]

**RC1.** For the indicated fixed \(m,s\), every \(\eta>0\),
and all sufficiently small fixed \(\lambda>0\),
the count is
\[
 O_{m,s,\eta,\lambda,\varepsilon}
 \left(X^{(\min\{s/k_M,,1/(k_M+2)\}+\eta)\lambda+\varepsilon}\right)
                              \qquad(\varepsilon>0).
 \tag{RC1}
\]
The unrestricted version on this sector takes \(s=1\).
It gives coefficients \(1/4,1/5,1/6\) for \(M=2,4,8\),
and \(1/6,1/8\) for \(M=6,12\), respectively. The same
statements apply when the individual literal orders of the
target and summand are multiples of the displayed fixed
\(M\). Their powers can first be regrouped to order \(M\).

Use a dyadic shell \(Y/2<c\le Y\) and put \(n=Y^\lambda\).
All raw masses in the proof are logarithms of radicals in
base \(n\). Bin the finitely many masses to a fixed small
mesh. Bounded factors and repetitions at fixed primes
then cost arbitrarily small normalized allowances at the
height limit. Let \(D\) denote the normalized exponent of a
single binned family, ignoring these allowances.

The required counting input is MCB2 with singleton valuation
coordinates. For a positive primitive additive row of height
between two fixed positive powers of \(Y\), and with all
three raw radical caps bounded in base \(n\), it bounds the
normalized exponent by half the sum of the costs of any
two specified roles, with an arbitrarily small allowance
at sufficiently small fixed \(\lambda\). The singleton
anatomic boxes cost at most their raw radical caps. A
bounded common gcd can first be divided out and recorded
among boundedly many integers. If the derived row has
bounded multiplicity as a map from the original family,
its bound also counts the original family. [The half-price estimate (MCB2)](ABC_MULTIPLE_SHORT_COFACTOR_BOUNDARY.md#bounded-group-dictionaries-with-arbitrary-fixed-coefficients)

Thus whenever two quantities \(A,B\) are the designated
roles of such a derived row, its bound supplies
\(2D\le \rho(A)+\rho(B)\) in limiting exponents, where
\(\rho(T)=\log_n\operatorname{rad}(T)\). The designated
roles may be fixed powers of \(A,B\) and may acquire
bounded coefficients or a bounded common divisor. The
third role must also have bounded raw radical mass, even
though its mass does not occur on the right of this inequality.

## The factors at the powers of two

Put
\[
 U=x-y,\qquad S_j=x^{2^j}+y^{2^j}\quad(0\le j\le m-1).
\]
All \(m+1\) quantities divide \(h\). The pair \((U,S_0)\)
occurs in \(U+2y=S_0\). Its gcd is at most two since
\(\gcd(x,y)=1\). The two values determine \(x,y\).

For \(0\le j\le m-2\), put \(a=x^{2^j}\) and
\(b=y^{2^j}\). The row
\[
                     S_j^2=S_{j+1}+2ab
\]
uses the pair \((S_j,S_{j+1})\). Its common gcd divides
two because \(\gcd(S_j,ab)=1\). The two entries determine
\(a+b=S_j\) and \(ab=(S_j^2-S_{j+1})/2\), and hence
\(x,y\) in their original order. If \(m\ge2\), the row
\[
                        S_1=U^2+2xy
\]
also gives the pair \((U,S_1)\) and has common gcd at
most two. Here \(x-y=U\) and \(xy=(S_1-U^2)/2\) determine
\(x,y\). Every third role in these rows has bounded raw
mass, since its primes divide \(2xy\). The largest value
in each row is at least \(x\), up to a bounded factor
after the primitive division.

These inequalities give a bound on the sum of the costs
of \(U,S_0,\ldots,S_{m-1}\). If \(m+1\) is even, use
\((U,S_0),(S_1,S_2),\ldots\), each with weight one. If
\(m+1\) is odd, use weight \(1/2\) on each of the three
pairs among \(U,S_0,S_1\). Pair the remaining entries
\((S_2,S_3),\ldots\) with weight one. Every entry gets
total weight one, so summing the pair inequalities gives
\[
                   (m+1)D\le \rho(U)+\sum_{j=0}^{m-1}\rho(S_j).
 \tag{RC2}
\]

## The additional cyclotomic factors when the order contains three

Suppose \(M=3\cdot2^m\). Add the \(m+1\) factors
\[
 R^+=x^2+xy+y^2,\qquad
 R_j=x^{2^{j+1}}-(xy)^{2^j}+y^{2^{j+1}}
                  \quad(0\le j\le m-1).
\]
These are \(\Phi_3(x,y)\) and \(\Phi_{3\cdot 2^{j+1}}(x,y)\),
respectively, with homogeneous cyclotomic polynomials.
For example, a primitive \(3\cdot 2^{j+1}\)-th root raised
to \(2^j\) is a primitive sixth root, giving the formula
for \(R_j\) by degree and monicity. The other group
contains \(\Phi_1(x,y)\) and the \(\Phi_{2^{j+1}}(x,y)\).
For either choice of \(M\), the groups used are exactly the
factors of \(x^M-y^M\).

The row \(R^+=R_0+2xy\) gives the pair \((R^+,R_0)\).
Both factors are odd for coprime \(x,y\), and each is
coprime to \(xy\), so the row is primitive. The sum of
the factors determines \(x^2+y^2\), and their difference
determines \(xy\). This recovers the ordered positive roots.

For \(0\le j\le m-2\), again put \(a=x^{2^j}, b=y^{2^j}\).
Use
\[
             R_{j+1}=R_j^2+2ab(a-b)^2
 \tag{RC3}
\]
for the pair \((R_j,R_{j+1})\). Expansion gives the
identity. The third role uses only primes of \(2xy\) and
of \(a-b\). The latter divides \(h\), since \(2^j\mid M\).
The row is primitive since \(R_j\) is odd and coprime to
\(ab(a-b)\). For the recovery, write \(P=ab\). From
the two designated entries,
\[
              P(R_j-P)=(R_{j+1}-R_j^2)/2 .
\]
There are at most two possibilities for \(P\). Then
\(a^2+b^2=R_j+P\) determines \(a,b\) and the original
roots up to this bounded multiplicity.

If \(m\ge2\), the row
\[
                     (R^+)^2=R_1+2xy(x+y)^2
\]
gives the additional pair \((R^+,R_1)\).
Its third role uses primes of \(2xyS_0\). The row is
primitive since \(R^+\) is odd and coprime to \(xy(x+y)\).
With \(P=xy\), the two designated entries give
\[
                    P(R^++P)=((R^+)^2-R_1)/2 .
\]
This again recovers the positive roots in bounded multiplicity.
The largest entry in each of these rows has numerical height
at least a fixed power of \(Y\). For instance,
\(a^2-ab+b^2\ge3a^2/4\) for \(a\ge b>0\).
All entries are of bounded polynomial height in \(Y\).

There is now the same pair arrangement on
\(R^+,R_0,R_1,\ldots,R_{m-1}\). Use disjoint pairs starting
with \((R^+,R_0)\) if the number of entries is even.
If it is odd, use the half weights on the three pairs
among \(R^+,R_0,R_1\) and pair the remaining consecutive
entries. This gives
\[
                   (m+1)D\le\rho(R^+)+\sum_{j=0}^{m-1}\rho(R_j).
 \tag{RC4}
\]

## Combining the radical caps

The shift factors are pairwise coprime outside primes dividing
\(M\). Indeed if \(p\mid h\), then both roots are units at
\(p\). For \(p\nmid M\), the polynomial \(T^M-1\) is
squarefree over \(\mathbb F_p\), so the common ratio \(x/y\)
can be a root of at most one of its reduced cyclotomic factors.
Consequently the sum of the masses of all the shift factors
is at most \(\rho(h)+o(1)\) at the height limit, uniformly
in the shell. Equations (RC2) and (RC4) give
\(k_MD\le\rho(h)+o(1)\).

Also apply MCB2 to the two roles \(x^M,y^M\) in the original
row. This gives \(2D\le\rho(x)+\rho(y)\), up to the same
type of allowance. Since \(x,y,h\) are pairwise coprime,
adding the last two estimates gives
\[
          (k_M+2)D\le \rho(x)+\rho(y)+\rho(h)+o(1)\le1+o(1).
\]
There are finitely many derived rows for a fixed \(m\). Each
has a fixed positive numerical height floor in the shell
and bounded raw caps. Choose the common small-budget
threshold, count allowances and cap mesh for this finite
family before taking the height limit. The resulting
finite-shell pair bounds and their sums give (RC1), with
the shell sum absorbed by the stated height allowance.

This argument uses the rational factorization of the original
power difference. For two weighted powers with a nonsplit
ratio of coefficients the analogous factors need not be
rational. A decomposition over an extension does not give
the same radical accounting on rational additive rows.

## Arbitrarily many groups indexed by powers of three

**RC5.** Fix \(m\ge1\) and \(b\ge0\), and take
\(M=2^m3^b\). Under the same count and cap hypotheses as
RC1, the same estimate holds with
\[
                        k_M=(m+1)(b+1).
 \tag{RC5}
\]
For example \(M=36\) gives an unrestricted sector
coefficient of \(1/11\).

Keep the group in RC2. For each \(0\le j<b\), apply
the additional group of RC4 with
\[
               X_j=x^{3^j},\qquad Y_j=y^{3^j}
\]
in place of \(x,y\). The resulting \(m+1\) factors
have orders \(3^{j+1},2\cdot3^{j+1},\ldots,2^m3^{j+1}\).
Indeed the plus quadratic at these bases is
\(\Phi_{3^{j+1}}(x,y)\), and the successive minus
quadratics are \(\Phi_{2^{i+1}3^{j+1}}(x,y)\) for
\(0\le i<m\), by the same primitive-root and degree
argument. The factors over all groups have distinct
indices and account for the whole shift.

The third roles in the scaled identities continue to
use only primes of the original roots, the original
shift and a bounded integer. For a minus edge, the
difference on the scaled bases is
\(x^{2^i3^j}-y^{2^i3^j}\), which divides the original
shift. For the additional triangle edge,
\(X_j+Y_j=x^{3^j}+y^{3^j}\) also divides the original
shift since \(2\cdot3^j\mid M\). Coprimality of the
scaled bases, the height floor and bounded recovery
hold just as for RC4. Thus every added group yields
its own \((m+1)D\) inequality. Summing the \(b+1\)
groups and then using the original root pair proves
the assertion.
