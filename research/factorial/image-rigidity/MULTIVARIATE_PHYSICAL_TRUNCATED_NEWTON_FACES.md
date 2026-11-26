# Physical comparisons of truncated Newton faces

Consider the ordinary adjacent factorial moments \(P_n,P_N\) with
\(N=n+1\) and \(n\ge2\), in their original common variable. A shared
Newton slope supplies a possible root valuation. Its two face
polynomials must also be compatible at the precision supplied by
both adjacent coefficient arrays. This comparison can remove slopes
that pass the Newton width and residual gcd tests.

## Common-root widths at higher precision

Fix a prime \(p\) and a candidate root valuation \(s=a/e\) in lowest
terms, with \(e>0\). Let \(\pi^e=p\), let \(\mathcal O\) be the
valuation ring of \(\mathbb Q_p(\pi)\), and set \(w=e v_p\). This
ring has uniformizer \(\pi\) by the Eisenstein criterion, including
the case \(e=1\). Write \(P_j(z)=\sum_k c_{j,k}z^k\), put
\[
        m_j=\min_k \bigl(w(c_{j,k})+ak\bigr),\qquad
        F_j(X)=\pi^{-m_j}P_j(\pi^a X),
\]
and fix a positive integer \(L\). In \(F_j\) discard every term
whose coefficient has \(w\)-valuation at least \(L\) and divide by
the largest remaining monomial power of \(X\), meaning \(X\) to
the minimum retained index. Call the resulting polynomial \(A_j\).

**Theorem.** Let \(t_1\le\cdots\le t_D\) be the valuations of the Smith
entries of the Sylvester matrix of \(A_n,A_N\) over
\(\mathcal O\), with \(w(0)=+\infty\). Every rational common factor
has at most \(\#\{i:t_i\ge L\}\) roots of valuation \(s\), counted
with multiplicity. This applies at each \(L\). A constant \(A_j\)
can be interpreted as giving width zero, since its coefficient
is a unit.

For the proof, if a rational common factor has \(r_s\) roots of
valuation \(s\), the monic polynomial \(H_s(X)\) on their images
under \(z\mapsto z/\pi^a\) is invariant under the Galois group over
\(\mathbb Q_p(\pi)\). Its coefficients belong to \(\mathcal O\)
because every root is a unit. Its constant coefficient is a unit.
Monic division makes the divisibilities \(H_s\mid F_j\) valid over
\(\mathcal O\). Reduction modulo \(\pi^L\) and invertibility of \(X\)
in the quotient by \(H_s\) then give
\[
                H_s\mid A_n,\ A_N \pmod{\pi^L}.
\]
In this ring the degree bound \(r_s\le\min(\deg A_n,\deg A_N)\)
follows from monic division and the unit coefficient retained in
each \(A_j\). Put \(R=\mathcal O/(\pi^L)\) and
\(D=\deg A_n+\deg A_N\). The cokernel of the Sylvester map onto
\(R[X]_{<D}\) surjects by reduction on \(R[X]/(H_s)\cong R^{r_s}\).
The cokernel is \(\bigoplus_i R/(\pi^{t_i})\), where exponents
\(\ge L\) give copies of \(R\). Multiplying the surjection by
\(\pi^{L-1}\) gives the stated bound over the residue field.

The ordinary Newton polygon product rule assigns root valuations
with multiplicity to the negatives of the unit slopes.
[Filaseta and Rocks, Theorem 2.1](https://doi.org/10.4064/cm6761-3-2016)
state the polygon form of this rule. Combining the bounds over
shared root valuations and using positivity of the even factorial
moment gives a coprimality certificate whenever the total remaining
width is at most one. The case \(L=1\) is residual gcd matching.
At larger \(L\) the argument allows any prime and all residual
multiplicities.

## Binomial faces and denominator bounds

Suppose the two minima for root valuation \(s\) are attained
exactly at \(l_j,u_j\) in each order. Write
\[
 \begin{aligned}
   b_j&=u_j-l_j,&
   R_j&=-c_{j,l_j}/c_{j,u_j},\\
   \delta_j&=\min_{k\notin\{l_j,u_j\}}\!
        \left(v_p(c_{j,k})+ks-v_p(c_{j,l_j})-l_js\right),&
   g&=\gcd(b_n,b_N).
 \end{aligned}
\]
The numbers \(v_p(R_j)=b_j s\) are integers and \(\delta_j>0\).
A common root of valuation \(s\) forces
\[
       v_p\!\left(\frac{R_n^{\,b_N/g}}{R_N^{\,b_n/g}}-1\right)
                       \ \ge\ \min(\delta_n,\delta_N).           \tag{1}
\]
If this test passes, its common width is at most \(g\).

Indeed take \(L=e\min(\delta_n,\delta_N)\).
The two retained binomials, made monic by unit division, are
\(X^{b_n}-B_n,X^{b_N}-B_N\), with \(B_j=p^{-b_js}R_j\) units.
Set \(u=b_n/g,\ v=b_N/g\), choose integers \(a',b'\) with
\(a'u+b'v=1\), and put \(E=B_n^{a'}B_N^{b'}\).
In the Laurent polynomial ring the two binomials generate
the same ideal as \(X^g-E\) and \(B_n^v-B_N^u\).
The forward implication follows from taking the corresponding
powers of \(X^{b_j}=B_j\). For the reverse implication,
\(B_n^v=B_N^u\) gives \(E^u=B_n\) and \(E^v=B_N\).
The quotient is unchanged when \(X\) is inverted because all
of the constant terms in the binomials are units.
The Sylvester cokernel is consequently a free module of rank
\(g\) over \(\mathcal O/(B_n^v-B_N^u)\). To see the cokernel
identification in this case, divide an arbitrary polynomial
by the first monic binomial. An element of their ideal of
degree less than \(b_n+b_N\) can be written with multiplier
degrees less than \(b_N,b_n\), respectively, by reducing the
second multiplier modulo the first binomial and comparing the
remaining degrees. Thus the positive Smith valuations are
read off from \(g\) copies of \(w(B_n^v-B_N^u)\), with the
other entries at valuation zero. The common
normalization powers cancel in \(B_n^v/B_N^u\).
This proves (1) and the width bound. Equivalently, (1) follows
by comparing the two approximate equations for a shared root.

At a prime, delete failed slopes and reduce each surviving width
by the strongest available bound. The
[exact denominator formulas](MULTIVARIATE_ORIENTED_BLOCK_DENOMINATOR_GCD.md#means-of-consecutive-blocks-of-valuations)
for \(L_H,D_H\) may then be maximized over this smaller multiset
instead of the full shared Newton multiset. This applies with no bound on the number of row types, the gaps,
or the primes. For instance, if only one root valuation survives at a
prime, that prime occurs in neither denominator. A prime with too
little remaining width cannot support a common factor of the
proposed degree at any valuation.

## A high-prime terminal-strip comparison

Fix a prime \(p>nN\max_i d_i\), where \(d_i>0\) is the absolute
gap of row \(i\) and \(a_i\) its smaller exponent. Suppose that
every hit in either stripped ramp \(ja_i+1,\ldots,ja_i+jd_i\)
has offset \(j t\) for an integer \(t\) with
\(Nt>(N-1)d_i\). Such a hit is synchronous in both orders and
lies in the last rising block of its row. It has valuation
\(v_p(a_i+t)\) in both orders. Define, for each orientation,
\[
 V^\pm=\prod_{i:\,\pm}\ \prod_{(N-1)d_i<Nt\le Nd_i}(a_i+t),
 \qquad E_\pm=v_p(V^\pm).
\]
After stripping the row factorials the ordinates of each parent
are \((E_-,0,\ldots,0,E_+)\). Hence a negative root valuation can
only be \(-E_+\), of width one, when \(E_+>0\). Set
\[
 W_j=\frac{c_{j,j}}{c_{j,j-1}V^+}.
\]
These are \(p\)-adic units. The off-face depth at root valuation
\(-E_+\) is \(E_+\) when \(j\ge3\), and it is
\(E_++E_-\) when \(j=2\). Equation (1) gives the necessary
condition \(v_p(W_n/W_N-1)\ge E_+\). The positive root valuation
has the reversed condition obtained by interchanging the
columns. Each failed condition deletes its valuation group
from the denominator budget even if the parent block gcds
contain the entire shared strip. The calculation applies to
arbitrarily many moving row types.

There is a direct global gcd bound. Let \(\mathcal P\) be any
set of primes satisfying these terminal-strip hypotheses for the
same profile. Write \(\Theta^+\) for the reduced numerator of
\((c_{n,n}/c_{n,n-1})/(c_{N,N}/c_{N,N-1})-1\) and
\(\Theta^-\) for the numerator after reversing both columns.
If \(H\) has degree \(d_H\), then
\[
        (D_H)_{\mathcal P}\ \mid\
          \bigl(\gcd(V^+,\Theta^+)\gcd(V^-,\Theta^-)\bigr)^{d_H-1}.
\]
Here \((D_H)_{\mathcal P}\) is the part supported on \(\mathcal P\)
and \(\gcd(V,0)=V\). At a prime in \(\mathcal P\) the factors
\(V^\pm\) cancel in the cross ratio, and its denominator is a
unit. A surviving nonzero valuation of each sign therefore has
magnitude at most the corresponding gcd valuation. Each sign
occurs at most once. The exact formula
\(v_p(D_H)=\sum_{i<j}|s_i-s_j|\) for root valuations \(s_i\)
then gives at most \(d_H-1\) times the sum of the two magnitudes.

## An exact precision example

At \(n=3,N=4,p=17\) use seventeen increasing unit-gap rows
with smaller exponents
\[
       3371,\quad 22108,\quad 730808,\quad 17^7-1,\quad 15,
       \quad\underbrace{0,\ldots,0}_{12},
\]
and seventeen decreasing unit-gap rows of smaller exponent zero.
The adjacent coefficient valuation differences in the lower order
are \((0,2,7)\). In the upper order they are \((0,3,4,7)\).
The shared slopes are exactly \(0,7\), both of width one.
At slope zero the unit first-block ratios agree modulo \(17\).
At slope \(7\) the ratio of the two last-block ratios is
\((3/4)^{16}\), which is one modulo \(17\) and differs from one
at valuation exactly one. The off-face depths of this second
pair are \(5\) and \(3\). Thus both residual matches pass, while
(1) removes slope \(7\) and leaves width at most one.

The numbers in this check require only short products. The
three nonsynchronous active ramps have
\(v_{17}(3\cdot3371+2)=2\),
\(v_{17}(4\cdot22108+2)=3\), and
\(v_{17}(4\cdot730808+3)=4\).
Their residues \(5,8,12\) are distinct from the other order's
active ramp residues. The synchronous base is \(-1\) and
its last factors have valuation seven. The neutral residues
\(15,0\) meet neither ramp. In the first-block ratio comparison
\((c_{3,1}/c_{3,0})/(c_{4,1}/c_{4,0})\), the factor from the
binomial coefficients and negative rows is \((3/4)^{-16}\).
The product of the five nonzero-base factors
\((3a+1)/(4a+1)\) is one modulo \(17\).
Finally \((3^{16}-4^{16})/17\equiv1\pmod{17}\).

## All-degree congruence classes with arbitrarily many types

Fix any finite ordinary core and any integer signed gap \(h\) for
the additional rows. There are finite residue prescriptions for
a finite labelled list of increasing and decreasing unit-gap
rows and a positive period \(B\), with the following properties.
All labelled heights may vary independently in their prescribed
classes. Any number of further increasing and decreasing unit-gap
rows may be added independently at smaller exponents divisible by
\(B\), provided the two extra counts are equal and are multiples
of an explicitly specified fixed integer. Every profile in this
family is coprime in every degree. For any \(\tau>0\) and
\(|\lambda|<\tau\), profiles in the family can have maximum tail
relative gap tending to zero, with curvature and signed bias
tending to \(\tau,\lambda\), and with the number of distinct tail
types tending to infinity.

To construct the classes, let \(J\) be the
[absolute angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md).
It suffices to handle \(3\le N<J\). The order-one root is real
and does not meet the strictly positive even moment. Choose a
number \(r_+\) of labelled increasing rows with \(r_+\ge 2J\)
and \(r_-=r_+-h\ge0\). At each pair \(n,N\), let \(K_j\) be the
product of the core's last adjacent factorial ratios, with no
binomial multiplier. Choose \(r_+\) so that
\[
                 Z_{n,N}=(n/N)^{r_+-1}K_n/K_N\ne1
                         \quad(3\le N<J).
\]
Each pair excludes at most one integer \(r_+\), because \(n/N<1\).

For each pair choose a separate prime \(p>nN\) also greater than
all ramp factors of the fixed core at that pair. The core then
makes no contribution to the stripped ordinates. Put
\(\nu=v_p(Z_{n,N}-1)\ge0\). In every prescribed class for this
prime use a modulus \(p^{E+1}\) with \(E>2N+1+\nu\).
At each order \(j\), give one increasing labelled row to each
offset \(2\le k<j\). Prescribe its base at this prime by
\[
        v_p(ja+k)=
         \begin{cases}2k,&j=n,\\2k+1,&j=N.\end{cases}
\]
There is no other hit from that row, as the residue classes
\(-k/j\) are distinct except for the shared terminal
\(-n/n=-N/N\). To see this across the orders, clear denominators
and use \(p>nN\) together with \(\gcd(n,N)=1\).
Give one increasing labelled row the shared terminal
\(v_p(a+1)=E\), and make all other labelled rows have bases
zero modulo \(p^{E+1}\). The distinct assigned labels fit
inside \(r_+\). All constraints are congruences at the stated
modulus, and all selected valuations are stable there.

The lower Newton unit slopes are the successive differences
\(0,2k\ (2\le k<n),E\), in order, and the upper ones are
\(0,2k+1\ (2\le k<N),E\). Thus only \(0,E\) are shared and
both have width one. Let \(V=\prod_+(a_i+1)\) for the full
unit-gap tail, including the extra rows. At the shared terminal
slope, the ratio of the two last adjacent coefficient ratios
is
\[
  \frac{c_{n,n}/c_{n,n-1}}{c_{N,N}/c_{N,N-1}}
      = Z_{n,N}(n/N)^t\prod_-\frac{Na_i+1}{na_i+1},
\]
where \(t\) is the extra increasing count. The common factor
\(V\) has canceled. If all extra bases are zero modulo
\(p^{E+1}\), and \(t\) is a multiple of \(\varphi(p^{E+1})\),
this ratio differs from one at valuation \(\nu\). The off-face
depth in either parent is \(E\) minus the largest preceding
slope. It exceeds \(\nu\), so (1) deletes the terminal slope.
Taking \(B\) to be the product of the prime power moduli and
using the Chinese remainder theorem gives the independent
height prescriptions. The extra count condition is divisibility
by the least common multiple of their Euler totients. The
cutoff completes the all-degree statement.

Finally let the common extra count \(t\) tend to infinity through
these multiples and move the labelled bases to infinity in their
classes. Put \(T_\pm=(\tau\pm\lambda)/2>0\), and choose
\(A_\pm=(\exp(BT_\pm)-1)^{-1}\).
In orientation \(\pm\), take the \(t\) distinct extra bases
\(B(\lfloor A_\pm t\rfloor+i)\) for \(1\le i\le t\).
The sum of their unit-gap curvatures tends by a Riemann sum to
\[
       B^{-1}\log\bigl((A_\pm+1)/A_\pm\bigr)=T_\pm .
\]
Their signed bias is the corresponding signed curvature, and
the difference \(h\) of total degrees is unchanged by the
balanced extras. Thus these exactly coprime classes can
approximate each interior diffuse parameter, including any
positive exceptional reference parameter of the chosen core.
The congruence construction supplies sufficient conditions.
The width and denominator tests alone give no universal supply
of favorable primes for an arbitrary diffuse sequence.
