# Gamma quotients and coprimality on the binary factorial boundary

For every integer \(u\ge0\) and every \(m\ge1\), the factorial moments
of \(Y^u(aX+bY)\) in orders \(m,m+1\) have no common projective
coefficient zero over any characteristic-zero field. Interchanging
\(X,Y\) proves the two coordinate axes of
`GAP-FACTORIAL-BINARY-LINKED-BETA-COPRIMALITY`. The general two-positive,
unequal-exponent case remains open.

## 1. Gamma-quotient determinant signs

In the [external proof, Sections 2–7](https://github.com/1337isnot1337/solved-conjectures/blob/758066055ef49c2c66cb1087918ae62ad75e7ac1/results/kvz-gamma-quotient/PROOF.md),
put \(\delta=c-d>0\). Translate each finite increasing grid so that the
denominator is \(x_i+y_j\) and all \(x_i,y_j>0\). Such a split of the
translation \(d\) exists exactly when \(d+x_1+y_1>0\).
With \(G(a,t)=t^{a-1}e^{-t}/\Gamma(a)\), the change of variables
\(r=u+v,\ z=u/(u+v)\), of Jacobian \(r\), gives the decomposition
\[
 \frac{\Gamma(x+y+\delta)}{\Gamma(x+y)}
 =\int_0^\infty\int_0^\infty G(x,u)(u+v)^\delta G(y,v)\,du\,dv .
\]
A direct verification of the signs uses the following facts.

- For increasing \(a_i\), the initial Wronskians of \(e^{a_i t}\) are a
  positive exponential times the positive Vandermonde in the \(a_i\).
  Consequently the gamma-density determinant on ordered grids is
  strictly positive after extracting positive row and column factors.
- For fixed increasing \(v_j\), the order \(h\) initial Wronskian of
  \((u+v_j)^\delta\) has sign
  \(\operatorname{sgn}\prod_{j=0}^{h-1}(-\delta)_j\).
  Differentiation extracts \(\delta^{\underline j}\) from each row,
  \((u+v_j)^\delta\) from each column, and leaves the Vandermonde
  in the decreasing nodes \(1/(u+v_j)\).
- A system with all initial Wronskians nonzero has nonzero ordered
  evaluation determinants of the same signs. Dividing by the first
  function and differentiating reduces the system size, because the
  new Wronskians are the old Wronskians divided by the corresponding
  powers of the first function. Rolle's theorem gives the zero bound
  by induction. Coalescing ordered nodes determines the sign.
- Expand the three determinants in the two continuous Cauchy–Binet
  formulas. All terms are absolutely integrable by
  \((u+v)^\delta\le(1+u)^\delta(1+v)^\delta\). On the strictly ordered
  integration regions the determinant factors have fixed signs, and
  the product is strictly of that sign except in the finite-rank case.

Here \((z)_j=z(z+1)\cdots(z+j-1)\) with \((z)_0=1\). Thus for noninteger
\(\delta>0\) the order-\(h\) sign is precisely the sign of
\(\prod_{j<h}(-\delta)_j\). For an integer \(L=\delta\ge1\) the argument
works through order \(L+1\). Higher minors vanish because the quotient
is a polynomial of degree \(L\) in \(x+y\). In particular, order three
has sign \(\operatorname{sgn}(1-\delta)\), with zero at \(\delta=1\).
For \(c<d\) the Euler beta integral has positive measure and both
alternants \(t^{x_i}\), \(t^{y_i}\) have positive determinants by the
same exponential argument. This gives strict total positivity in that
regime whenever \(c+x_1+y_1>0\).

The external [literal Lean theorem](https://github.com/1337isnot1337/solved-conjectures/blob/758066055ef49c2c66cb1087918ae62ad75e7ac1/lean/Conjectures/KVZGamma.lean)
and [verification record](https://github.com/1337isnot1337/solved-conjectures/blob/758066055ef49c2c66cb1087918ae62ad75e7ac1/docs/verification/lean-kvz-20261121/receipt.json)
are pinned to the same revision. The written proof above rechecks the
strictness, rank and sign independently. The external verification
record reports a strict replay of the frozen closure. This audit uses
the written argument and does not claim a new Lean build.

### A product formula on unit-spaced columns

For increasing \(x_i\), put \(a_i=c+x_i+y_0,\ b_i=d+x_i+y_0>0\), also
assuming \(a_i>0\). Then
\[
 \det_{1\le i\le h,\ 0\le j<h}
       \frac{\Gamma(a_i+j)}{\Gamma(b_i+j)}
 =
 \prod_i\frac{\Gamma(a_i)}{\Gamma(b_i+h-1)}
 \prod_{i<j}(x_j-x_i)
 \prod_{j=0}^{h-1}(d-c)_j.                                         \tag{1}
\]
Indeed, scale each row to obtain evaluation of the polynomials
\((c+x+y_0)_j(d+x+y_0+j)_{h-1-j}\). Their determinant is a scalar
times the Vandermonde since their degrees are at most \(h-1\).
Evaluate at \(x_i=-c-y_0-(i-1)\) to make the polynomial matrix lower
triangular, with diagonal \((-1)^{i-1}(i-1)!(d-c)_{h-i}\). The
signs and factorials cancel those in the Vandermonde at these nodes.
This computation uses only polynomial evaluation after clearing the
Gamma factors. It gives the actual determinant value, including its
size and zeros, on unit-spaced columns.

## 2. An obstruction on the actual linked coefficient matrix

The linked orbit in the
[canonical-density audit](canonical-density/moment_pair_obstructions.md)
has row \(m\) given by the moments of
\(\operatorname{Beta}(mr+1,mu+1)\). After removing the
signed binomial factors, the raw moment matrix is
\[
 M_{m,k}=\frac{(mr+1)_k}{(m(r+u)+2)_k},\qquad
 \frac{M_{m,k+1}}{M_{m,k}}=\frac{mr+1+k}{m(r+u)+2+k}.
\]
For rows \(m_1<m_2\), subtract the consecutive-column ratios. The
numerator is \((m_2-m_1)(r-u(k+1))\) after putting over a common
denominator, so the column-\((k,k+1)\) minor has that sign.
Already \(r=2,u=1\) and \(m_1=3,m_2=4\) give a positive minor for
\(k=0\) and a negative one for \(k=2\), both in the actual range
\(k+1\le m_1\). Positive diagonal rescaling preserves both signs,
and uniform reversal of either index changes both together.

A composition of kernels having a fixed sign per minor order against
a positive measure has a fixed minor sign by Cauchy–Binet, when the
integrals are defined. This excludes such a composition for the
whole linked raw moment matrix in its native orders. It leaves
particular subquotients and the nonlinear Sylvester determinant
as possible receivers. The [shifted-resultant audit](SFC_SHIFTED_RESULTANT_POSITIVITY_AUDIT.md)
poses its positivity question for that different determinant.


The same polynomial determinant argument as in (1) yields a
stronger test on the linked matrix itself. For \(r>0,\ u\ge0\) and
increasing rows \(m_1,\ldots,m_h\), the consecutive-column minor is
\[
\begin{aligned}
 \det_{i,\ 0\le j<h} M_{m_i,k+j}
 ={}&
 \prod_i\frac{M_{m_i,k}}{((r+u)m_i+k+2)_{h-1}}\,
 \prod_{i<j}(m_j-m_i)\,r^{\binom h2}\\
 &{}\times\prod_{i=0}^{h-1}
        \left(1-\frac ur(k+i+1)\right)_{h-1-i}.
\end{aligned}                                                       \tag{2a}
\]
After clearing the indicated row factors the column polynomials
are \((r m+k+1)_j((r+u)m+k+2+j)_{h-1-j}\). Their determinant is a
constant times the Vandermonde. Evaluate at
\(m=-(k+1+i)/r\) to obtain the triangular determinant. Its
\(i\)-th diagonal entry is \((-1)^i i!\) times
\((1-\frac ur(k+i+1))_{h-1-i}\). Dividing by the evaluation
Vandermonde leaves exactly the claimed scalar.

In particular, if \(u/r=\eta>1\) is nonintegral, use
\(h=\lfloor\eta\rfloor+2\). At \(k=0\) exactly one of the
\(\binom h2\) linear factors of the final product in (2a) is positive.
It is the last factor for \(i=0\). Every other factor is strictly
negative, since \(h-1<\eta+1\) and \(i\ge1\) starts at
\(1-2\eta\) or lower. For sufficiently large \(k\) every factor is
negative, so the two minors have opposite signs. Take the same
rows in both minors with \(m_1\ge k+h-1\) to stay inside the actual
coefficient ranges. Reflection swaps \(r,u\). When the ratio of the
larger exponent to the smaller is nonintegral, the reversed
orientation with \(r>u>0\) already fails by the order-two calculation.
Thus neither orientation has a sign-regular raw matrix on these
profiles. For example, at \(r=2,u=3\) and rows \(3,4,5\), the
contiguous minors beginning at columns \(0,1\) are respectively
\(1/14632002\) and \(-7/139605768\).


The order \(h=\lfloor\eta\rfloor+2\) is the first failure of
sign regularity in this orientation. For \(s<h\) even the
largest linear factor in the last product of (2a), obtained at
\(i=k=0\), is at most \(s-1-\eta<0\). Hence every contiguous
minor through order \(h-1\) has the strict sign
\((-1)^{\binom s2}\). The discrete argument below, with all
thresholds zero through that order, gives the same strict signs
on arbitrary columns. In the reflected orientation \(r>u\),
entries are positive and the first failure is at order two.



There is an exact sign and zero pattern in the orientation \(u=Lr\),
where \(L\ge1\) is an integer. Formula (2a) gives the order-\(h\)
consecutive-column minor beginning at column \(k\) sign
\((-1)^{\binom h2}\) if \(h\le L(k+1)\) and zero otherwise.
When \(h\le L(k+1)\), every linear factor in the last product is
strictly negative. If \(h>L(k+1)\), the factor for \(i=0\)
contains zero. For **arbitrary** increasing columns
\(k_1<\cdots<k_h\), the minor has the same weak sign and is
nonzero exactly when
\[
                    k_j\ge j-2+\lceil j/L\rceil
                    \qquad(1\le j\le h).                         \tag{2b}
\]
The row indices only need to be strictly increasing. For the
physical array, take \(k_h\le m_1\).

Here is the general discrete argument underlying (2b). Suppose
the order-\(s\) determinant using **any** \(s\) ordered rows and
\(s\) consecutive columns starting at \(k\) is zero exactly for
\(k<K_s\) and has sign \(\epsilon_s\in\{1,-1\}\) otherwise.
Assume \(K_s\) is nondecreasing. Entries in columns \(k<K_1\)
are all zero, and the criterion is immediate at order one.
Restrict to \(k\ge K_1\), fix \(h\) rows, and divide each column by
the absolute value of its first entry. The new first row is
constant \(\epsilon_1\). Form the one-step horizontal differences
\(D_i(k)\) of the remaining normalized rows. For any \(s\)
of these rows, their consecutive \(s\)-column determinant at
\(k\), multiplied by \(\epsilon_1\), equals the normalized
former consecutive determinant starting at \(k\) on those rows
plus the first row. This is column subtraction from right to
left followed by expansion along the constant first row.
In particular, every entry of \(D(k)\) is zero for \(k<K_2\).
On the surviving column interval, the difference matrix has
thresholds \(K_{s+1}\) and signs \(\epsilon_1\epsilon_{s+1}\).

For skipped columns \(k_j\), use the same subtraction. Each
remaining column is now a sum of the \(D(\ell)\) over
\(k_j\le\ell<k_{j+1}\). The intervals are disjoint and ordered,
so determinant expansion chooses \(h-1\) increasing indices.
The induction hypothesis gives every nonzero summand sign
\(\epsilon_1\epsilon_h\), and its criterion is
\(\ell_j\ge K_{j+1}+j-1\). Such a choice exists exactly when the
largest index in each interval, \(k_{j+1}-1\), meets its bound.
Expansion along the first row gives the original sign
\(\epsilon_h\). Thus the original minor is nonzero exactly when
\(k_j\ge K_j+j-1\) for all \(j\). Taking \(K_s=\lceil s/L\rceil-1\)
proves (2b).

Consequently, for \(r,u>0\) the raw matrix in increasing row
and column order is weakly sign regular exactly when \(u/r\)
is a positive integer. The orientation with \(r>u\) fails at
order two even when \(r/u\) is an integer. For comparison,
the \(u=0<r\) axis is strictly totally positive by (2a).
On the \(r=0<u\) axis, \(M_{m,k}=k!/(mu+2)_k\). Remove the
positive column factorials and clear row denominators to
compute each consecutive-column determinant from polynomials
\((mu+2+k+j)_{h-1-j}\) of strictly decreasing degrees.
The determinant has sign \((-1)^{\binom h2}\), and the above
argument extends strictness to arbitrary ordered columns.
At \(r=u=0\) all rows coincide, giving rank one.
These sign and rank results concern the raw moment matrix.
The coefficient binomial factors and the two consecutive
polynomial degrees still enter the coprimality problem.



## 3. All orders when one monomial exponent is zero

Normalize as in [the shifted-window identity](SFC_INITIAL_RESULTANT_AND_SHIFTED_WINDOWS.md).
For \(a=1\) and \(s=b/a\), the actual moment of \(Y^u(X+sY)\) at
order \(m\) is \(m!(mu)!\) times
\[
 P_{m,mu}(s),\qquad P_{k,q}(s)=\sum_{j=0}^k\binom{q+j}{j}s^j.       \tag{2}
\]
This follows by expanding \([z^m](1-z)^{-1}(1-sz)^{-mu-1}\).
The chart \(a=0\) has a nonzero pure monomial moment for \(b\ne0\).
It suffices to prove coprimality over \(\mathbb C\) since the moment
polynomials are rational. When \(u=0\), the difference of the two
polynomials is \(s^{m+1}\), and both have constant term one.

For real \(s\) the even member of the consecutive pair is strictly
positive. The factorial functional is integration against
\(e^{-X-Y}\,dX\,dY\) on the positive quadrant, where the even power
of \(Y^u(X+sY)\) is positive off its zero set. We only need to
exclude common nonreal roots when \(u\ge1\).

For \(q\ge0\) put \(A_{k,q}=(1-s)^{q+1}P_{k,q}\). One has
\[
 \begin{aligned}
 A_{k+1,q}-A_{k,q}
   &=\binom{k+q+1}{k+1}s^{k+1}(1-s)^{q+1},\\
 A_{k,q+1}-A_{k,q}
   &=-\binom{k+q+1}{k}s^{k+1}(1-s)^{q+1}.
 \end{aligned}                                                     \tag{3}
\]
The first adds the terminal coefficient and the second follows by
multiplying \(P_{k,q+1}\) by \(1-s\) and using Pascal's identity.
First advance \(k\) from \(m\) to \(m+1\), then advance \(q\) by \(u\).
With \(q=mu,\ v=1-s\), the result is
\[
 \frac{A_{m+1,q+u}-A_{m,q}}
      {\binom{m+q+1}{m+1}s^{m+1}(1-s)^{q+1}}
  =R(v)=1-(1-v)\sum_{\ell=1}^u a_\ell v^{\ell-1}
  =a_uv^u-\sum_{i=0}^{u-1}(a_{i+1}-a_i)v^i,
 \qquad a_\ell=\frac{(m+q+2)_\ell}{(q+1)_\ell}.                     \tag{4}
\]
All denominators are nonzero at a common nonreal root. The quotient
coefficients are governed by the same unit-shift relation as the
Gamma factors in (1), namely
\(a_{i+1}-a_i=(m+1)a_i/(q+i+1)\).

Set \(\rho=u/(u+1),\ L=a_u\rho^u,\ c_i=(a_{i+1}-a_i)\rho^i\).
Then \(R(\rho w)=Lw^u-\sum_{i=0}^{u-1}c_iw^i\), and
\[
 \frac{c_i}{c_{i-1}}
  =1-\frac{i+1}{(u+1)(q+i+1)}\quad(1\le i<u),\qquad
 (1-1/m)c_0<c_i\le c_0\le(1+1/m)/u.                              \tag{5}
\]
For the lower bound each factor in the ratio product is at least
\(1-1/(m(u+1))\), giving
\(c_i/c_0\ge1-i/(m(u+1))>1-1/m\).
Also \(a_i\rho^i\) is nondecreasing for \(0\le i\le u\), since the
one-step ratio is \(1+(u-i-1)/((u+1)(q+i+1))\) for \(i<u\).
It follows that \(L\ge1\) and
\(R(\rho)=1-u^{-1}\sum_{\ell=1}^u a_\ell\rho^\ell\le0\).
There is exactly one positive root \(w_0\) of \(R(\rho w)\).
The equation \(L=\sum_i c_i w^{-(u-i)}\) has a strictly decreasing
right side from infinity to zero, and gives \(w_0\ge1\).

For \(m\ge4\), take \(x=3/(mu)\). For \(1\le j\le u\),
\((1+x)^{-j}\le(1+jx)^{-1}\le1-jx/(1+ux)\). Sum this bound,
using (5) and \(\sum_{j=1}^u j\ge u^2/2\), to obtain
\[
 \sum_{i=0}^{u-1} c_i(1+x)^{-(u-i)}
 \le (1+1/m)\left(1-\frac{3}{2(m+3)}\right)\le1\le L.
\]
Thus \(1\le w_0\le1+3/(mu)\). Divide \(R(\rho w)\) by \(w-w_0\).
The quotient coefficients \(b_i\) satisfy
\(b_0=c_0/w_0,\ b_{i+1}=(b_i+c_{i+1})/w_0\), so \(b_i>0\) and
\(b_i\le(i+1)c_0\) by induction. For \(i\le u-2\),
\[
 (w_0-1)b_i\le \frac{3(i+1)}{mu}c_0
   <\frac3m c_0\le (1-1/m)c_0<c_{i+1}.
\]
Consequently \(b_{i+1}>b_i\). The quotient is constant for \(u=1\).
Otherwise all its roots satisfy \(|w|<1\) by the
Eneström–Kakeya bound. In the form used here, multiplying a
polynomial with nondecreasing positive coefficients by \(1-w\)
shows that the leading term dominates all the others for \(|w|>1\).
Scaling the argument by the largest ratio of consecutive
coefficients gives a strict disk when the coefficients increase strictly.

For \(P_{m,mu}\) the ratio of a coefficient to the next is
\((j+1)/(q+j+1)\le 1/(u+1)\) for \(0\le j<m\). The same bound
puts every lower root in \(|s|\le1/(u+1)\), hence
\(|v|=|1-s|\ge \rho\). By (4), a nonreal common root would be a
nonreal zero of \(R(v)\) with \(|v|\ge\rho\). All zeros except
the positive \(v=\rho w_0\) satisfy \(|v|<\rho\), a contradiction.

The exact Sylvester resultants for \(m=1,2,3\), respectively, are
\[
 \begin{split}
 &(u+1)^2,\\
 &\frac{(u+1)^3}{4}(17u^3+37u^2+22u+4),\\
 &\frac{(u+1)^4}{432}(44519u^8+235223u^7+518466u^6+622063u^5\\
 &\hspace{38mm}+444511u^4+193890u^3+50496u^2+7200u+432).
 \end{split}                                                       \tag{6}
\]
They are positive for \(u\ge0\). Each is directly checkable from (2).
The Sylvester products use coefficient indices summing to
\(m(m+1)\) and the coefficient of \(s^j\) has degree at most \(j\)
in \(u\). Thus the \(m(m+1)+1\) exact evaluations at successive
integers starting at zero give a full interpolation certificate
for each formula. The [shifted-resultant audit](SFC_SHIFTED_RESULTANT_POSITIVITY_AUDIT.md)
also certifies the full bivariate identities at these orders.
Equations (4)–(6) prove the boundary theorem for every starting order.

## 4. A shifted Gamma shape on the same boundary

Let \(X,Y\) be independent Gamma variables of shapes \(1,\beta\),
respectively, at rate one, with any real \(\beta>0\). For every
integer \(u\ge0\) and every \(m\ge1\), the two moments of
\(Y^u(aX+bY)\) at orders \(m,m+1\) have no common projective
zero over \(\mathbb C\). Interchanging the variables gives the
reflected statement with the shifted shape on the variable
carrying the common monomial power. This includes the ordinary
exponent-axis claim at \(\beta=1\). For positive real \(\beta\)
the polynomial coprimality also holds over any extension of the
coefficient field \(\mathbb Q(\beta)\), by the Bézout identity.

Indeed, on the chart \(a=1\), the order-\(m\) moment divided by
\(m!(\beta)_{mu}\) is \(P_{m,q}(s)\) of (2) with \(q=mu+\beta-1\).
The Gamma moments \(\mathbb E[X^j]=j!\) and
\(\mathbb E[Y^j]=(\beta)_j\) give this directly from the
binomial expansion. Pure monomials and real candidate roots
are handled as in Section 3, since the product density is
positive on the positive quadrant. For \(u=0\), the upper
polynomial adds a single nonzero monomial to the lower one
with the same \(q=\beta-1\).

Suppose \(u\ge1\), so \(q>0\). The lower root bound is
\[
 |s|\le \frac{m}{m+q}<1,\qquad
 |1-s|\ge \rho=\frac q{q+m}.
\]
The successive coefficient ratios in (2) establish the first
inequality. On \(|s|<1\) define \((1-s)^{q+1}\) by the analytic
logarithm real at \(s=0\). The same Pascal calculation proves
(3) for real \(q\), and the successive \(u\) unit shifts preserve
the physical parameters in (4). In particular its expression
\(R(v)\) is polynomial in \(v=1-s\) even for noninteger \(q\).
For \(u=1\) it is linear with real coefficients and nonzero
leading coefficient.

For \(u\ge2,\ m\ge4\), use \(L,c_i\) from Section 3 with this
new value of \(\rho\). The ratio formulas now read
\[
 \frac{a_{i+1}\rho^{i+1}}{a_i\rho^i}
 =1+\frac{q-m(i+1)}{(q+m)(q+i+1)},\qquad
 \frac{c_i}{c_{i-1}}
 =1-\frac{m(i+1)}{(q+m)(q+i+1)}.
\]
In the second formula \(1\le i<u\), so \(q+i+1>mu\) and each
loss is below \(1/(q+m)<1/(m(u+1)-1)\). Consequently
\(c_i/c_0\ge1-i/(m(u+1)-1)>1-1/m\), while
\(c_0=(m+1)/(q+1)\le(1+1/m)/u\). For the first formula,
all factors preceding the last are at least one. The last is
also at least one if \(q\ge mu\). Otherwise the penultimate
factor is at least \(1+1/D_{u-2}\) and the final factor is
greater than \(1-1/D_{u-1}\), where
\(D_i=(q+m)(q+i+1)\). Since \(D_{u-1}-D_{u-2}=q+m>1\),
their product is greater than one. Thus \(L\ge1\) in every case.

The positive root \(w_0\) of \(R(\rho w)\) is still unique.
If \(w_0\le1\), every root \(w\) satisfies \(|w|\le w_0\):
the root equation and the triangle inequality give
\(L\le\sum_i c_i/|w|^{u-i}\), a strictly decreasing function of
\(|w|>0\). At equality the terms of degrees \(0\) and \(1\)
must have the same argument, which forces \(w>0\).
Hence each nonreal root has \(|w|<1\) in this case.
If \(w_0>1\), the bounds on the \(c_i\) and \(L\) give
\(w_0\le1+3/(mu)\) by the same estimate after (5).
Dividing by \(w-w_0\) gives strictly increasing positive
quotient coefficients by the same recurrence and inequalities
in Section 3. Thus every nonreal zero of \(R(v)\) satisfies
\(|v|<\rho\). This excludes a common nonreal zero of the
physical moments for \(m\ge4\).

For the remaining starts, the resultant of the two shifted
polynomials has nonnegative ordinary coefficients in \(u,\beta\)
and strictly positive coefficient of \(\beta^{m(m+1)}\). The
complete finite certificates are as follows.

| Starting order \(m\) | Total degree bound | Positive coefficients | Coefficient of \(\beta^{m(m+1)}\) |
|---:|---:|---:|---:|
| 1 | 2 | 5 | \(1/2\) |
| 2 | 6 | 25 | \(1/36\) |
| 3 | 12 | 85 | \(1/13824\) |

To verify the certificate, the coefficient of \(s^j\) is
\((mu+\beta)_j/j!\), of total degree at most \(j\).
Every Sylvester product uses indices summing to \(m(m+1)\),
as in (6), so evaluating the fixed-size Sylvester determinant
at integer pairs \((i,j)\), \(i,j\ge0,\ i+j\le m(m+1)\), and
interpolating in the product binomial basis determines the
whole polynomial. The standard-library replay below performs
this calculation directly. These determinants are strictly
positive for \(u\ge0,\beta>0\), completing the assertion.
This extends the boundary result for
`GAP-FACTORIAL-BINARY-LINKED-BETA-COPRIMALITY` to this shifted
shape. The actual univariate Gamma families and the
two-positive unequal-exponent linked beta problem remain
separate.



## 5. Two Gaussian blocks, fixed markers and sharp nilpotence

Take independent standard, variance-one Gaussian coordinates
\(x_1,x_2,y_1,\ldots,y_d\), \(d\ge1\), and put
\[
 S=(x_1^2+x_2^2)/2,\qquad T=(y_1^2+\cdots+y_d^2)/2,\qquad
 H=T^u(A S+B T),\qquad u\in\mathbb N_0.
\]
Then \(S,T\) are independent Gamma variables of shapes \(1,d/2\)
at rate one. In particular, Section 4 rules out an ordinary
Gaussian zero pair in every two consecutive positive orders
for \((A,B)\ne(0,0)\). The same result holds with any fixed
homogeneous polynomial marker \(C(y)\) of even degree \(2j\)
and nonzero Gaussian mean. The restriction to the \(y\)-block
preserves the shape-one \(S\)-variable.

The exact moment identity for this marker is
\[
 \mathbb E[C(y)H^m]
   =\mathbb E[C(y)]\,\mathbb E_{1,\ d/2+j}
                  [T^{um}(A S+B T)^m].                           \tag{7}
\]
On the right, the subscripts are the shapes of \(S,T\). Indeed,
Gaussian rescaling of the \(y\)-block against \(e^{tT}\) gives
\(\mathbb E[C(y)e^{tT}]=\mathbb E[C(y)](1-t)^{-d/2-j}\)
for real \(t<1\); differentiation at zero gives the identity
for all polynomial powers of \(T\). For the monomial marker
\(C=y^{2\nu}\), the scalar is
\(\prod_{i=1}^{d}(2\nu_i-1)!!\), with \((-1)!!=1\), and
\(j=|\nu|\). A marker with an odd \(y_i\) exponent is killed
by the corresponding sign change. This is the same shape shift as
`FACT-PRIME-POWER-MARKED-RADIAL-SFC-LEM-RADIAL-TILT`, used here
in a phase depending on two independent blocks.

### Universal coefficient statement

Let \(\mathcal R\) be any unital commutative \(\mathbb Q\)-algebra
and \(C\) have rational coefficients with \(\mathbb E[C]\ne0\).
Use the Gaussian moment map coefficientwise over \(\mathcal R\).
If the literal marked moments in orders \(m,m+1\) vanish, then
\[
                  (A,B)^{2m}=0,\qquad H^{2m}=0.                    \tag{8}
\]
The coefficient-ideal exponent is sharp for each fixed order
pair \(m,m+1\), fixed \(u,d\) and eligible marker.

To prove this, take universal symbols \(a,b\) over \(\mathbb Q\).
After removing the rational units
\(\mathbb E[C]\,\ell!(d/2+j)_{u\ell}\) at each order \(\ell\),
the moment forms are
\[
        F_\ell(a,b)
         =\sum_{k=0}^{\ell}
                \frac{(u\ell+d/2+j)_k}{k!}\,a^{\ell-k}b^k .
\]
Section 4 shows that \(F_m,F_{m+1}\) are coprime in
\(\mathbb Q[a,b]\), since they have no common projective complex
zero. Hence they form a regular sequence. The standard
two-form Hilbert calculation gives
\[
 \operatorname{Hilb}\bigl(\mathbb Q[a,b]/(F_m,F_{m+1}),t\bigr)
       =\frac{(1-t^m)(1-t^{m+1})}{(1-t)^2}.
\]
Its length is \(m(m+1)\), and its top degree is \(2m-1\)
with dimension one. Thus \((a,b)^{2m}\) lies in the moment
ideal and \((a,b)^{2m-1}\) does not. Evaluate the universal
forms at \(a=A,b=B\) to get the first part of (8).
For every \(N\ge1\), the coefficients of \(H^N\) generate
exactly \((A,B)^N\). For each \(k=0,\ldots,N\), the monomial
\(x_1^{2(N-k)}y_1^{2(uN+k)}\) has coefficient
\(2^{-(u+1)N}\binom Nk A^{N-k}B^k\), a rational unit times
that coefficient monomial. The first-block degree separates
these terms, and the reverse inclusion follows by expansion.
This proves literal nilpotence and its sharpness by taking
\(\mathcal R\) to be the universal moment quotient itself.
In particular every polynomial-marked moment of \(H^n\)
vanishes for \(n\ge2m\) if the two selected moments vanish.
Over a reduced coefficient ring, the two equations force
\(A=B=0\).

The regular-sequence argument is the same fixed-pair mechanism
as in [the retained coefficient theorem, Section 1](SUPPORT_UNIFORM_COEFFICIENT_NILPOTENCE_AND_RADIAL_GAUSSIAN.md),
`FACTORIAL-SUPPORT-UNIFORM-COEFFICIENT-NILPOTENCE-AND-RADIAL-GAUSSIAN`.
The receiver here is an all-order class with two radial blocks.
Its two moments are different from the ordinary univariate
binomial pair used in that node. For arbitrary polynomial
markers the boundary documented by
`FACT-MARKED-RADIAL-UNIVERSAL-POLYNOMIAL-MARKER-EXTENSION`
persists here. At \(u=0,\ d=2,\ H=T\), the marker
\(T^2-6T+6\) has mean \(2\) and kills the moments at orders
one and two, by the factorial moments of the Gamma1 law.

### One-sided additions and a power of the Laplacian

Use first-block Witt coordinates \(z,w\), so \(S=zw\), and
assign weights \(1,-1\) to \(z,w\) and weight zero to all
\(y_i\). The Gaussian map kills nonzero weights.
If \(J\) has all its weights at least \(\gamma\ge1\), then
\(P=H+J\) has exactly the same moments as \(H\) with any
marker \(C(y)\) in (7): every term containing \(J\) has
positive weight. This is the vacuum and weight-zero case of
`GAUSSIAN-ANGULAR-CORE-ACTUAL-MODULE-INVARIANCE` in
[the retained angular-core theorem, (A3)–(A5)](../../moments/vanishing/GAUSSIAN_ANGULAR_CORE_FULL_IMAGE_INVARIANCE.md).
Over \(\mathbb C\), purity of \(P\) and eventual vanishing for
every fixed polynomial marker are each equivalent to \(H=0\).
For nonzero \(H\) the unmarked sequence has a nonzero moment
in every adjacent pair. If \(H=0\), then \(QJ^n\) has strictly
positive weight as soon as \(n\gamma>-\min\mathrm{wt}(Q)\),
for each fixed nonzero marker \(Q\).

More generally, the two vanishing moments in (7) for \(P\)
imply (8) for the core over any \(\mathbb Q\)-algebra in
Wick coordinates. In \(QP^n\), the terms with more than
\(\max(0,\lfloor-\min\mathrm{wt}(Q)/\gamma\rfloor)\) copies
of \(J\) have positive weight. All other terms contain
\(H^{2m}\) once \(n\) is at least \(2m\) plus that bound.
This gives a literal fixed-marker cutoff for \(P\); the
one-sided summand itself need not be nilpotent.

If \(J\) is also homogeneous of spatial degree \(2(u+1)\),
let \(\mathscr D=2\partial_z\partial_w+\sum_i\partial_{y_i}^2\).
The homogeneous Gaussian contraction formula
\(\mathbb E[P^n]=\mathscr D^{(u+1)n}P^n/
 (2^{(u+1)n}((u+1)n)!)\) gives a scoped generalized-vanishing
consequence for the operator \(\mathscr D^{u+1}\).
Indeed, zero contractions of \(P^m,P^{m+1}\) force (8) for the
core. For a polynomial marker \(Q\) of degree \(b\) and
minimum weight \(q\), the contraction
\(\mathscr D^{(u+1)n}(QP^n)\) is zero for
\(n\ge2m+\lfloor(b-q)/\gamma\rfloor\).
Each surviving monomial after differentiation has degree
at most \(b\); \(\mathscr D\) preserves angular weight, so
terms with \(q+j\gamma>b\), where \(j\) is the number of
copies of \(J\), are killed. For the other terms the given
bound ensures at least \(2m\) copies of \(H\). If \(J=0\) the
bound is simply \(n\ge2m\). This conclusion applies to the
specified homogeneous phases and does not use the general
dimension classification `MOMENTS-GVC-CLASSIFICATION`.



## 6. Regular Euclidean chains on diagonal profiles

For every integer \(r=u\ge0\) and every \(n\ge1\), the actual
linked pair \(J_n,J_{n+1}\) has a regular monic Euclidean chain:
its remainders have degrees \(n-1,n-2,\ldots,0\).
This proves the diagonal restriction of
`GAP-FACTORIAL-LINKED-BETA-REGULAR-CHAIN`.

Set \(\lambda=nr+1\), and define \(C_j^\lambda(x)\) by
\[
 (1-2xv+v^2)^{-\lambda}=\sum_{j\ge0}C_j^\lambda(x)v^j.
\]
The physical coefficient polynomial is
\(g_n(s)=[v^n](1-v)^{-\lambda}(1-sv)^{-\lambda}\).
For \(s=e^{2i\theta}\) with \(0<\theta<\pi\), substitution of
\(e^{i\theta}v\) in the generating series gives
\(g_n(s)=e^{in\theta}C_n^\lambda(\cos\theta)\).
The upper order uses \(\lambda+r\). We show strict
interlacing in \(\theta\) despite this change of parameter.

Differentiating the series in \(v\) gives
\[
 (j+1)C_{j+1}^\lambda
 =2(j+\lambda)x C_j^\lambda-(j+2\lambda-1)C_{j-1}^\lambda.
\]
Taking leading terms gives the monic Jacobi recurrence with
zero diagonal and offdiagonals \(a_j(\lambda)\) where
\[
 a_j(\lambda)^2=\frac14\left(1-\frac{\lambda(\lambda-1)}
                  {(j+\lambda)(j+\lambda-1)}\right).
\]
For \(\lambda>1\), these are positive and below \(1/2\).
The row sums put the Jacobi spectrum strictly in \((-1,1)\).
The recurrence determines an eigenvector from its first
entry, so the symmetric matrix has simple eigenvalues.
For \(r=0\), the recurrence gives
\(C_n^1(\cos\theta)=\sin((n+1)\theta)/\sin\theta\). The
consecutive lists of equally spaced interior zeros then
strictly interlace.

Assume \(r\ge1\), and put \(N=n+\lambda\).
The Gegenbauer differential equation follows by applying
\((1-x^2)\partial_x^2-(2\lambda+1)x\partial_x+
(v\partial_v)^2+2\lambda v\partial_v\) to the generating
series, which gives zero. In angle coordinates, set
\(U_n=(\sin\theta)^\lambda C_n^\lambda(\cos\theta)\). Then
\[
 U_n''+Q_nU_n=0,\quad
 Q_n=N^2-\lambda(\lambda-1)\csc^2\theta,\quad
 Q_{n+1}-Q_n=(r+1)(2N+r+1)-r(2\lambda+r-1)\csc^2\theta.
\]
Let \(\sin^2\theta_0=z_0=r(2\lambda+r-1)/
((r+1)(2N+r+1))\), \(0<\theta_0<\pi/2\), and take the
turning point \(\theta_t\in(0,\pi/2)\) of \(Q_n\).
For \(d_0=2nr+r+1\), exact substitution gives
\[
 Q_0=Q_n(\theta_0)=
 \frac{n^2(r^2-1)+n(r^2+2r-1)+r+1}{d_0},\qquad
 (N-Q_0)d_0=n^2(r+1)^2+2n(r+1).
\]
Thus \(0<Q_0<N\) and \(\theta_t<\theta_0\), including
\(n=r=1\), where \(Q_0=1\) and \(N=3\). Also
\(\sin\theta_t=\sqrt{\lambda(\lambda-1)}/N
\ge r/(r+2)\ge1/3\), and \(z_0\le r/(r+1)\) gives
\(\cos\theta_0\ge1/\sqrt{r+1}\).
Since \(z_0-\sin^2\theta_t=z_0 Q_0/N^2\), the mean-value
theorem applied to \(\sin^2\theta\) gives
\[
 (\theta_0-\theta_t)\sqrt{Q_0}
 \le \frac32\sqrt{r+1}\,\frac{z_0 Q_0^{3/2}}{N^2}
 <\frac32\sqrt{(r+1)/N}<\frac\pi2 .
\]
This bound holds uniformly as \(\lambda\) grows.

Near zero, \(U_n\) is a positive constant times
\(\theta^\lambda(1+O(\theta^2))\), since
\(C_n^\lambda(1)=(2\lambda)_n/n!>0\).
It is positive with positive derivative up to \(\theta_t\):
while \(U_n>0\), the equation gives \(U_n''\ge0\) there.
Start at \(\theta_t\) with this positive value and derivative.
The solution for the constant potential \(Q_0\) stays
positive up to \(\theta_0\) by the strict angle bound.
The variation-of-constants integral for
\(U_n''+Q_0U_n=(Q_0-Q_n)U_n\) has a nonnegative sine
kernel and nonnegative integrand up to any putative first
zero. Hence \(U_n\) stays positive through \(\theta_0\).

By symmetry, all lower gaps lie in
\((\theta_0,\pi-\theta_0)\), where \(Q_{n+1}>Q_n\).
Each consecutive pair of lower zeros contains an upper zero.
Indeed, if the upper solution has no interior zero, choose
both solutions positive in the gap. Then
\((U_n' U_{n+1}-U_n U_{n+1}')'
=(Q_{n+1}-Q_n)U_n U_{n+1}>0\),
whereas this Wronskian is nonnegative at the left lower
zero and nonpositive at the right.

The extreme upper zeros lie strictly outside the lower ones.
For \(1\le j\le n-1\), the trailing \(n\times n\) principal
block of the upper Jacobi matrix has offdiagonal
\(a_{j+1}(\lambda+r)>a_j(\lambda)\). In the formula for
\(a_j^2\), the denominator in the subtracted term equals
\((1+j/\lambda)(1+j/(\lambda-1))\) after normalization.
Both factors strictly increase when \(j\) becomes \(j+1\)
and \(\lambda\) becomes \(\lambda+r\), since
\(\lambda>jr\) and \(\lambda-1>jr\).
For \(n\ge2\), the Rayleigh quotient on the positive top
eigenvector of the lower matrix shows that the trailing
block has a strictly larger largest eigenvalue. A positive
top vector exists because replacing entries by absolute
values does not decrease the quotient, and its eigenvector
equation prevents any zero entry. The full upper matrix
has largest eigenvalue at least that of its trailing block.
For \(n=1\) the upper eigenvalue is already positive and
the lower root is zero. Symmetry gives the other endpoint.
Together with the interior comparisons this proves strict
interlacing.

Finally, the [actual linked beta formula](canonical-density/moment_pair_obstructions.md)
takes \(s=z/(z-1)\) from the coefficients \(z-1,z\).
The point \(s=1\) is not a root since the coefficients of
\(g_n\) and \(g_{n+1}\) are positive. For \(s=e^{2i\theta}\),
\(z=(1-i\cot\theta)/2\), so the \(J\)-roots interlace on
the line \(\Re z=1/2\). The fixed affine change
\(z=1/2+ix\) and monic rescaling give real monic
polynomials and preserve Euclidean degrees. At the lower
roots the first remainder takes alternating nonzero
values of the upper polynomial. Thus it has degree
\(n-1\) and interlaces the lower polynomial, and this
argument repeats down to a constant. Unequal profiles
retain a separate regular-chain obligation.



## 7. Other factorial obligations

The fixed-shape auxiliary production matrix in the
[cyclic Wronskian note](gap_three_hurwitz_structure.tex) already
has a totally nonnegative bidiagonal factorization on its stated
range. The actual \(H_3\) construction uses cyclic sections of a
rotated Wronskian, multiplying auxiliary minors by the signed
periodic weights \(1,1,0,-1,-1,0\). The positive-measure composition in
Section 1 gives no comparison between the positive and negative
terms. Thus `FACTORIAL-GAP3-HURWITZ-TN`,
`FACTORIAL-GAP3-TERMINAL-MINOR-STRICT` and
`FACTORIAL-LOWER-ONE-TERMINAL-DETERMINANT` stay open.
The explicit eventual cutoff for the ordinary univariate binomial
concerns a different support family.

## 8. Cumulant tails and zero-geometry tests

The [Nuttall proof, Sections 2–6](https://github.com/1337isnot1337/solved-conjectures/blob/758066055ef49c2c66cb1087918ae62ad75e7ac1/results/nuttall-cumulants/PROOF.md)
derives signs for the theta cumulants for all sufficiently large
integration orders, at each fixed minor order. Its transferable
step is the following determinant identity. For a smooth
rapidly decaying function \(\phi\) and its \(m\)-fold tail primitive
\(\Psi_m\), with \(m\ge h\),
\[
 \det[\Psi_m^{(i+j)}(x)]_{i,j=0}^{h-1}
 =
 \frac{\int_{0<t_1<\cdots<t_h}
  \prod_i t_i^{m-h}\,\prod_{i<j}(t_j-t_i)\,
  \det[\phi^{(j-1)}(x+t_i)]_{i,j=1}^{h}\,d\boldsymbol t}
 {\prod_{i=0}^{h-1}(m-i-1)!}.
\]
Integrate by parts \(i\) times in row \(i\). The boundary at zero
vanishes before each integration by parts. The resulting row signs
\((-1)^i\) cancel exactly the sign reversing \(h\) descending
powers of the integration variable in Cauchy–Binet. Both determinants
change sign under the same permutation, leaving a symmetric product.
The \(h!\) cancels on restriction to ordered variables.

In the original theta application, the first theta summand makes
\(\phi(t)\) asymptotic to a positive constant times
\(y^{9/4}e^{-y}\), \(y=\pi e^{4t}\), with a differentiable relative
error \(O(y^{-1})\). The Hankel-Wronskian identity
\(\tau_{h+1}\tau_{h-1}=\tau_h\tau_h''-(\tau_h')^2\) then gives
alternating tail signs: its leading logarithmic second derivative
is \(-16h y\). A determinant factorization evaluates a separated
box \(t_i=\tfrac18\log m+O_h(1)\) with the same sign and size bounded
below by a positive multiple of
\((\tfrac18\log m)^{h(m-h)}e^{-C\sqrt m}\) before the positive scalar.
The possibly opposite-sign region has some \(t_i\) in a fixed
interval and contributes at most
\(C U^m(2\log m)^{m(h-1)}\). The ratio tends to zero because a full
factor \((\log m)^m\) outweighs every fixed-base exponential in
\(m\). These statements use the theta tail and its derivative bounds,
and the compact residual range in \(x\) is uniform because a common
\(x\) multiplies every \(y_i\) by the same positive factor.
The method would need an actual tail-primitive representation with
these quantitative controls for a factorial determinant. The
ordinary-binomial explicit-cutoff obligation has no such identified
representation, and the identity alone gives no cutoff there.

The [Csordas–Vishnyakova proof](https://github.com/1337isnot1337/solved-conjectures/blob/758066055ef49c2c66cb1087918ae62ad75e7ac1/results/csordas-vishnyakova/PROOF.md)
has a directly verifiable two-parameter extension. Let
\(p(z)=z^4+A z^2+B\), where \(B>0\) and \(A\ge0\), and put
\(\varphi(z)=e^{-z^4}\). The roots of \(p\), equivalently their
negatives, have second and fourth power sums \(-2A\) and
\(2A^2-4B\), with odd power sums zero by evenness. Therefore the
product over their shifts is
\[
 \prod_{p(\alpha)=0}\varphi(z+\alpha t)
   =e^{-4z^4}\exp(12A z^2t^2)\exp((4B-2A^2)t^4).
\]
All its Taylor coefficients in \(t\) are nonnegative on the real
\(z\)-axis exactly when \(A^2\le2B\), the necessity following from
the \(t^4\) coefficient at \(z=0\). Yet
\(\varphi(it)\varphi(-it)=e^{-2t^4}\) has a negative coefficient.
A real-rooted approximating polynomial gives a product on the
imaginary axis with only nonnegative coefficients, since each real
root \(r\) contributes \(t^2+r^2\). Compact convergence preserves
each coefficient by Cauchy's formula. This verifies the failed
implication for a full range of test polynomials and excludes a
route from those coefficient inequalities alone to global
Laguerre–Pólya membership.

The [Cardon et al. counterexample](https://github.com/1337isnot1337/solved-conjectures/blob/758066055ef49c2c66cb1087918ae62ad75e7ac1/results/cardon-zero-sector/PROOF.md)
transfers its infinitesimal cubic test to the **actual** factorial
multipliers in the [Grace strip question](ITERATED_GRACE_POLYNOMIAL_ROOT_STRIP.md).
For positive weights \(\gamma_k\), let \(D_t(z^k)=\gamma_k^t z^k\).
On \(z^j(z-1)^3\) first remove the \(z^j\) factor and translate the
cubic to the origin. Its coefficients from degree three to zero
are \(1+O(t),O(t),O(t),t\Delta^3\log\gamma_j+O(t^2)\).
The cubic discriminant is
\[
 -27t^2(\Delta^3\log\gamma_j)^2+O(t^3).
\]
Hence any nonzero third difference forces a nonreal conjugate
pair at all sufficiently small positive \(t\). For the original
Cardon weights \(\log\gamma_k=-k^p\), \(p>2\), the third difference
at zero is strictly negative by integrating the strictly increasing
second derivative of \(x^p\) over translated unit squares. At
\(p=3\) and \(t=\log(16/15)/6\), scaling the variable as in the
external proof gives \(r^4z^3-3rz^2+3z-1\), \(r=15/16\), whose
discriminant is \(-72517275/4294967296\).

For the actual degree-\(n\) factorial factor take
\(\gamma_k=(an+1)_{dk}\), integers \(d\ge1\), \(a\ge0\), and
\(n\ge3\). For \(j+3\le n\) the third difference is strictly negative
by finite products alone:
\[
 \Delta^3\log\gamma_j
 =
 \sum_{\ell=1}^d
 \log\frac{(an+dj+\ell)(an+dj+\ell+2d)}
               {(an+dj+\ell+d)^2}<0 .
\]
Take \(j=n-3\) to make the witness degree exactly \(n\).
The input zeros are \(0,1\), in every positive-width sector about
the nonnegative ray. A nonreal output zero excludes any uniform
sector reduction, by taking initial width below its argument.
The cubic gives the same failure for the inverse fractional
multipliers, since changing the sign of \(\log\gamma\) does not
change the leading squared coefficient. Thus continuous
interpolation by positive fractional powers of the physical
multiplier cannot preserve real-rootedness or reduce sectors
throughout the interpolation. This does not address the direct
finite-order Grace strip or its actual common-root criterion
`GAP-FACTORIAL-POLYNOMIAL-GRACE-STRIP-NONCOLLISION`.



## Verification

The pinned KVZ diagnostic passed all 115 exact rational checks,
and all 20 recorded source hashes matched the pinned revision.
Here, exact rational arithmetic matched (1) on 32 test grids and
(2a) on 120 test grids, including positive, negative and zero
determinants. Rational elimination also matched the skipped-column
classification on 23,792 minors, using rows from 7 through 12,
columns from 0 through 7, orders 1 through 5, and the profiles
((1,1),(1,2),(2,4),(2,6),(1,4),(4,0),(0,4),(0,0)).
The exact first-failure orders also agreed with 1,915 rational
determinant controls on the profiles (2,3), (2,5), (3,10), (4,19),
(3,2), and (5,2).
Direct polynomial expansion matched (4) on 49 pairs
of integer parameters. The determinant proofs above give the
universal identities.

Direct Wick expansion matched 72 two-block marker identities.
Exact quotient calculations confirmed sharp nilpotence on 36 pairs
and checked the angular and Laplacian cutoffs in a nonreduced example.
Rational division verified 96 diagonal Euclidean chains at
r=0,1,2,5,16,1000 and n=1 through 16. Symbolic expansion checked
the Gegenbauer generating differential identity and the displayed
identities for the potential at the crossover angle.

The finite interpolation in (6) was checked at all 23 points by
rational elimination on Sylvester matrices. An independent short
replay using the repository's SymPy environment is

```python
from math import comb
from sympy import Poly, Symbol, resultant

s = Symbol("s")
quotients = (
    (1, [1]),
    (4, [4, 22, 37, 17]),
    (432, [432, 7200, 50496, 193890, 444511, 622063, 518466, 235223, 44519]),
)
for m, (denominator, coefficients) in enumerate(quotients, start=1):
    for u in range(m * (m + 1) + 1):
        lower = Poly(sum(comb(m * u + j, j) * s**j for j in range(m + 1)), s)
        upper = Poly(
            sum(comb((m + 1) * u + j, j) * s**j for j in range(m + 2)), s
        )
        expected = (u + 1)**(m + 1) * sum(
            coefficient * u**j for j, coefficient in enumerate(coefficients)
        )
        assert denominator * resultant(lower, upper, s) == expected
```

The shifted-shape certificates use 125 exact bivariate interpolation
values in total. The following replay uses rational elimination and
ordinary-coefficient reconstruction by two-variable finite differences.

```python
from fractions import Fraction
from math import comb, prod

def rising(x, length):
    return prod(x + k for k in range(length))

def sylvester_at(m, u, beta):
    n = m + 1
    def row(order):
        return [Fraction(rising(order*u + beta, j), prod(range(1, j+1)))
                for j in range(order, -1, -1)]
    lower, upper = row(m), row(n)
    matrix = ([[0]*i + lower + [0]*(n-1-i) for i in range(n)] +
              [[0]*i + upper + [0]*(m-1-i) for i in range(m)])
    determinant = Fraction(1)
    for j in range(m+n):
        pivot = next((i for i in range(j, m+n) if matrix[i][j]), None)
        if pivot is None:
            return Fraction(0)
        if pivot != j:
            matrix[j], matrix[pivot] = matrix[pivot], matrix[j]
            determinant = -determinant
        determinant *= matrix[j][j]
        for i in range(j+1, m+n):
            factor = matrix[i][j] / matrix[j][j]
            for k in range(j+1, m+n):
                matrix[i][k] -= factor * matrix[j][k]
    return determinant

binomial_polynomials = [[Fraction(1)]]
for j in range(12):
    p = binomial_polynomials[-1]
    binomial_polynomials.append(
        [(a - j*b)/(j+1) for a,b in zip([0]+p,p+[0])]
    )

for m, count, denominator in [(1,5,2), (2,25,36), (3,85,13824)]:
    degree = m*(m+1)
    indices = [(i,j) for i in range(degree+1) for j in range(degree+1-i)]
    values = {(i,j): sylvester_at(m,i,j) for i,j in indices}
    coefficients = dict.fromkeys(indices, Fraction(0))
    for i,j in indices:
        difference = sum(
            (-1)**(i+j-a-b)*comb(i,a)*comb(j,b)*values[a,b]
            for a in range(i+1) for b in range(j+1)
        )
        for a,c in enumerate(binomial_polynomials[i]):
            for b,d in enumerate(binomial_polynomials[j]):
                coefficients[a,b] += difference*c*d
    assert all(value >= 0 for value in coefficients.values())
    assert sum(value > 0 for value in coefficients.values()) == count
    assert coefficients[0,degree] == Fraction(1,denominator)
    print(m, len(values), count, coefficients[0,degree])
```
