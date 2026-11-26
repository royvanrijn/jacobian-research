# Common shifted numerators over a prescribed prime pool

The quantitative Subspace Theorem gives a uniform curve cover for two
original abc ratios whose shifted numerators have a large common divisor
and whose outside-pool numerical heights are small. It also gives a
one-row count from the combined numerical height on all three letters.
The pool can move with the height. The proofs below use Evertse's theorem
directly.

## Notation and quantitative input

Let \(S\) consist of the real place and a finite set of primes, and
write \(s=|S|\). For a positive integer \(N\), put
\(N_{\overline S}=\prod_{p\notin S}p^{v_p(N)}\).
For nonzero rationals \(u,v\), define their outside-pool numerical height by
\[
 E_S(u,v)=\sum_{p\notin S}
 \left(\max(0,\log|u|_p,\log|v|_p)+
       \max(0,-\log|u|_p,-\log|v|_p)\right).
\]
Write \(E_S(u)\) for the one-coordinate version. For reduced positive
\(u=a/b\), it is exactly \(\log (ab)_{\overline S}\).
For logarithmic heights use \(h(u)=\log\max(|a|,|b|)\) and
\(h([1:u:v])=\sum_\nu\log\max(1,|u|_\nu,|v|_\nu)\).

For a rational vector \(P\), let \(\|P\|_\nu\) be the sup norm at
finite places and Euclidean at the real place. Put
\(h_E(P)=\sum_\nu\log\|P\|_\nu\). Use the following special case of
Evertse's quantitative Subspace Theorem. In fixed dimension \(r\), if
at every \(\nu\in S\) the \(r\) rational linear forms \(L_{i,\nu}\)
have bounded coefficient height and coefficient determinant of absolute
value 1 at that place, the solutions of
\[
 \sum_{\nu\in S}\sum_i\log\frac{|L_{i,\nu}(P)|_\nu}{\|P\|_\nu}
       < -(r+\delta)h_E(P)
 \qquad (0<\delta<1)
 \tag{AP0}
\]
lie in \(\exp(O_{r,\delta}(s))\) proper rational linear subspaces, with
the bounded coefficient height absorbed in the constant. Both the
large-height and small-height parts are available, so there is no
height cutoff. Here all forms have coefficients in \(\{0,\pm1\}\),
and all coefficient determinants are \(\pm1\).
[Evertse, An improvement of the quantitative Subspace theorem, main theorem (i)–(ii), p. 228](https://numdam.org/item/CM_1996__101_3_225_0.pdf#page=5)

## Numerical capture in one original row

**AP1.** Fix \(0<w<1\). Uniformly in \(S\), for all \(X\) sufficiently
large in terms of \(w\), the number of positive primitive rows
with \(a+b=c\), \(X/2<c\le X\), and
\[
                   \log(abc)_{\overline S}\le (1-w)\log X
 \tag{AP1}
\]
is at most \(\exp(O_w(s))\).

Use \(P=(a,b)\) in (AP0). At the real place choose the forms \(a,b\).
At each \(p\in S\) choose the smallest in \(p\)-adic absolute value
among \(a,b,a+b\) and any other of the three. The two letters other
than a nonunit are units because the row is primitive. Each pair of
distinct forms is unimodular. The product of the numerators in (AP0)
is therefore
\[
                \frac{ab}{(abc)_S}
                =\frac{(abc)_{\overline S}}{c}
                \le 2X^{-w}.
\]
All finite norms of \(P\) are 1. Since \(h_E(P)\le\log X+O(1)\),
(AP0) follows with \(\delta=w/2\) for large \(X\). The at most \(3^s\)
form patterns only affect the constant in the exponential. A proper
linear subspace of \(\mathbb Q^2\) contains at most one positive
primitive pair on its ray.

## Common gcd of the third letters

**AP2.** Fix \(0<w<1\) and set
\[
 m=\lceil64/w\rceil,\qquad \rho=\frac{w}{48m}.
\]
For every \(S\) and every \(X\) sufficiently large in terms of \(w\),
the pairs of positive primitive rows
\[
 u=a/b,\quad v=a'/b',\qquad
 a+b=c,\quad a'+b'=c',\qquad X/2<c,c'\le X,
\]
with
\[
       E_S(u,v)\le \rho\,h([1:u:v]),\qquad
       \gcd(c,c')\ge X^w
 \tag{AP2}
\]
lie in a union of at most \(\exp(O_w(s))\) zero sets of nonzero
polynomials in \(u,v\) of degrees at most \(m\). The polynomials and
\(S\) may depend on \(X\). The estimate is uniform in the primes of
\(S\).

Put \(T=\log X\) and \(E=E_S(u,v)\). Then \(h(u),h(v)\le T\) and
\(E\le 2\rho T=wT/(24m)\). Since \(u+1=c/b\) and \(v+1=c'/b'\)
are in lowest terms, the local expression
\[
 g_p=-\min(0,\log\max(|u+1|_p,|v+1|_p))
        =\min(v_p(c),v_p(c'))\log p
 \tag{AP3}
\]
is exact. The corresponding real term is zero since both shifted
values exceed 1.

If \(\sum_{p\in S}g_p\ge wT/2\), apply (AP0) to \(P=(u,1)\).
At \(p\in S\) with \(|u+1|_p<1\), replace the coordinate form \(u\)
by \(u+1\). It was a unit there and the replacement gains
\(-\log|u+1|_p\). Use standard coordinate forms at all other places
of \(S\). The product formula for \(u\) shows that the value
of the left side in (AP0) before the replacement is
\[
        -2h_E(P)+E_S(u).
\]
Indeed, the contribution gained by removing the places outside \(S\)
is \(\sum_{p\notin S}(2\log\max(1,|u|_p)-\log|u|_p)\).
The replacement gains at least \(wT/2\), and \(E_S(u)\le E\).
Thus (AP0) holds with \(\delta=w/8\) for large \(X\).
After considering at most \(2^s\) choices of forms, \(u\) lies in a
set of size \(\exp(O_w(s))\), giving vertical lines.

Otherwise \(G=\sum_{p\notin S}g_p\ge wT/2\). Let \(M_i\) be the
\(r+1=\binom{m+2}{2}\) monomials of total degree at most \(m\) in
two variables, including \(M_0=1\). Put \(s_i=M_i(-1,-1)\) and
\(P=(M_i(u,v)-s_i)_{i\ne0}\). The polynomials forming these
coordinates are a basis for the degree-at-most-\(m\) polynomials
vanishing at \((-1,-1)\). In particular \(P\ne0\). At each place
\(\nu\in S\), choose a monomial \(M_z\) of smallest absolute value
at \((u,v)\), and for each \(i\ne z\) use the linear form on \(P\)
corresponding to
\[
                       M_i-s_i s_z M_z .
\]
These form an integral basis of the same kernel of evaluation.
For \(z\ne0\), the form indexed by \(i=0\) recovers the \(z\)-th
old coordinate up to sign, and adding its multiples recovers the
remaining old coordinates. This proves determinant \(\pm1\)
and checks the coefficient bound.

Let \(\sigma_\nu=\max(|\log |u|_\nu|,|\log |v|_\nu|)\). With a factor
of 2 at the real place only, the \(i\)-th linear form is at most
\(|M_i(u,v)|_\nu\) in absolute value. The full monomial product
formula, with the chosen monomial omitted, implies
\[
 \sum_{\nu\in S}\sum_{i\ne z(\nu)}\log|L_{i,\nu}(P)|_\nu
 \le r\log2+(r+1)m E+4mT.
\]
In detail, the omitted term costs at most
\(m\sum_{\nu\in S}\sigma_\nu\le 2m(h(u)+h(v))\le4mT\).
The monomial logs outside \(S\) have absolute sum at most
\((r+1)mE\). If a form vanishes, the corresponding nonzero
polynomial itself gives one of at most \((r+1)^2\) exceptional
curves, and the remaining calculation may use finite logarithms.

For \(p\notin S\), each \(M_i-s_i\) can be written as
\((u+1)A_i+(v+1)B_i\) with integral polynomials of degree at most
\(m-1\), by first substituting \(-1\) for \(u\) and then for \(v\).
If \(g_p>0\), both \(u,v\) are \(p\)-adic units and that
identity gives \(\|P\|_p\le e^{-g_p}\).
If \(g_p=0\), use
\(\|P\|_p\le\max(1,|u|_p,|v|_p)^m\).
Thus \(\sum_{p\notin S}\log\|P\|_p\le mE-G\).
Combining these bounds and subtracting the norm terms in (AP0)
gives at most
\[
              -r h_E(P)+r\log2+(2r+1)mE+4mT-rG.
\]
Here \((2r+1)mE\le rwT/8\) and \(4m/r\le 8/m\le w/8\).
For large \(X\) the bound is at most \(-r h_E(P)-rwT/8\).
Also \(h_E(P)\le m h([1:u:v])+\log(2\sqrt r)\le3mT\),
so (AP0) holds with \(\delta=\min(1/2,rw/(32m))\).
The number of choices of local form bases is at most \((r+1)^s\).
A nonzero form defining a hyperplane containing any proper
exceptional subspace pulls back to a nonzero polynomial of degree
at most \(m\), by independence of the original coordinates of \(P\).

The same conclusion holds for the gcd of the reduced numerators
of \(u-\alpha,v-\beta\) for fixed nonzero rational shifts, including
distinct shifts, with constants and the lower threshold for \(X\)
allowed to depend on the shifts and with \(\rho\) decreased by an
absolute factor. Omit exact zeros of the shifted values as whole
vertical or horizontal lines. Substitute \(U=-u/\alpha\) and
\(V=-v/\beta\), and add the prime supports of the shifts to \(S\).
The outside-pool height weakly decreases. The scalar and joint
heights change by \(O_{\alpha,\beta}(1)\). The logarithmic gcd of
reduced numerators changes by at most a constant because each
of the positive valuations in the definition changes by at most
the absolute valuation of the fixed scalar. The preceding proof
works for signed \(U,V\) with the same finite-place calculations
and with a slightly weakened \(w\) to absorb these constants.
For arbitrary nonzero rational shifted values, the all-place
logarithmic gcd also has the real term
\(\max(0,-\log\max(|u-\alpha|,|v-\beta|))\).
AP2 only counts the integer numerator gcd.

## Counting and the endpoint hypothesis

If \(\mathcal F\) is any family of original rows in the shell on which
\(E_S(a/b)\le(\rho/3)\log X\), all pairs in \(\mathcal F\) satisfy
the height condition of AP2 for sufficiently large \(X\).
Indeed \(E_S(u,v)\le E_S(u)+E_S(v)\) and
\(h([1:u:v])\ge h(u)\ge\log(X/4)\).
A nonzero polynomial of degree at most \(m\) vanishes on at most
\(2m|\mathcal F|\) ordered pairs of ratios from \(\mathcal F\).
At most \(m\) first coordinates can make the specialization
identically zero. Every other fibre has at most \(m\) points.
The ratio determines the positive primitive row. Consequently
\[
 \#\{(x,y)\in\mathcal F^2:\gcd(c_x,c_y)\ge X^w\}
                   \le \exp(O_w(s))|\mathcal F|.
 \tag{AP4}
\]
The estimates are slicewise in \(X\), with no restriction on the order
of taking \(X\) and a radical parameter. If \(S\) is supported on any
prescribed integer at most \(X^A\) for fixed \(A\), then
\(s=O_A(\log X/\log\log X)\) by the product bound for distinct primes.
For fixed \(w,A\) and each \(\varepsilon>0\), the factor in (AP4)
is \(O_{w,A,\varepsilon}(X^\varepsilon)\). Equivalently the number
of pairs in (AP4) is \(O_{w,A,\lambda,\eta}(X^{\eta\lambda}|\mathcal F|)\)
for every fixed \(\lambda,\eta>0\).

For a fixed integer \(M\ge X^w\) and \(\tau\ge0\), there is also
an elementary bound for the number of original rows with \(M\mid c\)
and \(E_S(a/b)\le\tau\log X\). Write \(a=a_S a_{\overline S}\)
and \(b=b_S b_{\overline S}\). The two parts \(a_S,b_S\le X\)
have \(\exp(o(\log X))\) possibilities uniformly in pools supported
on integers at most \(X^A\). For instance, for any fixed \(0<t<1\)
the number of possible \(a_S\) is at most
\(X^t\prod_{p\in S}(1-p^{-t})^{-1}=X^{t+o(1)}\).
For each fixed \(a_S,b_S\), the parts of \(a,b\) are units modulo
\(M\), and \(a_{\overline S}/b_{\overline S}\equiv-b_S/a_S\bmod M\).
Group \(a_{\overline S}\) in dyadic intervals \([H,2H)\).
The inverse residual ratios lie in \((0,X^\tau/H^2]\), since
\(a_{\overline S}b_{\overline S}\le X^\tau\).
Distinct inverse ratios differ by at least \(M/(4H^2)\), by
divisibility of their cross-difference. The residual ratios
are reduced, so they recover the original rows. Summing over
the dyadic intervals and the supported parts gives
\[
                O_{A,\varepsilon}\bigl(X^\varepsilon(1+X^\tau/M)\bigr)
 \tag{AP5}
\]
rows uniformly in \(M\) and \(\tau\).

There is a corresponding anchored count for a moving pool. At each
fixed anchor \(x\), choose any \(S_x\) supported on an integer of
height at most \(X^A\). For each \(\varepsilon>0\) there are
\(O_{A,\varepsilon}(X^\varepsilon(1+X^{\tau-w}))\) query rows \(y\) with
\(E_{S_x}(a_y/b_y)\le\tau\log X\) and \(\gcd(c_x,c_y)\ge X^w\).
Indeed, this gcd is a divisor \(M\) of the complete anchor target.
There are \(O_{\varepsilon}(X^{\varepsilon/2})\) divisors to consider,
and (AP5) applies to each. The anchor and query can have different
outside numerical heights. AP2 provides a bounded-degree geometric
cover in addition to these cardinality bounds.

The numerical threshold for subpower growth is sharp without a
radical cap. For \(0<w<\tau<1\) set \(b=3^k\),
\(X=3b/2\), and \(M=5^{\lceil w\log_5 X\rceil}\).
With \(S=\{\infty,2,3\}\), take
\[
  a=2^{\lfloor(1-\tau)\log_2X\rfloor-2}j,\qquad
  1\le j\le X^\tau,\qquad 3\nmid j,\qquad a\equiv-b\bmod M.
\]
There are \(\gg X^{\tau-w}\) choices because the last congruence
sets a unit residue class in \(j\bmod M\) and only one in three
of its members are divisible by 3. We have \(a\le X/4\), so all
the rows are in the shell. They are primitive, \(M\mid c\),
and \(\log (ab)_{\overline S}\le\log j\le\tau\log X\).

The existing fixed-pool saving discounts radical mass measured on the
pool. AP1 uses the numerical mass of the full prime powers. For polynomial-height pools
its height threshold can also be handled by recording the three full supported
parts at subpower cost and using the three-factor determinant bound in dyadic
height boxes.
AP4 is a pair count using the outside numerical heights of two original
letters and the common divisor of the third. It holds without capping
radicals. An application to a saturated abc family still needs a
pool per anchor capturing the required numerical height on sufficiently
many queries having large third-letter common gcds with that anchor.
Allowing a separate pool for every query can remove the counting implication.
[Prescribed-pool radical saving (EB5)](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)
[Moving curves and secant gcds (BC2–BC5)](ABC_BOUNDED_DEGREE_CURVES_AND_SECANT_GCD.md)

## Scope of the preprint bound with moving coefficients

Noubissie's Theorem 9 and Proposition 4 in arXiv:2611.11341v1 give
bounds polynomial in the normalized logarithmic height of the joint
coefficient vector \([f,g,1]\), at fixed degrees, field and \(S\).
Theorem 9 expressly calls this a logarithmic height. The displayed
definition omits a logarithm, in which case its sum would diverge
for a vector containing 1. With the standard logarithmic height,
there are counterexamples to both uniform bounds.
[Noubissie, Theorem 9, Proposition 4, and Lemma 5](https://arxiv.org/html/2611.11341v1)

Fix \(S=\{\infty,2\}\), \(\epsilon=1/4\), any fixed allowed
\(\eta>0\), and \(0<\theta<\min(\eta/3,1/8)\).
Put \(R=2^k\), and take \(u=Ra,v=Rb\) with \(a,b\) odd positive
integers at most \(R^\theta\). Then \(E_S(u,v)\le2\theta\log R\),
\(h([1:u:v])\ge\log R\), and
\(\epsilon(h(u)+h(v))\le(1+\theta)\log R/2\).
For the coprime monic primitive polynomials \(f=x+R\), \(g=y+2R\),
the finite-place gcd is at least \(R\). Their shifts are distinct
\(S\)-units and \([f,g,1]\) has logarithmic height \(\log(2R)\).
This meets the full-gcd condition. For the outside-\(S\) assertion
take instead \(f=3^k(x+1)\), \(g=3^k(y+1)\) on the same grid.
They are coprime over \(\mathbb Q[x,y]\), their outside-\(S\)
gcd is at least \(3^k\), and the joint logarithmic coefficient
height is \(k\log3\).

For both examples the grid has \(\asymp R^{2\theta}\) points.
The fibre argument used for (AP4) bounds each nonzero
degree-at-most-\(m\) plane curve by \(O_m(R^\theta)\) grid points.
The number of curves must therefore grow at least as \(R^\theta\).
The claimed bounds grow only polynomially in \(k\) for these examples.
The counterexample with monic primitive polynomials also shows that
allowing only \(S\)-unit shifts does not repair the full-gcd estimate.
The outside-\(S\) example uses polynomial contents supported away from
\(S\), which the all-place gcd can detect. Separately, Lemma 5 asserts
a polynomial count in a logarithmic height cutoff, already false for
\(\mathbb Q\) with \(\theta=1\) in that lemma. The step from the height
cutoff to equation (48) uses this count. AP1–AP4 use the fixed-form
Subspace Theorem and are independent of these preprint assertions.

Qualitative fixed-\(S\) statements do not provide the uniform
curve number in AP2. The moving-target theorem of Grieve and Wang
assumes coefficient heights small relative to a sequence of point
heights and uses a fixed \(S\). The examples above have coefficient
and point logarithmic heights of comparable size. The explicit
specialization here concerns bounded-degree covers and counting for
moving prime pools with fixed shifts.
[Xiao, Theorem 3.6](https://arxiv.org/html/2110.01751)
[Grieve and Wang, Theorem 1.2](https://arxiv.org/html/1902.09109)
