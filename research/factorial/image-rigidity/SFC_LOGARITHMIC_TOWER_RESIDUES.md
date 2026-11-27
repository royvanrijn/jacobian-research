# Logarithmic residues for factorial moments at prime-power orders

Work over a finite extension of \(\mathbb Q_p\), with \(v(p)=1\).
Write the phase as \(F=c\,G\), with
\[
 G=X^a+P,\qquad P=\sum_{\nu\ne a}t_\nu X^\nu,\qquad
 \delta=\min_{\nu\ne a}
        \left(v(t_\nu)+\frac{|\nu|-|a|}{p-1}\right)>0 .
 \tag{1}
\]
Thus \(X^a\) is the unique least-slope term. The dimension and the
finite support are arbitrary. Put \(H_n=c^n\prod_i(a_i n)!\).
Let \(h\) be the largest number of coordinates in which a convex
combination of the other columns is strictly below \(a\), and set
\(m=\max(0,h-1)\). The monomial-only case has normalized moment one
and is omitted below.

For every fixed nonnegative integer \(r\) and every nonnegative
integer \(J\), there are coefficients \(B_0(r),\ldots,B_J(r)\)
such that, as \(\ell\) tends to infinity,
\[
 v\left((n-r)^m\frac{\Lambda(F^n)}{H_n}
             -\sum_{j=0}^J B_j(r)(n-r)^j\right)
       \ \ge\ (J+1)\ell-C_J ,
 \qquad n=r+p^\ell u_\ell .
 \tag{2}
\]
Here \(u_\ell\) is any positive integer prime to \(p\) with
\(\log_p(n+1)\le A\ell+B\) for fixed \(A,B\).
The constant in (2) may depend on this size bound. Each \(B_j(r)\)
is given by a convergent series over the perturbation counts and
admits an effective finite truncation at each fixed precision.
If some \(B_j(r)\ne0\), then the moments along these orders are
eventually nonzero. In particular, a single nonzero coefficient
among a fixed set of adjacent offsets certifies eventual
nonvanishing in the corresponding windows.

For \(h\ge2\) the leading coefficient across all offsets has a
common linear-functional description. For \(I\subseteq\{i:a_i>0\}\)
define \(\mathcal L_I\) on Laurent monomials \(X^\alpha\) by
\[
 \mathcal L_I(X^\alpha)=
 \prod_{i\in I}\frac{(-1)^{-\alpha_i-1}}{(-\alpha_i-1)!}
 \prod_{i\notin I}\alpha_i!
\]
if \(\alpha_i<0\) exactly on \(I\), and zero otherwise.
Then
\[
 B_0(r)=\frac{1}{\prod_i(a_i r)!}
 \sum_{\substack{I\subseteq\{i:a_i>0\}\\ |I|=h}}
 \frac{\mathcal L_I\bigl(G^r\log(1+P/X^a)\bigr)}
      {\prod_{i\in I}a_i}.
 \tag{3}
\]
The logarithm is expanded by degree in the perturbation coefficients,
and the functional on it means the convergent sum in this order.
For two variables with \(h=2\), the single functional is the residue
of \(e^{-x-y}G^r\log(1+P/X^a)\,dx\,dy\).
This is a weighted logarithmic functional. The theorem on unweighted
constant terms of Laurent powers of
[Duistermaat and van der Kallen, Sections 1 and 3](https://webspace.science.uu.nl/~kalle101/powers.pdf)
does not give the consecutive-window bound for (3).

## 1. Count expansion and the pole bound

For perturbation counts \(k=(k_\nu)\), let \(K=\sum k_\nu\) and
\(\Delta_i=\sum_\nu k_\nu(\nu_i-a_i)\). At a usual integer order
\(n\ge K\), the contribution to the normalized moment is
\[
 R_k(n)=\frac{(n)_{\underline K}}{\prod_\nu k_\nu!}
             \prod_\nu t_\nu^{k_\nu}
             \prod_i \frac{(a_i n+\Delta_i)!}{(a_i n)!}.
 \tag{4}
\]
Interpret each factorial ratio as a finite product or reciprocal
product of linear factors, so (4) is a rational function of \(n\).
The exact moment uses only \(K\le n\). Terms with \(K>n\) are
used below only in the Laurent expansions at \(r\).

At \(n=r\) a pole factor in row \(i\) occurs exactly when
\(a_i r+\Delta_i<0\). That row has \(a_i>0\) and its pole is simple.
If any row has a pole then \(K>r\), since
\(a_i r+\Delta_i=a_i(r-K)+\sum_\nu k_\nu\nu_i\).
The falling factorial supplies one zero. There are at most \(h\)
pole rows and hence \(z^mR_k(r+z)\) is regular at \(z=0\) for every
\(k\). If \(h\ge2\), a nonzero constant term in that series
requires exactly \(h\) pole rows.

## 2. Uniform control of the count tail

First let \(n\ge K\). Legendre's digit formula applied to (4) gives
\[
                       v(R_k(n))\ge \delta K-C(1+\log_p(n+1))
 \tag{5}
\]
for a constant depending only on the support and prime. Indeed the
falling factorial has linear valuation part \(K/(p-1)\), which
cancels the linear part of \(\prod k_\nu!\). The linear contribution
of each row ratio is \(\Delta_i/(p-1)\). All remaining factorial digit
sums in (4) have size \(O(1+\log_p(n+1))\), since \(T_i=a_i n+\Delta_i\)
and the counts are bounded by a fixed multiple of \(n\).
This uses the same termwise valuation as the
[relative valuation automata](FACTORIAL_RELATIVE_VALUATION_AUTOMATA.md),
equation (1).

For the expansions at \(r\), the nonzero integer offsets of the
linear factors have absolute value at most \(D_r(K+1)\) for a fixed
\(D_r\). If an integer interval has length \(q\) and its endpoints
are in that bound, the sum of valuations of its nonzero entries
is \(q/(p-1)+O(\log(K+2))\). To see this, count its multiples of
\(p^j\) with an error at most two per \(j\), including the removal
of zero. Beyond the logarithm of the largest absolute entry there
are no more terms, and truncating the sum of \(q/p^j\) costs only
\(O(1)\). Apply this interval estimate also to the falling
factorial. The zero factors themselves contribute only their fixed
integer slopes after removing their powers of \(z\). Thus the first
coefficient in \(z^mR_k(r+z)\) has valuation at least
\(\delta K-C_r\log(K+2)\).

Upon removing that first power and coefficient, the factors have
form \(1+\alpha z/b\) or its inverse, where \(b\ne0\) is one of the
bounded integer offsets. A coefficient of degree \(j\) in their
product loses at most \(j\log_p(D_r(K+1))\). The inverse factors
use geometric expansions. Consequently every fixed-degree
coefficient of \(z^mR_k(r+z)\) has valuation at least
\[
                         \delta K-C_{r,j}\log(K+2).
 \tag{6}
\]
The coefficient sum over \(k\) converges, defining \(B_j(r)\).
The same estimates yield a cutoff at each fixed precision because
the linear term in \(K\) dominates the logarithmic error.

For (2), choose a sufficiently large fixed \(M_J\).
When \(K>M_J\ell\), (5) and the size bound on \(n\) discard all
actual terms to precision \((J+1)\ell-O(1)\). The same holds for
the Taylor coefficients of degrees at most \(J\) by (6).
For \(K\le M_J\ell\) all nonzero offsets have valuation less than
\(\ell\) for sufficiently large \(\ell\), so the geometric expansions
converge at \(z=n-r\). The Taylor remainder of degrees at least
\(J+1\) has valuation at least
\((J+1)\ell+\delta K-C'_{r,J}\log(K+2)\).
The second and third terms have a fixed lower bound. These
\(K\) are all at most \(n\) for large \(\ell\). This proves (2),
using the nonarchimedean inequality in both tails.
If the first nonzero coefficient is \(B_j(r)\), use (2) with
\(J=j\). Its term eventually has strictly lower valuation than
the remainder.

## 3. The same residue functional for adjacent offsets

At a pole row put \(\alpha_i=a_i r+\Delta_i=-E_i<0\).
After stripping its factor \(z^{-1}\), the row ratio at \(z=0\)
is \((-1)^{E_i-1}/(a_i(E_i-1)!(a_i r)!)\).
This follows directly by multiplying the offsets from
\(-E_i+1\) through \(a_i r\), omitting zero.
At any other row the value is \(\alpha_i!/(a_i r)!\).
The numerator zero contributes the derivative of
\(\binom nK\) at \(r\), times \(K!/\prod k_\nu!\).
For \(K>r\) the derivative is
\((-1)^{K-r-1}r!(K-r-1)!/K!\).
In the formal binomial expansion,
\[
 \left.\partial_n(1+Z)^n\right|_{n=r}=(1+Z)^r\log(1+Z),
 \qquad Z=P/X^a .
\]
A term of degree \(K\) in \(Z\) multiplied by \(X^{ar}\) has
exponent \(ar+\Delta\) after the multinomial expansion. If \(K\le r\)
it has no negative coordinate and gives zero to every contributing
\(\mathcal L_I\). Combining the formulas proves (3) coefficient
by coefficient. Convergence follows from (6).

Thus the search for a consecutive block of vanishing leading
coefficients is a single weighted Laurent moment problem for
\(G^r\), with a common logarithmic factor. It is stronger to
control that common functional on a window than to check the
individual valuation minima separately. The reduction keeps the
actual signs and units even when many counts have the same valuation.

The common functional extends to every Taylor coefficient in (2).
Use \(\phi_b(T)\) for the finite factorial ratio \((T+b)!/T!\),
and set
\[
 \mathcal K_z(X^\alpha)=\prod_i\phi_{\alpha_i}(a_i z),
 \qquad Q_r(z)=\prod_i\phi_{a_i r}(a_i z),\qquad
 L=\log(1+P/X^a).
\]
The functional is defined only on exponents occurring in these
expansions. If \(a_i=0\), those exponents are nonnegative in
coordinate \(i\), so \(\phi_{\alpha_i}(0)\) is defined.
Since \(Q_r(0)=\prod_i(a_i r)!\ne0\), one has the formal identity
\[
 \sum_{j\ge0} B_j(r)z^j
       = Q_r(z)^{-1}\,z^m\mathcal K_z\!\left(G^r\exp(zL)\right).
 \tag{7}
\]
First read the right side by degree in perturbation coefficients,
then sum each Laurent coefficient in \(z\) by (6). Indeed
\(\phi_{\Delta_i}(a_i(r+z))
 =\phi_{a_i r+\Delta_i}(a_i z)/\phi_{a_i r}(a_i z)\).
For each count, the binomial factor in
\(G^r\exp(zL)=X^{ar}(1+P/X^a)^{r+z}\) delivers exactly (4).

For a given \(z\)-degree \(j\) in (7), use only the powers
\(L^k\) with \(k\le j+h-m\). Every Laurent monomial arising in
these powers times \(G^r\) can have poles in at most \(h\)
coordinates. This follows from the convex condition on
\(\Delta\), and gives the asserted bound on \(k\).
In the two-variable case with \(h=2\), for example, the next
coefficient after (3) uses only \(G^r\), \(G^r L\), and
\(G^r L^2\). Thus a vanished leading layer leaves another common
functional calculation on the same adjacent powers.

Explicitly in the two-variable case with \(h=2\), let \(K_b\)
be the coefficient of \(z^b\) in the monomial functional
\(\mathcal K_z\) and write \(H_j=\sum_{\ell=1}^j1/\ell\)
with \(H_0=0\). Then (7) gives
\[
 \begin{aligned}
 Q_r(0)B_1(r)
  ={}&\Lambda(G^r)+K_{-1}(G^rL)
             +\tfrac12 K_{-2}(G^rL^2)\\
     &-\left(\sum_i a_iH_{a_i r}\right)K_{-2}(G^rL).
 \end{aligned}
\]
For \(\alpha_i\ge0\), the row factor begins
\(\alpha_i!(1+a_i H_{\alpha_i}z)\). For \(\alpha_i=-k<0\)
it begins
\((-1)^{k-1}(1+a_i H_{k-1}z)/(a_i(k-1)!z)\).
These first two coefficients suffice for both displayed Laurent
functionals. If \(B_0(r)=0\), the final term drops out.
In that case the relative moment converges to \(B_1(r)\) along
the specified approaches to \(r\). In the two-row factored
example of Sections 5–6 under (10), that limit is one modulo
the maximal ideal because all the ordinary normalized moments
have the same unit residue.

## 4. Noninteger rational centers and a common offset functional

Let \(\rho\in\mathbb Q_{>0}\cap\mathbb Z_p\) be nonintegral.
In (2) allow positive integer sequences \(n_\ell\) with
\(v(n_\ell-\rho)=\ell\) and the same size bound. Use
\[
 I_\rho=\{i:a_i>0,\ a_i\rho\in\mathbb Z\},\qquad
 h_\rho=\max\#\{i\in I_\rho:b_i<a_i\},
\]
where the maximum is over convex combinations \(b\) of the
other columns. For this center the power of \(n_\ell-\rho\)
in (2) is \(m=h_\rho\). There is no falling-factorial zero
to subtract from the number of pole rows. In particular a
single descent can now give a genuine pole.
For \(a_i\rho\notin\mathbb Z\), a row cannot have a pole.
For \(a_i\rho\in\mathbb Z\), the pole criterion is still
\(a_i\rho+\Delta_i<0\), with at most one pole factor per row.

The estimates of Section 2 hold after this change.
Write \(\rho=b/q\) with \(p\nmid q\). Numerators of the
offsets in a row form consecutive terms of an arithmetic
progression of step \(q\). At precision \(p^j\) that step
is a unit, so the count of multiples differs from its
length divided by \(p^j\) by at most one. Removal of a
possible zero accounts for one further error. The numerators
remain bounded by a fixed multiple of \(K+1\), and dividing
by \(q\) leaves their valuations unchanged. This proves
the same interval bound and the same Taylor estimates.

There is again a common functional across adjacent
centers. Define \(\phi_b(T)\) for integral \(b\) to be the
rational factorial ratio \((T+b)!/T!\) in the sense of (4).
For a Laurent monomial \(X^\alpha\) use the rational
function in \(z\)
\[
                  \prod_i\phi_{\alpha_i}(a_i(\rho+z)).
\]
Define \(\mathcal T_\rho(X^\alpha)\) to be the coefficient
of \(z^{-h_\rho}\) in its Laurent expansion at zero.
For the exponents appearing below the pole order is at most
\(h_\rho\). For each fixed \(j\ge0\) put
\(C_{\rho,j}=\prod_i\phi_{a_i j}(a_i\rho)\ne0\).
Then the leading coefficient at the center \(\rho+j\) is
\[
 B_0(\rho+j)
       = C_{\rho,j}^{-1}\mathcal T_\rho
                          \bigl(G^j(1+P/X^a)^\rho\bigr).
 \tag{8}
\]
Use the formal binomial series in the perturbation coefficients.
For its term of count \(k\), the exponent in \(G^j(1+P/X^a)^\rho\)
is \(aj+\Delta\) and its scalar coefficient uses
\(\binom{\rho+j}{K}\). The ratio for that term at \(n=\rho+j+z\)
is the ratio with shift \(aj+\Delta\) at \(\rho+z\), divided
by the ratio with shift \(aj\) at \(\rho+z\). The denominator
is analytic and nonzero at zero. Since the pole order is at
most \(h_\rho\), its value at zero is all that enters the
top coefficient. This gives (8) and the series converges
by (6). In particular, nonintegral limiting orders produce
a binomial factor in the common functional, while integer
centers with canceled multinomial factors produce the
logarithmic factor in (3).

## 5. Systematic cancellation of the highest pole

For \(h\ge2\), any Laurent monomial with a nonnegative exponent
in one of the coordinates of \(I\), \(|I|=h\), is killed by
\(\mathcal L_I\) even after multiplication by \(G^r\).
Hence \(B_0(r)\) vanishes for every \(r\) if the logarithm \(L\)
has no monomials negative in a full set of \(h\) coordinates.
This condition can hold even when individual moment terms have
poles of the full order \(h-1\).

In particular, suppose \(G=\prod_q G_q\) and the anchor
\(X^a\) is the product of chosen unit-coefficient monomial anchors
\(M_q\) of the factors. Assume each quotient \(G_q/M_q=1+Z_q\),
where \(Z_q\) has zero constant term in a common system of
perturbation parameters. Then \(L=\sum_q\log(1+Z_q)\) formally.
For each size-\(h\) set \(I\), if every \(Z_q\) is nonnegative
in at least one coordinate of \(I\), every monomial of the
factor logarithm is also nonnegative in that coordinate.
All leading coefficients (3) then vanish.

For a mixed bivariate example take
\[
                 G=(x^2y+u y^2)(xy^2+v x^2),
                 \qquad v_p(u),v_p(v)>1/(p-1).
\]
The anchor \(x^3y^3\) has a unique least slope. The term \(uvx^2y^2\)
is lower in both coordinates, so \(h=2\). Yet the logarithm is
\(\log(1+u y/x^2)+\log(1+v x/y^2)\). Each summand has no
negative exponent in one coordinate, so every \(B_0(r)=0\).
This rules out a universal sharp-window proof from the leading
pole coefficients alone, even for finite supports with strict
slope isolation and genuine two-row termwise poles.
The next Laurent coefficient is retained by (7).

## 6. Factorwise expansion can recover an all-order window

There is an all-order certificate for a class of these cancellations.
Write a phase as
\[
 F=c\prod_{j=1}^{b}
       \left(X^{a_j}+\sum_{\nu\ne a_j}t_{j,\nu}X^\nu\right),
 \qquad A=\sum_j a_j .
 \tag{9}
\]
No upper bound on exponents, dimension, support or factor count is
imposed. Suppose the following hold at a prime \(p\).

- In each factor no convex combination of its other columns is
  strictly below \(a_j\) in two coordinates.
- For each coordinate in which any factor can descend from its
  anchor, the corresponding total anchor exponent \(A_i\) is a
  power of \(p\), allowing \(A_i=1\).

Let \(D_j\) be the possible descending coordinates in factor \(j\).
Put \(D_j^*=\max_{i\in D_j}(A_i-1)/(p-1)\), taking zero for an
empty set. If
\[
                 v(t_{j,\nu})\ge D_j^*+\frac1{p-1}
                 \quad\text{for every perturbation coefficient},
 \tag{10}
\]
then \(\Lambda(F^n)/(c^n\prod_i(A_i n)!)\) is one modulo the
positive-valuation ideal at every order. In particular all the
ordinary moments are nonzero. This applies over every
characteristic-zero field whose coefficient ratios admit a
valuation with these bounds, without requiring completeness.

Expand each factor of \(F^n\) separately. For a group of counts
in factor \(j\), let \(K_j\) be the sum of its nonanchor counts
and \(\Delta_{j,i}\) its displacement from the pure anchor in
coordinate \(i\). The convex condition allows at most one
negative \(\Delta_{j,i}\) in any factor for the chosen counts.
If \(T_i=A_i n+\sum_j\Delta_{j,i}<A_i n\), let \(J_i\) be the
factors with \(\Delta_{j,i}<0\), and choose \(j(i)\in J_i\)
with largest \(K_{j(i)}\). The sets \(J_i\) are disjoint
for the given counts, even if different factors can descend
in the same coordinate for different counts. No factor is
assigned twice. Since exponents are nonnegative,
\[
 T_i\ge A_i n-\sum_{j\in J_i}a_{j,i}K_j
                         \ge A_i(n-K_{j(i)}).
\]
The factorial ratio in row \(i\) therefore has valuation at least
\[
 -\frac{A_i-1}{p-1}K_{j(i)}
       -v\left(\frac{n!}{(n-K_{j(i)})!}\right)
\]
by the prime-power scaling identity. The independently expanded
factor \(j(i)\) supplies this exact falling factorial in its
multinomial, so it cancels the corresponding loss. All unused
falling factorials and nondecreasing row ratios are integral.
The remaining relative valuation is bounded below by
\[
 \sum_{j,\nu} k_{j,\nu}
   \left(v(t_{j,\nu})-D_j^*-\frac1{p-1}\right)
             +\frac{1}{p-1}\sum_{j,\nu}s_p(k_{j,\nu}).
\]
This is strictly positive if any perturbation is used. The
all-anchor term gives the stated unit congruence.

This combines the
[pure-anchor support test](FACTORIAL_RELATIVE_VALUATION_AUTOMATA.md#5-which-supports-allow-a-uniformly-dominant-pure-anchor)
with independent groupings. In particular the mixed example in
Section 5 has all moments nonzero at \(p=3\) when
\(v_3(u),v_3(v)\ge3/2\). The two total anchor exponents are
both \(3\), and each factor decreases in its own coordinate.
The expanded support would fail the uniform pure-anchor test
because the monomial \(uvx^2y^2\) is lower in both coordinates.
Its logarithmic leading residue also vanishes identically.
The independent falling factorials give the all-order result.

Overlapping possible descents are allowed. At \(p=2\), for example,
\[
      F=(x^2y^2+u y^2)(x^2y^2+v xy^3+w x^3y)
\]
has total anchor \((4,4)\) and both factors can descend in \(x\).
The second factor's two displacements sum to zero coordinatewise,
so no combination descends in two rows within that factor.
The support expansion includes terms strictly lower in both rows.
If the three coefficient valuations are at least four,
(10) gives all-order nonvanishing for this family as well.

## 7. When meromorphic interpolation fails

Suppose \(h\ge2\) and keep the support and anchor fixed. In any
sufficiently small closed coefficient polydisc about the origin,
each \(B_0(r)\) is a nonzero restricted power series in the \(t_\nu\).
To see it is nonzero, realize a rational convex combination
descending in \(h\) coordinates with integer counts. Scale those
counts until precisely those \(h\) rows cross below zero at \(r\).
The corresponding monomial in the independent \(t_\nu\) has a
nonzero coefficient from exactly one \(I\) in (3).
Estimate (6), with a positive slope margin throughout the
polydisc, gives restricted convergence. A nonzero restricted
series has a zero set with empty interior. In one variable this
follows by translating the center within the disk and comparing
the first nonzero Taylor coefficient at sufficiently small inputs.
Translation of restricted series on the unit disk is invertible.
In several variables apply the one-variable argument successively
on a product of subdisks. Thus the complement of the zeros of all
\(B_0(r)\), \(r\ge0\), is dense by the Baire category theorem on
the closed polydisc. One can also exclude the coordinate
hyperplanes to keep the specified support.

At any such coefficient vector, the normalized moments are
unbounded on orders \(r+p^\ell\) tending to each nonnegative
integer \(r\), by (2) with \(m=h-1>0\) and \(B_0(r)\ne0\).
Consequently the pure-normalized sequence has no meromorphic
interpolation on any open disk in \(\mathbb Z_p\), even if agreement
is required only at sufficiently large ordinary integers in the
disk. A meromorphic function on a small disk about any point
is a quotient of analytic functions with nonzero denominator
series. Its possible poles cannot accumulate at that point,
by the first nonzero denominator coefficient. The required
integer poles are dense. This does not concern interpolation
after a different normalization such as the
[\(n!\)-normalized series](FACTORIAL_LINEAR_PRECISION_ANALYTIC_INTERPOLATION.md),
which can lose the residue detected here. The dense coefficient
statement is over the local field and does not assert algebraicity.

There is an explicit rational example of the obstruction near
zero, as well as a contiguous window result, in the next section.

## 8. A test of the window transfer on a product plus a constant

Take \(G=X_1\cdots X_d+t\), \(d\ge2\), \(t\ne0\), with
\(v(t)>d/(p-1)\). Put \(\xi=(-1)^{d-1}t\) and, for integer \(s\),
\[
                 Z_s(\xi)=\sum_{j\ge0}\frac{\xi^j}{(j!)^{d-1}(s+j)!},
\]
where the reciprocal of a factorial at a negative integer is
taken as zero. Formula (3) gives
\[
                      B_0(r)=\frac{t^{r+1}}{(r!)^{d-1}}Z_{r+1}(\xi).
 \tag{11}
\]
Indeed, with \(K=r+j+1\) the \(d\) exponents at \(r\) are
all \(-j-1\). The factorial derivative in Section 3 and the
\(d\) reciprocal factorial factors give (7).

Use \(\theta=\xi\,d/d\xi\). Coefficient comparison gives
\[
 (s+\theta)Z_s=Z_{s-1},\qquad \theta^{d-1}Z_{s-1}=\xi Z_s .
\]
On rewriting each Euler operation with the first identity, the
second gives a recurrence with endpoint terms \(Z_{s-d}\) with
coefficient one and \(\xi Z_s\) on the other side. Its remaining
terms have indices between \(s-d+1\) and \(s-1\). A positive
starting index with \(d\) consecutive zeros would therefore
propagate backwards to \(Z_1=0\). But \(Z_1\ne0\): its
constant term has valuation zero and all other terms have
positive valuation, by \(v(t)>d/(p-1)\) and Legendre's formula.
Thus limit residues cannot have a zero \(d\)-window. By (2),
at each fixed nonnegative offset every sufficiently late
\(d\)-window of the specified prime-power orders contains a
nonzero moment. For \(d=2\) the length is the support size two.
This checks the transfer in a family whose Strong Factorial
conclusion in that dimension is already known at all orders.
See [the constant plus product theorem](CONSTANT_PLUS_PRODUCT_OF_TWO_VARIABLES_ALL_ORDERS.md).

Also for \(d=2\), put \(r=p^k\) with \(k\ge1\). For \(j\ge1\),
\[
 v(j!)+v((r+j+1)!/(r+1)!)\le 2j/(p-1)
\]
because \(s_p(r+1)=2\) and every positive integer has digit sum
at least one. Therefore \((r+1)!Z_{r+1}(\xi)\) has unit leading
term and all other terms have positive valuation. At every
\(r=p^k\) the pure-normalized moments are unbounded along
\(r+p^\ell\). These centers accumulate at zero. In particular
any rational \(t\ne0\) with \(v(t)>2/(p-1)\) gives a normalized
sequence with no meromorphic interpolation on a neighborhood
of zero.

The tower route to unrestricted SFC(2) still needs sharp
support-size separation for the full family (7) on arbitrary
supports and control of general starts and tied least slopes. For a bivariate phase with one least slope, the concrete
next test is whether at least one of the full series (7) is
nonzero in each window of \(s\) integer centers, where \(s\)
is the actual support size. The highest logarithmic layer alone
fails by Section 5. Section 6 recovers all orders for grouped
factors whose active descents can be assigned separately, leaving the analogous
higher-layer separation for general factors open.
