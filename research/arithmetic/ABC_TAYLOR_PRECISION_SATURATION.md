# Taylor precision up to the full modulus

The first residue determinant in (QD2) can retain its
Taylor gain after the row count exceeds the available
lifting order. The gain saturates instead of stopping
at half the backed full-modulus weight.

Fix a modulus \(N=\prod_p p^{e_p}\) on one letter and an
exact unit congruence \(v z^q\equiv K\bmod N\) for a
grouped ratio on the other two. At every modulus prime
choose an arbitrary shallow depth \(b_p\ge0\). A chart
records \(z\bmod p^{b_p}\) where \(b_p>0\). Depth beyond \(e_p\) gives the same valuation guarantee
in the capped estimate. Choices of depths, if adaptive, cost at most
divisor-many records for fixed \(N\).

Consider a square evaluation minor of size \(D_0\) on
any selected polynomial monomial columns. On a chart
at \(p\), replace \(v\) by \(K z^{-q}\); entries agree
modulo \(p^{e_p}\), and row denominator clearing is a
unit. Expand the exact monomials about a unit chart
center. The Taylor vector of order \(j\) has
\(p\)-valuation at least \(j b_p\). Keep the Taylor
orders \(j<e_p/b_p\), and collect every higher order
and the original correction into a vector divisible
by \(p^{e_p}\). A surviving multilinear term in the
determinant cannot use a Taylor vector twice, whereas
arbitrarily many residual columns are allowed.
Consequently
\[
 v_p(\det)\ge
       \sum_{j=0}^{D_0-1}\min(j b_p,e_p).
 \tag{ST1}
\]
For \(b_p=0\) the displayed contribution is zero.
This argument uses the individual \(e_p\)'s directly;
no common Taylor truncation bound is needed.

Write \(\gamma_p=\log_X(p)/L\) for the normalized mass
of a prime, \(r_p=e_pL\) for its scaled exponent,
\(b=\sum_p\gamma_p b_p\) for the shallow record cost,
and \(w=\sum_p\gamma_p r_p\) for the full lock weight.
In the limit \(D_0L\to y\), the per-row divisibility
from (ST1) is at least
\[
 \sum_p\gamma_p\Phi(b_p y,r_p)+o(1),\qquad
 \Phi(u,r)=\int_0^1\min(ut,r)\,dt
   =\begin{cases}
      u/2,&0\le u\le r,\\[2pt]
      r-r^2/(2u),&u\ge r>0.
     \end{cases}
 \tag{ST2}
\]
The left Riemann sum differs from the integral by at
most \(w/D_0\), after summing over primes. In particular
if \(r_p\ge\kappa b_p\) at all retained shallow places,
the gain is at least \(b\Phi(y,\kappa)\).
This is \(by/2\) for \(y\le\kappa\), and
\(b(\kappa-\kappa^2/(2y))\) for \(y\ge\kappa\).
The ordinary rank defect is forced whenever this exceeds
the per-row Archimedean bound. It can replace the first
Taylor inequalities in (DNR1), with the same Newton recovery
of supported root values after the auxiliary is produced.
The terminal lift in (QD2) has its own height and precision
requirements.
[First determinant and recovery bounds](ABC_QUANTITATIVE_DEGREE_WINDOWS.md#direct-newton-recovery-at-arbitrary-nonatomic-mass)

For an upper bound on this particular determinant criterion,
\(\Phi\) is jointly concave and positively homogeneous.
Hence for all adaptive depth assignments,
\[
 \sum_p\gamma_p\Phi(b_p y,r_p)\le\Phi(b y,w).
 \tag{ST3}
\]
This bounds the gain certified by (ST1), not additional
divisibility a determinant might have for other reasons.

## An exact new wholly diffuse sector

Take all three formal degree measures equal to
\[
  \tfrac13\operatorname{Unif}
       [\,3-117/1000,\,3+117/1000\,],
\]
with full numerical heights one. Use the entire free
pair at common group base \(a=3\) with multiplier one.
The narrow height cap is \(h=1/3\), and the wider cap
is \(E=39/2000\). Indeed the positive and negative
numerical deviations on a single letter are each
\(117/(12000)\), and the ratio adds one from each
letter on each side. Lock normalized mass
\(\mu=33/200\) of the supplying letter in a continuity
band strictly above three. The available mass above
three is \(1/6\). The bounded-degree roots are squarefree,
and depth one gives \(b=33/200\) and mass depth
\(\kappa=3\). Set
\[
 M=7,\qquad \alpha=99/100,\qquad r_0=35/66,\qquad
 y=\alpha(M+1)r_0=21/5.
\]
Then
\[
 B=h r_0+ME=\frac{35}{198}+\frac{273}{2000}
  < \frac{891}{2800}
  = b\left(\kappa-\frac{\kappa^2}{2y}\right).
\]
The wider polynomial degree is fixed, and the nonatomic
atomic charge in (DNR2) is zero. Thus this sector has a
strict counted bound \(d\le\mu+b=33/100<1/3\).
All margins are strict.

For comparison, even the first old Taylor determinant for
this same full free grouping cannot reach a saving with
any independent rectangular parameters. Indeed its
conditions and \(b\kappa\le R\mu\), where
\(R=3117/1000\), imply
\[
 E M < R\mu\left(\frac12-\frac{h}{b(M+1)}\right).
\]
Under \(\mu+b<d=h=1/3\) the largest right side for
fixed \(M\) is at most
\(\frac12 R d(1-\sqrt{2/(M+1)})^2\), by optimizing in
\(b\). Even a fixed amplification has the same upper
bound, because in this derivation the effective chart
cost becomes \(b/\ell\), while an extra factor
\(1/\ell\) on the final right side can be dropped.
Put \(x=\sqrt{2/(M+1)}\). Then
\[
 \frac{(1-\sqrt{2/(M+1)})^2}{2M}
   = \frac{x^2(1-x)^2}{2(2-x^2)}
   \le\frac{25}{1358}.
\]
For \(x\le11/20\) the numerator is at most \(1/16\)
and the denominator is at least \(679/200\).
For \(x>11/20\) the expression is decreasing, since
the derivative has the sign of \(x^3-4x+2<0\).
Consequently the old criterion would require
\(E<(3117/3000)(25/1358)<39/2000\).
The comparison fixes the full grouping and permits arbitrary
old rectangular depths, ranks, and amplification.

## The broad sparse one-root obstruction survives

Now take the common formal profile
\(\nu=(1/3)\operatorname{Unif}[2,4]\).
Allow arbitrary retained submeasures on the free letters,
possibly different nonzero signed integer multipliers, and
arbitrary sparse horizontal slices in a one-root
projective determinant. Also allow all adaptive
per-prime shallow depths on the chosen full modulus.

Let \(s\in(1/6,1/3]\) be half the total retained free
mass, and \(g=2s-1/3\).
A strict count improvement from the first determinant
and (DNR2) requires \(\mu+b<g\). Write \(J\) for the gap
between extreme wider slices, \(R\) for the largest
within-slice width, \(\mathcal A\) for the summed slice
length, and \(B\) for the row height.
Then
\[
 \mathcal A\le(J+1)R,\quad B\ge hR,\quad B\ge JE\ (J>0).
\]
Here \(h\ge s\) is the average retained integral of the
absolute multipliers. For a shifted base \(a'\) set by two
monomials in the extreme slices, \(E\) is the average
of their two free integrals of \(|r-a'i(r)|\).
These bounds hold for asymmetric assignments by viewing
two evaluation monomials as a primitive two-coordinate
subrow, whose height is at most that of the full row.
The logarithmic height of their quotient is at least
half the sum of absolute prime valuations.
[Sparse width bounds](ABC_QUANTITATIVE_DEGREE_WINDOWS.md#a-broad-all-diffuse-profile-beyond-the-one-root-criterion)

Write \(y=\alpha\mathcal A\) for the Taylor size and
\(t=by/w\). A successful saturated determinant requires
\(B<w\phi(t)\) by (ST3), where
\[
 \phi(t)=
 \begin{cases}t/2,&0<t\le1,\\1-1/(2t),&t\ge1.\end{cases}
\]
Also \(w\le4\mu\). It follows that
\[
 h < b(J+1)\phi(t)/t,\qquad
 E < \frac{4\mu}{J}\phi(t).
\]
The first bound uses \(R\ge\mathcal A/(J+1)\).
Because \(\phi(t)/t\le1/2\) and \(b<g\le s\le h\),
we must have \(J\ge2\) and \(b>2s/(J+1)\).
The two bounds give
\[
 Eh < 4\mu b(1+1/J)\frac{\phi(t)^2}{t}
     \le \frac{32}{81}g^2.                       \tag{ST4}
\]
Indeed \(\phi^2/t\le8/27\), by optimizing the elementary
one-variable function. For \(J\ge3\) use
\(\mu b\le g^2/4\). For \(J=2\) use
\(b>2g/3\), which gives \(\mu b\le2g^2/9\).

A second consequence uses both budgets before multiplying:
\[
 \frac{E}{s}
 < \frac4J\left(\phi(t)-\frac{t}{J+1}\right)
 \le \frac4J\left(1-\sqrt{\frac2{J+1}}\right)
 \le \frac25.                                      \tag{ST5}
\]
The middle inequality follows by maximizing in \(t\)
on the two pieces of \(\phi\). The last holds by direct
comparison for \(2\le J\le7\), and for \(J\ge8\) by
bounding the middle term by \(4(J-1)/(J(J+1))\).

For completeness these inequalities contradict the true
folded-distance bounds. Put \(u=3s\in(1/2,1]\), and
let \(\nu'\) be the average retained measure, of mass \(s\)
and dominated by \(\nu\).
For \(a'=0\) one has \(E\ge2s\), contradicting (ST5). Put \(a=|a'|\).
For \(a>0\), each integrand is at least the distance
from \(r\) to \(a\mathbb Z\), and
\[
 ah+E\ge\int r\,d\nu'(r)\ge(2+u)s.
\]
Therefore (ST5) gives \(ah>(8/5+u)s\).
If \(0<a\le2\), folded phase density is at most
\((2+a)/3\), so
\(E\ge 3as^2/[2(2+a)]\).
This yields
\[
 Eh\ge \frac38(8/5+u)s^3\quad(0<a\le2).
\]
The ratio \(g^2/(3s^3)=(2u-1)^2/u^3\)
increases on \([1/2,1]\).
For \(u\le3/4\) it is at most \(16/27\), and for
\(3/4\le u\le4/5\) it is at most \(45/64\).
The last lower bound on \(Eh/s^3\) is respectively
at least \(63/80\) and \(141/160\), both strictly
above the respective right sides of (ST4).
For \(u\ge4/5\) and \(a\le1\) the same density bound
improves its \(3/8\) to \(1/2\), giving
\(Eh\ge(6/5)s^3\).

It remains to check \(u\ge4/5\) and \(1<a\le2\).
Split the interval \([2,4]\) into one full period
of length \(a\) and the remaining interval of length
\(2-a\). The folded numerical-distance density is
at most \(1/3\) on each separately.
For a fixed split of the retained masses \(s_1,s_2\),
this gives \(E\ge(3/2)(s_1^2+s_2^2)\).
Since the latter interval contains mass at most
\((2-a)/6\), optimization over the split gives
\[
 E\ge\frac32s^2((1-z)^2+z^2),\qquad
 z=\min\left(\frac12,\frac{2-a}{2u}\right).
\]
If \(z=1/2\), then \(a\le 2-u\le6/5\), so
\(((1-z)^2+z^2)/a\ge5/12\).
Otherwise \(a=2-2uz\le2-(8/5)z\), and
\[
 (1-z)^2+z^2\ge \frac{2-(8/5)z}{3}
\]
because the difference is
\(1/3-(22/15)z+2z^2>0\).
In particular \(Eh\ge(6/5)s^3\) here too.
For \(a\ge2\) the whole interval occupies at most
one period, and folded distance density at most
\(1/3\) gives \(Eh\ge(3/2)s^3\) from \(h\ge s\).
Finally \(g^2\le3s^3\), so these bounds contradict
(ST4), whose right side is at most \((32/27)s^3\).
Thus the broad profile still evades this saturated,
sparse one-root determinant even with different free
assignments and adaptive per-prime depth choices.

## Anisotropic shallow records on rectangular multigroup supports

For \(m\ge2\) narrow groups, let \(b_{i,p}\) be the
shallow record depth of the \(i\)-th coordinate at prime
\(p\). Use the primes at which all these depths are
positive. The corresponding determinant calculation
takes each Taylor multiindex \(\mathbf j\) with local
weight \(\sum_i b_{i,p}j_i\), capped at \(e_p\) as in
(ST1). For a minor of size \(N\) it certifies the sum
of the \(N\) smallest capped weights, counting all
multiindices and using arbitrary residual columns
at weight \(e_p\).

Fix \(m\). The weight of the \(k\)-th cheapest
uncapped index is at most
\((m! k\prod_i b_{i,p})^{1/m}\). Indeed the simplex up
to that weight has volume \(k\), and taking coordinate
floors covers it by unit boxes whose integer bases lie
in the same simplex. Hence, with \(NL^m\to\alpha V\)
and \(y=(m!\alpha V)^{1/m}\), the limiting per-row
gain certified by the selected Taylor weights is at most
\[
 \sum_p\gamma_p \Phi_m(y c_p,r_p)
       \le \Phi_m(y b,w),\qquad
 \Phi_m(u,r)=\int_0^1\min(u t^{1/m},r)\,dt .
 \tag{ST6}
\]
Here \(c_p=(\prod_i b_{i,p})^{1/m}\),
\(B_i=\sum_p\gamma_p b_{i,p}\), \(b=(\prod_i B_i)^{1/m}\),
and \(w=\sum_p\gamma_p r_p\), with \(\gamma_p,r_p\) as in
(ST2). The Riemann-sum error after summing over primes
is at most \(w/N\), and the second inequality uses
homogeneity and joint concavity of \(\Phi_m\) and
Hölder on \(\sum_p\gamma_p c_p\). In the uniform
depth case, if all \(r_p\ge\kappa b_p\), the gain
given by the same calculation is at least
\(b\Phi_m(y,\kappa)+o(1)\), using the asymptotic count
of multiindices of each total order. Explicitly,
\[
 \Phi_m(y,\kappa)=
 \begin{cases}
   m y/(m+1), & y\le\kappa,\\
   \kappa-\kappa^{m+1}/((m+1)y^m), & y\ge\kappa .
 \end{cases}
\]
This is the extension of (TN4) that uses higher rank
than the first uncapped Taylor range.

Even this criterion has a broad-profile obstruction.
Give all three letters mass \(1/3\) uniformly on \([2,4]\),
and use the same rectangular hypothesis as in (BG1).
Allow arbitrary retained submeasures, bounded signed
multipliers and distinct depth records for each group.
Write \(s_i\) for the average retained mass of group \(i\)
across the two free letters, \(S=\sum_i s_i\), \(g=2S-1/3\),
and \(\mathfrak b=\sum_i B_i\). Suppose each entire group
dictionary costs at least its own shallow chart, \(2s_i\ge B_i\).
If \(\mu\) is the full lock cost, a strict saving would require
\(\mu+\mathfrak b<g\). This implies \(1/6<S\le1/3\) and
\(g\le S\). All \(B_i\) may be taken positive, since otherwise
(ST6) certifies no gain.

Let \(h_i\) be the narrow height caps and \(E\) the wider
height cap. The folding estimate (BG1) applies. In addition,
since \(s_i\ge B_i/2\), one ratio \(2s_i/B_i\) is at least
\(2S/\mathfrak b\), so
\[
 \prod_i h_i\ge\prod_i s_i
 \ge \frac{\prod_i B_i}{2^m}\frac{2S}{\mathfrak b},
 \qquad
 E\prod_i h_i\ge
    \frac{3}{4m2^m}\frac{\prod_i B_i\,S^3}{\mathfrak b}.
 \tag{ST7}
\]
The second bound follows from \(\sum_i s_i^2\ge S^2/m\).
The rectangular row height is \(B=\sum_i h_i r_i+ME\),
with \(V=(M+1)\prod_i r_i\). Put \(\kappa=w/b\) and
\[
 x=\frac{m+1}{b}
       \left(\frac{\prod_i h_i}{m!\alpha(M+1)}\right)^{1/m}.
\]
By the arithmetic-geometric mean inequality,
\(B\ge bmxy/(m+1)+ME\), so (ST6) forces \(x<1\).
For fixed \(x\), the largest possible
\(\Phi_m(y,\kappa)-mxy/(m+1)\) over \(y>0\) is
\(\kappa(1-x^{m/(m+1)})\), attained at
\(y=\kappa x^{-1/(m+1)}\). Using \(w\le4\mu\) gives
\(EM<4\mu(1-x^{m/(m+1)})\). By the definition of \(x\)
and (ST7),
\[
 \frac1M\le
 \frac{m!\,b^m x^m}{(m+1)^m\prod_i h_i}
   \frac{1}{1-u_m\mathfrak b/(2S)},\qquad
 u_m=\frac{2^m m!}{(m+1)^m}.
\]
In particular \(M\ge1\), since \(1/(M+1)\le
u_m\mathfrak b/(2S)<1\). Optimize
\(x^m(1-x^{m/(m+1)})\) by putting \(z=x^{m/(m+1)}\);
its maximum is \((m+1)^{m+1}/(m+2)^{m+2}\).
Consequently
\[
 E\prod_i h_i<
  \frac{4m!(m+1)}{(m+2)^{m+2}}\,
  \frac{\mu\prod_i B_i}{1-u_m\mathfrak b/(2S)}.
 \tag{ST8}
\]
For \(m\ge2\), \(u_m\le8/9\) by the successive ratio
\(2((m+1)/(m+2))^{m+1}<1\). Write \(v=\mathfrak b/g\).
Using \(\mu+\mathfrak b<g\), \(g\le S\) and
\(g^2\le3S^3\) gives
\[
 \frac{\mu\mathfrak b}{1-u_m\mathfrak b/(2S)}
 \le g^2\,\frac{v(1-v)}{1-4v/9}
 \le\frac{g^2}{3}\le S^3 .
\]
The middle inequality follows from
\(v^2-31v/27+1/3>0\), whose discriminant is \(-11/729\).
Combining (ST7)–(ST8) would require
\[
 1<
 \frac{16m\,2^m m!(m+1)}{3(m+2)^{m+2}}
 \le 1 .
\]
The coefficient on the right equals one at \(m=2\),
decreases at \(m=3\) by direct comparison, and thereafter
has successive ratio
\(2(1+1/m)((m+2)/(m+3))^{m+3}<1\).
This proves the claimed obstruction. Fixed amplification
does not evade it: divide the amplified heights and
chart costs by its order and use the resulting full
weight cap \(w/\ell\le4\mu\). For the same reason the
sparse one-root obstruction above holds under fixed
amplification.

The conclusion concerns the first determinant forced
by capped Taylor weights. In particular, the broad profile
still requires at least one group whose whole dictionary
costs less than its own chart, even if shallow depths
are assigned separately in all coordinates and at all
primes. The count on the surviving group configurations
requires arithmetic information beyond that Taylor bound.
[Original rectangular setup and scalar folding bound](ABC_QUANTITATIVE_DEGREE_WINDOWS.md#a-broad-diffuse-obstruction-for-rectangular-multigroup-determinants)


## Expanding the exact congruence error as a second variable

The original congruence gives more than a cap if the
auxiliary support has bounded degree \(M\) in the wider
variable. At a retained prime write
\[
 z=z_0+p^{b_p}t,\qquad
 v=Kz^{-q}+p^{e_p}u,\qquad t,u\in\mathbb Z_p .
\]
The choice of shallow record is the same as before.
The substituted monomial evaluations are power series in
\(t\) with powers \(u^j\), \(0\le j\le M\).
The vector of index \((i,j)\) has divisibility
\(p^{b_pi+e_pj}\). A multilinear expansion of an evaluation
minor of size \(N\) uses distinct pairs, so it gives the
sum of the \(N\) smallest such weights. It now uses the
structure of the error above the congruence modulus.
No record of \(u\) is needed.

For uniform depth ratio \(\kappa=L e_p/b_p\) and
\(y=LN\), the normalized mean order is
\[
 \Psi_M(y,\kappa)=\frac{1}{2y}
      \sum_{0\le j\le M}(\tau^2-\kappa^2 j^2)
                   \mathbf1_{\kappa j<\tau},\qquad
 y=\sum_{j=0}^M(\tau-\kappa j)_+ .
 \tag{ST9}
\]
This follows by counting the strips of indices with
\(Li+\kappa j\le\tau\) and summing their weights.
For a lower bound it suffices to have \(Le_p/b_p\ge\kappa\).
For arbitrary shallow depth assignments, the sum of the
\(N\) smallest local weights is an infimum of linear
functions of \((b_p,e_p)\), so it is jointly concave
and homogeneous. Aggregate as in (ST3). The criterion
certifies at most \(b\Psi_M(y,w/b)\) for total shallow
cost \(b\) and full numerical weight \(w\).
A prime without shallow information can be given zero
gain in this comparison.

### A sector beyond the capped Taylor criterion

Use three equal mass-\(1/3\) measures uniformly distributed
on \([3-3/20,3+3/20]\). The full free grouping at scale
three gives \(h=1/3\) and \(E=1/40\).
On the supplying letter lock mass \(\mu=b=33/200\)
above three, using depth one, so \(\kappa=3\) is valid.
Take \(\alpha=99/100,\ M=7,\ r_0=25/33\), giving \(y=6\).
Only the first two \(j\)-strips are active in (ST9).
In this range \(\Psi_M(y,\kappa)=y/4+\kappa/2-
\kappa^2/(4y)\). Thus
\[
 B=\frac{25}{99}+\frac7{40}
   <\frac{693}{1600}=b\Psi_M(6,3).
\]
The first determinant and (DNR2) count this sector with
normalized exponent at most \(\mu+b=33/100\). Strict
continuity margins permit nearby profiles.

For comparison, the capped criterion on this full rectangle
would require
\[
 E< \frac{w b x^2(1-x)}{2h-bx^2},\qquad
 x^2=\frac{2h}{\alpha b(M+1)}, \quad 0<x<1.
\]
This follows by optimizing \(b\Phi(y,\kappa)-h y/
(\alpha(M+1))\) in \(y\), with \(w=b\kappa\), and substituting
for \(M\). On the displayed profile \(w\le(63/20)\mu\) and
\(\mu+b<d=1/3=h\). If \(b\le2d/3\), use
\(\mu b\le d^2/4\) and
\(x^2(1-x)/(2-(2/3)x^2)\le32/351\).
Otherwise use \(\mu b\le2d^2/9\) and
\(x^2(1-x)/(2-x^2)\le64/621\).
For \(x\le3/4\) both rational bounds follow from the
numerator bound \(4/27\) and the minimum denominator
on that range. For \(x\ge3/4\) the fractions decrease,
since the derivative of \(x^2(1-x)/(2-a x^2)\) has
the sign of \(4-6x+a x^3<0\) for \(0<a\le1\).
The upper bounds on \(E\) are respectively
\((21/20)(32/351)/4\) and \((21/20)(2/9)(64/621)\),
both below \(1/40\). The same upper bounds allow
fixed amplification by dropping its helpful full
weight loss.

### A broad-profile barrier for the exact error on rectangles

Consider the broad \([2,4]\) profile with one rectangular
group, arbitrary retained submeasures and signed
multipliers as above. Use \(s,g,h,E\) as in (ST4).
For a putative saving, \(\mu+b<g\le s\le h\), and
\(w\le4\mu\). Set \(\kappa=w/b\) and
\(a=h/(b\alpha(M+1))\). The aggregate error criterion
requires
\[
 ME < b\bigl(\Psi_M(y,\kappa)-a y\bigr).
 \tag{ST10}
\]
The \(j=0\) strip alone shows \(\Psi_M\le y/2\), so
\(a<1/2\) is necessary. In particular \(M\ge2\).

If exactly \(k\) strips are active in (ST9), solving for
\(\tau\) and summing their squared offsets give
\[
 \Psi_M(y,\kappa)=
   \frac{y}{2k}+\frac{\kappa(k-1)}2
      -\frac{\kappa^2 k(k^2-1)}{24y}.
\]
When the last strip becomes active,
\(y_0=\kappa M(M+1)/2\), the derivative equals
\((2M+1)/(3M(M+1))\le1/(M+1)\). Thereafter it is
decreasing. Since \(a>1/(M+1)\), the supremum in
(ST10) occurs with \(y\le y_0\), in the range where
allowing all \(j\ge0\) leaves (ST9) unchanged.
In that range the strip count satisfies
\(\sum_{j\ge0}(\tau-\kappa j)_+\ge\tau^2/(2\kappa)\)
by the integral comparison. Its inverse mean order
is at most \((2/3)\sqrt{2\kappa y}\), by integrating
the corresponding inverse bound up to \(y\).
Optimizing this square root minus \(a y\) gives
\(ME<2b\kappa/(9a)\). Hence
\[
 Eh<\frac29 b w(1+1/M)\le\frac13g^2,\qquad E<s/3 .
\]
Here \(w b\le4\mu b\le g^2\), \(M\ge2\), and \(g\le s\le h\).
The folded distance and retained moment inequalities
following (DNR3) contradict these two bounds. In applying
them to a rectangle \(E\) refers to its actual grouping
base, which is the unshifted instance of those inequalities.
Thus even the exact error expansion on a single rooted
rectangle cannot produce a saving on this broad law.
The mixed area argument below extends the obstruction
to sparse supports.


### A mixed area barrier for sparse exact-error supports

Keep the broad uniform \([2,4]\) profile, and now use any
sparse one-root projective support in the nonnegative
horizontal degree and the integer wider degree, with
bounded nonzero signed multipliers on the free letters. They
need not agree on the two sides. Let \(P\) be the convex
hull of the scaled support, \(\mathcal A\) the sum of its
prescribed horizontal slice lengths before taking the hull,
and \(J\) the gap of its extreme wider slices. Factoring
out the common \(v\)-power at the bottom of the support
does not change projective row heights or the \(p\)-adic
minors at the unit places. The exact-error expansion therefore
has the finite strip cap \(M=J\).

We first bound the projective row height from below.
On each free side the direction \(n(r,i)=(i,r-a i)\)
is the valuation direction of a root at scaled degree
\(r\) assigned multiplier \(i\), with the chosen base \(a\).
The normalized finite-place valuation vectors for the two
sides have opposite signs. The Archimedean vector closes
their sum, by the product formula. The projective row
height is consequently the sum of the support function
of \(P\) on these weighted vectors, including the closing
one. For finite vectors, rotate the outward normal vectors
by a right angle in the positive boundary orientation
and concatenate them in angle order. Their sum is zero.
They are the successive edges of a convex polygon \(K\).
The height functional is \(2V(P,K)\), with \(V\) the mixed
area, so it is at least \(2\sqrt{\operatorname{area}(P)
\operatorname{area}(K)}\) by the mixed area inequality.
The formulas pass to continuous marked degree measures
by approximation of their bounded supports.
[Schneider, Mixed volumes and mixed area measures](https://doi.org/10.1017/CBO9781139003858.008)
· [Minkowski inequality, §7.2](https://doi.org/10.1017/CBO9781139003858.010)

The finite-prime vectors of each side lie in one open
semicircle. Indeed the linear functional \((x,y)\mapsto ax+y\)
is \(r>0\) on \(n(r,i)\) for either sign of \(i\); the two
opposite signs of the letters give opposite semicircles.
The two groups of consecutive edges delimit two caps on
opposite sides of their joining chord. Include the closing
edge in whichever semicircle contains it. The cap area
for each group is half the sum of the nonnegative pairwise
edge determinants within that group. In particular it
is at least half the area of its directional zonotope.
To bound that area, fix \(r_0,i_0\) on the side. For
\(0\le v\le1\), the condition
\[
        |i_0 r-i r_0|\le v
\]
places \(r\in[2,4]\) in intervals of length \(2v/|i_0|\)
with spacing \(r_0/|i_0|\). These intervals have no interior
overlap, since \(2v\le r_0\). An interval of length two
meets at most \(|i_0|\) periods in the sense of total
covered length, since \(r_0\ge2\), so the total covered
length is at most \(2v\). Thus a retained marked measure
of mass \(s_W\) on this side has at most \(v/3\) mass in
this condition. Integrating the complementary mass from
zero to \(3s_W\le1\) shows that the determinant integral
against this fixed vector is at least \(3s_W^2/2\).
The area of the side's zonotope is half the two-vector
integral, hence at least \(3s_W^3/4\). Summing the two caps
gives \(\operatorname{area}(K)\ge3s^3/4\), where \(s\) is
the average retained mass across the two free letters.
The same lower bound holds if the closing vector lies on
a boundary of the semicircles, by taking limits.

For \(J>0\), let \(I\) be the segment obtained by scaling
a segment between extreme-slice points by \(1/J\).
It has vertical displacement one. Sweeping each integer
slice through \(I\) shows
\(\operatorname{area}(P+I)\ge\mathcal A\).
The height functional is translation invariant and
Minkowski additive. As \(B(I)\le B(P)/J\), it follows that
\[
 B(P)\ge \frac{J}{J+1}\sqrt{3s^3\mathcal A}.
 \tag{ST11}
\]
This thickening accounts for the discrete slices.

Let \(g=2s-1/3\), \(v=b/g\), \(t=by/w\), \(k=J+1\),
and \(\psi=\Psi_J(t,1)\). A successful determinant has
\(B(P)<w\psi\). The largest within-slice horizontal width
gives \(B(P)\ge h\mathcal A/k\), with \(h\ge s\);
also \(y\le\mathcal A\). Since \(\psi\le t/2\) and
\(b<g\le s\), this forces \(J\ge2\). Using
\(wb\le4g^2 v(1-v)\) and \(g^2\le3s^3\), (ST11) and
the width inequality would require
\[
       4(k/J)^2v(1-v)\,\psi^2/t>1,\qquad
       \psi/t>1/(kv).
 \tag{ST12}
\]
We check both regimes of the finite strip count.

If \(t\le Jk/2\), allow all nonnegative strips without
changing \(\Psi_J(t,1)\). The strip count satisfies
\(\sum_{j\ge0}(\tau-j)_+\ge(\tau^2+\tau)/2\).
Integrating the inverse upper bound
\(\sqrt{2t+1/4}-1/2\) yields, for
\(z=\sqrt{2t+1/4}-1/2\),
\[
    t=z(z+1)/2,\qquad
    \psi\le z\left(\frac23-\frac1{6(z+1)}\right).
\]
The second condition in (ST12) implies \(z+1<4k/3\),
so \(\psi\le C z\) with \(C=2/3-1/(8k)\).
Consequently \(\psi^2/t\le2C^2-C\psi/t
<2C^2-C/(kv)\). The left side of (ST12) is then at most
\(2(C+(C-1/2)/J)^2=2(2/3+1/(24J))^2\le121/128<1\),
by maximizing the resulting quadratic in \(v\).

If \(t\ge Jk/2\), put \(\zeta=2t/(Jk)\ge1\) and
\(u=k\psi/t=1/2+1/\zeta-c/\zeta^2\), where
\(c=(J+2)/(6J)\le1/3\). The second condition gives
\(v>1/u\) with \(1<u\le3/2\), so
\(u^2v(1-v)\le u-1\). The left side of (ST12) is
at most \(2(k/J)\zeta(u-1)\). Since
\(\zeta(u-1)=1-\zeta/2-c/\zeta\) decreases for
\(\zeta\ge1\), this bound is at most
\(2(J^2-1)/(3J^2)<1\). This also contradicts (ST12).

The exact-error Taylor certificate with the primitive
maximum row bound therefore fails on this broad formal
profile for arbitrary nonzero signed
one-root assignments and sparse slice shapes, including
adaptive per-prime depths and fixed amplification. In the
amplified case the same argument uses the effective shallow
cost and full weight divided by the amplification order.


## Archimedean column weights on an anatomic box

At fixed box scales and multiplier assignments, the
primitive row estimates can also retain the Archimedean
weight of each selected column. In the same notation
as the mixed area argument, let \(P\) be the convex
hull of the full monomial support and let
\(\Delta\) be the Archimedean vector. Write
\[
 C_{\rm fin}=\sum_{\text{finite directions }n} h_P(n)
\]
with the finite-place vectors weighted by their normalized
logarithmic sizes as above. All these quantities are fixed
up to vanishing inner-height errors on an anatomic box.
This follows because the dyadic root sizes are fixed,
and the total contribution of a coordinate's primes is
its logarithmic size. Rebox any subproducts assigned
separately before applying this estimate. The original exceptional constants
are removed from the grouped variables.

For a column with exponent vector \(w\), its primitive
integer entry in any row has archimedean logarithmic
size in base \(X\) at most
\[
                 C_{\rm fin}+\langle w,\Delta\rangle+o_X(1).
 \tag{ST13}
\]
Indeed row primitivization adds to every column the same
negative sum of the local minimum valuations.
This is \(C_{\rm fin}\) in the stated scaling.
The absolute magnitude of the original monomial itself
contributes \(\langle w,\Delta\rangle\).
At a prime used by the Taylor determinant the narrow
and wider variables are units, so the row scaling
is a local unit and does not change those divisibilities.
Using (ST13) in the Leibniz expansion of a minor with
column set \(I\) and size \(N\) gives an average row
upper bound
\[
                   C_{\rm fin}+
                   \frac1N\sum_{w\in I}\langle w,\Delta\rangle .
\]
To force rank below \(N\), use the largest such mean
over the size-\(N\) column sets. When \(N\) is the
entire column count this is the usual mean of the
column weights. All estimates precede the outer
budget limit at fixed truncation and monomial degrees.

The mixed area bound controls
\(C_{\rm fin}+h_P(\Delta)\), the primitive maximum
row expression. The column mean can be smaller when
the two free sides have different archimedean exponent
vectors. It agrees with the maximum expression when
\(\Delta=0\), including equal retained profiles with
identical multiplier assignments. Thus the broad
profile obstruction for maximum row height applies
directly there, while estimates for asymmetric
column means pose a separate optimization problem.
