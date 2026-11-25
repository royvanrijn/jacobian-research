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

## Column-mean bounds for the broad uniform profile

### A one-sided mixed-area inequality at the sample centroid

Let \(S\subset\mathbb R^2\) be a bounded measurable set of
positive area \(A\), with centroid \(m\), and let \(P\)
be its closed convex hull. One free letter contributes a finite
positive measure of normal vectors \(\eta\) in an open
semicircle, including the common sign appropriate to that
letter. Write \(Z=\int[0,\mathcal R\eta]\,d\nu(\eta)\) for
its directional zonotope, with \(\mathcal R\) rotation by
a positive right angle. Its Archimedean mean contribution
is
\[
 B_{\rm side}=\int [h_P(\eta)-\eta\cdot m]\,d\nu(\eta).
\]
Then
\[
                  B_{\rm side}\ge
             \frac{2\sqrt2}{3}\sqrt{A\,|Z|}.
 \tag{CM1}
\]

Here is a proof retaining the missing closing edge.
First take an atomic measure. Sort the edges
\(\mathcal R\eta\) with their positive measure weights,
and append their negative sum as the closing base.
The resulting convex cap \(K\) has area \(|Z|/2\):
its signed polygon area is half the sum of determinants
of the ordered pairs from the one-sided chain.
Position the base on \(y=0\), with \(K\) above it,
and translate \(P\) by \(-m\). Call the base length
\(l\) and the downward support height of \(P\) equal to
\(a\). The mixed-area expression is
\[
 B_{\rm side}=2V(P,K)-la.
\]
As \(t\downarrow0\), the portion of \(K+tP\) below
\(y=0\) has area \(lat+o(t)\).
For the upper bound, it is confined to depth \(ta\)
and its horizontal projection tends to that of the
base, by compactness. For the lower bound, translate
\(K\) by \(t\) times a point of \(P\) at depth \(a\),
and use that the horizontal sections immediately above
the base have length tending to \(l\).
Thus \(B_{\rm side}\) is the derivative at zero of
the area of \((K+tP)\cap\{y\ge0\}\). This set contains
\(K+t(P\cap\{y\ge0\})\), so Brunn–Minkowski gives
\[
 B_{\rm side}\ \ge\
 2\sqrt{|K|\,|P\cap\{y\ge0\}|}.
\]

Even for nonconvex \(S\), the latter cap has area at
least \(4A/9\). Indeed choose the upper cap of \(P\)
of area \(A\). Its \(y\)-centroid is at least the
\(y\)-centroid of \(S\), namely zero, since among
subsets of area \(A\) it maximizes the first moment
of \(y\). The halfspace above the centroid of this
convex cap has at least \(4/9\) of its area by planar
Grünbaum, and lies within \(P\cap\{y\ge0\}\).
This proves (CM1). Approximate the normal measure by
finite atomic measures to obtain the general bounded
case. [Grünbaum, Partitions of mass-distributions and
of convex bodies by hyperplanes](https://msp.org/pjm/1960/10-4/pjm-v10-n4-p18-s.pdf)

### Application when the area is continuous

For each free letter use the corresponding signed
normal vectors of (ST13). The vectors
\(\eta=(i,r-ai)\) lie in the open semicircle
\(a\eta_x+\eta_y=r>0\), even for signed \(i\).
The sign of the chosen letter only reverses this
semicircle and does not affect its zonoid area.
The periodic-window proof preceding (ST11) gives
\[
 |Z_W|\ge \frac34s_W^3
\]
for any retained mass \(s_W\) and any signed
nonzero integer multipliers. Let \(s=(s_A+s_B)/2\) and
\(g=2s-1/3\) as before. Adding the two applications of
(CM1) and using convexity of \(x^{3/2}\) gives
\[
 B_{\rm mean}\ge\frac{2\sqrt2}{3}\sqrt A
      \left(\sqrt{|Z_A|}+\sqrt{|Z_B|}\right)
 \ge\sqrt{\frac83 A s^3}.
 \tag{CM2}
\]
This accounts for the Archimedean centroid shift without a maximum-row bound.

In the continuous-area model the rank size satisfies
\(y\le A\). The unrestricted two-dimensional error-order
mean is at most \((2/3)\sqrt{2bwy}\). Thus under the
strict root budget \(\mu+b<g\), and \(w\le4\mu\), its
square is strictly less than
\[
       \frac89 A g^2\le\frac83 A s^3.
\]
The finite error strip at large rank is not covered by
this unrestricted mean, so this comparison is made
in the range before its last strip is crossed. The
following same calculation excludes the capped range
as well in the continuous-area model. Put \(n=J+1\),
\(t=by/w\), \(r_*=\phi_J(t)/\sqrt t\), and \(v=b/g\).
Necessity of the horizontal-width lower bound, in a
normalization where the support uses at most \(n\)
horizontal levels, gives \(v>t/(n\phi_J(t))\).
In the actual discretized setting this inequality holds
with its stated summed slice area, independently of
(CM2): on a slice of horizontal length \(\ell\), integration
of each one-sided horizontal distance costs at least
\(\ell^2/2\). Summing and using Cauchy gives the mean
height at least \(h\mathcal A/n\), where \(h\ge s\).
In the finite-strip range \(t\ge Jn/2\), the formula for \(u=n\phi_J(t)/t\) in the proof of (ST12) then
gives
\[
 [2\sqrt{v(1-v)}\,r_*]^2
    \le 4\bigl(\phi_J(t)/n-t/n^2\bigr)
    \le 2/3-4/(3n) <8/9.
\]
The same estimate applies in a continuous model that
retains this width bound. In the unlimited model (CM2)
alone suffices at every rank.

### A quantitative discrete correction

For actual integer horizontal slices, let
\(\mathcal A\) be the sum of their individual lengths,
and take \(m\) to be the centroid of their length measure.
For each side separately, choose any vector \(v_W\)
with vertical coordinate one and sweep every slice by
\([-v_W/2,v_W/2]\). These bands have disjoint interiors,
total area \(\mathcal A\), and the same centroid \(m\).
Their convex hull is \(P+[-v_W/2,v_W/2]\).
Apply (CM1) to this thickened support.
Its additional one-sided height is
\(\frac12\int|\eta_W\cdot v_W|\,d\nu_W'\).
Consequently, with the optimal tilt error
\[
 e_W=\inf_{a'}\int|r-a'i_W(r)|\,d\nu_W'(r),
\]
one obtains the exact perturbative bound
\[
 B_{\rm mean}+\frac{e_A+e_B}{2}
 \ \ge\ \frac{2\sqrt2}{3}\sqrt{\mathcal A}
            \left(\sqrt{|Z_A|}+\sqrt{|Z_B|}\right).
 \tag{CM3}
\]
The error is at most \(2(s_A+s_B)\) by taking \(a'=0\).
This comparison holds for arbitrary signed nonzero
multipliers and arbitrary sparse slices.

There is a strict uniform numerical gap in the zonotope bound preceding (ST11) for
retained masses bounded away from zero, even when
integer multipliers are allowed to grow. To prove this,
suppose masses \(s_j\to s>0\) and the zonoid area tends
to \(3s^3/4\). Fix \(s_0>0\) below the masses.
For anchored \(r_0\ge2+3s_0\) and sufficiently large
\(\lvert i_0\rvert\), the periodic-window bound improves
the coefficient \(1/3\) in the mass bound \(t/3\)
by a uniform factor. In fact the interval length
two is at most \(2|i_0|/r_0+1\) periods after rounding
up, so the covered length is at most
\(2t(2/r_0+1/|i_0|)\).
Every such anchor then contributes a fixed excess above
\(3s_j^2/2\) in its determinant integral. Its mass
must tend to zero. The set \(r_0\ge2+3s_0\) itself
has mass at least \(s_0/2\) by the original density
cap, so this gives a bounded-\(|i_0|\) anchor subset
of fixed positive mass. Integrating
\(|i_0r-i r_0|\ge 2|i|-4|i_0|\) against those anchors
shows that \(\int |i|\) stays bounded.

The marked measures in \(r,i\) are tight. Their limit
preserves the density cap, and lower semicontinuity of
the nonnegative determinant integral would give
equality in the zonotope bound preceding (ST11). At equality the conditional
integral at almost every anchor equals \(3s^2/2\),
so its sublevel mass is exactly \(t/3\) on
\(0\le t\le3s\), and its support has determinant
at most \(3s\). The same holds at all support anchors
by continuity on each integer multiplier stratum,
using the finite first moment of \(|i|\).
For any one fixed anchor the resulting bound on
\(|i_0r-i r_0|\) bounds every \(|i|\), and opposite
signs cannot coexist because then the determinant is
at least four. Choose the common sign positive and an
anchor minimizing \(r/i\) on this compact support.
Its determinants all have one sign. The periodic
sublevel windows now have only half their former
length, giving conditional mass at most \(t/6\).
This contradicts equality and proves the uniform gap.

In an endpoint test at least one free mass is at least
\(s>1/6\), so this uniform gap adds a fixed positive
multiple of \(\sqrt{\mathcal A}\) to the right side
of (CM3), beyond the continuous constant in (CM2).
It absorbs the absolute tilt error in (CM3) whenever
\(\mathcal A\to\infty\). In the unlimited Taylor
range this excludes the one-root column-mean
criterion on the broad profile. In the finite-strip
range the horizontal-width constraint and the
calculation above already give a strict upper
coefficient below the continuous constant. Thus a
possible escape for the actual sparse determinant
must have bounded summed slice area and an
asymmetric column mean at this stage. The next argument excludes
that range as well.

## The full discrete column-mean obstruction

Keep the broad profile with density \(1/6\) on \([2,4]\) for
each letter. Let \(A=\mathcal A\) be the sum of the horizontal
slice lengths, \(s=(s_A+s_B)/2\), \(g=2s-1/3\), and \(B\) the
primitive column mean. The whole-support mean is a lower
bound for the largest mean over subsets of any specified size,
so this argument also applies when the determinant uses the
worst such subset.

The generic exact-error Taylor bound (ST9) combined with the
column-mean height (ST13) cannot save on this profile, even
with arbitrary asymmetric retained submeasures, signed nonzero
integer multipliers, sparse integer slices, and the finite
strip cap. It remains possible to seek further valuation from
the actual shape of the support or from a different determinant
basis. The obstruction concerns the stated order estimate.

### Optimizing the centered thickening

Write \(Z_W\) for the one-sided zonotope area and \(E_W\) for the error
of the best tilt on side \(W\). The bounds used below are
\[
 Z_W\ge\frac34s_W^3,\qquad
 Z_W\ge\frac12s_W E_W.
\]
The first is the periodic-window estimate preceding ST11. For the
second, write \(q=r/i\) and use the measure \(|i|d\nu_W'\), of mass
\(h_W\ge s_W\). The determinant integral giving the zonotope area
is half the two-fold mean of \(|q-q'|\) against that measure. Minimizing
the inner integral in \(q'\) yields \(Z_W\ge h_W E_W/2\). This applies
also to signed multipliers.

For any \(0<\delta\le1\), sweep each support slice by
\([-\delta v_W/2,\delta v_W/2]\) for a best unit-vertical tilt \(v_W\).
The swept pieces have disjoint interiors, combined area \(\delta A\),
and the same centroid as the original sample. The convex hull is the
Minkowski sum with the sweep segment. By CM1,
\[
 B_W\ge\frac{2\sqrt2}{3}\sqrt{\delta A\max
          \left(\frac34s_W^3,\frac12s_W E_W\right)}
          -\frac{\delta E_W}{2}.
\]
Consequently
\[
 B_W\ge L(A,s_W)=
 \begin{cases}
 (2/9)A s_W,& A\le(27/8)s_W,\\
 \sqrt{2A/3}\, s_W^{3/2}-(3/4)s_W^2,& A\ge(27/8)s_W.
 \end{cases}                                                    \tag{1}
\]
Indeed take \(E_*=(3/2)s_W^2\). For \(E_W\le E_*\) the optimized
lower bound with \(Z_W\ge3s_W^3/4\) is worst at \(E_W=E_*\).
The stationary value of \(\delta\) there is \(8A/(27s_W)\).
For \(E_W\ge E_*\) the bound with \(Z_W\ge s_W E_W/2\) has interior
maximum \((2/9)A s_W\) whenever \(E_W\ge(4/9)A s_W\).
In the remaining range the endpoint \(\delta=1\) gives a concave
function of \(E_W\), whose endpoint values are the second branch
of (1) and \((2/9)A s_W\). The former is at most the latter by
completing the square. Zero mass is interpreted by continuity.

Suppose first that \(A\le9/8\), so \(s_*=8A/27\le d\). The difference
between the first and second formulas in (1) is
\((2/9)A s_W(1-\sqrt{s_W/s_*})^2\), used only for \(s_W\le s_*\).
As \(u(1-\sqrt u)^2\le1-u\) on \([0,1]\), this difference is at
most \((2/9)A(d-s_W)\). Summing yields
\[
                 B\ge (4/9)A g.                                \tag{2}
\]
For \(A\ge9/8\), both sides use the second line of (1). Factoring
\(s_W^{3/2}\), replacing \(\sqrt{s_W}\) in the negative term by
\(\sqrt d\), and applying Jensen yields
\[
 B\ge 2s^{3/2}\left(\sqrt{2A/3}-(3/4)\sqrt d\right).
\]
Since \(g^2\le3s^3\), this gives
\[
                 B/g\ge c\sqrt A-1/2,\qquad c=2\sqrt2/3.          \tag{3}
\]

### Taylor comparisons

A potential success has \(g>0\), \(w>0\), \(b>0\), \(b+\mu<g\),
and \(w\le4\mu\). Taking rank-size parameter \(y=A\) only increases
the upper bound on the Taylor gain \(G\), since the original parameter
satisfies \(y\le A\). Let \(t=bA/w\), let \(J\) denote the vertical
range of the support, put \(n=J+1\), and let \(\phi_J(t)=\Psi_J(t,1)\)
as in (ST9). Then
\[
 \frac{G}{g}\le\frac{4A\phi_J(t)}{A+4t}.                         \tag{4}
\]
This follows from \(b+w/4\le g\). Independently \(G\le bA/2\) by
\(\phi_J(t)\le t/2\). The horizontal bound \(B\ge hA/n\), where
\(h=(h_A+h_B)/2\ge s\ge g\), follows by averaging the one-sided
width cost on each slice and applying Cauchy. Success would therefore
have \(n\ge3\) and \(B<w\phi_J(t)\).

For \(t\le1\), \(\phi_J(t)=t/2\); for \(1\le t\le3\),
\(\phi_J(t)=t/4+1/2-1/(4t)\). The function is concave on its
whole domain, as follows either from the strip formula or by
integrating its inverse count. On the second branch
\(9\phi_J(t)-4t\) is maximized at \(t=3/\sqrt7\), with value \(A_0\).
On the first branch it is at most \(1/2<A_0\). For \(t\ge3\) it
decreases because the derivative at \(t=3\) is \(5/18<4/9\).
Thus (4) is at most \((4/9)A\) whenever \(A\ge A_0\).
Together with (2) this excludes \(A_0\le A\le9/8\).

Now assume \(A\ge9/8\), and write \(Y=c\sqrt A-1/2\).
For \(t\le1\), (4) is at most \(2A/(A+4)\). This is at most
\(1/2\le Y\) when \(A\le4/3\). For \(A\ge4/3\) it is at most
\(\sqrt A/2\le Y\). For \(1\le t\le3\) and \(A\le9/4\), we have
\(1/2\le Y\le1\). That (4) is at most \(Y\) is equivalent to
\((4Y-A)t+A/t\ge A(2-Y)\). By AM--GM it suffices that
\(16Y\ge A(4+(2-Y)^2)\), which also implies \(4Y>A\). The last
inequality follows from \(A=(9/8)(Y+1/2)^2\),
\((Y+1/2)^2/Y\le9/4\), and \(4+(2-Y)^2\le25/4\).

Before the finite strip cap is crossed, the unlimited strip
count is at least \((\tau^2+\tau)/2\). Integrating the resulting
upper bound for the inverse count, and then using
\(\sqrt{2u+1/4}\le \sqrt{2u}+1/(8\sqrt{2u})\) for \(u>0\),
gives \(\phi_J(t)\le c\sqrt t-1/2+1/(4\sqrt{2t})\).
For \(t\ge1\) this is at most \(c\sqrt t-8/25\), and for
\(t\ge3\) it is at most \(c\sqrt t-19/48\).
If \(\phi_J(t)\le c\sqrt t-\lambda\), the maximum over real
\(\sqrt t\) of the right side of (4) with that bound is
\(\sqrt{c^2A+4\lambda^2}-2\lambda\). For \(\lambda=8/25\)
and \(A\ge9/4\) this is at most \(Y\), since
\(c\sqrt A\ge\sqrt2\) and \(39/100\le(7/25)\sqrt2\).
For \(\lambda=19/48\) and \(A\ge9/8\), the same follows from
\(c\sqrt A\ge1\) and \(13/24\le7/12\).
This completes the comparison in the uncapped range.

For the capped range \(t\ge Jn/2\), write \(v=b/g\). A potential
success and the horizontal bound give \(v>v_0=t/(n\phi_J(t))\).
The formula
\[
 \phi_J(t)=\frac{t}{2n}+\frac{n-1}{2}
                  -\frac{n(n^2-1)}{24t}
\]
shows \(\phi_J(t)/t\le (4J-1)/(3Jn)<4/(3n)\), since this ratio
decreases by concavity and vanishes at zero before division.
In particular \(v_0>1/2\). From \(wb\le4g^2v(1-v)\) and its
decrease for \(v\ge v_0\), we get
\[
 \frac{G^2}{A g^2}
    \le4\left(\frac{\phi_J(t)}n-\frac{t}{n^2}\right)
    \le\frac23-\frac4{3n}.
\]
The last expression before the upper bound decreases for
\(t\ge Jn/2\) and has its maximum there. Also \(G/g\le4A/(3n)\)
by \(b\le g\) and the bound for \(\phi_J(t)/t\).
Set \(x=\sqrt A\), \(c_n=\sqrt{2/3-4/(3n)}\), and
\(X_n=1/(2(c-c_n))\). If \(x\ge X_n\) the first bound is at
most \(cx-1/2\). From \(\sqrt{9/8}\) to \(X_n\) use the
second bound, by concavity of \(cx-1/2-4x^2/(3n)\) and the
endpoint checks. At the left endpoint the difference is
\(1/2-3/(2n)\ge0\). At the right endpoint one needs
\(c_n(c-c_n)\ge2/(3n)\). Squaring the equivalent comparison
\(c c_n\ge2(n-1)/(3n)\) has difference
\(4(n-3)(n+1)/(27n^2)\ge0\).
This proves \(G/g\le Y\) also in the capped range, finishing
the exclusion by (3).

### The remaining small-area range

The common centroid supplies a two-sided tilted segment even when
the free directions and the wider slice shapes are asymmetric. This
completes the exclusion for \(A<A_0\).

Suppose the criterion succeeds and let \(Q=\sum_j p_j^2\), where
\(p_j=\ell_j/A\) is the fraction of the sample on wider level \(j\).
The horizontal argument used above retains the stronger bound
\(B\ge h A Q\): on a slice the direction's horizontal contribution
to the average difference between its support and a sample point
is at least \(|i_W|\ell_j/2\). This holds on each side separately
and summing gives the assertion. Since \(G\le bA/2\) and
\(b<g\le s\le h\), success implies \(Q<1/2\).

Choose an integer median \(j_0\) of the sample levels and split the
sample into two sets of equal length, respectively below and above
the median, splitting its own level if needed. Put
\[
 D=\sum_j p_j|j-j_0|.
\]
The two barycenters lie in \(P\), are opposite about the full sample
centroid, and have vertical coordinates differing by \(2D\).
Thus their segment, after translating to the sample centroid and
shearing back to the raw normal coordinates \((i,r)\), has endpoints
\(\pm D(-a,1)\) for some real \(a\) whenever \(D>0\). On either
side this gives the same absolute tilt error, so
\[
 B\ge 2 D E,\qquad
 E=\frac12\sum_{W=A,B}\int|r-a i_W(r)|\,d\nu_W'(r).               \tag{5}
\]

The discrete distribution has the exact estimates
\[
                       D\ge1/3,\qquad DQ\ge1/6.                 \tag{6}
\]
To prove them translate \(j_0\) to zero. If \(D\le1/2\), the
three-point probability \(q_0=1-D,\ q_{\pm1}=D/2\) is a
pointwise majorant of \(1-D-(1-3D/2)|j|\). Thus
\(\sum p_jq_j\ge1-2D+3D^2/2=\sum q_j^2\).
Cauchy yields \(Q\ge1-2D+3D^2/2\). With \(Q<1/2\), this
implies \(D>1/3\), and \(D(1-2D+3D^2/2)\) is increasing,
with value \(1/6\) at \(D=1/3\). If \(D>1/2\) and
\(Q\ge9/25\), the desired bound for the product follows
immediately. Otherwise convolve the integer distribution with a
uniform variable on \([-1/2,1/2]\). Its density has squared
integral \(Q\) and its absolute first moment is \(D+p_0/4\).
For any real density with squared integral \(Q\) and absolute
moment \(m\), Cauchy applied to \((1-|x|/R)_+\) gives
\(1-m/R\le\sqrt{2QR/3}\). Taking \(R=3m\) gives
\(Qm\ge2/9\). Using \(p_0\le\sqrt Q\), we get in our case
\[
              DQ\ge2/9-Q^{3/2}/4
                  \ge2/9-27/500>1/6.
\]
This establishes (6).

We next show \(s>3/10\). Let \(s_*=8A/27\), which is less
than \(d/2\) since \(A<A_0<8/15\). If some mass \(s_W\le s_*\),
write that mass as \(v\), the other mass as \(u\), and
\(x=v/s_*\). The other mass \(u\) exceeds \(s_*\), since
\(u+v>d\). By (1), a necessary condition for
\(B<gA/2\) is
\[
       \frac29\{u+v(2\sqrt x-x)\}<\frac{u+v-d}{2}.
\]
The right minus left side increases with \(u\le d\). Even at
\(u=d\) the inequality fails because
\[
 v(9/4-2\sqrt x+x)
   =s_* x(9/4-2\sqrt x+x)\le (5/4)s_*<d.
\]
The function of \(x\) in the middle has derivative
\(9/4-3\sqrt x+2x>0\). Both masses must therefore exceed
\(s_*\). Now (1) gives \(B\ge(4/9)A s\), and
\(B<gA/2\) forces \(s>3/10\).

By (5), \(B<G\le sA/2\) and (6) give
\(E/s<A/(4D)\le3A/4<2/5\).
Here \(A_0<8/15\) follows from \(119^2<7\cdot45^2\).
We may now apply the folded-distance argument of ST4--ST5
to the common tilt \(a\), independently of the original
chart tilt. Its hypotheses are \(3s\ge4/5\) and
\(E/s<2/5\), and it gives
\[
                             Eh\ge(6/5)s^3.                     \tag{7}
\]
For clarity, let \(a_*=|a|\). The retained moment of the average
of the two side measures is at least \((2+3s)s\), so
\(a_* h+E\ge(2+3s)s\) and \(a_*h>(12/5)s\).
For \(a_*\ge2\), folding \([2,4]\) into one period gives
\(E\ge(3/2)s^2\), and \(h\ge s\) suffices. For \(0<a_*\le1\),
folded numerical-distance density at most \((2+a_*)/(3a_*)\)
gives \(E\ge a_*s^2/2\). For \(1<a_*<2\) split \([2,4]\)
into one period and a remainder of length \(2-a_*\).
Each part's folded numerical-distance density is at most
\(1/3\). Optimizing the retained split gives
\(E\ge(3/2)s^2((1-z)^2+z^2)\), where
\(z=\min(1/2,(2-a_*)/(6s))\). In the case \(z=1/2\),
\(a_*\le2-3s\le6/5\). Otherwise \(a_*\le2-(8/5)z\) and
\((1-z)^2+z^2\ge(2-(8/5)z)/3\) by completing the square.
Thus \(E\ge a_*s^2/2\) in this last range as well. These density
bounds apply to signed marks because \(|r-ai|\) is at least
the folded distance to \(a_*\mathbb Z\). The case \(a_*=0\)
violates the retained moment and the strict bound on \(E/s\).

Finally, multiplying \(B\ge h A Q\) by (5), and comparing
with \(B^2<b^2A^2/4\), forces
\[
 A>\frac{8 DQ Eh}{b^2}
   \ge\frac{48 DQ s^3}{5 b^2}
   \ge\frac8{15}.
\]
The last step uses (6) and \(b^2<g^2\le3s^3\). This contradicts
\(A<A_0<8/15\). Thus the full sparse ST13 criterion cannot
save on this broad profile under the stated one-root exact-error
hypotheses.

[The weighted interpolation refinement](ABC_WEIGHTED_HILBERT_TAYLOR.md) uses the actual exponent set to compute the best universal local Taylor weight. It rules out two-level single-root supports and nested toric supports with any fixed number of roots on the same broad profile.
