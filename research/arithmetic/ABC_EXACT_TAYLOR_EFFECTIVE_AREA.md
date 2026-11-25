# Effective horizontal area for arbitrary one-root exact Taylor supports

Take the broad profile with density \(1/6\) on \([2,4]\)
for each letter, and permit arbitrary signed nonzero multipliers
and asymmetry between the two free letters. The support is an
arbitrary finite list of one-root monomials, with no restriction
on the number, location, or nesting of its wider levels.

Let \(L\) be the small budget scale, \(x=I-qj\) the torically
sheared horizontal exponent, \(m_x\) the number of occupied
wider levels above \(x\), and
\[
 n_k=\#\{x:m_x>k\},\qquad N=\sum_k n_k,\qquad A=LN,\qquad
 K=\frac L N\sum_{k\ge0} n_k^2,\qquad
 \bar k=\frac1N\sum_{k\ge0} k n_k .
\]
In the limiting normalized first determinant, a strict endpoint
saving by the exact universal support-dependent Taylor weights
requires \(K>4/7\). In particular the range of sufficiently small
supports excluded in this way is genuinely nonnested. Since
\(K\le A\), this contains the entire region \(A\le4/7\), and
also larger supports with sufficiently high projection overlap. The shifted-median refinement below excludes \(K\le1\).

## Absorbing the toric overlap cost

Let \(s_W\le d=1/3\) be the two retained masses, \(s\) their
average, \(g=2s-d>0\), \(h\) the average absolute multiplier
integral, and \(R\) the sum of their retained numerical moments.
Let \(b,\mu,w\) be the shallow root record cost, full lock cost,
and full modulus weight. The contemplated saving requires
\[
       \mu+b<g,\quad w\le4\mu,\quad
       R\ge4s+6s^2>4g,\quad h\ge s\ge g .
\]
The moment lower bound follows by filling the bottom of
\([2,4]\) separately on each side at density \(1/6\), then
applying convexity. Its difference from \(4g\) is at least
\(6(s-1/3)^2+2/3\).

The vertical-first interpolation basis from the finite weighted
Hilbert calculation consists of \(x^a j^k\) for \(0\le a<n_k\).
It is independent by descending in the \(j\)-degree, using that
a polynomial of degree \(k\) vanishing on a fiber of more than
\(k\) levels has every coefficient zero there. Hence the exact
universal Taylor mean \(G\) is bounded above by
\[
 G\le \frac{bL}{2N}\sum_k n_k(n_k-1)+w\bar k
      \le bK/2+w\bar k.
\]
The Archimedean column mean \(B\) independently satisfies
\(B\ge R\bar k\). For each direction and each actual projection
fiber, compare the evaluation with the endpoint of that
fiber favored by the vertical sign. The horizontal coordinate
is unchanged. The average one-sided distance along \(m_x\)
distinct integer levels is at least \((m_x-1)/2\), and its
average over columns is precisely \(\bar k\).

Consequently \(B<G\) implies
\[
 B< b_{\mathrm{eff}}K/2\le gK/2,\qquad
 b_{\mathrm{eff}}=\frac{bR}{R-w}.
\]
The last inequality is equivalent to \(g w\le R(g-b)\) and
follows from \(w\le4(g-b)\) and \(R\ge4g\). This eliminates the
overlap premium from the comparison without a loss in the
root budget.

## Excluding effective area at most \(4/7\)

Use the geometric lower bounds from the full discrete
column-mean obstruction for ST13. For a finite lattice support,
give each column a horizontal interval of length \(L\) centered
on it. Intervals at the same wider level have disjoint interiors.
The sample then has exactly summed area \(A=LN\), the same
centroid, and the support increases \(B\) by \(hL\). This is
\(o(1)\) in the limiting normalized comparison at bounded
multipliers. The formulas below are written on that limiting
scale. At finite \(L\), one can replace \(K\) in the upper bound
\(B<gK/2\) by \(K+2hL/g\) for the thickened support.

For \(A\ge9/8\), the large-area geometric estimate is
\(B/g\ge(2\sqrt2/3)\sqrt A-1/2\ge1/2\). This rules out
\(K\le4/7\). Suppose \(A<9/8\), and put \(s_*=8A/27\le d\).
The sidewise centered thickening lower bound in ST13
starts at \((2/9)A s_W\) and can lose at most
\[
 (2/9)A s_W(1-\sqrt{s_W/s_*})^2
     \le (2/9)A s_*/16
\]
on the side where \(s_W\le s_*\). Indeed
\(x(1-\sqrt x)^2\le1/16\) on \([0,1]\).
Thus
\[
 B\ge(2/9)A(2s-s_*/8).
\]
Together with \(B<gK/2\le gA/2\) this forces
\(s>3/10-s_*/20\ge17/60>4/15\).

Let \(p_j\) be the fraction of support length at wider level
\(j\), \(Q=\sum_j p_j^2\), and \(D=\sum_j p_j|j-j_0|\) at an
integer vertical median \(j_0\). The horizontal inequality
\(B\ge h A Q\) gives \(Q<1/2\). The discrete median estimate in
ST13 then gives
\[
                        D\ge1/3,\qquad DQ\ge1/6.
\]
The barycenters of the two equal halves of the sample are
opposite about the full sample centroid, with vertical
displacements \(\pm D\). Their chord gives \(B\ge2DE\),
where \(E=\frac12\sum_W\int|r-a i_W|\,d\nu'_W\) for a
common tilt \(a\). Thus \(E/s<(3/4)K\).

Write \(u=3s\in(4/5,1]\). For \(0<|a|\le2\), the
folded-distance argument of ST4--ST5 gives
\(E\ge |a|s^2/2\). This estimate only requires \(u\ge4/5\).
Indeed for \(|a|\le1\), use the folded numerical density
bound \((2+|a|)/(3|a|)\le1/|a|\). For \(1<|a|\le2\),
split the interval into one full period and its remainder,
each with folded density at most \(1/3\). With
\(z=\min(1/2,(2-|a|)/(2u))\), it gives
\(E\ge(3/2)s^2((1-z)^2+z^2)\). If \(z=1/2\), then
\(|a|\le2-u\le6/5\). Otherwise
\(|a|\le2-(8/5)z\) and
\(3((1-z)^2+z^2)\ge2-(8/5)z\), completing the proof.
The retained moment gives \(|a|h+E\ge(2+u)s\). Combining,
\[
 \frac{Eh}{g^2}>
   \frac{(2+u-3K/4)u^3}{6(2u-1)^2}
   \ge \frac12-\frac K8 .
\]
For the last comparison write \(C=2-3K/4\). The logarithmic
derivative of \((C+u)u^3/(2u-1)^2\) is
\(1/(C+u)+(2u-3)/(u(2u-1))<0\) on this range when
\(K\le4/7\). The last summand is at most \(-1\), and
\(C+u>1\). Evaluation at \(u=1\) gives the stated bound.
For \(|a|\ge2\) one full period bounds the folded density
by \(1/3\), so \(Eh\ge(3/2)s^3\), which gives
\(Eh/g^2\ge1/2\). The case \(a=0\) violates the retained
moment and \(E/s<3K/4\).

Finally multiply \(B\ge h A Q\) and \(B\ge2 D E\) and
compare with \(B^2<g^2 K^2/4\). Then
\[
     K\ \ge\ \frac{K^2}{A}\
       >\ 8 DQ\frac{Eh}{g^2}\
       \ge\ \frac23-\frac K6 ,
\]
which forces \(K>4/7\). A fixed positive slack in the intended
root-budget saving also preserves the exclusion at the boundary
under horizontal lattice thickening. If \(\mu+b\le g-\epsilon\),
then \(b_{\mathrm{eff}}\le g-\epsilon/3\), since \(4g/R\le2/3\)
on the broad profile. In general the finite-scale correction
in \(K\) is \(o(1)\) as described above.

## A second exact orientation

For completeness there is an independent interpolation upper
bound controlled by wider slice sizes even when the horizontal
point sets do not overlap. Let \(N_j\) be the count on wider
level \(j\). Transpose the previous fiber interpolation, pairing
\(x^a\) with the polynomial basis in \(j\) on levels with
\(N_j>a\). The exact gain obeys
\[
 G\le \frac{bL}{2N}\sum_j N_j(N_j-1)+wV,
 \qquad V=\sum_{j<\ell}\min(p_j,p_\ell).
\]
Indeed \(\sum_a \binom{\#\{j:N_j>a\}}2
=\sum_{j<\ell}\min(N_j,N_\ell)\). Both orientations can be
used as weight-by-weight rank lower bounds when computing
the exact local filtration. There is also \(V\le2D\). At
each superlevel \(u\) of the masses \(p_j\), the number of
pairs is \(\binom{t}{2}\), while the sum of distances from the
median over \(t\) distinct integers is at least
\(\lfloor t^2/4\rfloor\ge\frac12\binom{t}{2}\).
Integrate in \(u\).

## The shifted centroid estimate

Keep the effective-area notation from above. Write
\(A=LN\), let \(B\) be the primitive column mean, let \(s\) be
the mean free retained mass, and let \(h\) be the mean absolute
multiplier integral. A possible success already gives
\[
 B<gK/2,\qquad B\ge h A Q,\qquad h\ge s\ge g,\qquad K\le A,
\]
where \(g=2s-1/3\), and \(Q=\sum_j p_j^2\) for the slice
fractions \(p_j\). In particular \(Q<1/2\).
The horizontal thickening correction for an actual lattice
sample is \(hL\), as in the \(4/7\) proof.

Let \(j_0\) be an integer median, let
\(D=\sum_j p_j |j-j_0|\), and split the \(p_j\) into upper and
lower parts \(w_{j,+},w_{j,-}\) of equal total mass \(1/2\).
All nonmedian levels are assigned to their respective sides.
Set
\[
 U=\sum_j p_j w_{j,+},\quad V=\sum_j p_j w_{j,-},\quad
 \sigma=U-V,\quad F=Q-|\sigma|=2\min(U,V).
\]
The ordinary median chord gives \(B\ge2DE\), where \(E\) is the
mean absolute multiplier error at a single common tilt \(a\).
There is also the additive inequality
\[
                          B\ge h A F+2DE.                       \tag{1}
\]

To prove it, let \(\ell_j=A p_j\), and let \(\bar x_j\) be the
horizontal mean position on slice \(j\). The two points
\((\bar x_j\pm\ell_j/2,j)\) belong to the sliced support hull:
the mean distance of a length-\(\ell_j\) subset of a line to
either of its extremes is at least \(\ell_j/2\).
For a normal with horizontal sign \(\epsilon\), replace the
slice center by \(\bar x_j+\epsilon\ell_j/2\), and form the
full mean and two half barycenters of the resulting points.
The full mean moves horizontally by \(\epsilon A Q/2\).
The upper half barycenter relative to the full mean moves by
\(\epsilon A\sigma/2\). The two half barycenters are still
opposite about their full mean and have vertical displacements
\(\pm D\). Comparing the support function to the larger of
these two points supplies the original median-chord absolute
value, minus at most \(|i|A|\sigma|/2\), after adding the
full-mean displacement contribution \(|i|A Q/2\). Integrate
over the two free sides to obtain (1). This selects a virtual
barycenter separately for each direction and uses the same
original chord tilt \(a\) throughout.

## Three bounds for the median split

The discrete median inequality already gives \(D\ge1/3\)
when \(Q\le1/2\). We use two refinements and one split-energy
bound:
\[
 |\sigma|\le1/4,\qquad
 DQ\ge(2-Q)/9,\qquad
 DQ\ge\frac{1}{18t(1-t)}-\frac Q2,\quad
                   t=\frac{\min(U,V)}Q.                         \tag{2}
\]

For the first inequality, if \(p_{j_0}\le1/2\), every \(p_j\)
is at most \(1/2\) by the median property, and \(U,V\le1/4\).
If \(p_{j_0}>1/2\), write \(p_0=p_{j_0}\). Equal half masses
express \(\sigma\), up to its sign, as the difference of the
two nonnegative tail sums
\(\sum p_j(p_0-p_j)\). Each sum is at most
\(p_0(1-p_0)\le1/4\).

For the second, when \(1/3\le D\le1/2\), the three-point
triangular comparison gives \(Q\ge1-2D+3D^2/2\). Multiplying
by \(9D+1\), subtracting 2, and factoring gives
\((D-1/3)(27D^2/2-12D+3)\ge0\).
For \(D\ge1/2\), translate the median to zero, use the
discrete tent \(f(j)=(1-|j|/R)_+\) with \(R=3D\), and apply
Cauchy. Since \(\sum_j p_j f(j)\ge2/3\), it suffices to
bound its squared norm by \(2R/3+2/9\).
In fact, for \(\theta=R-\lfloor R\rfloor\), direct summation is
\[
 \sum_{j\in\mathbf Z} f(j)^2=
 \frac{2R}{3}+\frac{1}{3R}
   +\frac{-2\theta^3+3\theta^2-\theta}{3R^2}.
\]
The last numerator has absolute value at most \(1/4\).
For \(R\ge2\) this gives the stated bound. For
\(R=3/2+z\in[3/2,2]\), the last numerator is
\(z/2-2z^3\le z/2\), whereas the allowable numerator
\(2R^2/3-R\) is \(z+2z^2/3\). Thus \(Q(9D+1)\ge2\)
also in this range.

For the third inequality, separately smooth each half's
discrete distances \(k\ge0\) from the median into \([k,k+1]\).
The resulting density has mass \(1/2\), absolute moment
\(D_\pm+1/4\), and squared integral at most \(U\) or \(V\).
For a nonnegative density of mass \(M\), moment \(m\) and
squared integral \(Q_*\) on \([0,\infty)\), the one-sided
tent and Cauchy give \(mQ_*\ge4 M^3/9\). In detail its tent
integral is at least \(M-m/R\) and at most
\(\sqrt{Q_*R/3}\), and take \(R=3m/M\). Hence
\(D\ge(1/U+1/V)/18-1/2\). Multiplying by \(Q\) gives (2).

## Excluding effective area at most \(16/19\)

Assume \(K\le16/19\). If \(A\ge9/8\), the large-area
column-mean bound gives \(B/g\ge1/2\), already impossible.
For \(A<9/8\), the centered-thickening estimate loses at most
\((2/9)A s_*/16\) per side from \((2/9)A s_W\), with
\(s_*=8A/27\le1/3\). As in the \(4/7\) proof this gives
\(s>17/60>4/15\). Also \(E/s<K/(4D)\) from the plain
median chord.

The same folded estimate as in the \(4/7\) proof, but
retaining \(D\) in this last inequality, gives
\[
                          \frac{Eh}{g^2}\ge
                              \frac12-\frac{K}{24D}.            \tag{3}
\]
Indeed for \(0<|a|\le2\), it gives \(E\ge |a|s^2/2\);
combine this with \(|a|h+E\ge(2+3s)s\).
Writing \(u=3s\), the necessary comparison reduces to
minimizing \((2+u-K/(4D))u^3/(2u-1)^2\) on
\(4/5<u\le1\). Its logarithmic derivative is negative:
the positive summand \(1/(2+u-K/(4D))<1\), since
\(D\ge1/3\) and \(K\le16/19\), and the other summand
\((2u-3)/(u(2u-1))\le-1\). The endpoint \(u=1\) gives
(3). For \(|a|\ge2\), the one-period fold gives
\(Eh/g^2\ge1/2\). The zero tilt contradicts the small
absolute error.

Now \(B<gK/2\) and \(B\ge hAQ\) give
\(h/g<K/(2AQ)\). Combining (1) and (3), and using \(A\ge K\),
forces
\[
             F+\frac{2}{K}DQ-\frac Q6<\frac12.                  \tag{4}
\]
This contradicts a universal consequence of (2):
\[
                    F+\frac{19}{8} DQ-\frac Q6\ge\frac12.       \tag{5}
\]
To see it, put \(p=19/8\) and assume failure of (5).
Using \(DQ\ge(2-Q)/9\) and \(F\ge Q-1/4\) yields
\(Q<16/41\). Using the same bound for \(DQ\) and \(F=2Qt\)
then yields \(t<23/128\). The remaining bound in (2)
gives
\[
 F+pDQ-Q/6
 \ge \frac{p}{18t(1-t)}
      -(p/2+1/6-2t)\frac{16}{41}.
\]
The right side decreases up to \(t=23/128\). For an exact
derivative check, dropping the factor \((1-t)^2\) in the
denominator of the negative derivative, its magnitude there
still exceeds \(32/41\) by \(358672/195201\).
At \(t=23/128\) the displayed lower bound exceeds \(1/2\)
by \(24719/3564540\). This proves (5), completing the
contradiction to (4).

The argument is on the limiting normalized scale. A fixed
strict root-record budget margin leaves positive slack after
toric overlap absorption, so the \(o(1)\) horizontal lattice
thickening shift does not alter the boundary assertion.

## The unit-area strengthening

The same argument actually excludes \(K\le1\). Suppose such
a support saves. The previous overlap reduction and
large-area comparison show \(A\le9/8\), \(u=3s>17/20\),
\(B/g<1/2\), \(h/g<1/(2Q)\), and \(E/s<v=K/(4D)\le3/4\).
In what follows \(a\ge0\) denotes the absolute value of the
common median tilt. The case \(a=0\) contradicts the retained
moment, so assume \(a>0\).

### Tilts of size at most one

The retained moment gives \(h/g\ge(3-v)/a\), since
\((2+u-v)u/(2u-1)\) is decreasing in \(u\) on the range
under consideration. The small-period folded density estimate
gives \(E/g\ge a/(2(2+a))\), using
\(u^2/(2u-1)\ge1\). The horizontal upper constraint implies
\(Q\le a/(2(3-v))\le2/9<1/4\).
The last estimate in (2) yields
\[
 DF=2t\,DQ\ge\frac{1}{9(1-t)}-tQ\ge1/12.
\]
The final lower bound follows by minimizing the right side
with \(Q=1/4\), at \(t=1/3\).

More precisely \(h A Q\le B<gK/2\) gives
\(h/g<K/(2 A Q)\), and we use \(A\ge K\). Inequality (1)
therefore gives
\[
 \frac{B}{gK}\ge
   \frac{(3-v)F}{a}+\frac{Da}{3K}
   \ge2\sqrt{\frac{(3-v)DF}{3K}}\ge\frac12.
\]
Here \(2/(2+a)\ge2/3\), \(3-v\ge9/4\), and \(K\le1\).
This contradicts \(B<gK/2\).

### Tilts of size at least two

The one-period fold gives \(E/g\ge u^2/(2(2u-1))\ge1/2\).
Then \(B/(gK)\ge F+D/K\ge F+D\) by (1).
A possible success would have \(D\le1/2\). The triangular
bound in this range and \(F\ge Q-1/4\) imply
\[
 F+D\ge3/4-D+(3/2)D^2\ge7/12,
\]
again impossible.

### Intermediate tilts

For \(1<a<2\), put \(z=\min(1/2,(2-a)/(2u))\). The
period-and-remainder folded estimate gives \(E\ge c(a,u)a s^2\),
where \(c(a,u)=(3/(2a))((1-z)^2+z^2)\). Since
\(a h\ge(2+u-v)s\),
\[
 \frac{Eh}{g^2}\ge
 c(a,u)\frac{(2+u-v)u^3}{3(2u-1)^2}
 \ge\frac{31}{50}\left(1-\frac{K}{12D}\right).                  \tag{6}
\]
For the second inequality, the logarithmic derivative of
\(c(a,u)\) in \(u\) is zero on the capped branch and on the
remaining branch is at most \(1/(2u)\). This follows from
\(4z(1-2z)\le1-2z+2z^2\), whose difference is a positive
quadratic. The logarithmic derivative of the whole middle
expression is thus at most
\[
 \frac{7}{2u}-\frac4{2u-1}+\frac{1}{2+u-v}
 \le-\frac12+\frac1{2+u-v}<0.
\]
The inequality before the strict sign follows from \(u\le1\),
and \(2+u-v>2\) follows from \(u>17/20\) and \(v\le3/4\).
Thus it suffices to evaluate at \(u=1\), where
\(c(a,1)\ge31/50\). In fact on writing \(a=2-2z\), this
is equivalent to \(150z^2-88z+13\ge0\), which has negative
discriminant.

As in (4), the horizontal constraint, (1), and (6) reduce
success to
\[
 F+\frac{62}{25}DQ-\frac{31}{150}Q<\frac12 .
\]
This inequality is ruled out by (2). In detail put
\(p=62/25\) and \(q=31/150\). Applying \(DQ\ge(2-Q)/9\)
and \(F\ge Q-1/4\) first yields \(Q<179/466\).
Applying the same bound for \(DQ\) and \(F=2Qt\)
then yields \(t<125/716\). The half-density estimate gives
\[
 F+pDQ-qQ\ge
 \frac{p}{18t(1-t)}-\left(\frac p2+q-2t\right)\frac{179}{466}.
\]
This expression decreases up to \(t=125/716\). As an
exact derivative bound, \(p(1-2t)/(18t^2)-2(179/466)\)
at that endpoint is \(1780700413/819140625>0\), so
retaining the additional factor \(1/(1-t)^2\) only
strengthens the sign. Its value at the endpoint exceeds
\(1/2\) by \(267508951/7745793750\). This completes the
exclusion at effective area \(K\le1\).
