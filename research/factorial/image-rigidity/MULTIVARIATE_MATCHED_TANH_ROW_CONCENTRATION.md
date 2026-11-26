# Matched tanh concentration for one factorial row

Let \(n\ge2\), \(N=n+1\), and consider a shared unequal row with
smaller exponent \(\ell\) and absolute gap \(d\).
Write \(v=d/(2\ell+d)\), \(c_h=2h/d-1\), and
\[
       T_{j,h}=\frac{(j-1)v}{j+c_hv},\qquad 1\le h\le d .
\]
These are the
[exact matched Gauss slots](MULTIVARIATE_MATCHED_TANH_INTEGER_MULTIPLICITY.md#recovering-the-gauss-grid).
The following comparisons concern unsigned moments, so either
orientation of the row gives the same result.

## The sign at the common leading scale

Put \(R=n^2/(n^2-1)\) and \(V_{p,j}=\sum_h T_{j,h}^p\).
For \(p>0\) set
\[
     F_p(w)=\sum_h(1+c_hw)^{-p},\qquad w=v/n,\qquad r=n/N .
\]
The sign of \(V_{p,N}-R^pV_{p,n}\) is the sign of
\(F_p(rw)-F_p(w)\). For \(w>0\) in the domain of both terms,
\[
       \frac{F_p(rw)-F_p(w)}w=-\int_r^1 F_p'(tw)\,dt .
\]
Since \(F_p''>0\) and \(F_p'(0)=-p\), this expression decreases
strictly with \(w\) from a positive value. If \(d\le2\) then
all \(c_h\ge0\) and the difference is always positive.
If \(d\ge3\), the abstract function has one sign change before
the first singularity of \(F_p\), so there is at most one in
the ordinary range of \(v\). At any fixed ordinary \(v>0\),
the difference is negative for sufficiently large \(p\):
the negative extreme \(c_1\) gives a larger leading value at
\(w\) than at \(rw\), and dominates the two power sums as
\(p\) tends to infinity.

In degrees \(2,3\), the ordinary row \((1,2)\) gives
\[
              V_{2,3}-R^2 V_{2,2}=\frac{41}{11025}>0,
\]
whereas \((1,6)\) gives
\[
              V_{2,3}-R^2 V_{2,2}
                    =-\frac{5509361483}{307213015824}<0 .
\]
Appending the reflection of either row doubles its unsigned
difference, giving ordinary mirrored profiles of opposite signs.

There is a uniform positive sign at fixed moment order and
sufficiently large \(n\) when the total absolute curvature \(S\)
is bounded. For \(d\ge3\), the exact sum
\(\sum_h c_h^2=(d^2+2)/(3d)\le(11/27)d\) yields
\[
 F_p'(w)\le -p+
      p(p+1)w(1-w)^{-p-2}(11/27)d .
\]
The row contributes \(2dv/(1+v)\ge dv\) to \(S\). Thus
\(n>(11/27)(p+1)S+(p+2)\) ensures \(F_p'(v/n)<0\)
for each row, by \((1-1/n)^{p+2}\ge1-(p+2)/n\).
Summing gives \(V_{p,N}>R^p V_{p,n}\) for every nonempty such
profile. In particular \(S\le12\) suffices at \(n\ge13,19,25\)
for \(p=1,2,3\), respectively.

## Matching the actual unsigned first moment

Normalize the two slot lists of a single row to sum to one.
The ratio \(T_{N,h}/T_{n,h}\) is strictly increasing in \(h\),
while each list is decreasing in \(h\). Normalization preserves
the ratio order. Since the sums agree, the upper normalized
list is at most the lower one on an initial segment, and at
least the lower one thereafter. Its leading partial sums are
therefore at most those of the lower list. Thus the lower
list majorizes the upper one, strictly for \(d\ge2\).
Every strictly convex moment of the normalized list, in particular
each power \(p>1\), decreases. The same argument works across
all rows at once if their values of \(v\) agree, by ordering
all slots by \(c_h\).

The comparison also holds for fixed positive weights on the slots
in this common-\(v\) class, using the weighted first moment for
normalization. At the single crossing of the two lists choose
a level separating the initial and final segments in both lists.
A subgradient of any convex function at that level bounds each
convex difference below by the subgradient times the corresponding
slot difference. The weighted linear differences sum to zero.
The [weighted ordinary-row-shaped collision](MULTIVARIATE_WEIGHTED_TANH_ANGLE_OBSTRUCTION.md)
uses reflections of the same row with the same positive real
weight. It therefore satisfies this stronger moment comparison and
still has coincident consecutive arguments. The example uses
fractional row multiplicities.

For general profiles, matching the first moment of the entire
array gives both signs, even with exact ordinary mirroring.
At orders \(2,3\), the two unit-gap rows \((1,2),(2,3)\)
have normalized slot lists
\[
                  (11/18,7/18),\qquad (8/13,5/13).
\]
The second power sum increases by \(53/27378\).
For the single row \((1,3)\), the lists are
\[
                  (5/9,4/9),\qquad (7/13,6/13),
\]
and the second power sum decreases by \(44/13689\).
Reflecting each row duplicates its slots and halves each
difference of normalized second moments.

Both strict signs persist at arbitrarily small relative gaps.
For an increasing sign take two unit-gap rows with distinct,
arbitrarily large bases. Their single-slot \(T_{n,1}\) and
the ratio \(T_{N,1}/T_{n,1}\) both increase strictly with \(v\),
so the larger slot gains normalized mass. For a decreasing
sign use a gap-two row and the within-row argument.
In each construction reflect every row and repeat the resulting
profile an integer number of times. Repetition divides the second
moment difference by that number. As the base tends to infinity
the curvature of one copy tends to zero, so one can choose the
copy count to approach any prescribed positive total curvature.
This gives both signs for arbitrarily diffuse ordinary mirrored
profiles in every positive open curvature interval.

These comparisons concern the matched moment arrays. The
[angular certificate](MULTIVARIATE_MATCHED_TANH_FINITE_ANGULAR_CERTIFICATE.md)
requires a positive gap between root arguments in its finite
window, and these moment comparisons give no such gap. They
specify which one-sided moment inequalities are available for
an angular argument using the complete factorial rows.
