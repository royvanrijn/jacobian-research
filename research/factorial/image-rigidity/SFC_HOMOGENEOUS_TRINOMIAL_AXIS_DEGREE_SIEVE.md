# A fractional degree test for the initial homogeneous trinomial window

Fix three nonnegative integers \(u_1,u_2,u_3\) and assign each
column to one of the two coordinate axes. In degree \(D\),
a column has exponents \((u_j,D-u_j)\) or \((D-u_j,u_j)\)
according to the assignment. Exponents assigned to the same
side must be distinct. Consider ordinary integer degrees
for which the three columns are nonnegative and distinct.

**Theorem.** For each such profile, there is a computable finite
list of degrees that can fail the initial factorial window
\(1,2,3\), over any characteristic-zero field. Together with the
[uniform interior exclusion](SFC_UNIFORM_HOMOGENEOUS_TRINOMIAL_INITIAL_CORNER.md),
this gives an absolute bound on the degree of a failing
homogeneous bivariate trinomial. The per-profile bound below
does not require that the limiting gamma system be nonvanishing.

**Completion target and priority.** Finish the initial homogeneous
bivariate trinomial theorem: exclude simultaneous zeros at orders
\(1,2,3\) in every degree, over every characteristic-zero field.
The absolute degree bound makes this a finite global reduction in
the integer supports. Each remaining support must be checked exactly
over all projective coefficient directions; sampling coefficients
does not discharge that test. The bound has not been numerically
evaluated, the complete residual has not been checked, and the
practical cost is unknown. The finite reduction is proved; the
nonvanishing theorem remains open.

There are two completion routes. Evaluate or sharpen the nested
cutoffs in Section 3 and run the exact residual tests, or prove the
uniform [signed remainder inequality](SFC_HOMOGENEOUS_TRINOMIAL_SHORT_SPANS.md)
\(R(L,H;a,d)<0\) for integer gaps and nonnegative heights.
The latter would give a direct finish without evaluating the
absolute bound. This initial-window target has priority over the
separate [arbitrary-dimensional two-monomial project](MULTIVARIATE_ALL_DEGREE_REDUCTION.md),
whose finitely many remaining terminal orders still admit infinitely
many exponent profiles. Neither route here asserts windows at
arbitrary positive starting orders.

## 1. Isolating the highest exponent at a fractional degree

Let \(M=\max_j u_j\), which is at least one because some side
has two columns. Write \(F_m(D,C)\) for the order-\(m\) moment
of the profile divided by \((mD)!\) at ordinary degrees.
On a term \(C^\alpha\) supported in just one side,
its coefficient, including the multinomial coefficient, is
\[
               \binom m\alpha\frac{U!}{(mD)_{\underline U}},
           \qquad U=\sum_j\alpha_j u_j ,
 \tag{1}
\]
where \((x)_{\underline U}=\prod_{j=0}^{U-1}(x-j)\).
The empty product is one. Define \(P_m(D,C)\) over
\(\mathbb Q(D)\) by keeping exactly these terms, with the
same coefficients. In particular \(P_1=F_1\), and if
all columns are assigned to one side then \(P_m=F_m\)
at every ordinary degree.

Set \(D_*=M-1/3\). For \(m=1,2\) the factors in (1) are
positive at \(D_*\), since \(mD_*>(mM)-1\).
In order three the factor \(3D-(3M-1)\) gives a simple
pole in exactly the terms with \(U=3M\). Within each side
only one exponent can be \(M\), so those terms are the
pure cubes of the columns with exponent \(M\). They
all have the same nonzero residue. After multiplying
\(P_3\) by \(D-D_*\), it specializes, up to a common
nonzero factor, to the sum of those pure cubes.

If only one column is maximal, its coefficient must
vanish. The other two columns have a nonzero order-one
linear equation with positive real coefficients. Thus
any common complex direction is proportional to a
nonzero real vector. The second form on these columns
is strictly positive definite over the reals. For two
columns on one side, it is a positive multiple of the
square integral on \((0,1)\) of their linear combination
of \(t^{u_j}(1-t)^{D_*-u_j}\). This follows from the beta
integral, and their quotient is nonconstant. The
integrals converge because \(2(D_*-u_j)>-1\).
For columns on opposite sides the second form is
diagonal with positive entries. Hence there is no
common direction in the unique-maximum case.

If the maximum occurs twice, there is a pair with
coefficients \(X,Y\) on one side, and a singleton \(Z\)
on the other, with \(Y,Z\) the equal maximal columns.
At \(D_*\) the first two forms have the shape
\[
             A_1X+B_1(Y+Z),\qquad
             A_2X^2+2B_{12}XY+B_2(Y^2+Z^2)
\]
with all displayed constants positive. The specialized
third form is a nonzero multiple of \(Y^3+Z^3\).
If \(Y=0\), both \(Z\) and then \(X\) vanish. Otherwise
put \(Y=1,\ Z=\xi,\ r=B_1/A_1>0\); the first equation
gives \(X=-r(1+\xi)\). At \(\xi=-1\) the quadratic
equals \(2B_2>0\). For the other roots of
\(\xi^3=-1\), we have \(\xi^2-\xi+1=0\).
The quadratic becomes
\[
       (3A_2r^2+B_2)\xi-2B_{12}r(1+\xi).
\]
This is a real linear expression in a nonreal \(\xi\)
with nonzero constant term, so it cannot vanish.

The specialization therefore has no common nonzero
direction. The homogeneous resultant of \(P_1,P_2,P_3\)
is consequently a nonzero rational function of \(D\).
The multiplication by \(D-D_*\) used to specialize the
third form does not affect this conclusion over
\(\mathbb Q(D)\).

## 2. Returning to the ordinary mixed-axis moments

A mixed-side term of \(F_m\), with \(m=2\) or \(3\), has
coefficient
\[
  \binom m\alpha \binom{mD}{E_\alpha(D)}^{-1},
\]
where both \(E_\alpha(D)\) and \(mD-E_\alpha(D)\)
are at least \(D-mM\). At \(D\ge6M\) set
\(k=\min(E_\alpha,mD-E_\alpha)\ge\lfloor D/2\rfloor\).
Then \(\binom{mD}{k}\ge\binom{2k}{k}\ge2^k\).
Every coefficient of \(F_m-P_m\) is therefore bounded
by \(m!2^{-\lfloor D/2\rfloor}\). All coefficients of
both families are bounded by \(m!\), since the
within-side coefficient is also a multinomial
coefficient divided by a binomial coefficient at
ordinary degrees.

One can implement the resultant comparison using a
binary Sylvester determinant. Write the common first
form as \(c_1(D)X+c_2(D)Y+c_3(D)Z\) and parameterize
its zero plane by
\[
        X=-c_2 U-c_3 V,\qquad Y=c_1 U,\qquad Z=c_1 V .
\]
The parameters span the plane injectively both at
ordinary degrees and at \(D_*\), since \(c_1>0\) there.
The binary degree-two and degree-three resultant
formed from \(P_2,P_3\) after this substitution is a
nonzero rational function \(A(D)/B(D)\), by Section 1.
By the leading term of \(A\) and the coefficient
sum of \(B\), one can compute positive \(c\), an
integer \(j\ge0\), and a threshold such that
\[
                              |A(D)/B(D)|\ge c D^{-j}
\]
above the threshold. The same \(5\)-by-\(5\) Sylvester
determinant for \(F_2,F_3\) differs by at most a
computable constant times \(2^{-\lfloor D/2\rfloor}\).
This follows by expanding the determinant, using
the bounded coefficients and the mixed-term estimate.
For all sufficiently large \(D\), this difference
is less than \(cD^{-j}\). An integer cutoff is
computable by increasing \(D\) until this inequality
holds on both parities in the range where
\(D^j2^{-\lfloor D/2\rfloor}\) decreases separately on
each parity. The finitely many smaller degrees can
then be tested by their exact original forms, including
all coefficient directions over \(\mathbb C\). A
nonzero characteristic-zero direction would also
exist over \(\overline{\mathbb Q}\), so exclusions
over \(\mathbb C\) hold over every characteristic-zero
field.

The uniform interior theorem forces every exponent pair
in a failing profile to have its smaller entry below one
absolute cutoff. There are only finitely many assignments
of these bounded entries to the two sides. Applying the
per-profile conclusion to each gives the asserted absolute
degree bound.

The [exact profile checker](verify_sfc_homogeneous_trinomial_axis_degrees.py)
implements this comparison and tests every smaller admissible degree.
It uses an upper bound of \(500\) on the absolute Sylvester
entries after the linear substitution: each coefficient is bounded
by \(10\cdot6\cdot2^3\), using \(0<c_j(D)\le1\) and at most ten
terms before substitution. The entry differences have the same
bound times \(2^{-\lfloor D/2\rfloor}\). Expanding the
five-row determinant gives the error constant \(600\cdot500^5\).
For example

~~~sh
python3 research/factorial/image-rigidity/verify_sfc_homogeneous_trinomial_axis_degrees.py 0 1 1 --opposite-third
~~~

checks the tied-maximum profile with fixed exponents \(0,1\)
on one side and \(1\) on the other. Replays for fixed exponents
\(0,1,2\) on one side and for \(0,1\) opposite each of
\(0,1,2\) find no failed degrees.

## 3. Making the absolute bound effective

**Theorem.** An absolute upper bound on the degree of any
failed initial homogeneous bivariate trinomial window is
computable. Thus the unrestricted initial window for this
support class reduces to one finite exact calculation.
The algorithm below specifies a cutoff without evaluating it.

We make the interior cutoff effective. Use the notation
and the local moments of the
[uniform interior proof](SFC_UNIFORM_HOMOGENEOUS_TRINOMIAL_INITIAL_CORNER.md).
In that proof the following estimates can be used with
computable constants.

- The normalized local kernel on a bounded \(\eta\)-interval
  converges with every fixed finite number of derivatives
  at a computable rate as the smaller coordinate grows.
  For the Taylor estimates there, \(s(c)\ge\sqrt{B/2}\)
  if both coordinates are at least \(B\), and
  \(s(c)/c,s(c)/(D-c)\le B^{-1/2}\).
  The derivatives of \(D\sin^2(\eta/(2\sqrt D))\) of order
  at least two are bounded by one for \(D\ge1\).
  These give explicit bounds in the displayed Taylor
  estimates of the interior proof.
- For separated arguments the mixed derivative bound
  has computable \(c_0>0,C\) and, using the modified
  width at the boundary, decay \((1+T)^3e^{-c_0T^2+CT}\).
  One can use \(c_0=1/16\), since
  \(-2\log\cos z\ge z^2\) for \(|z|<\pi/2\)
  and at least one extreme is at \(\eta\)-distance \(T/2\)
  from the mean parameter. All derivative bounds
  used there are finite product and sum estimates
  involving orders at most six.
- A source of explicit gamma remainder bounds for these
  estimates is
  [Binet's Laplace integral](https://dlmf.nist.gov/5.9.E10_2).
  Its kernel
  \(\bigl(1/(e^t-1)-1/t+1/2\bigr)/t\) is bounded in
  absolute value by three on \(t>0\).
  For \(t\le1\), expand \(e^t-1=t+t^2/2+r\) with
  \(0\le r\le t^3\) and use \(e^t-1\ge t\).
  For \(t\ge1\) bound the three terms separately.
  Differentiating the Laplace integral gives the bound
  \(3k!/x^{k+1}\) for the \(k\)-th Stirling remainder
  derivative at \(x>0\).
- For bounded boundary entries the pure rescaled angular
  moments converge at a computable rate as \(D\to\infty\).
  For \(U=\sum_i a_i\) on the left boundary, the exact
  rescaled product integral is
  \[
    D K_{D,m}=
       \left(\prod_i (1-a_i/D)^{-(D-a_i)}a_i^{-a_i}\right)
       \frac{U!D^{U+1}}{(mD+1)_{\underline{U+1}}},
  \]
  with \(0^0=1\). For bounded integers \(a_i\), a
  geometric-series bound for the logarithm of each
  factor supplies a computable error bound. The
  right boundary is the same after reflection.

We detail the cutoff selection so that the unbounded
coefficient ratios and the angular and width thresholds
can be handled independently.

### A uniform margin in the separated interior model

Fix a bound \(L\) on the span of each local node cluster.
Allow all partitions of at most three nodes into clusters
and any corresponding consecutive subwindow of \(1,2,3\).
The local models are the Gaussian systems in the Newton
bases, including coincident nodes. Give each cluster a
positive scalar width \(\lambda_g\), and let its coefficient
vector have maximum norm \(N_g\). Normalize the summed
pure moments in each row by
\[
                           W(m)=\max_g \lambda_g N_g^m .
\]
At least one coefficient vector is nonzero. There is
a computable positive margin for these normalized models,
uniform even in the ratios of the widths. Here the
margin is a lower bound on the largest absolute moment
over the prescribed subwindow.

For one cluster this is compactness, using the Gaussian
exclusions of the interior proof. For several clusters,
induct on the number of clusters. Choose \(R>1\) so that
at a row where one group has weight less than \(1/R\)
of the other, its entire normalized pure contribution
is less than half the smaller-model margin.
Such a choice is computable, since each local pure form
has a uniform upper bound on unit coefficient vectors
in the compact node range. Sort the widths in decreasing
order. If a ratio between consecutive widths exceeds
\(R^8\), split at that gap. A group weight is the maximum
of its cluster weights. If the later group is within
\(1/R\) of the earlier at row \(m\), and cluster \(g\)
attains the earlier maximum at row \(m'>m\), compare
that \(g\) with a maximizer of the later group at row
\(m\). The weight ratio in row \(m'\) is bounded below
by \((\lambda_g/\lambda_{\rm late})^{m'/m-1}R^{-m'/m}>R\).
Here \(m'/m-1\ge1/2\) and \(m'/m\le3\).
Thus the active sets of the two groups are ordered.
The [assignment lemma](SFC_SEPARATED_DEGREE_CLUSTER_WINDOWS.md#1-the-ordered-assignment-of-moment-rows)
assigns one group its own full subwindow with contamination
less than half its margin.

If no width ratio is that large, normalize the maximum
width and maximum coefficient norm to one. All widths
belong to a compact interval bounded away from zero.
The separated Gaussian model on this compact set has
a nonzero moment by the singleton, pair and triple
exclusions of the interior proof. For a pair and a separate
singleton the subwindow can only be \(1,2,3\), and
the exclusion there permits every positive weight ratio.
The minimum margin is computably positive.
For example, subdivide the compact parameter boxes using
interval bounds for the exponential kernel and its divided
differences, and refine until every box has a positive
lower bound in at least one row. The simplex derivative
integral for the divided differences gives uniform
derivative bounds at coincident nodes as well.
The normalization of maximum coefficient norm one can
be parameterized by finitely many choices of a unit
coefficient, with other coordinates in the closed unit
disk. This yields the required margin by induction.

### Selecting the angle and exponent cutoffs

Fix a cap \(B\ge1\) for all columns treated as boundary
columns. Suppose there are at most two of these and that
each other column has both exponents at least \(C_B\),
to be chosen. For a prescribed local interior span \(L\),
use the computable interior margin just proved. For the
boundary part use the minimum of the positive margins
over its finitely many patterns \(u_j\le B\) at the
last one or two orders, respectively. Such a margin is
computable from the rational gamma forms of the interior
proof by compactness, and their positivity there is
proved for a same-side pair and for opposite singletons.

First treat the interior and boundary parts as two groups,
with weights \(W_{\rm in}(m)=\max_g h_g N_g^m\) and
\(W_{\rm bd}(m)=D^{-1}N_{\rm bd}^m\).
If the boundary part is absent, use the interior margin
on the full window. Otherwise choose \(R>1\) so that an
inactive group's normalized model contribution is less
than one quarter of the smaller of the two margins.
For \(C_B\) large enough, every interior width exceeds
\(D^{-1}R^8\), since \(D h(a)=s(a)\ge\sqrt{C_B/2}\).
The same ordered comparison as above now assigns one
group its full subwindow. The boundary group, if assigned,
gets the last one or two orders.

To ensure that the interior clusters have span at most a
chosen \(L\) while mixed terms are small, choose
\(L_0<L_1<L_2<L_3\) as follows. Take \(L_0\ge 8(B+1)\).
For each \(i=0,1,2\), use span bound \(2L_i\), compute
the preceding margins and ratio cutoff, and choose
\(L_{i+1}\) large enough to make all mixed errors smaller
than one quarter of the smaller margin whenever the
separation is at least \(L_{i+1}\). This is possible
computably by the exponential bound. When bounding
a kernel for this step, widths within an interior cluster
of span at most \(2L_i\) differ by at most a computable
factor from its reference width, by the log-width
derivative bound. Modified boundary widths are at most a
constant depending on \(B\) times \(D^{-1}\). After
inserting the coefficient sizes, the product
of the row weights to the proportions of selected factors
is bounded by their maximum, exactly as in the interior
proof. Thus the mixed-error threshold does not depend on
the large width ratios or the degree.

Among these three adjacent cutoff bands at least one
contains no gap between successive interior nodes, as
there are at most two such gaps. Group across the gaps
below that band. The spans within clusters are at most
twice the lower cutoff, and all separations across clusters
are at least the upper cutoff. Boundary columns on the
same side stay together. Increase \(C_B\) to ensure
the required separations to the boundary and between its
sides: the boundary angular span on one side is bounded
by \(4\sqrt{B+1}\) for \(D\ge2(B+1)\), whereas interior
nodes have distance at least \(2\sqrt{C_B}\) from
either endpoint of the angular interval. Also increase
\(C_B\) to meet the width-ratio cutoffs and to make the
local convergence errors less than one quarter of the
chosen margins, for all three possible bands.
All choices are effective by the kernel estimates above.
A selected group then retains a nonzero moment in its
assigned window. This proves the interior exclusion
whenever each column is either bounded by \(B\) in its
smaller entry or has both entries at least \(C_B\),
with at least one column of the latter kind.

Finally take \(B_0=1\). For \(i=0,1,2\) set
\(B_{i+1}>\max(B_i,C_{B_i})\). If some column has both
entries at least \(B_3\), at most two smaller entries
lie below \(B_3\). One of the three half-open intervals
\([B_i,B_{i+1})\) therefore contains no smaller entry.
Apply the preceding criterion with \(B=B_i\) and
\(C_B\le B_{i+1}\). Thus \(B_3\) is a computable uniform
interior cutoff. Only the finitely many axis profiles
with \(u_j<B_3\) remain. The maximum of their
computable cutoffs from Section 2 is a computable
absolute bound on the degree of a failed window.
