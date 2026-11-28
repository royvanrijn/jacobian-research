# A PSD face with exact mixed-radius comparison

Let \(F\) satisfy the range, balance and triangle conditions (AT1)
of the three-label agreement relaxation, and suppose its integral
operator on \(L^2(\gamma_2)\) is positive semidefinite. We may
average under common rotations and reflections. Write \(F(r,s,t)\)
for the resulting even angular kernel, \(t\in[-\pi,\pi]\) modulo
\(2\pi\), and let \(d\mu(r)=r e^{-r^2/2}\,dr\).
For the comparison below there is a canonical same-radius
representative \(F(r,r,t)\), defined by its angular Fourier
Gram coefficients. Values of an arbitrary measurable version
on the radial diagonal play no role.
[The circle certificate and the mixed-radius counterexample](GAUSSIAN_AGREEMENT_TRIANGLE_RELAXATION.md)

The simplex comparison is exact at every \(0<\rho<1\) on the
face where the angular mean of \(F(r,s,\cdot)\) is \(1/3\) for
almost every radial pair. It is also exact if
\(F(r,r,2\pi/3)=0\) almost everywhere in the representative
above. Equality in either case gives the simplex kernel
almost everywhere. The following decomposition identifies
the term left over for a general PSD kernel.

## Angular Gram coefficients and same-radius triangle inequalities

Use angular measure \(dt/(2\pi)\). Put
\(K_m(r,s)=\widehat F(r,s,m)\) for \(m\in\mathbb Z\).
Testing the PSD kernel on products of a radial function and
the \(m\)-th angular mode gives a positive semidefinite
radial integral operator \(K_m\). The kernels are real and
\(K_{-m}=K_m\) after reflection averaging. They admit vectors
\(v_m(r)\) in real Hilbert spaces with
\[
 K_m(r,s)=\langle v_m(r),v_m(s)\rangle,\qquad
                \sum_m\|v_m(r)\|^2\le1\quad\text{a.e.}         \tag{AF1}
\]
Here one may take \(v_{-m}=v_m\). To justify (AF1), average
the radial variables independently on cells of a finite
radial partition \(\mathcal P\). For every cell \(i\), the
same-cell angle function \(F_{ii}\) has nonnegative Fourier
coefficients by PSD and has values between zero and one
almost everywhere. Fejér averaging at zero shows that the
sum of its Fourier coefficients is at most one. Weight by
the cell masses and sum. It follows that
\(\sum_m\operatorname{Tr}(P_\mathcal P K_m P_\mathcal P)\le1\),
where \(P_\mathcal P\) is radial conditional expectation.
On a generating refinement the trace bound passes to
\(\sum_m\operatorname{Tr}K_m\le1\). The spectral
representations of these positive trace-class operators
give the Gram vectors. Their cell averages have squared
norm sum at most one by the preceding finite bound.
Convergence of conditional expectations in the direct-sum
Hilbert space now gives the pointwise norm bound, using
almost-everywhere convergence along a subsequence.

For almost all \(r,s\), the continuous function
\(\sum_m K_m(r,s)e^{imt}\) agrees almost everywhere in angle
with the original kernel. Absolute uniform convergence in
angle follows from Cauchy–Schwarz and (AF1). Define
\(F_r(t)=\sum_m\|v_m(r)\|^2e^{imt}\). It obeys the circular
range and triangle inequalities. Indeed each averaged
same-cell angle function obeys them, by taking three
independent radii from that cell in the original triangle
inequalities. The continuous angular extensions preserve
the inequalities, and convergence of the cell-averaged
Gram vectors at \(r\) gives uniform convergence of their
same-cell angle functions to \(F_r\). This proves the
assertion without assuming regularity of the radial
variables.

Write
\[
 p(r)=\|v_0(r)\|^2,\quad v_*=\mathbb E_\mu v_0,\quad
 \Delta(r)=p(r)-\tfrac13,\quad
 H(r,s)=\langle v_0(r)-v_*,v_0(s)-v_*\rangle .
\]
Balance gives \(\langle v_0(r),v_*\rangle=\|v_*\|^2=1/3\),
so \(\Delta=\|v_0-v_*\|^2\ge0\), \(H(r,s)=K_0(r,s)-1/3\),
and \(H\) has zero row mean. In particular if
\(\int\Delta\,d\mu=0\) then \(K_0(r,s)=1/3\) almost everywhere.
Putting \(a=2\pi/3\) and summing over the three equal
angles yields
\[
 1+2F_r(a)\ \ge\ F_r(0)+2F_r(a)
           =3\!\sum_{m:\,3\mid m}\|v_m(r)\|^2
           \ \ge\ 3p(r).
\]
Thus
\[
                  0\le\Delta(r)\le \tfrac23 F_r(a).            \tag{AF2}
\]

## Mixed-radius comparison and its remainder

Fix \(0<\rho<1\) and, for each \(r,s\), put
\(k_{rs}(t)=M_\rho(r,s,t)\), the relative Mehler density
in radial and angular coordinates. Its Fourier coefficients
\(\widehat k_{rs}(m)\) are positive when \(r,s>0\), as follows
by expanding \(e^{\kappa\cos t}\) into the product of the two
exponential series in \(e^{it}\), \(e^{-it}\), with
\(\kappa=\rho rs/(1-\rho^2)>0\). Write \(c_{rs}\) for \(c_k\)
in (AT3) with this angular weight. Let
\(\mathcal D_r(k_{rs})\) be \(1/\pi\) times the right-hand
side of (AT3) with \(F_r\) in place of \(F\), even when
\(F_r\) has mean different from \(1/3\).
It is nonnegative, and the same coefficient calculation
gives the arbitrary-mean identity
\[
 \int_{\mathbb T} k_{rs}(t)(F_T-F_r)(t)\,\frac{dt}{2\pi}
          =-c_{rs}\Delta(r)+\mathcal D_r(k_{rs}).
\]
The Gram formula also gives, for every mode,
\[
  \frac{\widehat F_r(m)+\widehat F_s(m)}2-K_m(r,s)
                   =\frac12\|v_m(r)-v_m(s)\|^2.
\]
Combining the two equalities and integrating over the
radii gives the exact comparison
\[
\begin{split}
 \mathcal S_\rho(F)-S_\rho(T)
     ={}& \iint \widehat k_{rs}(0) H(r,s)\,d\mu^2\\
        &+\iint [c_{rs}-\widehat k_{rs}(0)]
                        \frac{\Delta(r)+\Delta(s)}2\,d\mu^2\\
        &-\frac12\iint
                       [\mathcal D_r(k_{rs})+\mathcal D_s(k_{rs})]\,d\mu^2\\
        &-\frac12\sum_{m\ne0}\iint\widehat k_{rs}(m)
                          \|v_m(r)-v_m(s)\|^2\,d\mu^2 .       \tag{AF3}
\end{split}
\]
All sums converge by (AF1), the inequality
\(|\widehat k_{rs}(m)|\le\widehat k_{rs}(0)\), and the unit
row integral of the radial Mehler kernel
\(\widehat k_{rs}(0)\). In particular the last two terms
are nonpositive. The zeroth-frequency gain is fully
displayed in the first two terms.

If \(\int \Delta\,d\mu=0\), then \(H=0\) and (AF3)
proves the simplex upper bound. By (AF2) the angle-\(a\)
zero condition suffices. For equality, the positive
diagonal triangle certificate with its exact circular
mean forces \(F_r=F_T\) for almost every \(r\).
The last term in (AF3) then makes every \(v_m\) for
\(m\ne0\) independent of the radius almost everywhere,
since the Mehler coefficients are strictly positive.
The zero-frequency kernel is already constant, so
\(F(r,s,\cdot)=F_T\) almost everywhere. If \(F\) is the
common-rotation average of a genuine balanced agreement,
this equality determines the underlying partition.
[Recovering the ordinary partition](GAUSSIAN_AGREEMENT_TRIANGLE_DISCRETIZATION.md#recovering-a-partition-from-its-averaged-agreement)

There is also a dimension-specific uniform defect bound
for the planar radial measure:
\[
       \mathcal S_\rho(F)-S_\rho(T)
       \le (\rho^2+\tfrac12)\int\Delta(r)\,d\mu(r)
       \le \frac{1+2\rho^2}{3}\int F_r(a)\,d\mu(r).           \tag{AF4}
\]
For the first term in (AF3), the radial mean-zero
functions \(v_0-v_*\) have no first Gaussian chaos,
so the radial OU operator has norm at most \(\rho^2\)
on them. For the second term,
\(c_{rs}\le(\pi/a)\widehat k_{rs}(0)
=\tfrac32\widehat k_{rs}(0)\). In fact \(c_{rs}\) averages
the decreasing weight \(k_{rs}\) against the density
\(2t/a^2\) on \([0,a]\). The cumulative distribution
of this density is bounded by \(\pi/a\) times the
uniform cumulative distribution on \([0,\pi]\).
Integration against a nonnegative decreasing weight
gives the comparison. Integrate the radial unit-row
identity and use (AF2). The planar radial setting is
sufficient here because ordinary three-cell maximizers
admit dimension reduction to a plane. For an ordinary
partition with circular label masses \(p_i(r)\), one has
\(\Delta(r)=\sum_i(p_i(r)-1/3)^2\) and \(F_r(a)\) is the
angular probability of agreement after a rotation by \(a\).
These identities follow by using the ordinary circular
indicator coefficients as Gram vectors in (AF1).

The general PSD comparison requires an estimate for the
positive radial contribution in (AF3). Neither the
circle certificate nor the angle-\(a\) PSD packing
forces \(F_r(a)=0\) for arbitrary feasible kernels.

## Four-point constraints on mixed angular means

Suppose in addition that for almost every four observations
their six agreements sum to at least one. Fix three points
on radius \(r\) with angles \(0,a,2a\), and let the fourth
point have radius \(s\) and angle \(t\). Then
\[
      \sum_{j=0}^{2} F(r,s,t+ja)\ \ge\ 1-3F_r(a).                \tag{AF5}
\]
For a measurable kernel, take the first three radii
independently from the radial averaging cell containing
\(r\), and the fourth from the cell containing \(s\).
The averaged inequality holds almost everywhere in the
four angles and then everywhere by continuity of the
Fourier representatives. Use the convergence in (AF1)
to take both radial cell sizes to zero.

Integrating (AF5) in angle and interchanging the roles of
the two radii gives
\[
 H(r,s)\ge-\min\{F_r(a),F_s(a)\},\qquad
       \int |H(r,s)|\,d\mu(s)\le 2F_r(a).
\]
The last inequality uses the zero row mean of \(H\).
The pointwise lower bound in (AF5) also yields
\[
 \int\!\!\int_{\mathbb T}
    \left|\sum_{j=0}^{2}F(r,s,t+ja)-1\right|
                  \frac{dt}{2\pi}\,d\mu(s)\le 6F_r(a)
\]
because the expression inside the absolute value has
zero integral by balance. Thus a small same-radius overlap
at \(120^\circ\) controls both the radial mean defect and
the full error in the three-angle sum on each row. These
bounds still allow a positive variance of the radial
means in (AF3).

A three-facet perturbation gives a strict Mehler gain near
correlation one even with every four-point label test. The
full PSD relaxation therefore permits positive total excess
in (AF3). An application to all correlations needs stronger
consistency or restrictions on the agreement of an ordinary
partition.
[The strengthened relaxation counterexample](GAUSSIAN_AGREEMENT_MULTI_FACET_COUNTEREXAMPLE.md)
