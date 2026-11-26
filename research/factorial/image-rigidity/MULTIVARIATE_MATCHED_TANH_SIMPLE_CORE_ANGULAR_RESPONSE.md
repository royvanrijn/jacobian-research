# Matched tanh angular response at simple roots of an arbitrary core

Fix a finite ordinary two-monomial exponent profile with distinct columns,
called the core. Denote its raw order-\(m\) moment polynomial by \(Q_m\).
Append any nonempty finite collection of unequal shared rows. Their
orientations and integer total degree difference are unrestricted.

For an appended row, let \(\ell\ge1\) be its smaller exponent,
\(d\ge1\) its gap, and \(\epsilon\in\{-1,1\}\) its orientation. In the
[matched tanh grid](MULTIVARIATE_MATCHED_TANH_INTEGER_MULTIPLICITY.md#recovering-the-gauss-grid)
the slot \(1\le h\le d\) has coordinates
\[
 v=\left(\frac{n}{T_{n+1}}-\frac{n-1}{T_n}\right)^{-1}
       =\frac{d}{2\ell+d},\qquad
 c=\frac{2h}{d}-1.
\]
The value \(v\) is the same across the row and independent of the
comparison order \(n\). Define three exact tail statistics, summing
over slots,
\[
 W=2\sum v,\qquad P=\sum\epsilon v^2,\qquad
 E=4\sum c v^2=4\sum_{\rm rows}v^2 .
 \tag{1}
\]
The last equality uses \(\sum_{h=1}^d (2h/d-1)=1\).

There is a computable all-degree small-tail angular test on the
exact core. It applies even if the core has coincident arguments,
provided the roots involved are simple and the first or second angular
response specified below separates each coincidence. The tolerance is
uniform over tails of sufficiently small \(W>0\), including tails
with just one added row. The calculation uses both signed and
unsigned slot statistics and retains the restrictions on their joint
contributions from complete integer rows.

## The matched coefficient expansion

For each degree \(m\ge2\), strip the tail's positive constant
coefficient and positively rescale the variable to normalize both
tail endpoints to one. Dividing the resulting coefficient by the
corresponding coefficient of \(Q_m\) gives the normalized factorial
multiplier of the tail. Its logarithm at index \(k\) is
\[
 -\frac{W}{2m}k(m-k)
 +\frac{2P}{3m}\left(\frac{k}{m}-\frac12\right)k(m-k)
 +\frac{E}{4m^2}k(m-k)
 +O(mW^3).
 \tag{2}
\]
The error is uniform in the number, gaps, and orientations of the
appended rows.

To see this, put \(M=\ell+d/2\) on one row and \(s=k/m\). Write
\(x=M+\epsilon d(s-1/2)\) and let \(L_m(s)\) be the sum of the
endpoint-normalized log-Gamma chords of the tail per unit order.
Its second derivative is a sum of \(d^2 m\psi'(1+mx)\).
By the
[trigamma integral and midpoint bound](MULTIVARIATE_HOLOMORPHIC_MATCHED_LOG_TRANSPORT.md#the-exact-curvature-kernel),
\[
 m\psi'(1+mx)=\int_0^\infty e^{-xu}\frac{u/m}{e^{u/m}-1}\,du
             =\frac1x-\frac1{2m x^2}+R,
 \qquad 0\le R\le\frac1{6m^2x^3}.
\]
Indeed \(0\le t/(e^t-1)-1+t/2\le t^2/12\), which follows by
integrating the cited bound on the derivative. Let \(b=\max v\) over
the tail and assume \(b\le1/4\). With \(y=2\epsilon v(s-1/2)\),
one has \(x=M(1+y)\) and \(|y|\le v\). Replacing \((1+y)^{-1}\)
by \(1-y\) costs at most \(v^2/(1-b)\). Replacing \((1+y)^{-2}\)
by \(1\) costs at most \(2v/(1-b)^3\). Also \(1/M=2v/d\le2v\).
After multiplying by \(d^2\), the combined error in this row is at
most \(5(d^2/M)v^2\), since \(m\ge2\) and
\[
 \frac1{1-b}+\frac{2}{m(1-b)^3}
                 +\frac{2}{3m^2(1-b)^3}\le5 .
\]
The retained terms sum to
\[
 L_m''(s)=W-4P(s-1/2)-E/(2m)+\mathcal E(s),
 \qquad |\mathcal E(s)|\le5Wb^2 .
\]
Here \(d^2/M=2dv\), \(d^2/M^2=4v^2\), and the linear term
of \(d^2/x\) is \(-4\epsilon d v^2(s-1/2)\).
Each term is integrated with zero endpoint values. The Green
kernel has absolute integral \(s(1-s)/2\), and the function
\(-(s-1/2)s(1-s)/6\) has second derivative \(s-1/2\).
Multiplying the resulting \(L_m\) by \(m\) gives (2), with
absolute remainder at most \(5mWb^2/8\). The bound \(b\le W/2\)
follows from (1) and makes the error cubic in \(W\).

## The simple-root response

Use \(D=z\,d/dz\), and define three diagonal operators on coefficients
at degree \(m\):
\[
 A_m=-\frac{D(m-D)}{2m},\qquad
 C_m=\frac{D(m-D)}{4m^2},\qquad
 B_m=\frac{2}{3m}\left(\frac Dm-\frac12\right)D(m-D).
\]
For one simple core root \(\zeta\), write \(G=(DQ_m)(\zeta)\ne0\).
All operator expressions in the following definitions are evaluated
at \(\zeta\):
\[
 \begin{aligned}
 \alpha&=-A_mQ_m/G,&
 p&=-B_mQ_m/G,&
 q&=-C_mQ_m/G=-\alpha/(2m),\\
 \beta&=-\bigl(\tfrac12 A_m^2 Q_m+
                 \alpha D A_m Q_m+
                 \tfrac12\alpha^2 D^2Q_m\bigr)/G .
 \end{aligned}
 \tag{3}
\]
If \(\widetilde\zeta\) is the nearby root after the stated tail
normalization, the branch of the complex logarithm near one gives
\[
 \log(\widetilde\zeta/\zeta)
        =W\alpha+W^2\beta+Pp+Eq+O_{Q_m,m}(W^3).
 \tag{4}
\]
The constant is uniform over the tails.

In fact \(|P|\le W^2/4\) and \(0<E\le W^2\).
Exponentiating (2) gives the polynomial
\[
 Q_m+\bigl(W A_m+\tfrac12W^2 A_m^2+P B_m+E C_m\bigr)Q_m
                        +O_{Q_m,m}(W^3)
\]
coefficientwise. On a fixed small disk about \(\zeta\), Rouché's
theorem gives the unique perturbed root for small \(W\), and local
division by \(Q_m(z)/(z-\zeta)\) gives displacement \(O(W)\).
Taylor expand in \(w=\log(z/\zeta)\) at \(w=0\). The coefficient of
\(w\) from \(Q_m(\zeta e^w)\) is \(G\), its coefficient of \(w^2\)
is \((D^2Q_m)(\zeta)/2\), and the coefficient of \(w\) from
\((A_m Q_m)(\zeta e^w)\) is \((D A_m Q_m)(\zeta)\).
The degree-one comparison first gives \(w=W\alpha+O(W^2)\).
Substitute this into the terms already quadratic and divide by \(G\)
to obtain (4). The angular response is the imaginary part of (4),
since the discarded tail normalizations are positive.

## A sufficient test in every degree

Take an integer \(J\ge4\) at least the
[effective absolute ordinary angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md).
For each \(2\le n<N=n+1<J\), consider every pair
\((\zeta_n,\zeta_N)\) of core roots with equal argument. Assume the
two roots in each such pair are simple. Use (3) at each, and write
\(\Delta\alpha=\alpha_N-\alpha_n\), and likewise for the other response
coefficients. Each pair must satisfy at least one of the following
tests.

- \(\Im\Delta\alpha\ne0\).
- \(\Im\Delta\alpha=0\) and the real quadratic
  \[
       H(u)=\Im\Delta\beta+\tfrac14 u\,\Im\Delta p
                                  +u^2\Im\Delta q
  \]
  has the same strict sign at \(u=0\) and at every \(u=1/d,-1/d\)
  for positive integers \(d\).

The first test can also be evaluated as
\[
 \Im\Delta\alpha
   =\frac1{2nN}\Im\left[
       N\,\frac{D^2Q_n}{DQ_n}(\zeta_n)
       -n\,\frac{D^2Q_N}{DQ_N}(\zeta_N)\right].
\]
The second test follows from a joint integer restriction. A row
of gap \(d_i\) and rate \(v_i\) contributes weight
\(\lambda_i=(2d_iv_i/W)^2\), and \(\sum\lambda_i\le1\). The
normalized pair of quadratic statistics therefore has the form
\[
  (P/W^2,E/W^2)=\sum_i\lambda_i
                             \left(\frac{\epsilon_i}{4d_i},\frac1{d_i^2}\right).
\]
The remaining weight \(1-\sum_i\lambda_i\) can be placed at \((0,0)\).
Thus the second-order angular coefficient in the difference of (4)
is a convex combination of the tested values \(H(u)\).
They have a positive minimum in absolute value, since the only
limit point of the test set is \(0\). The cubic error in (4) is
smaller than this second-order margin for small \(W>0\).

The tests are effective. The core roots and all coefficients in (3)
are algebraic, and equality of two rays is decidable by real algebraic
arithmetic. If \(H(0)\ne0\), for all sufficiently large \(d\) the
explicit upper bound
\[
       |H(\pm1/d)-H(0)|
              \le |\Im\Delta p|/(4d)+|\Im\Delta q|/d^2
\]
is smaller than \(|H(0)|/2\), leaving only finitely many sign checks.
The coefficient and root Taylor bounds can be made rational by
algebraic root isolation and Rouché disks. Pairs with distinct core
arguments stay separated by continuity, including pairs involving a
multiple root. There are finitely many such comparisons below \(J\).
The angular cutoff handles \(N\ge J\) for every appended ordinary
profile. The pair of orders one and two is separated because the
first-order root is real and the second factorial moment is the
expectation of an even power of a nonzero real polynomial on
independent unit exponentials. This is strictly positive on the
real axis. These observations give a computable threshold \(W_0>0\)
such that every appended tail with \(0<W<W_0\) has disjoint
consecutive arguments and hence satisfies SFC(2) in every degree.

For a real reciprocal core after positive endpoint normalization,
every simple root on the unit circle has \(p\) real. In fact \(B_mQ_m\)
is anti-reciprocal, since reciprocity interchanges \(D\) with \(m-D\).
For the perturbation \(Q_m+t B_m Q_m\), real \(t\), the root branch
at the opposite parameter is \(1/\overline{\zeta(t)}\). Differentiating
the logarithmic displacement relation \(w(-t)=-\overline{w(t)}\)
gives the assertion. Consequently, when the common first angular
speed is \(a=\Im\alpha_n=\Im\alpha_N\), the second test for two such
roots reduces to checking that the endpoints
\(\Im\Delta\beta\) and \(\Im\Delta\beta+a/(2nN)\) have the same
strict sign. This simplified test still permits arbitrary asymmetry
in the appended rows.

## Sharpness of the second-order test

Fix one simple core ray coincidence for which the first angular
response vanishes. The strict common sign condition is equivalent
to a uniform positive lower bound of order \(W^2\) for the angular
splitting of these two roots over all sufficiently small ordinary
tails.

For necessity, fix a gap \(d\ge1\), an orientation \(\epsilon\), and
a desired mass fraction \(r\in[0,1]\). For each positive integer
\(L\), if \(r>0\) take one row of this gap and orientation at
smaller base \(\lceil d^2 L/r\rceil\). Add
\(\lfloor(1-r)L\rfloor\) unit-gap rows of smaller base \(L^2\),
with arbitrary orientations. For \(r=0\) use just the unit rows.
The total \(W\) is asymptotic to \(1/L\). The distinguished row,
when present, has asymptotic mass fraction \(r\). The unit rows'
total contributions to both quadratic statistics are \(O(L^{-3})\).
Therefore
\[
 (P/W^2,E/W^2)\ \longrightarrow\
        r^2\left(\frac{\epsilon}{4d},\frac1{d^2}\right).
\]
This realizes every segment from the origin to one of the joint
integer test points as a normalized-moment limit.

If \(H(0)=0\) the diffuse limit already gives splitting
\(o(W^2)\). Otherwise failure of the strict test supplies some
\(H(\epsilon/d)\) that is zero or has sign opposite to \(H(0)\).
Choose \(r^2\) so the linear combination
\((1-r^2)H(0)+r^2H(\epsilon/d)\) vanishes. The displayed
ordinary sequence then has splitting \(o(W^2)\) by (4). This
limits a quadratic-order angular margin and makes no assertion
of an exact zero at positive tail strength.

The theorem is a local sufficient condition. It leaves open the
cases of a multiple coincident core root and of vanishing or
opposing second-order responses. It does not assert an angular
gap for an arbitrary unmodified ordinary profile.
