# Sign changes required for exact factorial replacement contacts

In a sequence of ordinary profiles with bounded curvature, let
\(C\) be the countable retained core and let \(E^+,E^-\) be the
added and omitted multisets as in the
[signed replacement equations](MULTIVARIATE_SIGNED_CORE_REPLACEMENT_CONTACTS.md).
Consider a positive Gaussian reference strength \(W_*>0\) at
adjacent orders \(n,N=n+1\), \(n\ge2\), with the same fixed
integer degree difference as the ordinary profiles. Choose
\(r\) distinct nonreal reference root pairs in the upper
half-plane, simple in both orders.

High-order cancellation can make the signed replacement much
smaller than the absolute second-order errors of its two parts.
A bound on sign changes of the signed factorial slots controls
the exact error by the *signed* low moments instead. The control
holds when slot heights coalesce. It gives a rank test excluding
exact ordinary common roots continuing the selected pairs. No
row count or gap bound is needed for the slot test. A fixed-gap
version uses fewer clock directions and also gives a scalar test
when earlier moments match exactly.

## A signed interpolation bound

Let \(\mu\) be a finite real signed atomic measure supported on
\([-\eta,\eta]\), with no mass at zero and with zero the only
possible accumulation point of its atoms. Merge weights at the
same point and omit zero weights. A sign change is a transition
between consecutive sign blocks in the ordered nonzero weights.
The transition can be at zero if both sides accumulate there.
Suppose there are at most \(s\) sign changes and write
\(M_i=\int t^i\,d\mu(t)\).

If a vector-valued function \(A\) is analytic on a fixed disk
about zero, with Taylor coefficients \(A_i=[t^i]A\), then
\[
 \int A(t)\,d\mu(t)
  =\sum_{i=0}^s A_i M_i
                 +O_A\!\left(\sum_{i=0}^s\eta^{s+1-i}|M_i|\right).
                                                               \tag{1}
\]
The constant also depends on \(s\). In particular the error is
\(O_A(\eta\sum_{i=0}^s|M_i|)\). If all the moments in (1)
vanish, then \(\mu=0\).

To prove this, choose one knot between each pair of adjacent sign
blocks, with zero allowed at an accumulation boundary. Let \(Q\)
be the monic polynomial with these \(q\le s\) knots as roots.
Choose knots outside the atoms, so \(Q\mu\) has one sign and
vanishes only if \(\mu=0\). In the variable \(x=t/\eta\), the
scaled knots and the support lie in \([-1,1]\). Write
\(\widehat Q(x)=\eta^{-q}Q(\eta x)\) and let \(\widehat\mu\)
be the pushforward of \(\mu\). Its moments are
\(\widehat M_i=M_i/\eta^i\). The coefficients of \(\widehat Q\)
are bounded by \(2^q\), so
\[
       \int|\widehat Q|\,d|\widehat\mu|
             =\left|\int\widehat Q\,d\widehat\mu\right|
             \le 2^q\sum_{i=0}^q|\widehat M_i|.
\]
Apply interpolation at the scaled knots to
\(H(x)=A(\eta x)-\sum_{i=0}^s A_i\eta^i x^i\).
On \(|x|=2\) its norm is \(O_A(\eta^{s+1})\).
The degree-\((q-1)\) interpolation polynomial has coefficients
of that size, and the remaining part is \(\widehat Q(x)\) times
a function bounded by the same quantity on \([-1,1]\).
Indeed the divided differences in Newton's formula are Cauchy
integrals on \(|x|=2\), whose denominators have all factors
at distance at least one from that circle. This bound
requires no knot separation. If \(q=0\), use the zero
interpolation polynomial and \(\widehat Q=1\).
Integrating the two parts gives (1). Zero moments force
\(\int\widehat Q\,d\widehat\mu=0\), proving the last assertion.

## Exact signed slots and two Gaussian clocks

For an unequal row let \(d=b-a\), \(D=|d|\), \(m=(a+b)/2\),
\(\epsilon=\operatorname{sgn}(d)\), and \(\ell=\min(a,b)\).
At order \(j\) its Gauss slot \(h\), \(1\le h\le D\), has
base \(b_h=(j\ell+h)/D\). Define the signed inverse slot
\[
 t_h=\frac{\epsilon}{b_h+(j-1)/2}
              =\frac{d}{jm+h-D/2}.
\]
Use sign \(+1\) for the occurrences in \(E^+\) and \(-1\)
for those in \(E^-\). In order \(j\) put
\[
    \mu_j=\sum_{\text{differing slots}}\!
            (\text{occurrence sign})\,\epsilon\,t_h^2\delta_{t_h},
    \qquad W_j=j\int \frac1t\,d\mu_j(t).
\]
Every sum with this integrand converges absolutely. In fact for
\(\eta=\sup_{\text{differing rows}}|d|/m\to0\) the slot formula
gives \(|t_h|\le \eta/(j-\eta/2)\). At each row,
\[
      j\sum_{h=1}^D|t_h|=\frac{d^2}m+O_j(\eta d^2/m),
      \qquad \sum_{h=1}^D|t_h|^2=O_j(\eta d^2/m).
\]
The total unsigned midpoint curvature in the differences is
bounded. Its removed part tends to zero, and its added part
tends to \(W_*\). Hence \(W_j\to W_*\). At any stage,
only finitely many differing slots have \(|t|\) bounded away
from zero: such rows have \(|d|/m\) bounded below and each
costs at least that amount of curvature.

Let \(x_i=i-(j-1)/2\) for the step \(i\to i+1\). Up to a term
constant in \(i\), the logarithm of the binomial-stripped
successive coefficient ratio contributed by a slot is
\(\epsilon\log(1+x_i t_h)\). This is the exact factorial
product of the \(D\) factors between adjacent arguments,
with its reciprocal for a decreasing row. Subtract the
average of the \(j\) step values and sum the first \(k\) of
them to obtain the normalized log coefficient. For \(p\ge1\)
write
\[
 S_{j,p}(k)=\sum_{i=0}^{k-1}x_i^p
                                  -\frac{k}{j}\sum_{i=0}^{j-1}x_i^p.
\]
Here \(S_{j,1}(k)=jF_j(k)\) with \(F_j(k)=-k(j-k)/(2j)\).
The exact normalized log multiplier of the signed difference
is therefore
\[
       W_jF_j(k)+\int A_j(t)_k\,d\mu_j(t),\qquad
 A_j(t)_k=\sum_{h\ge0}
              \frac{(-1)^{h+3}}{h+2}S_{j,h+2}(k)t^h.        \tag{2}
\]
The series is analytic on a fixed disk for the chosen orders.
The absolute bounds above justify all uses of it for the
infinite omitted multiset. The physical comparison still
uses the actual ratio of the endpoint argument scales. Its
logarithm tends to the reference center by compactification.

## Contact rank tests

At the chosen reference, let \(\mathscr L(U_n,U_N,c)\in\mathbb C^r\)
be the first derivative of the selected logarithmic root ratios
(upper minus lower), after perturbing the normalized log
coefficient arrays by \(U_n,U_N\), and *subtracting* the
center perturbation \(c\mathbf1\). Explicitly the variation
of the log root of \(G_j\) at \(\zeta\) is
\[
                 -\frac{\sum_k [z^k]G_j(z)\,U_{j,k}\zeta^k}
                           {\zeta G_j'(\zeta)}.
\]
Every direction and span below is over \(\mathbb R\).
Put
\[
 a_n=\mathscr L(F_n,0,0),\qquad
 a_N=\mathscr L(0,F_N,0),\qquad e=\mathbf1/(2nN).
\]
For \(h\ge0\) let \(R_{j,h}\) be \(\mathscr L\) applied to
\([t^h]A_j\) in its own order, with the remaining two inputs
zero.

**Slot test.** Fix two integers \(s_n,s_N\ge0\). If the
following list of vectors is linearly independent,
\[
            a_n,\ a_N,\ e,\quad
            (R_{n,h})_{0\le h\le s_n},\
            (R_{N,h})_{0\le h\le s_N},                           \tag{3}
\]
then an exact sequence of ordinary common roots continuing all
the selected pairs must, at every sufficiently late stage,
have more than \(s_j\) sign changes in at least one of the
two \(\mu_j\). This applies even with unbounded signed gaps
and arbitrary row counts.

For the proof, suppose both budgets hold on a subsequence.
Set \(X_j=W_j-W_*\), and let \(c\) be the difference of the
actual center from its reference value. These tend to zero.
Write \(M_{j,h}\) for the moments in (1).
Use (1)–(2) to expand the two signed profile perturbations.
The selected root branches are analytic in these perturbations
because the roots are simple. With
\[
 A=|X_n|+|X_N|+|c|
                     +\sum_{j=n,N}\sum_{h=0}^{s_j}|M_{j,h}|,
\]
the joint physical root equations are
\[
   0=a_nX_n+a_NX_N-c\mathbf1+
       \sum_{j=n,N}\sum_{h=0}^{s_j}R_{j,h}M_{j,h}
                      +O(\eta A+A^2).
\]
The linear part has a fixed real left inverse by (3).
Since \(A\to0\), it gives \(A=0\) for sufficiently late
stages. The interpolation bound then forces \(\mu_n=\mu_N=0\),
hence \(W_n=W_N=0\). This contradicts \(X_n=X_N=0\) and
\(W_*>0\).

There is a different chart when the differing rows have a
fixed signed gap \(d\). Let \(u=1/m\) and let
\(\omega_d=\sum(\text{occurrence sign})u^2\delta_u\)
be the signed midpoint measure. To define its exact jet,
write \(J_{j,k}(u,d)\) for the log factorial chord of a row
and \(E(u,d)\) for its center correction, with the signed
\(d\log(N/n)\) removed. The vector
\[
 A_d(u)=u^{-2}\bigl(
       (J_{n,k}(u,d)-d^2u F_n(k))_k,\
       (J_{N,k}(u,d)-d^2u F_N(k))_k,\
                              E(u,d)+du/(2nN)\bigr)             \tag{4}
\]
is analytic at zero for fixed \(d\). To see this directly,
a step ratio is the product of the \(D\) factors
\(jm+d(k-j/2)+h\) if \(d>0\). If \(d<0\) it is the
reciprocal of the product of the \(D\) factors
\(jm+d(k+1-j/2)+h\). In both cases its log minus
\(d\log(jm)\) is a finite sum of analytic logarithms in \(u\).
The linear coefficient in \(u\) is
\(d^2(k-(j-1)/2)/j+d/(2j)\).
The normalized chord keeps the first term and removes
the second. The step averages in the order difference keep
the center correction \(-du/(2nN)\). This proves (4),
including its linear subtractions.

Define \(R_{d,h}=\mathscr L([u^h]A_d)\) and \(a=a_n+a_N\).
The uses of the midpoint chart are as follows.

- For a single fixed signed gap \(d\), the bias deviation is
  exactly \(1/d\) times the midpoint clock deviation \(X\):
  \(\lambda_* = W_*/d\) and
  \(\delta=\lambda_+-\lambda_- -\lambda_*=X/d\).
  If \(a+e/d,R_{d,0},\ldots,R_{d,s}\) are independent, every
  exact contact sequence as above must eventually exceed the
  sign change budget \(s\) for \(\omega_d\).
- For a finite menu of signed gaps \(\mathcal D\) with budgets
  \(s_d\), the same exclusion holds if \(a,e\), together with
  all \(R_{d,h}\) for \(d\in\mathcal D,\ 0\le h\le s_d\),
  are independent. At least one budget must eventually be
  exceeded. Here the signed bias deviation \(\delta\) is
  independent in the rank test.
- For one simple pair and a single fixed \(d\), suppose the
  initial moments of \(\omega_d\) up to order \(s-1\) are zero
  exactly, or their absolute sum is \(o(|M_{d,s}|)\).
  If \(a+e/d\) and \(R_{d,s}\) are independent, the budget
  \(s\) is again excluded. This test imposes no condition on
  the earlier response vectors.

Indeed the signed normalized coefficient and center residuals
after subtracting the Gaussian and bias terms are *exactly*
\(\sum_d\int A_d(u)\,d\omega_d(u)\).
With \(X=W_+-W_--W_*\), the linear physical equation is
\(aX+e\delta+\sum_d\sum_{h\le s_d}R_{d,h}M_{d,h}\).
The same left-inverse argument now applies using (1).
In the scalar test the earlier terms are absorbed in the
error \(o(|M_{d,s}|)\). For exact matching the error from (1)
is \(O(\eta|M_{d,s}|)\).
In all cases the forced vanishing of the measure or measures
would leave zero net clock, contrary to \(X=0\) and \(W_*>0\).

If a finite menu is used in the middle test and there are
also differing rows outside it, write \(V_{\rm out}\) for
their exact compensated coefficient and center vector,
with their leading Gaussian and bias terms included in
\(X,\delta\). The same proof gives the necessary inequality
\[
      |X|+|\delta|+\sum_{d\in\mathcal D}\sum_{h\le s_d}|M_{d,h}|
                           \ \le\ C\|V_{\rm out}\|             \tag{5}
\]
whenever the menu budgets hold. Its right hand side is at
most a constant times
\(\eta\sum_{\text{outside rows}}d^2/m\), by the exact matched
factorial chord and center bounds. For a single-gap test with
outliers, the quantity \(|\delta-X/d|\) must also be added
to the right side, since its clock and bias are no longer
tied exactly. No finite-menu rank test by itself discards
older or newly occurring gap groups. The slot test (3)
allows them directly because its two measures use the exact
slot data from every row at the tested orders.

## Restrictions on a whole sequence

Before imposing independence, the same estimate gives a kernel
condition. Whenever the alternation budgets hold, collect the
clock, center or bias, and moment coordinates of the appropriate
chart in a vector \(y\). In a chart without outliers, \(y=0\) is
impossible at positive strength, by (1). At any continuing exact
contact, every limit of \(y/\|y\|\) belongs to the kernel of the
displayed linear response map. Thus any other constraints that
keep these normalized directions at a positive distance from
the kernel suffice for exclusion as well.

If the rank tests pass for every allocation of budgets with
\(s_n+s_N\le S\), the total number of sign changes eventually
exceeds \(S\). Indeed there are only finitely many allocations
under this bound. The same assertion applies to the fixed-menu
budgets with their sum in place of \(s_n+s_N\).
At any fixed pair the usable ranks are bounded by twice the
number of selected simple pairs. Unbounded growth therefore
does not follow from the fixed-pair rank tests. A sequence
with coalescing limiting roots is outside these charts.
