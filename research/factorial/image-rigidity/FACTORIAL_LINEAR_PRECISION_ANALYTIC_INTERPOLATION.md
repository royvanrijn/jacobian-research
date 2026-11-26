# Linear precision and analytic interpolation for arbitrary factorial supports

Fix a finite support of nonzero exponent columns \(\nu\) in \(d\)
variables, a prime \(p\), and a finite extension \(K/\mathbb Q_p\)
with ring of integers \(\mathcal O\), uniformizer \(\pi\), and \(v(p)=1\).
Write
\[
 f=\sum_\nu \frac{u_\nu}{\nu!}X^\nu,\qquad
 P_n=\frac{\Lambda_d(f^n)}{n!},\qquad E(T)=\sum_{n\ge0}P_nT^n.
\]
The natural coefficients range over a fixed integral residue class
\(u_\nu=u_{\nu,0}+\pi y_\nu,\ y_\nu\in\mathcal O\). Write
\(A=\mathcal O\langle(y_\nu)_\nu\rangle\) for the ring of
restricted power series with integral coefficients.

In each variable \(i\), retain the unit columns \(p^b e_i\)
and let \(h_i\) be the largest retained position. Put
\[
 U_i(T)=1-\sum_{p^b e_i\ {\rm retained}}(u_{p^b e_i}T)^{p^{h_i-b}},
 \qquad U=\prod_i U_i,
\]
taking \(U_i=1\) when the variable has no retained column.

**Theorem.** At each precision \(k\ge1\),
\[
                   E(T)\equiv N_k(T)/U(T)^{m_k}\pmod{p^kA[[T]]},
                                                               \tag{1}
\]
with \(N_k\in A[T]\), \(\deg N_k=O(k)\), and \(m_k=O(k)\).
The same conclusion holds with a constant term whose Gauss valuation
on the parameter polydisc is strictly above \(1/(p-1)\).

There are a fixed \(L\ge1\) and functions \(F_a(z,\mathbf y)\),
\(0\le a<L\), locally analytic in \(z\in\mathbb Z_p\) and analytic
in the coefficient polydisc, such that
\[
 P_n(\mathbf y)\equiv F_a((n-a)/L,\mathbf y)\pmod{p^k}
                 \quad(n\equiv a\bmod L,\ n\ge Ck+C_0).         \tag{2}
\]
The constants and finite-precision expressions are effective from
the support and residue class. The linear bounds are in the precision.
They can depend on the support, prime and ramification degree.

For a proposed \(w\)-window, put \(a+j=a_j+Lb_j\) with
\(0\le a_j<L\), for \(0\le j<w\). If the vectors
\[
                    \bigl(F_{a_j}(z+b_j,\mathbf y)\bigr)_{j=0}^{w-1}
                                                                \tag{3}
\]
are nonzero everywhere on the \(L\) domains
\(\mathbb Z_p\times\mathcal O^{\#\operatorname{supp}f}\), then
all sufficiently late \(w\)-windows contain a nonzero moment,
uniformly on the coefficient class. This hypothesis can be
certified at a finite precision. An unbounded sequence of exact
zero windows, allowing the coefficients to vary in this class,
forces a common zero of (3) on at least one domain.

## A fixed matrix with differential corrections

With zero constant term, use
\(W_N=N![X^N]\exp(Tf)\), \(W_0=1\), and \(W_N=0\)
for negative indices. For \(N\ne0\) let \(i\) be its first
positive coordinate. The
[partition recurrence](FACTORIAL_MOMENT_SERIES_PRIME_POWER_RATIONALITY.md#the-integral-partition-recurrence)
writes \(W_N\) as the sum of \(T W_{N-\nu}\) with weights
\[
               u_\nu\binom{N_i-1}{\nu_i-1}
                       \prod_{j\ne i}\binom{N_j}{\nu_j},
                 \qquad \nu_i>0.                              \tag{4}
\]

Choose one power \(K_0\) of \(p\) larger than every exponent
entry and such that \(v(K_0)>v_p(b!)\) for all \(b\) up
to the largest entry. This power is independent of precision.
Group vectors by residues \(r\bmod K_0\) and their set \(Z\)
of positive coordinates. For each state sum its \(W_N(T)\)
with weight \(\mathbf x^{(N-r)/K_0}\). Let \(S\) be the vector
of these series and \(e_0\) the vector of the zero exponent
state. At each \(T\)-degree every auxiliary series is a polynomial.

Adding \(\nu\) to a source vector \(r+K_0t\) determines the
target state from \(r,Z,\nu\). It also multiplies its auxiliary
weight by \(\prod_j x_j^{\lfloor(r_j+\nu_j)/K_0\rfloor}\).
The weight in (4) is an integral polynomial in \(t\) whose
nonconstant coefficients are divisible by \(p\). Indeed
\(\binom{c+K_0t}{b}\) has integral constant \(\binom cb\), and
every nonconstant numerator coefficient contains \(K_0\)
before division by \(b!\). Here each \(c\) is nonnegative,
including \(r_i+\nu_i-1\) in the distinguished coordinate.
The total degree in \(t\) is at most \(d_*-1\) where
\(d_*=\max_\nu|\nu|\).

Set \(\theta_j=x_j\partial_{x_j}\). The weight polynomial in
\(\theta_j\) acts on each coefficient of \(\mathbf x^t\)
by its value at \(t\). Multiplication by the carry monomial
after differentiation sends it to the correct target exponent.
The constant and nonconstant parts of the weights give the
integral matrix \(C(\mathbf x)\) and differential operator matrix
\(H(\mathbf x,\boldsymbol\theta)\) in
\[
                          S=e_0+TC S+pT H S .
\]
Both are independent of \(T\) and of precision. The first
positive-coordinate rule uses the target positivity set, which
is determined by \(Z\) and \(\nu\).

Let \(r_0\) be the number of states, \(G=(I-TC)^{-1}\), and
\(\Delta(T,\mathbf x)=\det(I-TC)\). This denominator has
constant coefficient one. The exact equation gives
\[
                   S\equiv\sum_{\ell=0}^{k-1}
                                  (pTG H)^\ell G e_0\pmod {p^k}.
                                                               \tag{5}
\]
Indeed the remainder \((pTG H)^k S\) has integral coefficients
divisible by \(p^k\). Each derivative in \(\mathbf x\) can
increase the power of \(\Delta\) in the denominator by one,
and increase numerator \(T\)-degree by at most \(r_0\) over
that next power. Multiplication by \(TG\) costs one denominator
power and at most \(r_0\) numerator degrees. Thus the \(\ell\)-th
term over \(\Delta^{1+d_*\ell}\) has numerator \(T\)-degree
at most \(r_0(1+d_*\ell)-1\). This bound also holds when
earlier terms are cleared against the last denominator.
Evaluating \(\mathbf x=\mathbf1\) and summing the states
gives an approximation to \(E\) with both degrees \(O(k)\)
and fixed denominator base \(\Delta_*(T)=\Delta(T,\mathbf1)\).

## Reducing the denominator to the singleton factors

The [fixed singleton pole theorem](FACTORIAL_MOMENT_SERIES_PRIME_POWER_RATIONALITY.md#which-poles-can-occur-at-higher-precision)
gives another expression at each precision with denominator a
power of \(U\). Its proof works over \(A\). It truncates the
nonsingleton counts by the native integer coefficient bound,
and the nonunit singleton counts by their positive Gauss
valuation. The resulting one-color matrices and polynomial
congruences apply directly over \(A\).

Reverse \(\Delta_*\) as the monic \(X^{r_0}\Delta_*(1/X)\).
Its reduction has constant coefficients because the natural
coefficients have fixed residues. Split that reduction into
three monic factors, according to whether a reciprocal root
is zero, is nonzero and outside the roots of the reversal
of \(\bar U\), or is among those roots. These factors are
pairwise coprime. Hensel factorization lifts them over the
complete ring \(A\): the residue Bézout identities solve
each successive coefficient correction, and the monic factors
converge coefficientwise. Reverse them to get
\(\Delta_*=Q_0Q_bQ_g\), each with constant coefficient one.
The reduction of \(Q_0\) is one, and \(Q_b,Q_g\) have
unit leading coefficients.

Set \(m=1+d_*(k-1)\) and \(e=1/v(\pi)\). In (5) after
specialization, remove \(Q_0^m\) from the denominator by
multiplying the numerator by the polynomial truncation
of \(Q_0^{-m}\) modulo \(p^k\). Since \(Q_0-1\) is divisible
by \(\pi\), only the first \(ek\) binomial terms contribute.
This increases the numerator degree by at most
\((ek-1)\deg Q_0\).

The factors \(Q_b\) and \(U\) have unit resultant by their
disjoint residue root sets and their unit leading coefficients.
Comparing \(N/(Q_b^mQ_g^m)\) to the second expression,
with denominator a power of \(U\), shows that \(N\) is
divisible by \(Q_b^m\) over \((A/p^k)[T]\): multiply across
denominators and apply the Bézout identity. Division by its
unit leading coefficient does not increase numerator degree.

Choose \(b_0\ge1\) with \(\bar Q_g\mid\bar U^{b_0}\).
Division over \(A\) gives \(U^{b_0}=Q_g B+\pi R\).
Thus
\[
                     Q_g^{-m}=\frac{B^m}{U^{b_0m}}
                                 (1-\pi R/U^{b_0})^{-m}.
\]
The binomial truncation at \(ek\) terms gives denominator
\(U^{b_0(m+ek-1)}\) and numerator degree \(O(k)\).
When \(Q_g=1\) there is no need for this step.
This proves (1) with computable linear bounds. A constant term
of valuation bounded below strictly above \(1/(p-1)\)
adds the factor \(\exp(cT)\), whose truncation modulo \(p^k\)
has degree \(O(k)\) by \(v_p(j!)\le j/(p-1)\).

## Compatible unit recurrences on the coefficient class

Let \(V(X)=X^{\deg U}U(1/X)\), a monic polynomial with
unit constant coefficient. Its reduction has constant
coefficients on the coefficient polydisc. If \(U=1\),
(1) already gives \(F_a=0\) and the required comparison.

Otherwise we may take \(m_k\ge1\), multiplying numerator
and denominator by \(U\) if needed. Equation (1) gives the
forward-shift recurrence \(V(\mathscr S)^{m_k}P=0\bmod p^k\)
beyond index \(O(k)\), by reading coefficients above
\(\deg N_k\) in \(U^{m_k}E\).
Its unit constant coefficient extends this tail uniquely
backwards to all nonnegative orders, giving \(u_k\) with
entries in \(A/p^k\). Adjacent precisions agree on a
sufficiently late tail modulo the lower precision. Their
difference satisfies the product of their two unit
recurrences, so agreement propagates backwards to order
zero. The compatible \(u_k\) define a limit \(u\) over
\(A\) agreeing with \(P_n\) modulo \(p^k\) for \(n\ge Ck+C_0\).

Choose \(L\) divisible by the orders of the roots of
\(\bar V\) in a finite residue extension. Write \(q=\deg V\)
and let \(J(Y)\) be the characteristic polynomial of
multiplication by \(X^L-1\) on \(A[X]/(V)\).
Its residue eigenvalues are all zero, also when
\(\bar V\) has repeated roots, so \(J\equiv Y^q\bmod\pi\).
Cayley--Hamilton gives \(V(X)\mid J(X^L-1)\); consequently
\(J(\mathscr S^L-1)^{m_k}\) annihilates \(u_k\).

The Gauss valuation of the coefficient of \(Y^i\) in
\(J\), for \(i<q\), is at least \(1/e\), hence at least
\((q-i)/(eq)\). Products give that same slope for the
coefficients of \(J^{m_k}\), measured from leading
degree \(q m_k\). Reduction of \(Y^r\) by \(J^{m_k}\)
then gives valuation at least \((r-i)/(eq)\) for each
remainder coefficient of degree \(i<qm_k\). This follows
by induction on \(r\), substituting the monic relation.
In particular
\[
                (\mathscr S^L-1)^r u(a)\equiv0\bmod p^k A
                                     \quad(r\ge qm_k+eqk).    \tag{6}
\]
This is a linear bound in the precision on the order of
the finite difference.

Define \(F_a(z,\mathbf y)=\sum_{r\ge0}
(\mathscr S^L-1)^r u(a)\binom zr\). At nonnegative integer
\(z\) it equals \(u(a+Lz)\) by binomial inversion.
On \(z=b+p^\ell t\) with integral \(b\), the numerator
factors of \(\binom zr\) have total Gauss valuation at least
\(\sum_{j=1}^{\ell}\lfloor r/p^j\rfloor\), by counting
multiples of \(p^j\) in the consecutive integers
\(b,\ldots,b-r+1\). After division by \(r!\) the loss is
at most \(r/(p^\ell(p-1))\).
Choose \(\ell\) so this slope is smaller than the positive
slope in (6). The series then converges in \(A\langle t\rangle\)
on these order disks. This proves joint analyticity and (2).
The [quadratic interpolation theorem](QUADRATIC_DYADIC_ANALYTIC_UNIT_ENVELOPE.md)
gives sharper support-specific bounds for degree two.
The use of binomial finite differences has the same
interpolation basis as Poonen's theorem for iterates
close to the identity; the recurrence construction and
the convergence bound here give a direct proof.
[Poonen, *p-adic interpolation of iterates*, Theorem 1](https://math.mit.edu/~poonen/papers/p-iteration.pdf)

Compactness of \(\mathbb Z_p\times\mathcal O^{\#\operatorname{supp}f}\)
gives the assertions about zero windows. If the vectors in (3)
are nonzero everywhere, their sup norms have a positive lower
bound. At precision \(k\) one can truncate the binomial series
using (6) and test on a finite product of residue sets.
Evaluation of elements of \(A\) modulo \(p^k\) is fixed by
the coefficient inputs modulo \(p^k\); each binomial term
requires at most a further \(v_p(r!)\) input digits in the
order parameter. Conversely, in an unbounded sequence of
exact zero windows one can choose a convergent subsequence
of both order and coefficient parameters in a fixed class.
Equation (2) and continuity make the limit a common zero
of (3). The argument leaves open the exclusion of these
common zeros for arbitrary supports.

## Spacing of exact zero windows at a fixed phase

Fix the coefficients, and assume that on each of the finitely many
order disks of analyticity at least one component of each vector
(3) is not identically zero. An analytic component on such a disk
has finitely many zeros there unless it vanishes identically.
Indeed at any zero its first nonzero Taylor coefficient makes
that zero isolated, and infinitely many zeros would accumulate
in the compact disk. Thus the common zero set in each class
\(a\) is finite.

A zero window beginning at a sufficiently large \(n\) in class \(a\)
makes every component of (3) have valuation at least
\(\lfloor(n-C_0)/C\rfloor\) at \(z=(n-a)/L\).
Outside neighborhoods of the finitely many common zeros the vectors
have a positive lower bound in sup norm. At a common zero
\(\zeta\), choose a component that is not identically zero on
its disk. On a smaller neighborhood it has the form
\((z-\zeta)^r g(z)\), \(r\ge1\), with \(v(g(z))\) constant,
by its first nonzero Taylor coefficient. Therefore every sufficiently
large zero window has an associated common zero \(\zeta\) with
\[
                         v((n-a)/L-\zeta)\ge c_1 n-c_0          \tag{8}
\]
for fixed \(c_1>0,c_0\). One can use a common pair of constants
over the finitely many classes and common zeros.

If two successive starts \(n<n'\) have the same \(a,\zeta\),
(8) gives \(v_p(n'-n)\ge c_1 n-c_0+v_p(L)\).
Since \(n'-n\) is a positive integer, these starts grow at least
exponentially in the earlier start. Thus the late zero starts
lie in finitely many sequences with exponential spacing along
each sequence. If all the relevant \(\zeta\) are algebraic over
\(\mathbb Q\), there are only finitely many exact zero windows.
In fact, for an integer \(z\) and a fixed algebraic \(\zeta\),
an integer polynomial \(g_0\) having root \(\zeta\) gives
\(v_p(z-\zeta)\le \log_p|g_0(z)|+O(1)=O(\log(z+1))\)
whenever \(g_0(z)\ne0\). The other factors of \(g_0(z)\)
have valuations bounded below, since \(z\) is integral.
The finitely many integer zeros of \(g_0\) can be omitted
in this estimate. This contradicts (8) for an unbounded
sequence. In particular an infinite zero-window sequence
under the nonidentity hypothesis requires a common zero
that is transcendental over \(\mathbb Q\) and exponentially
well approximated by positive integers at that prime.

## Exact replay

The [finite checker](verify_factorial_linear_precision_interpolation.py)
compares the partition transfer to 126 literal factorial coefficients
and checks 168 coefficients beyond the predicted numerator bound.
It uses SymPy for the exact characteristic polynomial and takes
the remaining operations over integers or residues. From the
repository root run:

    .venv/bin/python research/factorial/image-rigidity/verify_factorial_linear_precision_interpolation.py

These controls exercise the matrix transcription and the degree
bound. The differential and Hensel arguments above give the
all-precision and all-order conclusions.
