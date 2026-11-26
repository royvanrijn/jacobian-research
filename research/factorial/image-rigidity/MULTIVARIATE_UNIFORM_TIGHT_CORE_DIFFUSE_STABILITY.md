# Uniform diffuse angular stability for cores with tight exponent tails

For an ordinary two-monomial profile \(C\), ignore equal rows.
For each remaining row \(w=(a,b)\) put
\[
       r_w=\max(a,b),\qquad d_w=b-a,\qquad s_w=d_w^2/r_w .
\]
Call \(C\) an admissible core if it has at most \(R\) unequal
rows, \(\sum_w s_w\le K\), and for every positive integer \(t\)
it satisfies
\[
                       \sum_{r_w>B_t}s_w\le 2^{-t}.            \tag{1}
\]
Here \(R\) is fixed, \(K\) is a fixed positive rational cap, and
\(B_t\) is a computable nondecreasing sequence of integers tending
to infinity. The empty core is allowed. Membership of an individual
core in a smaller prescribed family need not be decidable, provided
each core in that family satisfies these bounds.

**Theorem.** Fix a positive real algebraic number \(\tau_0\).
There are computable \(\delta>0,\eta_0>0\), independent of the
admissible core, with the following property. Append any finite
set of unequal shared ordinary rows \(E\) to any admissible core.
Write \(\tau=\sum_{w\in E}s_w\) and
\(\eta=\max_{w\in E}|d_w|/r_w\). If
\[
                      |\tau-\tau_0|\le\delta,\qquad \eta\le\eta_0 ,
\]
then the consecutive moment polynomials on the extended support
have disjoint root arguments in every degree, including after
independent positive rescalings. The conclusion implies physical
SFC(2) over any characteristic-zero field.

In particular, it suffices to put a fixed bound \(L\) on the total
absolute gap \(\sum_w|d_w|\) of the core, with arbitrary ordinary
baselines and arbitrary shared or private orientations. For
\(L\ge1\), take \(R=K=L\) and \(B_t=\lceil 2^t L^2\rceil\).
Indeed \(s_w\le |d_w|\) and the omitted curvature is at most
\(L^2/(B_t+1)\). The case \(L=0\) is the empty core.
The appended rows have no bound on their number or total absolute
gap. The effective procedure uses the
[absolute ordinary angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md),
whose numerical value is unevaluated.

## Uniform gaps for the Gaussian references

Fix an integer \(J\ge4\) at least the absolute ordinary cutoff.
For each pair \(2\le n<N=n+1<J\), endpoint-normalize the
core polynomials at degrees \(j=n,N\) and multiply coefficient
\(k\) by \(e^{-\tau_0 k(j-k)/(2j)}\). Call these references
\(G_j(C,z)\). Their endpoints are one and their coefficients
are no larger than the binomial coefficients. The moduli of
their roots therefore lie between \(2^{-j}\) and \(2^j\).
For the upper bound, when \(|z|\ge1\), the root equation gives
\(|z|^j\le(2^j-1)|z|^{j-1}\); use the reciprocal polynomial
for the lower bound.

Use the [ordered quotient resultant](MULTIVARIATE_ARBITRARY_CORE_DIFFUSE_EXTENSION.md#an-angular-eliminant-for-the-exact-core)
of the two references, including the factor for the real ray
in the odd order. Denote its value by \(\mathcal R_n(C)\).
It is nonzero for every finite ordinary core. Indeed the
corresponding common-clock polynomial in
\(q=e^{-\tau_0/(2nN)}\) has algebraic coefficients and a nonzero
specialization at \(q=0\). The binomials at \(q=0\) have distinct
nontrivial root-of-unity quotient sets because \(n,N\) are
coprime. The real-ray factor enters only in the odd degree.
Hermite–Lindemann makes the actual \(q\) transcendental,
so the resultant cannot vanish. The same argument applies
to the empty core. The resultant is a continuous function of
the two normalized coefficient arrays: the root products giving
the quotient polynomials are symmetric, with denominator a power
of their fixed nonzero endpoint coefficient.
[Hermite–Lindemann and the algebraic-strength specialization](MULTIVARIATE_ARBITRARY_CORE_DIFFUSE_EXTENSION.md#an-angular-eliminant-for-the-exact-core)

For any sequence of cores with the row cap, fix a subsequence
on which its number of rows is constant, and label the rows.
After a further subsequence, each row either is constant or
has \(r_w\to\infty\). For a diverging row, (1) implies
\(s_w\to0\), because for every fixed \(t\) the row eventually
lies above \(B_t\). Its normalized coefficient contribution
tends to one. To see this, eventually \(s_w\le1/2\) and
\(|d_w|/r_w\le s_w\), so the row is shared and
\(\min(a,b)\ge r_w/2\). The second derivative of
\(\log\Gamma(1+x)\) is at most \(1/x\) for \(x>0\).
Its endpoint chord at index \(k\) consequently has log modulus
bounded by \(s_w k(j-k)/j\), which tends to zero.
Thus every coefficient limit is the reference on a finite
ordinary subcore, keeping the eventual constant rows. By
continuity and nonvanishing on these subcores, the quantities
\(|\mathcal R_n(C)|\) have a positive uniform lower bound.
The number of pairs below \(J\) is finite.

A positive resultant bound gives an explicit lower bound
for the distance between unit projections of roots of the
two references. Indeed, the quotient resultant is the product
of differences of the root quotients, including the odd factor.
Each quotient has modulus at most \(2^{2N}\). For any two
reference roots \(z,w\) the product includes
\(z/\bar z-w/\bar w\), using the dummy quotient \(1\) for an
odd real root. Its modulus is at most twice the distance between
\(z/|z|\) and \(w/|w|\). If the two quotient polynomials have
degrees \(M_n,M_N\), this gives a chord lower bound
\[
           \frac{|\mathcal R_n(C)|}
             {2(2^{2N+1})^{M_nM_N-1}} .                      \tag{2}
\]
The even reference has no real root at positive strength
by lognormal averaging of an ordinary nonnegative even
moment. Hence the real ray observation used in (2) holds
even when the core is empty.

## A terminating calculation without testing core membership

For a fixed \(t\), enumerate all finite candidate cores
with at most \(R\) rows, \(r_w\le B_t\), total curvature at most
\(K\), and satisfying the cutoff inequalities (1) at indices
\(1,\ldots,t\). All these tests are rational. This list
contains the truncations of every admissible core to \(r_w\le B_t\).
It may also contain profiles outside the prescribed family.

At each candidate, compute a positive lower bound on the
reference gap (2) by first isolating the algebraic coefficients,
then enclosing \(q\) with convergent rational exponential
bounds. Refine the enclosure of each nonzero
\(\mathcal R_n(C)\) until its lower bound in modulus is
at least half its upper bound. Take the minimum of the
computed bounds in (2), over all candidates and pairs, and
call it \(g_t>0\).

As \(t\to\infty\), these finite candidate lists still have
a common positive resultant margin for all sufficiently
large \(t\). Otherwise a sequence of near-zero candidates
admits the same constant-or-diverging row extraction as
before. For every fixed \(u\), all candidates with \(t\ge u\)
satisfy the required tail bound at \(B_u\), so every diverging
row again gives zero curvature in the limit. The limit core
is finite ordinary and its resultant is nonzero, a contradiction.
Thus \(g_t\) is bounded below by a positive number for all
large \(t\).

The omitted rows of an admissible core have total curvature
at most \(2^{-t}\). The chord estimate used above shows that
deleting them changes the log of each normalized coefficient
in degree \(j\) by at most \(j2^{-t}/4\). We use the following
root transfer bound. If two monic normalized polynomials of
degree \(j\) have coefficients bounded by \(\binom jk\)
and log coefficient differences at most \(\varepsilon\),
a root \(z\) of either has some root \(\zeta\) of the other
with
\[
     |z-\zeta|\le(1+2^j)\varepsilon^{1/j},\qquad
 \left|\frac z{|z|}-\frac\zeta{|\zeta|}\right|
           \le 2^{j+1}(1+2^j)\varepsilon^{1/j}.                \tag{3}
\]
Indeed each coefficient difference is at most
\(\binom jk\varepsilon\) by the exponential mean-value bound
on the normalized stripped coefficients. Thus the absolute
value of the other polynomial at \(z\) is at most
\((1+2^j)^j\varepsilon\). Factor it and use the annular
bounds for the unit projection inequality. No simplicity
assumption is needed.

Let \(A_J=\max_{2\le j<J}2^{j+1}(1+2^j)\). Increase \(t\)
until \(\varepsilon_t=J2^{-t}/4\le1\) and
\(2A_J\varepsilon_t^{1/(J-1)}<g_t/2\).
This search terminates by the uniform lower bound for the
candidate reference gaps. Each of the two roots on an
admissible core is then within the transfer bound of a
reference root on its truncated core. Hence \(g_t/2\) is
a certified common chord gap for the complete admissible
core references at every pair below \(J\).

## Transfer to the appended ordinary rows

Take \(\eta\le1/2\). For an appended shared row and \(j\ge2\),
its trigamma along the index segment is between
\(1/(j r_w+1)\) and \(1/[j r_w(1-\eta)]\).
Since \(1/r_w\le\eta\), these differ from \(1/(j r_w)\)
by at most \(2\eta/(j r_w)\). Comparing the chord with its
constant-curvature quadratic and summing the rows shows
that its total error relative to Gaussian strength \(\tau\)
has log modulus at most \(j\tau\eta/4\). Replacing
\(\tau\) by \(\tau_0\) contributes at most
\(j|\tau-\tau_0|/8\). This applies independently of the
core since endpoint normalization factors rowwise.
[Diffuse coefficient transfer](MULTIVARIATE_ARBITRARY_CORE_DIFFUSE_EXTENSION.md#uniform-transfer-to-the-added-ordinary-rows)

Choose positive rational \(\delta,\eta_0\) small enough that
\(\delta<\tau_0/2\), \(\eta_0\le1/2\), and the expression
\(\varepsilon=J((\tau_0+\delta)\eta_0/4+\delta/8)\) satisfies
\(\varepsilon\le1\) and \(2A_J\varepsilon^{1/(J-1)}<g_t/2\).
These choices are effective by rational enclosures.
Equation (3) now transfers the uniform no-common-ray
bound from each core reference to its ordinary extension.
The two ordinary low orders \(1,2\) separate by real
positivity, and the ordinary cutoff handles \(N\ge J\).
Bézout over \(\mathbb Q\) gives the field assertion.

The row cap is used to ensure that every coefficient limit
is a finite ordinary core, whose algebraic coefficients
permit the transcendence argument. This theorem does not
assert uniform stability over all countably retained cores.

